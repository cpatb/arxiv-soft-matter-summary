"""Batch analysis over every frame of every clip.

Design note: all geometric quantities are stored in *pixels* and converted to
body lengths afterwards, in analyse_results.py. Doing it this way avoids a
circular dependency -- the body-length estimate is itself a measurement drawn
from the same frames -- and means the normalisation can be revised without
re-running the expensive pass.

Filter scales still need some length to be set in, so a per-clip working scale
is estimated in a short pre-scan. Results are insensitive to it: it only sets
smoothing kernel sizes, not any reported number.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from catfish import calibrate, flow, frames, segment, topology  # noqa: E402

RESULTS = frames.ROOT / "results"

#: Isolated-silhouette measurement is expensive and body length changes only
#: with slow camera drift, so it is sampled rather than computed every frame.
BL_EVERY = 5

#: The radial packing profile is likewise sampled.
PROFILE_EVERY = 15

#: Frames used in the pre-scan that sets the working filter scale.
PRESCAN_FRAMES = 30


def prescan_body_px(clip: frames.Clip) -> float:
    """Rough per-clip body length in px, used only to set filter scales."""
    idx = np.linspace(0, len(clip) - 1, PRESCAN_FRAMES).astype(int)
    vals = []
    for i in idx:
        df = calibrate.isolated_silhouettes(clip.rgb(int(i)))
        if len(df):
            vals.extend(df["major_px"].tolist())
    if len(vals) < 20:
        return segment.NOMINAL_BODY_PX
    return float(np.median(vals))


def process_clip(clip: frames.Clip, limit: int | None = None) -> dict:
    t_start = time.perf_counter()
    work_px = prescan_body_px(clip)
    print(f"[{clip.name}] working scale {work_px:.1f} px "
          f"({time.perf_counter()-t_start:.0f}s prescan)", flush=True)

    n = len(clip) if limit is None else min(limit, len(clip))

    rows: list[dict] = []
    bl_samples: list[dict] = []
    corr_curves: list[dict] = []
    profiles: list[dict] = []
    events: list[dict] = []

    # Fission/fusion is detected incrementally against the previous frame only.
    # Buffering every frame's label image would need ~8 MB x n_frames (19 GB on
    # the longest clip), so only the previous one is retained.
    prev_labels: np.ndarray | None = None
    prev_gray = None
    for i in range(n):
        rgb = clip.rgb(i)
        # Derive grey from the already-loaded frame rather than re-reading and
        # re-converting it: L* is already computed inside the shared fields.
        fields = segment.compute_fields(rgb, body_px=work_px)
        gray = np.clip(fields.L * 255.0 / 100.0, 0, 255).astype(np.uint8)
        env = segment.school_envelope(fields=fields)

        row: dict = {"frame": i, "t_s": clip.time_s(i), "present": env.present}

        if env.present:
            # Shape metrics are computed in px here (body_px=1) and normalised
            # later; the dimensionless ones are already final.
            shapes = topology.component_shapes(env, body_px=1.0)
            shapes.sort(key=lambda s: s.area_bl2, reverse=True)
            big = shapes[0]
            row.update(
                {
                    "n_components": env.n_components,
                    "total_area_px2": sum(s.area_bl2 for s in shapes),
                    "area_px2": big.area_bl2,
                    "perimeter_px": big.perimeter_bl,
                    "circularity": big.circularity,
                    "solidity": big.solidity,
                    "aspect_ratio": big.aspect_ratio,
                    "orientation_rad": big.orientation_rad,
                    "gyration_px": big.gyration_bl,
                    "boundary_roughness": big.boundary_roughness,
                    "centroid_y": big.centroid_y,
                    "centroid_x": big.centroid_x,
                    "largest_fraction": big.area_bl2 / max(sum(s.area_bl2 for s in shapes), 1e-9),
                    "packing_fraction": segment.packing_fraction(None, env, fields=fields),
                }
            )
            if prev_labels is not None:
                events.extend(
                    topology.detect_events([prev_labels, env.labels], [i - 1, i])
                )
            prev_labels = env.labels.copy()

            if i % PROFILE_EVERY == 0:
                prof = topology.packing_profile(env, fields.fishness, body_px=1.0)
                if prof["depth_bl"]:
                    profiles.append({"frame": i, **prof})

        if i % BL_EVERY == 0 and env.present:
            df = calibrate.isolated_silhouettes(fields=fields, env=env)
            for mj in df["major_px"].to_numpy():
                bl_samples.append({"frame": i, "t_s": clip.time_s(i), "major_px": float(mj)})
            row["n_isolated"] = len(df)

        # --- kinematics: needs the previous frame -------------------------
        if prev_gray is not None and env.present:
            cam = flow.estimate_camera_motion(prev_gray, gray, exclude=env.mask)
            fl = flow.dense_flow(prev_gray, gray)
            s = flow.sample_school_velocity(fl, env.mask, work_px, camera=cam)

            row["camera_ok"] = cam.ok
            row["camera_inliers"] = cam.n_inliers
            if cam.ok:
                # Magnitude of the removed camera drift, for QA.
                cm = cam.velocity_at(np.array([gray.shape[0] / 2]),
                                     np.array([gray.shape[1] / 2]))
                row["camera_shift_px"] = float(np.hypot(cm[0, 0], cm[0, 1]))

            row.update(flow.order_parameters(s["u"], s["v"]))

            if len(s["u"]) >= 60:
                c = flow.correlation_function(
                    s["y"], s["x"], s["u"], s["v"], body_px=work_px
                )
                row["xi_bl_work"] = c["xi_bl"]
                row["n_flow_points"] = c["n_points"]
                if i % PROFILE_EVERY == 0 and len(c["r_bl"]):
                    corr_curves.append(
                        {
                            "frame": i,
                            "r_bl": c["r_bl"].tolist(),
                            "C": c["C"].tolist(),
                            "counts": c["counts"].tolist(),
                            "xi_bl": c["xi_bl"],
                        }
                    )

        rows.append(row)
        prev_gray = gray

        if i % 100 == 0:
            el = time.perf_counter() - t_start
            print(f"[{clip.name}] {i}/{n}  {el:.0f}s  "
                  f"({el/max(i,1):.2f} s/frame)", flush=True)

    RESULTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(RESULTS / f"{clip.name}_frames.csv", index=False)
    pd.DataFrame(bl_samples).to_csv(RESULTS / f"{clip.name}_bodylength.csv", index=False)

    with open(RESULTS / f"{clip.name}_events.json", "w") as fh:
        json.dump(events, fh, indent=2)
    with open(RESULTS / f"{clip.name}_correlation.json", "w") as fh:
        json.dump(corr_curves, fh, indent=2)
    with open(RESULTS / f"{clip.name}_profiles.json", "w") as fh:
        json.dump(profiles, fh, indent=2)

    meta = {
        "clip": clip.name,
        "n_frames": n,
        "fps": clip.fps,
        "working_body_px": work_px,
        "n_present": int(sum(r["present"] for r in rows)),
        "n_bl_samples": len(bl_samples),
        "n_events": len(events),
        "elapsed_s": time.perf_counter() - t_start,
    }
    with open(RESULTS / f"{clip.name}_meta.json", "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"[{clip.name}] done in {meta['elapsed_s']:.0f}s: {meta}", flush=True)
    return meta


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clips", nargs="*", default=list(frames.CLIPS))
    ap.add_argument("--limit", type=int, default=None,
                    help="process only the first N frames (for timing tests)")
    args = ap.parse_args()

    for name in args.clips:
        process_clip(frames.load_clip(name), limit=args.limit)


if __name__ == "__main__":
    main()
