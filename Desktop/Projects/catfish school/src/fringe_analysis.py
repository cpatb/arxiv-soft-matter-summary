"""Individual-level analysis at the school fringe, and a check on the flow speeds.

Two purposes.

1. **Cross-validation.** All the kinematics in the main pipeline come from dense
   optical flow, a single method. Tracking individual fish at the fringe gives a
   second, methodologically independent estimate of swimming speed. If the two
   agree the flow field can be believed; if they disagree, one of them is wrong
   and saying so is more useful than picking a favourite.

2. **Nearest-neighbour spacing**, reported with its bias stated. Detection is
   incomplete -- verified at roughly half of individuals in moderately dense
   regions -- and missed neighbours can only push a nearest-neighbour distance
   *upward*. The values here are therefore an upper bound on true spacing, and
   the packing fraction from the main pipeline (measured from area, needing no
   individuation) is the trustworthy density statistic.

Two scales, and they are not interchangeable
--------------------------------------------
``working_body_px`` is the scale the segmentation filters were built at in the
batch pass. It is a *filter* parameter and it is 1.4 to 1.6x larger than the fry
actually are; passing anything else to ``segment`` would change what is detected.

``body_length.median_px`` is the *measured* fry length, and it is the unit every
number in this project is divided by. An earlier version of this module fell back
to the working scale whenever the measured value was missing from the meta file
-- which it always was -- so every speed and every spacing here came out 1.4-1.6x
too small while every ratio still looked right. There is no fallback now: a
missing measurement raises.

Camera motion
-------------
The footage is handheld and the pan is not small: the median frame-to-frame
camera shift is 1.0 to 2.8 px against a fish step of order 2 px. Detections are
therefore mapped back into the coordinate frame of each window's first frame,
using the same ``flow.estimate_camera_motion`` similarity transform the rest of
the pipeline uses, before anything is linked or differenced. Nearest-neighbour
distances are measured within a single frame and are untouched by this.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import trackpy as tp
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))

from catfish import flow, frames, plotstyle, segment  # noqa: E402

RESULTS = frames.ROOT / "results"
FIGS = frames.ROOT / "figures"

#: Length of each tracking window, in frames. Short windows are deliberate:
#: fringe fish leave the fringe, and a long window would keep only the atypical
#: individuals that happen to stay put.
WINDOW_FRAMES = 45

#: Number of windows sampled per clip.
N_WINDOWS = 6

#: Maximum frame-to-frame displacement for linking, as a fraction of the
#: *working* scale -- the same scale the detector runs at, so the radius stays
#: matched to the detections it links. 0.55 x working is about 0.85 measured body
#: lengths, which is generous headroom over the measured swimming speeds.
LINK_RANGE_WORK = 0.55

#: Frames a track may vanish for and still be continued.
LINK_MEMORY = 2

#: Tracks shorter than this are dropped before any statistic is computed.
MIN_TRACK_FRAMES = 12

tp.quiet()


def measured_body_px(meta: dict, name: str) -> float:
    """The measured fry length in px: the unit every result here is divided by.

    Deliberately has no fallback. Substituting the working scale when this is
    absent is exactly the failure this function exists to prevent -- it produced
    a whole results file that was internally consistent, plausible, and out by a
    factor of 1.4 to 1.6. Written by ``analyse_results.py`` from the isolated
    silhouettes measured in ``measure_bodylength.py``.
    """
    bl = meta.get("body_length", {}).get("median_px")
    if bl is None:
        raise KeyError(
            f"{name}_meta.json has no body_length.median_px. Run "
            f"src/analyse_results.py first -- it measures the body length and "
            f"writes it back into the meta file. Refusing to substitute "
            f"working_body_px ({meta.get('working_body_px')!r} px), which is a "
            f"segmentation filter scale and not a fish."
        )
    bl = float(bl)
    if not np.isfinite(bl) or bl <= 0:
        raise ValueError(f"{name}: body_length.median_px is {bl!r}, not a usable length")
    return bl


def frame_state(clip: frames.Clip, index: int, work_px: float):
    """Fringe detections, greyscale, and school mask for one frame.

    Returned together because all three come from the same ``compute_fields``
    call, which is the expensive part.
    """
    fields = segment.compute_fields(clip.rgb(index), body_px=work_px)
    gray = np.clip(fields.L * 255 / 100, 0, 255).astype(np.uint8)
    env = segment.school_envelope(fields=fields)
    empty = pd.DataFrame(columns=["y", "x", "frame"])
    if not env.present:
        return empty, gray, None
    det = segment.detect_individuals(fields=fields, envelope=env)
    if det.empty:
        return empty, gray, env.mask
    det = det[~det["in_core"]].copy()
    det["frame"] = index
    return det[["y", "x", "frame"]], gray, env.mask


def fringe_detections(clip: frames.Clip, index: int, work_px: float) -> pd.DataFrame:
    """Detections outside the dense core, where they can actually be believed."""
    return frame_state(clip, index, work_px)[0]


def track_window(clip: frames.Clip, start: int, work_px: float) -> tuple[pd.DataFrame, pd.DataFrame, bool]:
    """Link fringe detections across a short window, in a camera-fixed frame.

    Returns (tracks, first-frame detections in raw image coordinates, truncated).

    Every detection is mapped into the coordinate frame of ``start`` by the
    accumulated inverse of the frame-to-frame camera transform, so what the
    linker sees -- and what is later differenced into a speed -- is fish motion
    with the photographer removed. If the camera estimate fails on some pair the
    window is truncated there rather than continued with an assumed pan: the
    alternative is silently mixing corrected and uncorrected displacements.
    """
    feats: list[pd.DataFrame] = []
    det_first = pd.DataFrame(columns=["y", "x", "frame"])
    prev_gray, prev_mask = None, None
    to_start = np.eye(3)  # maps this frame's coordinates into `start`'s
    truncated = False

    for i in range(start, min(start + WINDOW_FRAMES, len(clip))):
        det, gray, mask = frame_state(clip, i, work_px)
        if i == start:
            det_first = det
        else:
            cam = flow.estimate_camera_motion(prev_gray, gray, exclude=prev_mask)
            if not cam.ok:
                truncated = True
                break
            M = np.vstack([cam.matrix, [0.0, 0.0, 1.0]])  # previous -> this frame
            to_start = to_start @ np.linalg.inv(M)
        prev_gray, prev_mask = gray, mask

        if len(det):
            xy = det[["x", "y"]].to_numpy(dtype=float)
            fixed = xy @ to_start[:2, :2].T + to_start[:2, 2]
            d = det.copy()
            d["x"], d["y"] = fixed[:, 0], fixed[:, 1]
            feats.append(d)

    if not feats:
        return pd.DataFrame(), det_first, truncated
    f = pd.concat(feats, ignore_index=True)
    if f["frame"].nunique() < MIN_TRACK_FRAMES:
        return pd.DataFrame(), det_first, truncated

    search = max(2.0, LINK_RANGE_WORK * work_px)
    linked = tp.link(f, search_range=search, memory=LINK_MEMORY,
                     pos_columns=["y", "x"], t_column="frame")
    tracks = tp.filter_stubs(linked, MIN_TRACK_FRAMES).reset_index(drop=True)
    return tracks, det_first, truncated


def track_speeds(tracks: pd.DataFrame, body_px: float, fps: float) -> np.ndarray:
    """Per-track median speed in body lengths per second.

    ``body_px`` is the *measured* fry length. Positions are already in a
    camera-fixed frame, so the differences below are fish motion.

    The per-track median, rather than the mean of all steps, keeps a single bad
    link from dominating: a mislink produces one huge displacement, which the
    median ignores and the mean does not.
    """
    out = []
    for _, g in tracks.groupby("particle"):
        g = g.sort_values("frame")
        dy = np.diff(g["y"].to_numpy())
        dx = np.diff(g["x"].to_numpy())
        dt = np.diff(g["frame"].to_numpy())
        ok = dt > 0
        if ok.sum() < 5:
            continue
        step = np.hypot(dy[ok], dx[ok]) / dt[ok]  # px per frame
        out.append(np.median(step) / body_px * fps)
    return np.array(out)


def nnd_bl(det: pd.DataFrame, body_px: float) -> np.ndarray:
    """Nearest-neighbour distances among fringe detections, in body lengths.

    Measured within one frame, so camera motion cannot enter.
    """
    if len(det) < 5:
        return np.array([])
    pts = det[["y", "x"]].to_numpy()
    d, _ = cKDTree(pts).query(pts, k=2)
    return d[:, 1] / body_px


def main() -> None:
    plotstyle.apply()
    FIGS.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20260808)  # recorded seed; window choice is random

    out: dict = {}
    speed_by_clip: dict[str, np.ndarray] = {}
    nnd_by_clip: dict[str, np.ndarray] = {}

    for name in frames.CLIPS:
        meta_p = RESULTS / f"{name}_meta.json"
        frames_p = RESULTS / f"{name}_frames.csv"
        if not (meta_p.exists() and frames_p.exists()):
            print(f"skip {name}: pipeline output missing")
            continue
        with open(meta_p) as fh:
            meta = json.load(fh)
        work_px = meta["working_body_px"]      # segmentation scale
        bl_px = measured_body_px(meta, name)   # the unit; no fallback
        clip = frames.load_clip(name)
        print(f"{name}: working scale {work_px:.2f} px, measured BL {bl_px:.2f} px "
              f"(ratio {work_px / bl_px:.3f})", flush=True)

        fr = pd.read_csv(frames_p)
        present = fr.loc[fr["present"], "frame"].to_numpy()
        usable = present[present < len(clip) - WINDOW_FRAMES]
        if len(usable) < WINDOW_FRAMES:
            print(f"skip {name}: too few frames with a school present")
            continue

        starts = rng.choice(usable, size=min(N_WINDOWS, len(usable)), replace=False)

        speeds, nnds, n_tracks, n_det, n_trunc = [], [], 0, [], 0
        for s in sorted(starts):
            tr, det0, truncated = track_window(clip, int(s), work_px)
            n_trunc += int(truncated)
            if len(tr):
                n_tracks += tr["particle"].nunique()
                speeds.append(track_speeds(tr, bl_px, clip.fps))
            n_det.append(len(det0))
            v = nnd_bl(det0, bl_px)
            if len(v):
                nnds.append(v)

        sp = np.concatenate(speeds) if speeds else np.array([])
        nd = np.concatenate(nnds) if nnds else np.array([])
        speed_by_clip[name] = sp
        nnd_by_clip[name] = nd

        # The flow-based speed over the same clip, for comparison. Divided by the
        # same measured body length, so the ratio below compares like with like.
        flow_speed = np.nan
        if "median_speed_px" in fr:
            d = fr[fr["present"]]
            flow_speed = float(np.nanmedian(d["median_speed_px"]) / bl_px * clip.fps)

        out[name] = {
            "body_px_measured": bl_px,
            "working_body_px": float(work_px),
            "n_windows": int(len(starts)),
            "n_windows_truncated_camera": int(n_trunc),
            "n_tracks": int(n_tracks),
            "n_fringe_detections_median": float(np.median(n_det)) if n_det else 0.0,
            "track_speed_bl_s_median": float(np.median(sp)) if len(sp) else np.nan,
            "track_speed_bl_s_iqr": [float(np.percentile(sp, 25)),
                                     float(np.percentile(sp, 75))] if len(sp) else None,
            "flow_speed_bl_s_median": flow_speed,
            "speed_ratio_track_over_flow": float(np.median(sp) / flow_speed)
            if len(sp) and np.isfinite(flow_speed) and flow_speed > 0 else np.nan,
            "nnd_bl_median": float(np.median(nd)) if len(nd) else np.nan,
            "nnd_bl_iqr": [float(np.percentile(nd, 25)),
                           float(np.percentile(nd, 75))] if len(nd) else None,
        }
        print(f"{name}: {out[name]}", flush=True)

    # --- figure ----------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(8.4, 2.9))
    for name in speed_by_clip:
        c = plotstyle.CLIP_COLOUR[name]
        sp, nd = speed_by_clip[name], nnd_by_clip[name]
        if len(sp) > 5:
            axes[0].hist(sp, bins=25, histtype="step", density=True, color=c,
                         label=name.replace("IMG_", "clip "))
        if len(nd) > 5:
            axes[1].hist(nd, bins=30, histtype="step", density=True, color=c)
    axes[0].set_xlabel(r"track speed  (BL s$^{-1}$)")
    axes[0].set_ylabel("density")
    axes[0].legend(loc="upper right")
    axes[1].set_xlabel(r"fringe NND  (BL)")
    axes[1].set_ylabel("density")

    # Agreement between the two independent speed estimates.
    xs = [v["flow_speed_bl_s_median"] for v in out.values()]
    ys = [v["track_speed_bl_s_median"] for v in out.values()]
    cs = [plotstyle.CLIP_COLOUR[k] for k in out]
    axes[2].scatter(xs, ys, c=cs, s=28)
    lim = [0, max([v for v in xs + ys if np.isfinite(v)] + [0.1]) * 1.25]
    axes[2].plot(lim, lim, color="0.5", ls="--", lw=1.0)
    axes[2].set_xlim(lim)
    axes[2].set_ylim(lim)
    axes[2].set_xlabel(r"optical flow  (BL s$^{-1}$)")
    axes[2].set_ylabel(r"individual tracks  (BL s$^{-1}$)")
    for ax, L in zip(axes, "abc"):
        plotstyle.panel_label(ax, L)
    fig.tight_layout()
    fig.savefig(FIGS / "fig5_fringe_individuals.png")
    plt.close(fig)

    with open(RESULTS / "fringe_analysis.json", "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    print("\nwrote", RESULTS / "fringe_analysis.json")


if __name__ == "__main__":
    main()
