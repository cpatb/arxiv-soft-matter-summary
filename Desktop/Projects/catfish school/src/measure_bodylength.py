"""Dedicated body-length pass, with its systematic uncertainty measured.

Run separately from the main pipeline for two reasons.

First, the batch pass collected silhouettes using an acceptance window scaled to
the *assumed* body length, which made the estimate circular and self-reinforcing
(see calibrate.MIN_ABS_AREA_PX for what went wrong). The criteria are now
scale-free and the measurement has to be redone.

Second, body length is a property of the population and changes only through slow
camera drift, so it does not need every fifth frame. A sample of ~120 frames per
clip is ample and costs minutes rather than hours.

The output includes a threshold sweep. The measured length is the major axis of
the *thresholded dark silhouette*, and a stricter threshold clips more of the
faint tail, so the value depends on where the threshold is set. Sweeping it
converts that from a hidden bias into a stated uncertainty.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import measure, morphology

sys.path.insert(0, str(Path(__file__).resolve().parent))

from catfish import calibrate, frames, segment  # noqa: E402

RESULTS = frames.ROOT / "results"

#: Frames sampled per clip.
N_FRAMES = 120

#: Background percentile defining the fish threshold. The reported value.
THRESHOLD_PCT = 99.5

#: Thresholds swept to quantify the systematic.
SWEEP_PCT = (98.0, 99.0, 99.5, 99.8, 99.9)

SEED = 20260808


def silhouettes_at(fields, env, pct: float) -> list[dict]:
    """Accepted isolated silhouettes at a given background threshold."""
    f = fields.fishness
    outside = ~env.mask
    if outside.sum() < 10_000:
        return []
    thr = float(np.percentile(f[outside], pct))
    grown = ndi.binary_dilation(env.mask, iterations=int(round(fields.body_px)))
    fg = (f > thr) & grown
    fg = morphology.remove_small_objects(fg, max_size=int(calibrate.MIN_ABS_AREA_PX) - 1)
    if not fg.any():
        return []
    labels, _ = ndi.label(fg)

    out = []
    for p in measure.regionprops(labels):
        if not (calibrate.MIN_ABS_AREA_PX <= p.area <= calibrate.MAX_ABS_AREA_PX):
            continue
        if p.axis_minor_length <= 0:
            continue
        aspect = p.axis_major_length / p.axis_minor_length
        if not (calibrate.ASPECT_RANGE[0] <= aspect <= calibrate.ASPECT_RANGE[1]):
            continue
        if p.solidity < calibrate.ISOLATED_MIN_SOLIDITY:
            continue
        if p.eccentricity < calibrate.ISOLATED_MIN_ECCENTRICITY:
            continue

        clearance = max(2, int(round(calibrate.ISOLATION_CLEARANCE_BL * p.axis_major_length)))
        r0, c0, r1, c1 = p.bbox
        pad = clearance + 2
        R0, C0 = max(0, r0 - pad), max(0, c0 - pad)
        R1, C1 = min(labels.shape[0], r1 + pad), min(labels.shape[1], c1 + pad)
        sub = labels[R0:R1, C0:C1]
        me = sub == p.label
        if (ndi.binary_dilation(me, iterations=clearance) & ((sub > 0) & (sub != p.label))).any():
            continue

        out.append(
            {
                "y": float(p.centroid[0]), "x": float(p.centroid[1]),
                "area_px": float(p.area),
                "major_px": float(p.axis_major_length),
                "minor_px": float(p.axis_minor_length),
                "aspect": float(aspect),
                "solidity": float(p.solidity),
                "eccentricity": float(p.eccentricity),
            }
        )
    return out


def main() -> None:
    rng = np.random.default_rng(SEED)
    summary = {}

    for name in frames.CLIPS:
        mp, fp = RESULTS / f"{name}_meta.json", RESULTS / f"{name}_frames.csv"
        if not (mp.exists() and fp.exists()):
            print(f"skip {name}: pipeline output missing")
            continue
        with open(mp) as fh:
            work = json.load(fh)["working_body_px"]
        clip = frames.load_clip(name)
        fr = pd.read_csv(fp)
        present = fr.loc[fr["present"], "frame"].to_numpy()
        if len(present) < 5:
            continue
        picks = sorted(rng.choice(present, size=min(N_FRAMES, len(present)),
                                  replace=False).tolist())

        rows = []
        sweep = {p: [] for p in SWEEP_PCT}
        for i in picks:
            i = int(i)
            fields = segment.compute_fields(clip.rgb(i), body_px=work)
            env = segment.school_envelope(fields=fields)
            if not env.present:
                continue
            for p in SWEEP_PCT:
                sil = silhouettes_at(fields, env, p)
                sweep[p].extend(s["major_px"] for s in sil)
                if p == THRESHOLD_PCT:
                    for s in sil:
                        rows.append({"frame": i, "t_s": clip.time_s(i), **s})

        df = pd.DataFrame(rows)
        df.to_csv(RESULTS / f"{name}_bodylength_v2.csv", index=False)

        if df.empty:
            print(f"{name}: no silhouettes accepted")
            continue

        med = float(df["major_px"].median())
        sweep_med = {str(p): float(np.median(v)) if len(v) else np.nan
                     for p, v in sweep.items()}
        finite = [v for v in sweep_med.values() if np.isfinite(v)]
        summary[name] = {
            "n_frames_sampled": len(picks),
            "n_silhouettes": int(len(df)),
            "body_px": med,
            "iqr_px": [float(df["major_px"].quantile(.25)),
                       float(df["major_px"].quantile(.75))],
            "aspect_median": float(df["aspect"].median()),
            "threshold_sweep_px": sweep_med,
            "systematic_range_px": [float(min(finite)), float(max(finite))] if finite else None,
            "systematic_pct": float(100 * (max(finite) - min(finite)) / med)
            if finite and med > 0 else np.nan,
        }
        print(f"{name}: BL = {med:.1f} px  (n={len(df)}, "
              f"threshold sweep {min(finite):.1f}-{max(finite):.1f} px)", flush=True)

    with open(RESULTS / "bodylength_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    print("\nwrote", RESULTS / "bodylength_summary.json")


if __name__ == "__main__":
    main()
