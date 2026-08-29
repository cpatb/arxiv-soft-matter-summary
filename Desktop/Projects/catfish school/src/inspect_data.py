"""Inspect the raw frame data before choosing any analysis method.

Prints shape, dtype, value range, and channel statistics for sample frames from
each clip, and characterises the fish-vs-background contrast in several colour
spaces so the segmentation choice is made from data rather than assumption.

No units are available in these recordings (no ruler, no calibration target), so
everything here is reported in pixels; physical normalisation happens later via
the measured body length.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from skimage import color, io

ROOT = Path(__file__).resolve().parent.parent
FRAMES = ROOT / "data" / "frames"
CLIPS = ["IMG_0544", "IMG_0545", "IMG_0546", "IMG_0547"]


def frame_paths(clip: str) -> list[Path]:
    return sorted((FRAMES / clip).glob("*.jpg"))


def describe(clip: str) -> dict:
    paths = frame_paths(clip)
    n = len(paths)
    # Sample at 15%, 45%, 75% of the clip -- same instants as the visual previews.
    idx = [int(n * f) for f in (0.15, 0.45, 0.75)]
    out = {"clip": clip, "n_frames": n, "samples": []}

    for i in idx:
        img = io.imread(paths[i])
        assert img.ndim == 3 and img.shape[2] == 3, f"expected RGB, got {img.shape}"
        assert img.dtype == np.uint8, f"expected uint8, got {img.dtype}"

        # skimage reads RGB; float conversion keeps us in [0, 1] for all downstream work.
        rgb = img.astype(np.float64) / 255.0
        lab = color.rgb2lab(rgb)          # L in [0, 100]
        gray = color.rgb2gray(rgb)        # luminance in [0, 1]

        rec = {
            "frame_index": i,
            "shape": list(img.shape),
            "dtype": str(img.dtype),
            "rgb_min": [float(rgb[..., c].min()) for c in range(3)],
            "rgb_max": [float(rgb[..., c].max()) for c in range(3)],
            "rgb_mean": [float(rgb[..., c].mean()) for c in range(3)],
            "gray_percentiles": {
                str(p): float(np.percentile(gray, p))
                for p in (0.1, 1, 5, 25, 50, 75, 95, 99, 99.9)
            },
            "L_percentiles": {
                str(p): float(np.percentile(lab[..., 0], p))
                for p in (0.1, 1, 5, 25, 50, 75, 95, 99, 99.9)
            },
            "n_nan": int(np.isnan(rgb).sum()),
            "n_inf": int(np.isinf(rgb).sum()),
        }

        # The fish are dark. Estimate how separable "dark" is by looking at the
        # gap between the low tail and the bulk of the luminance histogram.
        lo, med = rec["gray_percentiles"]["5"], rec["gray_percentiles"]["50"]
        rec["dark_contrast_ratio"] = float(med / lo) if lo > 0 else float("inf")
        out["samples"].append(rec)

    return out


def main() -> None:
    results = [describe(c) for c in CLIPS]
    for r in results:
        print(f"\n=== {r['clip']}  ({r['n_frames']} frames) ===")
        for s in r["samples"]:
            g = s["gray_percentiles"]
            print(
                f"  frame {s['frame_index']:5d}  shape={s['shape']}  {s['dtype']}  "
                f"NaN={s['n_nan']} Inf={s['n_inf']}"
            )
            print(
                f"    gray pct  0.1%={g['0.1']:.3f}  1%={g['1']:.3f}  5%={g['5']:.3f}  "
                f"50%={g['50']:.3f}  95%={g['95']:.3f}  99.9%={g['99.9']:.3f}"
            )
            print(
                f"    RGB mean={[round(v,3) for v in s['rgb_mean']]}  "
                f"dark contrast (med/5%)={s['dark_contrast_ratio']:.2f}"
            )

    (ROOT / "data" / "meta").mkdir(parents=True, exist_ok=True)
    with open(ROOT / "data" / "meta" / "data_inspection.json", "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nwrote {ROOT / 'data' / 'meta' / 'data_inspection.json'}")


if __name__ == "__main__":
    main()
