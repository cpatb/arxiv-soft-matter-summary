"""
Wrinkle analysis: per-segment measurements and batch processing.

Usage:
    python analyze.py <root_folder> [px_per_um]
"""

import re
import sys
import numpy as np
import pandas as pd
import tifffile
from pathlib import Path
from skimage.morphology import skeletonize
from skimage.measure import regionprops
from scipy.ndimage import label, convolve


_EIGHT_CONN = np.ones((3, 3), dtype=int)
MIN_SKEL_PX = 5
_NEIGHBOR_KERN = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]], dtype=np.uint8)


def count_forks(skel: np.ndarray) -> int:
    """Branch points: skeleton pixels with 3+ skeleton neighbors (8-connected)."""
    neighbors = convolve(skel.astype(np.uint8), _NEIGHBOR_KERN, mode="constant")
    return int(np.sum((neighbors >= 3) & (skel > 0)))


def analyze_mask(
    binary_mask: np.ndarray, px_per_um: float | None = None
) -> pd.DataFrame:
    """Return one row per connected wrinkle segment measured on the skeleton."""
    skel = skeletonize(binary_mask)
    labeled, _ = label(skel, structure=_EIGHT_CONN)  # type: ignore[misc]
    rows = []
    for p in regionprops(labeled):
        if p.area < MIN_SKEL_PX:
            continue
        cy, cx = p.centroid
        rows.append({
            "label": p.label,
            "length_px": p.area,
            "length_um": p.area / px_per_um if px_per_um is not None else None,
            "orientation_deg": (np.degrees(p.orientation) + 180) % 180,
            "centroid_x": cx,
            "centroid_y": cy,
            "area_px": p.area,
        })
    return pd.DataFrame(rows)


def summarize_image(df: pd.DataFrame, image_shape: tuple) -> dict:
    """Aggregate statistics for a single image."""
    total_px = int(df["length_px"].sum()) if len(df) > 0 else 0
    image_area = image_shape[0] * image_shape[1]
    return {
        "n_wrinkles": len(df),
        "total_length_px": total_px,
        "mean_length_px": float(df["length_px"].mean()) if len(df) > 0 else 0.0,
        "std_length_px": float(df["length_px"].std()) if len(df) > 0 else 0.0,
        "mean_orientation_deg": float(df["orientation_deg"].mean()) if len(df) > 0 else None,
        "areal_density": total_px / image_area,
    }


def run_batch(
    root_folder: str | Path, px_per_um: float | None = None
) -> pd.DataFrame:
    """Walk root_folder for *_roi.tif files and return a concatenated DataFrame."""
    root = Path(root_folder)
    all_dfs: list[pd.DataFrame] = []

    for roi_path in sorted(root.rglob("*_roi.tif")):
        stem = roi_path.stem[:-4]  # strip '_roi'
        m = re.search(r"(\d+)min", stem)
        if m is None:
            print(f"Skipping {roi_path.name}: no timepoint found")
            continue
        timepoint_min = int(m.group(1))

        # Expected layout: <root>/<condition>/<sample_N>/<face>/<files>
        face_dir = roi_path.parent
        sample_dir = face_dir.parent
        condition_dir = sample_dir.parent

        face = face_dir.name
        sm = re.search(r"(\d+)", sample_dir.name)
        sample = int(sm.group(1)) if sm else None
        condition = condition_dir.name

        roi_map = tifffile.imread(roi_path)
        h, w = roi_map.shape[-2], roi_map.shape[-1]
        skel_path = roi_path.parent / f"{stem}_skeleton.tif"
        skel = (tifffile.imread(skel_path) > 0) if skel_path.exists() else skeletonize(roi_map > 0)
        n_forks = count_forks(skel)

        df = analyze_mask(roi_map > 0, px_per_um=px_per_um)
        df["condition"] = condition
        df["sample"] = sample
        df["face"] = face
        df["timepoint_min"] = timepoint_min
        df["image_h"] = h
        df["image_w"] = w
        df["n_forks"] = n_forks
        df["roi_file"] = str(roi_path.relative_to(root))

        all_dfs.append(df)
        print(f"  {roi_path.relative_to(root)}: {len(df)} segments, {n_forks} forks")

    if not all_dfs:
        print("No _roi.tif files found.")
        return pd.DataFrame()

    result = pd.concat(all_dfs, ignore_index=True)
    out = root / "results.csv"
    result.to_csv(out, index=False)
    print(f"Saved {out}  ({len(result)} rows, {len(all_dfs)} images)")
    return result


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: python analyze.py <root_folder> [px_per_um]")
    px = float(sys.argv[2]) if len(sys.argv) > 2 else None
    run_batch(sys.argv[1], px)
