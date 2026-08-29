"""
Interactive ROI alignment using napari.

Usage:
    python align.py <face_folder>

Opens all *_roi.tif files as separate napari layers (color-coded by timepoint).
Move each layer using either:
  - the "translate" boxes in the layer-controls panel, or
  - Transform mode (the square-with-arrows icon) + drag.
Close the viewer — offsets are written to <face_folder>/offsets.json and read
automatically by visualize.py.
"""

import sys
import json
import re
import numpy as np
import tifffile
import napari
from napari.layers import Image
from pathlib import Path
from typing import cast

COLORMAPS = ["red", "green", "cyan", "yellow", "magenta", "bop blue", "bop orange"]


def main(face_folder: str) -> None:
    folder    = Path(face_folder)
    roi_paths = sorted(p for p in folder.glob("*_roi.tif"))
    if not roi_paths:
        sys.exit(f"No _roi.tif files found in {folder}")

    offsets_file = folder / "offsets.json"
    existing_offsets = {}
    if offsets_file.exists():
        with open(offsets_file) as f:
            existing_offsets = json.load(f)
        print(f"Loaded existing offsets from {offsets_file}")
        for fname, off in existing_offsets.items():
            print(f"  {fname}: dy={off['dy']:.1f}, dx={off['dx']:.1f}")

    viewer = napari.Viewer(title="ROI Alignment — shift layers, then close to save offsets.json")

    for i, roi_path in enumerate(roi_paths):
        stem = roi_path.stem[:-4]
        m    = re.search(r"(\d+)min", stem)
        name = m.group(0) if m else roi_path.stem
        img  = (tifffile.imread(roi_path) > 0).astype(np.float32)
        layer = cast(Image, viewer.add_image(
            img,
            name=name,
            colormap=COLORMAPS[i % len(COLORMAPS)],
            blending="additive",
            opacity=0.6,
            metadata={"roi_filename": roi_path.name},
        ))

        if roi_path.name in existing_offsets:
            off = existing_offsets[roi_path.name]
            translate = [0.0] * layer.ndim
            translate[-2] = float(off["dy"])
            translate[-1] = float(off["dx"])
            layer.translate = translate

    print("Tip: select a layer, then either edit the translate boxes in the layer")
    print("     controls panel, or click the Transform icon (square with arrows)")
    print("     and drag.  Close the viewer when done.")
    napari.run()

    # Preserve existing entries; only overwrite layers that are still present.
    offsets = dict(existing_offsets)
    for layer in viewer.layers:
        fname = layer.metadata.get("roi_filename") if hasattr(layer, "metadata") else None
        if not fname:
            continue
        # data_to_world captures combined translate + affine + scale offset,
        # so it works whether the user used the translate boxes or transform-mode drag.
        origin = np.asarray(layer.data_to_world([0] * layer.ndim))
        dy = float(origin[-2])
        dx = float(origin[-1])
        offsets[fname] = {"dx": dx, "dy": dy}

    out = folder / "offsets.json"
    with open(out, "w") as f:
        json.dump(offsets, f, indent=2)
    print(f"Saved {out}")
    for fname, off in offsets.items():
        print(f"  {fname}: dy={off['dy']:.1f}, dx={off['dx']:.1f}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: python align.py <face_folder>")
    main(sys.argv[1])
