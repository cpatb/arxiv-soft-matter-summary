"""
Wrinkle annotation navigator.

Usage:
    python annotate.py <folder_of_tifs>

Navigate with ← / → arrow keys or Prev/Next buttons.
Annotations auto-save (_roi.tif and _skeleton.tif) on advance and on close.
"""

import sys
import numpy as np
import tifffile
import napari
from pathlib import Path
from skimage.morphology import skeletonize
from scipy.ndimage import label


def _load_img(path: Path) -> np.ndarray:
    img = tifffile.imread(path)
    if img.ndim == 3 and img.shape[0] < 10:
        return img.max(axis=0)
    return img


def _to_display(img: np.ndarray) -> np.ndarray:
    arr = img.astype(np.float32)
    m = arr.max()
    return arr / m if m > 0 else arr


def _load_annotation(path: Path, shape: tuple) -> np.ndarray:
    roi_path = path.parent / f"{path.stem}_roi.tif"
    if roi_path.exists():
        return (tifffile.imread(roi_path) > 0).astype(np.uint8)
    return np.zeros(shape, dtype=np.uint8)


def _save(path: Path, ann_data: np.ndarray) -> None:
    binary = ann_data > 0
    labeled, _ = label(binary)  # type: ignore[misc]
    tifffile.imwrite(path.parent / f"{path.stem}_roi.tif", labeled.astype(np.int32))
    tifffile.imwrite(path.parent / f"{path.stem}_skeleton.tif", skeletonize(binary).astype(np.uint8))
    print(f"Saved {path.stem}")


def main(folder: str) -> None:
    images = sorted(
        p for p in Path(folder).glob("*.tif")
        if not p.stem.endswith("_roi") and not p.stem.endswith("_skeleton")
    )
    if not images:
        sys.exit(f"No source TIFs found in {folder}")

    idx = [0]

    first = _load_img(images[0])
    first_shape = first.shape[:2] if first.ndim == 3 else first.shape

    viewer = napari.Viewer()
    image_layer = viewer.add_image(_to_display(first), name="image", colormap="gray")
    ann_layer = viewer.add_labels(
        _load_annotation(images[0], first_shape), name="annotation"
    )
    ann_layer.brush_size = 5
    ann_layer.mode = "paint"

    def _update_title() -> None:
        viewer.title = f"{images[idx[0]].name}  [{idx[0] + 1} / {len(images)}]"

    _update_title()

    def _go(delta: int) -> None:
        new = idx[0] + delta
        if not (0 <= new < len(images)):
            return
        _save(images[idx[0]], np.asarray(ann_layer.data))
        idx[0] = new
        img = _load_img(images[idx[0]])
        shape = img.shape[:2] if img.ndim == 3 else img.shape
        image_layer.data = _to_display(img)
        image_layer.reset_contrast_limits()
        ann_layer.data = _load_annotation(images[idx[0]], shape)
        _update_title()

    @viewer.bind_key("Left")
    def _prev(_v: napari.Viewer) -> None:
        _go(-1)

    @viewer.bind_key("Right")
    def _next(_v: napari.Viewer) -> None:
        _go(1)

    from qtpy.QtWidgets import QHBoxLayout, QPushButton, QWidget

    class NavWidget(QWidget):
        def __init__(self) -> None:
            super().__init__()
            layout = QHBoxLayout(self)
            btn_prev = QPushButton("← Prev")
            btn_next = QPushButton("Next →")
            btn_prev.clicked.connect(lambda: _go(-1))
            btn_next.clicked.connect(lambda: _go(1))
            layout.addWidget(btn_prev)
            layout.addWidget(btn_next)

    viewer.window.add_dock_widget(NavWidget(), name="Navigation", area="bottom")

    napari.run()
    _save(images[idx[0]], np.asarray(ann_layer.data))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: python annotate.py <folder>")
    main(sys.argv[1])
