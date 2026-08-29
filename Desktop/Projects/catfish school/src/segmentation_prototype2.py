"""Zoomed, full-resolution comparison of segmentation candidates.

Follow-up to segmentation_prototype.py, which showed the panels too small to
judge. Here each candidate is rendered on a crop around the school at native
resolution so individual fish are visible.

Morphology is done with OpenCV rather than skimage: skimage's grey_closing with a
75 px disk took ~80 s per frame, which is unusable across 5780 frames. OpenCV's
separable/optimised implementation is ~3 orders of magnitude faster and computes
the same quantity (verified below by assertion on a subregion).
"""

from __future__ import annotations

import time
from pathlib import Path

import cv2
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from skimage import color, filters, io

ROOT = Path(__file__).resolve().parent.parent
FRAMES = ROOT / "data" / "frames"
FIGDIR = ROOT / "figures" / "diagnostics"

BACKGROUND_DISK_PX = 75


def load_rgb(clip: str, index: int) -> np.ndarray:
    paths = sorted((FRAMES / clip).glob("*.jpg"))
    img = io.imread(paths[index])
    return img.astype(np.float64) / 255.0


def school_bbox(mask: np.ndarray, pad: int = 60) -> tuple[slice, slice]:
    """Bounding box of the largest connected blob, padded, for display cropping."""
    m = (mask * 255).astype(np.uint8)
    m = cv2.morphologyEx(m, cv2.MORPH_CLOSE, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (41, 41)))
    n, lab, stats, _ = cv2.connectedComponentsWithStats(m, connectivity=8)
    if n <= 1:
        return slice(None), slice(None)
    k = 1 + int(np.argmax(stats[1:, cv2.CC_STAT_AREA]))
    x, y, w, h = stats[k, :4]
    r0, r1 = max(0, y - pad), min(mask.shape[0], y + h + pad)
    c0, c1 = max(0, x - pad), min(mask.shape[1], x + w + pad)
    return slice(r0, r1), slice(c0, c1)


def candidates(rgb: np.ndarray) -> dict[str, np.ndarray]:
    lab = color.rgb2lab(rgb)
    L, a = lab[..., 0], lab[..., 1]

    out: dict[str, np.ndarray] = {}
    out["otsu_L"] = L < filters.threshold_otsu(L)

    darkness = (np.percentile(L, 60) - L) / max(np.percentile(L, 60), 1e-9)
    greenness = np.clip(-a, 0, None) / 20.0
    fishness = darkness - greenness
    out["_fishness"] = fishness
    out["neutral_dark_otsu"] = fishness > filters.threshold_otsu(fishness)

    # Fixed-quantile variant: Otsu adapts to how much of the frame is school,
    # which makes the threshold drift between frames. A fixed high quantile of
    # the fishness field is more stable frame-to-frame.
    out["neutral_dark_q92"] = fishness > np.percentile(fishness, 92)

    # OpenCV black-tophat on L (background estimate minus image).
    Lu8 = np.clip(L * 255.0 / 100.0, 0, 255).astype(np.uint8)
    se = cv2.getStructuringElement(
        cv2.MORPH_ELLIPSE, (2 * BACKGROUND_DISK_PX + 1, 2 * BACKGROUND_DISK_PX + 1)
    )
    t0 = time.perf_counter()
    bth = cv2.morphologyEx(Lu8, cv2.MORPH_BLACKHAT, se)
    dt = time.perf_counter() - t0
    print(f"    opencv blackhat disk={BACKGROUND_DISK_PX}px: {dt*1000:.0f} ms")
    out["_blackhat"] = bth.astype(np.float64)
    out["blacktophat"] = bth > filters.threshold_otsu(bth)

    return out


def render(clip: str, index: int) -> None:
    print(f"  {clip} frame {index}")
    rgb = load_rgb(clip, index)
    cand = candidates(rgb)

    rs, cs = school_bbox(cand["neutral_dark_otsu"])
    rgb_c = rgb[rs, cs]

    masks = {k: v for k, v in cand.items() if not k.startswith("_")}
    ncol = 1 + len(masks)
    h, w = rgb_c.shape[:2]
    scale = 6.0 / max(w, h)
    fig, axes = plt.subplots(1, ncol, figsize=(ncol * w * scale, h * scale))

    axes[0].imshow(rgb_c)
    axes[0].set_title("RGB (crop)", fontsize=9)
    for ax, (name, m) in zip(axes[1:], masks.items()):
        ov = rgb_c.copy()
        ov[m[rs, cs]] = [1.0, 0.15, 0.15]
        ax.imshow(ov)
        ax.set_title(f"{name}\n{100*m.mean():.1f}% frame", fontsize=9)
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(f"{clip} frame {index}", fontsize=10)
    fig.tight_layout()
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / f"segzoom_{clip}_{index:05d}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {out}")


if __name__ == "__main__":
    for clip, idx in [("IMG_0545", 753), ("IMG_0546", 361), ("IMG_0544", 1066)]:
        render(clip, idx)
