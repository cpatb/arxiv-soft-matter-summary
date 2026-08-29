"""Compare candidate segmentation methods for separating fish from substrate.

The fish are dark and roughly colour-neutral; the confounders are (a) dark green
filamentous algae, (b) shaded substrate, (c) surface glint. This script renders
several candidate discriminators side by side so the choice is made by looking at
the result, not by assuming a threshold works.

Colour convention: skimage.io reads RGB. CIELAB is used because it separates
lightness (L*) from the green-red opponent axis (a*), and the algae are strongly
green (a* < 0) while the fish are near-neutral (a* ~ 0).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage as ndi
from skimage import color, filters, io, morphology

ROOT = Path(__file__).resolve().parent.parent
FRAMES = ROOT / "data" / "frames"
FIGDIR = ROOT / "figures" / "diagnostics"

# Sampled from the preview frames: a fish is roughly 30-45 px long at these
# camera distances. The background scale must exceed that to avoid eating fish.
BACKGROUND_DISK_PX = 75


def load_rgb(clip: str, index: int) -> np.ndarray:
    paths = sorted((FRAMES / clip).glob("*.jpg"))
    img = io.imread(paths[index])
    assert img.dtype == np.uint8 and img.ndim == 3
    return img.astype(np.float64) / 255.0


def candidates(rgb: np.ndarray) -> dict[str, np.ndarray]:
    """Return several candidate binary masks plus the scalar fields behind them."""
    lab = color.rgb2lab(rgb)
    L, a = lab[..., 0], lab[..., 1]

    out: dict[str, np.ndarray] = {}

    # (1) Global Otsu on lightness. Simplest possible baseline.
    out["otsu_L"] = L < filters.threshold_otsu(L)

    # (2) Local (adaptive) threshold: handles slow illumination gradients but the
    #     block must be much larger than a fish or it erases the school interior.
    block = 401  # odd, ~10x a fish length
    out["local_L"] = L < filters.threshold_local(L, block_size=block, offset=4.0)

    # (3) Neutrality-weighted darkness. Algae are green (a* strongly negative);
    #     fish are near-neutral. Score high = dark AND not green.
    #     Normalise both terms to comparable scale before combining.
    darkness = (np.percentile(L, 60) - L) / max(np.percentile(L, 60), 1e-9)
    greenness = np.clip(-a, 0, None) / 20.0  # a* ~ -20 is strongly green
    fishness = darkness - greenness
    out["_fishness"] = fishness
    out["neutral_dark"] = fishness > filters.threshold_otsu(fishness)

    # (4) Morphological background residual: grey-closing with a disk larger than a
    #     fish estimates "the scene without dark objects"; the residual isolates
    #     fish-sized dark features. This is top-hat filtering (black top-hat).
    #     Note: it removes the *interior* of a school larger than the disk, so it
    #     detects individuals, not the school envelope.
    se = morphology.disk(BACKGROUND_DISK_PX)
    closed = ndi.grey_closing(L, footprint=se)
    residual = closed - L
    out["_residual"] = residual
    out["blacktophat"] = residual > filters.threshold_otsu(residual)

    return out


def render(clip: str, index: int) -> None:
    rgb = load_rgb(clip, index)
    cand = candidates(rgb)

    masks = {k: v for k, v in cand.items() if not k.startswith("_")}
    fields = {k: v for k, v in cand.items() if k.startswith("_")}

    ncol = 1 + len(masks) + len(fields)
    fig, axes = plt.subplots(1, ncol, figsize=(4.2 * ncol, 7.5))

    axes[0].imshow(rgb)
    axes[0].set_title("RGB")

    for ax, (name, field) in zip(axes[1 : 1 + len(fields)], fields.items()):
        im = ax.imshow(field, cmap="magma")
        ax.set_title(name.lstrip("_"))
        plt.colorbar(im, ax=ax, fraction=0.03)

    for ax, (name, m) in zip(axes[1 + len(fields) :], masks.items()):
        ov = rgb.copy()
        ov[m] = [1.0, 0.0, 0.0]
        ax.imshow(ov)
        ax.set_title(f"{name}\n{100 * m.mean():.1f}% of frame")

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(f"{clip} frame {index}")
    fig.tight_layout()
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / f"segcompare_{clip}_{index:05d}.png"
    fig.savefig(out, dpi=85, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    for clip, idx in [("IMG_0545", 753), ("IMG_0546", 361), ("IMG_0544", 1066)]:
        render(clip, idx)
