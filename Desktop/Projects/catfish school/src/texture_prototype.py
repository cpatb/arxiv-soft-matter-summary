"""Test a texture cue to separate the school from its own cast shadow.

Motivation: the 3-sigma-above-background envelope threshold failed on
IMG_0547 frame 559 because that frame's background is genuinely heterogeneous
(bright sand, dark shadow, dark water), which inflates the MAD until the
threshold exceeds the field maximum. Separately, the school casts a diffuse
shadow on the lake bed that is dark and colour-neutral, so the existing
"dark and not green" score cannot reject it.

Both problems have the same fix: the school is dark AND finely textured (many
resolvable silhouettes at fish scale), whereas a shadow is dark and smooth. This
script renders the candidate fields so the discriminator is chosen by looking.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage as ndi
from skimage import color, filters

sys.path.insert(0, str(Path(__file__).resolve().parent))
from catfish import frames, segment  # noqa: E402

FIGDIR = frames.ROOT / "figures" / "diagnostics"
BODY_PX = 32.0


def local_texture(L: np.ndarray, sigma: float) -> np.ndarray:
    """Local standard deviation of lightness at scale `sigma`.

    Computed as sqrt(E[L^2] - E[L]^2) with the variance clipped at zero: the
    subtraction is not guaranteed non-negative in floating point and a bare
    sqrt would emit NaNs on flat regions.
    """
    mu = ndi.gaussian_filter(L, sigma, mode="nearest")
    mu2 = ndi.gaussian_filter(L * L, sigma, mode="nearest")
    var = np.clip(mu2 - mu * mu, 0.0, None)
    return np.sqrt(var)


def fields(rgb: np.ndarray) -> dict[str, np.ndarray]:
    lab = color.rgb2lab(rgb)
    L = lab[..., 0]

    fish = segment.fishness_field(rgb)
    # Texture at half a body length: responds to the alternation between fish
    # silhouette and interstitial gap, which is the signature of a packed school.
    tex = local_texture(L, 0.5 * BODY_PX)
    # Smooth both to the same scale so their product is a coherent field.
    fish_s = ndi.gaussian_filter(np.clip(fish, 0, None), 0.5 * BODY_PX, mode="nearest")
    tex_s = ndi.gaussian_filter(tex, 0.5 * BODY_PX, mode="nearest")

    # Normalise texture by its own high percentile so the product is scale-free.
    tex_n = tex_s / max(float(np.percentile(tex_s, 99)), 1e-9)

    return {
        "fishness (dark & not green)": fish_s,
        "texture (local SD of L*)": tex_n,
        "product = school score": fish_s * tex_n,
    }


def render(clip_name: str, index: int) -> None:
    clip = frames.load_clip(clip_name)
    rgb = clip.rgb(index)
    fs = fields(rgb)

    ncol = 1 + len(fs) + 1
    fig, axes = plt.subplots(1, ncol, figsize=(4.6 * ncol, 8))
    axes[0].imshow(rgb)
    axes[0].set_title("RGB", fontsize=9)

    for ax, (name, f) in zip(axes[1:], fs.items()):
        im = ax.imshow(f, cmap="magma")
        ax.set_title(name, fontsize=9)
        plt.colorbar(im, ax=ax, fraction=0.035)

    score = fs["product = school score"]
    t = filters.threshold_otsu(score)
    m = score > t
    ov = rgb.copy()
    ov[m] = [1.0, 0.1, 0.1]
    axes[-1].imshow(ov)
    axes[-1].set_title(f"Otsu on score\nthr={t:.4f}, {100*m.mean():.1f}% of frame", fontsize=9)

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(f"{clip_name} frame {index}", fontsize=10)
    fig.tight_layout()
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / f"texture_{clip_name}_{index:05d}.png"
    fig.savefig(out, dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}  (Otsu thr={t:.4f}, area={100*m.mean():.2f}%)")


if __name__ == "__main__":
    # 0547/559 is the failure case (school + strong cast shadow);
    # 0546/482 is a known-good case; 0545/1004 has no school at all.
    for c, i in [("IMG_0547", 559), ("IMG_0546", 482), ("IMG_0545", 1004)]:
        render(c, i)
