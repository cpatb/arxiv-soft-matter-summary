"""Visual and numerical validation of the segmentation, before any analysis runs.

Produces, for sample frames from every clip:
  - the school envelope outline over the raw frame
  - individual detections marked, coloured by edge/core
  - a zoomed crop so individual silhouettes can actually be checked by eye
and prints the numbers that would reveal a silent failure (component count,
area fraction, detection count, body-length estimate, completeness).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from catfish import frames, segment  # noqa: E402

FIGDIR = frames.ROOT / "figures" / "diagnostics"


def validate_frame(clip: frames.Clip, index: int) -> dict:
    rgb = clip.rgb(index)
    env = segment.school_envelope(rgb)
    det = segment.detect_individuals(rgb, envelope=env)
    comp = segment.detection_completeness(rgb, det, env, segment.NOMINAL_BODY_PX)

    clean = det[
        (det["eccentricity"] >= segment.FISH_MIN_ECCENTRICITY) & (~det["in_core"])
    ] if len(det) else det
    body_est = float(np.median(clean["major_px"])) if len(clean) else np.nan

    info = {
        "clip": clip.name,
        "frame": index,
        "t_s": round(clip.time_s(index), 2),
        "n_components": env.n_components,
        "env_area_frac": round(float(env.mask.mean()), 4),
        "n_detections": len(det),
        "n_edge": int((~det["in_core"]).sum()) if len(det) else 0,
        "n_core": int(det["in_core"].sum()) if len(det) else 0,
        "body_px_est": round(body_est, 1),
        "packing_fraction": round(comp.get("packing_fraction", np.nan), 3),
        "areal_completeness": round(comp.get("areal_completeness", np.nan), 3),
    }

    # --- figure -----------------------------------------------------------
    ys, xs = np.where(env.mask)
    if len(ys) == 0:
        return info
    pad = 80
    r0, r1 = max(0, ys.min() - pad), min(rgb.shape[0], ys.max() + pad)
    c0, c1 = max(0, xs.min() - pad), min(rgb.shape[1], xs.max() + pad)

    fig, axes = plt.subplots(1, 3, figsize=(19, 7))

    axes[0].imshow(rgb)
    axes[0].contour(env.mask, levels=[0.5], colors="red", linewidths=1.2)
    axes[0].set_title(f"envelope: {env.n_components} component(s), "
                      f"{100*env.mask.mean():.1f}% of frame", fontsize=9)

    axes[1].imshow(rgb[r0:r1, c0:c1])
    if len(det):
        edge = det[~det["in_core"]]
        core = det[det["in_core"]]
        axes[1].scatter(edge["x"] - c0, edge["y"] - r0, s=5, c="lime",
                        marker="o", linewidths=0, label=f"edge n={len(edge)}")
        axes[1].scatter(core["x"] - c0, core["y"] - r0, s=5, c="magenta",
                        marker="o", linewidths=0, label=f"core n={len(core)}")
        axes[1].legend(fontsize=7, loc="upper right", framealpha=0.8)
    axes[1].set_title(f"detections (n={len(det)}), BL={body_est:.1f} px", fontsize=9)

    # Zoom on a 300 px box at the school edge, where detection should be cleanest.
    zh = 300
    zr = int(np.clip(ys.mean() - zh / 2, 0, rgb.shape[0] - zh))
    zc = int(np.clip(xs.min() + 40, 0, rgb.shape[1] - zh))
    axes[2].imshow(rgb[zr:zr + zh, zc:zc + zh])
    if len(det):
        m = ((det["y"] >= zr) & (det["y"] < zr + zh)
             & (det["x"] >= zc) & (det["x"] < zc + zh))
        d = det[m]
        axes[2].scatter(d["x"] - zc, d["y"] - zr, s=42, facecolors="none",
                        edgecolors="red", linewidths=0.9)
    axes[2].set_title(f"zoom {zh}x{zh} px at school edge", fontsize=9)

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(f"{clip.name}  frame {index}  (t={clip.time_s(index):.1f} s)", fontsize=10)
    fig.tight_layout()
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / f"validate_{clip.name}_{index:05d}.png"
    fig.savefig(out, dpi=110, bbox_inches="tight")
    plt.close(fig)
    info["figure"] = str(out.relative_to(frames.ROOT))
    return info


def main() -> None:
    results = []
    for clip in frames.all_clips():
        n = len(clip)
        for frac in (0.25, 0.60):
            results.append(validate_frame(clip, int(n * frac)))

    hdr = ["clip", "frame", "t_s", "n_components", "env_area_frac", "n_detections",
           "n_edge", "n_core", "body_px_est", "packing_fraction", "areal_completeness"]
    print("  ".join(f"{h:>17s}" for h in hdr))
    for r in results:
        print("  ".join(f"{str(r.get(h, '')):>17s}" for h in hdr))


if __name__ == "__main__":
    main()
