"""Image-based figures: the schools themselves, and a fission sequence.

These are the figures a reader looks at first, and the only ones that let
someone check the segmentation by eye rather than trusting a number.
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

sys.path.insert(0, str(Path(__file__).resolve().parent))

from catfish import flow, frames, plotstyle, segment  # noqa: E402

RESULTS = frames.ROOT / "results"
FIGS = frames.ROOT / "figures"


def _pick_representative(name: str) -> int:
    """A frame with the school present and near-median area."""
    df = pd.read_csv(RESULTS / f"{name}_frames.csv")
    d = df[df["present"] & df["area_px2"].notna()]
    if d.empty:
        return len(frames.load_clip(name)) // 2
    target = d["area_px2"].median()
    return int(d.iloc[(d["area_px2"] - target).abs().argmin()]["frame"])


def _crop_to_school(rgb: np.ndarray, mask: np.ndarray, pad: int = 90):
    ys, xs = np.where(mask)
    if len(ys) == 0:
        return slice(None), slice(None)
    r0, r1 = max(0, ys.min() - pad), min(rgb.shape[0], ys.max() + pad)
    c0, c1 = max(0, xs.min() - pad), min(rgb.shape[1], xs.max() + pad)
    return slice(r0, r1), slice(c0, c1)


def available_clips() -> list[str]:
    """Clips whose pipeline output exists, so partial runs still produce figures."""
    return [c for c in frames.CLIPS if (RESULTS / f"{c}_meta.json").exists()]


def fig_gallery() -> None:
    """One representative frame per clip, with the fitted envelope outlined."""
    names = available_clips()
    if not names:
        print("no clips available; skipping gallery")
        return
    fig, axes = plt.subplots(1, len(names), figsize=(2.9 * len(names), 3.6),
                             squeeze=False)
    axes = axes[0]
    for ax, name in zip(axes, names):
        clip = frames.load_clip(name)
        with open(RESULTS / f"{name}_meta.json") as fh:
            work = json.load(fh)["working_body_px"]
        i = _pick_representative(name)
        rgb = clip.rgb(i)
        env = segment.school_envelope(rgb, body_px=work)
        rs, cs = _crop_to_school(rgb, env.mask)
        ax.imshow(rgb[rs, cs])
        if env.present:
            ax.contour(env.mask[rs, cs], levels=[0.5],
                       colors=[plotstyle.CLIP_COLOUR[name]], linewidths=1.1)
        # Scale bar: 10 body lengths, the only annotation on the image.
        h, w = rgb[rs, cs].shape[:2]
        bar = 10 * work
        ax.plot([0.06 * w, 0.06 * w + bar], [0.94 * h, 0.94 * h],
                color="white", lw=2.4, solid_capstyle="butt")
        ax.plot([0.06 * w, 0.06 * w + bar], [0.94 * h, 0.94 * h],
                color="black", lw=1.2, solid_capstyle="butt")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel(name.replace("IMG_", "clip "))
    fig.tight_layout()
    fig.savefig(FIGS / "fig0_gallery.png")
    plt.close(fig)
    print("wrote fig0_gallery.png")


def _find_multicomponent_run(name: str) -> tuple[int, int, str] | None:
    """Longest run of frames with more than one component, and its direction.

    A run of separated components is not automatically a split. It may equally
    be two schools that were already apart and then merged -- which is what the
    longest run in IMG_0546 turned out to be. The direction is therefore read
    off the ends of the run rather than assumed:

        one component before  -> the school came apart      (fission)
        one component after   -> the pieces came together   (fusion)

    Labelling every multi-component episode "fission" would have described a
    merger as a split.
    """
    df = pd.read_csv(RESULTS / f"{name}_frames.csv")
    if "n_components" not in df:
        return None
    d = df[df["present"]].reset_index(drop=True)
    split = (d["n_components"] > 1).to_numpy()
    if not split.any():
        return None
    fr = d["frame"].to_numpy()

    runs, start = [], None
    for k, s in enumerate(split):
        if s and start is None:
            start = k
        elif not s and start is not None:
            runs.append((start, k - 1))
            start = None
    if start is not None:
        runs.append((start, len(split) - 1))
    if not runs:
        return None

    i0, i1 = max(runs, key=lambda r: r[1] - r[0])
    before_one = i0 > 0 and not split[i0 - 1]
    after_one = i1 < len(split) - 1 and not split[i1 + 1]
    if before_one and not after_one:
        kind = "fission"
    elif after_one and not before_one:
        kind = "fusion"
    elif before_one and after_one:
        kind = "fission and refusion"
    else:
        kind = "separated"
    return int(fr[i0]), int(fr[i1]), kind


def fig_fission_sequence() -> None:
    """The longest multi-component episode, shown as a strip through the event."""
    candidates = []
    for name in available_clips():
        r = _find_multicomponent_run(name)
        if r:
            candidates.append((name, r, r[1] - r[0]))
    if not candidates:
        print("no fission/fusion episode found; skipping fig6")
        return
    name, (f0, f1, kind), _ = max(candidates, key=lambda t: t[2])

    clip = frames.load_clip(name)
    with open(RESULTS / f"{name}_meta.json") as fh:
        work = json.load(fh)["working_body_px"]

    span = f1 - f0
    lead = max(15, span // 2)
    idx = np.linspace(max(0, f0 - lead), min(len(clip) - 1, f1 + lead), 6).astype(int)

    fig, axes = plt.subplots(1, len(idx), figsize=(2.0 * len(idx), 3.4))
    for ax, i in zip(axes, idx):
        rgb = clip.rgb(int(i))
        env = segment.school_envelope(rgb, body_px=work)
        rs, cs = _crop_to_school(rgb, env.mask, pad=70)
        ax.imshow(rgb[rs, cs])
        if env.present:
            for lab in range(1, env.n_components + 1):
                ax.contour((env.labels == lab)[rs, cs], levels=[0.5],
                           colors=[plotstyle.PALETTE[(lab - 1) % len(plotstyle.PALETTE)]],
                           linewidths=1.0)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel(f"{clip.time_s(int(i)):.1f} s")
    fig.tight_layout()
    fig.savefig(FIGS / "fig6_fission_fusion.png")
    plt.close(fig)
    print(f"wrote fig6_fission_fusion.png -- {kind.upper()} "
          f"({name}, frames {f0}-{f1}, {(f1-f0)/29.97:.1f} s)")


def fig_velocity_field() -> None:
    """The measured velocity field inside one school, camera motion removed."""
    avail = available_clips()
    if not avail:
        print("no clips available; skipping velocity field")
        return
    name = "IMG_0545" if "IMG_0545" in avail else avail[0]
    with open(RESULTS / f"{name}_meta.json") as fh:
        work = json.load(fh)["working_body_px"]
    clip = frames.load_clip(name)
    i = _pick_representative(name)

    rgb = clip.rgb(i)
    f1 = segment.compute_fields(rgb, body_px=work)
    g1 = np.clip(f1.L * 255 / 100, 0, 255).astype(np.uint8)
    f2 = segment.compute_fields(clip.rgb(i + 1), body_px=work)
    g2 = np.clip(f2.L * 255 / 100, 0, 255).astype(np.uint8)

    env = segment.school_envelope(fields=f1)
    cam = flow.estimate_camera_motion(g1, g2, exclude=env.mask)
    fl = flow.dense_flow(g1, g2)
    s = flow.sample_school_velocity(fl, env.mask, work, camera=cam)

    u, v = s["u"], s["v"]
    du, dv = u - u.mean(), v - v.mean()

    rs, cs = _crop_to_school(rgb, env.mask)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.2))
    for ax, (a, b, ttl) in zip(
        axes, [(u, v, "velocity"), (du, dv, "fluctuation")]
    ):
        ax.imshow(rgb[rs, cs])
        ax.quiver(
            s["x"] - cs.start, s["y"] - rs.start, a, -b,  # -b: image y runs down
            color="white", scale_units="xy", angles="xy",
            scale=0.05, width=0.0035, headwidth=3,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel(ttl)
    fig.tight_layout()
    fig.savefig(FIGS / "fig7_velocity_field.png")
    plt.close(fig)
    print("wrote fig7_velocity_field.png")


def main() -> None:
    plotstyle.apply()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig_gallery()
    fig_fission_sequence()
    fig_velocity_field()


if __name__ == "__main__":
    main()
