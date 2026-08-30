"""Map the disc-band turbulence over the pitch, and compare designs against it.

Three views of the same field, because the zone this study is about is a
three-dimensional object and any one projection of it can be argued with:

  PLAN      sigma_u at z = 2 m over the site, with the pitch outlined. Shows
            where on the field the zone is.
  SECTION   sigma_u on the pitch centreline, x-z, with the disc band marked.
            Shows how far up it reaches and whether a barrier has put a shear
            layer inside the band.
  PROFILE   sigma_u along the pitch, y-averaged over the playing surface, in
            the disc band. This is the curve the gates in `windbreak_gates.py`
            are written against.

Passing `--designs W1 W2 W3` overlays each design's profile on the baseline's
and draws one section per design, so the comparison is on one figure rather
than across three.

COLOUR AND AXES. The plan and section share one colour scale across every
panel drawn in a run -- a per-panel scale would make a design that halved the
turbulence look identical to the baseline. The profile's y axis starts at zero.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

from windbreak_gates import DISC_BAND, PITCH_X, PITCH_Y, DESIGNS  # noqa: E402


def load(path, heading):
    F = np.load(path, allow_pickle=False)
    k = F[f"dir_{heading:03d}__k"].astype(float)
    solid = F[f"dir_{heading:03d}__solid"]
    sig = np.sqrt(np.maximum(2.0 / 3.0 * k, 0.0))
    sig[solid] = np.nan
    return F["x"], F["y"], F["z"], sig


def profile(x, y, z, sig):
    """sigma_u along the pitch, averaged over the pitch width and disc band."""
    iy = (y >= PITCH_Y[0]) & (y <= PITCH_Y[1])
    iz = (z >= DISC_BAND[0]) & (z <= DISC_BAND[1])
    sub = sig[:, iy, :][:, :, iz]
    return np.nanmean(sub, axis=(1, 2))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cfd", default=str(ROOT / "P3_CFD"))
    ap.add_argument("--heading", type=int, default=270)
    ap.add_argument("--designs", nargs="*", default=[])
    ap.add_argument("--figdir", default=str(HERE / "figs_windbreak"))
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cfd = Path(args.cfd)
    cases = [("baseline", cfd / "vertical_fields.npz")]
    for t in args.designs:
        p = cfd / f"vertical_fields_dir_{args.heading:03d}_{t}.npz"
        if p.exists():
            cases.append((t, p))
        else:
            print(f"  {t}: {p.name} not written yet -- skipping")

    data = [(tag, *load(p, args.heading)) for tag, p in cases]
    # One colour scale for every panel in the run.
    vmax = float(np.nanpercentile(
        np.concatenate([d[4][np.isfinite(d[4])].ravel() for d in data]), 99))
    print(f"shared colour scale 0 to {vmax:.3f} m/s over {len(data)} case(s)")

    outdir = Path(args.figdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- plan + section, one row per case --------------------------------
    fig, axes = plt.subplots(len(data), 2, figsize=(11.5, 3.1 * len(data)),
                             squeeze=False)
    for row, (tag, x, y, z, sig) in enumerate(data):
        kz = int(np.argmin(np.abs(z - 2.0)))
        ax = axes[row][0]
        im = ax.pcolormesh(x, y, sig[:, :, kz].T, cmap="magma", vmin=0,
                           vmax=vmax, shading="nearest")
        ax.plot([PITCH_X[0], PITCH_X[1], PITCH_X[1], PITCH_X[0], PITCH_X[0]],
                [PITCH_Y[0], PITCH_Y[0], PITCH_Y[1], PITCH_Y[1], PITCH_Y[0]],
                color="#7ef", lw=1.4)
        if tag in DESIGNS:
            ax.axvline(DESIGNS[tag]["x"], color="#4f8", lw=2.0, ls="--")
        ax.set_xlim(-60, 90)
        ax.set_ylim(-40, 45)
        ax.set_aspect("equal")
        ax.set_ylabel(f"{tag}\ny (m)", fontsize=9)
        if row == len(data) - 1:
            ax.set_xlabel("x (m)   east is +", fontsize=9)
        ax.tick_params(labelsize=8)
        if row == 0:
            ax.set_title(r"plan: $\sigma_u$ at z = 2 m", fontsize=10)

        jy = int(np.argmin(np.abs(y - 2.5)))
        ax = axes[row][1]
        iz = z <= 30.0
        ax.pcolormesh(x, z[iz], sig[:, jy, :][:, iz].T, cmap="magma", vmin=0,
                      vmax=vmax, shading="nearest")
        ax.axhspan(DISC_BAND[0], DISC_BAND[1], color="#7ef", alpha=0.18, lw=0)
        for zz in DISC_BAND:
            ax.axhline(zz, color="#7ef", lw=1.0)
        if tag in DESIGNS:
            d = DESIGNS[tag]
            ax.plot([d["x"], d["x"]], [0, d["height"]], color="#4f8", lw=3.0)
        ax.set_xlim(-60, 90)
        ax.set_ylim(0, 30)
        ax.set_ylabel("z (m)", fontsize=9)
        if row == len(data) - 1:
            ax.set_xlabel("x (m)", fontsize=9)
        ax.tick_params(labelsize=8)
        if row == 0:
            ax.set_title(r"centreline section: $\sigma_u$", fontsize=10)
    cb = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.02)
    cb.set_label(r"$\sigma_u$  (m/s)", fontsize=9)
    fig.savefig(outdir / f"fig_zone_{args.heading:03d}.png", dpi=190,
                facecolor="white", bbox_inches="tight")
    plt.close(fig)

    # ---- profile along the pitch -----------------------------------------
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    inp = None
    for tag, x, y, z, sig in data:
        p = profile(x, y, z, sig)
        inp = (x >= PITCH_X[0]) & (x <= PITCH_X[1])
        ax.plot(x[inp], p[inp], lw=1.8,
                label=f"{tag}" + (f"  ({DESIGNS[tag]['name']})"
                                  if tag in DESIGNS else ""))
        if tag in DESIGNS:
            ax.axvline(DESIGNS[tag]["x"], color="#888", lw=0.8, ls=":")
    for a, b, nm in ((-31.0, -5.0, "west"), (-5.0, 21.0, "middle"),
                     (21.0, 47.0, "east")):
        ax.axvline(b, color="#ccc", lw=0.7)
        ax.text(0.5 * (a + b), 0.02, nm, ha="center", fontsize=8,
                color="#888", transform=ax.get_xaxis_transform())
    ax.set_xlabel("x along the pitch (m)   east is +")
    ax.set_ylabel(r"$\sigma_u$ in the disc band  (m/s)")
    ax.set_xlim(PITCH_X)
    ax.set_ylim(0, None)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / f"fig_profile_{args.heading:03d}.png", dpi=200,
                facecolor="white")
    plt.close(fig)

    print(f"wrote figures in {outdir}")
    for tag, x, y, z, sig in data:
        p = profile(x, y, z, sig)
        inp = (x >= PITCH_X[0]) & (x <= PITCH_X[1])
        xe = x[inp]
        pe = p[inp]
        east = (xe >= 21.0)
        west = (xe <= -5.0)
        print(f"  {tag:<9} disc-band sigma_u: west {np.nanmean(pe[west]):.3f}, "
              f"east {np.nanmean(pe[east]):.3f} m/s, "
              f"ratio {np.nanmean(pe[east]) / np.nanmean(pe[west]):.2f}x")


if __name__ == "__main__":
    main()
