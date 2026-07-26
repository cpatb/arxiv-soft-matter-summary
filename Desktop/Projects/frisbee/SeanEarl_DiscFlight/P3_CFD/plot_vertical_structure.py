"""
Figures for the vertical structure of the flow over Sean Earl Field.

Fig A  pitch-averaged vertical profiles: |U|, sigma_u, TI, dU/dz, each with the
       empty-fetch counterfactual overlaid and the disc band shaded. This is
       the figure that answers "how high can I safely throw".

Fig B  vertical (x-z) slice along the pitch centreline for the easterly case,
       showing mean speed and turbulence intensity with the buildings masked
       out. This is the figure that shows WHY: the Campus Towers wake and the
       shear layer that caps it.

Building interiors are masked from every panel using the 3-D fluid mask, not
drawn as zero-velocity flow.
"""

from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402

from vertical_structure import (PITCH_X, PITCH_Y, DISC_BAND,   # noqa: E402
                                pitch_mask, profiles)

COLORS = {90: "#c1272d", 210: "#0b7a3b", 270: "#1f4e9c"}
LABEL = {90: r"$90^\circ$ easterly", 210: r"$210^\circ$ SSW",
         270: r"$270^\circ$ westerly"}


def fig_profiles(F, C, outpath, zmax=60.0):
    headings = [int(h) for h in F["headings"]]
    P = {h: profiles(F, h) for h in headings}
    cp = profiles(C, int(C["headings"][0])) if C is not None else None

    panels = [("U", r"mean speed $|U|$  [m s$^{-1}$]", False),
              ("sigma", r"fluctuation $\sigma_u=\sqrt{2k/3}$  [m s$^{-1}$]",
               False),
              ("TI", r"turbulence intensity $\sigma_u/|U|$", False),
              ("dUdz", r"vertical shear $\partial|U|/\partial z$  [s$^{-1}$]",
               True)]

    fig, axes = plt.subplots(1, 4, figsize=(15, 5.6), sharey=True)
    for ax, (key, xlabel, zeroline) in zip(axes, panels):
        for h in headings:
            p = P[h]
            ax.plot(p[key], p["z"], "-o", ms=3, lw=1.8, color=COLORS[h],
                    label=LABEL[h])
        if cp is not None:
            ax.plot(cp[key], cp["z"], "--", lw=2.0, color="0.35",
                    label="no buildings (control)")
        ax.axhspan(DISC_BAND[0], DISC_BAND[1], color="#ffd54f", alpha=0.30,
                   zorder=0)
        if zeroline:
            ax.axvline(0.0, color="0.6", lw=0.8, zorder=0)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.grid(alpha=0.3)
        ax.set_ylim(0, zmax)
    axes[0].set_ylabel("height above ground  $z$  [m]")
    axes[0].annotate("disc band\n1–6 m", xy=(0.97, 3.5),
                     xycoords=("axes fraction", "data"), ha="right",
                     va="center", fontsize=8, color="#7a5c00")
    # Campus Towers roof height, the source of the shear layer.
    for ax in axes:
        ax.axhline(44.8, color="0.45", lw=1.0, ls=":", zorder=0)
    axes[-1].annotate("Campus Towers roof (44.8 m)", xy=(0.98, 44.8),
                      xycoords=("axes fraction", "data"), ha="right",
                      va="bottom", fontsize=8, color="0.35")
    axes[0].legend(fontsize=8, loc="upper right", framealpha=0.9)
    fig.suptitle("Vertical structure over the playing surface "
                 "(pitch-averaged, building interiors masked)", fontsize=11)
    fig.tight_layout()
    fig.savefig(outpath, dpi=160)
    plt.close(fig)
    print(f"wrote {outpath}")


def fig_slice(F, theta, outpath, zmax=80.0):
    """x-z slice along the pitch centreline: speed and TI, buildings masked."""
    x, y, z = F["x"], F["y"], F["z"]
    j = int(np.argmin(np.abs(y - 0.5 * (PITCH_Y[0] + PITCH_Y[1]))))
    U = F[f"dir_{theta:03d}__U"][:, j]                 # (nx,nz,3)
    k = F[f"dir_{theta:03d}__k"][:, j]                 # (nx,nz)
    ok = F[f"dir_{theta:03d}__valid3d"][:, j]          # (nx,nz)

    mag = np.linalg.norm(U, axis=-1)
    sig = np.sqrt(2.0 * np.maximum(k, 0.0) / 3.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        TI = sig / np.maximum(mag, 1e-6)
    mag = np.where(ok, mag, np.nan)
    TI = np.where(ok, TI, np.nan)

    Z, X = np.meshgrid(z, x)
    fig, axes = plt.subplots(2, 1, figsize=(11, 7.2), sharex=True)
    for ax, (fldd, ttl, cmap, vmax) in zip(
            axes,
            [(mag, r"mean speed $|U|$  [m s$^{-1}$]", "viridis", None),
             (TI, r"turbulence intensity $\sigma_u/|U|$", "magma", 1.0)]):
        pc = ax.pcolormesh(X, Z, fldd, cmap=cmap, shading="nearest",
                           vmin=0, vmax=vmax)
        fig.colorbar(pc, ax=ax, pad=0.01).set_label(ttl, fontsize=9)
        # Buildings as solid grey where the mask says there is no fluid.
        ax.pcolormesh(X, Z, np.where(ok, np.nan, 1.0), cmap="Greys",
                      vmin=0, vmax=1.6, shading="nearest")
        ax.axhspan(DISC_BAND[0], DISC_BAND[1], color="#ffd54f", alpha=0.35,
                   zorder=3)
        ax.axvspan(PITCH_X[0], PITCH_X[1], color="none", ec="w", lw=1.2,
                   zorder=4)
        ax.set_ylim(0, zmax)
        ax.set_ylabel("z  [m]")
        ax.grid(alpha=0.2)
    axes[0].set_title(f"Centreline slice, "
                      r"$\theta_{\rm wind}=$" + f"{theta}"
                      r"$^\circ$ — wind from the right; "
                      "the pitch (white box) sits inside the tower wake",
                      fontsize=10)
    axes[-1].set_xlabel("x  [m]   (east positive; playing surface "
                        f"{PITCH_X[0]:g} to {PITCH_X[1]:g} m)")
    fig.tight_layout()
    fig.savefig(outpath, dpi=160)
    plt.close(fig)
    print(f"wrote {outpath}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fields", default="vertical_fields.npz")
    ap.add_argument("--control", default="abl_control.npz")
    ap.add_argument("--outdir", default="../P5_Synthesis/figs")
    args = ap.parse_args()

    F = np.load(args.fields, allow_pickle=False)
    try:
        C = np.load(args.control, allow_pickle=False)
    except FileNotFoundError:
        C = None
        print(f"!! no control {args.control}; profiles drawn without it")

    fig_profiles(F, C, f"{args.outdir}/fig5_vertical_profiles.png")
    fig_slice(F, 90, f"{args.outdir}/fig6_centreline_slice.png")


if __name__ == "__main__":
    main()
