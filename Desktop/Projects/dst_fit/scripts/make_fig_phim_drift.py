"""Figure: the corn-starch series, its Wyart-Cates fits, and the phi_m that will not hold still.

Deliberately spare. Axis labels, a phi legend, panel letters, nothing else. What
the panels mean belongs in the sibling fig_phim_drift.md, not printed on the
axes.

Panel (b) uses a linear y axis starting at zero rather than cropping to the data
range. That makes the phi_m drift look small, which is honest: it is a change of
0.09 on a quantity bounded near 0.64. Panel (c) is there because small and
insignificant are different things, and the drift is 32 reduced chi-squared.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dst_fit import fit_model, read_series, thickening_branch  # noqa: E402
from dst_fit.models import WYART_CATES  # noqa: E402

DATA_DIR = Path(r"C:\Users\compl\Desktop\Projects\Corn Starch Bulk Rheology AM")
OUT = Path(__file__).resolve().parents[1] / "figures" / "fig_phim_drift.png"

mpl.rcParams.update({
    "font.size": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.7,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "lines.linewidth": 1.0,
    "savefig.dpi": 400,
    "figure.dpi": 130,
})


def main() -> None:
    sweeps = read_series(sorted(DATA_DIR.glob("*_e1.xls")))
    colors = plt.cm.viridis(np.linspace(0.05, 0.9, len(sweeps)))

    fits, phim, err, nominal = [], [], [], []
    for s in sweeps:
        br = thickening_branch(s.viscosity)
        f = fit_model("wyart_cates", s.stress, s.viscosity, branch=br, n_starts=12)
        j = WYART_CATES.param_names.index("phi_m")
        fits.append((s, br, f))
        phim.append(f.params[j])
        err.append(f.stderr[j])
        nominal.append(s.phi)

    phim, err, nominal = np.array(phim), np.array(err), np.array(nominal)
    w = 1.0 / err**2
    mean = float(np.sum(w * phim) / np.sum(w))
    dev = (phim - mean) / err

    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.35), constrained_layout=True)

    # (a) flow curves and the fits
    for (s, br, f), c in zip(fits, colors):
        ax[0].loglog(s.stress, s.viscosity, "o", ms=2.2, mfc="none", mew=0.6, color=c,
                     label=f"{s.phi:.2f}")
        grid = np.logspace(np.log10(s.stress[br].min()), np.log10(s.stress[br].max()), 200)
        ax[0].loglog(grid, WYART_CATES(grid, *f.params), "-", color=c, lw=1.0)
    ax[0].set_xlabel(r"$\sigma$  (Pa)")
    ax[0].set_ylabel(r"$\eta$  (Pa$\,$s)")
    ax[0].legend(title=r"$\phi_{\rm nom}$", frameon=False, ncol=2,
                 handlelength=1.0, labelspacing=0.25, columnspacing=0.8,
                 borderpad=0.2, title_fontsize=6.5, loc="upper left")

    # (b) phi_m against the nominal fraction, zero-based
    ax[1].axhline(mean, color="0.6", lw=0.8, zorder=0)
    ax[1].errorbar(nominal, phim, yerr=err, fmt="s", ms=3.4, lw=0.9,
                   capsize=2, color="C3", zorder=3)
    ax[1].plot(nominal, nominal, "--", color="0.35", lw=0.8, zorder=1)
    ax[1].set_ylim(0.0, 0.68)
    ax[1].set_xlim(0.0, 0.68)
    ax[1].set_xlabel(r"$\phi_{\rm nom}$")
    ax[1].set_ylabel(r"$\phi_m$")

    # (c) the same points in units of their own uncertainty
    ax[2].axhline(0.0, color="0.6", lw=0.8, zorder=0)
    ax[2].axhspan(-2, 2, color="0.88", zorder=0)
    ax[2].bar(nominal, dev, width=0.012, color=colors, zorder=3)
    ax[2].set_xlim(0.0, 0.68)
    ax[2].set_xlabel(r"$\phi_{\rm nom}$")
    ax[2].set_ylabel(r"$(\phi_m-\overline{\phi_m})/\sigma$")

    for a, letter in zip(ax, "abc"):
        a.text(0.03, 0.97, f"({letter})", transform=a.transAxes,
               va="top", ha="left", fontsize=8)

    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    print(f"wrote {OUT}")
    print(f"weighted mean phi_m = {mean:.4f}")
    print("deviations (sigma): " + ", ".join(f"{d:+.1f}" for d in dev))


if __name__ == "__main__":
    main()
