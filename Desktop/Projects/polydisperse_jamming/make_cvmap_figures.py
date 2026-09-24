"""make_cvmap_figures.py -- the CV-at-fixed-weighed-phi drop figure (PREREG_CVMAP_2026-09-24).

C1: Rg_r/D0 against t/tau at x0 = 0.932 (x at CV 0), U0 = 4 m/s, for CV 0 / 0.2 / 0.5 and the
    window-only arm W; (a) default EPS_MIN 0.01, (b) EPS_MIN 0.03 (Amendment 2). Healthy cells only.
Explanations in FIGURES_CVMAP.md; house style from make_figures.py.

    python make_cvmap_figures.py
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

from make_figures import COL, save   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
X0, U0 = "0.932", "4.0"
ARMS = (("c0.0", "CV 0", COL[0], "-"), ("c0.2", "CV 0.2", COL[1], "-"),
        ("c0.5", "CV 0.5", COL[2], "-"), ("W", "window only", COL[3], ":"))


def fig_c1():
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.6), sharey=True)
    rows = []
    for ax, camp, tag in ((axs[0], "cvmap", "EPS_MIN 0.01"), (axs[1], "cvmapE03", "EPS_MIN 0.03")):
        for arm, lab, col, ls in ARMS:
            p = os.path.join(OUT, camp, f"{arm}_x{X0}_U{U0}", "spread.dat")
            a = np.loadtxt(p, comments="#", ndmin=2)
            assert a[-1, 1] >= 5.9, (camp, arm, a[-1, 1])       # healthy cells only
            ax.plot(a[:, 1], a[:, 11], ls, color=col, lw=1.2, label=lab)
            rows += [dict(cap=tag, arm=arm, t_over_tau=float(t), rg_r=float(r)) for t, r in a[:, [1, 11]]]
        ax.set_title(tag, fontsize=9)
        ax.set_xlim(0, 6)
        ax.set_xlabel(r"$t/\tau$")
    axs[0].set_ylim(0, 0.7)          # one explicit limit for both shared panels
    axs[0].set_ylabel(r"$R_{g,r}/D_0$")
    axs[1].legend(frameon=False, loc="lower right", fontsize=8)
    fig.tight_layout()
    save(fig, "C1_cvmap_spread", rows, f"PREREG_CVMAP: Rg_r(t), x0 {X0}, U0 {U0}, L9")


def main():
    fig_c1()
    return 0


if __name__ == "__main__":
    sys.exit(main())
