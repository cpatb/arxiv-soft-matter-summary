"""make_type_figures.py -- figure T1, the thickening-type map (derivations/thickening_type_2026-09-25.md).

T1 (a) the CST / DST / SJ regions in the (chi, x) plane, x = phi/phi_0, chi = phi_m/phi_0, with the DST boundary
       x_DST(chi, n*delta) drawn for the n*delta values of the cvmap inputs (2.147 CV 0.5 delta, 2.634 CV 0 delta)
       and for n = 5 (5 * 1.3169); the SJ edge is x = chi.
   (b) the 24 cvmap DNS drops (EPS_MIN 0.03) at their (chi, x), marker = type by their own n*delta, colour = peak
       spread Rg_max. Inputs: out/cvmap_cells.txt, derivations/cvmap_by_type_2026-09-26.txt (itself read from
       out/grade_cvmapE03.txt).
Explanations in FIGURES_TYPE.md; house style from make_figures.py.

    python make_type_figures.py [--check]
"""
import argparse
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402
from scipy.optimize import brentq  # noqa: E402

from make_figures import COL, FIG, md5, save   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def hmax(x, chi, nd):
    uu = np.logspace(-6, 3, 20001)
    J = 1 - (1 - chi) * np.exp(-uu)
    return float((nd * x * (1 - chi) * uu * np.exp(-uu) / (J * (J - x))).max())


def x_dst(chi, nd):
    f = lambda x: hmax(x, chi, nd) - 1.0  # noqa: E731
    if f(chi * (1 - 1e-9)) < 0:
        return np.nan
    return brentq(f, 1e-6, chi * (1 - 1e-9), xtol=1e-9)


def build():
    rows = []
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.0))
    chis = np.linspace(0.70, 0.99, 59)
    ax[0].plot(chis, chis, "-", color=COL[0], lw=1.0)
    for nd, ls, c in ((2.147, "--", COL[1]), (2.634, "-.", COL[2]), (5 * 1.3169, ":", COL[3])):
        xd = [x_dst(c_, nd) for c_ in chis]
        ax[0].plot(chis, xd, ls, color=c, lw=0.9)
        rows += [{"panel": "a", "n_delta": nd, "chi": float(ch), "x_dst": float(x)} for ch, x in zip(chis, xd)]
    # (b) cvmap drops: peak spread against the distance from each drop's OWN DST boundary
    txt = open(os.path.join(HERE, "derivations", "cvmap_by_type_2026-09-26.txt"), encoding="utf-8").read()
    pts = re.findall(r"^(\S+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(SJ|DST|CST)\s+\|\s+([\d.]+)\s+([\d.]+)",
                     txt, re.M)
    mk = {"SJ": ("s", COL[0]), "DST": ("D", COL[2]), "CST": ("o", COL[1])}
    for t in ("SJ", "DST", "CST"):
        sel = [p for p in pts if p[5] == t]
        m, c = mk[t]
        ax[1].plot([float(p[1]) - float(p[4]) for p in sel], [float(p[7]) for p in sel], m, color=c, ms=4,
                   mfc=c, mec="k", mew=0.4, lw=0)
        rows += [{"panel": "b", "label": p[0], "chi": float(p[2]), "x": float(p[1]), "type": p[5],
                  "Rg_max": float(p[7])} for p in sel]
    ax[1].axvline(0.0, color="0.5", lw=0.7, ls="--")
    ax[1].set_ylim(0, 0.75)
    ax[1].set_xlabel(r"$x - x_{DST}$")
    ax[1].set_ylabel(r"peak spread $R_g/D_0$")
    ax[0].set_xlim(0.70, 0.99)
    ax[0].set_ylim(0.0, 1.0)
    ax[0].set_xlabel(r"$\chi = \phi_m/\phi_0$")
    ax[0].set_ylabel(r"$x = \phi/\phi_0$")
    for a_, lab in zip(ax, "ab"):
        a_.text(0.02, 0.97, f"({lab})", transform=a_.transAxes, va="top")
    fig.tight_layout()
    keys = ["panel", "n_delta", "chi", "x_dst", "label", "x", "type", "Rg_max"]
    rows = [{k: r.get(k, "") for k in keys} for r in rows]
    save(fig, "T1_thickening_type", rows,
         "make_type_figures.py; x_DST from the closed form (T); cvmap drops from out/grade_cvmapE03.txt via "
         "derivations/cvmap_by_type_2026-09-26.txt; phi_0 0.62, n 2")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    build()
    if a.check:
        h1 = md5(os.path.join(FIG, "T1_thickening_type.png"))
        build()
        h2 = md5(os.path.join(FIG, "T1_thickening_type.png"))
        print(f"T1 png md5 {h1} / {h2}: {'IDENTICAL' if h1 == h2 else 'DIFFERS'}")


if __name__ == "__main__":
    main()
