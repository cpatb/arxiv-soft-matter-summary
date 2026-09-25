"""make_lfsig_figures.py -- figure L1 for PREREG_LFSIG_2026-09-24.md: flow curves from LF_DEM.

    python lfdem/make_lfsig_figures.py   -> figures/L1_lfsig_flowcurves.{pdf,png,csv}

(a) relative viscosity eta_r vs stress sigma/sigma_0 (log-log), (b) contact number per particle Z_c vs
sigma (log x, linear y from 0). CV 0.2 blue, CV 0.5 orange; x_m 0.90 open, 0.95 filled. Points are the
block means over strain [1.5, 4] of HEALTHY cells (grade_lfsig's registered rule), bars 1 block s.e.
Explanations live in FIGURES_LFSIG.md, not on the axes.
"""
from __future__ import annotations

import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import grade_lfsig as G  # noqa: E402

COL = {0.2: "#1f77b4", 0.5: "#ff7f0e"}


def main() -> int:
    import glob
    cells = [c for c in (G.load_cell(d) for d in sorted(glob.glob(os.path.join(G.OUT, "lfsig", "c*"))))
             if c is not None and c["healthy_reg"]]
    fig, ax = plt.subplots(1, 2, figsize=(6.8, 2.8), constrained_layout=True)
    rows = []
    for cv in G.CVS:
        for xm in G.XMS:
            cs = sorted([c for c in cells if c["cv"] == cv and c["xm"] == xm], key=lambda c: c["sig"])
            if not cs:
                continue
            s = [c["sig"] for c in cs]
            kw = dict(color=COL[cv], marker="o", ms=4, lw=1,
                      mfc=COL[cv] if xm == 0.95 else "white", label=f"CV {cv}, x {xm}")
            ax[0].errorbar(s, [c["eta"] for c in cs], yerr=[c["eta_se"] for c in cs], **kw)
            ax[1].errorbar(s, [c["zc"] for c in cs], yerr=[c["zc_se"] for c in cs], **kw)
            for c in cs:
                rows.append([cv, c["cv_real"], xm, c["phi"], c["sig"], c["eta"], c["eta_se"], c["zc"], c["zc_se"]])
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel(r"$\sigma/\sigma_0$"); ax[0].set_ylabel(r"$\eta_r$")
    ax[1].set_xscale("log"); ax[1].set_ylim(0, None)
    ax[1].set_xlabel(r"$\sigma/\sigma_0$"); ax[1].set_ylabel(r"$Z_c$")
    ax[0].legend(frameon=False, fontsize=7)
    os.makedirs(os.path.join(ROOT, "figures"), exist_ok=True)
    base = os.path.join(ROOT, "figures", "L1_lfsig_flowcurves")
    fig.savefig(base + ".pdf"); fig.savefig(base + ".png", dpi=200)
    with open(base + ".csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["cv_nominal", "cv_realised", "x_m", "phi", "sigma_over_sigma0", "eta_r", "eta_se", "Zc", "Zc_se"])
        w.writerows(rows)
    print(f"L1: {len(rows)} points from {len(cells)} healthy cells -> {base}.*")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
