"""Figure L2 for PREREG_LFSIG2_2026-09-27.md: the lubricated activation curves by CV, 3 seeds.

    python lfdem/make_lfsig2_figures.py [--check]

Panels: (a) frictional-contact fraction f = Z_c/Z_c(100) vs sigma/sigma_0 per arm (points = seeds, line = the
WC form exp(-(sigma*/sigma)^delta) with the arm-mean parameters); (b) relative viscosity eta_r vs sigma;
(c) contact share of the stress vs sigma. Minimal text on the axes (house style); the explanation lives in
FIGURES_LFSIG2.md. Data are re-read from out/lfsig2/ through the grader's loader; the plotted numbers are
written beside the figure as CSV. --check re-renders and compares the PNG md5 with the stored one.
"""
from __future__ import annotations

import csv
import hashlib
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import grade_lfsig2 as G  # noqa: E402

FIG = os.path.join(os.path.dirname(HERE), "figures")
OUT = os.path.join(os.path.dirname(HERE), "out")
COL = {0.2: "#1f77b4", 0.35: "#2ca02c", 0.5: "#d62728"}
MK = {1: "o", 2: "s", 3: "^"}


def main() -> int:
    check = "--check" in sys.argv
    import glob
    dirs = [x for x in sorted(glob.glob(os.path.join(OUT, "lfsig2", "c*"))) if ".killed_" not in x]
    cells = [c for c in (G.load_cell(d, False) for d in dirs) if c and c["healthy"]]
    by = {(c["cv"], c["seed"], c["sig"]): c for c in cells}
    arms = {(cv, sd): G.arm_readouts(by, cv, sd, G.SIGMAS) for cv in G.CVS for sd in G.SEEDS}
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.3))
    rows = []
    s = np.array(G.SIGMAS)
    sf = np.logspace(0, 2, 200)
    for cv in G.CVS:
        pars = []
        for sd in G.SEEDS:
            r = arms[(cv, sd)]
            if not r["complete"]:
                continue
            ax[0].plot(s, r["f"], MK[sd], color=COL[cv], ms=4, mfc="none")
            ax[1].plot(s, r["eta"], MK[sd], color=COL[cv], ms=4, mfc="none")
            ax[2].plot(s, r["fcont"], MK[sd], color=COL[cv], ms=4, mfc="none")
            if r["r1_ok"]:
                pars.append((r["sstar"], r["delta"]))
            for k, sig in enumerate(G.SIGMAS):
                rows.append([cv, sd, sig, r["f"][k], r["z"][k], r["eta"][k], r["fcont"][k]])
        if pars:
            ss, dd = np.mean([p[0] for p in pars]), np.mean([p[1] for p in pars])
            ax[0].plot(sf, G.wc_form(sf, ss, dd), "-", color=COL[cv], lw=1.2, label=f"CV {cv}")
    # Amendment 4 CV 0 arm (five rungs), same readout; ORDERED cells are not excluded here -- the grader prints Q6
    m_dirs = [x for x in sorted(glob.glob(os.path.join(OUT, "lfsig2m", "c00*"))) if ".killed_" not in x]
    m_cells = {(0.0, c["seed"], c["sig"]): c for c in (G.load_cell(d, False) for d in m_dirs) if c and c["healthy"]}
    s5 = np.array(G.SIGMAS_X)
    for sd in G.SEEDS:
        r = G.arm_readouts(m_cells, 0.0, sd, G.SIGMAS_X)
        if r["complete"]:
            ax[0].plot(s5, r["f"], MK[sd], color="#7f7f7f", ms=4, mfc="none")
            ax[1].plot(s5, r["eta"], MK[sd], color="#7f7f7f", ms=4, mfc="none")
            ax[2].plot(s5, r["fcont"], MK[sd], color="#7f7f7f", ms=4, mfc="none")
            for k, sig in enumerate(G.SIGMAS_X):
                rows.append([0.0, sd, sig, r["f"][k], r["z"][k], r["eta"][k], r["fcont"][k]])
    for a in ax:
        a.set_xscale("log"); a.set_xlabel(r"$\sigma/\sigma_0$")
    ax[0].set_ylabel(r"$Z_c/Z_c(100)$"); ax[0].set_ylim(0, 1.05); ax[0].legend(frameon=False, fontsize=8)
    ax[1].set_ylabel(r"$\eta_r$"); ax[1].set_ylim(0, None)
    ax[2].set_ylabel("contact share"); ax[2].set_ylim(0, 1)
    for a, lab in zip(ax, "abc"):
        a.text(0.02, 0.95, f"({lab})", transform=a.transAxes, va="top")
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    png = os.path.join(FIG, "L2_lfsig2_activation.png")
    pdf = os.path.join(FIG, "L2_lfsig2_activation.pdf")
    csvp = os.path.join(FIG, "L2_lfsig2_activation.csv")
    old = None
    if check:
        old = hashlib.md5(open(png, "rb").read()).hexdigest() if os.path.exists(png) else None
    fig.savefig(png, dpi=200, metadata={"Software": None}); fig.savefig(pdf, metadata={"Creator": None, "Producer": None, "CreationDate": None})
    with open(csvp, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(["cv", "seed", "sigma", "f", "zc", "eta_r", "fcont"]); w.writerows(rows)
    new = hashlib.md5(open(png, "rb").read()).hexdigest()
    print(f"wrote {png} ({len(rows)} rows) md5 {new}" + (f"  check: {'IDENTICAL' if old == new else 'CHANGED'} vs {old}" if check else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
