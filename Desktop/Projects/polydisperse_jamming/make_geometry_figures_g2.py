"""make_geometry_figures_g2.py -- figure G2: out-of-sample tests of the contact-geometry laws on binary and
ternary packs (derivations/binary_contact_check.py; rattler_closure_2026-09-24.md §4d).

(a) per-class rattler fraction, predicted (Poisson, mean ~ a^2, no parameter) vs measured; binaries by size
    ratio, ternary 1:2:4 separately; the size-blind prediction in grey for contrast. Identity line.
(b) SL contact count measured / configuration-model prediction against size ratio (one point per record).
Explanations in FIGURES_GEOMETRY.md.

    python make_geometry_figures_g2.py
"""
import glob
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402
from scipy.stats import poisson   # noqa: E402

from make_figures import COL, save   # noqa: E402
from psd import BinaryPSD            # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")


def binaries():
    recs = {}
    for f in sorted(glob.glob(os.path.join(OUT, "phase12*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            bc = r.get("by_class")
            if (not bc or not r.get("converged") or "mix_ratio" not in r or "x_small_vol" not in r
                    or not (0.0 < r["x_small_vol"] < 1.0) or r["mix_ratio"] <= 1.0
                    or not bc.get("n_S") or not bc.get("n_L") or r.get("mu_p") != 0.3):
                continue
            recs[r["key"]] = r
    return list(recs.values())


def ternaries():
    recs = {}
    for f in sorted(glob.glob(os.path.join(OUT, "phase12_tern*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            bd = r.get("by_design_class")
            if r.get("converged") and bd and len(bd) == 3 and all(bd[c]["n"] for c in bd):
                recs[r["key"]] = r
    return list(recs.values())


def main():
    fig, axs = plt.subplots(1, 2, figsize=(6.8, 2.9))
    rows = []
    ax = axs[0]
    ratios = sorted({r["mix_ratio"] for r in binaries()})
    for i, ratio in enumerate(ratios):
        xs, ys, ns = [], [], []
        for r in binaries():
            if r["mix_ratio"] != ratio:
                continue
            bc = r["by_class"]
            psd = BinaryPSD(d10=1.0, ratio=ratio, x_small=r["x_small_vol"])
            n = np.array([bc["n_S"], bc["n_L"]], float)
            d = np.array([psd.d_small, psd.d_large])
            Z = (bc["n_S"] * bc["z_S"] + bc["n_L"] * bc["z_L"]) / n.sum()
            a2 = np.sum(n * d**2) / n.sum()
            meas = [bc["rat_S_k2"], bc["rat_L_k2"]]
            pred = poisson.cdf(1, Z * d**2 / a2)
            xs += meas; ys += list(pred); ns += [poisson.cdf(1, Z)] * 2
            for m, pr, cls in zip(meas, pred, "SL"):
                rows.append(dict(kind="binary", ratio=ratio, cls=cls, meas=m, pred_a2=pr, pred_blind=poisson.cdf(1, Z)))
        ax.plot(xs, ns, ".", color="0.8", ms=3, zorder=1)
        ax.plot(xs, ys, "o", color=COL[(i % 5) + 1], ms=3, label=f"binary {ratio:.3g}", zorder=3)
    xs, ys = [], []
    for r in ternaries():
        cl = list(r["by_design_class"].values())
        n = np.array([c["n"] for c in cl], float); d = np.array([c["d"] for c in cl], float)
        z = np.array([c["z_mean"] for c in cl], float); fr = [c["rattler_frac"] for c in cl]
        Z = np.sum(n * z) / n.sum(); a2 = np.sum(n * d**2) / n.sum()
        pr = poisson.cdf(1, Z * d**2 / a2)
        xs += fr; ys += list(pr)
        for m, p_ in zip(fr, pr):
            rows.append(dict(kind="ternary", ratio=4.0, cls="", meas=m, pred_a2=p_, pred_blind=poisson.cdf(1, Z)))
    ax.plot(xs, ys, "^", color=COL[0], ms=3, label="ternary 1:2:4", zorder=3)
    ax.plot([0, 1], [0, 1], "-", color="0.5", lw=0.7, zorder=0)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("measured rattler fraction (per class)"); ax.set_ylabel("predicted")
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    ax = axs[1]
    for r in binaries():
        bc = r["by_class"]
        DS, DL = bc["n_S"] * bc["z_S"], bc["n_L"] * bc["z_L"]
        pred = DS * DL / (DS + DL)
        v = bc["nc_SL"] / pred
        ax.plot(r["mix_ratio"], v, "o", color=COL[1], ms=3, alpha=0.6)
        rows.append(dict(kind="binary_SL", ratio=r["mix_ratio"], cls="SL", meas=bc["nc_SL"], pred_a2=pred, pred_blind=v))
    ax.axhline(1.0, color="0.5", lw=0.7)
    ax.set_xlim(0, 5.5); ax.set_ylim(0, 1.3)
    ax.set_xlabel(r"size ratio $d_L/d_S$"); ax.set_ylabel(r"$n_{SL}$ / configuration model")
    fig.tight_layout()
    save(fig, "G2_contact_laws_oos", rows,
         "(a) per-class rattler fraction pred (Poisson a^2) vs meas, grey = size-blind; (b) SL count / configuration model; "
         "for kind=binary_SL the pred_blind column holds the ratio")


if __name__ == "__main__":
    sys.exit(main())
