"""make_kmer_figures.py -- the singlet/doublet figure (PREREG_KMER).

K1 (a) chi(w2) at p_hat 3e-5 for mu_m = 3 and 0.3, seed mean +- sd, with the H_mix curve
       drawn from the campaign's own endpoints; (b) sigma*(w2)/sigma*(0) and delta(w2) at
       mu_p 0.3, from the per-contact dumps (R* = 1/4 for every contact).
Explanations in FIGURES_KMER.md; house style from make_figures.py.

    python make_kmer_figures.py [--check]
"""
import argparse
import glob
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

import actshape as A              # noqa: E402
from make_figures import COL, FIG, md5, save   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
WS = (0.0, 0.25, 0.5, 0.75, 1.0)
RUNG = 3e-5


def load():
    recs = []
    for f in glob.glob(os.path.join(OUT, "shape", "phase2_kmer_*.jsonl")):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            if r.get("failed") or not r.get("converged") or r.get("triplet_w", 0) != 0:
                continue
            recs.append(r)
    return recs


def fig_k1():
    recs = load()
    by = {}
    for r in recs:
        by[(r["doublet_w"], r["mu_p"], r["p_hat_target"], r["seed"])] = r
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.6))
    rows = []
    for k, mm in enumerate((3.0, 0.3)):
        pts = []
        for w in WS:
            c = [by[(w, mm, RUNG, s)]["phi"] / by[(w, 0.0, RUNG, s)]["phi"] for s in (1, 2, 3)
                 if (w, mm, RUNG, s) in by and (w, 0.0, RUNG, s) in by]
            if c:
                pts.append((w, float(np.mean(c)), float(np.std(c, ddof=1)) if len(c) > 1 else 0.0,
                            len(c)))
        if not pts:
            continue
        ww, mean, sd, n = zip(*pts)
        # mu_m 3 clump cells fail the post-hoc overlap gate K4 (PREREG_KMER Amendment 2):
        # drawn HOLLOW and labelled, never as clean data
        flagged = mm == 3.0
        ax[0].errorbar(ww, mean, yerr=sd, fmt="o", color=COL[k + 1], capsize=2,
                       mfc="white" if flagged else COL[k + 1],
                       label=f"μm = {mm:g}" + (" (K4-flagged)" if flagged else ""))
        # H_mix from the campaign's own endpoint means
        e = {}
        for w in (0.0, 1.0):
            for m_ in (mm, 0.0):
                v = [by[(w, m_, RUNG, s)]["phi"] for s in (1, 2, 3) if (w, m_, RUNG, s) in by]
                e[(w, m_)] = np.mean(v) if v else np.nan
        wg = np.linspace(0, 1, 101)
        hm = (1 / ((1 - wg) / e[(0.0, mm)] + wg / e[(1.0, mm)])) / \
             (1 / ((1 - wg) / e[(0.0, 0.0)] + wg / e[(1.0, 0.0)]))
        ax[0].plot(wg, hm, color=COL[k + 1], lw=1.0, ls="--")
        rows += [{"panel": "a", "quantity": "chi", "mu_m": mm, "w2": a, "value": b, "sd": c, "n": d}
                 for a, b, c, d in pts]
    ax[0].set_xlabel(r"doublet volume fraction $w_2$")
    ax[0].set_ylabel(r"$\chi=\phi_m/\phi_0$")
    ax[0].set_xlim(-0.05, 1.05)
    ax[0].set_ylim(0, 1)
    ax[0].legend(frameon=False, fontsize=8, loc="lower left")
    # (b) onset and sharpness from the dumps, mu_p 0.3
    sst, dl = {}, {}
    for w in WS:
        for s in (1, 2, 3):
            r = by.get((w, 0.3, RUNG, s))
            if r is None:
                continue
            p = os.path.join(OUT, "kmer_dump", r["key"] + ".npz")
            if not os.path.exists(p):
                continue
            dd = np.load(p)
            fn = dd["fn"][dd["fn"] > 0]
            g = fn / (0.25 * float(dd["p_hat"]))
            sst.setdefault(w, {})[s] = A.sigma_star_over_kappa(g)
            dl.setdefault(w, {})[s] = A.delta_from_quantiles(g, 0.1, 0.9)
    ax2 = ax[1].twinx()
    for w in WS:
        if w in sst and 0.0 in sst:
            rr = [sst[w][s] / sst[0.0][s] for s in sst[w] if s in sst[0.0]]
            if rr:
                ax[1].errorbar(w, np.mean(rr), yerr=np.std(rr, ddof=1) if len(rr) > 1 else 0,
                               fmt="s", color=COL[1], capsize=2)
                rows.append({"panel": "b", "quantity": "sigma_star_ratio", "mu_m": 0.3, "w2": w, "value": np.mean(rr),
                             "sd": 0.0, "n": len(rr)})
        if w in dl:
            v = list(dl[w].values())
            ax2.errorbar(w, np.mean(v), yerr=np.std(v, ddof=1) if len(v) > 1 else 0, fmt="^",
                         color=COL[2], capsize=2)
    ax[1].axhspan(0.55, 0.95, color=COL[1], alpha=0.10, lw=0)
    ax[1].set_ylim(0, 1.2)
    ax2.set_ylim(0, 1.6)
    ax[1].set_xlim(-0.05, 1.05)
    ax[1].set_xlabel(r"$w_2$")
    ax[1].set_ylabel(r"$\sigma^*(w_2)/\sigma^*(0)$", color=COL[1])
    ax2.set_ylabel(r"$\delta_{10/90}$", color=COL[2])
    fig.tight_layout()
    save(fig, "K1_doublet_window", rows,
         "PREREG_KMER p_hat 3e-5; (a) chi seed-paired, dashed = H_mix of the campaign's own "
         "endpoints; (b) mu_p 0.3 dumps, shaded = registered sigma* band [0.55, 0.95]")


def fig_k2():
    """K2: the DNS drop at fixed x = 0.85 across the chi(w2) ladder (PREREG_KMERDROP):
    jammed volume J(f_hat > f_c) and the maximum radial radius of gyration (continuous;
    beta_max is grid-quantised)."""
    labs = (("w0.00", 0.0), ("w0.25", 0.25), ("w0.50", 0.5), ("w0.75", 0.75), ("w1.00", 1.0))
    rows, ws, J, R = [], [], [], []
    for lab, w in labs:
        p = os.path.join(OUT, "kmerdrop", lab, "spread.dat")
        if not os.path.exists(p):
            continue
        a = np.loadtxt(p, comments="#", ndmin=2)
        ws.append(w)
        J.append(float(a[:, 7].max()))
        R.append(float(a[:, 11].max()))
        rows.append({"w2": w, "J_max": J[-1], "Rg_r_max_over_D0": R[-1],
                     "beta_max": float(a[:, 2].max())})
    if not ws:
        return
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.6))
    ax[0].plot(ws, J, "o-", color=COL[2])
    ax[0].set_ylim(0, 1)
    ax[0].set_xlabel(r"$w_2$")
    ax[0].set_ylabel(r"jammed volume fraction")
    ax[1].plot(ws, R, "s-", color=COL[1])
    ax[1].set_ylim(0, 0.5)
    ax[1].set_xlabel(r"$w_2$")
    ax[1].set_ylabel(r"max $R_{g,r}/D_0$")
    for a_ in ax:
        a_.set_xlim(-0.05, 1.05)
    fig.tight_layout()
    save(fig, "K2_doublet_drop", rows,
         "PREREG_KMERDROP: x = 0.85, sigma* 300 Pa, delta 1.3169, L9; chi(w2) from the DEM")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    fig_k1()
    fig_k2()
    if a.check:
        h1 = md5(os.path.join(FIG, "K1_doublet_window.png"))
        fig_k1()
        print("K1 deterministic:", h1 == md5(os.path.join(FIG, "K1_doublet_window.png")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
