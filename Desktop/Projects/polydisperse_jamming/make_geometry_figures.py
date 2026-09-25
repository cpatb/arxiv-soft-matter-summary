"""make_geometry_figures.py -- figures for the 2026-09-24 geometric results.

G1 (a) sigma*(CV)/sigma*(0) at fixed d10 split into its exact factors (Love-Weber identity,
       derivations/sigma_star_geometry_2026-09-24.md eq. 7): Z phi, <a_i a_j>_c/<a^3>, C, 1/q, and the
       measured ratio; mu_p 0.3 at p_hat 3e-4 and 1e-3 (the p_hat 1e-4 cells unjam intermittently).
   (b) per-grain rattler probability against d/<d> (out/zbysize, isotropic, mu_p 0.3) with the Poisson
       predictions for lambda ~ a^2 and closure C (derivations/rattler_closure_2026-09-24.md).
Numbers are parsed from out/sigma_star_geometry.txt and recomputed from out/zbysize/*.npz.
Explanations in FIGURES_GEOMETRY.md.

    python make_geometry_figures.py
"""
import glob
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402
from scipy.stats import poisson   # noqa: E402

from make_figures import COL, save   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")


def parse_factors():
    txt = open(os.path.join(OUT, "sigma_star_geometry.txt"), encoding="utf-8").read()
    rows = {}
    for m in re.finditer(r"^\s+(0\.3) (\de-0\d)\s+(0\.[25])\s+(\d+)\s+(\S+)\s+=\s+(\S+)\s+x\s+(\S+)\s+x\s+(\S+)", txt, re.M):
        mu, p, cv, n, rm, g, c, iq = m.groups()
        rows[(p, cv)] = dict(R=float(rm), G=float(g), C=float(c), iq=float(iq), n=int(n))
    for m in re.finditer(r"^\s+(0\.3) (\de-0\d)\s+(0\.[25])\s+(\d+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\|", txt, re.M):
        mu, p, cv, n, zp, aa = m.groups()
        if (p, cv) in rows:
            rows[(p, cv)].update(Zphi=float(zp), aa=float(aa))
    return rows


def fig_g1():
    rows = parse_factors()
    fig, axs = plt.subplots(1, 2, figsize=(6.8, 2.7))
    ax = axs[0]
    keys = [("3e-04", "0.2"), ("3e-04", "0.5"), ("1e-03", "0.2"), ("1e-03", "0.5")]
    names = [("Zphi", r"$Z\phi$"), ("aa", r"$\langle a_ia_j\rangle_c/\langle a^3\rangle$"),
             ("C", r"$C$"), ("iq", r"$1/q$"), ("R", r"$\sigma^*$ (measured)")]
    w = 0.16
    data = []
    for j, (k, lab) in enumerate(names):
        vals = [rows[kk][k] for kk in keys]
        ax.bar(np.arange(len(keys)) + (j - 2) * w, vals, w, color=COL[j + 1] if j < 4 else COL[0], label=lab)
        for kk, v in zip(keys, vals):
            data.append(dict(p_hat=kk[0], cv=kk[1], factor=k, ratio_vs_cv0=v))
    ax.axhline(1.0, color="0.5", lw=0.6)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([f"CV {c}\n$\\hat p$ {p}" for p, c in keys], fontsize=7)
    ax.set_ylim(0, 1.2)
    ax.set_ylabel("ratio to CV 0")
    ax.legend(frameon=False, fontsize=6, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    ax = axs[1]
    for i, f in enumerate(sorted(glob.glob(os.path.join(OUT, "zbysize", "zb_cv*_mu0.3_*.npz")))):
        z = np.load(f)
        d, k, rat = z["d"], z["k"].astype(float), z["rattler"]
        a = d / 2
        Z = k.mean()
        edges = np.unique(np.quantile(d, np.linspace(0, 1, 9)))
        idx = np.clip(np.searchsorted(edges, d, side="right") - 1, 0, len(edges) - 2)
        m1, m2 = a.mean(), (a * a).mean()
        sC = a * (a * m1 + m2); sC /= sC.mean()
        s2 = a * a / (a * a).mean()
        x, y, pC, p2 = [], [], [], []
        for b in range(len(edges) - 1):
            m = idx == b
            x.append(d[m].mean() / d.mean()); y.append(rat[m].mean())
            pC.append(np.mean(poisson.cdf(1, Z * sC[m]))); p2.append(np.mean(poisson.cdf(1, Z * s2[m])))
        cv = re.search(r"zb_cv([\d.]+)_", f).group(1)
        col = COL[1 + i]
        ax.plot(x, y, "o", color=col, ms=4, label=f"CV {cv}")
        ax.plot(x, p2, "-", color=col, lw=1)
        ax.plot(x, pC, ":", color=col, lw=1)
        for xi, yi, a2i, cci in zip(x, y, p2, pC):
            data.append(dict(cv=cv, d_over_mean=xi, p_rattler=yi, pred_a2=a2i, pred_C=cci))
    ax.set_xlim(0, None); ax.set_ylim(0, 1)
    ax.set_xlabel(r"$d/\langle d\rangle$"); ax.set_ylabel("P(rattler)")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    cols = ["p_hat", "cv", "factor", "ratio_vs_cv0", "d_over_mean", "p_rattler", "pred_a2", "pred_C"]
    data = [{c: r.get(c, "") for c in cols} for r in data]
    save(fig, "G1_sigma_star_geometry", data,
         "(a) factor ratios CV/CV0 (mu 0.3) from out/sigma_star_geometry.txt; (b) per-grain rattler prob. and Poisson preds (solid a^2, dotted closure C)")


if __name__ == "__main__":
    fig_g1()
    print("G1 written")
