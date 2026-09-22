"""make_activation_figures.py -- figures for the activation-shape work (2026-09-22, block 3).

House style from make_figures.py (imported, not restated). Explanations live in
derivations/jam_stress_2026-09-22.md and FIGURES_ACTIVATION.md, not on the axes.

    python make_activation_figures.py            # build
    python make_activation_figures.py --check    # build twice, compare md5

J1  THE THRESHOLD ACROSS THE WINDOW, two ways. (a) stress-controlled sigma_jam/sigma* =
    [ln(1/f_c)]^(-1/delta); (b) the solver's rate-driven gd_c eta_s/sigma*. chi = 0.80,
    delta in {0.7, 1.3169, 2.6}. (a) reverses order at f_c = 1/e (dotted, x = 0.926): the
    whole band chi < x < 0.926 has the opposite order to x > 0.926. (b) reverses only where
    f_c > ~0.95, i.e. within ~0.01 of the lower edge (delta 0.7 crosses 1.32 at x ~ 0.809 and
    2.6 at x ~ 0.804), and is NON-monotone in delta there. Arm B of jamwin (x = 0.85) sits
    outside that sliver, which is Amendment 1's prediction drawn.
J2  (once actshape lands) the measured activation curves fhat(sigma/kappa) by CV.
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

import actshape as A              # noqa: E402
from make_figures import COL, FIG, md5, save   # noqa: E402

CHI = 0.80
DELTAS = (0.7, 1.3169, 2.6)
HERE = os.path.dirname(os.path.abspath(__file__))


def j1_data():
    x = np.linspace(CHI + 1e-3, 1 - 1e-3, 120)
    rows, curves = [], {}
    for d in DELTAS:
        s = np.array([A.weibull_sigma_jam_ratio(A.f_crit(v, CHI), d) for v in x])
        r = np.array([A.rate_threshold(lambda f, d=d: np.log(1 / f) ** (-1 / d), v, CHI,
                                       n=20_001)[0] for v in x])
        curves[d] = (s, r)
        rows += [{"delta": d, "x": float(a), "stress_sjam_over_sstar": float(b),
                  "rate_gdc_etas_over_sstar": float(c)} for a, b, c in zip(x, s, r)]
    return x, curves, rows


def fig_j1():
    x, curves, rows = j1_data()
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.6))
    xe = 1 - (1 - CHI) * math.exp(-1)          # f_c = 1/e
    for k, d in enumerate(DELTAS):
        s, r = curves[d]
        ax[0].plot(x, s, color=COL[k + 1], lw=1.4, label=f"δ = {d:g}")
        ax[1].plot(x, r, color=COL[k + 1], lw=1.4)
    ax[0].axvline(xe, color=COL[0], lw=0.7, ls=":")
    ax[0].set_yscale("log")
    ax[1].set_yscale("log")
    ax[0].set_ylabel(r"$\sigma_{jam}/\sigma^*$")
    ax[1].set_ylabel(r"$\dot\gamma_c\,\eta_s/\sigma^*$")
    for a, t in zip(ax, ("(a) stress-controlled", "(b) rate-driven (solver)")):
        a.set_xlabel(r"$x=\phi/\phi_0$")
        a.set_xlim(CHI, 1.0)
        a.set_title(t, fontsize=9)
    ax[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    save(fig, "J1_threshold_across_window", rows,
         f"chi={CHI}, EPS_MIN=0.01, MP_EXPO=2; actshape.weibull_sigma_jam_ratio / "
         "rate_threshold; dotted line f_c = 1/e")


def j2_data():
    """delta_10/90 (derjaguin) per (mu, CV) at p 3e-4, N = 1000, W4-healthy cells; and the
    seed-pooled activation curves at mu 0.3."""
    import grade_actshape as G
    recs = G.load_records(os.path.join(HERE, "out"), "phase9_actshape_c*.jsonl")
    per, pooled = {}, {}
    for key, r in recs.items():
        if "error" in r or r["p_hat_target"] != 3e-4:
            continue
        if not (r["steady"] and r["servo_held"] and r["quasistatic"]):
            continue
        d = np.load(os.path.join(HERE, "out", "actshape_dump", key + ".npz"))
        snaps = G.snapshots(d)
        g = np.concatenate([A.g_values(s["fn"], s["r1"], s["r2"], s["p"], "derjaguin")
                            for s in snaps])
        per.setdefault((r["mu_p"], r["cv"]), []).append(A.delta_from_quantiles(g, 0.1, 0.9))
        if r["mu_p"] == 0.3:
            pooled.setdefault(r["cv"], []).append(g)
    return per, {cv: np.concatenate(v) for cv, v in pooled.items()}


def fig_j2():
    per, pooled = j2_data()
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.6))
    rows = []
    for k, mu in enumerate((0.0, 0.3)):
        for cv in (0.0, 0.2, 0.5):
            v = per.get((mu, cv))
            if not v:
                continue
            m, sd = float(np.mean(v)), float(np.std(v, ddof=1)) if len(v) > 1 else 0.0
            flagged = cv == 0.5                      # finite-size flagged (Amendment 1)
            ax[0].errorbar(cv, m, yerr=sd, fmt="o", color=COL[k + 1], capsize=2,
                           mfc="white" if flagged else COL[k + 1],
                           label=f"μp = {mu:g}" if cv == 0.0 else None)
            rows.append({"mu_p": mu, "cv": cv, "delta_mean": m, "delta_sd": sd, "n": len(v),
                         "finite_size_flagged": int(flagged)})
    ax[0].axhspan(1.3169, 1.3852, color=COL[0], alpha=0.12, lw=0)
    ax[0].set_ylim(0, 1.7)
    ax[0].set_xlim(-0.05, 0.55)
    ax[0].set_xlabel("CV")
    ax[0].set_ylabel(r"$\delta_{10/90}$")
    ax[0].legend(frameon=False, fontsize=8, loc="lower left")
    s = np.logspace(-2, 1.5, 351)
    for k, cv in enumerate((0.0, 0.2, 0.5)):
        if cv not in pooled:
            continue
        f = A.activation_curve(pooled[cv], s)
        ax[1].plot(s, f, color=COL[k + 1], lw=1.4, ls="--" if cv == 0.5 else "-",
                   label=f"CV = {cv:g}")
        rows += [{"mu_p": 0.3, "cv": cv, "sigma_over_kappa": float(a), "fhat": float(b)}
                 for a, b in zip(s[::10], f[::10])]
    ax[1].set_xscale("log")
    ax[1].set_ylim(0, 1)
    ax[1].set_xlabel(r"$\sigma/\kappa$  ($d_{10}$ units)")
    ax[1].set_ylabel(r"$\hat f$")
    ax[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    keys = sorted({k for r in rows for k in r})
    rows = [{k: r.get(k, "") for k in keys} for r in rows]
    save(fig, "J2_activation_by_cv", rows,
         "actshape p=3e-4 N=1000 derjaguin; hollow/dashed = CV 0.5 finite-size flagged; "
         "grey band = solver DELTA cs_fit 1.3169 .. si_fit 1.3852")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    fig_j1()
    fig_j2()
    if a.check:
        h1 = md5(os.path.join(FIG, "J1_threshold_across_window.png"))
        fig_j1()
        h2 = md5(os.path.join(FIG, "J1_threshold_across_window.png"))
        print("J1 png deterministic:", h1 == h2, h1[:12])
    return 0


if __name__ == "__main__":
    sys.exit(main())
