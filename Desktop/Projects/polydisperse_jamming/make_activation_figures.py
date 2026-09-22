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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    fig_j1()
    if a.check:
        h1 = md5(os.path.join(FIG, "J1_threshold_across_window.png"))
        fig_j1()
        h2 = md5(os.path.join(FIG, "J1_threshold_across_window.png"))
        print("J1 png deterministic:", h1 == h2, h1[:12])
    return 0


if __name__ == "__main__":
    sys.exit(main())
