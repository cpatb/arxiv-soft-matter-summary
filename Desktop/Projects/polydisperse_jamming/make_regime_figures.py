"""make_regime_figures.py -- figure M1, the ball-vs-hybrid plane (REGIME_MAP_2026-09-26 Gate C).

Drainage Deborah number De against dilatant-strength ratio S, both from the project's poroelastic model
(derivations/ball_hybrid_scaling.groups; impact_stress/pore_relaxation). Paths (each a knob at fixed everything else):
  silica 800 nm, phi 0.44, lab impact: salt, lambda_D 30 -> 1 nm (K_p = sigma* a / 6 lambda_D, sigma* fixed);
  silica at fixed chemistry: size 0.4 -> 20 um (Derjaguin convention [D], as ball_hybrid_scaling section 3);
  cornstarch phi 0.42, [M26] impact: size fraction 5 -> 20 um (Hertz K_p recomputed);
  cornstarch 10 um: solvent viscosity 1e-3 -> 1 Pa s (tau_p with mu_f passed explicitly; groups() is inert to it).
A scaling estimate, not a measurement. Explanations in FIGURES_REGIME.md; house style from make_figures.py.

    python make_regime_figures.py [--check]
"""
import argparse
import dataclasses
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

from make_figures import COL, FIG, md5, save   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "derivations"))
import ball_hybrid_scaling as bhs  # noqa: E402

pr, idm = bhs.pr, bhs.idm


def silica_salt():
    base = pr.silica_poro(phi=0.44)
    out = []
    for lam in (30e-9, 10e-9, 3e-9, 1e-9):
        K = pr.network_modulus_edl(bhs.SIGMA_STAR_SILICA[1], base.sys.a, lam)
        g = bhs.groups(pr.PoroSystem(sys=base.sys, K_p=K, K_p_low=K, K_p_high=K, K_p_origin="salt"), **bhs.LAB["silica"])
        out.append((g["De"], g["S"], f"lambda_D {lam*1e9:g} nm"))
    return out


def silica_size():
    base = pr.silica_poro(phi=0.44)
    a0, s0 = base.sys.a, bhs.SIGMA_STAR_SILICA[1]
    out = []
    for a in (0.4e-6, 0.8e-6, 2e-6, 5e-6, 10e-6, 20e-6):
        sysd = idm.DilatantSystem(name="ladder", phi=base.sys.phi, phi_c=base.sys.phi_c, a=a, rho=base.sys.rho)
        K = pr.network_modulus_edl(s0 * (a / a0) ** -1, a)
        g = bhs.groups(pr.PoroSystem(sys=sysd, K_p=K, K_p_low=K, K_p_high=K, K_p_origin="size"), **bhs.LAB["silica"])
        out.append((g["De"], g["S"], f"a {a*1e6:g} um"))
    return out


def cs_poro(a, mu=None):
    b = pr.cornstarch_poro(phi=0.42)
    s = dataclasses.replace(b.sys, a=a) if mu is None else dataclasses.replace(b.sys, a=a, mu_w=mu)
    k = pr.network_modulus_affine(pr.contact_stiffness_hertz(pr.E_STARCH_PA, 0.3, s.a, s.phi, pr.P_CONF_CORNSTARCH),
                                  s.a, s.phi)
    return pr.PoroSystem(sys=s, K_p=k, K_p_low=k, K_p_high=k, K_p_origin="cs ladder")


def cs_size():
    out = []
    for a in (5e-6, 7e-6, 10e-6, 14e-6, 20e-6):
        g = bhs.groups(cs_poro(a), **bhs.LAB["cornstarch"])
        out.append((g["De"], g["S"], f"a {a*1e6:g} um"))
    return out


def cs_visc():
    out = []
    op = bhs.LAB["cornstarch"]
    for mu in (1e-3, 1e-2, 1e-1, 1.0):
        p = cs_poro(10e-6, mu)
        s = p.sys
        Ld = pr.drainage_lengths(op["D0"], s.a, s.phi)["contact"]
        De = pr.drainage_number(pr.pore_relaxation_time(Ld, s.kappa, p.K_p, s.phi, mu_f=mu), op["D0"], op["U0"])
        out.append((De, bhs.groups(p, **op)["S"], f"mu_f {mu:g}"))
    return out


def build():
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    rows = []
    for name, pts, c, m in (("silica: salt", silica_salt(), COL[1], "o"),
                            ("silica: size (D)", silica_size(), COL[4], "s"),
                            ("cornstarch: size fraction", cs_size(), COL[2], "^"),
                            ("cornstarch: solvent viscosity", cs_visc(), COL[3], "v")):
        de, s_, lab = zip(*pts)
        ax.plot(s_, de, m + "-", color=c, ms=4, lw=0.9)
        rows += [{"path": name, "point": l_, "De": d, "S": sv} for d, sv, l_ in pts]
    ax.axhline(1.0, color="0.6", lw=0.6, ls="--")
    ax.axvline(1.0, color="0.6", lw=0.6, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"dilatant-strength ratio $S = K_p\Delta\phi/(\phi\rho U_0^2)$")
    ax.set_ylabel(r"drainage Deborah number $De$")
    fig.tight_layout()
    save(fig, "M1_ball_hybrid_plane", rows,
         "make_regime_figures.py; ball_hybrid_scaling.groups + pore_relaxation; scaling estimate")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    build()
    if a.check:
        h1 = md5(os.path.join(FIG, "M1_ball_hybrid_plane.png"))
        build()
        h2 = md5(os.path.join(FIG, "M1_ball_hybrid_plane.png"))
        print(f"M1 png md5 {h1} / {h2}: {'IDENTICAL' if h1 == h2 else 'DIFFERS'}")


if __name__ == "__main__":
    main()
