"""make_mechanism_figures.py — figures for the dilatancy-sign mechanism.

House style is `make_figures.py`'s and is imported from it rather than restated:
minimal text on the axes (the explanation lives in FIGURES_MECHANISM.md), linear axes
start at zero, log only above a decade, and every figure saves .pdf + .png + the
PLOTTED DATA as .csv.

    python make_mechanism_figures.py            # build
    python make_mechanism_figures.py --check    # build twice, compare md5

Figures
-------
M1  the regime map: which side of the drained/undrained line and which side of the
    dilatancy-sign line the lab's silica sits on, under BOTH readings of the drainage
    length. The two readings straddle the boundary, which is the point.
M2  the transient, from poro_impact.py: the pore-pressure sign flip across phi_c, and
    the strain clock that shuts the source off.
M3  phi_c(CV) from phase 8, with the crossing. Skipped with a printed reason if the
    sweep has not returned.
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

import poro_impact as PI          # noqa: E402
from make_figures import COL, FIG, md5, save  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "impact_stress")))

PHI = 0.44
PHI_C0 = 0.42
D10 = 800e-9
U0, D0 = 3.0, 2.5e-3
RHO = 1440.0
N_SHEARBAND = 10.0


def fig_M1() -> None:
    """(a) the SIGNED pore pressure across the critical state; (b) where Pi = 1 sits.

    Panel (a) is the mechanism: the pore pressure changes SIGN at phi = phi_c, and the
    network-modulus ceiling flattens both branches. Plotting |Pi| against |Delta phi|
    would hide exactly the thing the figure exists to show, which is what the first
    version of this figure did.

    Panel (b) is the open question: Pi depends on the drainage length SQUARED and the
    length is unmeasured, so the Pi = 1 line is crossed somewhere between the grain and
    the drop.
    """
    import pore_relaxation as PR
    from inertial_dilatant import ID, DilatantSystem

    dphis = np.linspace(-0.05, 0.05, 401)
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.7), constrained_layout=True)

    # --- (a) signed pore pressure: the two Darcy branches and the one ceiling --------
    # The capped answer is the inner envelope, and it turns out to be the CEILING under
    # BOTH drainage lengths -- so the section-6 ambiguity, which spans five decades in
    # Pi, does not touch the pore-pressure magnitude at all. Drawing only the capped
    # curves would have hidden that by making them coincide with no explanation.
    series = [(None, "ceiling", COL[0], "-"),
              (D0, r"Darcy, $\ell=D_0$", COL[1], "--"),
              (N_SHEARBAND * D10, r"Darcy, $\ell=10\,d_{32}$", COL[2], "-.")]
    for ell, label, c, ls in series:
        vals = []
        for d in dphis:
            if d == 0.0:
                vals.append(0.0)
                continue
            ps = PR.silica_poro(phi=PHI, phi_c=PHI - abs(d))
            v = (PR.dilatant_pressure_ceiling(ps.sys, ps.K_p) if ell is None
                 else PR.dilatant_pressure_darcy(ps.sys, L_d=ell,
                                                 gammadot=U0 / D0))
            vals.append(-np.sign(d) * v)
        axes[0].plot(dphis, vals, color=c, lw=1.4, ls=ls, label=label)
        for d, v in zip(dphis, vals):
            rows.append({"figure": "M1a", "ell_m": ell if ell else float("nan"),
                         "dphi": d, "p_f_Pa": v, "Pi": float("nan")})
    axes[0].axhline(0.0, color=COL[0], lw=0.6, ls=":")
    axes[0].axvline(0.0, color=COL[0], lw=0.6, ls=":")
    axes[0].set_yscale("symlog", linthresh=1.0)   # spans 6 decades AND changes sign
    axes[0].set_xlabel(r"$\phi-\phi_c$")
    axes[0].set_ylabel(r"$\Delta p_f$  [Pa]")
    axes[0].set_xlim(-0.05, 0.05)
    axes[0].legend(frameon=False, fontsize=6.5, loc="lower left")

    # --- (b) Pi against the drainage length ------------------------------------------
    s = DilatantSystem(name="x", phi=PHI, phi_c=PHI_C0, a=D10, rho=RHO)
    pi_d0 = 1.0 / ID(s, U0, D0)
    ells = np.logspace(np.log10(D10), np.log10(D0), 300)
    pis = pi_d0 * (ells / D0) ** 2
    axes[1].plot(ells / D0, pis, color=COL[0], lw=1.4)
    axes[1].axhline(1.0, color=COL[0], lw=0.8, ls="--")
    for ell, c in ((N_SHEARBAND * D10, COL[2]), (D0, COL[1])):
        axes[1].plot([ell / D0], [pi_d0 * (ell / D0) ** 2], "o", color=c, ms=5)
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")      # both span decades
    axes[1].set_xlabel(r"$\ell/D_0$")
    axes[1].set_ylabel(r"$\Pi$")
    for e, p in zip(ells, pis):
        rows.append({"figure": "M1b", "ell_m": e, "dphi": PHI - PHI_C0, "Pi": p,
                     "p_f_Pa": float("nan")})

    save(fig, "M1_regime_map", rows,
         "(a) signed pore pressure min(Darcy Eq.(4), network ceiling Eq.(16)), phi=0.44, "
         "d10=800nm, U0=3 m/s, D0=2.5mm; positive = OVERpressure. Both axes linear and "
         "both include zero. NOTE the ceiling uses K_p from impact_stress' EDL route "
         "whose P_conf is the WITHDRAWN si_fit sigma*, so the Pa scale is provisional; "
         "the SIGN and the shape are not. (b) Pi at dphi=0.02 against drainage length; "
         "dashed Pi=1; markers at ell=10 d32 and ell=D0. Both axes LOG, each spanning "
         "more than a decade.")


def fig_M2() -> None:
    """(a) the sign, in the UNDRAINED limit; (b) why the pack locks there and not else.

    ⚠ REBUILT 2026-09-08 after a referee found that the phi-evolution law was missing
    its elastic term. The first version of this figure showed the pore pressure peaking
    near gamma = 0.55 and decaying, and that decay was an artefact of the incomplete
    equation. Corrected, an undrained pack LOCKS -- the elastic compression cancels the
    Reynolds dilation exactly -- and the suction grows linearly instead.
    """
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.6), constrained_layout=True)
    s_dil = 1.0
    lam_und, lam_dr = 1e-6, 2.0

    for dphi0, lab, c in ((+0.02, r"$\phi>\phi_c$", COL[1]),
                          (0.0, r"$\phi=\phi_c$", COL[0]),
                          (-0.02, r"$\phi<\phi_c$", COL[2])):
        r = PI.solve(lam=lam_und, s_dil=s_dil, dphi0=dphi0, t_end=4.0, n_nodes=201)
        axes[0].plot(r.t, r.p_base, color=c, lw=1.4, label=lab)
        for t, p, d in zip(r.t, r.p_base, r.dphi):
            rows.append({"figure": "M2a", "lam": lam_und, "dphi0": dphi0, "t": t,
                         "p_base": p, "dphi": d})

    for lam, lab, c, ls in ((lam_und, r"undrained, $\Lambda=10^{-6}$", COL[1], "-"),
                            (lam_dr, r"drained, $\Lambda=2$", COL[3], "--")):
        r = PI.solve(lam=lam, s_dil=s_dil, dphi0=+0.02, t_end=4.0, n_nodes=201)
        axes[1].plot(r.t, r.dphi, color=c, lw=1.4, ls=ls, label=lab)
        if lam != lam_und:
            for t, p, d in zip(r.t, r.p_base, r.dphi):
                rows.append({"figure": "M2b", "lam": lam, "dphi0": 0.02, "t": t,
                             "p_base": p, "dphi": d})

    axes[0].axhline(0.0, color=COL[0], lw=0.6, ls=":")
    axes[1].axhline(0.0, color=COL[0], lw=0.6, ls=":")
    axes[0].set_xlabel(r"$\gamma$")
    axes[0].set_ylabel(r"$\hat p$ at the substrate")
    axes[1].set_xlabel(r"$\gamma$")
    axes[1].set_ylabel(r"$\phi-\phi_c$")
    for a in axes:
        a.set_xlim(0.0, 4.0)
    axes[0].legend(frameon=False, fontsize=7, loc="upper left")
    axes[1].legend(frameon=False, fontsize=7, loc="center right")
    save(fig, "M2_transient", rows,
         "poro_impact.solve, S=1, 201 nodes, with the CORRECTED phi-evolution law. "
         "(a) Lambda=1e-6, the lab's regime: pore pressure at the impermeable "
         "substrate, gauge, positive = OVERpressure. It grows linearly and without "
         "bound -- dilatancy hardening -- because an undrained pack cannot change "
         "volume. (b) phi - phi_c against strain: LOCKED when undrained, relaxing "
         "toward the critical state only when the pack can drain. All axes linear and "
         "all include zero.")


def fig_M3() -> None:
    """phi_c(CV) from phase 8, if it has returned."""
    import csv as _csv

    path = os.path.join(HERE, "out", "phase8_summary.csv")
    if not os.path.exists(path):
        print("  M3 SKIPPED: out/phase8_summary.csv does not exist yet "
              "(run analyse_phase8.py after the sweep)")
        return
    with open(path, encoding="utf-8") as fh:
        recs = list(_csv.DictReader(fh))
    if not recs:
        print("  M3 SKIPPED: out/phase8_summary.csv is empty")
        return

    tags = ("E1", "E2", "E3", "E4")
    fig, ax = plt.subplots(figsize=(3.4, 2.8))
    rows = []
    for tag, c in zip(tags, COL[1:5]):
        by_cv: dict[float, list[float]] = {}
        for r in recs:
            try:
                v = float(r.get(f"phi_c_{tag}", "nan"))
            except ValueError:
                continue
            if np.isfinite(v):
                by_cv.setdefault(float(r["cv"]), []).append(v)
        if len(by_cv) < 2:
            continue
        cvs = sorted(by_cv)
        m = [float(np.mean(by_cv[c_])) for c_ in cvs]
        sd = [float(np.std(by_cv[c_], ddof=1)) if len(by_cv[c_]) > 1 else 0.0
              for c_ in cvs]
        ax.errorbar(cvs, m, yerr=sd, color=c, lw=1.3, marker="o", ms=3, capsize=2,
                    label=tag)
        for c_, mm, ss in zip(cvs, m, sd):
            rows.append({"figure": "M3", "estimator": tag, "cv": c_, "phi_c": mm,
                         "sd": ss, "n": len(by_cv[c_])})
    # NO reference line at the lab's phi = 0.44. A frictional-sphere DEM's phi_c (~0.58)
    # is not commensurable with the 0.37-0.42 that impact_stress uses for cornstarch and
    # silica, so drawing 0.44 here would invite reading CV* off an intersection that
    # means nothing. Only the SHIFT transfers, and the shift is what the figure shows.
    ax.set_xlabel("CV")
    ax.set_ylabel(r"$\phi_c$")
    ax.set_xlim(0.0, None)
    ax.legend(frameon=False, fontsize=7)
    save(fig, "M3_phi_c_vs_CV", rows,
         "phase8_criticalstate.py, constant-volume shear. Four estimators, see "
         "analyse_phase8.py. x starts at zero; y does not, because phi_c is a fraction "
         "near 0.6 and the differences are ~0.06 -- the suppressed range is stated in "
         "FIGURES_MECHANISM.md instead. NO phi = 0.44 line is drawn: the DEM's absolute "
         "phi_c is not commensurable with the experiment's, only its SHIFT is.")


def fig_M4() -> None:
    """p_ss(phi) per CV -- the ladder shifting right, and the matched-phi collapse.

    This is the whole phase-8 result in one panel: three ladders, each the constant-
    volume grain pressure against packing fraction, displaced to the right as CV rises.
    The vertical gap where two ladders overlap in phi IS the headline number.
    """
    import csv as _csv
    import json as _json
    import glob as _glob

    recs = []
    for f in _glob.glob(os.path.join(HERE, "out", "phase8_prod2*.jsonl")):
        for line in open(f, encoding="utf-8"):
            r = _json.loads(line)
            if "error" not in r:
                recs.append(r)
    if not recs:
        print("  M4 SKIPPED: no phase8_prod2 records")
        return
    del _csv
    rows = []
    fig, ax = plt.subplots(figsize=(3.6, 2.9), constrained_layout=True)
    cvs = sorted({round(float(r["cv"]), 3) for r in recs})
    for cv, c in zip(cvs, COL[1:]):
        by_phi = {}
        for r in recs:
            if round(float(r["cv"]), 3) != cv:
                continue
            by_phi.setdefault(round(r["phi"], 3), []).append(r["p_ss_hat"])
        xs = sorted(by_phi)
        ys = [float(np.mean(by_phi[x])) for x in xs]
        es = [float(np.std(by_phi[x], ddof=1)) if len(by_phi[x]) > 1 else 0.0
              for x in xs]
        ax.errorbar(xs, ys, yerr=es, color=c, lw=1.3, marker="o", ms=3.5, capsize=2,
                    label=f"CV = {cv:.1f}")
        for x, y, e in zip(xs, ys, es):
            rows.append({"figure": "M4", "cv": cv, "phi": x, "p_ss_hat": y, "sd": e,
                         "n": len(by_phi[x])})
    ax.set_yscale("log")          # p_ss spans four decades
    ax.set_xlabel(r"$\phi$")
    ax.set_ylabel(r"$\hat p_{ss}$")
    ax.legend(frameon=False, fontsize=7)
    save(fig, "M4_pss_vs_phi", rows,
         "phase8_criticalstate.py, protocol cv-shear-v2, mu_p=0.3, lognormal, 3 seeds. "
         "Constant-volume steady-state grain pressure against packing fraction, one "
         "ladder per CV. y is LOG (four decades); x is linear and does NOT start at "
         "zero because phi is a fraction near 0.6 and the ladders span 0.09 -- the "
         "range shown is the range occupied. The vertical gap where the CV=0 and "
         "CV=0.5 ladders overlap in phi is the matched-phi result (5.4x to 43x).")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    os.makedirs(FIG, exist_ok=True)
    names = ["M1_regime_map", "M2_transient"]
    fig_M1()
    fig_M2()
    fig_M3()
    fig_M4()
    if os.path.exists(os.path.join(FIG, "M4_pss_vs_phi.png")):
        names.append("M4_pss_vs_phi")
    if os.path.exists(os.path.join(FIG, "M3_phi_c_vs_CV.png")):
        names.append("M3_phi_c_vs_CV")
    if args.check:
        def _h(n):
            return md5(os.path.join(FIG, f"{n}.png"))
        before = {n: _h(n) for n in names}
        fig_M1()
        fig_M2()
        fig_M3()
        fig_M4()
        bad = [n for n in names if _h(n) != before[n]]
        print("BYTE-IDENTICAL REGENERATION: "
              + ("PASSED" if not bad else f"FAILED for {bad}"))
        return 1 if bad else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
