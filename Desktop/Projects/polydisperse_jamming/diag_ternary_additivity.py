"""diag_ternary_additivity.py -- can CALIBRATION explain the ternary residual, or only
NON-ADDITIVITY?

POST-HOC. Written after `grade_ternary.py` returned FAILURES ['T1', 'T2'] and is NOT a
gate: nothing here enters FAILURES, no bar is moved, and the T1 verdict stands on the
pre-registered bars whatever this file says. It exists because the pre-registration named
TWO explanations for the continuous-family under-prediction and the ternary was chosen to
separate them:

    (i)  CALIBRATION -- a(u) is simply wrong at the u values involved;
    (ii) ADDITIVITY  -- the pair sum is the wrong way to combine them.

T1 failing 6/6 positive is consistent with BOTH. This file asks the question that is not:

    hold the model's STRUCTURE fixed and ask what interaction value a(u) the ADDITIVE
    model would need to reproduce each measured ternary. Then compare that value with
    the one the BINARIES measured at the same u, WITH ITS PROFILE INTERVAL.

If the required a(u) lies inside the binary interval, calibration can explain the
ternary and additivity is not implicated. If it lies far outside -- or outside [0, 1],
where no interaction value exists at all -- then no recalibration of that interaction can
absorb the residual, and additivity is the remaining explanation.

TWO HONEST LIMITS OF THIS DIAGNOSTIC, stated before its numbers:
  * It inverts ONE interaction at a time. The residual could in principle be absorbed by
    a joint move of both, so a single-parameter inversion landing outside its interval
    is evidence, not proof. Both single-parameter inversions are reported.
  * The interval is the binary profile interval, which is a statistical interval on a
    fit; the ternary residual is a difference of two protocol means. They are not the
    same kind of uncertainty and the comparison is reported as such.

Run:  python diag_ternary_additivity.py     (writes out/diag_ternary_additivity.txt)
"""

from __future__ import annotations

import os

import numpy as np
from scipy.optimize import brentq

import grade_mixture_ext as gx
import grade_ternary as gt
import mixture_closure as mc
import mixture_model as mm
from phase12_ternary import COMPOSITIONS, make_ternary

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = mm.PHI_SCALE
MU = 0.3
LINES: list[str] = []


def say(s: str = "") -> None:
    LINES.append(s)
    print(s)


def main() -> int:
    meas3 = mm.measured(mu_p=MU)
    beta = meas3[("lognormal", 0.0)][0]
    pts23 = mc.binary_points(mu_p=MU)
    allpts = gx.ext_points(beta)
    cand, _, _ = gx.candidates(beta, pts23)
    paJ, pbJ = cand["J"](0)
    tab, fa, fb = gx.tabulated_closure(beta, allpts)

    say("=== diag_ternary_additivity.py -- POST-HOC, not a gate ===")
    say("    T1 stands on the pre-registered bars whatever this file says.")
    say()

    # measured ternary phi, through the grader's own loader and admission rule
    g, _, _ = gt.load()
    tphi = {}
    for comp in COMPOSITIONS:
        key = tuple(round(v, 4) for v in comp)
        ok = [c for c in g.get(key, []) if gt.admissible(c)[0]]
        v = np.array([c["phi"] for c in ok])
        tphi[comp] = (float(v.mean()), float(v.std(ddof=1)), v.size)

    # --- the binary interaction values with their profile intervals -------------------
    say("=== what the BINARIES measured, with their profile intervals ===")
    say(f"   {'r':>5}{'u':>8}{'a':>9}{'a interval':>20}{'b':>9}{'b interval':>20}")
    binab = {}
    for r in (2.0, 4.0):
        sub = [p for p in allpts if abs(p[0] - r) < 1e-9 and 0.0 < p[1] < 1.0]
        a_, b_, _ = mc.fit_ab_direct(sub, beta)
        alo, ahi = mc.profile_interval(sub, beta, 0)
        blo, bhi = mc.profile_interval(sub, beta, 1)
        binab[1.0 / r] = (a_, alo, ahi, b_, blo, bhi)
        say(f"   {r:5.1f}{1 / r:8.3f}{a_:9.4f}   [{alo:.4f}, {ahi:.4f}]"
            f"{b_:9.4f}   [{blo:.4f}, {bhi:.4f}]   n = {len(sub)}")
    say()

    # --- the inversion -----------------------------------------------------------------
    def phi_with(comp, u_target, a_override=None, b_override=None):
        """The ADDITIVE model with the interaction at ONE u replaced."""
        p = make_ternary(comp)

        def A(u):
            u = np.asarray(u, float)
            out = np.asarray(fa(u), float).copy()
            if a_override is not None:
                out = np.where(np.abs(u - u_target) < 1e-9, a_override, out)
            return out

        def B(u):
            u = np.asarray(u, float)
            out = np.asarray(fb(u), float).copy()
            if b_override is not None:
                out = np.where(np.abs(u - u_target) < 1e-9, b_override, out)
            return out

        return mc.packing_density_ab(p.diameters, p.volume_fractions, beta, A, B)

    say("=== INVERSION: what a(u) would the ADDITIVE model need? ===")
    say("   One interaction moved at a time; everything else held at its measured value.")
    say("   a is an interaction function and is only DEFINED on [0, 1]: a = 0 is perfect")
    say("   filling (the fines cost the coarse skeleton nothing), a = 1 is full")
    say("   interference. A required value below 0 means NO admissible interaction")
    say("   reproduces the measurement -- the additive model cannot get there at all.")
    say()
    say(f"   {'y_L':>6}{'y_M':>6}{'y_S':>6}{'phi meas':>10}{'phi ADD':>9}"
        f"{'req a(1/4)':>12}{'in [lo,hi]?':>13}{'req a(1/2)':>12}{'in [lo,hi]?':>13}")
    rows = []
    for comp in COMPOSITIONS:
        m, sd, n = tphi[comp]
        p = make_ternary(comp)
        padd = mc.packing_density_tabulated(p.diameters, p.volume_fractions, beta, fa, fb)
        out = {}
        for u in (0.25, 0.5):
            lo_, hi_ = -3.0, 1.0
            f = lambda a: phi_with(comp, u, a_override=a) - m           # noqa: E731
            try:
                if f(lo_) * f(hi_) > 0:
                    # no root in the admissible-plus-slack window: report the endpoint
                    # behaviour rather than a number a root finder invented
                    out[u] = (float("nan"),
                              f"no root in [{lo_}, {hi_}]: f({lo_}) {f(lo_):+.4f}, "
                              f"f({hi_}) {f(hi_):+.4f}")
                else:
                    ra = float(brentq(f, lo_, hi_, xtol=1e-10))
                    a_, alo, ahi, *_ = binab[u]
                    inside = alo <= ra <= ahi
                    out[u] = (ra, ("INSIDE" if inside else
                                   ("outside, ADMISSIBLE" if 0 <= ra <= 1 else
                                    "outside AND INADMISSIBLE (a < 0)")))
            except Exception as exc:                       # noqa: BLE001
                out[u] = (float("nan"), repr(exc)[:40])
        rows.append((comp, m, padd, out))

        def fmt(t):
            v, why = t
            return ("  NO ROOT" if not np.isfinite(v) else f"{v:8.4f}")

        def tag(t):
            v, why = t
            if not np.isfinite(v):
                return "no admissible a reproduces it"
            return why

        say(f"   {comp[0]:6.2f}{comp[1]:6.2f}{comp[2]:6.2f}{m:10.4f}{padd:9.4f}"
            f"  {fmt(out[0.25])}  {fmt(out[0.5])}   "
            f"u=1/4: {tag(out[0.25])} | u=1/2: {tag(out[0.5])}")
    say()
    say("   🔴 THE 'NO ROOT' ENTRIES ARE NOT EVIDENCE, AND AN EARLIER VERSION OF THIS FILE")
    say("      CALLED THEM 'the stronger half of the result'. THAT WAS BACKWARDS.")
    say("      A referee pass (2026-09-14) asked what phi actually does as a(1/4) moves.")
    say("      Measured below: it does NOTHING. phi is BIT-IDENTICAL across the whole")
    say("      admissible range, so 'no root' says a(1/4) is NOT IDENTIFIABLE from these")
    say("      compositions -- a statement about the branch structure, not about filling.")
    say("      A degenerate parameter cannot testify either way, and the honest evidence")
    say("      is the a(1/2) half alone, which IS responsive.")
    say()
    say("   sweep of a(1/4) at (0.60, 0.25, 0.15), everything else held:")
    sweep = [-3.0, -1.0, 0.0, 0.2, 0.3010, 0.4, 0.5, 0.6, 0.8, 1.0]
    base = phi_with(COMPOSITIONS[0], 0.25, a_override=sweep[0])
    for a in sweep:
        v = phi_with(COMPOSITIONS[0], 0.25, a_override=a)
        say(f"      a(1/4) = {a:+7.4f} -> phi = {v:.12f}   change from a = -3: "
            f"{v - base:+.3e}")
    say("   phi does not move until a ~ 0.6, where gamma_L finally becomes the binding")
    say("   branch. Below that the L-S pair term is simply absent from the minimum.")
    say()
    say("   the same sweep for a(1/2), which enters BOTH the L-M and M-S terms:")
    for a in (0.3, 0.5, 0.6, 0.7, 0.7854, 0.9):
        v = phi_with(COMPOSITIONS[0], 0.5, a_override=a)
        say(f"      a(1/2) = {a:+7.4f} -> phi = {v:.12f}")
    say("   THIS one is responsive over its whole range, so its root is a real")
    say("   constraint. It is the only half of the inversion that carries evidence.")
    say()
    a14, alo14, ahi14 = binab[0.25][0], binab[0.25][1], binab[0.25][2]
    a12, alo12, ahi12 = binab[0.5][0], binab[0.5][1], binab[0.5][2]
    say(f"   binary-measured a(1/4) = {a14:.4f} [{alo14:.4f}, {ahi14:.4f}]; "
        f"a(1/2) = {a12:.4f} [{alo12:.4f}, {ahi12:.4f}]")
    r14 = [r[3][0.25][0] for r in rows if np.isfinite(r[3][0.25][0])]
    r12 = [r[3][0.5][0] for r in rows if np.isfinite(r[3][0.5][0])]
    if r14:
        say(f"   required a(1/4): {min(r14):.4f} to {max(r14):.4f} over the six "
            f"compositions; {sum(1 for v in r14 if alo14 <= v <= ahi14)}/{len(r14)} inside "
            f"the binary interval, {sum(1 for v in r14 if v < 0)}/{len(r14)} INADMISSIBLE")
    if r12:
        say(f"   required a(1/2): {min(r12):.4f} to {max(r12):.4f}; "
            f"{sum(1 for v in r12 if alo12 <= v <= ahi12)}/{len(r12)} inside the binary "
            f"interval, {sum(1 for v in r12 if v < 0)}/{len(r12)} INADMISSIBLE")
    say()

    # --- THE BRACKET, which a referee pass identified as the cleanest statement --------
    say("=== BRACKET: the measurement lies BETWEEN the additive and sequential rules ===")
    say("   Raised by the physics-referee pass, 2026-09-14. Comparing rms residuals says")
    say("   which of two candidates is nearer; it does NOT say the winner is right. The")
    say("   SIGN structure says more, and it says something neither rms does:")
    say("   ADD is systematically BELOW the measurement and SEQd43 systematically ABOVE.")
    say(f"   {'y_L':>6}{'y_M':>6}{'y_S':>6}{'ADD':>9}{'meas':>9}{'SEQd32':>9}{'SEQd43':>9}"
        f"   ADD < meas < SEQ?")
    import grade_ternary as gt2
    nb32 = nb43 = 0
    for comp in COMPOSITIONS:
        m, _, _ = tphi[comp]
        p = make_ternary(comp)
        add = mc.packing_density_tabulated(p.diameters, p.volume_fractions, beta, fa, fb)
        s32 = gt2.sequential(comp, "d32", "TAB", beta, paJ, pbJ, fa, fb)
        s43 = gt2.sequential(comp, "d43", "TAB", beta, paJ, pbJ, fa, fb)
        b32, b43 = add < m < s32, add < m < s43
        nb32 += int(b32); nb43 += int(b43)
        say(f"   {comp[0]:6.2f}{comp[1]:6.2f}{comp[2]:6.2f}{add:9.4f}{m:9.4f}"
            f"{s32:9.4f}{s43:9.4f}   d32 {'yes' if b32 else 'NO ':>3}  "
            f"d43 {'yes' if b43 else 'NO ':>3}")
    say(f"   bracketed by (ADD, SEQd32): {nb32}/6;  by (ADD, SEQd43): {nb43}/6")
    say("   CONSEQUENCE for how this block is written up: a SMALLER rms for SEQd32 is")
    say("   fully consistent with the truth being an INTERMEDIATE, partial-synergy")
    say("   combination rule rather than a sequential one -- the rms winner still carries")
    say("   the OPPOSITE bias, of comparable size, to the rms loser. So 'sequential beats")
    say("   additive' must NOT be read as 'the rule is sequential'. The defensible")
    say("   statement is the bracket: the additive rule under-fills and the fully")
    say("   sequential rule over-fills, on every composition, and the truth is between.")
    say()

    # --- THE CONFOUND THIS DESIGN CANNOT BREAK -----------------------------------------
    say("=== CONFOUND: the residual is PERFECTLY rank-correlated with the coarse count ===")
    say("   Number fraction goes as y_i/d_i^3, so raising the fines VOLUME at fixed N")
    say("   necessarily strips grains from the coarse class. Composition and coarse-phase")
    say("   SAMPLE SIZE are therefore collinear by construction in this block, and a")
    say("   perfect rank correlation with a sample-size variable is what a finite-size")
    say("   artefact looks like. Printed here because it is the strongest argument")
    say("   AGAINST the block's own headline.")
    from scipy.stats import spearmanr
    nl, rs, ys = [], [], []
    for comp in COMPOSITIONS:
        m, _, _ = tphi[comp]
        p = make_ternary(comp)
        padd = mc.packing_density_tabulated(p.diameters, p.volume_fractions, beta, fa, fb)
        nl.append(int(p.class_counts(1500)[2]))
        rs.append((m - padd) / SCALE)
        ys.append(comp[2])
    say(f"   {'n_L':>6}{'y_S':>7}{'residual':>11}")
    for a_, b_, c_ in sorted(zip(nl, ys, rs)):
        say(f"   {a_:6d}{b_:7.2f}{c_:+11.2f}")
    for nm, v in (("n_L", nl), ("y_S", ys)):
        rho, pv = spearmanr(v, rs)
        say(f"   spearman(residual, {nm:>3}) = {rho:+.4f}  (p = {pv:.4f}, n = 6)")
    say("   -> The two cannot be separated at fixed N. The discriminator is to hold")
    say("      composition FIXED and change N, which moves n_L and moves nothing else:")
    say("      PREREG_TERNARY_N_2026-09-14.md, graded by grade_ternary_N.py.")
    say("      UNTIL THAT BLOCK IS GRADED, T1's FAILURE IS AMBIGUOUS between")
    say("      non-additivity and coarse-phase finite size.")
    say()

    # --- the ceiling check: is the measurement even below ideal filling? ----------------
    say("=== CEILING: the measurement must sit below perfect separation ===")
    say("   If a ternary packed ABOVE the a = b = 0 limit, the filling picture would be")
    say("   refuted outright rather than found non-additive, so this is checked, not")
    say("   assumed. It also sizes how much room the model had.")
    say(f"   {'y_L':>6}{'y_M':>6}{'y_S':>6}{'phi meas':>10}{'ADD':>9}{'ideal':>9}"
        f"{'meas as % of the ADD->ideal gap':>34}")
    for comp in COMPOSITIONS:
        m, _, _ = tphi[comp]
        p = make_ternary(comp)
        padd = mc.packing_density_tabulated(p.diameters, p.volume_fractions, beta, fa, fb)
        pid = mc.packing_density_ab(p.diameters, p.volume_fractions, beta, 0.0, 0.0)
        frac = 100.0 * (m - padd) / (pid - padd)
        say(f"   {comp[0]:6.2f}{comp[1]:6.2f}{comp[2]:6.2f}{m:10.4f}{padd:9.4f}"
            f"{pid:9.4f}{frac:33.1f} %")
    say("   All six sit between the additive prediction and the ideal ceiling, so the")
    say("   filling picture is INTACT and it is the COMBINATION RULE that is at issue.")
    say()

    # --- how big is the miss, against the effect it is a miss on? -----------------------
    say("=== SIZE: the residual against the effect it is a residual on ===")
    say("   A 7.5-scale miss means nothing without the size of what is being predicted.")
    for comp in COMPOSITIONS:
        m, _, _ = tphi[comp]
        p = make_ternary(comp)
        padd = mc.packing_density_tabulated(p.diameters, p.volume_fractions, beta, fa, fb)
        eff = (m - beta) / SCALE
        res = (m - padd) / SCALE
        say(f"   {comp[0]:.2f}/{comp[1]:.2f}/{comp[2]:.2f}: composition effect "
            f"{eff:+.1f} scale, additive model captures {100 * (padd - beta) / (m - beta):.0f} %, "
            f"residual {res:+.1f} scale = {100 * res / eff:.0f} % of the effect")
    say()

    with open(os.path.join(OUT, "diag_ternary_additivity.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LINES) + "\n")
    print(f"-> {os.path.join(OUT, 'diag_ternary_additivity.txt')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
