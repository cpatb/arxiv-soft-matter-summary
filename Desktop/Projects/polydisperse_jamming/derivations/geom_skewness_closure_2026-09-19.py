"""geom_skewness_closure_2026-09-19.py

Prints every number quoted in geom_skewness_closure_2026-09-19.md.

Question: is the family effect a SKEWNESS effect?  Tested against an EXTERNAL,
published, independently fitted closure that this project had no hand in.

  Desmond & Weeks, Phys. Rev. E 90, 022204 (2014), doi:10.1103/PhysRevE.90.022204
  (preprint arXiv:1303.4627), Eq. (4):

      phi_RCP = phi*_RCP + c1 delta + c2 S delta^2
      phi*_RCP = 0.634,  c1 = 0.0658,  c2 = 0.0857
      delta = sqrt(<dR^2>)/<R>            (NUMBER-weighted, scale-invariant)
      S     = <dR^3>/<dR^2>^(3/2)         (third standardised moment)

  !! c2 = 0.0857 is the PUBLISHED value.  arXiv:1303.4627v2 misprints it as 0.857
     in both the Fig. 2 caption and the body.  Resolved by API lookup 2026-09-19;
     the published value is the one that reproduces the paper's own stated residual
     and its own worked example (c2 * 0.4^2 * 1 = 0.0137 ~ "0.02"; 0.857 gives 0.137).

  D&W are FRICTIONLESS (inferred from their cited algorithm, Xu/Blawzdziewicz/O'Hern
  PRE 71, 061306 (2005); the word friction does not appear in their paper).  The
  like-for-like arm here is therefore mu_p = 0, which this project has for lognormal
  and uniform over the full CV ladder.

READ ONLY on out/.  Run:  python geom_skewness_closure_2026-09-19.py
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import geom_common_20260919 as G  # noqa: E402

CVS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
TC = G.T_CRIT_2DOF

PHI_STAR, C1, C2 = 0.634, 0.0658, 0.0857
DW_DELTA_MAX = 0.413          # their posted data's largest delta
DW_S_ABS_MAX = 2.50           # their posted data's largest |S|, binary packings only
DW_S_CONTINUOUS_MAX = 1.17    # largest |S| among their CONTINUOUS distributions


def sec(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


def dw(delta, S):
    return PHI_STAR + C1 * delta + C2 * S * delta**2


# ------------------------------------------------- 1. the skewness of each family
def skew_of(family: str, cv: float) -> float:
    """Third standardised moment of the NUMBER-weighted diameter law, closed form."""
    if cv == 0.0:
        return 0.0
    if family == "lognormal":
        return 3.0 * cv + cv**3            # exact for a lognormal with CV = cv
    if family in ("uniform", "bidisperse"):
        return 0.0                         # both are symmetric by construction
    raise ValueError(family)


def part0_c2():
    sec("0. SETTLING c2 = 0.0857 vs 0.857 HERE, not by trusting a lookup")
    print("   Two independent literature passes on 2026-09-19 came back with DIFFERENT")
    print("   values of c2 from the same paper: 0.0857 (reported as the PUBLISHED PRE")
    print("   value) and 0.857 (which is what arXiv:1303.4627v2 prints, in both the")
    print("   Fig. 2 caption and the body).  A factor of ten in the skewness term is")
    print("   not something to resolve by which agent spoke last, so it is settled")
    print("   here by arithmetic on quantities the paper states in words.")
    print("\n   (a) D&W's own worked example, quoted in their results section:")
    print("       'for delta = 0.40 and S = 0, phi_RCP is shifted upward by ~0.02.")
    print("        Fixing that value of delta, changing S to +-1 shifts phi_RCP by")
    print("        ~+-0.02.'")
    print(f"       delta term  c1 * 0.40                  = {C1*0.40:.5f}  vs '~0.02'")
    for c2 in (0.0857, 0.857):
        print(f"       skew term   c2 * 1 * 0.40^2, c2={c2:<7} = {c2*1*0.16:.5f}  "
              f"vs '~0.02'  -> {'CONSISTENT' if abs(c2*0.16-0.02) < 0.01 else 'OFF BY 7x'}")
    print("\n   (b) Admissibility.  phi_RCP cannot exceed the FCC limit pi/(3 sqrt 2) =")
    print(f"       {math.pi/(3*math.sqrt(2)):.6f}.  Evaluate Eq. (4) at their own largest")
    print("       posted |S| = 2.50 and largest delta = 0.413:")
    for c2 in (0.0857, 0.857):
        v = PHI_STAR + C1 * 0.413 + c2 * 2.50 * 0.413**2
        print(f"       c2={c2:<7} -> phi_RCP = {v:.5f}   "
              f"{'admissible' if v < math.pi/(3*math.sqrt(2)) else 'IMPOSSIBLE (above FCC)'}")
    print("\n   (c) This project's own data, section 5 below: at the matched-binary pair")
    print("       (delta = 0.5, S = 8/3) the two values predict a gap of")
    for c2 in (0.0857, 0.857):
        print(f"       c2={c2:<7} -> {c2*(8/3)*0.25:+.5f}   against a measured "
              f"+0.0635 +- 0.0032")
    print("\n   All three agree: c2 = 0.0857.  Used throughout.  Anyone citing the")
    print("   PREPRINT for this coefficient gets the skewness term ten times too large.")
    assert abs(C2 - 0.0857) < 1e-12


def part1():
    sec("1. THE THREE FAMILIES IN D&W's OWN COORDINATES (delta, S)")
    print("   delta and S are scale-invariant ratios, so radius vs diameter is")
    print("   immaterial; both are NUMBER-weighted, which these families are.")
    print("   Closed forms, verified below against a 10^6-sample draw:")
    print("     lognormal  S = 3 CV + CV^3")
    print("     uniform    S = 0 (symmetric on an interval)")
    print("     bidisperse S = 0 (symmetric two-point, 50:50 by number)")
    rng = np.random.default_rng(20260919)
    print("\n     family       CV    S closed    S sampled    |diff|   in D&W box?")
    for fam in ("lognormal", "uniform", "bidisperse"):
        for cv in (0.2, 0.5):
            n = 1_000_000
            if fam == "lognormal":
                s_ = math.sqrt(math.log(1 + cv * cv))
                d = np.exp(rng.normal(-0.5 * s_ * s_, s_, n))
            elif fam == "uniform":
                w = cv * math.sqrt(3.0)
                d = rng.uniform(1 - w, 1 + w, n)
            else:
                d = np.where(rng.random(n) < 0.5, 1 - cv, 1 + cv)
            m = d.mean()
            v = ((d - m) ** 2).mean()
            ssamp = ((d - m) ** 3).mean() / v**1.5
            scl = skew_of(fam, cv)
            inbox = (cv <= DW_DELTA_MAX) and (abs(scl) <= DW_S_ABS_MAX)
            print(f"     {fam:11s} {cv:.1f}   {scl:+8.4f}   {ssamp:+9.4f}   "
                  f"{abs(scl-ssamp):.4f}   {'yes' if inbox else 'NO -- extrapolation'}")
            assert abs(scl - ssamp) < 0.02
    print("\n   D&W's fitted box, from their posted 81-packing data set:")
    print(f"     delta in [0.049, {DW_DELTA_MAX}];  |S| <= {DW_S_ABS_MAX}, but |S| above")
    print(f"     {DW_S_CONTINUOUS_MAX} comes ONLY from their BINARY packings.  This")
    print("     project's lognormal reaches S = 1.625 at CV = 0.5, inside the binary")
    print("     range and above the continuous one.")


# ---------------------------------------- 2. the like-for-like frictionless test
def part2(fits):
    sec("2. THE FRICTIONLESS ARM -- D&W's increment against this project's, mu_p = 0")
    print("   D&W's absolute level cannot be compared (different protocol, different")
    print("   algorithm, and 0.634 is their own fitted intercept).  What IS comparable")
    print("   is the CV INCREMENT, because Eq. (4) is additive in the monodisperse")
    print("   intercept -- which is exactly the structure item A measured (additive")
    print("   separability in friction).  So compare")
    print("       D_measured(CV) = phi_J(CV) - phi_J(0)")
    print("   against")
    print("       D_DW(CV) = c1 delta + c2 S delta^2   with delta = CV.")
    out = {}
    for fam in ("lognormal", "uniform"):
        g = fits[(fits.family == fam) & (fits.mu_p == 0.0)]
        piv = g.pivot_table(index="seed", columns="cv_nom", values="phi_J_lin")
        print(f"\n   --- {fam}, mu_p = 0 (frictionless), n = {len(piv)} seeds ---")
        print("     CV      S      D_DW      D_meas     se       diff      t     "
              "ratio meas/DW")
        for cv in CVS[1:]:
            if cv not in piv.columns:
                continue
            d = (piv[cv] - piv[0.0]).to_numpy()
            m, se, t, n = G.tstat(d)
            S = skew_of(fam, cv)
            pred = C1 * cv + C2 * S * cv**2
            tt = (m - pred) / se if se > 0 else float("nan")
            out[(fam, cv)] = (pred, m, se, tt)
            print(f"     {cv:.1f}  {S:+6.3f}  {pred:+.5f}  {m:+.5f}  {se:.5f}  "
                  f"{m-pred:+.5f}  {tt:+6.2f}    {m/pred:.3f}")
        print(f"     (t_crit at 2 dof = {TC:.3f})")
    return out


# --------------------------------- 3. the discriminating half: does S do the work?
def part3(out):
    sec("3. IS THE FAMILY EFFECT A SKEWNESS EFFECT?  -- the discriminating half")
    print("   A closure that returns the same number for both families would 'agree'")
    print("   with neither family's shape, so the informative statistic is the RATIO")
    print("   lognormal/uniform, which in D&W is set ENTIRELY by the c2 S delta^2 term")
    print("   (uniform has S = 0 exactly).")
    for cv in (0.3, 0.4, 0.5):
        pl, ml, sel, _ = out[("lognormal", cv)]
        pu, mu_, seu, _ = out[("uniform", cv)]
        rp, rm = pl / pu, ml / mu_
        # se of the measured ratio by first order, treating the two families as
        # independent (they share the CV=0 cell, so this is CONSERVATIVE -- the
        # shared subtrahend cancels in a difference and inflates an independent var)
        ser = abs(rm) * math.sqrt((sel / ml) ** 2 + (seu / mu_) ** 2)
        print(f"     CV = {cv:.1f}:  D&W ratio {rp:.3f}   measured ratio {rm:.3f} "
              f"+- {ser:.3f}   t = {(rm-rp)/ser:+.2f}  -> "
              f"{'DISAGREE' if abs((rm-rp)/ser) > TC else 'agree'}")
    print("\n   AND THE HALF THAT CAN REFUTE D&W.  Uniform and bidisperse both have")
    print("   S = 0 EXACTLY and the same delta, so Eq. (4) gives them the SAME phi_RCP")
    print("   at every CV -- D&W's two parameters cannot tell them apart.  If the data")
    print("   separates them, the two-parameter closure is refuted as a complete")
    print("   description, whatever it gets right about lognormal-vs-uniform.")


def part4(fits):
    sec("4. UNIFORM vs BIDISPERSE -- the test D&W cannot pass")
    print("   Both S = 0, both delta = CV: D&W predicts D_unif = D_bidi at every CV.")
    print("   Only mu_p = 0.3 has both families.  They SHARE their CV = 0 cell")
    print("   (paper/main.tex:149-151), so the reference cancels EXACTLY in the")
    print("   difference and the comparison reduces to phi_J(CV) itself, seed-paired.")
    a = fits[(fits.family == "uniform") & (fits.mu_p == 0.3)].pivot_table(
        index="seed", columns="cv_nom", values="phi_J_lin")
    b = fits[(fits.family == "bidisperse") & (fits.mu_p == 0.3)].pivot_table(
        index="seed", columns="cv_nom", values="phi_J_lin")
    seeds = sorted(set(a.index) & set(b.index))
    print(f"   seeds: {seeds}")
    print("     CV     D_unif     D_bidi     difference   se       t      verdict")
    nres = 0
    for cv in CVS[1:]:
        du = (a.loc[seeds, cv] - a.loc[seeds, 0.0]).to_numpy()
        db = (b.loc[seeds, cv] - b.loc[seeds, 0.0]).to_numpy()
        m, se, t, n = G.tstat(du - db)
        res = abs(t) > TC
        nres += res
        print(f"     {cv:.1f}   {du.mean():+.5f}  {db.mean():+.5f}   {m:+.5f}   "
              f"{se:.5f}  {t:+6.2f}   "
              f"{'RESOLVED -> D&W REFUTED' if res else 'not resolved'}")
    print(f"\n   CVs at which the two S = 0 families are resolved apart: {nres} of 5.")
    print("   NOTE the friction mismatch: this arm is mu_p = 0.3 and D&W is")
    print("   frictionless.  Item A measured the friction interaction at CV = 0.5 to be")
    print("   +0.00404 in the lognormal family, 6.5% of the shift, so a difference")
    print("   between two families at matched friction is not made by friction.")
    return nres


# --------------------------------- 5. the matched-binary pair: D&W's hardest case
def part5():
    sec("5. THE MATCHED-BINARY PAIR -- D&W's prediction against the headline number")
    print("   The manuscript's pair: size ratio 3, number-CV 0.500, two compositions.")
    print("   From geom_psd_functional_class: q_large = 1/2 gives S = 0 and q_large =")
    print("   1/10 gives S = +8/3.  delta = 0.5 for both.  So D&W predicts a gap made")
    print("   entirely by the skewness term:")
    for q_l, label in ((0.5, "A (n_small = 0.500, X_s = 1/28)"),
                       (0.1, "B (n_small = 0.900, X_s = 0.25)")):
        n_s = 1 - q_l
        s_val = 1.0 / (n_s + q_l * 3.0)
        d = np.array([s_val, 3.0 * s_val])
        p = np.array([n_s, q_l])
        m = (p * d).sum()
        v = (p * (d - m) ** 2).sum()
        delta = math.sqrt(v) / m
        S = (p * (d - m) ** 3).sum() / v**1.5
        print(f"     {label}: delta = {delta:.6f}, S = {S:+.6f}, "
              f"phi_DW = {dw(delta, S):.6f}")
    gap_dw = dw(0.5, 8 / 3) - dw(0.5, 0.0)
    print(f"\n   D&W predicted gap = c2 * S_B * delta^2 = {C2:.4f} * {8/3:.6f} * 0.25 "
          f"= {gap_dw:+.6f}")
    print("   Manuscript measured gap (paper/numbers.tex \\pairDPhi, \\pairDPhiSe,")
    print("   at p_hat = 1e-3, mu_p = 0.3, N = 1500):  +0.0635 +- 0.0032")
    meas, se = 0.0635, 0.0032
    print(f"   difference = {meas-gap_dw:+.5f}  = {(meas-gap_dw)/se:+.2f} se   -> "
          f"{'DISAGREE' if abs(meas-gap_dw)/se > TC else 'AGREE at 2 dof'}")
    print("   The project's own two-exponent filling closure predicts +0.0604")
    print("   (paper/numbers.tex \\pairFillingPred), which is "
          f"{(0.0604-gap_dw):+.5f} from D&W.")
    print("\n   WHAT THIS DOES AND DOES NOT SHOW.")
    print("   (a) It is a genuine EXTERNAL prediction: D&W's three constants were fitted")
    print("       in 2013 to 81 frictionless packings of other distributions, with no")
    print("       knowledge of this pair, and the pair is the manuscript's headline.")
    print("   (b) S = +2.667 is OUTSIDE D&W's posted range (max |S| = 2.50) by 6.7%,")
    print("       and that end of their range is carried by BINARY packings only --")
    print("       which is, at least, the same kind of distribution as this pair.")
    print("   (c) It is a comparison of INCREMENTS, not levels: D&W's absolute numbers")
    print("       (0.6669 and 0.7240) are frictionless and the measurement is at")
    print("       mu_p = 0.3 (0.6042 and 0.6677).  The comparison is only licensed by")
    print("       item A's additive separability, and inherits its ~6% residual.")
    print("   (d) One agreement is not a validation.  The refuting test is section 4.")
    return gap_dw, meas, se


def main():
    fits = G.load_fits()
    part0_c2()
    part1()
    out = part2(fits)
    part3(out)
    nres = part4(fits)
    gap_dw, meas, se = part5()
    sec("VERDICT")
    ln5 = out[("lognormal", 0.5)]
    un5 = out[("uniform", 0.5)]
    print(f"   Frictionless arm, CV 0 -> 0.5:  lognormal D&W {ln5[0]:+.5f} vs measured "
          f"{ln5[1]:+.5f} +- {ln5[2]:.5f} (t = {ln5[3]:+.2f});")
    print(f"                                   uniform   D&W {un5[0]:+.5f} vs measured "
          f"{un5[1]:+.5f} +- {un5[2]:.5f} (t = {un5[3]:+.2f}).")
    print(f"   Matched-binary pair: D&W {gap_dw:+.5f} vs measured {meas:+.5f} +- {se:.5f}"
          f"  ({(meas-gap_dw)/se:+.2f} se).")
    print(f"   Uniform vs bidisperse, where D&W must predict ZERO difference: resolved "
          f"apart at {nres} of 5 CVs.")
    print("\n   VERDICT: PARTIALLY SUPPORTED, AND REFUTED AS A COMPLETE DESCRIPTION.")
    print("   Skewness at fixed width is the right leading variable -- an external")
    print("   two-parameter closure fitted to other data predicts the manuscript's")
    print("   headline gap to within its own quoted accuracy -- but (delta, S) is not")
    print("   sufficient, because two S = 0 families at identical delta are measured")
    print("   apart.")
    print("\n   ANCHOR SCOPE.  This anchor is external and could have failed: the")
    print("   coefficients are someone else's and the pair is not in their data set.")
    print("   It CANNOT decide the friction question (D&W is frictionless and the pair")
    print("   was measured at mu_p = 0.3), it CANNOT decide the level (only increments),")
    print("   and at S = 2.667 it is a 6.7% extrapolation past their own range.")


if __name__ == "__main__":
    main()
