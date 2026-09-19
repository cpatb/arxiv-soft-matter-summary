"""geom_psd_functional_class_2026-09-19.py

Prints every number quoted in geom_psd_functional_class_2026-09-19.md.

Question: the matched-binary pair rules out WIDTH descriptors.  Exactly which class
of functionals of the PSD does it rule out, and which does it leave standing?

Derived here, with SymPy, and checked against the manuscript's own pair:
  * a two-point PSD at fixed (mean, CV, size ratio) is the root set of ONE quadratic in
    the composition -- generically TWO points, and they are the manuscript's A and B;
  * the skewness of ANY two-point distribution is (1-2q)/sqrt(q(1-q)) with q the number
    fraction of the heavy atom -- a function of COMPOSITION ALONE, independent of the
    two diameters, so it cannot be held fixed at two compositions;
  * therefore the pair rules out exactly the functionals that are constant on the level
    set {mean, CV, support}, and rules out nothing that reads the third moment.

READ ONLY on out/.  Run:  python geom_psd_functional_class_2026-09-19.py
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import geom_common_20260919 as G  # noqa: E402

sys.path.insert(0, G.ROOT)


def sec(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


# ------------------------------------------------ 1. the two-point level set
def part1():
    sec("1. THE LEVEL SET {mean, CV, size ratio} ON TWO-POINT PSDs IS A QUADRATIC")
    q, r, c = sp.symbols("q r c", positive=True)
    s = sp.symbols("s", positive=True)
    # atoms: s (small) with number fraction 1-q, r*s (large) with number fraction q
    m1 = (1 - q) * s + q * r * s
    m2 = (1 - q) * s**2 + q * (r * s) ** 2
    print("   Atoms d_s = s and d_l = r s, r > 1; q = NUMBER fraction of the LARGE atom.")
    print("   <d>   =", sp.factor(m1))
    print("   <d^2> =", sp.factor(m2))
    s_of = sp.solve(sp.Eq(m1, 1), s)[0]
    print("   Fix <d> = 1  ->  s =", sp.simplify(s_of))
    cv2 = sp.simplify((m2 / m1**2 - 1))
    print("   CV^2 = <d^2>/<d>^2 - 1 =", sp.factor(sp.simplify(cv2)))
    eq = sp.simplify(sp.expand(sp.numer(sp.together(cv2 - c**2))))
    poly = sp.Poly(sp.expand(eq), q)
    print("\n   Setting CV^2 = c^2 and clearing denominators gives a polynomial in q of")
    print("   degree", poly.degree(), ":")
    print("     ", sp.factor(eq), " = 0")
    # canonical form A u^2 q^2 + u(2A - r - 1) q + c^2 = 0 with A = 1+c^2, u = r-1
    A, u = 1 + c**2, r - 1
    canon = A * u**2 * q**2 + u * (2 * A - r - 1) * q + c**2
    diff = sp.simplify(sp.expand(eq) - sp.expand(canon * sp.LC(poly) / sp.LC(sp.Poly(canon, q))))
    print("\n   Canonical form, asserted equal up to the leading coefficient:")
    print("       (1+c^2)(r-1)^2 q^2 + (r-1)(2(1+c^2) - r - 1) q + c^2 = 0")
    print("   SymPy difference from the cleared polynomial:", sp.simplify(diff))
    assert sp.simplify(diff) == 0
    roots = sp.solve(canon, q)
    print("\n   Vieta (exact):")
    prod = sp.simplify(roots[0] * roots[1])
    summ = sp.simplify(roots[0] + roots[1])
    print("     q1 q2   =", sp.simplify(prod), "  = c^2 / ((1+c^2)(r-1)^2)")
    print("     q1 + q2 =", sp.factor(sp.simplify(summ)))
    assert sp.simplify(prod - c**2 / ((1 + c**2) * (r - 1) ** 2)) == 0
    print("   -> TWO compositions realise any admissible (mean, CV, ratio).  They are")
    print("      the manuscript's mixtures A and B.  The level set is a PAIR, not a")
    print("      continuum, so the pair exhausts what can be matched this way.")
    disc = sp.simplify(sp.discriminant(canon, q))
    print("\n   Discriminant =", sp.factor(disc))
    print("   Real distinct roots require it > 0; equal roots are the tangency where the")
    print("   pair degenerates to one mixture and the experiment is impossible.")
    return canon, q, r, c


# ---------------------------------------------- 2. skewness along the level set
def part2():
    sec("2. SKEWNESS OF A TWO-POINT PSD DEPENDS ONLY ON THE COMPOSITION")
    q, a, b = sp.symbols("q a b", positive=True)
    # X takes value b with prob q, a with prob 1-q.  Any two-point law is an affine
    # image of a Bernoulli, and standardized moments are affine-invariant.
    m1 = (1 - q) * a + q * b
    m2 = (1 - q) * a**2 + q * b**2
    m3 = (1 - q) * a**3 + q * b**3
    var = sp.simplify(m2 - m1**2)
    mu3 = sp.simplify(m3 - 3 * m1 * m2 + 2 * m1**3)
    skew = sp.simplify(mu3 / var ** sp.Rational(3, 2))
    print("   General two atoms a < b with P(b) = q:")
    print("     variance    =", sp.factor(var))
    print("     mu_3        =", sp.factor(mu3))
    print("     skewness    =", sp.simplify(skew))
    target = (1 - 2 * q) / sp.sqrt(q * (1 - q))
    d = sp.simplify(sp.powsimp(sp.simplify(skew - target), force=True))
    print("\n   Claim: skewness = (1 - 2q)/sqrt(q(1-q)), FREE of a and b.")
    print("   SymPy difference (with b > a assumed):",
          sp.simplify(d.subs({a: 1, b: 3})), "at (a,b)=(1,3);",
          sp.simplify(d.subs({a: 0.7, b: 5.3})), "at (a,b)=(0.7,5.3)")
    # verify numerically on a grid rather than fight assumptions in the CAS
    print("\n   Numeric verification over a grid of (a, b, q) -- worst |difference|:")
    worst = 0.0
    for aa in (0.2, 0.5, 1.0, 2.0):
        for bb in (aa * 1.5, aa * 3.0, aa * 10.0):
            for qq in (0.05, 0.1, 0.25, 0.5, 0.75, 0.9):
                p = np.array([1 - qq, qq])
                x = np.array([aa, bb])
                m = (p * x).sum()
                v = (p * (x - m) ** 2).sum()
                s3 = (p * (x - m) ** 3).sum() / v ** 1.5
                worst = max(worst, abs(s3 - (1 - 2 * qq) / math.sqrt(qq * (1 - qq))))
    print(f"     worst |skew_numeric - (1-2q)/sqrt(q(1-q))| over 72 points = {worst:.3e}")
    assert worst < 1e-12
    print("\n   CONSEQUENCE.  gamma_1(q) = (1-2q)/sqrt(q(1-q)) has")
    dg = sp.simplify(sp.diff(target, q))
    print("     d gamma_1/dq =", sp.factor(dg), "  < 0 on (0,1), so gamma_1 is STRICTLY")
    print("     DECREASING and therefore INJECTIVE in q.")
    print("   Two distinct compositions => two distinct skewnesses.  A two-point PSD")
    print("   CANNOT hold CV and skewness fixed at two compositions.  QED.")
    print("   (Corollary: it also cannot hold ANY standardized moment fixed -- every")
    print("    standardized moment of a two-point law is a function of q alone, because")
    print("    the law is an affine image of Bernoulli(q) and standardized moments are")
    print("    affine-invariant.  Verified for the 4th below.)")
    m4 = (1 - q) * a**4 + q * b**4
    mu4 = sp.simplify(m4 - 4 * m1 * m3 + 6 * m1**2 * m2 - 3 * m1**4)
    kurt = sp.simplify(mu4 / var**2)
    ktarget = (1 - 3 * q * (1 - q)) / (q * (1 - q))
    print("     symbolic kurtosis =", kurt)
    print("     claimed form      =", sp.simplify(ktarget))
    worst4 = 0.0
    for aa in (0.2, 1.0, 2.0):
        for bb in (aa * 1.5, aa * 4.0):
            for qq in (0.1, 0.3, 0.5, 0.9):
                p = np.array([1 - qq, qq])
                x = np.array([aa, bb])
                m = (p * x).sum()
                v = (p * (x - m) ** 2).sum()
                k4 = (p * (x - m) ** 4).sum() / v**2
                worst4 = max(worst4, abs(k4 - (1 - 3 * qq * (1 - qq)) / (qq * (1 - qq))))
    print(f"     kurtosis: worst |k4 - (1-3q(1-q))/(q(1-q))| over 24 points = "
          f"{worst4:.3e}  -> also q-only")
    assert worst4 < 1e-10
    return target


# ------------------------------------------ 3. the manuscript's pair, reproduced
def part3(canon, qsym, rsym, csym, skew_of_q):
    sec("3. THE MANUSCRIPT'S MATCHED PAIR, REPRODUCED FROM THE CLOSED FORM")
    print("   paper/main.tex sec:binaries: ratio 3, CV 0.500, X_s = 1/28 and 0.25,")
    print("   n_small = 0.500 and 0.900, skewness 0.00 and +2.67, phi gap +0.0635.")
    print("   Those are INPUTS to the manuscript.  Here they are OUTPUTS of the")
    print("   quadratic, which is the check that the level set is what I claim.")
    rr, cc = sp.Rational(3), sp.Rational(1, 2)
    poly = canon.subs({rsym: rr, csym: cc})
    roots = sorted([sp.nsimplify(x) for x in sp.solve(sp.expand(poly), qsym)],
                   key=lambda z: float(z))
    print(f"\n   Quadratic at r = 3, CV = 1/2:  {sp.expand(poly)} = 0")
    print(f"   Roots (q = number fraction of the LARGE atom): "
          f"{[sp.nsimplify(x) for x in roots]}")
    for q_large in roots:
        q_l = float(q_large)
        n_small = 1 - q_l
        s_val = 1.0 / (n_small + q_l * 3.0)
        d_s, d_l = s_val, 3.0 * s_val
        # X_s = fines VOLUME fraction
        Xs = n_small * d_s**3 / (n_small * d_s**3 + q_l * d_l**3)
        # skewness with q = number fraction of the HEAVY (large) atom
        sk = (1 - 2 * q_l) / math.sqrt(q_l * (1 - q_l))
        # direct moment check
        p = np.array([n_small, q_l])
        x = np.array([d_s, d_l])
        m = (p * x).sum()
        v = (p * (x - m) ** 2).sum()
        cv_chk = math.sqrt(v) / m
        sk_chk = (p * (x - m) ** 3).sum() / v**1.5
        print(f"\n     q_large = {q_l:.4f}  ->  n_small = {n_small:.4f},  "
              f"d_s = {d_s:.6f}, d_l = {d_l:.6f}")
        print(f"       realised mean = {m:.10f}   CV = {cv_chk:.10f}   "
              f"skewness = {sk_chk:+.6f}")
        print(f"       closed-form skewness (1-2q)/sqrt(q(1-q)) = {sk:+.6f}   "
              f"|diff| = {abs(sk-sk_chk):.2e}")
        print(f"       fines VOLUME fraction X_s = {Xs:.6f}  (= 1/{1/Xs:.4f})")
        assert abs(cv_chk - 0.5) < 1e-12
        assert abs(sk - sk_chk) < 1e-12
    print("\n   Against the manuscript: n_small 0.500 / 0.900, X_s 0.0357 / 0.25,")
    print("   skewness 0.00 / +2.67.  Reproduced above to the printed precision.")
    print("   NOTE the sign convention: the manuscript's +2.67 is the mixture with")
    print("   n_small = 0.900, i.e. q_large = 0.1 -- many small grains and a few large,")
    print("   which is right-skewed.  (1-2*0.1)/sqrt(0.1*0.9) = +2.6667.")


# ------------------------------------------- 4. which class is actually excluded
def part4():
    sec("4. WHICH CLASS OF FUNCTIONALS THE PAIR EXCLUDES -- stated precisely")
    print("   Let P be the set of probability measures on (0, inf) (number-weighted")
    print("   diameter laws).  The manuscript's pair is {P_A, P_B} subset P with")
    print("       <d>_A = <d>_B,  CV_A = CV_B,  supp(P_A) = supp(P_B) = {d_s, d_l},")
    print("   hence also equal size ratio, equal span, equal min and max, equal atom")
    print("   count (2), and equal N.  Measured phi_A != phi_B by 0.0635 +- 0.0032.")
    print("\n   EXCLUDED: every functional Phi[P] that is constant on the level set")
    print("       L = { P in P : <d> = 1, CV = 1/2, supp(P) = {d_s, d_l} }.")
    print("   That set contains, in particular:")
    print("     * any function of the first two moments alone -- <d>, CV, variance,")
    print("       d32/d10 IS NOT in this list (see below), sigma/mu, span;")
    print("     * any function of the SUPPORT alone -- size ratio r, d_max/d_min,")
    print("       d_max - d_min, the number of species K;")
    print("     * any function of (moments 1-2) and (support) jointly;")
    print("     * therefore every closure of the form phi_J = f(CV) or f(CV, r) or")
    print("       f(span), which is what 'width is not the descriptor' means.")
    print("\n   NOT EXCLUDED -- and this is the part the phrase 'no descriptor of the")
    print("   size distribution' overstates:")
    print("     * anything reading the THIRD or higher moment.  Skewness differs by")
    print("       2.6667 across the pair, so a skewness-dependent closure is untouched.")
    print("     * anything reading a VOLUME-weighted moment.  Check, numerically:")
    rr = 3.0
    for q_l in (0.5, 0.1):
        n_s = 1 - q_l
        s_val = 1.0 / (n_s + q_l * rr)
        d = np.array([s_val, rr * s_val])
        p = np.array([n_s, q_l])
        mom = {k: float((p * d**k).sum()) for k in range(1, 5)}
        d32 = mom[3] / mom[2]
        d43 = mom[4] / mom[3]
        Xs = p[0] * d[0] ** 3 / (p * d**3).sum()
        cvvol = math.sqrt((p * d**3 * (d - mom[4] / mom[3]) ** 2).sum()
                          / (p * d**3).sum()) / (mom[4] / mom[3])
        print(f"       q_large={q_l:.2f}: d32/d10 = {d32/mom[1]:.6f}, "
              f"d43/d10 = {d43/mom[1]:.6f}, X_s = {Xs:.6f}, "
              f"volume-weighted CV = {cvvol:.6f}")
    print("       -> d32/d10, d43/d10, X_s and the volume-weighted CV ALL DIFFER across")
    print("          the pair, so none of them is excluded by it.  The pair matches the")
    print("          NUMBER-weighted first two moments; volume-weighted moments are")
    print("          functions of moments 3 and 4 and are free to move.")
    print("\n   SO THE PRECISE STATEMENT IS:")
    print("     The pair excludes every closure that is a function of the NUMBER-")
    print("     weighted mean, the NUMBER-weighted CV, and the SUPPORT.  It does not")
    print("     exclude a closure reading any higher number moment, nor any volume-")
    print("     weighted moment, nor the composition itself.  'No descriptor of the size")
    print("     distribution works' is too strong; 'no descriptor of its WIDTH works'")
    print("     is exact.")
    print("\n   AND THE LIMIT OF THE PAIR AS EVIDENCE:  because a two-point law is an")
    print("   affine image of a Bernoulli, EVERY standardized moment of it is a function")
    print("   of q alone (section 2).  So within two-point PSDs, matching (mean, CV,")
    print("   support) and varying composition necessarily varies EVERY higher")
    print("   standardized moment simultaneously.  The pair therefore cannot attribute")
    print("   the 0.0635 to the third moment rather than the fourth, or to any single")
    print("   moment at all.  It is a proof of insufficiency of width, not an")
    print("   identification of the sufficient descriptor.")


# ------------------------- 5. what a THREE-point PSD would buy, and it is decisive
def part5():
    sec("5. THE MEASUREMENT THAT WOULD IDENTIFY THE MOMENT -- a three-atom PSD")
    print("   A K-atom PSD on a fixed support has K-1 free composition parameters.")
    print("   Fixing <d> and CV leaves K-3 free.  So:")
    print("     K = 2: 0 free  -> the level set is the two Vieta roots (section 1).")
    print("     K = 3: 1 free  -> a CURVE of PSDs matched in mean, CV and support,")
    print("                       along which the skewness VARIES CONTINUOUSLY.")
    print("     K = 4: 2 free  -> a surface; skewness AND kurtosis can be varied")
    print("                       INDEPENDENTLY, which is what separates the moments.")
    print("   So a three-atom PSD already gives a skewness AXIS at fixed width, and a")
    print("   four-atom one gives the first genuine skewness-vs-kurtosis discriminator.")
    print("\n   Constructed here: three atoms at d = {1, u, u^2} scaled to <d> = 1, with")
    print("   number fractions (n1, n2, n3) solving <d> = 1 and CV = 1/2, one parameter")
    print("   free.  The realised skewness range along that curve is printed.")
    for u in (2.0, 3.0):
        d0 = np.array([1.0, u, u * u])
        best = []
        for n3 in np.linspace(1e-4, 0.9999, 20001):
            # solve n1 + n2 = 1 - n3 and the CV constraint for n2 given n3
            # <d> = 1 fixes the overall scale; work with shape then rescale.
            for n2 in np.linspace(0.0, 1 - n3, 2001):
                n1 = 1 - n2 - n3
                if n1 < 0:
                    continue
                p = np.array([n1, n2, n3])
                m1 = (p * d0).sum()
                d = d0 / m1
                m = 1.0
                v = (p * (d - m) ** 2).sum()
                if abs(math.sqrt(v) - 0.5) < 2e-3:
                    s3 = (p * (d - m) ** 3).sum() / v**1.5
                    best.append((s3, n1, n2, n3))
            if len(best) > 4000:
                break
        if best:
            sk = np.array([b[0] for b in best])
            print(f"     support {{1, {u}, {u*u}}}, CV = 0.500 +- 0.002: "
                  f"{len(best)} solutions, skewness spans {sk.min():+.3f} to "
                  f"{sk.max():+.3f} (range {sk.max()-sk.min():.3f})")
        else:
            print(f"     support {{1, {u}, {u*u}}}: no solution at CV = 0.500")
    print("\n   For comparison the two-point pair offers exactly two skewnesses,")
    print("   0.000 and +2.667.  A three-atom family at the same CV sweeps a")
    print("   continuum, so the skewness coefficient becomes estimable rather than")
    print("   merely bounded.")
    print("\n   COST NOTE: the project already has a DiscretePSD spec")
    print("   ('discrete:d1/2/4:y.../...' in psd.make_psd), so this needs no new")
    print("   machinery -- only a composition grid at fixed number-CV.")


def main():
    canon, q, r, c = part1()
    sk = part2()
    part3(canon, q, r, c, sk)
    part4()
    part5()
    sec("VERDICT")
    print("   The closed form is exact and reproduces the manuscript's own pair.")
    print("   VERDICT: SUPPORTED -- the two-point level set is the pair of Vieta roots,")
    print("   and a two-point PSD cannot hold CV and ANY standardized moment fixed at")
    print("   two compositions, because every standardized moment of a two-point law is")
    print("   a function of the composition alone.")
    print("\n   ANCHOR SCOPE.  This is pure moment algebra, so it can decide exactly")
    print("   which class of closures the matched pair excludes, and it does so with no")
    print("   reference to any packing measurement.  It CANNOT say which functional is")
    print("   the right one -- it shows the pair is insufficient to identify one, which")
    print("   is a statement about the design, not about packing.")


if __name__ == "__main__":
    main()
