"""three_knobs_check.py -- the symbolic checks `three_knobs.py` claims.

⚠ THIS FILE EXISTS BECAUSE `three_knobs.py` CITED IT BEFORE IT DID. The printed output
said "SymPy-verified in three_knobs_check" while no such file was on disk -- a claim
pointing at nothing, which is the defect this repository records as prose contradicting
its own artifact. Written 2026-09-22.

  C1  The Wyart-Cates closure is EXACTLY covariant under (phi, phi0, phim) -> c*(...).
  C2  It therefore depends only on x = phi/phi0 and chi = phim/phi0 (and fhat).
  C3  The phase-diagram boundaries: phiJ(0) = phi0, phiJ(1) = phim, so the shear-jamming
      window is chi <= x < 1.
  C4  The thickening ratio Theta = eta(fhat=1)/eta(fhat=0) = chi^2(1-x)^2/(chi-x)^2,
      and it DIVERGES exactly at x = chi -- so "shear-jams" and "Theta infinite" are the
      same statement rather than two.
  C5  Limiting cases, and a DISCRIMINATING negative: Theta must NOT be 1 when chi = 1.
  C6  A numeric spot-check against a direct evaluation of the unreduced closure.

Run:  python3 three_knobs_check.py
"""

from __future__ import annotations

import sys

import sympy as sp

FAIL: list[str] = []


def check(name: str, ok, detail: str = "") -> None:
    ok = bool(ok)
    print(f"   [{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


def main() -> int:
    print("=" * 92)
    print("three_knobs_check.py -- the reduction, symbolically")
    print("=" * 92)

    phi, phi0, phim, fhat, eta_s, c = sp.symbols(
        "phi phi0 phim fhat eta_s c", positive=True)
    x, chi = sp.symbols("x chi", positive=True)

    phiJ = phi0 - (phi0 - phim) * fhat
    eta = eta_s * (1 - phi / phiJ) ** -2

    print("\n=== C1  exact covariance ===")
    scaled = eta.subs({phi: c * phi, phi0: c * phi0, phim: c * phim})
    check("C1  eta is invariant under (phi, phi0, phim) -> c*(phi, phi0, phim)",
          sp.simplify(scaled - eta) == 0)

    print("\n=== C2  the two-parameter reduction ===")
    eta_red = eta_s * (1 - x / (1 - (1 - chi) * fhat)) ** -2
    check("C2  eta(phi, phi0, phim) == eta(x, chi) under phi = x*phi0, phim = chi*phi0",
          sp.simplify(eta.subs({phi: x * phi0, phim: chi * phi0}) - eta_red) == 0)

    print("\n=== C3  the phase-diagram boundaries ===")
    check("C3a phiJ(fhat = 0) = phi0", sp.simplify(phiJ.subs(fhat, 0) - phi0) == 0)
    check("C3b phiJ(fhat = 1) = phim", sp.simplify(phiJ.subs(fhat, 1) - phim) == 0)
    # phiJ is monotone in fhat, so phi >= phim is the condition to jam at SOME fhat
    check("C3c phiJ is monotone decreasing in fhat (so phim is its minimum)",
          sp.simplify(sp.diff(phiJ, fhat) + (phi0 - phim)) == 0)
    print("       => jams at some fhat iff phi >= phim  <=>  x >= chi")
    print("       => jammed already at rest iff phi >= phi0 <=> x >= 1")
    print("       => the SHEAR-JAMMING WINDOW is chi <= x < 1, width 1 - chi")

    print("\n=== C4  the thickening ratio ===")
    theta = sp.simplify(eta_red.subs(fhat, 1) / eta_red.subs(fhat, 0))
    closed = chi ** 2 * (1 - x) ** 2 / (chi - x) ** 2
    check("C4a Theta = chi^2 (1-x)^2 / (chi - x)^2",
          sp.simplify(theta - closed) == 0, f"Theta = {sp.simplify(theta)}")
    # ⚠ THE CONDITION `chi < 1` IS NOT OPTIONAL, and this check is what found that. The
    # first version asserted "Theta -> oo as x -> chi" unconditionally and FAILED: SymPy
    # returns `oo*sign((1-chi)^2)`, because at `chi = 1` the numerator vanishes too and
    # the limit is 1, not infinity. `chi = 1` is exactly `phim = phi0` -- no contrast, no
    # thickening (C5a) -- so the divergence is a statement about a material that HAS two
    # branches. Verified at concrete chi in (0, 1) rather than asserted generically.
    div = []
    for cval in (sp.Rational(1, 2), sp.Rational(3, 4), sp.Rational(9, 10)):
        div.append(sp.limit(closed.subs(chi, cval), x, cval, "-") == sp.oo)
    check("C4b Theta -> oo as x -> chi from below, FOR chi < 1", all(div),
          f"checked at chi = 1/2, 3/4, 9/10: {div}")
    check("C4b' and NOT at chi = 1, where Theta == 1 identically",
          sp.limit(closed.subs(chi, 1), x, 1, "-") == 1)
    check("C4c eta_s cancels out of Theta", eta_s not in theta.free_symbols)

    print("\n=== C5  LIMITING CASES, including a discriminating negative ===")
    # chi -> 1 means phim = phi0: the two branches coincide, there is NO thickening
    check("C5a chi = 1 (phim = phi0) gives Theta = 1 exactly -- no thickening at all",
          sp.simplify(closed.subs(chi, 1) - 1) == 0)
    # x -> 0 is a dilute suspension: eta -> eta_s on both branches, so Theta -> 1
    check("C5b x -> 0 (dilute) gives Theta -> 1",
          sp.simplify(sp.limit(closed, x, 0)) == 1)
    # DISCRIMINATING: Theta must NOT be 1 for a generic (x, chi) -- a formula that
    # returned 1 everywhere would pass C5a and C5b and be useless.
    generic = closed.subs({x: sp.Rational(1, 2), chi: sp.Rational(3, 4)})
    check("C5c DISCRIMINATING: Theta != 1 at a generic point (x=1/2, chi=3/4)",
          sp.simplify(generic - 1) != 0, f"Theta = {generic} = {float(generic):.4f}")
    # and Theta must exceed 1 wherever there is any contrast and any concentration
    check("C5d Theta > 1 whenever 0 < x < chi < 1",
          bool(sp.simplify(closed.subs({x: sp.Rational(1, 3), chi: sp.Rational(2, 3)})
                           > 1)))

    print("\n=== C6  NUMERIC SPOT-CHECK against the UNREDUCED closure ===")
    print(f"   {'phi':>7}{'phi0':>7}{'phim':>7}{'fhat':>6}"
          f"{'eta (unreduced)':>18}{'eta (x, chi)':>16}{'rel err':>11}")
    worst = 0.0
    for p, p0, pm, f in ((0.52, 0.62, 0.50, 0.3), (0.44, 0.70, 0.5176, 0.9),
                         (0.59, 0.64, 0.57, 0.0), (0.30, 0.64, 0.50, 1.0)):
        a = float(eta.subs({phi: p, phi0: p0, phim: pm, fhat: f, eta_s: 1}))
        b = float(eta_red.subs({x: p / p0, chi: pm / p0, fhat: f, eta_s: 1}))
        rel = abs(a - b) / abs(a)
        worst = max(worst, rel)
        print(f"   {p:7.3f}{p0:7.3f}{pm:7.4f}{f:6.1f}{a:18.8f}{b:16.8f}{rel:11.2e}")
    check("C6  the reduced and unreduced forms agree numerically", worst < 1e-12,
          f"worst relative error {worst:.2e}")

    print("\n" + "=" * 92)
    if FAIL:
        print(f"FAILED {len(FAIL)} check(s):")
        for f in FAIL:
            print(f"   - {f}")
    else:
        print("ALL CHECKS PASSED")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
