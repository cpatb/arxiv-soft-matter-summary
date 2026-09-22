"""multiplet_psd_2026-09-22_check.py -- the four attacks on the multiplet derivation.

Checks, in the order derive-and-check prescribes:

  C1  SymPy: the closed-form moment identities against a direct symbolic sum (K = 3).
  C2  Limiting cases: gamma = 1, monomers-only, all-k-mers, and the K -> 1 collapse.
  C3  Dimensions: every expression, by hand and with pint.
  C4  Numeric spot-check: Monte-Carlo realisations of an actual multiplet population
      against the closed forms, at >= 3 parameter points.
  C5  The bound Lambda >= 1 and the EXISTENCE of a critical gamma*, both proved in the
      .md and here tested numerically over a grid.
  C6  Anchor that does not pass through the code under test: the closed forms are
      re-derived by an INDEPENDENT route (explicit per-particle arrays, no moments).

Run:  python3 derivations/multiplet_psd_2026-09-22_check.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

FAIL: list[str] = []
OUT: list[str] = []


def say(s: str = "") -> None:
    print(s)
    OUT.append(s)


def check(name: str, ok: bool, detail: str = "") -> None:
    say(f"   [{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


# =============================================================================================
# The model, written once, in one place.
#
#   monomer diameter          d1
#   a k-mer contains          k monomers, so SOLID volume = k * v1,  v1 = (pi/6) d1^3
#   its EFFECTIVE (packing) sphere has diameter   d_k = k^(gamma/3) d1
#   so its EFFECTIVE volume is                    v_k = k^gamma v1
#   with the compaction exponent                  gamma = 3 / d_f
#
# gamma = 1  <=> d_f = 3 : compact, effective volume == solid volume, no trapped fluid.
# gamma > 1  <=> d_f < 3 : open aggregate, encloses fluid.
# =============================================================================================


def moments(n, k, s):
    """M_s = <k^s> under NUMBER fractions n. Plain sum; this is the definition."""
    n = np.asarray(n, float)
    k = np.asarray(k, float)
    return float(np.sum(n * k ** s))


def closed_forms(n, k, gamma):
    """The five closed forms the .md derives, eqs (3)-(7)."""
    g3 = gamma / 3.0
    M = lambda s: moments(n, k, s)
    cv2 = M(2 * g3) / M(g3) ** 2 - 1.0
    return dict(
        d10_over_d1=M(g3),
        cv=float(np.sqrt(max(cv2, 0.0))),
        d32_over_d10=M(gamma) / (M(2 * g3) * M(g3)),
        d43_over_d10=M(4 * g3) / (M(gamma) * M(g3)),
        Lambda=M(gamma) / M(1.0),
    )


def direct_from_particles(n, k, gamma, d1=1.0, npart=2_000_000, seed=20260922):
    """C6 ANCHOR: the same five numbers built from an explicit particle array.

    This never forms a moment. It draws k for each particle, converts to a diameter and
    a solid volume ONE PARTICLE AT A TIME, and reduces with plain means. If the closed
    forms are wrong, this disagrees.
    """
    rng = np.random.default_rng(seed)
    kk = rng.choice(np.asarray(k, float), size=npart, p=np.asarray(n, float))
    d = d1 * kk ** (gamma / 3.0)
    v_eff = (np.pi / 6.0) * d ** 3
    v_sol = kk * (np.pi / 6.0) * d1 ** 3
    return dict(
        d10_over_d1=float(d.mean() / d1),
        cv=float(d.std(ddof=0) / d.mean()),
        d32_over_d10=float((d ** 3).mean() / (d ** 2).mean() / d.mean()),
        d43_over_d10=float((d ** 4).mean() / (d ** 3).mean() / d.mean()),
        Lambda=float(v_eff.sum() / v_sol.sum()),
    )


def main() -> int:
    say("=" * 92)
    say("multiplet_psd_2026-09-22 -- derivation checks")
    say("=" * 92)

    # ---------------------------------------------------------------------------- C1
    say("\n=== C1  SYMPY: closed forms against a direct symbolic sum, K = 3 ===")
    n1, n2, n3, g = sp.symbols("n1 n2 n3 gamma", positive=True)
    d1s = sp.symbols("d_1", positive=True)
    ns = [n1, n2, n3]
    ks = [sp.Integer(1), sp.Integer(2), sp.Integer(3)]
    # normalisation is imposed at the end, so carry sum(n) symbolically
    S = sum(ns)
    d = [d1s * kk ** (g / 3) for kk in ks]

    def E(expr_list):
        return sum(ni * ei for ni, ei in zip(ns, expr_list)) / S

    def Msym(s):
        return sum(ni * kk ** s for ni, kk in zip(ns, ks)) / S

    # (3) d10
    lhs = E(d)
    rhs = d1s * Msym(g / 3)
    check("C1a  d10 = d1 * M_{gamma/3}", sp.simplify(lhs - rhs) == 0)

    # (4) CV^2
    lhs = E([x ** 2 for x in d]) / E(d) ** 2 - 1
    rhs = Msym(2 * g / 3) / Msym(g / 3) ** 2 - 1
    check("C1b  CV^2 = M_{2g/3}/M_{g/3}^2 - 1", sp.simplify(lhs - rhs) == 0)

    # (5) d32/d10
    lhs = (E([x ** 3 for x in d]) / E([x ** 2 for x in d])) / E(d)
    rhs = Msym(g) / (Msym(2 * g / 3) * Msym(g / 3))
    check("C1c  d32/d10 = M_g/(M_{2g/3} M_{g/3})", sp.simplify(lhs - rhs) == 0)

    # (6) d43/d10
    lhs = (E([x ** 4 for x in d]) / E([x ** 3 for x in d])) / E(d)
    rhs = Msym(4 * g / 3) / (Msym(g) * Msym(g / 3))
    check("C1d  d43/d10 = M_{4g/3}/(M_g M_{g/3})", sp.simplify(lhs - rhs) == 0)

    # (7) Lambda = total effective volume / total solid volume
    v_eff = sum(ni * (sp.pi / 6) * di ** 3 for ni, di in zip(ns, d))
    v_sol = sum(ni * kk * (sp.pi / 6) * d1s ** 3 for ni, kk in zip(ns, ks))
    check("C1e  Lambda = M_gamma / M_1",
          sp.simplify(v_eff / v_sol - Msym(g) / Msym(sp.Integer(1))) == 0)

    # the identity that makes gamma=1 exact for EVERY distribution
    check("C1f  Lambda == 1 identically at gamma = 1 (any n)",
          sp.simplify((v_eff / v_sol).subs(g, 1) - 1) == 0)

    # ---------------------------------------------------------------------------- C2
    say("\n=== C2  LIMITING CASES ===")
    k = [1, 2, 3, 4]
    n = [0.55, 0.25, 0.13, 0.07]

    r = closed_forms(n, k, 1.0)
    check("C2a  gamma = 1 (d_f = 3, compact): Lambda = 1 exactly",
          abs(r["Lambda"] - 1.0) < 1e-14, f"Lambda = {r['Lambda']:.16f}")

    r = closed_forms([1.0], [1], 1.7)
    check("C2b  monomers only: CV = 0 and Lambda = 1 at ANY gamma",
          abs(r["cv"]) < 1e-14 and abs(r["Lambda"] - 1) < 1e-14,
          f"CV = {r['cv']:.2e}, Lambda = {r['Lambda']:.16f}")

    r = closed_forms([1.0], [7], 1.7)
    exp_lam = 7.0 ** (1.7 - 1.0)
    check("C2c  ALL particles are 7-mers: CV = 0, Lambda = K^(gamma-1)",
          abs(r["cv"]) < 1e-12 and abs(r["Lambda"] - exp_lam) < 1e-12,
          f"Lambda = {r['Lambda']:.10f} vs K^(g-1) = {exp_lam:.10f}")

    r = closed_forms([0.5, 0.5], [1, 2], 3.0)
    # gamma = 3 <=> d_f = 1, a LINEAR chain: d_k = k d1
    check("C2d  gamma = 3 (d_f = 1, linear chain): d_k = k*d1 so d10/d1 = <k>",
          abs(r["d10_over_d1"] - 1.5) < 1e-14, f"d10/d1 = {r['d10_over_d1']}")

    # a monodisperse SPHERE pack must reproduce psd.py's own lognormal CV=0 identities
    from psd import moment_ratio_d32_d10, moment_ratio_d43_d10
    r = closed_forms([1.0], [1], 1.0)
    check("C2e  monodisperse: d32/d10 = d43/d10 = 1, matching psd.py at CV = 0",
          abs(r["d32_over_d10"] - moment_ratio_d32_d10(0.0)) < 1e-14
          and abs(r["d43_over_d10"] - moment_ratio_d43_d10(0.0)) < 1e-14)

    # ---------------------------------------------------------------------------- C3
    say("\n=== C3  DIMENSIONS ===")
    try:
        import pint
        u = pint.UnitRegistry()
        d1q = 0.8e-6 * u.m
        gam = 1.7
        kq = 3.0
        dk = (kq ** (gam / 3.0)) * d1q
        vk = (np.pi / 6.0) * dk ** 3
        vsol = kq * (np.pi / 6.0) * d1q ** 3
        lam = (vk / vsol).to_base_units()
        ok = (dk.dimensionality == u.m.dimensionality
              and vk.dimensionality == (u.m ** 3).dimensionality
              and lam.dimensionless)
        check("C3a  d_k is a LENGTH, v_k a VOLUME, Lambda DIMENSIONLESS (pint)", ok,
              f"d_k = {dk.to(u.um):~P}, Lambda = {lam:~P}")
    except ImportError:
        check("C3a  pint available", False, "pint not installed -- dimension check UNRUN")
    say("   by hand: k is a pure COUNT, gamma a pure EXPONENT, so k^(gamma/3) is")
    say("   dimensionless and every d_k carries d1's metres exactly once. M_s is a sum of")
    say("   (dimensionless number fraction) x (dimensionless k^s): dimensionless. CV, the")
    say("   two moment RATIOS and Lambda are ratios of like quantities: dimensionless.")

    # ---------------------------------------------------------------------------- C4
    say("\n=== C4  NUMERIC SPOT-CHECK vs an independent particle-array route (C6 anchor) ===")
    cases = [
        ("DLCA-like, mostly singlets", [0.70, 0.20, 0.07, 0.03], [1, 2, 3, 4], 1.667),
        ("RLCA-like, heavier tail   ", [0.45, 0.28, 0.17, 0.10], [1, 2, 3, 4], 1.429),
        ("compact doublets only     ", [0.60, 0.40, 0.00, 0.00], [1, 2, 3, 4], 1.000),
        ("strong aggregation        ", [0.25, 0.25, 0.25, 0.25], [1, 2, 3, 4], 2.000),
    ]
    say(f"   {'case':28s}{'gamma':>7}{'quantity':>16}{'closed form':>14}{'particles':>14}{'rel err':>11}")
    worst = 0.0
    for label, n, k, gam in cases:
        cf = closed_forms(n, k, gam)
        dp = direct_from_particles(n, k, gam)
        for q in ("d10_over_d1", "cv", "d32_over_d10", "d43_over_d10", "Lambda"):
            a, b = cf[q], dp[q]
            rel = abs(a - b) / max(abs(a), 1e-30)
            worst = max(worst, rel)
            say(f"   {label:28s}{gam:7.3f}{q:>16}{a:14.8f}{b:14.8f}{rel:11.2e}")
    check("C4  closed forms reproduce the particle-array route on 4 cases x 5 quantities",
          worst < 3e-3, f"worst relative error {worst:.2e} (Monte-Carlo, 2e6 particles)")

    # ---------------------------------------------------------------------------- C5
    say("\n=== C5  Lambda >= 1, and the EXISTENCE of a critical gamma* ===")
    rng = np.random.default_rng(7)
    bad = 0
    for _ in range(20000):
        K = rng.integers(1, 9)
        kk = np.arange(1, K + 1)
        nn = rng.dirichlet(np.ones(K))
        gam = float(rng.uniform(1.0, 3.0))
        if closed_forms(nn, kk, gam)["Lambda"] < 1.0 - 1e-12:
            bad += 1
    check("C5a  Lambda >= 1 on 20000 random (n, K, gamma in [1,3]) draws", bad == 0,
          f"{bad} violations")

    # G is bounded: phi_J < 1 always, so G = phi_J(CV)/phi_J_mono <= 1/phi_J_mono.
    PHI_J_MONO = 0.5748          # this project's own lognormal CV=0, mu_p=0.3, matched rung
    G_MAX = 1.0 / PHI_J_MONO
    n, k = [0.70, 0.20, 0.07, 0.03], [1, 2, 3, 4]
    gams = np.linspace(1.0, 3.0, 2001)
    lam = np.array([closed_forms(n, k, gm)["Lambda"] for gm in gams])
    crossed = np.any(lam > G_MAX)
    gstar = float(gams[np.argmax(lam > G_MAX)]) if crossed else float("nan")
    check("C5b  Lambda is monotone increasing in gamma", bool(np.all(np.diff(lam) > 0)))
    check("C5c  Lambda exceeds the ABSOLUTE ceiling on G before gamma = 3", crossed,
          f"G_max = 1/phi_J_mono = {G_MAX:.4f}, crossed at gamma = {gstar:.4f} "
          f"(d_f = {3.0/gstar:.4f})")
    say(f"   so a critical gamma* EXISTS for this population: below d_f = {3.0/gstar:.3f}")
    say("   aggregation must promote jamming, whatever phi_J(CV) does, because Lambda")
    say("   diverges while G cannot exceed 1/phi_J_mono.")

    # ---------------------------------------------------------------------------- done
    say("\n" + "=" * 92)
    if FAIL:
        say(f"FAILED {len(FAIL)} check(s):")
        for f in FAIL:
            say(f"   - {f}")
    else:
        say("ALL CHECKS PASSED")
    with open(os.path.join(HERE, "multiplet_psd_2026-09-22_check.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
