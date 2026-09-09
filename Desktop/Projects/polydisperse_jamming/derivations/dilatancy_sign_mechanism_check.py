"""dilatancy_sign_mechanism_check.py — the four checks for derivations/dilatancy_sign_mechanism.md

Run:  python derivations/dilatancy_sign_mechanism_check.py

Checks, in the order the derive-and-check skill requires:

  C1 SymPy    the closed form for the pore-pressure ratio Pi reproduces, term by term,
              the Reynolds-Darcy balance it was derived from, and equals 1/ID of
              [M26] Eq. 6 when the drainage length is the drop diameter.
  C2 Limits   seven named limits, each recovered or explained.
  C3 Dims     every symbol carried through `pint`; Pi must come out dimensionless.
  C4 Numeric  the closed form evaluated against the LAB'S OWN independent
              implementation (impact_stress.inertial_dilatant.ID) at 5 parameter
              points spanning the experimental range. That module was written for a
              different purpose by a different pass and is not imported by the
              derivation, so it is an anchor outside the code under test.
  C5 Channels the two CV channels (permeability and critical state) and their ratio,
              including the analytic statement that the permeability channel vanishes
              at CV -> 0 while the critical-state channel does not.

No number in this file is a fitted or plausible-looking value: every constant is
either an identity, a published value with its source named, or a parameter swept.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
IMPACT_STRESS = os.path.abspath(os.path.join(HERE, "..", "..", "impact_stress"))

FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


# =====================================================================================
# C1  SymPy: the closed form follows from the Reynolds-Darcy balance
# =====================================================================================
def c1_symbolic() -> None:
    print("\n=== C1  SymPy ===")
    (mu_w, alpha, gdot, dphi, ell, kappa, rho, U0, D0, A) = sp.symbols(
        "mu_w alpha gammadot Delta_phi ell kappa rho U0 D0 A", positive=True
    )

    # (1) Reynolds dilatancy: the volumetric strain rate the skeleton demands.
    eps_v = alpha * gdot * dphi                       # [1/s]
    # (2) volume flux per unit area the fluid must supply over the drainage length
    q = eps_v * ell                                   # [m/s]
    # (3) Darcy, with the pressure drop developing over the same length
    #     u_D = (kappa/mu_w) * (dp/ell)   ==>   dp = q * mu_w * ell / kappa
    dp = sp.solve(sp.Eq(kappa / mu_w * (sp.Symbol("dp", positive=True) / ell), q),
                  sp.Symbol("dp", positive=True))[0]
    dp_expected = alpha * gdot * dphi * mu_w * ell**2 / kappa
    check("C1.1 pore-pressure magnitude from Darcy + Reynolds",
          sp.simplify(dp - dp_expected) == 0, f"dp = {sp.simplify(dp)}")

    # (4) nondimensionalise on the inertial stress rho U0^2, with gdot ~ U0/D0
    Pi = (dp_expected.subs(gdot, U0 / D0)) / (rho * U0**2)
    Pi_expected = alpha * mu_w * dphi * ell**2 / (rho * kappa * U0 * D0)
    check("C1.2 Pi = alpha mu_w dphi ell^2 / (rho kappa U0 D0)",
          sp.simplify(Pi - Pi_expected) == 0, f"Pi = {sp.simplify(Pi)}")

    # (5) with ell = D0 and alpha -> A, Pi is exactly 1/ID of [M26] Eq. 6.
    ID_M26 = rho * kappa * U0 / (A * mu_w * dphi * D0)
    Pi_at_D0 = Pi_expected.subs({ell: D0, alpha: A})
    check("C1.3 Pi(ell = D0) == 1 / ID  ([M26] Eq. 6)",
          sp.simplify(Pi_at_D0 - 1 / ID_M26) == 0,
          f"Pi(ell=D0) = {sp.simplify(Pi_at_D0)}")

    # (6) the drainage-length rescaling, stated as an identity
    ratio = sp.simplify(Pi_expected / Pi_expected.subs(ell, D0))
    check("C1.4 Pi(ell)/Pi(D0) = (ell/D0)^2", sp.simplify(ratio - (ell / D0) ** 2) == 0,
          f"ratio = {ratio}")


# =====================================================================================
# C2  Limits
# =====================================================================================
def c2_limits() -> None:
    print("\n=== C2  limiting cases ===")
    (mu_w, alpha, dphi, ell, kappa, rho, U0, D0) = sp.symbols(
        "mu_w alpha Delta_phi ell kappa rho U0 D0", positive=True
    )
    Pi = alpha * mu_w * dphi * ell**2 / (rho * kappa * U0 * D0)

    lims = [
        ("mu_w -> 0 (inviscid pore fluid): no pore pressure, purely inertial",
         sp.limit(Pi, mu_w, 0), 0),
        ("kappa -> oo (freely draining): no pore pressure",
         sp.limit(Pi, kappa, sp.oo), 0),
        ("kappa -> 0 (impermeable): pore pressure unbounded -> the CEILING binds",
         sp.limit(Pi, kappa, 0), sp.oo),
        ("U0 -> oo (fast impact): dilatancy cannot keep up, inertia wins",
         sp.limit(Pi, U0, 0 + sp.oo), 0),
        ("U0 -> 0 (quasi-static): pore pressure dominates the inertial scale",
         sp.limit(Pi, U0, 0), sp.oo),
        ("Delta_phi -> 0 (at the critical state): NO dilatancy, no pore pressure",
         sp.limit(Pi, dphi, 0), 0),
        ("ell -> 0 (grain-scale drainage): pore pressure vanishes",
         sp.limit(Pi, ell, 0), 0),
    ]
    for name, got, want in lims:
        check(f"C2 {name}", sp.simplify(got - want) == 0 if want != sp.oo
              else got == sp.oo, f"-> {got}")

    # The sign statement is not a limit of Pi (which is written for |dphi|); it is the
    # sign of the pressure change, and it is the whole mechanism. Stated explicitly:
    print("     note: Pi is written with |Delta_phi|. The SIGN of the pore-pressure "
          "change is sign(phi - phi_c):")
    print("           phi > phi_c  dilates -> fluid drawn IN  -> SUCTION      -> "
          "effective stress RISES -> solid-like")
    print("           phi < phi_c  contracts -> fluid driven OUT -> OVERPRESSURE -> "
          "effective stress FALLS -> liquefied")


# =====================================================================================
# C3  Dimensions, with pint
# =====================================================================================
def c3_dimensions() -> None:
    print("\n=== C3  dimensions (pint) ===")
    import pint

    u = pint.UnitRegistry()
    mu_w = 1.0e-3 * u.Pa * u.s
    alpha = 1.0 * u.dimensionless
    dphi = 0.02 * u.dimensionless
    ell = 2.5e-3 * u.m
    kappa = 1.0e-15 * u.m**2
    rho = 1440.0 * u.kg / u.m**3
    U0 = 3.0 * u.m / u.s
    D0 = 2.5e-3 * u.m

    Pi = alpha * mu_w * dphi * ell**2 / (rho * kappa * U0 * D0)
    check("C3.1 Pi is dimensionless", Pi.dimensionality == u.dimensionless.dimensionality,
          f"Pi = {Pi.to_base_units():~P}")

    # Carman-Kozeny, and the Sauter-mean identity
    phi = 0.44
    d32 = 8.0e-7 * u.m
    kappa_ck = (1 - phi) ** 3 * d32**2 / (180 * phi**2)
    check("C3.2 Carman-Kozeny kappa has dimensions of area",
          kappa_ck.dimensionality == (u.m**2).dimensionality,
          f"kappa = {kappa_ck.to('m**2'):~.3e}")

    # the pore-pressure magnitude itself
    dp = alpha * mu_w * dphi * (U0 / D0) * ell**2 / kappa_ck
    check("C3.3 the pore pressure has dimensions of pressure",
          dp.dimensionality == (u.Pa).dimensionality, f"dp = {dp.to('kPa'):~.4g}")


# =====================================================================================
# C4  Numeric, against the lab's independent implementation
# =====================================================================================
def c4_numeric() -> None:
    print("\n=== C4  numeric, against impact_stress.inertial_dilatant (independent) ===")
    if IMPACT_STRESS not in sys.path:
        sys.path.insert(0, IMPACT_STRESS)
    try:
        from inertial_dilatant import ID, DilatantSystem  # noqa: N811
    except Exception as exc:                              # noqa: BLE001
        check("C4 import impact_stress.inertial_dilatant", False, repr(exc))
        return

    # Five points spanning the lab's range: two particle sizes, two velocities, and
    # [M26]'s own cornstarch. dphi is held at the module default so the comparison is
    # of the ALGEBRA, not of phi_c.
    pts = [
        ("silica 800 nm, U0=3",  800e-9, 0.44, 0.42, 3.0, 2.5e-3, 2000.0),
        ("silica 800 nm, U0=1",  800e-9, 0.44, 0.42, 1.0, 2.5e-3, 2000.0),
        ("silica 1.5 um, U0=3",  1.5e-6, 0.44, 0.42, 3.0, 2.5e-3, 2000.0),
        ("silica 800 nm, phi=0.46", 800e-9, 0.46, 0.42, 3.0, 2.5e-3, 2000.0),
        ("cornstarch 10 um, U0=3", 10e-6, 0.42, 0.37, 3.0, 4.0e-3, 1550.0),
    ]
    worst = 0.0
    for name, a, phi, phi_c, U0, D0, rho_p in pts:
        rho = rho_p * phi + 1000.0 * (1 - phi)
        s = DilatantSystem(name=name, phi=phi, phi_c=phi_c, a=a, rho=rho)
        theirs = ID(s, U0, D0)
        # ours: Pi = A mu_w dphi D0 / (rho kappa U0) with ell = D0, so ID = 1/Pi
        kappa = (1 - phi) ** 3 * a**2 / (180 * phi**2)
        Pi = s.A * s.mu_w * (phi - phi_c) * D0 / (rho * kappa * U0)
        ours = 1.0 / Pi
        rel = abs(ours - theirs) / theirs
        worst = max(worst, rel)
        print(f"     {name:26s} ID_theirs={theirs:.6e}  ID_ours={ours:.6e}  "
              f"rel={rel:.2e}")
    check("C4 closed form agrees with the independent implementation to 1e-12",
          worst < 1e-12, f"worst relative difference {worst:.2e}")

    # A DISCRIMINATING negative case: the same call with a DIFFERENT dphi must NOT
    # agree, otherwise the comparison above is insensitive to the quantity the whole
    # mechanism turns on.
    s = DilatantSystem(name="neg", phi=0.44, phi_c=0.42, a=800e-9, rho=1440.0)
    s2 = DilatantSystem(name="neg2", phi=0.44, phi_c=0.43, a=800e-9, rho=1440.0)
    r = ID(s2, 3.0, 2.5e-3) / ID(s, 3.0, 2.5e-3)
    check("C4-neg halving dphi changes ID by more than 1.5x (the test can fail)",
          r > 1.5, f"ID(dphi=0.01)/ID(dphi=0.02) = {r:.4f}")


# =====================================================================================
# C5  The two CV channels
# =====================================================================================
def c5_channels() -> None:
    print("\n=== C5  the two CV channels ===")
    CV, phi, s, phic0 = sp.symbols("CV phi s phi_c0", positive=True)

    # lognormal at fixed NUMBER mean: d32/d10 = (1 + CV^2)^2  (psd.moment_ratio_d32_d10)
    d32_ratio = (1 + CV**2) ** 2
    kappa_channel = sp.simplify(sp.diff(sp.log(d32_ratio**2), CV))
    check("C5.1 d ln(kappa)/dCV = 8 CV/(1+CV^2)",
          sp.simplify(kappa_channel - 8 * CV / (1 + CV**2)) == 0,
          f"= {kappa_channel}")

    # critical-state channel, with phi_c(CV) linearised as phi_c0 + s*CV
    dphi_expr = phi - (phic0 + s * CV)
    phic_channel = sp.simplify(-sp.diff(sp.log(dphi_expr), CV))
    check("C5.2 d ln(1/dphi)/dCV = s/(phi - phi_c)",
          sp.simplify(phic_channel - s / dphi_expr) == 0, f"= {phic_channel}")

    # THE ANALYTIC CORE: the permeability channel is O(CV) and the critical-state
    # channel is O(1) as CV -> 0, so the leading CV dependence of ID is entirely the
    # packing one, for ANY non-zero s.
    lim_kappa = sp.limit(kappa_channel, CV, 0)
    lim_phic = sp.limit(phic_channel, CV, 0)
    check("C5.3 permeability channel VANISHES at CV -> 0", lim_kappa == 0,
          f"-> {lim_kappa}")
    check("C5.4 critical-state channel does NOT vanish at CV -> 0",
          sp.simplify(lim_phic - s / (phi - phic0)) == 0, f"-> {lim_phic}")

    # numeric ratio at the lab's numbers, with s left as a parameter
    R = sp.simplify(phic_channel / kappa_channel)
    print(f"     channel ratio R(CV) = {R}")
    s_val = 0.1246   # PLACEHOLDER: the ISOTROPIC phi_J slope, +0.0623 over CV 0->0.5
                     # (out/phase2_summary.csv, lognormal, mu_p = 0.3, N = 1500).
                     # phase8_criticalstate.py replaces it with the CRITICAL-STATE
                     # slope, which is the quantity the mechanism actually needs.
    cv_star = float((0.44 - 0.42) / s_val)
    print(f"     CV* (where phi_c(CV) = phi, i.e. the pole) = {cv_star:.4f} "
          f"at phi = 0.44, phi_c0 = 0.42, s = {s_val}")
    for cv_val in (0.05, 0.10, 0.15, 0.20, 0.50):
        dphi_val = 0.44 - (0.42 + s_val * cv_val)
        expr = float(R.subs({CV: cv_val, phi: 0.44, phic0: 0.42, s: s_val}))
        if dphi_val > 0:
            print(f"     CV={cv_val:.2f}: dphi={dphi_val:+.5f}  R = {expr:8.3f}  "
                  f"(critical-state channel is {expr:.1f}x the permeability one)")
        else:
            print(f"     CV={cv_val:.2f}: dphi={dphi_val:+.5f}  PAST THE POLE -- the "
                  f"packing is CONTRACTIVE, [M26] Eq. 6 does not apply, and R is "
                  f"not a ratio of two same-signed channels")


# =====================================================================================
# C6  The third term of Eq. (11): d ln A/dCV, which the prose asserts is small
# =====================================================================================
def c6_A_channel() -> None:
    """[M26]'s A depends on CV through tan(beta) = tan(beta0) + alpha*dphi.

    Step 7 of the derivation asserts this term is smaller than the critical-state one by
    a factor of order dphi. That was ASSERTED and not checked until a referee computed
    it; it is checked here.
    """
    print("")
    print("=== C6  the A channel of Eq. (11) ===")
    if IMPACT_STRESS not in sys.path:
        sys.path.insert(0, IMPACT_STRESS)
    from inertial_dilatant import DilatantSystem

    phi, phic0, s_slope = 0.44, 0.42, 0.1246
    worst = 0.0
    for cv in (0.0, 0.05, 0.10, 0.15):
        dphi = phi - (phic0 + s_slope * cv)
        h = 1e-6
        def lnA(c):
            d = phi - (phic0 + s_slope * c)
            return float(np.log(DilatantSystem(name="x", phi=phi, phi_c=phi - d,
                                               a=800e-9, rho=1440.0).A))
        dlnA = (lnA(cv + h) - lnA(cv - h)) / (2 * h)
        crit = s_slope / dphi
        ratio = abs(dlnA) / crit
        worst = max(worst, ratio)
        print(f"     CV={cv:.2f}  |dlnA/dCV|={abs(dlnA):.4f}  crit_term={crit:8.3f}  "
              f"ratio={ratio:.4f}")
    check("C6 the A channel is under 5 % of the critical-state channel everywhere",
          worst < 0.05, f"worst ratio {worst:.4f}")


# =====================================================================================
# C7  EXTERNAL ANCHOR: the same machinery against a MEASURED pore suction
# =====================================================================================
def c7_external_anchor() -> None:
    """Jerome, Vandenberghe & Forterre measured -1e4 to -1e5 Pa in 170 um glass beads.

    This is the only check here that compares against DATA rather than against algebra.
    It exercises the Reynolds + Darcy + network-modulus chain that Eq. (4) and section 7
    rest on, through the repo's own `pore_relaxation` implementation, on a system the
    derivation has no freedom to tune.
    """
    print("")
    print("=== C7  external anchor: JVF's measured glass-bead suction ===")
    if IMPACT_STRESS not in sys.path:
        sys.path.insert(0, IMPACT_STRESS)
    try:
        import pore_relaxation as PR
    except Exception as exc:                          # noqa: BLE001
        check("C7 import pore_relaxation", False, repr(exc))
        return
    g = PR.glass_beads_JVF()
    ceil = PR.dilatant_pressure_ceiling(g.sys, g.K_p)
    lo, hi = 1.0e4, 1.0e5           # [JVF16]'s reported range, Pa
    print(f"     predicted ceiling K_p dphi/phi = {ceil:.3e} Pa "
          f"(K_p = {g.K_p:.3e} Pa, dphi = {g.sys.dphi:.4f})")
    print(f"     measured range                 = {lo:.0e} to {hi:.0e} Pa")
    check("C7 the prediction is within one decade of the measured range, and ABOVE it "
          "(the affine modulus is an upper bound)",
          hi <= ceil <= 10 * hi, f"ratio to the top of the measured range "
                                 f"{ceil / hi:.3f}")
    # A DISCRIMINATING negative case has to be one the anchor could fail. Predicting
    # the SILICA ceiling with the glass-bead modulus must NOT land in the glass-bead
    # range, or C7 would pass for a routine that ignores its inputs.
    si = PR.silica_poro(phi=0.44, phi_c=0.42)
    ceil_si = PR.dilatant_pressure_ceiling(si.sys, si.K_p)
    check("C7-neg the same formula on 800 nm silica lands far OUTSIDE the glass-bead "
          "range, so C7 is sensitive to its inputs",
          ceil_si < lo / 100.0,
          f"silica ceiling {ceil_si:.3e} Pa against a measured bead range starting at "
          f"{lo:.0e} Pa -- a factor {lo / ceil_si:.3g}")
    # NOT a check: the withdrawn 1-atm cavitation cap happens to sit AT the top of the
    # measured bead range, so this dataset cannot tell the two ceilings apart. It is
    # recorded because a reader would otherwise assume C7 also refutes the cap.
    from inertial_dilatant import P_ATM
    print(f"     note: the withdrawn cavitation cap is {P_ATM:.3e} Pa, which is inside "
          f"the measured range -- this anchor does NOT discriminate the two ceilings;")
    print("           what does is that a fluid constant carries no particle size, so "
          "it gives the SAME number for 170 um beads and 800 nm silica.")


if __name__ == "__main__":
    c1_symbolic()
    c2_limits()
    c3_dimensions()
    c4_numeric()
    c5_channels()
    c6_A_channel()
    c7_external_anchor()
    print("\n" + ("ALL CHECKS PASSED" if not FAILURES
                  else f"{len(FAILURES)} FAILED: {FAILURES}"))
    raise SystemExit(1 if FAILURES else 0)
