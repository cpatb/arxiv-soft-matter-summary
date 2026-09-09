"""drop_size_test.py — the drainage length, testable with a different nozzle.

THE PROBLEM THIS ATTACKS
------------------------
`derivations/dilatancy_sign_mechanism.md` §6 calls the drainage length `ell` the largest
single uncertainty in the lab's use of the inertial-dilatant number: `ID` depends on it
and nobody has measured it. §6 proposes a PARTICLE-SIZE ladder to settle it, which needs
new particles.

This module points out that **drop diameter does the same job and needs only a different
nozzle** -- and does it better, because the three readings of `ell` predict the boundary
velocity moving in three DIFFERENT DIRECTIONS, not merely by different amounts.

THE ALGEBRA, all three cases
----------------------------
The dilatant pressure comes from a Darcy balance, `grad^2 P = (mu_w/kappa) gammadot
tan(psi)`, so `P ~ (mu_w/kappa) * gammadot * dphi * L_p^2` where `L_p` is the length the
pressure gradient develops over, and `gammadot` is whatever rate the deformation imposes.
The dilatant FORCE is that pressure over the contact area, `~ D0^2`. The inertial force
is `~ rho U0^2 D0^2`. So

    ID = F_i/F_d ~ rho U0^2 / [ (mu_w/kappa) gammadot dphi L_p^2 ]

and everything turns on what `gammadot` and `L_p` are:

  (M) MACROSCOPIC  -- both set by the drop:  gammadot = U0/D0,  L_p = D0
      ID ~ rho kappa U0 / (A mu_w dphi D0)                      -> ID ~ 1/D0
      This is [M26] Eq. 6 verbatim, and it gives  U0* ~ D0^(+1).

  (B) BAND DRAINAGE -- shear at the drop scale, drainage in a grain-scale band:
      gammadot = U0/D0,  L_p = ell (fixed, ~10 grain diameters)
      ID ~ rho kappa U0 D0 / (A mu_w dphi ell^2)                -> ID ~ D0
      giving  U0* ~ D0^(-1).   OPPOSITE SIGN.

  (L) FULLY LOCALISED -- both set by the band:  gammadot = U0/ell,  L_p = ell
      ID ~ rho kappa U0 / (A mu_w dphi ell)                     -> ID independent of D0
      giving  U0* ~ D0^(0).    NO DEPENDENCE.

⚠ Case (B) is the one the derivation's §6 actually wrote down (it fixed `gammadot` at
`U0/D0` and let `ell` vary), and its `Pi ~ ell^2/D0` is reproduced here. Case (L) was not
distinguished from it there. They differ, and this module keeps them apart.

WHY THIS IS WORTH DOING
-----------------------
A factor-2 change in drop diameter moves the predicted boundary velocity by 2x UP, 2x
DOWN, or not at all. Three readings, three directions, one nozzle change. No competing
mechanism built on suspension VISCOSITY predicts any boundary in (U0, D0) at all, so a
null result is informative too.

    python drop_size_test.py
"""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "impact_stress")))

from inertial_dilatant import ID, DilatantSystem  # noqa: E402

ELL_GRAINS = 10.0     # the band thickness in grain diameters, for cases (B) and (L).
                      # O(10) is the standard granular shear-band width and it is used
                      # here as an ILLUSTRATION; the point of the test is the EXPONENT,
                      # which does not depend on this number at all.


def system(name, a, phi, phi_c, rho_p):
    rho = rho_p * phi + 1000.0 * (1 - phi)
    return DilatantSystem(name=name, phi=phi, phi_c=phi_c, a=a, rho=rho)


def id_macroscopic(s, u0, d0):
    """(M) [M26] Eq. 6 verbatim."""
    return ID(s, u0, d0)


def id_band(s, u0, d0, ell):
    """(B) shear at the drop scale, drainage over a fixed band: ID = ID_M * (D0/ell)^2."""
    return ID(s, u0, d0) * (d0 / ell) ** 2


def id_localised(s, u0, d0, ell):
    """(L) both at the band scale: ID = ID_M * (D0/ell), i.e. D0-independent."""
    return ID(s, u0, d0) * (d0 / ell)


def u_star(fn, s, d0, target, ell=None, lo=1e-4, hi=1e6):
    """Bisect for the U0 at which `fn` reaches `target`. ID is linear in U0 in all
    three cases, so this is exact up to the tolerance and needs no fitting."""
    def f(u):
        return (fn(s, u, d0) if ell is None else fn(s, u, d0, ell)) - target
    if f(lo) * f(hi) > 0:
        return float("nan")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def main() -> int:
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    # [M26]'s own cornstarch is the system where the boundary is REACHABLE, so it is the
    # one this test can actually be run on. Their transition sits inside 0.5-7.2 m/s at
    # D0 = 4.0 mm; ID there is the target.
    cs = system("cornstarch [M26]", 10.0e-6, 0.42, 0.37, 1550.0)
    si = system("silica (lab)", 0.8e-6, 0.44, 0.42, 2000.0)
    d0_ref = 4.0e-3
    target = id_macroscopic(cs, 7.2, d0_ref)

    p("=== THE DROP-SIZE TEST OF THE DRAINAGE LENGTH ===")
    p("  🔴 CALIBRATION, and the first version of this script got it wrong. Each reading")
    p("  is a DIFFERENT force ratio, so each has its OWN boundary value; applying (M)'s")
    p("  number to (B) and (L) made both look as though cornstarch were always inertial,")
    p("  which was an artefact of the mis-calibration and not a refutation of anything.")
    p("  Every reading is now pinned to the SAME observation -- [M26] see the transition")
    p(f"  at the top of their range, U0 = 7.2 m/s at D0 = {d0_ref * 1e3:.1f} mm -- so all")
    p("  three agree there BY CONSTRUCTION and the test is entirely about how they MOVE.")
    p()
    ell_cs = ELL_GRAINS * cs.a
    tgt = dict(
        M=id_macroscopic(cs, 7.2, d0_ref),
        B=id_band(cs, 7.2, d0_ref, ell_cs),
        L=id_localised(cs, 7.2, d0_ref, ell_cs),
    )
    p(f"  boundary value per reading: (M) {tgt['M']:.4e}   (B) {tgt['B']:.4e}   "
      f"(L) {tgt['L']:.4e}")
    p()
    p(f"  {'D0 [mm]':>8s} {'U0* (M)':>10s} {'U0* (B)':>10s} {'U0* (L)':>10s}   "
      f"{'(M)/ref':>8s} {'(B)/ref':>8s} {'(L)/ref':>8s}")
    for d0_mm in (1.0, 2.0, 4.0, 6.0, 8.0):
        d0 = d0_mm * 1e-3
        um = u_star(id_macroscopic, cs, d0, tgt["M"])
        ub = u_star(id_band, cs, d0, tgt["B"], ell=ell_cs)
        ul = u_star(id_localised, cs, d0, tgt["L"], ell=ell_cs)
        p(f"  {d0_mm:8.2f} {um:10.2f} {ub:10.2f} {ul:10.2f}   "
          f"{um / 7.2:8.3f} {ub / 7.2:8.3f} {ul / 7.2:8.3f}")
    p()
    p("  🔑 THE THREE READINGS MOVE IN THREE DIRECTIONS. Halving the drop diameter from")
    p("  4 to 2 mm HALVES the boundary velocity under (M), DOUBLES it under (B), and")
    p("  leaves it unchanged under (L) -- a 4x separation between the extremes over one")
    p("  nozzle change. The SIGN alone decides it, so the test does not depend on the")
    p("  O(1) band thickness, on the ID prefactor, or on where the boundary sits.")
    p()
    p("  ⚠ (M) and (B) both predict the transition LEAVES the accessible range at one")
    p("  end of the ladder -- (M) at small D0, (B) at large -- so the ladder should be")
    p("  centred on 4 mm and kept inside about 2-6 mm for the transition to stay")
    p("  observable under every reading.")
    p()
    p("  A competing mechanism built on the suspension VISCOSITY (phi_m rising with CV,")
    p("  eta ~ (phi_m - phi)^-2) predicts NO boundary in (U0, D0) at all, so a null")
    p("  result is informative rather than a failed experiment.")
    p()
    p("=== IS IT RUNNABLE ON THE LAB'S OWN SILICA? NO -- and that is worth knowing ===")
    for d0_mm in (2.5,):
        d0 = d0_mm * 1e-3
        um = u_star(id_macroscopic, si, d0, target)
        p(f"  silica, D0 = {d0_mm} mm: U0* = {um:.0f} m/s under reading (M) -- "
          f"unreachable on any bench tower.")
    p("  So the drop-size test belongs on CORNSTARCH (or on a larger particle), where")
    p("  the boundary is already inside the accessible velocity range. That is a")
    p("  feature: it tests [M26]'s construction on [M26]'s own material, and it is the")
    p("  cheapest available attack on the drainage length.")
    p()
    p("  ⚠ WHAT WOULD MAKE THIS INVALID: if changing D0 also changes the impact regime")
    p("  in some other way that moves the observed transition -- Weber and Reynolds both")
    p("  scale with D0, so the SPREADING behaviour changes even with no rheology change.")
    p("  The observable must therefore be the REGIME (solid-like vs re-entrant liquid),")
    p("  not the spread factor, and the control is a Newtonian liquid of matched")
    p("  viscosity and density run over the same D0 ladder.")

    out = os.path.join(HERE, "out", "drop_size_test.txt")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
