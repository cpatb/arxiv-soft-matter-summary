"""dilatancy_sign_mechanism_numbers.py — every number quoted in the derivation's prose.

The repo rule is that a number in a document is printed by a script, not typed into the
prose. This prints the §1 dimensionless-group table, the §4 channel table, the §6
drainage-length contrast and the §7 ceilings, to
`derivations/dilatancy_sign_mechanism_numbers.txt`.

    python derivations/dilatancy_sign_mechanism_numbers.py

Inputs are the lab's own values, each named where it comes from. `phi_c` and its CV
slope are the two that are NOT measurements yet and they are labelled PLACEHOLDER at
every point of use.
"""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "impact_stress")))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))

from inertial_dilatant import ID, P_ATM, DilatantSystem  # noqa: E402
from psd import moment_ratio_d32_d10  # noqa: E402

# --- the lab's operating point ------------------------------------------------------
PHI = 0.44           # Paper 1 working solid fraction
PHI_C0 = 0.42        # PLACEHOLDER, impact_stress.lab_silica_dilatant default
D10 = 800e-9         # 800 nm Stoeber silica
U0 = 3.0             # m/s, bench tower
D0 = 2.5e-3          # m
MU_W = 1.0e-3        # Pa s, water
GAMMA_LV = 0.072     # N/m
RHO_P = 2000.0       # kg/m^3, silica-in-water suspension density model in that module
N_SHEARBAND = 10.0   # shear-band thickness in grain diameters; O(10) is the standard
                     # granular value and it is used here ONLY as an illustration of
                     # the sensitivity, never as a measurement.

# PLACEHOLDER critical-state slope: the ISOTROPIC phi_J slope, +0.0623 over CV 0 -> 0.5
# (out/phase2_summary.csv, lognormal, mu_p = 0.3, N = 1500, three seeds).
S_PLACEHOLDER = 0.0623 / 0.5


def rho_susp(phi: float) -> float:
    return RHO_P * phi + 1000.0 * (1 - phi)


def main() -> int:
    lines: list[str] = []

    def p(msg: str = "") -> None:
        print(msg)
        lines.append(msg)

    rho = rho_susp(PHI)
    d32 = D10 * moment_ratio_d32_d10(0.0)
    sysm = DilatantSystem(name="lab silica", phi=PHI, phi_c=PHI_C0, a=d32, rho=rho)
    id_d0 = ID(sysm, U0, D0)

    p("=== SECTION 1: dimensionless groups at the lab's operating point ===")
    p(f"  phi = {PHI}, d10 = {D10 * 1e9:.0f} nm, U0 = {U0} m/s, D0 = {D0 * 1e3:.1f} mm")
    p(f"  rho_suspension            = {rho:.1f} kg/m^3")
    p(f"  Re  = rho U0 D0 / mu_w    = {rho * U0 * D0 / MU_W:.3e}")
    p(f"  We  = rho U0^2 D0 / gamma = {rho * U0**2 * D0 / GAMMA_LV:.3e}")
    p(f"  St  = rho_p a^2 U0/(9 mu_w D0) "
      f"= {RHO_P * d32**2 * U0 / (9 * MU_W * D0):.3e}")
    p(f"  kappa (Carman-Kozeny)     = {sysm.kappa:.3e} m^2")
    p(f"  ID  (ell = D0)            = {id_d0:.3e}")
    p(f"  Pi  (ell = D0) = 1/ID     = {1.0 / id_d0:.3e}   -> UNDRAINED")
    ell_grain = N_SHEARBAND * d32
    pi_grain = (1.0 / id_d0) * (ell_grain / D0) ** 2
    p(f"  Pi  (ell = {N_SHEARBAND:.0f} d32 = {ell_grain * 1e6:.1f} um) "
      f"= {pi_grain:.3e}   -> {'DRAINED' if pi_grain < 1 else 'UNDRAINED'}")
    p(f"  ratio of the two Pi       = {(1.0 / id_d0) / pi_grain:.3e} "
      f"= (D0/ell)^2 = {(D0 / ell_grain) ** 2:.3e}")
    p()

    p("=== SECTION 4: the two CV channels, at the PLACEHOLDER critical-state slope ===")
    p(f"  s = dphi_c/dCV = {S_PLACEHOLDER:.4f}  <-- PLACEHOLDER (isotropic phi_J slope)")
    cv_star = (PHI - PHI_C0) / S_PLACEHOLDER
    p(f"  CV* = (phi - phi_c0)/s = {cv_star:.4f}   <-- PROVISIONAL, do not quote")
    p("   CV     dphi        kappa-chan   phic-chan    ratio    d32/d10   ID/ID(CV=0)")
    id0 = None
    for cv in (0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50):
        r32 = moment_ratio_d32_d10(cv)
        dphi = PHI - (PHI_C0 + S_PLACEHOLDER * cv)
        kch = 8 * cv / (1 + cv**2)
        s_cv = DilatantSystem(name="s", phi=PHI, phi_c=PHI_C0 + S_PLACEHOLDER * cv,
                              a=D10 * r32, rho=rho)
        idv = ID(s_cv, U0, D0)
        if cv == 0.0:
            id0 = idv
        if dphi > 0:
            pch = S_PLACEHOLDER / dphi
            ratio = f"{pch / kch:8.2f}" if kch > 0 else "     inf"
            p(f"  {cv:4.2f}  {dphi:+9.5f}   {kch:9.4f}   {pch:9.4f}  {ratio}  "
              f"{r32:8.4f}  {idv / id0:11.4g}")
        else:
            p(f"  {cv:4.2f}  {dphi:+9.5f}   {kch:9.4f}         --        --  "
              f"{r32:8.4f}   CONTRACTIVE: Eq. 6 does not apply")
    p()

    p("=== SECTION 7: the ceiling on the dilatant pore pressure ===")
    p("  The binding ceiling is the NETWORK-MODULUS one of PORE_PHYSICS.md Eq. (4.5),")
    p("  p_dil_max = K_p dphi / phi -- a statement about how much dilation strain the")
    p("  packing HAS, not about what the fluid can survive. The cavitation cap that an")
    p("  earlier brief proposed was WITHDRAWN there (PORE_PHYSICS 4.4) for four reasons,")
    p("  the strongest being that it removes all particle-size dependence.")
    import pore_relaxation as PR
    ps = PR.silica_poro(phi=PHI, phi_c=PHI_C0)
    dp = PR.dilatant_pressure_darcy(ps.sys, L_d=D0, gammadot=U0 / D0)
    ceil = PR.dilatant_pressure_ceiling(ps.sys, ps.K_p)
    clo = PR.dilatant_pressure_ceiling(ps.sys, ps.K_p_low)
    chi = PR.dilatant_pressure_ceiling(ps.sys, ps.K_p_high)
    p(f"  uncapped Darcy, Eq. (4), ell = D0    = {dp:.4g} Pa")
    p(f"  network-modulus ceiling (4.5)        = {ceil:.4g} Pa "
      f"[bracket {clo:.4g} to {chi:.4g} Pa over lambda_D 30 nm -> 1 nm]")
    p(f"  overshoot of the uncapped expression = {dp / ceil:.4g}x  "
      f"(SIX decades, not the ~460x a cavitation cap would give)")
    p(f"  withdrawn cavitation cap (1 atm)     = {P_ATM:.4g} Pa  <-- NOT USED")
    p(f"  impact stress rho U0^2               = {rho * U0**2:.4g} Pa")
    p(f"  ceiling / impact stress              = {ceil / (rho * U0**2):.4g}")
    p("  ANCHOR-FREE RATIO (PORE_PHYSICS 9.4): p_dil_max/P_conf = (a/6 lambda_D)(dphi/phi)")
    p("  -- P_conf cancels, so this survives the withdrawal of si_fit's sigma*, which is")
    p("  the ONLY thing the absolute Pa figures above depend on.")
    for lam_nm, salt in ((30.0, "nominal 0 mM"), (10.0, "1 mM"), (1.0, "100 mM")):
        ratio = (D10 / (6 * lam_nm * 1e-9)) * ((PHI - PHI_C0) / PHI)
        p(f"    lambda_D = {lam_nm:4.0f} nm ({salt:11s}): "
          f"p_dil_max/P_conf = {ratio:8.3f}  "
          f"{'AVAILABLE as strengthening' if ratio > 1 else 'NOT available'}")
    p("  KEY STRUCTURAL POINT: the ceiling contains NO permeability. Where it binds --")
    p("  and for this system it binds by six decades -- the M1 channel is not merely")
    p("  sub-dominant, it is ABSENT, and the only PSD dependence left is dphi(CV) and")
    p("  whatever moment of the size distribution K_p reads (unmeasured; K_p ~ a for the")
    p("  double-layer contact, so it is a FIRST moment, not d32).")
    p()

    p("=== SECTION 6: the drainage-length ambiguity, as a size ladder ===")
    p("   a[um]   ID(ell=D0)    ID(ell=10 d32)   ratio to a=0.8um, each reading")
    base = {}
    for a_um in (0.8, 2.5, 8.0):
        a = a_um * 1e-6
        s_a = DilatantSystem(name="a", phi=PHI, phi_c=PHI_C0, a=a, rho=rho)
        id_macro = ID(s_a, U0, D0)
        id_local = id_macro * (D0 / (N_SHEARBAND * a)) ** 2
        if a_um == 0.8:
            base = {"m": id_macro, "l": id_local}
        p(f"  {a_um:5.1f}   {id_macro:.4e}    {id_local:.4e}    "
          f"macro {id_macro / base['m']:8.2f}x   local {id_local / base['l']:8.2f}x")
    p("  -> macroscopic drainage gives ID ~ a^2; grain-scale drainage gives ID ~ a^0.")

    out = os.path.join(HERE, "dilatancy_sign_mechanism_numbers.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
