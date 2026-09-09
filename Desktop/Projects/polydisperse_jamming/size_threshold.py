"""size_threshold.py — is there a critical PARTICLE SIZE for solid-like -> inertial?

THE HYPOTHESIS
--------------
The author's, 2026-09-08: *maybe there is simply a size threshold `a*`, above which an
impacting dense suspension stops behaving solid-like and starts behaving as the
inertial-dilatant (re-entrant liquid) mode.*

A size threshold is what a critical value of a size-dependent dimensionless group IS, so
the hypothesis is a family with one member per such group. This module enumerates them,
evaluates each **at each material's own operating point**, and reports how each moves
with the knobs the bench can turn.

⚠⚠ REWRITTEN 2026-09-08 after a referee pass. THE FIRST VERSION HAD A FATAL DEFECT and
its headline ("of seven candidates exactly one survives") is WITHDRAWN. It computed a
single `a*` at the LAB SILICA's operating point and then sorted CORNSTARCH with it --
but `a*` depends on `phi_c`, `D0` and `rho` and [M26]'s cornstarch has different values
for all three. Recomputed at cornstarch's own point, `a*` = 15.4 um and cornstarch falls
BELOW its own threshold, so the "survivor" does not separate the materials at all. This
is the project's own recurring lesson -- never grade a knob at one setting of the knob it
interacts with -- committed again.

WHAT REPLACES IT
----------------
Two statements, of which only the second is well posed:

  (A) CROSS-MATERIAL. Do not convert to a critical size at all. Evaluate each group at
      each material's OWN parameters and compare the NUMBERS. Whether one boundary value
      separates them is then a question about the numbers and carries no hidden
      operating point.

  (B) WITHIN ONE MATERIAL. A critical size IS well defined along a ladder that holds
      chemistry, `phi`, `phi_c`, `D0` and `U0` fixed and moves only `a` -- which is
      exactly the experiment being proposed (E-4). `a_star_within` computes it, and it
      is the only place `a*` is quoted.

    python size_threshold.py            # the group table and the boundary bracket
    python size_threshold.py --sweep    # how each threshold moves, and the CV routes
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from dataclasses import dataclass

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "impact_stress")))

from inertial_dilatant import ID, DilatantSystem  # noqa: E402

# --- constants, each with its origin -------------------------------------------------
MU_W = 1.0e-3          # Pa s, water
GAMMA_LV = 0.072       # N/m, suspension-air surface tension
K_B = 1.380649e-23     # J/K, exact (SI 2019)
T_LAB = 293.15         # K
G_ACCEL = 9.80665      # m/s^2
CHI_W = 4.55e-10       # 1/Pa, isothermal compressibility of water
                       # (impact_stress.inertial_dilatant.poroelastic_time_ratio)
R_THROAT_COEF = 0.155 / 2.0   # r_throat = 0.155 * R = 0.0775 * a  (diameter convention)
LAMBDA_D_LADDER = ((30.0e-9, "0 mM"), (10.0e-9, "1 mM"),
                   (3.0e-9, "10 mM"), (1.0e-9, "100 mM"))
LAMBDA_D = 30.0e-9     # nominal zero added salt

# ⚠⚠ PLACEHOLDER, and the sensitivity is large: the ISOTROPIC phi_J slope
# (+0.0623 over CV 0 -> 0.5, out/phase2_summary.csv, lognormal, mu_p = 0.3, N = 1500)
# standing in for the CRITICAL-STATE slope that phase8_criticalstate.py is measuring.
S_PHI_C_PLACEHOLDER = 0.0623 / 0.5


@dataclass(frozen=True)
class Material:
    """A material AT ITS OWN operating point. Every field is a parameter that at least
    one candidate threshold depends on, which is the whole reason this type exists."""
    name: str
    a: float            # particle DIAMETER [m]
    phi: float
    phi_c: float
    d0: float           # drop diameter [m]
    rho_p: float        # particle density [kg/m^3]
    contact: str        # "edl" or "hertz" -- decides whether K_p carries a size
    lam_d: float = LAMBDA_D

    @property
    def rho(self) -> float:
        return self.rho_p * self.phi + 1000.0 * (1 - self.phi)

    @property
    def dphi(self) -> float:
        return self.phi - self.phi_c

    def sys(self, a: float | None = None) -> DilatantSystem:
        return DilatantSystem(name=self.name, phi=self.phi, phi_c=self.phi_c,
                              a=self.a if a is None else a, rho=self.rho)

    def kappa(self, a: float | None = None) -> float:
        return self.sys(a).kappa

    def k_p(self, a: float | None = None) -> float:
        """Drained network modulus [Pa], with P_conf set to the impact stress.

        ⚠ The contact law decides whether this carries a particle size at all:
        double-layer gives K_p = P a /(6 lambda_D), Hertz gives K_p ~ E^(2/3) P^(1/3)
        with NO size in it (PORE_PHYSICS 6.3 vs 6.4). The first version of this module
        used the double-layer form for cornstarch, which is the wrong law for it.
        """
        aa = self.a if a is None else a
        p_conf = self.rho * 1.0  # placeholder scale; cancels in every ratio used here
        if self.contact == "edl":
            return p_conf * aa / (6.0 * self.lam_d)
        return float("nan")      # Hertz: no size dependence, so no size threshold


# [M26]'s cornstarch: phi_c = 0.37, D0 = 4.0 mm, rho_p = 1550, and it is HERTZIAN
# (PORE_PHYSICS 6.4). The lab silica: phi_c = 0.42 is a PLACEHOLDER, D0 = 2.5 mm,
# rho_p = 2000, double-layer contacts.
SILICA = Material("silica (lab)", 0.8e-6, 0.44, 0.42, 2.5e-3, 2000.0, "edl")
CORNSTARCH = Material("cornstarch [M26]", 10.0e-6, 0.42, 0.37, 4.0e-3, 1550.0, "hertz")
MATERIALS = (SILICA, CORNSTARCH)

U0_LAB = 3.0
U0_M26 = (0.5, 7.2)     # [M26]'s own velocity range


# =====================================================================================
# (A) the groups, evaluated at each material's OWN point
# =====================================================================================
def group_values(m: Material, u0: float) -> dict[str, float]:
    """Every size-dependent dimensionless group, at this material's own parameters."""
    gdot = u0 / m.d0
    a = m.a
    out: dict[str, float] = {}

    # G1 ID -- macroscopic drainage ([M26] Eq. 6). ID ~ a^2. Large = re-entrant liquid.
    out["ID (drainage, ell=D0)"] = ID(m.sys(), u0, m.d0)

    # G2 ID with a grain-scale drainage length: every factor of a cancels, so this
    # group is SIZE-BLIND by construction and predicts no threshold.
    out["ID (drainage, ell=10a)"] = ID(m.sys(), u0, m.d0) * (m.d0 / (10.0 * a)) ** 2

    # G3 dilatant strengthening capacity, PORE_PHYSICS Eq. 9.1. ~ a^1 for a
    # double-layer contact, SIZE-BLIND for a Hertzian one.
    out["p_dil_max/P_conf"] = ((a / (6.0 * m.lam_d)) * (m.dphi / m.phi)
                               if m.contact == "edl" else float("nan"))

    # G4 Peclet. CONVENTION: Pe = 6 pi mu R^3 gammadot / (kB T), R = RADIUS.
    out["Pe (Brownian)"] = (6.0 * math.pi * MU_W * (a / 2) ** 3 * gdot / (K_B * T_LAB))

    # G5 Stokes, in the convention used elsewhere in this repo.
    out["St (grain inertia in fluid)"] = m.rho_p * a**2 * u0 / (9.0 * MU_W * m.d0)

    # G6 capillary: impact stress against the pore-throat entry pressure.
    out["rho U0^2 / p_capillary"] = (m.rho * u0**2
                                     / (2.0 * GAMMA_LV / (R_THROAT_COEF * a)))

    # G7 inertial number at the impact stress scale.
    out["I (inertial number)"] = gdot * a * math.sqrt(m.rho_p / (m.rho * u0**2))

    # --- ADDED after the referee pass; these were missing ---------------------------
    # G8 De = tau_p/tau_imp, the poroelastic drainage time over the contact time. This
    # is the group PORE_PHYSICS section 8 says the document is built on, and it answers
    # the same drainage question as G1 with a different closure.
    kp = m.k_p()
    out["De (tau_p/tau_imp)"] = (MU_W * m.d0 * u0 / (m.kappa() * kp)
                                 if math.isfinite(kp) else float("nan"))

    # G9 tau_chi/tau_imp -- Jerome's own validity condition, storage in the FLUID's
    # compressibility rather than the network's. inertial_dilatant.poroelastic_time_ratio
    # is the repo's implementation of it.
    out["tau_chi/tau_imp"] = (MU_W * CHI_W * m.d0 * u0 / m.kappa())

    # G10 gravitational Peclet: does the pack sediment while the drop hangs?
    drho = m.rho_p - 1000.0
    out["Pe_grav"] = (drho * G_ACCEL * (math.pi / 6.0) * a**3 * (a / 2)
                      / (K_B * T_LAB))

    # G11 particle Reynolds number.
    out["Re_p"] = m.rho * gdot * a**2 / MU_W

    # G12 double layer against the pore throat (PORE_PHYSICS 5.2): above ~1 the
    # Carman-Kozeny permeability is not the bare one.
    out["lambda_D / r_throat"] = m.lam_d / (R_THROAT_COEF * a)
    return out


# =====================================================================================
# (B) a critical size, WITHIN one material's ladder -- the only well-posed version
# =====================================================================================
def a_star_within(m: Material, target_id: float, u0: float) -> float:
    """The diameter at which THIS material, at ITS OWN phi/phi_c/D0/U0, reaches
    `target_id` on the macroscopic-drainage group. ID ~ a^2, so it is a pure rescaling
    of the material's own current ID -- no other material's parameters enter."""
    return m.a * math.sqrt(target_id / ID(m.sys(), u0, m.d0))


def a_star_capacity(m: Material) -> float:
    """The diameter at which p_dil_max/P_conf = 1 for THIS material.

    Returns NaN for a Hertzian contact, where K_p carries no size and the criterion has
    no root -- which is the honest answer for cornstarch and was wrong in v1.
    """
    if m.contact != "edl":
        return float("nan")
    return 6.0 * m.lam_d * m.phi / m.dphi


# how each threshold moves. ⚠ the T1 rows were WRONG in v1 -- see the notes.
SCALINGS = {
    "G1 ID (ell = D0)": dict(
        U0="U0^-1/2", salt="⚠ NOT blind -- see note", dphi="dphi^+1/2", D0="D0^+1/2",
        CV="d32 AND phi_c -- see note"),
    "G2 ID (ell ~ a)": dict(U0="—", salt="—", dphi="—", D0="—",
                            CV="—  (size-blind: no threshold)"),
    "G3 capacity": dict(U0="NONE", salt="~ lambda_D", dphi="dphi^-1", D0="NONE",
                        CV="through dphi(CV)"),
    "G4 Pe": dict(U0="U0^-1/3", salt="—", dphi="—", D0="D0^+1/3", CV="—"),
    "G5 St": dict(U0="U0^-1/2", salt="—", dphi="—", D0="D0^+1/2", CV="—"),
    "G6 capillary": dict(U0="U0^-2", salt="—", dphi="—", D0="NONE", CV="—"),
    "G8 De": dict(U0="U0^-1/3", salt="~ lambda_D^1/3", dphi="—", D0="D0^-1/3",
                  CV="through d32"),
    "G9 tau_chi/tau_imp": dict(U0="U0^-1/2", salt="—", dphi="—", D0="D0^-1/2",
                               CV="through d32"),
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", action="store_true")
    args = ap.parse_args()
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    p("=== (A) EVERY SIZE-DEPENDENT GROUP, AT EACH MATERIAL'S OWN OPERATING POINT ===")
    p("  Not converted to a critical size: converting requires all the OTHER parameters")
    p("  to be shared, and they are not. This is the defect that made v1's headline")
    p("  wrong, and it is avoided by comparing the numbers directly.")
    p(f"  silica     : a = {SILICA.a * 1e6:.2f} um, phi = {SILICA.phi}, "
      f"phi_c = {SILICA.phi_c} (PLACEHOLDER), D0 = {SILICA.d0 * 1e3:.1f} mm, "
      f"U0 = {U0_LAB}")
    p(f"  cornstarch : a = {CORNSTARCH.a * 1e6:.2f} um, phi = {CORNSTARCH.phi}, "
      f"phi_c = {CORNSTARCH.phi_c}, D0 = {CORNSTARCH.d0 * 1e3:.1f} mm, "
      f"U0 = {U0_M26[0]} and {U0_M26[1]}")
    p()
    gs = group_values(SILICA, U0_LAB)
    gc_lo = group_values(CORNSTARCH, U0_M26[0])
    gc_hi = group_values(CORNSTARCH, U0_M26[1])
    p(f"  {'group':30s} {'silica @3':>12s} {'starch @0.5':>12s} "
      f"{'starch @7.2':>12s} {'separated?':>11s}")
    for k in gs:
        a_, b_, c_ = gs[k], gc_lo[k], gc_hi[k]
        if not (math.isfinite(a_) and math.isfinite(b_)):
            sep = "n/a"
        else:
            sep = "yes" if (a_ < min(b_, c_) or a_ > max(b_, c_)) else "no"
        p(f"  {k:30s} {a_:12.4g} {b_:12.4g} {c_:12.4g} {sep:>11s}")

    p()
    p("=== THE BOUNDARY BRACKET, done WITHOUT circularity ===")
    p("  v1 claimed the hypothesis 'independently requires' an ID boundary near 6.6e-3.")
    p("  🔴 WITHDRAWN: that bracket was algebraically ID(silica) to ID(cornstarch) at")
    p("  the lab's own parameters, i.e. an identity, and 6.6e-3 is itself just")
    p("  ID(cornstarch) evaluated at the top of [M26]'s velocity range:")
    p(f"     ID(cornstarch, M26 pars, U0 = 7.2) = {gc_hi['ID (drainage, ell=D0)']:.4e}")
    p("  so 'the bracket contains 6.6e-3' restates where 6.6e-3 came from.")
    p()
    p("  What IS non-circular: [M26] observe the transition somewhere inside their own")
    p("  velocity range, so the boundary lies between the ID their cornstarch has at")
    p("  the bottom and at the top of it:")
    p(f"     ID_boundary in [{gc_lo['ID (drainage, ell=D0)']:.3e}, "
      f"{gc_hi['ID (drainage, ell=D0)']:.3e}]")
    p(f"  and the lab's silica sits at {gs['ID (drainage, ell=D0)']:.3e}, which is")
    p(f"  {gc_lo['ID (drainage, ell=D0)'] / gs['ID (drainage, ell=D0)']:.1f}x below "
      f"even the BOTTOM of that interval.")
    p("  ⇒ On the drainage reading, the lab's silica cannot reach the inertial regime")
    p("    by velocity alone at the lab's drop size. That statement needs no a*.")

    p()
    p("=== (B) A CRITICAL SIZE, WITHIN ONE MATERIAL'S LADDER (the proposed experiment) ===")
    p("  Here every other parameter IS held fixed, so a* is well posed. This is the only")
    p("  place a critical size is quoted.")
    target = gc_hi["ID (drainage, ell=D0)"]
    for u in (0.5, 1.0, 3.0, 5.0):
        a_s = a_star_within(SILICA, target, u)
        p(f"  silica ladder, U0 = {u:4.1f} m/s: to reach ID = {target:.2e} needs "
          f"a = {a_s * 1e6:8.3f} um  ({a_s / SILICA.a:6.2f}x its present size)")
    p(f"  capacity edge for silica (p_dil_max/P_conf = 1): "
      f"a = {a_star_capacity(SILICA) * 1e6:.3f} um at 0 mM")
    p(f"  capacity edge for cornstarch: "
      f"{'n/a -- HERTZIAN contact, K_p carries no size' if not math.isfinite(a_star_capacity(CORNSTARCH)) else ''}")
    p("  ⚠ v1 built a size WINDOW from these two edges and quoted it as [3.96, 7.66] um,")
    p("  a factor 1.94. 🔴 WITHDRAWN as a quantitative claim: PORE_PHYSICS 6.3 says the")
    p("  affine K_p is an UPPER bound and real packings measure 2-5x softer, which moves")
    p("  the lower edge UP by the same factor and CLOSES the window at 2x. A 1.94x band")
    p("  cannot be asserted against a >=2x one-sided uncertainty on one of its edges.")
    p("  What survives is the STRUCTURE -- a capacity condition below and a rate")
    p("  condition above -- not the numbers, and not the sorting of any material by it.")

    if args.sweep:
        p()
        p("=== HOW EACH THRESHOLD MOVES ===")
        p(f"  {'group':22s} {'U0':>10s} {'salt':>24s} {'dphi':>12s} {'D0':>10s}")
        for k, s in SCALINGS.items():
            p(f"  {k:22s} {s['U0']:>10s} {s['salt']:>24s} {s['dphi']:>12s} "
              f"{s['D0']:>10s}")
        p()
        p("  🔴 NOTE, and it breaks v1's headline discriminator. v1 said G3 is the ONLY")
        p("  candidate that moves with salt, and built the experiment on that. It is")
        p("  false: PORE_PHYSICS 5.2 says the double layer occupies the pore throat at")
        p("  low salt, and at 800 nm")
        for lam, salt in LAMBDA_D_LADDER:
            m2 = Material(SILICA.name, SILICA.a, SILICA.phi, SILICA.phi_c, SILICA.d0,
                          SILICA.rho_p, SILICA.contact, lam)
            p(f"     lambda_D/r_throat = {m2.lam_d / (R_THROAT_COEF * m2.a):.3f} "
              f"at {salt}")
        p("  so the permeability that G1 and G8 are built on is itself salt-dependent")
        p("  at the small end, by an order-unity factor. A salt ladder therefore moves")
        p("  the drainage groups too, and -- separately -- moves sigma*, i.e. the")
        p("  frictional mechanism as well. The salt axis is NOT a clean discriminator.")
        p()
        p("=== WHAT CV DOES, and the two routes are not comparable ===")
        p("  Route 1, through d32: CV moves the SAMPLE along the size axis, because")
        p("  Carman-Kozeny's length is the Sauter mean and d32/d10 = (1+CV^2)^2 for a")
        p("  lognormal at fixed NUMBER mean.")
        p("  ⚠ PROVENANCE UNSTATED: whether the lab's '800 nm' is a number mean at all")
        p("  is not on record. psd.py notes laser diffraction reports d43; if 800 nm is")
        p("  a volume or intensity mean then d32 FALLS as CV rises and this route")
        p("  reverses sign. Q1 on the Status card.")
        p()
        p("  Route 2, through phi_c: CV raises phi_c, which shrinks dphi, which moves")
        p("  a* DOWN (a* ~ sqrt(dphi)). v1's SCALINGS said G1 depends on CV 'through")
        p("  d32 only', which contradicted this project's own measured +0.062 shift.")
        p()
        p(f"  {'phi_c':>7s} {'dphi':>8s} {'a*(silica ladder)':>18s} "
          f"{'CV needed via d32':>18s}")
        for phi_c in (0.420, 0.425, 0.430, 0.435):
            m2 = Material(SILICA.name, SILICA.a, SILICA.phi, phi_c, SILICA.d0,
                          SILICA.rho_p, SILICA.contact)
            a_s = a_star_within(m2, target, U0_LAB) * 1e6
            need = a_s / (SILICA.a * 1e6)
            cv_need = math.sqrt(math.sqrt(need) - 1.0) if need >= 1 else float("nan")
            p(f"  {phi_c:7.3f} {SILICA.phi - phi_c:8.4f} {a_s:18.3f} "
              f"{cv_need:18.2f}")
        p()
        p("  🔑 BOTH CV routes push the same way and they are NOT independent, so the")
        p("  'CV = 1.45 needed' figure of v1 is an UPPER bound that assumed phi_c fixed.")
        p("  A phi_c shift of +0.015 -- a quarter of what this project measured at")
        p("  CV = 0.5 -- already takes it to ~1.1. The honest statement is that the")
        p("  required CV is between about 1.0 and 1.5 and depends on a phi_c(CV) slope")
        p("  that phase8_criticalstate.py is measuring right now.")

    out = os.path.join(HERE, "out", "size_threshold.txt")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
