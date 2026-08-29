"""
inertial_dilatant.py — the inertial-dilatant number ID, applied to colloidal silica.

WHAT ID IS, AND WHERE IT COMES FROM
-----------------------------------
The number the PI asked about is defined, named, and used in:

  [M26] A. Mobaseri, L. Gordillo, C. Burton, S. Yoon, D. Lee, S. Kumar,
        M. M. Driscoll & X. Cheng, "Inertia-Dilatancy Interplay Governs
        Shear-Thickening Drop Impact", Phys. Rev. Lett. 136, 148201 (2026),
        doi:10.1103/fyx7-jb1d; preprint arXiv:2601.12642 (19 Jan 2026).
        Equations below are from the arXiv v1; the Supplemental Material is
        NOT on arXiv and may carry prefactors this module cannot see.

They "coin" it there (their word), for cornstarch-water drops on aluminium,
D0 = 4.0 mm, U0 = 0.5-7.2 m/s, phi = 0.30-0.43. Note the author list: this is
Driscoll's own paper with Cheng's group, so "the way Xiang Cheng and company do
it" means this construction specifically.

THE CONSTRUCTION, IN THEIR NOTATION
-----------------------------------
Two competing forces on the substrate.

INERTIAL (Wagner/Gordillo-Sun-Cheng, their Eq. 1):
    F_i(t) = c A0 sqrt(t) exp(-t/tau),    A0 = rho U0^(5/2) D0^(3/2)
    measured c = 3.99 +/- 0.75  vs  theory 3 sqrt6 / 2 = 3.67
    tau = a t_max^F with a = 1.88 +/- 0.05 (theory 2)
  scaling:  F_i ~ rho U0^2 D0^2

DILATANT (Reynolds-Darcy, after Jerome, Vandenberghe & Forterre PRL 2016):
    grad^2 P_f = (mu_w/kappa) gammadot tan(Psi),   tan(Psi) = alpha * dphi
    P_f      ~ (mu_w/kappa) dphi r_c |dh/dt|                      their Eq. 3
    F(t)     ~ A pi r_c^3 (mu_w/kappa) dphi |dh/dt|               their Eq. 4
    A        = [(1-sin B)/(1+sin B)] exp(pi tan B)                their Eq. 5
               with tan B = tan B0 + tan Psi
    kappa    = (1-phi)^3 a^2 / (180 phi^2)     [Carman-Kozeny; a = DIAMETER]
    dphi     = phi - phi_c,   alpha ~ 1,   phi_c ~ 0.37 for cornstarch
    B0       = 44 +/- 2 deg (fitted; = cornstarch's repose angle)
  scaling (setting r_c ~ h ~ D0 and |dh/dt| ~ U0):
    F_d ~ A mu_w U0 D0^3 dphi / kappa

THE NUMBER (their Eq. 6):

    ID = F_i / F_d = rho kappa U0 / (A mu_w dphi D0) = t_d / t_i             (6)
    t_i ~ D0/U0                (inertial time)
    t_d ~ rho kappa/(A mu_w dphi)   (dilatant time)

    ID << 1 : dilatancy wins early -> SOLID-LIKE impact          (their Regime 2)
    ID >> 1 : inertia wins early -> RE-ENTRANT LIQUID-LIKE       (their Regime 3)
    ID ~ 1  : the two forces trade off -> extended force plateau, whose
              measured duration peaks at ID ~ 1 (their Fig. S8)

WHAT MATTERS FOR SILICA, AND IT IS ONE THING
--------------------------------------------
    ID  proportional to  kappa  proportional to  a^2

Cornstarch a ~ 10 um; the lab's silica a = 0.8 um. That is a factor
(0.8/10)^2 = 1/156 in kappa, hence in ID, before anything else is considered.
Colloidal silica therefore sits VERY deep on the dilatant side, and the
re-entrant liquid regime that motivates the whole paper is essentially
unreachable for it at bench impact speeds.

THE CAVEAT THAT HAS TO TRAVEL WITH EVERY NUMBER BELOW
-----------------------------------------------------
[M26] write F_d with "~", after setting r_c ~ h ~ D0 and |dh/dt| ~ U0. Those
are order-unity substitutions for quantities that in reality vary by factors of
several during one impact (r_c reaches ~2 D0 while h falls to ~0.2 D0, and
F_d goes as r_c^3). So the ABSOLUTE value of ID carries an undetermined
prefactor, plausibly O(10)-O(100). Evaluating our implementation at [M26]'s own
cornstarch conditions returns ID ~ 1e-2 at the top of their velocity range,
not ID ~ 1, which is the size of that prefactor freedom.

CONSEQUENCE, AND IT IS THE MAIN METHODOLOGICAL POINT OF THIS MODULE:
    quote RATIOS of ID between systems, which are prefactor-free.
    do NOT quote ID against 1 as if the threshold were calibrated.
The ratio silica:cornstarch is the answer the question actually needs.

A SHARPER VERSION OF THE SAME OBJECTION, from the referee pass. ID depends on
U0 and D0 only through U0/D0, which is exactly the abscissa of [M26]'s Fig. 1(d).
So the ID = 1 contour is a determinate curve in the plotted plane, and Eq. 6 as
printed puts it at U0/D0 ~ 1.9e5 1/s while the experiment spans 125-1800 1/s --
two decades off the right edge of the figure. A curve is nevertheless drawn
inside the figure, and its SHAPE is right (monotone in phi, 1.9x over
phi = 0.41-0.43). So the discrepancy is a constant prefactor of ~150, not a
different functional form, and it must live somewhere this module cannot see.

WHAT IS AND IS NOT A RESOLUTION.
  * `ID_capped` was proposed here and is now LARGELY WITHDRAWN -- read its
    docstring before using it. The premise (Eq. 3 returns unreachable pressures
    for submicron particles) survives; the fix does not.
  * `poroelastic_time_ratio` is the size-dependent statement that DOES survive,
    and it is [M26]'s own source's criterion rather than an invention.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MU_W = 1.0e-3      # water viscosity [Pa.s], as fixed in [M26]
ALPHA_DIL = 1.0    # tan(Psi) = alpha * dphi, alpha ~ 1 [M26]
BETA0_CORNSTARCH_DEG = 44.0   # fitted friction angle [M26]; = repose angle
P_ATM = 1.013e5    # [Pa] — the suction a pore fluid can sustain before cavitating


@dataclass(frozen=True)
class DilatantSystem:
    """A suspension described the way [M26] describes one."""

    name: str
    phi: float
    phi_c: float          # dilatancy-onset (critical) volume fraction
    a: float              # particle DIAMETER [m]  (note: diameter, not radius)
    rho: float            # suspension density [kg/m^3]
    gamma: float = 0.072  # suspension-air surface tension [N/m]
    beta0_deg: float = BETA0_CORNSTARCH_DEG
    mu_w: float = MU_W

    @property
    def dphi(self) -> float:
        """phi - phi_c. NEGATIVE means the packing is below the dilatancy onset,
        in which case the Reynolds-Darcy branch does not exist at all."""
        return self.phi - self.phi_c

    @property
    def kappa(self) -> float:
        """Carman-Kozeny permeability (1-phi)^3 a^2 / (180 phi^2) [m^2]."""
        return (1.0 - self.phi) ** 3 * self.a**2 / (180.0 * self.phi**2)

    @property
    def A(self) -> float:
        """Effective friction coefficient, [M26] Eq. 5."""
        tan_psi = ALPHA_DIL * max(self.dphi, 0.0)
        tan_b = np.tan(np.radians(self.beta0_deg)) + tan_psi
        b = np.arctan(tan_b)
        return float((1 - np.sin(b)) / (1 + np.sin(b)) * np.exp(np.pi * tan_b))

    @property
    def t_d(self) -> float:
        """Dilatant time rho kappa / (A mu_w dphi) [s], [M26]."""
        if self.dphi <= 0:
            return float("inf")
        return self.rho * self.kappa / (self.A * self.mu_w * self.dphi)

    @property
    def r_throat(self) -> float:
        """Radius of the pore throat in a dense sphere pack [m].

        CORRECTED after review. The air-invasion (drainage) criterion is set by
        the largest sphere that fits through a triangular throat between three
        touching particles, r_throat = 0.155 R_particle, NOT by the particle
        size itself. With a = DIAMETER, R_particle = a/2, so

            r_throat = 0.0775 a

        Using a instead of r_throat understates the capillary ceiling by 12.9x,
        which is enough to change every conclusion that depends on it.
        """
        return 0.155 * (self.a / 2.0)

    @property
    def p_capillary(self) -> float:
        """Capillary entry pressure of the pack, 2 gamma / r_throat [Pa].

        Brown & Jaeger (Rep. Prog. Phys. 77, 046602 (2014)) argue that the
        largest impact stresses in a dense suspension come from "frustrated
        dilatancy against a stiff boundary", and that the stiff boundary is
        often the LIQUID-AIR INTERFACE: the interface holds the particles in
        until the suction can drive air through the pore throats, after which
        the drop granulates.

        For both systems here this comes out well above one atmosphere
        (186 kPa for 10 um cornstarch, 2.3 MPa for 800 nm silica), so the
        capillary branch NEVER binds and the ceiling is always cavitation.
        That is a real result and it is what removes the particle-size
        dependence from `ID_capped` - see the caveats there.
        """
        return 2.0 * self.gamma / self.r_throat

    @property
    def p_dilatant_ceiling(self) -> float:
        """The largest dilatant pore suction the system can actually sustain [Pa].

        Two independent ceilings, whichever binds first:
          * cavitation of the pore fluid at ~1 atm of suction
          * the interface being pulled into a pore, 2 gamma / a
        Below either, Darcy sets the pressure; above, the Reynolds-Darcy
        formula is describing a state the material cannot reach.
        """
        return min(P_ATM, self.p_capillary)


# --- the two systems -------------------------------------------------------
def cornstarch_M26(phi: float = 0.42) -> DilatantSystem:
    """[M26]'s own system: cornstarch-water, a ~ 10 um, phi_c ~ 0.37."""
    return DilatantSystem(name=f"cornstarch-water (M26), phi={phi:g}", phi=phi,
                          phi_c=0.37, a=10.0e-6,
                          rho=1550.0 * phi + 1000.0 * (1 - phi))


def lab_silica_dilatant(phi: float = 0.44, phi_c: float = 0.42,
                        beta0_deg: float = BETA0_CORNSTARCH_DEG) -> DilatantSystem:
    """The Driscoll-lab 800 nm colloidal silica, in [M26]'s variables.

    phi_c is the one genuinely unmeasured input. What is known:
      * Shah, Arora & Driscoll (Commun. Phys. 5, 222 (2022)) see the first
        localised solid patches at phi = 0.47 and bulk solidification at
        phi >= 0.49, for the same 830 nm Stoeber silica.
      * the repo's Wyart-Cates fit gives a frictional jamming fraction
        phi_M = 0.5176.
      * this lab's own drop-impact sessions see solidification fronts at
        phi = 0.42-0.45.
    Those three cannot all be right on one phi axis — which is the documented
    phi-scale offset (Paper 1 card, section 8). The default here, phi_c = 0.42,
    is the value implied by the lab's OWN observation that fronts appear at
    phi = 0.42; the sensitivity to it is reported explicitly, and it is large.

    beta0 (the friction angle) is unmeasured for silica; 44 deg is carried over
    from cornstarch. A depends on it exponentially, so it is the second-largest
    uncertainty after phi_c.
    """
    return DilatantSystem(name=f"lab silica 800 nm, phi={phi:g}", phi=phi,
                          phi_c=phi_c, a=800.0e-9,
                          rho=2000.0 * phi + 1000.0 * (1 - phi),
                          beta0_deg=beta0_deg)


# --- the number ------------------------------------------------------------
def ID(sys: DilatantSystem, U0: float, D0: float) -> float:
    """Inertial-dilatant number, [M26] Eq. 6, verbatim.

        ID = rho kappa U0 / (A mu_w dphi D0)

    Returns inf if dphi <= 0 (no dilatant branch: the packing is below the
    dilatancy onset, so there is no dilatant force for inertia to compete with).
    """
    if sys.dphi <= 0:
        return float("inf")
    return sys.rho * sys.kappa * U0 / (sys.A * sys.mu_w * sys.dphi * D0)


def ID_from_times(sys: DilatantSystem, U0: float, D0: float) -> float:
    """Same number as t_d/t_i — an independent route, asserted equal in tests."""
    return sys.t_d / (D0 / U0)


def dilatant_pressure(sys: DilatantSystem, U0: float, D0: float,
                      r_c: float | None = None, dhdt: float | None = None) -> float:
    """Reynolds-Darcy pore suction, [M26] Eq. 3, with their scaling substitutions.

        P_f ~ (mu_w/kappa) dphi r_c |dh/dt| ,   r_c ~ D0,  |dh/dt| ~ U0

    Returned as a POSITIVE magnitude of a suction. Compare it against
    `p_dilatant_ceiling` before believing it.
    """
    if sys.dphi <= 0:
        return 0.0
    r_c = D0 if r_c is None else r_c
    dhdt = U0 if dhdt is None else dhdt
    return (sys.mu_w / sys.kappa) * sys.dphi * r_c * dhdt


def ID_capped(sys: DilatantSystem, U0: float, D0: float) -> float:
    """ID with the dilatant pressure capped at the cavitation ceiling.

    *** LARGELY WITHDRAWN AFTER REVIEW. KEPT FOR THE RECORD, NOT FOR USE. ***

    THE PREMISE STILL STANDS: Eq. (3) has P_f ~ 1/kappa ~ 1/a^2 and returns
    pore suctions of 1e7-1e8 Pa for 800 nm particles, which is not a state a
    pore fluid can reach. Something must limit it.

    THE FIX DOES NOT. Four objections, all of which land:

    1. WRONG CAPILLARY LENGTH, and it was load-bearing. Air invasion is set by
       the PORE THROAT, r_throat = 0.155 R_particle, not by the particle size.
       Corrected (see `r_throat`), the capillary ceiling is 186 kPa for
       cornstarch and 2.3 MPa for silica, so it NEVER binds and the ceiling is
       always ~1 atm for both systems. Then kappa, a, dphi and D0 all cancel
       out of ID_cap, and the number can no longer tell 10 um cornstarch from
       800 nm silica - which was its entire purpose.

    2. THE SUPPORTING EVIDENCE WAS MISREAD. Jerome, Vandenberghe & Forterre
       (PRL 117, 098003 (2016)) measured -10 to -100 kPa in 170 um GLASS BEADS,
       where uncapped Darcy also predicts ~1e4-1e5 Pa. They demonstrate no
       saturation. Their own closing paragraph extrapolates to cornstarch and
       predicts "about 1-10 MPa ... in good agreement with impact stress
       estimations in cornstarch". The cited authority asserts the opposite of
       what was attributed to it.

    3. ONE ATMOSPHERE IS NOT A MATERIAL CONSTANT. Degassed water sustains -25
       to -140 MPa. The ~1 atm figure is the heterogeneous-nucleation threshold
       of gas-saturated water, i.e. a sample-preparation property. That makes it
       a falsifiable hypothesis (degas and repeat), not a ceiling.

    4. CAPPING A PRESSURE DOES NOT CAP A FORCE. r_c keeps growing as sqrt(t),
       so the cap only turns r_c^3 into r_c^2, and once the fluid does cavitate
       the pack granulates and carries a capillary-cohesive stress instead. The
       force changes microscopic origin, not magnitude.

    Note also that ID_cap = rho U0^2 / (A P_ceiling) is algebraically the
    rho U0^2 / (a stress) form that `dilatancy.py` argues against. The two
    proposals are the same construction with different denominators.

    USE `poroelastic_time_ratio` INSTEAD. It is the size-dependent statement
    that survives, it is the criterion Jerome et al. themselves state, and it
    reaches a defensible version of the same conclusion.
    """
    if sys.dphi <= 0:
        return float("inf")
    P_uncapped = dilatant_pressure(sys, U0, D0)
    P_eff = min(P_uncapped, sys.p_dilatant_ceiling)
    return sys.rho * U0**2 / (sys.A * P_eff)


def reentrant_velocity(sys: DilatantSystem, D0: float, capped: bool = True) -> float:
    """Impact speed at which ID reaches 1 — the re-entrant liquid threshold [m/s].

    Uncapped:  ID = 1  =>  U0 = A mu_w dphi D0 / (rho kappa)
    Capped:    ID = 1  =>  U0 = sqrt(A P_ceiling / rho)
    """
    if sys.dphi <= 0:
        return float("nan")
    if not capped:
        return sys.A * sys.mu_w * sys.dphi * D0 / (sys.rho * sys.kappa)
    return float(np.sqrt(sys.A * sys.p_dilatant_ceiling / sys.rho))


def compare(sys_a: DilatantSystem, sys_b: DilatantSystem, U0: float, D0: float) -> dict:
    """Ratio of ID between two systems — the prefactor-free comparison."""
    ida, idb = ID(sys_a, U0, D0), ID(sys_b, U0, D0)
    return {
        "system_a": sys_a.name, "system_b": sys_b.name,
        "ID_a": ida, "ID_b": idb, "ratio_a_over_b": ida / idb,
        "kappa_ratio": sys_a.kappa / sys_b.kappa,
        "kappa_a_m2": sys_a.kappa, "kappa_b_m2": sys_b.kappa,
    }


def report(sys: DilatantSystem, U0: float, D0: float) -> dict:
    """Everything about one system at one impact condition."""
    P_un = dilatant_pressure(sys, U0, D0)
    return {
        "name": sys.name,
        "phi": sys.phi, "phi_c": sys.phi_c, "dphi": sys.dphi,
        "a_m": sys.a, "rho": sys.rho,
        "kappa_m2": sys.kappa, "A": sys.A,
        "U0": U0, "D0": D0,
        "t_i_s": D0 / U0, "t_d_s": sys.t_d,
        "ID": ID(sys, U0, D0),
        "P_darcy_Pa": P_un,
        "P_capillary_Pa": sys.p_capillary,
        "P_ceiling_Pa": sys.p_dilatant_ceiling,
        "darcy_overshoot": P_un / sys.p_dilatant_ceiling if sys.dphi > 0 else 0.0,
        "ID_capped": ID_capped(sys, U0, D0),
        "U_reentrant_capped_m_s": reentrant_velocity(sys, D0, capped=True),
        "U_reentrant_uncapped_m_s": reentrant_velocity(sys, D0, capped=False),
        "rho_U0_sq_Pa": sys.rho * U0**2,
    }


# --------------------------------------------------------------------------
# The size-dependent statement that survives review
# --------------------------------------------------------------------------
CHI_WATER = 1.0 / 2.2e9   # isothermal compressibility of water [1/Pa]


def poroelastic_time_ratio(sys: DilatantSystem, U0: float, D0: float,
                           chi: float = CHI_WATER) -> dict:
    """tau_d / tau_i: does the pore pressure have time to equilibrate?

    THIS REPLACES `ID_capped` AS THE DEFENSIBLE SIZE-DEPENDENT RESULT, and it
    is not an invention: Jerome, Vandenberghe & Forterre (PRL 117, 098003
    (2016)) state the criterion themselves. Their Darcy-Reynolds pore-pressure
    equation assumes the pore fluid is incompressible, and they say that holds

        "as long as the diffusion timescale for gas expulsion
         tau_d ~ eta chi L^2 / kappa ... is small compared to the impact
         timescale tau_i ~ D/V"

    with chi the fluid compressibility. Evaluating it:

        tau_d = mu_w * chi * D0^2 / kappa ,      tau_i = D0 / U0
        tau_d / tau_i = mu_w * chi * D0 * U0 / kappa                      (19)

    and since kappa ~ a^2, this ratio goes as a^(-2): it is the SAME particle
    size sensitivity that makes ID small for colloids, but expressed as a
    validity condition on the model rather than as a fudge to its output.

    WHAT IT SAYS FOR THE TWO SYSTEMS (pure, fully-degassed water):
        cornstarch, a = 10 um, D0 = 4 mm, U0 = 5.1 m/s :  tau_d/tau_i = 0.014
        silica,     a = 800 nm, D0 = 2.5 mm, U0 = 2.3 m/s:  tau_d/tau_i = 0.81

    So for cornstarch the incompressible-Darcy assumption holds with a factor
    70 to spare, and the ID paper is on safe ground. For 800 nm silica it is
    marginal at best: the pore fluid's own compliance matters over the impact
    time, so Eq. 3 -- and therefore ID itself -- is being used at the edge of
    its derivation. Entrained gas makes this dramatically worse: 0.1% gas by
    volume raises chi by ~20x and pushes silica to tau_d/tau_i ~ 18.

    THE HONEST CONCLUSION about applying ID to submicron colloids: not "the
    Darcy pressure saturates at 1 atm" (that was wrong, see `ID_capped`), but
    "the Darcy limit itself stops applying". The saturated pore pressure
    P_f^sat = Delta_phi/(phi chi) is still enormous, so compressibility alone
    does not deliver a low ceiling -- only a phase change would, and that is a
    sample-preparation property, testable by degassing.
    """
    tau_i = D0 / U0
    tau_d = sys.mu_w * chi * D0**2 / sys.kappa
    return {
        "tau_i_s": tau_i,
        "tau_d_s": tau_d,
        "ratio": tau_d / tau_i,
        "darcy_incompressible_valid": bool(tau_d / tau_i < 0.1),
        "P_f_saturated_Pa": (sys.dphi / (sys.phi * chi)) if sys.dphi > 0 else 0.0,
        "chi_used_per_Pa": chi,
        "note": ("Jerome, Vandenberghe & Forterre PRL 117, 098003 (2016): the "
                 "Reynolds-Darcy pore-pressure equation assumes an incompressible "
                 "pore fluid, valid only for tau_d << tau_i."),
    }
