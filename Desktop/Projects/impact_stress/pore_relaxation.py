"""
pore_relaxation.py — how fast can pore pressure relax, and how big can the
Reynolds-dilatant pressure get, for 800 nm silica versus 10 um cornstarch.

WHY THIS MODULE EXISTS (read before using any number out of it)
--------------------------------------------------------------
`inertial_dilatant.py` already carries the Reynolds-Darcy machinery: the
Carman-Kozeny permeability, the inertial-dilatant number ID, and the suction
P_f ~ (mu/kappa) dphi r_c |dh/dt|. It also carries ONE pore-relaxation estimate,
`poroelastic_time_ratio`, and that estimate answers a narrower question than it
looks like it answers.

It puts the whole storage into the PORE FLUID'S OWN COMPRESSIBILITY,
chi_w = 4.55e-10 /Pa. That is the rigid-skeleton limit: grains rigid, network
rigid, and the only thing that can hold a pore-pressure increment is the water
squeezing. It is the right criterion for the question Jerome, Vandenberghe &
Forterre actually ask -- "may I treat the pore fluid as incompressible in the
Darcy-Reynolds equation?" -- and it gives tau_d/tau_i = 0.014 (cornstarch) and
0.81 (silica).

It is NOT the pore relaxation time. A saturated packing relaxes an imposed pore
pressure by draining, and the drainage is resisted by the PARTICLE NETWORK's
stiffness, not by water's bulk modulus. The two storage terms are in series and
differ by six to ten decades, so which one you use changes the answer by six to
ten decades. This module carries both, in one storage coefficient, and reduces
to the existing result in the appropriate limit (asserted in the tests).

THE PHYSICS, IN FIVE LINES
--------------------------
Mixture continuity + Darcy for a saturated deformable packing gives a diffusion
equation for the pore pressure (Biot 1941; Terzaghi consolidation):

    S dp/dt = div[ (kappa/mu_f) grad p ] + (dilation source)

    S      = (1 - phi) chi_f  +  1/K_p            storage coefficient   [1/Pa]
    D_p    = kappa / (mu_f S)                     pore-pressure diffusivity [m^2/s]
    tau_p  = L^2 / D_p                            relaxation over a length L [s]

with K_p the DRAINED bulk modulus of the particle network. The two limits:

    K_p -> inf   : S -> (1-phi) chi_f     the existing rigid-skeleton estimate
    K_p finite   : S -> 1/K_p             D_p = kappa K_p / mu_f, which is also
                                          exactly the colloid's COLLECTIVE
                                          DIFFUSION coefficient

For every case in this module the second term dominates by many decades, i.e.
the network is far more compliant than water. That is the whole point.

THE DILATANCY APPROXIMATION
---------------------------
Reynolds dilatancy (Pailha & Pouliquen JFM 633, 115 (2009); Jerome, Vandenberghe
& Forterre PRL 117, 098003 (2016)): a packing denser than phi_c must expand to
shear. The expansion has to be filled by solvent drawn in from a distance L_d,
Darcy resists that flow, and the resulting suction presses the particles
together. Scaling the Darcy balance (kappa/mu)(p/L_d^2) ~ dphi * gammadot,

    p_dil  =  (mu_f / kappa) * dphi * L_d^2 * gammadot                    (D1)

which is the same closure the DNS uses behind DIL_COEFF
(`wc_drop_impact.c`: p_dil = (eta_f/k) eps_dil L_d^2 gammadot) and the same as
[M26] Eq. 3 once |dh/dt| = L_d gammadot.

THE IDENTITY THAT TIES THE TWO HALVES TOGETHER. Divide (D1) by K_p and use
tau_p = mu_f L_d^2 / (kappa K_p):

    p_dil  =  K_p * dphi * (gammadot * tau_p)                             (D2)

The dilatant pressure is the network modulus times the dilation strain that
accumulates in one pore-relaxation time. Read that way it has an obvious
ceiling, because the strain cannot exceed the total dilation the packing has
available:

    p_dil_max  =  K_p * dphi / phi                                        (D3)

(with incompressible grains the particle volume phi*V is fixed, so the mixture's
volumetric strain is eps_v = -dphi/phi -- referenced to the SOLID fraction, and
consistent with the definition K_p = phi dP/dphi).

(D3) is a stiffness ceiling, and it is the honest replacement for the cavitation
cap that the 2026-08-06 brief proposed and then withdrew. Unlike the cavitation
cap it keeps a particle-size dependence, because K_p does, and it is validated
against two independent measurements in `test_pore_relaxation.py`:

  * Jerome et al. measure -10 to -100 kPa in 170 um glass beads; (D3) gives
    ~10^4-10^5 Pa there.
  * Jerome et al. extrapolate 1-10 MPa for cornstarch; (D3) gives ~10^6 Pa.

For 800 nm silica (D3) gives O(10^1-10^2) Pa -- five to six decades below the
uncapped Darcy value of 10^7-10^8 Pa, and, notably, the same order as the fitted
Wyart-Cates onset sigma* = 44 Pa.

THE LARGE UNCERTAINTY, STATED UP FRONT
--------------------------------------
K_p is not measured for either system. Everything here that depends on its
ABSOLUTE value is an order-of-magnitude statement, and the module makes K_p an
explicit argument everywhere rather than hiding a default. What is robust is
    (a) the SCALING, tau_p ~ mu L^2 / (kappa K_p) with kappa ~ a^2, and
    (b) the RATIO between the two systems, which is where the answer lives.
`network_modulus_affine` gives the micromechanical estimate and
`K_p_bracket` gives the range that should be carried through any conclusion.

CONVENTIONS
-----------
`a` is a particle DIAMETER everywhere in this file, matching Carman-Kozeny and
matching `inertial_dilatant.DilatantSystem`. (It is a RADIUS in the Wyart-Cates
literature. This has bitten the project before.)  SI units throughout.

SOURCES
-------
[BIO41] Biot, J. Appl. Phys. 12, 155 (1941).
[PP09]  Pailha & Pouliquen, J. Fluid Mech. 633, 115 (2009).
[JVF16] Jerome, Vandenberghe & Forterre, Phys. Rev. Lett. 117, 098003 (2016).
[M26]   Mobaseri et al. (incl. Driscoll & Cheng), Phys. Rev. Lett. 136,
        148201 (2026).
[ATH25] Athani, Metzger, Forterre & Mari, Phys. Rev. Fluids 10, 043301 (2025).
[MAK99] Makse, Gland, Johnson & Schwartz, Phys. Rev. Lett. 83, 5070 (1999) --
        affine bulk modulus of a sphere packing.
[CS69]  Carnahan & Starling, J. Chem. Phys. 51, 635 (1969).
[SAD22] Shah, Arora & Driscoll, Commun. Phys. 5, 222 (2022).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from inertial_dilatant import DilatantSystem, cornstarch_M26, lab_silica_dilatant

# --------------------------------------------------------------------------
# constants
# --------------------------------------------------------------------------
MU_W = 1.0e-3           # water dynamic viscosity [Pa.s] (20 C: 1.002e-3)
CHI_WATER = 1.0 / 2.2e9  # isothermal compressibility of water [1/Pa]
K_BOLTZ = 1.380649e-23  # [J/K]
T_LAB = 293.15          # [K], 20 C

# Coordination number of a dense random packing at/above isostaticity.
# Frictional isostatic value is 4; frictionless is 6; real dense packs sit at
# 5-6.  Enters the affine modulus linearly, so it is not a leading uncertainty.
Z_COORD = 6.0

# Contact-stiffness inputs, both order-of-magnitude and both flagged.
E_STARCH_PA = 1.0e8     # hydrated starch granule Young's modulus [Pa].
#   Literature for hydrated/gelatinised starch spans ~1e7-1e9 Pa; 1e8 is the
#   geometric middle. Enters K_p linearly -> a decade of spread each way.
DEBYE_DEFAULT_M = 10.0e-9
#   Debye length of the lab's silica suspension [m]. 10 nm corresponds to ~1 mM
#   1:1 electrolyte. The salt-screening series takes this from ~30 nm (0 mM,
#   set by dissolved CO2) to ~1 nm (100 mM), so K_p for silica is a function of
#   the very variable the experiment scans. See `debye_length`.
EPS_WATER = 78.4 * 8.8541878128e-12   # permittivity of water [F/m]
E_CHARGE = 1.602176634e-19            # [C]
N_AVOGADRO = 6.02214076e23            # [1/mol]


# --------------------------------------------------------------------------
# 1. network modulus K_p
# --------------------------------------------------------------------------
def osmotic_modulus_hard_sphere(phi: float, a: float, T: float = T_LAB) -> float:
    """Brownian hard-sphere osmotic bulk modulus K = phi dPi/dphi [Pa].

    Carnahan-Starling [CS69]:  Pi = n kT Z(phi),
        Z(phi) = (1 + phi + phi^2 - phi^3) / (1 - phi)^3
    with n = 6 phi / (pi a^3) the number density (a = DIAMETER).

    This is the correct K_p only BELOW the contact onset, where the particles
    are held apart by thermal/repulsive pressure alone. It is included because
    it is the one rigorous end of the bracket, and because it shows how small
    a colloid's equilibrium stiffness is: O(0.1-10) Pa at phi ~ 0.5 for 800 nm
    spheres, which is why colloidal consolidation takes days.

    Beyond phi ~ 0.494 the hard-sphere fluid branch is metastable; CS is used
    here as an extrapolated compressed-branch estimate, which is standard but
    is a real approximation. Do not read it past phi ~ 0.58 (glass).
    """
    _check_phi(phi)
    n1 = 6.0 / (np.pi * a**3)        # number density per unit volume fraction
    om = 1.0 - phi
    u = 1.0 + phi + phi**2 - phi**3
    Zc = u / om**3
    # dZ/dphi = [u'(1-phi) + 3u] / (1-phi)^4,  u' = 1 + 2phi - 3phi^2
    dZ = ((1.0 + 2.0 * phi - 3.0 * phi**2) * om + 3.0 * u) / om**4
    # Pi = n1 kT phi Z(phi)  =>  dPi/dphi = n1 kT (Z + phi Z')
    dPi_dphi = n1 * K_BOLTZ * T * (Zc + phi * dZ)
    return float(phi * dPi_dphi)


def network_modulus_affine(k_n: float, a: float, phi: float,
                           Z: float = Z_COORD) -> float:
    """Affine (Voigt) bulk modulus of a sphere packing with normal springs [Pa].

        K_p = Z phi k_n / (6 pi a)          (a = DIAMETER; [MAK99] uses radius,
                                             K = Z phi k_n / (12 pi R))

    The affine estimate is an UPPER bound: real packings relax non-affinely and
    measure 2-5x softer near isostaticity. Used here for an order of magnitude
    and a ratio, not for an absolute number.
    """
    _check_phi(phi)
    if k_n <= 0 or a <= 0:
        raise ValueError(f"k_n and a must be positive, got k_n={k_n}, a={a}")
    return float(Z * phi * k_n / (6.0 * np.pi * a))


def debye_length(c_salt_mM: float, T: float = T_LAB, valence: int = 1) -> float:
    """Debye screening length of a 1:1 electrolyte [m].

        lambda_D = sqrt( eps kT / (2 N_A e^2 z^2 c) ),  c in mol/m^3

    At 0 mM added salt the real screening is set by dissolved CO2 and the
    silica's own counterions, ~1e-5 M, i.e. ~100 nm; the function is not valid
    at c -> 0 and raises rather than returning infinity.
    """
    if c_salt_mM <= 0:
        raise ValueError("c_salt_mM must be > 0; at nominal 0 mM use the "
                         "residual ionic strength (~0.01 mM) explicitly.")
    c = c_salt_mM * 1.0e-3 * 1.0e3        # mM -> mol/m^3
    return float(np.sqrt(EPS_WATER * K_BOLTZ * T
                         / (2.0 * N_AVOGADRO * E_CHARGE**2 * valence**2 * c)))


def contact_stiffness_repulsive(sigma_star: float, a: float,
                                lam: float = DEBYE_DEFAULT_M) -> float:
    """Normal stiffness of an electrostatically stabilised colloidal contact [N/m].

        F*   = sigma_star * a^2        the Wyart-Cates onset force
        k_n ~ F* / lambda_D            the repulsive barrier has range lambda_D

    Rationale: in Wyart-Cates the onset stress is sigma* ~ F*/a^2, so the fitted
    sigma* = 44.2 Pa IS a statement about the repulsive force scale of these
    particles. Dividing by the range over which that force decays gives the
    stiffness of the "spring" the particles sit on. Both steps are scaling
    arguments, so this is an order of magnitude.

    Note the consequence: k_n ~ 1/lambda_D, so ADDING SALT STIFFENS THE NETWORK.
    That is a testable prediction of the salt-screening series and it is the
    same direction as sigma* rising with salt.
    """
    if lam <= 0:
        raise ValueError(f"lambda_D must be positive, got {lam}")
    return float(sigma_star * a**2 / lam)


def contact_stiffness_hertz_scale(E_particle: float, a: float) -> float:
    """Pressure-free scale of a solid-solid contact stiffness [N/m], k_n ~ E_p a.

    *** UPPER BOUND ONLY. DO NOT USE FOR A NUMBER. ***

    Kept because it is the estimate a scaling argument reaches for first, and
    because it is instructive how badly it fails. Checked against [JVF16]'s
    MEASURED pore suction in 170 um glass beads it overshoots by three decades:
    it corresponds to a contact squashed by delta ~ a, whereas a real pack at
    10 kPa confinement has delta/a ~ 5e-5. Use `contact_stiffness_hertz`, which
    carries the confining pressure and lands within a decade of the data.
    """
    return float(E_particle * a)


def mean_contact_force(P_conf: float, a: float, phi: float,
                       Z: float = Z_COORD) -> float:
    """Mean normal force per contact in a packing at isotropic pressure P [N].

        f = pi P a^2 / (Z phi)

    From the standard stress-from-contacts identity
    P = (Z phi / pi a^2) f  (e.g. Agnolin & Roux, PRE 76, 061302 (2007)),
    obtained by summing f*l over the Z n /2 contacts per unit volume with
    n = 6 phi/(pi a^3).
    """
    _check_phi(phi)
    if P_conf <= 0:
        raise ValueError(f"P_conf must be positive, got {P_conf}")
    return float(np.pi * P_conf * a**2 / (Z * phi))


def contact_stiffness_hertz(E_particle: float, nu: float, a: float, phi: float,
                            P_conf: float, Z: float = Z_COORD) -> float:
    """Hertzian normal stiffness at a stated confining pressure [N/m].

        f      = pi P a^2 / (Z phi)                mean contact force
        R*     = a / 4                             two spheres of radius a/2
        E*     = E / (2 (1 - nu^2))                two identical spheres
        delta  = (3 f / (4 E* sqrt(R*)))^(2/3)     Hertz overlap
        k_n    = 2 E* sqrt(R* delta)               = df/d(delta)

    k_n ~ P^(1/3), so the network modulus of a solid-grain packing is
    pressure-dependent and CANNOT be quoted without stating P. This is the
    correction that brings Eq. (D3) into agreement with [JVF16]'s measurement;
    the pressure-free k_n ~ E a is three decades too stiff at their confinement.

    Hertz assumes delta << a. The returned overlap ratio is not checked here;
    `hertz_overlap_ratio` exposes it so callers can verify.
    """
    if E_particle <= 0 or a <= 0:
        raise ValueError(f"need E>0, a>0; got E={E_particle}, a={a}")
    if not (-1.0 < nu < 0.5):
        raise ValueError(f"Poisson ratio out of range: {nu}")
    E_star = E_particle / (2.0 * (1.0 - nu**2))
    R_star = a / 4.0
    f = mean_contact_force(P_conf, a, phi, Z)
    delta = (3.0 * f / (4.0 * E_star * np.sqrt(R_star))) ** (2.0 / 3.0)
    return float(2.0 * E_star * np.sqrt(R_star * delta))


def hertz_overlap_ratio(E_particle: float, nu: float, a: float, phi: float,
                        P_conf: float, Z: float = Z_COORD) -> float:
    """delta / a for the Hertz contact above. Hertz needs this << 1."""
    E_star = E_particle / (2.0 * (1.0 - nu**2))
    R_star = a / 4.0
    f = mean_contact_force(P_conf, a, phi, Z)
    delta = (3.0 * f / (4.0 * E_star * np.sqrt(R_star))) ** (2.0 / 3.0)
    return float(delta / a)


def contact_stiffness_edl(P_conf: float, a: float, phi: float,
                          lam: float = DEBYE_DEFAULT_M,
                          Z: float = Z_COORD) -> float:
    """Normal stiffness of an electrostatically stabilised contact under load [N/m].

    An electric-double-layer repulsion decays exponentially with separation h,
    F = F* exp(-h/lambda_D), so its stiffness is set by the load it is already
    carrying rather than by a material modulus:

        k_n = -dF/dh = F / lambda_D = f / lambda_D,    f = pi P a^2/(Z phi)

    Feeding that into the affine modulus gives a strikingly simple result with
    no free parameter and no particle modulus in it at all:

        K_p = Z phi k_n / (6 pi a) = P a / (6 lambda_D)

    exposed as `network_modulus_edl`. Two consequences worth stating:
      * K_p / P = a / (6 lambda_D) ~ 13 for 800 nm silica at lambda_D = 10 nm,
        so the network is stiff RELATIVE TO ITS OWN CONFINEMENT but soft in
        absolute terms at the low pressures of a colloid.
      * K_p ~ 1/lambda_D: adding salt stiffens the network. This is a testable
        prediction of the lab's salt-screening series.
    """
    return float(mean_contact_force(P_conf, a, phi, Z) / lam)


def network_modulus_edl(P_conf: float, a: float,
                        lam: float = DEBYE_DEFAULT_M) -> float:
    """K_p = P a / (6 lambda_D) [Pa] -- the affine modulus of an EDL-stabilised
    packing, with Z and phi cancelling out exactly. See `contact_stiffness_edl`.
    """
    if lam <= 0 or a <= 0:
        raise ValueError(f"need lambda_D>0, a>0; got {lam}, {a}")
    if P_conf <= 0:
        raise ValueError(f"P_conf must be positive, got {P_conf}")
    return float(P_conf * a / (6.0 * lam))


# --------------------------------------------------------------------------
# 2. storage, diffusivity, relaxation time
# --------------------------------------------------------------------------
def storage_coefficient(phi: float, K_p: float,
                        chi_f: float = CHI_WATER) -> float:
    """Poroelastic storage S = (1-phi) chi_f + 1/K_p  [1/Pa].

    Two compliances in series-of-storage: the pore fluid squeezing inside the
    pore volume fraction (1-phi), and the skeleton compressing. K_p = inf
    recovers the rigid-skeleton limit used by
    `inertial_dilatant.poroelastic_time_ratio`.
    """
    _check_phi(phi)
    if K_p <= 0:
        raise ValueError(f"K_p must be positive, got {K_p}")
    return float((1.0 - phi) * chi_f + 1.0 / K_p)


def pore_diffusivity(kappa: float, K_p: float, phi: float,
                     mu_f: float = MU_W, chi_f: float = CHI_WATER) -> float:
    """Pore-pressure diffusivity D_p = kappa / (mu_f S)  [m^2/s].

    Identical to the colloid's collective diffusion coefficient in the
    skeleton-dominated limit, D_c = kappa K_p / mu_f.
    """
    if kappa <= 0:
        raise ValueError(f"kappa must be positive, got {kappa}")
    return float(kappa / (mu_f * storage_coefficient(phi, K_p, chi_f)))


def pore_relaxation_time(L: float, kappa: float, K_p: float, phi: float,
                         mu_f: float = MU_W, chi_f: float = CHI_WATER) -> float:
    """tau_p = L^2 / D_p  [s]: time to relax pore pressure over a length L.

    L is the DRAINAGE LENGTH -- the distance pore fluid must actually travel to
    equilibrate -- and choosing it is a modelling decision, not a material
    property. The three that matter for a drop impact are exposed as named
    lengths in `drainage_lengths`.
    """
    if L <= 0:
        raise ValueError(f"L must be positive, got {L}")
    return float(L**2 / pore_diffusivity(kappa, K_p, phi, mu_f, chi_f))


def slowest_mode_time(L: float, kappa: float, K_p: float, phi: float,
                      mu_f: float = MU_W, chi_f: float = CHI_WATER) -> float:
    """The exact e-folding time of the slowest drainage mode [s].

    tau_p = L^2/D_p is the SCALING. For the standard one-dimensional problem --
    a layer of thickness L drained at one face, sealed at the other -- the
    solution of dp/dt = D_p d2p/dx2 is a cosine series whose slowest mode decays
    as exp(-t/tau_1) with

        tau_1 = 4 L^2 / (pi^2 D_p) = 0.405 * tau_p

    so the true relaxation is about 2.5x FASTER than the scaling estimate. That
    factor is reported rather than folded in, because every conclusion in this
    module is a statement about decades and folding an O(1) factor into tau_p
    would make the scaling harder to check by hand.
    """
    return float((4.0 / np.pi**2)
                 * pore_relaxation_time(L, kappa, K_p, phi, mu_f, chi_f))


def drainage_lengths(D0: float, a: float, phi: float) -> dict[str, float]:
    """The three drainage lengths a drop impact offers, smallest first [m].

    'particle'  : a. Local pore-scale equilibration between neighbouring pores.
                  tau_p at this length is size-INDEPENDENT, because kappa ~ a^2
                  cancels L^2 = a^2. This is the fastest relaxation available
                  and it is the same for silica and cornstarch.
    'contact'   : the contact-patch radius at the force peak, 0.49 D0
                  (from the Wagner analysis in `stress.py`; the value is the
                  measured r_c at peak force, not an assumption).
                  This is the length that governs whether the SPREADING LAMELLA
                  can drain, and it is the one [M26]'s Eq. 3 uses.
    'drop'      : D0. Global drainage of the whole drop.

    The particle length is included precisely because it is the counterexample
    to "smaller particles always drain slower": at the pore scale they do not.
    The size dependence enters only through geometric lengths held fixed while
    a shrinks.
    """
    _check_phi(phi)
    return {"particle": a, "contact": 0.49 * D0, "drop": D0}


# --------------------------------------------------------------------------
# 2b. does the continuum Darcy picture survive at 800 nm?
# --------------------------------------------------------------------------
def debye_to_throat_ratio(a: float, lam: float = DEBYE_DEFAULT_M) -> float:
    """lambda_D / r_throat, with r_throat = 0.0775 a (a = DIAMETER).

    Carman-Kozeny treats the pore space as an inert plumbing network. That
    assumption needs the double layer to be thin compared with the pore throat.
    When this ratio approaches 1 the double layers of opposing pore walls
    overlap, the counterion cloud is dragged along with the flow, and the
    back-EMF it generates opposes the flow: the electroviscous effect. The
    consequence is an EFFECTIVE permeability BELOW the geometric one.

    For 800 nm silica the throat radius is 62 nm, so this ratio runs from 0.016
    at 100 mM to ~1.5 at the residual ionic strength of nominally salt-free
    water. At the low-salt end the continuum-Darcy permeability is an
    overestimate, which makes the undrained conclusion stronger, not weaker.
    No correction is applied here; the diagnostic is reported so the assumption
    is visible.
    """
    if a <= 0 or lam <= 0:
        raise ValueError(f"need a>0, lambda_D>0; got a={a}, lam={lam}")
    return float(lam / (0.0775 * a))


def peclet(gammadot: float, a: float, mu_f: float = MU_W,
           T: float = T_LAB) -> float:
    """Shear Peclet number, Pe = 6 pi mu_f (a/2)^3 gammadot / kT.

    a is a DIAMETER, so the radius a/2 appears cubed. Pe >> 1 means shear wins
    over Brownian relaxation and the packing can be treated as athermal, which
    is what every formula in this module assumes. Pe ~ 300 for the lab's silica
    at gammadot = 920 1/s: athermal, but only by two decades, and Brownian
    stress is not negligible in the osmotic modulus (see
    `osmotic_modulus_hard_sphere`). For 10 um cornstarch Pe is ~1e6 and the
    question does not arise.
    """
    if gammadot < 0:
        raise ValueError(f"gammadot must be >= 0, got {gammadot}")
    return float(6.0 * np.pi * mu_f * (a / 2.0) ** 3 * gammadot / (K_BOLTZ * T))


# --------------------------------------------------------------------------
# 3. the drainage (Deborah) number
# --------------------------------------------------------------------------
def drainage_number(tau_p: float, D0: float, U0: float) -> float:
    """De = tau_p / tau_impact, with tau_impact = D0/U0.

    De >> 1  UNDRAINED: pore pressure cannot relax during the impact, so an
             imposed dilation is fully frustrated and generates the largest
             particle pressure the network can carry.
    De << 1  DRAINED:   suction bleeds away as fast as dilation makes it, and
             there is no dilatant strengthening.
    De ~ 1   the interesting case, where the impact velocity selects the
             behaviour.
    """
    if U0 <= 0 or D0 <= 0:
        raise ValueError(f"D0 and U0 must be positive, got D0={D0}, U0={U0}")
    return float(tau_p / (D0 / U0))


# --------------------------------------------------------------------------
# 4. the dilatancy approximation
# --------------------------------------------------------------------------
def dilatant_pressure_darcy(sys: DilatantSystem, L_d: float,
                            gammadot: float) -> float:
    """Uncapped Reynolds-Darcy dilatant pressure, Eq. (D1)  [Pa].

        p_dil = (mu_f / kappa) * dphi * L_d^2 * gammadot

    Returns 0 below the dilatancy onset (dphi <= 0), where there is no dilatant
    branch at all. This is the quantity that returns 1e7-1e8 Pa for 800 nm
    silica; compare it against `dilatant_pressure_capped` before using it.
    """
    if sys.dphi <= 0:
        return 0.0
    if L_d <= 0 or gammadot < 0:
        raise ValueError(f"need L_d>0, gammadot>=0; got L_d={L_d}, gd={gammadot}")
    return float((sys.mu_w / sys.kappa) * sys.dphi * L_d**2 * gammadot)


def dilatant_pressure_ceiling(sys: DilatantSystem, K_p: float) -> float:
    """Stiffness ceiling on the dilatant pressure, Eq. (D3)  [Pa].

        p_dil_max = K_p * dphi / phi

    The packing can only generate pore suction while it is still dilating; once
    it has expanded from phi to phi_c there is no further dilation to frustrate.

    THE STRAIN IS dphi/phi, NOT dphi/(1-phi). Corrected 2026-08-10 while writing
    the derivation out in full. With incompressible grains the particle volume
    V_s = phi V is fixed, so V = V_s/phi and

        eps_v = dV/V = -dphi/phi,

    i.e. the mixture's volumetric strain is referenced to the SOLID fraction,
    not the porosity. The same phi appears in the definition of the modulus,
    K_p = phi dP/dphi, so the two are consistent by construction; using
    1/(1-phi) mixed two conventions. The error is O(1) -- exactly zero at
    phi = 0.5, 27% at phi = 0.44, 38% at phi = 0.58 -- and changes no
    conclusion, but it moves every printed ceiling.

    This replaces the cavitation cap the earlier brief proposed and withdrew.
    It is validated in the test suite against [JVF16]'s measured suction in
    170 um glass beads.
    """
    if sys.dphi <= 0:
        return 0.0
    return float(K_p * sys.dphi / sys.phi)


def dilatant_pressure_capped(sys: DilatantSystem, K_p: float, L_d: float,
                             gammadot: float) -> float:
    """min of the Darcy build-up (D1) and the stiffness ceiling (D3)  [Pa]."""
    return float(min(dilatant_pressure_darcy(sys, L_d, gammadot),
                     dilatant_pressure_ceiling(sys, K_p)))


def dilatant_pressure_from_identity(K_p: float, sys: DilatantSystem,
                                    gammadot: float, tau_p: float) -> float:
    """Eq. (D2): p_dil = K_p * dphi * gammadot * tau_p.

    Algebraically identical to (D1) when tau_p is computed at the SAME L_d.
    Exists so the tests can assert the two routes agree, which is the check
    that the derivation was transcribed correctly.
    """
    if sys.dphi <= 0:
        return 0.0
    return float(K_p * sys.dphi * gammadot * tau_p)


# --------------------------------------------------------------------------
# 5. packaged systems
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class PoroSystem:
    """A DilatantSystem plus the poroelastic inputs it needs."""

    sys: DilatantSystem
    K_p: float               # drained network bulk modulus [Pa]
    K_p_low: float           # bracket, low  [Pa]
    K_p_high: float          # bracket, high [Pa]
    K_p_origin: str          # provenance string, reprinted in every report

    @property
    def kappa(self) -> float:
        return self.sys.kappa

    def tau_p(self, L: float, K_p: float | None = None) -> float:
        return pore_relaxation_time(L, self.kappa,
                                    self.K_p if K_p is None else K_p,
                                    self.sys.phi)

    def D_p(self, K_p: float | None = None) -> float:
        return pore_diffusivity(self.kappa, self.K_p if K_p is None else K_p,
                                self.sys.phi)


P_CONF_SILICA = 44.2327
#   Confining particle pressure for the silica network [Pa]. Set equal to the
#   fitted Wyart-Cates onset stress, on the argument that the pressure at which
#   frictional contacts first form IS the pressure at which a contact network
#   first exists to have a modulus. Carries the same caveats as sigma* itself
#   (failed fit gate, phi-axis offset) -- see materials.py.
P_CONF_CORNSTARCH = 1.0e5
#   Confining particle pressure for cornstarch under impact [Pa], the scale of
#   the jammed-column stresses reported by Waitukaitis & Jaeger, Nature 487, 205
#   (2012). K_p ~ P^(1/3) here, so a decade of error in P is a factor 2.2 in K_p.
P_CONF_JVF = 1.0e4
#   [JVF16]'s own measured pore suctions are -1e4 to -1e5 Pa; 1e4 is the low end,
#   which is the conservative choice for a stiffness that grows with pressure.


def silica_poro(phi: float = 0.50, phi_c: float = 0.42,
                lam: float = DEBYE_DEFAULT_M,
                P_conf: float = P_CONF_SILICA) -> PoroSystem:
    """The lab's 800 nm Stoeber silica as a poroelastic medium.

    K_p from the electric-double-layer route, `network_modulus_edl`:

        K_p = P a / (6 lambda_D)

    which has no particle modulus in it -- correctly, since these contacts never
    touch. The bracket spans lambda_D = 1 nm (100 mM added salt, stiffest) to
    30 nm (nominal 0 mM, where screening is set by dissolved CO2 and the
    particles' own counterions, softest). That is exactly the range the lab's
    salt-screening series covers, so the bracket is a measurable axis rather
    than an ignorance interval.
    """
    s = lab_silica_dilatant(phi=phi, phi_c=phi_c)
    K_mid = network_modulus_edl(P_conf, s.a, lam)
    K_lo = network_modulus_edl(P_conf, s.a, 30.0e-9)
    K_hi = network_modulus_edl(P_conf, s.a, 1.0e-9)
    return PoroSystem(
        sys=s, K_p=K_mid, K_p_low=K_lo, K_p_high=K_hi,
        K_p_origin=(f"EDL affine modulus K_p = P({P_conf:g} Pa) a / (6 lambda_D),"
                    f" lambda_D = {lam*1e9:g} nm; bracket = lambda_D 30 nm (0 mM)"
                    " to 1 nm (100 mM). P_conf set to the fitted sigma*."))


def cornstarch_poro(phi: float = 0.42, E_particle: float = E_STARCH_PA,
                    nu: float = 0.3,
                    P_conf: float = P_CONF_CORNSTARCH) -> PoroSystem:
    """[M26]'s cornstarch-water as a poroelastic medium.

    K_p from the affine modulus with a HERTZIAN, pressure-dependent contact
    stiffness. The bracket is the literature spread on hydrated starch,
    E_p = 1e7 to 1e9 Pa, which because k_n ~ E^(2/3) is only a factor ~4.6 each
    way -- smaller than the spread on P_conf would suggest.
    """
    s = cornstarch_M26(phi=phi)

    def _K(E):
        return network_modulus_affine(
            contact_stiffness_hertz(E, nu, s.a, phi, P_conf), s.a, phi)

    return PoroSystem(
        sys=s, K_p=_K(E_particle), K_p_low=_K(1.0e7), K_p_high=_K(1.0e9),
        K_p_origin=(f"affine modulus with Hertz k_n at P_conf = {P_conf:g} Pa,"
                    f" E_p = {E_particle:.0e} Pa (bracket 1e7-1e9); "
                    f"delta/a = {hertz_overlap_ratio(E_particle, nu, s.a, phi, P_conf):.2e}"))


def glass_beads_JVF(phi: float = 0.58, phi_c: float = 0.575,
                    a: float = 170e-6, E_particle: float = 70.0e9,
                    nu: float = 0.25, P_conf: float = P_CONF_JVF) -> PoroSystem:
    """[JVF16]'s 170 um glass beads -- the one system with a MEASURED suction.

    This exists so the test suite can check Eq. (D3) against data rather than
    against another scaling argument. Jerome, Vandenberghe & Forterre report
    pore suctions of -1e4 to -1e5 Pa under impact.

    dphi is not reported directly; 0.005 is the dilation scale for a dense pack
    shearing from phi ~ 0.58 toward its critical state, and it is the single
    largest uncertainty in the comparison. The affine modulus is an upper bound
    (real packings relax non-affinely and measure 2-5x softer near
    isostaticity), so the predicted ceiling is expected to sit at or somewhat
    above the measured range, not below it.
    """
    s = DilatantSystem(name=f"glass beads 170 um (JVF16), phi={phi:g}", phi=phi,
                       phi_c=phi_c, a=a, rho=2500.0 * phi + 1000.0 * (1 - phi))
    K_mid = network_modulus_affine(
        contact_stiffness_hertz(E_particle, nu, a, phi, P_conf), a, phi)
    return PoroSystem(
        sys=s, K_p=K_mid, K_p_low=K_mid / 5.0, K_p_high=K_mid,
        K_p_origin=(f"affine modulus with Hertz k_n at P_conf = {P_conf:g} Pa,"
                    f" E_glass = {E_particle:.0e} Pa; low end = /5 for "
                    "non-affine relaxation. Affine is an upper bound."))


def K_p_bracket(p: PoroSystem) -> tuple[float, float]:
    """(low, high) network modulus for p, for propagating the uncertainty."""
    return (p.K_p_low, p.K_p_high)


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def _check_phi(phi: float) -> None:
    if not (0.0 < phi < 1.0):
        raise ValueError(f"phi must lie in (0,1), got {phi}")


def report(p: PoroSystem, D0: float, U0: float,
           K_p: float | None = None) -> dict:
    """Everything this module computes for one system, as a flat dict."""
    K = p.K_p if K_p is None else K_p
    s = p.sys
    L = drainage_lengths(D0, s.a, s.phi)
    gammadot = U0 / D0                    # convention (A), Shah/Arora/Driscoll
    tau_i = D0 / U0
    out = {
        "system": s.name,
        "phi": s.phi, "phi_c": s.phi_c, "dphi": s.dphi,
        "a_diameter_m": s.a,
        "kappa_m2": s.kappa,
        "K_p_Pa": K,
        "K_p_bracket_Pa": list(K_p_bracket(p)),
        "K_p_origin": p.K_p_origin,
        "K_osmotic_hard_sphere_Pa": osmotic_modulus_hard_sphere(s.phi, s.a),
        "storage_S_per_Pa": storage_coefficient(s.phi, K),
        "storage_fluid_fraction": ((1 - s.phi) * CHI_WATER
                                   / storage_coefficient(s.phi, K)),
        "D_p_m2_s": pore_diffusivity(s.kappa, K, s.phi),
        "tau_impact_s": tau_i,
        "gammadot_per_s": gammadot,
        "tau_p_s": {k: pore_relaxation_time(v, s.kappa, K, s.phi)
                    for k, v in L.items()},
        "De": {k: drainage_number(pore_relaxation_time(v, s.kappa, K, s.phi), D0, U0)
               for k, v in L.items()},
        "drainage_lengths_m": L,
        "p_dil_darcy_Pa": {k: dilatant_pressure_darcy(s, v, gammadot)
                           for k, v in L.items()},
        "p_dil_ceiling_Pa": dilatant_pressure_ceiling(s, K),
        "p_dil_capped_Pa": {k: dilatant_pressure_capped(s, K, v, gammadot)
                            for k, v in L.items()},
    }
    return out
