"""
materials.py — material property containers for drop-impact stress analysis.

Everything is SI:
    length  m
    time    s
    mass    kg
    stress / pressure   Pa
    density kg/m^3
    viscosity   Pa.s
    surface tension N/m

Design rule (see CLAUDE.md): no magic numbers buried in expressions. Every
material constant is a named field with a documented provenance string, so a
number can always be traced back to where it came from.

PROVENANCE OF THE SILICA SUSPENSION NUMBERS
-------------------------------------------
The Wyart-Cates parameters for the lab's 800 nm colloidal silica are NOT
invented here. They are the calibration that already lives in the drop-impact
repo, fitted to 8 measured steady flow curves:

    C:\\Users\\compl\\Github\\Drop Impact\\closure\\params_fitted.yaml
        eta_s      = 0.004736 Pa.s
        phi0       = 0.70      (frictionless jamming fraction, fit hit its cap)
        phiM       = 0.5176    (frictional jamming fraction)
        sigma_star = 44.2327 Pa
        delta      = 1.3852
        fit rms(log eta) = 0.3742

Two caveats travel with those numbers and must travel with every result
derived from them:

  (1) The fit FAILED its own acceptance gate (rms_ok requires < 0.35; it is
      0.374).  Status card: "Drop Impact - Simulation" section 6.
  (2) The rheometer dataset's volume-fraction axis is offset from this
      project's drop-sample volume fraction by an unquantified amount.
      Status card: "Drop Impact - Paper 1" section 8.  So sigma_star = 44.2 Pa
      is "the right order of magnitude for this material", not "the onset
      stress of the drops that were impacted".

Both caveats are carried in the `provenance` and `caveats` fields so that any
report generated from this module reprints them.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# --------------------------------------------------------------------------
# Physical constants
# --------------------------------------------------------------------------
G_ACCEL = 9.80665  # m/s^2, standard gravity


@dataclass(frozen=True)
class Fluid:
    """A simple Newtonian fluid (or a gas).

    Attributes
    ----------
    rho : density [kg/m^3]
    c   : speed of sound [m/s]
    eta : dynamic viscosity [Pa.s]
    gamma : surface tension against air [N/m]; None if not meaningful (gases)
    """

    name: str
    rho: float
    c: float
    eta: float
    gamma: float | None = None
    provenance: str = ""

    @property
    def Z(self) -> float:
        """Acoustic impedance Z = rho * c  [kg m^-2 s^-1 = Pa s/m].

        Z is the constant of proportionality between a pressure disturbance and
        the material velocity it drives:  p = Z u.  It is the single number
        that decides how much of an impact pressure a target *accepts* rather
        than reflects.
        """
        return self.rho * self.c

    @property
    def nu(self) -> float:
        """Kinematic viscosity [m^2/s]."""
        return self.eta / self.rho


# --- Reference fluids (20 degC unless stated) -----------------------------

WATER = Fluid(
    name="water (20 C)",
    rho=998.2,
    c=1482.0,
    eta=1.002e-3,
    gamma=72.8e-3,
    provenance="CRC Handbook / NIST, 20 degC, 1 atm",
)

AIR = Fluid(
    name="air (20 C, 1 atm)",
    rho=1.204,
    c=343.2,
    eta=1.82e-5,
    gamma=None,
    provenance="CRC Handbook, 20 degC, 1 atm",
)

# Amorphous / fused silica, bulk. Stoeber colloidal silica is porous and less
# dense than fused silica; the repo's calibration script uses 2.0 g/cm^3 for
# the lab's particles (closure/calibrate_from_rheometry.py:64).
SILICA_SOLID = Fluid(
    name="colloidal silica particle (porous)",
    rho=2000.0,
    c=5900.0,  # longitudinal sound speed of fused silica; used only in Wood's law
    eta=float("inf"),
    gamma=None,
    provenance=(
        "rho = 2.0 g/cm^3 from closure/calibrate_from_rheometry.py:64 "
        "('porous colloidal silica'); c = 5900 m/s is the fused-silica "
        "longitudinal wave speed (CRC) and enters only through Wood's law, "
        "where the answer is insensitive to it."
    ),
)


# --------------------------------------------------------------------------
# Wyart-Cates suspension
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class WyartCates:
    """Wyart-Cates shear-thickening suspension, in the repo's own parameterisation.

    Model
    -----
    Fraction of particle contacts that are frictional at steady stress sigma:

        f_ss(sigma) = exp[ -(sigma_star / sigma)^delta ]                    (1)

    Stress-dependent jamming fraction, linear in f:

        phi_J(f) = phi0 - (phi0 - phiM) f                                   (2)

    Effective viscosity (Krieger-Dougherty / Maron-Pierce form, exponent -2):

        eta(f) = eta_s [1 - phi/phi_J(f)]^(-2)                              (3)

    Steady state closes the loop:  sigma = eta(f_ss(sigma)) * gammadot.

    Physical reading
    ----------------
    sigma_star is the *dilatancy onset stress*: the shear stress at which the
    repulsive force holding particles apart is overcome, contacts become
    frictional, and the jamming fraction drops from phi0 toward phiM. Below
    sigma_star the suspension is a lubricated (frictionless) liquid; above it,
    a frictional one.  Everything in this project that says "did the impact
    thicken the drop?" is asking whether the local stress crossed sigma_star.

    delta controls how *sharp* that crossover is. delta -> infinity is a step
    function at sigma_star; delta = 1 is the original Wyart-Cates exponential.

    IMPORTANT CONSEQUENCE AT phi = 0.44
    -----------------------------------
    Shear jamming requires phi > phiM, i.e. phi_J(f=1) < phi.  With
    phiM = 0.5176 and phi = 0.44 the suspension can NEVER jam by this model,
    at any stress; it can only thicken continuously.  The drops in the
    experiment demonstrably do arrest.  That mismatch is real and is discussed
    in the report - it is either the phi-axis offset (caveat 2 above) or
    evidence that impact arrest is not bulk shear jamming.
    """

    name: str
    eta_s: float        # solvent (suspending fluid) viscosity [Pa.s]
    phi0: float         # frictionless jamming fraction [-]
    phiM: float         # frictional jamming fraction [-]
    sigma_star: float   # dilatancy onset stress [Pa]
    delta: float        # sharpness exponent [-]
    phi: float          # particle volume fraction [-]
    a: float            # particle RADIUS [m]
    rho_p: float        # particle density [kg/m^3]
    rho_f: float        # suspending-fluid density [kg/m^3]
    gamma: float        # surface tension of the suspension [N/m]
    provenance: str = ""
    caveats: tuple[str, ...] = ()

    # ---------------------------------------------------------------- checks
    def __post_init__(self):
        if not (0.0 < self.phi < self.phi0):
            raise ValueError(f"need 0 < phi < phi0; got phi={self.phi}, phi0={self.phi0}")
        if not (self.phiM < self.phi0):
            raise ValueError("need phiM < phi0 (frictional packs looser than frictionless)")
        if self.delta <= 0:
            raise ValueError("delta must be > 0")
        if self.sigma_star <= 0:
            raise ValueError("sigma_star must be > 0")

    # ------------------------------------------------------------ primitives
    def phi_J(self, f):
        """Jamming fraction at frictional-contact fraction f. Eq. (2)."""
        f = np.clip(np.asarray(f, dtype=float), 0.0, 1.0)
        return self.phi0 - (self.phi0 - self.phiM) * f

    def eta_of_f(self, f, eps_min: float = 0.01):
        """Effective viscosity at frictional fraction f [Pa.s]. Eq. (3).

        eps_min floors (1 - phi/phi_J) so the divergence is cut off at a finite
        jammed viscosity eta_s/eps_min^2 - this mirrors the repo's closure
        (finite contact stiffness), and only bites when phi -> phi_J.
        """
        gap = np.maximum(1.0 - self.phi / self.phi_J(f), eps_min)
        return self.eta_s * gap ** -2.0

    def f_ss(self, sigma):
        """Steady-state frictional fraction at shear stress sigma [Pa]. Eq. (1)."""
        sigma = np.asarray(sigma, dtype=float)
        out = np.zeros_like(sigma)
        pos = sigma > 0
        out[pos] = np.exp(-((self.sigma_star / sigma[pos]) ** self.delta))
        return out

    # --------------------------------------------------------- derived scales
    @property
    def eta_0(self) -> float:
        """Low-stress ("frictionless branch") viscosity of the suspension [Pa.s].

        This is the viscosity the drop has BEFORE the impact thickens it, and
        it is the correct viscosity to use in the impact Reynolds number.
        """
        return float(self.eta_of_f(0.0))

    @property
    def eta_frictional(self) -> float:
        """Fully-frictional (f = 1) viscosity [Pa.s] - the thickened branch."""
        return float(self.eta_of_f(1.0))

    @property
    def rho(self) -> float:
        """Suspension density by volume mixing [kg/m^3]."""
        return self.phi * self.rho_p + (1.0 - self.phi) * self.rho_f

    @property
    def c_wood(self) -> float:
        """Effective sound speed of the suspension by Wood's law [m/s].

            1 / (rho_eff c_eff^2) = phi/(rho_p c_p^2) + (1-phi)/(rho_f c_f^2)

        Wood's law is the low-frequency (wavelength >> particle) limit of a
        two-phase mixture: the phases share stress and add compliances. It is
        the standard estimate and is what sets the suspension's acoustic
        impedance for the water-hammer stage.

        NOTE this needs c_p and c_f, which are supplied by the caller through
        `sound_speed_wood`; this property uses fused silica + water defaults.
        """
        return sound_speed_wood(self.phi, SILICA_SOLID, WATER, self.rho)

    def sigma_at_f(self, f_target: float) -> float:
        """Shear stress at which the frictional fraction equals f_target [Pa].

        Invert Eq. (1):   sigma = sigma_star * [ -ln f ]^(-1/delta)
        """
        if not (0.0 < f_target < 1.0):
            raise ValueError("f_target must be strictly between 0 and 1")
        return self.sigma_star * (-np.log(f_target)) ** (-1.0 / self.delta)

    def gammadot_at_f(self, f_target: float) -> float:
        """Steady shear rate at which the frictional fraction equals f_target [1/s].

        Uses the steady-state closure sigma = eta(f) gammadot, evaluated
        self-consistently at f = f_target (so the viscosity used is the
        already-partly-thickened one, not the low-stress value).
        """
        sigma = self.sigma_at_f(f_target)
        return sigma / float(self.eta_of_f(f_target))

    @property
    def gammadot_c_naive(self) -> float:
        """Onset shear rate if the suspension had NOT yet thickened [1/s].

            gammadot_c = sigma_star / eta_0

        This is the number you get by asking "at what shear rate does the
        low-stress suspension first generate a stress equal to sigma_star?".
        It is an UPPER bound on the true onset rate in the sense that the real
        suspension thickens on the way up, so it reaches sigma_star at a
        somewhat lower applied rate. Both are reported in the results table.
        """
        return self.sigma_star / self.eta_0

    def flow_curve(self, sigma_grid=None):
        """Steady flow curve, parameterised by stress (single-valued in sigma).

        Returns dict of arrays: sigma [Pa], f [-], eta [Pa.s], gammadot [1/s].

        Parameterising by STRESS rather than by shear rate is deliberate: the
        Wyart-Cates flow curve can be S-shaped (multivalued in gammadot) but is
        always single-valued in sigma, so no branch selection or continuation
        is needed and the curve is exact rather than solver-dependent.
        """
        if sigma_grid is None:
            sigma_grid = np.logspace(-2, 4, 601)
        sigma = np.asarray(sigma_grid, dtype=float)
        f = self.f_ss(sigma)
        eta = self.eta_of_f(f)
        gammadot = sigma / eta
        return {"sigma": sigma, "f": f, "eta": eta, "gammadot": gammadot}


def sound_speed_wood(phi: float, particle: Fluid, fluid: Fluid, rho_mix: float) -> float:
    """Wood's-law effective sound speed of a two-phase mixture [m/s].

    1/(rho_mix c^2) = phi/(rho_p c_p^2) + (1-phi)/(rho_f c_f^2)

    Validity: wavelength >> particle size (true here by ~6 orders of
    magnitude), no resonant scattering, no gas phase.
    """
    if not (0.0 <= phi <= 1.0):
        raise ValueError("phi must be in [0, 1]")
    inv_K = phi / (particle.rho * particle.c**2) + (1.0 - phi) / (fluid.rho * fluid.c**2)
    return float(np.sqrt(1.0 / (rho_mix * inv_K)))


# --------------------------------------------------------------------------
# The lab's own system
# --------------------------------------------------------------------------
def lab_silica(phi: float = 0.44, a: float = 400e-9) -> WyartCates:
    """The Driscoll-lab 800 nm colloidal silica suspension.

    Parameters
    ----------
    phi : particle volume fraction. Default 0.44 = the salt-axis pivot.
    a   : particle RADIUS in m. Default 400 nm, i.e. 800 nm diameter.

    The suspending-fluid viscosity is the FITTED eta_s = 4.736 mPa.s, not
    1.0 mPa.s. That is ~4.7x water and is flagged in the simulation status
    card: either the solvent is a glycerol/water mix or the phi labels are
    offset. It is used as-fitted here so that the results are consistent with
    the repo's own calibration; the sensitivity to it is reported separately.
    """
    return WyartCates(
        name=f"lab silica, 800 nm, phi={phi:g}",
        eta_s=0.004736,
        phi0=0.70,
        phiM=0.5176,
        sigma_star=44.2327,
        delta=1.3852,
        phi=phi,
        a=a,
        rho_p=2000.0,
        rho_f=1000.0,
        gamma=72.0e-3,
        provenance=(
            "Wyart-Cates parameters from Github/Drop Impact/closure/"
            "params_fitted.yaml, fitted to 8 measured steady flow curves "
            "(2026-07-20). rho_p = 2.0 g/cm^3 from closure/"
            "calibrate_from_rheometry.py:64. Particle diameter 800 nm from "
            "Status/'Drop Impact - Paper 1' section 1."
        ),
        caveats=(
            "The flow-curve fit FAILED its own acceptance gate: "
            "rms(log eta) = 0.374 against a required < 0.35.",
            "The rheometer dataset's phi axis is offset from the drop "
            "samples' phi by an unquantified amount (Paper 1 card, section 8), "
            "so sigma_star = 44.2 Pa is an order-of-magnitude anchor for this "
            "material, not a measurement on the impacted drops.",
            "eta_s = 4.74 mPa.s is ~4.7x water; the solvent identity is "
            "unresolved (Simulation card, line 501).",
            "At phi = 0.44 < phiM = 0.5176 the model cannot shear-jam at any "
            "stress, yet the drops arrest. See report section on this.",
        ),
    )


def suspension_as_fluid(s: WyartCates) -> Fluid:
    """Wrap a WyartCates suspension as a Fluid using its LOW-STRESS viscosity.

    This is the drop as it arrives: unthickened, viscosity eta_0, density by
    volume mixing, sound speed by Wood's law. Use it for the inertial and
    acoustic stage estimates, which happen before (or while) the suspension
    thickens.
    """
    return Fluid(
        name=s.name + " (low-stress branch)",
        rho=s.rho,
        c=s.c_wood,
        eta=s.eta_0,
        gamma=s.gamma,
        provenance=(
            "density by volume mixing; sound speed by Wood's law; viscosity "
            "is the Wyart-Cates f=0 branch value eta_0 = eta_s (1-phi/phi0)^-2"
        ),
    )
