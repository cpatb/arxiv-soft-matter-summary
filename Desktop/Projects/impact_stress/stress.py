"""
stress.py — the analytical stress field under an impacting drop.

WHAT THIS MODULE IS
-------------------
An impacting drop does not apply "a stress". It applies a *sequence* of very
different stresses, each on its own length and time scale. This module computes
all of them for a given (drop, speed, target) and puts them on one axis.

Stages, in time order:

  Stage 0  GAS CUSHION    before contact; the air film carries the load
  Stage 1  WATER HAMMER   the compressible instant of first touch
  Stage 2  WAGNER         incompressible, expanding contact patch; the pressure
                          is self-similar and concentrated in a RING at the
                          moving turning line
  Stage 3  LAMELLA        the drop is a spreading sheet; the surviving stress
                          is viscous SHEAR, not pressure

PRIMARY SOURCES (equations below are taken from these, not invented)
--------------------------------------------------------------------
[PLJ16]  Philippi, Lagree & Antkowiak, "Drop impact on a solid surface:
         short-time self-similarity", J. Fluid Mech. 795, 96-135 (2016).
         arXiv:1504.05847.  (NOTE: the third author is Antkowiak, not
         Josserand -- a common mis-citation.)
         PDF on disk: Desktop/Papers/Soft Matter & Fluids/Drop Impact/
         _div_class__title__Drop_impact_on_a_solid_surface__short-time_self-similarity__div_.pdf
[CSG22]  Cheng, Sun & Gordillo, "Drop Impact Dynamics: Impact Force and Stress
         Distributions", Annu. Rev. Fluid Mech. 54, 57-81 (2022).
         doi:10.1146/annurev-fluid-030321-103941
         PDF on disk: .../Drop Impact/Drop Impact Dynamics_Impact Force and
         stress cheng review.pdf.  (Figure 1a of that review is Michelle
         Driscoll's own drop-impact image.)
Equation numbers in [square brackets] below are CSG22's numbering.

NOTATION (fixed throughout; note the two clocks)
-----------------------------------------------
    R    drop radius [m]          D = 2R is the diameter
    U    impact speed [m/s]       (CSG22 calls it V)
    rho  drop density [kg/m^3]    eta  drop viscosity [Pa.s]
    t    time since first contact [s]
    t*   = U t / R    dimensionless time on the RADIUS clock  (this module)
    ttil = U t / D    dimensionless time on the DIAMETER clock (CSG22)
    t*   = 2 * ttil.   Every factor-of-2 confusion in this subject is this.
    r    radial coordinate on the substrate [m]
    d(t) turning-line (contact) radius [m]
    xi   = r / d(t), the self-similar radial coordinate
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from materials import AIR, Fluid

# --------------------------------------------------------------------------
# Literature constants — each traceable to a printed equation
# --------------------------------------------------------------------------

WAGNER_PREFACTOR = 3.0
"""d(t) = sqrt(3 R U t).

CSG22 write the turning-line radius as r0 = sqrt(6 D V t)/2 (their Figure 2a,
verified experimentally by Gordillo et al. 2018). With D = 2R:
    sqrt(6 * 2R * U t)/2 = sqrt(12 R U t)/2 = sqrt(3 R U t).
Pure geometry (intersecting a falling sphere with a plane) would give
sqrt(2 R U t); the extra factor comes from Wagner's free-surface pile-up.
"""

PLJ_CUTOFF = 0.02
"""Screening cut-off Delta = 0.02 sqrt(V D t) near the turning line.

CSG22 Eq. [14] / PLJ16: the 1/x shear-stress singularity is screened within
Delta of the turning line, with Delta chosen from simulations at several Re.
In units of d: Delta/d = 0.02 sqrt(2 R U t)/sqrt(3 R U t) = 0.02 sqrt(2/3)
= 0.01633.  So the singular formulae are trusted only for xi < 1 - 0.0163.
"""

DELTA_OVER_D = PLJ_CUTOFF * np.sqrt(2.0 / 3.0)  # = 0.016330

F_MAX_TILDE = 0.85
T_MAX_TILDE = 0.17
"""Universal peak force and peak time, CSG22 Section 2.1:
    Fmax/(rho D^2 V^2) = 0.85 +/- 0.05   at   t V/D = 0.17 +/- 0.02
valid at high Re (> 200). Below that, use `peak_force_finite_Re`.
"""

MEAN_SHEAR_COEFF = 10.7 / (3.0 * np.pi)
"""1.13525 = (CSG22 Eq. [14] drag prefactor 5.35) x 2 / (3 pi).
The single source of truth for the contact-area-averaged shear
stress coefficient; imported by dilatancy.py rather than retyped."""

MEAN_SHEAR_COEFF_BAND = (0.76, 1.14)
"""Honest spread on MEAN_SHEAR_COEFF, from PLJ16 Table 2's own factor-2 range in
the screening length Delta (0.03-0.06 sqrt(ttil)) plus the +3.3% area-
normalisation inconsistency in the published construction. The default 1.135 is
the literature value and sits at the top of this band. U_c goes as C^(-2/3), so
the band is +8% to +31% on every onset velocity."""

RE_CORRECTION_COEFF = 8.0 * np.sqrt(6.0) / (3.0 * np.pi**1.5)
"""8 sqrt(6) / (3 pi^{3/2}) = 1.1734, the first-order viscous correction in
CSG22 Eq. [2]:   Ftilde = (3 sqrt6 / 2) (1 + 1.1734 Re^{-1/2}) ttil^{1/2}.
Re here is CSG22's Re = D V / nu (diameter-based).
"""


# --------------------------------------------------------------------------
# Targets
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Target:
    """What the drop hits.

    rho   : density [kg/m^3]. Sets the target's INERTIA in the incompressible
            stage. np.inf = immovable.
    Z     : acoustic impedance rho*c [Pa s/m]. Sets how much water-hammer
            pressure the target accepts. np.inf = perfectly rigid.
    E     : effective elastic modulus E* = E/(1-nu^2) [Pa]. np.inf = rigid,
            0 = a liquid (no shear rigidity).
    tau_y : yield stress [Pa] for yield-stress baths (Carbopol). 0 otherwise.
    """

    name: str
    rho: float
    Z: float
    E: float
    tau_y: float = 0.0
    provenance: str = ""


def rigid_target(name: str = "rigid solid (perfectly rigid)") -> Target:
    """The idealisation. NOTE that real glass is NOT this in the acoustic
    stage: fused silica has Z = 1.3e7 against the suspension's 2.3e6, so
    Z2/(Z1+Z2) = 0.85 and a glass slide takes 15% off the water-hammer
    pressure. Use `glass_target()` for that. In the incompressible stage
    glass really is rigid (Cauchy number 8e-8)."""
    return Target(name=name, rho=np.inf, Z=np.inf, E=np.inf,
                  provenance="idealisation, exact only in the Wagner stage")


def glass_target(name: str = "fused-silica microscope slide") -> Target:
    """A real glass slide: finite acoustic impedance, effectively infinite
    stiffness. rho = 2200 kg/m^3, c_long = 5900 m/s, E = 70 GPa."""
    return Target(name=name, rho=np.inf, Z=2200.0 * 5900.0, E=70e9 / (1 - 0.17**2),
                  provenance="fused silica; rho c = 1.30e7 Pa s/m, E = 70 GPa")


def liquid_target(fluid: Fluid, tau_y: float = 0.0, name: str | None = None) -> Target:
    return Target(name=name or f"deep bath of {fluid.name}", rho=fluid.rho,
                  Z=fluid.Z, E=0.0, tau_y=tau_y,
                  provenance=f"deep-pool limit (depth >> R) of {fluid.name}")


def elastic_target(rho: float, E: float, nu: float = 0.5,
                   c_long: float = 1480.0, name: str = "elastic solid") -> Target:
    """Soft elastic solid. E = Young's modulus; E* = E/(1-nu^2).

    For a nu -> 0.5 hydrogel the LONGITUDINAL wave speed is essentially the
    solvent's (~1480 m/s) even when the shear modulus is only kPa, because
    bulk compression is stiff while shear is soft. So a soft gel has a
    water-like acoustic impedance but a tiny elastic modulus: it is 'rigid'
    to the water hammer and 'soft' to the Wagner stage. That split is real
    and is one of the more useful things in this analysis.
    """
    return Target(name=name, rho=rho, Z=rho * c_long, E=E / (1.0 - nu**2),
                  provenance="E* = E/(1-nu^2); Z from the longitudinal wave speed")


# --------------------------------------------------------------------------
# The impact
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Impact:
    """A single impact condition: which drop, how fast, onto what."""

    drop: Fluid
    U: float
    R: float
    target: Target
    gas: Fluid = AIR

    # ------------------------------------------------------------- scales
    @property
    def D(self) -> float:
        """Drop diameter [m]."""
        return 2.0 * self.R

    @property
    def tau(self) -> float:
        """Inertial time R/U [s] — the radius clock this module runs on."""
        return self.R / self.U

    @property
    def tau_D(self) -> float:
        """Impact time D/U [s] — the diameter clock CSG22 uses."""
        return self.D / self.U

    @property
    def p_dyn(self) -> float:
        """rho U^2 [Pa] — the reference stress. CSG22 nondimensionalise by it."""
        return self.drop.rho * self.U**2

    # ------------------------------------------------- dimensionless groups
    @property
    def Re(self) -> float:
        """CSG22's Reynolds number Re = D U / nu = rho U D / eta."""
        return self.drop.rho * self.U * self.D / self.drop.eta

    @property
    def Re_R(self) -> float:
        """Radius-based Reynolds number rho U R / eta = Re/2."""
        return 0.5 * self.Re

    @property
    def We(self) -> float:
        """Weber number rho D U^2 / gamma."""
        if self.drop.gamma is None:
            return float("nan")
        return self.drop.rho * self.D * self.U**2 / self.drop.gamma

    @property
    def Oh(self) -> float:
        """Ohnesorge number eta / sqrt(rho gamma D)."""
        if self.drop.gamma is None:
            return float("nan")
        return self.drop.eta / np.sqrt(self.drop.rho * self.drop.gamma * self.D)

    @property
    def Fr(self) -> float:
        """Froude number U^2/(g D). Fr >> 1 means gravity is negligible."""
        return self.U**2 / (9.80665 * self.D)

    @property
    def Ma(self) -> float:
        """Mach number U/c."""
        return self.U / self.drop.c

    @property
    def St_gas(self) -> float:
        """Gas Stokes number, CSG22 Section 4.2: St = eta_gas/(rho U D)."""
        return self.gas.eta / (self.drop.rho * self.U * self.D)

    def t_of_tstar(self, t_star):
        """Convert dimensionless radius-clock time to seconds."""
        return np.asarray(t_star, dtype=float) * self.tau


# ==========================================================================
# STAGE 1 — WATER HAMMER
# ==========================================================================
def water_hammer_pressure(imp: Impact) -> float:
    """Peak contact pressure of the compressible instant [Pa].

    DERIVATION
    ----------
    At first contact, no signal has yet crossed the contact patch, so the
    problem is one-dimensional: two half-spaces approaching at relative speed
    U. Velocity and stress are continuous across the interface, and each side
    responds acoustically (p = Z * velocity change on that side):

        p = Z1 u1 = Z2 u2,   u1 + u2 = U   =>   p = U Z1 Z2 / (Z1 + Z2)   (WH)

    Limits (all in the test suite):
        rigid  Z2 -> inf : p = rho1 c1 U      <- the textbook water hammer,
                                                 CSG22 quote it as "~rho V c"
        matched Z2 = Z1  : p = (1/2) rho c U  <- a bath HALVES it
        free surface     : p = 0

    So "the water-hammer pressure" is really the RIGID-WALL SPECIAL CASE of an
    impedance-sharing law. That generalisation is the first of the two results
    this analysis needs for non-rigid substrates.
    """
    Z1, Z2 = imp.drop.Z, imp.target.Z
    return Z1 * imp.U if np.isinf(Z2) else imp.U * Z1 * Z2 / (Z1 + Z2)


def shock_detachment(imp: Impact) -> dict:
    """When and where the compressible stage ends.

    Geometry: a sphere descending at U wets a circle of radius
    sqrt(2 R U t) at short time, so the contact edge sweeps outward at
    sqrt(R U / 2t) — infinite at t = 0, decreasing. While that exceeds the
    sound speed the compressed liquid is trapped and the water-hammer
    pressure is realised. Setting it equal to c:

        t_d = R U / (2 c^2)        d_d = R U / c

    CSG22 (Section 1) quote exactly these: duration "~R V/c^2", contact radius
    "~R V/c ~ 10 um" — our t_d carries the explicit factor 1/2.

    The affected volume fraction is (3/8)(U/c)^3 — parts per billion. The
    water-hammer pressure is the largest number in the whole problem and it
    touches essentially none of the drop. CSG22 add that air cushioning
    screens it further, so it is "hard to detect experimentally".
    """
    c, U, R = imp.drop.c, imp.U, imp.R
    t_d = R * U / (2.0 * c**2)
    d_d = R * U / c
    return {
        "t_d": t_d,
        "t_d_star": t_d / imp.tau,
        "d_d": d_d,
        "d_d_over_R": d_d / R,
        "p": water_hammer_pressure(imp),
        "volume_fraction_affected": 3.0 / 8.0 * (U / c) ** 3,
    }


# ==========================================================================
# STAGE 2 — WAGNER: pressure
# ==========================================================================
def contact_radius(imp: Impact, t):
    """Turning-line radius d(t) = sqrt(3 R U t) [m].  CSG22 Fig. 2a."""
    return np.sqrt(WAGNER_PREFACTOR * imp.R * imp.U * np.asarray(t, dtype=float))


def contact_line_speed(imp: Impact, t):
    """ddot = (1/2) sqrt(3 R U / t) [m/s]."""
    t = np.asarray(t, dtype=float)
    return 0.5 * np.sqrt(WAGNER_PREFACTOR * imp.R * imp.U / t)


def wagner_pressure(imp: Impact, t, r):
    """Pressure on the wall, p(r,t) [Pa].  CSG22 Eq. [10] (from PLJ16).

    CSG22 print it as

        p(r,t) = (3 sqrt2 / 2 pi) rho V^2 D / sqrt(3 V D t - 2 r^2)      [10]

    which, with D = 2R and d^2 = 3 R U t (so 3VDt = 2 d^2), is identically

        p(r,t) = (3/pi) rho U^2 R / sqrt(d(t)^2 - r^2)
               = (sqrt3 / pi) rho U^2  t*^(-1/2)  (1 - xi^2)^(-1/2)       (5)

    HOW TO READ EQUATION (5) — three factors, three separate ideas:
      rho U^2            the natural stress scale
      t*^(-1/2)          the amplitude: the earlier you look, the larger it is
      (1-xi^2)^(-1/2)    the SHAPE: flat over the contact patch, divergent at
                         its rim. THIS is the stress localisation.

    The field is self-similar — one universal shape times a decaying
    amplitude — which is what PLJ16's title means by "short-time
    self-similarity", and what makes the pressure maximum sit OFF-AXIS, in a
    ring that races outward with the turning line. CSG22 stress that this is
    the opposite of a steady impinging jet, whose maximum is on the axis.

    Returns NaN where r >= d(t): outside the contact patch this solution does
    not apply.
    """
    t = np.asarray(t, dtype=float)
    r = np.asarray(r, dtype=float)
    d2 = WAGNER_PREFACTOR * imp.R * imp.U * t
    arg = d2 - r**2
    with np.errstate(invalid="ignore", divide="ignore"):
        p = (3.0 / np.pi) * imp.drop.rho * imp.U**2 * imp.R / np.sqrt(arg)
    return np.where(arg > 0, p, np.nan)


def wagner_pressure_selfsimilar(imp: Impact, t_star, xi):
    """Eq. (5) in self-similar variables. Same numbers, clearer structure."""
    t_star = np.asarray(t_star, dtype=float)
    xi = np.asarray(xi, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = ((np.sqrt(3.0) / np.pi) * imp.drop.rho * imp.U**2
               * t_star**-0.5 / np.sqrt(1.0 - xi**2))
    return np.where(np.abs(xi) < 1.0, out, np.nan)


def wagner_pressure_axis(imp: Impact, t_star, include_stagnation: bool = False):
    """On-axis pressure p(0,t) [Pa].

        singular part:  (sqrt3/pi) rho U^2 t*^(-1/2) = 0.5513 rho U^2 t*^(-1/2)
        CSG22 Eq. [9]:  p(0,t) = (1/2) rho U^2 + [the singular part]

    CSG22 call the constant (1/2) rho U^2 "a constant dynamic pressure that
    can be absorbed by redefining the reference pressure", and Eq. [10] omits
    it. That is fine for the SHAPE and for anything deviatoric, and it is
    what the default does. It is NOT fine if you want the gauge pressure on
    the wall late in the impact: the omitted term is 9% of the total at
    t* = 0.01, 29% at t* = 0.1 and 50% at t* = 0.3. Pass
    include_stagnation=True there.

    The singular part alone equals rho U^2 at t* = 3/pi^2 = 0.304.
    """
    p = ((np.sqrt(3.0) / np.pi) * imp.drop.rho * imp.U**2
         * np.asarray(t_star, dtype=float) ** -0.5)
    return p + 0.5 * imp.p_dyn if include_stagnation else p


def wagner_pressure_peak(imp: Impact, t_star):
    """Pressure at the screened rim [Pa] — the actual peak stress.

    Eq. (5) diverges at xi = 1. PLJ16 find the singularity is screened within
    Delta = 0.02 sqrt(V D t) of the turning line (CSG22 Eq. [14]), i.e.
    1 - xi_cut = DELTA_OVER_D = 0.01633. Evaluating (5) there:

        p_peak = 3.063 rho U^2 t*^(-1/2) = 5.556 x the on-axis value

    THE STRESS IS IN A RING. But the amplification factor is not a constant of
    nature, and three things limit how far it can be quoted:

    1. It is time-independent ONLY because Delta is taken proportional to
       sqrt(t). That is an empirical fit; PLJ16 say the origin of Delta "needs
       to be clearly identified". Their two candidate mechanisms disagree
       qualitatively -- a viscous length gives an amplification that GROWS with
       time and depends on Re, an inertial Wagner inner layer gives one that is
       UNBOUNDED as t -> 0 -- and both are excluded by PLJ16's own Table 2,
       which shows Delta/sqrt(ttil) varying by less than 2x across Re 250-5000
       and ttil 5e-3 to 0.1. So the sqrt(t) fit survives, unexplained, over
       that window only.
    2. Inside PLJ16's own Delta range (0.03-0.06 sqrt(ttil)) the ratio is
       3.8-5.4, not 5.56. Quote 4-6.
    3. Below t* ~ 2e-5 the rim is capped by compressibility. The AXIS crossover
       with the water hammer is at t* = 6e-7; the rim, being 5.6x higher, hits
       the same ceiling ~31x later in t*, and between those two times the
       amplification falls toward 1.

    Correct statement: p_peak/p_axis = 4-6 for 5e-3 < t* < 0.1, set by an
    empirical screening length of unknown physical origin, factor-2 uncertain,
    and not established outside that window.
    """
    return wagner_pressure_selfsimilar(imp, t_star, 1.0 - DELTA_OVER_D)


def wagner_hammer_crossover(imp: Impact) -> dict:
    """Where the Wagner field's blow-up is cut off by compressibility.

    Setting the on-axis Wagner pressure equal to the water-hammer pressure:

        (sqrt3/pi) rho U^2 t*^(-1/2) = rho c U  =>  t*_x = (3/pi^2) Ma^2

    Compare the independent geometric estimate t*_d = Ma^2/2 from
    `shock_detachment`. Two unrelated arguments — one about wave speeds, one
    about matching pressure formulae — agree to a factor 1.64.

    STRUCTURAL POINT: the compressible stage is not a separate phenomenon
    bolted on at t = 0. It is what makes the incompressible solution finite.
    """
    Ma = imp.Ma
    return {
        "Ma": Ma,
        "t_star_pressure_match": (3.0 / np.pi**2) * Ma**2,
        "t_star_shock_detach": 0.5 * Ma**2,
        "p_ceiling": water_hammer_pressure(imp),
    }


# ==========================================================================
# STAGE 2/3 — SHEAR STRESS (the one that matters for a suspension)
# ==========================================================================
def wall_shear_stress(imp: Impact, t, r):
    """Wall shear stress tau(r,t) [Pa].  CSG22 Eq. [13] (from PLJ16).

    CSG22 print it as

        tau(r,t) = [sqrt6 rho^(1/2) V^(3/2) D^(1/2) eta^(1/2) / pi^(3/2)]
                   * 2r / (3 V D t - 2 r^2)                              [13]

    With D = 2R and 3VDt - 2r^2 = 2(d^2 - r^2) this is identically

        tau(r,t) = (2/pi^(3/2)) rho U^2 Re_R^(-1/2) t*^(-1/2) xi/(1-xi^2)  (13')

    with Re_R = rho U R / eta.  WHY THIS FORM IS THE USEFUL ONE:

      * PREFACTOR rho U^2 Re_R^(-1/2).  The shear stress is smaller than the
        pressure by ~1/sqrt(Re). At the lab's Re_R = 121, the area-averaged
        shear stress is 5.3x below the on-axis pressure and 9.7x below
        rho U^2 at t* = 1 (5.5x at t* = 0.315). CORRECTED AFTER REVIEW: an
        earlier version of this docstring said '~100x', which is wrong by an
        order of magnitude and contradicted dilatancy.py. The right
        statement is that rho U^2 overstates the stress that drives
        thickening by roughly one decade, which is the difference between a
        gauge whose threshold is reachable and one whose is not.

      * SHAPE xi/(1-xi^2).  A 1/(distance from the rim) divergence — STRONGER
        than the pressure's 1/sqrt. The shear stress is more sharply localised
        than the pressure, and unlike the pressure it is exactly ZERO on the
        axis (by axisymmetry). CSG22 note the shear maximum stays near the
        drop's edge for the entire impact, while the pressure maximum migrates
        to the centre.

      * tau ~ sqrt(eta).  Not eta. So a factor-5 uncertainty in the drop's
        viscosity is only a factor-2.2 uncertainty in the stress — which is
        what makes this gauge usable despite the lab's unresolved eta_0.

    Validity: high Re, early time, Newtonian eta. For a shear-thickening drop
    this is the stress the impact WOULD impose on the unthickened fluid — i.e.
    exactly the right thing to compare against the dilatancy onset stress to
    ask "does the impact trigger thickening at all?". Once it does, eta rises
    and the true stress rises with sqrt(eta): the feedback is positive but
    weak, so the ONSET criterion is self-consistent even though the
    post-onset magnitude is not.
    """
    t = np.asarray(t, dtype=float)
    r = np.asarray(r, dtype=float)
    d2 = WAGNER_PREFACTOR * imp.R * imp.U * t
    arg = d2 - r**2
    pref = (np.sqrt(6.0) / np.pi**1.5) * np.sqrt(
        imp.drop.rho * imp.U**3 * imp.D * imp.drop.eta)
    with np.errstate(invalid="ignore", divide="ignore"):
        tau = pref * 2.0 * r / (2.0 * arg)
    return np.where(arg > 0, tau, np.nan)


def wall_shear_selfsimilar(imp: Impact, t_star, xi):
    """Eq. (13') in self-similar variables [Pa].

        tau = (2/pi^(3/2)) rho U^2 Re_R^(-1/2) t*^(-1/2) xi/(1-xi^2)
    """
    t_star = np.asarray(t_star, dtype=float)
    xi = np.asarray(xi, dtype=float)
    amp = ((2.0 / np.pi**1.5) * imp.drop.rho * imp.U**2
           * imp.Re_R**-0.5 * t_star**-0.5)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = amp * xi / (1.0 - xi**2)
    return np.where((np.abs(xi) < 1.0) & (xi >= 0), out, np.nan)


def wall_shear_peak(imp: Impact, t_star):
    """Shear stress at the screened rim [Pa] — the largest shear in the drop.

    Evaluating (13') at xi = 1 - DELTA_OVER_D:
        xi/(1-xi^2) = 0.9837/0.03240 = 30.36
        tau_peak = 10.90 rho U^2 Re_R^(-1/2) t*^(-1/2)

    The rim runs 30x the amplitude — a much stronger localisation than the
    pressure's 5.5x, because the shear singularity is 1/x rather than
    1/sqrt(x). It is also the part most sensitive to the screening length,
    so quote it as an upper bound, not a measurement.
    """
    return wall_shear_selfsimilar(imp, t_star, 1.0 - DELTA_OVER_D)


def wall_shear_mean(imp: Impact, t_star):
    """Contact-area-averaged shear stress [Pa] — the ROBUST gauge.

    Rather than pick a point, integrate. PLJ16 integrate Eq. [13] over the
    contact area out to the screening cut-off to get the drag force
    (CSG22 Eq. [14]):

        F_d(t) = 5.35 eta^(1/2) rho^(1/2) V^2 D sqrt(t)                  [14]

    Dividing by the contact area pi d^2 = pi * 3 R U t:

        <tau> = 5.35 sqrt(eta rho) U^2 (2R) sqrt(t) / (3 pi R U t)
              = (10.7 / 3 pi) sqrt(eta rho) U / sqrt(t)
              = 1.135 rho U^2 Re_R^(-1/2) t*^(-1/2)                      (14')

    HOW MUCH TO TRUST THE PREFACTOR. Not as much as "5.35" suggests.
      * The integral DOES NOT converge. int xi^2/(1-xi^2) dxi is
        logarithmically divergent -- PLJ16 say so in as many words -- so
        <tau> depends on the screening length, logarithmically. A factor 1.4
        in Delta is 11% in F_d; PLJ16 Table 2's own factor-2 spread in
        Delta is 24%.
      * The published 5.35 is internally inconsistent with the Delta that
        PLJ16 report. Their stated Delta/sqrt(ttil) ~ 0.03 (Table 2:
        0.03-0.06 over Re 250-5000) integrates to a coefficient of
        1.02 here, not 1.135; recovering 1.135 needs Delta = 0.0202
        sqrt(ttil), below their whole quoted range.
      * Integrating F_d out to d-Delta but dividing by the full area pi d^2
        (as the published construction does) is a further +3.3%.
    So: use C = 1.0 with a band 0.76-1.14, not 1.1353 to five figures. The
    published value is kept as the default so the number is traceable, and
    `MEAN_SHEAR_COEFF_BAND` carries the honest spread. U_c goes as
    C^(-2/3), so this is +8% to +31% on every onset velocity.
    """
    coeff = MEAN_SHEAR_COEFF
    return (coeff * imp.drop.rho * imp.U**2 * imp.Re_R**-0.5
            * np.asarray(t_star, dtype=float) ** -0.5)


def shear_rate(imp: Impact, t_star, xi):
    """Local shear rate gammadot = tau/eta [1/s] from Eq. (13')."""
    return wall_shear_selfsimilar(imp, t_star, xi) / imp.drop.eta


def shear_rate_mean(imp: Impact, t_star):
    """Contact-averaged shear rate <tau>/eta [1/s] from Eq. (14')."""
    return wall_shear_mean(imp, t_star) / imp.drop.eta


def naive_impact_shear_rate(imp: Impact) -> float:
    """gammadot_impact = U/D — the convention of Shah, Arora & Driscoll (2022).

    Their argument (Comm. Phys. 5, 222 (2022), Methods): at the instant of
    impact the bottom of the drop is at rest while the apex still falls at U,
    so U/D is 'an estimate of the maximum shear rate at the moment of impact'.
    They quote 233 - 1333 1/s for D = 3 mm. It is the abscissa of the
    front-speed collapse and therefore of the lab's fitted gammadot_c.

    Kept explicit so that every comparison states its convention. The same
    paper ALSO uses gammadot_spr = U/l with l ~ 100 um (giving ~1e4 1/s) for
    the effective Reynolds number — a factor ~30 apart. Section 4 of the
    report shows these are the same Eq. (13') shear rate evaluated at two
    different times, which reconciles them.
    """
    return imp.U / imp.D


# ==========================================================================
# FORCE (whole-drop, measurable with a piezo sensor)
# ==========================================================================
def impact_force(imp: Impact, t_star):
    """Early-time impact force F(t) [N].  CSG22 Eq. [1] with the [2] correction.

        Ftilde = (3 sqrt6 / 2) (1 + 1.1734 Re^(-1/2)) ttil^(1/2)
        F = Ftilde * rho D^2 U^2 ,   ttil = t U / D = t*/2

    The Re^(-1/2) term is Gordillo et al.'s (2018) matched-asymptotic
    correction: the viscous boundary layer speeds up the propagation of the
    self-similar pressure field, so a viscous drop hits HARDER, not softer.
    At Re = 200 it is +8%; at Re = 20, +26%.
    """
    ttil = 0.5 * np.asarray(t_star, dtype=float)
    Ft = (3.0 * np.sqrt(6.0) / 2.0) * (1.0 + RE_CORRECTION_COEFF * imp.Re**-0.5) * np.sqrt(ttil)
    return Ft * imp.drop.rho * imp.D**2 * imp.U**2


def peak_force(imp: Impact) -> dict:
    """Peak impact force and when it happens.  CSG22 Eqs. [3], [4].

    High-Re limit: Fmax = 0.85 rho D^2 U^2 at ttil = 0.17.
    Finite Re: solve  (1 - ttil_max)^4 - A^4 (1 + Utl Re^(-1/2))^2 ttil_max = 0
    with A = 1.29 and Utl = 8 sqrt6/(3 pi^{3/2}), then

        Fmax(Re) = Fmax_inf * ttil_max_inf / ttil_max(Re)                [4]

    NOTE THE SENSE OF EQ. [4] — it is easy to get upside-down, and the physics
    fixes it three ways:
      * CSG22's text: "When Re is below 200, Fmax INCREASES and tmax DECREASES
        with decreasing Re." So Fmax and tmax move OPPOSITELY.
      * The force-time curve's shape is invariant and the area under it is the
        drop's fixed initial momentum, so Fmax * tmax = const.
      * A viscous drop hits harder, not softer (the same Re^(-1/2) correction
        that raises `impact_force`).
    """
    A = 1.29
    k = (1.0 + RE_CORRECTION_COEFF * imp.Re**-0.5) ** 2

    # solve the quartic on (0,1) by bisection — monotone decreasing there
    def g(x):
        return (1.0 - x) ** 4 - A**4 * k * x

    lo, hi = 1e-12, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(mid) > 0:
            lo = mid
        else:
            hi = mid
    ttil_max = 0.5 * (lo + hi)
    F_tilde_max = F_MAX_TILDE * T_MAX_TILDE / ttil_max
    return {
        "t_tilde_max": ttil_max,
        "t_star_max": 2.0 * ttil_max,
        "t_max_s": ttil_max * imp.tau_D,
        "F_tilde_max": F_tilde_max,
        "F_max_N": F_tilde_max * imp.drop.rho * imp.D**2 * imp.U**2,
        "mean_pressure_at_peak_Pa": (
            F_tilde_max * imp.drop.rho * imp.D**2 * imp.U**2
            / (np.pi * (contact_radius(imp, ttil_max * imp.tau_D)) ** 2)),
    }


# ==========================================================================
# STAGE 0 — GAS CUSHION
# ==========================================================================
def gas_film(imp: Impact) -> dict:
    """Air-film scales just before touchdown.

    CSG22 Section 4.2 (citing Josserand et al. 2016) give the draining-regime
    gas pressure directly:

        p_gas ~ rho U^2 St^(-1/3),   St = eta_gas / (rho U D)

    The companion geometric scalings (Mandre/Mani/Brenner; Bouwhuis et al.)
    are h_dimple/R ~ St^(2/3), l_dimple/R ~ St^(1/3). These are SCALINGS with
    O(1) prefactors omitted — enough to answer "does the gas film matter for
    the suspension?", not enough to quote to two figures.

    CSG22's structural point: air cushioning SCREENS both the water-hammer
    pressure and the 1/sqrt(t) pressure singularity, and moves the gas-pressure
    maximum off-axis to the contact ring — the same ring geometry as the
    liquid-phase result.
    TWO STOKES NUMBERS, AND THEY ARE NOT INTERCHANGEABLE. CSG22 define
    St_D = eta_g/(rho U D) and quote p_gas ~ rho U^2 St_D^(-1/3) with it.
    Mandre/Mani/Brenner and Bouwhuis define St_R = eta_g/(rho U R), and the
    geometric scalings h ~ R St^(2/3), l ~ R St^(1/3) belong to THAT one.
    Using the wrong one shifts h by 2^(2/3) = 1.59 and l by 1.26. Both are
    returned, each used only with its own source's scaling.
    """
    St_D = imp.St_gas                    # CSG22 convention, diameter-based
    St_R = 2.0 * St_D                    # MMB convention, radius-based
    h_d = imp.R * St_R ** (2.0 / 3.0)
    l_d = imp.R * St_R ** (1.0 / 3.0)
    return {
        "St_gas_D_based_CSG22": St_D,
        "St_gas_R_based_MMB": St_R,
        "h_dimple": h_d,
        "l_dimple": l_d,
        "p_gas": imp.drop.rho * imp.U**2 * St_D ** (-1.0 / 3.0),
        "t_dimple_star": (h_d / imp.U) / imp.tau,
        "gas_shear_stress": imp.gas.eta * imp.U / h_d,
        "_precision": ("SCALINGS ONLY - order-unity prefactors omitted. "
                       "Do not quote p_gas to more than one significant figure."),
    }


# ==========================================================================
# TARGET COMPLIANCE
# ==========================================================================
def compliance_factors(imp: Impact) -> dict:
    """How much a non-rigid target reduces the peak stress.

    (a) COMPRESSIBLE STAGE — acoustic impedance sharing, from Eq. (WH):

            p / p_rigid = Z2 / (Z1 + Z2)                                  (12)

    (b) INCOMPRESSIBLE (WAGNER) STAGE — inertia sharing.
        Redo the Wagner derivation with a deformable target. The approach
        speed is split, u1 + u2 = U, and each side develops the SAME
        interfacial pressure through its own potential flow, so

            rho1 u1 = rho2 u2  =>  u1 = U rho2/(rho1+rho2)
            p / p_rigid = rho2 / (rho1 + rho2)                            (13)

        Equations (12) and (13) are the SAME LAW with Z -> rho: two media
        sharing an interface share the load in proportion to their inertia,
        with 'inertia' meaning acoustic impedance when waves matter and mass
        density when they do not. For a drop hitting a bath of ITSELF, both
        give exactly 1/2.

        THE 3 DOES NOT CHANGE, and an earlier version of this note wrongly
        said it did. In the linearised Wagner problem each side's free
        surface responds with the same shape function at the same d, so the
        Wagner condition involves only the RELATIVE approach u1 + u2 = U and
        is identical to the one-sided problem. Cimpeanu & Moore, JFM 856,
        764 (2018) Eq. (3.17) confirm it: a drop on a rigid wall and a drop
        on a deep flat pool give an identical d(t). For a CURVED target it
        becomes a reduced radius, 1/R_eff = 1/R1 + 1/R2.

        The matching condition is continuity of PRESSURE IMPULSE, rho*phi,
        not of the potential. Both sides solve the same flat-disc mixed
        boundary-value problem and therefore carry the same universal 2/pi,
        which is what makes rho1 u1 = rho2 u2 follow rather than be assumed.
        Cimpeanu & Moore Eq. (3.13) give the interface coefficient (1+V)/2,
        equal to 1 at a rigid wall and 1/2 at a deep pool: a published
        confirmation of the factor 1/2 at rho1 = rho2. The generalisation to
        rho1 != rho2 is not published there and is derived here.

        DO NOT CONFLATE THIS WITH Ve = V/2. Langley et al. and Bisighini et
        al. report the interface receding at V/2, which is the same u1 -- but
        substituting V -> V/2 into the Wagner pressure would give 2^(-3/2) =
        0.354, or 0.25 on a rho V^2 gauge. The answer is 1/2 because d(t) is
        set by the RELATIVE approach U while the amplitude is set by u1, so
        only one velocity factor is halved.

        AND THE SHEAR STRESS DOES SOMETHING ELSE ENTIRELY -- see
        `shear_sharing_factor`. It does not halve; for a bath of the same
        liquid it VANISHES.

    (c) ELASTIC SOLID — a deflection criterion, not a sharing law.
        The target deflects w ~ p d / E*, which matters once w is comparable
        to the gap the drop is closing at the same radius, h_gap ~ d^2/(2R).
        With the Wagner pressure this reduces to

            t*  <~  (2/pi) C ,      C = rho U^2 / E*                      (15)

        C is the CAUCHY NUMBER (Yokoyama et al., Nat. Commun. 17, 607 (2026)
        name it that). The 2/pi comes out of the algebra and must not be
        dropped: with p = (3/pi) rho U^2 R/d and d = R sqrt(3 t*), the
        deflection w = p d/E* = 3 rho U^2 R/(pi E*) is TIME-INDEPENDENT while
        the gap being closed is h_gap = (3/2) R t*; equating them gives
        t* = (2/pi) C.

        C >~ 1: soft for the whole impact.  C << 1: rigid except for the first
        (2/pi)C-th of it. For this impact: fused silica C = 8e-8, a 10 kPa PDMS
        gel C = 0.6.

        HOWLAND'S CRITERION IS THIS ONE. Howland et al. PRL 117, 184502 (2016)
        write p_peak >~ E with p_peak = rho R U/(t_ej - t_i), which is
        algebraically rho U^2/t*_ej -- i.e. "C >~ t*_ej". Their worked example
        (R = 1.5 mm, U = 2.36 m/s, t_ej - t_i = 30 us, rho = 789 kg/m^3)
        reproduces their quoted E_c = 93 kPa exactly. Note they use E, not
        E* = E/(1-nu^2), which differ by 1.33 at nu = 0.5.

        C IS NECESSARY, NOT SUFFICIENT. Yokoyama et al. show the peak force
        does NOT collapse on Ca alone: the inertial-to-Hertzian crossover needs
        a second group theta = eta/sqrt(rho E R^2) built from the DROP's
        viscosity, through Z = theta^0.28 Ca, and they classify the problem as
        self-similarity of the second kind, so the 0.28 is not obtainable by
        dimensional analysis.

    (d) YIELD-STRESS BATH — the lab's own index.

            J = tau_y / (rho U^2)                                         (16)

        J << 1: the bath yields everywhere the impact touches, so it is
        liquid-like and Eq. (13) applies.  J >~ 1: it never yields and behaves
        as an elastic solid, so Eq. (15) applies with E* the elastic modulus.
        J is therefore the switch BETWEEN the two laws, which is a cleaner
        statement of what the Carbopol ladder is scanning than 'stiffness'.
    """
    Z1, Z2 = imp.drop.Z, imp.target.Z
    rho1, rho2 = imp.drop.rho, imp.target.rho

    # WHICH LAW APPLIES. The inertia-sharing law (13) is derived for a target
    # that FLOWS, so it is meaningful only for a liquid (E* = 0). A solid does
    # not share inertia with the drop; it deflects elastically, and its
    # response is set by the Cauchy number C, not by its density. Reporting
    # rho2/(rho1+rho2) for a solid would be a silently wrong number, so it is
    # labelled here rather than left for a reader to misuse.
    if np.isinf(imp.target.E):
        law = "rigid: no reduction"
    elif imp.target.E == 0.0:
        law = "liquid target: use inertial_factor (Eq. 13); Cauchy number not meaningful"
    else:
        law = "elastic solid: use compliance_number_C (Eq. 15); inertial_factor NOT applicable"

    # The inertia-sharing law is derived for a target that FLOWS. Emitting a
    # number for a solid would be a silent 2.5x error (a stiff PDMS block would
    # come back as 0.40 when the answer is 1.0), so the switch is implemented,
    # not merely documented. NaN, not a plausible-looking float.
    if np.isinf(rho2):
        f_inertial = 1.0
    elif imp.target.E == 0.0:
        f_inertial = rho2 / (rho1 + rho2)
    else:
        f_inertial = float("nan")

    return {
        "applicable_law": law,
        "acoustic_factor": 1.0 if np.isinf(Z2) else Z2 / (Z1 + Z2),
        "inertial_factor": f_inertial,
        # C = rho U^2 / E*, the Cauchy number. E* = inf (rigid) -> C = 0;
        # E* = 0 (a liquid, no shear rigidity) -> C = inf, i.e. infinitely
        # compliant, which is the correct statement and must not divide by zero.
        "compliance_number_C": (
            0.0 if np.isinf(imp.target.E)
            else np.inf if imp.target.E == 0.0
            else imp.p_dyn / imp.target.E),
        "yield_number_J": imp.target.tau_y / imp.p_dyn,
        "p_water_hammer": water_hammer_pressure(imp),
    }


def shear_sharing_factor(imp: Impact) -> dict:
    """How the WALL SHEAR STRESS changes on a non-rigid target.

    THIS IS NOT THE PRESSURE LAW, AND THE DIFFERENCE IS LARGE.
    Wall shear comes from SLIP between liquid and surface. On a rigid wall the
    slip is the full outer radial velocity. On a liquid target both sides move,
    with tangential velocities u_{r,i} = (2/pi) u_i r / sqrt(d^2 - r^2) and
    u_1 = U rho2/(rho1+rho2), u_2 = U rho1/(rho1+rho2). The slip driving the
    Stokes layer is their DIFFERENCE:

        |u_1 - u_2| / U = |rho2 - rho1| / (rho1 + rho2)                  (18)

    CONSEQUENCE, and it is the sharpest prediction in this module:

        a drop impacting a deep bath OF ITSELF feels NO wall shear at all.

    At matched density the interface is a symmetry plane carrying a slip
    boundary condition (CSG22's Fig. 7a caption says exactly that), so the
    Stokes layer never forms. The pressure halves; the shear goes to zero.
    Anything that depends on wall shear -- including the dilatancy gauge in
    dilatancy.py -- is switched off, not merely reduced.

    Returned values are the SLIP ratio. The true tau ratio carries a further
    two-sided boundary-layer factor <= 1 that is not computed here, so treat
    these as upper bounds on tau/tau_rigid.
    """
    rho1, rho2 = imp.drop.rho, imp.target.rho
    slip = 1.0 if np.isinf(rho2) else abs(rho2 - rho1) / (rho1 + rho2)
    return {
        "slip_ratio": slip,
        "tau_over_tau_rigid_upper_bound": slip,
        "note": ("exactly zero for a bath of the same liquid: the interface is a "
                 "symmetry plane and no Stokes layer forms"),
    }
