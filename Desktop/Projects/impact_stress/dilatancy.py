"""
dilatancy.py — does an impact push a dense suspension past its dilatancy onset?

THE QUESTION
------------
A dense suspension shear-thickens when the local stress exceeds an onset stress
sigma*, at which the repulsive force holding particles apart is overcome and
contacts turn frictional (Wyart & Cates 2014). Drop impact imposes stress. So:

    does the impact cross sigma*, where, when, and over how much of the drop?

A dimensionless group answering that is what "an inertio-dilatant number" has
to be. This module builds three candidates, shows that two of them cannot work,
keeps the third, and evaluates it for the lab's 800 nm silica.

THE THREE CANDIDATES
--------------------
    Lambda_P  = rho U^2 / sigma*                     (inertial / pressure gauge)
    Lambda_F  = <tau_wall> / sigma*                  (shear gauge, Eq. 14')
    Phi_f(t)  = fraction of the contact area with tau > sigma*   (field version)

Lambda_P is the obvious one and it is the one the lab's own analysis used. It
DOES NOT WORK, and the project already established why (referee_objections_
reframe.md): rho U^2 = 1.5-37 kPa over the accessible velocity range while
every well-characterised shear-thickening suspension has sigma* = O(1-100) Pa.
Every drop of every suspension sits 15x to 37000x above sigma*, at the lowest
velocity anyone can drop from. A gauge with no reachable threshold cannot
measure a threshold.

Lambda_F fixes exactly that:

    <tau> / (rho U^2) = 1.135 Re_R^(-1/2) t*^(-1/2)

so the shear gauge is smaller than the inertial gauge by ~1/sqrt(Re). At the
lab's Re_R ~ 60-200 that is a factor 8-14, which moves a 40-800x overshoot into
a range where a threshold is crossed near the bottom of the experiment's
velocity window.

WHY THIS IS PHYSICS AND NOT JUST A SMALLER NUMBER -- CORRECTED AFTER REVIEW.
An earlier version of this module argued "dilatancy is a shear phenomenon;
sigma* is a shear stress; rho U^2 is an isotropic pressure". That is BACKWARDS
relative to Wyart & Cates, who write their model in PARTICLE PRESSURE: their
P* ~ F*/R^2, their control parameter is p = P/P*, and their Fig. 1 caption says
"onset pressure p ~ 1". A flow-curve fit cannot tell the two apart, because
sigma/P = mu(I_v) = O(1) in steady simple shear, so the repo's sigma* = 44.2 Pa
is really mu * P* and P* is 44-150 Pa.

The argument that survives is about TOTAL versus EFFECTIVE stress. rho U^2 is a
total stress applied to a saturated two-phase medium on a timescale far shorter
than pore-pressure equilibration. Under undrained loading the pore fluid carries
essentially all of an isotropic increment (Skempton B -> 1), so the EFFECTIVE
particle pressure -- the quantity that actually sets p = P/P* -- barely moves.
What does generate effective particle pressure is shear, through Reynolds
dilatancy: the packing must expand to shear, the expansion is resisted by Darcy
backflow, and the resulting suction presses the particles together.

So the shear stress is the right gauge not because pressure is irrelevant, but
because in an undrained impact the deviatoric part is what converts applied load
into particle-scale normal force. Lambda_F is a proxy for a pressure criterion,
not a rival to one. (This is also the mechanism `inertial_dilatant.py` uses,
which is why the two modules must not be allowed to contradict each other.)

Phi_f is the version with spatial structure, and it is the one that connects to
what the lab actually measures — a solidification FRONT, i.e. a moving boundary
between thickened and unthickened material. In this framework that front is the
tau = sigma* contour.

SOURCES
-------
[WC14]   Wyart & Cates, Phys. Rev. Lett. 112, 098302 (2014).
[PLJ16]  Philippi, Lagree & Antkowiak, J. Fluid Mech. 795, 96 (2016).
[CSG22]  Cheng, Sun & Gordillo, Annu. Rev. Fluid Mech. 54, 57 (2022).
[SAD22]  Shah, Arora & Driscoll, Commun. Phys. 5, 222 (2022),
         doi:10.1038/s42005-022-00998-w.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

import stress as st
from materials import WyartCates
from stress import DELTA_OVER_D, Impact

T_STAR_EVAL = 0.315
"""Default dimensionless time at which the gauges are evaluated.

WHY NOT t* = 1.  The Wagner self-similar solution is derived for short
times and PLJ16 report departure from similarity once t* = U t/R is of
order 1. At t* = 1 the contact radius is d = sqrt(3) R = 1.73 R, larger
than the drop itself, so the solution is being used outside its domain.

WHY 0.315.  That is the PEAK-FORCE time for the lab condition (CSG22's
ttil_max = 0.157 on the diameter clock, doubled), i.e. the instant of
maximum loading, and there d = 0.97 R -- just inside validity.

The choice matters: moving from t* = 1 to t* = 0.315 raises every stress
by 1/sqrt(0.315) = 1.78 and lowers every onset velocity by 0.315^(1/3) =
0.68. That factor is part of the honest uncertainty and is reported as a
sensitivity rather than buried in a default."""

# --------------------------------------------------------------------------
# Published onset stresses for ~800 nm charge-stabilised silica
# --------------------------------------------------------------------------
SIGMA_STAR_LITERATURE = {
    "Ishibashi (2026), fitted law with measured F* = 0.3 nN": 19.1,
    "Lin et al., PRL 115, 228304 (2015), scaled by a^-2": 62.5,
    "Shah, Arora & Driscoll (2022), the group's own 830 nm particles": 100.0,
    "repo flow-curve fit, params_fitted.yaml (provisional)": 44.2327,
}
"""Independent estimates of sigma* for this material, Pa.

The first three are external to this project and do not depend on the
retracted rheology file. Their geometric mean is 49.2 Pa; the repo's own
flow-curve fit (44.2 Pa) sits inside their range, which is the main reason to
trust it as an order-of-magnitude anchor despite its failed fit gate.
"""


def sigma_star_bracket() -> tuple[float, float, float]:
    """(low, geometric-mean, high) of the published external estimates [Pa]."""
    vals = [v for k, v in SIGMA_STAR_LITERATURE.items() if "repo" not in k]
    return min(vals), float(np.exp(np.mean(np.log(vals)))), max(vals)


def sigma_star_from_repulsion(F_star: float, a: float) -> float:
    """Onset stress from a measured repulsive force scale: sigma* = F*/(6 pi a^2).

    The standard scaling argument: contacts turn frictional when the shear
    stress supplies a force ~F* across a particle cross-section ~a^2. The
    6 pi is the usual Stokes-drag convention; different papers differ by O(1),
    so use this for SCALING WITH PARTICLE SIZE (sigma* ~ a^-2), not for
    absolute values.
    """
    return F_star / (6.0 * np.pi * a**2)


# --------------------------------------------------------------------------
# The three candidate numbers
# --------------------------------------------------------------------------
def lambda_inertial(imp: Impact, sigma_star: float) -> float:
    """Lambda_P = rho U^2 / sigma*.  The gauge that does not work — see above."""
    return imp.p_dyn / sigma_star


def lambda_shear(imp: Impact, sigma_star: float, t_star: float = 1.0) -> float:
    """Lambda_F = <tau_wall>(t*) / sigma*.  THE PROPOSED INERTIO-DILATANT NUMBER.

        Lambda_F(t*) = 1.135 (rho U^2 / sigma*) Re_R^(-1/2) t*^(-1/2)

    Written out in primitive variables (the form worth memorising):

        Lambda_F(t*) = 1.135 t*^(-1/2) sqrt(rho eta U^3 / R) / sigma*

    Note what it depends on and how weakly:
        U^(3/2)   strong — this is what makes an onset velocity sharp
        eta^(1/2) weak — a 5x uncertainty in eta_0 is only 2.2x in Lambda
        R^(-1/2)  weak
        sigma*^-1 direct

    Lambda_F = 1 defines an onset impact velocity

        U_c = [ sigma*^2 R / (1.29 rho eta) ]^(1/3)  * t*^(1/3)

    i.e. U_c ~ sigma*^(2/3): a suspension whose onset stress is lowered 8x by
    salt should start jamming at a 4x lower impact speed. That is a falsifiable
    prediction and it is the main experimental output of this analysis.
    """
    return float(st.wall_shear_mean(imp, t_star)) / sigma_star


def lambda_shear_peak(imp: Impact, sigma_star: float, t_star: float = 1.0) -> float:
    """Same, using the screened RIM stress instead of the area average.

    Larger by 10.90/1.135 = 9.6x. Quote as an upper bound: it depends on
    PLJ16's screening length, which is an empirical cut-off, whereas the area
    average is a convergent integral.
    """
    return float(st.wall_shear_peak(imp, t_star)) / sigma_star


# --------------------------------------------------------------------------
# The field version: where is the suspension frictional?
# --------------------------------------------------------------------------
def frictional_contour(imp: Impact, sigma_star: float, t_star):
    """Radius xi_c(t*) of the tau = sigma* contour, in units of the contact radius.

    From Eq. (13'),  tau = A t*^(-1/2) xi/(1-xi^2)  with
    A = (2/pi^(3/2)) rho U^2 Re_R^(-1/2), the contour tau = sigma* is

        xi/(1-xi^2) = K,   K = sigma* sqrt(t*) / A
        =>  xi_c = [ -1 + sqrt(1 + 4 K^2) ] / (2K)                        (17)

    Material at xi > xi_c is above the onset stress and turning frictional;
    material at xi < xi_c is still lubricated. So xi_c is an INNER boundary:
    the thickened region is an ANNULUS hugging the rim, not a disc at the
    centre. That is the qualitative prediction to check against the lab's
    under-view footage.

    xi_c -> 1 (nothing thickens) when K is large, i.e. weak impact or large
    sigma*; xi_c -> 0 (everything thickens) when K -> 0.
    """
    t_star = np.asarray(t_star, dtype=float)
    A = (2.0 / np.pi**1.5) * imp.p_dyn * imp.Re_R**-0.5
    K = sigma_star * np.sqrt(t_star) / A
    # algebraically (-1+sqrt(1+4K^2))/(2K), rearranged to avoid the
    # catastrophic cancellation of that form as K -> 0.
    return 2.0 * K / (1.0 + np.sqrt(1.0 + 4.0 * K**2))


def frictional_area_fraction(imp: Impact, sigma_star: float, t_star):
    """Fraction of the wetted contact area above sigma*.

        Phi_f(t*) = [ (1-Delta/d)^2 - xi_c^2 ] / (1-Delta/d)^2 ,  clipped to [0,1]

    The outer limit is PLJ16's screening cut-off rather than xi = 1, because
    the formula is not trusted inside it. That choice changes Phi_f by ~3%.

    Phi_f is the natural order parameter here: it is dimensionless, bounded,
    resolution-independent, and it is what an under-view image of a partly
    solidified drop actually shows.
    """
    xi_out = 1.0 - DELTA_OVER_D
    xi_c = frictional_contour(imp, sigma_star, t_star)
    frac = (xi_out**2 - np.clip(xi_c, 0.0, xi_out) ** 2) / xi_out**2
    return np.clip(frac, 0.0, 1.0)


# --------------------------------------------------------------------------
# Onset velocity
# --------------------------------------------------------------------------
def onset_velocity(susp: WyartCates, R: float, sigma_star: float | None = None,
                   t_star: float = 1.0, criterion: str = "mean",
                   area_threshold: float = 0.5, eta: float | None = None,
                   U_bracket: tuple[float, float] = (0.01, 50.0)) -> float:
    """Impact speed at which the chosen dilatancy criterion is first met [m/s].

    criterion:
      "mean"  : <tau>(t*) = sigma*                    (Lambda_F = 1)
      "peak"  : rim tau(t*) = sigma*                  (upper bound)
      "area"  : Phi_f(t*) = area_threshold            (the field criterion)

    eta defaults to the suspension's low-stress viscosity eta_0 — the right
    choice, because the question is whether the UNTHICKENED drop gets pushed
    over the onset.

    The three criteria differ by a factor of a few in U_c. That spread is the
    honest uncertainty of the method and is reported alongside every number
    rather than being hidden by picking one.
    """
    from materials import suspension_as_fluid

    sigma_star = susp.sigma_star if sigma_star is None else sigma_star
    base = suspension_as_fluid(susp)
    if eta is not None:
        base = type(base)(name=base.name, rho=base.rho, c=base.c, eta=eta,
                          gamma=base.gamma, provenance=base.provenance + " [eta overridden]")

    def f(U):
        imp = Impact(drop=base, U=U, R=R, target=st.rigid_target())
        if criterion == "mean":
            return lambda_shear(imp, sigma_star, t_star) - 1.0
        if criterion == "peak":
            return lambda_shear_peak(imp, sigma_star, t_star) - 1.0
        if criterion == "area":
            return float(frictional_area_fraction(imp, sigma_star, t_star)) - area_threshold
        raise ValueError(f"unknown criterion {criterion!r}")

    lo, hi = U_bracket
    if f(lo) * f(hi) > 0:
        raise ValueError(
            f"onset criterion {criterion!r} is not crossed anywhere in "
            f"U = [{lo}, {hi}] m/s for sigma* = {sigma_star:g} Pa. Returning "
            f"NaN here would put an unflagged NaN into the results file.")
    return float(brentq(f, lo, hi, xtol=1e-9, rtol=1e-12))


def onset_velocity_closed_form(sigma_star: float, R: float, rho: float,
                               eta: float, t_star: float = 1.0) -> float:
    """Closed form of the "mean" criterion — no root-finding needed.

        1.135 sqrt(rho eta U^3 / R) t*^(-1/2) = sigma*
        =>  U_c = [ sigma*^2 R t* / (1.2889 rho eta) ]^(1/3)

    with 1.2889 = (10.7/3pi)^2. Provided as an independent check on
    `onset_velocity(criterion="mean")`; the test suite asserts they agree.
    """
    coeff = st.MEAN_SHEAR_COEFF
    return (sigma_star**2 * R * t_star / (coeff**2 * rho * eta)) ** (1.0 / 3.0)


# --------------------------------------------------------------------------
# Reconciling the lab's two shear-rate conventions
# --------------------------------------------------------------------------
def tstar_matching_shear_rate(imp: Impact, gammadot_target: float) -> float:
    """Which t* makes the contact-averaged shear rate equal a given value.

        <gammadot>(t*) = 1.135 (rho U^2/eta) Re_R^(-1/2) t*^(-1/2)
        =>  t* = [ 1.135 rho U^2 / (eta Re_R^(1/2) gammadot) ]^2

    Used to show that SAD22's two conventions — gammadot_impact = U/D
    (~1e3 1/s) and gammadot_spr = U/l with l ~ 100 um (~1e4 1/s) — are the
    same field evaluated at two different instants, rather than two
    incompatible estimates.
    """
    amp = st.MEAN_SHEAR_COEFF * imp.p_dyn / imp.drop.eta * imp.Re_R**-0.5
    return (amp / gammadot_target) ** 2


# --------------------------------------------------------------------------
# Convenience: everything at once for one impact
# --------------------------------------------------------------------------
def dilatancy_report(imp: Impact, susp: WyartCates, sigma_star: float | None = None,
                     t_star: float = 1.0) -> dict:
    """All three candidate numbers plus the field diagnostics, for one impact."""
    s = susp.sigma_star if sigma_star is None else sigma_star
    return {
        "U": imp.U,
        "sigma_star_Pa": s,
        "rho_U2_Pa": imp.p_dyn,
        "Re_R": imp.Re_R,
        "Lambda_P": lambda_inertial(imp, s),
        "tau_mean_Pa": float(st.wall_shear_mean(imp, t_star)),
        "Lambda_F": lambda_shear(imp, s, t_star),
        "tau_peak_Pa": float(st.wall_shear_peak(imp, t_star)),
        "Lambda_F_peak": lambda_shear_peak(imp, s, t_star),
        "xi_contour": float(frictional_contour(imp, s, t_star)),
        "frictional_area_fraction": float(frictional_area_fraction(imp, s, t_star)),
        "gammadot_mean_per_s": float(st.shear_rate_mean(imp, t_star)),
        "gammadot_naive_per_s": st.naive_impact_shear_rate(imp),
        "t_star": t_star,
    }
