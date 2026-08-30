"""Push the calibration into the drop impact, and ask which half of it matters.

The point of calibrating this closure was never the flow curve. It was to put
numbers into a drop-impact solver. So the calibration is only finished when it
has been asked the question the solver asks.

The question
------------
Wyart-Cates has two halves:

- a **stress-activation** half, ``f(sigma) = exp(-(sigma*/sigma)^beta)``, which
  says *when* frictional contacts form. Its parameter ``sigma*`` is the one
  quantity in the model that a single flow curve measures cleanly: it comes from
  the curvature of the curve, it is untouched by the packing gauge, and it
  checks out to 0.73x and 1.83x against an independent shear-reversal
  measurement (``dst_fit`` gate G10).
- a **packing** half, ``phi_J(f) = (1-f) phi_0 + f phi_m``, which says *what
  happens* once they have. Its parameter ``phi_m`` is the one a flow curve does
  not measure: one material's own curves return a range 0.023 to 0.060 wide,
  against a per-curve error near 0.001.

In a rheometer both halves are in play, because the sweep straddles ``sigma*``.
In a drop impact they are not. The impact stress scale is

    Sigma = rho k U0^2

with ``k`` the added-mass coefficient of the jamming front, measured at 2 to 12
(Waitukaitis & Jaeger 2012; Peters, Majumdar & Jaeger 2016; Han et al. 2018).
At the drop-impact rig's own onset velocity that is kilopascals, one to three
decades above every ``sigma*`` fitted here. So ``f`` is pinned at 1 for the
whole event, the activation half is saturated, and everything the prediction
depends on sits in the half that was not measured.

That is the finding this module computes rather than asserts.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from wyart_cates import frictional_fraction_ss, relative_viscosity  # noqa: E402

# Added-mass coefficient of an impact-driven jamming front. Measured, not
# assumed: k = 2 to 12 across Waitukaitis & Jaeger, Nature 487, 205 (2012),
# Peters, Majumdar & Jaeger, Nature 532, 214 (2016), and Han, Peters & Jaeger,
# Nat. Commun. 7, 12243 (2016) / PRE 2018. Sigma = rho k U0^2.
K_ADDED_MASS = (2.0, 12.0)

# Suspension density at phi = 0.44 silica (rho_p = 2000 kg/m^3 for Stober
# silica) in water: rho = phi rho_p + (1-phi) rho_w.
RHO_SUSPENSION = 0.44 * 2000.0 + 0.56 * 1000.0  # 1440 kg/m^3

# The drop-impact rig's measured jamming onset, from the programme's own data:
# u0* ~ 1.15 m/s, quoted three ways at four salt levels.
U0_ONSET = 1.15  # m/s

# Volume-fraction tolerance the retake campaign holds itself to. Used only as
# the yardstick the calibrated phi_m band is compared against.
PHI_TOLERANCE = 0.002


# The deviatoric gauge. rho k U0^2 is a TOTAL stress, and under undrained
# loading -- which this silica is, by five to six decades -- the pore fluid
# carries the isotropic part, so only the deviatoric part generates particle
# pressure. impact_stress/STRESS_AND_DILATANCY_BRIEF.md establishes this for
# this geometry and measures the ratio: the area-averaged Wagner wall shear
# stress is 5 to 10x below rho U0^2. Both gauges are carried here, because
# which one you use changes whether the activation half is saturated, and the
# earlier version of this module used only the one the lab's own brief refutes.
DEVIATORIC_OVER_INERTIAL = (1.0 / 10.0, 1.0 / 5.0)


def impact_stress(u0=U0_ONSET, rho=RHO_SUSPENSION, k=K_ADDED_MASS):
    """Stress scale of an impact, ``rho k U0^2``, as a (low, high) pair over k.

    This is the INERTIAL (total) gauge. For the stress that actually reaches
    the particle network see :func:`deviatoric_stress`.
    """
    return tuple(float(rho * ki * u0**2) for ki in k)


def deviatoric_stress(u0=U0_ONSET, rho=RHO_SUSPENSION):
    """Area-averaged wall shear stress, as a (low, high) pair [Pa].

    Taken as ``rho U0^2`` times the 1/10 to 1/5 band the lab's own Wagner
    analysis measures, rather than recomputing the Wagner integral here: the
    band is what that analysis reports and the conclusions below are quoted
    across its whole width.
    """
    return tuple(float(rho * u0**2 * c) for c in DEVIATORIC_OVER_INERTIAL)


@dataclass(frozen=True)
class SaturationCheck:
    """Is the stress-activation half of the model saturated during an impact?

    "Saturated" is judged against the calibration's own resolution, not against
    an arbitrary closeness to 1. What the residual activation ``1 - f`` buys is a
    shift in the jamming point,

        phi_J(f) - phi_m = (1 - f) (phi_0 - phi_m),

    and that only matters if it is comparable to the width of the ``phi_m`` band
    the calibration actually returns. Comparing a threshold with the smallest
    change the measurement can express is the check this programme has had to
    add four separate times; it is built in here.
    """

    sigma_star: float
    sigma_impact_lo: float
    sigma_impact_hi: float
    f_lo: float
    f_hi: float
    beta: float
    phi_0: float
    phi_m: float
    band_width: float

    @property
    def decades_above(self) -> float:
        return float(np.log10(self.sigma_impact_lo / self.sigma_star))

    @property
    def phi_j_residual(self) -> float:
        """How far above phi_m the jamming point still sits at the weakest impact."""
        return float((1.0 - self.f_lo) * (self.phi_0 - self.phi_m))

    @property
    def saturated(self) -> bool:
        """True when the residual activation moves phi_J by under a tenth of the band."""
        if self.band_width <= 0:
            raise ValueError("band_width must be positive to judge saturation")
        return bool(self.phi_j_residual < 0.1 * self.band_width)

    def __str__(self) -> str:
        return (f"sigma* = {self.sigma_star:8.3g} Pa  ({self.decades_above:4.1f} dec below "
                f"impact)  f = {self.f_lo:.4f}-{self.f_hi:.4f}  "
                f"phi_J - phi_m = {self.phi_j_residual:.5f} vs band {self.band_width:.4f} "
                f"({self.phi_j_residual / self.band_width:6.1%})  "
                f"{'SATURATED' if self.saturated else 'NOT saturated'}")


def saturation(sigma_star, beta, phi_0, phi_m, band_width,
               u0=U0_ONSET, rho=RHO_SUSPENSION, k=K_ADDED_MASS, gauge="inertial"):
    """How close to 1 the frictional fraction sits at impact stresses, judged
    against the width of the phi_m band the same calibration returned.

    ``gauge`` selects ``"inertial"`` (rho k U0^2) or ``"deviatoric"`` (the wall
    shear). The verdict differs between them, so the gauge is named rather than
    assumed.
    """
    if gauge == "deviatoric":
        lo, hi = deviatoric_stress(u0, rho)
    elif gauge == "inertial":
        lo, hi = impact_stress(u0, rho, k)
    else:
        raise ValueError(f"gauge must be 'inertial' or 'deviatoric', got {gauge!r}")
    f_lo = float(frictional_fraction_ss(np.array([lo]), sigma_star, beta)[0])
    f_hi = float(frictional_fraction_ss(np.array([hi]), sigma_star, beta)[0])
    return SaturationCheck(float(sigma_star), lo, hi, f_lo, f_hi, float(beta),
                           float(phi_0), float(phi_m), float(band_width))


def viscosity_at_impact(phi, phi_m, phi_0, sigma_star, beta, n=2.0,
                        u0=U0_ONSET, rho=RHO_SUSPENSION, k=K_ADDED_MASS):
    """Relative viscosity the solver would use during the impact, (low, high) over k.

    ``inf`` means the closure has jammed: ``phi >= phi_J(f)``.

    ``n`` must be the divergence exponent the supplied ``phi_m`` was FITTED at.
    A phi_m obtained at n = 2.2 evaluated at n = 2 understates eta_r by about
    2.4x at the gaps these calibrations sit in; Guy et al.'s PMMA is n = 2.2 and
    Richards' corn starch is n = 2.
    """
    lo, hi = impact_stress(u0, rho, k)
    out = []
    for s in (lo, hi):
        f = frictional_fraction_ss(np.array([s]), sigma_star, beta)
        out.append(float(relative_viscosity(phi, f, phi_0, phi_m, n=n)[0]))
    return tuple(out)


@dataclass(frozen=True)
class BandConsequence:
    """What one material's phi_m band does to the impact prediction at a given phi."""

    phi: float
    phi_m_lo: float
    phi_m_hi: float
    eta_r_at_lo: float   # phi_m at the bottom of the band
    eta_r_at_hi: float   # phi_m at the top
    jams_at_lo: bool
    jams_at_hi: bool

    @property
    def verdict_flips(self) -> bool:
        return self.jams_at_lo != self.jams_at_hi

    @property
    def eta_ratio(self) -> float:
        """Ratio of the two predicted impact viscosities; inf if one has jammed."""
        if not np.isfinite(self.eta_r_at_lo) or not np.isfinite(self.eta_r_at_hi):
            return float("inf")
        lo, hi = sorted((self.eta_r_at_lo, self.eta_r_at_hi))
        return float(hi / lo) if lo > 0 else float("inf")


def band_consequence(phi, phi_m_lo, phi_m_hi, phi_0, sigma_star, beta, n=2.0, **kw):
    """Evaluate the impact prediction at both ends of a calibrated phi_m band.

    The two ends are two calibrations of the *same material from the same paper*,
    differing only in which of its own flow curves was used. If the verdict flips
    between them, no solver run using either is defensible on its own.

    ``n`` is the divergence exponent the band was fitted at; see
    :func:`viscosity_at_impact`.
    """
    e_lo = viscosity_at_impact(phi, phi_m_lo, phi_0, sigma_star, beta, n=n, **kw)[0]
    e_hi = viscosity_at_impact(phi, phi_m_hi, phi_0, sigma_star, beta, n=n, **kw)[0]
    # Judge jamming at the same f the viscosity was evaluated at, not at f -> 1.
    # The earlier version mixed the two and printed a finite eta_r labelled
    # "jammed" wherever f fell short of 1 (Guy's 404 nm, where sigma* is 514 Pa
    # and f is only 0.87 at the impact stress).
    return BandConsequence(
        phi=float(phi), phi_m_lo=float(phi_m_lo), phi_m_hi=float(phi_m_hi),
        eta_r_at_lo=e_lo, eta_r_at_hi=e_hi,
        jams_at_lo=not np.isfinite(e_lo),
        jams_at_hi=not np.isfinite(e_hi),
    )


def flip_window(phi_m_lo, phi_m_hi):
    """Range of packings over which the band's two ends disagree about jamming.

    Exactly ``(phi_m_lo, phi_m_hi)``: below the lower end both say flowing, above
    the upper end both say jammed, and in between they disagree. Returned as a
    function so the statement is computed rather than asserted, and so the width
    can be compared with the tolerance an experiment actually holds phi to.
    """
    lo, hi = sorted((float(phi_m_lo), float(phi_m_hi)))
    return lo, hi, hi - lo, (hi - lo) / PHI_TOLERANCE
