"""Is the lowest-stress end of a flow curve the frictionless branch, or Brownian?

Why this module exists
----------------------
Wyart-Cates has no Brownian physics in it. Its ``f = 0`` limit is the
*hard-sphere, high-Peclet* viscosity ``eta_r = (1 - phi/phi_0)^-n``, so the
low-stress plateau of a flow curve is only that branch when the suspension is
already out of the Brownian regime at the lowest stress sampled.

That is not a detail. Reading the volume fraction off a low-stress viscosity is
the only fit-free handle on ``phi`` a flow curve offers, and it is the handle the
tier-0 calibration needs, because a drop-impact solver advects the *real* volume
fraction and cannot use a fitted one. On Guy et al.'s d = 3770 nm PMMA that
handle works. On their d = 404 nm PMMA -- same chemistry, same lab, same figure,
9.3x smaller -- it is wrong by a factor of 260 in viscosity, and the Peclet
number at the lowest sampled stress says why: 2e-4.

So the admission rule below is not a convenience. It is the difference between a
calibration and a number.

Definition
----------
The shear Peclet number for a sphere of radius ``a`` in a solvent of viscosity
``eta_s`` at temperature ``T``, sheared at rate ``gamma_dot``:

    Pe = 6 pi eta_s a^3 gamma_dot / (k_B T)

which is the ratio of the time a particle takes to diffuse its own radius,
``a^2/D_0`` with ``D_0 = k_B T / (6 pi eta_s a)``, to the time the flow takes to
strain by one, ``1/gamma_dot``. ``Pe >> 1`` is the non-Brownian limit Wyart-Cates
is written in.

The suspension's own Peclet number is larger than the single-particle one by the
crowding factor eta_r, but the single-particle value is the conservative choice:
it is the one that does not depend on the model being calibrated.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# CODATA 2018, exact by definition of the kelvin.
K_B = 1.380649e-23  # J/K

# Guy, Hermes & Poon, PRL 115, 088304 (2015) state their solvent viscosities
# "at 19 degC". Nothing else in the deposit records a temperature.
T_GUY = 292.15  # K

# The admission threshold. Pe = 1 is where the two timescales cross, so it is a
# definition rather than a tuned cut. It is left at exactly 1 deliberately: the
# two Guy series sit about four decades apart in Pe (1e-4 against 1e0 to 1e2),
# so every threshold between 1e-3 and 1e2 returns the same verdict, and the
# calibration reports the margin rather than resting on the number.
PE_NON_BROWNIAN = 1.0


def peclet(shear_rate, radius, eta_s, temperature=T_GUY):
    """Single-particle shear Peclet number, ``6 pi eta_s a^3 gamma_dot / (k_B T)``.

    All arguments SI: rate in 1/s, radius in m, viscosity in Pa.s, temperature
    in K. Vectorised over ``shear_rate``.
    """
    shear_rate = np.asarray(shear_rate, dtype=float)
    if radius <= 0 or eta_s <= 0 or temperature <= 0:
        raise ValueError("radius, eta_s and temperature must be positive")
    return 6.0 * np.pi * eta_s * radius**3 * shear_rate / (K_B * temperature)


@dataclass(frozen=True)
class PlateauVerdict:
    """Whether one flow curve's lowest-stress point is on the frictionless branch."""

    stress: float          # Pa, the lowest stress sampled
    shear_rate: float      # 1/s there
    eta_r: float           # relative viscosity there
    pe: float              # single-particle Peclet number there
    admitted: bool
    reason: str

    def __str__(self) -> str:
        return (
            f"sigma_lo = {self.stress:.4g} Pa, gamma_dot = {self.shear_rate:.4g} 1/s, "
            f"eta_r = {self.eta_r:.4g}, Pe = {self.pe:.3g} -> "
            f"{'admitted' if self.admitted else 'REJECTED'} ({self.reason})"
        )


def plateau_admissible(
    stress,
    shear_rate,
    viscosity,
    *,
    radius,
    eta_s,
    temperature=T_GUY,
    pe_min=PE_NON_BROWNIAN,
) -> PlateauVerdict:
    """Judge the lowest-stress point of one flow curve.

    Admitted only if the suspension is non-Brownian there (``Pe >= pe_min``) and
    the relative viscosity is above 1, i.e. there is a suspension effect at all.
    The verdict carries the numbers so a rejection can be reported rather than
    silently dropping a curve.
    """
    stress = np.asarray(stress, dtype=float)
    shear_rate = np.asarray(shear_rate, dtype=float)
    viscosity = np.asarray(viscosity, dtype=float)
    if not (stress.shape == shear_rate.shape == viscosity.shape):
        raise ValueError("stress, shear_rate and viscosity must have the same shape")
    i = int(np.argmin(stress))
    pe = float(peclet(shear_rate[i], radius, eta_s, temperature))
    eta_r = float(viscosity[i] / eta_s)
    if not np.isfinite(eta_r) or eta_r <= 1.0:
        return PlateauVerdict(float(stress[i]), float(shear_rate[i]), eta_r, pe,
                              False, "relative viscosity <= 1, no suspension effect")
    if pe < pe_min:
        return PlateauVerdict(float(stress[i]), float(shear_rate[i]), eta_r, pe,
                              False, f"Brownian, Pe < {pe_min:g}")
    return PlateauVerdict(float(stress[i]), float(shear_rate[i]), eta_r, pe,
                          True, f"non-Brownian, Pe >= {pe_min:g}")


def invert_plateau(eta_r, phi_0, n):
    """Volume fraction implied by a frictionless-branch relative viscosity.

    Inverts ``eta_r = (1 - phi/phi_0)^-n`` to ``phi = phi_0 (1 - eta_r^(-1/n))``.
    Fit-free: no optimiser, no starting guess, no branch cut. This is the only
    place in the calibration where ``phi`` is obtained without an optimiser, and
    it is what the fitted ``phi`` is checked against.

    Returns NaN where ``eta_r <= 1``, which is outside the model's range rather
    than an error in it.
    """
    eta_r = np.asarray(eta_r, dtype=float)
    out = np.full(eta_r.shape, np.nan, dtype=float)
    ok = eta_r > 1.0
    out[ok] = phi_0 * (1.0 - eta_r[ok] ** (-1.0 / n))
    return out if out.shape else float(out)
