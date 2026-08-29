"""
Elastic modulus E from gravity-sphere indentation (camera-measured depth).

A dense sphere of known buoyant weight rests on the gel; the camera measures the
indentation depth delta. Hertzian contact for a rigid sphere on an incompressible
(nu=0.5) half-space gives

    F = (16/9) E sqrt(R) delta^{3/2}   ->   E = 9 F / (16 sqrt(R) delta^{3/2}),

where F = (rho_sphere - rho_fluid) V g is the buoyant weight.

Finite-thickness correction (the artifact is up to ~600% if ignored; MOC section 5):
for a bonded incompressible layer of thickness h, Dimitriadis et al. (Biophys. J.
82, 2798, 2002) give a polynomial stiffening factor in chi = sqrt(R delta)/h:

    F = (16/9) E sqrt(R) delta^{3/2} * P(chi),
    P(chi) = 1 + 1.133 chi + 1.497 chi^2 + 1.469 chi^3 + 0.755 chi^4   (bonded).

Because chi is computable directly from measured (R, delta, h), inverting for E
needs no iteration. Keep delta/h ~ 0.1-0.4 and read at swelling equilibrium (this
is the *drained* modulus, ~10-40% below a rheometer G'; reconcile when cross-checking).

Reference: Dimitriadis, Horkay, Maresca, Kachar & Chadwick, Biophys. J. 82, 2798 (2002).
"""
from __future__ import annotations

import numpy as np

G_ACCEL = 9.80665  # m / s^2

# Dimitriadis bonded-incompressible finite-thickness polynomial coefficients.
_DIMITRIADIS_BONDED = (1.0, 1.133, 1.497, 1.469, 0.755)


def buoyant_weight(radius_m: float, rho_sphere: float, rho_fluid: float) -> float:
    """F = (rho_sphere - rho_fluid) * V * g for a sphere of given radius."""
    if radius_m <= 0:
        raise ValueError("radius must be positive")
    V = 4.0 / 3.0 * np.pi * radius_m**3
    return (rho_sphere - rho_fluid) * V * G_ACCEL


def finite_thickness_factor(R, delta, h) -> float:
    """Dimitriadis bonded-incompressible stiffening factor P(chi), chi=sqrt(R delta)/h."""
    if h <= 0:
        raise ValueError("h must be positive")
    chi = np.sqrt(R * delta) / h
    c0, c1, c2, c3, c4 = _DIMITRIADIS_BONDED
    return c0 + c1 * chi + c2 * chi**2 + c3 * chi**3 + c4 * chi**4


def modulus_from_indentation(F, R, delta, h=None) -> float:
    """Young's modulus E (Pa) from indentation force F, sphere radius R, depth delta.

    If layer thickness h is given, applies the Dimitriadis finite-thickness
    correction; otherwise returns the half-space (infinite-thickness) Hertz value.
    All lengths in metres, F in newtons.
    """
    if not (F > 0 and R > 0 and delta > 0):
        raise ValueError(f"F,R,delta must be positive; got F={F},R={R},delta={delta}")
    E_inf = 9.0 * F / (16.0 * np.sqrt(R) * delta**1.5)
    if h is None:
        return float(E_inf)
    return float(E_inf / finite_thickness_factor(R, delta, h))


def forward_indentation_force(E, R, delta, h=None) -> float:
    """Forward model: predicted force for a given E (used for validation)."""
    F_inf = (16.0 / 9.0) * E * np.sqrt(R) * delta**1.5
    if h is None:
        return float(F_inf)
    return float(F_inf * finite_thickness_factor(R, delta, h))
