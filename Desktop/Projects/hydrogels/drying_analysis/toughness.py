"""
Fracture toughness Gamma (pure-shear / planar tension, dead-load) and the
fractocohesive length l_fc (notch-series knee).

--- Pure-shear toughness (Rivlin-Thomas) ---
In the pure-shear (planar-tension) geometry, the fracture energy is

    Gamma = W(lambda_c) * H0,

where H0 is the unstretched gauge height and W(lambda_c) is the strain-energy
density of the *unnotched* sheet at the critical stretch lambda_c at which the
notched sheet's crack propagates. W is obtained by integrating the measured
nominal-stress vs stretch curve:  W = integral_1^{lambda_c} sigma_nominal dlambda.
A dead-load ladder gives sigma at each hung weight; the camera gives lambda.

--- Fractocohesive length ---
From a notch-length series, the critical stretch (or critical nominal stress) is
flat for short notches and falls for long notches; the knee is at c ~ l_fc
(Chen-Wang-Suo; Yang-Yin-Suo). We locate the knee with a continuous two-segment
(plateau-then-line) fit and return the breakpoint.

References: Rivlin & Thomas, J. Polym. Sci. 10, 291 (1953); Chen, Wang & Suo,
Extreme Mech. Lett. 10, 50 (2017).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares


# --------------------------------------------------------------------------
# Strain-energy density and pure-shear toughness
# --------------------------------------------------------------------------
def strain_energy_density(stretch, sigma_nominal) -> float:
    """W = integral of nominal stress d(stretch), from the first sample to the last.

    `stretch` should start at ~1 (unstretched) and increase to lambda_c; sigma is
    the nominal (engineering) stress. Trapezoidal integration.
    """
    lam = np.asarray(stretch, dtype=float)
    sig = np.asarray(sigma_nominal, dtype=float)
    if lam.shape != sig.shape:
        raise ValueError("stretch and sigma_nominal must have equal shape")
    if lam.size < 2:
        raise ValueError("need >=2 points to integrate")
    if np.any(np.diff(lam) < 0):
        raise ValueError(
            "stretch must be monotonically non-decreasing (a loading ramp); got "
            "a decreasing segment -- pass only the loading branch, not unloading "
            "or a mis-ordered trace.")
    if not np.isclose(lam[0], 1.0, atol=0.05):
        import warnings
        warnings.warn(
            f"strain-energy integral starts at lambda={lam[0]:.3f}, not ~1.0; "
            "the segment from 1 to lambda[0] is omitted from W.", stacklevel=2)
    return float(np.trapezoid(sig, lam))


def toughness_pure_shear(W_c: float, H0: float) -> float:
    """Gamma = W(lambda_c) * H0. W_c in J/m^3 (=Pa), H0 in m -> Gamma in J/m^2."""
    if H0 <= 0:
        raise ValueError("H0 must be positive")
    if W_c < 0:
        raise ValueError("W_c must be non-negative")
    return float(W_c * H0)


def neo_hookean_W_pure_shear(G: float, lam: float) -> float:
    """Analytic strain-energy density for an incompressible neo-Hookean solid in
    pure shear (lambda1=lam, lambda2=1, lambda3=1/lam):
        W = (G/2)(lam^2 + 1/lam^2 - 2).  (Used for validation and quick estimates.)
    """
    return float(0.5 * G * (lam**2 + 1.0 / lam**2 - 2.0))


# --------------------------------------------------------------------------
# Fractocohesive length: knee of critical-quantity vs notch-length
# --------------------------------------------------------------------------
@dataclass
class KneeFit:
    l_fc: float               # breakpoint (same units as notch length c)
    plateau_value: float      # critical quantity on the short-notch plateau
    slope_after: float        # slope of the falling branch
    rms_residual: float


def _piecewise_plateau_line(params, c):
    """y = plateau for c <= b; plateau + slope*(c-b) for c > b."""
    plateau, slope, b = params
    return np.where(c <= b, plateau, plateau + slope * (c - b))


def fractocohesive_knee(notch_length, critical_quantity) -> KneeFit:
    """Locate l_fc as the breakpoint of a plateau-then-line fit to
    critical_quantity(notch_length). `critical_quantity` is lambda_c or the
    critical nominal stress; both flatten below l_fc and fall above it.
    """
    c = np.asarray(notch_length, dtype=float)
    y = np.asarray(critical_quantity, dtype=float)
    if c.size < 4:
        raise ValueError(f"need >=4 notch lengths to fit a knee; got {c.size}")
    order = np.argsort(c)
    c, y = c[order], y[order]

    plateau0 = float(np.median(y[: max(2, c.size // 3)]))
    slope0 = float((y[-1] - y[0]) / (c[-1] - c[0]))
    b0 = float(np.median(c))

    def resid(p):
        return _piecewise_plateau_line(p, c) - y

    lo = [-np.inf, -np.inf, c.min()]
    hi = [np.inf, np.inf, c.max()]
    res = least_squares(resid, [plateau0, slope0, b0], bounds=(lo, hi),
                        max_nfev=2000)
    plateau, slope, b = res.x
    rms = float(np.sqrt(np.mean(res.fun**2)))
    return KneeFit(l_fc=float(b), plateau_value=float(plateau),
                   slope_after=float(slope), rms_residual=rms)


# --------------------------------------------------------------------------
# Network length ell = Gamma / E, with the neo-Hookean l_fc cross-check
# --------------------------------------------------------------------------
def network_length(Gamma: float, E: float) -> float:
    """ell = Gamma / E (units of length)."""
    if E <= 0:
        raise ValueError("E must be positive")
    return float(Gamma / E)


def ell_from_lfc(l_fc: float, lam_c: float) -> float:
    """Cross-check: ell = l_fc * f(lam_c) with the PURE-SHEAR neo-Hookean
    f = W/E = (lam_c^2 + 1/lam_c^2 - 2)/6 (consistent with the pure-shear notch
    series that produces lam_c and with `neo_hookean_W_pure_shear`, since E=3G).
    For brittle low-stretch gels (lam_c ~ 1.1-1.5), f -> ~1 so ell ~ l_fc.
    """
    f = (lam_c**2 + 1.0 / lam_c**2 - 2.0) / 6.0
    return float(l_fc * f)
