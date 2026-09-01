"""Layer 2: a shear-thickening element with a finite strain clock.

Why this module exists
----------------------
Every shear-thickening damper model in the device literature is quasi-steady: it
reads a viscosity off a measured steady flow curve at the instantaneous shear
rate. That is only legitimate when the microstructure reaches its steady state
much faster than the mechanical event the device is built to resist.

In a shear-thickening suspension the microstructure does not relax on a *time*,
it relaxes on a *strain*. The frictional contact fraction ``f`` follows

    d f / d gamma = (f_ss(sigma) - f) / gamma_0

with ``gamma_0`` of order one strain unit. During a gait event of duration
``tau`` at shear rate ``gamma_dot`` the material accumulates ``gamma_dot * tau``
of strain, so whether it thickens at all is set by

    A = gamma_dot * tau / gamma_0                     (:func:`activation_number`)

and the quasi-steady assumption is the limit ``A -> infinity``. Section 3 of
``SCALING_2026-08-31.md`` works out that ``A`` is order one across the published
design band, which is why this module exists rather than a lookup on a flow
curve.

Conventions, and one that is a trap
-----------------------------------
SI throughout. Stress-controlled Wyart-Cates, reusing ``dst_fit.models`` rather
than restating it: ``frictional_fraction``, ``jamming_fraction`` and
``divergence`` are imported, not copied, so a correction there propagates here.

The trap: Wyart-Cates is exactly covariant under
``(phi, phi_m, phi_0) -> c (phi, phi_m, phi_0)``, so ``phi_0`` is a gauge choice
and only the ratios ``phi/phi_0`` and ``phi_m/phi_0`` are ever measured. Do not
present ``phi_0`` as a design knob and do not sweep it expecting information;
sweeping it re-measures algebra. See the ``wyart-cates-phi0-is-gauge`` note.
:func:`reduced_packing` returns the two combinations that do carry meaning.

Jamming is reported, not clipped
--------------------------------
When ``phi`` reaches ``phi_J(f)`` the viscosity diverges. That is shear jamming,
a physical event, and :func:`transient_response` terminates and flags it rather
than returning a large finite number a fit could chase.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

from dst_fit.models import divergence, frictional_fraction, jamming_fraction

__all__ = [
    "STFMaterial",
    "activation_number",
    "reduced_packing",
    "reduced_rate",
    "steady_curve",
    "steady_branch_count",
    "steady_stress",
    "quasi_steady_ratio",
    "transient_response",
    "TransientResult",
]

# Integration stops when phi/phi_J reaches this, i.e. just short of the pole.
# Not a physical constant: it is the numerical stand-off from a divergence, and
# results must not depend on it. test_activation.py::test_jam_margin_is_inert
# checks that halving it moves nothing that is reported.
JAM_STANDOFF = 1.0 - 1e-6


@dataclass(frozen=True)
class STFMaterial:
    """A shear-thickening suspension, Wyart-Cates plus a strain clock.

    Parameters are SI. ``gamma_0`` is the strain clock: the strain over which
    the frictional fraction relaxes towards its steady value. Measured at 0.5 to
    1.6 on silica and PMMA from shear-reversal transients
    (``dst_fit/STRAIN_CLOCK_2026-08-28.md``), against the 0.2 that drop-impact
    solvers borrow from corn starch.

    ``phi_0`` is a gauge choice, not a measurement. See the module docstring.
    """

    phi: float
    phi_m: float
    sigma_star: float
    eta_s: float
    gamma_0: float
    phi_0: float = 0.64
    n: float = 2.0
    beta: float = 1.0

    def __post_init__(self) -> None:
        if not 0.0 < self.phi_m < self.phi_0:
            raise ValueError(f"need 0 < phi_m ({self.phi_m}) < phi_0 ({self.phi_0})")
        if not 0.0 < self.phi < self.phi_0:
            raise ValueError(f"need 0 < phi ({self.phi}) < phi_0 ({self.phi_0})")
        for name in ("sigma_star", "eta_s", "gamma_0", "n"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive, got {getattr(self, name)}")

    def viscosity(self, f: float | np.ndarray) -> np.ndarray:
        """Suspension viscosity at frictional fraction ``f``. inf when jammed."""
        phi_j = jamming_fraction(np.asarray(f, dtype=float), self.phi_0, self.phi_m)
        return divergence(self.phi, phi_j, eta_s=self.eta_s, n=self.n)

    def f_steady(self, sigma: float | np.ndarray) -> np.ndarray:
        """Steady frictional fraction at stress ``sigma``."""
        return frictional_fraction(sigma, self.sigma_star, beta=self.beta)

    def is_jammed(self, f: float) -> bool:
        """True when phi has reached the jamming fraction at this ``f``."""
        return self.phi / jamming_fraction(f, self.phi_0, self.phi_m) >= JAM_STANDOFF

    def f_jam(self) -> float:
        """Frictional fraction at which this phi jams, or nan if it never does.

        From phi = (1-f) phi_0 + f phi_m, solved for f. Values outside [0, 1]
        mean the packing cannot shear jam at any stress.
        """
        f = (self.phi_0 - self.phi) / (self.phi_0 - self.phi_m)
        return float(f) if 0.0 <= f <= 1.0 else float("nan")


def activation_number(gamma_dot: float | np.ndarray, tau_event: float, gamma_0: float) -> np.ndarray:
    """A = gamma_dot * tau_event / gamma_0, the strain available in units of the clock.

    A >> 1 is the quasi-steady limit the device literature assumes. A <~ 1 means
    the suspension does not finish thickening during the event.
    """
    if tau_event <= 0 or gamma_0 <= 0:
        raise ValueError("tau_event and gamma_0 must be positive")
    return np.asarray(gamma_dot, dtype=float) * tau_event / gamma_0


def reduced_packing(mat: STFMaterial) -> tuple[float, float]:
    """The two gauge-invariant packing combinations, ``(phi/phi_0, phi_m/phi_0)``.

    These, not the three fractions separately, are what a flow curve measures.
    """
    return mat.phi / mat.phi_0, mat.phi_m / mat.phi_0


def reduced_rate(gamma_dot: float | np.ndarray, mat: STFMaterial) -> np.ndarray:
    """Lambda = eta_s gamma_dot / sigma_star: solvent stress in units of the onset stress.

    The stress is ``sigma = sigma_star * Lambda * (eta/eta_s)``, so Lambda and the
    reduced viscosity carry the whole rate dependence.
    """
    return mat.eta_s * np.asarray(gamma_dot, dtype=float) / mat.sigma_star


def steady_curve(sigma: np.ndarray, mat: STFMaterial) -> tuple[np.ndarray, np.ndarray]:
    """Steady flow curve, parameterised by stress.

    Returns ``(gamma_dot, eta)``. Parameterised by stress rather than rate
    because for a discontinuously thickening suspension ``sigma(gamma_dot)`` is
    multivalued while ``gamma_dot(sigma)`` is not. This mirrors the choice in
    ``dst_fit.models`` and avoids a root find on an S-shaped branch.
    """
    sigma = np.asarray(sigma, dtype=float)
    eta = mat.viscosity(mat.f_steady(sigma))
    with np.errstate(divide="ignore", invalid="ignore"):
        gamma_dot = np.where(np.isfinite(eta), sigma / eta, 0.0)
    return gamma_dot, eta


def steady_stress(gamma_dot: float, mat: STFMaterial, *, n_grid: int = 20001) -> float:
    """Steady stress at an imposed shear rate, or inf if the rate is unreachable.

    Found by interpolating the parametric curve rather than root finding, so a
    multivalued region is detected instead of silently returning one branch. If
    the requested rate is bracketed more than once the curve is S-shaped there
    and the *lowest* stress solution is returned, with the ambiguity reported by
    :func:`steady_branch_count`.
    """
    sigma_grid = np.geomspace(1e-6 * mat.sigma_star, 1e6 * mat.sigma_star, n_grid)
    rates, _ = steady_curve(sigma_grid, mat)
    finite = np.isfinite(rates) & (rates > 0)
    if not finite.any():
        return float("inf")
    rates, sigma_grid = rates[finite], sigma_grid[finite]
    crossings = np.flatnonzero(np.diff(np.sign(rates - gamma_dot)) != 0)
    if crossings.size == 0:
        return float("inf")
    i = crossings[0]
    lo, hi = rates[i], rates[i + 1]
    w = 0.0 if hi == lo else (gamma_dot - lo) / (hi - lo)
    return float(np.exp(np.log(sigma_grid[i]) * (1 - w) + np.log(sigma_grid[i + 1]) * w))


def steady_branch_count(gamma_dot: float, mat: STFMaterial, *, n_grid: int = 20001) -> int:
    """How many steady stresses satisfy this rate. >1 means the curve is S-shaped here."""
    sigma_grid = np.geomspace(1e-6 * mat.sigma_star, 1e6 * mat.sigma_star, n_grid)
    rates, _ = steady_curve(sigma_grid, mat)
    finite = np.isfinite(rates) & (rates > 0)
    rates = rates[finite]
    return int(np.count_nonzero(np.diff(np.sign(rates - gamma_dot)) != 0))


@dataclass(frozen=True)
class TransientResult:
    """Trajectory of a strain-clocked element under an imposed shear-rate history."""

    t: np.ndarray
    f: np.ndarray
    sigma: np.ndarray
    eta: np.ndarray
    gamma_dot: np.ndarray
    strain: np.ndarray
    jammed: bool
    t_jam: float

    @property
    def sigma_peak(self) -> float:
        return float(np.nanmax(self.sigma))


def transient_response(
    gamma_dot_of_t: Callable[[float], float],
    mat: STFMaterial,
    t_end: float,
    *,
    f0: float = 0.0,
    n_out: int = 2001,
    rtol: float = 1e-8,
    atol: float = 1e-10,
) -> TransientResult:
    """Integrate the strain clock under an imposed shear-rate history.

    State is the frictional fraction. At each instant the stress is algebraic,
    ``sigma = eta(f) * gamma_dot``, so no implicit solve is needed: ``f`` is the
    only thing with memory.

        d f / d t = gamma_dot * (f_ss(sigma) - f) / gamma_0

    ``f0`` is the frictional fraction at rest. Zero is the fully lubricated
    start, appropriate for a device that has been still between steps.

    Terminates and sets ``jammed`` if phi reaches phi_J. The viscosity is not
    clipped, so a jam is reported rather than absorbed.
    """
    if t_end <= 0:
        raise ValueError("t_end must be positive")
    if not 0.0 <= f0 <= 1.0:
        raise ValueError(f"f0 must lie in [0, 1], got {f0}")

    def rhs(t: float, y: np.ndarray) -> list[float]:
        f = float(np.clip(y[0], 0.0, 1.0))
        gd = float(gamma_dot_of_t(t))
        eta = float(mat.viscosity(f))
        if not np.isfinite(eta):
            return [0.0]
        sigma = eta * abs(gd)
        f_ss = float(mat.f_steady(sigma))
        return [abs(gd) * (f_ss - f) / mat.gamma_0]

    def jam_event(t: float, y: np.ndarray) -> float:
        phi_j = jamming_fraction(float(np.clip(y[0], 0.0, 1.0)), mat.phi_0, mat.phi_m)
        return JAM_STANDOFF - mat.phi / phi_j

    jam_event.terminal = True  # type: ignore[attr-defined]
    jam_event.direction = -1.0  # type: ignore[attr-defined]

    t_eval = np.linspace(0.0, t_end, n_out)
    sol = solve_ivp(
        rhs, (0.0, t_end), [f0], t_eval=t_eval, events=jam_event,
        rtol=rtol, atol=atol, method="LSODA", max_step=t_end / 200.0,
    )
    if not sol.success and sol.status != 1:
        raise RuntimeError(f"strain-clock integration failed: {sol.message}")

    t = sol.t
    f = np.clip(sol.y[0], 0.0, 1.0)
    gd = np.array([float(gamma_dot_of_t(ti)) for ti in t])
    eta = mat.viscosity(f)
    sigma = eta * np.abs(gd)
    strain = np.concatenate([[0.0], np.cumsum(np.abs(gd[:-1]) * np.diff(t))]) if t.size > 1 else np.zeros(1)

    jammed = bool(sol.t_events[0].size)
    t_jam = float(sol.t_events[0][0]) if jammed else float("nan")

    assert t.shape == f.shape == sigma.shape, "trajectory arrays must align"
    return TransientResult(
        t=t, f=f, sigma=sigma, eta=eta, gamma_dot=gd, strain=strain,
        jammed=jammed, t_jam=t_jam,
    )


def quasi_steady_ratio(
    gamma_dot_peak: float,
    tau_event: float,
    mat: STFMaterial,
    *,
    f0: float = 0.0,
) -> tuple[float, float, float]:
    """Peak stress delivered against what quasi-steady design would predict.

    Drives a constant shear rate ``gamma_dot_peak`` for ``tau_event`` (the
    harshest case for the quasi-steady assumption, since any softer ramp
    accumulates strain more slowly) and compares the peak stress reached with
    the steady stress at the same rate.

    Returns ``(ratio, sigma_transient, sigma_steady)``. A ratio below one is the
    shortfall the device literature does not model. ``nan`` if the steady state
    is jammed, where the ratio is not defined.
    """
    res = transient_response(lambda _t: gamma_dot_peak, mat, tau_event, f0=f0)
    sigma_qs = steady_stress(gamma_dot_peak, mat)
    if not np.isfinite(sigma_qs) or sigma_qs <= 0:
        return float("nan"), res.sigma_peak, sigma_qs
    return res.sigma_peak / sigma_qs, res.sigma_peak, sigma_qs
