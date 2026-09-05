"""Structural (thixotropic) family: a scalar structure parameter with first-order
kinetics, Moore / Coussot lineage.

Ported from ``Desktop/Projects/bio_fluids_bridge/blood_rheology.py`` (read, not
re-derived; SymPy checks in ``derivations/structural_kinetics.py``). The source writes

    dlam/dt = k_build (1 - lam) - k_break gd lam,     eta = eta_inf (1 + (eta_0/eta_inf - 1) lam)

which is written here with ``t_r = 1/k_build`` (rebuild time, s), ``alpha = k_break``
(breakdown per unit strain, dimensionless) and ``beta = eta_0/eta_inf - 1`` (structural
viscosity gain, dimensionless):

    dlam/dt = (1 - lam)/t_r - alpha gd lam,           eta = eta_inf (1 + beta lam).

Three facts, each derived in ``derivations/structural_kinetics.md`` and each a test:

1. Steady state ``lam_ss = 1/(1 + alpha t_r gd)``. ``alpha`` and ``t_r`` enter only as
   the product, so a steady flow curve identifies ``alpha_tr = alpha t_r`` and nothing
   else about the kinetics. :data:`STRUCTURAL_STEADY` is therefore parametrised by the
   product. :data:`STRUCTURAL_STEADY_SEPARATE` keeps ``alpha`` and ``t_r`` apart and is
   NOT in the registry: it exists so ``tests/test_structural.py`` can show the identify
   layer returning ``degenerate`` for the pair, which is the gate for the
   ``parameters-that-enter-as-a-ratio`` lesson.
2. The steady viscosity is the Cross model with ``m = 1``, ``eta_0 = eta_inf (1 + beta)``,
   ``lam_cross = alpha t_r``. A steady flow curve cannot tell "structural" from "Cross":
   only a transient can, and only a transient carries ``t_r`` on its own.
3. At constant rate ``lam(t) = lam_ss + (lam0 - lam_ss) exp(-t/tau)`` with
   ``1/tau = 1/t_r + alpha gd``; the rest relaxation time is ``t_r``.

Against human blood (2026-09-05, ``data/horner2020/step_response.md``): after a
300 -> 1 1/s step the measured stress rises to a maximum at 4-6 s and then decays by
6-31 % on all 20 Horner 2020 donors. Fact 3 makes ``lam(t)`` and hence the stress
monotone at constant rate for any parameter values, so this family can describe the rise
and not the decay. The thesis attributes the decay to syneresis (a wall depletion layer,
``data/horner2020/SOURCE.md``), not to bulk rheology; on that reading the family is not
refuted by the sheet, and no model should be fitted to the decay. A first reading of the
sheet (commit 39b0988) called the family "refuted on shape"; withdrawn once the thesis
was read. Fitted to the rise (``scripts/structural_step_horner.py``, pinned in
``tests/test_structural_step.py``): the timescale is right on 20/20 donors (``tau`` =
1.1-2.2 s, 95 % approach at 3.4-6.6 s, where the peak is) and the SHAPE is not -- a
single exponential misses the rise by 1-12 % of its amplitude (median 3.4 %) against a
sample noise of 0.04-1 mPa, with the same residual sign pattern on 15/20, and anchoring
both endpoints to the steady sheet instead of to this fit does not change that. The
rise is less concave than one exponential; one structure variable with first-order
kinetics cannot give that. What it would take (a lag, a second timescale) is untested.

The source also added a Casson-like yield offset ``sigma_y`` to the steady stress. It is
left out here: it is a fourth steady parameter with its own identifiability question, and
the registry composes a yield stress explicitly rather than hiding one inside a "structural"
model. The source's constants are ballpark, not fitted (its own README says so), and are
not shipped as a record.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

from .base import RateModel

__all__ = [
    "STRUCTURAL_STEADY", "STRUCTURAL_STEADY_SEPARATE",
    "lambda_steady", "relaxation_time", "lambda_step", "lambda_transient",
    "eta_of_lambda", "RampLoop", "ramp_loop",
]


# --- steady state ------------------------------------------------------------------

def lambda_steady(gd, alpha_tr: float):
    """``1/(1 + alpha_tr gd)``, the steady structure fraction; in (0, 1]."""
    return 1.0 / (1.0 + alpha_tr * np.asarray(gd, dtype=float))


def eta_of_lambda(lam, eta_inf: float, beta: float):
    """``eta_inf (1 + beta lam)``: ``eta_inf`` at lam = 0, ``eta_inf (1 + beta)`` at lam = 1."""
    return eta_inf * (1.0 + beta * np.asarray(lam, dtype=float))


def _structural_steady(gd: np.ndarray, eta_inf: float, beta: float, alpha_tr: float) -> np.ndarray:
    return eta_of_lambda(lambda_steady(gd, alpha_tr), eta_inf, beta)


def _guess_structural(gd, eta):
    eta_0, eta_inf = np.max(eta), np.min(eta)
    beta = max(eta_0 / eta_inf - 1.0, 1e-6)
    # alpha_tr gd = 1 where eta sits halfway between the plateaux (limit 2 above)
    i = int(np.argmin(np.abs(eta - 0.5 * (eta_0 + eta_inf))))
    return np.array([eta_inf, beta, 1.0 / gd[i]])


STRUCTURAL_STEADY = RateModel(
    name="structural_steady",
    param_names=("eta_inf", "beta", "alpha_tr"),
    func=_structural_steady,
    bounds=((1e-6, 0.0, 1e-6), (1e6, 1e6, 1e6)),
    describes_thickening=False,
    notes=("eta = eta_inf (1 + beta lam_ss), lam_ss = 1/(1 + alpha_tr gd); the steady state "
           "of first-order structure kinetics, identical to Cross with m = 1"),
    guess_fn=_guess_structural,
    limits={
        "gd -> 0": "eta -> eta_inf (1 + beta)",
        "gd -> inf": "eta -> eta_inf",
        "alpha_tr gd = 1": "eta = eta_inf (1 + beta/2), the plateau mean",
        "any gd": "equals Cross(eta_0 = eta_inf (1 + beta), eta_inf, lam = alpha_tr, m = 1)",
        "beta = 0": "Newtonian eta_inf at every rate (alpha_tr drops out)",
    },
    units={"eta_inf": "Pa*s", "beta": "dimensionless", "alpha_tr": "s"},
)


def _structural_steady_separate(gd, eta_inf, beta, alpha, t_r):
    return _structural_steady(gd, eta_inf, beta, alpha * t_r)


def _guess_structural_separate(gd, eta):
    eta_inf, beta, alpha_tr = _guess_structural(gd, eta)
    return np.array([eta_inf, beta, np.sqrt(alpha_tr), np.sqrt(alpha_tr)])


STRUCTURAL_STEADY_SEPARATE = RateModel(
    name="structural_steady_separate",
    param_names=("eta_inf", "beta", "alpha", "t_r"),
    func=_structural_steady_separate,
    bounds=((1e-6, 0.0, 1e-6, 1e-6), (1e6, 1e6, 1e6, 1e6)),
    describes_thickening=False,
    notes=("DIAGNOSTIC ONLY, not in MODELS: the steady structural model with alpha and t_r "
           "kept separate although only their product enters; fitting it to a flow curve "
           "must return degenerate for the pair"),
    guess_fn=_guess_structural_separate,
    limits={"any gd": "equals structural_steady with alpha_tr = alpha t_r"},
    units={"eta_inf": "Pa*s", "beta": "dimensionless", "alpha": "dimensionless", "t_r": "s"},
)


# --- transients --------------------------------------------------------------------

def relaxation_time(gd: float, alpha: float, t_r: float) -> float:
    """``1/(1/t_r + alpha gd)``: the structure's relaxation time at constant rate."""
    return 1.0 / (1.0 / t_r + alpha * gd)


def lambda_step(t, gd: float, alpha: float, t_r: float, lam0: float):
    """Closed-form structure after a step to constant ``gd`` from ``lam0`` at t = 0."""
    t = np.asarray(t, dtype=float)
    lam_ss = lambda_steady(gd, alpha * t_r)
    return lam_ss + (lam0 - lam_ss) * np.exp(-t / relaxation_time(gd, alpha, t_r))


def lambda_transient(t_eval, gd_of_t: Callable[[float], float], alpha: float, t_r: float,
                     lam0: float, *, rtol: float = 1e-8, atol: float = 1e-10) -> np.ndarray:
    """Integrate ``dlam/dt = (1 - lam)/t_r - alpha gd(t) lam`` from ``lam0`` on ``t_eval``.

    ``gd_of_t`` is any rate history (a step, a ramp, a measured protocol). The ODE is
    linear in ``lam`` and stiff when ``alpha gd t_r >> 1``, hence LSODA, as in the source.
    """
    t_eval = np.asarray(t_eval, dtype=float)
    if t_eval.ndim != 1 or t_eval.size < 2 or np.any(np.diff(t_eval) <= 0):
        raise ValueError("t_eval must be a strictly increasing 1-d array of >= 2 times")
    if not 0.0 <= lam0 <= 1.0:
        raise ValueError(f"lam0 = {lam0} is outside [0, 1]")

    def rhs(t, y):
        return [(1.0 - y[0]) / t_r - alpha * float(gd_of_t(t)) * y[0]]

    sol = solve_ivp(rhs, (t_eval[0], t_eval[-1]), [lam0], t_eval=t_eval,
                    method="LSODA", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"lambda_transient: {sol.message}")
    lam = sol.y[0]
    # the exact solution stays in [0, 1] for lam0 in [0, 1]; allow only solver rounding
    if lam.min() < -1e-6 or lam.max() > 1.0 + 1e-6:
        raise RuntimeError(f"lambda_transient left [0, 1]: [{lam.min():.3g}, {lam.max():.3g}]")
    return np.clip(lam, 0.0, 1.0)


@dataclass(frozen=True)
class RampLoop:
    """Up-then-down linear rate ramp: ``t``, ``gd``, ``lam``, ``eta``, ``sigma`` and the
    signed loop area ``area = oint sigma d(gd)`` (positive when the up branch carries
    the higher stress, which is thixotropic hysteresis)."""
    t: np.ndarray
    gd: np.ndarray
    lam: np.ndarray
    eta: np.ndarray
    sigma: np.ndarray
    area: float
    n_up: int  # samples on the up branch; the down branch is the rest

    @property
    def up(self) -> tuple[np.ndarray, np.ndarray]:
        return self.gd[: self.n_up], self.sigma[: self.n_up]

    @property
    def down(self) -> tuple[np.ndarray, np.ndarray]:
        return self.gd[self.n_up:], self.sigma[self.n_up:]


def ramp_loop(eta_inf: float, beta: float, alpha: float, t_r: float, *,
              gd_max: float, t_ramp: float, n: int = 600, lam0: float = 1.0) -> RampLoop:
    """Shear-rate ramp 0 -> ``gd_max`` -> 0, each leg ``t_ramp`` long, from ``lam0``.

    The loop area is a measurement of ``t_r`` that the steady curve cannot make: it
    vanishes as ``t_ramp/t_r -> inf`` (quasi-static, both branches on the steady curve)
    and is positive for any finite ramp.
    """
    if n < 3:
        raise ValueError("n must be >= 3")
    t_up = np.linspace(0.0, t_ramp, n)
    t_dn = np.linspace(t_ramp, 2 * t_ramp, n + 1)[1:]  # no duplicate apex sample
    t = np.concatenate([t_up, t_dn])

    def gd_of_t(tt: float) -> float:
        return gd_max * (tt / t_ramp if tt <= t_ramp else 1.0 - (tt - t_ramp) / t_ramp)

    lam = lambda_transient(t, gd_of_t, alpha, t_r, lam0)
    gd = np.array([gd_of_t(tt) for tt in t])
    eta = eta_of_lambda(lam, eta_inf, beta)
    sigma = eta * gd
    # signed area by the trapezoid rule around the closed path (last point is at gd = 0
    # again, so the path closes on the first point up to solver tolerance)
    area = float(np.trapezoid(sigma, gd))
    return RampLoop(t=t, gd=gd, lam=lam, eta=eta, sigma=sigma, area=area, n_up=n)
