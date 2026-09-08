"""Device design-space questions built on :mod:`biorheo.apps.injectability`.

Everything here answers one of the three questions a device or formulation group asks
after a viscosity curve has been fitted, and each is a *derivative or a threshold* of the
force the injectability layer already computes, not a new physical model.

1. **How much does a finer needle cost me?**  :func:`bore_exponent` measures

       b = d ln F / d ln R                (fixed Q, L, syringe; friction excluded)

   by central difference in ``ln R``. Two closed forms bound it and neither passes
   through the force solver, which is what makes them usable as anchors
   (:func:`bore_exponent_closed_form`):

   - Newtonian, ``dp = 8 mu L Q / (pi R^4)``  ->  ``b = -4`` exactly, for any ``mu, Q, L``.
   - Power law, ``sigma = k gd^n``: ``Q = (pi n/(3n+1)) R^3 (tau_w/k)^(1/n)`` with
     ``tau_w = dp R/(2L)``, so ``Q ~ R^(3+1/n) dp^(1/n)`` and ``dp ~ Q^n R^-(3n+1)``,
     i.e. ``b = -(3n+1)`` exactly.

   For a shear-thinning fluid ``n < 1`` and ``b > -4``: shrinking the bore raises the wall
   shear rate, which thins the fluid, so the force rises *less steeply* than the Newtonian
   ``R^-4`` rule. :func:`newtonian_rule_error` reports the size of that over-prediction,
   which is the quantity the Roche 2018 injection-time paper describes qualitatively.

   ``b`` is a property of the fluid AND the operating point. It is not a constant of the
   material: a fluid sitting on its infinite-shear plateau at the wall is Newtonian there
   and returns ``-4`` however shear-thinning it is at rheometer rates.

2. **What flow rate fits inside a force budget?**  :func:`max_flow_rate` inverts
   ``F(Q) = budget`` on a logarithmic bracket. ``F`` is monotone increasing in ``Q`` for
   every generalised-Newtonian model in this registry, which is asserted rather than
   assumed (the bracket search fails loudly if the sign does not change).

3. **How much of that budget does my parameter uncertainty eat?**  :func:`derating`
   computes ``Q_max`` twice: at the fitted parameters, and against the one-at-a-time
   profile push-through *upper* force envelope. The ratio is a derating factor, and it is
   the only number here that a designer can act on directly.

What a budget number is and is not
----------------------------------
``max_flow_rate`` crosses a *model* against a *number the caller supplies*, with every
omission in :mod:`biorheo.apps.injectability` still omitted -- no barrel pressure drop, no
entrance loss, no wall slip, no needle-hub friction, no plunger stiction. The **ratio** of
two such crossings is meaningful because the omissions largely cancel; the **absolute**
crossing is not, and callers are expected to say so wherever they print one.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from ..constitutive.base import RateModel
from .injectability import Needle, Syringe, injection_force

__all__ = [
    "BORE_STEP", "bore_exponent", "bore_exponent_closed_form",
    "newtonian_rule_error", "max_flow_rate", "DeratingResult", "derating",
]

BORE_STEP = 1e-3        # relative step in ln R for the central difference
_Q_BRACKET = (1e-12, 1e-3)   # m^3/s; 6e-5 to 6e4 mL/min


def bore_exponent(Q: float, needle: Needle, syringe: Syringe, model: RateModel, params,
                  *, step: float = BORE_STEP) -> float:
    """``d ln F / d ln R`` at fixed ``Q``, by central difference in ``ln R``.

    ``friction`` is deliberately not a parameter: it is a property of the syringe, not the
    fluid, and including it would make ``b`` depend on a constant the fluid cannot change.
    """
    params = np.asarray(params, dtype=float)

    def lnF(scale: float) -> float:
        n = Needle(radius=needle.radius * scale, length=needle.length)
        return float(np.log(injection_force(Q, n, syringe, model, params).force))

    return float((lnF(np.exp(step)) - lnF(np.exp(-step))) / (2.0 * step))


def bore_exponent_closed_form(model: RateModel, params) -> float | None:
    """The analytic ``b`` where one exists, computed from the parameters alone.

    Returns ``None`` for every model that has no closed form -- this function must never
    fall back on the numerical route, because its whole purpose is to be an independent
    check of it.
    """
    pd = model.params_dict(np.asarray(params, dtype=float))
    if model.name == "newtonian":
        return -4.0
    if model.name == "power_law":
        return -(3.0 * float(pd["n"]) + 1.0)
    return None


def newtonian_rule_error(Q: float, needle: Needle, syringe: Syringe, model: RateModel,
                         params, *, bore_ratio: float) -> dict[str, float]:
    """How wrong a Newtonian ``R^-4`` sizing rule is at a finer bore.

    The rule takes the force already known at ``needle`` and scales it by
    ``bore_ratio ** -4``; the truth is the force the fitted law gives at the finer needle.
    Returns both, their ratio, and the relative over-prediction
    ``F_rule / F_true - 1`` (positive = the Newtonian rule is pessimistic).
    """
    if not 0.0 < bore_ratio <= 1.0:
        raise ValueError(f"bore_ratio must be in (0, 1], got {bore_ratio}")
    params = np.asarray(params, dtype=float)
    f_ref = injection_force(Q, needle, syringe, model, params).force
    fine = Needle(radius=needle.radius * bore_ratio, length=needle.length)
    f_true = injection_force(Q, fine, syringe, model, params).force
    f_rule = f_ref * bore_ratio ** -4.0
    return {"bore_ratio": float(bore_ratio), "F_ref": float(f_ref),
            "F_true": float(f_true), "F_newtonian_rule": float(f_rule),
            "over_prediction": float(f_rule / f_true - 1.0)}


def max_flow_rate(budget: float, needle: Needle, syringe: Syringe, model: RateModel,
                  params, *, friction: float = 0.0,
                  bracket: tuple[float, float] = _Q_BRACKET) -> float:
    """The ``Q`` at which the plunger force reaches ``budget`` newtons.

    Raises if the bracket does not straddle the budget, rather than returning an endpoint:
    a device that cannot reach the budget anywhere in ``bracket`` is a finding, not a
    number to be clipped.
    """
    params = np.asarray(params, dtype=float)

    def g(logQ: float) -> float:
        F = injection_force(float(np.exp(logQ)), needle, syringe, model, params,
                            friction=friction).force
        return F - budget

    lo, hi = (float(np.log(b)) for b in bracket)
    g_lo, g_hi = g(lo), g(hi)
    if g_lo > 0:
        raise ValueError(f"force already exceeds the {budget} N budget at Q = {bracket[0]:g} m^3/s "
                         f"(F = {g_lo + budget:.4g} N): the budget is unreachable from below")
    if g_hi < 0:
        raise ValueError(f"force never reaches the {budget} N budget by Q = {bracket[1]:g} m^3/s "
                         f"(F = {g_hi + budget:.4g} N)")
    return float(np.exp(float(brentq(g, lo, hi, xtol=1e-10))))


@dataclass(frozen=True)
class DeratingResult:
    """``Q_max`` at the fit and against the upper force envelope, and their ratio."""
    budget: float
    q_nominal: float            # m^3/s
    q_upper_bound: float        # m^3/s, using the highest force any profile end gives
    factor: float               # q_upper_bound / q_nominal, <= 1
    binding_parameter: str      # the parameter whose profile end sets the upper envelope
    binding_end: str            # "lo" or "hi"
    open_parameters: tuple[str, ...]   # profile ends that hit a box bound

    def __str__(self) -> str:
        ml_min = 6e7            # m^3/s -> mL/min
        note = ("" if not self.open_parameters
                else f"  [OPEN in {', '.join(self.open_parameters)}: this is a bound, not a limit]")
        return (f"Q_max at the fit {self.q_nominal * ml_min:.4g} mL/min; "
                f"at the 95 % upper force envelope {self.q_upper_bound * ml_min:.4g} mL/min; "
                f"derating factor {self.factor:.4f} "
                f"(set by {self.binding_parameter} at its {self.binding_end} end){note}")


def derating(budget: float, needle: Needle, syringe: Syringe, model: RateModel, params,
             intervals: dict[str, tuple[float, float]],
             *, friction: float = 0.0, open_parameters: tuple[str, ...] = (),
             bracket: tuple[float, float] = _Q_BRACKET) -> DeratingResult:
    """``Q_max`` at the fit, and ``Q_max`` against the one-at-a-time profile push-through
    *upper* force envelope.

    ``intervals`` maps each parameter name to its ``(lo, hi)`` profile interval. The upper
    force envelope is the largest force any single-parameter push gives -- the same
    one-at-a-time hull :func:`biorheo.apps.injectability.force_with_uq` reports, which is
    WIDER than the joint 95 % region it is drawn from. The factor is therefore
    conservative in the direction a designer wants, and the caller must say which interval
    it is quoting.
    """
    params = np.asarray(params, dtype=float)
    q_nom = max_flow_rate(budget, needle, syringe, model, params,
                          friction=friction, bracket=bracket)

    # find the single-parameter push that gives the largest force at q_nom, then invert
    # the budget against that parameter vector. F is monotone in Q for every model here,
    # so the largest force at one Q is the smallest Q_max.
    best, best_F, best_name, best_end = None, -np.inf, "", ""
    for j, name in enumerate(model.param_names):
        if name not in intervals:
            continue
        for end, val in zip(("lo", "hi"), intervals[name]):
            if not np.isfinite(val):
                continue
            p = params.copy()
            p[j] = val
            try:
                F = injection_force(q_nom, needle, syringe, model, p, friction=friction).force
            except (ValueError, RuntimeError):
                continue
            if F > best_F:
                best, best_F, best_name, best_end = p, F, name, end
    if best is None:
        raise ValueError("no usable profile end: every interval was missing or non-finite")

    q_hi = max_flow_rate(budget, needle, syringe, model, best,
                         friction=friction, bracket=bracket)
    return DeratingResult(budget=float(budget), q_nominal=float(q_nom),
                          q_upper_bound=float(q_hi), factor=float(q_hi / q_nom),
                          binding_parameter=best_name, binding_end=best_end,
                          open_parameters=tuple(open_parameters))
