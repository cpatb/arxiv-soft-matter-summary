"""Regularised viscosity laws for solvers, with the knob and its error made explicit.

A yield-stress law has eta -> inf as gammadot -> 0, and no CFD solver can take that.
Every practical solver therefore evaluates a regularised law instead, and the choice of
scheme and knob is a modelling decision that changes the answer (the
``yfrac-is-a-clamp-readout`` episode: an observable defined through the regularised law
measured the knob, not the fluid). ``regularise`` refuses to build a model without the
knob, names the scheme in the model's ``name`` and ``fixed`` so a solver export cannot
lose it, and ``induced_error`` reports the relative viscosity error against the exact
law at any rate, so a study can sweep the knob before quoting a number.

Schemes (yield term written as it replaces ``tau_y/gd`` in ``eta``):

papanastasiou (knob ``m``, s)
    ``tau_y (1 - exp(-m gd)) / gd``. Finite at gd = 0, where it equals ``tau_y m``.
    Papanastasiou, J. Rheol. 31, 385 (1987) [UNVERIFIED: not DOI-resolved in this
    session; the formula is the one every CFD textbook prints and the tests check it
    against the exact law, not against the paper].
bercovier_engelman (knob ``eps``, 1/s)
    ``tau_y / sqrt(gd^2 + eps^2)``. Finite at gd = 0, where it equals ``tau_y/eps``.
    Bercovier & Engelman, J. Comput. Phys. 36, 313 (1980) [UNVERIFIED, same caveat].
clamp (knob ``eta_max``, Pa s)
    ``min(eta_exact, eta_max)``. This is the bi-viscosity model of the literature: the
    Newtonian branch ``eta_max gd`` meets the exact stress ``tau_y + k gd^n`` at exactly
    the rate where ``eta_exact = eta_max``, so the two definitions are the same function.
    ``tests/test_regularise.py`` checks this identity numerically rather than asserting
    it.

For Herschel-Bulkley with n < 1 the power-law term ``k gd^(n-1)`` still diverges at
gd -> 0 under the first two schemes; only ``clamp`` bounds the whole viscosity. The
tests make that distinction explicit because a solver author will otherwise assume
"regularised" means "bounded".
"""

from __future__ import annotations

from dataclasses import replace
from typing import Callable

import numpy as np

from .base import RateModel

__all__ = ["SCHEMES", "regularise", "induced_error", "critical_rate"]

# scheme -> (knob name, knob unit)
SCHEMES: dict[str, tuple[str, str]] = {
    "papanastasiou": ("m", "s"),
    "bercovier_engelman": ("eps", "1/s"),
    "clamp": ("eta_max", "Pa*s"),
}


# --- yield-term replacements --------------------------------------------------------
# Each takes (tau_y, gd, knob) and returns the regularised value of tau_y/gd.

def _yield_papanastasiou(tau_y: float, gd: np.ndarray, m: float) -> np.ndarray:
    # -expm1(-m gd) = 1 - exp(-m gd), stable for small m gd; the gd = 0 limit is tau_y m.
    with np.errstate(divide="ignore", invalid="ignore"):
        out = tau_y * (-np.expm1(-m * gd)) / gd
    return np.where(gd > 0, out, tau_y * m)


def _yield_bercovier_engelman(tau_y: float, gd: np.ndarray, eps: float) -> np.ndarray:
    return tau_y / np.sqrt(gd * gd + eps * eps)


_YIELD: dict[str, Callable[[float, np.ndarray, float], np.ndarray]] = {
    "papanastasiou": _yield_papanastasiou,
    "bercovier_engelman": _yield_bercovier_engelman,
}


# --- per-model regularised functions -------------------------------------------------
# Keyword-only knob so that dst_fit's Model.fixed / accepts_fixed machinery carries it.

def _hb_regularised(scheme: str):
    yield_term = _YIELD[scheme]
    knob = SCHEMES[scheme][0]

    if knob == "m":
        def f(gd, tau_y, k, n, *, m):
            with np.errstate(divide="ignore"):
                return yield_term(tau_y, gd, m) + k * gd ** (n - 1.0)
    else:
        def f(gd, tau_y, k, n, *, eps):
            with np.errstate(divide="ignore"):
                return yield_term(tau_y, gd, eps) + k * gd ** (n - 1.0)
    return f


def _casson_regularised(scheme: str):
    # sqrt(sigma) = sqrt(tau_y) + sqrt(eta_c gd) => eta = (sqrt(tau_y/gd) + sqrt(eta_c))^2,
    # and the scheme replaces tau_y/gd inside the first square root.
    yield_term = _YIELD[scheme]
    knob = SCHEMES[scheme][0]

    if knob == "m":
        def f(gd, tau_y, eta_c, *, m):
            return (np.sqrt(yield_term(tau_y, gd, m)) + np.sqrt(eta_c)) ** 2
    else:
        def f(gd, tau_y, eta_c, *, eps):
            return (np.sqrt(yield_term(tau_y, gd, eps)) + np.sqrt(eta_c)) ** 2
    return f


def _clamped(base: RateModel):
    def f(gd, *params, eta_max):
        with np.errstate(divide="ignore", invalid="ignore"):
            eta = base.func(gd, *params, **base.fixed)
        return np.minimum(eta, eta_max)
    return f


_YIELD_MODELS = {"herschel_bulkley": _hb_regularised, "casson": _casson_regularised}


def regularise(model: RateModel, scheme: str, **knob: float) -> RateModel:
    """Return ``model`` regularised by ``scheme`` with the knob value given.

    ``knob`` must be exactly the scheme's knob (``m``, ``eps`` or ``eta_max``), positive.
    The result keeps the base model's parameter names and bounds, so it fits and
    profiles like the exact model, and records ``base`` for ``induced_error``.
    """
    if model.is_regularised:
        raise ValueError(f"{model.name!r} is already regularised")
    if scheme not in SCHEMES:
        raise ValueError(f"unknown scheme {scheme!r}; known: {sorted(SCHEMES)}")
    knob_name, _ = SCHEMES[scheme]
    if set(knob) != {knob_name}:
        raise ValueError(
            f"scheme {scheme!r} takes exactly one knob, {knob_name!r}; got {sorted(knob)}"
        )
    value = float(knob[knob_name])
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{knob_name} must be finite and positive, got {value}")

    if scheme == "clamp":
        func = _clamped(model)
    else:
        # The yield-term schemes below rebuild the law from the model's PARAMETERS only;
        # unlike the clamp, which wraps base.func, they never see base.fixed, so a fixed
        # constant on the base model would be silently dropped. No shipped model has one,
        # and this refuses rather than assuming that stays true (NEXT_STEPS item 5).
        if model.fixed:
            raise ValueError(
                f"scheme {scheme!r} rebuilds {model.name!r} from its parameters and does "
                f"not carry fixed constants, but {model.name!r} fixes {sorted(model.fixed)}; "
                "regularising it would silently drop them. Use 'clamp', which wraps the "
                "model's own func and does receive them.")
        if model.name not in _YIELD_MODELS:
            raise ValueError(
                f"{scheme!r} regularises a yield term; {model.name!r} has none "
                f"(use 'clamp', or one of {sorted(_YIELD_MODELS)})"
            )
        func = _YIELD_MODELS[model.name](scheme)

    return replace(
        model,
        name=f"{model.name}+{scheme}",
        func=func,
        # The clamp wraps base.func, which already receives base.fixed; the yield-term
        # schemes are written for models that have no fixed constants.
        fixed={knob_name: value},
        notes=f"{model.notes}; regularised: {scheme} with {knob_name} = {value:g}",
        regularisation={scheme: value},
        base=model,
        limits={
            **model.limits,
            "gd -> 0 (regularised)": {
                "papanastasiou": "yield term -> tau_y m; power-law term still diverges if n < 1",
                "bercovier_engelman": "yield term -> tau_y/eps; power-law term still diverges if n < 1",
                "clamp": "eta -> eta_max, bounded",
            }[scheme],
        },
    )


def induced_error(reg: RateModel, gammadot, params) -> np.ndarray:
    """Relative viscosity error ``(eta_reg - eta_exact)/eta_exact`` at each rate.

    Negative where the regularisation under-predicts the viscosity (always, for the
    schemes here: each replaces ``tau_y/gd`` by something smaller). Where the exact
    viscosity is infinite (gd = 0) the error is ``-1`` by convention: the regularised
    law is finite there, so the whole viscosity is missing.
    """
    if not reg.is_regularised or reg.base is None:
        raise ValueError(f"{reg.name!r} is not a regularised model")
    gd = np.asarray(gammadot, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        exact = reg.base(gd, *params)
        approx = reg(gd, *params)
        err = (approx - exact) / exact
    return np.where(np.isfinite(exact), err, -1.0)


def critical_rate(reg: RateModel, params, rtol: float) -> float:
    """Smallest rate above which ``|induced_error| < rtol`` holds, found by bisection.

    Walks a log grid from 1e-6 to 1e6 1/s, locates the last crossing of ``rtol`` and
    refines it with ``scipy.optimize.brentq``. Raises if the error never drops below
    ``rtol`` on the grid (a knob too coarse for the requested tolerance) or is below it
    everywhere (nothing to find).
    """
    from scipy.optimize import brentq

    # (eta_reg - eta)/eta is a difference of two doubles; below ~1e-12 it is rounding.
    if rtol < 1e-12:
        raise ValueError(f"rtol={rtol} is below what a relative difference can resolve")
    grid = np.logspace(-6, 6, 241)
    err = np.abs(induced_error(reg, grid, params))
    above = err >= rtol
    if not above.any():
        raise ValueError(f"induced error is below rtol={rtol} on the whole grid")
    if above.all() or above[-1]:
        raise ValueError(f"induced error never drops below rtol={rtol} up to 1e6 1/s")
    i = int(np.flatnonzero(above)[-1])  # last grid point still above tolerance
    g = lambda lg: abs(float(induced_error(reg, 10.0 ** lg, params))) - rtol
    return float(10.0 ** brentq(g, np.log10(grid[i]), np.log10(grid[i + 1])))
