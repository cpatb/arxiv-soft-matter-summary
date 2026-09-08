"""A temperature axis for a record that has none, without inventing one.

Derived and pre-registered in ``data/temperature_shift_prereg.md`` before this file
existed. **No record in this registry carries a temperature axis and this module does not
add one.** It carries a shift the CALLER supplies, exactly, and refuses to guess.

The assumption, stated so it can be refused
-------------------------------------------
Take the fluid to be *thermorheologically simple* over the range in question: changing the
temperature rescales the viscosity and the rate axis together and does not change the shape
of the flow curve. Then, with a dimensionless shift factor ``a > 0``,

    eta(gd; T)  =  a * eta_ref(a * gd)

``a > 1`` is the colder, more viscous fluid. The whole content of the assumption is that
ONE number describes the entire curve's move -- which is falsifiable by the client's own
data, and :func:`collapse_shift` returns the residual so it is checked rather than asserted.

The parameter maps
------------------
Substituting ``a*gd`` and multiplying by ``a`` gives a map on the fitted parameters,
derived in the pre-registration's §2 and checked numerically by ``tests/test_shift.py``:

===================  ====================================================================
``newtonian``        ``mu -> a mu``
``power_law``        ``k -> a^n k``; ``n`` fixed
``herschel_bulkley`` ``tau_y`` UNCHANGED; ``k -> a^n k``; ``n`` fixed
``casson``           ``tau_y`` UNCHANGED; ``eta_c -> a eta_c``
``cross``            ``eta_0, eta_inf -> a *``; ``lam -> a lam``; ``m`` fixed
``carreau_yasuda``   ``eta_0, eta_inf -> a *``; ``lam -> a lam``; ``a``, ``n`` fixed
===================  ====================================================================

**The yield stress comes out invariant in both yield-stress families and that was not
imposed.** For Casson the ``sqrt(a)`` outside cancels the ``sqrt(1/a)`` inside
(``sqrt(a) sqrt(tau_y/(a gd)) = sqrt(tau_y/gd)``); for Herschel-Bulkley the ``a`` cancels
the ``1/a`` in ``tau_y/(a gd)``. A yield stress is a stress scale set by structure, not a
viscosity, so this is the right answer -- but here it is a consequence of the
superposition assumption rather than a choice, and the opposite result would have been a
reason to distrust the construction.

The envelope moves, and forgetting that is the trap
---------------------------------------------------
A record's envelope is where its curve was MEASURED. The shifted curve needs ``a*gd``
inside that envelope, so it is valid only on ``[gd_min/a, gd_max/a]``. A shift that left
the envelope alone would silently widen the record's claim -- the same failure this package
already records for the OpenFOAM export, where a dictionary has no envelope at all.

Where ``a`` comes from, and it is never from here
-------------------------------------------------
:func:`shift_from_two_viscosities` (the ratio of two measured viscosities at the same
shear rate) or :func:`arrhenius_shift` (requiring the caller's own activation energy;
``R`` comes from :mod:`scipy.constants`, not from memory). **There is no third route and
no default activation energy is shipped for any fluid.**
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
from scipy.constants import R as GAS_CONSTANT      # J mol^-1 K^-1, exact since 2019
from scipy.optimize import minimize_scalar

from ..registry.schema import Envelope, FluidRecord, Parameter
from .base import RateModel

__all__ = [
    "GAS_CONSTANT", "ThermalShiftNotRepresentable", "SHIFTABLE",
    "arrhenius_shift", "shift_from_two_viscosities",
    "shift_parameters", "shift_envelope", "shift_record", "collapse_shift",
]


class ThermalShiftNotRepresentable(NotImplementedError):
    """The family's temperature dependence is not a one-parameter curve shift."""


# exponent of ``a`` applied to each parameter, per family. A parameter absent from a
# family's map is UNCHANGED, and that is a claim, not a default: see the module docstring.
_MAP: dict[str, dict[str, float]] = {
    "newtonian": {"mu": 1.0},
    "power_law": {"k": float("nan")},          # a^n, filled from the record's own n
    "herschel_bulkley": {"k": float("nan")},   # a^n; tau_y unchanged
    "casson": {"eta_c": 1.0},                  # tau_y unchanged
    "cross": {"eta_0": 1.0, "eta_inf": 1.0, "lam": 1.0},
    "carreau_yasuda": {"eta_0": 1.0, "eta_inf": 1.0, "lam": 1.0},
}
SHIFTABLE = tuple(_MAP)


# Below this, an argument said to be in kelvin is a celsius value that was not converted.
# 100 K is far under any temperature this package's fluids are measured at and far over any
# celsius value one would pass, so it separates the two cases without excluding real data.
MIN_PLAUSIBLE_KELVIN = 100.0


def arrhenius_shift(E_a: float, T: float, T_ref: float) -> float:
    """``exp(E_a/R (1/T - 1/T_ref))``, temperatures in KELVIN.

    ``E_a`` is the CALLER's activation energy in J/mol. No default exists for any fluid in
    this package, and none will: an activation energy is a measurement, and putting a
    plausible one here would be exactly the kind of invented constant this repo forbids.
    """
    if not (T > 0 and T_ref > 0):
        raise ValueError(f"temperatures must be in kelvin and positive, got {T}, {T_ref}")
    # `T > 0` does not distinguish kelvin from celsius, and the only per-record temperature
    # field in the schema is `envelope.temperature_C` -- so the obvious wrong call,
    # `arrhenius_shift(E_a, rec.envelope.temperature_C, ref.envelope.temperature_C)`, would
    # have passed silently and returned a shift factor wrong by many orders of magnitude.
    # No fluid this package is for is studied below 100 K, so a value under that is a unit
    # error rather than a cryogenic measurement.
    if T < MIN_PLAUSIBLE_KELVIN or T_ref < MIN_PLAUSIBLE_KELVIN:
        raise ValueError(
            f"temperatures must be in KELVIN; got {T} and {T_ref}, at least one below "
            f"{MIN_PLAUSIBLE_KELVIN} K. If these are celsius, add 273.15. No fluid this "
            "package is intended for is measured below that, so this is read as a unit "
            "error rather than a cryogenic measurement.")
    return float(np.exp(E_a / GAS_CONSTANT * (1.0 / T - 1.0 / T_ref)))


def shift_from_two_viscosities(eta_at_T: float, eta_at_T_ref: float) -> float:
    """``a = eta(T) / eta(T_ref)``, both measured at the SAME shear rate.

    This is the route a client can actually supply. It is exact only if the two curves
    really do collapse under one factor; :func:`collapse_shift` is how that is checked over
    a whole curve rather than at one point.
    """
    if not (eta_at_T > 0 and eta_at_T_ref > 0):
        raise ValueError(f"viscosities must be positive, got {eta_at_T}, {eta_at_T_ref}")
    return float(eta_at_T / eta_at_T_ref)


def _exponents(model: RateModel, params: np.ndarray) -> dict[str, float]:
    if model.name not in _MAP:
        raise ThermalShiftNotRepresentable(
            f"{model.name}: a temperature shift is not representable for this family. "
            "The one-parameter superposition eta -> a*eta(a*gd) rescales a FLOW CURVE; "
            "this family carries a STRUCTURE EQUATION with its own relaxation time, and "
            "that time does not follow from the curve's shift factor. Fit the family at "
            "each temperature, or use its steady-curve counterpart and state that the "
            "transient is not covered."
        )
    out = dict(_MAP[model.name])
    if "n" in model.param_names:
        n = float(params[model.param_names.index("n")])
        for name, e in out.items():
            if np.isnan(e):
                out[name] = n
    if any(np.isnan(v) for v in out.values()):
        raise ValueError(f"{model.name}: exponent map left unfilled: {out}")
    return out


def shift_parameters(model: RateModel, params, a: float) -> np.ndarray:
    """Parameters of the curve ``a * eta_ref(a * gd)``, in the model's own order.

    At ``a == 1`` the array is returned unchanged bit for bit, deliberately: an identity
    that is only approximate would hide a defect in the map behind a tolerance.
    """
    if not a > 0:
        raise ValueError(f"shift factor must be positive, got {a}")
    p = np.asarray(params, dtype=float).copy()
    exps = _exponents(model, p)
    if a == 1.0:
        return p
    for name, e in exps.items():
        p[model.param_names.index(name)] *= a ** e
    return p


def shift_envelope(env: Envelope, a: float, *, temperature_C: float | None = None) -> Envelope:
    """The shifted curve is valid on ``[gd_min/a, gd_max/a]`` and nowhere else."""
    if not a > 0:
        raise ValueError(f"shift factor must be positive, got {a}")
    return replace(env,
                   gammadot_min=env.gammadot_min / a,
                   gammadot_max=env.gammadot_max / a,
                   temperature_C=env.temperature_C if temperature_C is None else temperature_C,
                   geometry=env.geometry + f"; time-temperature shifted by a = {a:.6g}")


def shift_record(record: FluidRecord, a: float, *,
                 temperature_C: float | None = None) -> FluidRecord:
    """A record for the same fluid at the shifted temperature.

    The result is an ordinary :class:`FluidRecord` -- it exports, evaluates and reports
    like any other -- with three differences that are all stated on it: the parameters are
    mapped, the envelope is divided by ``a``, and the identifiability sidecar is DROPPED.

    **The sidecar is dropped and not mapped.** A profile-likelihood interval belongs to a
    fit against measured data; there is no data at the shifted temperature, so there is no
    interval, and carrying the old one across would be a confidence statement about
    measurements nobody made. The resulting record therefore reports as a POINT ESTIMATE in
    `biorheo.apps.report`, which is the honest description.
    """
    p_new = shift_parameters(record.model, record.params_si, a)
    exps = _exponents(record.model, record.params_si)
    params = {}
    for i, name in enumerate(record.model.param_names):
        old = record.parameters[name]
        scale = a ** exps[name] if name in exps else 1.0
        params[name] = Parameter(name=name, value=old.value * scale, unit=old.unit,
                                 si=float(p_new[i]))
    note = (f"{record.notes}\n\nTIME-TEMPERATURE SHIFTED from `{record.id}` by a caller-"
            f"supplied factor a = {a:.6g}, under eta(gd) = a*eta_ref(a*gd). No data was "
            f"measured at this temperature. The identifiability sidecar is deliberately "
            f"NOT carried across: an interval belongs to a fit, and there is no fit here.")
    return replace(record,
                   id=f"{record.id}__shift_a{a:.6g}",
                   parameters=params,
                   envelope=shift_envelope(record.envelope, a, temperature_C=temperature_C),
                   identifiability=None,
                   notes=note,
                   path=None)


def collapse_shift(gammadot, eta_ref, eta_other, *,
                   bracket: tuple[float, float] = (1e-3, 1e3)) -> dict:
    """Best single shift factor taking ``eta_ref`` onto ``eta_other``, AND its residual.

    Both curves are given at the SAME rates. The objective is the rms of
    ``log(eta_other(gd)) - log(a * eta_ref_interp(a * gd))`` over the rates where
    ``a*gd`` stays inside the reference curve's own range -- extrapolating the reference to
    score the fit would be scoring an extrapolation.

    **The residual is returned beside the factor and callers are expected to look at it.**
    Time-temperature superposition is an assumption; two curves that differ in SHAPE still
    have a best ``a``, and it means nothing. ``n_used`` says how many rates survived the
    range restriction, because a factor fitted on three points is not a collapse either.
    """
    gd = np.asarray(gammadot, dtype=float)
    e_ref = np.asarray(eta_ref, dtype=float)
    e_oth = np.asarray(eta_other, dtype=float)
    if not (gd.shape == e_ref.shape == e_oth.shape):
        raise ValueError(f"shapes must match, got {gd.shape}, {e_ref.shape}, {e_oth.shape}")
    if gd.size < 4:
        raise ValueError(f"need at least 4 rates to fit a collapse, got {gd.size}")
    if np.any(gd <= 0) or np.any(e_ref <= 0) or np.any(e_oth <= 0):
        raise ValueError("rates and viscosities must be positive (the fit is in log space)")

    order = np.argsort(gd)
    lg, lref = np.log(gd[order]), np.log(e_ref[order])
    lother = np.log(e_oth[order])

    def residual(log_a: float) -> float:
        a = float(np.exp(log_a))
        target = lg + np.log(a)
        inside = (target >= lg[0]) & (target <= lg[-1])
        if inside.sum() < 4:
            return np.inf
        pred = np.log(a) + np.interp(target[inside], lg, lref)
        return float(np.sqrt(np.mean((lother[inside] - pred) ** 2)))

    res = minimize_scalar(residual, bounds=(np.log(bracket[0]), np.log(bracket[1])),
                          method="bounded", options={"xatol": 1e-10})
    a = float(np.exp(res.x))
    target = lg + np.log(a)
    n_used = int(((target >= lg[0]) & (target <= lg[-1])).sum())
    return {"a": a, "rms_log_residual": float(res.fun), "n_used": n_used,
            "n_total": int(gd.size), "converged": bool(res.success),
            "note": "rms residual in log(eta). A large value means the two curves do not "
                    "collapse under any single factor, i.e. the fluid is not "
                    "thermorheologically simple over this range, and the returned `a` is "
                    "meaningless rather than approximate."}
