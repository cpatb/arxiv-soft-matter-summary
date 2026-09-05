"""Fit a rate model to a flow curve and say, per parameter, whether it was measured.

The verdict rules are fixed here, before any data is looked at, and applied in code
(the ``enforce-gates-in-code`` lesson: a rule applied by hand at read time gets
skipped). They are read off ``dst_fit.profiling.profile_likelihood``, which walks the
profile likelihood outward from the best fit until the F statistic crosses its
critical value or the parameter's bound is hit:

measurable
    the profile interval closes on both sides, and neither end sits on a bound;
bounded_only
    the interval closes on one side only (the other end ran into a bound), so the data
    gives a one-sided bound and not a value;
not_identifiable
    both ends are at bounds: the profile never crossed the critical value, the data
    says nothing about this parameter on its own;
degenerate
    at least one end is open AND the parameter's fit correlation with another one has
    ``|corr| > 0.99``: it is not identifiable because it is traded against that other
    parameter (a ratio or product enters the model, the ``parameters-that-enter-as-a-
    ratio`` lesson). Checked before ``bounded_only``/``not_identifiable``, so a
    degenerate pair is named as such rather than reported as two unrelated bounds.

"Bounded" here means the likelihood-ratio interval at ``level``, not the Wald error;
``dst_fit`` reports both, and the profile is the one that is honest near a bound. A
``measurable`` verdict is a statement about closure, not precision; ``ParamEstimate.decades``
carries the width and must be quoted beside it.

This module is the only writer of a record's ``identifiability`` block.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from dst_fit import FitResult, fit_model
from dst_fit.profiling import Profile, profile_likelihood

from ..constitutive import RateModel
from ..registry.schema import Envelope

__all__ = ["VERDICTS", "ParamEstimate", "IdentifiedModel", "fit_rate_curve", "verdict_for"]

VERDICTS = ("measurable", "bounded_only", "degenerate", "not_identifiable")
DEGENERATE_CORR = 0.99  # |corr| above which two parameters count as one


@dataclass(frozen=True)
class ParamEstimate:
    name: str
    value: float
    lo: float  # profile-likelihood interval; equals the bound where at_bound
    hi: float
    verdict: str
    lo_at_bound: bool
    hi_at_bound: bool
    wald_stderr: float
    partner: str | None = None  # the parameter it is degenerate with, if any

    def __post_init__(self):
        if self.verdict not in VERDICTS:
            raise ValueError(f"verdict {self.verdict!r} not in {VERDICTS}")

    @property
    def decades(self) -> float:
        """Width of the interval in decades, ``log10(hi/lo)``; inf if ``lo <= 0``.

        "measurable" says the interval CLOSES, not that it is narrow: on a curve with
        no yield plateau the fitter has returned a closed tau_y interval spanning
        [0.002, 9] Pa (3.7 decades) around a truth of 1 Pa. Quote the width with the
        verdict.
        """
        if self.lo <= 0 or not np.isfinite(self.hi):
            return float("inf")
        return float(np.log10(self.hi / self.lo))


@dataclass(frozen=True)
class IdentifiedModel:
    model: RateModel
    params: dict[str, ParamEstimate]
    fit: FitResult
    envelope: Envelope
    level: float

    @property
    def values(self) -> np.ndarray:
        return np.array([self.params[p].value for p in self.model.param_names])

    @property
    def verdicts(self) -> dict[str, str]:
        return {k: v.verdict for k, v in self.params.items()}

    def identifiability_block(self) -> dict[str, str]:
        """What the registry stores under ``identifiability``; nothing else writes it."""
        return {k: (v.verdict if v.partner is None else f"{v.verdict} (with {v.partner})")
                for k, v in self.params.items()}


def verdict_for(profile: Profile, corr_row: np.ndarray, index: int,
                param_names: tuple[str, ...]) -> tuple[str, str | None]:
    """Apply the rules in the module docstring to one parameter's profile."""
    open_side = profile.lower_at_bound or profile.upper_at_bound
    if open_side:
        others = [(abs(float(corr_row[j])), param_names[j])
                  for j in range(len(param_names)) if j != index]
        if others:
            c, partner = max(others)
            if np.isfinite(c) and c > DEGENERATE_CORR:
                return "degenerate", partner
    if profile.lower_at_bound and profile.upper_at_bound:
        return "not_identifiable", None
    if open_side:
        return "bounded_only", None
    return "measurable", None


def fit_rate_curve(model: RateModel, gammadot, eta, *, level: float = 0.95,
                   seed: int = 0, n_starts: int = 8, temperature_C: float = float("nan"),
                   geometry: str = "unspecified") -> IdentifiedModel:
    """Fit ``model`` to a measured ``(gammadot, eta)`` curve and profile every parameter.

    The envelope is the measured rate range of the data passed in; ``eta_at`` on a
    record built from this result will refuse to extrapolate beyond it.
    """
    gd = np.asarray(gammadot, dtype=float)
    et = np.asarray(eta, dtype=float)
    if gd.shape != et.shape or gd.ndim != 1:
        raise ValueError(f"gammadot and eta must be equal-length 1-d, got {gd.shape} and {et.shape}")
    good = np.isfinite(gd) & np.isfinite(et) & (gd > 0) & (et > 0)
    n_dropped = int((~good).sum())
    if n_dropped:
        print(f"fit_rate_curve: dropping {n_dropped} of {gd.size} points (non-finite or <= 0)")
    gd, et = gd[good], et[good]
    if gd.size <= model.n_params:
        raise ValueError(f"{gd.size} points cannot determine {model.n_params} parameters")

    # dst_fit's fitter is generic in the independent variable; its own guess is keyed on
    # model name and is not this model's, so the guess is passed explicitly.
    fit = fit_model(model, gd, et, n_starts=n_starts, seed=seed, start=model.guess(gd, et))
    if not fit.success:
        raise RuntimeError(f"fit of {model.name} failed: {fit.message}")

    params: dict[str, ParamEstimate] = {}
    for i, name in enumerate(model.param_names):
        prof = profile_likelihood(model, gd, et, name, level=level, fit=fit, seed=seed)
        verdict, partner = verdict_for(prof, fit.corr[i], i, model.param_names)
        params[name] = ParamEstimate(
            name=name, value=float(fit.params[i]), lo=float(prof.lower), hi=float(prof.upper),
            verdict=verdict, lo_at_bound=bool(prof.lower_at_bound),
            hi_at_bound=bool(prof.upper_at_bound), wald_stderr=float(prof.wald_stderr),
            partner=partner,
        )
    envelope = Envelope(gammadot_min=float(gd.min()), gammadot_max=float(gd.max()),
                        temperature_C=temperature_C, geometry=geometry)
    return IdentifiedModel(model=model, params=params, fit=fit, envelope=envelope, level=level)
