"""Synthetic flow curves with known truth, and the recovery study that gates the fitter.

A round trip ``fit(synthesise(truth)) == truth`` certifies the optimiser, not the
physics; here it is used for exactly that, to check that the *verdicts* and *intervals*
of :mod:`biorheo.identify.fit` mean what they claim. Noise is lognormal (multiplicative,
``eta_meas = eta exp(cv * z)``), matching the log-viscosity residual that ``dst_fit``
minimises, so the interval's nominal coverage is the thing under test.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..constitutive import RateModel
from .fit import IdentifiedModel, fit_rate_curve

__all__ = ["synthetic_rate_curve", "RecoveryReport", "recovery_study"]


def synthetic_rate_curve(model: RateModel, truth, gammadot, noise_cv: float,
                         seed: int) -> np.ndarray:
    """``eta(gammadot; truth) * exp(noise_cv * N(0,1))`` with a recorded seed."""
    gd = np.asarray(gammadot, dtype=float)
    truth = np.asarray(truth, dtype=float)
    model.check_params(truth)
    rng = np.random.default_rng(seed)
    eta = model(gd, *truth)
    if not np.all(np.isfinite(eta) & (eta > 0)):
        raise ValueError("truth gives a non-finite or non-positive viscosity on this grid")
    return eta * np.exp(noise_cv * rng.standard_normal(gd.shape))


@dataclass(frozen=True)
class RecoveryReport:
    model_name: str
    truth: dict[str, float]
    n_seeds: int
    within_interval: dict[str, int]  # seeds whose profile interval contains the truth (any interval)
    n_closed: dict[str, int]         # seeds whose interval closed on both sides
    within_closed: dict[str, int]    # of those, how many contain the truth: the coverage that can fail
    verdicts: dict[str, dict[str, int]]  # param -> verdict -> count
    fits: tuple[IdentifiedModel, ...]

    def summary(self) -> str:
        lines = [f"{self.model_name}: {self.n_seeds} seeds"]
        for p in self.truth:
            v = ", ".join(f"{k} {n}" for k, n in sorted(self.verdicts[p].items()))
            lines.append(f"  {p:>8s} = {self.truth[p]:<10.4g} in interval {self.within_interval[p]:2d}/{self.n_seeds}"
                         f"  closed {self.within_closed[p]}/{self.n_closed[p]}  ({v})")
        return "\n".join(lines)


def recovery_study(model: RateModel, truth, gammadot, *, n_seeds: int = 20,
                   noise_cv: float = 0.03, level: float = 0.95) -> RecoveryReport:
    """Fit ``n_seeds`` synthetic curves and count, per parameter, how often the profile
    interval contains the truth and which verdict was returned."""
    truth = np.asarray(truth, dtype=float)
    within = {p: 0 for p in model.param_names}
    n_closed = {p: 0 for p in model.param_names}
    within_closed = {p: 0 for p in model.param_names}
    verdicts: dict[str, dict[str, int]] = {p: {} for p in model.param_names}
    fits = []
    for seed in range(n_seeds):
        eta = synthetic_rate_curve(model, truth, gammadot, noise_cv, seed)
        ident = fit_rate_curve(model, gammadot, eta, level=level, seed=seed)
        fits.append(ident)
        for i, p in enumerate(model.param_names):
            est = ident.params[p]
            if est.lo <= truth[i] <= est.hi:
                within[p] += 1
            # an interval that ran to a bound contains the truth for free (tau_y's lower
            # bound is 0); only closed intervals make a coverage count that can fail
            if not (est.lo_at_bound or est.hi_at_bound):
                n_closed[p] += 1
                if est.lo <= truth[i] <= est.hi:
                    within_closed[p] += 1
            verdicts[p][est.verdict] = verdicts[p].get(est.verdict, 0) + 1
    return RecoveryReport(model_name=model.name, truth=dict(zip(model.param_names, truth.tolist())),
                          n_seeds=n_seeds, within_interval=within, n_closed=n_closed,
                          within_closed=within_closed, verdicts=verdicts,
                          fits=tuple(fits))
