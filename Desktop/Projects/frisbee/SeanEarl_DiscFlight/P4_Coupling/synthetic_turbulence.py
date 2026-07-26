"""
Simple synthetic-turbulence sampler.

For the first-pass Monte Carlo we use a **white-in-time** Gaussian fluctuation
with magnitude sigma_i = sqrt(2*k/3) per component (isotropic assumption).
This captures the local TKE-driven scatter of landing points but does NOT
model temporal/spatial gust correlation. Replace with a Mann (1998)-based
Fourier sampler for production runs (stub in `mann_gusts` below).
"""

from __future__ import annotations
import numpy as np


class WhiteNoiseGusts:
    """Zero-time-correlation Gaussian gusts scaled by local TKE."""

    def __init__(self, field, rng: np.random.Generator):
        self._field = field
        self._rng = rng

    def __call__(self, r, t):
        k_loc = max(self._field.k(r, t), 0.0)
        sigma = np.sqrt(2.0 * k_loc / 3.0)
        return self._rng.normal(0.0, sigma, size=3)


def mann_gusts(field, rng, L_int=20.0, dt=0.01, T=8.0):
    """Stub for a Mann-1998 / von Kármán Fourier-sampler gust generator.

    For each Monte Carlo throw, pre-generate a time series u'(t) of length
    T/dt with spectrum matching local k and integral length L_int, then
    interpolate during integration. Implementation deferred to v2.
    """
    raise NotImplementedError(
        "Mann-spectrum gust generator not implemented; use WhiteNoiseGusts."
    )


def wind_with_gusts(field, gust_gen):
    """Compose a `wind(r, t)` callable: mean from CFD + fluctuation from generator."""
    def _w(r, t):
        return np.asarray(field.U(r, t)) + np.asarray(gust_gen(r, t))
    return _w
