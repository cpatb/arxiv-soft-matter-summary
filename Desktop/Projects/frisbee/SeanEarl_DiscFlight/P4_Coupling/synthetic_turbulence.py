"""
Synthetic-turbulence gust generation for the coupled trajectory model.

WHY NOT WHITE NOISE
-------------------
The obvious first pass -- draw a fresh Gaussian fluctuation each time the wind
is queried -- is not usable with an adaptive ODE integrator, and the reason is
worth recording because it is easy to reintroduce.

An adaptive solver evaluates the right-hand side several times per step at
different stage points and compares the results to estimate local truncation
error. If the RHS returns a fresh random number on every call, that comparison
measures the noise rather than the discretisation error, so the controller
rejects the step, halves it, and tries again -- forever. Measured cost of one
throw with a white-noise gust term at rtol=1e-4: >100 s, versus ~0.1 s without.
A 10^4-throw Monte Carlo sweep would never finish. Worse, the resulting
"solution" is not a solution of any well-posed equation: it is grid-dependent
noise, so the landing scatter it produces is an artefact of the step size.

WHAT WE DO INSTEAD
------------------
The gust field must be a deterministic function of position and time for the
duration of one throw. We use

    u'(r, t) = sigma(r) * xi(t)

where `xi(t)` is a pre-generated, unit-variance, three-component
Ornstein-Uhlenbeck process sampled on a fixed time grid and linearly
interpolated, and `sigma(r) = sqrt(2 k(r) / 3)` is the local RMS fluctuation
from the RANS turbulent kinetic energy (isotropic assumption).

The OU process has autocorrelation exp(-|dt| / T_L), i.e. an exponential
correlation with integral time scale T_L. That is the standard first-order
approximation to a von Karman / Mann spectrum: it gets the variance and the
integral time scale right and rolls off as f^-2 at high frequency rather than
the true f^-5/3. For landing-point scatter, which is dominated by the
low-frequency, large-eddy content, this is the right leading-order model, and
unlike white noise it is a well-posed forcing.

T_L is taken from the RANS as k/epsilon, the standard Lagrangian time-scale
estimate. For the site fields here that is O(1 s), comparable to the flight
time, which is precisely the regime in which temporal correlation matters and
white noise is worst: a gust that persists for a large fraction of the flight
displaces the disc systematically, whereas white noise averages out.

LIMITATIONS (state these in the paper)
--------------------------------------
  * Frozen in space: xi(t) depends on time only, so the model imposes a gust
    that is spatially uniform over the disc's path at any instant. It captures
    temporal gust coherence but not the spatial decorrelation of eddies smaller
    than the trajectory. A full Mann (1998) box would be needed for that.
  * Isotropic: sigma is the same in all three components, since RANS k carries
    no anisotropy information. Real ABL turbulence has sigma_u > sigma_v >
    sigma_w.
  * Gaussian: no gust intermittency or extreme-value tail, so peak gusts are
    under-represented. Combined with steady RANS, results should be read as
    mean-condition playability, not worst-case.
"""

from __future__ import annotations

import numpy as np

# Bounds on the Lagrangian integral time scale inferred from k/epsilon. The
# lower bound keeps the OU process from degenerating toward white noise (and
# reintroducing the stiffness problem above) where epsilon is large; the upper
# bound stops a near-laminar cell from producing an essentially constant offset
# over the whole flight.
T_L_MIN = 0.05      # s
T_L_MAX = 20.0      # s

# Sampling interval of the pre-generated gust series. The disc travels ~0.3 m
# in this time at release speed, well below the CFD grid spacing, and it
# resolves the shortest correlation times we admit.
GUST_DT = 0.02      # s


class CorrelatedGusts:
    """Smooth band-limited gusts with a prescribed integral time scale.

    Uses the Shinozuka-Deodatis spectral representation: the fluctuation is a
    sum of N cosines with frequencies spanning the resolved band, amplitudes
    set by the target spectrum, and independent uniform random phases,

        xi(t) = sum_j A_j cos(2 pi f_j t + phi_j),   A_j = sqrt(2 S(f_j) df)

    normalised to unit variance. The target is the exponential-correlation
    (Ornstein-Uhlenbeck) spectrum

        S(f) = 4 sigma^2 T_L / (1 + (2 pi f T_L)^2)

    whose autocorrelation is exp(-|dt| / T_L).

    This is preferred over integrating an OU difference equation and
    interpolating it. A sampled-and-interpolated series is only piecewise
    linear, so the RHS has a slope discontinuity at every node; the adaptive
    integrator detects those kinks and cuts its step, which measured ~5x slower
    than the smooth form here. A cosine sum is analytic, so the solver steps
    freely, and it is the same construction family as the Mann/von Karman
    sampler the method promises -- upgrading amounts to changing `_spectrum`.
    """

    # Number of spectral modes. 32 is enough to avoid audible periodicity over
    # a flight while keeping evaluation cheap; the series repeats only on
    # 1/df = 2 * t_max, longer than any throw.
    N_MODES = 32

    def __init__(self, field, rng, t_max=8.0, T_L=None,
                 ref_pos=(0.0, 0.0, 2.0)):
        self._field = field

        if T_L is None:
            k_ref = max(float(field.k(np.asarray(ref_pos), 0.0)), 1e-12)
            eps_ref = max(float(field.eps(np.asarray(ref_pos), 0.0)), 1e-12)
            T_L = k_ref / eps_ref
        self.T_L = float(np.clip(T_L, T_L_MIN, T_L_MAX))

        # Frequencies are LOG-spaced over a wide band. Linear spacing starting
        # at 1/(2 t_max) silently discards the sub-band [0, f_min), which for
        # this spectrum holds ~24% of the variance; the missing low-frequency
        # content made the autocorrelation decay too fast and go negative
        # (a band-truncation artefact) instead of relaxing as exp(-dt/T_L).
        # Log spacing reaches two decades below the spectral corner at modest
        # mode count. Very low modes act as a slowly-varying mean-wind offset
        # over a flight, which is the physically correct behaviour.
        f_corner = 1.0 / (2.0 * np.pi * self.T_L)
        f_min = 0.01 * f_corner
        f_max = 20.0 * f_corner
        f = np.logspace(np.log10(f_min), np.log10(f_max), self.N_MODES)
        df = np.gradient(f)

        S = 4.0 * self.T_L / (1.0 + (2.0 * np.pi * f * self.T_L) ** 2)
        A = np.sqrt(2.0 * S * df)
        # Normalise so the process has exactly unit stationary variance: for a
        # sum of cosines with random phase, var = sum(A_j^2) / 2.
        A = A / np.sqrt(np.sum(A * A) / 2.0)

        self._omega = 2.0 * np.pi * f                      # (N,)
        self._A = A                                        # (N,)
        self._phase = rng.uniform(0.0, 2.0 * np.pi, size=(3, self.N_MODES))

    def _spectrum_check_variance(self, n=200000, t_max=8.0, rng=None):
        """Sample variance of xi over a long record; should be ~1."""
        rng = rng or np.random.default_rng(0)
        t = np.linspace(0.0, 200.0, n)
        xi = (self._A * np.cos(self._omega * t[:, None]
                               + self._phase[0])).sum(axis=1)
        return float(xi.var())

    def autocorrelation(self, lags):
        """Analytic autocorrelation of the realised mode set, for testing."""
        lags = np.atleast_1d(lags)
        w = self._A ** 2 / 2.0
        return (w * np.cos(self._omega * lags[:, None])).sum(axis=1) / w.sum()

    def shape(self, t):
        """Unit-variance three-component fluctuation at time t."""
        # (3, N) -> (3,)
        return (self._A * np.cos(self._omega * t + self._phase)).sum(axis=1)

    def fluctuation(self, k_local, t):
        """Gust given an ALREADY-INTERPOLATED local TKE."""
        return np.sqrt(2.0 * max(float(k_local), 0.0) / 3.0) * self.shape(t)

    def __call__(self, r, t):
        return self.fluctuation(self._field.k(r, t), t)


class _OUGustsReference:
    """Ornstein-Uhlenbeck gusts by direct time-stepping, kept as a reference.

    Same statistics as `CorrelatedGusts` but built by integrating the OU
    difference equation and interpolating, which is easier to verify by eye.
    Used in the test suite to confirm the spectral construction reproduces the
    intended autocorrelation. Not used for production runs -- see the note in
    `CorrelatedGusts` about interpolation kinks slowing the integrator.

    One instance corresponds to ONE throw: the time series is drawn in the
    constructor, so repeated calls at the same `t` return the same value and
    the integrator sees a smooth, deterministic forcing.

    Parameters
    ----------
    field : object with `.k(r, t)` and `.eps(r, t)`
    rng : np.random.Generator
    t_max : float
        Length of the series to pre-generate, seconds.
    T_L : float, optional
        Integral time scale. If None it is taken from k/epsilon at `ref_pos`.
    ref_pos : array_like, optional
        Position at which to evaluate k/epsilon for the time scale.
    """

    def __init__(self, field, rng, t_max=8.0, T_L=None,
                 ref_pos=(0.0, 0.0, 2.0)):
        self._field = field

        if T_L is None:
            k_ref = max(float(field.k(np.asarray(ref_pos), 0.0)), 1e-12)
            eps_ref = max(float(field.eps(np.asarray(ref_pos), 0.0)), 1e-12)
            T_L = k_ref / eps_ref
        self.T_L = float(np.clip(T_L, T_L_MIN, T_L_MAX))

        # Exact discrete update for an OU process with unit stationary
        # variance: xi_{n+1} = a xi_n + sqrt(1 - a^2) z, a = exp(-dt / T_L).
        # Using the exact solution (rather than an Euler-Maruyama step) keeps
        # the stationary variance exactly 1 for any dt.
        n = int(np.ceil(t_max / GUST_DT)) + 2
        a = np.exp(-GUST_DT / self.T_L)
        b = np.sqrt(max(1.0 - a * a, 0.0))

        xi = np.empty((n, 3))
        xi[0] = rng.standard_normal(3)          # start from the stationary pdf
        noise = rng.standard_normal((n, 3))
        for i in range(1, n):
            xi[i] = a * xi[i - 1] + b * noise[i]

        self._t = np.arange(n) * GUST_DT
        self._xi = xi

    def __call__(self, r, t):
        # Local RMS fluctuation from the RANS TKE, isotropic.
        k_loc = max(float(self._field.k(r, t)), 0.0)
        sigma = np.sqrt(2.0 * k_loc / 3.0)
        xi = np.array([
            np.interp(t, self._t, self._xi[:, 0]),
            np.interp(t, self._t, self._xi[:, 1]),
            np.interp(t, self._t, self._xi[:, 2]),
        ])
        return sigma * xi


class WhiteNoiseGusts:
    """DEPRECATED -- kept only so the failure mode stays documented and tested.

    Draws an independent fluctuation on every call. See the module docstring:
    this makes the ODE right-hand side non-deterministic, which defeats
    adaptive step-size control (>1000x slowdown) and yields a step-size-
    dependent, non-convergent "solution". Do not use it for results.
    """

    def __init__(self, field, rng: np.random.Generator):
        self._field = field
        self._rng = rng

    def __call__(self, r, t):
        k_loc = max(self._field.k(r, t), 0.0)
        sigma = np.sqrt(2.0 * k_loc / 3.0)
        return self._rng.normal(0.0, sigma, size=3)


def mann_gusts(field, rng, L_int=20.0, dt=0.01, T=8.0):
    """Stub for a full Mann (1998) three-dimensional Fourier gust box.

    Would replace `CorrelatedGusts` for production runs: generates a frozen 3-D
    turbulence field with the correct spectral tensor, so that gusts decorrelate
    in SPACE as well as time. Deferred; see the limitations note in the module
    docstring.
    """
    raise NotImplementedError(
        "Mann-spectrum gust box not implemented; use CorrelatedGusts."
    )


def wind_with_gusts(field, gust_gen):
    """Compose a `wind(r, t)` callable: RANS mean plus gust fluctuation.

    Takes the fused `U_and_k` path when the field offers one, so that the mean
    and the gust scaling share a single interpolation per evaluation.
    """
    fused = getattr(field, "U_and_k", None)
    frac = getattr(gust_gen, "fluctuation", None)

    if fused is not None and frac is not None:
        def _w(r, t):
            U, k_loc = fused(r, t)
            return U + frac(k_loc, t)
        return _w

    def _w(r, t):
        return np.asarray(field.U(r, t)) + np.asarray(gust_gen(r, t))
    return _w
