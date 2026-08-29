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

The process has autocorrelation exp(-|ds| / L), i.e. an exponential correlation
with integral scale L. That is the standard first-order approximation to a von
Karman / Mann spectrum: it gets the variance and the integral scale right and
rolls off as f^-2 at high frequency rather than the true f^-5/3. For
landing-point scatter, which is dominated by the low-frequency, large-eddy
content, this is the right leading-order model, and unlike white noise it is a
well-posed forcing.

WHICH COORDINATE THE GUST DECORRELATES IN  (this is the load-bearing choice)
---------------------------------------------------------------------------
The disc supplies its own advection. Its airspeed is V_d ~ 14 m/s against a
site mean wind of ~1.2-1.4 m/s and sigma_u ~ 0.3-0.6 m/s, so V_d / sigma_u =
17-36: an eddy barely evolves while the disc traverses it. The field is
therefore FROZEN, and the disc decorrelates by flying out of one eddy and into
the next -- that is, over PATH ARCLENGTH, not over wall-clock time.

    tau_c = L_int / V_d  ~  6 m / 14 m/s  ~  0.43 s

The `temporal` model retained below instead sets the correlation time to the
RANS turbulence decay time T_L = k/eps, evaluated once at the launch point.
Measured in the disc band at this site that is 10.5-18.0 s, against a flight
time of ~2.5 s. A gust correlated over 10-18 s is not turbulence at all: it is
a per-throw constant. Since the wind enters only through V_rel = v - U and
Galilean invariance is exact in this code, a spatially uniform quasi-constant
u' IS a randomly drawn mean wind -- which is precisely the quantity the study
claims a thrower CAN aim off. Ratio of the two time scales: ~30x.

`sweeping` (the default) makes xi a function of arclength s = |r - r_launch|
with decorrelation length L_int, and is a deterministic function of position,
so it remains a well-posed forcing for the adaptive integrator. `temporal` is
kept, bit-reproducible, so the two can be compared on identical release draws.

ANISOTROPY
----------
RANS k carries no anisotropy information, so the split is imposed: sigma_v /
sigma_u = 0.75 and sigma_w / sigma_u = 0.5, the standard surface-layer values.
The ratios are applied at FIXED k, i.e. sigma_u is raised so that
k = (sigma_u^2 + sigma_v^2 + sigma_w^2) / 2 is preserved; k is what the RANS
actually provides and it must not be silently changed. The u axis is the local
mean-wind direction and w is vertical. This matters because w' drives the angle
of attack and hence the pitching moment, so isotropy inflates the moment
channel by ~2x in variance.

LIMITATIONS (state these in the paper)
--------------------------------------
  * The frozen field is one-dimensional: xi varies along the path but is
    uniform across it, so the model has no spanwise decorrelation and no
    rolling moment from a gust gradient across the disc. A full Mann (1998)
    box or a 3-D random-Fourier-mode field would be needed for that, and is
    the natural next upgrade -- it amounts to giving each mode a wavevector
    instead of a scalar wavenumber.
  * L_int is evaluated once per throw at the launch point and held fixed along
    the path. sigma IS updated along the path (it is interpolated from the
    local k at every RHS evaluation), so the height dependence of the
    fluctuation amplitude is represented, but the height dependence of the
    integral SCALE is not. Making L position-dependent requires a path
    integral of ds / L(r(s)) and therefore an extra ODE state; deferred.
  * The integral scale convention is not unique. This module uses
    L_int = sigma_u * k / eps = sqrt(2/3) k^1.5 / eps, which is exactly 4.97x
    the standard RANS mixing length C_mu^0.75 k^1.5 / eps. Report which one you
    mean; the resulting N = L_f / L_int differs by the same factor.
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

# Bounds on the integral LENGTH scale, the sweeping model's analogue of
# T_L_MIN/T_L_MAX. The lower bound keeps the spatial process from degenerating
# toward white noise along the path; the upper bound stops a near-laminar cell
# from producing a constant offset over the whole flight. At V_d ~ 14 m/s these
# correspond to correlation times of ~0.02 s and ~4.3 s.
L_INT_MIN = 0.25    # m
L_INT_MAX = 60.0    # m

# Surface-layer anisotropy, as (sigma_u, sigma_v, sigma_w) RATIOS. Applied at
# fixed k, so sigma_u is scaled up to compensate; see the module docstring.
# u is along the local mean wind, w is vertical.
SIGMA_RATIOS = (1.0, 0.75, 0.5)
ISOTROPIC_RATIOS = (1.0, 1.0, 1.0)

# Below this horizontal mean speed the wind direction carries no information
# and the turbulence axes are taken as (x, y, z) instead. The site contains
# cells at |U| = 0.02 m/s, where the interpolated direction is numerical noise.
U_DIR_MIN = 0.05    # m/s


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
    #
    # This count is also a load-bearing invariant of the PAIRED design: the
    # constructor draws exactly 3 * N_MODES uniforms whatever the model, field
    # or scale, so the gust and no-gust arms advance the RNG identically and
    # genuinely share release draws. Do not make it adaptive without revisiting
    # `test_paired_arms_consume_identical_rng_state`.
    N_MODES = 32

    def __init__(self, field, rng, t_max=8.0, T_L=None,
                 ref_pos=(0.0, 0.0, 2.0), model="sweeping",
                 L_int=None, sigma_ratios=SIGMA_RATIOS):
        self._field = field
        if model not in ("sweeping", "temporal"):
            raise ValueError(f"unknown gust model {model!r}; "
                             "expected 'sweeping' or 'temporal'")
        self.model = model
        ref = np.asarray(ref_pos, dtype=float)
        self._ref_pos = ref

        k_ref = max(float(field.k(ref, 0.0)), 1e-12)
        eps_ref = max(float(field.eps(ref, 0.0)), 1e-12)

        # Both models are the same spectral construction; they differ only in
        # WHICH coordinate the process decorrelates over, and hence in the
        # scale. `scale` is a time (s) for `temporal` and a length (m) for
        # `sweeping`, and `_omega` is correspondingly rad/s or rad/m.
        if model == "temporal":
            raw = k_ref / eps_ref if T_L is None else float(T_L)
            scale = float(np.clip(raw, T_L_MIN, T_L_MAX))
            lo, hi = T_L_MIN, T_L_MAX
            self.T_L = scale
            self.L_int = float("nan")
        else:
            # L_int = sigma_u * k / eps = sqrt(2/3) k^1.5 / eps.
            # See the module docstring: this is 4.97x the RANS mixing length,
            # and which convention is meant must be stated in the paper.
            raw = (np.sqrt(2.0 * k_ref / 3.0) * k_ref / eps_ref
                   if L_int is None else float(L_int))
            scale = float(np.clip(raw, L_INT_MIN, L_INT_MAX))
            lo, hi = L_INT_MIN, L_INT_MAX
            self.L_int = scale
            self.T_L = float("nan")

        self.scale = scale
        self.scale_raw = float(raw)
        # Recorded per throw so a run can be audited for silent clipping, which
        # previously happened at T_L_MAX with nothing written down.
        self.clipped = bool(raw < lo or raw > hi)

        # Modes are LOG-spaced over a wide band. Linear spacing starting at
        # 1/(2 t_max) silently discards the sub-band [0, f_min), which for this
        # spectrum holds ~24% of the variance; the missing low-frequency
        # content made the autocorrelation decay too fast and go negative
        # (a band-truncation artefact) instead of relaxing as exp(-ds/L).
        # Log spacing reaches two decades below the spectral corner at modest
        # mode count. Very low modes act as a slowly-varying offset over a
        # flight, which is the physically correct behaviour.
        f_corner = 1.0 / (2.0 * np.pi * scale)
        f = np.logspace(np.log10(0.01 * f_corner),
                        np.log10(20.0 * f_corner), self.N_MODES)
        df = np.gradient(f)

        S = 4.0 * scale / (1.0 + (2.0 * np.pi * f * scale) ** 2)
        A = np.sqrt(2.0 * S * df)
        # Normalise so the process has exactly unit stationary variance: for a
        # sum of cosines with random phase, var = sum(A_j^2) / 2.
        A = A / np.sqrt(np.sum(A * A) / 2.0)

        self._omega = 2.0 * np.pi * f                      # (N,)
        self._A = A                                        # (N,)
        self._phase = rng.uniform(0.0, 2.0 * np.pi, size=(3, self.N_MODES))

        # --- anisotropy -------------------------------------------------
        r_u, r_v, r_w = (float(x) for x in sigma_ratios)
        # Preserve k: k = (su^2 + sv^2 + sw^2)/2, so su = sqrt(2k / sum(r^2)).
        # The isotropic case r=(1,1,1) recovers su = sv = sw = sqrt(2k/3).
        self._ratios = np.array([r_u, r_v, r_w])
        self._sigma_gain = self._ratios / np.sqrt(np.sum(self._ratios ** 2))
        # sigma_i = _sigma_gain[i] * sqrt(2 k)

        # Turbulence axes, fixed for the throw. The u axis is the mean-wind
        # direction at the launch point, NOT the local one: a frame that
        # re-aligns with the local flow would spin discontinuously through the
        # near-calm cells this site contains, and surface-layer anisotropy is
        # set by the ground and the large-scale wind, not by the local eddy.
        U_ref = np.asarray(field.U(ref, 0.0), dtype=float)
        e_u = np.array([U_ref[0], U_ref[1], 0.0])
        n_u = float(np.linalg.norm(e_u))
        if n_u < U_DIR_MIN:
            e_u = np.array([1.0, 0.0, 0.0])
        else:
            e_u = e_u / n_u
        e_w = np.array([0.0, 0.0, 1.0])
        e_v = np.cross(e_w, e_u)
        self._axes = np.stack([e_u, e_v, e_w])              # (3, 3) rows

    def _spectrum_check_variance(self, n=200000, t_max=8.0, rng=None):
        """Sample variance of xi over a long record; should be ~1."""
        rng = rng or np.random.default_rng(0)
        t = np.linspace(0.0, 200.0, n)
        xi = (self._A * np.cos(self._omega * t[:, None]
                               + self._phase[0])).sum(axis=1)
        return float(xi.var())

    def autocorrelation(self, lags):
        """Analytic autocorrelation of the realised mode set, for testing.

        `lags` are seconds for the `temporal` model and metres for `sweeping`.
        """
        lags = np.atleast_1d(lags)
        w = self._A ** 2 / 2.0
        return (w * np.cos(self._omega * lags[:, None])).sum(axis=1) / w.sum()

    def coordinate(self, r, t):
        """The coordinate this realisation decorrelates over.

        `sweeping`: path arclength from the launch point, approximated by the
        chord |r - r_launch|. Exact for a straight path, and a disc trajectory
        is straight to within a few percent in plan; more importantly it is a
        smooth analytic function of position, so the frozen field stays a
        well-posed forcing for the adaptive integrator without adding an ODE
        state to carry true arclength.

        `temporal`: wall-clock time, the pre-existing behaviour.
        """
        if self.model == "temporal":
            return float(t)
        return float(np.linalg.norm(np.asarray(r, dtype=float) - self._ref_pos))

    def shape(self, coord):
        """Unit-variance 3-component fluctuation in TURBULENCE axes (u, v, w)."""
        # (3, N) -> (3,)
        return (self._A * np.cos(self._omega * coord + self._phase)).sum(axis=1)

    def fluctuation(self, k_local, r, t):
        """Gust in WORLD axes, given an ALREADY-INTERPOLATED local TKE.

        sigma is evaluated from the LOCAL k at every call, so the height (and
        horizontal) dependence of the fluctuation amplitude is carried along
        the path. Only the integral scale is frozen at the launch value.
        """
        xi = self.shape(self.coordinate(r, t))
        sigma = self._sigma_gain * np.sqrt(2.0 * max(float(k_local), 0.0))
        return (sigma * xi) @ self._axes

    def __call__(self, r, t):
        return self.fluctuation(self._field.k(r, t), r, t)


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
            return U + frac(k_loc, r, t)
        return _w

    def _w(r, t):
        return np.asarray(field.U(r, t)) + np.asarray(gust_gen(r, t))
    return _w


def wind_mean_only(field, gust_gen):
    """`wind(r, t)` using the RANS mean ONLY, with the gusts suppressed.

    The counterfactual half of a paired experiment: identical site, identical
    mean flow, identical release-parameter draws, zero fluctuation. Comparing
    landing scatter against `wind_with_gusts` isolates the part of the scatter
    caused by turbulence from the part caused by the mean flow and by flight
    time alone -- a disc released higher flies longer and would disperse more
    even in perfectly still air.

    `gust_gen` is accepted and deliberately ignored. The caller must still
    CONSTRUCT it, because building a gust realization consumes random numbers;
    dropping that draw would desynchronise the RNG stream and the two arms
    would no longer share release conditions, destroying the pairing.
    """
    del gust_gen

    fused = getattr(field, "U_and_k", None)
    if fused is not None:
        def _w(r, t):
            U, _ = fused(r, t)
            return U
        return _w

    def _w(r, t):
        return np.asarray(field.U(r, t))
    return _w
