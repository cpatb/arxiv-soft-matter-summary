"""The strain clock: how fast the frictional network rebuilds after a reversal.

Why this module exists
----------------------
Everything else in this package fits *steady* flow curves. Wyart-Cates infers a
frictional fraction ``f(sigma)`` from the shape of ``eta(sigma)`` and never
watches it move. But the closure that a transient solver actually integrates is
not the steady relation at all; it is a relaxation law. Richards et al. (2019),
Eq. 4, write it as

    df/dt = -(gamma_dot / gamma_0 + 1 / t_r) (f - f_ss)

with a *strain* clock ``gamma_0`` and a *time* clock ``t_r``. Those two constants
set the whole transient response, and neither is measurable from a steady flow
curve: at steady state the left-hand side is zero and both clocks drop out.

A shear-reversal experiment measures them. Reversing the flow destroys the
frictional contact network, and the instantaneous viscosity then climbs back to
its steady value over accumulated strain. Lin et al. (2015) deposited exactly
that: instantaneous relative viscosity against accumulated strain, at four rates,
for two materials. This module turns those records into ``gamma_0``.

The inversion, and why it is the identifiable half of the model
---------------------------------------------------------------
At constant shear rate, strain is the natural clock and the relaxation law
integrates in closed form. Writing ``gamma = gamma_dot t``,

    df/dgamma = -(1/gamma_0 + 1/(gamma_dot t_r)) (f - f_ss)

so ``f`` approaches ``f_ss`` exponentially in strain with a single constant

    1/gamma_eff = 1/gamma_0 + 1/(gamma_dot t_r).                          (*)

``gamma_eff`` is therefore measurable one rate at a time, and the two clocks are
separated by how ``gamma_eff`` varies across rates. Note the sign that (*)
predicts: ``gamma_eff`` *rises* with shear rate towards ``gamma_0`` from below,
because the time term matters less the faster you shear.

Getting ``f`` from the measured viscosity means inverting Wyart-Cates. That
inversion is analytic, and the useful part is what it depends on::

    eta_r = (1 - phi/phi_J)^-n,     phi_J = (1 - f) phi_0 + f phi_m

    =>  phi/phi_J = 1 - eta_r^(-1/n) =: u

    =>  f = (phi_0/phi - 1/u) / (phi_0/phi - phi_m/phi)

**Only the ratios ``phi_0/phi`` and ``phi_m/phi`` appear.** That matters more
than it looks. The central result of this package is that ``phi`` and ``phi_m``
enter Wyart-Cates only through ``phi/phi_J``, so a single flow curve determines
that ratio and neither number; the fitted ``phi`` and ``phi_m`` slide together
along a degenerate valley. The inversion above is a function of exactly the
ratios the valley preserves, so it is *invariant along the degenerate direction*.
``f``, and therefore ``gamma_0``, inherits the identifiability of ``sigma*``
rather than the degeneracy of ``phi_m``. :func:`frictional_fraction_from_viscosity`
is where that is used and ``tests/test_transient.py`` is where it is proved.

Where the inversion stops working, which is a result rather than a nuisance
--------------------------------------------------------------------------
Wyart-Cates bounds the viscosity from below at ``f = 0``: with no frictional
contacts the jamming point is the frictionless one and

    eta_r >= (1 - phi/phi_0)^-n =: eta_frictionless.

The measured post-reversal minimum sits *below* that floor on every record
tried here. No value of ``f`` reproduces it, and the inversion returns ``f < 0``
rather than pretending otherwise. This is not a fitting failure; it is the model
missing a variable. Reversal does not only break contacts, it also reverses the
microstructural anisotropy that carries hydrodynamic stress, and Wyart-Cates has
no state variable for microstructure. So the first fraction of a strain unit
after a reversal is outside the model's reach entirely, and any transient solver
running a Wyart-Cates closure will under-predict the softening there.
:func:`admissible_window` reports how much of a record that costs instead of
silently clipping it.

References
----------
Lin, Guy, Hermes, Ness, Sun, Poon & Cohen, PRL 115, 228304 (2015).
Richards, Royer, Liebchen, Guy & Poon, PRL 123, 038004 (2019), Eq. 4.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import uniform_filter1d

from dst_fit.models import WC_EXPONENT

# --- reduction of the raw records --------------------------------------------
# Log bins across the strain axis. The gate sweeps 50 to 800 and reports the
# spread, because that spread is the dominant systematic on gamma_eff.
N_BINS_DEFAULT = 200

# A bin with fewer than this many raw points has no usable scatter estimate.
MIN_BIN_COUNT = 3

# Consistency constant for a MAD-based standard deviation under normality, and
# the asymptotic efficiency penalty of a median relative to a mean, sqrt(pi/2).
MAD_TO_SIGMA = 1.4826
MEDIAN_SE_FACTOR = 1.2533141373155003

# Digitisation artefact. Several deposited records end with a handful of points
# that fall off a cliff at the right-hand frame edge (28 points at gamma_dot = 2,
# 48 at gamma_dot = 5 on Fig2b): the plot border traced as data. The rule is
# deliberately narrow so it cannot touch the interior, and the count is reported.
ARTEFACT_EDGE_FRACTION = 0.98  # only the last 2% of the strain axis is eligible
ARTEFACT_DROP_FRACTION = 0.70  # ... and only if it falls below 70% of the plateau

# --- the separable scan -------------------------------------------------------
# Log-spaced over four and a half decades of strain; 1200 points resolves lambda
# to 0.9%, comfortably inside every interval this returns.
LAMBDA_GRID = (1e-3, 5e1, 1200)

# chi-squared quantile for a 95% profile interval on one parameter.
CHI2_1DOF_95 = 3.841458820694124

# Boxcar width used to defeat single-point noise before landmarking, given as a
# *strain* window rather than a sample count. The Lin records are each sampled
# uniformly in strain but not at the same step: d_gamma runs 4.0e-4 to 1.0e-3
# across the deposit, so a fixed sample count would smooth a 2.5x different
# physical window on different records and the landmark would carry that
# difference. 0.05 strain units sits about 20x below a recovery clock of order 1
# strain unit: wide enough to kill the isolated spikes that sit at the end of
# several records, far too narrow to touch the quantity being measured.
SMOOTH_STRAIN = 0.05

# Floor on the boxcar, in samples. A window of one is no filter at all.
SMOOTH_MIN_SAMPLES = 3

# The post-reversal minimum is sought below this strain. Unrestricted, the global
# argmin of several records lands on a one-point dip at the very end of the
# record rather than on the reversal dip, which sits near gamma ~ 0.01 to 0.5.
INTERIOR_MAX_STRAIN = 2.0

# A record must retain at least this share of its post-minimum points inside the
# model's admissible range for a recovery fit to be attempted. Below it there is
# no exponential left to fit and a returned constant would be an artefact of the
# window rather than a measurement.
MIN_ADMISSIBLE_FRACTION = 0.5


def frictional_fraction_from_viscosity(
    eta_r: np.ndarray,
    phi: float,
    phi_m: float,
    phi_0: float,
    n: float = WC_EXPONENT,
) -> np.ndarray:
    """Invert Wyart-Cates for the frictional fraction ``f`` implied by ``eta_r``.

    The exact analytic inverse of :func:`dst_fit.models._wyart_cates`, in
    *relative* viscosity so that no solvent viscosity is needed:
    ``f = (phi_0/phi - 1/u) / (phi_0/phi - phi_m/phi)`` with
    ``u = 1 - eta_r^(-1/n)``.

    Values outside ``[0, 1]`` are returned as they come out, **not clipped**.
    ``f > 1`` means the measured viscosity exceeds what the model allows even
    with every contact frictional; ``f < 0`` means it sits below the frictionless
    branch. Both say the data is outside the model's range, which is information.
    Clipping would convert that into a plausible-looking number.

    ``eta_r <= 1`` returns NaN: a relative viscosity at or below the solvent's
    is not a suspension viscosity, and the inversion has no root there.

    Parameters
    ----------
    eta_r
        Relative viscosity, dimensionless.
    phi, phi_m, phi_0
        Volume fraction, frictional jamming fraction, frictionless jamming
        fraction. Only the ratios ``phi_0/phi`` and ``phi_m/phi`` affect the
        result; see the module docstring.
    n
        Divergence exponent, matching the fit that produced ``phi`` and ``phi_m``.
    """
    eta_r = np.asarray(eta_r, dtype=float)
    if not (0.0 < phi < phi_m < phi_0):
        raise ValueError(
            f"need 0 < phi < phi_m < phi_0; got phi={phi:g}, phi_m={phi_m:g}, "
            f"phi_0={phi_0:g}"
        )
    if not np.isfinite(n) or n <= 0:
        raise ValueError(f"exponent n must be finite and positive, got {n!r}")

    out = np.full(eta_r.shape, np.nan, dtype=float)
    ok = np.isfinite(eta_r) & (eta_r > 1.0)
    if not np.any(ok):
        return out

    u = 1.0 - eta_r[ok] ** (-1.0 / n)  # = phi / phi_J, strictly in (0, 1)
    r0 = phi_0 / phi
    rm = phi_m / phi
    out[ok] = (r0 - 1.0 / u) / (r0 - rm)
    return out


def frictionless_relative_viscosity(phi: float, phi_0: float, n: float = WC_EXPONENT) -> float:
    """The model's floor, ``eta_r`` at ``f = 0``. Nothing below this is invertible."""
    if not (0.0 < phi < phi_0):
        raise ValueError(f"need 0 < phi < phi_0; got phi={phi:g}, phi_0={phi_0:g}")
    return float((1.0 - phi / phi_0) ** -float(n))


def reduced_jamming_distance(eta_r: np.ndarray, n: float = WC_EXPONENT) -> np.ndarray:
    """``u = phi_J / phi``, the exact inverse of ``eta_r = (1 - phi/phi_J)^-n``.

    This is the *identified combination*. The central result of this package is
    that ``phi`` and ``phi_m`` enter Wyart-Cates only through ``phi/phi_J``, so a
    flow curve approaching jamming pins that ratio and neither number. Reading a
    transient in ``u`` therefore uses only what the steady data determines, and
    :func:`frictional_fraction_from_viscosity` is exactly ``A + B u`` on top of
    it. Defined and monotone for every ``eta_r > 1``.
    """
    eta = np.asarray(eta_r, dtype=float)
    return 1.0 / (1.0 - eta ** (-1.0 / float(n)))


def d_reduced_jamming_distance(eta_r: np.ndarray, n: float = WC_EXPONENT) -> np.ndarray:
    """``du / d eta_r``, analytic, for propagating measured scatter into ``u``.

    Analytic rather than a finite difference: ``u`` is order 1 while ``du`` over
    a plausible step is order 1e-8, so a central difference loses most of its
    significant figures to cancellation exactly where the weights matter most.
    """
    eta = np.asarray(eta_r, dtype=float)
    nn = float(n)
    s = eta ** (-1.0 / nn)
    return -(s / (nn * eta)) / (1.0 - s) ** 2


def trim_edge_artefact(strain: np.ndarray, eta_r: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    """Drop the digitisation cliff at the right-hand end of a deposited record.

    A point is removed only if it lies in the last :data:`ARTEFACT_EDGE_FRACTION`
    of the strain axis *and* falls below :data:`ARTEFACT_DROP_FRACTION` of the
    late-strain plateau. Nothing in the interior can satisfy both, so this cannot
    reshape the recovery being measured; it removes the traced plot border.
    """
    g = np.asarray(strain, dtype=float)
    e = np.asarray(eta_r, dtype=float)
    plateau = float(np.median(e[g >= 0.5 * g.max()]))
    bad = (g >= ARTEFACT_EDGE_FRACTION * g.max()) & (e < ARTEFACT_DROP_FRACTION * plateau)
    return g[~bad], e[~bad], int(bad.sum())


def log_bin(
    strain: np.ndarray, eta_r: np.ndarray, n_bins: int = N_BINS_DEFAULT
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reduce a record to log-spaced strain bins with a **measured** uncertainty.

    Returns ``(strain, eta_r, eta_err)``: the geometric-mean strain, the median
    viscosity, and the standard error of that median estimated from the
    within-bin MAD. Median and MAD rather than mean and standard deviation
    because the raw records carry single-point digitisation spikes that both
    non-robust estimators would propagate straight into the weights.

    Binning is not cosmetic. The raw records hold 10 000 to 25 000 samples that
    are strongly autocorrelated, and they differ between rates by 2.5x in count;
    fitting at native sampling therefore weights the rates against each other by
    sampling density rather than by information, and reports an interval far
    tighter than the analysis choices support.
    """
    g, e, _ = trim_edge_artefact(strain, eta_r)
    edges = np.logspace(np.log10(g.min()), np.log10(g.max()), int(n_bins) + 1)
    idx = np.clip(np.digitize(g, edges) - 1, 0, int(n_bins) - 1)

    gm, em, ee = [], [], []
    for b in range(int(n_bins)):
        m = idx == b
        k = int(m.sum())
        if k < MIN_BIN_COUNT:
            continue
        med = float(np.median(e[m]))
        sd = MAD_TO_SIGMA * float(np.median(np.abs(e[m] - med)))
        gm.append(float(np.exp(np.mean(np.log(g[m])))))
        em.append(med)
        ee.append(MEDIAN_SE_FACTOR * sd / np.sqrt(k))

    if len(gm) < 10:
        raise ValueError(
            f"log_bin: only {len(gm)} usable bins from {g.size} points at "
            f"n_bins={n_bins}; too few to fit a relaxation"
        )
    err = np.asarray(ee, dtype=float)
    # A bin whose points happen to be identical gives sd = 0. Floor it at the
    # smallest non-zero error the record itself produced, rather than at an
    # invented constant, so the weighting stays on the data's own scale.
    positive = err[err > 0]
    floor = float(positive.min()) if positive.size else 1e-9
    return (
        np.asarray(gm, dtype=float),
        np.asarray(em, dtype=float),
        np.where(err > 0, err, floor),
    )


@dataclass(frozen=True)
class RelaxationScan:
    """A globally-scanned profile of the relaxation strain.

    ``lo``/``hi`` bound ``lam`` at 95% by the F threshold on the cost *scaled by
    the achieved reduced chi-squared*. That scaling is the Birge correction and
    is exactly equivalent to inflating every uncertainty by ``sqrt(chi2_r)``, so
    a model whose shape is wrong widens rather than reporting a tightness it has
    not earned.
    """

    lam: float
    lo: float
    hi: float
    chi2_reduced: float
    y_inf: float
    y_0: float
    n_points: int
    grid: np.ndarray = field(repr=False)
    ssr: np.ndarray = field(repr=False)

    @property
    def at_edge(self) -> bool:
        """True when the optimum sits on the grid boundary: a bound, not a value."""
        return bool(self.lam <= self.grid[1] or self.lam >= self.grid[-2])


def fit_relaxation(
    strain: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    *,
    grid: tuple[float, float, int] = LAMBDA_GRID,
) -> RelaxationScan:
    """Fit ``y = y_inf + (y_0 - y_inf) exp(-(gamma - gamma_start)/lambda)``.

    The two amplitudes are eliminated exactly at each trial ``lambda`` (Golub and
    Pereyra 1973), so the reported optimum is the global one on the scanned grid
    and depends on no initial guess and no optimiser tolerance.

    Because the model is affine in ``y``, rescaling and offsetting ``y`` and
    ``sigma`` together changes nothing here. That invariance is precisely why
    ``lambda`` does not depend on ``phi``, ``phi_m`` or ``phi_0``: those enter
    the inversion only through the affine map ``f = A + B u``.
    """
    g = np.asarray(strain, dtype=float)
    yv = np.asarray(y, dtype=float)
    sv = np.asarray(sigma, dtype=float)
    if not (g.shape == yv.shape == sv.shape):
        raise ValueError(f"shape mismatch: {g.shape}, {yv.shape}, {sv.shape}")
    if g.size < 5:
        raise ValueError(f"fit_relaxation: {g.size} points is too few for 3 parameters")
    if np.any(sv <= 0) or not np.all(np.isfinite(sv)):
        raise ValueError("fit_relaxation: uncertainties must be finite and positive")

    g0 = float(g[0])
    w = 1.0 / sv
    y_w = yv * w
    lams = np.logspace(np.log10(grid[0]), np.log10(grid[1]), int(grid[2]))

    ssr = np.empty(lams.size)
    best_i, best_coef = 0, None
    for i, lam in enumerate(lams):
        basis = np.column_stack([np.ones_like(g), np.exp(-(g - g0) / lam)]) * w[:, None]
        coef, *_ = np.linalg.lstsq(basis, y_w, rcond=None)
        resid = y_w - basis @ coef
        ssr[i] = float(resid @ resid)
        if i == 0 or ssr[i] < ssr[best_i]:
            best_i, best_coef = i, coef

    dof = max(g.size - 3, 1)
    inside = np.flatnonzero(ssr <= ssr[best_i] * (1.0 + CHI2_1DOF_95 / dof))
    return RelaxationScan(
        lam=float(lams[best_i]),
        lo=float(lams[inside[0]]),
        hi=float(lams[inside[-1]]),
        chi2_reduced=float(ssr[best_i] / dof),
        y_inf=float(best_coef[0]),
        y_0=float(best_coef[0] + best_coef[1]),
        n_points=int(g.size),
        grid=lams,
        ssr=ssr,
    )


def smoothing_samples(strain: np.ndarray, window: float = SMOOTH_STRAIN) -> int:
    """Boxcar width in samples that spans ``window`` strain units on this record.

    Uses the median step rather than the mean so a single gap cannot set the
    width. Records here are uniformly sampled, so the two agree; the median is
    the safe choice if that ever stops being true.
    """
    strain = np.asarray(strain, dtype=float)
    if strain.size < 2:
        return SMOOTH_MIN_SAMPLES
    step = float(np.median(np.diff(strain)))
    if not np.isfinite(step) or step <= 0:
        raise ValueError("strain steps must be positive and finite")
    return max(SMOOTH_MIN_SAMPLES, int(round(window / step)))


def smooth(y: np.ndarray, size: int) -> np.ndarray:
    """Boxcar smoothing that does not restore the raw values at the edges.

    ``mode="nearest"`` extends the end samples rather than wrapping or zeroing,
    which matters because the artefact this is here to remove lives in the last
    handful of samples of several records. An earlier version of this smoothing
    copied the raw edge values back over the filtered ones, which reinstated
    exactly the spikes it had just removed.
    """
    y = np.asarray(y, dtype=float)
    if size < 1:
        raise ValueError(f"smoothing width must be >= 1, got {size}")
    if y.size < size:
        return y.copy()
    return uniform_filter1d(y, size=int(size), mode="nearest")


@dataclass(frozen=True)
class Window:
    """Which part of a reversal record the model can actually describe."""

    i_min: int
    gamma_min: float
    eta_min: float
    eta_frictionless: float
    i_start: int
    gamma_start: float
    n_post_min: int
    n_admissible: int

    @property
    def admissible_fraction(self) -> float:
        return self.n_admissible / self.n_post_min if self.n_post_min else 0.0

    @property
    def depth_below_floor(self) -> float:
        """``eta_min / eta_frictionless``. Below 1 the model cannot reach the dip."""
        return self.eta_min / self.eta_frictionless

    @property
    def usable(self) -> bool:
        return self.admissible_fraction >= MIN_ADMISSIBLE_FRACTION

    def __str__(self) -> str:
        return (
            f"minimum eta_r = {self.eta_min:.3f} at gamma = {self.gamma_min:.4f}, "
            f"{self.depth_below_floor:.3f}x the frictionless floor "
            f"{self.eta_frictionless:.3f}; the model becomes admissible at "
            f"gamma = {self.gamma_start:.4f} and covers "
            f"{100 * self.admissible_fraction:.1f}% of the post-minimum record"
        )


def admissible_window(
    strain: np.ndarray,
    eta_r: np.ndarray,
    phi: float,
    phi_m: float,
    phi_0: float,
    n: float = WC_EXPONENT,
    interior_max: float = INTERIOR_MAX_STRAIN,
    require_admissible: bool = True,
) -> Window:
    """Locate the reversal minimum and the first strain the model can describe.

    The minimum is taken on the smoothed record and restricted to
    ``gamma < interior_max``; see :data:`INTERIOR_MAX_STRAIN` for why. The start
    of the admissible window is the first sample at or after that minimum whose
    inverted ``f`` is non-negative, i.e. the first point that lies on or above
    the model's frictionless branch.

    ``require_admissible=False`` returns a window with ``n_admissible = 0``
    instead of raising when no part of the record is invertible. That is how
    :func:`fit_recovery` calls it: the clock is measured on ``u = phi_J/phi``,
    which is defined for every ``eta_r > 1``, so a record entirely below the
    floor still carries a clock even though it carries no ``f``. The raising
    default is kept because "the model cannot describe this record" is a real
    verdict and callers asking specifically about ``f`` should get it.
    """
    strain = np.asarray(strain, dtype=float)
    eta_r = np.asarray(eta_r, dtype=float)
    if strain.shape != eta_r.shape:
        raise ValueError(f"shape mismatch: strain {strain.shape} vs eta_r {eta_r.shape}")
    if strain.size < 3:
        raise ValueError("need at least three samples")
    if np.any(np.diff(strain) <= 0):
        raise ValueError("strain must be strictly increasing")

    es = smooth(eta_r, smoothing_samples(strain))
    interior = np.flatnonzero(strain < interior_max)
    if interior.size == 0:
        raise ValueError(
            f"no samples below gamma = {interior_max:g}; the reversal dip cannot be found"
        )
    i_min = int(interior[int(np.argmin(es[interior]))])

    f = frictional_fraction_from_viscosity(es, phi, phi_m, phi_0, n=n)
    post = np.arange(i_min, strain.size)
    adm = post[np.isfinite(f[post]) & (f[post] >= 0.0)]
    if adm.size == 0:
        if require_admissible:
            raise ValueError(
                "the whole post-minimum record sits below the frictionless branch; "
                "no part of it is invertible"
            )
        # No admissible sample: report the minimum as the start and let
        # n_admissible = 0 carry the verdict, rather than inventing an index.
        i_start, gamma_start = i_min, float(strain[i_min])
    else:
        i_start, gamma_start = int(adm[0]), float(strain[adm[0]])
    return Window(
        i_min=i_min,
        gamma_min=float(strain[i_min]),
        eta_min=float(es[i_min]),
        eta_frictionless=frictionless_relative_viscosity(phi, phi_0, n=n),
        i_start=i_start,
        gamma_start=gamma_start,
        n_post_min=int(post.size),
        n_admissible=int(adm.size),
    )


@dataclass(frozen=True)
class Recovery:
    """One rate's exponential rebuild of the frictional network, fitted in ``f``."""

    rate: float
    gamma_eff: float
    gamma_eff_err: float
    gamma_eff_lo: float
    gamma_eff_hi: float
    f_start: float
    f_plateau: float
    f_plateau_err: float
    window: Window
    gamma_fit_min: float
    gamma_fit_max: float
    n_fit: int
    rmse: float
    chi2_reduced: float
    success: bool

    @property
    def t_eff(self) -> float:
        """The same clock read as a time, ``gamma_eff / gamma_dot``, in seconds."""
        return self.gamma_eff / self.rate

    def __str__(self) -> str:
        return (
            f"gamma_dot = {self.rate:g} 1/s:  gamma_eff = {self.gamma_eff:.4f} "
            f"[{self.gamma_eff_lo:.4f}, {self.gamma_eff_hi:.4f}]  "
            f"(t_eff = {self.t_eff:.4f} s), f rebuilds {self.f_start:.3f} -> "
            f"{self.f_plateau:.3f}, chi2_r {self.chi2_reduced:.1f} on "
            f"{self.n_fit} bins over gamma {self.gamma_fit_min:.3f}-"
            f"{self.gamma_fit_max:.3f}"
        )


def fit_recovery(
    strain: np.ndarray,
    eta_r: np.ndarray,
    rate: float,
    phi: float,
    phi_m: float,
    phi_0: float,
    n: float = WC_EXPONENT,
    gamma_fit_max: float = np.inf,
) -> Recovery:
    """Fit ``f(gamma) = f_ss + (f_0 - f_ss) exp(-(gamma - gamma_start)/gamma_eff)``.

    Three free parameters, fitted by ordinary least squares on ``f`` rather than
    on ``eta_r``. Fitting on ``eta_r`` would measure the e-folding strain of the
    *viscosity*, and ``eta`` is strongly nonlinear in ``f`` near jamming, so the
    two differ by a factor that depends on how close the sample is to its
    jamming point. Only the ``f``-space constant is the ``gamma_0`` of the
    relaxation law.

    Two things about *how* it is fitted are load-bearing.

    **The fit starts at the reversal minimum, not where the model becomes
    admissible.** The admissible window (where the inverted ``f`` is
    non-negative) opens at ``gamma = 0.7`` to ``1.3`` on these records, which is
    at or past the recovery constant itself: a fit started there sees only the
    tail and returns a number set by the plateau. Starting at the minimum keeps
    the whole rebuild in the window, at the price of running the inversion below
    its own floor -- which is a stated extrapolation, not a hidden one, and
    which the affine invariance below makes harmless for ``gamma_eff``.

    **The amplitudes are eliminated exactly rather than optimised.** At fixed
    ``gamma_eff`` the model is linear in ``(f_0, f_ss)``, so they are solved by
    weighted linear least squares and the cost becomes a function of
    ``gamma_eff`` alone (Golub and Pereyra 1973). Scanning it on a log grid is
    globally reliable and the scan is the profile likelihood. The
    ``least_squares`` call this replaced returned ``gamma_eff`` = 0.587 or 25.3
    for the *same* record depending only on the convergence tolerance.

    The record is reduced to log-spaced bins first, one median per bin with the
    measured within-bin scatter as its uncertainty. The raw records hold 10 000
    to 25 000 strongly autocorrelated samples, so an unweighted fit at native
    sampling gives the plateau roughly 90% of the weight and reports an interval
    two orders of magnitude tighter than the analysis choices justify.

    ``gamma_fit_max`` truncates the upper end of the window, which is the
    stability check the gate applies.

    Raises
    ------
    ValueError
        If the window is too short to carry a recovery. Returning a number from
        a handful of plateau points would be reporting the window, not the
        material.
    """
    strain = np.asarray(strain, dtype=float)
    eta_r = np.asarray(eta_r, dtype=float)
    if not np.isfinite(rate) or rate <= 0:
        raise ValueError(f"shear rate must be finite and positive, got {rate!r}")

    # Reported for the diagnostics even though the fit does not start there,
    # and non-raising because the clock outlives the f-inversion (see above).
    win = admissible_window(
        strain, eta_r, phi, phi_m, phi_0, n=n, require_admissible=False
    )

    g_b, eta_b, err_b = log_bin(strain, eta_r)
    i0 = int(np.argmin(eta_b))
    keep = np.zeros(g_b.size, dtype=bool)
    keep[i0:] = True
    keep &= g_b <= gamma_fit_max
    g_fit = g_b[keep]
    if g_fit.size < 10:
        raise ValueError(
            f"gamma_dot = {rate:g}: {g_fit.size} bins inside the fit window, need 10"
        )
    eta_fit, eta_err = eta_b[keep], err_b[keep]

    # f is affine in eta only through u = phi_J/phi, so the uncertainty is
    # propagated analytically rather than by finite differences.
    f_fit = frictional_fraction_from_viscosity(eta_fit, phi, phi_m, phi_0, n=n)
    slope = phi / (phi_0 - phi_m)  # |df/du|
    f_err = slope * np.abs(d_reduced_jamming_distance(eta_fit, n=n)) * eta_err

    scan = fit_relaxation(g_fit, f_fit, f_err)
    return Recovery(
        rate=float(rate),
        gamma_eff=scan.lam,
        gamma_eff_err=0.5 * (scan.hi - scan.lo),
        gamma_eff_lo=scan.lo,
        gamma_eff_hi=scan.hi,
        f_start=scan.y_0,
        f_plateau=scan.y_inf,
        f_plateau_err=float("nan"),
        window=win,
        gamma_fit_min=float(g_fit[0]),
        gamma_fit_max=float(g_fit[-1]),
        n_fit=int(g_fit.size),
        rmse=float(np.sqrt(scan.ssr.min() / g_fit.size)),
        chi2_reduced=scan.chi2_reduced,
        success=not scan.at_edge,
    )


@dataclass(frozen=True)
class StrainClock:
    """Richards Eq. 4's two clocks, separated by how ``gamma_eff`` varies with rate.

    ``1/gamma_eff = 1/gamma_0 + 1/(gamma_dot t_r)`` is linear in ``1/gamma_dot``,
    so a straight line through ``(1/gamma_dot, 1/gamma_eff)`` gives ``1/gamma_0``
    as the intercept and ``1/t_r`` as the slope. The law requires **both to be
    positive**: a negative slope means ``gamma_eff`` falls as the rate rises,
    which the additive-time form cannot produce at any ``gamma_0`` or ``t_r``.
    """

    rates: np.ndarray
    gamma_eff: np.ndarray
    gamma_eff_err: np.ndarray
    gamma_0: float
    gamma_0_err: float
    t_r: float
    t_r_err: float
    slope: float
    slope_err: float
    intercept: float
    r_squared: float

    @property
    def form_holds(self) -> bool:
        """True when the fitted time clock is positive, i.e. the law can hold."""
        return bool(np.isfinite(self.t_r) and self.t_r > 0)

    @property
    def strain_spread(self) -> float:
        return float(np.max(self.gamma_eff) / np.min(self.gamma_eff))

    @property
    def time_spread(self) -> float:
        t = self.gamma_eff / self.rates
        return float(np.max(t) / np.min(t))

    @property
    def rate_span(self) -> float:
        return float(np.max(self.rates) / np.min(self.rates))

    def report(self) -> str:
        lines = [
            f"    rate span {self.rate_span:.0f}x over {self.rates.size} rates",
            f"    gamma_eff spread {self.strain_spread:.2f}x   "
            f"t_eff spread {self.time_spread:.2f}x   "
            f"({'STRAIN' if self.strain_spread < self.time_spread else 'TIME'} "
            f"is the tighter collapse)",
            f"    1/gamma_eff vs 1/gamma_dot: slope {self.slope:+.4f} "
            f"+/- {self.slope_err:.4f}, intercept {self.intercept:+.4f}, "
            f"r^2 = {self.r_squared:.3f}",
        ]
        if self.form_holds:
            lines.append(
                f"    gamma_0 = {self.gamma_0:.4f} +/- {self.gamma_0_err:.4f}, "
                f"t_r = {self.t_r:.4f} +/- {self.t_r_err:.4f} s"
            )
        else:
            lines.append(
                f"    t_r = {self.t_r:.4f} s is NOT positive: gamma_eff falls with rate, "
                f"the opposite of what an additive time clock produces. The additive "
                f"form of Richards Eq. 4 cannot describe this rate dependence."
            )
        return "\n".join(lines)


def fit_strain_clock(
    rates: np.ndarray, gamma_eff: np.ndarray, gamma_eff_err: np.ndarray | None = None
) -> StrainClock:
    """Weighted straight line through ``(1/gamma_dot, 1/gamma_eff)``.

    Weights are propagated from the ``gamma_eff`` errors by
    ``d(1/gamma_eff) = err / gamma_eff^2``. Uses ``numpy.polyfit`` rather than a
    hand-rolled normal equation.
    """
    rates = np.asarray(rates, dtype=float)
    gamma_eff = np.asarray(gamma_eff, dtype=float)
    if rates.shape != gamma_eff.shape:
        raise ValueError(f"shape mismatch: rates {rates.shape} vs gamma_eff {gamma_eff.shape}")
    if rates.size < 3:
        raise ValueError("separating two clocks needs at least three rates")
    if np.any(rates <= 0) or np.any(gamma_eff <= 0):
        raise ValueError("rates and gamma_eff must all be positive")

    x = 1.0 / rates
    y = 1.0 / gamma_eff
    if gamma_eff_err is None:
        w = None
    else:
        e = np.asarray(gamma_eff_err, dtype=float) / gamma_eff**2
        w = np.where(np.isfinite(e) & (e > 0), 1.0 / np.where(e > 0, e, 1.0), 0.0)
        if not np.any(w > 0):
            w = None

    coeffs, cov = np.polyfit(x, y, 1, w=w, cov="unscaled" if w is not None else True)
    slope, intercept = float(coeffs[0]), float(coeffs[1])
    slope_err, intercept_err = float(np.sqrt(cov[0, 0])), float(np.sqrt(cov[1, 1]))

    resid = y - np.polyval(coeffs, x)
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - float(resid @ resid) / ss_tot if ss_tot > 0 else float("nan")

    # gamma_0 = 1/intercept and t_r = 1/slope, with first-order error propagation.
    # Both are reported even when negative; the caller decides what that means,
    # and StrainClock.form_holds is the flag that says the law is violated.
    with np.errstate(divide="ignore", invalid="ignore"):
        gamma_0 = 1.0 / intercept if intercept != 0 else np.inf
        t_r = 1.0 / slope if slope != 0 else np.inf
        gamma_0_err = intercept_err / intercept**2 if intercept != 0 else np.inf
        t_r_err = slope_err / slope**2 if slope != 0 else np.inf

    return StrainClock(
        rates=rates,
        gamma_eff=gamma_eff,
        gamma_eff_err=(
            np.asarray(gamma_eff_err, dtype=float)
            if gamma_eff_err is not None
            else np.full(rates.shape, np.nan)
        ),
        gamma_0=float(gamma_0),
        gamma_0_err=float(gamma_0_err),
        t_r=float(t_r),
        t_r_err=float(t_r_err),
        slope=slope,
        slope_err=slope_err,
        intercept=intercept,
        r_squared=float(r2),
    )


def half_recovery_strain(
    strain: np.ndarray,
    eta_r: np.ndarray,
    phi: float,
    phi_m: float,
    phi_0: float,
    n: float = WC_EXPONENT,
    fraction: float = 0.5,
    plateau_from: float = 0.5,
) -> tuple[float, float]:
    """Strain from the reversal minimum to where ``f`` first reaches ``fraction``
    of its plateau. Returns ``(gamma_half, f_plateau)``.

    A landmark, not a fit. :func:`fit_recovery` returns an exponential constant
    measured over a window whose *start* depends on where the record climbs back
    above the model's frictionless floor, and that start moves with shear rate
    because a slower recovery crosses the floor later. So a rate trend in
    ``gamma_eff`` could in principle be a trend in the window rather than in the
    material. This landmark shares none of that: it is read off the record
    directly, needs no window, and asks only when the frictional fraction got
    half way back.

    The plateau is the median of ``f`` over the last ``plateau_from`` fraction of
    the record's strain range -- a median, not a mean, because several records
    carry a spike in their final samples.

    Returns NaN for the strain when the target is never reached inside the record.
    """
    strain = np.asarray(strain, dtype=float)
    eta_r = np.asarray(eta_r, dtype=float)
    if not (0.0 < fraction < 1.0):
        raise ValueError(f"fraction must lie in (0, 1), got {fraction!r}")
    if not (0.0 < plateau_from < 1.0):
        raise ValueError(f"plateau_from must lie in (0, 1), got {plateau_from!r}")

    win = admissible_window(strain, eta_r, phi, phi_m, phi_0, n=n)
    f = frictional_fraction_from_viscosity(
        smooth(eta_r, smoothing_samples(strain)), phi, phi_m, phi_0, n=n
    )

    g_lo = strain[-1] - plateau_from * (strain[-1] - strain[0])
    tail = f[strain >= g_lo]
    tail = tail[np.isfinite(tail)]
    if tail.size == 0:
        raise ValueError("no finite frictional fraction in the plateau window")
    f_plateau = float(np.median(tail))

    target = fraction * f_plateau
    after = np.arange(win.i_min, strain.size)
    hit = after[np.isfinite(f[after]) & (f[after] >= target)]
    if hit.size == 0:
        return float("nan"), f_plateau
    return float(strain[hit[0]] - strain[win.i_min]), f_plateau
