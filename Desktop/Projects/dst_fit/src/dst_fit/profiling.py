"""Profile-likelihood confidence intervals, for when the error bar is the lie.

Why the linearised error is not enough
--------------------------------------
:mod:`dst_fit.fitting` reports ``sqrt(diag(cov))`` from the Jacobian at the
optimum. That is a quadratic approximation to the cost surface, and it measures
how sharply the cost rises along each parameter axis *near the minimum*. Two
things break it here, both routinely:

1. **The surface is not quadratic.** Wyart-Cates diverges as ``phi -> phi_J``,
   so the cost is strongly asymmetric in ``phi_m``: moving phi_m down towards phi
   is catastrophic, moving it up is nearly free. A symmetric ``+/- sigma`` cannot
   express that, and reports the harmless direction and the catastrophic one as
   equally constrained.
2. **A near-degenerate direction is reported as sharp.** When two parameters
   trade off, the marginal width along each axis stays small while the joint
   region is long and thin. The correlation matrix shows this; the standard error
   does not.

The result on real data was a curve returning ``phi_m = 0.5876 +/- 0.0013`` from
a fit whose ``phi_m - phi`` was itself 0.0013. Nobody should believe that
interval, and the package had no way to say so.

What this module computes
-------------------------
For one parameter ``p``, hold it at a fixed value, re-minimise over all the
others, and record the resulting sum of squared residuals ``SSR(p)``. The
profile likelihood ratio

    F(p) = [SSR(p) - SSR_min] / [SSR_min / (n - k)]

is compared against the ``F(1, n-k)`` critical value. The set of ``p`` with
``F(p) <= F_crit`` is the confidence interval. This is the standard nonlinear
regression interval of Bates & Watts (*Nonlinear Regression Analysis and Its
Applications*, 1988, ch. 6); it is exact for a linear model and, unlike the Wald
interval, it follows the actual shape of the cost surface, so it comes out
asymmetric when the problem is asymmetric and unbounded when the data does not
bound the parameter.

What it still assumes
---------------------
That the residuals are independent draws from one distribution. On every real
flow curve tried here the runs test rejects that: the misfit is systematic. So
a profile interval computed on real data is still a **lower bound** on the true
uncertainty, and this module says so rather than implying otherwise. What it
fixes is the part that is fixable without knowing the true noise model.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import f as f_dist

from dst_fit.fitting import _PENALTY, FitResult, fit_model
from dst_fit.models import Model, get_model

# Below this rmse in log viscosity the residuals are float64 rounding, the
# denominator of the F statistic is meaningless, and the interval collapses to a
# width set by arithmetic rather than by data. Refuse rather than return it.
_MIN_RMSE = 1e-8

# Bisection refinements once the crossing is bracketed. Each one costs a single
# nuisance minimisation, and 10 halvings of the bracket is far finer than any
# quoted precision.
_REFINEMENTS = 10


@dataclass(frozen=True)
class Profile:
    """One parameter's profile likelihood and the interval read off it."""

    param: str
    best: float
    values: np.ndarray  # ascending, includes the optimum
    f_stat: np.ndarray  # profile F statistic at each value
    ssr_min: float
    n_points: int
    n_free: int
    level: float
    f_crit: float
    lower: float
    upper: float
    lower_at_bound: bool
    upper_at_bound: bool
    wald_stderr: float
    bounds: tuple[float, float]

    @property
    def bounded(self) -> bool:
        """True when the data closes the interval on both sides."""
        return not (self.lower_at_bound or self.upper_at_bound)

    @property
    def interval(self) -> tuple[float, float]:
        return (self.lower, self.upper)

    @property
    def width(self) -> float:
        return self.upper - self.lower

    @property
    def wald_width(self) -> float:
        """Width of the two-sided Wald interval at the same confidence level."""
        z = np.sqrt(self.f_crit)
        return float(2.0 * z * self.wald_stderr)

    @property
    def widening(self) -> float:
        """How many times wider the profile interval is than the Wald one."""
        return float(self.width / self.wald_width) if self.wald_width > 0 else np.inf

    @property
    def asymmetry(self) -> float:
        """Ratio of the longer arm to the shorter one. 1.0 is symmetric."""
        up, down = self.upper - self.best, self.best - self.lower
        lo, hi = sorted((abs(down), abs(up)))
        return float(hi / lo) if lo > 0 else np.inf

    def contains(self, value: float) -> bool:
        return bool(self.lower <= value <= self.upper)

    def summary(self) -> str:
        arrow = ""
        if self.lower_at_bound and self.upper_at_bound:
            arrow = "   UNBOUNDED both ways: the data does not constrain this parameter"
        elif self.lower_at_bound:
            arrow = "   unbounded below (hit the parameter bound)"
        elif self.upper_at_bound:
            arrow = "   unbounded above (hit the parameter bound)"
        return (
            f"{self.param} = {self.best:.4f}  "
            f"[{self.lower:.4f}, {self.upper:.4f}] at {100 * self.level:.0f}%"
            f"{arrow}\n"
            f"    Wald +/- {self.wald_stderr:.4f} would give "
            f"[{self.best - np.sqrt(self.f_crit) * self.wald_stderr:.4f}, "
            f"{self.best + np.sqrt(self.f_crit) * self.wald_stderr:.4f}]; "
            f"the profile is {self.widening:.1f}x wider and "
            f"{self.asymmetry:.1f}x asymmetric"
        )


def _nuisance_minimiser(model: Model, sigma, log_eta, index: int):
    """Return f(fixed_value, warm_start) -> (ssr, params) with ``index`` held fixed."""
    lo = np.array(model.bounds[0], dtype=float)
    hi = np.array(model.bounds[1], dtype=float)
    keep = [j for j in range(model.n_params) if j != index]

    def run(value: float, warm: np.ndarray) -> tuple[float, np.ndarray]:
        full = np.array(warm, dtype=float)
        full[index] = value

        def residual(free):
            p = full.copy()
            p[keep] = free
            with np.errstate(all="ignore"):
                pred = np.asarray(model(sigma, *p), dtype=float)
            out = np.full(pred.shape, _PENALTY)
            ok = np.isfinite(pred) & (pred > 0)
            out[ok] = np.log(pred[ok]) - log_eta[ok]
            return out

        x0 = np.clip(full[keep], lo[keep] + 1e-12, hi[keep] - 1e-12)
        sol = least_squares(
            residual, x0, bounds=(lo[keep], hi[keep]), method="trf", max_nfev=3000
        )
        out = full.copy()
        out[keep] = sol.x
        return float(sol.fun @ sol.fun), out

    return run


def profile_likelihood(
    model: Model | str,
    stress: np.ndarray,
    viscosity: np.ndarray,
    param: str,
    *,
    branch: slice | None = None,
    level: float = 0.95,
    fit: FitResult | None = None,
    max_steps: int = 120,
    step_growth: float = 1.18,
    n_starts: int = 12,
    seed: int = 0,
) -> Profile:
    """Profile ``param`` and return its likelihood-ratio confidence interval.

    Walks outward from the best fit in both directions, re-minimising over the
    remaining parameters at each step and warm-starting from the previous
    solution so the profile follows one continuous branch of the cost surface
    rather than hopping between local minima. Stops at the parameter's bound or
    once the statistic is well past its critical value, then bisects to locate
    the crossing.

    The step starts at a quarter of the Wald standard error and grows
    geometrically. That matters precisely in the pathological case: where the
    Wald error is absurdly small, a fixed step sized from it would need
    thousands of iterations to escape, and a fixed step sized from the bounds
    would step straight over the minimum.

    Parameters
    ----------
    fit
        A completed fit to profile around. Recomputed if omitted; pass it in
        when profiling several parameters of the same fit.
    """
    mdl: Model = get_model(model) if isinstance(model, str) else model
    if param not in mdl.param_names:
        raise ValueError(f"{mdl.name} has no parameter {param!r}; has {mdl.param_names}")
    index = mdl.param_names.index(param)

    if fit is None:
        fit = fit_model(mdl, stress, viscosity, branch=branch, n_starts=n_starts, seed=seed)

    sigma = np.asarray(stress, dtype=float)
    eta = np.asarray(viscosity, dtype=float)
    if branch is not None:
        sigma, eta = sigma[branch], eta[branch]
    good = np.isfinite(sigma) & np.isfinite(eta) & (sigma > 0) & (eta > 0)
    sigma, eta = sigma[good], eta[good]
    log_eta = np.log(eta)

    n, k = int(sigma.size), mdl.n_params
    dof = n - k
    if dof < 1:
        raise ValueError(
            f"profile_likelihood: {n} points and {k} parameters leaves {dof} degrees "
            "of freedom; there is no scale to compare the profile against"
        )
    f_crit = float(f_dist.ppf(level, 1, dof))
    ssr_min = float(fit.residuals @ fit.residuals)
    rmse = float(np.sqrt(ssr_min / n))
    if rmse < _MIN_RMSE:
        raise ValueError(
            f"profile_likelihood: rmse in log viscosity is {rmse:.2e}, at the level of "
            "float64 rounding. There is no residual scale to compare the profile "
            "against, so the interval would be a statement about arithmetic precision. "
            "This is a synthetic noise-free curve, not a measurement."
        )

    minimise = _nuisance_minimiser(mdl, sigma, log_eta, index)
    lo_b, hi_b = float(mdl.bounds[0][index]), float(mdl.bounds[1][index])
    best = float(fit.params[index])
    wald = float(fit.stderr[index])

    def f_of(ssr: float) -> float:
        return (ssr - ssr_min) / (ssr_min / dof)

    def walk(direction: int) -> tuple[list[float], list[float], float, bool]:
        """Step outward until the statistic is past overshoot or the bound is hit."""
        limit = hi_b if direction > 0 else lo_b
        span = abs(limit - best)
        if span <= 0:
            return [], [], limit, True
        step = max(wald / 4.0 if np.isfinite(wald) and wald > 0 else span / 200.0,
                   span / 5000.0)
        warm = np.array(fit.params, dtype=float)
        vals: list[float] = []
        stats: list[float] = []
        prev_val = best
        at_bound = False

        for _ in range(max_steps):
            nxt = prev_val + direction * step
            if (direction > 0 and nxt >= limit) or (direction < 0 and nxt <= limit):
                nxt, at_bound = limit, True
            ssr, warm = minimise(nxt, warm)
            stat = f_of(ssr)
            vals.append(nxt)
            stats.append(stat)
            if stat >= f_crit:
                # bracketed: refine by bisection on the true profile
                a, b = prev_val, nxt
                for _ in range(_REFINEMENTS):
                    mid = 0.5 * (a + b)
                    s_mid, _ = minimise(mid, warm)
                    if f_of(s_mid) < f_crit:
                        a = mid
                    else:
                        b = mid
                return vals, stats, float(b), False
            if at_bound:
                return vals, stats, limit, True
            prev_val = nxt
            step *= step_growth
        # Ran out of steps without either crossing or reaching the bound. With a
        # geometrically growing step this is unreachable in practice, but if it
        # happens the interval was not closed by the data and must not be
        # reported as though it were.
        return vals, stats, prev_val, True

    up_vals, up_stats, upper, upper_at_bound = walk(+1)
    dn_vals, dn_stats, lower, lower_at_bound = walk(-1)

    values = np.array(dn_vals[::-1] + [best] + up_vals, dtype=float)
    stats = np.array(dn_stats[::-1] + [0.0] + up_stats, dtype=float)

    return Profile(
        param=param,
        best=best,
        values=values,
        f_stat=stats,
        ssr_min=ssr_min,
        n_points=n,
        n_free=k,
        level=level,
        f_crit=f_crit,
        lower=float(lower),
        upper=float(upper),
        lower_at_bound=bool(lower_at_bound),
        upper_at_bound=bool(upper_at_bound),
        wald_stderr=wald,
        bounds=(lo_b, hi_b),
    )


@dataclass(frozen=True)
class IntervalIntersection:
    """Is there one value of a shared constant compatible with every curve?

    The right question to ask of N intervals, and it needs no Gaussian
    assumption and no weighted mean. Either a single value sits inside all of
    them or it does not, and if it does not, the largest subset that does agree
    says how bad the disagreement is.
    """

    lower: float
    upper: float
    n_total: int
    n_overlapping: int  # curves containing the best common value
    non_empty: bool
    limiting_low: str  # which curve sets the lower edge
    limiting_high: str

    @property
    def width(self) -> float:
        return max(self.upper - self.lower, 0.0)

    def report(self) -> str:
        if self.non_empty:
            return (
                f"one common value is possible: [{self.lower:.4f}, {self.upper:.4f}], "
                f"width {self.width:.4f}. Set below by {self.limiting_low} and above by "
                f"{self.limiting_high}"
            )
        return (
            f"no single value fits all {self.n_total} curves: the highest lower limit "
            f"is {self.lower:.4f} ({self.limiting_low}) and the lowest upper limit is "
            f"{self.upper:.4f} ({self.limiting_high}), a gap of "
            f"{self.lower - self.upper:.4f}. At most {self.n_overlapping} of "
            f"{self.n_total} curves can share a value"
        )


def interval_intersection(
    profiles: list[Profile], labels: list[str] | None = None
) -> IntervalIntersection:
    """Intersect N confidence intervals for the same constant.

    This is the consistency test that survives asymmetric and unbounded
    intervals, which the 2-sigma-from-the-weighted-mean rule does not: a
    weighted mean of asymmetric intervals is not defined, and an unbounded one
    has no sigma at all.
    """
    if len(profiles) < 2:
        raise ValueError("intersecting needs at least two intervals")
    labels = labels or [p.param + f"[{i}]" for i, p in enumerate(profiles)]

    lows = np.array([p.lower for p in profiles], dtype=float)
    highs = np.array([p.upper for p in profiles], dtype=float)
    i_lo, i_hi = int(np.argmax(lows)), int(np.argmin(highs))
    lower, upper = float(lows[i_lo]), float(highs[i_hi])

    if lower <= upper:
        return IntervalIntersection(
            lower=lower, upper=upper, n_total=len(profiles),
            n_overlapping=len(profiles), non_empty=True,
            limiting_low=labels[i_lo], limiting_high=labels[i_hi],
        )

    # Empty. The largest subset sharing a value is the maximum overlap depth,
    # found by sweeping the interval endpoints.
    edges = np.unique(np.concatenate([lows, highs]))
    depth = max(
        int(np.count_nonzero((lows <= e) & (highs >= e))) for e in edges
    )
    return IntervalIntersection(
        lower=lower, upper=upper, n_total=len(profiles),
        n_overlapping=depth, non_empty=False,
        limiting_low=labels[i_lo], limiting_high=labels[i_hi],
    )


# --------------------------------------------------------------------------
# the correction the profile does not make
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class BootstrapInterval:
    """A parameter interval from resampling *blocks* of residuals.

    The profile interval and the Wald interval disagree about the shape of the
    cost surface but agree about the noise model: both assume the residuals are
    independent draws. On every real flow curve in this project the runs test
    rejects that, and rejects it hard. Neither interval can see the difference,
    and on the Guy et al. series the two come out equal to within a per cent,
    which is how the assumption was caught. Two methods that share an assumption
    agreeing is not corroboration of the assumption.

    A moving-block bootstrap keeps the assumption that is defensible (the
    residual pattern is a repeatable feature of this measurement) and drops the
    one that is not (that consecutive residuals are unrelated). Blocks of length
    ``n / n_runs`` reproduce the observed run structure by construction.
    """

    param: str
    best: float
    lower: float
    upper: float
    level: float
    values: np.ndarray  # the bootstrap distribution, for a histogram
    n_resamples: int
    n_converged: int
    block_length: int
    n_runs: int
    n_points: int
    wald_stderr: float

    @property
    def width(self) -> float:
        return self.upper - self.lower

    @property
    def wald_width(self) -> float:
        return float(2.0 * 1.959963984540054 * self.wald_stderr)

    @property
    def widening(self) -> float:
        """How many times wider than the Wald interval at the same level."""
        return float(self.width / self.wald_width) if self.wald_width > 0 else np.inf

    @property
    def asymmetry(self) -> float:
        up, down = self.upper - self.best, self.best - self.lower
        lo, hi = sorted((abs(down), abs(up)))
        return float(hi / lo) if lo > 0 else np.inf

    def contains(self, value: float) -> bool:
        return bool(self.lower <= value <= self.upper)

    def summary(self) -> str:
        return (
            f"{self.param} = {self.best:.4f}  [{self.lower:.4f}, {self.upper:.4f}] "
            f"at {100 * self.level:.0f}%, moving-block bootstrap\n"
            f"    blocks of {self.block_length} from {self.n_points} points "
            f"({self.n_runs} sign runs), {self.n_converged}/{self.n_resamples} "
            f"resamples converged; {self.widening:.1f}x the Wald width, "
            f"{self.asymmetry:.1f}x asymmetric"
        )


def block_length_from_runs(residuals: np.ndarray) -> tuple[int, int]:
    """Block length and run count for residuals ordered by the x axis.

    A sequence of ``n`` residuals with ``R`` sign runs carries, for the purpose
    of counting independent information, roughly ``R`` independent blocks of
    length ``n/R``. Under independence ``R ~ n/2``, giving blocks of 2; a badly
    misspecified fit leaves a handful of long runs and blocks to match.

    Using the runs count rather than a fitted autocorrelation time is
    deliberate: it is the same statistic
    :func:`dst_fit.identifiability.runs_test` already reports, so the block
    length and the misspecification verdict cannot disagree with each other.
    """
    from dst_fit.identifiability import runs_test

    n_runs, _, _ = runs_test(residuals)
    n = int(np.asarray(residuals).size)
    return max(1, int(round(n / max(n_runs, 1)))), int(n_runs)


def block_bootstrap(
    model: Model | str,
    stress: np.ndarray,
    viscosity: np.ndarray,
    param: str,
    *,
    branch: slice | None = None,
    level: float = 0.95,
    fit: FitResult | None = None,
    n_resamples: int = 400,
    block_length: int | None = None,
    seed: int = 0,
    n_starts: int = 12,
) -> BootstrapInterval:
    """Percentile interval for ``param`` from a moving-block residual bootstrap.

    Residuals, not data points, are resampled: the stresses are fixed by the
    experiment's design, so resampling ``(sigma, eta)`` pairs would throw that
    design away and widen the interval for the wrong reason.

    Each resample draws overlapping blocks of length ``block_length`` from the
    fitted residuals, lays them end to end, and adds them back to the fitted
    curve in log space. The residuals are log-viscosity, so this is
    multiplicative noise on eta, which is how the data actually scatters.

    ``block_length`` defaults to ``n / n_runs`` from
    :func:`block_length_from_runs`. Setting it to 1 recovers the ordinary
    residual bootstrap, i.e. the independence assumption, and is the control
    worth running: the gap between the two is the price of that assumption.

    A percentile interval is used rather than basic or BCa. It is the
    conventional default, it respects the parameter's bounds automatically, and
    the bias correction BCa applies needs a jackknife costing another ``n`` fits
    for a refinement smaller than the effect being measured here.
    """
    mdl: Model = get_model(model) if isinstance(model, str) else model
    if param not in mdl.param_names:
        raise ValueError(f"{mdl.name} has no parameter {param!r}; has {mdl.param_names}")
    index = mdl.param_names.index(param)

    if fit is None:
        fit = fit_model(mdl, stress, viscosity, branch=branch, n_starts=n_starts, seed=seed)

    sigma = np.asarray(stress, dtype=float)
    eta = np.asarray(viscosity, dtype=float)
    if branch is not None:
        sigma, eta = sigma[branch], eta[branch]
    good = np.isfinite(sigma) & np.isfinite(eta) & (sigma > 0) & (eta > 0)
    sigma, eta = sigma[good], eta[good]

    resid = np.asarray(fit.residuals, dtype=float)
    n = int(sigma.size)
    if resid.size != n:
        raise ValueError(
            f"block_bootstrap: {resid.size} residuals for {n} points; the fit passed "
            "here was not made on this data"
        )
    if n < 8:
        raise ValueError(f"block_bootstrap: {n} points is too few to resample blocks from")

    auto_len, n_runs = block_length_from_runs(resid)
    L = int(block_length) if block_length is not None else auto_len
    L = max(1, min(L, n))

    # The fitted curve the residuals are added back to. A non-finite prediction
    # cannot host a residual, so those points leave every resample rather than
    # being given a fabricated value.
    with np.errstate(all="ignore"):
        pred = np.asarray(mdl(sigma, *fit.params), dtype=float)
    usable = np.isfinite(pred) & (pred > 0)
    if int(usable.sum()) <= mdl.n_params:
        raise ValueError("block_bootstrap: the fitted curve is jammed at almost every point")

    n_blocks = int(np.ceil(n / L))
    n_offsets = n - L + 1
    rng = np.random.default_rng(seed)
    lo = np.array(mdl.bounds[0], dtype=float)
    hi = np.array(mdl.bounds[1], dtype=float)

    draws: list[float] = []
    for _ in range(n_resamples):
        starts = rng.integers(0, n_offsets, size=n_blocks)
        r_star = np.concatenate([resid[s : s + L] for s in starts])[:n]
        eta_star = pred * np.exp(r_star)
        try:
            boot = fit_model(mdl, sigma[usable], eta_star[usable], n_starts=1, seed=0)
        except (RuntimeError, ValueError):
            continue
        # A resample that pins the parameter on a bound is kept, not dropped:
        # dropping it would censor exactly the tail that is trying to escape and
        # narrow the interval for the wrong reason.
        draws.append(float(np.clip(boot.params[index], lo[index], hi[index])))

    values = np.asarray(draws, dtype=float)
    if values.size < 20:
        raise RuntimeError(
            f"block_bootstrap: only {values.size} of {n_resamples} resamples converged; "
            "the interval would be noise"
        )
    a = 0.5 * (1.0 - level)
    return BootstrapInterval(
        param=param,
        best=float(fit.params[index]),
        lower=float(np.quantile(values, a)),
        upper=float(np.quantile(values, 1.0 - a)),
        level=level,
        values=values,
        n_resamples=n_resamples,
        n_converged=int(values.size),
        block_length=L,
        n_runs=n_runs,
        n_points=n,
        wald_stderr=float(fit.stderr[index]),
    )
