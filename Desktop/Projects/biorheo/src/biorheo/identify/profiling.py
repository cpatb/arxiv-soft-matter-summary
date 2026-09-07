"""Profile-likelihood intervals that end where the F statistic crosses, not on a grid.

Why this module exists
----------------------
``dst_fit.profiling.profile_likelihood`` walks outward from the best fit with

    step = max(wald/4, span/5000),   span = |box_bound - best|

and, once the statistic passes ``f_crit``, runs ``_REFINEMENTS = 10`` bisections on the
bracket and returns its **outer** end. Both halves of that are sized from the parameter's
declared BOX BOUND, which is a guard-rail so the fitter cannot wander -- not a statement
about the data. For a viscosity bounded at ``1e6`` the step floor is ``200``, on a
parameter whose value is ``0.0038``: the first step overshoots by five orders of
magnitude and ten bisections leave a resolution of ``200/1024 = 0.195``.

The returned endpoint is then a point on the grid ``best + m * span / 5.12e6`` for an
integer ``m <= 1024``. Measured over this registry (``scripts/profile_interval_diagnosis.py``):
**203 of 296 closed upper ends** sit exactly on that grid, worst 1639x the best-fit value
above it.

``dst_fit`` is consumed unchanged (``tests/test_gate_0_1_dstfit_unchanged.py`` exists to
detect the absence of change), so the correction lives here. Gate and predictions were
pre-registered in ``data/profile_interval_prereg.md`` before this file was written.

What is different
-----------------
1. **The step is sized from the parameter, not from its box.** ``wald/4`` where the Wald
   error is usable, else a fixed fraction of ``|best|``, and only if both are unusable
   does it fall back to the box. Which rule fired is recorded in ``step_rule``.
2. **The endpoint is located by Brent's method on ``F(p) - f_crit``**, not by a fixed
   count of bisections, so its accuracy is set by a tolerance rather than by how wide the
   first step happened to be.
3. **The result carries its own check.** ``f_at_upper`` / ``f_at_lower`` are the statistic
   re-evaluated at the returned endpoints. A correct closed endpoint has
   ``|F - f_crit| / f_crit`` at the level of the root-finder tolerance; gate G13.3 reads
   exactly this field, so the object cannot claim a crossing it did not find.

What is deliberately NOT different
----------------------------------
The inner solve. ``dst_fit.profiling._nuisance_minimiser`` is imported and used unchanged,
so the two routines re-minimise the nuisance parameters identically and a difference
between them is a difference in the WALK, not in the objective. That is a private name and
therefore a coupling; ``tests/test_profiling_wrapper.py`` pins its presence and its
contract, and fails loudly rather than silently falling back if it moves.

One-parameter models have no nuisance parameters at all (``keep`` is empty, and
``least_squares`` has nothing to vary). Those are handled here by evaluating the sum of
squares directly, which is also what makes the ``newtonian`` closed form an exact anchor
for this module -- see ``data/profile_interval_prereg.md`` section 3, G13.4.

What this does NOT fix
----------------------
Nothing about the FIT. A profile interval is a statement about the cost surface of one
model on one dataset; if the model misfits, a correct interval on a wrong model is still
a wrong uncertainty. `residual_structure` is the field that says so.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq
from scipy.stats import f as f_dist

from dst_fit.fitting import _PENALTY, FitResult
from dst_fit.models import Model

# A private name from dst_fit, imported deliberately: using the SAME nuisance minimiser is
# what makes this routine's interval comparable with dst_fit's. Pinned by
# tests/test_profiling_wrapper.py. Import errors are not caught -- a silent fallback to a
# private copy of the objective is exactly how two implementations drift apart.
from dst_fit.profiling import _nuisance_minimiser

# Growth of the outward step. dst_fit uses 1.18; a larger ratio is used here because the
# step now STARTS at the parameter's own scale instead of at a fraction of its box, so the
# walk would otherwise need ~40 nuisance solves to reach a crossing 1000 Wald errors out.
# Brent's method is superlinear, so a wider bracket costs far less than the extra steps.
STEP_GROWTH = 1.35

# Initial step, in units of the Wald standard error. Same as dst_fit's.
WALD_FRACTION = 0.25

# Fallback initial step as a fraction of |best|, used when the Wald error is zero or
# non-finite (a flat or degenerate direction). Deliberately small: overshooting is what
# this module exists to stop.
VALUE_FRACTION = 1e-3

MAX_STEPS = 300

# Relative x-tolerance for Brent. The nuisance re-minimisation is itself only good to
# roughly 1e-10 in SSR, so tightening past this measures the inner optimiser's noise.
ROOT_RTOL = 1e-8

# Below this rmse in log viscosity there is no residual scale to compare a profile
# against; the same guard dst_fit applies, restated here because this routine can be
# called without going through dst_fit's.
_MIN_RMSE = 1e-8

# A located endpoint is only reported as a confidence limit if it survives this check:
# |F - f_crit| / f_crit within VERIFY_TOL, and F genuinely below f_crit just inside it and
# at or above just outside. On a degenerate surface the nuisance re-minimisation lands on
# different local optima at neighbouring parameter values, so F is not a function there --
# measured on two Carreau-Yasuda blood records, F moves 0.239 -> 10.06 between points
# 0.1 % apart. A root finder will still return SOMETHING from such a sign change; it is a
# discontinuity, not a crossing, and this check is what stops it being written down as an
# uncertainty. VERIFY_TOL is deliberately the same number as gate G13.3's bar.
VERIFY_TOL = 0.01
VERIFY_DELTA = 1e-3   # relative offset for the flanking points


@dataclass(frozen=True)
class ProfileInterval:
    """A profile-likelihood interval whose ends are F-crossings.

    Field names match ``dst_fit.profiling.Profile`` where they mean the same thing, so
    ``biorheo.identify.verdict_for`` accepts either without a shim.
    """

    param: str
    best: float
    lower: float
    upper: float
    lower_at_bound: bool
    upper_at_bound: bool
    wald_stderr: float
    f_crit: float
    n_points: int
    n_params: int
    bounds: tuple[float, float]
    # the self-check: the statistic re-evaluated at the returned endpoints
    f_at_lower: float
    f_at_upper: float
    step_rule: str
    n_evals: int
    # How far below zero the statistic went. F < 0 means the nuisance re-minimisation
    # found a LOWER SSR than the fit this interval is centred on, i.e. the reported
    # optimum is not the optimum and the interval is drawn around the wrong point.
    # Recorded rather than raised: it is a statement about the FIT, and the caller who
    # owns the fit is the one who can act on it.
    worst_negative_f: float = 0.0
    # values where a warm-started nuisance solve claimed F >= f_crit and a second
    # start refuted it; a nonzero count means multiple local optima along the profile
    false_crossings: tuple[float, ...] = ()
    # An end is UNVERIFIED when a sign change was found but it is a discontinuity of the
    # profile rather than a crossing of it. Such an end is reported as open (the endpoint
    # is set to the box bound and ``*_at_bound`` is True), because the data did not close
    # it; the value the root finder returned and the F there are kept so the reader can
    # see what was rejected and why.
    lower_unverified: tuple[float, float] | None = None
    upper_unverified: tuple[float, float] | None = None
    values: np.ndarray = field(repr=False, default_factory=lambda: np.empty(0))
    stats: np.ndarray = field(repr=False, default_factory=lambda: np.empty(0))

    @property
    def dof(self) -> int:
        return self.n_points - self.n_params

    @property
    def bounded(self) -> bool:
        return not (self.lower_at_bound or self.upper_at_bound)

    @property
    def interval(self) -> tuple[float, float]:
        return (self.lower, self.upper)

    @property
    def width(self) -> float:
        return self.upper - self.lower

    @property
    def wald_width(self) -> float:
        return 2.0 * np.sqrt(self.f_crit) * self.wald_stderr

    @property
    def widening(self) -> float:
        """How many times wider the profile interval is than the Wald one."""
        w = self.wald_width
        return float(self.width / w) if w > 0 and np.isfinite(w) else float("inf")

    @property
    def asymmetry(self) -> float:
        up, dn = self.upper - self.best, self.best - self.lower
        if dn <= 0 or not np.isfinite(dn) or not np.isfinite(up):
            return float("inf")
        return float(max(up / dn, dn / up))

    def contains(self, value: float) -> bool:
        return bool(self.lower <= value <= self.upper)

    def crossing_error(self) -> float:
        """Worst ``|F - f_crit| / f_crit`` over the CLOSED ends. G13.3 reads this.

        An end at a box bound is not a crossing and contributes nothing; a profile with
        both ends open returns 0.0, which is correct and not a pass by default -- such a
        profile is labelled ``not_identifiable`` and carries no interval to quote.
        """
        errs = []
        if not self.upper_at_bound:
            errs.append(abs(self.f_at_upper - self.f_crit) / self.f_crit)
        if not self.lower_at_bound:
            errs.append(abs(self.f_at_lower - self.f_crit) / self.f_crit)
        return float(max(errs)) if errs else 0.0

    def rejected_error(self) -> float:
        """Worst ``|F - f_crit| / f_crit`` over the ends REJECTED as unverified.

        Reported beside :meth:`crossing_error` rather than folded into it, because the
        two answer different questions. ``crossing_error`` grades the ends this routine
        reports as confidence limits, which is the population gate G13.3 was written over.
        This one says how far the rejected candidates were from being limits -- a small
        value here would mean the verification is throwing away good endpoints, and a
        large one means it is catching genuine discontinuities. Measured: 0.69 to 0.95.
        """
        errs = [abs(u[1] - self.f_crit) / self.f_crit
                for u in (self.lower_unverified, self.upper_unverified) if u is not None]
        return float(max(errs)) if errs else 0.0

    @property
    def verified(self) -> bool:
        """True when no end was rejected as a discontinuity rather than a crossing."""
        return self.lower_unverified is None and self.upper_unverified is None


def _direct_ssr(model: Model, x: np.ndarray, log_eta: np.ndarray):
    """Nuisance minimiser for a model with ONE parameter: nothing to minimise over."""

    def run(value: float, warm: np.ndarray) -> tuple[float, np.ndarray]:
        p = np.array(warm, dtype=float)
        p[0] = value
        with np.errstate(all="ignore"):
            pred = np.asarray(model(x, *p), dtype=float)
        # dst_fit.fitting._PENALTY, imported rather than copied: a residual penalty
        # that drifts between the two would make the one-parameter path disagree with
        # the nuisance path on exactly the models where predictions go non-positive.
        out = np.full(pred.shape, _PENALTY)
        ok = np.isfinite(pred) & (pred > 0)
        out[ok] = np.log(pred[ok]) - log_eta[ok]
        return float(out @ out), p

    return run


def profile_interval(
    model: Model,
    x: np.ndarray,
    eta: np.ndarray,
    param: str,
    *,
    fit: FitResult,
    level: float = 0.95,
) -> ProfileInterval:
    """Profile ``param`` and return the interval whose ends are ``F = f_crit``.

    ``x`` is the model's independent variable -- shear rate for a ``RateModel``, stress
    for a ``dst_fit`` model. ``fit`` is required rather than recomputed: profiling every
    parameter of one curve must profile around ONE optimum, and letting this function
    refit would let two parameters of the same record be profiled around two different
    ones.
    """
    if param not in model.param_names:
        raise ValueError(f"{model.name} has no parameter {param!r}; has {model.param_names}")
    index = model.param_names.index(param)

    x = np.asarray(x, dtype=float)
    eta = np.asarray(eta, dtype=float)
    good = np.isfinite(x) & np.isfinite(eta) & (x > 0) & (eta > 0)
    x, eta = x[good], eta[good]
    log_eta = np.log(eta)

    n, k = int(x.size), model.n_params
    dof = n - k
    if dof < 1:
        raise ValueError(
            f"profile_interval: {n} points and {k} parameters leaves {dof} degrees of "
            "freedom; there is no scale to compare the profile against"
        )
    f_crit = float(f_dist.ppf(level, 1, dof))
    ssr_min = float(fit.residuals @ fit.residuals)
    rmse = float(np.sqrt(ssr_min / n))
    if rmse < _MIN_RMSE:
        raise ValueError(
            f"profile_interval: rmse in log viscosity is {rmse:.2e}, at the level of "
            "float64 rounding. There is no residual scale to compare the profile "
            "against, so the interval would be a statement about arithmetic precision. "
            "This is a synthetic noise-free curve, not a measurement."
        )

    minimise = (_direct_ssr(model, x, log_eta) if k == 1
                else _nuisance_minimiser(model, x, log_eta, index))
    lo_b, hi_b = float(model.bounds[0][index]), float(model.bounds[1][index])
    best = float(fit.params[index])
    wald = float(fit.stderr[index])

    n_evals = 0
    # Points where a warm-started solve reported F >= f_crit and the second nuisance
    # start did not. Reported, not swallowed: a nonzero count says the cost surface has
    # more than one local optimum along this parameter's profile.
    false_crossings: list[float] = []
    best_params = np.array(fit.params, dtype=float)

    def f_of(value: float, warm: np.ndarray) -> tuple[float, np.ndarray]:
        nonlocal n_evals
        ssr, params = minimise(value, warm)
        n_evals += 1
        return (ssr - ssr_min) / (ssr_min / dof), params

    def f_best_of(value: float, warm: np.ndarray) -> tuple[float, np.ndarray]:
        """``F`` at ``value`` taking the BETTER of two nuisance starts.

        The profile is defined by the *minimum* of SSR over the nuisance parameters, so
        any single warm-started solve returns an upper bound on it and a value that
        depends on the path taken to get there. Two starts are used -- the rolling one
        (which follows the branch) and the best fit (which does not) -- and the lower SSR
        wins, because lower is unambiguously closer to the minimum being profiled.

        This is not a refinement of taste. Without it the walk reports a crossing wherever
        the warm-started solve happens to land on a worse optimum, which on two of the 100
        blood records manufactured a closed interval for a parameter whose profile in fact
        saturates at ``F = 0.19`` and never reaches ``f_crit`` at all -- i.e. it turned an
        UNBOUNDED parameter into a measured one.
        """
        nonlocal n_evals
        a_ssr, a_p = minimise(value, warm)
        b_ssr, b_p = minimise(value, best_params)
        n_evals += 2
        ssr, params = (a_ssr, a_p) if a_ssr <= b_ssr else (b_ssr, b_p)
        return (ssr - ssr_min) / (ssr_min / dof), params

    # --- the initial step, sized from the parameter and not from its box ----------------
    from_wald = WALD_FRACTION * wald if (np.isfinite(wald) and wald > 0) else 0.0
    from_value = VALUE_FRACTION * abs(best)
    if from_wald > 0 or from_value > 0:
        step0 = max(from_wald, from_value)
        step_rule = "wald" if from_wald >= from_value else "value"
    else:
        # best is 0 and the Wald error is unusable: nothing about the parameter gives a
        # scale, so the box is all that is left. This is dst_fit's rule, kept only as the
        # last resort and RECORDED so a reader can see it fired.
        step0 = max(abs(hi_b - best), abs(best - lo_b)) / 5000.0
        step_rule = "box_fallback"

    def walk(direction: int):
        """Step outward until F >= f_crit (bracketed) or the box bound is reached."""
        limit = hi_b if direction > 0 else lo_b
        span = abs(limit - best)
        if span <= 0:
            return [], [], limit, True, float("nan"), None
        step = min(step0, span)
        warm = np.array(fit.params, dtype=float)
        vals: list[float] = []
        stats: list[float] = []
        prev_val, prev_warm = best, warm
        for _ in range(MAX_STEPS):
            nxt = prev_val + direction * step
            at_bound = (direction > 0 and nxt >= limit) or (direction < 0 and nxt <= limit)
            if at_bound:
                nxt = limit
            # ``warm`` is the solution at ``prev_val``, so this is the first branch of
            # ``f_best_of(nxt, prev_warm)``; the second branch is only paid for when it
            # can change the answer. Keeping the two identical matters: the root finder
            # below evaluates ``f_best_of(., prev_warm)``, and a bracket found under one
            # objective and rooted under another is not a bracket.
            stat, warm = f_of(nxt, warm)
            if stat >= f_crit:
                # An apparent crossing is CONFIRMED against the second nuisance start
                # before it is believed. A warm-started solve that lands on a worse
                # optimum reports a large SSR, which is indistinguishable from a real
                # rise in the profile -- and manufactures a confidence limit for a
                # parameter that has none. Measured: without this, two of the 100 blood
                # records closed an interval for a Carreau-Yasuda ``eta_0`` whose profile
                # in fact saturates at F = 0.19 and never reaches f_crit anywhere.
                stat_b, warm_b = f_of(nxt, best_params)
                if stat_b < stat:
                    if stat_b < f_crit:
                        false_crossings.append(float(nxt))
                    stat, warm = stat_b, warm_b
            vals.append(nxt)
            stats.append(stat)
            if stat >= f_crit:
                # Bracketed on [prev_val, nxt]. ``g`` warm-starts from the INNER end,
                # which is on the branch the walk followed, and takes the better of that
                # and the best fit -- so it is a function of ``p`` alone. A root finder
                # requires that; a rolling warm start makes the objective path-dependent
                # and brentq then fails on a bracket the walk had genuinely found.
                inner = np.array(prev_warm, dtype=float)

                def g(p: float) -> float:
                    s, _ = f_best_of(p, inner)
                    return s - f_crit

                fa, fb = g(prev_val), g(nxt)
                if not (fa < 0.0 <= fb):
                    # Not a bracket under the deterministic objective, so there is nothing
                    # for a root finder to find. Rejected as unverified, not returned as
                    # an endpoint.
                    return vals, stats, limit, True, stat, (float(nxt), float(fb + f_crit))
                root = float(brentq(g, prev_val, nxt, rtol=ROOT_RTOL, maxiter=200))
                f_root, _ = f_best_of(root, inner)

                # VERIFY. A sign change is not a crossing if the profile jumps across
                # f_crit instead of passing through it, which is what a nuisance
                # minimisation landing on different local optima at neighbouring
                # parameter values produces. Three conditions, all cheap:
                #   (a) F at the located point is f_crit, to VERIFY_TOL;
                #   (b) F is below f_crit just INSIDE it;
                #   (c) F is at or above f_crit just OUTSIDE it.
                # (b) and (c) are what a discontinuity fails: brentq will happily return a
                # point where F is 0.05 f_crit on one side and 10 f_crit on the other.
                # The flanking tests carry the SAME tolerance band as (a). Without it
                # they reject the good case: just outside a genuine crossing,
                # F - f_crit is (dF/dp) * delta, which for a small delta is below the
                # nuisance solve's own noise, so a bare `f_out >= f_crit` fails on
                # rounding. Measured before the band was added: donor B's `eta_inf` was
                # rejected at F = f_crit to 15 digits. The band is far too narrow to
                # admit a discontinuity, which misses by a factor of 4 to 20.
                delta = VERIFY_DELTA * max(abs(root), abs(best))
                f_in, _ = f_best_of(root - direction * delta, inner)
                f_out, _ = f_best_of(root + direction * delta, inner)
                verified = (abs(f_root - f_crit) / f_crit <= VERIFY_TOL
                            and f_in <= f_crit * (1.0 + VERIFY_TOL)
                            and f_out >= f_crit * (1.0 - VERIFY_TOL))
                if not verified:
                    return vals, stats, limit, True, stat, (float(root), float(f_root))
                return vals, stats, root, False, f_root, None
            if at_bound:
                return vals, stats, limit, True, stat, None
            prev_val, prev_warm = nxt, np.array(warm, dtype=float)
            step *= STEP_GROWTH
        # Never crossed and never reached the bound. With a geometrically growing step
        # this needs a surface that is flat over 300 steps; it is reported as OPEN rather
        # than as an interval, because the data did not close it.
        return vals, stats, prev_val, True, (stats[-1] if stats else float("nan")), None

    up_vals, up_stats, upper, upper_open, f_up, up_unver = walk(+1)
    dn_vals, dn_stats, lower, lower_open, f_dn, dn_unver = walk(-1)

    values = np.array(dn_vals[::-1] + [best] + up_vals, dtype=float)
    stats = np.array(dn_stats[::-1] + [0.0] + up_stats, dtype=float)

    return ProfileInterval(
        param=param, best=best, lower=float(lower), upper=float(upper),
        lower_at_bound=bool(lower_open), upper_at_bound=bool(upper_open),
        wald_stderr=wald, f_crit=f_crit, n_points=n, n_params=k,
        bounds=(lo_b, hi_b), f_at_lower=float(f_dn), f_at_upper=float(f_up),
        step_rule=step_rule, n_evals=n_evals,
        worst_negative_f=float(min(0.0, stats.min())) if stats.size else 0.0,
        false_crossings=tuple(false_crossings),
        lower_unverified=dn_unver, upper_unverified=up_unver,
        values=values, stats=stats,
    )


def profile_statistic(model: Model, x: np.ndarray, eta: np.ndarray, param: str,
                      value: float, *, fit: FitResult, level: float = 0.95,
                      warm: np.ndarray | None = None) -> tuple[float, float]:
    """Return ``(F(value), f_crit)`` for one parameter held at ``value``.

    This is what makes an endpoint checkable without trusting the routine that produced
    it: a profile confidence limit is BY DEFINITION the value where ``F = f_crit``, so
    any claimed endpoint -- this module's, ``dst_fit``'s, or one read out of a committed
    sidecar -- can be graded by evaluating the statistic there. Gate G13.2 uses it on the
    shipped sidecar values.

    ``warm`` starts the nuisance minimisation. Left at ``None`` it starts from the best
    fit, which is the honest default for a point supplied from outside: no branch has been
    followed to it, so none can be assumed.
    """
    if param not in model.param_names:
        raise ValueError(f"{model.name} has no parameter {param!r}; has {model.param_names}")
    index = model.param_names.index(param)
    x = np.asarray(x, dtype=float)
    eta = np.asarray(eta, dtype=float)
    good = np.isfinite(x) & np.isfinite(eta) & (x > 0) & (eta > 0)
    x, eta = x[good], eta[good]
    n, k = int(x.size), model.n_params
    dof = n - k
    if dof < 1:
        raise ValueError(f"profile_statistic: {n} points, {k} parameters, dof {dof}")
    minimise = (_direct_ssr(model, x, np.log(eta)) if k == 1
                else _nuisance_minimiser(model, x, np.log(eta), index))
    ssr_min = float(fit.residuals @ fit.residuals)
    start = np.array(fit.params if warm is None else warm, dtype=float)
    ssr, _ = minimise(float(value), start)
    return (ssr - ssr_min) / (ssr_min / dof), float(f_dist.ppf(level, 1, dof))


def newtonian_closed_form(eta: np.ndarray, level: float = 0.95) -> tuple[float, float, float]:
    """The exact profile interval for a one-parameter constant-viscosity fit.

    For ``eta(gd) = eta_0`` with squared residuals in ``log eta`` there are no nuisance
    parameters, so with ``y_i = log eta_i`` and ``u = log eta_0``

        SSR(u) = sum_i (y_i - u)^2,     u_hat = ybar,
        SSR(u) - SSR_min = n (u - ybar)^2

    exactly. With ``s^2 = SSR_min / (n - 1)`` the statistic is ``F = n (u - ybar)^2 / s^2``
    and ``F = f_crit`` at ``u = ybar +/- s sqrt(f_crit / n)``.

    Returns ``(best, lower, upper)`` in viscosity units. This passes through neither the
    outward walk nor the root finder, which is what makes it an anchor (gate G13.4) rather
    than a round trip.
    """
    e = np.asarray(eta, dtype=float)
    e = e[np.isfinite(e) & (e > 0)]
    n = e.size
    if n < 2:
        raise ValueError(f"newtonian_closed_form needs at least 2 points, got {n}")
    y = np.log(e)
    ybar = float(y.mean())
    ssr_min = float(((y - ybar) ** 2).sum())
    s = float(np.sqrt(ssr_min / (n - 1)))
    f_crit = float(f_dist.ppf(level, 1, n - 1))
    half = s * np.sqrt(f_crit / n)
    return float(np.exp(ybar)), float(np.exp(ybar - half)), float(np.exp(ybar + half))
