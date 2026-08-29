"""Does the data actually determine the parameters, and is the model even right?

Two separate questions, and a fit can fail either one while looking fine.

**Identifiability.** A small standard error means "the optimum is sharp along
this axis," not "this parameter is pinned." When two parameters trade off, the
likelihood has a long flat valley and both marginal errors can still come back
small. The correlation matrix is where that shows up, so it is read directly
rather than being summarised into a single number and lost.

**Specification.** A model that cannot make the shape of the data still returns
parameters, and a least-squares fit will report them with error bars. The tell
is not the size of the residuals but their *structure*: a correct model leaves
residuals that flip sign like coin tosses, a wrong one leaves long runs of the
same sign. The Wald-Wolfowitz runs test picks that up without needing to know
the measurement noise, which is the point, because the noise is usually not
known.

Both verdicts are reported as text a person can act on, not as a score.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy.stats import norm

from dst_fit.fitting import FitResult
from dst_fit.models import WC_BETA, frictional_fraction

# Above this, two parameters are trading off badly enough that their individual
# values should not be quoted. 0.95 corresponds to the joint uncertainty region
# being about three times longer than wide.
CORR_DEGENERATE = 0.95

# One-sided p below this means the residuals have significantly too few sign
# changes, i.e. systematic structure the model failed to capture.
RUNS_ALPHA = 0.01

# Amplitude, in log viscosity, below which a *detected* systematic misfit is not
# worth calling a failure of the model.
#
# The runs test answers "is there structure", and its power grows with the number
# of points, so on a precise 30-point sweep it resolves structure far below
# anything that matters. Run on this project's real data it fired on fourteen of
# fifteen curves, which makes it a statement about data quality rather than a
# verdict on the model. What a user needs alongside it is *how big* the
# systematic part is, which is what SYSTEMATIC_TOL is compared against.
#
# 0.05 is five per cent in viscosity, set by the instrument rather than by the
# fits: a rheometer flow curve repeats to a few per cent between loadings, so a
# smooth trend smaller than that is inside the reproducibility of the
# measurement and cannot be attributed to the constitutive model. It is not
# tuned -- see the values it produces in README.md.
SYSTEMATIC_TOL = 0.05

# Degree of the polynomial used to separate a smooth trend from point-to-point
# scatter. Three is the lowest degree that can express the shapes a constitutive
# model actually gets wrong across a sweep -- a tilt, a bow, and an S -- while
# still being far too stiff to chase noise.
SYSTEMATIC_DEGREE = 3


@dataclass(frozen=True)
class Degeneracy:
    a: str
    b: str
    corr: float

    def __str__(self) -> str:
        return f"{self.a} and {self.b} (r = {self.corr:+.3f})"


@dataclass
class Assessment:
    fit: FitResult
    degeneracies: tuple[Degeneracy, ...]
    condition_number: float
    n_runs: int
    runs_z: float
    runs_p: float
    rmse_log: float
    systematic_amplitude: float = 0.0
    systematic_fraction: float = 0.0
    profiles: tuple = ()

    @property
    def identifiable(self) -> bool:
        return len(self.degeneracies) == 0

    @property
    def well_specified(self) -> bool:
        return self.runs_p >= RUNS_ALPHA

    @property
    def structurally_wrong(self) -> bool:
        """Detectably systematic *and* big enough to blame the model for.

        :attr:`well_specified` asks only whether structure is detectable, which
        on a precise sweep it almost always is. This asks the question a user
        has: is the model getting the shape wrong by more than the instrument
        repeats to?
        """
        return not self.well_specified and self.systematic_amplitude > SYSTEMATIC_TOL

    @property
    def usable(self) -> bool:
        """Both questions answered yes. Anything else needs the verdict read."""
        return self.identifiable and self.well_specified

    def verdict(self) -> str:
        name = self.fit.model.name
        lines: list[str] = []

        if self.structurally_wrong:
            lines.append(
                f"{name} does not describe this data. The residuals change sign "
                f"{self.n_runs} times where {self._expected_runs():.0f} would be expected "
                f"by chance (p = {self.runs_p:.2g}), and the systematic part of the "
                f"misfit is {self.systematic_amplitude:.3f} in log viscosity, about "
                f"{100 * self.systematic_amplitude:.0f}% -- far above what the "
                f"measurement repeats to. The parameters below are what least squares "
                f"returned; they are not measurements of anything."
            )
        elif not self.well_specified:
            lines.append(
                f"{name} is slightly off, and this data is precise enough to see it. "
                f"The residuals are systematic (p = {self.runs_p:.2g}, {self.n_runs} sign "
                f"changes against {self._expected_runs():.0f} expected), but the smooth "
                f"part of the misfit is only {self.systematic_amplitude:.3f} in log "
                f"viscosity, roughly {100 * self.systematic_amplitude:.1f}%, which is "
                f"inside the reproducibility of a flow curve. The structure is real; it "
                f"is too small to reject the model over. Treat the parameters as usable "
                f"and the residual shape as worth a look."
            )
        elif self.rmse_log > 0.5:
            lines.append(
                f"{name} fits without obvious structure in the residuals, but loosely "
                f"(rmse {self.rmse_log:.2f} in log viscosity, roughly {100 * self.rmse_log:.0f}% "
                f"scatter). Treat the parameters as indicative."
            )
        else:
            lines.append(
                f"{name} describes this data (rmse {self.rmse_log:.3f} in log viscosity, "
                f"residuals consistent with noise, p = {self.runs_p:.2f})."
            )

        if self.degeneracies:
            worst = self.degeneracies[0]
            lines.append(
                f"The data does not separate {worst}. Their individual values should not "
                f"be quoted; only the combination the data constrains is measured. "
                f"Correlation matrix condition number {self.condition_number:.3g}."
            )
            for d in self.degeneracies[1:]:
                lines.append(f"Also degenerate: {d}.")
        else:
            lines.append(
                f"No parameter pair is degenerate at |r| > {CORR_DEGENERATE:.2f} "
                f"(worst {abs(self.fit.worst_correlation()[2]):.3f}), so the individual "
                f"values carry meaning."
            )

        # The profile is the only thing here that can say a parameter is bounded
        # on one side and not the other. A Wald +/- cannot express that shape at
        # all, and on this model the asymmetry is the normal case rather than the
        # exception: moving phi_m down towards phi is catastrophic and moving it
        # up is nearly free.
        for prof in self.profiles:
            if prof.lower_at_bound and prof.upper_at_bound:
                lines.append(
                    f"{prof.param} is not bounded by this data in either direction: "
                    f"the profile likelihood runs to both parameter bounds without "
                    f"crossing the {100 * prof.level:.0f}% threshold. The fitted "
                    f"{prof.param} = {prof.best:.4g} and its +/- {prof.wald_stderr:.4g} "
                    f"are what least squares returned, not a measurement."
                )
            elif prof.upper_at_bound:
                lines.append(
                    f"{prof.param} is bounded below by {prof.lower:.4g} and unbounded "
                    f"above: the data rules out smaller values and says nothing about "
                    f"larger ones. Quote it as a lower limit, not as "
                    f"{prof.best:.4g} +/- {prof.wald_stderr:.4g}."
                )
            elif prof.lower_at_bound:
                lines.append(
                    f"{prof.param} is bounded above by {prof.upper:.4g} and unbounded "
                    f"below. Quote it as an upper limit, not as "
                    f"{prof.best:.4g} +/- {prof.wald_stderr:.4g}."
                )
            else:
                extra = ""
                if prof.asymmetry >= 1.5:
                    extra = (
                        f" The interval is {prof.asymmetry:.1f}x asymmetric, so the "
                        f"+/- {prof.wald_stderr:.4g} misstates it in one direction."
                    )
                elif prof.widening >= 1.5:
                    extra = (
                        f" That is {prof.widening:.1f}x wider than the linearised "
                        f"+/- {prof.wald_stderr:.4g}."
                    )
                lines.append(
                    f"{prof.param} = {prof.best:.4g}, profile-likelihood interval "
                    f"[{prof.lower:.4g}, {prof.upper:.4g}] at "
                    f"{100 * prof.level:.0f}%.{extra}"
                )
        return "\n".join(lines)

    def _expected_runs(self) -> float:
        r = np.asarray(self.fit.residuals)
        n_pos = int(np.count_nonzero(r > 0))
        n_neg = int(r.size - n_pos)
        if n_pos == 0 or n_neg == 0:
            return 1.0
        return 2.0 * n_pos * n_neg / r.size + 1.0


def systematic_misfit(
    residuals: np.ndarray, degree: int = SYSTEMATIC_DEGREE
) -> tuple[float, float]:
    """Split the residuals into a smooth trend and point-to-point scatter.

    Returns ``(amplitude, fraction)``: the rms of a least-squares polynomial of
    ``degree`` fitted through the residuals in order, and the share of the total
    residual variance that polynomial accounts for.

    This is the *size* of a systematic misfit, where :func:`runs_test` gives only
    its *detectability*. The two answer different questions and a verdict needs
    both. A model that cannot make the shape of the data leaves a large smooth
    trend; a model that is very slightly wrong, measured very precisely, leaves a
    small one that the runs test still rejects.

    The abscissa is the residual index, not the stress: the residuals are already
    ordered by stress and only their ordering matters here, so using the index
    keeps this independent of how the stress axis is spaced.

    A pure-noise sequence still projects onto the polynomial -- with ``n`` points
    and ``degree + 1`` coefficients it keeps about ``(degree + 1) / n`` of the
    variance -- so ``fraction`` has a floor near 0.13 on a 30-point sweep rather
    than near zero. The amplitude, being absolute, does not need that caveat and
    is the number the verdict uses.
    """
    r = np.asarray(residuals, dtype=float)
    r = r[np.isfinite(r)]
    if r.size <= degree + 1:
        return 0.0, 0.0
    total = float(np.sqrt(np.mean(r**2)))
    if total <= 0.0:
        return 0.0, 0.0

    x = np.linspace(-1.0, 1.0, r.size)
    trend = np.polynomial.polynomial.Polynomial.fit(x, r, degree)(x)
    amplitude = float(np.sqrt(np.mean(trend**2)))
    return amplitude, float(min(amplitude**2 / total**2, 1.0))


def runs_test(residuals: np.ndarray) -> tuple[int, float, float]:
    """Wald-Wolfowitz runs test on residual signs, one-sided for *too few* runs.

    Returns ``(n_runs, z, p)``. A large negative z means long same-sign runs,
    which is what a model that cannot make the shape of the data produces. Too
    *many* runs is not tested: that indicates over-fitting or correlated noise,
    which is a different problem and not what this gate is for.

    Residuals must already be ordered by the independent variable, otherwise the
    test measures nothing.
    """
    r = np.asarray(residuals, dtype=float)
    signs = np.sign(r)
    signs = signs[signs != 0]
    n = signs.size
    n_pos = int(np.count_nonzero(signs > 0))
    n_neg = n - n_pos

    if n < 3 or n_pos == 0 or n_neg == 0:
        # No sign variation at all, or too little data to say anything. All one
        # sign is maximally systematic, so report it as such rather than as a
        # pass, but a short array genuinely carries no information.
        if n >= 3 and (n_pos == 0 or n_neg == 0):
            return 1, -np.inf, 0.0
        return max(n, 1), 0.0, 1.0

    n_runs = 1 + int(np.count_nonzero(np.diff(signs) != 0))
    mean = 2.0 * n_pos * n_neg / n + 1.0
    var = (2.0 * n_pos * n_neg * (2.0 * n_pos * n_neg - n)) / (n * n * (n - 1.0))
    if var <= 0:
        return n_runs, 0.0, 1.0
    z = (n_runs - mean) / np.sqrt(var)
    return n_runs, float(z), float(norm.cdf(z))


@dataclass(frozen=True)
class TransitionCoverage:
    """Where the fitted sigma* sits relative to the stresses actually measured."""

    sigma_star: float
    sigma_min: float
    sigma_max: float
    f_min: float
    f_max: float
    decades_below: float
    decades_above: float

    @property
    def brackets(self) -> bool:
        """True when the data spans the frictional transition rather than one side of it."""
        return self.f_min < 0.25 and self.f_max > 0.75

    def __str__(self) -> str:
        if self.brackets:
            return (
                f"the measured stresses bracket the transition "
                f"(frictional fraction runs {self.f_min:.2f} to {self.f_max:.2f})"
            )
        # Deliberately no "measure N more decades" advice. That number would be
        # computed from the fitted sigma*, which in this branch is exactly what
        # the data does not determine, so it would be circular: the estimate
        # setting the scale for the measurement meant to constrain the estimate.
        #
        # Note also what the fit does here rather than what one might expect.
        # sigma* does not run off to a bound; it migrates *into* the measured
        # window. On synthetic data with a true sigma* of 5 Pa, a sweep over
        # 1e-4 to 1e-2 Pa returns 0.006 and a sweep over 1e3 to 1e5 Pa returns
        # 5.3e4. So a fitted sigma* sitting comfortably inside the measured range
        # is not evidence of anything; it is what happens either way.
        side = "below" if self.f_max <= 0.75 else "above"
        return (
            f"the measured stresses sit {side} the transition rather than across it: "
            f"the frictional fraction only runs {self.f_min:.2f} to {self.f_max:.2f} over "
            f"{self.sigma_min:.3g} to {self.sigma_max:.3g} Pa. sigma* is not measured here. "
            f"The fit will still place it inside this window, because that is what "
            f"least squares does with an unconstrained parameter, so its apparent "
            f"plausibility carries no information. Extend the sweep to "
            f"{'higher' if side == 'below' else 'lower'} stress"
        )


def transition_coverage(
    sigma_star: float, stress: np.ndarray, beta: float = WC_BETA
) -> TransitionCoverage:
    """Does the measured stress range actually span the frictional transition?

    Wyart-Cates puts the whole of its stress dependence in ``f(sigma) =
    exp(-sigma*/sigma)``, which goes from 0 to 1 over roughly two decades around
    ``sigma*``. If every measured point sits on one side of that, the fit still
    returns a ``sigma*`` with a standard error, and the number is a bound rather
    than a measurement: nothing in the data distinguishes it from any other value
    further out.

    This is a property of the *experiment*, not of the fit, so it is worth saying
    separately. It is also actionable, which a covariance matrix is not: the
    answer is to measure further, and this says how much further.
    """
    s = np.asarray(stress, dtype=float)
    s = s[np.isfinite(s) & (s > 0)]
    if s.size == 0:
        raise ValueError("no positive finite stresses")
    lo, hi = float(s.min()), float(s.max())
    f_lo = float(frictional_fraction(np.array([lo]), sigma_star, beta=beta)[0])
    f_hi = float(frictional_fraction(np.array([hi]), sigma_star, beta=beta)[0])

    # f = 0.25 at sigma = sigma*/ln(4)^(1/beta); f = 0.75 at sigma*/ln(4/3)^(1/beta).
    s_25 = sigma_star / np.log(4.0) ** (1.0 / beta)
    s_75 = sigma_star / np.log(4.0 / 3.0) ** (1.0 / beta)
    return TransitionCoverage(
        sigma_star=float(sigma_star),
        sigma_min=lo,
        sigma_max=hi,
        f_min=f_lo,
        f_max=f_hi,
        decades_below=float(max(np.log10(lo / s_25), 0.0)),
        decades_above=float(max(np.log10(s_75 / hi), 0.0)),
    )


@dataclass(frozen=True)
class ConsistencyCheck:
    """Does one parameter hold still across a series of independent fits?

    ``phi_m`` is a property of the particles, so N fits of one material at N
    volume fractions have to return one value. This is the statistic that asks,
    and it is here rather than in a gate script because two gates must apply an
    identical rule to two datasets for the comparison between them to mean
    anything. A rule copied between scripts drifts.
    """

    values: np.ndarray
    errors: np.ndarray
    mean: float
    mean_err: float
    deviations: np.ndarray  # |value - mean| / error, in sigma
    n_sigma: float
    n_consistent: int
    n_total: int
    chi2: float
    dof: int
    n_pairs_consistent: int
    n_pairs: int

    @property
    def reduced_chi2(self) -> float:
        return self.chi2 / self.dof if self.dof > 0 else float("nan")

    @property
    def unweighted_mean(self) -> float:
        return float(np.mean(self.values))

    @property
    def median(self) -> float:
        return float(np.median(self.values))

    def passes(self, n_required: int) -> bool:
        return self.n_consistent >= n_required

    def report(self, labels: list[str] | None = None) -> str:
        labels = labels or [f"curve {i + 1}" for i in range(self.n_total)]
        lines = [
            f"against the inverse-variance weighted mean {self.mean:.4f} "
            f"+/- {self.mean_err:.4f}"
        ]
        for lab, v, e, d in zip(labels, self.values, self.errors, self.deviations):
            mark = "ok" if d <= self.n_sigma else "OUT"
            lines.append(f"    {lab:>16s}   {v:.4f} +/- {e:.4f}   {d:7.1f} sigma   {mark}")
        lines.append(
            f"    {self.n_consistent} of {self.n_total} within {self.n_sigma:.0f} sigma"
        )
        lines.append(
            f"    chi2 for a single common value: {self.chi2:.1f} on {self.dof} dof "
            f"(reduced {self.reduced_chi2:.1f})"
        )
        lines.append(
            f"    pairwise: {self.n_pairs_consistent} of {self.n_pairs} pairs agree "
            f"within {self.n_sigma:.0f} sigma"
        )
        # The weighted mean is pulled hard by whichever curve has the smallest
        # error bar, and near jamming those errors are linearised and optimistic.
        # The unweighted centre says where the estimates actually sit.
        lines.append(
            f"    unweighted mean {self.unweighted_mean:.4f}, median {self.median:.4f}, "
            f"range {float(np.min(self.values)):.4f} to {float(np.max(self.values)):.4f}"
        )
        return "\n".join(lines)


def weighted_mean(values: np.ndarray, errors: np.ndarray) -> tuple[float, float]:
    """Inverse-variance weighted mean and its standard error."""
    values = np.asarray(values, dtype=float)
    errors = np.asarray(errors, dtype=float)
    ok = np.isfinite(values) & np.isfinite(errors) & (errors > 0)
    if not ok.any():
        return float("nan"), float("nan")
    w = 1.0 / errors[ok] ** 2
    mean = float(np.sum(w * values[ok]) / np.sum(w))
    return mean, float(np.sqrt(1.0 / np.sum(w)))


def consistency_across_curves(
    values: np.ndarray, errors: np.ndarray, n_sigma: float = 2.0
) -> ConsistencyCheck:
    """Test whether N independent estimates of one constant agree.

    Consistency is measured against the inverse-variance weighted mean of all N
    estimates, and an estimate counts as consistent when it sits within
    ``n_sigma`` of that mean *in units of its own error*. The alternative
    reading, all N(N-1)/2 pairwise comparisons, is stricter and is computed
    alongside for information rather than being chosen after the fact.

    Both readings are reported because they can disagree, and which one a gate
    uses has to be fixed before the data is seen.
    """
    values = np.asarray(values, dtype=float)
    errors = np.asarray(errors, dtype=float)
    if values.shape != errors.shape:
        raise ValueError(f"shape mismatch: values {values.shape} vs errors {errors.shape}")
    if values.size < 2:
        raise ValueError("consistency needs at least two estimates")

    mean, mean_err = weighted_mean(values, errors)
    with np.errstate(divide="ignore", invalid="ignore"):
        dev = np.abs(values - mean) / errors
    n_consistent = int(np.count_nonzero(dev <= n_sigma))
    chi2 = float(np.nansum(((values - mean) / errors) ** 2))

    pairs = [
        abs(values[i] - values[j]) / np.hypot(errors[i], errors[j])
        for i in range(values.size)
        for j in range(i + 1, values.size)
    ]
    pairs_arr = np.asarray(pairs, dtype=float)
    return ConsistencyCheck(
        values=values,
        errors=errors,
        mean=mean,
        mean_err=mean_err,
        deviations=dev,
        n_sigma=n_sigma,
        n_consistent=n_consistent,
        n_total=int(values.size),
        chi2=chi2,
        dof=int(values.size - 1),
        n_pairs_consistent=int(np.count_nonzero(pairs_arr <= n_sigma)),
        n_pairs=int(pairs_arr.size),
    )


def pairwise_degeneracies(
    fit: FitResult, corr_max: float = CORR_DEGENERATE
) -> tuple[Degeneracy, ...]:
    """Parameter pairs correlated beyond ``corr_max``, worst first."""
    corr = np.asarray(fit.corr)
    names = fit.model.param_names
    found = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            r = corr[i, j]
            if np.isfinite(r) and abs(r) > corr_max:
                found.append(Degeneracy(names[i], names[j], float(r)))
    return tuple(sorted(found, key=lambda d: -abs(d.corr)))


def assess(
    fit: FitResult,
    corr_max: float = CORR_DEGENERATE,
    *,
    stress: np.ndarray | None = None,
    viscosity: np.ndarray | None = None,
    profile: Sequence[str] = (),
    level: float = 0.95,
) -> Assessment:
    """Full verdict on one fit: is it identifiable, and is the model right?

    Parameters
    ----------
    profile
        Names of parameters to profile. Empty by default because a profile costs
        one nuisance minimisation per step and ``assess`` is called in loops;
        naming a parameter here is what turns the verdict from "these two trade
        off" into "this one is bounded below and not above", which is the thing
        a user can act on. Requires ``stress`` and ``viscosity``: a
        :class:`~dst_fit.fitting.FitResult` carries its residuals but not the
        data they were computed from.
    """
    corr = np.asarray(fit.corr)
    if np.all(np.isfinite(corr)):
        cond = float(np.linalg.cond(corr))
    else:
        cond = float("inf")

    n_runs, z, p = runs_test(fit.residuals)
    amplitude, fraction = systematic_misfit(fit.residuals)

    profiles: list = []
    if profile:
        if stress is None or viscosity is None:
            raise ValueError(
                "profiling needs the data: pass stress= and viscosity= alongside "
                f"profile={list(profile)!r}. A FitResult carries its residuals but "
                "not the curve they came from."
            )
        # Imported here rather than at module scope: profiling imports fitting,
        # and keeping this local documents that the profile is an opt-in extra
        # rather than part of the cheap path every caller takes.
        from dst_fit.profiling import profile_likelihood

        for name in profile:
            profiles.append(
                profile_likelihood(
                    fit.model,
                    stress,
                    viscosity,
                    name,
                    branch=fit.branch,
                    level=level,
                    fit=fit,
                )
            )

    return Assessment(
        fit=fit,
        degeneracies=pairwise_degeneracies(fit, corr_max),
        condition_number=cond,
        n_runs=n_runs,
        runs_z=z,
        runs_p=p,
        rmse_log=fit.rmse_log,
        systematic_amplitude=amplitude,
        systematic_fraction=fraction,
        profiles=tuple(profiles),
    )
