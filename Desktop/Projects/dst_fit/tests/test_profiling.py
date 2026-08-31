"""Profile likelihood, block bootstrap, and the interval-intersection test.

The three exist to answer one question -- were the error bars too small? -- so
the tests check the properties that would make the answer wrong: that the
machinery reproduces the Wald interval where the Wald interval is right, widens
where it is wrong, reports unboundedness rather than inventing a limit, and
intersects intervals without a Gaussian assumption.
"""

from __future__ import annotations

import numpy as np
import pytest

from dst_fit import fit_model
from dst_fit.fitting import synthetic_curve
from dst_fit.models import WYART_CATES
from dst_fit.profiling import (
    block_bootstrap,
    block_length_from_runs,
    interval_intersection,
    profile_likelihood,
    wald_width,
)

TRUTH = (0.50, 0.58, 5.0)
GRID = np.logspace(-1, 3, 50)


def _curve(seed=3, noise=0.03, truth=TRUTH):
    return synthetic_curve(WYART_CATES, truth, GRID, noise=noise, seed=seed)


# --------------------------------------------------------------------------
# profile likelihood
# --------------------------------------------------------------------------
def test_profile_brackets_the_optimum_and_contains_the_truth():
    s, eta = _curve()
    fit = fit_model(WYART_CATES, s, eta, n_starts=8)
    pr = profile_likelihood(WYART_CATES, s, eta, "phi_m", fit=fit)
    assert pr.lower < pr.best < pr.upper
    assert pr.contains(TRUTH[1])
    assert pr.bounded


def test_profile_agrees_with_wald_on_a_quadratic_surface():
    """Where the linearisation is right, the profile must not invent width.

    This is the control that makes a *disagreement* meaningful. A method that
    always came out wider would flag every fit and mean nothing.
    """
    s, eta = _curve()
    fit = fit_model(WYART_CATES, s, eta, n_starts=8)
    pr = profile_likelihood(WYART_CATES, s, eta, "phi_m", fit=fit)
    assert 0.9 <= pr.widening <= 1.3
    assert pr.asymmetry < 1.5


def test_profile_statistic_is_zero_at_the_optimum_and_rises_away_from_it():
    s, eta = _curve()
    pr = profile_likelihood(WYART_CATES, s, eta, "phi_m")
    j = int(np.argmin(np.abs(pr.values - pr.best)))
    assert pr.f_stat[j] == pytest.approx(0.0, abs=1e-9)
    assert np.all(pr.f_stat >= -1e-9), "the profiled SSR cannot beat the full optimum"
    # the endpoints of the walk are where the statistic reached its critical value
    assert pr.f_stat.max() >= pr.f_crit * 0.5


def test_profile_reports_unboundedness_instead_of_inventing_a_limit():
    """A sweep entirely below the transition cannot bound sigma* from above.

    G5 established that the *fit* still returns a plausible-looking sigma* in
    this case. The profile is what says the number is a bound, not a
    measurement, and the flag is the product feature.
    """
    s, eta = synthetic_curve(WYART_CATES, TRUTH, np.logspace(-4, -2, 40),
                             noise=0.02, seed=1)
    fit = fit_model(WYART_CATES, s, eta, n_starts=12)
    pr = profile_likelihood(WYART_CATES, s, eta, "sigma_star", fit=fit)
    assert not pr.bounded
    assert "unbounded" in pr.summary().lower() or "UNBOUNDED" in pr.summary()


def test_profile_refuses_an_unknown_parameter():
    s, eta = _curve()
    with pytest.raises(ValueError, match="no parameter"):
        profile_likelihood(WYART_CATES, s, eta, "phi_max")


def test_profile_refuses_a_noiseless_curve():
    """With residuals at float64 rounding there is no scale to compare against.

    The noiseless fit reaches rmse ~7e-16, not exactly zero, so the guard has to
    be on the residual *scale* rather than on SSR being identically zero.
    """
    s, eta = synthetic_curve(WYART_CATES, TRUTH, GRID, noise=0.0)
    with pytest.raises(ValueError, match="float64 rounding"):
        profile_likelihood(WYART_CATES, s, eta, "phi_m")


def test_profile_level_widens_the_interval():
    s, eta = _curve()
    fit = fit_model(WYART_CATES, s, eta, n_starts=8)
    narrow = profile_likelihood(WYART_CATES, s, eta, "phi_m", fit=fit, level=0.68)
    wide = profile_likelihood(WYART_CATES, s, eta, "phi_m", fit=fit, level=0.99)
    assert wide.width > narrow.width


# --------------------------------------------------------------------------
# block length
# --------------------------------------------------------------------------
def test_block_length_is_two_for_independent_residuals():
    """Independent signs give ~n/2 runs, hence blocks of 2."""
    rng = np.random.default_rng(0)
    r = rng.standard_normal(400)
    L, n_runs = block_length_from_runs(r)
    assert 180 < n_runs < 220
    assert L == 2


def test_block_length_grows_with_the_length_of_the_runs():
    """Four long runs in 40 points must give blocks of about 10."""
    r = np.concatenate([np.full(10, 1.0), np.full(10, -1.0),
                        np.full(10, 1.0), np.full(10, -1.0)])
    L, n_runs = block_length_from_runs(r)
    assert n_runs == 4
    assert L == 10


# --------------------------------------------------------------------------
# block bootstrap
# --------------------------------------------------------------------------
def test_bootstrap_contains_the_truth_on_correctly_specified_data():
    s, eta = _curve()
    bs = block_bootstrap(WYART_CATES, s, eta, "phi_m", n_resamples=200, seed=1)
    assert bs.contains(TRUTH[1])
    assert bs.n_converged >= 150


def test_bootstrap_is_not_wide_by_construction():
    """It must reproduce the Wald width when the residuals really are independent.

    Otherwise 'the bootstrap interval is wider' would carry no information: the
    method has to be able to agree.
    """
    s, eta = _curve()
    fit = fit_model(WYART_CATES, s, eta, n_starts=8)
    iid = block_bootstrap(WYART_CATES, s, eta, "phi_m", fit=fit,
                          n_resamples=300, block_length=1, seed=1)
    assert 0.7 <= iid.widening <= 1.4


def test_correlated_residuals_widen_the_bootstrap():
    """The claim the method exists to make, on a curve with a planted misfit.

    Wyart-Cates is fit to data that also carries a smooth systematic bump, so
    the residuals come out in long runs. Blocks sized from those runs must give
    a wider interval than blocks of 1.
    """
    s, eta_clean = synthetic_curve(WYART_CATES, TRUTH, GRID, noise=0.005, seed=2)
    bump = 1.0 + 0.10 * np.sin(np.pi * np.log10(s / s.min()) / np.log10(s.max() / s.min()))
    eta = eta_clean * bump

    fit = fit_model(WYART_CATES, s, eta, n_starts=12)
    blocked = block_bootstrap(WYART_CATES, s, eta, "phi_m", fit=fit,
                              n_resamples=300, seed=1)
    iid = block_bootstrap(WYART_CATES, s, eta, "phi_m", fit=fit,
                          n_resamples=300, block_length=1, seed=1)
    assert blocked.block_length > 2, "the planted bump must show up as long runs"
    assert blocked.width > iid.width


def test_bootstrap_refuses_a_fit_made_on_other_data():
    """Silently pairing a fit with the wrong data would resample the wrong residuals."""
    s, eta = _curve()
    fit = fit_model(WYART_CATES, s[:30], eta[:30], n_starts=4)
    with pytest.raises(ValueError, match="not the same fit|was not made on this data"):
        block_bootstrap(WYART_CATES, s, eta, "phi_m", fit=fit, n_resamples=30)


def test_bootstrap_refuses_too_few_points():
    s, eta = _curve()
    with pytest.raises(ValueError, match="too few"):
        block_bootstrap(WYART_CATES, s[:6], eta[:6], "phi_m", n_resamples=30)


# --------------------------------------------------------------------------
# intersecting intervals
# --------------------------------------------------------------------------
class _FakeInterval:
    def __init__(self, lower, upper):
        self.lower, self.upper, self.param = lower, upper, "phi_m"


def test_overlapping_intervals_admit_a_common_value():
    ivs = [_FakeInterval(0.55, 0.57), _FakeInterval(0.56, 0.59), _FakeInterval(0.54, 0.565)]
    out = interval_intersection(ivs, ["a", "b", "c"])
    assert out.non_empty
    assert out.lower == pytest.approx(0.56)
    assert out.upper == pytest.approx(0.565)
    assert out.limiting_low == "b" and out.limiting_high == "c"
    assert "one common value is possible" in out.report()


def test_disjoint_intervals_report_the_largest_agreeing_subset():
    """Three that agree and two that do not must report 3, not 5 and not 0."""
    ivs = [_FakeInterval(0.50, 0.52), _FakeInterval(0.505, 0.525),
           _FakeInterval(0.51, 0.53), _FakeInterval(0.58, 0.60),
           _FakeInterval(0.585, 0.605)]
    out = interval_intersection(ivs, list("abcde"))
    assert not out.non_empty
    assert out.n_overlapping == 3
    assert out.n_total == 5
    assert "no single value fits all 5" in out.report()


def test_intersection_needs_two_intervals():
    with pytest.raises(ValueError, match="at least two"):
        interval_intersection([_FakeInterval(0.5, 0.6)])


def test_intersection_of_real_profiles_from_one_curve_is_trivially_non_empty():
    """Sanity: profiling the same curve twice must not manufacture disagreement."""
    s, eta = _curve()
    fit = fit_model(WYART_CATES, s, eta, n_starts=8)
    a = profile_likelihood(WYART_CATES, s, eta, "phi_m", fit=fit)
    b = profile_likelihood(WYART_CATES, s, eta, "phi_m", fit=fit)
    assert interval_intersection([a, b], ["a", "b"]).non_empty


# --------------------------------------------------------------------------- #
# One Wald ruler.  Gate N7's argument is a comparison between the profile's
# widening and the bootstrap's, so the two had to stop measuring themselves
# against different Wald intervals.
# --------------------------------------------------------------------------- #
def test_wald_width_matches_the_t_quantile_by_hand():
    """2 * t_{27, 0.975} * se, computed from a table value, not from the code."""
    # t_{27, 0.975} = 2.051831, standard table.  se chosen as an exact binary
    # fraction so the comparison is not about float formatting.
    assert wald_width(0.25, dof=27, level=0.95) == pytest.approx(
        2.0 * 2.0518305 * 0.25, rel=1e-6
    )


def test_wald_width_tends_to_the_normal_quantile_at_large_dof():
    """The finite-sample ruler must reduce to the asymptotic one, not sit beside it."""
    assert wald_width(1.0, dof=10_000_000, level=0.95) == pytest.approx(
        2.0 * 1.959963984540054, rel=1e-4
    )
    # ...and be strictly wider at every finite dof, which is the whole reason
    # the two rulers disagreed.
    for dof in (5, 10, 26, 27, 100):
        assert wald_width(1.0, dof=dof, level=0.95) > 2.0 * 1.959963984540054


def test_wald_width_refuses_zero_degrees_of_freedom():
    with pytest.raises(ValueError, match="dof"):
        wald_width(1.0, dof=0)


def test_profile_and_bootstrap_use_the_same_ruler():
    """Same curve, same parameter: the two widenings must share a denominator."""
    s, eta = _curve()
    fit = fit_model(WYART_CATES, s, eta, n_starts=8)
    prof = profile_likelihood(WYART_CATES, s, eta, "phi_m", fit=fit)
    boot = block_bootstrap(WYART_CATES, s, eta, "phi_m", fit=fit, n_resamples=60)
    assert prof.dof == boot.dof
    assert prof.wald_width == pytest.approx(boot.wald_width, rel=1e-12)


# --------------------------------------------------------------------------- #
# Fitted residuals are shrunk, and a resample drawn from them inherits it.
# --------------------------------------------------------------------------- #
def test_rescaling_widens_the_bootstrap_by_the_predicted_factor():
    """sqrt(n / (n - k)) -- an analytic factor, not a tuned one.

    Discriminating: the *unrescaled* spread must be narrower by that factor, so
    a version of the code that rescaled twice, or not at all, fails this.
    """
    s, eta = _curve()
    fit = fit_model(WYART_CATES, s, eta, n_starts=8)
    kw = dict(fit=fit, block_length=1, n_resamples=400, seed=11)
    on = block_bootstrap(WYART_CATES, s, eta, "phi_m", **kw)
    off = block_bootstrap(WYART_CATES, s, eta, "phi_m", rescale_residuals=False, **kw)
    n, k = on.n_points, on.n_free
    predicted = np.sqrt(n / (n - k))
    ratio = np.std(on.values, ddof=1) / np.std(off.values, ddof=1)
    assert ratio == pytest.approx(predicted, rel=0.05)
    assert predicted > 1.0  # guards against a degenerate n == k slipping through


def test_rescaled_bootstrap_dispersion_matches_the_wald_standard_error():
    """Under iid noise the two estimate the same sigma; before rescaling they did not."""
    s, eta = _curve()
    fit = fit_model(WYART_CATES, s, eta, n_starts=8)
    i = WYART_CATES.param_names.index("phi_m")
    se = float(fit.stderr[i])
    boot = block_bootstrap(
        WYART_CATES, s, eta, "phi_m", fit=fit, block_length=1, n_resamples=600, seed=5
    )
    assert np.std(boot.values, ddof=1) / se == pytest.approx(1.0, abs=0.10)


def test_warm_start_is_not_what_sets_the_interval_width():
    """H_cold, refuted.  Kept as a test so the claim cannot quietly reverse.

    The session plan blamed the iid control's sub-Wald width on resamples being
    refitted single-start from a generic guess.  Warm-starting them changes the
    interval by well under the 5% the shortfall was; the shortfall is arithmetic
    (see scripts/probe_bootstrap_calibration.py).
    """
    s, eta = _curve()
    fit = fit_model(WYART_CATES, s, eta, n_starts=8)
    kw = dict(fit=fit, block_length=1, n_resamples=300, seed=7)
    warm = block_bootstrap(WYART_CATES, s, eta, "phi_m", **kw)
    cold = block_bootstrap(
        WYART_CATES, s, eta, "phi_m", warm_start=False, boot_starts=1, **kw
    )
    assert warm.width == pytest.approx(cold.width, rel=0.05)
