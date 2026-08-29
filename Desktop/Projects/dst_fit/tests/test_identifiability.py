"""Tests for the identifiability and specification verdicts."""

from __future__ import annotations

import numpy as np
import pytest

from dst_fit.fitting import FitResult, fit_model, synthetic_curve
from dst_fit.identifiability import (
    transition_coverage,
    CORR_DEGENERATE,
    assess,
    pairwise_degeneracies,
    runs_test,
    systematic_misfit,
)
from dst_fit.models import WYART_CATES

SIGMA = np.logspace(-2, 3, 40)
TRUTH_WC = (0.52, 0.58, 5.0)
TRUTH_HB = (2.0, 1.5, 0.6)


# --------------------------------------------------------------------------
# runs test
# --------------------------------------------------------------------------
def test_alternating_residuals_look_like_noise():
    r = np.array([1.0, -1.0] * 20)
    n_runs, z, p = runs_test(r)
    assert n_runs == 40
    assert z > 0 and p > 0.99  # too many runs is not what this test flags


def test_one_long_block_of_each_sign_is_flagged():
    r = np.concatenate([np.ones(20), -np.ones(20)])
    n_runs, z, p = runs_test(r)
    assert n_runs == 2
    assert z < -3 and p < 1e-3


def test_all_one_sign_is_maximally_systematic():
    n_runs, z, p = runs_test(np.ones(20))
    assert n_runs == 1 and p == 0.0


def test_random_signs_are_not_flagged():
    rng = np.random.default_rng(0)
    fired = sum(runs_test(rng.normal(size=60))[2] < 0.01 for _ in range(200))
    assert fired <= 6, f"{fired}/200 false positives at alpha 0.01"


def test_too_short_to_say_anything_returns_neutral():
    _, z, p = runs_test(np.array([1.0, -1.0]))
    assert z == 0.0 and p == 1.0


def test_zeros_are_dropped_not_counted_as_a_sign():
    n_runs, _, _ = runs_test(np.array([1.0, 0.0, 1.0, -1.0]))
    assert n_runs == 2  # the zero does not create a spurious run


# --------------------------------------------------------------------------
# degeneracy
# --------------------------------------------------------------------------
def test_degeneracies_are_sorted_worst_first():
    s, e = synthetic_curve("ellis", (10.0, 0.05, 3.0, 1.4), SIGMA, noise=0.01, seed=0)
    fit = fit_model("ellis", s, e)
    fit.corr = np.array(
        [
            [1.00, 0.99, 0.10, 0.20],
            [0.99, 1.00, -0.97, 0.30],
            [0.10, -0.97, 1.00, 0.40],
            [0.20, 0.30, 0.40, 1.00],
        ]
    )
    d = pairwise_degeneracies(fit)
    assert [abs(x.corr) for x in d] == sorted([abs(x.corr) for x in d], reverse=True)
    assert len(d) == 2


def test_threshold_is_respected():
    s, e = synthetic_curve(WYART_CATES, TRUTH_WC, SIGMA, noise=0.02, seed=1)
    fit = fit_model(WYART_CATES, s, e)
    _, _, worst = fit.worst_correlation()
    assert pairwise_degeneracies(fit, corr_max=abs(worst) + 0.01) == ()
    assert len(pairwise_degeneracies(fit, corr_max=abs(worst) - 0.01)) >= 1


def test_default_threshold_is_the_documented_one():
    assert CORR_DEGENERATE == 0.95


# --------------------------------------------------------------------------
# assess
# --------------------------------------------------------------------------
def test_correct_model_is_usable():
    s, e = synthetic_curve(WYART_CATES, TRUTH_WC, SIGMA, noise=0.02, seed=0)
    a = assess(fit_model(WYART_CATES, s, e))
    assert a.well_specified and a.identifiable and a.usable
    assert "describes this data" in a.verdict()


def test_wrong_model_is_caught_and_says_so_plainly():
    s, e = synthetic_curve(WYART_CATES, TRUTH_WC, SIGMA, noise=0.02, seed=0)
    a = assess(fit_model("herschel_bulkley", s, e, n_starts=6))
    assert not a.well_specified
    assert not a.usable
    verdict = a.verdict()
    assert "does not describe this data" in verdict
    assert "not measurements of anything" in verdict


def test_wrong_model_in_the_other_direction_is_also_caught():
    s, e = synthetic_curve("herschel_bulkley", TRUTH_HB, SIGMA, noise=0.02, seed=0)
    a = assess(fit_model(WYART_CATES, s, e, n_starts=6))
    assert not a.usable


def test_condition_number_is_finite_for_a_healthy_fit():
    s, e = synthetic_curve(WYART_CATES, TRUTH_WC, SIGMA, noise=0.02, seed=0)
    a = assess(fit_model(WYART_CATES, s, e))
    assert np.isfinite(a.condition_number) and a.condition_number >= 1.0


def test_verdict_names_the_degenerate_pair():
    s, e = synthetic_curve(WYART_CATES, TRUTH_WC, SIGMA, noise=0.02, seed=0)
    fit = fit_model(WYART_CATES, s, e)
    fit.corr = np.array([[1.0, 0.999, 0.1], [0.999, 1.0, 0.1], [0.1, 0.1, 1.0]])
    a = assess(fit)
    assert not a.identifiable
    v = a.verdict()
    assert "does not separate" in v and "phi" in v and "phi_m" in v
    assert "should not be quoted" in v


def test_non_finite_correlation_reports_infinite_condition_number():
    s, e = synthetic_curve(WYART_CATES, TRUTH_WC, SIGMA, noise=0.02, seed=0)
    fit = fit_model(WYART_CATES, s, e)
    fit.corr = np.full((3, 3), np.nan)
    assert assess(fit).condition_number == pytest.approx(float("inf"))


# --------------------------------------------------------------------------
# transition coverage
# --------------------------------------------------------------------------
def test_bracketing_window_is_reported_as_bracketing():
    cov = transition_coverage(5.0, np.logspace(-1, 3, 25))
    assert cov.brackets
    assert cov.f_min < 0.25 and cov.f_max > 0.75
    assert "bracket the transition" in str(cov)


def test_window_below_the_transition_is_refused():
    cov = transition_coverage(5.0, np.logspace(-4, -2, 25))
    assert not cov.brackets
    assert "sit below the transition" in str(cov)
    assert "Extend the sweep to higher stress" in str(cov)


def test_window_above_the_transition_is_refused():
    cov = transition_coverage(5.0, np.logspace(4, 6, 25))
    assert not cov.brackets
    assert "sit above the transition" in str(cov)
    assert "Extend the sweep to lower stress" in str(cov)


def test_message_never_quotes_a_decade_count_when_unconstrained():
    """That number would come from the sigma* the data failed to determine."""
    for stress in (np.logspace(-4, -2, 25), np.logspace(4, 6, 25)):
        msg = str(transition_coverage(5.0, stress))
        assert "more decades" not in msg
        assert "carries no information" in msg


def test_coverage_uses_only_positive_finite_stresses():
    cov = transition_coverage(5.0, np.array([np.nan, -1.0, 0.0, 1.0, 100.0]))
    assert cov.sigma_min == 1.0 and cov.sigma_max == 100.0


def test_coverage_rejects_an_empty_window():
    with pytest.raises(ValueError, match="no positive finite"):
        transition_coverage(5.0, np.array([0.0, -1.0, np.nan]))


def test_sigma_star_migrates_into_the_measured_window():
    """The failure mode G5 exposed: an unconstrained sigma* lands inside the sweep.

    Regression guard. If a future change made sigma* run to a bound instead, the
    coverage check would still be right but the README's explanation would not.
    """
    for lo, hi in ((-4.0, -2.0), (3.0, 5.0)):
        window = np.logspace(lo, hi, 25)
        s, e = synthetic_curve(WYART_CATES, TRUTH_WC, window, noise=0.02, seed=0)
        fit = fit_model(WYART_CATES, s, e, n_starts=8, seed=0)
        ss = fit.params[WYART_CATES.param_names.index("sigma_star")]
        assert s.min() <= ss <= s.max(), f"sigma*={ss:g} outside [{s.min():g}, {s.max():g}]"
        assert not transition_coverage(ss, s).brackets


# --------------------------------------------------------------------------
# separating a detectable misfit from a big one
# --------------------------------------------------------------------------
def test_systematic_misfit_recovers_a_planted_trend():
    """A smooth bow plus small noise: the amplitude must be the bow's."""
    rng = np.random.default_rng(0)
    x = np.linspace(-1.0, 1.0, 40)
    trend = 0.10 * (1.0 - x**2)
    resid = trend + 0.005 * rng.standard_normal(x.size)
    amp, frac = systematic_misfit(resid)
    np.testing.assert_allclose(amp, float(np.sqrt(np.mean(trend**2))), rtol=0.05)
    assert frac > 0.95


def test_systematic_misfit_of_pure_noise_is_small_but_not_zero():
    """A cubic through noise keeps about (degree + 1) / n of the variance."""
    rng = np.random.default_rng(1)
    resid = rng.standard_normal(40)
    amp, frac = systematic_misfit(resid)
    assert amp < 0.6  # against a unit-sigma sequence
    assert 0.02 < frac < 0.35


def test_systematic_misfit_handles_degenerate_input():
    assert systematic_misfit(np.zeros(20)) == (0.0, 0.0)
    assert systematic_misfit(np.array([1.0, 2.0])) == (0.0, 0.0)


def test_systematic_misfit_is_scale_equivariant():
    """Doubling the residuals doubles the amplitude and leaves the fraction alone."""
    rng = np.random.default_rng(2)
    resid = 0.02 * rng.standard_normal(50) + np.linspace(-0.05, 0.05, 50)
    a1, f1 = systematic_misfit(resid)
    a2, f2 = systematic_misfit(2.0 * resid)
    np.testing.assert_allclose(a2, 2.0 * a1, rtol=1e-12)
    np.testing.assert_allclose(f2, f1, rtol=1e-12)


def _fit_with_residuals(resid):
    """A minimal FitResult standing in for a real fit, to drive Assessment."""
    n = resid.size
    return FitResult(
        model=WYART_CATES,
        params=np.array([0.40, 0.52, 5.0]),
        stderr=np.array([0.01, 0.01, 0.5]),
        cov=np.eye(3) * 1e-4,
        corr=np.eye(3),
        residuals=np.asarray(resid, dtype=float),
        n_points=n,
        success=True,
        message="synthetic",
    )


def test_verdict_separates_a_big_misfit_from_a_small_detectable_one():
    """The 14-of-15 problem: both of these have systematic residuals."""
    x = np.linspace(-1.0, 1.0, 40)
    big = assess(_fit_with_residuals(0.30 * (1.0 - x**2)))
    small = assess(_fit_with_residuals(0.010 * (1.0 - x**2)))

    assert not big.well_specified and not small.well_specified
    assert big.structurally_wrong
    assert not small.structurally_wrong
    assert "does not describe this data" in big.verdict()
    assert "slightly off" in small.verdict()
    assert "too small to reject the model over" in small.verdict()


def test_verdict_still_passes_a_clean_fit():
    rng = np.random.default_rng(3)
    a = assess(_fit_with_residuals(0.01 * rng.standard_normal(40)))
    assert a.well_specified and not a.structurally_wrong
    assert "describes this data" in a.verdict()


def test_assess_refuses_to_profile_without_the_data():
    fit = _fit_with_residuals(np.full(20, 0.01))
    with pytest.raises(ValueError, match="profiling needs the data"):
        assess(fit, profile=("phi_m",))


def test_verdict_reports_a_one_sided_profile_as_a_limit():
    """The point of wiring the profile in: a bound is stated as a bound.

    A sweep entirely below the transition does not bound sigma* from above, and
    the Wald +/- cannot say so. The verdict must.
    """
    model = WYART_CATES.with_fixed(eta_s=1.0, phi_0=0.64)
    sigma = np.logspace(-4, -2, 25)  # far below any plausible sigma*
    rng = np.random.default_rng(4)
    eta = model(sigma, 0.40, 0.52, 5.0) * (1.0 + 0.01 * rng.standard_normal(sigma.size))
    fit = fit_model(model, sigma, eta, n_starts=8)
    a = assess(fit, stress=sigma, viscosity=eta, profile=("sigma_star",))

    assert len(a.profiles) == 1
    text = a.verdict()
    prof = a.profiles[0]
    if prof.bounded:
        assert "profile-likelihood interval" in text
    else:
        assert "unbounded" in text or "not bounded" in text
        assert "not a measurement" in text or "limit, not as" in text
