"""Gates for the identification layer (plan Task 4, rules written before the runs).

Gate A (passes): HERSCHEL_BULKLEY on gd = logspace(-2, 2, 25) at 3 % lognormal noise,
20 seeds -> truth inside the 95 % profile interval for >= 17/20 per parameter.
Measured 2026-09-05: tau_y 19, k 19, n 18.

Gate B as pre-registered (FAILS, kept as a strict xfail so the miss stays visible):
gd = logspace(1, 2, 8) with the same truth (12, 4, 0.45) -> tau_y returns
bounded_only or not_identifiable in >= 18/20. Measured: 16/20 non-measurable
(bounded_only 3, degenerate 13, measurable 4). Diagnosis: the design was meant to have
"no yield plateau", but with truth tau_y = 12 Pa the yield stress is 52 % of the total
stress at gd = 10 and 27 % at gd = 100 -- far above 3 % noise -- so tau_y IS partly
determinable there and the four "measurable" verdicts are not wrong (their intervals
contained the truth in 19/20 seeds). The test was the wrong test, not a failed fitter;
the bar is not loosened. Gate B' below is the sharper discriminator: tau_y = 0.3 Pa, so
its share of the stress is 2.6 % -> 0.9 %, at or below the noise.
"""

from __future__ import annotations

import numpy as np
import pytest

from biorheo.constitutive import get_model
from biorheo.identify import (
    VERDICTS,
    ParamEstimate,
    fit_rate_curve,
    recovery_study,
    synthetic_rate_curve,
    verdict_for,
)

HB = get_model("herschel_bulkley")
TRUTH = np.array([12.0, 4.0, 0.45])
GRID_FULL = np.logspace(-2, 2, 25)
GRID_HIGH = np.logspace(1, 2, 8)


class _FakeProfile:
    def __init__(self, lower_at_bound, upper_at_bound):
        self.lower_at_bound, self.upper_at_bound = lower_at_bound, upper_at_bound


# --- verdict rules, checked on synthetic profiles (no fitting) ---------------------

@pytest.mark.parametrize("lo_b, hi_b, corr, expected", [
    (False, False, 0.999, "measurable"),        # closed both sides: correlation irrelevant
    (True, False, 0.5, "bounded_only"),
    (False, True, 0.5, "bounded_only"),
    (True, True, 0.5, "not_identifiable"),
    (True, False, 0.995, "degenerate"),         # open + |corr| > 0.99
    (True, True, -0.995, "degenerate"),         # sign does not matter
    (True, False, 0.99, "bounded_only"),        # exactly 0.99 is not > 0.99
])
def test_verdict_rules(lo_b, hi_b, corr, expected):
    corr_row = np.array([1.0, corr, 0.1])
    verdict, partner = verdict_for(_FakeProfile(lo_b, hi_b), corr_row, 0, ("a", "b", "c"))
    assert verdict == expected
    assert (partner == "b") if expected == "degenerate" else (partner is None)


@pytest.mark.parametrize("lo_b, hi_b, cond, vec, expected, partner", [
    # rank signature: open side + null direction shared by a and b
    (True, False, 1e-12, [0.7071, -0.7071, 0.0], "degenerate", "b"),
    (True, True, 1e-12, [0.7071, -0.7071, 0.0], "degenerate", "b"),
    # closed both sides: the rank signature is not consulted
    (False, False, 1e-12, [0.7071, -0.7071, 0.0], "measurable", None),
    # null direction is this parameter alone (inert, no partner): existing rules apply
    (True, True, 0.0, [1.0, 0.0, 0.0], "not_identifiable", None),
    # null direction does not involve this parameter
    (True, False, 1e-12, [0.0, 0.7071, -0.7071], "bounded_only", None),
    # well-conditioned Jacobian: no degeneracy
    (True, False, 0.1, [0.7071, -0.7071, 0.0], "bounded_only", None),
    # NaN covariance (rank-deficient fit) must not hide the rank signature
    (True, False, 1e-15, [0.7071, -0.7071, 0.0], "degenerate", "b"),
])
def test_verdict_rank_signature(lo_b, hi_b, cond, vec, expected, partner):
    corr_row = np.array([1.0, np.nan, np.nan])  # what a singular covariance looks like
    verdict, got = verdict_for(_FakeProfile(lo_b, hi_b), corr_row, 0, ("a", "b", "c"),
                               cond=cond, null_vec=np.array(vec))
    assert (verdict, got) == (expected, partner)


def test_null_direction_separates_degenerate_from_well_posed():
    """The condition ratio must separate a product-only pair from a well-posed fit by
    orders of magnitude, not by a tuned threshold: measured 1e-11..1e-17 against 0.09..0.11
    on 2026-09-05, i.e. >= 8 decades either side of DEGENERATE_COND = 1e-8."""
    from biorheo.constitutive import STRUCTURAL_STEADY, STRUCTURAL_STEADY_SEPARATE
    from biorheo.identify.fit import DEGENERATE_COND, null_direction

    gd = np.logspace(-2, 3, 30)
    cond, v = null_direction(STRUCTURAL_STEADY_SEPARATE, gd, [3.5e-3, 56.14, 0.05, 4.0])
    assert cond < DEGENERATE_COND * 1e-2
    np.testing.assert_allclose(np.abs(v), [0, 0, 2 ** -0.5, 2 ** -0.5], atol=1e-6)
    assert np.sign(v[2]) == -np.sign(v[3])  # alpha up, t_r down: the product is fixed

    cond, _ = null_direction(STRUCTURAL_STEADY, gd, [3.5e-3, 56.14, 0.2])
    assert cond > DEGENERATE_COND * 1e6
    cond, _ = null_direction(HB, GRID_FULL, TRUTH)
    assert cond > DEGENERATE_COND * 1e6

    # beta = 0 makes alpha_tr inert: the null direction is alpha_tr alone, no partner
    cond, v = null_direction(STRUCTURAL_STEADY, gd, [3.5e-3, 0.0, 0.2])
    assert cond < DEGENERATE_COND
    np.testing.assert_allclose(np.abs(v), [0, 0, 1], atol=1e-9)


def test_single_parameter_model_cannot_be_degenerate():
    verdict, partner = verdict_for(_FakeProfile(True, True), np.array([1.0]), 0, ("mu",))
    assert verdict == "not_identifiable" and partner is None


def test_param_estimate_rejects_unknown_verdict():
    with pytest.raises(ValueError):
        ParamEstimate("x", 1.0, 0.5, 2.0, "plausible", False, False, 0.1)
    assert set(VERDICTS) == {"measurable", "bounded_only", "degenerate", "not_identifiable"}


def test_decades_is_interval_width():
    e = ParamEstimate("x", 1.0, 0.01, 100.0, "measurable", False, False, 0.1)
    assert e.decades == pytest.approx(4.0)
    assert ParamEstimate("x", 1.0, 0.0, 9.0, "measurable", False, False, 0.1).decades == np.inf


# --- one fit ----------------------------------------------------------------------

def test_synthetic_curve_is_seeded_and_lognormal():
    a = synthetic_rate_curve(HB, TRUTH, GRID_FULL, 0.03, 0)
    b = synthetic_rate_curve(HB, TRUTH, GRID_FULL, 0.03, 0)
    c = synthetic_rate_curve(HB, TRUTH, GRID_FULL, 0.03, 1)
    np.testing.assert_array_equal(a, b)
    assert not np.allclose(a, c)
    exact = HB(GRID_FULL, *TRUTH)
    # log-residual is N(0, 0.03): all 25 points inside 5 sigma, none exactly zero
    assert np.all(np.abs(np.log(a / exact)) < 0.15) and np.all(a != exact)


def test_noise_free_curve_is_refused():
    # dst_fit refuses to profile a residual at the level of float64 rounding; the
    # identify layer must surface that rather than return a "measurable" of nothing.
    eta = HB(GRID_FULL, *TRUTH)
    with pytest.raises(ValueError, match="synthetic noise-free"):
        fit_rate_curve(HB, GRID_FULL, eta)


def test_envelope_is_the_data_range_and_dropped_points_are_excluded(capsys):
    gd = np.concatenate([GRID_FULL, [0.0, np.nan]])
    eta = np.concatenate([synthetic_rate_curve(HB, TRUTH, GRID_FULL, 0.03, 0), [1.0, 1.0]])
    ident = fit_rate_curve(HB, gd, eta, temperature_C=20.0, geometry="test")
    assert "dropping 2 of 27" in capsys.readouterr().out
    assert ident.envelope.gammadot_min == pytest.approx(1e-2)
    assert ident.envelope.gammadot_max == pytest.approx(1e2)
    assert ident.fit.n_points == 25
    assert set(ident.identifiability_block()) == {"tau_y", "k", "n"}


def test_full_range_fit_is_measurable_and_tight():
    # Truth-in-interval is a 95 % statement and belongs to the 20-seed gate below (seed 0
    # is in fact one of the 2/20 misses for n). One fit checks structure only: every
    # verdict measurable, intervals narrow, value within 5 Wald errors of the truth.
    eta = synthetic_rate_curve(HB, TRUTH, GRID_FULL, 0.03, 0)
    ident = fit_rate_curve(HB, GRID_FULL, eta)
    for i, p in enumerate(HB.param_names):
        e = ident.params[p]
        assert e.verdict == "measurable"
        assert e.lo < e.value < e.hi
        assert abs(e.value - TRUTH[i]) < 5 * e.wald_stderr
        assert e.decades < 0.2  # tight: 25 points over 4 decades at 3 % noise


# --- the pre-registered recovery gates -------------------------------------------

def test_gate_A_full_range_coverage():
    r = recovery_study(HB, TRUTH, GRID_FULL, n_seeds=20, noise_cv=0.03)
    print(r.summary())
    for p in HB.param_names:
        assert r.within_interval[p] >= 17, r.summary()
        assert r.verdicts[p] == {"measurable": 20}, r.summary()


def _non_measurable(report, p="tau_y"):
    return sum(n for v, n in report.verdicts[p].items() if v != "measurable")


@pytest.mark.xfail(strict=True, reason=(
    "pre-registered Gate B missed on 2026-09-05: 16/20 non-measurable against a bar of 18. "
    "With tau_y = 12 Pa the yield term is 27-52 % of the stress on gd = 10-100 1/s, so the "
    "design was not 'no plateau'; see module docstring and Gate B'."))
def test_gate_B_as_preregistered():
    r = recovery_study(HB, TRUTH, GRID_HIGH, n_seeds=20, noise_cv=0.03)
    print(r.summary())
    assert _non_measurable(r) >= 18, r.summary()


def test_gate_B_prime_no_plateau_design():
    """tau_y share of stress 2.6 % -> 0.9 % on gd = 10-100 1/s, at or below 3 % noise."""
    truth = np.array([0.3, 4.0, 0.45])
    share = truth[0] / (truth[0] + truth[1] * GRID_HIGH ** truth[2])
    assert share.max() < 0.03, share
    r = recovery_study(HB, truth, GRID_HIGH, n_seeds=20, noise_cv=0.03)
    print(r.summary())
    assert _non_measurable(r) >= 18, r.summary()
    # k and n stay measurable even where tau_y is not
    assert r.verdicts["k"] == {"measurable": 20} and r.verdicts["n"] == {"measurable": 20}
    # whatever tau_y verdict came back, the interval kept the truth. Reviewer
    # 2026-09-05: this count is vacuous for tau_y (19/20 intervals run to the bound 0,
    # which contains 0.3 for free), so the coverage that CAN fail is asserted on the
    # closed intervals of k and n, which are 20/20 measurable.
    assert r.within_interval["tau_y"] >= 17
    for p in ("k", "n"):
        assert r.n_closed[p] == 20 and r.within_closed[p] >= 17, r.summary()


def test_casson_recovery():
    r = recovery_study(get_model("casson"), [4.0, 0.01], np.logspace(-1, 3, 25), n_seeds=20)
    for p in ("tau_y", "eta_c"):
        assert r.within_interval[p] >= 17 and r.verdicts[p] == {"measurable": 20}, r.summary()


# --- residual structure: the model-form flag beside the verdicts ----------------------

def test_runs_test_hand_anchor():
    """16 signs --+++++++++----- : 3 runs; n1 = 9, n2 = 7; mu = 1 + 2*63/16 = 8.875;
    var = 2*63*(126-16)/(256*15) = 13860/3840 = 3.609375; z = (3 - 8.875)/1.899836 = -3.0924."""
    from biorheo.identify import runs_test
    r = runs_test(np.array([-1, -1] + [1] * 9 + [-1] * 5, dtype=float))
    assert r.signs == "--+++++++++-----" and r.n_runs == 3
    assert r.expected_runs == pytest.approx(8.875)
    assert r.z == pytest.approx(-3.09237, abs=2e-5)
    assert r.structured
    # alternating signs: the most runs possible, z > 0, not structured
    alt = runs_test(np.array([1, -1] * 8, dtype=float))
    assert alt.n_runs == 16 and alt.z > 0 and not alt.structured
    # one-signed residuals: a bias, which the runs test cannot see; nan and unflagged
    one = runs_test(np.ones(10))
    assert np.isnan(one.z) and not one.structured and one.n_runs == 1
    # exact zeros are dropped, not assigned a sign
    assert runs_test(np.array([1.0, 0.0, -1.0, 0.0, 1.0])).n == 3


def test_synthetic_hb_is_not_flagged_structured():
    """A curve generated by the model being fitted has exchangeable residual signs; the
    flag must stay off on the same designs the verdict gates use (pre-registered z
    threshold -2, one-sided, so about 1 seed in 40 may trip it by chance)."""
    truth = np.array([1.0, 2.0, 0.5])
    tripped = 0
    for seed in range(20):
        eta = synthetic_rate_curve(HB, truth, GRID_FULL, 0.03, seed)
        ident = fit_rate_curve(HB, GRID_FULL, eta, seed=seed)
        assert ident.residual_structure is not None
        tripped += ident.residual_structure.structured
    assert tripped <= 2, tripped


def test_wrong_model_is_flagged_structured():
    """Discriminating case: a Newtonian fit to a shear-thinning curve leaves its
    residuals in two long runs (positive at low rate, negative at high), and the flag
    must come on."""
    from biorheo.constitutive import get_model
    eta = synthetic_rate_curve(HB, [1.0, 2.0, 0.5], GRID_FULL, 0.03, 0)
    ident = fit_rate_curve(get_model("newtonian"), GRID_FULL, eta)
    assert ident.residual_structure.structured, ident.residual_structure
