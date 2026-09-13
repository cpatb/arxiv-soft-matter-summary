"""test_firstpassage.py — the first-passage estimator, on series whose answer is known.

The estimator was adopted AFTER the registered one failed, so it carries more weight
than usual and is anchored on closed forms before any campaign number is read.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from firstpassage import (
    FRAC_ABOVE_MIN, P_PREP_MAX, admitted, bootstrap_exponent, first_passage,
    powerlaw_fit,
)


def _cell(cv, x, dphi, t, frac=1.0, p_prep=1e-11, seed=1):
    """A minimal row of the shape cell_table produces, for the fit tests."""
    return {"cv": cv, "x": x, "dphi": dphi, "seed": seed, "p_prep": p_prep,
            "t_0.0001": t, "frac_0.0001": frac}


# ------------------------------------------------------------------- the estimator
def test_first_passage_on_a_linear_ramp_is_exact():
    g = np.linspace(0, 1.0, 1001)
    p = 2e-3 * g                                   # reaches 1e-4 at gamma = 0.05
    t, frac = first_passage(g, p, 1e-4)
    assert t == pytest.approx(0.05, rel=1e-6)
    assert frac == pytest.approx(1.0)


def test_first_passage_needs_no_p_ss_and_no_monotonicity():
    """The whole point: a series that wanders after the crossing still has a crossing.
    gamma_c (half of p_ss) could not read this; t_fp can."""
    g = np.linspace(0, 1.0, 1001)
    rng = np.random.default_rng(0)
    p = 2e-3 * g + 3e-4 * rng.normal(size=g.size) * (g > 0.2)
    t, frac = first_passage(g, p, 1e-4)
    assert 0.04 < t < 0.07
    assert frac > 0.9


def test_frac_above_catches_a_transient_excursion():
    """A spike through the threshold that falls back must NOT be admitted."""
    g = np.linspace(0, 1.0, 1001)
    p = np.full_like(g, 1e-5)
    p[100:110] = 5e-4                              # a brief poke through 1e-4
    t, frac = first_passage(g, p, 1e-4)
    assert math.isfinite(t)
    assert frac < 0.05
    assert not admitted(_cell(0.0, 0.1, 0.004, t, frac), 1e-4)


def test_never_reached_is_nan_not_an_error():
    g = np.linspace(0, 3.0, 301)
    t, frac = first_passage(g, np.full_like(g, 1e-9), 1e-4)
    assert math.isnan(t) and math.isnan(frac)


def test_a_series_already_above_at_gamma_zero_returns_the_first_sample():
    g = np.linspace(0, 1.0, 101)
    p = np.full_like(g, 1e-3)
    t, frac = first_passage(g, p, 1e-4)
    assert t == pytest.approx(0.0)
    assert frac == pytest.approx(1.0)


def test_a_higher_threshold_is_reached_later():
    g = np.linspace(0, 1.0, 2001)
    p = 2e-3 * g
    assert first_passage(g, p, 3e-4)[0] > first_passage(g, p, 1e-4)[0]


# ------------------------------------------------------------- the admission rule
def test_admission_refuses_a_jammed_preparation():
    r = _cell(0.5, 0.85, 0.039, 0.01, frac=1.0, p_prep=2e-4)
    assert r["p_prep"] > P_PREP_MAX
    assert not admitted(r, 1e-4)


def test_admission_boundary_is_inclusive_at_the_registered_value():
    assert admitted(_cell(0.0, 0.1, 0.004, 0.05, frac=FRAC_ABOVE_MIN), 1e-4)
    assert not admitted(_cell(0.0, 0.1, 0.004, 0.05, frac=FRAC_ABOVE_MIN - 1e-9), 1e-4)


# --------------------------------------------------------------------- the fitting
def test_powerlaw_fit_recovers_a_known_exponent():
    """An anchor outside the code under test: synthetic cells on an exact power law."""
    a_true, A = 1.25, 3e-5
    rows = []
    for x, dphi in ((0.35, 0.0145), (0.60, 0.0248), (0.85, 0.0351)):
        for seed in (1, 2, 3):
            rows.append(_cell(0.0, x, dphi, A * dphi ** (-a_true), seed=seed))
    a, lnA, n = powerlaw_fit(rows, 0.0, 1e-4)
    assert a == pytest.approx(a_true, rel=1e-9)
    assert math.exp(lnA) == pytest.approx(A, rel=1e-9)
    assert n == 9


def test_powerlaw_fit_ignores_the_rungs_it_is_told_not_to_see():
    """F2 is only a test if the fit never touches x = 0.05. A rung at a wildly wrong
    value must not move the answer."""
    a_true, A = 1.25, 3e-5
    rows = []
    for x, dphi in ((0.35, 0.0145), (0.60, 0.0248), (0.85, 0.0351)):
        for seed in (1, 2, 3):
            rows.append(_cell(0.0, x, dphi, A * dphi ** (-a_true), seed=seed))
    clean = powerlaw_fit(rows, 0.0, 1e-4)[0]
    rows += [_cell(0.0, 0.05, 0.0021, 999.0, seed=s) for s in (1, 2, 3)]
    assert powerlaw_fit(rows, 0.0, 1e-4)[0] == pytest.approx(clean, rel=1e-12)


def test_powerlaw_fit_excludes_inadmissible_cells():
    a_true, A = 1.25, 3e-5
    rows = []
    for x, dphi in ((0.35, 0.0145), (0.60, 0.0248), (0.85, 0.0351)):
        for seed in (1, 2, 3):
            rows.append(_cell(0.0, x, dphi, A * dphi ** (-a_true), seed=seed))
    rows.append(_cell(0.0, 0.85, 0.0351, 50.0, frac=0.2, seed=4))    # transient
    rows.append(_cell(0.0, 0.85, 0.0351, 50.0, p_prep=1e-3, seed=5))  # prep jammed
    assert powerlaw_fit(rows, 0.0, 1e-4)[0] == pytest.approx(a_true, rel=1e-9)


def test_bootstrap_exponent_is_zero_on_noiseless_data_and_grows_with_scatter():
    """A se that does not respond to the scatter is not a se."""
    a_true, A = 1.25, 3e-5
    rungs = ((0.35, 0.0145), (0.60, 0.0248), (0.85, 0.0351))
    exact = [_cell(0.0, x, d, A * d ** (-a_true), seed=s)
             for x, d in rungs for s in (1, 2, 3)]
    assert bootstrap_exponent(exact, 0.0, 1e-4, n_boot=200) == pytest.approx(0.0, abs=1e-9)
    rng = np.random.default_rng(1)
    noisy = [_cell(0.0, x, d, A * d ** (-a_true) * math.exp(0.3 * rng.normal()), seed=s)
             for x, d in rungs for s in (1, 2, 3)]
    assert bootstrap_exponent(noisy, 0.0, 1e-4, n_boot=500) > 0.05


def test_bootstrap_is_reproducible_from_its_seed():
    rng = np.random.default_rng(2)
    rows = [_cell(0.0, x, d, 3e-5 * d ** (-1.25) * math.exp(0.3 * rng.normal()), seed=s)
            for x, d in ((0.35, 0.0145), (0.60, 0.0248), (0.85, 0.0351))
            for s in (1, 2, 3)]
    a = bootstrap_exponent(rows, 0.0, 1e-4, n_boot=300, seed=7)
    b = bootstrap_exponent(rows, 0.0, 1e-4, n_boot=300, seed=7)
    assert a == b


# ------------------------------------------------------------ F6: the ratio a/b
from firstpassage import bootstrap_ratio, pss_exponent      # noqa: E402


def _cell_p(cv, x, dphi, t, pss, frac=1.0, p_prep=1e-11, seed=1):
    r = _cell(cv, x, dphi, t, frac, p_prep, seed)
    r["p_ss"] = pss
    return r


def test_pss_exponent_recovers_a_known_exponent():
    b_true, C = 1.10, 0.21
    rows = [_cell_p(0.0, x, d, 1.0, C * d ** b_true, seed=s)
            for x, d in ((0.35, 0.0145), (0.60, 0.0248), (0.85, 0.0351))
            for s in (1, 2, 3)]
    b, n = pss_exponent(rows, 0.0, 1e-4)
    assert b == pytest.approx(b_true, rel=1e-9)
    assert n == 9


def test_ratio_is_exactly_one_on_a_universal_build_curve():
    """The F6 null, constructed rather than asserted: if t_fp = K / p_ss exactly (the
    small-threshold limit of p(gamma) = p_ss f(gamma/gamma*) with gamma* fixed), then
    a/b = 1 whatever b is. A routine that returned a/b = 1 by construction would pass
    this AND the next test, which is why both are here."""
    for b_true in (0.8, 1.0, 1.3):
        rows = [_cell_p(0.0, x, d, 1e-4 / (0.21 * d ** b_true), 0.21 * d ** b_true,
                        seed=s)
                for x, d in ((0.35, 0.0145), (0.60, 0.0248), (0.85, 0.0351))
                for s in (1, 2, 3)]
        a = powerlaw_fit(rows, 0.0, 1e-4)[0]
        b = pss_exponent(rows, 0.0, 1e-4)[0]
        assert b == pytest.approx(b_true, rel=1e-9)
        assert a / b == pytest.approx(1.0, rel=1e-9), b_true


def test_ratio_departs_from_one_when_the_clock_is_not_just_the_pressure_scale():
    """The discriminating case: a build whose gamma* itself depends on dphi must give
    a/b != 1, or the statistic cannot see a kinetic effect at all."""
    b_true, k = 1.0, 0.4
    rows = [_cell_p(0.0, x, d, 1e-4 / (0.21 * d ** b_true) * d ** (-k), 0.21 * d ** b_true,
                    seed=s)
            for x, d in ((0.35, 0.0145), (0.60, 0.0248), (0.85, 0.0351))
            for s in (1, 2, 3)]
    a = powerlaw_fit(rows, 0.0, 1e-4)[0]
    b = pss_exponent(rows, 0.0, 1e-4)[0]
    assert a / b == pytest.approx(1.0 + k, rel=1e-6)


def test_ratio_bootstrap_responds_to_scatter():
    rungs = ((0.35, 0.0145), (0.60, 0.0248), (0.85, 0.0351))
    exact = [_cell_p(0.0, x, d, 1e-4 / (0.21 * d), 0.21 * d, seed=s)
             for x, d in rungs for s in (1, 2, 3)]
    assert bootstrap_ratio(exact, 0.0, 1e-4, n_boot=200) == pytest.approx(0.0, abs=1e-9)
    rng = np.random.default_rng(3)
    noisy = [_cell_p(0.0, x, d, 1e-4 / (0.21 * d) * math.exp(0.3 * rng.normal()),
                     0.21 * d, seed=s) for x, d in rungs for s in (1, 2, 3)]
    assert bootstrap_ratio(noisy, 0.0, 1e-4, n_boot=500) > 0.05
