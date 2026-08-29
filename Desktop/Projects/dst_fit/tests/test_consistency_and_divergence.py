"""The cross-curve consistency statistic, the phi-space divergence fit, and
the known-phi variant of Wyart-Cates.

These three are what turned N4's single failure into a comparison across three
datasets, so each is tested against something other than a previous run of
itself: closed-form statistics for the consistency check, a generated curve with
a known pole for the divergence fit, and an identifiability property for the
known-phi model.
"""

from __future__ import annotations

import numpy as np
import pytest

from dst_fit import consistency_across_curves, fit_divergence, fit_model, weighted_mean
from dst_fit.fitting import synthetic_curve
from dst_fit.models import WYART_CATES, WYART_CATES_KNOWN_PHI


# --------------------------------------------------------------------------
# weighted mean and consistency
# --------------------------------------------------------------------------
def test_weighted_mean_matches_the_closed_form():
    v = np.array([1.0, 2.0, 3.0])
    e = np.array([1.0, 1.0, 2.0])
    w = 1.0 / e**2
    mean, err = weighted_mean(v, e)
    assert mean == pytest.approx(float(np.sum(w * v) / np.sum(w)))
    assert err == pytest.approx(float(np.sqrt(1.0 / np.sum(w))))


def test_weighted_mean_ignores_useless_errors():
    v = np.array([1.0, 2.0, 99.0])
    e = np.array([1.0, 1.0, 0.0])  # a zero error carries infinite weight; drop it
    mean, _ = weighted_mean(v, e)
    assert mean == pytest.approx(1.5)


def test_consistent_estimates_all_pass():
    """Same value, honest errors: everything agrees and chi2 is zero."""
    v = np.full(6, 0.558)
    e = np.full(6, 0.005)
    c = consistency_across_curves(v, e)
    assert c.n_consistent == 6
    assert c.n_pairs_consistent == c.n_pairs == 15
    assert c.chi2 == pytest.approx(0.0)
    assert c.mean == pytest.approx(0.558)


def test_a_spread_larger_than_the_errors_is_caught():
    v = np.array([0.49, 0.53, 0.54, 0.55, 0.56, 0.58])
    e = np.full(6, 0.002)
    c = consistency_across_curves(v, e)
    assert c.n_consistent < 4
    assert c.reduced_chi2 > 10


def test_chi2_is_the_sum_of_squared_deviations():
    v = np.array([0.50, 0.52, 0.54])
    e = np.array([0.01, 0.02, 0.01])
    c = consistency_across_curves(v, e)
    assert c.chi2 == pytest.approx(float(np.sum(((v - c.mean) / e) ** 2)))
    assert c.dof == 2
    np.testing.assert_allclose(c.deviations, np.abs(v - c.mean) / e)


def test_pairwise_is_the_stricter_reading():
    """It must never call more pairs consistent than a mean-based count implies."""
    rng = np.random.default_rng(0)
    for _ in range(20):
        v = 0.55 + 0.02 * rng.standard_normal(6)
        e = np.full(6, 0.01)
        c = consistency_across_curves(v, e)
        assert 0 <= c.n_pairs_consistent <= c.n_pairs == 15


def test_unweighted_centre_is_reported_alongside():
    """The weighted mean follows the smallest error bar; the median does not."""
    v = np.array([0.50, 0.51, 0.52, 0.90])
    e = np.array([0.05, 0.05, 0.05, 0.0001])
    c = consistency_across_curves(v, e)
    assert c.mean == pytest.approx(0.90, abs=1e-3), "weighted mean is captured"
    assert c.median == pytest.approx(0.515)
    assert "median" in c.report()


def test_consistency_refuses_a_single_estimate():
    with pytest.raises(ValueError, match="at least two"):
        consistency_across_curves(np.array([0.5]), np.array([0.01]))


def test_consistency_refuses_mismatched_shapes():
    with pytest.raises(ValueError, match="shape mismatch"):
        consistency_across_curves(np.zeros(3), np.zeros(4))


# --------------------------------------------------------------------------
# fit_divergence
# --------------------------------------------------------------------------
def test_divergence_recovers_a_known_pole():
    """Generate from A (1 - phi/phi_c)^-n and get all three back."""
    phi = np.linspace(0.30, 0.54, 12)
    A, phi_c, n = 0.7, 0.558, 2.2
    eta = A * (1.0 - phi / phi_c) ** -n
    fit = fit_divergence(phi, eta, free_exponent=True, n_starts=24)
    assert fit.phi_c == pytest.approx(phi_c, abs=1e-4)
    assert fit.n == pytest.approx(n, abs=1e-3)
    assert fit.A == pytest.approx(A, rel=1e-3)
    assert fit.rmse_log < 1e-6


def test_divergence_recovers_the_pole_under_noise():
    rng = np.random.default_rng(7)
    phi = np.linspace(0.30, 0.54, 12)
    eta = 0.7 * (1.0 - phi / 0.558) ** -2.2 * np.exp(rng.normal(0, 0.05, phi.size))
    fit = fit_divergence(phi, eta, free_exponent=True, n_starts=24)
    assert abs(fit.phi_c - 0.558) < 3 * fit.stderr["phi_c"]


def test_pinning_the_exponent_moves_the_pole():
    """Why the exponent has to be settled before phi_m values are compared."""
    phi = np.linspace(0.30, 0.54, 12)
    eta = 0.7 * (1.0 - phi / 0.558) ** -2.2
    at_2 = fit_divergence(phi, eta, free_exponent=False, n=2.0, n_starts=24)
    at_22 = fit_divergence(phi, eta, free_exponent=False, n=2.2, n_starts=24)
    assert at_22.phi_c == pytest.approx(0.558, abs=1e-4)
    # The bias is small in absolute terms (~0.004 here, ~0.005 on the real Guy
    # branches) but large against the error bar the fit reports, which is the
    # form in which it corrupts a comparison against a published phi_m.
    bias = abs(at_2.phi_c - 0.558)
    assert bias > 5 * at_2.stderr["phi_c"], (
        f"exponent bias {bias:.4f} should dominate the quoted error "
        f"{at_2.stderr['phi_c']:.4f}"
    )
    assert bias > 0.5 * 0.005, "and should be comparable with a published +/-0.005"
    assert at_2.free_exponent is False and at_2.n == 2.0
    assert at_2.stderr["n"] == 0.0


def test_divergence_refuses_a_pole_inside_the_data():
    """phi_c below the largest measured phi would mean the branch contains its own pole."""
    phi = np.linspace(0.30, 0.96, 12)
    eta = np.linspace(2.0, 50.0, 12)
    with pytest.raises(ValueError, match="no room for phi_c"):
        fit_divergence(phi, eta)


def test_divergence_refuses_an_underdetermined_fit():
    with pytest.raises(ValueError, match="underdetermined"):
        fit_divergence(np.array([0.3, 0.4, 0.5]), np.array([2.0, 3.0, 5.0]))


def test_divergence_summary_names_the_exponent_convention():
    phi = np.linspace(0.30, 0.54, 12)
    eta = 0.7 * (1.0 - phi / 0.558) ** -2.2
    assert "fixed" in fit_divergence(phi, eta, free_exponent=False, n=2.0).summary()
    assert "+/-" in fit_divergence(phi, eta, free_exponent=True).summary()


# --------------------------------------------------------------------------
# Wyart-Cates with phi supplied
# --------------------------------------------------------------------------
def test_known_phi_model_matches_the_free_one_at_the_same_parameters():
    """Same physics; the only difference is which numbers are fitted."""
    sigma = np.logspace(-2, 4, 80)
    free = WYART_CATES(sigma, 0.50, 0.58, 5.0)
    known = WYART_CATES_KNOWN_PHI.with_fixed(phi=0.50)(sigma, 0.58, 5.0)
    np.testing.assert_allclose(known, free, rtol=1e-12)


def test_known_phi_must_be_supplied():
    sigma = np.logspace(-2, 4, 10)
    with pytest.raises(TypeError):
        WYART_CATES_KNOWN_PHI(sigma, 0.58, 5.0)


def test_pinning_phi_removes_the_phi_phi_m_trade_off():
    """The identifiability claim, stated as a test.

    Generate a curve at a phi close to jamming, where the free fit is known to
    trade phi against phi_m. With phi supplied, phi_m must come back at the
    truth and with a smaller relative error than the free fit achieves.
    """
    truth_phi, truth_phi_m, truth_ss = 0.55, 0.58, 5.0
    sigma = np.logspace(-1, 3, 60)
    s, eta = synthetic_curve(WYART_CATES, (truth_phi, truth_phi_m, truth_ss),
                             sigma, noise=0.02, seed=11)

    free = fit_model(WYART_CATES, s, eta, n_starts=16)
    i = free.model.param_names.index("phi_m")
    j = free.model.param_names.index("phi")
    assert abs(free.corr[i, j]) > 0.9, "the free fit should show the trade-off"

    known = fit_model(WYART_CATES_KNOWN_PHI.with_fixed(phi=truth_phi), s, eta, n_starts=16)
    assert known.params[0] == pytest.approx(truth_phi_m, abs=0.005)
    assert known.stderr[0] < free.stderr[i]


def test_known_phi_bounds_track_phi_0():
    m = WYART_CATES_KNOWN_PHI.with_fixed(phi=0.50, phi_0=0.546)
    assert m.bounds[1][0] == pytest.approx(0.545)
    assert m.bounds[1][1] == WYART_CATES_KNOWN_PHI.bounds[1][1]
