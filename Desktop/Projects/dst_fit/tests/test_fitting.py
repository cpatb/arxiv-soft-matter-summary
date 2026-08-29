"""Tests for the fitter, including the properties gate N2 checks statistically."""

from __future__ import annotations

import numpy as np
import pytest

from dst_fit.fitting import fit_model, synthetic_curve
from dst_fit.models import WYART_CATES, get_model, thickening_branch

SIGMA = np.logspace(-2, 3, 40)


def test_recovers_herschel_bulkley_from_noiseless_data():
    truth = (2.0, 1.5, 0.6)
    s, e = synthetic_curve("herschel_bulkley", truth, SIGMA)
    r = fit_model("herschel_bulkley", s, e)
    np.testing.assert_allclose(r.params, truth, rtol=1e-4)
    assert r.rmse_log < 1e-6


def test_recovers_wyart_cates_within_its_own_error_bar():
    truth = (0.52, 0.58, 5.0)
    s, e = synthetic_curve(WYART_CATES, truth, SIGMA, noise=0.02, seed=3)
    r = fit_model(WYART_CATES, s, e)
    assert np.all(np.abs(r.params - np.array(truth)) <= 3.0 * r.stderr)


def test_phi_and_phi_m_are_reported_as_correlated():
    """They trade off against each other. A fit that hid that would be lying."""
    s, e = synthetic_curve(WYART_CATES, (0.52, 0.58, 5.0), SIGMA, noise=0.02, seed=1)
    r = fit_model(WYART_CATES, s, e)
    a, b, rho = r.worst_correlation()
    assert {a, b} == {"phi", "phi_m"}
    assert abs(rho) > 0.5


def test_truncation_costs_phi_m_more_than_phi():
    """The structural claim gate N2 tests: phi_m is read off the plateau."""
    s, e = synthetic_curve(WYART_CATES, (0.52, 0.58, 5.0), SIGMA, noise=0.02, seed=0)
    full = fit_model(WYART_CATES, s, e, n_starts=4)
    keep = slice(0, int(0.6 * s.size))
    cut = fit_model(WYART_CATES, s[keep], e[keep], n_starts=4)
    i_phi = WYART_CATES.param_names.index("phi")
    i_phim = WYART_CATES.param_names.index("phi_m")
    assert cut.stderr[i_phim] / full.stderr[i_phim] > cut.stderr[i_phi] / full.stderr[i_phi]


def test_correlation_matrix_has_unit_diagonal():
    s, e = synthetic_curve("ellis", (10.0, 0.05, 3.0, 1.4), SIGMA, noise=0.01, seed=0)
    r = fit_model("ellis", s, e)
    np.testing.assert_allclose(np.diag(r.corr), 1.0, atol=1e-9)


def test_branch_is_recorded_on_the_result():
    s, e = synthetic_curve(WYART_CATES, (0.52, 0.58, 5.0), SIGMA, noise=0.01, seed=0)
    br = slice(5, 40)
    r = fit_model(WYART_CATES, s, e, branch=br)
    assert r.branch == br
    assert r.n_points == 35


def test_refuses_an_underdetermined_fit():
    s, e = synthetic_curve(WYART_CATES, (0.52, 0.58, 5.0), np.logspace(0, 1, 3))
    with pytest.raises(ValueError, match="underdetermined"):
        fit_model(WYART_CATES, s, e)


def test_shape_mismatch_raises():
    with pytest.raises(ValueError, match="shape mismatch"):
        fit_model(WYART_CATES, np.ones(5), np.ones(4))


def test_synthetic_curve_rejects_out_of_bounds_truth():
    with pytest.raises(ValueError, match="outside bounds"):
        synthetic_curve(WYART_CATES, (0.52, 0.58, -3.0), SIGMA)


def test_synthetic_curve_drops_jammed_points_rather_than_clipping():
    """phi above phi_m jams at high stress; those points must vanish, not saturate."""
    s, e = synthetic_curve(WYART_CATES, (0.60, 0.58, 5.0), SIGMA)
    assert s.size < SIGMA.size
    assert np.all(np.isfinite(e))


def test_fit_is_reproducible_across_calls():
    s, e = synthetic_curve(WYART_CATES, (0.52, 0.58, 5.0), SIGMA, noise=0.02, seed=7)
    a = fit_model(WYART_CATES, s, e, seed=11)
    b = fit_model(WYART_CATES, s, e, seed=11)
    np.testing.assert_array_equal(a.params, b.params)


def test_thickening_branch_composes_with_the_fitter():
    """The intended real-data path: cut the branch, then fit only that."""
    gdot = np.logspace(-1, 3, 40)
    eta = 0.1 * gdot**-0.5 + 1e-4 * gdot  # thins then thickens
    sigma = eta * gdot
    br = thickening_branch(eta)
    r = fit_model(get_model("ellis"), sigma, eta, branch=br)
    assert r.n_points == eta.size - br.start
