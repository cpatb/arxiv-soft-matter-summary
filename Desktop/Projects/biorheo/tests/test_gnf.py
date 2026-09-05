"""Unit 0.3a: rate-space GNF models against hand values and analytic limits.

Every numeric anchor below was computed by hand from the formula (arithmetic shown)
rather than by running the module. The limits are those listed in each model's
``limits`` mapping and derived in ``derivations/gnf_limits.md``.
"""

from __future__ import annotations

import numpy as np
import pytest
from dst_fit import fit_model

from biorheo.constitutive import (
    CARREAU_YASUDA,
    CASSON,
    CROSS,
    HERSCHEL_BULKLEY,
    MODELS,
    NEWTONIAN,
    POWER_LAW,
    RateModel,
    get_model,
)


# --- registry shape -----------------------------------------------------------------

def test_registry_names_and_units_cover_every_parameter():
    assert set(MODELS) == {
        "newtonian", "power_law", "herschel_bulkley", "casson", "carreau_yasuda", "cross",
        "structural_steady",
    }
    for m in MODELS.values():
        assert isinstance(m, RateModel)
        assert m.control == "rate"
        assert set(m.units) == set(m.param_names), m.name
        assert m.limits, f"{m.name} documents no limiting cases"
        assert not m.is_regularised
        lo, hi = m.bounds
        assert len(lo) == len(hi) == m.n_params


def test_get_model_unknown_name_lists_options():
    with pytest.raises(KeyError, match="herschel_bulkley"):
        get_model("bingham")


# --- hand values --------------------------------------------------------------------

def test_power_law_hand_value():
    # k = 2, n = 0.5, gd = 4: eta = 2 * 4^(-0.5) = 2 * 0.5 = 1.0 exactly.
    assert POWER_LAW(4.0, 2.0, 0.5) == pytest.approx(1.0)


def test_herschel_bulkley_hand_value():
    # tau_y = 23.756, k = 8.903, n = 0.4526, gd = 10:
    # 23.756/10 = 2.3756; 10^(0.4526 - 1) = 10^(-0.5474) = 0.28356; 8.903 * 0.28356 = 2.5246
    # total 4.9002 Pa s.
    assert HERSCHEL_BULKLEY(10.0, 23.756, 8.903, 0.4526) == pytest.approx(4.9002, rel=2e-4)


def test_casson_hand_value():
    # tau_y = 4, eta_c = 0.01, gd = 100: sqrt(4/100) = 0.2, sqrt(0.01) = 0.1, (0.3)^2 = 0.09.
    assert CASSON(100.0, 4.0, 0.01) == pytest.approx(0.09)


def test_carreau_yasuda_hand_value():
    # eta_0 = 1, eta_inf = 0, lam = 1, a = 2, n = 0.5, gd = sqrt(3):
    # (1 + 3)^((0.5 - 1)/2) = 4^(-0.25) = 1/sqrt(2) = 0.70710678.
    assert CARREAU_YASUDA(np.sqrt(3.0), 1.0, 0.0, 1.0, 2.0, 0.5) == pytest.approx(
        0.70710678, rel=1e-7)


def test_cross_hand_value():
    # eta_0 = 1, eta_inf = 0.1, lam = 1, m = 1, gd = 9: 0.1 + 0.9/(1 + 9) = 0.19.
    assert CROSS(9.0, 1.0, 0.1, 1.0, 1.0) == pytest.approx(0.19)


def test_newtonian_is_flat():
    gd = np.logspace(-3, 3, 7)
    np.testing.assert_allclose(NEWTONIAN(gd, 0.0035), 0.0035)


# --- analytic limits ----------------------------------------------------------------

GD = np.logspace(-4, 4, 33)


def test_hb_tau_y_zero_is_power_law():
    np.testing.assert_allclose(HERSCHEL_BULKLEY(GD, 0.0, 3.0, 0.7), POWER_LAW(GD, 3.0, 0.7))


def test_hb_n_one_is_bingham():
    # eta = tau_y/gd + k
    np.testing.assert_allclose(HERSCHEL_BULKLEY(GD, 5.0, 2.0, 1.0), 5.0 / GD + 2.0)


def test_hb_stress_tends_to_tau_y_at_low_rate():
    sigma = HERSCHEL_BULKLEY.stress(np.array([1e-6, 1e-8]), 5.0, 2.0, 0.5)
    # sigma = 5 + 2 gd^0.5: at 1e-8 that is 5 + 2e-4, so relative excess 4e-5.
    assert sigma[1] == pytest.approx(5.0, rel=1e-4)
    assert sigma[0] > sigma[1] > 5.0


def test_hb_diverges_at_zero_rate():
    assert np.isinf(HERSCHEL_BULKLEY(0.0, 5.0, 2.0, 0.5))
    assert np.isnan(HERSCHEL_BULKLEY.stress(0.0, 5.0, 2.0, 0.5))


def test_power_law_n_one_is_newtonian():
    np.testing.assert_allclose(POWER_LAW(GD, 0.7, 1.0), NEWTONIAN(GD, 0.7))


def test_casson_tau_y_zero_is_newtonian():
    np.testing.assert_allclose(CASSON(GD, 0.0, 0.004), 0.004)


def test_casson_high_rate_expansion():
    # eta = eta_c + 2 sqrt(tau_y eta_c/gd) + tau_y/gd, exactly (it is a binomial square).
    tau_y, eta_c = 0.005, 0.004
    np.testing.assert_allclose(
        CASSON(GD, tau_y, eta_c),
        eta_c + 2 * np.sqrt(tau_y * eta_c / GD) + tau_y / GD, rtol=1e-12)


def test_carreau_yasuda_plateaux():
    p = (0.056, 0.00345, 3.313, 2.0, 0.3568)  # arbitrary blood-like magnitudes
    assert CARREAU_YASUDA(1e-9, *p) == pytest.approx(0.056, rel=1e-6)
    # at gd = 1e9 the residual is (eta_0 - eta_inf)(lam gd)^(n-1) = 0.05255 * (3.313e9)^(-0.6432)
    # = 0.05255 * 7.5e-7 = 4e-8, i.e. 1.2e-5 of eta_inf, so the plateau holds to 1e-4
    assert CARREAU_YASUDA(1e9, *p) == pytest.approx(0.00345, rel=1e-4)


def test_carreau_yasuda_high_rate_is_power_law_when_eta_inf_zero():
    # (lam gd)^a >> 1: eta ~ eta_0 (lam gd)^(n-1) = eta_0 lam^(n-1) gd^(n-1)
    eta_0, lam, a, n = 2.0, 0.5, 2.0, 0.4
    gd = np.array([1e3, 1e4])
    k = eta_0 * lam ** (n - 1)
    np.testing.assert_allclose(CARREAU_YASUDA(gd, eta_0, 0.0, lam, a, n),
                               POWER_LAW(gd, k, n), rtol=1e-5)


def test_carreau_yasuda_n_one_is_flat():
    np.testing.assert_allclose(CARREAU_YASUDA(GD, 2.0, 0.1, 3.0, 2.0, 1.0), 2.0)


def test_cross_half_way_rate():
    eta_0, eta_inf, lam, m = 3.0, 1.0, 0.25, 0.8
    assert CROSS(1.0 / lam, eta_0, eta_inf, lam, m) == pytest.approx(2.0)


def test_cross_high_rate_is_power_law_when_eta_inf_zero():
    eta_0, lam, m = 2.0, 0.5, 0.7
    # 1/(1 + x) differs from 1/x by 1/x; at gd = 1e6, x = (5e5)^0.7 = 9.8e3, so 1e-4
    gd = np.array([1e6, 1e7])
    np.testing.assert_allclose(CROSS(gd, eta_0, 0.0, lam, m),
                               POWER_LAW(gd, eta_0 * lam ** (-m), 1.0 - m), rtol=2e-4)


# --- discriminating negative ----------------------------------------------------------

def test_models_disagree_where_they_should():
    # A Bingham fluid and a Casson fluid with the same tau_y and the same high-rate
    # viscosity are different functions in the transition; a test suite in which
    # every model agrees with every other has tested nothing.
    gd = np.array([0.1, 1.0, 10.0])
    bing = HERSCHEL_BULKLEY(gd, 1.0, 0.01, 1.0)
    cass = CASSON(gd, 1.0, 0.01)
    assert np.all(cass > bing)  # the cross term 2 sqrt(tau_y eta_c/gd) is positive
    assert np.max(cass / bing) > 1.05


# --- the dst_fit fitter runs in rate space -------------------------------------------

def test_dst_fit_recovers_hb_in_rate_space():
    # Not an anchor for the physics (a round trip certifies the optimiser only); it
    # certifies that a RateModel plugs into dst_fit unchanged and that ``guess`` lands.
    true = np.array([12.0, 4.0, 0.45])
    gd = np.logspace(-2, 3, 40)
    eta = HERSCHEL_BULKLEY(gd, *true)
    fit = fit_model(HERSCHEL_BULKLEY, gd, eta, start=HERSCHEL_BULKLEY.guess(gd, eta),
                    n_starts=2, seed=0)
    assert fit.success
    np.testing.assert_allclose(fit.params, true, rtol=1e-6)
    assert fit.n_points == 40


def test_guess_lands_inside_bounds_for_every_model():
    gd = np.logspace(-2, 3, 25)
    eta = 0.05 + 0.5 / (1 + gd)  # a generic thinning curve
    for m in MODELS.values():
        g = m.guess(gd, eta)
        lo, hi = np.array(m.bounds[0]), np.array(m.bounds[1])
        assert g.shape == (m.n_params,), m.name
        assert np.all(g > lo) and np.all(g < hi), (m.name, g)
