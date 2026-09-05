"""Unit 0.5b: the structural (thixotropic) family.

Anchors outside the module: the closed forms in ``derivations/structural_kinetics.md``
(SymPy), hand arithmetic shown beside each number, and the Cross model already tested in
``test_gnf.py``. The identify gate is pre-registered in the plan (Task 5): fitting the
steady curve with ``alpha`` and ``t_r`` separate must return ``degenerate`` for the pair,
and ``measurable`` for the product.

Gate as written before the run (2026-09-05): on gd = logspace(-2, 3, 30), truth
(eta_inf, beta, alpha_tr) = (3.5e-3, 56.14, 0.2) at 3 % lognormal noise, 20 seeds,
STRUCTURAL_STEADY returns measurable for all three in >= 18/20 with coverage >= 17/20;
STRUCTURAL_STEADY_SEPARATE returns degenerate (partner = the other of the pair) for
alpha and t_r in >= 18/20 while eta_inf and beta stay measurable in >= 18/20.

First run (corr-only degeneracy rule): 16/20 degenerate, 4/20 bounded_only, 20/20 in
interval, eta_inf and beta 20/20 measurable. The 4 misses (seeds 0, 13, 14, 16) were
fits that never left the symmetric start alpha = t_r = sqrt(alpha_tr); there the
Jacobian is exactly rank-deficient, the covariance is all-NaN, and the corr rule read
NaN as "no evidence". The bar stayed at 18; the detector gained an independent
Fisher-rank channel (``null_direction``, condition ratio 1e-11..1e-17 on these fits
against 0.09..0.11 for well-posed ones). Second run: 20/20.
"""

from __future__ import annotations

import numpy as np
import pytest

from biorheo.constitutive import (
    CROSS,
    MODELS,
    STRUCTURAL_STEADY,
    STRUCTURAL_STEADY_SEPARATE,
    get_model,
    lambda_step,
    lambda_steady,
    lambda_transient,
    ramp_loop,
    relaxation_time,
)
from biorheo.identify import recovery_study

# the source's steady constants: eta_inf = 3.5 mPa s, eta_0 = 200 mPa s, gamma_c = 5 1/s
ETA_INF, BETA, ALPHA_TR = 3.5e-3, 0.2 / 3.5e-3 - 1.0, 0.2
TRUTH = np.array([ETA_INF, BETA, ALPHA_TR])
GRID = np.logspace(-2, 3, 30)


# --- steady state, hand anchors -----------------------------------------------------

def test_registry_has_steady_but_not_the_diagnostic_model():
    assert get_model("structural_steady") is STRUCTURAL_STEADY
    assert "structural_steady_separate" not in MODELS
    assert set(STRUCTURAL_STEADY.units) == set(STRUCTURAL_STEADY.param_names)


def test_steady_hand_values():
    # beta = 0.2/3.5e-3 - 1 = 57.142857 - 1 = 56.142857
    assert BETA == pytest.approx(56.142857142857, rel=1e-9)
    # gd = 5 = 1/alpha_tr: lam_ss = 1/2, eta = 3.5e-3 (1 + 28.0714286) = 3.5e-3 * 29.0714286 = 0.10175
    assert lambda_steady(5.0, ALPHA_TR) == pytest.approx(0.5)
    assert STRUCTURAL_STEADY(5.0, *TRUTH) == pytest.approx(0.10175, rel=1e-9)
    # gd = 45: alpha_tr gd = 9, lam_ss = 0.1, eta = 3.5e-3 (1 + 5.6142857) = 0.02315
    assert STRUCTURAL_STEADY(45.0, *TRUTH) == pytest.approx(0.02315, rel=1e-9)
    # plateaux
    assert STRUCTURAL_STEADY(1e-9, *TRUTH) == pytest.approx(0.2, rel=1e-6)
    assert STRUCTURAL_STEADY(1e9, *TRUTH) == pytest.approx(ETA_INF, rel=1e-6)


def test_steady_is_cross_with_m_equal_one():
    gd = np.logspace(-3, 4, 50)
    cross = CROSS(gd, ETA_INF * (1.0 + BETA), ETA_INF, ALPHA_TR, 1.0)
    np.testing.assert_allclose(STRUCTURAL_STEADY(gd, *TRUTH), cross, rtol=1e-12)
    # and the separate-parameter form is the same curve for any split of the product
    for alpha in (0.01, 0.2, 7.0):
        np.testing.assert_allclose(STRUCTURAL_STEADY_SEPARATE(gd, ETA_INF, BETA, alpha, ALPHA_TR / alpha),
                                   STRUCTURAL_STEADY(gd, *TRUTH), rtol=1e-12)


def test_beta_zero_is_newtonian():
    gd = np.logspace(-2, 3, 20)
    np.testing.assert_allclose(STRUCTURAL_STEADY(gd, ETA_INF, 0.0, ALPHA_TR), ETA_INF)


def test_steady_is_monotone_thinning_and_lambda_in_unit_interval():
    gd = np.logspace(-4, 5, 200)
    eta = STRUCTURAL_STEADY(gd, *TRUTH)
    assert np.all(np.diff(eta) < 0)
    lam = lambda_steady(gd, ALPHA_TR)
    assert np.all((lam > 0) & (lam <= 1)) and lam[0] > 0.999 and lam[-1] < 1e-3


# --- transients ---------------------------------------------------------------------

def test_relaxation_time_hand_values():
    # 1/tau = 1/t_r + alpha gd: t_r = 1, alpha = 0.05, gd = 20 -> 1/tau = 2 -> tau = 0.5
    assert relaxation_time(20.0, 0.05, 1.0) == pytest.approx(0.5)
    assert relaxation_time(0.0, 0.05, 1.0) == pytest.approx(1.0)  # rest: t_r


def test_step_response_hand_value():
    # from lam0 = 1 to gd = 20 (alpha = 0.05, t_r = 1): lam_ss = 0.5, tau = 0.5,
    # lam(0.5) = 0.5 + 0.5 exp(-1) = 0.5 + 0.1839397 = 0.6839397
    assert lambda_step(0.5, 20.0, 0.05, 1.0, 1.0) == pytest.approx(0.68393972058572, rel=1e-9)


def test_ode_integrator_matches_the_closed_form_step():
    t = np.linspace(0.0, 5.0, 201)
    for gd, alpha, t_r, lam0 in [(20.0, 0.05, 1.0, 1.0), (0.5, 2.0, 0.1, 0.0), (300.0, 0.01, 3.0, 0.4)]:
        lam = lambda_transient(t, lambda _: gd, alpha, t_r, lam0)
        np.testing.assert_allclose(lam, lambda_step(t, gd, alpha, t_r, lam0), rtol=1e-6, atol=1e-8)


def test_transient_rejects_bad_inputs():
    with pytest.raises(ValueError, match="strictly increasing"):
        lambda_transient(np.array([0.0, 1.0, 0.5]), lambda _: 1.0, 0.05, 1.0, 1.0)
    with pytest.raises(ValueError, match="outside"):
        lambda_transient(np.linspace(0, 1, 5), lambda _: 1.0, 0.05, 1.0, 1.2)


def test_ramp_loop_area_is_positive_and_vanishes_quasi_statically():
    """Thixotropic hysteresis: the up branch (structure not yet broken) carries the higher
    stress, so the signed area is positive; a ramp slow against t_r collapses both
    branches onto the steady curve and the area -> 0. Checked as a monotone ladder."""
    alpha, t_r = 0.05, 1.0
    areas = []
    for t_ramp in (0.5, 5.0, 50.0, 500.0):
        loop = ramp_loop(ETA_INF, BETA, alpha, t_r, gd_max=100.0, t_ramp=t_ramp)
        assert loop.gd[0] == 0.0 and loop.gd[-1] == pytest.approx(0.0, abs=1e-12)
        assert np.all((loop.lam >= 0) & (loop.lam <= 1))
        areas.append(loop.area)
    assert all(a > 0 for a in areas), areas
    assert all(a1 > a2 for a1, a2 in zip(areas, areas[1:])), areas
    # quasi-static: the down branch sits on the steady curve to < 1 %
    loop = ramp_loop(ETA_INF, BETA, alpha, t_r, gd_max=100.0, t_ramp=500.0)
    gd_dn, sig_dn = loop.down
    mask = gd_dn > 1.0
    steady = STRUCTURAL_STEADY(gd_dn[mask], ETA_INF, BETA, alpha * t_r) * gd_dn[mask]
    np.testing.assert_allclose(sig_dn[mask], steady, rtol=1e-2)
    # and the fast ramp's up branch sits ABOVE it (structure still intact)
    fast = ramp_loop(ETA_INF, BETA, alpha, t_r, gd_max=100.0, t_ramp=0.5)
    gd_up, sig_up = fast.up
    m = gd_up > 1.0
    assert np.all(sig_up[m] > STRUCTURAL_STEADY(gd_up[m], ETA_INF, BETA, alpha * t_r) * gd_up[m])


# --- the pre-registered identify gate ----------------------------------------------

def test_gate_product_is_measurable_and_the_pair_is_degenerate():
    r = recovery_study(STRUCTURAL_STEADY, TRUTH, GRID, n_seeds=20, noise_cv=0.03)
    print(r.summary())
    for p in STRUCTURAL_STEADY.param_names:
        assert r.verdicts[p].get("measurable", 0) >= 18, r.summary()
        assert r.within_interval[p] >= 17, r.summary()

    alpha, t_r = 0.05, ALPHA_TR / 0.05
    sep = recovery_study(STRUCTURAL_STEADY_SEPARATE, [ETA_INF, BETA, alpha, t_r], GRID,
                         n_seeds=20, noise_cv=0.03)
    print(sep.summary())
    for p in ("eta_inf", "beta"):
        assert sep.verdicts[p].get("measurable", 0) >= 18, sep.summary()
    for p, partner in (("alpha", "t_r"), ("t_r", "alpha")):
        assert sep.verdicts[p].get("degenerate", 0) >= 18, sep.summary()
        partners = {f.params[p].partner for f in sep.fits if f.params[p].verdict == "degenerate"}
        assert partners == {partner}, partners
    # the product is still recovered by the separate fit, whatever split it landed on
    products = np.array([f.params["alpha"].value * f.params["t_r"].value for f in sep.fits])
    assert np.median(np.abs(products / ALPHA_TR - 1.0)) < 0.05, products
