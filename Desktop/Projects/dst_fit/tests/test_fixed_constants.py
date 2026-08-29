"""The four fixed constants of Wyart-Cates: defaults, binding, and refusal.

Two jobs here.

**Regression.** ``eta_s``, ``phi_0``, the divergence exponent ``n`` and the
activation sharpness ``beta`` were module-level constants until 2026-08-27. They
are now per-dataset arguments. Every gate result recorded before that change was
obtained at the old values, so the defaults are pinned by test: if a default
moves, those results stop reproducing and the change has to be deliberate.

**Behaviour.** Binding a constant must actually change the prediction (otherwise
the argument is decorative), must move the bounds that depend on it, and must
refuse a name the model does not have (otherwise a typo is silently ignored,
which is the failure mode the whole exercise exists to remove).
"""

from __future__ import annotations

import numpy as np
import pytest

from dst_fit.fitting import fit_model, synthetic_curve
from dst_fit.models import (
    ETA_S_WATER,
    HERSCHEL_BULKLEY,
    MARON_PIERCE,
    PHI_0_RCP,
    PHI_BOUND_MARGIN,
    WC_BETA,
    WC_EXPONENT,
    WYART_CATES,
    frictional_fraction,
)

SIGMA = np.logspace(-2, 4, 120)


# --------------------------------------------------------------------------
# regression: the defaults are the values every earlier gate ran at
# --------------------------------------------------------------------------
def test_module_defaults_are_unchanged():
    """Water at 25 degC, RCP, Maron-Pierce exponent, original Wyart-Cates form.

    Failing this test is not a bug in the test. It means a default moved, and
    every number in README.md was produced at the old one.
    """
    assert ETA_S_WATER == 8.9e-4
    assert PHI_0_RCP == 0.64
    assert WC_EXPONENT == 2.0
    assert WC_BETA == 1.0


def test_registry_models_carry_no_bound_constants():
    """Straight out of the registry, every model runs on the module defaults."""
    assert WYART_CATES.fixed == {}
    assert MARON_PIERCE.fixed == {}
    assert set(WYART_CATES.accepts_fixed()) == {"eta_s", "phi_0", "n", "beta"}
    assert "default" in WYART_CATES.describe_fixed()


def test_default_wyart_cates_matches_the_closed_form_it_had_before():
    """Independent re-derivation of eta_s (1 - phi/phi_J)^-2 with phi_0 = 0.64.

    Written out longhand rather than compared against a stored array, so this
    tests the formula rather than a previous run of the same code.
    """
    phi, phi_m, s_star = 0.50, 0.58, 5.0
    f = np.exp(-(s_star / SIGMA))
    phi_j = (1.0 - f) * 0.64 + f * 0.58
    expected = 8.9e-4 * (1.0 - phi / phi_j) ** -2.0
    np.testing.assert_allclose(WYART_CATES(SIGMA, phi, phi_m, s_star), expected, rtol=1e-12)


def test_default_bounds_are_unchanged():
    """0.639 is what the phi and phi_m caps were before phi_0 became settable."""
    assert WYART_CATES.bounds[0] == (0.10, 0.30, 1e-4)
    assert WYART_CATES.bounds[1] == pytest.approx((0.639, 0.639, 1e5))


# --------------------------------------------------------------------------
# binding a constant changes the answer
# --------------------------------------------------------------------------
def test_solvent_viscosity_scales_the_whole_curve():
    """eta_s is a prefactor, so tripling it triples eta everywhere. Exactly."""
    base = WYART_CATES(SIGMA, 0.50, 0.58, 5.0)
    tripled = WYART_CATES.with_fixed(eta_s=3 * ETA_S_WATER)(SIGMA, 0.50, 0.58, 5.0)
    np.testing.assert_allclose(tripled, 3.0 * base, rtol=1e-12)


def test_lowering_phi_0_raises_the_low_stress_viscosity():
    """The lubricated limit is set by phi_0 alone, so this is a one-line check."""
    lowered = WYART_CATES.with_fixed(phi_0=0.611)
    sigma = np.array([1e-6])  # f -> 0, purely frictionless
    assert lowered(sigma, 0.50, 0.55, 5.0)[0] > WYART_CATES(sigma, 0.50, 0.55, 5.0)[0]
    np.testing.assert_allclose(
        lowered(sigma, 0.50, 0.55, 5.0)[0],
        ETA_S_WATER * (1.0 - 0.50 / 0.611) ** -2.0,
        rtol=1e-9,
    )


def test_exponent_changes_the_slope_of_the_divergence():
    """eta ratio at two phi values goes as the exponent. n = 2.2 is Guy et al."""
    n = 2.2
    m = WYART_CATES.with_fixed(n=n)
    sigma = np.array([1e-6])
    a = m(sigma, 0.40, 0.55, 5.0)[0]
    b = m(sigma, 0.50, 0.55, 5.0)[0]
    ratio_expected = ((1 - 0.50 / PHI_0_RCP) / (1 - 0.40 / PHI_0_RCP)) ** -n
    assert b / a == pytest.approx(ratio_expected, rel=1e-9)


def test_beta_sharpens_the_transition_without_moving_its_centre():
    """f(sigma*) = 1/e for any beta; away from sigma* a larger beta is steeper."""
    for beta in (0.94, 1.0, 1.5):
        assert frictional_fraction(np.array([5.0]), 5.0, beta=beta)[0] == pytest.approx(
            np.exp(-1.0)
        )
    sharp = frictional_fraction(SIGMA, 5.0, beta=1.5)
    soft = frictional_fraction(SIGMA, 5.0, beta=0.7)
    below = SIGMA < 5.0
    assert np.all(sharp[below] <= soft[below] + 1e-12)


def test_binding_phi_0_moves_the_phi_and_phi_m_caps():
    """Both must stay below the frictionless packing or thickening runs backwards."""
    m = WYART_CATES.with_fixed(phi_0=0.546)
    cap = 0.546 - PHI_BOUND_MARGIN
    assert m.bounds[1][0] == pytest.approx(cap)
    assert m.bounds[1][1] == pytest.approx(cap)
    assert m.bounds[1][2] == WYART_CATES.bounds[1][2], "sigma* bound must not move"
    assert m.bounds[0] == WYART_CATES.bounds[0], "lower bounds must not move"


def test_binding_is_cumulative_and_leaves_the_original_alone():
    m = WYART_CATES.with_fixed(eta_s=2.4e-3).with_fixed(n=2.2, phi_0=0.611)
    assert m.fixed == {"eta_s": 2.4e-3, "n": 2.2, "phi_0": 0.611}
    assert WYART_CATES.fixed == {}
    assert WYART_CATES.bounds[1][1] == pytest.approx(0.639)


# --------------------------------------------------------------------------
# refusal
# --------------------------------------------------------------------------
def test_unknown_constant_is_refused_rather_than_ignored():
    """A silently ignored typo would put the wrong number in phi_m and say nothing."""
    with pytest.raises(ValueError, match="cannot fix"):
        WYART_CATES.with_fixed(eta_f=2.4e-3)
    with pytest.raises(ValueError, match="cannot fix"):
        MARON_PIERCE.with_fixed(phi_0=0.611)  # no phi_0 in the null model
    with pytest.raises(ValueError, match="no fixed constants"):
        HERSCHEL_BULKLEY.with_fixed(eta_s=1e-3)


def test_nonsense_constant_values_are_refused():
    for bad in (0.0, -1.0, np.nan, np.inf):
        with pytest.raises(ValueError, match="finite and > 0"):
            WYART_CATES.with_fixed(eta_s=bad)


def test_phi_0_below_the_lower_bound_is_refused():
    """phi_0 = 0.2 would leave phi_m no admissible range at all."""
    with pytest.raises(ValueError, match="no room"):
        WYART_CATES.with_fixed(phi_0=0.25)


# --------------------------------------------------------------------------
# the constants travel through the fitter
# --------------------------------------------------------------------------
def test_fitting_recovers_parameters_under_bound_constants():
    """Round trip at Guy et al.'s constants: generate, fit, get the truth back.

    The point is that ``fit_model`` uses the bound constants for both the
    prediction and the bounds. If it silently fell back to the defaults, phi and
    phi_m would come back shifted by roughly the phi_0 difference.
    """
    m = WYART_CATES.with_fixed(eta_s=2.4e-3, phi_0=0.611, n=2.2)
    truth = (0.52, 0.558, 200.0)
    sigma, eta = synthetic_curve(m, truth, np.logspace(0, 4, 60), noise=0.01, seed=3)
    fit = fit_model(m, sigma, eta, n_starts=12)
    got = dict(zip(fit.model.param_names, fit.params))
    assert got["phi"] == pytest.approx(0.52, rel=0.02)
    assert got["phi_m"] == pytest.approx(0.558, rel=0.02)
    assert got["sigma_star"] == pytest.approx(200.0, rel=0.10)
    assert fit.model.fixed == m.fixed


def test_fitting_under_default_constants_misses_when_the_truth_is_not_default():
    """The reason this refactor was necessary, stated as a test.

    Data generated at Guy et al.'s constants and fit with water/RCP defaults
    returns a phi_m that is wrong by far more than its own error bar. Before the
    refactor this was the only thing the package could do.
    """
    m = WYART_CATES.with_fixed(eta_s=2.4e-3, phi_0=0.611, n=2.2)
    truth = (0.52, 0.558, 200.0)
    sigma, eta = synthetic_curve(m, truth, np.logspace(0, 4, 60), noise=0.01, seed=3)
    wrong = fit_model(WYART_CATES, sigma, eta, n_starts=12)
    phi_m = wrong.params[wrong.model.param_names.index("phi_m")]
    err = wrong.stderr[wrong.model.param_names.index("phi_m")]
    assert abs(phi_m - 0.558) > 5 * err, (
        f"expected the default-constant fit to be biased well beyond its error bar, "
        f"got phi_m = {phi_m:.4f} +/- {err:.4f}"
    )


def test_synthetic_curve_respects_the_bound_constants():
    m = WYART_CATES.with_fixed(eta_s=2.4e-3)
    _, eta = synthetic_curve(m, (0.50, 0.58, 5.0), np.logspace(-1, 3, 20))
    _, eta_default = synthetic_curve(WYART_CATES, (0.50, 0.58, 5.0), np.logspace(-1, 3, 20))
    np.testing.assert_allclose(eta / eta_default, 2.4e-3 / ETA_S_WATER, rtol=1e-12)


def test_transition_coverage_landmarks_are_exact_for_any_beta():
    """f = 0.25 and 0.75 must land where the closed form says, whatever beta is.

    The landmarks are sigma*/ln(4)^(1/beta) and sigma*/ln(4/3)^(1/beta). Getting
    the 1/beta exponent wrong would silently make the bracketing verdict a
    function of a constant it should not depend on.
    """
    from dst_fit import transition_coverage

    for beta in (0.7, 0.94, 1.0, 1.5):
        edges = np.array([5.0 / np.log(4.0) ** (1 / beta),
                          5.0 / np.log(4.0 / 3.0) ** (1 / beta)])
        cov = transition_coverage(5.0, edges, beta=beta)
        assert cov.f_min == pytest.approx(0.25)
        assert cov.f_max == pytest.approx(0.75)
        assert not cov.brackets, "exactly at the edges is not bracketing"
        wider = transition_coverage(5.0, edges * np.array([0.5, 2.0]), beta=beta)
        assert wider.brackets


def test_scalar_and_zero_dim_input_survives_every_thickening_model():
    """A python float must work as well as an array, jammed or not."""
    from dst_fit.models import MARON_PIERCE, WYART_CATES_KNOWN_PHI

    assert np.isfinite(WYART_CATES(5.0, 0.50, 0.58, 5.0))
    assert not np.isfinite(WYART_CATES(1e9, 0.60, 0.58, 5.0))
    assert WYART_CATES_KNOWN_PHI.with_fixed(phi=0.50)(5.0, 0.58, 5.0) == pytest.approx(
        float(WYART_CATES(5.0, 0.50, 0.58, 5.0))
    )
    assert np.isfinite(MARON_PIERCE(5.0, 0.50, 0.58))
    np.testing.assert_allclose(
        WYART_CATES([1.0, 5.0], 0.50, 0.58, 5.0),
        WYART_CATES(np.array([1.0, 5.0]), 0.50, 0.58, 5.0),
    )
