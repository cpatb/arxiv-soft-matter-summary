"""Invariant tests for the constitutive models.

Each test states a property the physics requires, not a number a previous run
happened to produce.
"""

from __future__ import annotations

import numpy as np
import pytest

from dst_fit.models import (
    ELLIS,
    HERSCHEL_BULKLEY,
    MARON_PIERCE,
    WYART_CATES,
    frictional_fraction,
    get_model,
    jamming_fraction,
    thickening_branch,
)

SIGMA = np.logspace(-3, 4, 200)


# --------------------------------------------------------------------------
# Wyart-Cates
# --------------------------------------------------------------------------
def test_frictional_fraction_spans_zero_to_one_monotonically():
    """Non-decreasing everywhere, strictly increasing wherever float64 can say so.

    exp(-(sigma*/sigma)) underflows to exactly 0 below sigma ~ sigma*/745, so
    consecutive values tie there. Asserting strict monotonicity across the whole
    range would be demanding a difference smaller than the type can represent.
    """
    f = frictional_fraction(SIGMA, sigma_star=5.0)
    assert f[0] < 1e-6 and f[-1] > 0.999
    assert np.all((f >= 0) & (f <= 1))
    assert np.all(np.diff(f) >= 0)
    resolvable = f > 0
    assert np.all(np.diff(f[resolvable]) > 0)


def test_jamming_fraction_interpolates_the_two_endpoints():
    assert jamming_fraction(0.0, phi_0=0.64, phi_m=0.58) == pytest.approx(0.64)
    assert jamming_fraction(1.0, phi_0=0.64, phi_m=0.58) == pytest.approx(0.58)


def test_wyart_cates_has_no_thinning_branch():
    """The reason models.py is stress-parameterised and the branch is cut.

    Wyart-Cates viscosity increases monotonically with stress. Real corn starch
    thins first, so a whole-sweep fit would be fitting a shape the model cannot
    make.
    """
    eta = WYART_CATES(SIGMA, 0.50, 0.58, 5.0)[np.isfinite(WYART_CATES(SIGMA, 0.50, 0.58, 5.0))]
    # The physics claim: it never thins.
    assert np.all(np.diff(eta) >= 0)
    # And it really does thicken, rather than being flat to within rounding.
    assert eta[-1] / eta[0] > 2.0
    # Ties exist at the low-stress end because eta's dependence on sigma there
    # is below the spacing of float64, not because the model thins. Asserting
    # strict monotonicity pointwise would be a claim about the dtype.
    assert np.count_nonzero(np.diff(eta) > 0) > eta.size // 2


def test_wyart_cates_jams_above_the_frictional_packing():
    """phi above phi_m must jam at high stress, and not be a large finite number."""
    eta = WYART_CATES(SIGMA, 0.60, 0.58, 5.0)
    assert np.isfinite(eta[0]), "should still flow in the lubricated limit"
    assert not np.isfinite(eta[-1]), "must diverge once phi exceeds phi_J"


def test_wyart_cates_viscosity_rises_with_volume_fraction():
    a = WYART_CATES(SIGMA, 0.45, 0.58, 5.0)
    b = WYART_CATES(SIGMA, 0.50, 0.58, 5.0)
    assert np.all(b >= a)


def test_sigma_star_sets_where_the_transition_sits():
    """Doubling sigma* moves the half-transition stress by the same factor."""
    for s_star in (2.0, 20.0):
        f = frictional_fraction(SIGMA, s_star)
        half = SIGMA[int(np.argmin(np.abs(f - np.exp(-1.0))))]
        assert half == pytest.approx(s_star, rel=0.05)


# --------------------------------------------------------------------------
# Herschel-Bulkley
# --------------------------------------------------------------------------
def test_herschel_bulkley_is_the_exact_inverse_of_its_own_definition():
    tau_y, k, n = 2.0, 1.5, 0.6
    gdot = np.logspace(-2, 3, 50)
    sigma = tau_y + k * gdot**n
    eta = HERSCHEL_BULKLEY(sigma, tau_y, k, n)
    np.testing.assert_allclose(sigma / eta, gdot, rtol=1e-10)


def test_herschel_bulkley_does_not_flow_below_yield():
    eta = HERSCHEL_BULKLEY(np.array([0.5, 1.0, 1.99, 2.01]), 2.0, 1.5, 0.6)
    assert np.all(~np.isfinite(eta[:3]))
    assert np.isfinite(eta[3])


# --------------------------------------------------------------------------
# Ellis and Maron-Pierce
# --------------------------------------------------------------------------
def test_ellis_thins_between_its_two_plateaus():
    eta = ELLIS(SIGMA, 10.0, 0.05, 3.0, 1.4)
    assert np.all(np.diff(eta) < 0)
    assert eta[0] == pytest.approx(10.0, rel=1e-3)
    assert eta[-1] == pytest.approx(0.05, rel=1e-2)


def test_maron_pierce_is_rate_independent():
    eta = MARON_PIERCE(SIGMA, 0.50, 0.58)
    assert np.ptp(eta) == 0.0


def test_maron_pierce_jams_at_phi_m():
    assert not np.isfinite(MARON_PIERCE(SIGMA, 0.58, 0.58)[0])


# --------------------------------------------------------------------------
# plumbing
# --------------------------------------------------------------------------
def test_thickening_branch_starts_at_the_viscosity_minimum():
    eta = np.array([5.0, 3.0, 1.0, 2.0, 8.0])
    assert thickening_branch(eta) == slice(2, 5)
    assert thickening_branch(eta, pad=1) == slice(3, 5)
    assert thickening_branch(eta, pad=-10) == slice(0, 5)


def test_thickening_branch_rejects_empty():
    with pytest.raises(ValueError, match="empty"):
        thickening_branch(np.array([]))


def test_check_params_rejects_out_of_bounds():
    with pytest.raises(ValueError, match="outside bounds"):
        WYART_CATES.check_params((0.50, 0.58, -1.0))
    with pytest.raises(ValueError, match="expected 3 params"):
        WYART_CATES.check_params((0.50, 0.58))


def test_get_model_names_the_alternatives_when_wrong():
    with pytest.raises(KeyError, match="available"):
        get_model("carreau")
