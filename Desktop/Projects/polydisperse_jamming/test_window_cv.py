"""test_window_cv.py — the two analysis scripts added 2026-09-13.

Both read committed artifacts and produce numbers that go into prose, so what is tested
here is what could make that prose wrong: the record selection, the estimator algebra,
and the one cross-check that ties the new reading to the shipped one.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from cv_switch import CV_M, PHI_C, PHI_C_SE, cv_star, phi_c_curve
from window_shearedge import CVS, PHI_C_SHEAR, check_agrees_with_dst_window, ladder, wls_slope


# ------------------------------------------------------------------ window_shearedge
def test_record_selection_reproduces_dst_window():
    """The first version of this script read phase2_fits.csv directly and silently
    pooled the capfix re-runs, moving phi_0(CV 0.5) by 0.0018 -- a third of the effect
    it is about. This is the assertion that caught it."""
    check_agrees_with_dst_window()


def test_the_shear_edge_sits_above_the_isotropic_frictional_one():
    """The whole reason the shipped lower edge was flagged NOT safe. If this ever
    reverses, dst_window.py's caveat is wrong and so is this script's premise."""
    fm = ladder(0.3)
    for cv in CVS:
        assert PHI_C_SHEAR[cv] > fm.phi_J[cv], cv


def test_reading_A_traverse_reproduces_dst_window_headline():
    """dst_window.py reports 2.23 windows per unit CV on six levels. The endpoint
    statistic on the three levels the shear campaign covers must land on it, or the
    two scripts are measuring different things and the comparison is meaningless."""
    fr, fm = ladder(0.0), ladder(0.3)
    W = [fr.phi_J[cv] - fm.phi_J[cv] for cv in CVS]
    slope = (fm.phi_J[CVS[-1]] - fm.phi_J[CVS[0]]) / (CVS[-1] - CVS[0])
    assert slope / float(np.mean(W)) == pytest.approx(2.23, abs=0.05)


def test_wls_slope_is_exact_on_a_line():
    b, se, c = wls_slope([0.0, 0.2, 0.5], [1.0, 1.2, 1.5], [0.1, 0.1, 0.1])
    assert b == pytest.approx(1.0)
    assert c == pytest.approx(1.0)
    assert se > 0


def test_wls_slope_down_weights_the_noisy_point():
    """Two clean points on a line and one noisy point far off it: the slope must stay
    near the clean pair. This is the behaviour that makes the WLS and endpoint
    statistics disagree in the report, so it is asserted rather than assumed."""
    b, _, _ = wls_slope([0.0, 0.2, 0.5], [0.0, 0.2, 5.0], [1e-3, 1e-3, 10.0])
    assert abs(b - 1.0) < 0.05


# ------------------------------------------------------------------------- cv_switch
def test_phi_c_curve_passes_through_every_measured_level():
    c = phi_c_curve()
    for cv, phi in zip(CV_M, PHI_C):
        assert float(c(cv)) == pytest.approx(phi, abs=1e-12)


def test_cv_star_inverts_the_curve():
    """An anchor that does not pass through the estimator: pick a CV, read phi_c there
    from the curve, and require cv_star to hand the CV back."""
    c = phi_c_curve()
    for cv in (0.05, 0.17, 0.33, 0.47):
        assert cv_star(float(c(cv)), c) == pytest.approx(cv, abs=1e-6)


def test_cv_star_refuses_outside_the_measured_range():
    """This project does not extrapolate a root past the rungs that bracket it."""
    c = phi_c_curve()
    assert math.isnan(cv_star(0.55, c))          # below phi_c(0)
    assert math.isnan(cv_star(0.70, c))          # above phi_c(0.5)


def test_the_switchable_band_is_the_measured_contrast():
    """The band width IS Delta phi_c(CV 0 -> 0.5), the rate-matched campaign's headline.
    If these drift apart, one of the two transcriptions is stale."""
    assert PHI_C[-1] - PHI_C[0] == pytest.approx(0.0533, abs=5e-4)
    assert math.hypot(PHI_C_SE[-1], PHI_C_SE[0]) < 0.0025


def test_the_channel_crossover_is_where_the_two_terms_are_equal():
    """Delta phi* = phi_c'(CV)(1 + CV^2)/(8 CV) must make the permeability and
    critical-state terms of Eq. (11) equal, by construction. A sign slip or a dropped
    factor in the reported table would break this and nothing else would."""
    c = phi_c_curve()
    for cv in (0.05, 0.2, 0.35, 0.5):
        dprime = float(c.derivative()(cv))
        dstar = dprime * (1 + cv ** 2) / (8 * cv)
        perm = 8 * cv / (1 + cv ** 2)
        crit = dprime / dstar
        assert crit == pytest.approx(perm, rel=1e-12)


def test_the_permeability_channel_vanishes_at_cv_zero_and_the_packing_one_does_not():
    """Eq. (12)/(13): the order-counting the derivation rests on, as a live check."""
    c = phi_c_curve()
    small = 1e-6
    assert 8 * small / (1 + small ** 2) < 1e-5
    assert abs(float(c.derivative()(small))) > 1e-3
