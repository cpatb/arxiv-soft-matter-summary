"""Tests for the tier-0 calibration layer.

Every physical claim here is anchored on something outside the code under test:
a hand-computed Peclet number, an analytic inverse, the model's own exact gauge
invariance, or a discriminating negative case. A round trip through the fitter
would certify the optimiser, not the transfer, so there are none.

Run:  python -m pytest calibration/test_calibration.py -q
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from brownian import K_B, invert_plateau, peclet, plateau_admissible  # noqa: E402
from propagate_impact import (  # noqa: E402
    band_consequence,
    flip_window,
    impact_stress,
    saturation,
)
from transfer import (  # noqa: E402
    jamming_stress,
    jams_under_stress,
    transfer_map,
    verdict_band,
)
from wyart_cates import jamming_fraction, relative_viscosity  # noqa: E402


# --------------------------------------------------------------------------
# Peclet: anchored on a hand computation, not on the function's own output
# --------------------------------------------------------------------------
def test_peclet_matches_a_hand_computation():
    # 6 pi (2.83e-3)(1.885e-6)^3 (0.0505) / (1.380649e-23 * 292.15)
    #   numerator   = 6 pi * 2.83e-3 * 6.6976e-18 * 0.0505 = 1.8054e-20
    #   denominator = 1.380649e-23 * 292.15                = 4.0341e-21
    #   ratio       = 4.475
    got = peclet(0.0505, 1.885e-6, 2.83e-3, 292.15)
    assert got == pytest.approx(4.475, rel=2e-3)


def test_peclet_is_cubic_in_radius():
    """Doubling the radius must multiply Pe by exactly 8."""
    a = peclet(1.0, 1e-6, 1e-3)
    b = peclet(1.0, 2e-6, 1e-3)
    assert b / a == pytest.approx(8.0, rel=1e-12)


def test_peclet_rejects_unphysical_inputs():
    for bad in ({"radius": 0.0}, {"radius": -1e-6}, {"eta_s": 0.0},
                {"temperature": -1.0}):
        kw = {"radius": 1e-6, "eta_s": 1e-3, "temperature": 300.0} | bad
        with pytest.raises(ValueError):
            peclet(1.0, **kw)


def test_boltzmann_constant_is_the_codata_value():
    assert K_B == 1.380649e-23


def test_guy_404_is_brownian_and_guy_3770_is_not():
    """The discriminating negative case the admission rule exists for.

    Same chemistry, same lab, same figure; 9.3x apart in size. The rule must
    separate them, and by many decades, or it is not measuring anything.
    """
    small = peclet(2.48e-3, 202e-9, 2.4e-3, 292.15)
    large = peclet(0.0505, 1.885e-6, 2.83e-3, 292.15)
    assert small < 1e-3
    assert large > 1.0
    assert large / small > 1e3


# --------------------------------------------------------------------------
# plateau inversion: analytic, so it can be checked against the forward model
# --------------------------------------------------------------------------
def test_invert_plateau_is_the_exact_inverse_of_the_frictionless_divergence():
    phi_0 = 0.6123
    for phi in (0.20, 0.35, 0.50, 0.58):
        eta_r = relative_viscosity(phi, 0.0, phi_0, 0.5)   # f = 0 -> phi_J = phi_0
        back = invert_plateau(np.array([float(eta_r)]), phi_0, n=2.0)
        # relative_viscosity is hardcoded at n = 2, so invert at n = 2 to match
        assert back[0] == pytest.approx(phi, abs=1e-12)


def test_invert_plateau_refuses_below_the_solvent_viscosity():
    out = invert_plateau(np.array([0.5, 1.0, 4.0]), 0.64, 2.0)
    assert np.isnan(out[0]) and np.isnan(out[1])
    assert np.isfinite(out[2])


def test_plateau_admissible_reports_rather_than_dropping():
    stress = np.array([0.01, 0.1, 1.0])
    rate = np.array([1e-3, 1e-2, 1e-1])
    visc = stress / rate
    v = plateau_admissible(stress, rate, visc, radius=1e-9, eta_s=1e-3,
                           temperature=300.0)
    assert not v.admitted
    assert "Brownian" in v.reason
    assert v.pe > 0                      # the number is carried, not discarded


def test_plateau_admissible_checks_shapes():
    with pytest.raises(ValueError):
        plateau_admissible(np.zeros(3), np.zeros(4), np.zeros(3),
                           radius=1e-6, eta_s=1e-3)


# --------------------------------------------------------------------------
# the transfer map
# --------------------------------------------------------------------------
def test_transfer_map_returns_the_identity_on_identical_inputs():
    phi = np.array([0.40, 0.45, 0.50, 0.55])
    tm = transfer_map(phi, phi, n_boot=200)
    assert tm.slope == pytest.approx(1.0, abs=1e-12)
    assert tm.intercept == pytest.approx(0.0, abs=1e-12)
    assert tm.is_identity


def test_transfer_map_recovers_a_known_affine_law():
    phi = np.array([0.40, 0.45, 0.50, 0.55, 0.60])
    tm = transfer_map(phi, 1.2 * phi - 0.13, n_boot=200)
    assert tm.slope == pytest.approx(1.2, abs=1e-10)
    assert tm.intercept == pytest.approx(-0.13, abs=1e-10)
    assert tm.phi_x == pytest.approx(0.13 / 1.2, abs=1e-10)
    assert not tm.is_identity


def test_phi_x_is_invariant_under_the_models_gauge():
    """phi_x = -b/a must not move when the assumed phi_0 is rescaled.

    Wyart-Cates is homogeneous of degree zero in (phi, phi_m, phi_0), so
    assuming a phi_0 that is c times too large scales every inferred phi by c.
    Both the slope and the intercept scale with it; their ratio must not.
    """
    phi = np.array([0.40, 0.45, 0.50, 0.55, 0.60])
    inferred = 1.2 * phi - 0.13
    base = transfer_map(phi, inferred, n_boot=200)
    for c in (0.85, 1.0, 1.17):
        scaled = transfer_map(phi, c * inferred, n_boot=200)
        assert scaled.slope == pytest.approx(c * base.slope, rel=1e-10)
        assert scaled.intercept == pytest.approx(c * base.intercept, rel=1e-10)
        assert scaled.phi_x == pytest.approx(base.phi_x, rel=1e-10)


def test_transfer_map_refuses_too_few_curves_and_a_degenerate_axis():
    with pytest.raises(ValueError):
        transfer_map([0.4, 0.5], [0.4, 0.5])
    with pytest.raises(ValueError):
        transfer_map([0.5, 0.5, 0.5], [0.4, 0.5, 0.6])


# --------------------------------------------------------------------------
# jamming: the yes/no a solver asks
# --------------------------------------------------------------------------
def test_jamming_stress_inverts_the_activation_by_hand():
    # phi_0 = 0.61, phi_m = 0.55, phi = 0.58 -> f_jam = 0.03/0.06 = 0.5
    # sigma_jam = sigma* / (-ln 0.5) = 10 / 0.693147 = 14.4270
    assert jamming_stress(0.58, 0.55, 0.61, 10.0, 1.0) == pytest.approx(14.4270, rel=1e-5)


def test_jamming_stress_is_gauge_invariant():
    """Rescale the whole packing triple; sigma_jam must not move at all."""
    base = jamming_stress(0.58, 0.55, 0.61, 10.0, 1.0)
    for c in (0.8, 1.0, 1.05):
        assert jamming_stress(0.58 * c, 0.55 * c, 0.61 * c, 10.0, 1.0) == \
            pytest.approx(base, rel=1e-12)


def test_a_material_below_phi_m_never_jams():
    assert jamming_stress(0.50, 0.55, 0.61, 10.0, 1.0) == float("inf")
    assert not jams_under_stress(0.50, 0.55)


def test_a_material_already_above_phi_0_is_jammed_at_rest():
    assert jamming_stress(0.62, 0.55, 0.61, 10.0, 1.0) == 0.0


def test_jamming_stress_rejects_phi_0_below_phi_m():
    with pytest.raises(ValueError):
        jamming_stress(0.50, 0.62, 0.61, 10.0, 1.0)


def test_jamming_fraction_at_the_jamming_stress_reproduces_phi():
    """Cross-check against the closure itself: at sigma_jam, phi_J must equal phi."""
    phi, phi_m, phi_0, ss, beta = 0.58, 0.55, 0.61, 10.0, 1.0
    s = jamming_stress(phi, phi_m, phi_0, ss, beta)
    f = np.exp(-((ss / s) ** beta))
    assert float(jamming_fraction(f, phi_0, phi_m)) == pytest.approx(phi, abs=1e-12)


# --------------------------------------------------------------------------
# impact propagation
# --------------------------------------------------------------------------
def test_impact_stress_is_rho_k_u_squared():
    lo, hi = impact_stress(u0=2.0, rho=1000.0, k=(2.0, 12.0))
    assert lo == pytest.approx(8000.0)
    assert hi == pytest.approx(48000.0)


def test_saturation_is_judged_against_the_band_not_against_one():
    """A residual that is small in absolute terms can still be large in band units."""
    wide = saturation(500.0, 1.0, 0.61, 0.56, band_width=0.30)
    narrow = saturation(500.0, 1.0, 0.61, 0.56, band_width=0.003)
    assert wide.phi_j_residual == pytest.approx(narrow.phi_j_residual, rel=1e-12)
    assert wide.saturated
    assert not narrow.saturated


def test_saturation_requires_a_positive_band():
    with pytest.raises(ValueError):
        _ = saturation(10.0, 1.0, 0.61, 0.56, band_width=0.0).saturated


def test_saturation_increases_as_sigma_star_falls():
    residuals = [saturation(ss, 1.0, 0.61, 0.56, 0.03).phi_j_residual
                 for ss in (1.0, 10.0, 100.0, 1000.0)]
    assert residuals == sorted(residuals)


def test_band_consequence_flips_only_inside_the_band():
    lo, hi = 0.53, 0.59
    inside = band_consequence(0.55, lo, hi, 0.6123, 10.0, 1.0)
    below = band_consequence(0.50, lo, hi, 0.6123, 10.0, 1.0)
    above = band_consequence(0.605, lo, hi, 0.6123, 10.0, 1.0)
    assert inside.verdict_flips
    assert not below.verdict_flips
    assert not above.verdict_flips
    assert not np.isfinite(inside.eta_r_at_lo)      # jammed on the low-phi_m end
    assert np.isfinite(inside.eta_r_at_hi)


def test_flip_window_is_the_band_and_is_reported_in_tolerance_units():
    lo, hi, width, xtol = flip_window(0.5876, 0.5297)   # deliberately reversed
    assert (lo, hi) == (0.5297, 0.5876)
    assert width == pytest.approx(0.0579, abs=1e-12)
    assert xtol == pytest.approx(0.0579 / 0.002, rel=1e-12)


def test_verdict_band_needs_two_values():
    with pytest.raises(ValueError):
        verdict_band([0.55])


# --------------------------------------------------------------------------
# the source labs' steadiness rules
# --------------------------------------------------------------------------
def test_richards_red_set_is_still_what_gate_n8_found():
    """Parsed from the deposit, so it breaks loudly if the deposit changes."""
    from steady import richards_jammed_labels
    assert richards_jammed_labels() == ("0.46", "0.48", "0.50")


def test_guy_steadiness_threshold_is_applied_to_the_true_volume_fraction():
    """The sheet labelled '56%' is phi = 0.5562, which the rule admits."""
    from steady import is_steady
    assert is_steady("guy_3770", 0.5363)
    assert is_steady("guy_404", 0.5562)      # label says 56%, phi says 0.5562
    assert not is_steady("guy_404", 0.5696)
    assert not is_steady("guy_3770", 0.5613)


def test_cornstarch_steadiness_comes_from_the_parsed_colours():
    from steady import is_steady
    assert is_steady("richards_cornstarch", 0.44)
    assert not is_steady("richards_cornstarch", 0.46)


def test_an_unknown_dataset_raises_rather_than_defaulting_to_admitting():
    from steady import is_steady
    with pytest.raises(KeyError):
        is_steady("some_other_paper", 0.5)
