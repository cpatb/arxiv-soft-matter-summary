"""
test_dilatancy_design.py -- validation of the beach-sand design calculator.

Every check is either a closed-form limit, a reproduction of a number in PROTOCOL.md, or
a conservation/scaling law.  Run: python -m pytest test_dilatancy_design.py -q
"""
import numpy as np
import pytest

import dilatancy_design as dd


# --------------------------------------------------------------------------- #
# Grain-size descriptors: analytic limits.
# --------------------------------------------------------------------------- #
def test_monodisperse_collapses_every_diameter_onto_d50():
    s = dd.size_descriptors(0.3e-3, Cu=1.0)
    for key in ("d10", "d60", "d_sauter"):
        np.testing.assert_allclose(s[key], 0.3e-3, rtol=1e-9)


def test_uniformity_coefficient_round_trips():
    for Cu in (1.2, 1.5, 2.0, 4.0):
        s = dd.sigma_ln_from_Cu(Cu)
        np.testing.assert_allclose(dd.uniformity_coefficient(s), Cu, rtol=1e-12)


def test_d60_over_d10_equals_Cu_by_construction():
    for Cu in (1.5, 2.0, 3.0):
        s = dd.size_descriptors(0.3e-3, Cu)
        np.testing.assert_allclose(s["d60"] / s["d10"], Cu, rtol=1e-10)


def test_sauter_is_between_d10_and_d50():
    """Surface-area weighting pulls the effective diameter below the median but nowhere
    near as far as the 10th percentile."""
    s = dd.size_descriptors(0.3e-3, Cu=2.0)
    assert s["d10"] < s["d_sauter"] < s["d50"]


def test_sauter_formula_matches_numerical_integration():
    """d_32 = 1 / E_vol[1/d] for a lognormal volume distribution, checked by quadrature
    rather than trusting the closed form."""
    d50, s = 0.3e-3, 0.45
    x = np.linspace(-6, 6, 200001)
    d = d50 * np.exp(s * x)
    pdf = np.exp(-0.5 * x**2) / np.sqrt(2 * np.pi)
    inv_mean = np.trapezoid(pdf / d, x)
    np.testing.assert_allclose(dd.sauter_diameter(d50, s), 1.0 / inv_mean, rtol=1e-6)


def test_bad_sorting_raises():
    with pytest.raises(ValueError):
        dd.sigma_ln_from_Cu(0.5)


# --------------------------------------------------------------------------- #
# Permeability: reproduce PROTOCOL.md numbers.
# --------------------------------------------------------------------------- #
def test_kozeny_carman_reproduces_the_protocol_shorthand():
    """k ~ 1e-3 d^2 at porosity 0.40, and 4e-11 m^2 for 0.2 mm sand."""
    np.testing.assert_allclose(dd.permeability_kozeny_carman(1.0) / 1.0**2,
                               9.876e-4, rtol=1e-3)
    np.testing.assert_allclose(dd.permeability_kozeny_carman(0.2e-3), 3.95e-11, rtol=0.02)


def test_permeability_scales_as_d_squared():
    a = dd.permeability_kozeny_carman(0.2e-3)
    b = dd.permeability_kozeny_carman(0.4e-3)
    np.testing.assert_allclose(b / a, 4.0, rtol=1e-12)


def test_porosity_and_solid_fraction_are_not_interchangeable():
    """Regression for the notation collision that cost a factor of 7.6 here.  Passing the
    solid fraction where the porosity belongs must give a visibly different answer."""
    as_porosity = dd.permeability_kozeny_carman(0.2e-3, porosity=0.40)
    as_solid = dd.permeability_kozeny_carman(0.2e-3, porosity=0.60)
    assert as_solid / as_porosity > 5.0


def test_out_of_range_porosity_raises():
    for bad in (0.0, 1.0, -0.1, 1.5):
        with pytest.raises(ValueError):
            dd.permeability_kozeny_carman(0.2e-3, porosity=bad)


# --------------------------------------------------------------------------- #
# Suction, hardening and cavitation.
# --------------------------------------------------------------------------- #
def test_suction_reproduces_the_protocol_table():
    """~100 kPa for 0.2 mm and ~4 kPa for 1.0 mm at V = 2 m/s, D = 20 mm, psi = 0.1."""
    fine = dd.darcy_reynolds_suction(2.0, 0.2e-3, cap_at_cavitation=False)
    coarse = dd.darcy_reynolds_suction(2.0, 1.0e-3, cap_at_cavitation=False)
    np.testing.assert_allclose(fine / 1e3, 101.0, rtol=0.05)
    np.testing.assert_allclose(coarse / 1e3, 4.0, rtol=0.05)


def test_cavitation_cap_is_applied_and_can_be_switched_off():
    capped = dd.darcy_reynolds_suction(10.0, 0.1e-3, cap_at_cavitation=True)
    free = dd.darcy_reynolds_suction(10.0, 0.1e-3, cap_at_cavitation=False)
    assert capped < free
    np.testing.assert_allclose(capped, dd.P_ATM - dd.P_VAP, rtol=1e-12)


def test_cavitation_speed_scales_as_d_squared():
    a = dd.cavitation_speed(0.1e-3)
    b = dd.cavitation_speed(0.2e-3)
    np.testing.assert_allclose(b / a, 4.0, rtol=1e-12)


def test_cavitation_speed_is_where_the_uncapped_suction_hits_the_ceiling():
    d = 0.2e-3
    V = dd.cavitation_speed(d)
    np.testing.assert_allclose(dd.darcy_reynolds_suction(V, d, cap_at_cavitation=False),
                               dd.P_ATM - dd.P_VAP, rtol=1e-9)


def test_fine_sand_hardens_and_coarse_sand_drains():
    assert dd.hardening_ratio(2.0, 0.2e-3, cap_at_cavitation=False) > 1.0
    assert dd.hardening_ratio(2.0, 1.0e-3, cap_at_cavitation=False) < 1.0


# --------------------------------------------------------------------------- #
# The two design findings.
# --------------------------------------------------------------------------- #
def test_q1_diameters_are_distinguishable_at_realistic_sorting():
    """The design question for Q1: at C_u = 2 the d10 and d50 predictions must differ by
    far more than measurement scatter, or the experiment cannot discriminate."""
    q = dd.q1_discrimination(0.3e-3, Cu=2.0)
    assert q["ratio_d10_over_d50"] > 1.5
    assert dd.minimum_Cu_to_discriminate(0.3e-3, 0.15) < 1.5


def test_q1_ratio_is_unity_for_monodisperse():
    q = dd.q1_discrimination(0.3e-3, Cu=1.0)
    np.testing.assert_allclose(q["ratio_d10_over_d50"], 1.0, rtol=1e-9)
    np.testing.assert_allclose(q["ratio_sauter_over_d50"], 1.0, rtol=1e-9)


def test_strader_bounds_psi_well_below_the_protocol_value():
    """The literature constraint on the loosest number in the chain."""
    b = dd.psi_bound_from_strader()
    assert b["psi_max"] < 0.02
    assert b["factor_below_protocol"] > 5.0


def test_consistency_window_across_straders_range_is_empty():
    """The finding: no single psi*D fits both ends of a 16x grain-size sweep."""
    w = dd.consistency_window()
    assert not w["consistent"]
    assert w["shortfall_factor"] > 10.0


def test_consistency_window_closes_for_a_narrow_size_range():
    """Sanity: the contradiction is caused by the width of the sweep, not by a coding
    error, so a narrow range must be consistent."""
    w = dd.consistency_window(d_fine=200e-6, d_coarse=210e-6, V_max=0.5)
    assert w["consistent"]


def test_regime_match_depends_only_on_D_over_V():
    a = dd.match_field_regime(0.2e-3, D_lab=0.02, V_lab=2.0)
    b = dd.match_field_regime(0.2e-3, D_lab=0.04, V_lab=4.0)
    np.testing.assert_allclose(a["mismatch_factor"], b["mismatch_factor"], rtol=1e-12)


def test_proposed_rig_is_an_order_of_magnitude_off_the_footstep_regime():
    m = dd.match_field_regime(0.2e-3, D_lab=0.020, V_lab=2.0)
    np.testing.assert_allclose(m["mismatch_factor"], 10.0, rtol=1e-9)
    np.testing.assert_allclose(m["V_lab_required"], 0.2, rtol=1e-9)
