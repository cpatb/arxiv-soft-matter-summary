"""Unit 0.6: needle injectability.

Anchors outside the module: the closed forms in ``derivations/buckingham_reiner.md``
(SymPy integral vs closed form, three named limits) and the hand arithmetic shown beside
every number here. The quadrature route and the closed form are independent
implementations and are checked against each other.

Pre-registered gate (2026-09-05, before the run): on a synthetic Herschel-Bulkley curve
(truth tau_y = 20 Pa, k = 2 Pa s^n, n = 0.5 on gd = logspace(-1, 3, 30), 3 % lognormal
noise, seed 0), ``force_with_uq`` at Q = 1e-8 m^3/s through a 0.1 mm x 20 mm needle and
a 2.3 mm-radius barrel must (i) contain the true-parameter force in its profile interval,
(ii) contain it in its sampled interval, (iii) report the wall shear rate INSIDE the
fitted envelope when it is, and EXTRAPOLATED for a Q that pushes gd_w past 1e3 1/s, and
(iv) give sensitivities to tau_y and k that sum to 1 within 1e-6 (homogeneity of the
stress balance).

First run: (iii) failed because the pre-registered Q was wrong, not the flag: at
Q = 1e-8 the wall rate is 1.65e4 1/s, sixteen times the envelope top, so the flag
correctly said EXTRAPOLATED. I had written "inside" without computing gd_w. The gate is
unchanged (inside when inside, extrapolated when not); the inside case now uses
Q = 2e-10 (gd_w = 403 1/s, computed below) and the extrapolated case Q = 1e-8.
"""

from __future__ import annotations

import numpy as np
import pytest

from biorheo.apps import (
    Needle,
    Syringe,
    flow_rate,
    flow_rate_gnf,
    flow_rate_hb,
    force_with_uq,
    injection_force,
    plug_radius,
    pressure_for_flow_rate,
    sensitivities,
    wall_shear_rate,
    wall_stress,
)
from biorheo.constitutive import get_model
from biorheo.identify import fit_rate_curve
from biorheo.identify.synthetic import synthetic_rate_curve

HB = get_model("herschel_bulkley")
PL = get_model("power_law")
NEWT = get_model("newtonian")
CROSS = get_model("cross")
CASSON = get_model("casson")

NEEDLE = Needle(radius=1e-4, length=0.02)
SYRINGE = Syringe(barrel_radius=2.3e-3)


# --- closed form against hand arithmetic --------------------------------------------

@pytest.mark.parametrize("dp, R, L, mu, expected", [
    # pi R^4 dp / (8 mu L), each worked in derivations/buckingham_reiner.py
    (1e5, 1e-4, 0.02, 1.0, 1.9634954e-10),
    (5e4, 2e-4, 0.04, 0.1, 7.8539816e-9),
    (1e6, 5e-5, 0.01, 10.0, 2.4543693e-11),
])
def test_newtonian_poiseuille_three_points(dp, R, L, mu, expected):
    assert flow_rate_hb(dp, R, L, 0.0, mu, 1.0) == pytest.approx(expected, rel=1e-7)
    assert flow_rate_gnf(dp, R, L, NEWT, [mu]) == pytest.approx(expected, rel=1e-7)


def test_power_law_hand_value():
    # tau_w = 250 Pa, (tau_w/k)^(1/n) = 125^2 = 15625, n/(3n+1) = 0.2:
    # Q = pi * 0.2 * 1e-12 * 15625 = 9.8174770e-9
    assert flow_rate_hb(1e5, 1e-4, 0.02, 0.0, 2.0, 0.5) == pytest.approx(9.8174770e-9, rel=1e-7)
    assert flow_rate_gnf(1e5, 1e-4, 0.02, PL, [2.0, 0.5]) == pytest.approx(9.8174770e-9, rel=1e-7)


def test_bingham_buckingham_reiner_hand_value():
    # tau_w = 250, tau_y = 50 -> phi = 0.2; bracket 1 - 4/3 0.2 + 1/3 0.2^4 = 0.7338667
    # Q = 1.9634954e-10 * 0.7338667 = 1.4409440e-10
    assert flow_rate_hb(1e5, 1e-4, 0.02, 50.0, 1.0, 1.0) == pytest.approx(1.4409440e-10, rel=1e-7)
    assert plug_radius(1e5, 0.02, 50.0) / 1e-4 == pytest.approx(0.2)


def test_herschel_bulkley_hand_value_and_quadrature_agree():
    # worked in the derivation: 5.5627134e-9
    q_closed = flow_rate_hb(1e5, 1e-4, 0.02, 50.0, 2.0, 0.5)
    assert q_closed == pytest.approx(5.5627134e-9, rel=1e-7)
    assert flow_rate_gnf(1e5, 1e-4, 0.02, HB, [50.0, 2.0, 0.5]) == pytest.approx(q_closed, rel=1e-8)


@pytest.mark.parametrize("params", [(0.0, 0.5, 0.8), (5.0, 0.3, 0.4), (200.0, 10.0, 1.3)])
def test_quadrature_matches_closed_form_across_n(params):
    dp = 4e5
    q = flow_rate_hb(dp, 2e-4, 0.03, *params)
    assert q > 0
    assert flow_rate_gnf(dp, 2e-4, 0.03, HB, params) == pytest.approx(q, rel=1e-8)


def test_no_flow_below_yield_and_the_onset_is_continuous():
    # tau_w = dp R/(2L) = dp * 0.0025; tau_y = 50 -> dp* = 20000 Pa
    assert flow_rate_hb(19999.0, 1e-4, 0.02, 50.0, 1.0, 1.0) == 0.0
    assert flow_rate_gnf(19999.0, 1e-4, 0.02, HB, [50.0, 1.0, 1.0]) == 0.0
    assert 0 < flow_rate_hb(20001.0, 1e-4, 0.02, 50.0, 1.0, 1.0) < 1e-16


def test_flow_rate_dispatch_and_pressure_inversion_round_trip():
    q = flow_rate(1e5, NEEDLE, HB, [50.0, 2.0, 0.5])
    assert q == pytest.approx(flow_rate_hb(1e5, 1e-4, 0.02, 50.0, 2.0, 0.5))
    assert pressure_for_flow_rate(q, NEEDLE, HB, [50.0, 2.0, 0.5]) == pytest.approx(1e5, rel=1e-9)
    # and for a model with no closed form
    q_c = flow_rate(1e5, NEEDLE, CROSS, [10.0, 0.01, 0.1, 0.8])
    assert pressure_for_flow_rate(q_c, NEEDLE, CROSS, [10.0, 0.01, 0.1, 0.8]) == pytest.approx(1e5, rel=1e-8)


def test_wall_quantities():
    assert wall_stress(1e5, 1e-4, 0.02) == pytest.approx(250.0)
    # HB: gd_w = ((250 - 50)/2)^(1/0.5) = 100^2 = 1e4
    assert wall_shear_rate(1e5, NEEDLE, HB, [50.0, 2.0, 0.5]) == pytest.approx(1e4, rel=1e-9)
    # Newtonian: 4 Q / (pi R^3) = dp R / (2 L mu) = 250
    assert wall_shear_rate(1e5, NEEDLE, NEWT, [1.0]) == pytest.approx(250.0, rel=1e-9)


def test_casson_below_and_above_yield():
    """Renamed 2026-09-05: this was `test_casson_goes_through_the_quadrature_route`, and
    its comment said Casson had "no closed form here". Both became false when
    `flow_rate_casson` was added and `flow_rate` was pointed at it, and the test kept
    passing because its assertions never touched the route. Nothing was silently pinned
    to the old value -- no Casson result was asserted by value anywhere -- but a test
    whose name certifies the opposite of what it checks is worse than no test.
    Route coverage now lives in `tests/test_casson_needle.py`."""
    assert flow_rate(1e4, NEEDLE, CASSON, [50.0, 0.01]) == 0.0   # tau_w = 25 < 50
    q = flow_rate(1e5, NEEDLE, CASSON, [50.0, 0.01])
    assert 0 < q < flow_rate_hb(1e5, 1e-4, 0.02, 0.0, 0.01, 1.0)


def test_non_monotone_flow_curve_is_refused():
    # Cross with m = 3 has a stress maximum; the inversion is not defined past it
    with pytest.raises(ValueError, match="not monotone"):
        flow_rate_gnf(1e5, 1e-4, 0.02, CROSS, [100.0, 0.001, 1.0, 3.0])


def test_bad_inputs_are_rejected():
    with pytest.raises(ValueError):
        Needle(0.0, 0.02)
    with pytest.raises(ValueError):
        Syringe(-1.0)
    with pytest.raises(ValueError):
        flow_rate_hb(-1.0, 1e-4, 0.02, 0.0, 1.0, 1.0)
    with pytest.raises(ValueError):
        pressure_for_flow_rate(0.0, NEEDLE, HB, [1.0, 1.0, 1.0])


# --- force ----------------------------------------------------------------------------

def test_force_is_pressure_times_plunger_area_plus_friction():
    q = 5.5627134e-9  # gives dp = 1e5 for the HB point above
    f = injection_force(q, NEEDLE, SYRINGE, HB, [50.0, 2.0, 0.5], friction=0.3, rho=1000.0)
    assert f.dp == pytest.approx(1e5, rel=1e-8)
    # pi * 2.3e-3^2 = 1.6619025e-5 m^2 -> 1.6619025 N, + 0.3 N
    assert f.force == pytest.approx(1.6619025 + 0.3, rel=1e-6)
    assert f.wall_shear_rate == pytest.approx(1e4, rel=1e-8)
    assert f.plug_fraction == pytest.approx(0.2, rel=1e-8)
    # Re = rho U D / eta_w: U = Q/(pi R^2) = 5.5627134e-9 / 3.1415927e-8 = 0.17706 m/s,
    # eta_w = (50 + 2 * 100) / 1e4 = 0.025 Pa s -> Re = 1000 * 0.17706 * 2e-4 / 0.025 = 1.4165
    assert f.reynolds == pytest.approx(1.4165, rel=1e-3)
    assert np.isnan(injection_force(q, NEEDLE, SYRINGE, HB, [50.0, 2.0, 0.5]).reynolds)


def test_sensitivities_hand_values():
    q = 5.5627134e-9
    s = sensitivities(q, NEEDLE, SYRINGE, HB, [50.0, 2.0, 0.5])
    # dp is homogeneous of degree 1 in (tau_y, k) at fixed Q and n, so the two
    # log-sensitivities sum to exactly 1 (Euler); the split is tau_y's share of the
    # wall stress weighted by the profile, here 0.245 / 0.755
    assert s["tau_y"] + s["k"] == pytest.approx(1.0, abs=1e-6)
    assert 0.2 < s["tau_y"] < 0.3
    # power law at fixed Q: ln dp = const + n ln(Q (3n+1) / (pi n R^3)) so
    # d ln dp / d ln n = n [ln gd_w + 3n/(3n+1) - 1]; gd_w = 15625 -> 0.5 (9.6566 + 0.6 - 1) = 4.628
    q_pl = 9.8174770e-9
    s_pl = sensitivities(q_pl, NEEDLE, SYRINGE, PL, [2.0, 0.5])
    assert s_pl["k"] == pytest.approx(1.0, abs=1e-6)
    assert s_pl["n"] == pytest.approx(0.5 * (np.log(15625.0) + 0.6 - 1.0), rel=1e-4)


# --- the pre-registered UQ gate --------------------------------------------------------

TRUTH = np.array([20.0, 2.0, 0.5])
GRID = np.logspace(-1, 3, 30)


@pytest.fixture(scope="module")
def identified():
    eta = synthetic_rate_curve(HB, TRUTH, GRID, noise_cv=0.03, seed=0)
    return fit_rate_curve(HB, GRID, eta)


def test_gate_uq_contains_truth_and_flags_extrapolation(identified):
    q = 2e-10  # gd_w = 403 1/s at the truth (Newtonian estimate 4Q/(pi R^3) = 255)
    truth_force = injection_force(q, NEEDLE, SYRINGE, HB, TRUTH).force
    uq = force_with_uq(identified, q, NEEDLE, SYRINGE)
    print(uq.summary())
    assert uq.profile_lo <= truth_force <= uq.profile_hi, uq.summary()
    assert uq.sampled_lo <= truth_force <= uq.sampled_hi, uq.summary()
    assert uq.sampled_lo < uq.nominal.force < uq.sampled_hi
    assert uq.profile_lo < uq.nominal.force < uq.profile_hi
    assert uq.n_samples == 2000 and uq.seed == 0
    assert set(uq.sensitivity) == {"tau_y", "k", "n"}
    assert uq.sensitivity["tau_y"] + uq.sensitivity["k"] == pytest.approx(1.0, abs=1e-6)
    # gd_w at this Q is inside the fitted envelope [0.1, 1000]
    assert uq.envelope == pytest.approx((0.1, 1000.0), rel=1e-6)
    assert 0.1 <= uq.nominal.wall_shear_rate <= 1000.0 and not uq.extrapolated
    # a Q large enough to push the wall rate past the envelope is flagged
    big = force_with_uq(identified, 1e-8, NEEDLE, SYRINGE, n_samples=200)
    assert big.nominal.wall_shear_rate > 1000.0 and big.extrapolated


def test_uq_is_reproducible_and_seed_dependent(identified):
    a = force_with_uq(identified, 2e-10, NEEDLE, SYRINGE, n_samples=300)
    b = force_with_uq(identified, 2e-10, NEEDLE, SYRINGE, n_samples=300)
    c = force_with_uq(identified, 2e-10, NEEDLE, SYRINGE, n_samples=300, seed=1)
    assert (a.sampled_lo, a.sampled_hi) == (b.sampled_lo, b.sampled_hi)
    assert (a.sampled_lo, a.sampled_hi) != (c.sampled_lo, c.sampled_hi)
    assert a.profile_lo == b.profile_lo == c.profile_lo  # the profile route has no RNG


def test_uq_without_a_covariance_still_reports_the_profile_interval(identified):
    from dataclasses import replace
    broken = replace(identified.fit, cov=np.full((3, 3), np.nan))
    ident = replace(identified, fit=broken)
    uq = force_with_uq(ident, 2e-10, NEEDLE, SYRINGE, n_samples=50)
    assert np.isnan(uq.sampled_lo) and np.isnan(uq.sampled_hi)
    assert uq.profile_lo < uq.nominal.force < uq.profile_hi
    assert any("singular" in n for n in uq.notes)
