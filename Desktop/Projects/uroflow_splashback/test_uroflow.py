"""
test_uroflow.py -- validation of the uroflow splashback model against analytic limits,
published constants, and conservation laws.

v1 of this project had no tests.  Every check below either reproduces a closed-form
limit, reproduces a number from a cited paper, or asserts a conservation law.  Two of
them (test_wall_normal_is_horizontal, test_pool_normal_is_vertical) are regression tests
for the geometry inversion that v1 shipped.

Run: python -m pytest test_uroflow.py -q
"""
import numpy as np
import pytest

import uroflow as u


# --------------------------------------------------------------------------- #
# 1. Jet stability: analytic limits.
# --------------------------------------------------------------------------- #
def test_inviscid_limit_recovers_long_wave_rayleigh():
    """With mu -> 0 the dispersion relation must give omega_max = sqrt(sigma/(8 rho r^3))
    at x = 1/sqrt(2).  That is the exact maximum of the LONG-WAVE relation implemented
    here (coefficient 0.35355)."""
    d = 4e-3
    r = d / 2.0
    omega, x = u.max_growth_rate(d, mu=1e-12)
    expected = np.sqrt(u.SIGMA / (8.0 * u.RHO * r**3))
    np.testing.assert_allclose(omega, expected, rtol=1e-4)
    np.testing.assert_allclose(x, 1.0 / np.sqrt(2.0), rtol=1e-3)


def test_long_wave_constant_is_within_4pc_of_exact_rayleigh():
    """The exact Rayleigh (Bessel) result is 0.3433 sqrt(sigma/(rho r^3)).  The long-wave
    approximation gives 0.35355, i.e. 3.0% high.  Assert the gap is where we think it is
    rather than silently inheriting a 3% bias."""
    d = 4e-3
    r = d / 2.0
    omega, _ = u.max_growth_rate(d, mu=1e-12)
    coefficient = omega / np.sqrt(u.SIGMA / (u.RHO * r**3))
    assert 0.343 < coefficient < 0.343 * 1.04


def test_viscosity_slows_growth():
    """Weber's viscous term can only reduce the growth rate."""
    d = 4e-3
    inviscid, _ = u.max_growth_rate(d, mu=1e-12)
    viscous, _ = u.max_growth_rate(d, mu=1.0)      # 1000x water
    assert viscous < inviscid
    assert viscous > 0.0


def test_growth_rate_scales_as_r_minus_three_halves():
    """omega ~ sqrt(sigma/(rho r^3)), so halving the diameter multiplies omega by 2^1.5."""
    o1, _ = u.max_growth_rate(4e-3, mu=1e-12)
    o2, _ = u.max_growth_rate(2e-3, mu=1e-12)
    np.testing.assert_allclose(o2 / o1, 2.0**1.5, rtol=1e-3)


def test_rayleigh_drop_diameter_is_1p89_jet_diameters():
    """One wavelength of jet collapses into one drop: D/d = (3 pi / x_opt)^(1/3) with
    x_opt = 1/sqrt(2) gives 1.89."""
    d = 4e-3
    D = u.rayleigh_drop_diameter(d, x_opt=1.0 / np.sqrt(2.0))
    np.testing.assert_allclose(D / d, 1.89, rtol=0.01)


# --------------------------------------------------------------------------- #
# 2. Falling-jet kinematics: conservation laws.
# --------------------------------------------------------------------------- #
def test_speed_after_drop_is_torricelli():
    v0 = 1.67
    fall = 0.75
    np.testing.assert_allclose(u.jet_speed_after_drop(v0, fall),
                               np.sqrt(v0**2 + 2 * u.G * fall), rtol=1e-12)


def test_jet_radius_conserves_volume_flux():
    """Q = pi r^2 v must be invariant along the jet."""
    r0, v0 = 2e-3, 1.67
    Q0 = np.pi * r0**2 * v0
    for fall in [0.0, 0.1, 0.5, 1.0]:
        v = u.jet_speed_after_drop(v0, fall)
        r = u.jet_radius(r0, v0, v)
        np.testing.assert_allclose(np.pi * r**2 * v, Q0, rtol=1e-12)


def test_falling_jet_thins():
    r0, v0 = 2e-3, 1.67
    v = u.jet_speed_after_drop(v0, 0.75)
    assert u.jet_radius(r0, v0, v) < r0


def test_disturbance_integral_reduces_to_textbook_parallel_jet():
    """With no vertical drop the accumulated-growth criterion must reproduce
    L_b = v/omega_max * Lambda exactly."""
    Q, d, lam = 21e-6, 4e-3, 12.0
    v0 = u.exit_velocity(Q, d)
    r0 = d / 2.0
    s_break = u.breakup_arclength(v0, r0, lambda s: 0.0, lam, s_max=5.0, n=8000)
    np.testing.assert_allclose(s_break, u.breakup_length_parallel(Q, d, lam), rtol=2e-3)


def test_breakup_returns_inf_when_jet_survives():
    v0, r0 = 1.67, 2e-3
    assert np.isinf(u.breakup_arclength(v0, r0, lambda s: 0.0, lam=50.0, s_max=0.05))


# --------------------------------------------------------------------------- #
# 3. Target geometry.  These are the regression tests for the v1 inversion.
# --------------------------------------------------------------------------- #
def test_wall_normal_is_horizontal():
    """A horizontally-travelling stream hits a VERTICAL wall head-on: all normal, no
    tangential.  v1 scored exactly this case as the driest."""
    wall = u.Target("wall", plane=0.3)
    v_n, v_t = u.resolve_impact(vx=2.0, vz=0.0, target=wall)
    np.testing.assert_allclose(v_n, 2.0, rtol=1e-12)
    np.testing.assert_allclose(v_t, 0.0, atol=1e-12)


def test_pool_normal_is_vertical():
    """A vertically-falling stream hits a POOL head-on."""
    pool = u.Target("pool", plane=0.0)
    v_n, v_t = u.resolve_impact(vx=0.0, vz=-3.0, target=pool)
    np.testing.assert_allclose(v_n, 3.0, rtol=1e-12)
    np.testing.assert_allclose(v_t, 0.0, atol=1e-12)


def test_grazing_impact_has_small_normal_component():
    wall = u.Target("wall", plane=0.3)
    v_n, v_t = u.resolve_impact(vx=0.1, vz=-3.0, target=wall)
    assert v_n < 0.2
    assert v_t > 2.9


def test_impact_speed_is_preserved_by_the_decomposition():
    for tgt in (u.Target("wall", 0.3), u.Target("pool", 0.0)):
        v_n, v_t = u.resolve_impact(1.3, -2.7, tgt)
        np.testing.assert_allclose(np.hypot(v_n, v_t), np.hypot(1.3, 2.7), rtol=1e-12)


def test_bad_target_kind_raises():
    with pytest.raises(ValueError):
        u.Target("ceiling", plane=1.0)


def test_negative_film_raises():
    with pytest.raises(ValueError):
        u.Target("pool", plane=0.0, film=-0.01)


# --------------------------------------------------------------------------- #
# 4. Trajectory.
# --------------------------------------------------------------------------- #
def test_ballistic_conserves_energy():
    """speed_impact^2 = v0^2 + 2 g * fall, independent of launch angle."""
    v0, h = 1.67, 0.75
    for theta in np.radians([0.0, -15.0, -40.0]):
        traj = u.ballistic_to_target(v0, theta, h, u.Target("pool", 0.0))
        np.testing.assert_allclose(
            traj["speed_impact"]**2, v0**2 + 2 * u.G * traj["fall_at_impact"], rtol=1e-6)


def test_arclength_is_at_least_the_straight_line():
    traj = u.ballistic_to_target(1.67, np.radians(-20), 0.75, u.Target("wall", 0.3))
    straight = np.hypot(traj["x"][-1], traj["z"][-1] - 0.75)
    assert traj["s"] >= straight * (1.0 - 1e-9)


def test_unreachable_pool_raises():
    """A stream launched steeply upward from below a high plane never gets there."""
    with pytest.raises(u.StreamMissesTarget):
        u.ballistic_to_target(0.5, np.radians(80), 0.0, u.Target("pool", plane=5.0))


def test_steep_aim_that_lands_short_is_rejected():
    """Regression: v2.0 of ballistic_to_target had no floor check, so a stream aimed
    steeply downward was reported as reaching a distant wall from underground."""
    with pytest.raises(u.StreamMissesTarget):
        u.ballistic_to_target(1.67, np.radians(-85), 0.75, u.Target("wall", plane=0.60))


def test_shallow_aim_from_the_same_height_does_reach_the_wall():
    """The floor guard must not reject legitimate configurations."""
    traj = u.ballistic_to_target(1.67, np.radians(-10), 0.75, u.Target("wall", plane=0.30))
    assert traj["s"] > 0.30
    assert traj["z"][-1] > 0.0


def test_missing_target_is_a_valueerror_subclass():
    """Callers written against ValueError must keep working."""
    assert issubclass(u.StreamMissesTarget, ValueError)


# --------------------------------------------------------------------------- #
# 5. Drag.
# --------------------------------------------------------------------------- #
def test_drag_deficit_is_small_for_a_large_drop_short_flight():
    """The v1 assertion 'drag is negligible' is now a computed number.  A ~5 mm drop
    over ~0.3 m should lose only a few per cent of its speed."""
    pool = u.Target("pool", plane=0.0)
    out = u.drop_flight_with_drag(vx=1.6, vz=-2.0, x0=0.0, z0=0.3, D=5e-3, target=pool)
    assert 0.0 < out["drag_speed_deficit"] < 0.05


def test_drag_deficit_grows_as_drops_shrink():
    pool = u.Target("pool", plane=0.0)
    big = u.drop_flight_with_drag(1.6, -2.0, 0.0, 0.6, 5e-3, pool)["drag_speed_deficit"]
    small = u.drop_flight_with_drag(1.6, -2.0, 0.0, 0.6, 0.5e-3, pool)["drag_speed_deficit"]
    assert small > big


# --------------------------------------------------------------------------- #
# 6. Splash criteria against published numbers.
# --------------------------------------------------------------------------- #
def test_dry_threshold_brackets_the_published_constant():
    """Find the speed at which K crosses 57.7 and check the criterion flips there."""
    D = 3e-3
    lo = u.splash_dry(0.3, D)
    hi = u.splash_dry(6.0, D)
    assert not lo["splashes"] and lo["K"] < u.K_DRY_CRITICAL
    assert hi["splashes"] and hi["K"] > u.K_DRY_CRITICAL


def test_dry_K_has_the_right_scaling():
    """K = We^(1/2) Re^(1/4).  We ~ v^2 so We^(1/2) ~ v; Re ~ v so Re^(1/4) ~ v^(1/4).
    K therefore scales as v^(5/4), and doubling v multiplies K by 2^1.25 = 2.378."""
    D = 3e-3
    k1 = u.splash_dry(2.0, D)["K"]
    k2 = u.splash_dry(4.0, D)["K"]
    np.testing.assert_allclose(k2 / k1, 2.0**1.25, rtol=1e-9)


def test_cossali_K_has_the_right_scaling():
    """K = We^0.8 Re^0.4 ~ v^(1.6+0.4) = v^2, so doubling v quadruples K."""
    D = 3e-3
    k1 = u.splash_film(2.0, D, film=0.5 * D)["K"]
    k2 = u.splash_film(4.0, D, film=0.5 * D)["K"]
    np.testing.assert_allclose(k2 / k1, 4.0, rtol=1e-9)


def test_cossali_threshold_matches_the_published_correlation():
    """K_c = 2100 + 5880 delta^1.44.  Check the two ends of the correlated range."""
    D = 3e-3
    r_thin = u.splash_film(3.0, D, film=0.1 * D)
    r_thick = u.splash_film(3.0, D, film=1.0 * D)
    np.testing.assert_allclose(r_thin["K_c"], 2100.0 + 5880.0 * 0.1**1.44, rtol=1e-9)
    np.testing.assert_allclose(r_thick["K_c"], 2100.0 + 5880.0, rtol=1e-9)


def test_cossali_flags_extrapolation_outside_the_correlated_range():
    D = 3e-3
    assert u.splash_film(3.0, D, film=0.01 * D)["extrapolated"]
    assert u.splash_film(3.0, D, film=3.0 * D)["extrapolated"]
    assert not u.splash_film(3.0, D, film=0.5 * D)["extrapolated"]


def test_deeper_film_raises_the_splash_threshold():
    """Cossali's K_c increases with delta, so a deeper pool is harder to splash at fixed
    impact conditions.  This is the standing-water lever that v1 did not have at all."""
    D = 3e-3
    shallow = u.splash_film(3.0, D, film=0.15 * D)
    deep = u.splash_film(3.0, D, film=0.9 * D)
    assert deep["K_c"] > shallow["K_c"]
    assert deep["excess"] < shallow["excess"]


def test_splash_excess_is_monotone_in_impact_speed():
    D = 3e-3
    e = [u.splash_film(v, D, film=0.5 * D)["excess"] for v in [1.0, 2.0, 3.0, 4.0]]
    assert all(b > a for a, b in zip(e, e[1:]))


def test_film_criterion_rejects_a_dry_surface():
    with pytest.raises(ValueError):
        u.splash_film(3.0, 3e-3, film=0.0)


# --------------------------------------------------------------------------- #
# 7. Watson reuse and the coherent-column branch.
# --------------------------------------------------------------------------- #
def test_watson_radius_matches_the_bidet_model_formula():
    a, v = 1.5e-3, 3.0
    nu = u.MU / u.RHO
    Re_a = v * a / nu
    expected = a * (Re_a / (4.0 * u.C_DELTA**2))**(1.0 / 3.0)
    np.testing.assert_allclose(u.watson_transition_radius(a, v), expected, rtol=1e-12)


def test_coherent_film_conserves_flux():
    """h(r) = Q/(2 pi r v) must give back Q when integrated round the annulus."""
    Q, r, v = 21e-6, 0.01, 3.0
    h = u.coherent_film_thickness(Q, r, v)
    np.testing.assert_allclose(2.0 * np.pi * r * h * v, Q, rtol=1e-12)


# --------------------------------------------------------------------------- #
# 8. Input validation.
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("bad", [0.0, -1.0, np.nan, np.inf])
def test_non_physical_diameter_raises(bad):
    with pytest.raises(ValueError):
        u.exit_velocity(21e-6, bad)


def test_zero_flow_raises():
    with pytest.raises(ValueError):
        u.exit_velocity(0.0, 4e-3)


# --------------------------------------------------------------------------- #
# 9. End-to-end sanity on the composite evaluation.
# --------------------------------------------------------------------------- #
def test_evaluate_runs_both_regimes_and_labels_them():
    Q, d = 21e-6, 4e-3
    pool = u.Target("pool", plane=0.0, film=1.5e-3)
    near = u.evaluate(Q, d, np.radians(-10), h=0.15, target=pool, lam=6.0)
    far = u.evaluate(Q, d, np.radians(-10), h=0.90, target=pool, lam=6.0)
    assert near["regime"] in ("coherent", "fragmented")
    assert far["regime"] == "fragmented"
    assert far["path"] > near["path"]


def test_impact_diameter_shrinks_with_fall_height():
    """[FIX 2] regression: the jet must thin, so the local diameter at impact falls."""
    Q, d = 21e-6, 4e-3
    pool = u.Target("pool", plane=0.0, film=1.5e-3)
    lo = u.evaluate(Q, d, np.radians(-10), h=0.15, target=pool, lam=6.0)
    hi = u.evaluate(Q, d, np.radians(-10), h=0.90, target=pool, lam=6.0)
    assert hi["d_at_impact"] < lo["d_at_impact"] < d


def test_splash_excess_is_monotone_in_aim_angle_for_a_wall():
    """The critical-aim-angle bisection in run_uroflow.py assumes monotonicity.  Assert
    it here rather than in a comment, over the reachable range."""
    Q, d, lam = 21e-6, 4e-3, 3.1
    wall = u.Target("wall", plane=0.30, film=0.5e-3)
    vals = []
    for th_deg in np.linspace(0.0, -55.0, 12):
        try:
            r = u.evaluate(Q, d, np.radians(th_deg), 0.75, wall, lam)
        except u.StreamMissesTarget:
            break
        if np.isfinite(r["excess"]):
            vals.append(r["excess"])
    assert len(vals) >= 6
    assert all(b <= a + 1e-9 for a, b in zip(vals, vals[1:])), vals


def test_steeper_arrival_lowers_the_normal_component_at_a_wall():
    """The mechanism behind the corrected advice, asserted directly."""
    Q, d, lam = 21e-6, 4e-3, 3.1
    wall = u.Target("wall", plane=0.30, film=0.5e-3)
    shallow = u.evaluate(Q, d, np.radians(0.0), 0.75, wall, lam)
    steep = u.evaluate(Q, d, np.radians(-45.0), 0.75, wall, lam)
    assert steep["v_normal"] < shallow["v_normal"]
    assert steep["incidence_deg"] > shallow["incidence_deg"]


def test_calibrated_lambda_reproduces_the_observed_breakup_band():
    Q, d = 21e-6, 4e-3
    lam_lo, lam_hi, dr_lo, dr_hi = u.calibrate_lambda(Q, d, launch_height=0.75)
    assert 0.0 < lam_lo <= lam_hi
    # A perturbed meatus jet must be noisier than a quiet laboratory jet.
    assert lam_hi < u.LAMBDA_IDEAL
    assert 0.0 < dr_lo <= dr_hi < 1.0
