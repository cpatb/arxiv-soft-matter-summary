"""
Unit and integration tests for the 6-DOF disc model.

Run with:  pytest test_disc_6dof.py -v

Three tiers:
  1. Exact mathematical invariants (quaternions, Galilean invariance).
  2. Quantitative agreement with Hummel (2003) published flights -- these are
     the tests that back the numbers quoted in the paper.
  3. Physical-plausibility bounds (range, headwind, spin decay).
"""

import numpy as np
import pytest

import disc_6dof as d6
from coefficients import (
    DEG, DEFAULT_COEFFS, HUMMEL_SHORT, POTTS_CROWTHER, TabulatedCoefficients,
    ULTRASTAR,
)
from trajectories import right_handed_backhand
import validation_hummel as vh


# ---------------------------------------------------------------------------
# Tier 1 -- exact invariants
# ---------------------------------------------------------------------------

def test_quaternion_identity_roundtrip():
    q = np.array([1.0, 0.0, 0.0, 0.0])
    R = d6.q_to_R(q)
    assert np.allclose(R, np.eye(3))


def test_quaternion_90deg_z():
    # 90-degree rotation about z: x_B -> y_N
    q = d6.axis_angle_to_q(np.array([0, 0, 1]), np.pi / 2)
    R = d6.q_to_R(q)
    e_x = np.array([1.0, 0.0, 0.0])
    assert np.allclose(R @ e_x, np.array([0.0, 1.0, 0.0]), atol=1e-12)


def test_rotation_matrix_is_orthonormal():
    """Any quaternion must map to a proper rotation (R^T R = I, det = +1)."""
    rng = np.random.default_rng(20260725)
    for _ in range(25):
        q = d6.q_normalize(rng.normal(size=4))
        R = d6.q_to_R(q)
        assert np.allclose(R.T @ R, np.eye(3), atol=1e-12)
        assert np.isclose(np.linalg.det(R), 1.0, atol=1e-12)


@pytest.mark.parametrize("U", [
    [5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [-3.0, 4.0, 0.0],
])
def test_uniform_wind_is_galilean_equivalent(U):
    """A uniform wind U must be exactly equivalent to shifting v0 by -U.

    This is a physical identity, not an approximation: in a frame translating
    with the air, still-air dynamics apply. It is the strongest available check
    that the wind enters the aerodynamics correctly and only through the
    relative velocity. A sign error or a missed wind term in the moment
    calculation would break it immediately.
    """
    U = np.asarray(U, dtype=float)
    r0, v0, q0, w0 = right_handed_backhand(release_alpha_deg=8.0, roll_deg=0.0)

    with_wind = d6.simulate(r0, v0, q0, w0,
                            wind=d6.uniform_wind(U), t_max=8.0)
    shifted = d6.simulate(r0, v0 - U, q0, w0, wind=d6.still_air, t_max=8.0)

    assert np.isclose(with_wind.flight_time, shifted.flight_time, rtol=1e-6)
    expected = shifted.r[-1] + U * shifted.flight_time
    np.testing.assert_allclose(with_wind.r[-1], expected, atol=1e-6)


def test_still_air_wind_callable_is_a_noop():
    """Passing an explicit zero wind must match the still-air default."""
    r0, v0, q0, w0 = right_handed_backhand()
    a = d6.simulate(r0, v0, q0, w0, wind=d6.still_air, t_max=8.0)
    b = d6.simulate(r0, v0, q0, w0, wind=d6.uniform_wind([0, 0, 0]), t_max=8.0)
    np.testing.assert_allclose(a.r[-1], b.r[-1], atol=1e-9)


# ---------------------------------------------------------------------------
# Tier 2 -- quantitative agreement with Hummel (2003)
# ---------------------------------------------------------------------------

def test_geometry_matches_hummel_table_3_1():
    """Inertia and geometry must be the MEASURED values, not thin-disc ones."""
    g = ULTRASTAR
    assert np.isclose(g.mass, 0.175)
    assert np.isclose(g.diameter, 0.269)
    assert np.isclose(g.area, 0.057)
    assert np.isclose(g.Is, 0.00235)
    assert np.isclose(g.Ip, 0.00122)
    # Guard against a regression to the thin-uniform-disc idealisation, which
    # underestimates both inertias by ~1/3 and so overstates precession.
    assert g.Is > 1.2 * (0.5 * g.mass * g.radius ** 2)
    assert g.Ip > 1.2 * (0.25 * g.mass * g.radius ** 2)


def test_zero_pitching_moment_angle_matches_published_9_deg():
    """Hummel Sec. 2.4.2: the pitching moment vanishes at alpha = 9 deg.

    This pins the single most consequential number in the model -- the trim
    angle of attack sets where turn flips to fade.
    """
    trim = vh.trim_glide(TabulatedCoefficients(), ULTRASTAR)
    assert trim is not None
    assert np.isclose(trim["alpha_trim_deg"], 9.0, atol=0.5), trim


def test_trim_glide_matches_published_glide_angle_and_CL():
    """Hummel Sec. 2.4.2: gamma = -19.3 deg and CL = 0.67 at that trim."""
    trim = vh.trim_glide(TabulatedCoefficients(), ULTRASTAR)
    assert np.isclose(trim["gamma_deg"], vh.REF_GLIDE_ANGLE_DEG, atol=1.0)
    assert np.isclose(trim["CL"], vh.REF_TRIM_CL, rtol=0.05)


def test_hummel_common_throw_reproduces_published_flight():
    """Hummel Sec. 2.4.3 / Fig. 2-9: 25 m, apex 2.55 m, 3.5 s, v_end ~ v0/3.

    These four numbers are the quantitative validation the paper cites.
    """
    res = vh.hummel_common_throw(TabulatedCoefficients(), ULTRASTAR)
    s = vh.summarize(res, TabulatedCoefficients())

    assert np.isclose(s["range_m"], vh.REF_RANGE_M, rtol=vh.TOL_FRAC), s
    assert np.isclose(s["apex_m"], vh.REF_APEX_M, rtol=vh.TOL_FRAC), s
    assert np.isclose(s["time_s"], vh.REF_TIME_S, rtol=vh.TOL_FRAC), s
    v_ref = s["initial_speed"] / 3.0
    assert np.isclose(s["terminal_speed"], v_ref, rtol=vh.TOL_FRAC), s


def test_hummel_common_throw_turns_right_then_fades_left():
    """Fig. 2-9a: curves right, then reverses back left near t ~ 3 s.

    In our z-up frame with the throw along +x, +y is LEFT, so 'right' is a
    decreasing y followed by a minimum before touchdown.
    """
    res = vh.hummel_common_throw(TabulatedCoefficients(), ULTRASTAR)
    s = vh.summarize(res, TabulatedCoefficients())
    assert s["lateral_excursion_m"] > 0.2, "no rightward turn at all"
    assert s["turn_reverses"], "disc never fades back left"
    assert 2.0 < s["t_reversal"] < res.flight_time


# ---------------------------------------------------------------------------
# Tier 3 -- physical plausibility
# ---------------------------------------------------------------------------

def test_still_air_range_reasonable():
    """A canonical flat backhand should travel a plausible distance."""
    r0, v0, q0, w0 = right_handed_backhand(release_alpha_deg=8.0, roll_deg=0.0)
    res = d6.simulate(r0, v0, q0, w0, wind=d6.still_air, t_max=8.0)
    assert 20.0 < res.range_m < 80.0, f"range = {res.range_m:.2f} m"


def test_headwind_reduces_range():
    """Range under a 5 m/s headwind must be shorter than in still air."""
    r0, v0, q0, w0 = right_handed_backhand(release_alpha_deg=8.0, roll_deg=0.0)
    R_still = d6.simulate(r0, v0, q0, w0, wind=d6.still_air, t_max=8.0).range_m
    R_head = d6.simulate(r0, v0, q0, w0,
                         wind=d6.uniform_wind([-5.0, 0.0, 0.0]),
                         t_max=8.0).range_m
    assert R_head < 0.95 * R_still


def test_crosswind_drift_is_downwind_and_bounded():
    """Crosswind drift must be downwind, and strictly below the U_c*T bound.

    NOTE ON THE BOUND. A common coaching heuristic says the disc drifts by
    about U_c * T, i.e. that it is fully advected by the crosswind. That is the
    behaviour of a passive particle, and this model does NOT reproduce it: it
    predicts roughly a quarter of U_c * T for a 5 m/s crosswind.

    That is a consequence of the model's axisymmetry, not a bug -- the Galilean
    test above shows the wind coupling is exact. An axisymmetric disc has no
    fin and therefore no yaw stiffness, and at the high angles of attack
    reached late in a flight the aerodynamic resultant is dominated by the
    disc-normal component, which carries almost no lateral force. Launch
    sideslip is therefore shed only slowly. Adding the rim side-force
    coefficient CS (P1 README item 4) would increase the drift; until then,
    crosswind drift from this model should be read as a LOWER bound.
    """
    r0, v0, q0, w0 = right_handed_backhand(release_alpha_deg=8.0, roll_deg=0.0)
    Uc = 5.0
    res = d6.simulate(r0, v0, q0, w0,
                      wind=d6.uniform_wind([0.0, Uc, 0.0]), t_max=8.0)
    drift = float(res.r[-1, 1] - r0[1])
    assert drift > 0.0, f"drift {drift:.2f} m is not downwind"
    assert drift < Uc * res.flight_time, "drift exceeds the full-advection bound"


def test_crosswind_drift_increases_with_wind_speed():
    """Monotonicity: a stronger crosswind must push the disc further downwind."""
    r0, v0, q0, w0 = right_handed_backhand(release_alpha_deg=8.0, roll_deg=0.0)
    drifts = []
    for Uc in (0.0, 2.5, 5.0, 7.5):
        res = d6.simulate(r0, v0, q0, w0,
                          wind=d6.uniform_wind([0.0, Uc, 0.0]), t_max=8.0)
        drifts.append(float(res.r[-1, 1] - r0[1]))
    assert np.all(np.diff(drifts) > 0), drifts


def test_spin_decay_slow():
    """Over a full flight, spin should decay by only a few percent."""
    r0, v0, q0, w0 = right_handed_backhand(release_alpha_deg=8.0, roll_deg=0.0)
    res = d6.simulate(r0, v0, q0, w0, wind=d6.still_air, t_max=8.0)
    ratio = abs(res.omega_B[-1, 2]) / abs(w0[2])
    assert 0.7 < ratio < 1.05, (
        f"|Omega_final/Omega_initial| = {ratio:.3f}, expected near 1.0"
    )


def test_lift_and_drag_signs():
    """Drag must oppose the relative wind; lift must have an upward component
    for a level disc at positive angle of attack."""
    R = np.eye(3)                      # disc level, body frame = inertial
    v_N = np.array([14.0, 0.0, -1.0])  # forward and slightly descending
    omega = np.array([0.0, 0.0, -50.0])
    F_N, M_B, diag = d6.aero_force_moment(
        v_N, R, omega, np.zeros(3), DEFAULT_COEFFS, ULTRASTAR)
    assert diag["alpha_deg"] > 0.0
    assert F_N[0] < 0.0, "drag must decelerate a disc moving in +x"
    assert F_N[2] > 0.0, "lift must push up at positive alpha"


def test_zero_relative_velocity_is_finite():
    """Edge case: no relative wind must give zero force, not a NaN."""
    F_N, M_B, diag = d6.aero_force_moment(
        np.zeros(3), np.eye(3), np.array([0.0, 0.0, -50.0]),
        np.zeros(3), DEFAULT_COEFFS, ULTRASTAR)
    assert np.all(np.isfinite(F_N)) and np.allclose(F_N, 0.0)
    assert np.all(np.isfinite(M_B))


def test_coefficient_sets_expose_a_consistent_interface():
    """Every published set must be usable interchangeably by the integrator."""
    for cf in (HUMMEL_SHORT, POTTS_CROWTHER, TabulatedCoefficients()):
        assert np.isfinite(cf.CL(5 * DEG))
        assert cf.CD(5 * DEG) > 0.0
        assert np.isfinite(cf.CM(5 * DEG, 0.0))
        assert np.isfinite(cf.CR(0.0, -50.0, ULTRASTAR.diameter))
        assert np.isfinite(cf.spin_down_moment(-50.0))


def test_drag_is_minimised_at_the_zero_lift_angle():
    """CD is a parabola with its vertex at alpha0 = -CL0/CLa (Hummel eq. 3.4)."""
    cf = HUMMEL_SHORT
    a0 = cf.alpha0
    assert cf.CD(a0) < cf.CD(a0 + 5 * DEG)
    assert cf.CD(a0) < cf.CD(a0 - 5 * DEG)
    assert np.isclose(cf.CL(a0), 0.0, atol=1e-12), "alpha0 must be zero-lift"
