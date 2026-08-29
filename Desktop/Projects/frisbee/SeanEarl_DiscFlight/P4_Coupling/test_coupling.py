"""
Tests for the P4 coupling layer: wind-field interpolation, gust generation,
and the Monte Carlo cell driver.

Run with:  pytest test_coupling.py -v
"""

import sys
import time
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "P1_AnalyticalModel"))
import disc_6dof as d6                              # noqa: E402

from cfd_field import StubField, CFDField           # noqa: E402
from synthetic_turbulence import (                  # noqa: E402
    CorrelatedGusts, WhiteNoiseGusts, wind_with_gusts,
)
from monte_carlo import sweep_one_launch_point      # noqa: E402
import stub_flow_fields as sff                      # noqa: E402


class ConstantTurbField(StubField):
    """Stub with prescribed, uniform k and eps."""
    def k(self, r, t):
        return 3.0

    def eps(self, r, t):
        return 3.0


# ---------------------------------------------------------------------------
# Gust generator
# ---------------------------------------------------------------------------

def test_gusts_are_deterministic_in_time():
    """The gust field must be a function of (r, t), not a fresh random draw.

    This is the property that makes the coupled ODE well-posed. If it fails,
    the adaptive integrator cannot control its step and the landing scatter
    becomes an artefact of the step size rather than of the turbulence.
    """
    f = ConstantTurbField(np.zeros(3))
    g = CorrelatedGusts(f, np.random.default_rng(0), t_max=8.0)
    r = np.array([1.0, 2.0, 3.0])
    for t in (0.0, 0.37, 1.5, 7.9):
        a = g(r, t)
        b = g(r, t)
        np.testing.assert_array_equal(a, b)


def test_isotropic_gusts_have_unit_normalised_variance():
    """sigma must come out as sqrt(2k/3) per component when isotropic."""
    import synthetic_turbulence as st
    f = ConstantTurbField(np.zeros(3))
    g = CorrelatedGusts(f, np.random.default_rng(1), t_max=8.0,
                        model="temporal",
                        sigma_ratios=st.ISOTROPIC_RATIOS)
    assert 0.9 < g._spectrum_check_variance() < 1.1

    r = np.array([0.0, 0.0, 2.0])
    ts = np.linspace(0.0, 200.0, 40000)
    samples = np.array([g(r, t) for t in ts])
    sigma_target = np.sqrt(2.0 * 3.0 / 3.0)
    for comp in range(3):
        assert np.isclose(samples[:, comp].std(), sigma_target, rtol=0.2)


def test_anisotropy_preserves_k_and_splits_the_components():
    """sigma_v/sigma_u = 0.75, sigma_w/sigma_u = 0.5, at UNCHANGED k.

    k is what the RANS supplies. Imposing anisotropy must redistribute the
    energy between components, never change the total, or the gust amplitude
    silently stops matching the flow field it was derived from.
    """
    import synthetic_turbulence as st
    k_val = 3.0
    f = ConstantTurbField(np.array([4.0, 0.0, 0.0]))   # mean wind along +x
    g = CorrelatedGusts(f, np.random.default_rng(5), t_max=8.0,
                        model="temporal", sigma_ratios=st.SIGMA_RATIOS)
    r = np.array([0.0, 0.0, 2.0])
    ts = np.linspace(0.0, 400.0, 60000)
    s = np.array([g(r, t) for t in ts])
    su, sv, sw = s.std(axis=0)          # mean wind is +x, so u=x, v=y, w=z
    assert np.isclose(sv / su, 0.75, rtol=0.05), (su, sv, sw)
    assert np.isclose(sw / su, 0.50, rtol=0.05), (su, sv, sw)
    # k = (su^2 + sv^2 + sw^2) / 2 must be recovered
    assert np.isclose((su**2 + sv**2 + sw**2) / 2.0, k_val, rtol=0.08)


def test_anisotropy_axes_follow_the_mean_wind():
    """The w axis is vertical and the u axis lies along the mean wind."""
    import synthetic_turbulence as st
    r = np.array([0.0, 0.0, 2.0])
    ts = np.linspace(0.0, 400.0, 40000)
    # mean wind along +y: the LARGEST component must now be y, not x
    f = ConstantTurbField(np.array([0.0, 4.0, 0.0]))
    g = CorrelatedGusts(f, np.random.default_rng(6), t_max=8.0,
                        model="temporal", sigma_ratios=st.SIGMA_RATIOS)
    sd = np.array([g(r, t) for t in ts]).std(axis=0)
    assert np.argmax(sd) == 1, sd
    assert np.argmin(sd) == 2, sd     # w (vertical) is always the smallest


def test_temporal_integral_timescale_comes_from_k_over_eps():
    class F(StubField):
        def k(self, r, t):
            return 2.0

        def eps(self, r, t):
            return 0.5
    g = CorrelatedGusts(F(np.zeros(3)), np.random.default_rng(0), t_max=8.0,
                        model="temporal")
    assert np.isclose(g.T_L, 4.0)
    assert np.isnan(g.L_int)


def test_sweeping_integral_length_is_sigma_u_times_k_over_eps():
    """L_int = sqrt(2k/3) * k/eps, the convention stated in the module."""
    class F(StubField):
        def k(self, r, t):
            return 2.0

        def eps(self, r, t):
            return 0.5
    g = CorrelatedGusts(F(np.zeros(3)), np.random.default_rng(0), t_max=8.0,
                        model="sweeping")
    assert np.isclose(g.L_int, np.sqrt(2 * 2.0 / 3.0) * 4.0)
    assert np.isnan(g.T_L)
    assert not g.clipped


def test_gust_autocorrelation_matches_target_scale():
    """Autocorrelation relaxes as exp(-lag / scale) in EITHER coordinate."""
    f = ConstantTurbField(np.zeros(3))
    for model in ("temporal", "sweeping"):
        g = CorrelatedGusts(f, np.random.default_rng(2), t_max=8.0,
                            model=model)
        lags = np.array([0.0, 0.25, 0.5, 1.0, 2.0]) * g.scale
        got = g.autocorrelation(lags)
        want = np.exp(-lags / g.scale)
        assert np.max(np.abs(got - want)) < 0.06, (model, got, want)


def test_sweeping_gust_decorrelates_over_distance_not_time():
    """The defining property of the fix.

    Under `sweeping` the gust at a FIXED point must not change with time (the
    field is frozen), and must decorrelate as the disc moves L_int downrange.
    Under `temporal` it is the other way round -- which, at this site's
    T_L = 10-18 s against a 2.5 s flight, means it barely changes at all.
    """
    f = ConstantTurbField(np.zeros(3))
    r0 = np.array([0.0, 0.0, 2.0])

    sw = CorrelatedGusts(f, np.random.default_rng(3), t_max=8.0,
                         model="sweeping", L_int=6.0, ref_pos=r0)
    # frozen: same place, different time -> identical
    np.testing.assert_allclose(sw(r0, 0.0), sw(r0, 2.4), atol=1e-12)
    # swept: one integral length downrange -> substantially decorrelated
    far = sw(r0 + np.array([6.0, 0.0, 0.0]), 0.0)
    assert np.linalg.norm(far - sw(r0, 0.0)) > 0.1 * np.linalg.norm(sw(r0, 0.0))

    tp = CorrelatedGusts(f, np.random.default_rng(3), t_max=8.0,
                         model="temporal", T_L=14.0, ref_pos=r0)
    # a 2.4 s flight at T_L = 14 s: the "gust" is nearly a constant offset
    rho = float(tp.autocorrelation(np.array([2.4]))[0])
    assert rho > 0.8, rho


def test_correlated_gusts_do_not_slow_the_integrator():
    """Regression guard on the white-noise pathology.

    A non-deterministic RHS made a single throw take >100 s. The smooth gust
    field must keep it in the same order of magnitude as the gust-free case.
    """
    from trajectories import right_handed_backhand
    f = ConstantTurbField(np.zeros(3))
    r0, v0, q0, w0 = right_handed_backhand()

    t0 = time.perf_counter()
    plain = d6.simulate(r0, v0, q0, w0, wind=f.U, rtol=1e-4, atol=1e-6,
                        compute_diag=False)
    dt_plain = time.perf_counter() - t0

    g = CorrelatedGusts(f, np.random.default_rng(0), t_max=8.0)
    t0 = time.perf_counter()
    gusty = d6.simulate(r0, v0, q0, w0, wind=wind_with_gusts(f, g),
                        rtol=1e-4, atol=1e-6, compute_diag=False)
    dt_gusty = time.perf_counter() - t0

    assert dt_gusty < 10.0 * dt_plain + 1.0, (dt_plain, dt_gusty)
    assert len(gusty.t) < 20 * len(plain.t)


def test_white_noise_gusts_are_documented_as_nondeterministic():
    """Pin the deprecated generator's failure mode so it is not reintroduced."""
    f = ConstantTurbField(np.zeros(3))
    g = WhiteNoiseGusts(f, np.random.default_rng(0))
    r = np.array([0.0, 0.0, 2.0])
    a, b = g(r, 1.0), g(r, 1.0)
    assert not np.allclose(a, b), (
        "WhiteNoiseGusts is supposed to be non-deterministic; if this now "
        "passes, the deprecation note in synthetic_turbulence.py is stale"
    )


def test_zero_tke_gives_zero_gust():
    f = StubField(np.zeros(3))          # k = 0
    g = CorrelatedGusts(f, np.random.default_rng(0), t_max=8.0)
    np.testing.assert_allclose(g(np.zeros(3), 1.0), np.zeros(3), atol=1e-15)


# ---------------------------------------------------------------------------
# Field interpolation and the P3 -> P4 schema contract
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def stub_npz(tmp_path_factory):
    path = tmp_path_factory.mktemp("fields") / "stub.npz"
    x = np.arange(-60.0, 60.1, 4.0)
    y = np.arange(-60.0, 60.1, 4.0)
    z = np.array([1.0, 2.0, 3.0, 4.0])
    out = dict(x=x, y=y, z=z, is_stub=True, headings=np.array([90]))
    U, k, eps = sff.build_field(x, y, z, 90)
    out["dir_090__U"] = U
    out["dir_090__k"] = k
    out["dir_090__eps"] = eps
    np.savez_compressed(path, **out)
    return str(path)


def test_cfd_field_reads_the_documented_schema(stub_npz):
    f = CFDField(stub_npz, heading=90)
    U = f.U(np.array([0.0, 0.0, 2.0]), 0.0)
    assert U.shape == (3,) and np.all(np.isfinite(U))
    assert f.k(np.array([0.0, 0.0, 2.0]), 0.0) >= 0.0
    assert f.eps(np.array([0.0, 0.0, 2.0]), 0.0) > 0.0
    assert f.heading == 90


def test_field_interpolation_matches_grid_values_at_nodes(stub_npz):
    """At a grid node the interpolant must return the stored value exactly."""
    d = np.load(stub_npz)
    f = CFDField(stub_npz, heading=90)
    ix, iy, iz = 7, 11, 1
    r = np.array([d["x"][ix], d["y"][iy], d["z"][iz]])
    np.testing.assert_allclose(f.U(r, 0.0), d["dir_090__U"][ix, iy, iz],
                               rtol=1e-9, atol=1e-9)
    assert np.isclose(f.k(r, 0.0), d["dir_090__k"][ix, iy, iz], rtol=1e-9)


def test_field_below_lowest_sample_is_clamped_not_extrapolated(stub_npz):
    """Below the lowest sampled altitude the field must clamp, not run away."""
    d = np.load(stub_npz)
    f = CFDField(stub_npz, heading=90)
    z0 = float(d["z"][0])
    at_floor = f.U(np.array([0.0, 0.0, z0]), 0.0)
    below = f.U(np.array([0.0, 0.0, 0.05]), 0.0)
    np.testing.assert_allclose(below, at_floor, rtol=1e-9)


def test_stub_field_wake_reduces_speed_downstream():
    """The stub's obstacle must actually shelter the region behind it."""
    # Obstacle spans x in [-80, -50], centred at x = -65. With the wind blowing
    # toward +x, x = 0 is 65 m downstream (sheltered) and x = -110 is upstream
    # (undisturbed). y = 120 is far off the wake centreline at the same
    # streamwise station.
    y = np.array([0.0])
    z = np.array([2.0])
    up, _, k_up = None, None, None
    U, k, eps = sff.build_field(np.array([-110.0, 0.0]), y, z, 0)
    speed = np.linalg.norm(U, axis=-1)
    assert speed[1, 0, 0] < speed[0, 0, 0], "no wake deficit downstream"
    assert k[1, 0, 0] > k[0, 0, 0], "wake should raise TKE"

    # Lateral: at the same downstream station, off-centreline is less sheltered.
    U2, _, _ = sff.build_field(np.array([0.0]), np.array([0.0, 120.0]), z, 0)
    s2 = np.linalg.norm(U2, axis=-1)
    assert s2[0, 1, 0] > s2[0, 0, 0], "wake should be confined laterally"


# ---------------------------------------------------------------------------
# Monte Carlo cell
# ---------------------------------------------------------------------------

def test_sweep_cell_is_reproducible():
    f = StubField(np.array([3.0, 0.0, 0.0]))
    kw = dict(launch=(0, 0, 1.5), heading_deg=0,
              throw_type="right_backhand", n_throws=12, field=f)
    a = sweep_one_launch_point(seed=123, **kw)
    b = sweep_one_launch_point(seed=123, **kw)
    np.testing.assert_allclose(a["landings"], b["landings"])
    assert a["success_rate"] == b["success_rate"]


def test_sweep_cell_different_seeds_differ():
    f = StubField(np.array([3.0, 0.0, 0.0]))
    kw = dict(launch=(0, 0, 1.5), heading_deg=0,
              throw_type="right_backhand", n_throws=12, field=f)
    a = sweep_one_launch_point(seed=1, **kw)
    b = sweep_one_launch_point(seed=2, **kw)
    assert not np.allclose(a["landings"], b["landings"])


def test_success_rate_is_a_fraction_and_target_is_downrange():
    f = StubField(np.zeros(3))
    s = sweep_one_launch_point(launch=(10.0, -5.0, 1.5), heading_deg=0,
                               throw_type="right_backhand", bearing_deg=90.0,
                               n_throws=20, field=f, seed=0)
    assert 0.0 <= s["success_rate"] <= 1.0
    # bearing 90 deg -> target is 25 m in +y from the launch point
    np.testing.assert_allclose(s["target_xy"], [10.0, 20.0], atol=1e-9)


def test_headwind_cell_scores_worse_than_still_air():
    """Physical sanity: a stiff headwind must reduce delivery accuracy."""
    still = sweep_one_launch_point(
        launch=(0, 0, 1.5), heading_deg=0, throw_type="right_backhand",
        n_throws=60, field=StubField(np.zeros(3)), seed=5)
    head = sweep_one_launch_point(
        launch=(0, 0, 1.5), heading_deg=0, throw_type="right_backhand",
        n_throws=60, field=StubField(np.array([-6.0, 0.0, 0.0])), seed=5)
    assert head["success_rate"] <= still["success_rate"]
    assert head["mean"][0] < still["mean"][0], "headwind must shorten the throw"


def test_fast_trilinear_matches_scipy(stub_npz):
    """The hand-rolled trilinear must reproduce RegularGridInterpolator.

    Guards the speed optimisation in cfd_field.py: five scipy interpolator
    calls per RHS evaluation dominated the Monte Carlo sweep, so the eight
    corners are now combined by hand with shared cell location and weights.
    """
    from scipy.interpolate import RegularGridInterpolator
    d = np.load(stub_npz)
    f = CFDField(stub_npz, heading=90)
    axes = (d["x"], d["y"], d["z"])
    ref_U = [RegularGridInterpolator(axes, d["dir_090__U"][..., i],
                                     bounds_error=False, fill_value=None)
             for i in range(3)]
    ref_k = RegularGridInterpolator(axes, d["dir_090__k"],
                                    bounds_error=False, fill_value=None)

    rng = np.random.default_rng(3)
    for _ in range(300):
        r = np.array([rng.uniform(d["x"][0], d["x"][-1]),
                      rng.uniform(d["y"][0], d["y"][-1]),
                      rng.uniform(d["z"][0], d["z"][-1])])
        got = f.U(r, 0.0)
        want = np.array([float(ri(r[None, :])[0]) for ri in ref_U])
        np.testing.assert_allclose(got, want, rtol=1e-10, atol=1e-12)
        assert np.isclose(f.k(r, 0.0), float(ref_k(r[None, :])[0]),
                          rtol=1e-10)


def test_out_of_bounds_queries_are_clamped_and_counted(stub_npz):
    f = CFDField(stub_npz, heading=90)
    (x0, x1), (y0, y1), (z0, z1) = f.bounds
    inside = f.U(np.array([x1, y1, z1]), 0.0)
    outside = f.U(np.array([x1 + 500.0, y1 + 500.0, z1 + 50.0]), 0.0)
    np.testing.assert_allclose(inside, outside)
    assert np.all(np.isfinite(outside))
    assert f.fraction_clamped > 0.0


def test_missing_heading_raises_a_helpful_error(stub_npz):
    with pytest.raises(KeyError, match="available directions"):
        CFDField(stub_npz, heading=45)


# ---------------------------------------------------------------------------
# Common random numbers and reference frames  (WP0.2 / WP0.3 / CRN)
# ---------------------------------------------------------------------------

def test_heading_index_is_absent_from_the_crn_spawn_key():
    """Cells differing only in heading must share release draws.

    The between-heading contrast is the comparison this study exists to make.
    With the heading index in the spawn key it was the one comparison getting
    no variance reduction at all.
    """
    from monte_carlo import crn_spawn_key
    base = (2, 1, 0, 1, 3)
    for i_h in range(4):
        idx = (base[0], base[1], i_h, base[3], base[4])
        assert crn_spawn_key(idx) == crn_spawn_key(base)
    # ...but every other axis must still separate cells.
    assert crn_spawn_key((3, 1, 0, 1, 3)) != crn_spawn_key(base)
    assert crn_spawn_key((2, 0, 0, 1, 3)) != crn_spawn_key(base)
    assert crn_spawn_key((2, 1, 0, 0, 3)) != crn_spawn_key(base)
    assert crn_spawn_key((2, 1, 0, 1, 2)) != crn_spawn_key(base)


def test_paired_arms_consume_identical_rng_state():
    """Gust and no-gust arms must advance the RNG identically.

    `CorrelatedGusts.__init__` draws exactly 3 * N_MODES uniforms whatever the
    field or time scale, which is what makes the paired design work. Nothing
    defended that invariant, and it breaks silently if N_MODES ever becomes
    adaptive or field-dependent.
    """
    from synthetic_turbulence import CorrelatedGusts

    def draws_consumed(field, **kw):
        rng = np.random.default_rng(12345)
        before = rng.bit_generator.state
        CorrelatedGusts(field, rng, t_max=8.0, **kw)
        after = rng.bit_generator.state
        # count by replaying: how many uniforms reproduce the same end state?
        probe = np.random.default_rng(12345)
        for n in range(1, 4096):
            probe.uniform()
            if probe.bit_generator.state == after:
                return n
        raise AssertionError("RNG advance exceeded 4096 draws")

    quiet = StubField(np.zeros(3))
    windy = StubField(np.array([8.0, 2.0, 0.0]))
    n_quiet = draws_consumed(quiet)
    n_windy = draws_consumed(windy)
    n_forced = draws_consumed(quiet, T_L=0.05)
    assert n_quiet == n_windy == n_forced, (
        f"gust construction consumed {n_quiet}/{n_windy}/{n_forced} draws; "
        "the paired arms no longer share release conditions")
    assert n_quiet == 3 * CorrelatedGusts.N_MODES


def test_landings_are_actually_on_the_ground():
    """Valid landings must have z ~ 0, not 'finite'.

    `solve_ivp` returns the state at t_max when the ground event never fires,
    and that state is finite, so a finiteness test alone scores a disc still in
    the air as having landed at its current altitude.
    """
    from monte_carlo import GROUND_TOL_M
    s = sweep_one_launch_point(
        launch=(0, 0, 1.5), heading_deg=0, throw_type="right_backhand",
        n_throws=25, field=StubField(np.zeros(3)), seed=7)
    land = s["landings"]
    finite = np.isfinite(land[:, 0])
    assert np.all(np.abs(land[finite, 2]) <= GROUND_TOL_M)
    assert s["n_airborne"] == 0
    assert s["n_valid"] == int(finite.sum())


def test_landing_bias_is_in_the_throw_frame_for_every_bearing():
    """The frame defect regression test.

    A rotationally symmetric field must give a bias that is independent of
    bearing in the THROW frame. Reading `landing_mean` as if it were target-
    relative instead produces an apparent bias of order 2R sin(b/2) -- up to
    50 m at bearing 180 -- which is what was previously mistaken for a defect
    in the simulation.
    """
    f = StubField(np.zeros(3))          # still air: rotationally symmetric
    biases, means = [], []
    for bearing in (0.0, 60.0, 120.0, 180.0, 240.0, 300.0):
        s = sweep_one_launch_point(
            launch=(0, 0, 1.5), heading_deg=0, throw_type="right_backhand",
            bearing_deg=bearing, n_throws=40, field=f, seed=11)
        biases.append(s["landing_bias"])
        means.append(s["mean"][:2])
        # target_xy must be the ROTATED target, in world coordinates
        b = np.radians(bearing)
        np.testing.assert_allclose(
            s["target_xy"], 25.0 * np.array([np.cos(b), np.sin(b)]), atol=1e-9)

    biases = np.array(biases)
    # Same physical throw in every direction => same throw-frame bias.
    assert np.ptp(biases[:, 0]) < 0.5, "downrange bias must not depend on bearing"
    assert np.ptp(biases[:, 1]) < 0.5, "crossrange bias must not depend on bearing"
    # ...and the naive target-relative reading of landing_mean must NOT be
    # bearing-independent, which is exactly why the distinction matters.
    naive = np.linalg.norm(np.array(means) - np.array([25.0, 0.0]), axis=1)
    assert np.ptp(naive) > 20.0


def test_dispersion_and_delivery_use_stated_populations():
    """Blocked throws: excluded from dispersion, counted as failed deliveries."""

    class HalfBlocked(StubField):
        """Everything with y > 0 is inside a building."""
        def is_solid(self, x, y):
            return bool(y > 0.0)

    s = sweep_one_launch_point(
        launch=(0, 0, 1.5), heading_deg=0, throw_type="right_backhand",
        n_throws=40, field=HalfBlocked(np.zeros(3)), seed=3)
    land = s["landings"]
    valid = np.isfinite(land[:, 0]) & (np.abs(land[:, 2]) <= 1e-3)
    blocked = valid & (land[:, 1] > 0.0)
    assert s["n_valid"] == int(valid.sum())
    assert s["n_delivered"] == int((valid & ~blocked).sum())
    if s["n_delivered"] >= 2:
        # the centroid must come from the delivered population only
        np.testing.assert_allclose(
            s["mean"], land[valid & ~blocked].mean(axis=0), atol=1e-12)
    # a blocked throw can never be scored a hit
    assert s["success_rate"] <= s["n_delivered"] / max(s["n_valid"], 1) + 1e-12
