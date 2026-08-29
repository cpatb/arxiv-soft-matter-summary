"""
test_uncertainty.py -- checks for the BPA error-propagation layer.

Run: python -m pytest test_uncertainty.py -q
The posterior tests use small draw counts so the suite stays fast; the production
numbers come from run_uncertainty.py.
"""
import numpy as np
import pytest

import bpa
import uncertainty as unc


@pytest.fixture(scope="module")
def stains():
    return bpa.emit_and_land(np.array([0.0, 0.0, 1.15]), n_drops=14,
                             speed_range=(3.0, 9.0), d0_range=(1.0e-3, 2.0e-3),
                             cone_half_angle=np.radians(35), seed=7)


# --------------------------------------------------------------------------- #
# The error amplifier.
# --------------------------------------------------------------------------- #
def test_amplifier_is_greater_than_one():
    """If it were <= 1 the stain-size route would not amplify error and none of this
    module would be needed.  It is ~3-4."""
    for v in (1.0, 4.0, 9.0, 16.0):
        k = unc.spread_sensitivity(v, 1.5e-3)
        assert 2.5 < k < 4.5, (v, k)


def test_amplifier_matches_a_direct_finite_difference():
    """Cross-check the reported slope against an independent perturbation of the
    forward model, so the sensitivity is not an artefact of the inverse solver."""
    d0, v = 1.5e-3, 4.0
    we, re = bpa._we_re_from_vn(v, d0)
    xi = float(bpa.spread_factor(we, re))
    eps = 0.02
    v_hi = bpa.vn_from_spread(xi * (1 + eps), d0)
    k_direct = (v_hi / v - 1.0) / eps
    k_module = unc.spread_sensitivity(v, d0)
    assert abs(k_direct - k_module) / k_module < 0.10


def test_error_budget_combines_in_quadrature():
    b = unc.speed_error_budget(4.0, 1.5e-3, sigma_W=0.03, sigma_d0=0.04)
    np.testing.assert_allclose(b["sigma_xi"], 0.05, rtol=1e-12)
    np.testing.assert_allclose(b["sigma_v_rel"], b["k"] * 0.05, rtol=1e-12)


def test_zero_speed_sensitivity_raises():
    with pytest.raises(ValueError):
        unc.spread_sensitivity(0.0, 1.5e-3)


# --------------------------------------------------------------------------- #
# The noise model.
# --------------------------------------------------------------------------- #
def test_only_zeroes_every_other_source():
    n = unc.MeasurementNoise().only("rel_d0")
    assert n.rel_d0 > 0.0
    assert n.rel_W == 0.0 and n.rel_L == 0.0 and n.az_deg == 0.0


def test_only_rejects_unknown_field():
    with pytest.raises(ValueError):
        unc.MeasurementNoise().only("rel_pressure")


def test_zero_noise_is_detected():
    assert unc.MeasurementNoise(0.0, 0.0, 0.0, 0.0).is_zero
    assert not unc.MeasurementNoise().is_zero


def test_perturbation_keeps_minor_axis_below_major(stains):
    """Regression: independent noise on W and L can invert the axes, after which
    arcsin(W/L) saturates at 90 degrees and the direction information is silently
    destroyed.  perturb_stains re-orders them."""
    rng = np.random.default_rng(0)
    big = unc.MeasurementNoise(rel_W=0.30, rel_L=0.30, az_deg=0.0, rel_d0=0.0)
    for _ in range(40):
        for s in unc.perturb_stains(stains, big, rng):
            assert s.W <= s.L
            assert 0.0 < s.alpha <= np.pi / 2 + 1e-12


def test_perturbation_preserves_positivity(stains):
    rng = np.random.default_rng(1)
    noisy = unc.perturb_stains(stains, unc.MeasurementNoise(), rng)
    for s in noisy:
        assert s.W > 0 and s.L > 0 and s.d0 > 0


def test_zero_noise_leaves_stains_unchanged(stains):
    rng = np.random.default_rng(2)
    zero = unc.MeasurementNoise(0.0, 0.0, 0.0, 0.0)
    for a, b in zip(stains, unc.perturb_stains(stains, zero, rng)):
        np.testing.assert_allclose([a.W, a.L, a.d0], [b.W, b.L, b.d0], rtol=1e-12)


# --------------------------------------------------------------------------- #
# The posterior.
# --------------------------------------------------------------------------- #
def test_posterior_with_no_noise_collapses_to_a_point(stains):
    zero = unc.MeasurementNoise(0.0, 0.0, 0.0, 0.0)
    post = unc.origin_posterior(stains, zero, n_draws=4, seed=0)
    assert post["cloud"].shape[1] == 3
    assert np.std(post["cloud"][:, 2]) < 1e-6


def test_posterior_spread_grows_with_noise(stains):
    small = unc.MeasurementNoise(0.01, 0.01, 0.5, 0.01)
    large = unc.MeasurementNoise(0.08, 0.08, 4.0, 0.16)
    s_small = np.std(unc.origin_posterior(stains, small, n_draws=12, seed=3)["cloud"][:, 2])
    s_large = np.std(unc.origin_posterior(stains, large, n_draws=12, seed=3)["cloud"][:, 2])
    assert s_large > s_small


def test_posterior_returns_a_usable_credible_interval(stains):
    post = unc.origin_posterior(stains, unc.MeasurementNoise(), n_draws=16, seed=4)
    assert post["n_ok"] >= 8
    assert np.all(post["p2p5"] <= post["p97p5"])
    assert post["cov"].shape == (3, 3)


def test_bias_curve_pushes_the_height_down_as_d0_uncertainty_grows(stains):
    """The central empirical result: more parent-drop uncertainty pushes the
    reconstructed height down.

    Compared at the well-separated endpoints and against the Monte Carlo standard
    error, because adjacent points on this curve are not separable at small draw
    counts.  An earlier version of this test compared 10% against 20% with 10 draws,
    where the standard error is ~10 cm and the true difference is smaller than that;
    it failed for want of power, not because the trend is absent.
    """
    n_draws = 30
    curve = unc.height_bias_curve(stains, 1.15, [0.0, 0.20], n_draws=n_draws, seed=5)
    clean, noisy = curve[0], curve[1]
    sem = noisy["sd"] / np.sqrt(n_draws)
    assert clean["bias"] - noisy["bias"] > 2.0 * sem, (clean, noisy, sem)


def test_zero_uncertainty_point_of_the_curve_is_deterministic(stains):
    curve = unc.height_bias_curve(stains, 1.15, [0.0], n_draws=5, seed=5)
    assert curve[0]["sd"] == 0.0
    np.testing.assert_allclose(curve[0]["lo"], curve[0]["hi"], rtol=1e-12)


def test_attribution_covers_every_source(stains):
    att = unc.attribute_height_bias(stains, 1.15, unc.MeasurementNoise(),
                                    n_draws=6, seed=6)
    assert set(att) == {"rel_W", "rel_L", "az_deg", "rel_d0", "all"}
    for b, s in att.values():
        assert np.isfinite(b) and np.isfinite(s)


def test_posterior_raises_rather_than_returning_an_empty_cloud(stains):
    """Fewer than three stains cannot constrain a 3-D origin, so every draw is skipped.
    A silent empty return would be far worse than an exception here.

    Note what this test does NOT claim: enormous noise does not make the reconstruction
    fail, it makes it wrong.  An earlier version of this test assumed 10^9 relative
    noise would raise; it does not, because perturb_stains re-orders the axes and the
    fit still converges to a meaningless point.  Garbage in, confident garbage out is
    exactly the failure mode this module exists to expose.
    """
    with pytest.raises(RuntimeError):
        unc.origin_posterior(stains[:2], unc.MeasurementNoise(), n_draws=3, seed=7)
