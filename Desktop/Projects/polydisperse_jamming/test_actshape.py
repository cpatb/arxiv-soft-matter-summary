"""Tests for actshape.py. Anchors are closed forms and known distributions, never the
module's own output."""

import math

import numpy as np
import pytest

import actshape as A

RNG = np.random.default_rng(20260922)


def test_weibull_constant_closed_form():
    # ln(ln 0.1 / ln 0.9) and ln(ln 0.25 / ln 0.75), by hand
    assert math.isclose(A.weibull_const(0.1, 0.9), math.log(2.302585093 / 0.105360516),
                        rel_tol=1e-8)
    assert math.isclose(A.weibull_const(0.25, 0.75), math.log(1.386294361 / 0.287682072),
                        rel_tol=1e-8)


@pytest.mark.parametrize("delta", [0.7, 1.0, 1.4, 2.5])
def test_delta_recovered_from_weibull_sample(delta):
    u = RNG.weibull(delta, 400_000) * 3.7          # scale must not matter
    for lo, hi in A.Q_PAIRS:
        assert abs(A.delta_from_quantiles(u, lo, hi) / delta - 1) < 0.02


def test_exponential_is_delta_one_and_normal_is_not():
    u = RNG.exponential(1.0, 400_000)
    for lo, hi in A.Q_PAIRS:
        assert abs(A.delta_from_quantiles(u, lo, hi) - 1.0) < 0.02
    # discriminating negative: a narrow distribution peaked at its mean is NOT Weibull-1,
    # and its two quantile pairs must disagree in the direction H5 predicts
    v = np.abs(RNG.normal(1.0, 0.3, 400_000))
    assert A.delta_from_quantiles(v, 0.25, 0.75) > A.delta_from_quantiles(v, 0.1, 0.9)


def test_derjaguin_is_blind_at_equal_radii_and_not_otherwise():
    fn = RNG.exponential(1.0, 5000)
    r = np.full(5000, 0.5)
    ds = [A.delta_from_quantiles(A.g_values(fn, r, r, 1.0, n), 0.1, 0.9)
          for n in ("constant", "derjaguin", "area")]
    assert max(ds) - min(ds) < 1e-12
    r1 = RNG.lognormal(0.0, 0.5, 5000)
    r2 = RNG.lognormal(0.0, 0.5, 5000)
    dc = A.delta_from_quantiles(A.g_values(fn, r1, r2, 1.0, "constant"), 0.1, 0.9)
    dd = A.delta_from_quantiles(A.g_values(fn, r1, r2, 1.0, "derjaguin"), 0.1, 0.9)
    assert abs(dc - dd) > 1e-3


def test_reduced_radius():
    assert math.isclose(float(A.reduced_radius(np.array([1.0]), np.array([1.0]))[0]), 0.5)
    assert math.isclose(float(A.reduced_radius(np.array([1.0]), np.array([3.0]))[0]), 0.75)


def test_f_crit_edges():
    assert A.f_crit(0.9, 0.9) == 1.0            # lower window edge: everything must activate
    assert A.f_crit(0.95, 0.9) == pytest.approx(0.5)
    assert A.f_crit(1.0, 0.9) == 0.0
    assert A.f_crit(0.85, 0.9) == math.inf
    with pytest.raises(ValueError):
        A.f_crit(0.9, 1.0)


def test_sigma_jam_matches_weibull_closed_form():
    # for a Weibull g with scale lam and shape d, Q_{1-fc} = lam*(ln 1/fc)^(1/d), so
    # sigma_jam/kappa = 1/Q and sigma_jam/sigma* = (ln 1/fc)^(-1/d)
    d, lam = 1.4, 2.0
    g = RNG.weibull(d, 1_000_000) * lam
    for fc in (0.1, 0.25, 0.5, 0.75, 0.9):
        want = 1.0 / (lam * math.log(1 / fc) ** (1 / d))
        assert A.sigma_jam_over_kappa(g, fc) == pytest.approx(want, rel=0.01)
        ratio = A.sigma_jam_over_kappa(g, fc) / A.sigma_star_over_kappa(g)
        assert ratio == pytest.approx(A.weibull_sigma_jam_ratio(fc, d), rel=0.01)
    # at fc = 1/e the ratio is exactly 1 for ANY delta: that is the sigma* convention
    for dd in (0.5, 1.0, 3.0):
        assert A.weibull_sigma_jam_ratio(A.INV_E, dd) == pytest.approx(1.0, abs=1e-12)


def test_activation_curve_is_survival():
    g = np.array([1.0, 2.0, 3.0, 4.0])
    # s = sigma/kappa; contact frictional iff g > 1/s
    assert A.activation_curve(g, np.array([1 / 2.5]))[0] == pytest.approx(0.5)
    assert A.activation_curve(g, np.array([1 / 0.5]))[0] == pytest.approx(1.0)
    assert A.activation_curve(g, np.array([1 / 10.0]))[0] == pytest.approx(0.0)


def _relax_reaches(stress_at, x, chi, gd, fc):
    """Brute force: df/dt = F(f) - f from 0, F(f) = fraction active at stress eta(f)*gd."""
    f = 0.0
    for _ in range(40_000):
        s = A.eta_rel(f, x, chi) * gd
        # invert stress_at on a grid: largest f' with stress_at(f') <= s
        grid = np.linspace(1e-6, 1 - 1e-6, 4001)
        ok = stress_at(grid) <= s
        F = grid[ok].max() if ok.any() else 0.0
        f += 0.05 * (F - f)
    return f >= fc


@pytest.mark.parametrize("delta", [0.7, 2.6])
def test_rate_threshold_matches_brute_relaxation(delta):
    x, chi = 0.85, 0.8
    fc = A.f_crit(x, chi)

    def weib(f):
        return np.log(1 / np.asarray(f)) ** (-1 / delta)

    gdc, fstar = A.rate_threshold(weib, x, chi)
    assert 0 < fstar < fc
    assert not _relax_reaches(weib, x, chi, 0.97 * gdc, fc)
    assert _relax_reaches(weib, x, chi, 1.03 * gdc, fc)


def test_rate_threshold_refuses_outside_window():
    with pytest.raises(ValueError):
        A.rate_threshold(lambda f: f, 0.7, 0.8)
