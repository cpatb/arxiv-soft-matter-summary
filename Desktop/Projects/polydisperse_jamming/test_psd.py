"""test_psd.py — anchors for the size-distribution layer.

Engine-independent: runs on Windows without Yade.

Every test here is anchored on something that does NOT pass through the code under
test: a closed-form moment identity derived by hand, an independent Monte-Carlo
estimate, an exactly-known special case, or a discriminating negative. There are no
round-trip tests -- `inverse(forward(x)) == x` would certify the arithmetic and say
nothing about whether the moments are right.
"""

import numpy as np
import pytest

from psd import (
    LogNormalPSD,
    carman_kozeny_kappa,
    moment_ratio_d32_d10,
    moment_ratio_d43_d10,
)

# --------------------------------------------------------------------------
# A. Closed-form identities, hand-derived
# --------------------------------------------------------------------------


def test_monodisperse_is_exact_not_a_limit():
    """CV = 0 must give every moment ratio exactly 1.0, not 1+eps."""
    assert moment_ratio_d32_d10(0.0) == 1.0
    assert moment_ratio_d43_d10(0.0) == 1.0
    p = LogNormalPSD(d10=1e-6, cv=0.0)
    assert p.sigma_l == 0.0
    assert p.d32 == 1e-6
    assert p.d43 == 1e-6
    assert p.d50 == pytest.approx(1e-6, rel=1e-15)


def test_d32_ratio_matches_hand_derivation():
    """d32/d10 = (1+CV^2)^2, derived by hand:

        <d^n> = exp(n mu + n^2 s^2 / 2)          [lognormal]
        d32   = <d^3>/<d^2> = exp(mu + 5 s^2/2)
        d10   = <d>         = exp(mu +   s^2/2)
        ratio = exp(2 s^2),  and  s^2 = ln(1+CV^2)
              = (1+CV^2)^2
    """
    for cv in (0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5):
        s2 = np.log1p(cv**2)
        assert moment_ratio_d32_d10(cv) == pytest.approx(np.exp(2 * s2), rel=1e-14)
        assert moment_ratio_d43_d10(cv) == pytest.approx(np.exp(3 * s2), rel=1e-14)


def test_specific_hand_computed_points():
    """Two points computed by hand, independently of the implementation.

    CV = 0.5:  (1 + 0.25)^2 = 1.25^2 = 1.5625          exactly
    CV = 1.0:  (1 + 1.00)^2 = 2.00^2 = 4.0             exactly
    kappa ratio at CV = 0.5: 1.5625^2 = 2.44140625     exactly
    """
    assert moment_ratio_d32_d10(0.5) == pytest.approx(1.5625, rel=1e-15)
    assert moment_ratio_d32_d10(1.0) == pytest.approx(4.0, rel=1e-15)
    p0 = LogNormalPSD(d10=1e-6, cv=0.0)
    p5 = LogNormalPSD(d10=1e-6, cv=0.5)
    k0 = carman_kozeny_kappa(p0.d32, 0.44)
    k5 = carman_kozeny_kappa(p5.d32, 0.44)
    assert k5 / k0 == pytest.approx(2.44140625, rel=1e-12)


def test_analytic_moment_agrees_with_direct_integration():
    """<d^n> from the closed form vs numerical quadrature of the lognormal pdf.

    Integrated in u = ln(d), where the pdf is an ordinary Gaussian:

        <d^n> = INT e^(n u) * N(u; mu, s) du

    In the linear variable the integrand is a spike of width ~1e-6 sitting inside a
    range of 1e-3, and `quad`'s adaptive subdivision steps straight over it -- the
    n = 2 moment came back 3.1x low. That was a QUADRATURE failure, not a physics
    one, and the fix is a better-conditioned integral rather than a looser bound.

    Measured agreement in the log variable, n = 1..4:
        relerr = 4.9e-15, 4.7e-09, 6.6e-09, 9.3e-09
    so the bar is set at 1e-7 -- ten times looser than the worst observed point, and
    ten times TIGHTER than the 1e-6 the linear-variable version failed at.

    `quad`'s own `abserr` is NOT used as a gate: it is wildly conservative for an
    integrand spanning this many decades (up to 32% relative at n = 4) and asserting
    on it fails while the value is correct to nine digits. A self-reported error
    estimate is not a measurement of the error.

    This is still an independent path: the closed form exp(n mu + n^2 s^2/2) is the
    analytic evaluation of exactly this integral, and the quadrature never sees it.
    """
    from scipy.integrate import quad

    p = LogNormalPSD(d10=2e-6, cv=0.4)
    mu, s = p.mu_l, p.sigma_l

    def integrand(u, n):
        return np.exp(n * u) * np.exp(-((u - mu) ** 2) / (2 * s**2)) / (
            s * np.sqrt(2 * np.pi)
        )

    for n in (1, 2, 3, 4):
        num, _abserr = quad(integrand, mu - 14 * s, mu + 14 * s, args=(n,), limit=500)
        assert num == pytest.approx(p.moment(n), rel=1e-7), f"moment {n}"


# --------------------------------------------------------------------------
# B. Independent Monte-Carlo estimate
# --------------------------------------------------------------------------


def test_monte_carlo_reproduces_d32_ratio():
    """A large untruncated MC sample must reproduce the closed-form d32/d10.

    Tolerance widens with CV because d32 is a third/second moment ratio and the
    lognormal tail converges slowly -- at CV = 1 a 2e7-draw sample lands near 3.96
    against the exact 4.00. That is a property of the estimator, not an error, so
    the test asserts a CV-dependent tolerance and DOCUMENTS the drift rather than
    hiding it behind a loose global bound.
    """
    rng = np.random.default_rng(20260831)
    for cv, tol in ((0.1, 2e-3), (0.3, 8e-3), (0.5, 3e-2), (0.8, 1.2e-1)):
        p = LogNormalPSD(d10=1.0, cv=cv)
        d = rng.lognormal(p.mu_l, p.sigma_l, 4_000_000)
        d32_mc = np.mean(d**3) / np.mean(d**2)
        d10_mc = np.mean(d)
        assert d32_mc / d10_mc == pytest.approx(moment_ratio_d32_d10(cv), rel=tol), cv


# --------------------------------------------------------------------------
# C. Truncation must be reported, never assumed away
# --------------------------------------------------------------------------


def test_truncation_lowers_d32_and_realised_moments_say_so():
    """The whole point of `realised_moments`: after truncation the ANALYTIC d32 is
    wrong, and quoting it would be the exact error this module exists to prevent.

    Truncation removes the upper tail, and d32 is dominated by the upper tail, so the
    realised d32 must be strictly BELOW the analytic one -- a directional prediction,
    not just 'different'.
    """
    rng = np.random.default_rng(7)
    p = LogNormalPSD(d10=1.0, cv=0.8)
    d = p.sample(200_000, rng, truncate=(0.2, 5.0))
    m = p.realised_moments(d)
    assert m["d32"] < p.d32, "truncated d32 must fall below the analytic value"
    assert m["d_max"] <= 5.0
    assert m["d_min"] >= 0.2
    # and the untruncated draw must recover it
    d_full = p.sample(200_000, rng, truncate=None)
    m_full = p.realised_moments(d_full)
    assert m_full["d32"] == pytest.approx(p.d32, rel=0.10)


def test_realised_cv_tracks_requested_cv_when_untruncated():
    rng = np.random.default_rng(11)
    for cv in (0.1, 0.3, 0.5):
        p = LogNormalPSD(d10=1.0, cv=cv)
        d = p.sample(400_000, rng, truncate=None)
        assert p.realised_moments(d)["cv"] == pytest.approx(cv, rel=0.05)


# --------------------------------------------------------------------------
# D. Discriminating negatives
# --------------------------------------------------------------------------


def test_d32_and_d43_are_ordered_and_distinct():
    """A discriminating negative against the commonest confusion in this project:
    d10 < d32 < d43 strictly, for any CV > 0. Code that returned 'the diameter' for
    all three would pass any single-moment test and fail this one.
    """
    for cv in (0.05, 0.2, 0.5, 1.0):
        p = LogNormalPSD(d10=1e-6, cv=cv)
        assert p.d10 < p.d32 < p.d43, f"ordering broken at cv={cv}"
        assert p.d50 < p.d10, "median must sit below the number mean for a lognormal"


def test_carman_kozeny_is_not_merely_proportional_to_phi():
    """Discriminating negative on the permeability: kappa must carry the FULL
    (1-phi)^3/phi^2 porosity factor, not just d^2. Doubling phi from 0.30 to 0.60
    changes the porosity factor by a hand-computable amount:

        f(phi) = (1-phi)^3 / phi^2
        f(0.30) = 0.343 / 0.09   = 3.811111...
        f(0.60) = 0.064 / 0.36   = 0.177777...
        ratio   = 21.4375
    """
    d32 = 1e-6
    k30 = carman_kozeny_kappa(d32, 0.30)
    k60 = carman_kozeny_kappa(d32, 0.60)
    assert k30 / k60 == pytest.approx(21.4375, rel=1e-10)


def test_rejects_bad_input():
    with pytest.raises(ValueError):
        LogNormalPSD(d10=-1.0, cv=0.2)
    with pytest.raises(ValueError):
        LogNormalPSD(d10=1e-6, cv=-0.1)
    with pytest.raises(ValueError):
        moment_ratio_d32_d10(-0.5)
    with pytest.raises(ValueError):
        carman_kozeny_kappa(1e-6, 0.0)
    with pytest.raises(ValueError):
        carman_kozeny_kappa(1e-6, 1.0)
    with pytest.raises(ValueError):
        carman_kozeny_kappa(-1e-6, 0.5)


def test_radii_conversion_is_exactly_half_the_diameters():
    """The single diameter->radius conversion in the project. A factor-of-2 error
    here rescales every packing fraction by 8x and is the units-auditor's canonical
    finding, so it gets its own test.
    """
    rng = np.random.default_rng(3)
    p = LogNormalPSD(d10=2e-6, cv=0.0)
    r = p.radii_for_yade(100, rng)
    assert np.allclose(r, 1e-6)
    p2 = LogNormalPSD(d10=2e-6, cv=0.4)
    rng2a = np.random.default_rng(99)
    rng2b = np.random.default_rng(99)
    d = p2.sample(1000, rng2a)
    r2 = p2.radii_for_yade(1000, rng2b)
    assert np.allclose(r2, d / 2.0)


# --------------------------------------------------------------------------
# E. The M1 table in the spec must regenerate from this module
# --------------------------------------------------------------------------


def test_spec_m1_table_regenerates():
    """SPEC_2026-08-31.md §1.2 quotes d32/d10 at six CVs. If this module and the spec
    ever disagree, one of them is wrong; this test makes that a failure rather than a
    silent drift.
    """
    expected = {0.00: 1.0000, 0.10: 1.0201, 0.20: 1.0816,
                0.30: 1.1881, 0.50: 1.5625, 0.80: 2.6896}
    for cv, ratio in expected.items():
        assert moment_ratio_d32_d10(cv) == pytest.approx(ratio, abs=5e-5), cv
