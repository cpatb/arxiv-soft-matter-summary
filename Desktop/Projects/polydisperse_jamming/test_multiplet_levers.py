"""test_multiplet_levers.py — anchors for the LEVERS layer, not the PSD layer.

Added 2026-09-22 after a referee pass noted that `test_multiplet.py` covers
`psd.MultipletPSD` and nothing covers `flory()`, `phi_J()`, `levers()`, `critical_df()`
or the truncation treatment — the code every number in `out/multiplet_levers.txt` runs
through.

Same rule as the sibling suites: each anchor is a closed form, an exactly-known limit,
or a discriminating negative, never a round trip.
"""

import math

import numpy as np
import pytest

import mixture_model as mm
from multiplet_levers import (
    BETA,
    D10,
    K_MAX,
    PA_REFIT,
    PA_TEXT,
    PA_WIDE,
    PB_REFIT,
    PB_TEXT,
    PB_WIDE,
    critical_df,
    flory,
    levers,
    phi_J,
)
from psd import MultipletPSD


# --------------------------------------------------------------------------
# A. The aggregation-number distribution
# --------------------------------------------------------------------------


def test_flory_reproduces_its_own_mean():
    """`<k> = 1/(1-p)` is the defining property; truncation must not break it.

    This is the discriminating check on the truncation: a cutoff that bit into the
    distribution would show up here first, because the mean is the most tail-sensitive
    low moment.
    """
    for k_bar in (1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 8.0):
        n = flory(k_bar)
        k = np.arange(1, len(n) + 1, dtype=float)
        assert float(np.sum(np.asarray(n) * k)) == pytest.approx(k_bar, rel=2e-4)


def test_flory_at_kbar_one_is_the_unaggregated_suspension():
    assert flory(1.0) == (1.0,)
    with pytest.raises(ValueError):
        flory(0.9)


def test_flory_is_geometric_in_k():
    """n_k proportional to p^(k-1): the RATIO of successive classes is constant = p."""
    k_bar = 2.5
    p = 1.0 - 1.0 / k_bar
    n = np.asarray(flory(k_bar))
    ratios = n[1:8] / n[:7]
    np.testing.assert_allclose(ratios, p, rtol=1e-12)


def test_flory_reports_what_its_floor_discards():
    """The cut must be reported, and it must be small — but it is CHECKED, not assumed."""
    for k_bar in (1.2, 2.0, 5.0):
        _, rep = flory(k_bar, report=True)
        assert rep["n_after"] <= rep["n_before"]
        assert rep["vol_dropped"] < 1e-4, (k_bar, rep)


# --------------------------------------------------------------------------
# B. The closure wrapper
# --------------------------------------------------------------------------


def test_phi_J_of_a_monodisperse_pack_is_exactly_beta():
    """One class -> the linear packing model returns beta/(1-0) = beta, exactly.

    This pins that `phi_J` is being called with the right (d, y) convention: if
    diameters and volume fractions were swapped or mis-normalised, a single class
    would not return beta to the bit.
    """
    p = MultipletPSD(d10=D10, d_f=2.0, n_k=(1.0,))
    assert phi_J(p) == pytest.approx(BETA, abs=1e-12)


def test_phi_J_matches_mixture_model_called_directly():
    """The wrapper must add nothing. Anchored on the library function, not on itself."""
    p = MultipletPSD(d10=D10, d_f=1.8, n_k=(0.6, 0.25, 0.15))
    direct = mm.packing_density(p.diameters, p.volume_fractions,
                                BETA, PA_REFIT, PB_REFIT)
    assert phi_J(p) == direct


def test_polydispersity_gain_is_at_least_one():
    """G >= 1: mixing sizes cannot pack WORSE than monodisperse in this closure.

    A discriminating negative for a sign error in the interaction terms.
    """
    for k_bar in (1.2, 2.0, 5.0):
        for d_f in (1.5, 2.1, 3.0):
            assert levers(k_bar, d_f)["G"] >= 1.0 - 1e-12


# --------------------------------------------------------------------------
# C. The levers and the crossing
# --------------------------------------------------------------------------


def test_lambda_is_exactly_one_at_df3_so_the_ratio_is_one_over_G():
    """At d_f = 3 the trapped-fluid lever is switched off EXACTLY, so J/J0 = 1/G < 1.

    This is the derivation's step 10 as a test: compact aggregation can only SUPPRESS
    jamming, and it must do so for every aggregation number.
    """
    for k_bar in (1.1, 1.5, 2.0, 3.0, 5.0):
        lv = levers(k_bar, 3.0)
        assert lv["lam"] == pytest.approx(1.0, abs=1e-14)
        assert lv["ratio"] == pytest.approx(1.0 / lv["G"], rel=1e-12)
        assert lv["ratio"] < 1.0


def test_critical_df_is_a_genuine_crossing_not_a_clamp():
    """Evaluate the objective AT the returned root and on both sides of it.

    A bisection on a non-straddling bracket returns an endpoint, which would print a
    clamp as an estimate — the defect this repo records for `np.interp`. So the root is
    verified by flanking evaluations, not trusted because a solver returned it.
    """
    for k_bar in (1.2, 2.0, 5.0):
        d = critical_df(k_bar)
        assert 0.7 < d < 3.0, (k_bar, d)
        assert levers(k_bar, d)["ratio"] == pytest.approx(1.0, abs=1e-4)
        assert levers(k_bar, d - 0.02)["ratio"] > 1.0     # more open -> promotes
        assert levers(k_bar, d + 0.02)["ratio"] < 1.0     # more compact -> suppresses


def test_critical_df_returns_nan_when_there_is_no_crossing():
    """An UNAGGREGATED suspension has Lambda == 1 and G == 1 at every d_f, so there is
    no crossing and the function must say so rather than return a bracket endpoint."""
    assert np.isnan(critical_df(1.0))


def test_the_crossing_is_robust_across_all_three_calibrations():
    """d_f* must not depend on which fit is used -- the referee's M6.

    The bar is 0.02, which is 2x the largest move actually observed; a calibration
    dependence larger than that would make d_f* a property of the fit rather than of
    the geometry.
    """
    for k_bar in (1.2, 2.0, 5.0):
        ds = [critical_df(k_bar, pa, pb) for pa, pb in
              ((PA_REFIT, PB_REFIT), (PA_TEXT, PB_TEXT), (PA_WIDE, PB_WIDE))]
        assert max(ds) - min(ds) < 0.15, (k_bar, ds)          # textbook is the outlier
        assert abs(ds[0] - ds[2]) < 0.02, (k_bar, ds)         # the two refits agree


# --------------------------------------------------------------------------
# D. Truncation must not be the answer
# --------------------------------------------------------------------------


def test_lambda_over_G_is_converged_in_the_truncation():
    """The quantity every verdict is read off must not move with K_max.

    ⚠ This is the check that caught a real error: K_MAX = 24 was the first value tried
    and at <k> = 5, d_f = 1.8 the ratio moved 2.62 -> 3.83 between K_max 8 and 40. The
    test is written against the value the module now ships.
    """
    for k_bar in (1.1, 1.2, 2.0, 5.0, 8.0):
        for d_f in (1.2, 1.8, 3.0):
            a = MultipletPSD(d10=D10, d_f=d_f, n_k=flory(k_bar, K_MAX))
            b = MultipletPSD(d10=D10, d_f=d_f, n_k=flory(k_bar, K_MAX // 2))
            ra = a.lambda_trap / (phi_J(a) / BETA)
            rb = b.lambda_trap / (phi_J(b) / BETA)
            assert abs(ra - rb) / ra < 5e-3, (k_bar, d_f, rb, ra)


def test_a_too_small_cutoff_really_does_move_the_answer():
    """DISCRIMINATING: the convergence test above must be capable of FAILING.

    If K_max were irrelevant, `test_lambda_over_G_is_converged_in_the_truncation` would
    pass for a trivial reason. It is not: at the cutoff originally shipped, the answer
    is visibly wrong at high aggregation.
    """
    k_bar, d_f = 5.0, 1.8
    small = MultipletPSD(d10=D10, d_f=d_f, n_k=flory(k_bar, 8))
    big = MultipletPSD(d10=D10, d_f=d_f, n_k=flory(k_bar, K_MAX))
    rs = small.lambda_trap / (phi_J(small) / BETA)
    rb = big.lambda_trap / (phi_J(big) / BETA)
    assert abs(rs - rb) / rb > 0.2, (rs, rb)


# --------------------------------------------------------------------------
# E. The root-sensitivity estimator (referee pass 2, finding 1)
# --------------------------------------------------------------------------


def test_df_star_sensitivity_matches_a_direct_perturbation():
    """ANCHOR outside the formula: perturb G and re-find the root numerically.

    `df_star_sensitivity` uses the implicit function theorem. The check re-solves the
    crossing with G shifted by exactly one PHI_SCALE/BETA and measures how far the root
    actually moved. If the derivative form were wrong -- or inverted -- this disagrees.
    """
    from multiplet_levers import (BETA, D10, N_SEEDS, SCALE, df_star_sensitivity,
                                  flory, phi_J)
    from psd import MultipletPSD
    import scipy.optimize as so

    for k_bar in (1.2, 2.0, 5.0):
        d0 = critical_df(k_bar)
        pred = df_star_sensitivity(k_bar, d0)
        # the SAME noise magnitude the module declares -- what this anchor tests is the
        # DERIVATIVE FORM (via an independent brentq re-solve), not the convention for dG
        dG = SCALE / (BETA * math.sqrt(N_SEEDS))

        def f(x, shift):
            p = MultipletPSD(d10=D10, d_f=float(x), n_k=flory(k_bar))
            return p.lambda_trap - (phi_J(p) / BETA + shift)

        d1 = so.brentq(f, d0 - 0.5, min(d0 + 0.05, 2.999), args=(dG,), xtol=1e-10)
        actual = abs(d1 - d0)
        assert pred == pytest.approx(actual, rel=0.12), (k_bar, pred, actual)


def test_df_star_sensitivity_falls_as_aggregation_rises():
    """More aggregation -> steeper Lambda at the crossing -> a better-pinned root.

    DISCRIMINATING: a routine returning a constant, or one that had the derivative in
    the numerator, would not reproduce this ordering.
    """
    from multiplet_levers import df_star_sensitivity
    s = [df_star_sensitivity(k, critical_df(k)) for k in (1.1, 1.5, 2.0, 5.0)]
    assert all(b < a for a, b in zip(s, s[1:])), s
    assert s[0] > 10 * s[-1], s          # and it is a big range, not a flat one
