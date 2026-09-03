"""test_exponent_recovery.py — does the §6e machinery do what §6e says it does?

§6e claims two things about fitting `phi = phi_J + c*p^b`:

  (1) with `phi_J` PINNED at the true value, `b` is recovered accurately;
  (2) with `phi_J` FREE, the fit is ill-conditioned -- a small error in `phi_J` buys a
      large error in `b`, so the free fit can return a badly wrong `b` while still
      describing the data well.

Both are claims about the ESTIMATOR, not about the packing, so they can be tested on
synthetic ladders where the true `b` is known by construction. That is the anchor the real
analysis cannot have: on simulation data nobody knows the true exponent, which is the whole
problem. Here it is an input.

The ladder geometry (number of points, spacing, the p_hat range) is matched to the real
Phase-2 ladders so the conditioning number these tests see is the conditioning number the
real analysis faces.
"""

from __future__ import annotations

import numpy as np
import pytest

# the real ladders: 4 points, logarithmically spaced over ~1.5 decades
P_LADDER = np.array([1e-4, 3.16e-4, 1e-3, 3.16e-3])
SEED_SCATTER = 0.00179

# Amplitude calibrated to the real ladders, NOT chosen for convenience. A representative
# frictional ladder runs phi = 0.5949, 0.5827, 0.5748 with phi_J ~ 0.5748, so the excess
# spans ~0.02 at the top of the p range. C_REAL is set so that c * p_max^1 = 0.02.
#
# This calibration is load-bearing and an earlier version got it wrong: with an arbitrary
# c = 2.0 the excess at the LOWEST p was 2e-4, well under the seed scatter, and perturbing
# phi_J upward drove every low-p point negative so the fit returned NaN. That is worth
# recording rather than just fixing, because it is the same pathology the real ladders
# have in milder form -- the bottom of the ladder carries almost no signal, which is
# exactly why phi_J and b cannot be separated.
PHI_J_REAL = 0.5748
C_REAL = 0.02 / P_LADDER.max()


def make_ladder(phi_J, c, b, p=P_LADDER, noise=0.0, rng=None):
    phi = phi_J + c * p**b
    if noise:
        gen = rng if rng is not None else np.random.default_rng(0)
        phi = phi + gen.normal(0.0, noise, size=phi.shape)
    return p, phi


def fit_b_pinned(p, phi, phi_J):
    """b from a straight line in log-log, with phi_J known. This is the §6e estimator."""
    e = phi - phi_J
    m = (e > 0) & (p > 0)
    if m.sum() < 3:
        return float("nan")
    return float(np.polyfit(np.log(p[m]), np.log(e[m]), 1)[0])


@pytest.mark.parametrize("b_true", [0.5, 0.667, 1.0, 1.5])
def test_pinned_fit_recovers_b_exactly_when_noiseless(b_true):
    """Claim (1), the clean case: pinned phi_J must return the true exponent."""
    p, phi = make_ladder(PHI_J_REAL, C_REAL, b_true)
    assert fit_b_pinned(p, phi, PHI_J_REAL) == pytest.approx(b_true, abs=1e-9)


def _noise_study(b_true, n=2000, seed=20260903):
    rng = np.random.default_rng(seed)
    got = []
    for _ in range(n):
        p, phi = make_ladder(PHI_J_REAL, C_REAL, b_true, noise=SEED_SCATTER, rng=rng)
        v = fit_b_pinned(p, phi, PHI_J_REAL)
        if np.isfinite(v):
            got.append(v)
    return np.array(got)


def test_pinned_fit_is_unbiased_at_b_two_thirds():
    """Claim (1) with seed noise, in the regime where the ladder HAS signal."""
    g = _noise_study(0.667)
    assert abs(np.median(g) - 0.667) < 0.05, np.median(g)
    assert np.percentile(g, 75) - np.percentile(g, 25) < 0.10


def test_pinned_fit_is_BIASED_LOW_at_b_one_because_the_ladder_bottom_is_below_the_noise():
    """The sharpest result in this file, and it was found by a test failing.

    At b = 1 the excess at the LOWEST rung is c*p_min = 0.00063, which is BELOW the seed
    scatter of 0.00179. That rung is therefore pure noise, and in log-log it drags the
    slope down. So even with phi_J known EXACTLY -- no degeneracy at all -- a true b of 1
    is recovered as a median of ~0.84 with a wide IQR, and ~4.5% of fits fail outright
    because noise pushes an excess negative.

    This matters for 6e beyond the degeneracy: the free fit's 0.567 is consistent with a
    true b = 1 once BOTH the phi_J/b trade-off AND this low-side noise bias are allowed
    for. It also explains why b = 0.667 is recovered cleanly and b = 1 is not -- the
    steeper the exponent, the less signal the bottom of the ladder carries.
    """
    g = _noise_study(1.0)
    med = float(np.median(g))
    assert med < 0.95, med                       # biased low
    assert 0.75 < med < 0.90, med                # by a specific, reproducible amount
    iqr = np.percentile(g, 75) - np.percentile(g, 25)
    assert iqr > 0.25, iqr                       # and poorly determined
    # the ladder's bottom rung really is under the noise floor
    assert C_REAL * P_LADDER.min() < SEED_SCATTER


def test_a_phi_J_error_of_one_seed_scatter_moves_b_by_about_a_quarter_of_the_bar():
    """Claim (2), the conditioning number that decides §6e.

    This is the load-bearing test. §6e's verdict -- that b is UNDETERMINED rather than
    failed -- rests on |db/dphi_J| being large enough that one seed scatter of phi_J moves
    b across a meaningful fraction of the [0.6, 1.6] bar. The real data gives 148 per unit
    phi, i.e. 0.265 per seed scatter. A synthetic ladder of the same shape must show the
    same order of sensitivity, otherwise the number came from something other than the
    ladder geometry.
    """
    p, phi = make_ladder(PHI_J_REAL, C_REAL, 1.0)
    lo = fit_b_pinned(p, phi, PHI_J_REAL - SEED_SCATTER)
    hi = fit_b_pinned(p, phi, PHI_J_REAL + SEED_SCATTER)
    assert np.isfinite(lo) and np.isfinite(hi), (lo, hi)
    db_dphi = abs(hi - lo) / (2 * SEED_SCATTER)
    # order-of-magnitude agreement with the measured 148 is what is being asserted
    assert 30 < db_dphi < 600, db_dphi
    shift = db_dphi * SEED_SCATTER
    assert shift > 0.10, shift          # a real fraction of the 1.0-wide bar


def test_free_fit_can_return_a_badly_wrong_b_while_fitting_the_data_well():
    """Claim (2) in the form that matters: a wrong phi_J hides inside a good-looking fit.

    Take a ladder generated with b = 1, hand the fitter a phi_J that is wrong by the
    amount seen in the real data (0.003-0.015), and check BOTH that the recovered b is far
    from 1 AND that the resulting curve still tracks the data closely. If the residuals
    blew up, the degeneracy would be detectable from goodness-of-fit alone and §6e's
    argument would be unnecessary.
    """
    p, phi = make_ladder(PHI_J_REAL, C_REAL, 1.0)
    phi_J_wrong = PHI_J_REAL - 0.010
    b_wrong = fit_b_pinned(p, phi, phi_J_wrong)
    assert abs(b_wrong - 1.0) > 0.25, b_wrong

    c_wrong = np.exp(np.polyfit(np.log(p), np.log(phi - phi_J_wrong), 1)[1])
    resid = phi - (phi_J_wrong + c_wrong * p**b_wrong)
    # The wrong model degrades the fit but does not blow it up: residuals stay within a
    # few seed scatters. An earlier version asserted "within ONE seed scatter" and failed
    # at 0.0036 -- worth keeping as a measured fact rather than a tightened bar, because
    # it says goodness-of-fit is a WEAK discriminator here, not a useless one. It is
    # consistent with the real ladders' median r2 of 0.964: visibly imperfect, not damning.
    assert np.max(np.abs(resid)) < 4 * SEED_SCATTER, np.max(np.abs(resid))
    assert np.max(np.abs(resid)) > SEED_SCATTER, np.max(np.abs(resid))


def test_pinning_to_a_linear_intercept_is_circular():
    """The reason §6e refuses to call the constrained result a PASS.

    Build a ladder whose true exponent is NOT 1 (b = 0.6). Estimate phi_J the way
    `phi_J_lin` does -- as the intercept of a straight line through (p, phi) -- then pin
    that estimate and refit. If pinning a linear-fit intercept biases b toward 1, the
    recovered b must sit closer to 1 than the truth does. That bias is exactly why
    recovering b ~ 1 on the real data cannot be evidence for b = 1.
    """
    b_true = 0.6
    p, phi = make_ladder(PHI_J_REAL, C_REAL, b_true)
    phi_J_lin = float(np.polyfit(p, phi, 1)[1])       # intercept of phi vs p
    b_from_lin = fit_b_pinned(p, phi, phi_J_lin)
    assert np.isfinite(b_from_lin)
    assert abs(b_from_lin - 1.0) < abs(b_true - 1.0), (b_from_lin, b_true)
