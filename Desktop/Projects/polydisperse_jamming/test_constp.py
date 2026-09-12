"""test_constp.py — the constant-pressure servo and the dilatancy fit, without Yade.

Yade does not import on the Windows side where this suite runs, so the two pieces that
can silently invert a result -- the servo SIGN and the `alpha` extraction -- were written
as pure functions specifically so they could be tested here.

Anchors outside the code under test: the sign convention is fixed by Yade's documented
`velGrad` (diagonal positive in EXPANSION) together with `shear._stress_hat` (p positive
in COMPRESSION), and `alpha` is recovered from a synthetic exponential whose rate is
known exactly.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import shear as S

GDOT = 1.0e-4
P_STAR = 3.0e-4


# =====================================================================================
# the servo sign, which is the thing that runs away if it is backwards
# =====================================================================================
def test_over_pressured_cell_EXPANDS():
    """p > p* must give edot_v > 0. Yade expands on a positive velGrad diagonal."""
    assert S.servo_edot(2 * P_STAR, P_STAR, GDOT) > 0


def test_under_pressured_cell_CONTRACTS():
    assert S.servo_edot(0.5 * P_STAR, P_STAR, GDOT) < 0


def test_on_target_asks_for_nothing():
    assert S.servo_edot(P_STAR, P_STAR, GDOT) == 0.0


def test_the_gain_is_dimensionless_and_strain_clocked():
    """edot_v scales with gdot, so the gain cannot carry a hidden time unit."""
    a = S.servo_edot(2 * P_STAR, P_STAR, GDOT)
    b = S.servo_edot(2 * P_STAR, P_STAR, 10 * GDOT)
    assert b == pytest.approx(10 * a, rel=1e-12)


def test_the_response_is_linear_in_the_relative_error_below_the_clamp():
    small = S.servo_edot(1.01 * P_STAR, P_STAR, GDOT, gain=0.5)
    assert small == pytest.approx(0.5 * GDOT * 0.01, rel=1e-12)


def test_the_clamp_binds_and_binds_both_ways():
    lim = S.SERVO_EDOT_CLAMP * GDOT
    assert S.servo_edot(1e6 * P_STAR, P_STAR, GDOT) == pytest.approx(lim)
    assert S.servo_edot(0.0, P_STAR, GDOT) == pytest.approx(-lim)


def test_an_unjammed_cell_asks_for_maximum_contraction():
    """p_hat <= 0 is an unjammed pack; contracting is the right response, and a
    NEGATIVE p must not be allowed to produce a larger-than-clamp demand."""
    assert S.servo_edot(-5.0, P_STAR, GDOT) == pytest.approx(
        -S.SERVO_EDOT_CLAMP * GDOT)


def test_a_nonfinite_pressure_freezes_the_servo_rather_than_guessing():
    assert S.servo_edot(float("nan"), P_STAR, GDOT) == 0.0


def test_a_nonpositive_set_point_is_refused():
    for bad in (0.0, -1e-4):
        with pytest.raises(ValueError):
            S.servo_edot(P_STAR, bad, GDOT)


# =====================================================================================
# alpha: recovered from a synthetic decay whose rate is known exactly
# =====================================================================================
def _decay(rate: float, phi_c: float, phi0: float, n: int = 200, gmax: float = 3.0):
    """phi(gamma) = phi_c + (phi0 - phi_c) exp(rate * gamma), rate < 0."""
    g = np.linspace(0.0, gmax, n)
    return g, phi_c + (phi0 - phi_c) * np.exp(rate * g)


@pytest.mark.parametrize("phi0", [0.55, 0.64])
def test_recovers_a_known_dilatancy_rate_from_either_side(phi0):
    """alpha = -rate/phi_c, and it must come out the same whether the pack started
    above phi_c (dilating) or below it (contracting)."""
    rate, phi_c = -2.5, 0.59
    g, phi = _decay(rate, phi_c, phi0)
    out = S.dilatancy_rate(g, phi, phi_c)
    assert out["alpha_rate"] == pytest.approx(rate, rel=1e-8)
    assert out["alpha"] == pytest.approx(-rate / phi_c, rel=1e-8)
    assert out["alpha_sign_ok"]
    assert out["alpha_note"] == ""


def test_a_different_rate_gives_a_different_alpha():
    """The discriminating negative case: the fit must not return a constant."""
    phi_c = 0.59
    a = S.dilatancy_rate(*_decay(-2.5, phi_c, 0.64), phi_c)["alpha"]
    b = S.dilatancy_rate(*_decay(-5.0, phi_c, 0.64), phi_c)["alpha"]
    assert b == pytest.approx(2 * a, rel=1e-6)


def test_refuses_when_the_signal_never_moved():
    """A rate fitted to a flat signal is the 2026-09-01 lesson: gate on the SHAPE the
    functional form requires, not only on signal-to-noise."""
    rate, phi_c = -0.01, 0.59      # over gamma <= 1.2 this moves 0.005 decades
    g, phi = _decay(rate, phi_c, 0.64)
    out = S.dilatancy_rate(g, phi, phi_c)
    assert not math.isfinite(out["alpha"])
    assert "decades" in out["alpha_note"]
    assert out["alpha_decades"] < S.ALPHA_MIN_DECADES


def test_refuses_and_says_so_when_the_state_runs_AWAY():
    """A positive rate means the pack is leaving its own steady value: that is a broken
    servo or a mislabelled phi_ss, and it must not be reported as an alpha."""
    g, phi = _decay(+2.5, 0.59, 0.60)
    out = S.dilatancy_rate(g, phi, 0.59)
    assert not math.isfinite(out["alpha"])
    assert out["alpha_rate"] > 0
    assert not out["alpha_sign_ok"]
    assert "GREW" in out["alpha_note"]


def test_fit_window_does_not_reach_into_the_steady_window():
    """The alpha fit uses the first ALPHA_FIT_FRAC of the strain and the steady average
    the last STEADY_FRAC. If those overlapped, alpha would be fitted partly to the
    plateau it is measured against."""
    assert S.ALPHA_FIT_FRAC + S.STEADY_FRAC <= 1.0


def test_refuses_too_few_samples_and_a_nonpositive_phi_ss():
    g, phi = _decay(-2.5, 0.59, 0.64, n=3)
    assert not math.isfinite(S.dilatancy_rate(g, phi, 0.59)["alpha"])
    g, phi = _decay(-2.5, 0.59, 0.64)
    assert not math.isfinite(S.dilatancy_rate(g, phi, 0.0)["alpha"])


def test_the_constants_are_what_the_prereg_quotes():
    """PREREG_CONSTP_2026-09-12.md §3 and §5 quote these; a silent change here would
    re-decide gates that are already written down."""
    assert S.SERVO_GAIN == 0.5
    assert S.SERVO_EDOT_CLAMP == 0.5
    assert S.SERVO_P_TOL == 0.05
    assert S.ALPHA_MIN_DECADES == 0.3
