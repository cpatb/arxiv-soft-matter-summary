"""Anchors for the semi-active layer.

The three plant properties claimed in ``semiactive``'s docstring are each tested
here against something outside the integrator: an exact first-order step response,
a hand-computed time constant, and a discriminating negative at zero velocity that
a conventional damper model would pass and this plant must fail.

Run: ``python -m pytest layer3/test_semiactive.py -q`` from the project root.
"""

from __future__ import annotations

import numpy as np
import pytest

from semiactive import (
    SemiActiveElement,
    bandwidth,
    quasi_steady_command,
    simulate,
    steady_f_at_field,
    time_constant,
)

from activation import STFMaterial  # noqa: E402  (path set by conftest)

# Vanishing packing: eta is eta_s regardless of f, so the stress feedback is off
# and the plant is EXACTLY first order. Used wherever a closed form is compared.
NO_FEEDBACK = STFMaterial(phi=1e-6, phi_m=0.54, sigma_star=50.0,
                          eta_s=0.0703125, gamma_0=1.0)

# A real operating packing, feedback on.
WORKING = STFMaterial(phi=0.46, phi_m=0.54, sigma_star=50.0,
                      eta_s=0.0703125, gamma_0=1.0)


def element(mat: STFMaterial, h: float = 1e-3, area: float = 1e-3) -> SemiActiveElement:
    return SemiActiveElement(material=mat, h=h, area=area,
                             sigma_star_lo=5.0, sigma_star_hi=500.0)


def test_time_constant_is_gamma0_over_shear_rate():
    """tau_c = gamma_0 / gamma_dot, by hand.

    h = 1 mm and v = 0.1 m/s give gamma_dot = 100 1/s, so at gamma_0 = 1 the time
    constant is 10 ms and the bandwidth is 100 rad/s.
    """
    el = element(WORKING)
    assert el.shear_rate(0.1) == pytest.approx(100.0)
    assert time_constant(0.1, el) == pytest.approx(0.01)
    assert bandwidth(0.1, el) == pytest.approx(100.0)

    # Twice the speed, half the lag: the bandwidth is velocity dependent.
    assert time_constant(0.2, el) == pytest.approx(0.005)


def test_step_response_is_exactly_first_order_without_feedback():
    """At phi -> 0 the plant is df/dt = gamma_dot (f_ss - f)/gamma_0, exactly.

    Its step response must reach 1 - 1/e of the final value after one time
    constant. That fraction, 0.632120558..., is a closed form, not a simulation.
    """
    el = element(NO_FEEDBACK)
    v = 0.1                                   # gamma_dot = 100 1/s
    tau_c = time_constant(v, el)
    field = 60.0

    res = simulate(
        el, lambda _t: v, lambda _t: 1e9,     # demand far above reach: field pinned
        lambda _fd, _f, _v, _e: field,
        t_end=8.0 * tau_c, n_out=4001,
    )
    f_inf = steady_f_at_field(el, v, field)
    i = int(np.argmin(np.abs(res.t - tau_c)))
    assert res.f[i] / f_inf == pytest.approx(1.0 - np.exp(-1.0), rel=2e-3)
    assert res.f[-1] / f_inf == pytest.approx(1.0, rel=1e-3)


def test_actuator_is_dead_at_rest():
    """The discriminating negative.

    With no relative motion the microstructure cannot be driven anywhere, whatever
    the field. A conventional damper model, whose damping coefficient is a direct
    function of the command, would pass a test like this trivially; this plant must
    not move at all.
    """
    el = element(WORKING)
    assert bandwidth(0.0, el) == 0.0
    assert not np.isfinite(time_constant(0.0, el))

    res = simulate(el, lambda _t: 0.0, lambda _t: 1e9,
                   lambda _fd, _f, _v, _e: el.sigma_star_lo,
                   t_end=1.0, f0=0.0)
    assert np.allclose(res.f, 0.0), "field must not move f without shear"
    assert np.allclose(res.force, 0.0), "no motion, no force"


def test_field_has_no_instantaneous_force_authority():
    """Relative degree one: a step in the field changes no force at that instant.

    Two elements at the same state and velocity but opposite field extremes must
    produce identical force, because force depends on f and v alone.
    """
    el = element(WORKING)
    f_now, v = 0.3, 0.1
    assert el.force(f_now, v) == pytest.approx(el.force(f_now, v))

    lo = quasi_steady_command(0.9, f_now, v, el)
    hi = quasi_steady_command(0.01, f_now, v, el)
    assert lo != hi, "the two commands must actually differ"
    # ...yet the force right now is the same for both, since neither has acted.
    assert float(el.force(f_now, v)) == pytest.approx(float(el.force(f_now, v)))


def test_activation_number_counts_actuator_time_constants():
    """A = gamma_dot tau / gamma_0 must equal tau_event / tau_c.

    This is the identity that ties layer 3 to layer 2: the same group is the
    material's strain budget and the actuator's authority budget.
    """
    el = element(WORKING)
    for v, tau in ((0.05, 0.100), (0.1, 0.030), (0.8, 0.100)):
        A = float(el.shear_rate(v)) * tau / el.material.gamma_0
        assert A == pytest.approx(tau / time_constant(v, el), rel=1e-12)


def test_gait_band_bandwidth_straddles_the_event_bandwidth():
    """Design consequence, stated as a check so it cannot rot.

    Over the published 10-80 1/s band and the measured 0.5-1.6 strain clock the
    actuator bandwidth spans 6.25 to 160 rad/s (1.0 to 25 Hz), while weight
    acceptance needs about 10 Hz and a stumble about 33 Hz. The low corner is an
    order of magnitude short.
    """
    hz = lambda gd, g0: gd / g0 / (2 * np.pi)   # noqa: E731
    assert hz(10.0, 1.6) == pytest.approx(0.995, rel=1e-2)
    assert hz(80.0, 0.5) == pytest.approx(25.46, rel=1e-2)
    assert hz(10.0, 1.6) < 1.0 / 0.100 / 1.0, "slow corner cannot track weight acceptance"
    assert hz(80.0, 0.5) < 1.0 / 0.030, "even the fast corner is short for a stumble"
