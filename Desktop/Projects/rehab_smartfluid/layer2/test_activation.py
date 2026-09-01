"""Anchors for the strain-clocked shear-thickening element.

Every check here is either a closed form derived by hand, an analytic limit, an
exact symmetry of the model, or a discriminating negative. None of them invert
:mod:`activation` against itself: a round-trip test would certify the integrator
and say nothing about the physics.

Run: ``python -m pytest layer2/test_activation.py -q`` from the project root.
"""

from __future__ import annotations

import numpy as np
import pytest

from activation import (
    JAM_STANDOFF,
    STFMaterial,
    activation_number,
    quasi_steady_ratio,
    reduced_packing,
    reduced_rate,
    steady_branch_count,
    steady_curve,
    steady_stress,
    transient_response,
)

# Reference shear rate for the fixtures, inside the 10-80 1/s band measured for
# the published STF prosthetic foot (Smart Mater. Struct. 2020).
GAMMA_DOT_REF = 100.0

# The carrier viscosity is set BY CONSTRUCTION so that the lubricated stress at
# GAMMA_DOT_REF equals sigma_star exactly, i.e. the fixture sits at the onset:
#
#   eta_s * (1 - phi/phi_0)^-n * gamma_dot = sigma_star
#   eta_s = 50 / (7.111... * 100) = 0.0703125 Pa.s
#
# This is not a tuned number, it is the solution of that equation, and
# test_fixture_sits_at_the_onset re-derives it. The first version of this file
# used a water-like 1e-3 Pa.s, which put sigma at 1.4 Pa against sigma_star = 50
# so f_ss was 6e-16 and NOTHING thickened: every test passed and the
# discriminating negative below was the only one that noticed. The corroboration
# is that 0.07 Pa.s is PEG-like rather than water-like, which is what the
# published STF devices actually use.
ETA_S_ONSET = 0.0703125

# A packing well below jamming, so the trajectory stays finite and the shortfall
# is not masked by a divergence. Not a fit to any material.
DILUTE = STFMaterial(phi=0.40, phi_m=0.54, sigma_star=50.0,
                     eta_s=ETA_S_ONSET, gamma_0=1.0)

# Close enough to phi_m that the divergence amplifies the frictional fraction.
# Used for the amplification and jamming checks.
NEAR_JAM = STFMaterial(phi=0.56, phi_m=0.54, sigma_star=50.0,
                       eta_s=ETA_S_ONSET, gamma_0=1.0)


def test_fixture_sits_at_the_onset():
    """Guard on the fixture itself, so it cannot drift back into a dead regime.

    If the lubricated stress at the reference rate is far below sigma_star the
    material never thickens and every thickening test passes vacuously.
    """
    eta_lub = DILUTE.eta_s * (1.0 - DILUTE.phi / DILUTE.phi_0) ** -DILUTE.n
    sigma_lub = eta_lub * GAMMA_DOT_REF
    assert sigma_lub == pytest.approx(DILUTE.sigma_star, rel=1e-9)
    assert float(DILUTE.f_steady(sigma_lub)) == pytest.approx(np.exp(-1.0), rel=1e-9)


def test_f_jam_is_hand_computed():
    """phi_J(f) = (1-f) phi_0 + f phi_m solved for f, by hand.

    phi_0 = 0.64, phi_m = 0.54, phi = 0.60 gives f = 0.04/0.10 = 0.4 exactly.
    """
    mat = STFMaterial(phi=0.60, phi_m=0.54, sigma_star=50.0, eta_s=1e-3,
                      gamma_0=1.0, phi_0=0.64)
    assert mat.f_jam() == pytest.approx(0.4, abs=1e-12)

    # A packing below phi_m can never jam: no f in [0, 1] solves it.
    assert np.isnan(DILUTE.f_jam())


def test_strain_clock_matches_closed_form_at_vanishing_packing():
    """As phi -> 0 the viscosity is eta_s regardless of f, so the ODE is exactly solvable.

    Then sigma = eta_s * gamma_dot is constant, f_ss is constant, and

        f(t) = f_ss * (1 - exp(-gamma_dot t / gamma_0))     with f(0) = 0

    which is a hand integration, not a second call to the integrator.
    """
    mat = STFMaterial(phi=1e-6, phi_m=0.54, sigma_star=50.0, eta_s=1e-3, gamma_0=1.0)
    gamma_dot, t_end = 1.0e5, 3.0e-5           # 3 strain units at gamma_0 = 1

    res = transient_response(lambda _t: gamma_dot, mat, t_end)

    sigma_const = mat.eta_s * gamma_dot        # 100 Pa, so f_ss = exp(-0.5)
    f_ss = np.exp(-mat.sigma_star / sigma_const)
    expected = f_ss * (1.0 - np.exp(-gamma_dot * res.t / mat.gamma_0))

    assert np.max(np.abs(res.f - expected)) < 1e-6
    assert f_ss == pytest.approx(0.6065306597, rel=1e-9)   # exp(-1/2), by hand


def test_quasi_steady_limit_is_recovered():
    """A -> infinity must reproduce the steady flow curve. This is the limit the
    device literature assumes always holds."""
    gamma_dot = GAMMA_DOT_REF
    tau_long = 4000.0 * DILUTE.gamma_0 / gamma_dot        # A = 4000
    ratio, _, _ = quasi_steady_ratio(gamma_dot, tau_long, DILUTE)
    assert ratio == pytest.approx(1.0, rel=1e-3)


def test_frozen_limit_is_the_lubricated_viscosity():
    """A -> 0 must give the f = 0 viscosity, eta_s (1 - phi/phi_0)^-n, by hand.

    phi = 0.40, phi_0 = 0.64, n = 2: (1 - 0.625)^-2 = 1/0.140625 = 7.1111...
    """
    gamma_dot = GAMMA_DOT_REF
    tau_tiny = 1e-8 * DILUTE.gamma_0 / gamma_dot          # A = 1e-8
    _, sigma_tr, _ = quasi_steady_ratio(gamma_dot, tau_tiny, DILUTE)

    eta_lub_hand = DILUTE.eta_s * (1.0 - 0.40 / 0.64) ** -2.0
    assert eta_lub_hand == pytest.approx(0.5, rel=1e-9)   # 0.0703125 * 64/9
    assert sigma_tr == pytest.approx(eta_lub_hand * gamma_dot, rel=1e-6)


def test_shortfall_is_strict_at_order_one_activation():
    """The discriminating negative.

    If this module had quietly been written quasi-steady the ratio would be 1
    everywhere and every other test here would still pass. At A of order one it
    must be strictly below 1, and it must fall as A falls.
    """
    gamma_dot = GAMMA_DOT_REF
    ratios = []
    for A in (0.3, 1.0, 3.0, 10.0):
        tau = A * DILUTE.gamma_0 / gamma_dot
        r, _, _ = quasi_steady_ratio(gamma_dot, tau, DILUTE)
        ratios.append(r)

    assert all(r < 1.0 for r in ratios), ratios
    assert np.all(np.diff(ratios) > 0), f"shortfall must ease as A grows: {ratios}"
    assert ratios[0] < 0.95, f"A = 0.3 should be clearly short, got {ratios[0]}"


def test_gauge_covariance_holds_to_machine_precision():
    """Wyart-Cates is exactly covariant under (phi, phi_m, phi_0) -> c (phi, phi_m, phi_0).

    dst-fit verified this to 2e-15 and it invalidated three of its own claims. If
    this module broke the covariance it would have introduced a phi_0 dependence
    that does not exist in the model, so this is a check on the port, not on the
    theory.
    """
    c = 1.1
    scaled = STFMaterial(
        phi=DILUTE.phi * c, phi_m=DILUTE.phi_m * c, phi_0=DILUTE.phi_0 * c,
        sigma_star=DILUTE.sigma_star, eta_s=DILUTE.eta_s, gamma_0=DILUTE.gamma_0,
    )
    kw = dict(t_end=0.05)
    a = transient_response(lambda _t: GAMMA_DOT_REF, DILUTE, **kw)
    b = transient_response(lambda _t: GAMMA_DOT_REF, scaled, **kw)

    assert np.max(np.abs(a.sigma - b.sigma)) / np.max(a.sigma) < 1e-12
    assert reduced_packing(DILUTE) == pytest.approx(reduced_packing(scaled), rel=1e-12)


def test_jam_standoff_does_not_move_a_reported_result():
    """The stand-off from the pole is numerical, so an unjammed result must not feel it."""
    import activation as mod

    gamma_dot, tau = 200.0, 0.01
    before = quasi_steady_ratio(gamma_dot, tau, DILUTE)[0]
    original = mod.JAM_STANDOFF
    try:
        mod.JAM_STANDOFF = 1.0 - 1e-7
        after = quasi_steady_ratio(gamma_dot, tau, DILUTE)[0]
    finally:
        mod.JAM_STANDOFF = original
    assert before == pytest.approx(after, rel=1e-9)
    assert JAM_STANDOFF < 1.0


def test_jamming_is_reported_not_absorbed():
    """A packing above phi_m must jam under enough stress, and say so."""
    assert 0.0 < NEAR_JAM.f_jam() < 1.0
    res = transient_response(lambda _t: 5.0e4, NEAR_JAM, t_end=1.0)
    assert res.jammed, "near-jamming packing driven hard should reach phi_J"
    assert np.isfinite(res.t_jam)


def test_activation_number_and_reduced_rate_are_dimensionless():
    """pint check: both groups must come out dimensionless."""
    pint = pytest.importorskip("pint")
    u = pint.UnitRegistry()

    A = (200.0 / u.s) * (0.1 * u.s) / 1.0
    assert A.dimensionless
    assert A.to("dimensionless").magnitude == pytest.approx(
        float(activation_number(200.0, 0.1, 1.0))
    )

    lam = (DILUTE.eta_s * u.Pa * u.s) * (200.0 / u.s) / (DILUTE.sigma_star * u.Pa)
    assert lam.dimensionless
    assert lam.to("dimensionless").magnitude == pytest.approx(
        float(reduced_rate(200.0, DILUTE))
    )


def test_steady_curve_is_s_shaped_for_a_thickening_packing():
    """Structural check: near jamming the flow curve must be multivalued in rate.

    If it were single-valued everywhere the material would not be discontinuously
    thickening and the whole premise would be wrong.
    """
    # The grid must lie inside the rates this packing can actually reach. Its
    # steady curve turns over at 4.14 1/s (sigma = 31.9 Pa); probing 1e2 to 1e6,
    # as the first version of this test did, samples only rates the material
    # cannot sustain and finds no branches at all.
    branches = [steady_branch_count(g, NEAR_JAM) for g in np.geomspace(0.1, 4.0, 60)]
    assert max(branches) > 1, "expected an S-shaped region for a DST packing"
    assert np.isfinite(steady_stress(1.0, NEAR_JAM))


def test_jamming_caps_the_reachable_shear_rate():
    """Design constraint, and a discriminating negative against a phi-blind model.

    A packing at or above phi_m shear jams before it reaches gait shear rates, so
    the accessible rate window closes as the thickening gets stronger. Measured
    here: 4917 1/s at phi = 0.44 falling to 7.8 1/s at phi = phi_m = 0.54, against
    a gait band of 10 to 80 1/s. Any model without a phi term (Bouc-Wen, Bingham)
    cannot express this at all.
    """
    sigma = np.geomspace(1e-2, 1e4, 200001)
    rate_max = {}
    for phi in (0.44, 0.50, 0.52, 0.54, 0.56):
        mat = STFMaterial(phi=phi, phi_m=0.54, sigma_star=50.0,
                          eta_s=ETA_S_ONSET, gamma_0=1.0)
        g, e = steady_curve(sigma, mat)
        ok = np.isfinite(e) & (g > 0)
        rate_max[phi] = float(g[ok].max())

    values = [rate_max[p] for p in (0.44, 0.50, 0.52, 0.54, 0.56)]
    assert np.all(np.diff(values) < 0), f"reachable rate must fall with phi: {rate_max}"
    assert rate_max[0.44] > 80.0, "a dilute packing must reach the gait band"
    assert rate_max[0.54] < 10.0, "at phi_m the material jams below the gait band"
