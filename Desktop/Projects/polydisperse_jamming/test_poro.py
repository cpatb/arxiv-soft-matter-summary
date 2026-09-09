"""test_poro.py — anchors and gates for the 1D consolidation solver.

Every check either compares against something OUTSIDE this code (a published analytic
series, a closed form derived by hand, an independent scalar ODE solve) or is a
discriminating negative case. A round trip through `solve` is never an anchor.

    python -m pytest test_poro.py -q
"""

from __future__ import annotations

import numpy as np
import pytest

import poro_impact as PI


# =====================================================================================
# ANCHOR 1 — the published Terzaghi consolidation series
# =====================================================================================
def test_terzaghi_series_is_reproduced_with_the_source_off():
    """Source off, uniform initial p, drained top, sealed base: the classic solution.

    This is an anchor outside the code under test: `terzaghi_series` is the published
    Fourier series (Terzaghi 1943) and `solve` never calls it.
    """
    lam = 1.0
    r = PI.solve(lam=lam, s_dil=0.0, dphi0=0.0, t_end=0.30, p_init=1.0,
                 freeze_phi=True, n_nodes=401, rtol=1e-10, atol=1e-14)
    assert r.ok
    x = np.linspace(0.0, 1.0, 401)
    want = PI.terzaghi_series(x, T=lam * 0.30)
    got = np.asarray(r.p_final_field)
    err = float(np.max(np.abs(got - want)))
    assert err < 1e-5, f"max |numeric - analytic| = {err:.3e}"


def test_terzaghi_orientation_is_not_mirrored():
    """The series must be SEALED at x=0 and DRAINED at x=1, matching the solver.

    The first version of `terzaghi_series` was the mirror image and disagreed with the
    solver by 0.61 in a field of order 1. Pinned here so the orientation cannot flip
    back silently: at early time the sealed face still carries p0 and the drained face
    is already at zero.
    """
    assert PI.terzaghi_series(np.array([0.0]), 1e-6)[0] == pytest.approx(1.0, abs=2e-3)
    assert PI.terzaghi_series(np.array([1.0]), 0.10)[0] == pytest.approx(0.0, abs=1e-12)
    # and the mirrored series must NOT satisfy that, or the test is vacuous
    x = np.array([0.0, 1.0])
    mirrored = PI.terzaghi_series(1.0 - x, 0.10)
    assert mirrored[0] != pytest.approx(PI.terzaghi_series(x, 0.10)[0], abs=1e-3)


def test_terzaghi_anchor_is_not_vacuous():
    """The same comparison at a DIFFERENT time must disagree, or the anchor is empty."""
    x = np.linspace(0.0, 1.0, 401)
    a = PI.terzaghi_series(x, T=0.30)
    b = PI.terzaghi_series(x, T=0.60)
    assert float(np.max(np.abs(a - b))) > 0.15


# =====================================================================================
# ANCHOR 2 — the quasi-steady limit reproduces the closed form, WITH its prefactor
# =====================================================================================
def test_quasistatic_limit_recovers_Pi_over_two():
    """With phi frozen and enough time, p~(0) -> -S alpha dphi/(2 Lambda).

    The 1/2 is the geometric prefactor that Eq. (6)'s "~" cannot supply, and it is
    derived by hand in `quasistatic_peak`, not fitted here.
    """
    lam, s_dil, dphi0 = 2.0, 1.0, 0.02
    r = PI.solve(lam=lam, s_dil=s_dil, dphi0=dphi0, t_end=6.0, freeze_phi=True,
                 n_nodes=401, rtol=1e-10, atol=1e-14)
    assert r.ok
    want = PI.quasistatic_peak(lam, s_dil, dphi0)
    assert r.p_base_final == pytest.approx(want, rel=2e-3), (
        f"steady base pressure {r.p_base_final:.6e} vs closed form {want:.6e}"
    )
    # and the steady profile must be the parabola, not merely the right endpoint
    x = np.linspace(0.0, 1.0, 401)
    parab = want * (1.0 - x**2)
    assert float(np.max(np.abs(np.asarray(r.p_final_field) - parab))) < 2e-3 * abs(want)


def test_quasistatic_prefactor_test_can_fail():
    """A prefactor of 1 instead of 1/2 must be rejected by the same comparison."""
    lam, s_dil, dphi0 = 2.0, 1.0, 0.02
    r = PI.solve(lam=lam, s_dil=s_dil, dphi0=dphi0, t_end=6.0, freeze_phi=True,
                 n_nodes=401, rtol=1e-10, atol=1e-14)
    wrong = 2.0 * PI.quasistatic_peak(lam, s_dil, dphi0)
    assert r.p_base_final != pytest.approx(wrong, rel=0.05)


# =====================================================================================
# THE SIGN — the whole mechanism
# =====================================================================================
def test_dilating_pack_gives_SUCTION_and_contracting_gives_OVERPRESSURE():
    kw = dict(lam=2.0, s_dil=1.0, t_end=4.0, freeze_phi=True, n_nodes=201)
    dil = PI.solve(dphi0=+0.02, **kw)
    con = PI.solve(dphi0=-0.02, **kw)
    assert dil.p_base_final < 0, "phi > phi_c must give suction"
    assert con.p_base_final > 0, "phi < phi_c must give overpressure"
    # and the magnitudes must match: the mechanism is a sign flip, not an asymmetry,
    # at this level of the model
    assert abs(dil.p_base_final) == pytest.approx(abs(con.p_base_final), rel=1e-6)


def test_zero_dphi_gives_no_pore_pressure_at_all():
    r = PI.solve(lam=2.0, s_dil=1.0, dphi0=0.0, t_end=4.0, freeze_phi=True)
    assert abs(r.p_base_final) < 1e-12


# =====================================================================================
# THE STRAIN CLOCK — checked against an independent scalar integration
# =====================================================================================
def test_phi_relaxes_to_phi_c_on_the_strain_the_scalar_ODE_predicts():
    """dphi/dt = -phi alpha (phi - phi_c) has nothing to do with the PDE; solving it
    separately is an independent check of the source's shut-off clock."""
    from scipy.integrate import solve_ivp

    phi0, dphi0, alpha = 0.44, 0.02, 1.0
    phi_c = phi0 - dphi0
    sol = solve_ivp(lambda t, y: -y * alpha * (y - phi_c), (0, 6), [phi0],
                    rtol=1e-11, atol=1e-14, dense_output=True)
    d = lambda t: float(sol.sol(t)[0]) - phi_c          # noqa: E731
    lo, hi = 0.0, 6.0
    for _ in range(200):                                # bisect the half-crossing
        mid = 0.5 * (lo + hi)
        if d(mid) > 0.5 * dphi0:
            lo = mid
        else:
            hi = mid
    t_half_ref = 0.5 * (lo + hi)

    r = PI.solve(lam=2.0, s_dil=1.0, dphi0=dphi0, phi0=phi0, t_end=6.0,
                 n_nodes=101, rtol=1e-10, atol=1e-14)
    assert r.ok
    # solve_ivp's output grid limits the resolution of the crossing; require the
    # crossing to be inside one output step of the reference.
    step = float(np.max(np.diff(r.t)))
    assert abs(r.strain_to_half_dphi - t_half_ref) <= step, (
        f"strain to half dphi {r.strain_to_half_dphi:.4f} vs reference "
        f"{t_half_ref:.4f}, output step {step:.4f}"
    )
    assert 1.0 < t_half_ref < 2.5, f"unexpected reference clock {t_half_ref}"


def test_free_phi_never_reaches_the_frozen_steady_state():
    """The source shuts itself off, so the transient must UNDERSHOOT the closed form.

    This is the quantitative content of assumption A5, and it is the reason the
    derivation's Eq. (6) is an upper bound rather than a prediction.
    """
    lam, s_dil, dphi0 = 2.0, 1.0, 0.02
    free = PI.solve(lam=lam, s_dil=s_dil, dphi0=dphi0, t_end=2.0, n_nodes=201)
    frozen = PI.solve(lam=lam, s_dil=s_dil, dphi0=dphi0, t_end=2.0, freeze_phi=True,
                      n_nodes=201)
    assert abs(free.p_peak_hat) < abs(frozen.p_peak_hat)
    assert 0.0 < abs(free.steady_ratio) < 1.0


# =====================================================================================
# CAVITATION
# =====================================================================================
def test_cavitation_ceiling_binds_and_is_reported():
    kw = dict(lam=0.05, s_dil=1.0, dphi0=0.02, t_end=2.0, freeze_phi=True, n_nodes=201)
    free = PI.solve(**kw)
    capped = PI.solve(p_cav_hat=0.01, **kw)
    assert abs(free.p_base_final) > 0.01
    assert capped.cavitated
    assert capped.p_base_final >= -0.01 - 1e-9
    assert capped.frac_time_cavitated > 0.0


def test_cavitation_does_not_bind_when_the_ceiling_is_high():
    r = PI.solve(lam=2.0, s_dil=1.0, dphi0=0.02, t_end=2.0, freeze_phi=True,
                 p_cav_hat=10.0)
    assert not r.cavitated
    assert r.frac_time_cavitated == 0.0


# =====================================================================================
# DISCRETISATION
# =====================================================================================
def test_grid_convergence_is_second_order_or_better():
    kw = dict(lam=1.0, s_dil=0.0, dphi0=0.0, t_end=0.30, p_init=1.0, freeze_phi=True,
              rtol=1e-11, atol=1e-14)
    errs = []
    for n in (101, 201, 401):
        r = PI.solve(n_nodes=n, **kw)
        x = np.linspace(0.0, 1.0, n)
        errs.append(float(np.max(np.abs(np.asarray(r.p_final_field)
                                        - PI.terzaghi_series(x, T=0.30)))))
    p1 = np.log2(errs[0] / errs[1])
    p2 = np.log2(errs[1] / errs[2])
    assert p1 > 1.7 and p2 > 1.7, f"observed orders {p1:.2f}, {p2:.2f}; errs {errs}"


# =====================================================================================
# DIMENSIONS
# =====================================================================================
def test_physical_numbers_are_dimensionally_consistent():
    import pint

    u = pint.UnitRegistry()
    kappa = 3.225e-15 * u.m**2
    mu_w = 1.0e-3 * u.Pa * u.s
    kdr = 1.0e6 * u.Pa
    rho = 1440.0 * u.kg / u.m**3
    u0 = 3.0 * u.m / u.s
    d0 = 2.5e-3 * u.m
    ell = 2.5e-3 * u.m
    c_v = kappa * kdr / mu_w
    assert c_v.dimensionality == (u.m**2 / u.s).dimensionality
    lam = c_v * d0 / (u0 * ell**2)
    s = kdr * 1.0 / (rho * u0**2)
    assert lam.dimensionality == u.dimensionless.dimensionality
    assert s.dimensionality == u.dimensionless.dimensionality

    got = PI.physical_numbers(kappa=3.225e-15, mu_w=1.0e-3, k_drained=1.0e6,
                              rho=1440.0, u0=3.0, d0=2.5e-3, ell=2.5e-3)
    assert got["Lambda"] == pytest.approx(float(lam.magnitude), rel=1e-12)
    assert got["S"] == pytest.approx(float(s.magnitude), rel=1e-12)


def test_quasisteady_peak_equals_Pi_over_two_in_physical_variables():
    """S*alpha*dphi/(2 Lambda) must equal Pi/2 with Pi from the derivation's Eq. (6).

    K' cancels, which is the point: the steady pore pressure does not depend on the
    skeleton stiffness, only the time to reach it does.
    """
    kappa, mu_w, kdr, rho, u0, d0, ell = 3.225e-15, 1e-3, 1e6, 1440.0, 3.0, 2.5e-3, 2.5e-3
    n = PI.physical_numbers(kappa, mu_w, kdr, rho, u0, d0, ell)
    dphi = 0.02
    from_solver = abs(PI.quasistatic_peak(n["Lambda"], n["S"], dphi))
    pi_closed = PI.ALPHA_DIL * mu_w * dphi * ell**2 / (rho * kappa * u0 * d0)
    assert from_solver == pytest.approx(pi_closed / 2.0, rel=1e-12)
