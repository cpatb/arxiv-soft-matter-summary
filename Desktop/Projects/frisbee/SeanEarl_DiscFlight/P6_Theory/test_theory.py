"""Tests for the wind-speed theory layer.

Every physical assertion here is anchored on something that does NOT pass
through the code under test: a sympy limit, a hand-constructed field with a
known answer, or a deliberately corrupted input that must be caught. Round
trips are used only where the thing being checked really is a round trip
(the ScaledField wrapper at s = 1).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "P4_Coupling"))

import dispersion_law as dl                                # noqa: E402
import sweep_theory as sw                                  # noqa: E402
import wind_scaling as ws                                  # noqa: E402


# --------------------------------------------------------------------------
# Phi(N): anchors derived symbolically, not from the numeric implementation
# --------------------------------------------------------------------------

def test_phi_limits_match_sympy():
    d = dl.check_exponents_symbolically(4.0)
    assert d["Phi_at_zero"] == pytest.approx(1.0, abs=1e-12)
    assert d["N_Phi_at_infinity"] == pytest.approx(8.0 / 3.0, rel=1e-12)


def test_phi_numeric_matches_its_own_asymptotes():
    # Small N: the series and the closed form must agree where both are valid.
    assert dl.phi(1e-4).item() == pytest.approx(1.0 - 4e-4 / 15.0, rel=1e-8)
    # Large N: N*Phi -> 8/3 from the numeric branch.
    assert dl.phi(4000.0).item() * 4000.0 == pytest.approx(8.0 / 3.0, rel=2e-3)


def test_phi_is_monotone_decreasing():
    N = np.geomspace(1e-3, 1e3, 400)
    p = dl.phi(N)
    assert np.all(np.diff(p) < 0)


def test_ballistic_V_d_exponent_is_exactly_minus_three():
    """Hand check, independent of the sympy path: halve V_d at fixed L_f and
    the predicted sigma must rise by exactly 2^3."""
    L_f, L, su = 30.0, 5.0, 0.4
    a = dl.sigma_ballistic(su, L_f / 12.0, 12.0, L_f, L)
    b = dl.sigma_ballistic(su, L_f / 6.0, 6.0, L_f, L)
    assert b / a == pytest.approx(8.0, rel=1e-12)
    # and a scalar call must return a scalar, not a length-1 array
    assert isinstance(a, float)


def test_the_two_laws_are_actually_distinguishable():
    """Discriminating negative case. If the L_f exponents of the two laws came
    out equal, every exponent test in `analyse_theory` would be vacuous."""
    L, su, V = 5.0, 0.4, 12.0
    Lf = np.array([15.0, 27.0, 50.0, 65.0])
    bal = dl.sigma_ballistic(su, Lf / V, V, Lf, L)
    tra = dl.sigma_tracer(su, Lf / V, Lf, L)
    pb, *_ = dl.fit_exponent(Lf, bal)
    pt, *_ = dl.fit_exponent(Lf, tra)
    assert pt == pytest.approx(0.5, abs=1e-9)
    assert pb > 1.4                       # near the brief's 1.7 over this range
    assert abs(pb - pt) > 0.8             # and separated by far more than noise


def test_fit_exponent_recovers_a_known_power_law():
    x = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
    y = 3.7 * x ** 1.85
    p, lo, hi, r2 = dl.fit_exponent(x, y)
    assert p == pytest.approx(1.85, abs=1e-9)
    assert r2 == pytest.approx(1.0, abs=1e-12)


def test_local_exponent_on_a_pure_power_law():
    x = np.geomspace(1.0, 100.0, 40)
    y = 2.0 * x ** -0.75
    p = dl.local_exponent(x, y)
    assert np.allclose(p, -0.75, atol=1e-6)


# --------------------------------------------------------------------------
# ScaledField
# --------------------------------------------------------------------------

class _FakeField:
    """Constant field with known U, k, eps -- so the wrapper's arithmetic is
    checked against numbers written down here, not against a CFD file."""

    def __init__(self, U=(2.0, -1.0, 0.5), k=0.6, eps=0.05):
        self._U, self._k, self._eps = np.array(U, float), k, eps
        self.heading = 90

    def U(self, r, t=0.0):
        return self._U.copy()

    def U_and_k(self, r, t=0.0):
        return self._U.copy(), self._k

    def k(self, r, t=0.0):
        return self._k

    def eps(self, r, t=0.0):
        return self._eps

    def is_solid(self, x, y):
        return False

    @property
    def fraction_clamped(self):
        return 0.0

    def bounds(self):
        return ((-1, 1), (-1, 1), (-1, 1))


def test_scaled_field_is_identity_at_the_base_speed():
    f = _FakeField()
    g = ws.ScaledField(f, ws.BASE_UREF)
    r = np.zeros(3)
    assert np.allclose(g.U(r), f.U(r))
    assert g.k(r) == pytest.approx(f.k(r))
    assert g.eps(r) == pytest.approx(f.eps(r))


def test_scaled_field_applies_the_declared_powers():
    f = _FakeField()
    g = ws.ScaledField(f, 2 * ws.BASE_UREF)
    r = np.zeros(3)
    assert np.allclose(g.U(r), 2.0 * f.U(r))
    assert g.k(r) == pytest.approx(4.0 * f.k(r))
    assert g.eps(r) == pytest.approx(8.0 * f.eps(r))
    U2, k2 = g.U_and_k(r)
    assert np.allclose(U2, 2.0 * f.U(r))
    assert k2 == pytest.approx(4.0 * f.k(r))


def test_the_integral_length_is_the_invariant_the_theory_needs():
    """sigma_u ~ s and T_L ~ 1/s must cancel EXACTLY in L = sigma_u * T_L.
    This is the single property that makes one CFD solve cover every speed."""
    f = _FakeField()
    r = np.zeros(3)
    Ls = []
    for u in (1.0, 5.0, 20.0):
        g = ws.ScaledField(f, u)
        su = np.sqrt(2.0 / 3.0 * g.k(r))
        Ls.append(su * g.k(r) / g.eps(r))
    assert Ls[0] == pytest.approx(Ls[1], rel=1e-12)
    assert Ls[1] == pytest.approx(Ls[2], rel=1e-12)


# --------------------------------------------------------------------------
# scale_invariance_report on constructed fields
# --------------------------------------------------------------------------

def _write_field(path, s, corrupt_k=1.0, seed=0):
    """Write a small npz in the P3 schema, exactly scaled by `s`.

    `corrupt_k` multiplies k by an extra factor, so a test can check the
    report NOTICES a violation rather than only that it passes on clean input.
    """
    rng = np.random.default_rng(seed)
    nx, ny, nz = 12, 9, 4
    x = np.linspace(-40, 40, nx)
    y = np.linspace(-20, 20, ny)
    z = np.array([1.0, 2.0, 4.0, 6.0])
    U0 = rng.normal(size=(nx, ny, nz, 3))
    k0 = rng.uniform(0.05, 0.5, size=(nx, ny, nz))
    e0 = rng.uniform(0.01, 0.2, size=(nx, ny, nz))
    solid = np.zeros((nx, ny), bool)
    np.savez(path, x=x, y=y, z=z, is_stub=False,
             **{"dir_090__U": s * U0,
                "dir_090__k": corrupt_k * s ** 2 * k0,
                "dir_090__eps": s ** 3 * e0,
                "dir_090__solid": solid},
             headings=np.array([90]), meta_json="{}")


def test_report_is_zero_on_an_exactly_scaled_field(tmp_path):
    a, b = tmp_path / "u5.npz", tmp_path / "u10.npz"
    _write_field(a, 1.0)
    _write_field(b, 2.0)
    rep = ws.scale_invariance_report({5.0: str(a), 10.0: str(b)}, 90,
                                     pitch_x=(-40, 40), pitch_y=(-20, 20))
    for f in rep["cases"]["10.0"]["fields"].values():
        assert f["pitch_max"] == pytest.approx(0.0, abs=1e-12)


def test_report_catches_a_violation_it_should_not_miss(tmp_path):
    """Discriminating negative case: k scaled by s^2 * 1.10 must show up as a
    ~10% departure, not be absorbed."""
    a, b = tmp_path / "u5.npz", tmp_path / "u10.npz"
    _write_field(a, 1.0)
    _write_field(b, 2.0, corrupt_k=1.10)
    rep = ws.scale_invariance_report({5.0: str(a), 10.0: str(b)}, 90,
                                     pitch_x=(-40, 40), pitch_y=(-20, 20))
    f = rep["cases"]["10.0"]["fields"]
    assert f["U"]["pitch_max"] == pytest.approx(0.0, abs=1e-12)
    assert f["eps"]["pitch_max"] == pytest.approx(0.0, abs=1e-12)
    # k is uniform(0.05, 0.5); a 10% error relative to the field RMS is order
    # 10% * mean/rms, comfortably above any rounding.
    assert f["k"]["pitch_median"] > 0.03


def test_invariance_verdict_fails_when_the_exponent_is_wrong():
    good = [dict(u_ref=u, U_mean=u, sigma_u=0.05 * u, T_L=10.0 / u,
                 L_int=0.5, TI=0.05) for u in (2.5, 5.0, 10.0)]
    assert ws.invariance_verdict(good)["passes"]
    # sigma_u growing as u^1.5 must be rejected at the +/- 0.10 tolerance.
    bad = [dict(u_ref=u, U_mean=u, sigma_u=0.05 * u ** 1.5, T_L=10.0 / u,
                L_int=0.5, TI=0.05) for u in (2.5, 5.0, 10.0)]
    v = ws.invariance_verdict(bad)
    assert not v["passes"]
    assert not v["quantities"]["sigma_u"]["passes"]


# --------------------------------------------------------------------------
# the paired estimator
# --------------------------------------------------------------------------

def test_paired_effect_separates_bias_from_spread():
    """Hand-constructed displacement: a fixed 2 m downrange shift plus a
    +/- 1 m crossrange alternation. Bias and sigma are both known exactly."""
    n = 8
    base = np.zeros((n, 3))
    d = np.zeros((n, 3))
    d[:, 0] = 2.0
    d[:, 1] = np.where(np.arange(n) % 2 == 0, 1.0, -1.0)
    e_down, e_cross = sw.throw_frame(0.0)
    got = sw.paired_effect({"land": base + d}, {"land": base}, e_down, e_cross)
    assert got["bias_down"] == pytest.approx(2.0)
    assert got["bias_cross"] == pytest.approx(0.0)
    assert got["bias"] == pytest.approx(2.0)
    assert got["sigma"] == pytest.approx(1.0)


def test_paired_effect_respects_the_throw_frame():
    """The same displacement read at bearing 90 deg must swap the two
    components. A frame bug here would silently relabel downrange as
    crossrange, which is exactly the mistake `monte_carlo.py`'s docstring
    warns about."""
    n = 6
    base = np.zeros((n, 3))
    d = np.zeros((n, 3))
    d[:, 0] = 3.0                      # 3 m east
    e0 = sw.throw_frame(0.0)
    e90 = sw.throw_frame(90.0)
    a = sw.paired_effect({"land": base + d}, {"land": base}, *e0)
    b = sw.paired_effect({"land": base + d}, {"land": base}, *e90)
    assert a["bias_down"] == pytest.approx(3.0)
    assert a["bias_cross"] == pytest.approx(0.0, abs=1e-12)
    assert b["bias_down"] == pytest.approx(0.0, abs=1e-12)
    assert b["bias_cross"] == pytest.approx(-3.0)


def test_path_field_stats_weight_by_arc_length():
    """Hand-computed anchor. Two zones, k = 1.0 for x > 0 and k = 0.01 for
    x <= 0, equal segment lengths, so the answer is a weighted mean whose
    weights can be counted off by hand.

    The sampled points are r[idx] for idx = 0..9, i.e. x = -1, 0, 1, ..., 8.
    Eight of the ten are in the high-k zone, so at Uref = 5 (s = 1)

        sigma_u = 0.8*sqrt(2/3 * 1.0) + 0.2*sqrt(2/3 * 0.01) = 0.6695272

    An unweighted or off-by-one implementation gives 0.75 or 0.60 of the high
    value, both well outside the tolerance below."""
    class TwoZone:
        def k(self, p, t=0.0):
            return 1.0 if p[0] > 0 else 0.01

        def eps(self, p, t=0.0):
            return 0.1
    r = np.column_stack([np.linspace(-1.0, 9.0, 11), np.zeros(11), np.zeros(11)])
    seg = np.ones(10)
    su, L = sw._path_field_stats(TwoZone(), 5.0, r, seg, stride=1)
    expect = 0.8 * np.sqrt(2 / 3 * 1.0) + 0.2 * np.sqrt(2 / 3 * 0.01)
    assert su == pytest.approx(expect, rel=1e-12)
    # and it must be a WEIGHTED mean, not the plain average of the two zones
    assert su != pytest.approx(0.5 * (np.sqrt(2 / 3) + np.sqrt(2 / 3 * 0.01)),
                               rel=1e-3)
