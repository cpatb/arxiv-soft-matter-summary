"""test_shear.py — engine-independent tests for the phase-8 critical-state layer.

Runs on Windows without Yade: everything here tests the bookkeeping, the estimators and
the gates, on synthetic ladders whose phi_c is KNOWN. The Yade half is checked by the
smoke cell in `phase8_criticalstate.py --smoke`.

    python -m pytest test_shear.py -q
"""

from __future__ import annotations

import csv
import os

import numpy as np
import pytest

import analyse_phase8 as A
import phase8_criticalstate as P8

HERE = os.path.dirname(os.path.abspath(__file__))


# =====================================================================================
# the ladder and its centres
# =====================================================================================
def test_ladder_centres_match_phase2_summary():
    """PHI_J_ISO must be the project's OWN published phi_J, not a remembered number.

    This is the check that stops the ladder drifting away from the measurement it is
    centred on. It reads out/phase2_summary.csv, which is a committed artifact.
    """
    path = os.path.join(HERE, "out", "phase2_summary.csv")
    assert os.path.exists(path), "out/phase2_summary.csv is missing"
    table = {}
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if int(row["n_requested"]) != 1500:
                continue
            key = (row["family"], round(float(row["cv_nom"]), 3),
                   round(float(row["mu_p"]), 3))
            table[key] = float(row["phi_J"])
    missing = [k for k in P8.PHI_J_ISO if k not in table]
    assert not missing, f"PHI_J_ISO keys absent from phase2_summary.csv: {missing}"
    for k, v in P8.PHI_J_ISO.items():
        assert v == pytest.approx(table[k], abs=1e-12), (
            f"{k}: PHI_J_ISO has {v!r}, phase2_summary.csv has {table[k]!r}"
        )


def test_ladder_rungs_are_distinct_and_wide_enough():
    """Distinct rungs, uniform spacing, and wide enough to resolve the seed scatter.

    The largest seed scatter in the published phi_J table is 0.00475 (lognormal,
    cv 0.5, mu 0.3). A ladder spacing below that cannot resolve rungs at all.
    """
    grid = P8.build_grid("lognormal", (0.0,), (0.3,), (1,), 1500)
    phis = sorted(g[4] for g in grid)
    assert len(phis) == len(P8.PHI_OFFSETS)
    gaps = np.diff(phis)
    # one uniform spacing, whatever it is, and it must exceed the largest seed scatter
    # in the published phi_J table (0.00475) or the rungs cannot be told apart
    assert np.allclose(gaps, gaps[0], atol=1e-9)
    assert gaps.min() > 0.00475
    # the ladder must cover at least +/- one spacing around its centre, or a root
    # cannot be bracketed at all
    assert max(P8.PHI_OFFSETS) - min(P8.PHI_OFFSETS) >= 2 * gaps[0]


def test_ladder_centre_refuses_an_unknown_cell():
    """A missing centre must RAISE, not fall back on a guess."""
    with pytest.raises(KeyError):
        P8.ladder_centre("lognormal", 0.35, 0.3)


def test_cell_key_separates_every_swept_knob():
    """Two cells differing in ANY swept knob must not share a resume key.

    The N-pooling defect of 2026-09-08 was exactly a key that did not encode a knob
    the analysis grouped by.
    """
    base = dict(family="lognormal", cv=0.0, mu=0.3, seed=1, phi=0.6,
                n=1500, gamma_max=2.0, shear_steps=300000)
    k0 = P8.cell_key(**base)
    for knob, other in [("family", "uniform"), ("cv", 0.5), ("mu", 0.0), ("seed", 2),
                        ("phi", 0.615), ("n", 3000), ("gamma_max", 1.0),
                        ("shear_steps", 200000)]:
        alt = dict(base)
        alt[knob] = other
        assert P8.cell_key(**alt) != k0, f"key does not encode {knob}"


# =====================================================================================
# the gates stay in sync
# =====================================================================================
def test_i_max_agrees_between_producer_and_analyser():
    import shear
    assert shear.I_MAX == A.I_FIT_MAX


# =====================================================================================
# the four estimators, on synthetic ladders with a KNOWN phi_c
# =====================================================================================
def _synthetic(phi_c: float, b: float = 1.0, c: float = 0.5, gdot: float = 2.6e-5,
               d50: float = 1.0, offsets=(-0.045, -0.03, -0.015, 0, .015, .03, .045),
               centre: float | None = None):
    """A ladder obeying p = c (phi - phi_c)^b exactly, with I derived from it."""
    centre = phi_c + 0.02 if centre is None else centre
    phi = np.array([centre + o for o in offsets])
    p = np.where(phi > phi_c, c * np.maximum(phi - phi_c, 0) ** b, 0.0)
    with np.errstate(divide="ignore"):
        ivals = np.where(p > 0, gdot * d50 / np.sqrt(p), np.inf)
    return phi, p, ivals


def test_E2_recovers_a_known_root_and_exponent():
    phi, p, _ = _synthetic(0.5800, b=1.0)
    out = A._pressure_root(phi, p, None)
    assert out["phi_c_E2"] == pytest.approx(0.5800, abs=2e-4)
    assert out["b_E2"] == pytest.approx(1.0, abs=0.02)


def test_E2_recovers_a_DIFFERENT_root_discriminating_case():
    """The estimator must move when the truth moves -- otherwise it measures nothing."""
    a = A._pressure_root(*_synthetic(0.5800)[:2], None)["phi_c_E2"]
    b = A._pressure_root(*_synthetic(0.6200)[:2], None)["phi_c_E2"]
    assert abs(b - a) == pytest.approx(0.04, abs=1e-3)


def test_E2_recovers_a_non_unit_exponent():
    phi, p, _ = _synthetic(0.5800, b=2.5)
    out = A._pressure_root(phi, p, None)
    assert out["b_E2"] == pytest.approx(2.5, abs=0.05)
    assert out["phi_c_E2"] == pytest.approx(0.5800, abs=5e-4)


def test_E3_pinned_exponent_is_biased_when_the_truth_is_not_unity():
    """E3 pins b = 1 by convention; on b = 2.5 data it MUST be wrong.

    Recorded so that E3 is never read as an independent measurement of phi_c when the
    contact law is not harmonic -- the failure is the point of the test.
    """
    phi, p, _ = _synthetic(0.5800, b=2.5)
    out = A._pressure_root(phi, p, 1.0)
    assert abs(out["phi_c_E3"] - 0.5800) > 1e-3


def test_E1_recovers_phi_c_from_the_phi_of_I_curve():
    """phi(I) with a linear dilatancy law extrapolates to the right intercept."""
    ivals = np.linspace(1e-4, 8e-3, 7)
    phi_c, a = 0.5800, 1.5
    phi = phi_c - a * ivals
    out = A.e1_phi_of_I(phi, ivals)
    assert out["phi_c_E1"] == pytest.approx(phi_c, abs=1e-9)
    assert out["slope_E1"] == pytest.approx(-a, rel=1e-9)
    assert out["n_E1"] == 7


def test_E1_ignores_rungs_above_the_quasistatic_bar():
    ivals = np.array([1e-4, 1e-3, 5e-3, 5e-2, 2e-1])
    phi = 0.58 - 1.5 * ivals
    phi[3:] += 0.05                       # corrupt only the non-quasi-static rungs
    out = A.e1_phi_of_I(phi, ivals)
    assert out["n_E1"] == 3
    assert out["phi_c_E1"] == pytest.approx(0.58, abs=1e-9)


def test_E4_finds_the_isostatic_crossing_and_returns_nan_when_uncrossed():
    phi = np.array([0.55, 0.57, 0.59, 0.61])
    z = np.array([2.0, 3.0, 5.0, 6.0])
    got = A.e4_z_crossing(phi, z)["phi_c_E4"]
    assert got == pytest.approx(0.58, abs=1e-9)
    assert np.isnan(A.e4_z_crossing(phi, np.array([5.0, 5.5, 6.0, 6.5]))["phi_c_E4"])
    assert np.isnan(A.e4_z_crossing(phi, np.array([1.0, 1.5, 2.0, 2.5]))["phi_c_E4"])


# =====================================================================================
# the refusals
# =====================================================================================
def _recs_from(phi, p, ivals, z, mu=0.36, steady=True, qs=True):
    return [{"phi": float(a), "p_ss_hat": float(b), "i_ss": float(c),
             "z_nr_ss": float(d), "mu_macro": mu, "steady": steady,
             "quasistatic": qs} for a, b, c, d in zip(phi, p, ivals, z)]


def test_a_root_far_below_the_ladder_is_REFUSED():
    """A fit whose root sits below the ladder is an extrapolation, not a measurement."""
    phi = np.linspace(0.60, 0.69, 7)
    p = 0.5 * (phi - 0.40)                       # root 0.20 below the ladder
    ivals = np.full_like(phi, 1e-3)
    z = np.linspace(5.0, 6.5, 7)
    res = A.analyse_group(_recs_from(phi, p, ivals, z))
    assert np.isnan(res["phi_c_E2"])
    assert res["phi_c_E2_refused"] == pytest.approx(0.40, abs=5e-3)


def test_duplicate_realised_phi_is_detectable():
    """The 2026-09-08 compression defect: two rungs landing on the same state."""
    phis = [0.5595, 0.5745, 0.5745]           # the middle pair collapsed
    dupes = sum(1 for a, b in zip(sorted(phis), sorted(phis)[1:]) if abs(a - b) < 1e-4)
    assert dupes == 1
    ok = [0.5595, 0.5745, 0.5895]
    assert sum(1 for a, b in zip(ok, ok[1:]) if abs(a - b) < 1e-4) == 0


def test_analyse_group_counts_the_gates_it_was_given():
    phi = np.linspace(0.56, 0.62, 7)
    p = 0.5 * np.maximum(phi - 0.58, 0) + 1e-12
    ivals = np.full_like(phi, 1e-3)
    z = np.linspace(2.0, 6.0, 7)
    recs = _recs_from(phi, p, ivals, z)
    recs[0]["steady"] = False
    recs[1]["quasistatic"] = False
    res = A.analyse_group(recs)
    assert res["n_cells"] == 7
    assert res["n_steady"] == 6
    assert res["n_quasistatic"] == 6
