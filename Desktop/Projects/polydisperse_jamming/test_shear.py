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


def test_group_key_separates_protocols_and_run_parameters():
    """Two records from DIFFERENT protocols must never land in the same group.

    The default glob `phase8_prod*.jsonl` matches both the v1 and the v2 production
    runs, which differ in the periodic-cell flip fix, the ladder, gamma_max and the step
    count. This is the N-pooling defect of 2026-09-08 in a new place, and the key is
    where it is stopped.
    """
    base = {"family": "lognormal", "cv": 0.2, "mu_p": 0.3, "seed": 1,
            "n_requested": 1500, "protocol": "cv-shear-v2", "gamma_max": 1.5,
            "shear_steps": 200000}
    k0 = A.group_key(base)
    for knob, other in [("protocol", "cv-shear-v1"), ("gamma_max", 2.0),
                        ("shear_steps", 300000), ("cv", 0.5), ("seed", 2),
                        ("n_requested", 3000), ("mu_p", 0.0), ("family", "uniform")]:
        alt = dict(base)
        alt[knob] = other
        assert A.group_key(alt) != k0, f"group_key does not encode {knob}"


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


def test_E1_returns_the_DENSE_end_intercept_not_phi_c():
    """E1 is arithmetically a linear intercept at I = 0 -- and that is the DENSE end.

    Renamed and re-scoped 2026-09-08 when E1 was retracted as a phi_c estimator. The
    arithmetic is unchanged and still checked; what changed is the NAME of what it
    returns, so that no caller can read it as phi_c again.
    """
    ivals = np.linspace(1e-4, 8e-3, 7)
    intercept, a = 0.5800, 1.5
    phi = intercept - a * ivals
    out = A.e1_phi_of_I(phi, ivals)
    assert "phi_c_E1" not in out, "E1 must not expose a key that reads as phi_c"
    assert out["phi_dense_E1"] == pytest.approx(intercept, abs=1e-9)
    assert out["slope_E1"] == pytest.approx(-a, rel=1e-9)
    assert out["n_E1"] == 7


def test_on_a_FIXED_RATE_ladder_small_I_is_the_DENSE_end():
    """The fact that retired E1, pinned so it cannot be silently reverted.

    I = gdot * d50 * sqrt(rho/p). On this protocol gdot is set once per cell from
    gamma_max and the step count, so it is IDENTICAL across a ladder and I varies only
    through p. Since p RISES with phi, I FALLS with phi -- so extrapolating any quantity
    to I = 0 targets the densest rung, not the critical state.
    """
    gdot, d50 = 3.0e-5, 1.0
    phi = np.array([0.5745, 0.5945, 0.6145, 0.6345, 0.6545])
    p = np.array([3.40e-6, 2.15e-4, 5.81e-3, 1.06e-2, 1.65e-2])   # measured, cv = 0
    ivals = gdot * d50 / np.sqrt(p)
    assert np.all(np.diff(ivals) < 0), "I must FALL as phi rises on a fixed-rate ladder"
    assert ivals.argmin() == len(phi) - 1, "smallest I is the DENSEST rung"
    # and the intercept is therefore ABOVE the top rung, not below the bottom one
    out = A.e1_phi_of_I(phi, ivals)
    assert out["phi_dense_E1"] > phi.max() - 0.02


def test_E1_ignores_rungs_above_the_quasistatic_bar():
    ivals = np.array([1e-4, 1e-3, 5e-3, 5e-2, 2e-1])
    phi = 0.58 - 1.5 * ivals
    phi[3:] += 0.05                       # corrupt only the non-quasi-static rungs
    out = A.e1_phi_of_I(phi, ivals)
    assert out["n_E1"] == 3
    assert out["phi_dense_E1"] == pytest.approx(0.58, abs=1e-9)


def test_E4_isostatic_target_depends_on_FRICTION():
    """6 for frictionless spheres, 4 for frictional ones -- not one constant.

    Hardcoding 4 was correct for the mu_p = 0.3 production arm and silently wrong for the
    mu_p = 0 anchor run, where the pack is isostatic at 6. A constant that is right for
    one arm and wrong for another is the trap this project keeps meeting, so the choice
    is pinned here.
    """
    phi = np.array([0.613, 0.623, 0.633, 0.643, 0.653])
    z = np.array([4.4, 5.0, 5.8, 6.4, 6.8])
    frictionless = A.e4_z_crossing(phi, z, mu_p=0.0)
    frictional = A.e4_z_crossing(phi, z, mu_p=0.3)
    assert frictionless["z_iso_used"] == 6.0
    assert frictional["z_iso_used"] == 4.0
    # on this ladder Z never falls to 4, so the frictional target must find NO crossing
    assert np.isfinite(frictionless["phi_c_E4"])
    assert np.isnan(frictional["phi_c_E4"])
    # and the frictionless crossing must land where Z = 6, not somewhere else
    assert frictionless["phi_c_E4"] == pytest.approx(0.6363, abs=1e-3)


def test_E4_finds_the_isostatic_crossing_and_returns_nan_when_uncrossed():
    phi = np.array([0.55, 0.57, 0.59, 0.61])
    z = np.array([2.0, 3.0, 5.0, 6.0])
    got = A.e4_z_crossing(phi, z, mu_p=0.3)["phi_c_E4"]
    assert got == pytest.approx(0.58, abs=1e-9)
    assert np.isnan(A.e4_z_crossing(phi, np.array([5.0, 5.5, 6.0, 6.5]),
                                    mu_p=0.3)["phi_c_E4"])
    assert np.isnan(A.e4_z_crossing(phi, np.array([1.0, 1.5, 2.0, 2.5]),
                                    mu_p=0.3)["phi_c_E4"])


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


# =====================================================================================
# the quality gates -- defined AND applied (the 2026-09-09 defect)
# =====================================================================================
def test_every_defined_gate_is_APPLIED_not_merely_counted():
    """Until 2026-09-09 the gates were counted into `n_steady`/`n_quasistatic` and
    never used to exclude, so the primary estimator saw every cell whatever the flags
    said. One gate at a time, on a ladder built so that the flagged rung is the ONLY
    thing moving the root: a gate that is counted but not applied returns the UNGATED
    root and this fails.

    The second assertion is the discriminating half. A routine that dropped the rung
    but did not refit -- or one that dropped every rung -- would pass a bare count.
    """
    for gate in sorted(A.GATES):
        phi = np.array([0.60, 0.62, 0.64, 0.66, 0.68])
        p = 0.5 * (phi - 0.59)                  # four honest rungs, root at 0.59
        # one corrupt rung, carrying 5x too LITTLE stress -- which is the direction the
        # real non-quasi-static bottom rungs sit in (cv = 0 rung 0 reads p = 3.4e-06
        # against ~2e-04 from the rungs above it)
        p[0] *= 0.2
        ivals = np.full_like(phi, 1e-3)
        z = np.linspace(4.5, 6.5, 5)
        recs = _recs_from(phi, p, ivals, z)
        recs[0][gate] = False

        gated = A.analyse_group(recs, (gate,))
        ungated = A.analyse_group(recs, ())
        assert gated["n_cells"] == 4, gate
        assert gated["n_dropped_by_gate"] == 1, gate
        assert ungated["n_cells"] == 5 and ungated["n_dropped_by_gate"] == 0, gate
        assert gated["phi_c_E2"] == pytest.approx(0.59, abs=3e-3), gate
        assert abs(ungated["phi_c_E2"] - gated["phi_c_E2"]) > 5e-3, (
            f"gate {gate!r} did not move the primary estimator, so it was not applied")


def test_GATE_SETS_reports_every_defined_gate_AND_the_ungated_number():
    """A gate nobody prints is a gate nobody sees, and the ungated number must always
    be printed beside the gated ones -- the whole point of the 2026-09-09 fix is that
    neither is quoted alone."""
    covered = {g for _label, gs in A.GATE_SETS for g in gs}
    assert covered == set(A.GATES), (
        f"gates defined but never reported: {set(A.GATES) - covered}")
    assert ("ungated", ()) in A.GATE_SETS
    assert ("both", tuple(sorted(A.GATES))) in A.GATE_SETS


def test_apply_gates_refuses_an_unknown_gate_name():
    """An unknown name must raise, not silently admit everything."""
    with pytest.raises(KeyError):
        A.apply_gates([], ("no_such_gate",))


def test_a_gate_that_empties_the_ladder_returns_nan_not_a_number():
    phi = np.linspace(0.60, 0.68, 5)
    recs = _recs_from(phi, 0.5 * (phi - 0.59), np.full_like(phi, 1e-3),
                      np.linspace(4.5, 6.5, 5), steady=False)
    res = A.analyse_group(recs, ("steady",))
    assert res["n_cells"] == 0 and res["n_dropped_by_gate"] == 5
    assert np.isnan(res["phi_c_E2"])


def test_the_prod2_shear_rate_is_NOT_matched_across_the_CV_ladder():
    """gdot = gamma_max/(shear_steps*dt) and dt is set by the SMALLEST particle, so
    widening the size distribution RAISES the shear rate. Measured on the committed
    prod2 records; pinned here so that a future protocol change which fixes it (or
    makes it worse) cannot pass unnoticed."""
    import glob
    import json
    paths = sorted(glob.glob(os.path.join(HERE, "out", "phase8_prod2_shard*.jsonl")))
    assert paths, "the committed prod2 records are missing"
    by_cv = {}
    for path in paths:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                if "error" in r or not np.isfinite(r.get("gdot", np.nan)):
                    continue
                by_cv.setdefault(round(float(r["cv"]), 3), []).append(float(r["gdot"]))
    assert set(by_cv) == {0.0, 0.2, 0.5}
    means = {cv: float(np.mean(v)) for cv, v in by_cv.items()}
    span = max(means.values()) / min(means.values())
    assert span == pytest.approx(4.78, abs=0.05), means
    assert span > A.GDOT_MATCH_TOL, (
        "the rate confound is gone -- if that is real, update this test and the "
        "analyse_phase8 docstring together")
    # within one CV level the seeds span at most 1.10 (dt is set by the smallest
    # particle DRAWN, which is seed-dependent), so the 4.78x across CV is a systematic
    # several times larger than the seed scatter, not an instance of it
    for cv, v in by_cv.items():
        assert max(v) / min(v) < 1.10, (cv, v)


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
