"""test_strainclock.py — the strain-clock grid and the gamma_c estimator.

The estimator is graded on series whose answer is known in closed form BEFORE the
campaign's own numbers are read (the 2026-09-08 lesson: a diagnostic must reproduce
something already known before its new numbers are trusted). The non-monotone cases are
the point: this project has twice quoted a crossing time whose monotonicity was assumed.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from grade_strainclock import clock
from phase10_strainclock import (
    PHI_0, PHI_MU, SHEAR_STEPS, X_RUNGS, build_grid, rungs_for,
)


# --------------------------------------------------------------------------- the grid
def test_rungs_bracket_the_critical_line():
    """One rung below phi_mu, four above, and every rung below phi_0."""
    for cv in PHI_MU:
        r = rungs_for(cv)
        assert len(r) == len(X_RUNGS)
        assert sum(v < PHI_MU[cv] for v in r) == 1, f"cv {cv}: {r}"
        assert all(v < PHI_0[cv] for v in r), f"cv {cv} has a rung at/above phi_0: {r}"
        assert r == sorted(r)


def test_offsets_are_matched_across_cv_to_11_percent():
    """Both collapse coordinates (x and phi - phi_mu) live on this one grid, which is
    only true because W moves little; if a future anchor update breaks that, the
    pre-registration's 'both will be reported' has to be rewritten, not this test."""
    offs = {cv: [round(v - PHI_MU[cv], 6) for v in rungs_for(cv)] for cv in PHI_MU}
    for i in range(len(X_RUNGS)):
        vals = [abs(offs[cv][i]) for cv in PHI_MU]
        assert max(vals) / min(vals) < 1.15, f"rung {i}: {offs}"


def test_grid_is_45_cells_and_shards_mix_the_cvs():
    grid = build_grid(sorted(PHI_MU), (1, 2, 3), 1000)
    assert len(grid) == 45
    for shard in range(15):
        mine = [g for i, g in enumerate(grid) if i % 15 == shard]
        assert len(mine) == 3
        assert len({g[0] for g in mine}) == 3, f"shard {shard} drew {mine}"


def test_step_counts_match_the_rate_matched_ladder():
    assert SHEAR_STEPS == {0.0: 102_000, 0.2: 200_000, 0.5: 480_000}


# ---------------------------------------------------------------- the clock estimator
def test_clock_on_an_exponential_rise_recovers_the_analytic_half_strain():
    """p = p_ss (1 - exp(-gamma/tau)) reaches half at gamma = tau ln 2."""
    tau = 0.05
    g = np.linspace(0, 3.0, 1501)
    p_ss = 3e-3
    p = p_ss * (1 - np.exp(-g / tau))
    first, last, n_down, span = clock(g, p, 0.0, p_ss, 0.5)
    assert first == pytest.approx(tau * math.log(2), rel=2e-3)
    assert last == pytest.approx(first)
    assert n_down == 0 and span == pytest.approx(0.0, abs=1e-12)


def test_clock_is_offset_by_p_prep_not_by_the_raw_level():
    """The target is p_prep + level*(p_ss - p_prep): a pack that starts with pressure
    must not be credited with having already run half its clock."""
    g = np.linspace(0, 1.0, 1001)
    p = 1e-3 + 3e-3 * g                      # p_prep = 1e-3, p_ss = 4e-3 at gamma = 1
    first, _, _, _ = clock(g, p, 1e-3, 4e-3, 0.5)
    assert first == pytest.approx(0.5, rel=1e-3)


def test_clock_reports_a_non_monotone_series_rather_than_hiding_it():
    """A series that crosses, falls back, and crosses again: first and last must
    DIFFER, which is exactly what the readability flag keys on."""
    g = np.linspace(0, 1.0, 1001)
    p_ss = 2e-3
    p = np.where(g < 0.1, 2e-2 * g,                       # up through half at 0.05
                 np.where(g < 0.3, 2e-3 - 8e-3 * (g - 0.1),   # back down below it
                          8e-4 + 3.3e-3 * (g - 0.3)))         # and up again
    first, last, n_down, span = clock(g, p, 0.0, p_ss, 0.5)
    assert first < 0.1
    # the second up-crossing is 8e-4 + 3.3e-3 (g - 0.3) = 1e-3, i.e. g = 0.3606
    assert last == pytest.approx(0.3606, abs=2e-3)
    assert span > 0.25
    assert n_down > 0


def test_clock_returns_nan_when_the_pressure_never_builds():
    g = np.linspace(0, 3.0, 301)
    p = np.full_like(g, 1e-9)
    first, last, n_down, span = clock(g, p, 0.0, 3e-3, 0.5)
    assert math.isnan(first) and math.isnan(last)


def test_clock_at_90_percent_is_later_than_at_50_percent():
    tau = 0.08
    g = np.linspace(0, 3.0, 1501)
    p = 3e-3 * (1 - np.exp(-g / tau))
    c50 = clock(g, p, 0.0, 3e-3, 0.5)[0]
    c90 = clock(g, p, 0.0, 3e-3, 0.9)[0]
    assert c90 > c50
    assert c90 / c50 == pytest.approx(math.log(10) / math.log(2), rel=5e-3)


def test_sampling_resolves_the_shortest_clock_the_campaign_expects():
    """Delta gamma = 0.002 must put at least five samples under the shortest plausible
    clock (U4's fastest relaxation was gamma_half = 0.030)."""
    dgamma = 3.0 / 1500
    assert 0.030 / dgamma >= 5


# ------------------------------------------------- the landed campaign, pinned
# These assert the numbers that PREREG_STRAINCLOCK_2026-09-13.md §7 states in prose.
# A .md that quotes an artifact must not be able to drift away from it silently.
import glob      # noqa: E402
import json      # noqa: E402
import os        # noqa: E402

from grade_strainclock import cell_row, gm      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SHARDS = sorted(glob.glob(os.path.join(HERE, "out", "phase10_sc_shard*.jsonl")))


def _rows():
    seen, out = set(), []
    for f in SHARDS:
        for line in open(f, encoding="utf-8"):
            d = json.loads(line)
            if d["key"] in seen or d.get("error"):
                continue
            seen.add(d["key"])
            out.append(cell_row(d))
    return out


def test_the_campaign_is_complete_and_clean():
    assert len(SHARDS) == 15, "the shard files are part of the record"
    r = _rows()
    assert len(r) == 45
    assert all(x["protocol"] == "cv-shear-v3-fprep" for x in r)
    assert all(x["fswitch"] for x in r)


def test_g4_collapse_ratios_are_what_the_prereg_states():
    r = [x for x in _rows() if x["readable"] and x["z_nr_ss"] >= 4.0 and x["x"] > 0]
    tab = {}
    for cv in (0.0, 0.2, 0.5):
        tab[cv] = {x: gm([q["gamma_c"] for q in r
                          if q["cv"] == cv and abs(q["x"] - x) < 1e-9])
                   for x in (0.35, 0.60, 0.85)}
    for cv, expect in ((0.2, 1.069), (0.5, 0.914)):
        R = gm([tab[cv][x] / tab[0.0][x] for x in (0.35, 0.60, 0.85)])
        assert R == pytest.approx(expect, abs=5e-3), (cv, R)


def test_g5_crossings_recover_the_constant_pressure_phi_c():
    """The cross-protocol anchor, and the contrast that is the campaign's headline."""
    from phase10_strainclock import PHI_MU
    r = _rows()
    cross = {}
    for cv in (0.0, 0.2, 0.5):
        rr = sorted({q["phi_target"] for q in r if q["cv"] == cv})
        ps = [float(np.mean([q["p_ss"] for q in r
                             if q["cv"] == cv and q["phi_target"] == p])) for p in rr]
        for a, b, pa, pb in zip(rr, rr[1:], ps, ps[1:]):
            if (pa - 3e-4) * (pb - 3e-4) <= 0:
                cross[cv] = a + (3e-4 - pa) * (b - a) / (pb - pa)
                break
    assert len(cross) == 3
    for cv in cross:
        assert abs(cross[cv] - PHI_MU[cv]) < 0.010, cv
    assert cross[0.5] - cross[0.0] == pytest.approx(0.0523, abs=1e-3)
    assert abs((cross[0.5] - cross[0.0]) - (PHI_MU[0.5] - PHI_MU[0.0])) < 0.002


def test_the_frictionless_prep_carries_nothing_on_every_graded_rung():
    graded = [x for x in _rows() if x["x"] <= 0.60 + 1e-9]
    assert len(graded) == 36
    assert all(x["p_prep"] < 1e-6 for x in graded)
    assert min(x["p_prep"] for x in graded) < 1e-10
    assert max(x["ratio"] for x in _rows() if x["p_prep"] > 0) > 1e8


def test_the_clock_is_unreadable_only_near_the_critical_line():
    """The G2 failure is a physical result and this is the shape of it."""
    r = _rows()
    for x, want in ((-0.15, False), (0.10, False), (0.35, True), (0.60, True),
                    (0.85, True)):
        rr = [q for q in r if abs(q["x"] - x) < 1e-9]
        assert len(rr) == 9
        assert all(q["readable"] == want for q in rr), x
