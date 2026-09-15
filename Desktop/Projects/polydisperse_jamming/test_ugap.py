"""test_ugap.py -- regression tests for the U-GAP block (PREREG_UGAP_2026-09-16.md).

Pins the DESIGN (which cannot change without changing the claim), the WIRING (the
branch algebra and the vectorised kernel every interval rests on), and the VERDICTS once
graded.

Scope guard, as `test_kladder.py` and `test_ternary.py` carry: this file must be
incapable of re-deciding a closed block. It reads `phase12_ug*` only, and asserts so
structurally rather than by comment.
"""

from __future__ import annotations

import csv
import fnmatch
import math
import os

import numpy as np
import pytest

import ugap as ug

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
GRADE = os.path.join(OUT, "grade_ugap.txt")
CELLS = os.path.join(OUT, "ugap_cells.csv")


# ------------------------------------------------------------------ scope / structure

def test_tag_cannot_collide_with_a_closed_block():
    for tag in ("ug", "ug3k"):
        for closed in ("mixr", "extr", "tern", "kl"):
            assert not tag.startswith(closed)
            assert not fnmatch.fnmatch(f"phase12_{closed}_shard00of12.jsonl",
                                       f"phase12_{tag}*_shard*.jsonl")


def test_grader_scope_assertion_runs_and_passes():
    import grade_ugap as gu
    msg = gu.assert_scope()
    assert "scope OK" in msg


def test_the_identifiability_bar_is_the_one_grade_mixture_ext_uses():
    """The bar is QUOTED, not chosen. If grade_mixture_ext ever moves it, this fails."""
    import inspect

    import grade_mixture_ext as gx
    src = inspect.getsource(gx.tabulated_closure)
    assert "0.25" in src, "the closed block's rule no longer mentions 0.25"
    assert ug.IDENT_WIDTH == 0.25


def test_bars_have_not_moved():
    import grade_ugap as gu
    assert gu.BAR_RMS == 2.02          # TAB in-sample rms
    assert gu.BAR_MAX == 5.74          # TAB in-sample max
    assert gu.U2_MIN_NEG == 7
    assert gu.U3_MIN_LEVELS == 2


# ------------------------------------------------------------------ the design

def test_the_ratio_is_the_one_K5_actually_reads():
    """Arm A's ratio must be K = 5's adjacent pair on the span-4 ladder, not 1.25."""
    import kladder as kl
    d = kl.ladder(5)
    adjacent = d[0] / d[1]
    assert math.isclose(adjacent, ug.U_STAR, rel_tol=0, abs_tol=1e-15)
    assert math.isclose(ug.RATIOS["A"], 4.0**0.25, rel_tol=0, abs_tol=1e-15)
    # and it really is the MOST NUMEROUS pair type at K = 5
    us = [round(d[i] / d[j], 12) for j in range(5) for i in range(j)]
    assert us.count(round(ug.U_STAR, 12)) == 4
    assert len(us) == 10


def test_every_cell_is_affordable_and_none_is_silently_dropped():
    cells = ug.cells()
    assert len(cells) == len(ug.RATIOS) * len(ug.X_GRID) == 27
    assert len(ug.admitted()) == 27, "a refused cell must be recorded, not dropped"
    for c in cells:
        assert c["n_min_at_judge"] >= ug.MIN_PHASE_COUNT


def test_the_grid_straddles_the_branch_crossover_on_every_arm():
    """The whole design rests on this: a is seen only by the coarse branch, b only by
    the fine one, so both need their own half of the composition curve."""
    beta = 0.5946865273380492
    # the shipped interaction values at the three arms' own u, rounded; recomputing them
    # here would need the whole closure and this test is meant to be cheap.
    for arm, (a, b) in (("A", (0.9398, 0.9177)), ("B", (0.9695, 0.9574)),
                        ("C", (0.9249, 0.8985))):
        xc = ug.crossover(beta, a, b)
        below = sum(1 for x in ug.X_GRID if x < xc)
        above = len(ug.X_GRID) - below
        assert below >= 3, f"arm {arm}: only {below} coarse-branch points (x* {xc:.3f})"
        assert above >= 3, f"arm {arm}: only {above} fine-branch points (x* {xc:.3f})"


def test_the_bitid_cells_are_in_the_grid():
    for r, x in ug.BITID:
        assert x in ug.X_GRID
        assert r in ug.RATIOS.values()


# ------------------------------------------------------------------ wiring

def test_branch_algebra_matches_the_shipped_closure():
    """`branch_tots` is hand-derived from packing_density_ab; it is an ANCHOR only if it
    actually reproduces it. Checked over random (r, a, b, x), not asserted."""
    import mixture_closure as mc
    beta = 0.5946865273380492
    rng = np.random.default_rng(7)
    worst = 0.0
    for _ in range(200):
        a, b, x = rng.uniform(1e-3, 1 - 1e-3, 3)
        r = float(rng.uniform(1.05, 6.0))
        tc, tf = ug.branch_tots(float(x), beta, float(a), float(b))
        mine = beta / (1.0 - min(tc, tf))
        ref = mc.packing_density_ab(np.array([r, 1.0]), np.array([1 - x, x]),
                                    beta, float(a), float(b))
        worst = max(worst, abs(mine - ref))
    assert worst < 1e-12, f"branch algebra disagrees by {worst:.3e}"


def test_vectorised_kernel_reproduces_the_shipped_one():
    beta = 0.5946865273380492
    assert ug.check_kernel(beta, np.random.default_rng(20260916)) < 1e-12


def test_profile_fast_reproduces_profile_interval_exactly():
    """A faster reimplementation of the thing under test is exactly the kind of code
    that reports a bug as a finding. This pins it to the shipped routine."""
    import mixture_closure as mc
    beta = 0.5946865273380492
    pts = [(1.5, x, float(mc.packing_density_ab(np.array([1.5, 1.0]),
                                                np.array([1 - x, x]), beta, 0.92, 0.90)))
           for x in (0.05, 0.12, 0.20, 0.28, 0.36, 0.45)]
    for which in (0, 1):
        fast = ug.profile_fast(pts, beta, which)
        ship = mc.profile_interval(pts, beta, which)
        assert fast == pytest.approx(ship, abs=1e-12)


def test_the_profile_bar_is_per_point_so_duplicating_points_changes_nothing():
    """The design's central claim: more compositions at the same sensitivity buy nothing
    under this rule. If this ever fails, the whole 'spend on the grid, not on seeds'
    argument in PREREG_UGAP section 1 is wrong."""
    import mixture_closure as mc
    beta = 0.5946865273380492
    pts = [(1.5, x, float(mc.packing_density_ab(np.array([1.5, 1.0]),
                                                np.array([1 - x, x]), beta, 0.92, 0.90)))
           for x in (0.05, 0.12, 0.20, 0.28, 0.36, 0.45, 0.55, 0.70)]
    w1 = ug.profile_fast(pts, beta, 0)
    w2 = ug.profile_fast(pts * 2, beta, 0)
    w4 = ug.profile_fast(pts * 4, beta, 0)
    assert w1 == w2 == w4


def test_the_escape_model_really_is_the_constant_beta():
    """a = 1 and b = 1 must BOTH collapse the closure to phi = beta, for every x.
    This is the mechanism C4 rests on, so it is run rather than asserted."""
    beta = 0.5946865273380492
    xs = np.linspace(0.01, 0.99, 50)
    for a, b in ((1.0, 0.3), (0.3, 1.0), (1.0, 1.0)):
        got = ug.phi_fast(a, b, xs, beta)
        assert np.allclose(got, beta, atol=1e-15)


def test_identifiable_rejects_a_degenerate_case_and_accepts_a_good_one():
    """A DISCRIMINATING negative: a rule that said YES to everything would pass a
    one-sided test. Flat data (phi = beta everywhere) must be UNDETERMINED."""
    import mixture_closure as mc
    beta = 0.5946865273380492
    flat = [(1.5, x, beta) for x in ug.X_GRID]
    assert not ug.identifiable(flat, beta)[0]
    good = [(1.5, x, float(mc.packing_density_ab(np.array([1.5, 1.0]),
                                                 np.array([1 - x, x]), beta, 0.92, 0.90)))
            for x in ug.X_GRID]
    assert ug.identifiable(good, beta)[0]


# ------------------------------------------------------------------ the verdicts

@pytest.mark.skipif(not os.path.exists(GRADE), reason="block not graded yet")
def test_graded_verdicts_are_pinned():
    txt = open(GRADE, encoding="utf-8").read()
    # U0 -- bit identity, exactly zero
    assert "12 of 12 bit-identical" in txt
    assert "DIFFERS" not in txt
    # U1 -- the headline: a(u*) IS measurable
    assert "U1 (arm A, u = 0.70711): IDENTIFIABLE" in txt
    # U1c -- and u = 0.8 is not
    assert "U1c (arm B, u = 0.80000): UNDETERMINED" in txt
    # U1b -- the attribution control held
    assert "U1b PASS" in txt
    # U4 / U5
    assert "U4 monotone: a True, b True -- PASS" in txt
    assert "U5 PASS" in txt
    # the registered failures, both of them
    assert "'U1c: arm B came back UNDETERMINED'" in txt
    assert "U3: 1/3 and 0/3 levels <= 0" in txt


@pytest.mark.skipif(not os.path.exists(GRADE), reason="block not graded yet")
def test_C6_keeps_the_gate_and_the_science_as_SEPARATE_sentences():
    """The point of C6. A future edit that collapses 'the gate failed' into 'saturation
    is refuted' must break a test, because the corner sweep says UNDETERMINED."""
    txt = open(GRADE, encoding="utf-8").read()
    assert "NOT RESOLVED" in txt
    assert "saturation SURVIVES on 110 of 169 corners" in txt
    claims = open(os.path.join(HERE, "CLAIMS_UGAP_2026-09-16.md"), encoding="utf-8").read()
    assert "UNDETERMINED, NOT CONFIRMED" in claims
    assert "110 of 169 corners" in claims
    assert "not that it is" in claims and "refuted" in claims


@pytest.mark.skipif(not os.path.exists(GRADE), reason="block not graded yet")
def test_the_blind_caveats_reached_the_claims_file():
    """PREREG section 6 and lessons.md 2026-09-15: a caveat printed in the BLIND file
    must be carried into the claims. The K-ladder shipped without one until a referee
    asked. Every caveat's leading tag must appear in the claims."""
    import predict_ugap as pu
    claims = open(os.path.join(HERE, "CLAIMS_UGAP_2026-09-16.md"), encoding="utf-8").read()
    tags = sorted({c.split()[0] for c in pu.CAVEATS if c.split()[0].startswith("C-")})
    assert tags == ["C-i", "C-ii", "C-iii", "C-iv"]
    for t in tags:
        assert f"**{t}**" in claims, f"caveat {t} did not reach the claims file"


@pytest.mark.skipif(not os.path.exists(CELLS), reason="block not graded yet")
def test_cells_artifact_is_complete_and_matches_the_design():
    with open(CELLS, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 27, f"{len(rows)} rows, expected 27"
    assert {r["arm"] for r in rows} == {"A", "B", "C"}
    for r in rows:
        assert int(r["n_seeds"]) == len(ug.SEEDS)
        assert 0.55 < float(r["phi"]) < 0.68


@pytest.mark.skipif(not os.path.exists(GRADE), reason="block not graded yet")
def test_the_arm_B_bound_is_labelled_as_an_inference_from_a_null():
    """C4 must never read as a measurement. If someone deletes the label, this fails."""
    txt = open(GRADE, encoding="utf-8").read()
    assert "inference FROM A NULL" in txt
    assert "not a measurement" in txt
    claims = open(os.path.join(HERE, "CLAIMS_UGAP_2026-09-16.md"), encoding="utf-8").read()
    assert "inference FROM A NULL" in claims
