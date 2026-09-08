"""System size must not be pooled into a ladder, a scatter estimate or a seed mean.

Found 2026-09-08: `analyse_phase2.py` de-duplicated on a record key that encodes N -- its
docstring said so -- while the three keys it GROUPED by did not. 18 of 95 ladders and 20
of 114 scatter cells were mixing N = 1500 and N = 3000, so a system-size difference was
entering as ladder scatter and as seed noise.

These tests assert the fix is in the source and that its two dependent gates were made
N-aware in the same change -- because adding N to the summary table immediately made
anchor A4 report FAIL on a sequence that is supposed to be monotone in mu_p, which is the
same failure mode one level down.
"""

from __future__ import annotations

import json
import os
import re

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
SRC = os.path.join(HERE, "analyse_phase2.py")
VERDICT = os.path.join(OUT, "phase2_verdict.txt")
DIAG = os.path.join(OUT, "n_pooling.json")


@pytest.fixture(scope="module")
def src() -> str:
    with open(SRC, encoding="utf-8") as f:
        return f.read()


@pytest.fixture(scope="module")
def verdict() -> str:
    with open(VERDICT, encoding="utf-8") as f:
        return f.read()


@pytest.fixture(scope="module")
def diag() -> dict:
    with open(DIAG, encoding="utf-8") as f:
        return json.load(f)


# --- the fix is in the source, on all three keys --------------------------------------

def test_ladder_key_contains_N(src):
    m = re.search(r"groups\.setdefault\(\((.*?)\)", src, re.S)
    assert m, "the ladder grouping moved; re-point this test"
    assert 'r["n_requested"]' in m.group(1)


def test_seed_scatter_key_contains_N(src):
    m = re.search(r"scat\.setdefault\(\((.*?)\), \[\]\)", src, re.S)
    assert m, "the scatter grouping moved; re-point this test"
    assert 'r["n_requested"]' in m.group(1)


def test_seed_aggregate_key_contains_N(src):
    m = re.search(r"agg\.setdefault\(\((.*?)\),\s*\[\]\)", src, re.S)
    assert m, "the aggregate grouping moved; re-point this test"
    assert 'r["n_requested"]' in m.group(1)


# --- the two gates that the fix could have silently re-decided -------------------------

def test_anchor_A4_is_graded_per_system_size(verdict):
    """A4 asks whether phi_J falls monotonically in mu_p. With two N in the table its
    sequence gained a second mu_p = 0.3 entry and it reported FAIL -- on a different
    question. It must now be labelled with its N."""
    lines = [ln for ln in verdict.splitlines() if ln.startswith("ANCHOR A4")]
    assert lines, "A4 is missing from the verdict"
    for ln in lines:
        assert re.search(r"N = \d+", ln), f"A4 line without an N: {ln}"


def test_anchor_A4_passes_at_the_size_every_claim_is_quoted_at(verdict):
    ln = [x for x in verdict.splitlines()
          if x.startswith("ANCHOR A4") and "N = 1500" in x]
    assert len(ln) == 1, ln
    assert ln[0].rstrip().endswith("PASS")
    # and it must be a real sequence, not one point
    assert ln[0].count(":") >= 4


def test_every_A5_block_names_its_system_size(verdict):
    blocks = [ln for ln in verdict.splitlines() if ln.startswith("--- family")]
    assert blocks
    for ln in blocks:
        assert re.search(r"N = \d+", ln), f"A5 block without an N: {ln}"


def test_mismatched_CV_spans_are_flagged_not_silently_compared(verdict):
    """Two N blocks are only a finite-size comparison if they span the same CV range.
    The uniform family has N = 3000 only at cv 0.4-0.5 against 0.0-0.5 at N = 1500."""
    assert "CV SPAN DIFFERS BETWEEN SYSTEM SIZES" in verdict
    assert "NOT a matched finite-size comparison" in verdict


def test_the_lognormal_M2_verdict_agrees_across_a_MATCHED_pair(verdict):
    """The one family whose two sizes DO span the same CV range. Both must reach the
    same M2 verdict, and their effect at phi = 0.44 must agree to better than 10 %."""
    got = {}
    cur = None
    for ln in verdict.splitlines():
        m = re.match(r"--- family = lognormal, mu_p = 0\.30, N = (\d+) ---", ln)
        if m:
            cur = int(m.group(1))
        elif cur and "at the lab's phi = 0.44 the effect is" in ln:
            got[cur] = float(re.search(r"is ([\d.]+)x", ln).group(1))
            cur = None
    assert set(got) == {1500, 3000}, got
    assert abs(got[1500] - got[3000]) / got[1500] < 0.10, got
    assert verdict.count("M2 VERDICT: INSUFFICIENT") >= 2


# --- the diagnosis itself --------------------------------------------------------------

def test_the_defect_was_real_and_is_measured(diag):
    assert diag["n_ladders_mixing_N"] > 0
    assert diag["n_scatter_cells_mixing_N"] > 0
    assert diag["sizes_present"] == [1500, 3000]


def test_the_pooling_moved_phi_J_by_more_than_a_seed_scatter(diag):
    """If it had not, the fix would be housekeeping rather than a correction."""
    assert diag["d_phi_J_max_in_seed_sd"] > 1.0


def test_the_pooling_did_NOT_rescue_the_exponent_anchor(diag):
    """The useful outcome: a defect that COULD have explained the anchor's failure,
    ruled out by measurement. If this ever flips, neither number may be quoted."""
    assert diag["pooled"]["verdict"] == "FAIL"
    assert diag["split"]["verdict"] == "FAIL"
    lo, hi = diag["B_BAND"]
    assert (lo, hi) == (0.6, 1.6), "the pre-registered band was changed"
    assert diag["split"]["b_median"] < lo


def test_splitting_N_costs_ladders_and_that_is_reported(diag):
    """Separating N takes ladders from 8 points to 4, and the free fit needs 4, so the
    split DISQUALIFIES ladders as well as changing them. The honest trade must be visible."""
    assert diag["split"]["n_no_free_fit"] > diag["pooled"]["n_no_free_fit"]
    with open(os.path.join(OUT, "n_pooling.txt"), encoding="utf-8") as f:
        assert "DISQUALIFIES ladders" in f.read()


# --- the four siblings that had the same blindness -------------------------------------

CONSUMERS = ("analyse_skewness.py", "analyse_families.py",
             "diag_convergence.py", "diag_ladder.py",
             "analyse_exponent_constrained.py", "analyse_finite_size.py")


@pytest.mark.parametrize("fname", CONSUMERS)
def test_every_phase2_consumer_is_explicit_about_system_size(fname):
    """Grep for the PROPERTY a new N changes, not for the shape of the last defect.

    Four of these had zero references to N and their printed numbers moved the moment
    N = 3000 cells for two more families were ingested.
    """
    with open(os.path.join(HERE, fname), encoding="utf-8") as f:
        src = f.read()
    assert ("n_requested" in src) or ("PRIMARY_N" in src) or \
           ("restrict_to_primary_N" in src), \
        f"{fname} reads out/phase2*.jsonl without mentioning system size"


def test_the_restriction_helper_prints_what_it_drops():
    """A silent drop is how a later N sweep stops being included without anyone noticing."""
    with open(os.path.join(HERE, "analyse_phase2.py"), encoding="utf-8") as f:
        src = f.read()
    i = src.index("def restrict_to_primary_N")
    body = src[i:i + 1200]
    assert "print(" in body
    assert "dropped" in body and "sizes present" in body


def test_PRIMARY_N_is_the_size_the_claims_are_quoted_at():
    import importlib
    import sys as _s
    _s.path.insert(0, HERE)
    m = importlib.import_module("analyse_phase2")
    assert m.PRIMARY_N == 1500
