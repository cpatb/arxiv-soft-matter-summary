"""Tests for item 14: uncertainty on the shipped injectability forces.

The gate is pre-registered in `case_studies/blood_injectability/uq_prereg.md`. These
tests pin the verdicts and, more importantly, the ONE THAT FAILED -- G14.5 -- so that a
later run cannot turn a recorded failure into a pass by moving a bar.

They also pin the thing that made this item worth doing: `results.md` and `uq.md` must
not disagree with each other. Two sections of one page contradicting 20 lines apart is
this repo's "prose contradicts its own artifact" failure mode, and it happened here --
the crosscheck section's ordinal sentence was a pairwise claim that item 14 refuted.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

CASE = Path(__file__).resolve().parent.parent / "case_studies" / "blood_injectability"
UQ_JSON = CASE / "uq.json"
UQ_MD = CASE / "uq.md"
RESULTS_MD = CASE / "results.md"
PREREG = CASE / "uq_prereg.md"


@pytest.fixture(scope="module")
def uq():
    assert UQ_JSON.exists(), "run case_studies/blood_injectability/uq.py"
    return json.loads(UQ_JSON.read_text(encoding="utf-8"))


def test_the_four_passing_gates_pass(uq):
    v = uq["verdict"]
    for gid in ("G14.2", "G14.3", "G14.4", "G14.6"):
        assert v[gid]["pass"], f"{gid}: {v[gid]}"


def test_G14_5_is_recorded_as_FAILED_at_its_unchanged_bar(uq):
    """The gate this item existed to test, and it failed. Pinned so the record cannot
    drift: 36 of 190 pairs against a bar of 50 %, which was NOT moved.

    The pre-registration predicted a failure at 35-50 %. The measured 18.9 % is worse
    than that, so the prediction was right in direction and wrong in size -- which is
    also pinned, because it is the kind of thing that gets quietly reread later as
    having been expected."""
    g = uq["verdict"]["G14.5"]
    assert g["pass"] is False
    assert g["bar"] == 0.50
    assert g["n_pairs"] == 190
    assert g["n_disjoint_hull"] == 36
    assert g["rate_hull"] == pytest.approx(36 / 190)
    assert g["min_rank_gap_always_separated"] == 17


def test_the_geometry_is_asserted_equal_not_toleranced(uq):
    """G14.3 compares four constants that `run.py` and `crosscheck.py` each declare as
    their own literal. Equality is EXACT: two pieces of code constructing the same
    quantity by different routes is a failure mode this repo has already paid for."""
    g = uq["verdict"]["G14.3"]
    assert g["mismatches"] == {}
    for _, (a, b) in g["geometry"].items():
        assert a == b


def test_the_borrowed_intervals_belong_to_the_shipped_forces(uq):
    """G14.2. The intervals come from a refit in `profile_interval_g139.json`; the
    forces they are attached to come from the shipped records. The residual is the
    record's `%.6g` parameter storage, i.e. a part in 1e-6, and nothing else."""
    g = uq["verdict"]["G14.2"]
    assert g["n"] == 20
    assert g["worst_force_rel_dev"] < 1e-5
    assert g["worst_force_rel_dev"] > 1e-8, (
        "a deviation of exactly zero would mean the two are the same calculation, "
        "which would make this gate vacuous rather than passing")


def test_the_yield_stress_is_insensitive_by_a_third_independent_route(uq):
    """G14.6. `results.md` argues the yield stress is invisible from `x = tau_y/tau_w`
    and from zeroing `tau_y`. This is a finite-difference `d ln F / d ln p` on the force
    solver, which shares neither construction, and it agrees on 20/20."""
    g = uq["verdict"]["G14.6"]
    assert g["n_tau_y_insensitive"] == 20
    assert g["worst_abs_sens_tau_y"] < g["INSENSITIVE"]
    # and the discriminating half: the force IS sensitive to eta_c, near-proportionally.
    # A route that called everything insensitive would pass the line above and mean
    # nothing.
    assert g["worst_abs_sens_eta_c"] > 0.9


def test_every_interval_belongs_to_a_closed_profile(uq):
    """G14.4. An end reported at a box bound is not a limit, and item 13 made MORE ends
    open, not fewer. Printing a number for one would quote a limit the data does not
    set."""
    assert uq["verdict"]["G14.4"]["records_with_an_open_end"] == []
    for r in uq["rows"]:
        assert r["open_parameters"] == []
        assert r["hull_lo"] < r["force_N"] < r["hull_hi"]


def test_results_md_does_not_contradict_the_uq_section():
    """The reason this test exists: the crosscheck section said *which patient the blood
    came from still matters more*, a PAIRWISE claim, while the UQ section 20 lines above
    said donors cannot be ranked pairwise. One page, two sections, opposite claims."""
    md = RESULTS_MD.read_text(encoding="utf-8")
    assert "cannot be ranked pairwise at 95 %" in md
    assert "The PAIRWISE reading of that does NOT survive" in md
    # the retracted phrasing may appear ONLY as a quoted retraction, never as a claim
    idx = md.find("which patient the blood came from still matters more")
    if idx != -1:
        assert "an earlier version of this sentence said" in md[max(0, idx - 300):idx]


def test_the_report_names_which_interval_object_it_is_quoting():
    """A push-through hull and a covariance sample are different objects and the page
    must never print one without saying which. Checked on both artifacts."""
    for path in (UQ_MD, RESULTS_MD):
        md = path.read_text(encoding="utf-8")
        assert "hull" in md and "sampled" in md.lower()
        assert "WIDER than the joint region" in md


def test_the_prereg_predicted_this_failure_before_it_was_run():
    """A failure predicted in advance and a failure explained afterwards are different
    things. The pre-registration is committed separately (3d267db) and says so."""
    p = PREREG.read_text(encoding="utf-8")
    assert "**G14.5 FAILS.**" in p
    assert "35–50 %" in p or "35-50 %" in p


def test_the_family_term_is_not_silently_added_into_the_intervals(uq):
    """The Casson-vs-Cross difference is a separate term. Adding it in would make the
    interval look like a total uncertainty when it is a parameter uncertainty for one
    family, and the page says so."""
    md = UQ_MD.read_text(encoding="utf-8")
    assert "not added in" in md
    assert "CASSON interval" in md
