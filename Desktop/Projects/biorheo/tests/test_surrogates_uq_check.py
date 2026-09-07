"""Tests for item 15: how far the item-13 interval correction reached the Carbopol study.

The result worth pinning is not a verdict, it is a MECHANISM. Item 13's step floor is
`box_bound/5000`, so the size of the defect scales with how wide a parameter's declared
box is. Herschel-Bulkley's three boxes differ by six orders of magnitude, and the
parameter the force depends on is the one with the narrow box -- so the correction moved
the parameters this force does not use.

That is the same shape as item 13's own G13.9 failure one item earlier (there it was the
wrong END of an interval, here the wrong PARAMETER), and neither was predicted. If a
future change makes `n`'s interval move, these assertions fail and somebody has to look.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from dst_fit.models import Model

from biorheo.constitutive import get_model

CASE = (Path(__file__).resolve().parent.parent / "case_studies" / "injectability_surrogates")
CHECK = CASE / "uq_check.json"
RESULTS = CASE / "results.md"
PREREG = CASE / "uq_prereg.md"


@pytest.fixture(scope="module")
def chk():
    assert CHECK.exists(), "run case_studies/injectability_surrogates/uq_check.py"
    return json.loads(CHECK.read_text(encoding="utf-8"))


def test_the_three_gates_pass(chk):
    v = chk["verdict"]
    for gid in ("G15.1", "G15.2", "G15.3"):
        assert v[gid]["pass"], f"{gid}: {v[gid]}"


def test_the_parameter_intervals_moved_but_the_force_did_not(chk):
    """The finding. G15.1 asks whether the correction reached this page; it did, at the
    parameter level, and did NOT at the level of the number the page ships."""
    worst_param = chk["verdict"]["G15.1"]["worst_interval_rel_move"]
    worst_hull = max(g["hull_rel_move"] for g in chk["gels"])
    assert worst_param > 0.01, worst_param            # > 1 %: the bar
    assert worst_hull < 1e-3, worst_hull              # < 0.1 %: the force is untouched
    assert worst_param / worst_hull > 100


def test_the_mechanism_is_the_box_width_not_the_fix(chk):
    """`tau_y` (box 1e5) and `k` (box 1e6) move; `n` (box O(1)) does not, by three orders
    of magnitude. This is asserted against the model's DECLARED bounds rather than
    against remembered numbers, so it stays true if the bounds are ever changed -- and
    fails loudly if they are changed in a way that breaks the argument."""
    hb: Model = get_model("herschel_bulkley")
    hi = dict(zip(hb.param_names, hb.bounds[1]))
    assert hi["n"] < 100 < hi["tau_y"] <= hi["k"], hi

    moves = {}
    for p in chk["parameters"]:
        moves.setdefault(p["param"], []).append(max(p["rel_move_lo"], p["rel_move_hi"]))
    assert max(moves["n"]) < 1e-4, moves["n"]
    assert max(moves["k"]) > 1e-2, moves["k"]
    assert max(moves["n"]) * 100 < max(moves["k"])


def test_the_force_is_set_by_the_parameter_the_bug_could_not_touch(chk):
    """Why the two facts above compose into "the forces did not move"."""
    for g in chk["gels"]:
        assert g["widest_parameter"] == "n", g


def test_both_gels_narrowed(chk):
    assert all(g["narrowed"] for g in chk["gels"])


def test_no_open_ends_so_that_branch_is_still_untested(chk):
    """Recorded as UNTESTED rather than counted as verified. The sibling study's
    "an open interval is never printed as a number" branch has still never executed, in
    either case study, so its body is unexercised code."""
    assert chk["verdict"]["G15.3"]["n_open_ends"] == 0
    for p in chk["parameters"]:
        assert p["new_at_bound"] == [False, False], p


def test_the_quoted_interval_is_compared_with_the_models_own_misfit(chk):
    """G15.4. `rmse_log` is a relative viscosity error and `d ln F / d ln k` is ~1, so
    the two are comparable as percentages. Neither gel flipped side, and on the 0.1 wt%
    gel the ratio is close enough to 1 that the page has to say the interval is not a
    total error."""
    ratios = {g["record"]: g["ratio_new"] for g in chk["gels"]}
    assert chk["verdict"]["G15.4"]["flipped"] == []
    assert all(r > 1.0 for r in ratios.values()), ratios
    assert min(ratios.values()) < 3.0, ratios


def test_every_corrected_endpoint_is_an_actual_crossing(chk):
    """The item-13 self-check, carried through to this study's own parameters."""
    for p in chk["parameters"]:
        assert p["crossing_error"] <= 0.01, p
        # The invariant that matters is NOT that every step came from the Wald error --
        # the step is `max(wald/4, |best| * 1e-3)` and the floor legitimately wins on 2 of
        # these 6. It is that the BOX never entered, which is the direct statement that
        # item 13's defect could not act on this study.
        assert p["step_rule"] in ("wald", "value"), p
        assert p["step_rule"] != "box_fallback", p
    rules = {p["step_rule"] for p in chk["parameters"]}
    assert rules == {"wald", "value"}, (
        "if this ever narrows to one rule, the sentence in uq_check.md counting both "
        f"is stale: {rules}")


def test_results_md_carries_the_provenance_and_does_not_overstate_it():
    md = RESULTS.read_text(encoding="utf-8")
    assert "Which routine drew these intervals" in md
    assert "The correction moved the parameters this force does not depend on" in md
    assert "PARAMETER interval, not a total error" in md


def test_the_prereg_recorded_that_the_items_premise_was_false():
    """Item 15 said this study "ships point forces too". It never did. Pinned because a
    false premise that gets quietly worked around is the failure this repo keeps
    hitting -- four times now."""
    p = PREREG.read_text(encoding="utf-8")
    assert "premise is FALSE" in p
    assert "It does not ship point forces" in p


def test_the_prereg_named_its_unfalsifiable_gates_in_advance():
    """Item 14's lesson, applied one item later: two of the obvious gates here could not
    fail, and saying so before the run is what stops them being quoted as evidence."""
    p = PREREG.read_text(encoding="utf-8")
    assert "CANNOT fail" in p
    assert "is not counted as a gate" in p
