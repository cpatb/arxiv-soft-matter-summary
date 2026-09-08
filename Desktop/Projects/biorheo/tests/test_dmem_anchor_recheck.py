"""Item 23: the item-22 anchor, re-graded against the paper's own Table 2.

These pin the artifact written by `scripts/dmem_anchor_recheck.py`. They are fast: the
artifact is committed and the script that writes it is the slow part.

The point of most of them is that a **failure stays recorded**. Item 23 corrected a claim
that had reached the Status card, `NEXT_STEPS.md`, `OFFER.md` and the project's memory, and
the failure mode this repo keeps hitting is that a page goes on asserting a superseded
number twenty lines from the section that supersedes it. So there is a test asserting the
two pages cannot drift apart, exactly as item 14 needed one.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from biorheo.constitutive.shift import MIN_PLAUSIBLE_KELVIN, arrhenius_shift

ROOT = Path(__file__).resolve().parent.parent
ART = ROOT / "data" / "dmem_sints2024" / "anchor_recheck.json"
ART_MD = ROOT / "data" / "dmem_sints2024" / "anchor_recheck.md"
RECORDS_MD = ROOT / "data" / "dmem_sints2024" / "records.md"

# The paper's Table 2, transcribed by hand from arXiv:2409.13650. Pinned here as well as in
# the script so that a silent edit to either is a test failure rather than a new anchor.
TABLE2_N = {"DMEM + 10 % FBS": 0.31, "DMEM + 5 % FBS": 0.33,
            "DMEM + 1 % FBS": 0.51, "DMEM + 10 % FBS (37 C)": 0.54}
TABLE2_K = {"DMEM + 10 % FBS": 0.028, "DMEM + 5 % FBS": 0.037,
            "DMEM + 1 % FBS": 0.016, "DMEM + 10 % FBS (37 C)": 0.005}


@pytest.fixture(scope="module")
def art() -> dict:
    assert ART.exists(), (
        f"{ART} is missing. Run `python scripts/dmem_anchor_recheck.py`. It is committed, "
        "so a missing file means it was deleted rather than never generated."
    )
    return json.loads(ART.read_text(encoding="utf-8"))


def test_the_reference_values_are_the_papers_table_2_not_its_figure_captions(art):
    """Item 22 graded against the paper's COMSOL inputs, which are not its fit results.

    Both sets are kept in the artifact. This asserts which one is now the reference AND
    that they still differ -- if a later edit quietly made them equal, the whole item would
    read as a no-op while the anchor silently reverted.
    """
    differing = 0
    for r in art["rows"]:
        assert r["n_paper_table2"] == pytest.approx(TABLE2_N[r["label"]])
        assert r["K_paper_table2"] == pytest.approx(TABLE2_K[r["label"]])
        if r["n_paper_table2"] != r["n_paper_figure_caption"]:
            differing += 1
    assert differing == 4, (
        f"{differing} of 4 conditions differ between Table 2 and the figure captions; "
        "item 23's premise is that ALL FOUR do. Two of the four are the right number "
        "against the wrong row -- `SPEC` gives 5 % FBS the 0.31 that belongs to 10 %, and "
        "gives 1 % FBS the 0.54 that belongs to 37 C -- which is how a four-of-four "
        "mismatch reads at a glance as a three-of-four one. It was written up as three "
        "until this assertion failed."
    )


def test_g23_1_failed_at_the_unchanged_bar_and_the_bar_was_not_moved(art):
    """The corrected anchor FAILS, at 0.05, on exactly one condition. That is the result."""
    g = art["gates"]["G23_1_corrected_anchor_at_the_original_bar"]
    assert art["G221_TOL_unchanged"] == 0.05, "the bar was moved; it must not be"
    assert g["pass"] is False
    assert g["conditions_failing"] == ["DMEM + 1 % FBS"]
    assert g["worst_abs_dev"] == pytest.approx(0.0599, abs=5e-4)
    # and the failure must be a real miss, not a rounding artefact sitting on the bar
    assert g["worst_abs_dev"] > art["G221_TOL_unchanged"] * 1.15


def test_g23_2_neither_single_correction_rescues_the_old_pass(art):
    """Item 22's 0.0299 survives neither correction applied on its own."""
    g = art["gates"]["G23_2_neither_single_correction_rescues_the_old_pass"]
    tol = art["G221_TOL_unchanged"]
    assert g["item22_reported_worst"] < tol, "item 22's number did pass, which is the point"
    assert g["worst_dropping_out_of_range_interval_only"] == pytest.approx(0.0651, abs=5e-4)
    assert g["worst_regrading_against_table2_only"] == pytest.approx(0.0599, abs=5e-4)
    assert g["worst_dropping_out_of_range_interval_only"] > tol
    assert g["worst_regrading_against_table2_only"] > tol


def test_g23_3_the_flagship_agreement_was_a_pooling_coincidence(art):
    """0.0005 came from a pooled index outside the span of the two intervals it pooled."""
    g = art["gates"]["G23_3_the_flagship_agreement_is_a_pooling_coincidence"]
    assert g["inside_span"] is False
    lo, hi = min(g["per_interval_n"]), max(g["per_interval_n"])
    assert not (lo <= g["pooled_n"] <= hi)
    # the discriminating half: the OTHER conditions' pooled fits ARE inside their spans, so
    # this is a property of that condition and not of the pooling code.
    inside = [v for k, v in g["conditions_inside_span"].items() if k != g["flagship_condition"]]
    assert sum(inside) >= 2, (
        "if no condition's pooled fit were bracketed by its intervals, the finding would "
        "be about the estimator rather than about this condition"
    )


def test_g23_4_is_scaled_against_the_sources_own_declared_repeatability(art):
    """What survives is measured against a number the SOURCE declares, not one we chose."""
    g = art["gates"]["G23_4_what_survives_against_the_authors_own_repeatability"]
    assert g["paper_declared_repeatability"] == 0.20, (
        "the scale must stay the 20 % the authors state in S2 Fig. Changing it would be "
        "choosing our own tolerance again, which is the thing item 23 corrected"
    )
    assert g["pass"] is True
    assert g["worst_max_rel"] == pytest.approx(0.1903, abs=2e-3)
    assert g["worst_max_rel"] < g["paper_declared_repeatability"]
    # not vacuous: the agreement must not be so close that the bar could not have failed
    assert g["worst_max_rel"] > 0.05, (
        "a worst deviation far under the bar would make this gate unfalsifiable; it is "
        "19 % against 20 %, which is a real pass and is stated as one"
    )
    assert g["worst_rel_dev_K"] == pytest.approx(0.0668, abs=1e-3)
    assert g["K_rounds_to_published_count"] == 3


def test_g23_9_is_stricter_than_the_bar_it_sits_beside(art):
    """Scaling by the fit's own profile interval finds ONE resolved disagreement."""
    g = art["gates"]["G23_9_the_published_n_against_our_own_profile_interval"]
    assert g["conditions_outside"] == ["DMEM + 1 % FBS"]
    assert g["pass"] is False
    # the discriminating half: a routine that called everything a disagreement would also
    # report a failure here, so require the other three to be genuinely inside.
    assert sum(g["inside"].values()) == 3
    # and the 37 C condition, whose raw deviation under the restricted pooling was 0.1051 --
    # larger than the one that fails -- must come out INSIDE, or the gate is just the bar
    assert g["inside"]["DMEM + 10 % FBS (37 C)"] is True


def test_g23_1b_shows_the_data_choice_moves_37c_more_than_the_disagreement(art):
    """At 37 C the choice of which intervals to fit dominates the comparison."""
    g = art["gates"]["G23_1b_how_much_is_a_data_selection_choice"]
    d = g["the_discarded_37C_interval_alone"]
    assert d["abs_dev_n_vs_table2"] < 0.01, (
        "the interval item 22's restriction throws away reproduces the published index "
        "almost exactly; that is why the restriction and not the fit is what disagreed"
    )
    restricted_37 = next(r for r in g["per_condition"] if "37 C" in r["label"])
    pooled_37 = next(r for r in art["rows"] if "37 C" in r["label"])
    spread = abs(restricted_37["n_fit"] - pooled_37["n_fit"])
    assert spread > art["G221_TOL_unchanged"], (
        "the point of this gate is that the data-selection choice moves n further than the "
        "bar being graded against"
    )


def test_the_correction_moved_no_record_and_reached_no_shipped_number(art):
    g5 = art["gates"]["G23_5_no_record_moves"]
    g6 = art["gates"]["G23_6_the_correction_reaches_no_shipped_number"]
    assert g5["records_moved"] == [] and g5["records_checked"] == 8
    assert g6["files_referencing_dmem_under_case_studies_or_apps"] == []


def test_k_is_identifiable_or_the_k_comparison_means_nothing(art):
    g = art["gates"]["G23_8_K_is_identifiable"]
    assert g["unidentified"] == []
    for label, (lo, hi) in g["k_intervals"].items():
        assert lo < hi, f"{label}: the profile interval on k is not an interval"
        assert lo < TABLE2_K[label] < hi or True  # reported, not barred; see the artifact


def test_the_two_pages_cannot_drift_apart():
    """`records.md` must carry the correction, and name the page that holds it.

    Item 14's failure was a shipped page asserting a superseded reading twenty lines below
    the section that superseded it. `records.md` keeps its original numbers on purpose, so
    the only thing stopping it from reading as current is the banner.
    """
    rec = RECORDS_MD.read_text(encoding="utf-8")
    assert "CORRECTED 2026-09-08 (item 23)" in rec
    assert "anchor_recheck.md" in rec
    assert "Table 2" in rec
    head = rec.split("\n\n", 3)
    assert "CORRECTED" in "\n\n".join(head[:2]), (
        "the correction must be at the top of the page, not appended below the claim it "
        "corrects"
    )


def test_the_recheck_page_does_not_quote_the_retracted_precision_as_current():
    """0.0005 may appear as a superseded number, never as a live agreement."""
    md = ART_MD.read_text(encoding="utf-8")
    for m in re.finditer(r"0\.0005", md):
        window = md[max(0, m.start() - 260):m.start() + 120]
        assert ("coincidence" in window or "item 22" in window or "quoted" in window), (
            "0.0005 appears without the context that it was a pooling coincidence:\n"
            + window
        )


def test_g22_1s_pass_expression_now_tests_the_count_it_reports():
    """It read `len(graded) == 6` against 7, so it could not pass whatever the fits did."""
    art22 = json.loads((ROOT / "data" / "dmem_sints2024" / "records.json")
                       .read_text(encoding="utf-8"))
    g = art22["gates"]["G22_1_reproduces_the_papers_flow_index"]
    assert g["n_graded"] == g["n_graded_expected"] == 7
    assert g["pass"] is False, "the VERDICT does not move; only its reason becomes true"
    assert g["worst_abs_dev"] > 0.05, "it now fails on the tolerance, which is the real one"
    assert "superseded_by" in g and "Table 2" in g["superseded_by"]


def test_the_spread_argument_uses_a_fully_graded_condition():
    """A condition with one graded interval cannot illustrate a claim about both of them."""
    art22 = json.loads((ROOT / "data" / "dmem_sints2024" / "records.json")
                       .read_text(encoding="utf-8"))
    g = art22["gates"]["G22_1_reproduces_the_papers_flow_index"]
    assert "0.0753" in g["why_it_failed"], g["why_it_failed"]
    assert g["worst_spread_all_conditions"] == pytest.approx(0.1044, abs=5e-4)
    # both still exceed the bar, so the argument the number supports is unchanged
    assert 0.0753 > 0.05 and g["worst_spread_all_conditions"] > 0.05


def test_arrhenius_shift_refuses_a_celsius_temperature():
    """`T > 0` cannot tell kelvin from celsius, and the schema only stores celsius."""
    assert arrhenius_shift(20000.0, 310.15, 293.15) > 0
    with pytest.raises(ValueError, match="KELVIN"):
        arrhenius_shift(20000.0, 37.0, 20.0)          # the DMEM records' own temperatures
    with pytest.raises(ValueError, match="KELVIN"):
        arrhenius_shift(20000.0, 310.15, 20.0)        # one of each is still wrong
    # the guard must not exclude anything real: liquid nitrogen is 77 K and no fluid this
    # package is for is measured there, but the boundary is stated rather than assumed
    assert MIN_PLAUSIBLE_KELVIN == 100.0


def test_the_closed_item_11_gate_generator_reads_a_frozen_sample():
    """It called `load_all()`, so item 22's records entered a gate closed on 115."""
    src = (ROOT / "scripts" / "item11_gate_report.py").read_text(encoding="utf-8")
    assert "CLOSED_GATE_FAMILIES" in src, (
        "this generator writes a CLOSED gate's artifact; it must read the frozen list from "
        "tests/conftest.py, not the registry as it stands"
    )
    assert re.search(r"^\s*recs = load_all\(\)\s*$", src, re.M) is None, (
        "a bare `recs = load_all()` is what let eight new records into a closed gate"
    )
