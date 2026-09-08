"""Item 27 -- the external tube-flow anchor. Fast tests over the committed artifact.

These pin the numbers the report and the Status card quote, and they check the two things
that are easy to get wrong in prose: that a FAILED gate stays failed, and that the page
does not claim the discriminator discriminated when it did not.

The expensive re-run of `scripts/foam_pipe_anchor.py` is not here; it is a slow test.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent.parent
ART = HERE / "data" / "foam_koponen2023"
JSON = ART / "foam_pipe_anchor.json"
MD = ART / "foam_pipe_anchor.md"
PREREG = ART / "PREREG_item27.md"
SOURCE = ART / "SOURCE.md"


@pytest.fixture(scope="module")
def art() -> dict:
    return json.loads(JSON.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def md() -> str:
    return MD.read_text(encoding="utf-8")


def test_artifacts_are_committed():
    for p in (JSON, MD, PREREG, SOURCE):
        assert p.exists(), p


def test_provenance(art):
    assert art["source_doi"] == "10.5281/zenodo.8334762"
    assert art["paper_doi"] == "10.1016/j.ces.2023.119282"
    assert art["licence"] == "cc-by-4.0"
    assert art["generated_by"] == "scripts/foam_pipe_anchor.py"


def test_the_dataset_shape_is_what_the_prereg_says(art):
    assert art["n_rows_set2"] == 660
    assert art["n_cells"] == 15
    assert art["radii_mm"] == [4.075, 5.95, 9.5]


def test_middle_radius_is_5_95_not_6_0(art):
    """Addendum A.1: the prereg said 6.0 mm and the column says 11.90 mm diameter."""
    assert art["radii_mm"][1] == pytest.approx(5.95)
    assert "addendum A.1" in PREREG.read_text(encoding="utf-8")


# --- the compressibility correction, and the check that forced it in -------------------

def test_compressibility_factor_is_large_enough_to_matter(art):
    e = art["expansion"]
    assert e["n_rows"] == 660
    assert e["median"] == pytest.approx(1.5738, abs=5e-4)
    assert e["min"] > 1.0 and e["max"] > 3.0


def test_the_slip_bound_discriminates_between_the_two_reductions(art):
    """A flow with a plug cannot slip faster than its own mean velocity.

    This is the whole argument for the correction, and it is a check that CAN fail: it
    fails on the raw column and passes on the corrected one. If both passed, the
    correction would be unsupported by this evidence however right it looked.
    """
    e = art["expansion"]
    assert e["n_slipfrac_ge1_raw"] == 11, "the raw column must violate the bound"
    assert e["n_slipfrac_ge1_expanded"] == 0, "the corrected one must not"
    assert e["worst_slipfrac_raw"] > 1.0
    assert e["worst_slipfrac_expanded"] < 1.0


# --- gate verdicts, including the failures -------------------------------------------

def test_g271_is_recorded_FAILED_and_the_bar_was_not_moved(art):
    g = art["G27.1"]
    assert g["verdict"] == "FAIL"
    assert g["n_pass"] == 2 and g["n_total"] == 5
    assert g["bar"] == "|dn| <= 3.0 sigma_paper on >= 4 of 5"


def test_g271_grades_against_the_papers_FIT_TABLE_not_a_caption(art):
    """Item 23's lesson, applied before the fact: Table 2, with the paper's own sigmas."""
    published = {0.70: 0.361, 0.75: 0.363, 0.80: 0.309, 0.85: 0.265, 0.90: 0.238}
    rows = {round(r["phi"], 2): r for r in art["G27.1"]["rows"]}
    assert set(rows) == set(published)
    for phi, n in published.items():
        assert rows[phi]["n_paper"] == pytest.approx(n, abs=1e-9)
        assert rows[phi]["sigma_paper"] > 0


def test_g271_prime_is_a_size_not_a_verdict(art):
    g = art["G27.1p"]
    assert "verdict" not in g, "no bar was pre-registered for G27.1'; none may be invented"
    assert g["n_points"] == 660
    assert g["median"] == pytest.approx(0.0807, abs=2e-3)


def test_g272_passes_at_its_own_bar(art):
    g = art["G27.2"]
    assert g["bar"] == 0.20
    assert g["verdict"] == "PASS"
    assert g["median"] == pytest.approx(0.1874, abs=2e-3)
    assert g["n_rows"] == 210


def test_g272_pass_is_declared_NOT_INFORMATIVE(art):
    """The registered discriminator refused it, and the page must not read otherwise."""
    g3 = art["G27.3"]
    assert g3["informative"] is False
    assert g3["separation_ratio"] < g3["separation_bar"]
    assert g3["separation_ratio"] == pytest.approx(0.565, abs=5e-3)


def test_the_no_slip_arm_is_BETTER_which_is_the_refutation(art):
    assert art["G27.3"]["median"] < art["G27.2"]["median"]
    assert art["G27.3"]["median"] == pytest.approx(0.1058, abs=2e-3)


def test_g273_prime_measures_why_and_the_slip_term_is_subdominant(art):
    g = art["G27.3p"]
    assert g["c_monotone_in_radius"] is False, "the non-monotonicity is the mechanism"
    assert g["one_minus_c_ratio_big_over_small"] == pytest.approx(1.068, abs=5e-3)
    surviving = abs(g["one_minus_c_ratio_big_over_small"] - 1.0)
    assert surviving < art["G27.2"]["median"], "slip must come out sub-dominant HERE"


def test_positive_control_is_far_better_than_the_prediction(art):
    """Without this, G27.2's error could be misfit rather than extrapolation."""
    assert art["G27.4"]["verdict"] == "PASS"
    assert art["G27.4"]["median"] == pytest.approx(0.0240, abs=2e-3)
    assert art["G27.4"]["median"] < 0.5 * art["G27.2"]["median"]


def test_reverse_direction_refuted_the_registered_prediction(art):
    assert art["G27.5"]["worse_than_G27.2"] is False
    assert art["G27.5"]["median"] < art["G27.2"]["median"]


def test_g276_is_recorded_FAILED_with_its_number(art):
    g = art["G27.6"]
    assert g["verdict"] == "FAIL"
    assert g["bar"] == 0.50
    assert g["median"] > 1.0, "Mooney vs imaged slip is worse than a factor of two"


def test_g276_prime_reports_the_sign_and_the_inadmissible_answers(art):
    """An unsigned |rel err| hides both. 25 negative slip velocities is the finding."""
    g = art["G27.6"]
    assert g["signed_median"] > 0, "Mooney runs HIGH; the sign is part of the claim"
    assert g["n_mooney_negative"] == 25
    assert g["median_ratio_mooney_over_video"] == pytest.approx(2.198, abs=1e-2)


def test_g277_is_reported_but_no_gate_is_graded_on_it(art):
    g = art["G27.7"]
    assert "verdict" not in g and "bar" not in g
    assert g["median"] == pytest.approx(0.0638, abs=2e-3)


# --- the prose cannot drift from the artifact ----------------------------------------

def test_page_says_NOT_INFORMATIVE_in_words(md):
    assert "NOT INFORMATIVE" in md
    assert "INFORMATIVE**" not in md.replace("NOT INFORMATIVE**", "")


def test_page_does_not_call_the_foam_a_biofluid(md):
    low = md.lower()
    for banned in ("blood", "antibody", "mab formulation", "patient"):
        # they may appear only inside the explicit "does not transfer" sentence
        for m in re.finditer(banned, low):
            window = low[max(0, m.start() - 200): m.start() + 60]
            assert ("no number here transfers" in window
                    or "does not transfer" in window
                    or "not a biofluid" in window), f"unqualified mention of {banned!r}"


def test_scorecard_records_five_refutations(md):
    assert md.count("**REFUTED**") == 5
    assert "Five of seven refuted" in md


def test_dropped_rows_are_counted_not_hidden(art, md):
    assert art["G27.2"]["n_rows_unpredicted"] == 15
    assert art["G27.2"]["n_cells_dropped"] == 1
    assert "counted rather than dropped" in md


def test_no_registry_record_was_created_by_this_item():
    """The reduced curve is DERIVED; the registry has never held a derived curve."""
    recs = HERE / "src" / "biorheo" / "registry" / "records"
    assert not list(recs.glob("*foam*")), "item 27 must not have added a record"
    assert "No record enters the registry from this item" in MD.read_text(encoding="utf-8")


# --- the slow re-run -----------------------------------------------------------------

@pytest.mark.slow
def test_rerunning_reproduces_every_gate():
    """The fits are seeded, so a re-run must be byte-identical. ~3 minutes.

    NON-DESTRUCTIVE: the script overwrites the committed artifacts, so both are saved and
    RESTORED before any assertion can fail. The raw CSVs are committed, so their absence
    is a real error rather than a reason to skip.
    """
    import subprocess
    import sys as _sys

    raw = ART / "raw" / "Set_2_Data_PVOH_foam_rheology.csv"
    assert raw.exists(), raw
    original_json = JSON.read_text(encoding="utf-8")
    original_md = MD.read_text(encoding="utf-8")
    try:
        r = subprocess.run([_sys.executable, str(HERE / "scripts" / "foam_pipe_anchor.py")],
                           capture_output=True, text=True, cwd=HERE)
        rc, err = r.returncode, r.stderr[-2000:]
        after = JSON.read_text(encoding="utf-8")
    finally:
        JSON.write_text(original_json, encoding="utf-8", newline="")
        MD.write_text(original_md, encoding="utf-8", newline="")
    assert rc == 0, err
    assert json.loads(after) == json.loads(original_json)
