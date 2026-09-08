"""Item 24: two licensable Herschel-Bulkley records, and the first real test of the shift.

These pin the artifact written by `scripts/records_from_carbopol_pg.py`. Pre-registration
and the addendum recording two refuted predictions: `data/carbopol_pg_prereg.md`.

Two of these tests pin a **failure**. G24.3 and G24.4 were pre-registered to PASS and did
not; neither bar was moved. A test that let one of them start passing silently would be a
test that lets a recorded result change without anyone noticing, so both are asserted as
failures with their numbers.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from biorheo.registry import RECORDS_DIR, load_record

ROOT = Path(__file__).resolve().parent.parent
ART = ROOT / "data" / "carbopol_pg_vinutha2024" / "records.json"
ART_MD = ROOT / "data" / "carbopol_pg_vinutha2024" / "records.md"
IDS = ("carbopol_pg_vinutha2024_10C", "carbopol_pg_vinutha2024_50C")

# Published by the deposit in named columns, and by the paper's Methods. Not computed here
# and not computable here -- which is what makes G24.1 an anchor rather than a self-check.
MEDIUM_VISCOSITY = {"10.0": 0.102, "20.0": 0.055, "50.0": 0.012}
PUBLISHED_YIELD_STRESS = 8.9


@pytest.fixture(scope="module")
def art() -> dict:
    assert ART.exists(), f"{ART} is missing; run scripts/records_from_carbopol_pg.py"
    return json.loads(ART.read_text(encoding="utf-8"))


def test_the_two_records_load_and_are_licensable():
    for rid in IDS:
        rec = load_record(Path(RECORDS_DIR) / f"{rid}.yaml")
        assert rec.model.name == "herschel_bulkley"
        assert rec.commercial_ok is True
        assert rec.licence.name == "CC BY 4.0"
        assert {s.role for s in rec.sources} == {"dataset", "paper"}
        # the notes must carry the caveat that decides how the record may be used
        assert "NOT AN INJECTABILITY RECORD" in rec.notes
        assert "propylene glycol" in rec.notes


def test_no_record_comes_from_the_20C_block():
    """It is refused on four recomputed signals; see G24.7."""
    assert not (Path(RECORDS_DIR) / "carbopol_pg_vinutha2024_20C.yaml").exists()
    assert len(IDS) == 2


def test_g24_1_the_fitted_shift_reproduces_a_number_this_repo_did_not_compute(art):
    """The result of the item: `collapse_shift` never sees the medium viscosity."""
    g = art["gates"]["G24_1_fitted_shift_reproduces_the_published_viscosity_ratio"]
    assert art["medium_viscosity_Pa_s"] == MEDIUM_VISCOSITY
    # a_pred is the published ratio, exactly
    assert g["a_pred"] == pytest.approx(0.012 / 0.102, rel=1e-12)
    assert g["pass"] is True
    assert g["ratio"] == pytest.approx(1.1887, abs=5e-3)
    assert g["ln_ratio"] < g["ln_bar"]
    # not vacuous: the bar must be one the fit could have missed. a_fit is picked out of
    # collapse_shift's [1e-3, 1e3] bracket, six decades wide, against a 30 % bar.
    assert g["ln_bar"] < 0.3
    assert g["n_used"] >= 4


def test_g24_2_the_collapse_test_is_shown_informative_not_assumed(art):
    """A residual that does not move when `a` moves has discriminated nothing."""
    g = art["gates"]["G24_2_the_collapse_test_is_informative"]
    assert g["pass"] is True
    assert g["separation_factor"] > 3.0
    # and on a common point set, so the ratio is not partly a change of sample
    assert g["separation_factor_common_point_set"] > 3.0


def test_g24_3_FAILED_as_recorded_and_the_bar_was_not_moved(art):
    """Pre-registered to PASS at 10 %; the low-rate residual moved 69 %."""
    g = art["gates"]["G24_3_the_plateau_is_blind_to_the_shift_factor"]
    assert g["pass"] is False, (
        "this gate is a RECORDED REFUTATION. If it now passes, the result changed and the "
        "prereg addendum A.2 must be revisited rather than this assertion relaxed"
    )
    assert g["relative_move"] > 0.10
    assert "recorded_refutation" in g


def test_g24_3pp_explains_the_failure_by_measurement_not_by_argument(art):
    """The plateau is entirely outside the window a shift can act on."""
    g = art["gates"]["G24_3pp_is_the_plateau_even_inside_the_shiftable_window"]
    assert g["plateau_entirely_outside_the_shiftable_window"] is True
    assert g["plateau_upper_rate"] < g["shiftable_window_lower_rate"]
    # the split G24.3 used must be shown to be far above the plateau, or the diagnosis is
    # an assertion rather than a measurement
    assert g["g243_split_rate"] > 10 * g["plateau_upper_rate"]
    assert g["n_points_on_plateau"] >= 1


def test_g24_4_FAILED_as_recorded_with_the_scale_that_makes_it_narrow(art):
    """`tau_y` is not invariant -- by 0.022 Pa, against a replicate scatter ten times that."""
    g = art["gates"]["G24_4_yield_stress_invariant"]
    g8 = art["gates"]["G24_8_replicate_scatter_against_parameter_uncertainty"]
    assert g["pass"] is False, "a RECORDED REFUTATION; see prereg addendum A.2"
    assert g["interval_gap_Pa"] > 0
    assert g["relative_difference"] == pytest.approx(0.0356, abs=3e-3)
    # the point of reporting the scale: the statistical interval excludes the systematic,
    # and the systematic is the same size as the difference being called significant
    assert g8["worst_replicate_median_scatter"] > 0.5 * g["relative_difference"]


def test_g24_5_and_g24_6_confirm_the_rest_of_the_predicted_parameter_map(art):
    """`n` invariant and `k -> a^n k`, both predicted by item 21 before this data existed."""
    g5 = art["gates"]["G24_5_flow_index_invariant"]
    g6 = art["gates"]["G24_6_consistency_index_moves_as_a_to_the_n"]
    assert g5["pass"] is True and g6["pass"] is True
    assert g6["ln_ratio"] < g6["ln_bar"]
    # G24.6 must be a real constraint: k moves by a factor of ~3 between the two records,
    # so a routine returning k_ratio = 1 would fail it
    assert not (0.9 < g6["k_ratio"] < 1.1)


def test_g24_7_recomputes_the_refusal_rather_than_asserting_it(art):
    """All four signals reproduce, and they separate the 20 C block from the other two."""
    g = art["gates"]["G24_7_the_20C_block_is_refused"]
    assert g["pass"] is True
    assert all(g["reproduces"].values())
    s20, s10, s50 = (g["signals"]["20.0"], g["signals"]["10.0"], g["signals"]["50.0"])
    # a rate that goes BACKWARDS inside a decreasing programmed ramp is the discriminating
    # signal: the instrument cannot do it, and the two admitted blocks have none
    assert s20["n_rate_reversals"] == 2
    assert s10["n_rate_reversals"] == 0 and s50["n_rate_reversals"] == 0
    assert s20["grid_fraction"] < 0.1
    assert s10["grid_fraction"] > 0.7 and s50["grid_fraction"] > 0.7
    assert "20.0" not in [str(float(t)) for t in g["admitted_temperatures"]]


def test_g24_8_is_the_projects_first_independent_replicate(art):
    """A second flow curve of the same material under a different protocol."""
    g = art["gates"]["G24_8_replicate_scatter_against_parameter_uncertainty"]
    for T in ("10.0", "50.0"):
        r = g["replicate"][T]
        assert r["n_matched"] >= 4, "too few matched rates to call this a replicate"
        assert 0.9 < r["median_ratio"] < 1.2
    assert g["larger_term"] in ("replicate scatter", "parameter interval")


def test_g24_9_does_not_grade_a_hb_tau_y_against_a_three_component_yield_stress(art):
    """Two model forms; comparing them would call the form difference an error."""
    g = art["gates"]["G24_9_the_published_yield_stress_is_not_graded"]
    assert g["graded"] is False
    assert g["published_yield_stress_Pa"] == PUBLISHED_YIELD_STRESS
    assert "three-component" in g["published_model"]
    md = ART_MD.read_text(encoding="utf-8")
    assert "model-form comparison and not a disagreement" in md


def test_the_page_states_both_refutations_where_a_reader_meets_them():
    """A shipped page must not read as if every gate passed."""
    md = ART_MD.read_text(encoding="utf-8")
    assert "G24.3 FAILED" in md
    assert "G24.4 FAILED, and the bar is not moved" in md
    assert "inoperative on this dataset" in md
