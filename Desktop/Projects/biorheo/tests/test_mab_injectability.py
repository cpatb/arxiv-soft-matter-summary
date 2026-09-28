"""Item 25: the wedge fluid. 58 high-concentration mAb formulations through a needle.

Pins `case_studies/mab_injectability/results.json`. Pre-registration:
`data/mab_injectability_prereg.md`; provenance: `data/mab_prasnikar2025/SOURCE.md`.

The tests that matter most here are the ones asserting what the study does NOT do: it must
not quote an injectable-force limit, it must not use the source column that would not
reproduce, and it must keep using the blood study's device rather than a device of its own.
Each of those is a way this page could quietly become wrong while every number in it stayed
correct.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
ART = ROOT / "case_studies" / "mab_injectability" / "results.json"
ART_MD = ROOT / "case_studies" / "mab_injectability" / "results.md"
BLOOD = ROOT / "case_studies" / "blood_injectability" / "results.json"

# Quoted verbatim from the paper's Methods 2.2.4 and pinned here so a silent edit to the
# script's constant is a test failure rather than a new claim about the source.
NEWTONIAN_RANGE = [960.0, 3329.0]


@pytest.fixture(scope="module")
def art() -> dict:
    assert ART.exists(), f"{ART} is missing; run case_studies/mab_injectability/run.py"
    return json.loads(ART.read_text(encoding="utf-8"))


def test_all_58_formulations_are_present_and_none_was_dropped(art):
    assert art["n_formulations"] == 58
    assert len(art["rows"]) == 58
    c = [r["c_mab_mg_per_mL"] for r in art["rows"]]
    e = [r["eta_mPa_s"] for r in art["rows"]]
    assert min(c) == pytest.approx(159.3, abs=0.1) and max(c) == pytest.approx(212.2, abs=0.1)
    assert min(e) == pytest.approx(13.57, abs=0.05) and max(e) == pytest.approx(47.07, abs=0.05)
    assert all(np.isfinite([r["force_N"] for r in art["rows"]]))


def test_the_device_is_the_blood_studys_and_has_not_drifted(art):
    """Reusing it is what makes the two fluids comparable; a private copy would rot."""
    d = art["gates"]["G25_3_force_against_the_blood_study_same_device"]["device"]
    blood = json.loads(BLOOD.read_text(encoding="utf-8"))
    assert d["barrel_radius_m"] == blood["barrel_radius_m"]
    assert d["Q_mL_per_min"] == blood["Q_mL_per_min"]
    # the blood artifact stores the needle as a nested object; compare the length through it
    assert d["needle_length_m"] == pytest.approx(blood["needle"]["length_m"], rel=1e-12)


def test_the_source_verified_newtonian_over_a_stated_range(art):
    """The one-parameter model is the source's verdict, not a convenience."""
    assert art["newtonian_range_s"] == NEWTONIAN_RANGE
    assert art["temperature_C"] == 25.0
    src = (ROOT / "data" / "mab_prasnikar2025" / "SOURCE.md").read_text(encoding="utf-8")
    # the quote is line-wrapped AND inside a markdown blockquote, so the ">" markers have
    # to come off before flattening or they land in the middle of the sentence
    flat = " ".join(ln.lstrip("> ").strip() for ln in src.splitlines())
    flat = " ".join(flat.split())
    assert "960" in flat and "3329" in flat
    assert "categorized them as Newtonian" in flat


def test_g25_4_this_is_the_least_extrapolated_fluid_in_the_project(art):
    """0.74 decades against blood's 1.42 -- the reason this dataset was worth chasing."""
    g = art["gates"]["G25_4_how_far_the_device_extrapolates"]
    assert g["pass"] is True
    assert g["decades_past_measured_max"] == pytest.approx(0.741, abs=0.02)
    assert g["decades_past_measured_max"] < 1.42, (
        "the claim that this is the least-extrapolated fluid here is relative to blood's "
        "1.42 decades; if it is no longer smaller, the SOURCE.md headline must be withdrawn"
    )
    # a Newtonian wall rate is 4Q/(pi R^3) and does NOT depend on viscosity, so every
    # formulation must land on the same rate. If they differ, the device model is wrong.
    # 1e-12 was tighter than the force -> pressure -> model round trip that recovers this
    # rate; 1e-8 still discriminates, because a viscosity dependence would show at 1e-1.
    assert g["wall_rate_min_s"] == pytest.approx(g["wall_rate_max_s"], rel=1e-8)


def test_g25_5_the_confound_control_changed_the_answer(art):
    """The raw spread mixes the excipient with the concentration; that is the whole gate."""
    g = art["gates"]["G25_5_which_lever_excipient_or_concentration"]
    raw, norm = g["raw_viscosity_ratio"], g["concentration_normalised_viscosity_ratio"]
    assert norm > raw, (
        "at matched concentration the excipient lever is LARGER than the raw spread "
        "suggests; if that ever reverses the page's framing must change"
    )
    assert norm / raw > 1.3, "the confound control must matter, or the gate was decoration"
    assert g["excipient_lever_in_mg_per_mL_normalised"] == pytest.approx(62.6, abs=2.0)
    # and it must be arithmetically consistent with the Mooney slope it is built on
    assert (math.log(norm) / g["mooney_slope_mL_per_mg"]
            == pytest.approx(g["excipient_lever_in_mg_per_mL_normalised"], rel=1e-9))


def test_g25_1_failed_and_the_reduction_column_is_therefore_unused(art):
    """A number that cannot be reproduced cannot be extended."""
    g = art["gates"]["G25_1_the_reduction_column_does_not_reproduce"]
    assert g["pass"] is False and g.get("predicted_failure") is True
    assert g["r2_main_series"] < 0.70
    # the reduction column may be CARRIED in the artifact for audit, but no gate may use it
    src = (ROOT / "case_studies" / "mab_injectability" / "run.py").read_text(encoding="utf-8")
    body = src.split("def main(")[1]
    assert "reduction_pct" not in body, (
        "the reduction column reappeared inside main(); G25.1 says it does not reproduce"
    )


def test_g25_2_the_slope_survives_even_though_the_column_does_not(art):
    g = art["gates"]["G25_2_the_mooney_slope_is_robust"]
    assert g["pass"] is True
    assert g["relative_difference"] < 0.20
    assert g["doubling_concentration_mg_per_mL"] == pytest.approx(21.8, abs=1.0)


def test_g25_6_these_formulations_are_separable_unlike_the_blood_donors(art):
    """Item 14 got 36 of 190 pairs on blood; a formulation decision here is supported."""
    g = art["gates"]["G25_6_measurement_against_formulation"]
    assert g["pass"] is True
    assert g["fraction_disjoint"] > 0.5
    assert g["blood_comparison"]["n_disjoint"] == 36
    # the discriminating half: a routine that called every pair separable would also pass
    # the first assertion, so require the measurement to be a real, non-zero term
    assert g["median_relative_force_sd"] > 0.0
    assert g["n_disjoint"] < g["n_pairs"], "not every pair should separate"


def test_g25_7_shows_the_schema_boundary_with_its_exception(art):
    g = art["gates"]["G25_7_these_cannot_be_registry_records"]
    assert g["outcome"] == "FAILED"
    assert "cannot determine" in g["exception"]
    assert g["n_points"] == 1


def test_no_registry_record_was_created_from_this_deposit():
    """G25.7's consequence. The registry must be untouched by item 25."""
    from biorheo.registry import load_all
    recs = load_all()
    assert not any("prasnikar" in k or "mab_" in k for k in recs), sorted(recs)[:5]
    assert len(recs) == 125, len(recs)


def test_the_page_quotes_no_injectable_force_limit():
    """The retracted ISO 9626 comparison is the precedent; a human-factors limit is the
    same mistake in a different standard."""
    md = ART_MD.read_text(encoding="utf-8")
    assert "No injectable-force limit is quoted" in md
    for banned in ("ISO 11608", "ISO 7886", "20 N limit", "maximum injection force of"):
        assert banned not in md, f"the page quotes {banned!r}"
    # the scale it DOES use must be named
    assert "blood study" in md


def test_the_page_says_what_it_does_not_say():
    md = ART_MD.read_text(encoding="utf-8")
    for required in ("Nothing about whether any of these is injectable by a human",
                     "Nothing about stability",
                     "Nothing that uses the source's `viscosity reduction` column",
                     "Nothing about shear-thinning at the needle wall"):
        assert required in md, required


def test_the_page_says_the_calibration_exponent_is_unpublished():
    """B3 (item 30) row 5: the source fits eta = a*exp(b*c) but never prints a, b or R2, so
    the Mooney slope this page is built on has no external anchor. The page must say so,
    and must not claim the slope 'reproduces' or 'agrees with' the authors' value."""
    md = ART_MD.read_text(encoding="utf-8")
    assert "That slope has no external check" in md
    assert "printed nowhere" in md
    for banned in ("reproduces the authors' exponent", "agrees with the published slope",
                   "the published exponent"):
        assert banned not in md, banned
