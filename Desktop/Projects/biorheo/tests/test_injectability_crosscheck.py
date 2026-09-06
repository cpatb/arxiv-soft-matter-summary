"""NEXT_STEPS item 12: the Cross-vs-Casson needle-force cross-check.

Every bar is pre-registered in `case_studies/blood_injectability/crosscheck_prereg.md`,
committed at `1cd46de` BEFORE any of the comparison was computed. The bars are read out of
the artifact and checked against the case study's own `results.json`, so a loosened bar
fails here rather than passing quietly.

This file READS the committed `crosscheck.json`; it does not re-run the script, so a
failure cannot overwrite the artifact it is checking (three tests in this repo used to do
exactly that). The slow re-run lives at the foot of the file.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
CASE = ROOT / "case_studies" / "blood_injectability"
JSON = CASE / "crosscheck.json"
MD = CASE / "crosscheck.md"
PREREG = CASE / "crosscheck_prereg.md"


@pytest.fixture(scope="module")
def cc():
    assert JSON.exists(), "run case_studies/blood_injectability/crosscheck.py"
    return json.loads(JSON.read_text())


@pytest.fixture(scope="module")
def case():
    return json.loads((CASE / "results.json").read_text())


def test_the_scales_are_the_case_studys_own_not_retyped(cc, case):
    """The three comparison scales are READ from the sibling study's artifact at run time.
    If one of them is ever retyped into the script, this fails: a bar that can drift away
    from its source is a bar that can be tuned."""
    sc = cc["scales"]
    assert sc["source"] == "case_studies/blood_injectability/results.json"
    assert sc["resolution"] == max(case["entrance_loss"]["fraction_of_dp_min"],
                                  case["entrance_loss"]["fraction_of_dp_max"])
    assert sc["yield_term"] == case["yield_share_max"]
    assert sc["donor_spread"] == case["geometry"]["spread_all_20"] - 1.0
    # measured 2026-09-06: 0.987 %, 2.94 %, 128 %
    assert sc["resolution"] == pytest.approx(0.009872, rel=1e-3)
    assert sc["yield_term"] == pytest.approx(0.029415, rel=1e-3)
    assert sc["donor_spread"] == pytest.approx(1.279189, rel=1e-3)


def test_the_reference_is_the_case_study_itself(cc):
    """The Casson arm must reproduce the committed case study, or this cross-checks a
    different calculation than the one it claims to. Measured: EXACTLY 0."""
    assert cc["reference_reproduces_case_study_to"] == 0.0
    assert cc["reference"] == "casson"
    assert len(cc["donors"]) == 20


def test_g12_1_item_12s_own_bar_passes_and_is_labelled_weak(cc):
    """Item 12 proposed the donor-to-donor spread as the tolerance. It passes -- and the
    pre-registration called it WEAK in advance, because almost any pair of fits to the
    same curve would pass a 128 % bar. Recorded so the pass is not read as support."""
    g = cc["G12_1_item12_own_bar"]
    assert g["pass"] and g["worst"] < g["bar"]
    assert "establishes very little" in g["note"]


def test_g12_2_cross_and_casson_disagree_by_more_than_the_neglected_yield_term(cc):
    """**The finding, and it is not the one item 12 expected.**

    Item 12's premise was that the two families describe the same curve over the same
    envelope and the yield term contributes under 3 %, so the forces should agree. They do
    not: the worst Cross-vs-Casson force difference is 30.2 % and the median is -11.2 %,
    i.e. TEN TIMES the yield term the case study calls negligible and thirty times the
    calculation's own resolution. Carreau-Yasuda is 33.0 %.

    The case study's ORDINAL claim survives -- 30 % is still well inside the 128 % donor
    spread, so which patient the blood came from matters more than which family was
    fitted -- but its precision does not, and nothing on that page said so."""
    sc, s = cc["scales"], cc["summary"]
    for fam in ("cross", "carreau_yasuda"):
        w = s[fam]["worst_abs_rel_dev"]
        assert w > sc["yield_term"], (fam, w)
        assert w < sc["donor_spread"], (fam, w)
        assert s[fam]["band"] == "above the neglected yield term, below the donor spread"
    assert s["cross"]["worst_abs_rel_dev"] == pytest.approx(0.30179, rel=1e-3)
    assert s["cross"]["median_rel_dev"] == pytest.approx(-0.11192, rel=1e-3)
    assert s["carreau_yasuda"]["worst_abs_rel_dev"] == pytest.approx(0.33024, rel=1e-3)
    # the sign is part of the finding: both no-yield families predict a LOWER force
    assert s["cross"]["median_rel_dev"] < 0 and s["carreau_yasuda"]["median_rel_dev"] < 0


def test_g12_3_the_nulls_make_the_comparison_informative(cc):
    """The gate that decides whether G12.2 means anything. If a constant viscosity fitted
    to the whole 0.1-700 curve reproduced the Casson force as well as Cross does, the
    needle would have washed out every difference between every model and the comparison
    would say nothing about the fitted plateaux. It does not: the Newtonian null is 271 %
    off, a separation of 8.98x against the 1.6x floor this programme has used since
    2026-08-30, and it lands in a different band."""
    g = cc["G12_3_informative"]
    assert g["floor"] == 1.6
    assert g["informative"] and g["nulls_in_a_different_band"]
    assert g["separation"] >= g["floor"]
    assert g["separation"] == pytest.approx(8.9832, rel=1e-3)
    assert cc["summary"]["newtonian"]["worst_abs_rel_dev"] == pytest.approx(2.71107, rel=1e-3)
    # and the power-law null goes the OTHER way -- it keeps thinning past the envelope
    assert cc["summary"]["power_law"]["median_rel_dev"] < -0.5


def test_g12_4_the_verdict_does_not_depend_on_the_identifiability_split(cc):
    """Decided before it was looked at: only 11 of 20 Cross and 5 of 20 Carreau-Yasuda
    records are all-`measurable`. The band must be the same on both subsets, or the
    headline is the all-`measurable` one and the difference is the finding."""
    g = cc["G12_4_identifiability_split"]
    for fam, v in g.items():
        assert v["same"], (fam, v)
    assert g["cross"]["n_all_measurable"] == 11
    assert g["carreau_yasuda"]["n_all_measurable"] == 5
    assert g["newtonian"]["n_all_measurable"] == 20


def test_g12_5_no_donor_was_silently_dropped(cc):
    """A quadrature refusal is data about the fit, not a nuisance. All four families
    evaluated all 20 donors, which is itself worth pinning: `flow_rate_gnf` refuses a
    non-monotone flow curve, so 20/20 says no shipped Cross record has one."""
    for fam, v in cc["G12_5_accounting"].items():
        assert v["n_in"] == 20 and v["n_out"] == 20, (fam, v)
        assert v["failed"] == [], (fam, v)


def test_the_route_control_shows_the_effect_is_not_the_quadrature(cc):
    """Casson comes from the closed form and every other family from the quadrature, so a
    route difference would masquerade as a family difference -- and a 30 % effect is far
    too large to assume that away. Recomputing the Casson flow rate through the quadrature
    at the same pressure agrees to 1.0e-12, ten orders of magnitude below the smallest
    scale on the page."""
    rc = cc["route_control"]
    assert rc["worst_rel_dev"] < 1e-9, rc["worst_rel_dev"]
    assert rc["worst_rel_dev"] == pytest.approx(1.0152e-12, rel=1e-2)
    assert len(rc["rows"]) == 20


def test_the_plateau_control_identifies_the_mechanism_outside_the_flow_solver(cc):
    """At the needle's wall rate every family is effectively a constant viscosity, and at
    fixed Q the force is proportional to it -- so the force ratio should reproduce the
    ratio of the fits' OWN high-rate plateaux, read from the records and touching no part
    of `injectability`. It does, to 3.2 % (Cross) and 4.9 % (Carreau-Yasuda) against a
    30-33 % effect. So the families disagree BECAUSE their fitted plateaux disagree.

    The residual is not noise and is not swept up: it is the same size as the yield term
    (2.94 %), which is what it should be -- Casson has a yield stress in this geometry and
    Cross does not, so the two ratios cannot coincide exactly."""
    pc, sc, s = cc["plateau_control"], cc["scales"], cc["summary"]
    for fam in ("cross", "carreau_yasuda"):
        gap = pc[fam]["worst_gap"]
        effect = s[fam]["worst_abs_rel_dev"]
        assert gap < 0.25 * effect, (fam, gap, effect)
        assert len(pc[fam]["rows"]) == 20
    assert pc["cross"]["worst_gap"] == pytest.approx(0.031889, rel=1e-3)
    assert pc["carreau_yasuda"]["worst_gap"] == pytest.approx(0.049441, rel=1e-3)
    # the residual is the size of the yield term, not of the effect
    assert pc["cross"]["worst_gap"] < 2.0 * sc["yield_term"]


def test_the_prereg_predates_the_result_and_names_what_it_cannot_settle(cc):
    text = PREREG.read_text(encoding="utf-8")
    for token in ("G12.1", "G12.2", "G12.3", "G12.4", "G12.5", "1.6x"):
        assert token in text, token
    # the premise it checked FIRST, and found true
    assert "The premise is TRUE." in text
    # the scale it added that item 12 did not mention
    assert "unresolved systematic is of order 1 % in force" in text
    # and the limit it fixed in advance
    assert "It cannot validate the extrapolation." in text
    assert cc["prereg"].endswith("(commit 1cd46de)")


def test_the_md_quotes_its_own_json(cc):
    """Every headline number in the .md must be the JSON's."""
    md = MD.read_text(encoding="utf-8")
    s = cc["summary"]
    for fam in ("cross", "carreau_yasuda", "newtonian", "power_law"):
        assert f"{100 * s[fam]['worst_abs_rel_dev']:.4g} %" in md, fam
    assert f"{cc['G12_3_informative']['separation']:.4g}x" in md
    assert f"{cc['route_control']['worst_rel_dev']:.2e}" in md
    assert "UNINFORMATIVE" not in md.replace("reported UNINFORMATIVE", "")


@pytest.mark.slow
def test_rerunning_the_crosscheck_reproduces_the_verdicts(cc, tmp_path):
    """Re-run the whole comparison and require the same bands and the same worst values.

    RESTORES both artifacts before asserting, so a failure cannot destroy the reference
    it was checking against."""
    import shutil
    import subprocess
    import sys
    saved = {p: tmp_path / p.name for p in (JSON, MD)}
    for src, dst in saved.items():
        shutil.copy2(src, dst)
    try:
        r = subprocess.run([sys.executable, str(CASE / "crosscheck.py")],
                           capture_output=True, text=True, timeout=3600)
        assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-2000:]
        new = json.loads(JSON.read_text())
    finally:
        for src, dst in saved.items():
            shutil.copy2(dst, src)
    for fam in ("cross", "carreau_yasuda", "newtonian", "power_law"):
        assert new["summary"][fam]["band"] == cc["summary"][fam]["band"], fam
        assert new["summary"][fam]["worst_abs_rel_dev"] == pytest.approx(
            cc["summary"][fam]["worst_abs_rel_dev"], rel=1e-9), fam
    assert new["G12_3_informative"]["informative"] == cc["G12_3_informative"]["informative"]
    assert np.isclose(new["route_control"]["worst_rel_dev"],
                      cc["route_control"]["worst_rel_dev"], rtol=1e-6)
