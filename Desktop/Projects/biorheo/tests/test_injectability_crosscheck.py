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



def test_the_plateau_null_is_not_independent_evidence(cc, case):
    """A control that shares an identity with its target cannot check it.

    The Newtonian-at-`eta_c` null was added to separate "does the needle read the
    plateau?" from "do the plateaux disagree?". It cannot: Casson at `tau_y = 0` IS
    Newtonian at `eta_c`, so it recomputes the sibling study's `yield_share_max` by a
    different route. This pins that it agrees with that number -- which makes it a
    cross-implementation check worth having (one route through `flow_rate_casson`, one
    through `flow_rate_gnf`) -- and pins that the .md says so rather than claiming a
    second anchor."""
    pn = cc["plateau_null"]
    assert pn["worst_rel_dev"] == pytest.approx(case["yield_share_max"], rel=1e-4)
    # different routes, so NOT bit-identical: that is the content of the check
    assert pn["worst_rel_dev"] != case["yield_share_max"]
    md = MD.read_text(encoding="utf-8")
    assert "NOT independent evidence" in md
    assert "the same number" in md
    assert "MATCHED-PLATEAU null" in md


def test_g12_3_does_not_rest_on_which_null_was_picked(cc):
    """`max` over the two nulls was a free choice. The milder null must clear the floor
    too, or the verdict is a property of the choice."""
    g = cc["G12_3_informative"]
    assert g["separation_least_extreme_null"] >= g["floor"]
    assert g["informative_on_either_null"]
    assert g["separation_least_extreme_null"] == pytest.approx(2.705, rel=1e-3)
    assert g["separation_worst_null"] == pytest.approx(8.9832, rel=1e-3)



def test_the_yield_explanation_of_the_plateau_residual_is_withdrawn(cc):
    """**A claim this study made and its own artifact refutes.**

    "The residual is the size of the yield term, which is what it should be" is signed and
    quantitative: `F_casson = F_newt(eta_c)(1+y)` puts the force ratio BELOW the plateau
    ratio on EVERY donor. It does not. Six Cross donors and seven Carreau-Yasuda ones have
    the wrong sign, the worst gap in both families is a wrong-sign row, and CY's gap is
    1.68x the yield term rather than equal to it. The first version tested it with an
    unsigned bound and a factor-2 slop, which cannot see any of that."""
    pc, sc = cc["plateau_control"], cc["scales"]
    assert pc["cross"]["n_wrong_sign"] == 6
    assert pc["carreau_yasuda"]["n_wrong_sign"] == 7
    assert pc["cross"]["worst_gap_donor"] == "A"
    assert pc["cross"]["worst_gap_donor"] in pc["cross"]["wrong_sign_donors"]
    assert pc["carreau_yasuda"]["worst_gap"] > sc["yield_term"]
    md = MD.read_text(encoding="utf-8")
    assert "The sentence is withdrawn" in md


def test_the_matched_plateau_null_is_the_discriminating_control(cc):
    """The only null here that goes through the solver AND could have come out the other
    way: Cross with `eta_inf` set to Casson's `eta_c`, keeping its whole low-rate
    structure. It collapses 30.2 % into the yield band."""
    mpn, sc, s = cc["matched_plateau_null"], cc["scales"], cc["summary"]
    assert mpn["worst_abs_rel_dev"] == pytest.approx(0.02719, rel=1e-3)
    assert mpn["worst_abs_rel_dev"] < sc["yield_term"]
    assert mpn["worst_abs_rel_dev"] < 0.1 * s["cross"]["worst_abs_rel_dev"]
    assert len(mpn["rows"]) == 20


def test_the_in_envelope_needle_shows_the_extrapolation_is_not_the_cause(cc):
    """Run the same comparison at a needle wide enough that the wall rate is in the
    measured range. 99.5 % of the median difference survives, so the extrapolation
    contributes ~nothing to the typical donor and about a third of the worst one."""
    ien, s = cc["in_envelope_needle"], cc["summary"]
    assert ien["id_mm"] == pytest.approx(0.6236, rel=1e-3)
    c = ien["families"]["cross"]
    assert c["worst_abs_rel_dev"] == pytest.approx(0.204, rel=1e-2)
    assert c["median_rel_dev"] == pytest.approx(-0.1113, rel=1e-2)
    # the median barely moves; the worst does
    assert abs(c["median_rel_dev"] / s["cross"]["median_rel_dev"]) > 0.99
    assert c["worst_abs_rel_dev"] < s["cross"]["worst_abs_rel_dev"]
    # and it is still above the yield term, so the finding survives with NO extrapolation
    assert c["worst_abs_rel_dev"] > cc["scales"]["yield_term"]
    # honest about the residual extrapolation rather than claiming none
    assert ien["worst_decades_past_high"] < 0.2


def test_the_sign_and_misfit_claims_are_the_ones_the_data_supports(cc):
    """"Both families predict a LOWER force, consistently" was false: 4 of 20 Cross donors
    are POSITIVE, up to +14.2 %. And the families are not equally good -- Casson is worse
    on 20/20 and the force difference TRACKS the misfit."""
    s, mc = cc["summary"], cc["misfit_correlation"]
    g = s["cross"]["signs"]
    assert g["n_negative"] == 16 and g["n_positive"] == 4
    assert g["most_positive"] == pytest.approx(0.1419, rel=1e-2)
    assert g["sign_test_p"] == pytest.approx(0.01182, rel=1e-2)
    assert s["carreau_yasuda"]["signs"]["n_negative"] == 18
    for fam in ("cross", "carreau_yasuda"):
        assert mc[fam]["casson_worse_on"] == 20
        assert mc[fam]["rmse_ratio_min"] > 1.0
        assert mc[fam]["spearman_rho"] > 0.5 and mc[fam]["p"] < 0.02
    md = MD.read_text(encoding="utf-8")
    assert "OPPOSITE directions" in md and "not a symmetric" in md

def test_the_disagreement_is_misfit_not_extrapolation(cc):
    """**This study's own first explanation, refuted by its own data.**

    Finding a 30 % difference at a rate 1.42 decades past the data, both pages first said
    "the cause is the extrapolation". The families in fact disagree by MORE inside the
    measured envelope (42 % at the bottom of it) than at the needle (30.5 %), so the
    needle is not where the disagreement is created. The mechanism is that Casson misfits
    the curve: median rmse_log 0.146 against Cross's 0.030, and `residual_structure`
    flags 19/20 donors on Casson against 2/20 on Cross.

    Both pages now say this, and this pins that they do."""
    vl = cc["viscosity_ladder"]
    assert vl["worst_inside_cross"] > vl["at_needle_cross"]
    assert not vl["needle_is_the_worst_rate"]
    assert vl["worst_inside_cross"] == pytest.approx(0.4196, rel=1e-3)
    assert vl["at_needle_cross"] == pytest.approx(0.3054, rel=1e-3)
    # the worst rate is the BOTTOM of the envelope, not the top
    inside = [r for r in vl["rows"] if r["inside_envelope"]]
    assert max(inside, key=lambda r: r["cross"])["gammadot"] == min(
        r["gammadot"] for r in inside)
    # the misfit, computed from the records' own `extra`, not quoted from another .md
    assert vl["median_rmse_log"]["casson"] == pytest.approx(0.1462, rel=1e-3)
    assert vl["median_rmse_log"]["cross"] == pytest.approx(0.0301, rel=1e-3)
    assert vl["median_rmse_log"]["casson"] > 4 * vl["median_rmse_log"]["cross"]
    assert vl["n_structured_of_20"] == {"casson": 19, "cross": 2, "carreau_yasuda": 1}
    for page in (MD, CASE / "results.md"):
        t = page.read_text(encoding="utf-8")
        assert "NOT the extrapolation" in t or "was WRONG" in t, page.name


def test_the_case_study_page_carries_the_caveat(cc):
    """A reader of `results.md` must not take its forces at face value. The section is
    printed by `run.py` from `crosscheck.json`, so it cannot go stale silently."""
    t = (CASE / "results.md").read_text(encoding="utf-8")
    assert "moves these forces by more than the yield stress does" in t
    s = cc["summary"]
    assert f"{100 * s['cross']['worst_abs_rel_dev']:.4g} %" in t
    assert f"{100 * s['newtonian']['worst_abs_rel_dev']:.4g} %" in t
    # it must say which way the asymmetry runs: Cross is the BETTER fit, not an equal one
    assert "BETTER-fitting family" in t
    assert "not two equally good descriptions" in t

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
