"""The extended-domain refit (2026-09-08): gates from PREREG_domain_refit_2026-09-08.md.

Reads the committed artifact rather than re-running the CFD. What these pin is the thing
the grid study could not settle: whether the design propagation is still evaluating the
fitted closure outside the band it was fitted on.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
ART = ROOT / "cfd_openfoam" / "domain_refit.json"
MD = ROOT / "DOMAIN_REFIT_2026-09-08.md"
PREREG = ROOT / "PREREG_domain_refit_2026-09-08.md"
GRID_STUDY = ROOT / "GRID_STUDY_2026-09-03.md"


@pytest.fixture(scope="module")
def art() -> dict:
    return json.loads(ART.read_text(encoding="utf-8"))


def test_artifacts_exist():
    for p in (ART, MD, PREREG):
        assert p.exists(), p


def test_gb1_positive_control_passes_and_is_not_vacuous(art):
    """The wide case must reproduce the narrow one on the MATCHED band.

    Also asserted: the difference is not exactly zero. Two runs on different meshes and
    different domains agreeing to machine precision would mean the same numbers were read
    twice, and the gate would be vacuous rather than passing.
    """
    g = art["GB1"]
    assert g["verdict"] == "PASS"
    assert g["rel_C_P"] <= g["bar"] and g["rel_n_t"] <= g["bar"]
    assert g["rel_C_P"] > 1e-6 and g["rel_n_t"] > 1e-6, "identical => not two solves"
    assert g["rel_C_P"] == pytest.approx(0.00098, abs=2e-4)
    assert g["rel_n_t"] == pytest.approx(0.00212, abs=3e-4)


def test_gb2_one_power_law_still_reaches_r_over_d_40(art):
    g = art["GB2"]
    assert g["verdict"] == "PASS"
    assert g["r2"] >= g["bar"]
    assert g["r2"] == pytest.approx(0.99756, abs=1e-4)
    # widening the band did not cost fit quality
    assert g["r2"] > art["narrow"]["r2"]


def test_gb31_the_registered_reason_holds(art):
    """The far field decays FASTER than the [2, 20] power law, so n_t must rise."""
    g = art["GB3.1"]
    assert g["n_t_rose"] is True
    assert g["n_t_extended"] > g["n_t_narrow"]
    assert g["extrapolation_high_by"] == pytest.approx(0.380, abs=0.01)


def test_gb33_the_section4_defect_was_real(art):
    """Under the refit the design point is STILL above r/d = 20, so the narrow band would
    have been extrapolating whichever C_P row was used -- the defect was not an artefact
    of the grid study's row choice."""
    g = art["GB3.3"]
    assert g["still_above_20_under_refit"] is True
    assert g["inside_narrow_band"] is False
    assert g["inside_extended_band"] is True
    assert g["r_max_over_d_star_extended"] == pytest.approx(21.47, abs=0.05)


def test_gb34_the_design_number_moves_and_by_how_much(art):
    g = art["GB3.4"]
    assert g["smaller_than_grid_effect"] is True
    assert g["rel_move"] == pytest.approx(0.1428, abs=2e-3)
    assert 0.0 < g["rel_move"] < g["compare_against_grid_C_P_effect"]


def test_the_extended_CP_coincidence_is_named_not_traded_on(art):
    """The extended-band C_P lands within 0.1 % of the RETRACTED unconverged refine = 2
    value 1.5178. They are unrelated, and the page has to say so or a reader will read one
    as corroborating the other."""
    cp = art["wide_extended"]["C_P"]
    assert abs(cp - 1.5178) / 1.5178 < 0.01, "if this stops being true, drop the warning"
    md = MD.read_text(encoding="utf-8")
    assert "coincidence" in md.lower()
    assert "1.5178" in md


def test_the_page_does_not_claim_a_convergence_order(art):
    """Domain extent and mesh refinement are different axes. The JFE blocker stands."""
    md = MD.read_text(encoding="utf-8")
    assert "establishes no order of convergence" in md.lower()
    assert "diverges under refinement" in md
    assert "policy is not satisfied" in md.lower() or "not satisfied" in md


def test_both_closures_are_reported_side_by_side(art):
    """GB4: no number from this item may be quoted without its band."""
    md = MD.read_text(encoding="utf-8")
    assert "[2, 20]" in md and "[2, 40]" in md
    for key in ("design_narrow", "design_extended"):
        d = art[key]
        assert f"{d['Q_star_Lmin']:.4f}" in md
        assert f"{d['r_max_over_d_star']:.2f}" in md


def test_the_grid_study_section4_is_marked_as_addressed():
    """The sentence this item exists to repair must point at the repair."""
    gs = GRID_STUDY.read_text(encoding="utf-8")
    assert "DOMAIN_REFIT_2026-09-08" in gs, \
        "GRID_STUDY section 4 must reference the item that removed its extrapolation"
