"""Every capability claim in `../biofluids_venture/OFFER.md` is checked against an artifact.

A sales document is the one page in this project that will be read by someone who cannot
check it, and it is the one most likely to go stale: it is edited for a call, not for a
commit, and nothing else in the repo imports it. This repo has shipped prose that
contradicted its own artifact more than once, and the fix every time was a test that makes
the two impossible to separate. This is that test, pointed at the document where being
wrong costs the most.

Each test names the sentence it is defending. If one fails, the OFFER is making a claim
the code no longer supports -- fix the OFFER, or fix the code, but do not delete the test.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
OFFER = REPO.parent / "biofluids_venture" / "OFFER.md"
EXAMPLES = REPO / "case_studies" / "deliverable_example"
NEEDLE_BORE = REPO / "data" / "needle_bore.json"
RECORDS = REPO / "src" / "biorheo" / "registry" / "records"


@pytest.fixture(scope="module")
def offer() -> str:
    assert OFFER.exists(), f"{OFFER} missing"
    return OFFER.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def examples() -> dict:
    return json.loads((EXAMPLES / "summary.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def bore() -> dict:
    return json.loads(NEEDLE_BORE.read_text(encoding="utf-8"))


def test_the_bore_exponents_quoted_are_the_measured_ones(offer, examples):
    """OFFER §2 row 4: 'Blood in a fine needle: -3.98 ... A 0.5 wt% microgel: -2.31'."""
    assert f"{examples['01_clean']['bore_exponent']:.3g}" == "-3.98"
    assert f"{examples['04_thinning']['bore_exponent']:.3g}" == "-2.31"
    assert "`-3.98`" in offer and "`-2.31`" in offer


def test_the_82_percent_over_prediction_is_the_measured_one(offer, bore):
    """OFFER §2 row 4: 'over-predicts the force at a 30 % finer bore by 82 %'."""
    row = next(r for r in bore["rows"] if r["id"] == "carbopol_etd2020_0p5wt_40uL")
    assert bore["geometry"]["bore_ratio"] == 0.70
    assert round(row["newtonian_rule_over_prediction"] * 100) == 82
    assert "**82 %**" in offer


def test_the_blood_worst_case_is_1_24_percent(offer, bore):
    """OFFER §2 row 4: 'right to within 1.24 % on all 20 donors tested'."""
    g = bore["gates"]["G18_5_newtonian_rule_error_split"]
    assert g["blood_n"] == 20
    assert f"{g['blood_worst'] * 100:.2f}" == "1.24"
    assert "**1.24 %**" in offer and "all 20 donors tested" in offer


def test_the_derating_range_quoted_is_the_measured_one(offer, bore):
    """OFFER §2 row 5: 'that ratio runs 0.12 to 0.97 ... a factor of 8'."""
    g = bore["gates"]["G18_6_derating"]
    assert f"{g['strongest_derating']:.2f}" == "0.12"
    assert f"{g['worst_factor']:.2f}" == "0.97"
    assert round(1.0 / g["strongest_derating"]) == 8
    assert "0.12 to 0.97" in offer and "a factor of 8" in offer


def test_the_extrapolation_distance_quoted_is_the_measured_one(offer, examples):
    """OFFER §2 row 3: 'an extrapolation of 1.42 decades'."""
    assert f"{examples['01_clean']['decades_past_envelope']:.2f}" == "1.42"
    assert "1.42 decades" in offer


def test_the_eta_c_exponent_quoted_is_the_measured_one(offer):
    """OFFER §2 row 2: 'the force goes as eta_c^0.995'."""
    uq = json.loads((REPO / "case_studies" / "blood_injectability" / "uq.json")
                    .read_text(encoding="utf-8"))
    worst = max(abs(r["sensitivities"]["eta_c"]) for r in uq["rows"])
    assert f"{worst:.3f}" == "0.995"
    assert "eta_c^0.995" in offer


def test_the_record_count_and_licence_claim_are_true(offer):
    """OFFER §6: 'All 115 records in the registry are commercial_use: false'."""
    from biorheo.registry.loader import load_all
    recs = load_all()
    assert len(recs) == 115
    assert all(not r.commercial_ok for r in recs.values())
    assert "All 115 records" in offer and "commercial_use: false" in offer


def test_the_test_count_claim_is_not_wildly_stale(offer):
    """OFFER §6 quotes a test count. It cannot be exact (this file is in it), but a claim
    that drifts by more than 10 % has stopped being a description."""
    import re
    m = re.search(r"\((\d+) automated tests", offer)
    assert m, "OFFER §6 no longer states a test count"
    claimed = int(m.group(1))
    actual = sum(1 for p in (REPO / "tests").glob("test_*.py")
                 for line in p.read_text(encoding="utf-8").splitlines()
                 if line.startswith("def test_"))
    # function count is a lower bound on the collected count (parametrised tests expand)
    assert actual <= claimed <= actual * 3, (claimed, actual)


def test_the_offer_points_at_files_that_exist(offer):
    """Every path the OFFER names must resolve, or a call is made against a broken link."""
    for rel in ("case_studies/deliverable_example/", "biorheo.apps.report",
                "biorheo.solver.openfoam"):
        assert rel in offer, rel
    assert (EXAMPLES / "README.md").exists()
    assert (REPO / "src" / "biorheo" / "apps" / "report.py").exists()
    assert (REPO / "src" / "biorheo" / "solver" / "openfoam.py").exists()


def test_the_offer_names_the_capabilities_that_do_not_exist(offer):
    """The section that makes the document usable. If a future edit deletes a gap because
    it reads badly, this fails -- the gaps are the reason the claims are believable."""
    for gap in ("No temperature axis", "No wall slip", "No bench",
                "No commercially-usable reference data", "No regulatory authoring"):
        assert gap in offer, gap


def test_the_offer_refuses_hourly_billing_and_free_pilots(offer):
    """A design rule adopted from the research, not a preference."""
    assert "No hourly billing, no free pilot" in offer


def test_the_offer_says_its_prices_are_inferred(offer):
    """The research found no public price for this service anywhere. A price table that
    does not say so reads as a market rate."""
    assert "inferred, not observed" in offer


def test_the_offer_does_not_name_a_needle_gauge(offer):
    """Standing repo rule: the ISO 9626 comparison was retracted and the standard is still
    not in hand. A gauge number in a sales document is exactly where it would do damage."""
    for banned in ("27G", "30G", "27 G", "30 G", "ISO 9626"):
        assert banned not in offer, banned
