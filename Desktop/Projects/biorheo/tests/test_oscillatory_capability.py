"""Pins the item-7 capability verdict: the Horner workbook's oscillatory sheets CANNOT
say whether the blood step rise is a lag or a second timescale.

The five criteria and the verdict rule are pre-registered in
`scripts/oscillatory_capability_horner.py`'s docstring. Four of the five FAILED and the
one that passed (C5) is the one NEXT_STEPS item 7 actually asked about -- the item's own
premise, that sampling density is the constraint, is wrong.

Fast: read the committed `data/horner2020/oscillatory_capability.json`. The workbook
itself is gitignored (`SOURCE.md` has the URL and md5), so the slow half re-runs the
script only where the file is present.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
JSON = ROOT / "data" / "horner2020" / "oscillatory_capability.json"
MD = ROOT / "data" / "horner2020" / "oscillatory_capability.md"
WORKBOOK = ROOT / "data" / "horner2020" / "Horner-Thesis-Data.xlsx"


@pytest.fixture(scope="module")
def cap():
    assert JSON.exists(), "run scripts/oscillatory_capability_horner.py"
    return json.loads(JSON.read_text())


def _c(cap, prefix):
    return next(c for c in cap["criteria"] if c["name"].startswith(prefix))


def test_the_verdict_is_closed_and_the_rule_produced_it(cap):
    """The verdict rule was fixed before grading: the linear route is closed if C2, C3
    or C4 fails. All three did, and C1 as well."""
    assert cap["verdict"] == "CLOSED"
    assert cap["linear_route_open"] is False
    assert not _c(cap, "C1")["pass"]
    assert not _c(cap, "C2")["pass"]
    assert not _c(cap, "C3")["pass"]
    assert not _c(cap, "C4")["pass"]


def test_the_criterion_that_passed_is_the_one_item_7_asked_about(cap):
    """C5 is item 7's own question -- 'do the LAOS sheets have enough samples per cycle'.
    It passes by two orders of magnitude, which is why the item cannot be closed on the
    ground it proposed. Recording this stops a later session re-opening item 7 by
    re-asking the question that was never the constraint."""
    c5 = _c(cap, "C5")
    assert c5["pass"]
    spc = [r["samples_per_cycle"] for r in cap["laos"]]
    assert min(spc) >= 10 and min(spc) > 19 and max(spc) > 1500
    assert len(cap["laos"]) == 52


def test_no_phase_channel_exists_anywhere_in_the_workbook(cap):
    """C4, the cleanest of the four failures: a lag and a second timescale differ in
    PHASE, and the workbook stores a storage modulus and nothing conjugate to it."""
    assert cap["channels"]["phase_like"] == {}
    labels = set(cap["channels"]["labels"])
    assert "Storage modulus" in labels
    for absent in ("Loss modulus", "Phase angle", "tan(delta)", "Complex viscosity"):
        assert absent not in labels


def test_c3_fails_under_the_reading_that_is_not_guessed(cap):
    """The Read Me's 'a strain amplitude of 10' has two defensible units and the source
    does not say which. Both are carried; the rule fixed in advance requires every
    reading to pass, and one does not. The supporting measurements hold under both."""
    c3 = _c(cap, "C3")
    assert len(c3["readings"]) == 2
    passes = [v["pass"] for v in c3["readings"].values()]
    assert sorted(passes) == [False, True], c3["readings"]
    lve = cap["lve_plateau"]
    # nothing below 9.95 % strain was measured at all, on either ramp -- the decisive
    # fact, and the one that needs no reading of the ambiguous unit
    assert lve["strain_min_pct"] == pytest.approx(9.95271, rel=1e-5)
    # measured 2026-09-05: still rising at that point in BOTH ramps on 9 of 15 donors
    assert lve["n_donors"] == 15
    assert lve["n_rising_at_first_point_in_every_repeat"] == 9
    # and at ~1000 % strain G' is down to a median 0.207 of its own per-ramp maximum
    assert lve["median_ratio_near_1000_pct"] == pytest.approx(0.2067, rel=1e-2)


def test_the_sweeps_are_two_ramps_and_are_not_treated_as_one(cap):
    """Both sweep sheets concatenate two ramps over the same axis. An argmax over the
    whole column would mix them; and their spread at the same strain is far too large to
    call them repeats, which is reported rather than assumed away."""
    assert cap["amplitude_sweep"]["points_per_repeat"] == [26, 26]
    assert cap["frequency_sweep"]["points_per_repeat"] == [20, 20]
    lve = cap["lve_plateau"]
    assert lve["median_repeat_spread_at_min_strain"] > 0.5, lve
    for d, v in lve["per_donor"].items():
        assert v["n_repeats"] == 2, d
        for r in v["repeats"]:
            assert r["n"] == 26, (d, r)


def test_the_population_is_not_the_step_rise_population(cap):
    """C1. Item 2's finding is on all 20 donors; no oscillatory sheet has more than 13,
    and eight donors have no LAOS data at all."""
    assert cap["donors_with_no_laos"] == "ABCDFQRT"
    assert len(cap["laos_donors_union"]) == 12
    assert len(cap["laos_donors_every_sheet"]) == 10
    assert cap["frequency_sweep"]["n_donors"] == 13


def test_the_md_quotes_its_own_json(cap):
    md = MD.read_text(encoding="utf-8")
    assert md.startswith("# Can the oscillatory sheets")
    assert "**CLOSED.**" in md
    assert cap["donors_with_no_laos"] in md
    for c in cap["criteria"]:
        assert c["measured"] in md, c["name"]


@pytest.mark.slow
def test_rerunning_reproduces_the_verdict():
    if not WORKBOOK.exists():
        pytest.skip("workbook is gitignored; see data/horner2020/SOURCE.md")
    before = JSON.read_text()
    r = subprocess.run([sys.executable,
                        str(ROOT / "scripts" / "oscillatory_capability_horner.py")],
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr[-2000:]
    assert json.loads(JSON.read_text()) == json.loads(before)
