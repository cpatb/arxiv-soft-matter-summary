"""Pins the multi-start seed-scatter studies (NEXT_STEPS items 9 and 10).

Three criteria, fixed in `scripts/multistart_seed_scatter.py`'s docstring before the
held-out families were looked at. **One passed and two failed, and no bar was moved.**

  G9.1 PASSED  -- 54 all-`measurable` records, worst seed spread 5.997e-05 (bar 1e-3).
  G9.2 FAILED  -- 1 of 80 (model, donor) pairs moves more than 1e-1 without carrying a
                  `degenerate` or `not_identifiable` verdict: Herschel-Bulkley on donor
                  G, 4.149 in `tau_y`, graded `bounded_only`.
  G9.3 FAILED over all 80 on the same record, and PASSES over the 43 export-accepted
                  records in the families that actually ship (worst 5.997e-05).

The load-bearing statement is the narrow one: **the identifiability verdicts predict seed
scatter for everything the registry ships, and there is one demonstrated counterexample in
a family that does not ship.** The counterexample is what stops `bounded_only` being read
as "determined, with a one-sided interval".

The objective is what makes this a statement about FLATNESS rather than about a failed
search: it moves by at most 6.5e-10 relative across all 80 pairs, so no seed is landing on
a worse optimum.

Item 10 closes the scope gap item 9 left: the 15 Carbopol records were never in the blood
study. On that independent set (different material, instrument, 54 points per curve, and a
DIFFERENT FITTER for the shipped values) all three criteria PASS -- worst spread 5.409e-06,
seed-stable verdicts on 15/15, and all 225 (seed, parameter, record) intervals still
contain the value the record ships. So the blood failure is one record in one family, not
a property of the identify layer.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
JSON = ROOT / "data" / "horner2020" / "multistart_scatter.json"
MD = ROOT / "data" / "horner2020" / "multistart_scatter.md"
WORKBOOK_CSV = ROOT / "data" / "horner2020" / "steady_shear_human.csv"


@pytest.fixture(scope="module")
def ms():
    assert JSON.exists(), "run scripts/multistart_seed_scatter.py"
    return json.loads(JSON.read_text())


def test_the_bars_are_the_ones_that_were_pre_registered(ms):
    """A loosened bar fails here rather than passing quietly."""
    assert ms["G91_bar"] == 1e-3
    assert ms["G92_bar"] == 1e-1
    assert ms["G93_bar"] == 1e-3
    assert ms["seeds"] == [0, 1, 2, 3, 4]
    assert ms["families"] == ["casson", "herschel_bulkley", "cross", "carreau_yasuda"]
    assert ms["exploratory"] == ["carreau_yasuda"]
    assert ms["shipped_families"] == ["casson", "cross", "carreau_yasuda"]
    assert sum(len(v) for v in ms["records"].values()) == 80


def test_g91_passed_no_all_measurable_record_wanders(ms):
    """The criterion that did hold: a record whose parameters are all `measurable` does
    not move under a change of seed."""
    assert ms["G91_violations"] == []
    assert ms["n_all_measurable"] == 54
    assert ms["worst_spread_all_measurable"] == pytest.approx(5.996843e-05, rel=1e-4)
    assert ms["worst_spread_all_measurable"] < ms["G91_bar"]


def test_g92_FAILED_and_is_recorded_as_failed(ms):
    """G9.2 said every record moving more than 1e-1 must carry a `degenerate` or
    `not_identifiable` verdict. It does not. The bar stays where it was written and the
    failing number is pinned, so a later run that quietly stops reproducing it fails
    here."""
    assert ms["G92_violations"] == [["herschel_bulkley", "G", pytest.approx(4.1486, rel=1e-4)]]
    g = ms["records"]["herschel_bulkley"]["G"]
    assert g["worst_param"] == "tau_y"
    assert g["verdicts"]["tau_y"] == "bounded_only"
    # neither of the verdicts an export refuses -- which is the whole point of the failure
    assert not g["has_bad_verdict"]
    assert g["param_spread"] > ms["G92_bar"]


def test_g93_fails_over_everything_but_passes_over_what_ships(ms):
    """The stricter replacement. It fails on the same single record, and that record is
    Herschel-Bulkley, which is fitted for `model_comparison.md` and never shipped as a
    blood record. Over the families that do ship it passes by two orders of magnitude."""
    assert ms["G93_violations"] == [["herschel_bulkley", "G", pytest.approx(4.1486, rel=1e-4)]]
    assert ms["G93_violations_shipped"] == []
    assert ms["n_export_accepted"] == 63
    assert ms["n_export_accepted_shipped"] == 43
    assert ms["worst_spread_export_accepted_shipped"] == pytest.approx(5.996843e-05, rel=1e-4)
    assert ms["worst_spread_export_accepted_shipped"] < ms["G93_bar"]


def test_the_objective_does_not_move_so_this_is_flatness_not_a_failed_search(ms):
    """If the cost moved, some seeds would be landing on worse optima and the finding
    would be about convergence instead. It does not move."""
    assert ms["worst_cost_spread"] < 1e-8, ms["worst_cost_spread"]
    for fam, donors in ms["records"].items():
        for d, r in donors.items():
            assert r["cost_spread"] < 1e-8, (fam, d, r["cost_spread"])


def test_the_movers_are_the_three_recorded_ones(ms):
    movers = [(f, d) for f, ds in ms["records"].items() for d, r in ds.items()
              if r["param_spread"] > ms["G92_bar"]]
    assert sorted(movers) == [("carreau_yasuda", "B"), ("carreau_yasuda", "L"),
                              ("herschel_bulkley", "G")]
    assert ms["n_movers"] == 3


def test_the_md_quotes_its_own_json(ms):
    md = MD.read_text(encoding="utf-8")
    assert "**FAILED**" in md, "the two failures must be visible on the page"
    assert f"{ms['worst_cost_spread']:.2e}" in md
    assert f"{ms['worst_spread_all_measurable']:.3e}" in md


@pytest.mark.slow
def test_rerunning_reproduces_the_verdicts():
    """The fits are seeded, so a re-run must be identical. ~20 minutes."""
    if not WORKBOOK_CSV.exists():
        pytest.skip("source CSV missing")
    before = json.loads(JSON.read_text())
    r = subprocess.run([sys.executable,
                        str(ROOT / "scripts" / "multistart_seed_scatter.py")],
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr[-2000:]
    assert json.loads(JSON.read_text()) == before


# ------------------------------------------------ NEXT_STEPS item 10: the Carbopol set

CARBOPOL_JSON = ROOT / "data" / "carbopol_multistart_scatter.json"
CARBOPOL_MD = ROOT / "data" / "carbopol_multistart_scatter.md"


@pytest.fixture(scope="module")
def cms():
    assert CARBOPOL_JSON.exists(), "run scripts/multistart_seed_scatter_carbopol.py"
    return json.loads(CARBOPOL_JSON.read_text())


def test_carbopol_bars_are_the_pre_registered_ones(cms):
    assert cms["G101_bar"] == 1e-3 and cms["G102_bar"] == 1e-3
    assert cms["seeds"] == [0, 1, 2, 3, 4]
    assert len(cms["records"]) == 15


def test_carbopol_all_three_criteria_pass(cms):
    """G10.1/G10.2/G10.3 on an independent set -- different material, instrument, 54
    points per curve instead of 16, and a different fitter for the shipped values. The
    blood set failed the equivalent of G10.2 on one record; this set does not."""
    assert cms["G101_violations"] == []
    assert cms["G102_violations"] == []
    assert cms["G103_violations"] == []
    assert cms["n_all_measurable"] == 15 and cms["n_export_accepted"] == 15
    assert cms["worst_spread_all_measurable"] == pytest.approx(5.409e-06, rel=1e-3)
    assert cms["worst_spread_all_measurable"] < cms["G101_bar"]
    assert cms["worst_cost_spread"] < 1e-9


def test_carbopol_containment_holds_at_every_seed_not_just_seed_zero(cms):
    """G10.3, the one specific to this set. `write_identifiability` refuses a sidecar
    whose interval does not contain the shipped csv value, and the loader then refuses a
    record outside its sidecar interval -- both evaluated only at seed 0 until now. All
    5 x 3 x 15 = 225 intervals contain it."""
    assert cms["total_interval_failures"] == 0
    n = 0
    for rid, r in cms["records"].items():
        assert r["n_interval_failures"] == 0, rid
        for key, c in r["containment"].items():
            assert c["lo"] <= c["shipped"] <= c["hi"], (rid, key, c)
            n += 1
    assert n == 5 * 3 * 15 == 225


def test_carbopol_verdicts_do_not_depend_on_the_seed(cms):
    """A verdict that changes with the seed would be a statement about the optimiser; the
    sidecar records only seed 0's, so this has to hold for the sidecar to mean anything."""
    assert cms["verdict_unstable_records"] == []
    for rid, r in cms["records"].items():
        sets = {tuple(sorted(v.items())) for v in r["verdicts_by_seed"].values()}
        assert len(sets) == 1, (rid, sets)


def test_carbopol_md_quotes_its_own_json(cms):
    md = CARBOPOL_MD.read_text(encoding="utf-8")
    assert "**FAILED**" not in md
    assert f"{cms['worst_spread_all_measurable']:.3e}" in md
    assert f"{cms['worst_cost_spread']:.2e}" in md
