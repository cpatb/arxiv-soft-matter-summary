"""Pins the multi-start seed-scatter study (NEXT_STEPS item 9).

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
