"""Pins the multi-start seed-scatter studies (NEXT_STEPS items 9 and 10).

Three criteria, fixed in `scripts/multistart_seed_scatter.py`'s docstring. **All three
PASS.** An earlier run reported two of them FAILED; that failure is RETRACTED and the
retraction is the more useful part of this file.

  G9.1 no all-`measurable` record moves more than 1e-3 of its own confidence interval.
  G9.2 nothing `check_exportable` accepts moves more than 1e-1.
  G9.3 the same predicate at 1e-3. G9.3 implies G9.2; they are one nested criterion with
       two bars, not two independent results.

**THE RETRACTION.** The first version divided the seed-to-seed span by `|mean|`. For
Herschel-Bulkley on donor G that gave 4.149 -- "a 415 % move" -- and a docstring claim, a
pinned test and a new API flag were built on it. The five seeds return `tau_y` =
1.61e-19, 9.81e-19, 5.00e-22, 3.27e-20, 5.94e-21 Pa against a smallest measured stress of
1.583e-03 Pa, with an identical profile interval on all five: the optimiser
converged the parameter to ZERO five times, and the ratio was dust over dust. Had
`least_squares` returned exactly 0.0 the criterion would have passed. The metric now
measures a span in units of the parameter's own 95 % profile interval -- a parameter that
moves less than its stated uncertainty has not moved -- and the retracted ratio is still
computed and printed so the number stays visible.

That is this repo's own lesson ("a relative sensitivity whose denominator is itself
shrinking measures the shrinkage") applied to its own instrument, one session after it was
written down.

The objective moves by at most 6.5e-10 in `rmse_log` across all 80 pairs -- 1.3e-9 in the
least-squares cost, which is what the optimiser minimises -- so where parameters do move,
it is flatness and not a failed search.

Item 10 asks the same question of the 15 Carbopol records, an independent set (different
material, instrument, 54 points per curve, and a DIFFERENT FITTER for the shipped values).
All three of its criteria pass too.
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
    assert sum(len(v) for v in ms["records"].values()) == 80
    # the metric is part of the criterion, so a silent change to it must fail here
    assert ms["metric"].startswith("seed span in units of the parameter")


def test_all_three_criteria_pass(ms):
    assert ms["G91_violations"] == []
    assert ms["G92_violations"] == []
    assert ms["G93_violations"] == []
    assert ms["n_all_measurable"] == 54
    assert ms["n_export_accepted"] == 63
    assert ms["worst_spread_all_measurable"] == pytest.approx(2.9724e-05, rel=1e-3)
    assert ms["worst_spread_all_measurable"] < ms["G91_bar"]
    assert ms["worst_spread_export_accepted"] < ms["G93_bar"]


def test_the_retracted_metric_is_still_visible_and_is_the_artefact_it_was(ms):
    """The retracted number stays computed so nobody has to take the retraction on
    trust. Donor G's `tau_y` is zero to eighteen decimal places on all five seeds; the
    4.149 was that span over that mean."""
    assert ms["worst_spread_over_mean"] == pytest.approx(4.1486, rel=1e-4)
    g = ms["records"]["herschel_bulkley"]["G"]
    assert g["param_spread_over_mean"] == pytest.approx(4.1486, rel=1e-4)
    vals = g["values_by_seed"]["tau_y"]
    assert len(vals) == 5 and max(abs(v) for v in vals) < 1e-18, vals
    # against the parameter's own interval it is nothing, which is the corrected metric
    # 5.488e-04 since NEXT_STEPS item 13. This pin used to read 0.01953, and that
    # number was itself a bisection-grid point: 0.01953125 = 1000 * 1e5/5.12e6 for
    # the Herschel-Bulkley tau_y box bound of 1e5. The retraction this test guards
    # does not depend on the width -- it depends on the interval being IDENTICAL
    # across the five seeds while the parameter sits at ~1e-19 -- so the argument
    # is unchanged and only the number moved.
    assert g["ci_width"]["tau_y"] == pytest.approx(5.4877e-04, rel=1e-3)
    # tau_y's OWN contribution, not the record max (which is set by k at ~3e-09)
    # 1.79e-15 since item 13, from 5.0e-17: the SAME 9.8e-19 span over an interval
    # that got 35.6x narrower. Still 5.6e11 times under G9.1's own 1e-3 bar, so the
    # retraction is untouched; the pin moves because its denominator did.
    assert g["per_param_spread"]["tau_y"] < 1e-14, g["per_param_spread"]["tau_y"]
    assert g["param_spread"] < 1e-6, g["param_spread"]
    # and it is NOT a violation of anything any more
    assert ("herschel_bulkley", "G") not in [tuple(v[:2]) for v in ms["G92_violations"]]


def test_the_objective_does_not_move_so_this_is_flatness_not_a_failed_search(ms):
    """If the cost moved, some seeds would be landing on worse optima and the finding
    would be about convergence instead. `rmse_log` is a monotone function of the
    least-squares cost the optimiser minimises, so a bound on one bounds the other."""
    assert ms["worst_cost_spread"] < 1e-8, ms["worst_cost_spread"]
    for fam, donors in ms["records"].items():
        for d, r in donors.items():
            assert r["cost_spread"] < 1e-8, (fam, d, r["cost_spread"])


def test_the_movers_are_the_two_degenerate_carreau_yasuda_records(ms):
    """Both survive the metric change, and both carry a `degenerate` /
    `not_identifiable` verdict -- which is the positive content G9.2 was written for."""
    movers = [(f, d) for f, ds in ms["records"].items() for d, r in ds.items()
              if r["param_spread"] > ms["G92_bar"]]
    assert sorted(movers) == [("carreau_yasuda", "B"), ("carreau_yasuda", "L")]
    assert ms["n_movers"] == 2
    for f, d in movers:
        assert ms["records"][f][d]["has_bad_verdict"], (f, d)


def test_verdicts_are_seed_dependent_on_two_records_and_are_read_over_all_seeds(ms):
    """A real finding the review predicted. The verdict SET changes with the seed on 3 of
    80 pairs, so deciding which criterion a record falls under from seed 0 alone would
    make the populations a function of seed 0's luck. Both `all_measurable` and
    `has_bad_verdict` are therefore taken over all five seeds."""
    # Was three pairs until NEXT_STEPS item 13. Carreau-Yasuda donor A became
    # seed-STABLE once its intervals stopped being bisection-grid points, so the
    # correction removed a seed dependence rather than creating one. B and L remain, and
    # they are the two records whose profile is genuinely discontinuous (item 13 rejects
    # 5 of their interval ends as jumps rather than crossings).
    assert sorted(tuple(v) for v in ms["verdict_unstable"]) == [
        ("carreau_yasuda", "B"), ("carreau_yasuda", "L")]
    for f, d in [tuple(v) for v in ms["verdict_unstable"]]:
        r = ms["records"][f][d]
        assert not r["verdict_seed_stable"]
        assert len({tuple(sorted(v.items())) for v in r["verdicts_by_seed"].values()}) > 1
    # the unstable ones all still land on the same side of the export check
    for f, d in [tuple(v) for v in ms["verdict_unstable"]]:
        assert ms["records"][f][d]["has_bad_verdict"], (f, d)


def test_the_md_quotes_its_own_json_and_carries_the_retraction(ms):
    md = MD.read_text(encoding="utf-8")
    assert "**FAILED**" not in md, "all three criteria pass now"
    assert "## Retraction" in md and "415 %" in md
    assert f"{ms['worst_cost_spread']:.2e}" in md
    assert f"{ms['worst_spread_all_measurable']:.3e}" in md
    assert f"{ms['worst_spread_over_mean']:.3e}" in md


@pytest.mark.slow
def test_rerunning_reproduces_the_verdicts():
    """The fits are seeded, so a re-run must be identical. ~20 minutes.

    NON-DESTRUCTIVE. The script writes over the committed JSON, so the original is saved
    and RESTORED before any assertion can fail -- an earlier version compared and then
    left the artifact overwritten, so a failure destroyed the reference it was checking
    against (code review 2026-09-06). The source CSV is committed, so its absence is a
    real error rather than a reason to skip: a silent skip on the only test that can
    fail is not a test."""
    assert WORKBOOK_CSV.exists(), WORKBOOK_CSV
    original_json = JSON.read_text()
    original_md = MD.read_text(encoding="utf-8")
    try:
        r = subprocess.run([sys.executable,
                            str(ROOT / "scripts" / "multistart_seed_scatter.py")],
                           capture_output=True, text=True, cwd=ROOT)
        rc, err = r.returncode, r.stderr[-2000:]
        after = JSON.read_text()
    finally:
        JSON.write_text(original_json, encoding="utf-8", newline="")
        MD.write_text(original_md, encoding="utf-8", newline="")
    assert rc == 0, err
    assert json.loads(after) == json.loads(original_json)


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


@pytest.mark.slow
def test_carbopol_rerunning_reproduces_the_verdicts():
    """Same shape, same non-destructive discipline. A few minutes."""
    original_json = CARBOPOL_JSON.read_text()
    original_md = CARBOPOL_MD.read_text(encoding="utf-8")
    try:
        r = subprocess.run(
            [sys.executable,
             str(ROOT / "scripts" / "multistart_seed_scatter_carbopol.py")],
            capture_output=True, text=True, cwd=ROOT)
        rc, err = r.returncode, r.stderr[-2000:]
        after = CARBOPOL_JSON.read_text()
    finally:
        CARBOPOL_JSON.write_text(original_json, encoding="utf-8", newline="")
        CARBOPOL_MD.write_text(original_md, encoding="utf-8", newline="")
    assert rc == 0, err
    assert json.loads(after) == json.loads(original_json)


def test_carbopol_md_quotes_its_own_json(cms):
    md = CARBOPOL_MD.read_text(encoding="utf-8")
    assert "**FAILED**" not in md
    assert f"{cms['worst_spread_all_measurable']:.3e}" in md
    assert f"{cms['worst_cost_spread']:.2e}" in md
