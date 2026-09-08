"""The retracted 12 % must not survive anywhere in the pipeline (2026-09-08).

`cfd_fit_summary.json`'s `grid_independence` block shipped `C_P` = 1.5178 and a 12.28 %
change from a solve that stopped at its iteration cap. That value reached `run_bidet.py`'s
printed output and the manuscript's grid-independence sentence. These tests assert it is
gone from all three, that the block is regenerable from the GCI artifact rather than
hand-written, and that what replaced it does not read as an error estimate.
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SUMMARY = ROOT / "cfd_openfoam" / "cfd_fit_summary.json"
GCI = ROOT / "cfd_openfoam" / "convergence" / "gci_Re10000_H6_fine_conv.csv"
TEX = ROOT / "bidet_optimization_1.JOURNAL.tex"
RUNNER = ROOT / "model" / "run_bidet.py"


@pytest.fixture(scope="module")
def gi() -> dict:
    return json.loads(SUMMARY.read_text(encoding="utf-8"))["grid_independence"]


def test_the_block_matches_the_gci_artifact_it_claims_to_come_from(gi):
    """Not hand-written: every value must be readable out of the GCI CSV."""
    rows = {r["quantity"]: r for r in csv.DictReader(GCI.open(encoding="utf-8"))}
    for key, q in (("n_t", "n_t"), ("C_P", "C_P"), ("tau15", "tau_at_15mm_Pa"),
                   ("tau_peak", "tau_peak_Pa")):
        assert gi["base_refine_1"][key] == pytest.approx(float(rows[q]["phi3"]))
        assert gi["refined_x2"][key] == pytest.approx(float(rows[q]["phi2"]))
        assert gi["refined_x4"][key] == pytest.approx(float(rows[q]["phi1"]))


def test_the_regenerator_reports_the_block_as_up_to_date():
    r = subprocess.run([sys.executable,
                        str(ROOT / "cfd_openfoam" / "refresh_grid_block.py"), "--check"],
                       capture_output=True, text=True, cwd=ROOT / "cfd_openfoam")
    assert r.returncode == 0, r.stdout + r.stderr


def test_the_retracted_value_is_recorded_not_erased(gi):
    """A retraction whose old number vanishes cannot be checked by the next reader."""
    r = gi["RETRACTED_2026-09-03"]
    assert r["old_C_P"] == 1.5178
    assert r["old_rel_change_C_P"] == pytest.approx(0.1228, abs=1e-4)
    assert gi["refined_x2"]["C_P"] != r["old_C_P"]


def test_the_current_x2_change_is_3_percent_not_12(gi):
    assert gi["rel_change"]["C_P"] == pytest.approx(0.0303, abs=5e-4)
    assert gi["rel_change_x4"]["C_P"] == pytest.approx(0.1695, abs=5e-4)


def test_the_block_carries_the_divergence_warning(gi):
    """These percentages are sensitivities. Without the warning they read as error bars."""
    w = gi["WARNING"]
    assert "DIVERGES" in w and "NOT a convergence study" in w
    assert "never as an error estimate" in w


def test_the_runner_no_longer_calls_it_grid_independence():
    src = RUNNER.read_text(encoding="utf-8")
    assert "grid-independence:" not in src, "the printed label was retired 2026-09-08"
    assert "SENSITIVITY" in src and "DIVERGENCE" in src


def test_the_manuscript_no_longer_quotes_12_percent():
    tex = TEX.read_text(encoding="utf-8")
    assert r"$C_P$ by $12\%$" not in tex
    assert r"+3.0\%" in tex and r"+17.0\%" in tex
    assert "diverges under" in tex
    # the retraction itself must be visible in the paper, not silently swapped
    assert r"quoted $12\%$ for $C_P$" in tex


def test_the_manuscript_cites_the_procedure_it_invokes():
    tex = TEX.read_text(encoding="utf-8")
    assert r"\citet{Celik2008}" in tex
    assert r"\bibitem[Celik et~al.(2008)]{Celik2008}" in tex
    assert "10.1115/1.2960953" in tex, "the DOI was resolved through Crossref before use"
