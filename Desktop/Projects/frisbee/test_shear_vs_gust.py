"""The shear-vs-gust moment ratio, and the two retractions behind it.

This quantity has been stated three ways and two of them are retracted:

    2026-09-03 (retired)  sigma_u / DeltaU_shear                    37-122x
    2026-09-03 (retired)  sigma_u (d/L) / DeltaU_shear               1.7-4.9x
    2026-09-04 (current)  sqrt(C2) sigma_u (d/L)^(1/3) / DeltaU_shear  19-59x

The tests below pin the current numbers, check the REGIME TEST that selects between the
last two forms, and assert that the script still prints both retired values -- a
retraction whose old number stops being visible is not a retraction.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "shear_vs_gust_moment.py"
DOC = HERE / "FRISBEE_GENERALIZATION.md"


@pytest.fixture(scope="module")
def out() -> str:
    r = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True,
                       cwd=HERE)
    assert r.returncode == 0, r.stderr[-2000:]
    return r.stdout


def _ratios(out: str) -> dict[str, tuple[float, float, float]]:
    """theta -> (ratio_used, ratio_smooth, ratio_pointrms) from the summary block."""
    rows = {}
    for m in re.finditer(r"^\s*(90|210|270)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s*$",
                         out, re.M):
        rows[m.group(1)] = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
    return rows


def test_all_three_headings_are_reported(out):
    assert set(_ratios(out)) == {"90", "210", "270"}


def test_the_regime_test_puts_the_disc_in_the_inertial_subrange(out):
    """The whole 2026-09-04 correction turns on d/eta, so it must be printed and large.

    This is the check the first version did not do: it asserted the smooth-flow scaling
    in prose instead of computing the number that decides between the two.
    """
    m = re.search(r"d/eta = (\d+) to (\d+), against the inertial-subrange "
                  r"threshold of (\d+)", out)
    assert m, "the script must print d/eta beside the verdict"
    lo, hi, thresh = int(m.group(1)), int(m.group(2)), int(m.group(3))
    assert lo == 315 and hi == 468
    assert lo > 5 * thresh, "the regime call must not be marginal"
    assert "FIRMLY inertial-subrange" in out


def test_the_current_ratio_is_19_to_59(out):
    used = [v[0] for v in _ratios(out).values()]
    assert min(used) == pytest.approx(18.7, abs=0.3)
    assert max(used) == pytest.approx(59.0, abs=0.5)


def test_both_retracted_forms_are_still_printed(out):
    """A retraction whose old number disappears is not checkable by the next reader."""
    r = _ratios(out)
    smooth = [v[1] for v in r.values()]
    pointrms = [v[2] for v in r.values()]
    assert min(smooth) == pytest.approx(1.7, abs=0.1)
    assert max(smooth) == pytest.approx(4.9, abs=0.1)
    assert min(pointrms) == pytest.approx(37, abs=1)
    assert max(pointrms) == pytest.approx(122, abs=1)


def test_the_smooth_limit_understates_by_about_L_over_d_to_the_two_thirds(out):
    """(L/d)^(2/3) with d/L ~ 0.043 is ~8-12; the printed factor must match that scaling.

    A discriminating check rather than a restatement: if the script had kept the smooth
    form and only relabelled it, this ratio would be exactly 1.
    """
    m = re.search(r"understates the inertial one by ([\d.]+) to ([\d.]+) x", out)
    assert m
    lo, hi = float(m.group(1)), float(m.group(2))
    assert 8.0 < lo < 13.0 and 8.0 < hi < 13.0
    assert hi > lo


def test_the_null_survives_and_the_script_says_so_from_the_number(out):
    """The verdict text is SELECTED BY the computed ratio, not written beside it."""
    assert "The null SURVIVES with room to spare" in out
    used = [v[0] for v in _ratios(out).values()]
    assert min(used) >= 10.0, "that verdict branch requires >= 10x at every heading"


def test_the_open_fetch_margin_is_stated_and_is_single_figure(out):
    """Caveat 1 is the one that constrains what may be reported; it must carry a number."""
    m = re.search(r"leaves gust ahead by only\n([\d.]+)-([\d.]+) x there", out)
    assert m, "caveat 1 must state the open-fetch margin"
    assert float(m.group(1)) < 10.0


def test_the_document_quotes_what_the_script_prints(out):
    """The 2026-09-03 failure mode here was prose that disagreed with its own artifact:
    the docstring described the inertial form while the code ran the smooth one."""
    doc = DOC.read_text(encoding="utf-8")
    assert "19 to 59×" in doc
    assert "315 to 468" in doc
    for theta, ratio in (("90°", "26×"), ("210°", "19×"), ("270°", "59×")):
        assert f"**{ratio}**" in doc, (theta, ratio)
    # and the retracted numbers must still be visible, struck through
    assert "~~1.7×~~" in doc and "~~122×~~" in doc


def test_the_hardening_brief_is_marked_as_corrected():
    """The brief's '~20x' lands near the right answer for two wrong reasons at once."""
    brief = (HERE / "CLAUDE_CODE_FRISBEE_HARDENING.md").read_text(encoding="utf-8")
    i = brief.index("~20× the shear-induced one.")
    window = brief[i:i + 700]
    assert "CORRECTED 2026-09-04" in window
    assert "19 to 59×" in window
