"""Unit 0.1 gate: dst_fit is consumed unchanged.

Its suite is green when run from here, and its one figure is byte-identical to the
hash recorded on 2026-09-05 (246 tests green that day). The gate is the *absence* of
change, so the numbers are pinned rather than recomputed.

One deliberate change is pinned rather than forbidden: on 2026-09-26 the founder
retired dst_fit's gate G4c (decision D5; biorheo item 29 / work order D1). The
expected output of `scripts/gate_g4.py` therefore moved from exit 1 ("G4c: 1 of 6
within 2 sigma, needed 4") to exit 0 with the same 1 of 6 printed and labelled
RETIRED. `test_gate_g4_retired_g4c_keeps_its_number` pins BOTH halves: the exit code
AND the retained number, so a later edit that quietly fixes the data or the bar to
make G4c pass would fail this test as surely as one that dropped the number.
"""

import hashlib
import subprocess
import sys
from pathlib import Path

import pytest

DST_FIT = Path(r"C:\Users\compl\Desktop\Projects\dst_fit")
FIG = DST_FIT / "figures" / "fig_phim_drift.png"
FIG_MD5 = "e9d3c0430ce5b8fd09ba463901062218"  # md5 on 2026-09-05
GATE_G4 = DST_FIT / "scripts" / "gate_g4.py"
GUY_DEPOSIT = DST_FIT / "Towards a unified description of the rheology of hard-particle suspensions"
# G4c's own reading on 2026-08-27 and again on 2026-09-28 (script-printed, both runs):
G4C_RETAINED = "reads 1 of 6 within 2 sigma, which would have been FAIL"


def test_dst_fit_importable_at_pinned_version():
    import dst_fit

    assert dst_fit.__version__ == "0.1.0"


def test_dst_fit_source_is_the_sibling_checkout():
    import dst_fit

    assert Path(dst_fit.__file__).resolve().is_relative_to(DST_FIT)


@pytest.mark.slow
def test_dst_fit_suite_green():
    r = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider"],
        cwd=DST_FIT, capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stdout[-2000:]


def test_figure_hash_unchanged():
    assert FIG.exists(), FIG
    assert hashlib.md5(FIG.read_bytes()).hexdigest() == FIG_MD5


def test_gate_g4_source_marks_g4c_retired_without_loosening_the_bar():
    """Static half: the retirement is in the source, and the bar constants are intact."""
    src = GATE_G4.read_text(encoding="utf-8")
    assert "G4c  RETIRED 2026-09-26" in src
    assert "G4C_RETIRED_REASON" in src
    assert "G4_MIN_CONSISTENT = 4" in src  # the bar, kept verbatim
    assert "G4_N_SIGMA = 2.0" in src
    assert "G4_TOLERANCE = 0.02" in src


def test_gate_g4_retired_g4c_keeps_its_number():
    """Dynamic half (~6 s): exit 0 on G4a + G4b, and G4c's 1 of 6 is still printed.

    Needs the Guy et al. 2015 DataShare deposit on disk; it is committed data on this
    machine, so its absence is a failure, not a skip.
    """
    assert GUY_DEPOSIT.is_dir(), GUY_DEPOSIT
    r = subprocess.run(
        [sys.executable, "-X", "utf8", str(GATE_G4)],
        cwd=DST_FIT, capture_output=True, text=True, encoding="utf-8",
    )
    assert r.returncode == 0, r.stdout[-1500:]
    assert "G4a PASS   G4b PASS   G4c RETIRED (would have read FAIL)" in r.stdout
    assert G4C_RETAINED in r.stdout, r.stdout[-1500:]
    # the retired gate must never become a failure again by a different name
    assert "GATE G4 FAIL" not in r.stdout
