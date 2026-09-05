"""Unit 0.1 gate: dst_fit is consumed unchanged.

Its suite is green when run from here, and its one figure is byte-identical to the
hash recorded on 2026-09-05 (246 tests green that day). The gate is the *absence* of
change, so the numbers are pinned rather than recomputed.
"""

import hashlib
import subprocess
import sys
from pathlib import Path

import pytest

DST_FIT = Path(r"C:\Users\compl\Desktop\Projects\dst_fit")
FIG = DST_FIT / "figures" / "fig_phim_drift.png"
FIG_MD5 = "e9d3c0430ce5b8fd09ba463901062218"  # md5 on 2026-09-05


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
