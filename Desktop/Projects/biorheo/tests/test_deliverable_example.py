"""The four worked deliverable pages in `case_studies/deliverable_example/`.

These are sales artifacts, so the thing worth testing is not that they were produced but
that they still **demonstrate four different things**. A set of examples that all say the
same thing is a brochure, not evidence, and it would drift into that state silently the
first time a record or a bar changed.

The other half is the drift guard this repo keeps having to add: every number the README
quotes must equal the one in `summary.json`, because a page whose prose contradicts its
own artifact has shipped here more than once.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
CASE = REPO / "case_studies" / "deliverable_example"
SUMMARY = CASE / "summary.json"
README = CASE / "README.md"
SLUGS = ("01_clean", "02_open", "03_unbounded", "04_thinning", "05_screen")
# 01-04 are built from registry records; 05 is SCREEN data and cannot be one.
RECORD_SLUGS = SLUGS[:4]
SCREEN_SLUG = "05_screen"


@pytest.fixture(scope="module")
def s() -> dict:
    assert SUMMARY.exists(), "run case_studies/deliverable_example/run.py"
    return json.loads(SUMMARY.read_text(encoding="utf-8"))


def test_all_five_pages_exist(s):
    assert tuple(s) == SLUGS
    for slug in SLUGS:
        assert (CASE / f"{slug}.md").exists(), slug


def test_the_four_pages_demonstrate_four_different_things(s):
    """The point of the set: no page may be redundant.

    **This test first asserted four distinct FLAG SETS and failed at three**, and the
    failure was in the criterion, not the pages. `01_clean` and `04_thinning` carry the
    same flags (`non_commercial, extrapolated`) and differ on a different axis entirely --
    the Newtonian `R^-4` verdict, safe on one and not on the other. The demonstration axis
    is therefore (flags, bore verdict) and the test says so rather than being weakened to
    three.
    """
    axes = [(tuple(sorted(v["flags"])), v["newtonian_rule_safe_within_10pct"])
            for v in s.values()]
    assert len(set(axes)) == 5, axes
    assert len({a[0] for a in axes}) == 4, "four distinct flag sets, by construction"


def test_each_named_behaviour_is_actually_present(s):
    assert "open_interval" not in s["01_clean"]["flags"]
    assert "open_interval" in s["02_open"]["flags"]
    assert "hull_unbounded_above" in s["03_unbounded"]["flags"]
    assert "hull_unbounded_above" not in s["02_open"]["flags"]


def test_the_newtonian_rule_verdict_splits(s):
    """Both verdicts must appear. If every page said the rule was safe, the section would
    be passed by a routine that always says so -- the failure mode this repo calls 'both
    candidates pass'."""
    safe = [k for k, v in s.items() if v["newtonian_rule_safe_within_10pct"]]
    unsafe = [k for k, v in s.items() if not v["newtonian_rule_safe_within_10pct"]]
    assert len(safe) == 4 and unsafe == ["04_thinning"]
    assert s["04_thinning"]["bore_exponent"] > -3.0
    assert all(s[k]["bore_exponent"] < -3.8 for k in safe)
    # the screen page is a Newtonian fluid, so its exponent must be EXACTLY -4 -- item 18's
    # closed form, and a check that the screen path did not quietly change the physics
    assert s[SCREEN_SLUG]["bore_exponent"] == pytest.approx(-4.0, abs=1e-9)


def test_the_interval_widths_span_four_orders(s):
    """A clean fit and a degenerate one must not produce similar-looking error bars, or
    the identifiability verdict on the page is decorative."""
    w = {k: v["hull_relative_width"] for k, v in s.items()}
    assert w["01_clean"] < 1.0
    assert w["03_unbounded"] > 1000.0
    assert w["03_unbounded"] / w["04_thinning"] > 1e4


def test_the_licence_split_is_exactly_where_it_should_be(s):
    """This test used to assert every page was non-commercial, with the note "this fails
    the day a licensable record is added, which is the day the demo set should be rebuilt
    on it". **2026-09-08 was that day** (item 25).

    The four record-backed pages are still built on non-commercial data and must still be
    stamped. The screen page is CC BY 4.0 and must NOT be -- and it must be the one on the
    wedge fluid, or the set has a licensable page that demonstrates nothing new.
    """
    for slug in RECORD_SLUGS:
        assert s[slug]["commercial_use"] is False, slug
        assert "non_commercial" in s[slug]["flags"], slug
        assert "NOT LICENSABLE" in (CASE / f"{slug}.md").read_text(encoding="utf-8"), slug
    assert s[SCREEN_SLUG]["commercial_use"] is True
    assert "non_commercial" not in s[SCREEN_SLUG]["flags"]
    screen = (CASE / f"{SCREEN_SLUG}.md").read_text(encoding="utf-8")
    assert "NOT LICENSABLE" not in screen
    assert "CC BY 4.0" in screen
    assert "monoclonal antibody" in screen


def test_every_page_is_an_extrapolation_and_says_so(s):
    """Every page extrapolates -- but no longer all by more than a decade.

    This asserted `> 1.0` on every page until item 25 added the screen page, whose fluid
    was measured at 960-3329 1/s on a microfluidic viscometer and therefore extrapolates
    only 0.74 decades. That is the SMALLEST in the set and it is the point of that dataset,
    so the bar is restated as "every page is an extrapolation" and the screen page's rank
    is asserted rather than the old blanket threshold being kept or dropped.
    """
    assert all(v["decades_past_envelope"] > 0.0 for v in s.values())
    for slug in SLUGS:
        assert "decades outside the measured range" in (CASE / f"{slug}.md").read_text(
            encoding="utf-8"), slug
    others = [v["decades_past_envelope"] for k, v in s.items() if k != SCREEN_SLUG]
    assert all(x > 1.0 for x in others), "the four record pages all still exceed a decade"
    assert s[SCREEN_SLUG]["decades_past_envelope"] < min(others), (
        "the screen page must be the least-extrapolated of the set; that is why that "
        "dataset was chosen"
    )
    assert s[SCREEN_SLUG]["decades_past_envelope"] == pytest.approx(0.74, abs=0.02)


def test_the_readme_numbers_match_the_artifact(s):
    """The drift guard. Each number the README quotes is recomputed here from
    `summary.json` and required to appear verbatim in the prose."""
    text = README.read_text(encoding="utf-8")
    lo = min(v["decades_past_envelope"] for v in s.values())
    hi = max(v["decades_past_envelope"] for v in s.values())
    n_safe = sum(1 for v in s.values() if v["newtonian_rule_safe_within_10pct"])
    assert f"{lo:.2f} to {hi:.2f} decades" in text
    assert f"It is on {n_safe} of {len(s)} pages" in text
    assert f"{s['04_thinning']['bore_exponent']:.4g}" in text


def test_the_readme_links_every_page(s):
    text = README.read_text(encoding="utf-8")
    for slug in SLUGS:
        assert f"[`{slug}.md`]({slug}.md)" in text, slug


@pytest.mark.slow
def test_regenerating_the_pages_reproduces_them(tmp_path):
    """RESTORES every artifact before asserting, so a failure cannot destroy the reference
    it was checking against."""
    import shutil
    import subprocess
    import sys
    files = [SUMMARY, README] + [CASE / f"{x}.md" for x in SLUGS]
    saved = {p: tmp_path / p.name for p in files}
    for src, dst in saved.items():
        shutil.copy2(src, dst)
    try:
        r = subprocess.run([sys.executable, str(CASE / "run.py")], capture_output=True,
                           text=True, timeout=3600)
        assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-2000:]
        new = json.loads(SUMMARY.read_text(encoding="utf-8"))
        old = json.loads((tmp_path / SUMMARY.name).read_text(encoding="utf-8"))
    finally:
        for src, dst in saved.items():
            shutil.copy2(dst, src)
    for slug in SLUGS:
        assert new[slug]["flags"] == old[slug]["flags"], slug
        assert new[slug]["force_N"] == pytest.approx(old[slug]["force_N"], rel=1e-9), slug
