"""test_sweep.py -- regression tests for the 2026-09-18 descriptive sweep:

  * sweep_descriptive.py / out/sweep_descriptive.txt
  * REPORT_SWEEP_2026-09-18.md

STRUCTURAL SCOPE GUARD. Imports no project module at all; every claim below is read from
the committed artifact and re-checked with arithmetic implemented here. Asserted on the
import lines, as `test_additivity.py` and `test_fabric.py` do.

WHAT THESE EXIST TO CATCH. The first version of this block was refuted by an adversarial
referee on four counts, and each one is pinned here so it cannot come back: a 2*se
significance guard on 3 seeds (an 18 % test), an eyeball verdict where a heterogeneity
test belongs, an endpoint-only change definition on a family that is not monotone, and a
six-way ordering quoted as "the part worth carrying" when it survives a third of
resamples.
"""

from __future__ import annotations

import os
import re

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
SW = os.path.join(OUT, "sweep_descriptive.txt")
REPORT = os.path.join(HERE, "REPORT_SWEEP_2026-09-18.md")


def test_this_file_imports_no_project_module():
    """Parsed with `ast`, not matched with a regex on raw lines.

    A line-regex guard reads PROSE: the sibling guards match any line beginning with
    `import`, and the phrase "import lines" in this very docstring made the first version
    of this test fail on a word. Parsing the module gives the real import roots.
    """
    import ast
    src = open(os.path.join(HERE, "test_sweep.py"), encoding="utf-8").read()
    roots = set()
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.Import):
            roots.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            roots.add(node.module.split(".")[0])
    allowed = {"os", "re", "pytest", "ast", "__future__"}
    assert roots <= allowed, f"scope guard: unexpected import roots {roots - allowed}"


def _txt(p):
    assert os.path.exists(p), f"committed artifact missing: {p}"
    return open(p, encoding="utf-8").read()


@pytest.fixture(scope="module")
def sw():
    return _txt(SW)


@pytest.fixture(scope="module")
def rep():
    return _txt(REPORT)


def test_the_reproduction_gate_ran_first_and_passed(sw):
    assert sw.index("1. REPRODUCTION GATE") < sw.index("4. THE SHAPE-OVER-CV CONTRAST")
    assert "REPRODUCTION: PASS" in sw


def test_the_significance_guard_is_a_t_value_not_two_se(sw):
    """3 seeds means 2 dof and t = 4.303. A flat 2*se bar is an 18 % test, and it
    admitted exactly the change that carried the first draft's headline."""
    assert "t(df = n-1), NOT 2*se" in sw
    assert "4.303" in sw
    assert re.search(r"UNDETERMINED .*does not clear t=4\.30", sw), \
        "no ratio was rejected by the corrected guard; the guard may have been loosened"


def test_every_ranked_channel_is_marked_resolved_or_not(rep):
    """The report ranks every channel, including ones that fail a gate. It must SAY so
    where it ranks them -- the prose-contradicts-its-own-table failure."""
    m = re.search(r"Ranked by across-family spread.*?\n", rep, flags=re.S)
    assert m, "no ordering line"
    line = m.group(0)
    assert line.count("*(resolved)*") + line.count("*(not resolved)*") == 6, line
    assert "not resolved" in line


def test_BOTH_gates_are_applied_not_just_the_denominator(sw):
    """The denominator test asks whether the smallest family change is non-zero; the
    ordering claim needs whether the families DIFFER FROM EACH OTHER. A channel that
    passes one and fails the other must not be quoted."""
    blk = re.search(r"FINAL STATUS PER CHANNEL.*?\n\n", sw, flags=re.S)
    assert blk, "the two-gate status block is missing"
    rows = re.findall(r"^\s{5}(\S[\w ()=]*?)\s{2,}(PASS|FAIL)\s+(PASS|FAIL)\s+(.*)$",
                      blk.group(0), flags=re.M)
    assert len(rows) == 6, rows
    for nm, denom, het, verdict in rows:
        resolved = "RESOLVED" in verdict and "NOT RESOLVED" not in verdict
        assert resolved == (denom == "PASS" and het == "PASS"), (nm, denom, het, verdict)
    # and the gates must actually bite: a run where everything passes both would mean
    # the guards had stopped discriminating.
    assert any(d == "FAIL" or h == "FAIL" for _, d, h, _ in rows), rows


def test_the_conservative_heterogeneity_test_is_COMPUTED_not_described(sw, rep):
    """Q treats the se's as exactly known and is anti-conservative at 2 dof. An earlier
    version printed Q alone and put the conservative end in prose, so its column read as
    six unanimous rejections of family-independence."""
    assert "Cochran's Q" in sw and "ANOVA F" in sw
    assert "p(F)" in rep and "ANOVA F" in rep
    rows = re.findall(r"^\s{5}(\S[\w ()=]*?)\s+(-?[0-9.]+)\s+([0-9.]+)\s+"
                      r"(-?[0-9.]+)\s+([0-9.]+)\s+(YES|no) / (YES|no)\s*$",
                      sw, flags=re.M)
    assert len(rows) == 6, rows
    # The two tests must actually disagree somewhere, or the conservative one is not
    # doing any work and the earlier prose would have been adequate after all.
    assert any(q != f for *_, q, f in rows), rows


def test_the_endpoint_definition_is_shown_to_do_work(sw):
    """CV 0 -> 0.5 assumes monotonicity and one family is not monotone. The peak-based
    ratio must be printed beside the endpoint one, and must actually differ somewhere."""
    blk = re.search(r"ENDPOINT vs PEAK.*?IS THE ORDERING", sw, flags=re.S)
    assert blk, "the endpoint-vs-peak block is missing"
    rows = re.findall(r"^\s{5}(\S[\w ()=]*?)\s+([0-9.]+)x\s+([0-9.]+)x",
                      blk.group(0), flags=re.M)
    assert len(rows) == 6, ("every channel must be checked, including phi_J, which the "
                            "first version skipped because it has no shard getter", rows)
    assert any(r[0].strip() == "phi_J" for r in rows), rows
    assert any(abs(float(a) - float(b)) > 0.05 for _, a, b in rows), \
        "endpoint and peak agree everywhere; the check has stopped discriminating"


def test_the_peak_ratio_is_disclosed_as_biased_high(sw, rep):
    """It is a maximum over five noisy estimates, each family at its own best rung --
    an upper bound, biased the opposite way from the endpoint definition."""
    assert "winner's-curse" in sw or "winner" in sw
    assert "upper bound" in sw.lower()
    assert "upper bounds biased" in rep or "biased" in rep


def test_the_ordering_stability_is_measured_and_the_report_does_not_oversell_it(sw, rep):
    m = re.search(r"the exact six-way ordering survives (\d+) % of (\d+) draws", sw)
    assert m, "the ordering stability was not measured"
    frac = int(m.group(1))
    assert int(m.group(2)) >= 1000
    # The claim in the report must match the measurement, whichever way it came out.
    m2 = re.search(r"survives only (\d+)% of", rep)
    assert m2 and int(m2.group(1)) == frac, (frac, m2.group(1) if m2 else None)
    if frac < 50:
        assert "It is not a result" in rep


def test_Z_all_is_disclosed_as_not_independent(sw, rep):
    """Z_all = (1 - rattler_frac) * Z_nr to better than 0.5 %, so it is not a second
    confirmation of the rattler channel."""
    assert "Z_all IS NOT AN INDEPENDENT CHANNEL" in sw
    m = re.search(r"worst \|difference\| over the rows printed: ([0-9.]+)", sw)
    assert m, "the identity check printed no residual"
    assert float(m.group(1)) < 0.05, m.group(1)
    assert "not an independent channel" in rep.lower()


def test_every_filter_prints_its_N_and_the_duplicate_conflicts_are_surfaced(sw):
    ns = re.findall(r"N (?:records|after [^:]+): (\d+)", sw)
    assert len(ns) >= 3, ns
    assert [int(x) for x in ns] == sorted((int(x) for x in ns), reverse=True), ns
    assert re.search(r"duplicate keys whose phi DISAGREED: \d+", sw)
