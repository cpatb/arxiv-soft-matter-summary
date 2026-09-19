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
    assert blk, "the gate status block is missing"
    rows = re.findall(
        r"^\s{5}(\S[\w ()=]*?)\s{2,}(PASS|FAIL)\s+(PASS|FAIL)\s+(PASS|FAIL)\s+(.*)$",
        blk.group(0), flags=re.M)
    assert len(rows) == 6, rows
    for nm, denom, het, pair, verdict in rows:
        resolved = "RESOLVED" in verdict and "NOT RESOLVED" not in verdict
        assert resolved == (denom == "PASS" and het == "PASS" and pair == "PASS"), \
            (nm, denom, het, pair, verdict)
    # and each gate must actually bite somewhere: a run where one of them never fails
    # is a gate that has stopped discriminating.
    for j, gate in ((1, "denominator"), (2, "heterogeneity"), (3, "quoted pair")):
        assert any(r[j] == "FAIL" for r in rows), (gate, rows)


def test_the_quoted_pair_is_tested_against_a_SELECTION_CALIBRATED_threshold(sw):
    """An omnibus F says SOME pair differs; the ratio quotes ONE pair. And that pair is
    the argmax/argmin OF THE DATA, so a nominal two-sample t_crit under-corrects -- the
    same defect as 2*se on 3 seeds, one level up. The measured false-positive rate of
    the nominal gate must be printed, and the verdict must be read off the CALIBRATED
    threshold (referee, fourth pass)."""
    blk = re.search(r"A SIGNIFICANT OMNIBUS F.*?FINAL STATUS", sw, flags=re.S)
    assert blk, "the pairwise block is missing"
    rows = re.findall(r"^\s{5}(\S[\w ()=]*?)\s{2,}(\w+)\s+(\w+)\s+([0-9.]+)\s+"
                      r"([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)%\s+(YES|no)\s*$",
                      blk.group(0), flags=re.M)
    assert len(rows) == 6, rows
    for nm, fmax, fmin, t, tnom, tcal, fpr, ok in rows:
        assert fmax != fmin, nm
        # the verdict is read off the CALIBRATED threshold, not the nominal one
        assert (float(t) > float(tcal)) == (ok == "YES"), (nm, t, tcal, ok)
        # the calibration must actually be doing something
        assert float(tcal) > float(tnom), (nm, tnom, tcal)
        assert float(fpr) > 5.0, (nm, "the nominal gate's FPR should exceed 5 %", fpr)
    # and it must bite: at least one channel whose nominal verdict the calibration flips
    assert any(float(t) > float(tnom) and float(t) <= float(tcal)
               for _, _, _, t, tnom, tcal, _, _ in rows), rows


def test_the_heterogeneity_F_is_recomputable_from_the_printed_digits(sw):
    """The property that actually matters, tested the way the referee tested it.

    `contact force mean` is O(1e-3). At four decimals its between-family differences sat
    inside the rounding quantum, so reconstructing its heterogeneity F from the printed
    table gave p = 0.118 where the full-precision value is 0.0248 -- opposite sides of
    5 %, on a row the report calls RESOLVED. A claims table must be recomputable from
    what it prints. The ANOVA is re-implemented here from scratch; nothing is imported.
    """
    chg = re.findall(
        r"^\s{3}(\S[\w ()=]*?)\s{2,}([-+0-9.eE]+)\+-([0-9.eE+-]+)\s+"
        r"([-+0-9.eE]+)\+-([0-9.eE+-]+)\s+([-+0-9.eE]+)\+-([0-9.eE+-]+)\s+n=(\d+)/(\d+)/(\d+)",
        sw, flags=re.M)
    assert len(chg) == 6, chg
    printed = dict(re.findall(
        r"^\s{5}(\S[\w ()=]*?)\s+-?[0-9.]+\s+[0-9.]+\s+(-?[0-9.]+)\s+[0-9.]+\s+"
        r"(?:YES|no) / (?:YES|no)\s*$", sw, flags=re.M))
    assert len(printed) == 6, printed

    def f_from(ms, ses, ns):
        sd = [s * (n ** 0.5) for s, n in zip(ses, ns)]
        ssw = sum((n - 1) * d * d for n, d in zip(ns, sd))
        dfw = sum(n - 1 for n in ns)
        mbar = sum(n * m for n, m in zip(ns, ms)) / sum(ns)
        ssb = sum(n * (m - mbar) ** 2 for n, m in zip(ns, ms))
        return (ssb / (len(ms) - 1)) / (ssw / dfw)

    for nm, m1, s1, m2, s2, m3, s3, n1, n2, n3 in chg:
        nm = nm.strip()
        f_hat = f_from([float(m1), float(m2), float(m3)],
                       [float(s1), float(s2), float(s3)],
                       [int(n1), int(n2), int(n3)])
        f_ref = float(printed[nm])
        assert abs(f_hat - f_ref) <= 0.05 * f_ref, (nm, f_hat, f_ref)
        # A 5 % band in F is a PROXY for "the verdict does not change", and near the
        # F-critical value it is not a safe one. Require the same side of it too: a
        # future borderline channel could pass the band and flip its verdict silently
        # (referee, fourth pass). F_crit(2, 6) at 5 % = 5.1433.
        fcrit = 5.1433
        assert (f_hat > fcrit) == (f_ref > fcrit), (nm, f_hat, f_ref, fcrit)


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
    """`Z_all = (1 - rattler_frac) * Z_nr`, so it is not a second confirmation of the
    rattler channel.

    ⚠ 2026-09-19: this docstring said "to better than 0.5 %" and the assertion bounded
    only the max over printed ROWS (means over seeds at one rung). A row max cannot bound
    a cell: per CELL the residual is ~2.4x larger, and it is the per-cell figure that
    `main.tex` asserts "on every cell". Both are checked now, each against its own bound.
    """
    assert "Z_all IS NOT AN INDEPENDENT CHANNEL" in sw
    m = re.search(r"worst \|difference\| over the rows printed: ([0-9.]+)", sw)
    assert m, "the identity check printed no residual"
    assert float(m.group(1)) < 0.05, m.group(1)
    assert "not an independent channel" in rep.lower()
    # the per-CELL bound, which is the one the manuscript quotes
    mc = re.search(r"to within ([0-9.]+) on every one of the (\d+) cells", rep)
    assert mc, "REPORT_SWEEP no longer prints the per-cell identity bound"
    cell_dev, n_cells = float(mc.group(1)), int(mc.group(2))
    assert cell_dev < 0.05, cell_dev
    # and it must be at least the row max -- a per-cell max below a max over means of
    # those same cells is impossible, so this catches a wrong field or a wrong subset.
    assert cell_dev >= float(m.group(1)), (cell_dev, m.group(1))
    assert n_cells > 100, n_cells


def test_every_filter_prints_its_N_and_the_duplicate_conflicts_are_surfaced(sw):
    ns = re.findall(r"N (?:records|after [^:]+): (\d+)", sw)
    assert len(ns) >= 3, ns
    assert [int(x) for x in ns] == sorted((int(x) for x in ns), reverse=True), ns
    assert re.search(r"duplicate keys whose phi DISAGREED: \d+", sw)
