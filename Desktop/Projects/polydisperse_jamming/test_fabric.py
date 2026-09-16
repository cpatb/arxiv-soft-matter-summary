"""test_fabric.py -- regression tests for the 2026-09-18 fabric block:

  * PREREG_FABRIC_2026-09-18.md  -- the registration
  * grade_fabric.py / out/grade_fabric.txt  -- the graded campaign
  * diag_fabric_null.py / out/diag_fabric_null.txt -- the sizing, closed form and floor

STRUCTURAL SCOPE GUARD. This file must not be able to pass by calling the code it is
checking, and it must not be able to re-decide a gate that is closed. It therefore
imports **no project module at all**: the closed form E[A^2] = 1/n, the T statistic and
the floor arithmetic are all re-implemented here from scratch and compared against the
committed artifacts. The guard is asserted on the import lines, as `test_additivity.py`
and `test_u3decide.py` do.

WHY ARTIFACT + REIMPLEMENTATION RATHER THAN CALLING THE GRADER. `grade_fabric.py` reads
324 shards and Monte-Carlos its floor constant; a slow test that regenerates the artifact
it validates destroys its own reference when it fails. The artifacts are committed, so
their presence is asserted rather than skipped over.
"""

from __future__ import annotations

import math
import os
import re

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
GF = os.path.join(OUT, "grade_fabric.txt")
DN = os.path.join(OUT, "diag_fabric_null.txt")
PREREG = os.path.join(HERE, "PREREG_FABRIC_2026-09-18.md")


# --------------------------------------------------------------------- scope guard
def test_this_file_imports_no_project_module():
    src = open(os.path.join(HERE, "test_fabric.py"), encoding="utf-8").read()
    imports = [ln.strip() for ln in src.splitlines()
               if re.match(r"^\s*(import|from)\s", ln)]
    banned = ("grade_", "diag_", "packing", "measure", "psd", "shape", "phase2",
              "analyse_", "sweep_")
    for ln in imports:
        for bad in banned:
            assert bad not in ln, f"scope guard: '{ln}' pulls in {bad}"
    allowed = {"math", "os", "re", "numpy", "pytest", "__future__"}
    for ln in imports:
        root = re.sub(r"^(from|import)\s+([A-Za-z_][\w.]*).*$", r"\2", ln).split(".")[0]
        assert root in allowed, f"scope guard: unexpected import root {root!r}"


# ------------------------------------------------------------------------ fixtures
def _txt(path):
    assert os.path.exists(path), f"committed artifact missing: {path}"
    return open(path, encoding="utf-8").read()


@pytest.fixture(scope="module")
def gf():
    return _txt(GF)


@pytest.fixture(scope="module")
def dn():
    return _txt(DN)


@pytest.fixture(scope="module")
def prereg():
    return _txt(PREREG)


# ----------------------------------------------- the closed form, re-derived from scratch
def test_the_floor_closed_form_is_one_over_n():
    """E[A^2] = 1/n for n independent uniform normals -- re-implemented here, not
    imported. This is the identity the whole observable rests on, and it is the reason
    T = n*A^2 is free of the contact-count confound."""
    rng = np.random.default_rng(1234)
    for n in (300, 1200):
        v = rng.normal(size=(800, n, 3))
        v /= np.linalg.norm(v, axis=2, keepdims=True)
        F = np.einsum("rni,rnj->rij", v, v) / n
        dev = F - np.eye(3) / 3.0
        a2 = 1.5 * np.einsum("rij,rij->r", dev, dev)
        assert abs(float((n * a2).mean()) - 1.0) < 0.08, n


def test_the_analytic_moments_behind_the_closed_form():
    """The derivation's own intermediate steps, checked numerically: <n_i^4> = 1/5 and
    <n_i^2 n_j^2> = 1/15. If either is wrong the 1/n is wrong."""
    rng = np.random.default_rng(7)
    v = rng.normal(size=(400_000, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    assert abs(float((v[:, 0] ** 4).mean()) - 1 / 5) < 0.003
    assert abs(float((v[:, 0] ** 2 * v[:, 1] ** 2).mean()) - 1 / 15) < 0.002


def test_the_sizing_artifact_reports_the_closed_form_and_the_predicted_scatter(dn):
    assert "PASS: the Monte Carlo reproduces E[T] = 1" in dn
    m = re.search(r"sd\(T\) under the null.*?:\s*([0-9.]+)", dn)
    assert m, "the sizing artifact does not print sd(T)"
    sd = float(m.group(1))
    # chi2_5/5 predicts sqrt(2/5) = 0.6325; a value far from it means T is not the
    # statistic the design assumed.
    assert abs(sd - math.sqrt(2 / 5)) < 0.05, sd


def test_the_confound_is_real_and_in_the_dangerous_direction(dn):
    """The floor must RISE with CV -- i.e. the artefact points the same way a
    restructuring signal would. If it fell, the raw channel would be conservative and
    the whole T construction would be unnecessary."""
    slopes = [float(x) for x in re.findall(r"->\s*[0-9.]+\s+([+-][0-9.]+)\s*$",
                                           dn, flags=re.M)]
    assert len(slopes) >= 3, slopes
    assert all(s > 0 for s in slopes), slopes


# ------------------------------------------------------------------ prereg <-> grader
def test_prereg_registers_twelve_seeds_and_the_effect_size(prereg):
    assert "12 seeds is not a round number, it is the design" in prereg
    assert "DELTA_T = 0.5" in prereg or "ΔT = 0.5" in prereg


def test_prereg_registers_the_verdict_grammar_including_uninformative(prereg):
    """A null that cannot come back UNINFORMATIVE is not a graded null."""
    assert "UNINFORMATIVE" in prereg
    # case-insensitive: the prereg writes "**UNINFORMATIVE**. Not a pass." The prereg is
    # REGISTERED, so the test bends to the document and not the other way round.
    assert "not a pass" in prereg.lower()


def test_prereg_says_this_is_not_a_phi_J_result(prereg):
    assert "not a phi_J result" in prereg or "not a φ_J result" in prereg


def test_the_grader_ran_F0_before_anything_else(gf):
    i0 = gf.index("F0  REPRODUCTION")
    assert i0 < gf.index("PRODUCTION SET")
    assert i0 < gf.index("F2  PRIMARY")


def test_F0_and_F1_pass_in_the_committed_artifact(gf):
    f0 = re.search(r"F0  REPRODUCTION.*?VERDICT: (\w+)", gf, flags=re.S)
    f1 = re.search(r"F1  BIT IDENTITY.*?VERDICT: (\w+)", gf, flags=re.S)
    assert f0 and f0.group(1) == "PASS", f0.group(1) if f0 else "F0 missing"
    assert f1 and f1.group(1) == "PASS", f1.group(1) if f1 else "F1 missing"


def test_the_bit_identity_gate_is_read_on_phi_and_says_so(gf):
    assert "BIT-IDENTICAL True" in gf
    assert gf.count("BIT-IDENTICAL True") >= 3
    assert "quantized" in gf and "printed" in gf


def test_the_tag_guard_admits_only_fabric_shards(gf):
    files = re.findall(r"^\s{4}(phase2_\S+\.jsonl)$", gf, flags=re.M)
    assert files, "no input files listed"
    for f in files:
        assert f.startswith("phase2_fab"), f
        for bad in ("mixr", "extr", "tern", "phase2_kl", "phase2_ug"):
            assert bad not in f, (f, bad)


def test_k4_is_recorded_as_degenerate_not_as_a_result(gf, prereg):
    """At mu_p = 0.3 the k = 4 backbone is EMPTY, and a verdict quoted from it would be
    a verdict on nothing."""
    assert "DEGENERATE" in gf
    assert "k = 4 carries no F2" in prereg or "k = 4 carries no F2 verdict" in prereg


def test_every_filter_prints_its_N(gf):
    ns = re.findall(r"N (?:all records|after [^:]+): (\d+)", gf)
    assert len(ns) >= 4, ns
    counts = [int(x) for x in ns]
    assert counts == sorted(counts, reverse=True), counts


# ------------------------------------------------- the claim itself, pinned to the grade
CLAIMS = os.path.join(HERE, "CLAIMS_FABRIC_2026-09-18.md")


@pytest.fixture(scope="module")
def claims():
    return _txt(CLAIMS)


def _f2_rows(gf):
    """(family, slope, se) from the PRIMARY k = 2 block of the grade."""
    blk = re.search(r"T_nr at k = 2  \(PRIMARY\)(.*?)\n\n", gf, flags=re.S)
    assert blk, "no primary F2 block in the grade"
    out = {}
    for ln in blk.group(1).splitlines():
        m = re.match(r"\s+(\w+)\s+([+-][0-9.]+)\s+([0-9.]+)\s+([0-9.]+)", ln)
        if m:
            out[m.group(1)] = (float(m.group(2)), float(m.group(3)))
    return out


def test_F2_returns_the_null_on_all_three_families(gf):
    """The verdict, allowing for F6's label, which is appended in brackets."""
    for fam in ("lognormal", "uniform", "bidisperse"):
        m = re.search(rf"F2 VERDICT {fam}\s*:\s*(.+)", gf)
        assert m, fam
        assert m.group(1).strip().startswith("NULL, INFORMATIVE"), (fam, m.group(1))


def test_F6_label_is_WIRED_and_not_merely_printed(gf):
    """The registered safeguard: a family with a wholly flagged CV rung must have its
    verdict LABELLED. The first grader printed the flag counts and nothing downstream
    read them, so the label never fired (referee, 2026-09-18). This test fails if the
    wiring is removed OR if the label fires on a family with no wholly flagged rung."""
    wholly = dict(re.findall(r"^\s{4}(\w+)\s+\d+/\s*\d+ flagged.*$\n"
                             r"\s+-> WHOLLY FLAGGED RUNGS: (\[[^\]]*\])",
                             gf, flags=re.M))
    labelled = set(re.findall(r"F2 VERDICT (\w+)\s*:.*\[RESTS ON FLAGGED CELLS\]", gf))
    assert set(wholly) == labelled, (wholly, labelled)
    assert labelled, "no family carries the label; this run is known to have one"


def test_the_flagged_exclusion_does_not_flip_any_verdict(gf):
    """F6's own bar: the verdict must survive dropping the flagged rungs AND dropping
    every flagged cell. The artifact prints a warning line if either would flip it."""
    assert "would FLIP the verdict to DRIFT" not in gf


def test_the_companion_ANOVA_ran_and_agrees(gf):
    """An OLS slope has no power against a non-monotone CV dependence. The companion is
    reported BESIDE F2, never instead of it, and must agree with it here."""
    blk = re.search(r"COMPANION TEST.*?F2 VERDICT", gf, flags=re.S)
    assert blk, "the companion ANOVA is missing"
    rows = re.findall(r"^\s{6}(\w+)\s+([0-9.]+)\s+([0-9.]+)\s+(YES|no)\s*$",
                      blk.group(0), flags=re.M)
    assert len(rows) == 3, rows
    for fam, f_, p_, verdict in rows:
        assert verdict == "no", (fam, p_)


def test_the_dropped_cells_are_reported_against_the_axis_under_test(gf):
    """Cells dropped by the converged filter must be broken down, not counted. Here all
    of them sit on the top CV rung, which is the axis being measured."""
    m = re.search(r"N after requiring converged: (\d+)\s+\(dropped (\d+)\)", gf)
    assert m, "no converged filter line"
    if int(m.group(2)) > 0:
        assert "ITERATION CAP" in gf
        assert "correlate with the axis under test" in gf


def test_every_slope_is_inside_its_own_two_se(gf):
    """The verdict grammar, re-applied here from the printed numbers rather than
    trusted from the printed verdict."""
    rows = _f2_rows(gf)
    assert set(rows) == {"lognormal", "uniform", "bidisperse"}, rows
    for fam, (s, se) in rows.items():
        assert abs(s) <= 2 * se, (fam, s, se)


def test_the_null_is_INFORMATIVE_and_beats_the_registered_resolving_power(gf):
    """A null is only a result if the design could have found the registered effect.
    Registered: DELTA_T = 0.5 across CV 0 -> 0.5, i.e. 2*se*0.5 <= 0.5."""
    for fam, (s, se) in _f2_rows(gf).items():
        assert 2 * se * 0.5 <= 0.5, (fam, se)


def test_the_raw_aniso_channel_IS_the_floor(gf):
    """F4: the paired (aniso - floor) slope must be inside 2 se on every family, and the
    raw slope must be NON-zero -- a routine returning zero for both would pass the first
    half of that and mean nothing."""
    blk = re.search(r"F4  THE CONFOUND.*?F4 VERDICT", gf, flags=re.S).group(0)
    rows = re.findall(
        r"\s+(\w+)\s+([+-][0-9.]+)\s+([+-][0-9.]+)\s+([+-][0-9.]+)\s+([0-9.]+)\s+(YES|NO)",
        blk)
    assert len(rows) == 3, rows
    for fam, raw, floor, diff, se, verdict in rows:
        assert verdict == "YES", (fam, verdict)
        assert abs(float(diff)) <= 2 * float(se), (fam, diff, se)
        # The raw slope must be REAL and POSITIVE, or "the rise is the floor" is a
        # statement about nothing: a routine returning zero for every slope would
        # otherwise satisfy the line above. `se` here is the se of the DIFFERENCE, so
        # it is not the right yardstick for the raw slope's own significance; the
        # defensible statement is that the floor accounts for more of the raw slope
        # than is left over.
        assert float(raw) > 0 and float(floor) > 0, (fam, raw, floor)
        assert abs(float(diff)) < abs(float(raw)), (fam, diff, raw)


def test_the_claims_file_agrees_with_the_grade_it_was_generated_from(gf, claims):
    """This project has shipped artifacts whose prose contradicted their own table. The
    claims page is GENERATED, so this can only fail if the generator drifts -- which is
    exactly the failure worth catching."""
    for fam, (s, se) in _f2_rows(gf).items():
        assert f"| {fam} | {s:+.4f} | {se:.4f} |" in claims, fam


def test_the_claims_file_refuses_the_readings_the_campaign_cannot_support(claims):
    low = claims.lower()
    assert "not a `phi_j` result" in low or "not a phi_j result" in low
    assert "sheared" in low
    assert "absolute value of `t` carries no claim" in low
