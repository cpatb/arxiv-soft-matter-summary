"""test_additivity.py -- regression tests for the 2026-09-17 additivity block:

  * grade_additivity.py / out/grade_additivity.txt  -- the one-parameter mixing repair
  * diag_k6_interp.py  / out/diag_k6_interp.txt     -- is "K = 6 attemptable" a
                                                       measurement of the signal or of
                                                       the interpolation form?

STRUCTURAL SCOPE GUARD. This file must not be able to re-decide any CLOSED gate, and it
must not be able to pass by calling the same code it is checking. It therefore imports
**no project module at all**: every arithmetic claim below is re-implemented here from
scratch and compared against the committed artifact. The guard asserts that on the import
lines, as `test_u3decide.py` does.

WHY ARTIFACT + REIMPLEMENTATION RATHER THAN CALLING THE GRADER. `grade_additivity.py`
takes minutes (it rebuilds the tabulated closure), and a slow test that regenerates the
artifact it validates destroys its own reference on failure. The artifacts are committed,
so their presence is asserted rather than skipped over.
"""

from __future__ import annotations

import os
import re

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
GA = os.path.join(OUT, "grade_additivity.txt")
K6 = os.path.join(OUT, "diag_k6_interp.txt")
PREREG = os.path.join(HERE, "PREREG_ADDITIVITY_2026-09-17.md")


# --------------------------------------------------------------------- scope guard
def test_this_file_imports_no_project_module():
    """A test that imports the grader can pass by agreeing with the bug it contains."""
    src = open(os.path.join(HERE, "test_additivity.py"), encoding="utf-8").read()
    imports = [ln.strip() for ln in src.splitlines()
               if re.match(r"^\s*(import|from)\s", ln)]
    banned = ("grade_", "phase12_", "ugap", "kladder", "mixture_closure",
              "mixture_model", "transfer_", "diag_", "psd")
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
def ga():
    return _txt(GA)


@pytest.fixture(scope="module")
def k6():
    return _txt(K6)


def _num(text, pattern, group=1):
    m = re.search(pattern, text)
    assert m, f"pattern not found in artifact: {pattern}"
    return float(m.group(group))


# =================================================== the mixing fit, re-implemented
def _lam_bruteforce(m, A, S, scale):
    """The minimiser of chi2, found by a dense grid + parabolic refine. Deliberately NOT
    the closed form the grader uses -- two routes to the same number, the second of which
    knows nothing about the first."""
    dm = (np.asarray(m, float) - np.asarray(A, float)) / scale
    ds = (np.asarray(S, float) - np.asarray(A, float)) / scale
    grid = np.linspace(-3, 3, 600001)
    chi = ((dm[:, None] - grid[None, :] * ds[:, None]) ** 2).sum(axis=0)
    return float(grid[int(np.argmin(chi))])


def test_closed_form_lambda_matches_an_independent_minimiser():
    rng = np.random.default_rng(20260917)
    scale = 2.1e-3
    for _ in range(20):
        n = int(rng.integers(5, 25))
        A = 0.60 + 0.02 * rng.standard_normal(n)
        S = A + 0.01 * rng.standard_normal(n)
        m = A + 0.4 * (S - A) + 1e-3 * rng.standard_normal(n)
        dm = (m - A) / scale
        ds = (S - A) / scale
        closed = float(np.sum(dm * ds) / np.sum(ds * ds))
        assert abs(closed - _lam_bruteforce(m, A, S, scale)) < 1e-4


def test_the_closed_form_is_not_trivially_right():
    """A gate must be able to fail: a WRONG lambda formula must be caught by the above."""
    rng = np.random.default_rng(7)
    scale = 2.1e-3
    A = 0.60 + 0.02 * rng.standard_normal(12)
    S = A + 0.01 * rng.standard_normal(12)
    m = A + 0.4 * (S - A) + 1e-3 * rng.standard_normal(12)
    dm, ds = (m - A) / scale, (S - A) / scale
    wrong = float(np.sum(dm * ds) / np.sum(dm * dm))     # denominator swapped
    assert abs(wrong - _lam_bruteforce(m, A, S, scale)) > 1e-3


# ====================================================== the artifact's own integrity
def test_reproduction_gate_passed_and_every_line_is_OK(ga):
    block = ga.split("REPRODUCTION GATE")[1].split("=== WHAT lambda")[0]
    rows = [ln for ln in block.splitlines() if re.match(r"\s+(OK|FAIL)\s", ln)]
    assert len(rows) >= 14, f"only {len(rows)} reproduction rows"
    assert not [r for r in rows if r.strip().startswith("FAIL")], "a reproduction FAILED"
    assert "reproduction PASSES" in ga


def test_the_standing_numbers_are_the_recorded_ones(ga):
    """108.1 and 34.9 are the numbers this whole item is graded against."""
    assert abs(_num(ga, r"held-out chi2, ADDtab\s+got\s+([\d.]+)") - 108.1) < 0.06
    assert abs(_num(ga, r"held-out chi2, ADDj \(J\)\s+got\s+([\d.]+)") - 34.9) < 0.06


def test_sign_convention_is_asserted_in_the_artifact(ga):
    """transfer_kladder uses (model - phi) and this file uses (phi - model). The project
    has already shipped one artifact that contradicted its own table over exactly this."""
    assert "residual = MEASURED - MODEL" in ga
    assert abs(_num(ga, r"transfer_kladder -1\.69 flips\s+got\s+(-?[\d.]+)")
               - (-1.6946)) < 1e-3


def test_K2_cells_carry_no_leverage(ga):
    """At K = 2 the nested rule IS the additive one, so lambda cannot move those cells.
    If this ever stops being true the mixing model is reading a different object."""
    held = re.search(r"held out\s+(\d+) of\s+(\d+) cells carry leverage at d43"
                     r"\s+\(K = 2 cells: (\d+)\)", ga)
    assert held, "leverage line missing"
    lev, tot, k2 = (int(g) for g in held.groups())
    assert tot == 15 and k2 == 5 and lev == tot - k2, (lev, tot, k2)
    assert _num(ga, r"chi2 FLOOR: ([\d.]+) of") > 0


def test_the_registered_failures_are_exactly_the_recorded_ones(ga):
    """The bars do not move. If a later edit makes one of these pass, that is a finding
    to be explained, not a green test to be enjoyed."""
    line = ga.strip().splitlines()[-1]
    assert line.startswith("FAILURES: ["), line
    assert "P0 calibration arms disagree" in line
    assert "P2b" in line and "P3" in line
    for p in ("'P1 ", "'P2 held-out", "'P4 "):
        assert p not in line, f"{p} was recorded as PASSING; it now fails"


def test_P0_fails_on_the_registered_sentence_and_the_sharper_stat_is_printed(ga):
    """A pre-registered rule is a SENTENCE. P0 misses by 0.0017 while the proper
    two-sample statistic says the arms agree at t = 0.69; BOTH must stay printed and the
    verdict must stay FAIL."""
    assert "P0 STAYS FAILED on the sentence that was registered." in ga
    t = _num(ga, r"se of the DIFFERENCE [\d.]+: t = (\d+\.\d+)")
    assert t < 1.0, f"t = {t}: the arms no longer agree even on the loose test"
    assert _num(ga, r"it misses by (\d+\.\d+)") < 0.01


def test_the_oracle_separates_form_failure_from_transfer_failure(ga):
    """The oracle uses HELD-OUT data and is not a result. It is here because it decides
    which of two very different statements the P3 failure supports."""
    assert "It is NOT a result" in ga
    assert "TRANSFER failure, not a form failure" in ga
    gap = _num(ga, r"a gap of ([\d.]+) half-widths")
    assert gap > 5, f"gap {gap}: the ladders and the families no longer disagree"


def test_the_two_ladder_arms_agree_far_better_than_they_agree_with_the_families(ga):
    """The actual result of the item, in one inequality."""
    t_arms = _num(ga, r"ternary vs K-ladder : \|d\| = [\d.]+, t = ([\d.]+)")
    t_fam = _num(ga, r"ladders vs families : \|d\| = [\d.]+, t = ([\d.]+)")
    assert t_arms < 1.0 < t_fam
    assert t_fam / t_arms > 4


def test_the_convention_spread_is_reported_and_did_not_shrink(ga):
    """Mixing was not allowed to hide the composite-diameter convention. It made it
    WORSE, and that has to stay visible."""
    sp = _num(ga, r"spread [\d.]+ to [\d.]+, a factor ([\d.]+)")
    assert sp > 11.9, (f"mixed spread {sp} is below the UNMIXED 11.90x -- if that is "
                       "real it is a finding, not a green test")
    assert "11.90x" in ga


def test_within_group_split_is_printed_beside_the_pooled_rho(ga):
    """A weight table says where phi is SENSITIVE, not where the model FAILS; and the
    pooled rho overstates independence."""
    assert abs(_num(ga, r"pooled  \(lognormal \+ uniform, n = 10\): rho = (-?[\d.]+)")
               - (-0.939)) < 0.002
    assert abs(_num(ga, r"lognormal  n = 5: rho = (-?[\d.]+)") - (-1.0)) < 0.002
    assert abs(_num(ga, r"uniform    n = 5: rho = (-?[\d.]+)") - (-0.9)) < 0.002


def test_prereg_and_grader_agree_on_the_primary_convention():
    """A primary chosen after the fact is not a primary."""
    pr = _txt(PREREG)
    assert "Primary convention: `c = d43`" in pr
    src = open(os.path.join(HERE, "grade_additivity.py"), encoding="utf-8").read()
    assert re.search(r'^PRIMARY = "d43"', src, re.M)


# ======================================================== the K = 6 interpolation
def test_k6_control_reproduces_and_the_document_number_is_recovered(k6):
    assert abs(_num(k6, r"computed here ([\d.]+) -> OK") - 2.2056) < 0.002
    assert abs(_num(k6, r"by its own form: ([\d.]+) vs printed") - 2.99) < 0.006


def test_k6_quadratic_reproduces_independently(k6):
    """Re-derive the referee's counter-number with a hand-rolled Lagrange quadratic --
    no polyfit, no project module -- from the three signals the artifact prints."""
    rows = re.findall(r"^\s+(0\.\d+)\s+[\d.]+\s+([\d.]+)\s*$", k6, re.M)
    assert len(rows) == 3, rows
    us = [float(a) for a, _ in rows]
    ss = [float(b) for _, b in rows]
    u = 0.25 ** 0.2
    tot = 0.0
    for i in range(3):
        term = ss[i]
        for j in range(3):
            if j != i:
                term *= (u - us[j]) / (us[i] - us[j])
        tot += term
    printed = _num(k6, r"quadratic in \(u, s\), 3 pts\s+3\s+no\s+([\d.]+)")
    assert abs(tot - printed) < 5e-3, (tot, printed)
    assert abs(tot - 2.558) < 0.01


def test_k6_is_unresolved_and_the_bar_is_inside_the_span(k6):
    """The point of the whole file: the verdict flips with an undeclared convention."""
    lo = _num(k6, r"K = 6 \(u = [\d.]+\): signal spans ([\d.]+) to")
    hi = _num(k6, r"K = 6 \(u = [\d.]+\): signal spans [\d.]+ to ([\d.]+)")
    assert lo < 2.57 < hi, (lo, hi)
    assert "K = 6: UNRESOLVED" in k6
    assert "NOT established by this estimate" in k6


def test_k7_stays_closed_and_robustly(k6):
    """K >= 7 is the half of the forward content that survives; it must not drift."""
    hi = _num(k6, r"K = 7 \(u = [\d.]+\): signal spans [\d.]+ to ([\d.]+)")
    assert hi < 2.57, hi
    assert "verdict 'attemptable' on 0 of 7 forms" in k6
    assert "K = 7: CLOSED, and robustly" in k6


def test_no_two_parameter_form_fits_the_three_points(k6):
    """If a form cannot hit the point it was not given, the interpolation is genuinely
    under-determined and no choice among the forms is privileged -- which is a stronger
    statement than 'the margin is thin'."""
    worst = _num(k6, r"worst miss on the held-out middle point: ([\d.]+)%")
    assert worst > 10.0, worst


def test_the_anchored_forms_use_a_free_constraint_not_a_fitted_one(k6):
    """s(u = 1) = 0 is exact: identical classes, monodisperse pack, phi = beta. If this
    ever reads as a fitted parameter the spread below is being understated."""
    assert "costs no fitted parameter" in k6
    assert k6.count("YES") >= 3


# ============================== the claims file may not carry a hand-typed number
CLAIMS = os.path.join(HERE, "CLAIMS_ADDITIVITY_2026-09-17.md")

# Every headline number the claims file states, and which artifact must print it.
# Hardcoding is what lets a number go stale; this list is the anti-hardcoding guard.
CLAIMED = {
    "grade_additivity.txt": [
        "0.5214", "0.4969", "0.5069", "0.3033", "0.0274", "0.0228", "0.0175", "0.0313",
        "0.0245", "0.2036", "0.69", "5.67", "0.4894", "0.5245", "56.5", "42.8", "34.9",
        "108.1", "178.9", "15.38", "11.90", "14.2", "11.6", "5.4", "50.3", "2.10",
        "1.45", "0.4817", "0.5611", "0.4639", "0.5299", "0.0017", "0.0258", "108.0717",
        "34.8649", "7.5292", "3.3931", "-3.2955", "-0.939", "-1.000", "-0.900",
        "-0.636", "84.4", "47.8", "18.3", "6.7", "9.8", "18.4", "3.3", "-2.35", "2.64",
        "5.80", "0.926", "0.898", "0.671", "0.507",
    ],
    "diag_k6_interp.txt": [
        "7.0470", "4.4273", "2.2056", "2.9937", "2.558", "3.253", "1.27", "2.189",
        "2.356", "1.08", "9.4", "26.0", "18.6", "2.994", "3.214", "2.858", "3.128",
        "2.758", "1.31951",
    ],
}


def test_every_headline_number_in_the_claims_file_is_printed_by_a_script():
    """A number typed into a claims file is a number that can go stale. Each one below
    must appear verbatim in the artifact that is supposed to produce it."""
    # The .md is written with the typographic minus U+2212 and the artifacts print
    # ASCII '-'. Normalising is required for the comparison to mean anything; without
    # it every negative number silently "passes" by never being looked for.
    claims = _txt(CLAIMS).replace("−", "-")
    missing_from_claims, missing_from_artifact = [], []
    for fname, nums in CLAIMED.items():
        art = _txt(os.path.join(OUT, fname))
        for v in nums:
            if v not in claims:
                missing_from_claims.append(f"{v} (claims file no longer states it)")
            if v not in art:
                missing_from_artifact.append(f"{v} (not printed by {fname})")
    assert not missing_from_artifact, missing_from_artifact
    assert not missing_from_claims, missing_from_claims


def test_the_claims_file_records_the_failures_not_a_success(ga):
    """Generate the narrative FROM the verdict, never alongside it."""
    claims = _txt(CLAIMS)
    assert "P2b FAIL" in claims and "P3 FAIL" in claims
    assert "UNRESOLVED, not established" in claims
    assert "TRANSFER failure" in claims
    # and it must not claim the mix beat J, whatever a later edit does to the numbers
    assert "beats J" not in claims
    assert ga.strip().splitlines()[-1].count("'") == 6, "FAILURES no longer has 3 entries"
