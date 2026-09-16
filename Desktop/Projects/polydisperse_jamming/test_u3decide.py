"""test_u3decide.py -- regression tests for the 2026-09-17 block:

  * derivations/k_axis_bound.md   -- the K-axis feasibility bound
  * diag_u3decide.py             -- the go/no-go gate on a U3-narrowing block

STRUCTURAL SCOPE GUARD. This file must not be able to re-decide any CLOSED gate. It
touches no `mixr*`, `extr*`, `tern*`, `kl*` or `ug*` loader and asserts below that it
imports no grader; the only shared module it uses is `kladder.ladder`, which is a pure
function of K and reads no data at all.
"""

from __future__ import annotations

import itertools
import math
import os
import re

import numpy as np
import pytest

import kladder as kl

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
DERIV = os.path.join(HERE, "derivations")

U_LO_TAB = 0.25                      # ADDtab's lowest identifiable knot (r = 4)
U_HI_TAB = 4.0 ** -0.25              # the U-GAP knot, 0.70711
U_LO_FORM = 1.0 / 3.0                # BENCH row C's window for the parametric closure


# --------------------------------------------------------------------- scope guard
def test_this_file_imports_no_grader_and_no_runner():
    """A test file that imports a grader can re-decide a closed gate by accident."""
    src = open(os.path.join(HERE, "test_u3decide.py"), encoding="utf-8").read()
    imports = [ln.strip() for ln in src.splitlines()
               if re.match(r"^\s*(import|from)\s", ln)]
    banned = ("grade_", "phase12_", "ugap", "mixture_closure", "mixture_model",
              "transfer_")
    for ln in imports:
        for bad in banned:
            assert bad not in ln, f"scope guard: '{ln}' pulls in {bad}"
    # and the only project module imported at all must be kladder
    assert [ln for ln in imports if "kladder" in ln], "kladder import vanished"


# ----------------------------------------------------------------------- the bound
def kappa(u_lo: float, u_hi: float) -> float:
    return math.log(u_lo) / math.log(u_hi)


def k_max(u_lo: float, u_hi: float) -> int:
    return 1 + int(math.floor(kappa(u_lo, u_hi) + 1e-12))


def pair_ratios(d):
    d = sorted(d)
    return [d[i] / d[j] for i, j in itertools.combinations(range(len(d)), 2)]


def test_bound_is_exactly_four_on_the_tabulated_window():
    """kappa = 4.0000 EXACTLY, so K_max = 5. Not 3.9999 and not 4.0001."""
    assert kappa(U_LO_TAB, U_HI_TAB) == pytest.approx(4.0, abs=1e-12)
    assert k_max(U_LO_TAB, U_HI_TAB) == 5


def test_the_honest_form_window_gives_four_not_five():
    """BENCH row C's window admits K = 4 only -- the finding must not drift back."""
    assert k_max(U_LO_FORM, U_HI_TAB) == 4
    assert kappa(U_LO_FORM, U_HI_TAB) == pytest.approx(3.16993, abs=1e-5)


def test_the_pre_ugap_window_already_excluded_K5():
    """Reproduces G1c independently: before the U-GAP knot, K = 5 was out of window."""
    assert k_max(0.25, 2.0 / 3.0) == 4


def test_span_is_FORCED_at_the_top_rung_not_merely_tight():
    """At integer kappa and K = K_max the span interval collapses to a point."""
    k = k_max(U_LO_TAB, U_HI_TAB)
    s_lo, s_hi = U_HI_TAB ** -(k - 1), 1.0 / U_LO_TAB
    assert s_lo == pytest.approx(s_hi, abs=1e-12)
    assert s_hi == pytest.approx(4.0, abs=1e-12)
    # and one rung DOWN the span is genuinely free, so the collapse is not automatic
    s_lo4, s_hi4 = U_HI_TAB ** -(k - 2), 1.0 / U_LO_TAB
    assert s_hi4 - s_lo4 > 0.1


def test_shipped_ladder_crosses_out_of_the_window_at_exactly_K6():
    """Anchor outside the derivation: kladder.ladder, which the derivation never saw."""
    for k in (2, 3, 4, 5):
        us = pair_ratios(kl.ladder(k))
        assert min(us) >= U_LO_TAB - 1e-12
        assert max(us) <= U_HI_TAB + 1e-12
    us6 = pair_ratios(kl.ladder(6))
    assert max(us6) > U_HI_TAB
    # and it crosses at exactly the bound's own lower span end, 0.25**(1/5)
    assert max(us6) == pytest.approx(0.25 ** 0.2, abs=1e-12)


def test_bound_is_convention_independent():
    """u = small/large and u = large/small must give the same kappa."""
    assert (math.log(1 / U_LO_TAB) / math.log(1 / U_HI_TAB)
            == pytest.approx(kappa(U_LO_TAB, U_HI_TAB), abs=1e-12))


def test_the_geometric_ladder_attains_the_bound_so_it_is_not_merely_necessary():
    """Sufficiency: the bound must be ATTAINED, or K_max is an over-estimate."""
    for u_lo, u_hi in ((0.25, U_HI_TAB), (1 / 3, U_HI_TAB), (0.20, 0.75)):
        k = k_max(u_lo, u_hi)
        d = (1.0 / u_lo) ** (np.arange(k) / (k - 1))
        us = pair_ratios(d)
        assert min(us) >= u_lo - 1e-12
        assert max(us) <= u_hi + 1e-12
        # and K+1 must be infeasible, else the bound is not tight
        d2 = (1.0 / u_lo) ** (np.arange(k + 1) / k)
        assert max(pair_ratios(d2)) > u_hi + 1e-12


# ------------------------------------------------- the artifacts, read not recomputed
def _txt(name):
    p = os.path.join(OUT, name) if name.endswith(".txt") and "k_axis" not in name \
        else os.path.join(DERIV, name)
    if not os.path.exists(p):
        pytest.fail(f"{p} is missing -- the block's own artifact is not on disk")
    return open(p, encoding="utf-8").read()


def test_k_axis_check_transcript_has_no_failures():
    t = _txt("k_axis_bound_check.txt")
    assert "FAILURES: []" in t, "the derivation's own checks did not pass"
    assert "96 of 96 cells agree" in t, "the brute-force feasibility check moved"


def test_u3decide_control_reproduces_the_graded_number():
    """The gate is worthless if it cannot reproduce grade_ugap's own 33 of 87."""
    t = _txt("diag_u3decide.txt")
    assert "CONTROL PASS -- exact match" in t
    m = re.search(r"JOINT region: (\d+) of (\d+)", t)
    assert m and (int(m.group(1)), int(m.group(2))) == (33, 87)


def test_u3decide_verdict_is_NOT_WORTH_and_the_gap_is_recorded():
    """Pins the go/no-go AND the two numbers it turns on, so a later edit cannot
    flip the verdict without moving one of them."""
    t = _txt("diag_u3decide.txt")
    assert "NOT WORTH cells" in t, "the go/no-go verdict moved"
    m = re.search(r"realistic contraction f = ([0-9.]+) the survival fraction is "
                  r"([0-9.]+)%", t)
    assert m, "the verdict's own numbers are no longer printed"
    f_real, surv = float(m.group(1)), float(m.group(2))
    assert 0.10 < surv / 100 < 0.90, "the verdict says middling; the number must be"
    assert f_real > 0.60, "the verdict says the design misses the needed f; it must"


def test_u3decide_records_that_the_briefs_premise_was_refuted():
    """Clustering at the |dphi/da| peak DESTROYS identifiability. If a later edit
    makes a clustered design win, this test must fail loudly rather than the
    document quietly disagreeing with its own table."""
    t = _txt("diag_u3decide.txt")
    assert "THE BRIEF'S PREMISE IS REFUTED" in t
    # at least one clustered design must be non-identifiable, and the SPAN design must be
    assert re.search(r"n9 TIGHT.*NO", t), "the tight 9-point design no longer fails"
    assert re.search(r"n9 SPAN.*YES", t), "the spread 9-point design no longer passes"


def test_k6_estimate_is_labelled_as_an_estimate():
    """An inference from a 3-point interpolation must never read as a measurement."""
    t = _txt("k_axis_bound_check.txt")
    assert "ESTIMATE from a 3-point log interpolation, not a" in t
    assert "computed here 2.21" in t, "the arm-B control for the signal units is gone"


def test_signal_control_is_in_the_same_units_as_the_bar():
    """The 2.57 bar is in max(phi-beta)/scale. Comparing it against the U1 'phi span'
    column (arm B: 2.6) instead of max(phi-beta) (2.21) would be a units error that
    reverses the K = 7 verdict. Pinned by the control line."""
    t = _txt("k_axis_bound_check.txt")
    m = re.search(r"0\.80000\s+([0-9.]+)\s+computed here", t)
    assert m and float(m.group(1)) == pytest.approx(2.21, abs=0.01)
