"""test_kladder.py -- regression tests for the K-ladder / y_S-widening block.

Pins the DESIGN (which cannot change without changing the claim), the WIRING
(every reduction the block compares against), and the VERDICTS once they are graded.

Scope guard, as `test_ternary.py` carries: this file must be incapable of re-deciding a
closed block. It reads `phase12_kl*` only, and asserts so structurally.
"""

from __future__ import annotations

import csv
import fnmatch
import math
import os

import numpy as np
import pytest

import kladder as kl
import phase12_ternary as pt
from psd import parse_discrete_family

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")


# ------------------------------------------------------------------ scope / structure

def test_grader_globs_cannot_touch_a_closed_block():
    import grade_kladder as gk
    for g in gk.GLOBS.values():
        assert g.startswith("phase12_kl")
        for closed in ("phase12_tern", "phase12_mixr", "phase12_extr"):
            assert not fnmatch.fnmatch(closed + "_shard00of12.jsonl", g)


def test_bars_are_the_ternary_blocks_and_have_not_moved():
    import grade_kladder as gk
    assert gk.BAR_RMS == 2.02
    assert gk.BAR_MAX == 5.74
    assert gk.SEP_BAR == 1.6
    assert kl.MIN_PHASE_COUNT == 50
    assert kl.PHI_SCALE == 2.1e-3


# ------------------------------------------------------------------------- the design

def test_span_is_four_and_endpoints_are_shared_by_every_rung():
    for k in kl.K_LEVELS:
        d = kl.ladder(k)
        assert len(d) == k
        assert d[0] == 1.0
        assert math.isclose(d[-1], kl.SPAN, rel_tol=0, abs_tol=1e-12)


def test_every_pair_on_every_rung_is_at_u_at_least_one_quarter():
    """The constraint that refused d = 1, 2, 8 for the ternary block."""
    for k in kl.K_LEVELS:
        for _, _, u in kl.pair_table(kl.ladder(k)):
            assert u >= 0.25 - 1e-12


def test_q_solves_the_fines_fraction_exactly():
    for k in kl.K_LEVELS:
        d = kl.ladder(k)
        for ys in kl.YS_LEVELS:
            y = kl.y_from_q(d, kl.q_for_ys(d, ys))
            assert abs(y[0] - ys) < 1e-14
            assert abs(sum(y) - 1.0) < 1e-12


def test_the_two_axes_are_exactly_uncrossed_on_the_balanced_block():
    """corr(K, ln y_S) = 0 BY CONSTRUCTION. This is the whole design."""
    cells = [c for c in kl.build_cells()[0]
             if c["arm"] == "K" and c["y_S_target"] in kl.YS_BALANCED]
    assert len(cells) == 12
    r = np.corrcoef([c["K"] for c in cells],
                    np.log([c["y_S_target"] for c in cells]))[0, 1]
    assert abs(r) < 1e-12


def test_the_registered_refusals_are_refused_and_for_the_registered_reason():
    _, refused = kl.build_cells()
    got = {r["label"]: min(r["counts_judged"]) for r in refused}
    assert got == {"K4_yS0.13": 47, "K4_yS0.18": 31,
                   "K5_yS0.13": 30, "K5_yS0.18": 19}


def test_every_admitted_cell_clears_min_phase_count_at_N1500():
    for c in kl.build_cells()[0]:
        assert min(c["counts_judged"]) >= kl.MIN_PHASE_COUNT, c["label"]


def test_the_two_N_blocks_contain_the_same_core_cells():
    """Otherwise the N-doubling is partly a change of sample."""
    a = {c["label"] for c in kl.build_cells(kl.N_PROD)[0]}
    b = {c["label"] for c in kl.build_cells(kl.N_BIG)[0]}
    assert a <= b
    assert (b - a) == {"KX4_yS0.13", "KX4_yS0.18", "KX5_yS0.13"}
    assert all(c.startswith("KX") for c in b - a)


def test_family_strings_round_trip_to_the_same_floats():
    for c in kl.build_cells(kl.N_BIG)[0]:
        got = parse_discrete_family(c["family"])
        assert got is not None
        d, y = got
        assert tuple(d) == tuple(c["d_rel"])
        assert tuple(y) == tuple(c["y"])


def test_mono_and_NA2_are_byte_identical_to_the_ternary_blocks_strings():
    """So at the same seed they are the SAME packing -- the cross-block gate."""
    ours = {c["label"]: c["family"] for c in kl.build_cells()[0]}
    assert ours["mono"] == pt.y_string((1.0, 0.0, 0.0))
    assert ours["NA2_x6"] == pt.y_string((0.55, 0.35, 0.10))


def test_the_split_control_preserves_the_fines_fraction_and_moves_K():
    d3 = kl.ladder(3)
    y3 = kl.y_from_q(d3, kl.q_for_ys(d3, 0.05))
    d4, y4 = kl._split_middle(d3, y3)
    assert len(d4) == 4
    assert y4[0] == y3[0]                       # the fines class is untouched
    assert abs(sum(y4) - 1.0) < 1e-12
    assert abs((y4[1] + y4[2]) - y3[1]) < 1e-12  # the middle's volume is conserved
    assert abs(y4[3] - y3[2]) < 1e-15


# ------------------------------------------------------------------------ the wiring

def _closures(beta):
    import grade_mixture_ext as gx
    import mixture_closure as mc
    pts23 = mc.binary_points(mu_p=kl.MU_P)
    allpts = gx.ext_points(beta)
    cand, _, _ = gx.candidates(beta, pts23)
    paJ, pbJ = cand["J"](0)
    _, fa, fb = gx.tabulated_closure(beta, allpts)
    return (paJ, pbJ), fa, fb


@pytest.fixture(scope="module")
def ctx():
    import mixture_model as mm
    beta = mm.measured(mu_p=kl.MU_P)[("lognormal", 0.0)][0]
    (paJ, pbJ), fa, fb = _closures(beta)
    return beta, paJ, pbJ, fa, fb


def test_ideal_ceiling_reduces_to_furnas_exactly(ctx):
    import mixture_model as mm
    beta = ctx[0]
    for r in (2.0, 3.0, 4.0):
        for x in (0.05, 0.15, 0.30, 0.50):
            assert abs(kl.ideal_ceiling((1.0, r), (x, 1.0 - x), beta)
                       - mm.furnas_ideal(x, beta)) < 1e-12


def test_ideal_ceiling_is_not_the_numerical_limit(ctx):
    """The discriminating negative. The first version drove `packing_density`'s
    exponents to 1e-13, which leaves a ~3e-7 and misses furnas_ideal by 3e-8. This
    test fails for that implementation and passes for the substituted-zero one."""
    import mixture_model as mm
    beta = ctx[0]
    limit_version = mm.packing_density(np.array([1.0, 4.0]), np.array([0.15, 0.85]),
                                       beta, 1e-13, 1e-13)
    exact = kl.ideal_ceiling((1.0, 4.0), (0.15, 0.85), beta)
    assert abs(limit_version - exact) > 1e-9
    assert abs(exact - mm.furnas_ideal(0.15, beta)) < 1e-15


def test_sequential_k_at_K2_is_the_additive_binary_closure(ctx):
    import mixture_closure as mc
    import mixture_model as mm
    beta, paJ, pbJ, fa, fb = ctx
    for r in (2.0, 3.0, 4.0, 5.0):
        for x in (0.05, 0.15, 0.30, 0.50):
            d, y = (1.0, r), (x, 1.0 - x)
            assert abs(kl.sequential_k(d, y, beta, fa, fb, "d32")
                       - mc.packing_density_tabulated(np.array(d), np.array(y),
                                                      beta, fa, fb)) < 1e-12
            aJ = lambda u: float(np.sqrt(1 - (1 - u) ** paJ))   # noqa: E731
            bJ = lambda u: float(1 - (1 - u) ** pbJ)            # noqa: E731
            assert abs(kl.sequential_k(d, y, beta, aJ, bJ, "d32")
                       - mm.packing_density(np.array(d), np.array(y),
                                            beta, paJ, pbJ)) < 1e-12


def test_sequential_k_at_K3_reproduces_grade_ternarys_rule(ctx):
    """The rule the closed block's T4 was graded on. Both conventions, both closures."""
    import grade_ternary as gt
    beta, paJ, pbJ, fa, fb = ctx
    aJ = lambda u: float(np.sqrt(1 - (1 - u) ** paJ))   # noqa: E731
    bJ = lambda u: float(1 - (1 - u) ** pbJ)            # noqa: E731
    for comp in pt.COMPOSITIONS:
        y_l, y_m, y_s = comp
        for dmode in ("d32", "d43"):
            for closure, am, bm in (("TAB", fa, fb), ("J", aJ, bJ)):
                assert abs(kl.sequential_k((1.0, 2.0, 4.0), (y_s, y_m, y_l),
                                           beta, am, bm, dmode)
                           - gt.sequential(comp, dmode, closure, beta,
                                           paJ, pbJ, fa, fb)) < 1e-12


def test_sequential_k_is_NOT_the_additive_rule_at_K3(ctx):
    """A comparison whose two candidates agree has decided nothing. This asserts the
    DISCRIMINATING half: on the registered grid the two rules genuinely differ."""
    import mixture_closure as mc
    beta, _, _, fa, fb = ctx
    seps = []
    for c in kl.build_cells()[0]:
        if c["K"] < 3:
            continue
        d, y = np.asarray(c["d_rel"]), np.asarray(c["y"])
        seps.append(abs(kl.sequential_k(d, y, beta, fa, fb, "d32")
                        - mc.packing_density_tabulated(d, y, beta, fa, fb)) / kl.PHI_SCALE)
    assert float(np.mean(seps)) > 1.6 * 2.02, float(np.mean(seps))


def test_sequential_k_recursion_is_genuine_at_K4(ctx):
    """At K = 4 the coarse composite is itself a 3-class pack. A one-level rule that
    evaluated it ADDITIVELY would give a different number; assert it does."""
    import mixture_closure as mc
    beta, _, _, fa, fb = ctx
    d = kl.ladder(4)
    y = kl.y_from_q(d, kl.q_for_ys(d, 0.08))
    rec = kl.sequential_k(d, y, beta, fa, fb, "d32")

    # the one-level version: composite density from the ADDITIVE closure
    dv, yv = np.asarray(d), np.asarray(y) / sum(y)
    w = float(yv[1:].sum())
    ysub, dsub = yv[1:] / w, dv[1:]
    g_sub = float(mc.packing_density_tabulated(dsub, ysub, beta, fa, fb))
    nsub = ysub / dsub**3
    nsub = nsub / nsub.sum()
    d_eff = float((nsub * dsub**3).sum() / (nsub * dsub**2).sum())
    u = float(dv[0]) / d_eff
    a_, b_ = float(fa(u)), float(fb(u))
    x = float(yv[0])
    one_level = min(g_sub / (1.0 - (1.0 - a_) * x),
                    beta / (1.0 - (1.0 - beta + b_ * beta * (1.0 - 1.0 / beta)) * (1.0 - x)))
    assert abs(rec - one_level) > 1e-9


# ---------------------------------------------------------------------- the artifacts

def test_the_prediction_file_exists_and_was_written_blind():
    p = os.path.join(OUT, "predict_kladder.txt")
    assert os.path.exists(p), "predict_kladder.txt must exist before any grading"
    txt = open(p, encoding="utf-8").read()
    assert "=== BLIND PREDICTIONS ===" in txt
    assert "WIRING: PASS" in txt
    assert "NOT SEPARABLE by this design" in txt   # the H_K/H_P limit, registered


# ------------------------------------------------------------------- the verdicts
# Pinned from out/grade_kladder.txt and out/kladder_cells.csv, written by
# grade_kladder.py on the 135-run campaign (63 at N = 1500, 72 at N = 3000, 135
# admitted, 0 errors, 0 non-converged). These read the CSV rather than re-fitting the
# closure, so a test can never disagree with the artifact by fitting it differently.


def _cells():
    p = os.path.join(OUT, "kladder_cells.csv")
    if not os.path.exists(p):
        pytest.skip(f"{p} missing -- run grade_kladder.py")
    out = {}
    for r in csv.DictReader(open(p, encoding="utf-8")):
        out[(int(r["N"]), r["label"])] = r
    return out


def _R(c, n, lbl):
    return float(c[(n, lbl)]["R"])


def test_C1_the_second_K_increment_is_never_positive_beyond_noise():
    """THE HEADLINE. R4-R3 is large and positive; R5-R4 is zero or negative on 6 of 6.
    The registered growth laws predict increment ratios 1.00 (per class), 1.33 (per
    pair) and 2.00 (per triple); measured, they run -0.22 to +0.14."""
    c = _cells()
    ratios = []
    for n in (1500, 3000):
        for ys in kl.YS_BALANCED:
            r3 = _R(c, n, f"K3_yS{ys:.2f}")
            r4 = _R(c, n, f"K4_yS{ys:.2f}")
            r5 = _R(c, n, f"K5_yS{ys:.2f}")
            assert r4 - r3 > 2.02, (n, ys, r4 - r3)       # the step is real
            assert r5 - r4 < 2.02, (n, ys, r5 - r4)       # and it does not repeat
            ratios.append((r5 - r4) / (r4 - r3))
    assert len(ratios) == 6
    assert max(ratios) < 0.5, ratios                      # far below H_K's 1.00
    assert sum(1 for v in ratios if v <= 0) >= 4, ratios


def test_C1_sub_clause_R5_minus_R3_clears_the_bar_on_all_three_levels():
    c = _cells()
    got = [_R(c, 1500, f"K5_yS{ys:.2f}") - _R(c, 1500, f"K3_yS{ys:.2f}")
           for ys in kl.YS_BALANCED]
    assert all(v >= 2.02 for v in got), got


def test_C1a_the_registered_ranking_is_degenerate_because_R3_is_inside_the_bar():
    """Every registered hypothesis is a MULTIPLE of R3. On 5 of 6 rows |R3| < 2.02, so
    the ranking predicts nothing and the grader must say so."""
    c = _cells()
    small = [(n, ys) for n in (1500, 3000) for ys in kl.YS_BALANCED
             if abs(_R(c, n, f"K3_yS{ys:.2f}")) < 2.02]
    assert len(small) == 5, small
    txt = open(os.path.join(OUT, "grade_kladder.txt"), encoding="utf-8").read()
    assert "ranking DEGENERATE" in txt
    assert "ABSOLUTE INCREMENTS" in txt


def test_C2_a_degenerate_fourth_class_buys_far_less_than_a_genuine_one():
    """The sharpest single contrast in the block: same K, same pair count, and the
    residual moves by a fraction of what a genuinely intermediate class moves it."""
    c = _cells()
    for n in (1500, 3000):
        for ys in (0.05, 0.08):
            split = _R(c, n, f"SPLIT_yS{ys:.2f}") - _R(c, n, f"K3_yS{ys:.2f}")
            genuine = _R(c, n, f"K4_yS{ys:.2f}") - _R(c, n, f"K3_yS{ys:.2f}")
            assert genuine > 2.02
            assert abs(split) / genuine < 0.5, (n, ys, split, genuine)


def test_C2_the_split_barely_moves_the_MODEL_so_the_control_is_clean(ctx):
    """If the additive prediction itself jumped under the split, the control would be
    measuring the model rather than the packing."""
    import mixture_closure as mc
    beta, _, _, fa, fb = ctx
    cells = {c["label"]: c for c in kl.build_cells()[0]}
    for ys in (0.05, 0.08):
        a = cells[f"SPLIT_yS{ys:.2f}"]
        b = cells[f"K3_yS{ys:.2f}"]
        da = float(mc.packing_density_tabulated(np.asarray(a["d_rel"]),
                                                np.asarray(a["y"]), beta, fa, fb))
        db = float(mc.packing_density_tabulated(np.asarray(b["d_rel"]),
                                                np.asarray(b["y"]), beta, fa, fb))
        assert abs(da - db) / kl.PHI_SCALE < 2.02, (ys, (da - db) / kl.PHI_SCALE)


def test_C3_the_yS_trend_is_present_at_K3_and_absent_at_K2():
    """The K = 2 row is the control: the closure is IN-SAMPLE there, and the trend is
    absent and of the opposite sign."""
    from scipy import stats as _st
    c = _cells()
    for n in (1500, 3000):
        ys = list(kl.YS_LEVELS)
        r3 = [_R(c, n, f"K3_yS{v:.2f}") for v in ys]
        r2 = [_R(c, n, f"K2_yS{v:.2f}") for v in ys]
        rho3 = float(_st.spearmanr(ys, r3).statistic)
        rho2 = float(_st.spearmanr(ys, r2).statistic)
        assert rho3 >= 0.70, (n, rho3)
        assert r3[-1] - r3[0] >= 2.02, (n, r3)
        assert rho2 < 0, (n, rho2)          # the control runs the OTHER way


def test_C4_the_additive_closure_has_the_wrong_sign_in_K():
    """Measured phi RISES with K at matched y_S; the closure says it FALLS. A sign test
    that needs no bar. One of the three levels (y_S 0.08, +0.80 vs se 0.85) is not
    resolved on its own and the claim says so; the other two are."""
    c = _cells()
    resolved = 0
    for ys in kl.YS_BALANCED:
        a, b = c[(1500, f"K3_yS{ys:.2f}")], c[(1500, f"K5_yS{ys:.2f}")]
        d = (float(b["phi"]) - float(a["phi"])) / kl.PHI_SCALE
        se = math.hypot(float(a["se"]), float(b["se"])) / kl.PHI_SCALE
        pred = (float(b["ADDtab"]) - float(a["ADDtab"])) / kl.PHI_SCALE
        assert d > 0, (ys, d)               # measured rises, 3 of 3
        assert pred < 0, (ys, pred)         # the closure falls, 3 of 3
        resolved += d > se
    assert resolved == 2, resolved          # honest: 2 of 3 clear their own se


def test_C5_the_residual_GROWS_with_N_which_an_artefact_cannot_do():
    c = _cells()
    both = [lbl for (n, lbl) in c if n == 1500 and (3000, lbl) in c
            and c[(1500, lbl)]["arm"] in ("K", "SPLIT")]
    shifts = [_R(c, 3000, lbl) - _R(c, 1500, lbl) for lbl in both]
    assert len(shifts) == 18
    assert all(v > 0 for v in shifts), shifts
    assert 1.5 < float(np.mean(shifts)) < 2.3, float(np.mean(shifts))


def test_C6_the_NA2_shift_is_real_but_sits_at_the_block_mean():
    """Both halves, because either alone misleads: the shift is 7.6 se (so the
    registered 'seed noise' prediction FAILED) and it is the block-wide norm."""
    c = _cells()
    both = [lbl for (n, lbl) in c if n == 1500 and (3000, lbl) in c
            and c[(1500, lbl)]["arm"] in ("K", "SPLIT")]
    mean_shift = float(np.mean([_R(c, 3000, lbl) - _R(c, 1500, lbl) for lbl in both]))
    na2 = _R(c, 3000, "NA2_x6") - _R(c, 1500, "NA2_x6")
    assert abs(na2 - mean_shift) < 1.0, (na2, mean_shift)
    assert int(c[(1500, "NA2_x6")]["n_seeds"]) == 6
    assert int(c[(3000, "NA2_x6")]["n_seeds"]) == 6


def test_C7_the_sequential_advantage_is_the_size_of_its_own_convention_spread():
    """1.93x over the additive rule, against a 1.96x spread across the four undeclared
    composite-diameter conventions. Pinned so nobody quotes the first without the
    second."""
    txt = open(os.path.join(OUT, "grade_kladder.txt"), encoding="utf-8").read()
    assert "ADDtab/SEQd21 = 1.93" in txt
    assert "a factor 1.96" in txt
    assert "BRACKET: ADDtab below the measurement on 74% of cells" in txt


def test_G6_bracket_is_recorded_as_a_FAILURE_and_the_bar_was_not_moved():
    """The registered sub-clause wanted >= 80 % each way and got 74 %/74 %. The bar
    stays at 80 and the miss is in FAILURES -- alongside G3's N = 3000 miss, which
    only appeared once the gate was widened to the N the registered SENTENCE implies."""
    import grade_kladder as gk
    txt = open(os.path.join(OUT, "grade_kladder.txt"), encoding="utf-8").read()
    assert "G6 bracket 74%/74%" in txt
    assert "G3 N=3000 y_S 0.08: split moved R by +2.22" in txt
    assert "registered: both >= 80%" in txt
    assert gk.BAR_RMS == 2.02          # nothing else moved either


def test_G3_is_graded_at_BOTH_N_because_the_registered_sentence_names_no_N():
    """An earlier version scoped G3's failure condition to N = 1500, which is a
    NARROWER question than the prereg's "|delta| < 2.02 on 2 of 2". Widened after the
    physics-referee pass; this asserts it can now fail at either N."""
    src = open(os.path.join(HERE, "grade_kladder.py"), encoding="utf-8").read()
    g3 = src.split("=== G3 --")[1].split("=== G4 --")[0]
    assert "if abs(d) >= BAR_RMS:" in g3
    assert "if n == 1500 and abs(d) >= BAR_RMS:" not in g3


def test_C2_all_four_split_ratios_are_printed_including_the_negative_one():
    """The claims file quoted "11-32 %" and silently dropped the N = 3000, y_S = 0.05
    value of -4 %. All four are now printed by the grader, so the range cannot be
    quoted from a subset again."""
    txt = open(os.path.join(OUT, "grade_kladder.txt"), encoding="utf-8").read()
    for frac in ("the split buys +11%", "the split buys +18%",
                 "the split buys -4%", "the split buys +32%"):
        assert frac in txt, frac
    claims = open(os.path.join(HERE, "CLAIMS_KLADDER_2026-09-15.md"),
                  encoding="utf-8").read()
    assert "**−4 %**" in claims


def test_G1c_records_that_K5_is_evaluated_where_the_closure_has_no_data():
    """The saturation rung leans on 4 of 10 pairs at u = 0.7071, above every measured
    knot. Computed here from the ladder and the knot set, not read from prose."""
    import grade_mixture_ext as gx
    import mixture_model as mm
    beta = mm.measured(mu_p=kl.MU_P)[("lognormal", 0.0)][0]
    tab, _, _ = gx.tabulated_closure(beta, gx.ext_points(beta))
    top = max(tab)
    assert abs(top - 2.0 / 3.0) < 1e-9
    for k in (2, 3, 4):
        assert not [u for _, _, u in kl.pair_table(kl.ladder(k)) if u > top + 1e-9], k
    above5 = [u for _, _, u in kl.pair_table(kl.ladder(5)) if u > top + 1e-9]
    assert len(above5) == 4, above5
    assert all(abs(u - 0.5**0.5) < 1e-9 for u in above5)
    txt = open(os.path.join(OUT, "grade_kladder.txt"), encoding="utf-8").read()
    assert "SATURATION IS NOT ATTRIBUTABLE" in txt


def test_G1c_records_that_the_argmin_branch_reassigns_across_K(ctx):
    """A step in the ADDITIVE prediction across K is not by itself a filling
    statement: the min-over-branches argmin moves on 3 of 3 y_S rows."""
    beta, _, _, fa, fb = ctx
    cells = {c["label"]: c for c in kl.build_cells(kl.N_BIG)[0]}

    def argmin_branch(d_rel, y):
        d = np.asarray(sorted(d_rel, reverse=True), float)
        yv = np.asarray([y[list(d_rel).index(v)] for v in d], float)
        yv = yv / yv.sum()
        vals = []
        for i in range(d.size):
            tot = sum((1.0 - beta + float(fb(d[i] / d[j])) * beta * (1.0 - 1.0 / beta))
                      * yv[j] for j in range(i))
            tot += sum((1.0 - float(fa(d[j] / d[i]))) * yv[j]
                       for j in range(i + 1, d.size))
            vals.append(beta / (1.0 - tot) if 1.0 - tot > 1e-12 else np.inf)
        return int(np.argmin(vals))

    nswitch = 0
    for ys in kl.YS_BALANCED:
        idx = [argmin_branch(cells[f"K{k}_yS{ys:.2f}"]["d_rel"],
                             cells[f"K{k}_yS{ys:.2f}"]["y"]) for k in kl.K_LEVELS]
        nswitch += len(set(idx)) > 1
    assert nswitch == 3, nswitch
    txt = open(os.path.join(OUT, "grade_kladder.txt"), encoding="utf-8").read()
    assert "BRANCH IDENTITY" in txt


def test_n_min_is_collinear_with_both_axes_so_it_cannot_be_a_null():
    """spearman(y_S, n_min) = -1.0000 at every K, by construction. The 0.6 "bar" the
    prereg quoted beside the pooled correlation was therefore meaningless, and the
    grader now prints the collinearity instead of the bar."""
    from scipy import stats as _st
    for k in kl.K_LEVELS:
        pts = []
        for ys in kl.YS_LEVELS:
            d = kl.ladder(k)
            y = kl.y_from_q(d, kl.q_for_ys(d, ys))
            if min(kl.class_counts(d, y, kl.N_PROD)) >= kl.MIN_PHASE_COUNT:
                pts.append((ys, int(min(kl.class_counts(d, y, kl.N_PROD)))))
        if len(pts) >= 3:
            rho = float(_st.spearmanr([a for a, _ in pts], [b for _, b in pts]).statistic)
            assert abs(rho + 1.0) < 1e-12, (k, rho, pts)
    txt = open(os.path.join(OUT, "grade_kladder.txt"), encoding="utf-8").read()
    assert "CONFOUNDED DIAGNOSTICS, NOT A GATE" in txt
    assert "0.6 BAR THE" in txt          # and it is explicitly withdrawn


def test_C8_the_transfer_reproduces_the_standing_numbers_and_does_not_settle():
    p = os.path.join(OUT, "transfer_kladder.txt")
    if not os.path.exists(p):
        pytest.skip("transfer_kladder.txt missing -- run transfer_kladder.py")
    txt = open(p, encoding="utf-8").read()
    assert "ADDj      34.9" in txt        # the shipped closure's standing chi2
    assert "ADDtab     108.1" in txt      # the tabulated closure's standing chi2
    assert "a factor 11.90" in txt        # the convention spread
    assert "does NOT close the gap to J" in txt


def test_C9_the_top_descriptors_are_one_quantity_seen_four_ways():
    """pair_weight = (1 - sum y_i^2)/2 exactly, so it is a measure of how evenly the
    volume is spread -- and CV, triple_weight and fill_weight are its relatives. The
    scan cannot separate them, and the claim file says so rather than naming a winner."""
    import diag_kladder_descriptor as dd
    for c in kl.build_cells()[0]:
        if c["K"] < 2:
            continue
        y = np.asarray(kl.psd_for(c["d_rel"], c["y"]).volume_fractions, float)
        got = dd.descriptors(c["d_rel"], c["y"])["pair_weight"]
        assert abs(got - (1.0 - float((y**2).sum())) / 2.0) < 1e-12, c["label"]
