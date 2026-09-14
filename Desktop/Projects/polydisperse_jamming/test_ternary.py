"""test_ternary.py -- pin PREREG_TERNARY_2026-09-14.md and PREREG_MIXTURE_EXT_2026-09-14.md
against their artifacts.

Fast tests only: everything here reads the shard files and the model modules. No DEM.

TWO LOADER RULES, both learned the hard way and both enforced here rather than trusted:

  * This file loads `phase12_tern_shard*.jsonl` for the T-gates and
    `phase12_extr*_shard*.jsonl` for the X-gates, and NEVER `phase12_mixr*`. M0-M6 and
    F1-F3 are closed on the mixr block and nothing here may re-decide them.
  * Every loader carries its own friction filter. A test helper without one is how two
    CLOSED gates went red when the mu_p = 0 block landed (2026-09-13): a loader is a
    loader wherever it lives, and `test_mixture.py`'s docstring records the same lesson.
"""

from __future__ import annotations

import collections
import functools
import glob
import json
import os

import numpy as np
import pytest

import mixture_closure as mc
import mixture_model as mm
from phase12_mixture import MIN_PHASE_COUNT, N_PROD
from phase12_ternary import COMPOSITIONS, MONO, make_ternary, y_string
from psd import BinaryPSD, DiscretePSD, make_psd, parse_discrete_family

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
SCALE = mm.PHI_SCALE
MU = 0.3

TERN_SHARDS = sorted(glob.glob(os.path.join(OUT, "phase12_tern_shard*.jsonl")))
EXT_SHARDS = sorted(glob.glob(os.path.join(OUT, "phase12_extr*_shard*.jsonl")))

# The committed inputs are on disk, so their absence is a FAILURE, not a reason to skip.
# A `pytest.skip` guard on the one test in a file that can actually fail is not a test
# (lessons.md 2026-09-06).
assert TERN_SHARDS, "no phase12_tern shard files: the ternary block is missing"
assert EXT_SHARDS, "no phase12_extr shard files: the extension block is missing"


def _load(shards, keyfn):
    g = collections.defaultdict(list)
    n = 0
    for f in shards:
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                n += 1
                if r.get("error") or r.get("n_requested") != N_PROD:
                    continue
                if abs(float(r.get("mu_p", -1)) - MU) > 1e-9:    # the friction filter
                    continue
                if not r.get("converged"):
                    continue
                g[keyfn(r)].append(r)
    return g, n


def tern_cells():
    return _load(TERN_SHARDS,
                 lambda r: (round(float(r["tern_y_L"]), 4), round(float(r["tern_y_M"]), 4),
                            round(float(r["tern_y_S"]), 4)))


def ext_cells():
    return _load(EXT_SHARDS,
                 lambda r: (float(r["mix_ratio"]), round(float(r["x_small_vol"]), 6)))


@functools.lru_cache(maxsize=1)
def closure():
    """Cached: the profile-interval scan behind `tabulated_closure` costs ~2 min, and
    recomputing it once per test made the suite unusable rather than slow. The cache is
    keyed on nothing because every caller wants the same production closure."""
    meas = mm.measured(mu_p=MU)
    beta = meas[("lognormal", 0.0)][0]
    import grade_mixture_ext as gx
    pts23 = mc.binary_points(mu_p=MU)
    allpts = gx.ext_points(beta)
    cand, per, _ = gx.candidates(beta, pts23)
    tab, fa, fb = gx.tabulated_closure(beta, allpts)
    return beta, cand["J"](0), (tab, fa, fb), per, allpts


# =====================================================================================
# DiscretePSD -- the code gate, which has to hold whether or not any cell ever ran
# =====================================================================================

@pytest.mark.parametrize("r,xs", [(2.0, 0.05), (2.0, 0.30), (3.0, 0.25), (4.0, 0.15),
                                  (5.0, 0.10), (3.0, 0.60), (2.0, 0.95)])
def test_discrete_reduces_to_binary_bit_identically(r, xs):
    """The same distribution through two code paths must give the SAME ARRAY.

    Not "agrees to 1e-12": `array_equal`. The first implementation agreed to 1e-16 and
    that is enough to make a DEM packing a different packing.
    """
    rb = BinaryPSD(d10=1.0, ratio=r, x_small=xs).radii_for_yade(
        1500, np.random.default_rng(20260914))
    rd = DiscretePSD(d10=1.0, d_rel=(1.0, r), y=(xs, 1.0 - xs)).radii_for_yade(
        1500, np.random.default_rng(20260914))
    assert np.array_equal(rb, rd)
    # and it is not vacuous: the array must actually contain two distinct diameters
    assert np.unique(rb).size == 2


def test_discrete_k3_moments_match_closed_form():
    """A three-class case against expectations computed from the DEFINITION, not from
    the class under test. Bit-identity at K = 2 is degenerate in every way K = 3 could
    be wrong."""
    w = np.array([0.15 / 1.0**3, 0.25 / 2.0**3, 0.60 / 4.0**3])
    nf = w / w.sum()
    d = np.array([1.0, 2.0, 4.0]) / float(np.sum(nf * np.array([1.0, 2.0, 4.0])))
    m = [float(np.sum(nf * d**k)) for k in (1, 2, 3, 4)]
    p = DiscretePSD(d10=1.0, d_rel=(1.0, 2.0, 4.0), y=(0.15, 0.25, 0.60))
    assert np.allclose(p.number_fractions, nf, rtol=1e-14, atol=0)
    assert np.allclose(p.diameters, d, rtol=1e-14, atol=0)
    assert abs(p.moment(1) - 1.0) < 1e-14
    assert abs(p.d32 - m[2] / m[1]) < 1e-14
    assert abs(p.d43 - m[3] / m[2]) < 1e-14
    assert abs(p.cv - float(np.sqrt(m[1] - m[0] ** 2) / m[0])) < 1e-14
    # the discriminating half: a K=2 implementation cannot produce three diameters
    assert np.unique(p.sample(1500, np.random.default_rng(1))).size == 3


def test_discrete_family_string_round_trips_and_refuses_malformed():
    fam = y_string((0.60, 0.25, 0.15))
    assert parse_discrete_family(fam) == ((1.0, 2.0, 4.0), (0.15, 0.25, 0.6))
    a = make_psd(fam, 1.0, 0.0).radii_for_yade(300, np.random.default_rng(5))
    b = make_ternary((0.60, 0.25, 0.15)).radii_for_yade(300, np.random.default_rng(5))
    assert np.array_equal(a, b)
    for bad in ("discrete:d1/2/4:y0.5/0.5", "discrete:d1/2/4",
                "discrete:d1/2/4:z0.1/0.2/0.7", "discrete:dA/2/4:y0.1/0.2/0.7"):
        with pytest.raises(ValueError):
            parse_discrete_family(bad)
    with pytest.raises(ValueError):
        make_psd(fam, 1.0, 0.3)          # cv is derived, never chosen
    assert parse_discrete_family("binary:r3.0:x0.25") is None


def test_y_string_puts_the_large_class_on_the_large_diameter():
    """The design tuple is (y_L, y_M, y_S) and the family string is ascending in d.
    Getting that reversal wrong would silently run the mirror-image campaign, and the
    check that catches it is the grain COUNT: the coarse class must be the rarest."""
    p = make_ternary((0.60, 0.25, 0.15))
    c = p.class_counts(N_PROD)
    assert c[0] > c[1] > c[2], f"counts {c} are not descending with diameter"
    assert p.diameters[-1] / p.diameters[0] == pytest.approx(4.0, abs=1e-12)
    # and the coarse class carries the most VOLUME while being the fewest grains
    assert p.volume_fractions[-1] == pytest.approx(0.60, abs=1e-12)


# =====================================================================================
# The T gates -- PREREG_TERNARY_2026-09-14.md
# =====================================================================================

def test_ternary_block_is_complete_and_clean():
    g, nraw = tern_cells()
    assert nraw == 21, f"expected 21 ternary records, found {nraw}"
    assert len(g) == 7, f"expected 7 compositions (6 + mono), found {sorted(g)}"
    for k, v in g.items():
        assert len(v) == 3, f"{k} has {len(v)} admitted seeds, not 3"


def test_every_registered_composition_meets_its_own_admission_rule():
    """MIN_PHASE_COUNT was registered BEFORE the run. The refused composition must stay
    refused, and it is checked by the same arithmetic that refused it."""
    for comp in COMPOSITIONS:
        assert int(make_ternary(comp).class_counts(N_PROD).min()) >= MIN_PHASE_COUNT
    assert int(make_ternary((0.50, 0.30, 0.20)).class_counts(N_PROD).min()) < MIN_PHASE_COUNT


def test_t0_monodisperse_control_lands_on_beta():
    beta, _, _, _, _ = closure()
    g, _ = tern_cells()
    v = np.array([c["phi"] for c in g[tuple(round(x, 4) for x in MONO)]])
    dev = (float(v.mean()) - beta) / SCALE
    assert abs(dev) <= 2.0, f"T0: monodisperse control {dev:+.2f} scale from beta"


def test_t1_the_additive_closure_fails_and_fails_UPWARD():
    """THE RESULT. Pinned as a FAILURE with its direction, not as a number to preserve.

    The bar (rms 2.8, max 5.74) is 1.4x the closure's own in-sample rms on the binaries
    and its in-sample maximum. If a future change makes this test pass, the closure has
    genuinely changed and the ternary verdict must be revisited -- which is exactly what
    a regression test on a failing gate is for.
    """
    beta, _, (tab, fa, fb), _, _ = closure()
    g, _ = tern_cells()
    res = []
    for comp in COMPOSITIONS:
        m = float(np.mean([c["phi"] for c in g[tuple(round(x, 4) for x in comp)]]))
        p = make_ternary(comp)
        pred = mc.packing_density_tabulated(p.diameters, p.volume_fractions, beta, fa, fb)
        res.append((m - pred) / SCALE)       # MEASURED - PREDICTED; + = under-predicted
    res = np.asarray(res)
    assert float(np.sqrt(np.mean(res**2))) > 2.8, "T1 no longer fails its registered rms bar"
    assert float(np.max(np.abs(res))) > 5.74, "T1 no longer fails its registered max bar"
    assert int(np.sum(res > 0)) == 6, f"the 6/6 positive sign structure moved: {res}"
    assert float(res.mean()) >= 2.8, f"T3's mean bar no longer clears: {res.mean():+.2f}"


def test_the_k3_model_path_reduces_to_the_k2_one():
    """The check that could have overturned T1: if the K=3 wiring were wrong, the
    residual would be plumbing. Emptying a class must reproduce the binary EXACTLY,
    and which class is emptied decides which ratio survives."""
    beta, (paJ, pbJ), (tab, fa, fb), _, _ = closure()
    for comp3, r, xs in (((0.85, 0.0, 0.15), 4.0, 0.15),
                         ((0.75, 0.0, 0.25), 4.0, 0.25),
                         ((0.0, 0.60, 0.40), 2.0, 0.40),
                         ((0.70, 0.30, 0.0), 2.0, 0.30)):
        p = make_ternary(comp3)
        assert abs(p.size_ratio - r) < 1e-12, "the TEST's ratio is wrong, not the model's"
        got = mm.packing_density(p.diameters, p.volume_fractions, beta, paJ, pbJ)
        assert abs(got - mm.binary_density(r, xs, beta, paJ, pbJ)) < 1e-12
        gotT = mc.packing_density_tabulated(p.diameters, p.volume_fractions, beta, fa, fb)
        wantT = mc.packing_density_ab(np.array([r, 1.0]), np.array([1.0 - xs, xs]),
                                      beta, fa, fb)
        assert abs(gotT - wantT) < 1e-12


def test_t4_the_sequential_alternative_beats_the_additive_one():
    """The non-additive arm wins under the d32 convention, above the 1.6 informativeness
    bar. Pinned WITH the bar, because a win below the bar decides nothing."""
    import grade_ternary as gt
    beta, (paJ, pbJ), (tab, fa, fb), _, _ = closure()
    g, _ = tern_cells()
    rT, rS = [], []
    for comp in COMPOSITIONS:
        m = float(np.mean([c["phi"] for c in g[tuple(round(x, 4) for x in comp)]]))
        p = make_ternary(comp)
        rT.append((m - mc.packing_density_tabulated(p.diameters, p.volume_fractions,
                                                    beta, fa, fb)) / SCALE)
        rS.append((m - gt.sequential(comp, "d32", "TAB", beta, paJ, pbJ, fa, fb)) / SCALE)
    a = float(np.sqrt(np.mean(np.square(rT))))
    b = float(np.sqrt(np.mean(np.square(rS))))
    assert b < a, "SEQUENTIAL no longer beats ADDITIVE"
    assert a / b >= 1.6, f"the win is below the registered informativeness bar: {a / b:.2f}"


def test_the_measurement_sits_below_the_ideal_filling_ceiling():
    """A ternary packing ABOVE the a = b = 0 limit would refute the filling picture
    outright rather than find it non-additive. It does not, on all six."""
    beta, _, (tab, fa, fb), _, _ = closure()
    g, _ = tern_cells()
    for comp in COMPOSITIONS:
        m = float(np.mean([c["phi"] for c in g[tuple(round(x, 4) for x in comp)]]))
        p = make_ternary(comp)
        ideal = mc.packing_density_ab(p.diameters, p.volume_fractions, beta, 0.0, 0.0)
        add = mc.packing_density_tabulated(p.diameters, p.volume_fractions, beta, fa, fb)
        assert add < m < ideal, f"{comp}: {m} not between additive {add} and ideal {ideal}"


def test_t6_the_middle_class_is_resolved_and_is_not_an_endpoint():
    """The structural channel a binary cannot have. The middle class must be genuinely
    intermediate in BOTH z and rattler fraction -- if it simply tracked one of the
    extremes, the ternary would be a binary in disguise and T1 would be about nothing."""
    g, _ = tern_cells()
    n_ok = 0
    for comp in COMPOSITIONS:
        cs = g[tuple(round(x, 4) for x in comp)]
        z = [float(np.mean([c["by_design_class"][k]["z_mean"] for c in cs]))
             for k in ("c0", "c1", "c2")]
        rat = [float(np.mean([c["by_design_class"][k]["rattler_frac"] for c in cs]))
               for k in ("c0", "c1", "c2")]
        assert z[0] < z[1] < z[2], f"{comp}: z not ordered by size: {z}"
        assert rat[0] > rat[1] > rat[2], f"{comp}: rattler fraction not ordered: {rat}"
        n_ok += 1
    assert n_ok == 6


# =====================================================================================
# The NA gates -- PREREG_TERNARY_N_2026-09-14.md, the finite-size discriminator
# =====================================================================================

def na_blocks():
    import grade_ternary_N as gn
    beta_iso, _, (tab, fa, fb), _, _ = closure()
    g15, _ = gn.load("phase12_tern_shard*.jsonl", 1500)
    g30, _ = gn.load("phase12_tern3k_shard*.jsonl", 3000)
    b15 = float(np.mean([c["phi"] for c in g15[tuple(round(v, 4) for v in MONO)]]))
    b30 = float(np.mean([c["phi"] for c in g30[tuple(round(v, 4) for v in MONO)]]))
    return (gn.block_stats(g15, b15, fa, fb, 1500),
            gn.block_stats(g30, b30, fa, fb, 3000), g15, g30)


def test_the_N3000_block_is_complete():
    _, _, _, g30 = na_blocks()
    assert sum(len(v) for v in g30.values()) == 21
    assert len(g30) == 7


def test_na1_the_residual_does_not_shrink_when_the_coarse_phase_doubles():
    """THE DISCRIMINATOR. A finite-size artifact MUST shrink when n_L doubles. Pinned
    with its registered bars, and with the direction, because the direction is the
    result: the residual grows slightly rather than falling."""
    s15, s30, _, _ = na_blocks()
    r15 = np.array([s15[c]["res"] for c in COMPOSITIONS])
    r30 = np.array([s30[c]["res"] for c in COMPOSITIONS])
    assert float(r30.mean()) >= 4.7, "NA1 no longer clears its non-additive bar"
    assert float(r30.mean()) > 2.8, "NA1 fell to the finite-size bar: T1 is withdrawn"
    assert int(np.sum(r30 > 0)) == 6, f"the 6/6 sign structure moved at N=3000: {r30}"
    assert float(r30.mean()) > float(r15.mean()) - 1.0,         "the residual shrank materially at doubled N, which is the artifact signature"
    # n_L really did double, so the test is not vacuous
    for c in COMPOSITIONS:
        assert s30[c]["n_L"] >= 2 * s15[c]["n_L"] - 2


def test_na3_the_sample_size_confound_is_excluded_only_when_POOLED():
    """The within-block correlation is +-1 in BOTH blocks whatever the answer, because
    n_L and y_S are collinear at fixed N. Only the pooled number separates them, and
    asserting BOTH halves is what makes this test discriminating rather than a
    restatement."""
    from scipy.stats import spearmanr
    s15, s30, g15, g30 = na_blocks()
    for s in (s15, s30):
        rho, _ = spearmanr([s[c]["n_L"] for c in COMPOSITIONS],
                           [s[c]["res"] for c in COMPOSITIONS])
        assert rho == pytest.approx(-1.0, abs=1e-9),             f"the within-block collinearity is gone: {rho}"
    nl = [s[c]["n_L"] for c in COMPOSITIONS for s in (s15, s30)]
    rs = [s[c]["res"] for c in COMPOSITIONS for s in (s15, s30)]
    rho_pool, _ = spearmanr(nl, rs)
    assert abs(rho_pool) <= 0.6, f"n_L survives as an explanation pooled: {rho_pool}"
    # the surviving candidate is composition
    rho_y, _ = spearmanr([c[2] for c in COMPOSITIONS],
                         [s30[c]["res"] for c in COMPOSITIONS])
    assert rho_y > 0.6, f"the residual no longer tracks the fines volume: {rho_y}"


def test_the_finite_size_diagnostic_stayed_under_its_flag_in_both_blocks():
    """Read from the SOLVER's record, not recomputed. No cell may carry a note."""
    _, _, g15, g30 = na_blocks()
    for g, nm in ((g15, "N=1500"), (g30, "N=3000")):
        cs = [c for v in g.values() for c in v]
        assert len(cs) == 21, f"{nm}: {len(cs)} cells"
        assert max(c["d_max_over_L"] for c in cs) < 0.25, nm
        assert not [c for c in cs if c.get("notes")], f"{nm} has flagged cells"


# =====================================================================================
# The X gates -- PREREG_MIXTURE_EXT_2026-09-14.md, pinned as the ternary's calibration
# =====================================================================================

def test_ext_block_is_complete():
    g, nraw = ext_cells()
    assert nraw == 39, f"expected 39 extension records, found {nraw}"
    assert sorted({k[0] for k in g}) == [4.0, 5.0]


def test_x0_monodisperse_endpoints_land_on_beta():
    beta, _, _, _, _ = closure()
    g, _ = ext_cells()
    for key in ((4.0, 0.0), (4.0, 1.0), (5.0, 0.0), (5.0, 1.0)):
        v = np.array([c["phi"] for c in g[key]])
        assert abs((float(v.mean()) - beta) / SCALE) <= 2.0, f"X0 at {key}"


def test_x1_J_under_predicts_every_blind_composition():
    """X1's verdict: J is calibrated at u >= 1/3 and is NOT a predictor below it.
    Pinned as 9/9 POSITIVE, which is the statement, rather than as an rms to preserve."""
    beta, (paJ, pbJ), _, _, _ = closure()
    g, _ = ext_cells()
    res = []
    for (r, x), cs in sorted(g.items()):
        if not (0.0 < x < 1.0):
            continue
        bad = min(c["by_class"].get("n_S", 0) for c in cs) < MIN_PHASE_COUNT
        if bad:
            continue
        m = float(np.mean([c["phi"] for c in cs]))
        res.append((m - mm.binary_density(r, x, beta, paJ, pbJ)) / SCALE)
    res = np.asarray(res)
    assert res.size == 9, f"expected 9 blind compositions, got {res.size}"
    assert int(np.sum(res > 0)) == 9, f"X1's 9/9 under-prediction moved: {res}"
    assert float(np.sqrt(np.mean(res**2))) > 4.0, "X1 no longer fails its registered bar"


def test_x3_ratio_4_beats_ratio_3():
    """X3's bench reading: the composition maximum keeps RISING with size ratio."""
    beta, _, _, _, allpts = closure()
    r3 = [p[2] for p in allpts if abs(p[0] - 3.0) < 1e-9]
    r4 = [p[2] for p in allpts if abs(p[0] - 4.0) < 1e-9]
    assert r3 and r4
    assert (max(r4) - max(r3)) / SCALE > 3.0, "r4's peak no longer beats r3's by the bar"


def test_x5_the_large_skeleton_carries_the_load():
    """rat_L stays below 0.10 across the whole coarse branch while rat_S collapses."""
    g, _ = ext_cells()
    ratl, rats = [], []
    for (r, x), cs in sorted(g.items()):
        if not (0.0 < x < 1.0) or abs(r - 4.0) > 1e-9:
            continue
        ratl.append(float(np.mean([c["by_class"]["rat_L_k2"] for c in cs])))
        rats.append(float(np.mean([c["by_class"]["rat_S_k2"] for c in cs])))
    assert ratl, "the by_class rattler channel is absent from the extension block"
    assert float(np.max(ratl)) <= 0.10, f"X5's rat_L bound moved: max {float(np.max(ratl))}"
    # the DISCRIMINATING half: a routine returning a constant, or a channel that simply
    # never varies, would pass the bound above. rat_S must actually move across the
    # coarse branch, and it must sit far above rat_L.
    assert float(np.max(rats)) - float(np.min(rats)) > 0.3,         f"rat_S does not vary across the branch: {rats}"
    assert float(np.min(rats)) > float(np.max(ratl)), "rat_S and rat_L are not separated"


def test_the_closed_mixr_gates_are_not_reachable_from_this_file():
    """A structural guard, not a physics test: this module must never glob the mixr
    block, because M0-M6 and F1-F3 are closed on it."""
    src = open(os.path.join(HERE, "test_ternary.py"), encoding="utf-8").read()
    assert "phase12_mixr" not in src.split('"""')[2], \
        "test_ternary.py must not load the closed mixr block"
