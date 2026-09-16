"""kladder.py -- the K-LADDER / y_S-WIDENING design, and the models graded on it.

ONE SOURCE OF TRUTH for the grid. `predict_kladder.py` (blind predictions, run BEFORE
any cell), `phase12_kladder.py` (the runner) and `grade_kladder.py` (the verdicts) all
import the grid from here and none of them rebuilds it. lessons.md 2026-09-03: two
pieces of code that CONSTRUCT the same identifier by different rules will eventually
disagree, and one of them was about to `rm -rf` a converged case.

PRE-REGISTRATION: PREREG_KLADDER_2026-09-15.md.

THE QUESTION
------------
`grade_ternary.py` established that the pairwise-filling functional is NOT additive
beyond two classes: six ternaries at d = 1, 2, 4 are under-predicted 6/6, mean +7.53
scale, rms 7.83 against a bar of 2.8, confirmed at N = 3000 (+8.67, 6/6) with `n_L` and
`d_max/L` excluded as confounds. Two things weaken it and this block is both of them:

  (i)  `y_S` took only TWO values (0.10, 0.15) across the six cells, with
       `y_L` in [0.55, 0.70]. The residual tracks `y_S` at spearman +0.878, so the
       headline rests on a two-level contrast in a narrow corner.
  (ii) Nothing says whether the miss GROWS with the number of classes or SATURATES at
       the K = 3 value. If it saturates the correction is a constant and can be
       absorbed; if it grows, every continuous-PSD prediction this project makes needs
       the new rule -- and a continuous PSD is K -> infinity.

THE DESIGN, AND WHY IT IS NOT CONFOUNDED
-----------------------------------------
A GEOMETRIC diameter ladder of FIXED SPAN 4, filled in with more classes:

    K = 2 : d = 1, 4                          1 pair,   u = 1/4
    K = 3 : d = 1, 2, 4                       3 pairs,  u = 1/2 (x2), 1/4
    K = 4 : d = 1, 4^(1/3), 4^(2/3), 4        6 pairs,  u = 0.630 (x3), 0.397 (x2), 1/4
    K = 5 : d = 1, 4^(1/4), ..., 4           10 pairs,  u = 0.707 (x4), 1/2 (x3),
                                                            0.354 (x2), 1/4

EVERY pair on EVERY rung has u >= 1/4 BY CONSTRUCTION, because the span is 4 and no two
classes are further apart than the extremes. That is the constraint that refused
d = 1, 2, 8 for the ternary block, and here it is satisfied identically at every K.

(!) THAT BUYS "UNCONFOUNDED WITH EXTRAPOLATION" FOR **ADDtab AND NOTHING ELSE**, AND THE
CLOSURE HAS TO BE NAMED WHENEVER A K-LADDER NUMBER IS QUOTED. `u >= 1/4` is the
measured-knot floor of the TABULATED closure, whose knots are u = 0.250, 0.333, 0.500,
0.667 (and 0.707 since the U-GAP block). It is NOT the honest window of the PARAMETRIC
closures: BENCH row C records that `J` / `EXT` fail blind at size ratios 4 and 5, so
their window starts at u = 1/3, and on [1/3, 0.70711] the K-axis bound
(`derivations/k_axis_bound.md`) gives K_max = 4 rather than 5. So **a K = 5 statement
graded against `J` or `EXT` is outside its own window**, and the earlier wording of this
paragraph -- which said "additivity stays unconfounded with extrapolation on the whole
ladder" without naming a closure -- would have let a later reader quote a K = 5 number
against the wrong one. The ADDtab-graded gates (G1, G1c, G7, U3) are inside their window
and stand.

The extremes
d_min and d_max are the SAME on every rung, so the K axis adds classes without widening
the distribution -- `d_max/L`, the size ratio and the span are all held fixed.

Composition is a ONE-PARAMETER family, y_i proportional to d_i^q, and the grid is
addressed by the FINES VOLUME FRACTION `y_S` = y(d_min) rather than by q: for each
(K, y_S) the exponent q is SOLVED so that the finest class carries exactly `y_S`.

  ==> `y_S` and K are EXACTLY UNCROSSED on the balanced sub-block: at y_S = 0.03, 0.05
      and 0.08 the K axis runs 2, 3, 4, 5 at IDENTICAL fines volume fraction, so
      corr(K, y_S) = 0 by construction rather than by luck.

That is the direct repair of the defect the ternary block had to break with a second N:
there, `n_L` and `y_S` were collinear by construction because n goes as y/d^3. Here the
two axes are orthogonal by construction, and BOTH N are run anyway.

WHAT SEPARATES THE THREE HYPOTHESES (registered; see the prereg for the arithmetic)
----------------------------------------------------------------------------------
With R(K) the mean residual in PHI_SCALE units, calibrated on the measured K = 3 value
R3 and anchored at R(2) = 0 (K = 2 is TAB's own calibration set):

    H_sat  "already saturated / composition only"   R4 = R5 = R3
    H_K    "one unit per extra CLASS"               R4 = 2 R3,   R5 = 3 R3
    H_P    "one unit per extra PAIR", P = K(K-1)/2  R4 = 2.5 R3, R5 = 4.5 R3
    H_T    "one unit per TRIPLE", T = C(K,3)        R4 = 4 R3,   R5 = 10 R3

H_sat / H_K / H_T separate by 2.0x or more at K = 4 and 3.3x or more at K = 5.
**H_K and H_P do NOT separate** -- 1.25x at K = 4, 1.5x at K = 5, both below the 1.6x
informativeness bar -- and that is registered IN ADVANCE as a limit of the design, not
discovered afterwards. And the y_S axis is what separates all four of them from
H_y ("it is a composition effect and K is irrelevant"), because `y_S` is held FIXED
along the K axis and swept 6x along the composition axis.

Usage:
    python3 predict_kladder.py          # blind, before any cell exists
    python3 phase12_kladder.py --shard 0 --nshards 20 --tag kl
    python3 grade_kladder.py
"""

from __future__ import annotations

import math

import numpy as np

from psd import DiscretePSD

# ---------------------------------------------------------------- protocol constants
P_HAT = 1e-3                 # the rung every family, binary and ternary on disk has
N_PROD = 1500
N_BIG = 3000
MU_P = 0.3
SEEDS = (1, 2, 3)
MIN_PHASE_COUNT = 50         # inherited verbatim from phase12_mixture.py
PHI_SCALE = 2.1e-3           # the project's scale unit, unchanged
SPAN = 4.0                   # d_max / d_min, HELD FIXED on every rung of the K ladder

# The registered composition axis, as FINES VOLUME FRACTION of the finest class.
# 0.18 widens past the ternary block's maximum of 0.15; 0.03 is 5x below its minimum.
YS_LEVELS = (0.03, 0.05, 0.08, 0.13, 0.18)
# The balanced sub-block on which K and y_S are exactly uncrossed.
YS_BALANCED = (0.03, 0.05, 0.08)
K_LEVELS = (2, 3, 4, 5)


def ladder(k: int) -> tuple[float, ...]:
    """The K-class geometric diameter ladder of span SPAN, ascending.

    Exact at the endpoints (1.0 and SPAN) for every K, which is what makes `d_max/d_min`,
    the span and the worst pair ratio identical across the whole K axis.
    """
    if k < 1:
        raise ValueError(f"k must be >= 1, got {k}")
    if k == 1:
        return (1.0,)
    return tuple(SPAN ** (i / (k - 1)) for i in range(k))


def y_from_q(d_rel, q: float) -> tuple[float, ...]:
    """Volume fractions y_i proportional to d_i^q, normalised to sum 1."""
    d = np.asarray(d_rel, dtype=float)
    w = d**q
    y = w / w.sum()
    return tuple(float(v) for v in y)


def q_for_ys(d_rel, ys: float) -> float:
    """The exponent q such that the FINEST class carries exactly `ys` of the volume.

    With d_min normalised to 1 this is y_S(q) = 1 / sum_i d_i^q with every d_i >= 1, so
    the denominator is strictly increasing in q and the root is unique. Solved by plain
    BISECTION rather than `scipy.optimize.brentq`: this module is imported by the runner
    inside WSL under yade's own interpreter, which is not guaranteed to carry scipy, and
    a grid builder that cannot import is a campaign that does not launch. 200 halvings of
    [-8, 12] is machine precision.
    """
    d = np.asarray(d_rel, dtype=float)
    if not 0.0 < ys < 1.0:
        raise ValueError(f"ys must lie strictly in (0, 1), got {ys}")
    if d.size == 1:
        raise ValueError("a monodisperse ladder has y_S = 1 identically")

    def f(q: float) -> float:
        w = d**q
        return float(w[0] / w.sum()) - ys

    lo, hi = -8.0, 12.0
    if f(lo) * f(hi) > 0:
        raise ValueError(f"y_S = {ys} is not reachable on the ladder {tuple(d)}")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def psd_for(d_rel, y) -> DiscretePSD:
    return DiscretePSD(d10=1.0, d_rel=tuple(float(v) for v in d_rel),
                       y=tuple(float(v) for v in y))


def family_string(d_rel, y) -> str:
    """The family string, ASCENDING in diameter, at full round-trip precision.

    `repr` of a float is the shortest string that round-trips, so the runner's PSD and
    the PSD `packing.build_and_compress` re-parses out of this string are the SAME
    floats. lessons.md 2026-09-14: 1e-16 in d_small makes a DEM packing a different
    packing, so a rounded family string is a different campaign.
    """
    ds = "/".join(repr(float(v)) for v in d_rel)
    ys = "/".join(repr(float(v)) for v in y)
    return f"discrete:d{ds}:y{ys}"


def class_counts(d_rel, y, n: int) -> np.ndarray:
    return psd_for(d_rel, y).class_counts(n)


def d_max_over_L(d_rel, y, n: int, phi: float) -> float:
    """d_max / L at the DENSE state -- the diagnostic packing.py flags above 0.25.

    L^3 = (pi/6) sum_i n_i d_i^3 / phi. Evaluating at a LOWER bound on phi gives the
    WORST case, which is what an admissibility check wants.
    """
    p = psd_for(d_rel, y)
    counts = p.class_counts(n)
    vol = float(np.pi / 6.0 * np.sum(counts * p.diameters**3))
    return float(p.diameters[-1]) / (vol / phi) ** (1.0 / 3.0)


# ---------------------------------------------------------------------------- the grid

def _split_middle(d_rel, y, frac: float = 0.5, eps: float = 0.05):
    """Split the MIDDLE class of a K = 3 spec into two nearly-degenerate classes.

    d_mid -> d_mid*(1 - eps), d_mid*(1 + eps), each taking `frac` / (1 - `frac`) of the
    middle class's volume. The DISTRIBUTION barely moves; K goes 3 -> 4 and the pair
    count goes 3 -> 6.

    This is the sharpest available discriminator between "the residual is a functional
    of the size distribution" (it must not move) and "the residual counts classes or
    pairs" (it must jump by the same factor the K axis shows). It costs three cells.
    """
    if len(d_rel) != 3:
        raise ValueError(f"_split_middle needs a 3-class spec, got {len(d_rel)}")
    d0, d1, d2 = (float(v) for v in d_rel)
    y0, y1, y2 = (float(v) for v in y)
    return ((d0, d1 * (1.0 - eps), d1 * (1.0 + eps), d2),
            (y0, y1 * frac, y1 * (1.0 - frac), y2))


# The monodisperse protocol control. Written as the 3-class spec with two empty classes
# -- exactly the string `phase12_ternary.py` used for T0 -- so that at the same seed it
# is the BIT-IDENTICAL packing and the control doubles as a cross-block reproduction
# gate. DiscretePSD drops zero-volume classes, so this really is K = 1.
MONO_D = (1.0, 2.0, 4.0)
MONO_Y = (0.0, 0.0, 1.0)

# The NA2 exception carried over from PREREG_TERNARY_N_2026-09-14.md: (y_L, y_M, y_S) =
# (0.55, 0.35, 0.10) moved +1.81 against an se of 0.48 between N = 1500 and N = 3000
# where the other five agreed. Re-run at SIX seeds. Seeds 1-3 must reproduce the `tern`
# block BIT-IDENTICALLY (same string, same class, same seed); seeds 4-6 are the new
# information. The string is byte-for-byte what `phase12_ternary.y_string` emitted.
NA2_D = (1.0, 2.0, 4.0)
NA2_Y = (0.1, 0.35, 0.55)
NA2_SEEDS = (1, 2, 3, 4, 5, 6)


def build_cells(n: int = N_PROD):
    """Every cell of the campaign, as dicts. ADMISSIBILITY IS APPLIED HERE.

    Returns (cells, refused). A cell is REFUSED, not run and dropped, when any class
    falls below MIN_PHASE_COUNT at `n` -- the same rule and the same reason the ternary
    block refused (0.50, 0.30, 0.20) at 48 grains.

    THE CORE GRID'S ADMISSIBILITY IS ALWAYS JUDGED AT `N_PROD`, whatever `n` is asked
    for, so that the N = 1500 and N = 3000 blocks contain the SAME CELLS and the
    N-doubling comparison is a comparison of the same points. lessons.md 2026-09-08:
    evaluating two probes of a sensitivity on different samples makes the difference
    partly a change of sample. Cells that only N = 3000 can afford are a SEPARATE,
    LABELLED arm (`KX`), run only at N_BIG and never pooled into the N comparison.
    """
    cells: list[dict] = []
    refused: list[dict] = []

    def add(arm, label, d_rel, y, seeds=SEEDS, ys=None, kk=None, qq=None,
            judge_at: int = N_PROD, only_at: int | None = None):
        counts = class_counts(d_rel, y, judge_at)
        rec = {"arm": arm, "label": label, "d_rel": tuple(float(v) for v in d_rel),
               "y": tuple(float(v) for v in y), "seeds": tuple(seeds),
               "K": int(kk if kk is not None else len(counts)),
               "y_S_target": ys, "q": qq, "judged_at": judge_at, "only_at": only_at,
               "counts": [int(v) for v in class_counts(d_rel, y, n)],
               "counts_judged": [int(v) for v in counts],
               "family": family_string(d_rel, y)}
        if only_at is not None and n != only_at:
            return
        if len(counts) > 1 and int(counts.min()) < MIN_PHASE_COUNT:
            rec["why_refused"] = (f"min class count {int(counts.min())} < "
                                  f"{MIN_PHASE_COUNT} at N = {judge_at}")
            refused.append(rec)
        else:
            cells.append(rec)

    add("T0", "mono", MONO_D, MONO_Y, kk=1)

    for k in K_LEVELS:
        d = ladder(k)
        for ys in YS_LEVELS:
            q = q_for_ys(d, ys)
            y = y_from_q(d, q)
            add("K", f"K{k}_yS{ys:.2f}", d, y, ys=ys, kk=k, qq=q)

    # The degenerate-split control, at the two y_S levels of the balanced block that
    # have the most grain headroom.
    for ys in (0.05, 0.08):
        d3 = ladder(3)
        q = q_for_ys(d3, ys)
        y3 = y_from_q(d3, q)
        d4d, y4d = _split_middle(d3, y3)
        add("SPLIT", f"SPLIT_yS{ys:.2f}", d4d, y4d, ys=ys, kk=4, qq=q)

    add("NA2", "NA2_x6", NA2_D, NA2_Y, seeds=NA2_SEEDS, ys=0.10, kk=3)

    # Arm KX: the three (K, y_S) corners the coarse class cannot afford at N = 1500 and
    # can at N = 3000. They EXTEND the y_S axis to K = 4 and 5 -- the K x y_S grid is
    # 3 x 4 at N = 1500 and 5 x 4 at N = 3000 -- and they are labelled `KX` so that no
    # N-comparison can silently include a cell that exists at only one N.
    for k in (4, 5):
        d = ladder(k)
        for ys in YS_LEVELS:
            if class_counts(d, y_from_q(d, q_for_ys(d, ys)), N_PROD).min() >= MIN_PHASE_COUNT:
                continue                      # already in the core grid
            q = q_for_ys(d, ys)
            add("KX", f"KX{k}_yS{ys:.2f}", d, y_from_q(d, q), ys=ys, kk=k, qq=q,
                judge_at=N_BIG, only_at=N_BIG)
    return cells, refused


def grid(n: int):
    """(d_rel, y, seed, n, label) per run, flattened. The runner's iteration order."""
    cells, _ = build_cells(n)
    return [(c["d_rel"], c["y"], s, n, c["label"]) for c in cells for s in c["seeds"]]


def cell_key(label: str, seed: int, n: int, mu_p: float) -> str:
    return f"kl_{label}_s{seed}_p{P_HAT:.0e}_N{n}_mu{mu_p:.2f}"


# ------------------------------------------------------------------------ the models

def pair_table(d_rel):
    """Every ORDERED pair of the ladder as (d_fine, d_coarse, u = d_fine/d_coarse)."""
    d = sorted(float(v) for v in d_rel)
    return [(d[i], d[j], d[i] / d[j]) for j in range(len(d)) for i in range(j)]


def ideal_ceiling(d_rel, y, beta: float) -> float:
    """Perfect size separation for K classes: the a = b = 0 case of the linear model.

    `mixture_model.furnas_ideal` is the binary closed form; this is the same statement
    for any K. Written with a and b SUBSTITUTED as exact zeros rather than by driving
    `packing_density`'s exponents toward zero: pa = 1e-13 gives
    a = sqrt(1 - (1-u)^pa) ~ 3e-7, not 0, and the first version of this function missed
    `furnas_ideal` by up to 3e-8 on 6 of 12 binaries -- a limit evaluated numerically is
    not the limit. lessons.md 2026-09-14: match the ARITHMETIC, not the algebra.

    With a = b = 0 the loosening term is (1 - a) = 1 for every finer class and the wall
    term is (1 - beta + 0) = (1 - beta) for every coarser one, so

        gamma_i = beta / (1 - (1 - beta) * sum_{j coarser} y_j - sum_{j finer} y_j)

    and the density is the minimum over i. Gated against `furnas_ideal` in
    `predict_kladder.py` at 1e-12 rather than claimed here.
    """
    d = np.asarray(d_rel, dtype=float)
    yv = np.asarray(y, dtype=float)
    yv = yv / yv.sum()
    order = np.argsort(-d)                     # DESCENDING: index 0 is the coarsest
    yv = yv[order]
    best = np.inf
    for i in range(yv.size):
        tot = (1.0 - beta) * float(np.sum(yv[:i])) + float(np.sum(yv[i + 1:]))
        denom = 1.0 - tot
        if denom <= 1e-12:
            continue
        best = min(best, beta / denom)
    return float(best)


def sequential_k(d_rel, y, beta: float, amap, bmap, dmode: str = "d32") -> float:
    """The RECURSIVE nested-binary rule, generalised from `grade_ternary.sequential`.

    The coarse sub-pack packs first and the finest class fills the composite. The
    sub-pack's own density is obtained by calling THIS function again, so at K = 4 the
    nesting is genuine rather than one level deep with an additive core.

    `dmode` is the COMPOSITE's effective diameter, as the moment ratio M_{k+1}/M_k:
    "d10" (k = 0), "d21", "d32", "d43". Swept, not chosen -- the K = 3 block picked d32
    and d43 by hand and the two disagreed by rms 3.80 vs 7.19.

    REDUCTIONS, wiring-gated in `predict_kladder.py` rather than asserted here:
      * K = 2 -> exactly the ADDITIVE binary closure. Not asserted by inspection: with
        one coarse class g_sub = beta and d_eff = d_coarse, so `coarse` is the linear
        model's large-class branch and `fine` is its small-class branch, term for term.
        The gate checks it numerically anyway, because that is the claim
        PREREG_TERNARY_2026-09-14 got wrong by asserting it;
      * K = 3 -> bit-for-bit `grade_ternary.sequential` with the same `dmode`.
    """
    d = np.asarray(sorted(float(v) for v in d_rel), dtype=float)
    yv = np.asarray(y, dtype=float)
    order = np.argsort(np.asarray(d_rel, dtype=float))
    yv = yv[order]
    yv = yv / yv.sum()
    if d.size < 2:
        return float(beta)
    w = float(np.sum(yv[1:]))
    if w <= 0:
        return float(beta)
    ysub = yv[1:] / w
    dsub = d[1:]
    if dsub.size == 1:
        g_sub = float(beta)
        d_eff = float(dsub[0])
    else:
        g_sub = sequential_k(dsub, ysub, beta, amap, bmap, dmode)
        nsub = ysub / dsub**3
        nsub = nsub / nsub.sum()
        k = {"d10": 0, "d21": 1, "d32": 2, "d43": 3}[dmode]
        mk = float(np.sum(nsub * dsub**k))
        mk1 = float(np.sum(nsub * dsub ** (k + 1)))
        d_eff = mk1 / mk
    u = float(d[0]) / d_eff
    a_ = float(amap(u))
    b_ = float(bmap(u))
    x = float(yv[0])
    coarse = g_sub / (1.0 - (1.0 - a_) * x)
    fine = beta / (1.0 - (1.0 - beta + b_ * beta * (1.0 - 1.0 / beta)) * (1.0 - x))
    return float(min(coarse, fine))


def hypothesis_predictions(r3: float) -> dict[str, dict[int, float]]:
    """R(K) under each registered hypothesis, calibrated on the measured K = 3 value.

    R(2) = 0 is an ANCHOR, not a measurement: K = 2 is TAB's own calibration set, so its
    residual is the in-sample noise floor by construction. Registered as such.
    """
    out = {
        "H_sat": {k: r3 for k in K_LEVELS},
        "H_K": {k: r3 * (k - 2) for k in K_LEVELS},
        "H_P": {k: r3 * (k * (k - 1) / 2 - 1) / 2.0 for k in K_LEVELS},
        "H_T": {k: r3 * math.comb(k, 3) for k in K_LEVELS},
    }
    return out
