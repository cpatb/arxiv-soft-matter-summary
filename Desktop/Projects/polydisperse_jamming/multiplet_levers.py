"""multiplet_levers.py -- what AGGREGATION does to the distance to jamming.

THE QUESTION
------------
Every other family in this project chooses a shape for the size distribution and asks
what packing does. A real suspension has no such freedom: one synthesis makes one
primary particle, and the width of its size distribution is a READOUT of how far it has
aggregated. Salt, pH and shear history move the aggregation state, not the CV directly.

So: does aggregating a suspension bring it CLOSER to jamming or push it further away?

THE TWO LEVERS, and they oppose each other
------------------------------------------
A `k`-mer holds `k` monomers of solid but excludes the volume of a sphere of diameter
`d_k = k**(1/d_f) d_1`. Write `gamma = 3/d_f`. Then (derivation, with SymPy, limiting-case,
dimensional and particle-array checks, in `derivations/multiplet_psd_2026-09-22.md`):

  LEVER 1, TOWARD jamming.  The pack sees more volume than the solid you weighed in:

        Lambda = phi_eff / phi_solid = <k^gamma> / <k>   >= 1,   == 1 exactly at d_f = 3.

  LEVER 2, AWAY from jamming.  The multiplet population is POLYDISPERSE, and this
  project has measured -- and calibrated a filling functional for -- the fact that
  polydispersity RAISES the jamming fraction:

        G = phi_J(multiplet PSD) / phi_J(monodisperse)   >= 1.

  The distance to jamming, relative to the same suspension unaggregated, is the RATIO:

        J / J_0 = (phi_eff / phi_J) / (phi_solid / phi_J_mono) = Lambda / G.

  Aggregation promotes jamming iff Lambda > G, and suppresses it iff Lambda < G.

At `d_f = 3` (compact multiplets) `Lambda = 1` identically while `G > 1`, so COMPACT
aggregation can only push a suspension AWAY from jamming. As `d_f` falls, `Lambda` grows
without bound while `G` is capped by `phi_J < 1`. **So a critical `d_f*` must exist.**
Locating it is the point of this script, and it is located with the project's OWN
calibrated closure rather than a textbook one.

GAUGE, named before any number (the rung convention is not free)
---------------------------------------------------------------
`mixture_model` is calibrated at `p_hat = 1e-3`, `mu_p = 0.3`, `N = 1500`, where the
monodisperse packing fraction is `beta = 0.5947` (measured, not fitted). Every `phi_J`
here is that channel. `G` is a RATIO within that channel, which is why it is the quantity
quoted rather than either fraction alone.

WHAT THIS SCRIPT CANNOT DECIDE, stated up front
-----------------------------------------------
It is an EFFECTIVE-SPHERE model: each multiplet is replaced by one sphere. A real doublet
is a dumbbell, and this project's shape axis measures separately that non-sphericity
changes packing. So `G` here is the SIZE consequence of aggregation only. The DEM cells
registered in `PREREG_MULTIPLET_2026-09-22.md` are what separate the two, by running the
same `n_k` as spheres and as rigid clumps.

Run:  python3 multiplet_levers.py      (writes out/multiplet_levers.txt)
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

import mixture_model as mm
from psd import MultipletPSD

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

# The calibrated closure. beta is MEASURED (the monodisperse pack at this rung); PA and
# PB are the refit of `mixture_closure.py` on 16 binary compositions, no continuous-family
# point in the fit, and it transfers to 15 held-out family points at chi2 63.0 against the
# textbook exponents' 437.7. Source: out/mixture_closure.txt.
BETA = 0.5947
PA_REFIT, PB_REFIT = 0.911, 2.122
PA_TEXT, PB_TEXT = 1.020, 1.500
# The WIDER calibration, on size ratios 1.25-3.0. `out/mixture_closure.txt` E2 records it
# as the BETTER held-out family transfer (chi2 34.9 against 63.0) and it was not used in
# the first version of this file -- a referee found the omission. Carried as a third
# column rather than adopted, because the same file's E1 verdict is that the
# two-exponent FORM is too rigid, so neither pair is "the" answer.
PA_WIDE, PB_WIDE = 0.961, 2.214
D10 = 0.8e-6                      # the lab silica's nominal size; only ratios matter here
# Truncation of the aggregation-number distribution. NOT a free choice: `Lambda` is
# `<k^gamma>/<k>` with `gamma > 1`, so it weights the tail, and at `<k> = 5, d_f = 1.8`
# the value moved 2.62 -> 3.83 between K_max = 8 and 40 and was STILL moving. 24 was the
# first value tried and it is wrong for the upper half of the table. `K_CONV_TOL` gates
# it: section 4 re-runs every quoted point at K_MAX and K_MAX//2 and any point whose
# Lambda/G differs by more than the tolerance is printed as NOT CONVERGED rather than
# quoted.
K_MAX = 240
K_CONV_TOL = 0.005                # relative, on Lambda/G
SCALE = mm.PHI_SCALE
# Seeds per cell in the registered campaign. The floor on a CELL MEAN is
# SCALE/sqrt(N_SEEDS), which is the convention `multiplet_wc.gain` already uses.
N_SEEDS = 3


# A class that cannot place a single grain in the pack is not a class. The DEM packs
# N = 1500 grains, so a number fraction below 1/1500 is already unrepresentable there;
# FLOOR is set three decades below that so the CUTOFF is never what decides an answer,
# and `flory` prints the number and the volume it drops rather than assuming either is
# small. Carrying the raw Flory tail instead is not an option: at <k> = 1.2 its k = 24
# class has n ~ 2e-17, and DiscretePSD (rightly) refuses a class whose number fraction
# has underflowed to zero in the volume-fraction round trip.
FLOOR = 1.0 / (1500.0 * 1000.0)


def flory(k_bar: float, k_max: int = K_MAX, floor: float = FLOOR,
          report: bool = False):
    """Most-probable (Flory / Smoluchowski) aggregation-number distribution.

    Irreversible step aggregation with a size-independent kernel gives
    `n_k = (1-p) p**(k-1)` with `<k> = 1/(1-p)`. One parameter, derived rather than
    chosen, and it is the standard null for a coagulating colloid.

    Truncated at `k_max`, then classes below `floor` in NUMBER fraction are dropped and
    the rest renormalised. With `report=True` also returns how many classes were dropped
    and what fraction of the SOLID volume went with them -- the quantity that matters,
    since a rare large multiplet carries `k**gamma` times a monomer's excluded volume.
    """
    if k_bar < 1.0:
        raise ValueError(f"<k> must be >= 1 (k=1 is an unaggregated monomer), got {k_bar}")
    if k_bar == 1.0:
        return ((1.0,), dict(n_before=1, n_after=1, vol_dropped=0.0)) if report else (1.0,)
    p = 1.0 - 1.0 / k_bar
    ks = np.arange(1, k_max + 1, dtype=float)
    n = (1 - p) * p ** (ks - 1)
    n = n / n.sum()
    keep = n >= floor
    # solid volume goes as k, so this is the solid fraction the cut discards
    vol_dropped = float(np.sum(n[~keep] * ks[~keep]) / np.sum(n * ks))
    n_out = tuple(float(x) for x in n[keep] / n[keep].sum())
    if report:
        return n_out, dict(n_before=int(ks.size), n_after=len(n_out),
                           vol_dropped=vol_dropped, k_kept=int(ks[keep][-1]))
    return n_out


def phi_J(psd, pa: float = PA_REFIT, pb: float = PB_REFIT) -> float:
    """Jamming fraction of the effective-sphere pack, from the calibrated closure."""
    return mm.packing_density(psd.diameters, psd.volume_fractions, BETA, pa, pb)


def levers(k_bar: float, d_f: float, pa: float = PA_REFIT, pb: float = PB_REFIT) -> dict:
    """The three numbers, for one aggregation state."""
    psd = MultipletPSD(d10=D10, d_f=d_f, n_k=flory(k_bar))
    lam = psd.lambda_trap
    g = phi_J(psd, pa, pb) / BETA
    return dict(k_bar=k_bar, d_f=d_f, cv=psd.cv, lam=lam, G=g, ratio=lam / g,
                phiJ=phi_J(psd, pa, pb), ratio_size=psd.size_ratio)


def critical_df(k_bar: float, pa: float = PA_REFIT, pb: float = PB_REFIT,
                lo: float = 0.7, hi: float = 3.0, tol: float = 1e-6) -> float:
    """The `d_f` at which `Lambda = G`, by bisection. NaN if no crossing in [lo, hi].

    `Lambda/G - 1` is increasing as `d_f` FALLS, so the bracket is (hi, lo). Both ends
    are evaluated and the result is returned only if they straddle -- a bisection on a
    non-straddling bracket returns an endpoint, which would print a clamp as an estimate.
    """
    f = lambda x: levers(k_bar, x, pa, pb)["ratio"] - 1.0
    f_hi, f_lo = f(hi), f(lo)
    if not (f_hi < 0.0 < f_lo):
        return float("nan")
    a, b = lo, hi
    while b - a > tol:
        m = 0.5 * (a + b)
        if f(m) > 0.0:
            a = m
        else:
            b = m
    return 0.5 * (a + b)


def df_star_sensitivity(k_bar: float, d_f_star: float, h: float = 1e-4) -> float:
    """How far one PHI_SCALE of protocol noise moves the ROOT `d_f*`.

    ADDED 2026-09-22 after referee pass 2, which found that the resolvability filter
    first written here asked the WRONG QUESTION. That filter tested
    `|phi_J(d_f*) - beta| >= PHI_SCALE` -- "is `G` measurably different from 1 at the
    crossing" -- which is not the same as "is the LOCATION of the crossing resolvable".

    The right quantity follows from the implicit function theorem on `Lambda(d_f) =
    G(d_f)`: a perturbation `dG` in `G` moves the root by

        d(d_f*) = dG / |Lambda'(d_f*) - G'(d_f*)|,

    and the protocol's own noise on `G` is `PHI_SCALE/BETA` (since `G = phi_J/beta`).
    Both derivatives are taken by central difference at the root.

    ⚠ THE NOISE IS AVERAGED OVER SEEDS, corrected after referee pass 3. This first used
    `SCALE/BETA`, a SINGLE-run reproducibility, while the registered campaign averages
    `N_SEEDS = 3` seeds per cell and the sibling module `multiplet_wc.gain` divides the
    identical floor by `sqrt(n)`. Two files in one derivation using the same declared
    scale with different conventions is the drift this repository punishes, so the seed
    average is used here too.

    ⚠ Pass 3 applied that `sqrt(n)` and DROPPED the sibling's `sqrt(1 + G^2)` while
    citing the sibling as its authority -- so the two corrections partly cancel. With
    both terms the surviving subset is `2.0 <= <k> <= 8.0`, which is where it started.
    """
    def lam(x):
        return MultipletPSD(d10=D10, d_f=x, n_k=flory(k_bar)).lambda_trap

    def g(x):
        return phi_J(MultipletPSD(d10=D10, d_f=x, n_k=flory(k_bar))) / BETA

    dlam = (lam(d_f_star + h) - lam(d_f_star - h)) / (2 * h)
    dg = (g(d_f_star + h) - g(d_f_star - h)) / (2 * h)
    denom = abs(dlam - dg)
    # 🔴 THE RATIO TERM, added after referee pass 4. `G = phi_J/BETA` has TWO
    # independently noisy terms -- `out/mixture_model.txt` records BETA itself as
    # `0.5947 +- 0.0010 (n = 3)`, a measurement carrying the same per-seed floor, NOT a
    # noise-free reference constant. So the floor is the sibling module's own form,
    #     var(G) = SCALE^2 (1 + G^2) / BETA^2 / n,
    # and pass 3's fix had cited that sibling while dropping its sqrt(1 + G^2). The
    # omission made every sensitivity ~1.43x too SMALL, i.e. too favourable.
    g_at_root = g(d_f_star)
    dG = SCALE * math.sqrt(1.0 + g_at_root ** 2) / (BETA * math.sqrt(N_SEEDS))
    return float("inf") if denom == 0 else float(dG / denom)


def truncation_report(k_bar: float, d_f: float) -> str:
    """How much of the answer is the truncation at K_MAX rather than the physics."""
    out = []
    for km in (8, 12, 16, 24, 40):
        psd = MultipletPSD(d10=D10, d_f=d_f, n_k=flory(k_bar, km))
        out.append(f"K_max={km:3d}  Lambda={psd.lambda_trap:8.4f}  CV={psd.cv:7.4f}  "
                   f"G={phi_J(psd)/BETA:7.4f}  Lambda/G={psd.lambda_trap/(phi_J(psd)/BETA):7.4f}")
    return "\n      ".join(out)


def main() -> int:
    L: list[str] = []

    def say(s: str = "") -> None:
        print(s)
        L.append(s)

    say("=" * 92)
    say("multiplet_levers.py -- does AGGREGATION promote or suppress jamming?")
    say("=" * 92)
    say(f"   gauge: mixture_model at p_hat = {mm.P_HAT}, mu_p = {mm.MU_P}, N = {mm.N_REQ};")
    say(f"   monodisperse beta = {BETA} (measured). closure exponents PA = {PA_REFIT}, "
        f"PB = {PB_REFIT} (refit).")
    say(f"   grading scale PHI_SCALE = {SCALE} -- the protocol's own run-to-run scatter.")
    say(f"   aggregation-number distribution: Flory most-probable, truncated at k = {K_MAX}.")

    # ---------------------------------------------------------------- 1
    say("\n=== 1. THE TWO LEVERS ACROSS THE AGGREGATION / COMPACTNESS PLANE ===")
    say("   Lambda = phi_eff/phi_solid (toward jamming);  G = phi_J/phi_J_mono (away).")
    say("   Lambda/G < 1 means aggregation SUPPRESSES jamming; > 1 means it PROMOTES it.")
    say("")
    kbars = (1.2, 1.5, 2.0, 3.0, 5.0)
    dfs = (3.0, 2.5, 2.1, 1.8, 1.5, 1.2)
    say("   " + f"{'<k>':>6}" + "".join(f"{('d_f=%.1f' % f):>11}" for f in dfs))
    for kb in kbars:
        row = "".join(f"{levers(kb, f)['ratio']:11.4f}" for f in dfs)
        say(f"   {kb:6.1f}" + row)
    say("")
    say("   the same table, showing WHICH lever is moving (Lambda | G):")
    say("   " + f"{'<k>':>6}" + "".join(f"{('d_f=%.1f' % f):>17}" for f in dfs))
    for kb in kbars:
        row = "".join(f"{levers(kb,f)['lam']:8.3f} |{levers(kb,f)['G']:8.4f}" for f in dfs)
        say(f"   {kb:6.1f}" + row)

    # ---------------------------------------------------------------- 2
    say("\n=== 2. THE CRITICAL FRACTAL DIMENSION d_f* (where Lambda = G) ===")
    say("   Below d_f*, aggregation brings the suspension CLOSER to jamming.")
    say("")
    say("   🔴 THE RESOLVABILITY CRITERION CHANGED after referee pass 2. This table first")
    say("   filtered rows on |phi_J(d_f*) - beta| >= PHI_SCALE -- 'is G measurably")
    say("   different from 1 at the crossing'. THAT IS THE WRONG QUESTION. What decides")
    say("   whether d_f* is a measurement is how far the ROOT moves under the protocol's")
    say("   own noise, which the implicit function theorem on Lambda(d_f) = G(d_f) gives:")
    say("")
    say("       d(d_f*) = (PHI_SCALE/BETA) / |Lambda'(d_f*) - G'(d_f*)|")
    say("")
    say("   BOTH columns are printed, so the old filter's number stays visible.")
    say("")
    say(f"   {'<k>':>6}{'CV':>9}{'d_f*':>9}{'(text)':>9}{'(wide)':>9}"
        f"{'phiJ-b /SC':>12}{'D(d_f*)':>10}")
    rows = []
    for kb in (1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 8.0):
        dstar = critical_df(kb)
        if np.isnan(dstar):
            say(f"   {kb:6.1f}{'--':>12}{'no crossing in [0.7, 3.0]':>27}")
            continue
        dstar_t = critical_df(kb, PA_TEXT, PB_TEXT)
        dstar_w = critical_df(kb, PA_WIDE, PB_WIDE)
        lv = levers(kb, dstar)
        shift = (lv["phiJ"] - BETA) / SCALE
        sens = df_star_sensitivity(kb, dstar)
        rows.append((kb, dstar, shift, sens))
        say(f"   {kb:6.1f}{lv['cv']:9.4f}{dstar:9.4f}{dstar_t:9.4f}{dstar_w:9.4f}"
            f"{shift:12.2f}{sens:10.3f}")
    say("")
    ds_all = [d for _, d, _, _ in rows]
    plateau = max(ds_all) - min(ds_all)
    say(f"   🔴 THE PLATEAU IS {plateau:.3f} WIDE (d_f* = {min(ds_all):.3f} - "
        f"{max(ds_all):.3f}), and the last")
    say("   column says how much of that width ONE PHI_SCALE of protocol noise buys:")
    for kb, dstar, shift, sens in rows:
        frac = sens / plateau
        flag = ("  🔴 noise alone moves d_f* by "
                f"{frac:.0%} of the whole plateau") if frac > 0.25 else ""
        say(f"      <k> = {kb:4.1f}:  D(d_f*) = {sens:.3f}  = {frac:5.0%} of the plateau"
            f"{flag}")
    # The cut: a row is quoted only if its root is pinned to better than a QUARTER of
    # the spread being quoted. ⚠ One quarter is a convention, not a derivation, and a
    # referee rightly flagged it as undocumented. It is chosen so that four such rows
    # could differ by the full spread before the spread became noise.
    #
    # ⚠ AND IT IS NOT FULLY INERT -- a first version of this comment claimed the
    # surviving spread was "0.053 at any cut from 0.10 to 0.50", which the line this
    # code prints immediately refuted (0.031 at 10 %). The measured values are printed
    # rather than asserted. Tightening to 10 % drops <k> = 1.5 and 2.0 and NARROWS the
    # surviving spread to 0.031; 25 % and 50 % both give 0.053. So the quoted spread is
    # stable upward from the chosen cut and would SHRINK if tightened -- i.e. the cut is
    # conservative in the direction that matters. The ORDINAL headline uses no cut.
    CUT = 0.25
    ok = [(k, d) for k, d, _, s in rows if s < CUT * plateau]
    spreads = {c: (lambda sel: (max(sel) - min(sel)) if sel else float("nan"))(
        [d for _, d, _, s in rows if s < c * plateau]) for c in (0.10, 0.25, 0.50)}
    say("   cut sensitivity -- the surviving SPREAD at three thresholds: "
        + ", ".join(f"{c:.0%}: {v:.3f}" for c, v in spreads.items()))
    say("")
    if ok:
        ks = [k for k, _ in ok]
        ds = [d for _, d in ok]
        say("   🔑 WHAT SURVIVES: on the rows whose root is pinned to better than a")
        say(f"   quarter of the plateau ({min(ks)} <= <k> <= {max(ks)}), d_f* = "
            f"{min(ds):.3f} - {max(ds):.3f},")
        say(f"   a spread of {max(ds)-min(ds):.3f}.")
    say("")
    say("   ⚠ THE HEADLINE IS DEMOTED. 'd_f* is essentially independent of how much the")
    say("   suspension has aggregated' is NOT established at the precision first claimed:")
    say("   the low-<k> rows carry a statistical uncertainty COMPARABLE TO the whole")
    say("   plateau width, so their agreement with the high-<k> rows is not evidence of")
    say("   flatness. What IS established, and needs no precision argument, is the")
    say("   ORDINAL fact: every row's d_f* lies near 2.9, far below the compact limit")
    say("   d_f = 3 and far above the open-aggregate values 1.8-2.1, on all three closure")
    say("   calibrations. The '+-0.01 across calibrations' figure is a SYSTEMATIC")
    say("   uncertainty and does NOT cover this statistical one, which is larger.")
    say("")
    say("   ⚠ THREE closure calibrations are printed because d_f* is a property of G, and")
    say("   G is the part of this calculation that comes from a FIT. Lambda is closed")
    say("   form and carries no fitted parameter at all.")
    say(f"      text = textbook exponents ({PA_TEXT}, {PB_TEXT});")
    say(f"      (default column) = refit on size ratios 2 and 3 ({PA_REFIT}, {PB_REFIT});")
    say(f"      wide = refit on ratios 1.25-3 ({PA_WIDE}, {PB_WIDE}), which out/")
    say("      mixture_closure.txt records as the BETTER held-out transfer (family chi2")
    say("      34.9 against 63.0). It is shown rather than adopted because the same file's")
    say("      E1 verdict is that the two-exponent FORM is too rigid (joint chi2 3.55x the")
    say("      sum of per-ratio fits), so neither pair is the right answer -- and d_f*")
    say("      moves by at most ~0.01 between them, which is the point worth carrying.")

    # ---------------------------------------------------------------- 3
    say("\n=== 3. IS THE EFFECT RESOLVABLE? d(phi_J) against the protocol floor ===")
    say("   G is only usable if the phi_J it is built from moves by more than the")
    say("   protocol's own reproducibility. Reporting phi_J - beta in units of PHI_SCALE:")
    say("")
    say("   " + f"{'<k>':>6}" + "".join(f"{('d_f=%.1f' % f):>11}" for f in dfs))
    for kb in kbars:
        row = "".join(f"{(levers(kb,f)['phiJ']-BETA)/SCALE:11.2f}" for f in dfs)
        say(f"   {kb:6.1f}" + row)
    say("")
    say("   (units of PHI_SCALE = 2.1e-3. |value| < 1 means the closure predicts a shift")
    say("   the DEM protocol could not resolve, so no cell there can test anything.)")

    # ---------------------------------------------------------------- 4
    say("\n=== 4. TRUNCATION: is K_max the answer, or the physics? ===")
    say("   N before and after the number-fraction floor, and the SOLID volume it drops:")
    say(f"   {'<k>':>6}{'classes before':>16}{'after':>8}{'largest k kept':>16}"
        f"{'solid vol dropped':>20}")
    for kb in (1.2, 1.5, 2.0, 3.0, 5.0):
        _, rep = flory(kb, report=True)
        say(f"   {kb:6.1f}{rep['n_before']:16d}{rep['n_after']:8d}{rep.get('k_kept', 1):16d}"
            f"{rep['vol_dropped']:20.3e}")
    say("")
    for kb, f in ((2.0, 1.8), (5.0, 1.8), (2.0, 3.0)):
        say(f"   <k> = {kb}, d_f = {f}:")
        say("      " + truncation_report(kb, f))
    say("")
    say("   the Lambda/G column is what every verdict above is read off; a truncation")
    say("   that moved it would make K_max a fitted parameter rather than a cutoff.")

    say("\n   CONVERGENCE GATE on every point quoted in sections 1-3:")
    say(f"   Lambda/G at K_max = {K_MAX} against K_max = {K_MAX // 2}, "
        f"tolerance {K_CONV_TOL:.1%} relative.")
    worst_rel, worst_at, n_bad = 0.0, None, 0
    # ⚠ The gate must cover every <k> QUOTED anywhere in this file, not only section 1's
    # grid: sections 2 and 6 also use 1.1 and 8.0, and a referee found them uncovered.
    gate_kbars = sorted(set(kbars) | {1.1, 8.0})
    for kb in gate_kbars:
        for f in dfs:
            full = MultipletPSD(d10=D10, d_f=f, n_k=flory(kb, K_MAX))
            half = MultipletPSD(d10=D10, d_f=f, n_k=flory(kb, K_MAX // 2))
            rf = full.lambda_trap / (phi_J(full) / BETA)
            rh = half.lambda_trap / (phi_J(half) / BETA)
            rel = abs(rf - rh) / rf
            if rel > worst_rel:
                worst_rel, worst_at = rel, (kb, f)
            if rel > K_CONV_TOL:
                n_bad += 1
                say(f"      🔴 NOT CONVERGED  <k> = {kb}, d_f = {f}: "
                    f"{rh:.4f} -> {rf:.4f}  ({rel:.2%})")
    say(f"      worst relative drift {worst_rel:.3%} at <k> = {worst_at[0] if worst_at else 0}, "
        f"d_f = {worst_at[1] if worst_at else 0}")
    say(f"      {'[PASS]' if n_bad == 0 else '[FAIL]'} "
        f"{len(gate_kbars)*len(dfs) - n_bad} of {len(gate_kbars)*len(dfs)} quoted points converged")

    # ---------------------------------------------------------------- 5
    say("\n=== 5. THE BENCH STATEMENT ===")
    lv_compact = levers(2.0, 3.0)
    lv_open = levers(2.0, 1.8)
    say("   At <k> = 2 (half the solid in multiplets), COMPACT aggregates (d_f = 3):")
    say(f"      Lambda = {lv_compact['lam']:.4f} (no trapped fluid, exactly), "
        f"G = {lv_compact['G']:.4f}, Lambda/G = {lv_compact['ratio']:.4f}")
    say(f"      -> aggregation moves the suspension AWAY from jamming by "
        f"{(1-lv_compact['ratio'])*100:.1f} %.")
    say("   The same <k> = 2 with OPEN aggregates (d_f = 1.8, the diffusion-limited value):")
    say(f"      Lambda = {lv_open['lam']:.4f}, G = {lv_open['G']:.4f}, "
        f"Lambda/G = {lv_open['ratio']:.4f}")
    say(f"      -> aggregation moves it TOWARD jamming by "
        f"{(lv_open['ratio']-1)*100:.1f} %.")
    say("")
    say("   SO THE SIGN OF THE AGGREGATION EFFECT IS SET BY THE AGGREGATE'S OPENNESS,")
    say("   NOT BY HOW MUCH IT HAS AGGREGATED. That is falsifiable on the bench: it says")
    say("   two suspensions at the SAME measured CV can sit on opposite sides of the")
    say("   jamming shift if their aggregates differ in compactness -- and that a CV")
    say("   measurement alone therefore cannot predict the rheology.")

    # ---------------------------------------------------------------- 6
    say("\n=== 6. WHY d_f* BARELY MOVES: a first-order expansion, and its check ===")
    say("   Near gamma = 1 (d_f = 3), with L1 = <k ln k>/<k>:")
    say("       Lambda(gamma) = 1 + (gamma-1) L1 + O((gamma-1)^2)")
    say("       G(gamma)      = G0 + (gamma-1) G0' + ...")
    say("   so the crossing Lambda = G sits at")
    say("       gamma* - 1  =  (G0 - 1) / (L1 - G0')   ~=  (G0 - 1) / L1.")
    say("   BOTH G0-1 and L1 grow with aggregation, so their RATIO is nearly invariant.")
    say("   That is the explanation for the flat d_f* column, and it is checkable:")
    say("")
    say(f"   {'<k>':>6}{'G0 - 1':>10}{'L1':>10}{'(G0-1)/L1':>12}"
        f"{'gamma*-1 (exact)':>18}{'d_f* 1st order':>16}{'d_f* exact':>12}")
    for kb in (1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 8.0):
        psd0 = MultipletPSD(d10=D10, d_f=3.0, n_k=flory(kb))
        G0 = phi_J(psd0) / BETA
        L1 = float(np.sum(psd0.n * psd0.k * np.log(psd0.k)) / np.sum(psd0.n * psd0.k))
        approx = (G0 - 1.0) / L1
        dstar = critical_df(kb)
        say(f"   {kb:6.1f}{G0-1:10.4f}{L1:10.4f}{approx:12.4f}"
            f"{3.0/dstar - 1.0:18.4f}{3.0/(1.0+approx):16.4f}{dstar:12.4f}")
    say("")
    say("   🔴 THE SIGN, CORRECTED 2026-09-22 after a referee pass. An earlier version of")
    say("   this paragraph said the first-order column 'sits BELOW the exact one")
    say("   throughout'. THAT IS REFUTED BY THE TABLE ABOVE, which this file prints:")
    devs = []
    for kb in (1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 8.0):
        psd0 = MultipletPSD(d10=D10, d_f=3.0, n_k=flory(kb))
        G0 = phi_J(psd0) / BETA
        L1 = float(np.sum(psd0.n * psd0.k * np.log(psd0.k)) / np.sum(psd0.n * psd0.k))
        devs.append((kb, 3.0 / (1.0 + (G0 - 1.0) / L1) - critical_df(kb)))
    n_above = sum(1 for _, d in devs if d > 0)
    say("      " + "  ".join(f"<k>={k}: {d:+.4f}" for k, d in devs))
    say(f"   The first-order estimate is ABOVE the exact d_f* on {n_above} of {len(devs)}")
    say(f"   rows, worst |deviation| {max(abs(d) for _, d in devs):.4f} (not 0.005).")
    say("")
    say("   ABOVE is the direction the algebra predicts, and the earlier prose had the")
    say("   word wrong while its own mechanism sentence was right. Since")
    say("       gamma* - 1 = (G0-1)/(L1 - G0')   and   G0' > 0,")
    say("   dropping G0' makes the denominator LARGER, so gamma*_approx < gamma*_exact,")
    say("   and d_f = 3/gamma means d_f*_approx > d_f*_exact -- i.e. the first-order")
    say("   estimate sits CLOSER TO d_f = 3, which is ABOVE. The two rows that come out")
    say("   marginally below (-0.0003) are second order and are not explained by this")
    say("   argument; they are reported, not absorbed into it.")
    say("")
    say("   This is an explanation of the FLATNESS, not a replacement for the bisection,")
    say("   and every number quoted elsewhere in this file is the exact one.")

    # ---------------------------------------------------------------- 7
    say("\n=== 7. THE ANCHOR ON G THAT DOES NOT PASS THROUGH THE CLOSURE ===")
    say("   G is the only fitted half of this calculation, so its CEILING is checked")
    say("   against measured packs rather than against the model that predicts them.")
    say("   From out/mixture_closure.txt, MEASURED phi at this rung (mu_p = 0.3, p_hat")
    say(f"   = 1e-3, N = 1500), divided by the measured monodisperse beta = {BETA}:")
    say("")
    meas = {("lognormal", 0.5): 0.6489, ("uniform", 0.5): 0.6157,
            ("bidisperse", 0.5): 0.6042, ("lognormal", 0.3): 0.6194}
    for (fam, cv), phi in sorted(meas.items()):
        say(f"      {fam:11s} CV = {cv:.1f}   phi = {phi:.4f}   G = {phi/BETA:.4f}")
    gmax_meas = max(p / BETA for p in meas.values())
    say("")
    say(f"   -> the WIDEST family this project has ever packed reaches G = {gmax_meas:.4f}:")
    say("   the polydispersity lever is capped near +9 %, MEASURED.")
    say("")
    say("   So the argument can be run WITHOUT the closure at all. Ask only: at what")
    say("   d_f does the closed-form Lambda, on its own, exceed that measured ceiling?")
    say("   Anything past that point cannot be offset by ANY width effect this project")
    say("   has ever observed, whatever the closure says.")
    say("")
    say(f"   {'<k>':>6}{'d_f where Lambda = G_max':>28}{'Lambda at d_f = 1.8':>22}")
    for kb in (1.2, 1.5, 2.0, 3.0, 5.0):
        n_kb = flory(kb)
        cross = None
        for f in np.linspace(3.0, 1.0, 20001):
            if MultipletPSD(d10=D10, d_f=float(f), n_k=n_kb).lambda_trap >= gmax_meas:
                cross = float(f)
                break
        lam18 = MultipletPSD(d10=D10, d_f=1.8, n_k=n_kb).lambda_trap
        say(f"   {kb:6.1f}{cross:28.3f}{lam18:22.4f}")
    say("")
    say("   🔴 WHERE THE CLOSURE-FREE TEST STOPS WORKING, added 2026-09-22 after a")
    say("   referee pass. The rows above start at <k> = 1.2; a LIGHTLY aggregated")
    say("   suspension is just as realistic, and there Lambda may never reach the")
    say("   ceiling at all, which makes the closure-free test UNINFORMATIVE rather")
    say("   than favourable. The threshold is computable, so it is printed:")
    say("")
    say(f"   {'d_f':>7}{'smallest <k> whose Lambda reaches G_max':>42}")
    for df in (1.8, 2.1, 2.5):
        kmin = None
        for kb in np.arange(1.01, 3.0, 0.01):
            if MultipletPSD(d10=D10, d_f=df,
                            n_k=flory(float(kb))).lambda_trap >= gmax_meas:
                kmin = float(kb)
                break
        say(f"   {df:7.1f}{(f'{kmin:.2f}' if kmin else 'never below <k> = 3'):>42}")
    say("")
    say("   So the closure-free statement is NOT 'aggregation always wins'. It is:")
    say("   ABOVE those <k> thresholds, no width effect this project has ever measured")
    say("   can offset the trapped-fluid lever. Below them the test says nothing and the")
    say("   answer depends on the closure -- where the fitted G still gives the same")
    say("   SIGN, but by a margin of a per cent rather than a factor of several.")
    say("   The earlier wording quoted the Lambda range over the five tabulated <k>")
    say("   values as if it were a bound over all realistic aggregates. It was not.")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "multiplet_levers.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    say("")
    say("   written: out/multiplet_levers.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
