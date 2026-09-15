"""predict_kladder.py -- BLIND predictions for PREREG_KLADDER_2026-09-15.md.

RUN BEFORE ANY CELL EXISTS. Writes out/predict_kladder.txt, which the prereg quotes and
`grade_kladder.py` reproduces from disk before grading anything against it.

Nothing here reads a `kl*` record. It reads only the binary calibration set
(`mixr*`, `extr*`) and the closed ternary block (`tern*`), both of which are on disk and
CLOSED. If a `phase12_kl*` file existed this script would still not look at it.

Sections
--------
0  WIRING GATES -- every model path this block compares is checked to reduce to the path
   it claims to reduce to, numerically, at 1e-12. PREREG_TERNARY_2026-09-14 asserted one
   such reduction "by inspection" and it raised IndexError; nothing here is asserted.
1  the design: ladders, pair tables, admissibility, orthogonality of the two axes
2  the in-sample residuals the bars come from, re-derived and asserted
3  BLIND PREDICTIONS, one row per cell
4  candidate separation, and the 1.6x informativeness rule
5  the hypothesis table: what each of H_sat / H_K / H_P / H_T predicts, and which pairs
   of them this design can and CANNOT tell apart
"""

from __future__ import annotations

import os

import numpy as np

import grade_mixture_ext as gx
import grade_ternary as gt
import kladder as kl
import mixture_closure as mc
import mixture_model as mm
import phase12_ternary as pt
from psd import parse_discrete_family

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = kl.PHI_SCALE
MU = kl.MU_P

# Declared in PREREG_TERNARY_2026-09-14 and re-derived here. If either moves, the
# calibration set on disk has changed and every bar in this block is stale.
INSAMPLE_EXPECT = {("J", "own"): (2.83, 6.69), ("TAB", "all"): (2.02, 5.74)}
TAB_RMS_BAR = 2.02          # the bar's ORIGIN: TAB's own in-sample rms on the binaries
TAB_MAX_BAR = 5.74
SEP_BAR = 1.6               # below this a comparison is declared NOT INFORMATIVE


def main() -> int:
    lines: list[str] = []
    failures: list[str] = []

    def say(s: str = "") -> None:
        lines.append(s)
        print(s)

    meas = mm.measured(mu_p=MU)
    beta = meas[("lognormal", 0.0)][0]
    pts23 = mc.binary_points(mu_p=MU)
    allpts = gx.ext_points(beta)
    cand, _per, _ = gx.candidates(beta, pts23)
    paJ, pbJ = cand["J"](0)
    tab, fa, fb = gx.tabulated_closure(beta, allpts)

    def aJ(u):
        return float(np.sqrt(1 - (1 - u) ** paJ))

    def bJ(u):
        return float(1 - (1 - u) ** pbJ)

    def dens_tab(d, y):
        return float(mc.packing_density_tabulated(np.asarray(d, float),
                                                  np.asarray(y, float), beta, fa, fb))

    def dens_J(d, y):
        return float(mm.packing_density(np.asarray(d, float), np.asarray(y, float),
                                        beta, paJ, pbJ))

    say("=== predict_kladder.py -- PREREG_KLADDER_2026-09-15.md (BLIND) ===")
    say(f"   beta = {beta:.6f} (measured, mu_p {MU}); PHI_SCALE = {SCALE:.1e}")
    say(f"   J  PA {paJ:.4f} PB {pbJ:.4f};  TAB knots at u = "
        f"{', '.join(f'{u:.3f}' for u in sorted(tab))}")
    say(f"   calibration set: {len(pts23)} binary points at ratios "
        f"{sorted({p[0] for p in pts23})}; TAB fitted on {len(allpts)} incl. r = 4, 5")
    say()

    # ============================================================ 0. WIRING GATES
    say("=== 0. WIRING GATES: every model path reduces to the path it claims to ===")
    say("   (PREREG_TERNARY_2026-09-14 asserted one of these by inspection and it raised")
    say("    IndexError. Nothing below is asserted; all of it is run.)")
    nw = ok = 0

    # W1 -- the K-class ideal ceiling reduces to the binary Furnas closed form.
    say("   W1  ideal_ceiling(K=2) == mixture_model.furnas_ideal")
    for r in (2.0, 3.0, 4.0):
        for x in (0.05, 0.15, 0.30, 0.50):
            a = kl.ideal_ceiling((1.0, r), (x, 1.0 - x), beta)
            b = mm.furnas_ideal(x, beta)
            nw += 1
            ok += abs(a - b) < 1e-12
            if abs(a - b) >= 1e-12:
                failures.append(f"W1 r={r} x={x}: {a!r} vs {b!r}")
    say(f"       {ok} of {nw} to 1e-12 over r = 2, 3, 4 and x = 0.05-0.50")

    # W2 -- sequential_k at K = 2 IS the additive binary closure (both closures).
    n2 = o2 = 0
    for r in (2.0, 3.0, 4.0):
        for x in (0.05, 0.15, 0.30, 0.50):
            for nm, dn, am, bm in (("TAB", dens_tab, fa, fb), ("J", dens_J, aJ, bJ)):
                s = kl.sequential_k((1.0, r), (x, 1.0 - x), beta, am, bm, "d32")
                a = dn((1.0, r), (x, 1.0 - x))
                n2 += 1
                o2 += abs(s - a) < 1e-12
                if abs(s - a) >= 1e-12:
                    failures.append(f"W2 {nm} r={r} x={x}: {s!r} vs {a!r}")
    say(f"   W2  sequential_k(K=2) == the additive binary closure: {o2} of {n2} to 1e-12")

    # W3 -- sequential_k at K = 3 is BIT-FOR-BIT grade_ternary.sequential, which is the
    # rule the closed ternary block's T4 was graded on. Same compositions, both
    # conventions, both closures.
    n3 = o3 = 0
    worst3 = 0.0
    for comp in pt.COMPOSITIONS:
        y_l, y_m, y_s = comp
        for dmode in ("d32", "d43"):
            for nm, closure, am, bm in (("TAB", "TAB", fa, fb), ("J", "J", aJ, bJ)):
                a = gt.sequential(comp, dmode, closure, beta, paJ, pbJ, fa, fb)
                s = kl.sequential_k((1.0, 2.0, 4.0), (y_s, y_m, y_l), beta, am, bm, dmode)
                n3 += 1
                worst3 = max(worst3, abs(s - a))
                o3 += abs(s - a) < 1e-12
                if abs(s - a) >= 1e-12:
                    failures.append(f"W3 {nm} {dmode} {comp}: {s!r} vs {a!r}")
    say(f"   W3  sequential_k(K=3) == grade_ternary.sequential: {o3} of {n3} to 1e-12 "
        f"(worst |diff| {worst3:.2e})")

    # W4 -- the family strings this block will write ROUND-TRIP to the same floats, and
    # the two carried over from the ternary block are BYTE-IDENTICAL to what
    # phase12_ternary.y_string emitted, so at the same seed they are the same packing.
    n4 = o4 = 0
    for c in kl.build_cells(kl.N_BIG)[0]:
        got = parse_discrete_family(c["family"])
        assert got is not None
        d_got, y_got = got
        n4 += 1
        o4 += (tuple(d_got) == tuple(c["d_rel"]) and tuple(y_got) == tuple(c["y"]))
    say(f"   W4  family string round-trips to the SAME floats: {o4} of {n4} exact")
    mono_s = pt.y_string((1.0, 0.0, 0.0))
    na2_s = pt.y_string((0.55, 0.35, 0.10))
    ours = {c["label"]: c["family"] for c in kl.build_cells()[0]}
    for lbl, want in (("mono", mono_s), ("NA2_x6", na2_s)):
        same = ours[lbl] == want
        say(f"   W5  {lbl:>7} string == phase12_ternary.y_string: {same}  {ours[lbl]}")
        if not same:
            failures.append(f"W5 {lbl}: {ours[lbl]!r} != {want!r}")
    say("       -> those two cells are the SAME packing as the closed block's at the same")
    say("          seed, so seeds 1-3 must reproduce them BIT-IDENTICALLY. That is the")
    say("          cross-block gate; it costs nothing because the cells were needed anyway.")
    say(f"   WIRING: {'PASS' if not failures else 'FAIL ' + str(failures)}")
    say()

    # ============================================================ 1. THE DESIGN
    say("=== 1. the design: a geometric ladder of FIXED SPAN 4, filled in ===")
    say(f"   {'K':>3}{'P=pairs':>9}{'T=triples':>11}  diameters                        "
        f"pair u values")
    for k in kl.K_LEVELS:
        d = kl.ladder(k)
        us = sorted({round(u, 4) for _, _, u in kl.pair_table(d)})
        npair = k * (k - 1) // 2
        ntrip = k * (k - 1) * (k - 2) // 6
        say(f"   {k:3d}{npair:9d}{ntrip:11d}  "
            f"{', '.join(f'{v:.4f}' for v in d):32s} "
            f"{', '.join(f'{v:.4f}' for v in us)}")
    umin = min(u for k in kl.K_LEVELS for _, _, u in kl.pair_table(kl.ladder(k)))
    say(f"   worst pair ratio anywhere on the ladder: u = {umin:.4f} >= 1/4. That is the")
    say("   constraint that refused d = 1, 2, 8 for the ternary block, and span 4 holds it")
    say("   at EVERY K by construction, so additivity stays unconfounded with extrapolation.")
    say()
    say("   interaction values at each pair u (INTERPOLATED, never extrapolated):")
    say(f"   {'u':>8}{'a_J':>9}{'b_J':>9}{'a_TAB':>9}{'b_TAB':>9}{'|da|':>8}{'|db|':>8}"
        "   nearest TAB knots")
    knots = sorted(tab)
    for u in sorted({round(u, 6) for k in kl.K_LEVELS
                     for _, _, u in kl.pair_table(kl.ladder(k))}):
        lo = max([q for q in [0.0] + knots if q <= u], default=0.0)
        hi = min([q for q in knots + [1.0] if q >= u], default=1.0)
        say(f"   {u:8.4f}{aJ(u):9.4f}{bJ(u):9.4f}{float(fa(u)):9.4f}{float(fb(u)):9.4f}"
            f"{abs(aJ(u) - float(fa(u))):8.4f}{abs(bJ(u) - float(fb(u))):8.4f}"
            f"   [{lo:.3f}, {hi:.3f}]"
            + ("   <-- WIDEST knot gap; K=5's adjacent pair leans on it" if hi - lo > 0.3
               else ""))
    say()

    say("=== 1b. admissibility, applied BEFORE the run (MIN_PHASE_COUNT = "
        f"{kl.MIN_PHASE_COUNT} on EVERY class) ===")
    for n in (kl.N_PROD, kl.N_BIG):
        cells, refused = kl.build_cells(n)
        runs = sum(len(c["seeds"]) for c in cells)
        say(f"   --- N = {n}: {len(cells)} cells, {runs} runs, {len(refused)} REFUSED")
        say(f"   {'arm':>6}{'label':>14}{'K':>3}{'y_S':>7}{'q':>8}"
            f"{'n_i per class':>34}{'d_max/L':>9}{'CV':>7}{'d32':>7}{'d43':>7}")
        for c in cells:
            p = kl.psd_for(c["d_rel"], c["y"])
            dm = kl.d_max_over_L(c["d_rel"], c["y"], n, beta)
            ys = c["y_S_target"]
            flag = "  <-- d_max/L > 0.25" if dm > 0.25 else ""
            say(f"   {c['arm']:>6}{c['label']:>14}{c['K']:3d}"
                f"{(float('nan') if ys is None else ys):7.2f}"
                f"{(float('nan') if c['q'] is None else c['q']):8.3f}"
                f"{str(c['counts']):>34}{dm:9.4f}{p.cv:7.4f}{p.d32:7.3f}{p.d43:7.3f}{flag}")
        for r in refused:
            say(f"   REFUSED {r['label']:>13} {str(r['counts_judged']):>30}   "
                f"{r['why_refused']}")
    say()

    say("=== 1c. the two axes are UNCROSSED (this is the whole point of the design) ===")
    cells = kl.build_cells()[0]
    bal = [c for c in cells if c["arm"] == "K" and c["y_S_target"] in kl.YS_BALANCED]
    kk = np.array([c["K"] for c in bal], float)
    yy = np.array([c["y_S_target"] for c in bal], float)
    r_bal = float(np.corrcoef(kk, np.log(yy))[0, 1])
    core = [c for c in cells if c["arm"] == "K"]
    r_core = float(np.corrcoef([c["K"] for c in core],
                               np.log([c["y_S_target"] for c in core]))[0, 1])
    say(f"   balanced block ({len(bal)} cells, K = 2,3,4,5 x y_S = "
        f"{kl.YS_BALANCED}): corr(K, ln y_S) = {r_bal:+.6f}")
    say(f"   whole core grid ({len(core)} cells, K = 2 and 3 also carry y_S = 0.13, 0.18):"
        f" corr(K, ln y_S) = {r_core:+.4f}")
    say(f"   y_S range: {min(kl.YS_LEVELS):.2f} to {max(kl.YS_LEVELS):.2f}, a factor "
        f"{max(kl.YS_LEVELS) / min(kl.YS_LEVELS):.1f}; the closed ternary block had TWO "
        "levels, 0.10 and 0.15 (1.5x).")
    say("   `n_L` is NOT held fixed and is not claimed to be -- it cannot be, since")
    say("   n goes as y/d^3. What breaks that confound here is the N = 3000 block, which")
    say("   runs the IDENTICAL cell list and moves every count by exactly 2x.")
    say()

    # ============================================================ 2. THE BARS
    say("=== 2. in-sample residuals on the binary calibration set -- where the bars "
        "come from ===")
    for name, fn, setname, pset, key in (
            ("J", dens_J, "its OWN calibration set, ratios <= 3", pts23, ("J", "own")),
            ("TAB", dens_tab, "all binaries (it is fitted per ratio)", allpts,
             ("TAB", "all"))):
        res = np.array([(fn((r, 1.0), (1.0 - x, x)) - phi) / SCALE for r, x, phi in pset])
        rms = float(np.sqrt(np.mean(res**2)))
        mx = float(np.max(np.abs(res)))
        want = INSAMPLE_EXPECT[key]
        good = abs(rms - want[0]) < 0.02 and abs(mx - want[1]) < 0.02
        say(f"   {name:>4} on {setname}: n = {res.size}, rms {rms:.2f} (prereg "
            f"{want[0]}), max {mx:.2f} (prereg {want[1]}) -> {'OK' if good else 'MOVED'}")
        if not good:
            failures.append(f"in-sample {key} moved: {rms:.3f}/{mx:.3f} vs {want}")
    say(f"   BARS, unchanged and taken from TAB's own in-sample scatter: rms "
        f"{TAB_RMS_BAR}, max {TAB_MAX_BAR}. Nothing in this block refits anything.")
    say()

    # ============================================================ 3. PREDICTIONS
    say("=== 3. BLIND PREDICTIONS (phi), one row per cell ===")
    say("   ADD = the additive closure, TAB primary and J control. SEQ = the recursive")
    say("   nested rule of kladder.sequential_k, its composite diameter convention SWEPT")
    say("   over all four moment ratios rather than chosen. IDEAL = perfect separation.")
    say(f"   {'label':>14}{'K':>3}{'y_S':>6}{'ADDtab':>9}{'ADDj':>8}"
        f"{'SEQd10':>8}{'SEQd21':>8}{'SEQd32':>8}{'SEQd43':>8}{'IDEAL':>8}"
        f"{'head':>7}{'sep':>7}")
    say("=== BLIND PREDICTIONS ===")
    preds: dict[str, dict[str, float]] = {}
    for c in kl.build_cells(kl.N_BIG)[0]:
        d, y = c["d_rel"], c["y"]
        row = {"ADDtab": dens_tab(d, y), "ADDj": dens_J(d, y),
               "IDEAL": kl.ideal_ceiling(d, y, beta)}
        for dm in ("d10", "d21", "d32", "d43"):
            row[f"SEQ{dm}"] = kl.sequential_k(d, y, beta, fa, fb, dm)
        row["head"] = (row["IDEAL"] - row["ADDtab"]) / SCALE
        row["sep"] = abs(row["SEQd32"] - row["ADDtab"]) / SCALE
        preds[c["label"]] = row
        say(f"   {c['label']:>14}{c['K']:3d}"
            f"{(float('nan') if c['y_S_target'] is None else c['y_S_target']):6.2f}"
            f"{row['ADDtab']:9.4f}{row['ADDj']:8.4f}{row['SEQd10']:8.4f}"
            f"{row['SEQd21']:8.4f}{row['SEQd32']:8.4f}{row['SEQd43']:8.4f}"
            f"{row['IDEAL']:8.4f}{row['head']:7.1f}{row['sep']:7.1f}")
    say("=== END BLIND PREDICTIONS ===")
    say()

    # ============================================================ 4. SEPARATION
    say("=== 4. candidate separation, and the 1.6x rule ===")
    say("   A comparison whose candidates predict nearly the same thing has decided")
    say("   nothing whichever way it lands (lessons.md 2026-08-30). Separations below,")
    say("   in PHI_SCALE units, over the cells each comparison is made on.")
    grp = {}
    for c in kl.build_cells(kl.N_BIG)[0]:
        if c["arm"] in ("K",):
            grp.setdefault(c["K"], []).append(c["label"])
    say(f"   {'comparison':>26}{'mean |sep|':>12}{'min':>8}{'max':>8}   verdict")
    for nm, ka, kb in (("ADDtab vs SEQd32", "ADDtab", "SEQd32"),
                       ("ADDtab vs SEQd43", "ADDtab", "SEQd43"),
                       ("ADDtab vs SEQd10", "ADDtab", "SEQd10"),
                       ("ADDtab vs ADDj", "ADDtab", "ADDj"),
                       ("ADDtab vs IDEAL", "ADDtab", "IDEAL")):
        v = np.array([abs(preds[l][ka] - preds[l][kb]) / SCALE
                      for ls in grp.values() for l in ls])
        verdict = ("informative" if float(v.mean()) >= TAB_RMS_BAR * SEP_BAR
                   else "NOT INFORMATIVE (below 1.6 x the bar)")
        say(f"   {nm:>26}{float(v.mean()):12.2f}{float(v.min()):8.2f}"
            f"{float(v.max()):8.2f}   {verdict}")
    say()
    say("   the SPLIT control -- a K = 4 spec that is physically a K = 3 one:")
    for ys in (0.05, 0.08):
        a = preds[f"SPLIT_yS{ys:.2f}"]["ADDtab"]
        b = preds[f"K3_yS{ys:.2f}"]["ADDtab"]
        say(f"      y_S {ys:.2f}: ADDtab moves {abs(a - b) / SCALE:5.2f} scale when the "
            f"middle class is split in two (K 3->4, pairs 3->6)")
    say("      So the MODEL barely moves. If the measured residual jumps with the split,")
    say("      the residual is counting classes; if it does not, the residual is a")
    say("      functional of the distribution. Three cells, and it cannot be argued.")
    say()

    # ============================================================ 5. HYPOTHESES
    say("=== 5. the hypotheses, and which pairs this design can tell apart ===")
    r3_prior = 7.53        # the closed ternary block's measured mean residual, N = 1500
    say(f"   R3 is calibrated from the closed block's measured mean residual, {r3_prior}")
    say("   scale (PREREG_TERNARY_2026-09-14 T1). R(2) = 0 is an ANCHOR, not a blind")
    say("   prediction: K = 2 IS TAB's calibration set, so its residual is the in-sample")
    say("   noise floor by construction. Stated here so nobody reads it as a success.")
    hyp = kl.hypothesis_predictions(r3_prior)
    say(f"   {'hypothesis':>8}" + "".join(f"{'K=' + str(k):>9}" for k in kl.K_LEVELS))
    for nm, dd in hyp.items():
        say(f"   {nm:>8}" + "".join(f"{dd[k]:9.2f}" for k in kl.K_LEVELS))
    say()
    say(f"   pairwise separation RATIOS (larger / smaller), bar {SEP_BAR}:")
    names = list(hyp)
    say(f"   {'pair':>16}{'at K=4':>9}{'at K=5':>9}   verdict")
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            r4 = max(hyp[a][4], hyp[b][4]) / max(min(hyp[a][4], hyp[b][4]), 1e-9)
            r5 = max(hyp[a][5], hyp[b][5]) / max(min(hyp[a][5], hyp[b][5]), 1e-9)
            v = "separable" if max(r4, r5) >= SEP_BAR else "NOT SEPARABLE by this design"
            say(f"   {a + ' vs ' + b:>16}{r4:9.2f}{r5:9.2f}   {v}")
    say()
    say("   REGISTERED IN ADVANCE: H_K and H_P are NOT separable here (1.25x at K = 4,")
    say("   1.50x at K = 5, both below 1.6). No verdict of this block may claim to have")
    say("   distinguished 'per class' from 'per pair'. What it CAN decide is")
    say("   saturation vs growth, and linear-in-K vs three-body -- and, through the y_S")
    say("   axis at fixed K, whether the residual is a composition effect at all.")
    say()
    say("   And the falsifiable prediction the y_S axis carries on its own: the closed")
    say("   block found spearman(res, y_S) = +0.878 over TWO levels. Registered: at K = 3")
    say("   over FIVE levels spanning 6x, the residual must be MONOTONE INCREASING in")
    say("   y_S (spearman >= +0.7, and the 0.18 cell above the 0.03 cell by at least the")
    say(f"   {TAB_RMS_BAR} bar). If it is flat, the two-level correlation was a fluke of")
    say("   the corner and the composition reading of T1 fails.")
    say()

    say("=== PREDICTION FILE COMPLETE -- no kl* record was read to produce it ===")
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "predict_kladder.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    if failures:
        print("\nWIRING/BAR FAILURES:", failures)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
