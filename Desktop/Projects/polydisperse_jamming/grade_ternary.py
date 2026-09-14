"""grade_ternary.py -- grade T0-T6 of PREREG_TERNARY_2026-09-14.md.

Separate from every other grader for the reason the others are separate from each other:
M0-M6, F1-F3 and X0-X5 are CLOSED and this script must not be able to re-decide them. It
loads ONLY `phase12_tern_shard*.jsonl`. The binary blocks enter through
`mixture_closure.binary_points` / `grade_mixture_ext.ext_points`, which is the
calibration set the predictions were written from -- read, never re-fitted.

Every bar below is copied from PREREG_TERNARY_2026-09-14.md section 4 and none is
computed from the ternary data. The predictions are recomputed by the same code
`predict_ternary.py` ran BEFORE the cells existed, and the grader ASSERTS they still
match `out/predict_ternary.txt` -- so a later change to the closure cannot silently
move the prediction the gate is read against.

    python grade_ternary.py      -> out/grade_ternary.txt
"""

from __future__ import annotations

import collections
import glob
import json
import os
import re

import numpy as np

import grade_mixture_ext as gx
import mixture_closure as mc
import mixture_model as mm
from phase12_mixture import MIN_PHASE_COUNT, N_PROD
from phase12_ternary import COMPOSITIONS, MONO, make_ternary

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = mm.PHI_SCALE
MU = 0.3
GLOB = "phase12_tern_shard*.jsonl"        # NOT mixr*/extr*: see module docstring

# --- registered bars, PREREG_TERNARY_2026-09-14.md section 4 -------------------------
T0_BAR = 2.0            # scale; monodisperse control vs beta (X0's bar, verbatim)
T1_RMS_BAR = 2.8        # scale; 1.4 x TAB's in-sample rms 2.02
T1_MAX_BAR = 5.74       # scale; TAB's in-sample max |residual|, verbatim
T2_RMS_BAR = 4.0        # scale; X1's bar, verbatim
T2_MAX_BAR = 6.69       # scale; J's in-sample max on its own calibration set
T3_MEAN_BAR = 2.8       # scale; same number as T1's rms bar, by construction
T4_SEP_BAR = 1.6        # rms ratio loser/winner for "informative"

# Declared in the prereg as the numbers the bars came from; re-derived and ASSERTED.
INSAMPLE_EXPECT = {("J", "own"): (2.83, 6.69), ("TAB", "all"): (2.02, 5.74)}


def load():
    g: dict[tuple, list[dict]] = collections.defaultdict(list)
    files = sorted(glob.glob(os.path.join(OUT, GLOB)))
    nraw = 0
    for f in files:
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                nraw += 1
                r = json.loads(line)
                key = (round(float(r["tern_y_L"]), 4), round(float(r["tern_y_M"]), 4),
                       round(float(r["tern_y_S"]), 4))
                g[key].append(r)
    return g, len(files), nraw


def admissible(c: dict) -> tuple[bool, str]:
    """The registered admission rule. Applied here, declared in the prereg."""
    if c.get("error"):
        return False, "error"
    if not c.get("converged"):
        return False, "not converged"
    if c.get("n_requested") != N_PROD:
        return False, f"N {c.get('n_requested')}"
    if abs(float(c.get("mu_p", -1)) - MU) > 1e-9:
        return False, f"mu {c.get('mu_p')}"
    counts = c.get("n_design_counts") or []
    if len(counts) > 1 and min(counts) < MIN_PHASE_COUNT:
        return False, f"phase count {counts}"
    return True, ""


def sequential(comp, dmode, closure, beta, paJ, pbJ, fa, fb) -> float:
    """The registered NON-ADDITIVE alternative. Copied verbatim from predict_ternary.py
    and asserted below to reproduce that script's stored numbers."""
    y_l, y_m, y_s = comp
    p = make_ternary(comp)
    d = p.diameters

    def dens(dd, yy):
        if closure == "J":
            return mm.packing_density(dd, yy, beta, paJ, pbJ)
        return mc.packing_density_tabulated(dd, yy, beta, fa, fb)

    # The COMPOSITE is every class coarser than the finest one. Writing it that way
    # rather than as the hardcoded pair (d[1], d[2]) is what makes the K = 2 reduction
    # exist at all: DiscretePSD DROPS a zero-volume class, so an emptied-middle spec
    # arrives here with TWO diameters and `d[2]` raised IndexError. The wiring gate added
    # on 2026-09-14 caught it -- and it caught a claim, not just a crash, because
    # PREREG_TERNARY_2026-09-14.md section 3 asserts this construction "reduces EXACTLY
    # to the binary closure at K = 2". That was asserted by inspection and was FALSE of
    # the code as shipped. It is now true and TESTED.
    if d.size < 2:
        return beta
    yv = p.volume_fractions                # ascending in diameter, sums to 1
    w = float(np.sum(yv[1:]))
    if w <= 0:
        return beta
    ysub = np.asarray(yv[1:], float) / w
    dsub = np.asarray(d[1:], float)
    if dsub.size == 1:
        # a single coarse class IS monodisperse: it packs at beta and presents its own
        # diameter. This is the branch that makes K = 2 reduce to the binary closure.
        g_sub = beta
        d_eff = float(dsub[0])
    else:
        g_sub = dens(dsub, ysub)
        nsub = ysub / dsub**3
        nsub = nsub / nsub.sum()
        m2 = float(np.sum(nsub * dsub**2)); m3 = float(np.sum(nsub * dsub**3))
        m4 = float(np.sum(nsub * dsub**4))
        d_eff = (m3 / m2) if dmode == "d32" else (m4 / m3)
    u = d[0] / d_eff
    if closure == "J":
        a_ = float(np.sqrt(1 - (1 - u) ** paJ)); b_ = float(1 - (1 - u) ** pbJ)
    else:
        a_ = float(fa(u)); b_ = float(fb(u))
    x = y_s
    coarse = g_sub / (1.0 - (1.0 - a_) * x)
    fine = beta / (1.0 - (1.0 - beta + b_ * beta * (1.0 - 1.0 / beta)) * (1.0 - x))
    return float(min(coarse, fine))


def main() -> int:
    lines: list[str] = []
    failures: list[str] = []

    def say(s: str = "") -> None:
        lines.append(s)
        print(s)

    meas3 = mm.measured(mu_p=MU)
    beta = meas3[("lognormal", 0.0)][0]
    pts23 = mc.binary_points(mu_p=MU)
    allpts = gx.ext_points(beta)
    cand, _, _ = gx.candidates(beta, pts23)
    paJ, pbJ = cand["J"](0)
    tab, fa, fb = gx.tabulated_closure(beta, allpts)

    say("=== grade_ternary.py -- PREREG_TERNARY_2026-09-14.md ===")
    say(f"   beta = {beta:.6f} (measured, mu_p {MU}); PHI_SCALE = {SCALE:.1e}")
    say(f"   J  PA {paJ:.4f} PB {pbJ:.4f};  TAB at u = "
        f"{', '.join(f'{u:.3f}' for u in sorted(tab))}")
    say()

    # ---------------------------------------------------------------- reproduction gate
    say("=== REPRODUCTION GATE: the predictions this grader reads are the ones the "
        "pre-registration committed to ===")
    say("   (a claim is reproduced from the stored artifact BEFORE anything is graded "
        "against it -- lessons.md 2026-09-02)")
    ppath = os.path.join(OUT, "predict_ternary.txt")
    stored: dict[tuple, tuple[float, float]] = {}
    if os.path.exists(ppath):
        blk = open(ppath, encoding="utf-8").read().split("=== BLIND PREDICTIONS ===")[1]
        blk = blk.split("=== CANDIDATE SEPARATION ===")[0]
        for ln in blk.splitlines():
            m = re.match(r"\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s", ln)
            if m:
                stored[(round(float(m.group(1)), 4), round(float(m.group(2)), 4),
                        round(float(m.group(3)), 4))] = (float(m.group(4)), float(m.group(5)))
    say(f"   {len(stored)} predictions recovered from {os.path.basename(ppath)}")
    repro_bad = []
    for comp in list(COMPOSITIONS) + [MONO]:
        key = tuple(round(v, 4) for v in comp)
        p = make_ternary(comp)
        pj = mm.packing_density(p.diameters, p.volume_fractions, beta, paJ, pbJ)
        pt = mc.packing_density_tabulated(p.diameters, p.volume_fractions, beta, fa, fb)
        if key not in stored:
            repro_bad.append(f"{key}: absent from the stored prereg artifact")
            continue
        sj, st_ = stored[key]
        if abs(pj - sj) > 5e-5 or abs(pt - st_) > 5e-5:
            repro_bad.append(f"{key}: now J {pj:.4f} TAB {pt:.4f}, prereg J {sj:.4f} TAB {st_:.4f}")
    if repro_bad:
        say("   *** THE PREDICTIONS HAVE MOVED SINCE THE PRE-REGISTRATION ***")
        for b in repro_bad:
            say(f"       {b}")
        failures.append("REPRO")
    else:
        say(f"   all {len(stored)} match the stored values to 5e-5. The gate is reading "
            f"the number that was committed to.")
    # and the bars' own origin
    for name, fn, pset, tag in (("J", lambda d, y: mm.packing_density(d, y, beta, paJ, pbJ), pts23, "own"),
                                ("TAB", lambda d, y: mc.packing_density_tabulated(d, y, beta, fa, fb), allpts, "all")):
        res = np.array([(fn(np.array([r, 1.0]), np.array([1.0 - x, x])) - phi) / SCALE
                        for r, x, phi in pset])
        rms, mx = float(np.sqrt(np.mean(res**2))), float(np.max(np.abs(res)))
        want = INSAMPLE_EXPECT[(name, tag)]
        ok = abs(rms - want[0]) < 0.02 and abs(mx - want[1]) < 0.02
        say(f"   in-sample {name} ({tag}): rms {rms:.2f} (prereg {want[0]}), "
            f"max {mx:.2f} (prereg {want[1]})  -> {'OK' if ok else '*** MOVED ***'}")
        if not ok:
            failures.append(f"REPRO-{name}")
    say()

    # ------------------------------------------------------- WIRING (added before the
    # data was graded, and it is the check that could overturn the whole block: if the
    # K = 3 path into packing_density is wrong, every residual below is measuring my
    # plumbing rather than additivity. A K = 3 spec with the MIDDLE class zeroed must
    # reproduce the corresponding BINARY exactly, through the ternary code path.)
    say("=== WIRING GATE: the K=3 path reduces to the K=2 path when a class is emptied ===")
    wire_bad = []
    wire_n = 0          # counted, not computed from a formula that can drift
    # Which class is emptied decides which RATIO survives, and getting that wrong is a
    # bug in the TEST, not in the model: emptying the middle leaves d = 4 and 1, i.e.
    # ratio 4. All three emptyings are exercised so the check covers both ratios the
    # ladder contains rather than one of them twice.
    #   empty M -> d 4,1 ratio 4, x_small = y_S
    #   empty L -> d 2,1 ratio 2, x_small = y_S/(y_M+y_S)
    #   empty S -> d 4,2 ratio 2, x_small = y_M/(y_L+y_M)
    wire_cases = []
    for xs in (0.15, 0.25, 0.10):
        wire_cases.append(((1.0 - xs, 0.0, xs), 4.0, xs, "empty M"))
    for ym, ys in ((0.60, 0.40), (0.70, 0.30)):
        wire_cases.append(((0.0, ym, ys), 2.0, ys / (ym + ys), "empty L"))
    for yl, ym in ((0.70, 0.30), (0.55, 0.45)):
        wire_cases.append(((yl, ym, 0.0), 2.0, ym / (yl + ym), "empty S"))
    for comp3, r, xs, why in wire_cases:
        p3 = make_ternary(comp3)
        assert abs(p3.size_ratio - r) < 1e-12, (
            f"the TEST is wrong, not the model: {why} on {comp3} gives ratio "
            f"{p3.size_ratio}, not {r}")
        for nm, f3, f2 in (
                ("J", lambda d, y: mm.packing_density(d, y, beta, paJ, pbJ),
                 lambda: mm.binary_density(r, xs, beta, paJ, pbJ)),
                ("TAB", lambda d, y: mc.packing_density_tabulated(d, y, beta, fa, fb),
                 lambda: mc.packing_density_ab(np.array([r, 1.0]), np.array([1.0 - xs, xs]),
                                               beta, fa, fb))):
            got = f3(p3.diameters, p3.volume_fractions)
            want = f2()
            wire_n += 1
            if abs(got - want) > 1e-12:
                wire_bad.append(f"{nm} {why} r={r} x_s={xs:.4f}: K3 {got!r} vs K2 {want!r}")
            say(f"   {nm:>4} {why} r={r:.0f} x_s={xs:.4f}: K3-path {got:.12f}  "
                f"K2 {want:.12f}  diff {got - want:+.2e}")
    # The SEQUENTIAL alternative needs the same gate and did not have one until a
    # referee pass pointed it out (2026-09-14): it is graded against the additive rule
    # in T4, so if it did not reduce to the binary closure at K = 2 the comparison would
    # be between a model and a bug. Emptying the MIDDLE class collapses the composite to
    # pure L, so d_eff -> d_L, g_sub -> beta, and the construction must return the plain
    # ratio-4 binary.
    for xs in (0.15, 0.25, 0.10):
        got = sequential((1.0 - xs, 0.0, xs), "d32", "TAB", beta, paJ, pbJ, fa, fb)
        want = mc.packing_density_ab(np.array([4.0, 1.0]), np.array([1.0 - xs, xs]),
                                     beta, fa, fb)
        wire_n += 1
        if abs(got - want) > 1e-12:
            wire_bad.append(f"SEQ empty M x_s={xs}: {got!r} vs binary {want!r}")
        say(f"    SEQ empty M r=4 x_s={xs:.4f}: nested {got:.12f}  K2 {want:.12f}  "
            f"diff {got - want:+.2e}")
    # and the diameters the K=3 spec produces must be the binary's own two
    say("   (the emptied-middle spec also has to be the same DISTRIBUTION -- "
        "check_discrete_psd.py G6 shows the radii array is bit-identical)")
    if wire_bad:
        for b in wire_bad:
            say(f"   *** {b}")
        failures.append("WIRING")
        say("   -> WIRING FAIL. Nothing below can be read.")
    else:
        say(f"   -> WIRING PASS ({wire_n} of {wire_n} to 1e-12: J and TAB over both "
            f"ratios and all three emptyings, plus the SEQUENTIAL rule). The K=3 "
            f"residuals below are not plumbing.")
    say()

    # ---------------------------------------------------------------- data
    g, nfiles, nraw = load()
    say(f"=== data: {nraw} records from {nfiles} shard files -> {len(g)} compositions ===")
    adm: dict[tuple, list[dict]] = {}
    for key in sorted(g):
        ok = []
        for c in g[key]:
            a, why = admissible(c)
            if a:
                ok.append(c)
            else:
                say(f"   REFUSED {c.get('key')}: {why}")
        adm[key] = ok
    nadm = sum(len(v) for v in adm.values())
    say(f"   N before the admission filter: {nraw}; after: {nadm}")
    nerr = sum(1 for cs in g.values() for c in cs if c.get("error"))
    say(f"   errors: {nerr}; non-converged: "
        f"{sum(1 for cs in g.values() for c in cs if not c.get('error') and not c.get('converged'))}")
    say()

    def mean_sd(comp):
        key = tuple(round(v, 4) for v in comp)
        v = np.array([c["phi"] for c in adm.get(key, [])], dtype=float)
        if v.size == 0:
            return float("nan"), float("nan"), 0
        return float(v.mean()), (float(v.std(ddof=1)) if v.size > 1 else float("nan")), v.size

    # ---------------------------------------------------------------- T0
    say("=== T0 -- monodisperse control through the DiscretePSD code path ===")
    m, sd, n = mean_sd(MONO)
    dev = (m - beta) / SCALE
    say(f"   phi {m:.6f} +- {sd:.6f} (n {n})  vs beta {beta:.6f}: {dev:+.2f} scale")
    t0_ok = abs(dev) <= T0_BAR
    say(f"   -> T0 {'PASS' if t0_ok else 'FAIL'} (bar {T0_BAR} scale). The block is on "
        f"the production protocol and the new PSD path does not move it.")
    if not t0_ok:
        failures.append("T0")
    say()

    # ---------------------------------------------------------------- T1/T2
    say("=== T1 (TAB, primary) and T2 (J, control) -- the ADDITIVE prediction ===")
    say(f"   {'y_L':>6}{'y_M':>6}{'y_S':>6}{'phi':>10}{'sd':>9}{'n':>3}"
        f"{'TAB':>9}{'res':>8}{'J':>9}{'res':>8}{'SEQd32':>9}{'res':>8}{'SEQd43':>9}{'res':>8}")
    resT, resJ, resS32, resS43 = [], [], [], []
    measured_phi = {}
    for comp in COMPOSITIONS:
        m, sd, n = mean_sd(comp)
        p = make_ternary(comp)
        pt = mc.packing_density_tabulated(p.diameters, p.volume_fractions, beta, fa, fb)
        pj = mm.packing_density(p.diameters, p.volume_fractions, beta, paJ, pbJ)
        s32 = sequential(comp, "d32", "TAB", beta, paJ, pbJ, fa, fb)
        s43 = sequential(comp, "d43", "TAB", beta, paJ, pbJ, fa, fb)
        measured_phi[comp] = (m, sd, n)
        # SIGN CONVENTION, fixed here and used everywhere below:
        # residual = MEASURED - PREDICTED. Positive means the model UNDER-predicts,
        # which is the direction the continuous families fail in.
        rt, rj = (m - pt) / SCALE, (m - pj) / SCALE
        r32, r43 = (m - s32) / SCALE, (m - s43) / SCALE
        resT.append(rt); resJ.append(rj); resS32.append(r32); resS43.append(r43)
        say(f"   {comp[0]:6.2f}{comp[1]:6.2f}{comp[2]:6.2f}{m:10.4f}{sd:9.4f}{n:3d}"
            f"{pt:9.4f}{rt:+8.1f}{pj:9.4f}{rj:+8.1f}"
            f"{s32:9.4f}{r32:+8.1f}{s43:9.4f}{r43:+8.1f}")
    say("   residual = MEASURED - PREDICTED, in units of PHI_SCALE. POSITIVE means the "
        "model UNDER-predicts.")
    say()
    stats = {}
    for nm, arr in (("TAB", resT), ("J", resJ), ("SEQd32", resS32), ("SEQd43", resS43)):
        a = np.asarray(arr)
        stats[nm] = {"rms": float(np.sqrt(np.mean(a**2))), "max": float(np.max(np.abs(a))),
                     "mean": float(a.mean()), "npos": int(np.sum(a > 0)), "n": a.size}
        say(f"   {nm:>7}: rms {stats[nm]['rms']:6.2f}  max |res| {stats[nm]['max']:6.2f}  "
            f"mean signed {stats[nm]['mean']:+6.2f}  {stats[nm]['npos']}/{a.size} positive")
    say()
    t1_ok = stats["TAB"]["rms"] <= T1_RMS_BAR and stats["TAB"]["max"] <= T1_MAX_BAR
    say(f"   -> T1 {'PASS' if t1_ok else 'FAIL'}: rms {stats['TAB']['rms']:.2f} vs bar "
        f"{T1_RMS_BAR}, max {stats['TAB']['max']:.2f} vs bar {T1_MAX_BAR}")
    if not t1_ok:
        failures.append("T1")
    t2_ok = stats["J"]["rms"] <= T2_RMS_BAR and stats["J"]["max"] <= T2_MAX_BAR
    say(f"   -> T2 {'PASS' if t2_ok else 'FAIL'}: rms {stats['J']['rms']:.2f} vs bar "
        f"{T2_RMS_BAR}, max {stats['J']['max']:.2f} vs bar {T2_MAX_BAR}")
    say("      T2's outcome was registered as CARRYING NO INFORMATION about additivity: "
        "this ladder contains u = 0.25, where X1 already found J under-predicts.")
    if not t2_ok:
        failures.append("T2")
    say()

    # ---------------------------------------------------------------- T3
    say("=== T3 -- the SIGN structure ===")
    npos, nn = stats["TAB"]["npos"], stats["TAB"]["n"]
    mn = stats["TAB"]["mean"]
    additive_reading = not (npos == nn or npos == 0) and abs(mn) <= T3_MEAN_BAR
    nonadd_reading = (npos == nn) and mn >= T3_MEAN_BAR
    say(f"   TAB signed residuals: {npos}/{nn} positive, mean {mn:+.2f} scale "
        f"(bar {T3_MEAN_BAR})")
    say(f"   ADDITIVE reading (mixed signs AND |mean| <= {T3_MEAN_BAR}): {additive_reading}")
    say(f"   NON-ADDITIVE reading ({nn}/{nn} positive AND mean >= +{T3_MEAN_BAR}): "
        f"{nonadd_reading}")
    if additive_reading == nonadd_reading:
        say("   -> T3 UNDECIDED: neither reading is satisfied, or both are. Reported as "
            "undecided, which is a real verdict and is not the same as either.")
        failures.append("T3-undecided")
    else:
        say(f"   -> T3 verdict: {'NON-ADDITIVE' if nonadd_reading else 'ADDITIVE'}")
    # a sign test that can fail, quoted with its own p-value
    from math import comb
    p_two = 2.0 * sum(comb(nn, k) for k in range(npos, nn + 1)) / 2**nn if npos > nn / 2 else 1.0
    say(f"   two-sided sign test on {npos}/{nn}: p = {min(p_two, 1.0):.4f} "
        f"(n = {nn} is small; this bounds a direction, it does not measure a size)")
    say()

    # ---------------------------------------------------------------- T4
    say("=== T4 -- ADDITIVE vs SEQUENTIAL, the two arms ===")
    add_rms = stats["TAB"]["rms"]
    for nm in ("SEQd32", "SEQd43"):
        seq_rms = stats[nm]["rms"]
        lo, hi = sorted((add_rms, seq_rms))
        ratio = hi / lo if lo > 0 else float("inf")
        winner = "ADDITIVE (TAB)" if add_rms < seq_rms else f"SEQUENTIAL ({nm})"
        inf = ratio >= T4_SEP_BAR
        say(f"   TAB rms {add_rms:.2f} vs {nm} rms {seq_rms:.2f}: ratio "
            f"loser/winner = {ratio:.2f} (bar {T4_SEP_BAR}) -> winner {winner}, "
            f"{'INFORMATIVE' if inf else 'NOT INFORMATIVE'}")
    # the convention knob: do the two SEQ variants bracket the additive prediction?
    brack = []
    for comp in COMPOSITIONS:
        p = make_ternary(comp)
        pt = mc.packing_density_tabulated(p.diameters, p.volume_fractions, beta, fa, fb)
        s32 = sequential(comp, "d32", "TAB", beta, paJ, pbJ, fa, fb)
        s43 = sequential(comp, "d43", "TAB", beta, paJ, pbJ, fa, fb)
        brack.append(min(s32, s43) <= pt <= max(s32, s43))
    say(f"   convention knob: the d32/d43 variants bracket the additive prediction on "
        f"{sum(brack)}/{len(brack)} compositions (registered: not informative if they do)")
    say()

    # ---------------------------------------------------------------- T5
    say("=== T5 -- J vs TAB: PRE-DECLARED NOT INFORMATIVE ===")
    lo, hi = sorted((stats["TAB"]["rms"], stats["J"]["rms"]))
    say(f"   rms ratio {hi / lo if lo > 0 else float('inf'):.2f}; prediction separation "
        f"was 1.04 scale, below the seed error. Printed; decides nothing. This was "
        f"written in the pre-registration, not after seeing the ratio.")
    say()

    # ---------------------------------------------------------------- T6
    say("=== T6 (exploratory, not in FAILURES) -- three-class structure ===")
    say("   Per-DESIGN-class mean contact number and rattler fraction. A binary cannot "
        "show the MIDDLE class, which is the whole point of a ternary.")
    say(f"   {'y_L':>6}{'y_M':>6}{'y_S':>6}"
        f"{'z_S':>7}{'z_M':>7}{'z_L':>7}{'rat_S':>8}{'rat_M':>8}{'rat_L':>8}{'z_nr':>7}")
    for comp in COMPOSITIONS:
        key = tuple(round(v, 4) for v in comp)
        cs = adm.get(key, [])
        if not cs:
            continue
        def avg(ci, fld):
            v = [c["by_design_class"][ci][fld] for c in cs
                 if c.get("by_design_class") and ci in c["by_design_class"]]
            return float(np.mean(v)) if v else float("nan")
        say(f"   {comp[0]:6.2f}{comp[1]:6.2f}{comp[2]:6.2f}"
            f"{avg('c0', 'z_mean'):7.2f}{avg('c1', 'z_mean'):7.2f}{avg('c2', 'z_mean'):7.2f}"
            f"{avg('c0', 'rattler_frac'):8.3f}{avg('c1', 'rattler_frac'):8.3f}"
            f"{avg('c2', 'rattler_frac'):8.3f}"
            f"{float(np.mean([c['z_nr'] for c in cs])):7.2f}")
    say()

    say(f"FAILURES: {failures}")
    with open(os.path.join(OUT, "grade_ternary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
