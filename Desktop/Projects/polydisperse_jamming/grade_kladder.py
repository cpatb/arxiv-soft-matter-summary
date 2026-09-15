"""grade_kladder.py -- the verdicts for PREREG_KLADDER_2026-09-15.md.

Reads `phase12_kl_shard*.jsonl` (N = 1500) and `phase12_kl3k_shard*.jsonl` (N = 3000).
The globs match NEITHER `phase12_tern_shard*` NOR `phase12_mixr*` NOR `phase12_extr*`,
so T0-T6, NA1-NA4, M0-M6, F1-F3 and X0-X5 are CLOSED against this grader; a structural
guard below ASSERTS that rather than trusting the string.

A FAILURE IS THE MODEL'S RANGE, NOT A REFIT. No bar in this file may be moved. The
grid, the models and the admission rule come from `kladder.py`; nothing is rebuilt here.

Run it and READ THE EXIT CODE. lessons.md 2026-09-14: a new gate made the last grader
exit rc = 1 and a stale out/*.txt still read PASS four greps running.
"""

from __future__ import annotations

import collections
import csv
import glob
import json
import math
import os
import re
import sys

import numpy as np
from scipy import stats

import grade_mixture_ext as gx
import grade_ternary as gt
import kladder as kl
import mixture_closure as mc
import mixture_model as mm

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = kl.PHI_SCALE
MU = kl.MU_P
GLOBS = {1500: "phase12_kl_shard*.jsonl", 3000: "phase12_kl3k_shard*.jsonl"}

BAR_RMS = 2.02           # TAB's own in-sample rms on the binaries. NOT MOVABLE.
BAR_MAX = 5.74
SEP_BAR = 1.6
BAR_MONO = 2.0           # G0, the protocol control, in scale units
BAR_SPEARMAN = 0.70      # G2

# Structural guard: this grader must be incapable of loading a closed block's records.
for _g in GLOBS.values():
    assert not any(_g.startswith(p) for p in ("phase12_tern", "phase12_mixr",
                                              "phase12_extr")), _g
    assert "kl" in _g


def load(n: int):
    """Records for one N, grouped by cell label. Prints N before and after the filter."""
    g: dict[str, list[dict]] = collections.defaultdict(list)
    files = sorted(glob.glob(os.path.join(OUT, GLOBS[n])))
    nraw = 0
    for f in files:
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue            # a truncated last line from a killed shard
                nraw += 1
                g[str(r.get("kl_label"))].append(r)
    return g, len(files), nraw


def admissible(c: dict, n: int) -> tuple[bool, str]:
    """The registered admission rule, applied in code (prereg section 5)."""
    if c.get("error"):
        return False, "error"
    if not c.get("converged"):
        return False, "not converged"
    if c.get("n_requested") != n:
        return False, f"N {c.get('n_requested')}"
    if abs(float(c.get("mu_p", -1)) - MU) > 1e-9:
        return False, f"mu {c.get('mu_p')}"
    counts = c.get("n_design_counts") or []
    if len(counts) > 1 and min(counts) < kl.MIN_PHASE_COUNT:
        return False, f"phase count {counts}"
    return True, ""


def spearman(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3:
        return float("nan"), float("nan")
    r = stats.spearmanr(a[m], b[m])
    return float(r.statistic), float(r.pvalue)


def main() -> int:
    lines: list[str] = []
    failures: list[str] = []
    notinfo: list[str] = []

    def say(s: str = "") -> None:
        lines.append(s)
        print(s, flush=True)

    meas = mm.measured(mu_p=MU)
    beta1500 = meas[("lognormal", 0.0)][0]
    pts23 = mc.binary_points(mu_p=MU)
    allpts = gx.ext_points(beta1500)
    cand, _, _ = gx.candidates(beta1500, pts23)
    paJ, pbJ = cand["J"](0)
    _tab, fa, fb = gx.tabulated_closure(beta1500, allpts)

    def aJ(u):
        return float(np.sqrt(1 - (1 - u) ** paJ))

    def bJ(u):
        return float(1 - (1 - u) ** pbJ)

    def models(d, y, beta):
        out = {"ADDtab": float(mc.packing_density_tabulated(
                   np.asarray(d, float), np.asarray(y, float), beta, fa, fb)),
               "ADDj": float(mm.packing_density(
                   np.asarray(d, float), np.asarray(y, float), beta, paJ, pbJ)),
               "IDEAL": kl.ideal_ceiling(d, y, beta)}
        for dm in ("d10", "d21", "d32", "d43"):
            out[f"SEQ{dm}"] = kl.sequential_k(d, y, beta, fa, fb, dm)
        return out

    say("=== grade_kladder.py -- PREREG_KLADDER_2026-09-15.md ===")
    say(f"   beta(N=1500) = {beta1500:.6f} (measured, mu_p {MU}); "
        f"PHI_SCALE = {SCALE:.1e}; bars rms {BAR_RMS} max {BAR_MAX} (NOT MOVABLE)")
    say()

    # ------------------------------------------------- reproduction of the predictions
    say("=== REPRODUCTION GATE: the predictions graded here are the ones committed to ===")
    ppath = os.path.join(OUT, "predict_kladder.txt")
    stored: dict[str, dict[str, float]] = {}
    if os.path.exists(ppath):
        txt = open(ppath, encoding="utf-8").read()
        blk = txt.split("=== BLIND PREDICTIONS ===")[1].split("=== END BLIND")[0]
        for ln in blk.splitlines():
            m = re.match(r"\s*(\S+)\s+(\d+)\s+(nan|[\d.]+)((?:\s+[-\d.]+){9})\s*$", ln)
            if m:
                v = [float(x) for x in m.group(4).split()]
                stored[m.group(1)] = dict(zip(
                    ("ADDtab", "ADDj", "SEQd10", "SEQd21", "SEQd32", "SEQd43",
                     "IDEAL", "head", "sep"), v))
    say(f"   {len(stored)} prediction rows recovered from {os.path.basename(ppath)}")
    cells_all = {c["label"]: c for c in kl.build_cells(kl.N_BIG)[0]}
    bad = []
    for lbl, c in cells_all.items():
        if lbl not in stored:
            bad.append(f"{lbl} missing from the prediction file")
            continue
        now = models(c["d_rel"], c["y"], beta1500)
        for k in ("ADDtab", "ADDj", "SEQd32", "SEQd43", "IDEAL"):
            if abs(now[k] - stored[lbl][k]) > 5e-5:
                bad.append(f"{lbl}.{k} {now[k]:.6f} vs stored {stored[lbl][k]:.6f}")
    if bad:
        failures.append(f"REPRODUCTION: {bad}")
        say(f"   🔴 {len(bad)} mismatches: {bad[:6]}")
    else:
        say(f"   all {len(cells_all)} cells reproduce to 5e-5 on ADDtab/ADDj/SEQ/IDEAL.")
        say("   The grader is reading the number the pre-registration committed to.")
    say()

    # --------------------------------------------------------------------- the data
    data: dict[int, dict] = {}
    for n in (1500, 3000):
        g, nfiles, nraw = load(n)
        if nraw == 0:
            say(f"=== N = {n}: NO RECORDS ON DISK -- every gate that needs it is "
                f"reported UNRUN, not passed ===")
            say()
            continue
        nerr = sum(1 for cs in g.values() for c in cs if c.get("error"))
        nnc = sum(1 for cs in g.values() for c in cs
                  if not c.get("error") and not c.get("converged"))
        ok = {lbl: [c for c in cs if admissible(c, n)[0]] for lbl, cs in g.items()}
        nafter = sum(len(v) for v in ok.values())
        say(f"=== data, N = {n}: {nraw} records from {nfiles} shard files -> "
            f"{len(g)} labels ===")
        say(f"   N before the admission filter: {nraw}; after: {nafter}")
        say(f"   errors: {nerr}; non-converged: {nnc}")
        dropped = {lbl: admissible(cs[0], n)[1]
                   for lbl, cs in g.items() if not ok.get(lbl)}
        if dropped:
            say(f"   labels with NO admitted cell: {dropped}")
        data[n] = ok
        say()

    if 1500 not in data:
        say("NO N = 1500 DATA. Nothing can be graded.")
        return 1

    # per-label summary at each N
    summ: dict[int, dict[str, dict]] = {}
    for n, ok in data.items():
        beta = beta1500
        monos = ok.get("mono", [])
        if monos:
            beta = float(np.mean([c["phi"] for c in monos])) if n != 1500 else beta1500
        s: dict[str, dict] = {}
        for lbl, cs in ok.items():
            if not cs:
                continue
            c0 = cells_all[lbl]
            phi = np.array([float(c["phi"]) for c in cs])
            md = models(c0["d_rel"], c0["y"], beta)
            s[lbl] = {"phi": float(phi.mean()), "sd": float(phi.std(ddof=1))
                      if phi.size > 1 else float("nan"), "n": int(phi.size),
                      "se": (float(phi.std(ddof=1) / math.sqrt(phi.size))
                             if phi.size > 1 else float("nan")),
                      "K": c0["K"], "yS": c0["y_S_target"], "arm": c0["arm"],
                      "beta_used": beta, "models": md,
                      "R": (float(phi.mean()) - md["ADDtab"]) / SCALE,
                      "seR": ((float(phi.std(ddof=1)) / math.sqrt(phi.size)) / SCALE
                              if phi.size > 1 else float("nan")),
                      "nmin": min(cs[0]["n_design_counts"]),
                      "dmoL": float(np.mean([c.get("d_max_over_L", np.nan) for c in cs]))}
        summ[n] = s
        say(f"=== per-cell summary, N = {n} (beta used {beta:.6f}) ===")
        say(f"   {'label':>14}{'arm':>6}{'K':>3}{'y_S':>6}{'phi':>9}{'sd':>8}{'n':>3}"
            f"{'ADDtab':>9}{'R':>8}{'seR':>7}{'SEQd32':>9}{'res':>8}{'n_min':>7}"
            f"{'dmax/L':>8}")
        for lbl in sorted(s, key=lambda x: (s[x]["arm"], s[x]["K"] or 0,
                                            s[x]["yS"] or 0)):
            v = s[lbl]
            say(f"   {lbl:>14}{v['arm']:>6}{v['K']:3d}"
                f"{(float('nan') if v['yS'] is None else v['yS']):6.2f}"
                f"{v['phi']:9.4f}{v['sd']:8.4f}{v['n']:3d}{v['models']['ADDtab']:9.4f}"
                f"{v['R']:+8.2f}{v['seR']:7.2f}{v['models']['SEQd32']:9.4f}"
                f"{(v['phi'] - v['models']['SEQd32']) / SCALE:+8.2f}"
                f"{v['nmin']:7d}{v['dmoL']:8.3f}")
        say()

    s15 = summ[1500]

    # ------------------------------------------------------------------------- G0
    say("=== G0 -- protocol control and CROSS-BLOCK reproduction ===")
    if "mono" in s15:
        dm = (s15["mono"]["phi"] - beta1500) / SCALE
        say(f"   mono phi {s15['mono']['phi']:.6f} +- {s15['mono']['sd']:.6f} "
            f"(n {s15['mono']['n']}) vs beta {beta1500:.6f}: {dm:+.2f} scale "
            f"(bar {BAR_MONO})")
        if abs(dm) > BAR_MONO:
            failures.append(f"G0 mono {dm:+.2f} scale")
    tg, _, _ = gt.load()
    tern_by_comp = {k: [c for c in v if not c.get("error")] for k, v in tg.items()}
    pairs = [("mono", (1.0, 0.0, 0.0)), ("NA2_x6", (0.55, 0.35, 0.10))]
    nbit = obit = 0
    for lbl, comp in pairs:
        old = {int(c["seed"]): float(c["phi"])
               for c in tern_by_comp.get(tuple(round(v, 4) for v in comp), [])
               if c.get("n_requested") == 1500}
        new = {}
        for c in data[1500].get(lbl, []):
            new[int(c["seed"])] = float(c["phi"])
        for sd in sorted(set(old) & set(new)):
            nbit += 1
            d = new[sd] - old[sd]
            obit += (d == 0.0)
            say(f"   BIT {lbl:>7} seed {sd}: kl {new[sd]!r}  tern {old[sd]!r}  "
                f"diff {d:+.3e}")
    say(f"   -> bit-identity {obit} of {nbit} EXACT")
    if nbit and obit != nbit:
        failures.append(f"G0 bit-identity {obit}/{nbit}")
    say(f"   G0 {'PASS' if not [f for f in failures if f.startswith('G0')] else 'FAIL'}")
    say()

    # ------------------------------------------------------------------------- G1
    say("=== G1 (PRIMARY) -- the K profile at matched y_S ===")
    say("   R(2) = 0 is an ANCHOR (K = 2 IS TAB's calibration set), R3 is taken from")
    say("   THIS block's own K = 3 cell at the same y_S. Hypotheses: prereg section 4.")
    g1_sub = []
    for n in sorted(summ):
        s = summ[n]
        say(f"   --- N = {n}")
        say(f"   {'y_S':>6}" + "".join(f"{'R(K=' + str(k) + ')':>11}" for k in kl.K_LEVELS)
            + f"{'R5-R3':>9}{'se':>7}   nearest hypothesis (rms over K = 4, 5)")
        for ys in kl.YS_BALANCED:
            lab = {k: f"K{k}_yS{ys:.2f}" for k in kl.K_LEVELS}
            if not all(v in s for v in lab.values()):
                say(f"   {ys:6.2f}   incomplete: {[k for k, v in lab.items() if v not in s]}")
                continue
            R = {k: s[lab[k]]["R"] for k in kl.K_LEVELS}
            se53 = math.hypot(s[lab[5]]["seR"], s[lab[3]]["seR"])
            hyp = kl.hypothesis_predictions(R[3])
            dist = {nm: math.sqrt(np.mean([(hyp[nm][k] - R[k]) ** 2 for k in (4, 5)]))
                    for nm in hyp}
            rank = sorted(dist, key=lambda x: dist[x])
            d53 = R[5] - R[3]
            g1_sub.append((ys, n, d53, se53, rank[0], dist))
            degen = abs(R[3]) < BAR_RMS
            say(f"   {ys:6.2f}" + "".join(f"{R[k]:+11.2f}" for k in kl.K_LEVELS)
                + f"{d53:+9.2f}{se53:7.2f}   "
                + ", ".join(f"{nm} {dist[nm]:.1f}" for nm in rank)
                + ("   <-- R3 INSIDE THE BAR: ranking DEGENERATE" if degen else ""))
            if degen:
                notinfo.append(f"G1 ranking at N={n} y_S={ys}: |R3| = "
                               f"{abs(R[3]):.2f} < {BAR_RMS}")
        say()
        say("   ABSOLUTE INCREMENTS -- reported BESIDE the ranking above, because every")
        say("   hypothesis in the registered family is a MULTIPLE of R3, and a multiple")
        say("   of a quantity inside its own bar of zero predicts nothing (lessons.md")
        say("   2026-09-06: a relative measure whose denominator is shrinking measures")
        say("   the shrinkage). NO BAR IS MOVED; these are extra numbers beside them.")
        say(f"   {'y_S':>6}{'R4-R3':>9}{'se':>7}{'R5-R4':>9}{'se':>7}{'ratio':>8}"
            "   predicted ratio: H_K 1.00, H_P 1.33, H_T 2.00")
        for ys in kl.YS_BALANCED:
            lab = {k: f"K{k}_yS{ys:.2f}" for k in kl.K_LEVELS}
            if not all(v in s for v in lab.values()):
                continue
            R = {k: s[lab[k]]["R"] for k in kl.K_LEVELS}
            i1, i2 = R[4] - R[3], R[5] - R[4]
            s1 = math.hypot(s[lab[4]]["seR"], s[lab[3]]["seR"])
            s2 = math.hypot(s[lab[5]]["seR"], s[lab[4]]["seR"])
            rat = i2 / i1 if abs(i1) > 1e-9 else float("nan")
            say(f"   {ys:6.2f}{i1:+9.2f}{s1:7.2f}{i2:+9.2f}{s2:7.2f}{rat:8.2f}")
        say()
        say("   DECOMPOSITION of each step into the part the MEASUREMENT moved and the")
        say("   part the MODEL moved: R = (phi - ADD)/scale, so dR = dphi - dADD. Printed")
        say("   because a step in a RESIDUAL is not evidence that the packing did")
        say("   anything -- the closure's own K dependence lands in the same number.")
        say(f"   {'y_S':>6}{'step':>9}{'dR':>9}{'dphi':>9}{'dADD':>9}{'model share':>13}")
        for ys in kl.YS_BALANCED:
            lab = {k: f"K{k}_yS{ys:.2f}" for k in kl.K_LEVELS}
            if not all(v in s for v in lab.values()):
                continue
            for ka, kb in ((3, 4), (4, 5)):
                a, b = s[lab[ka]], s[lab[kb]]
                dR = b["R"] - a["R"]
                dphi = (b["phi"] - a["phi"]) / SCALE
                dadd = (b["models"]["ADDtab"] - a["models"]["ADDtab"]) / SCALE
                share = abs(dadd) / max(abs(dphi) + abs(dadd), 1e-12)
                say(f"   {ys:6.2f}{f'K{ka}->K{kb}':>9}{dR:+9.2f}{dphi:+9.2f}"
                    f"{dadd:+9.2f}{share:12.0%}")
        say()
    sub15 = [t for t in g1_sub if t[1] == 1500]
    npass = sum(1 for t in sub15 if t[2] >= BAR_RMS)
    say(f"   SUB-CLAUSE  R(5) - R(3) >= {BAR_RMS} on at least 2 of 3 y_S levels "
        f"(N = 1500): {npass} of {len(sub15)}")
    if len(sub15) == 3:
        if npass >= 2:
            say("   -> the residual GROWS with K, as registered.")
        else:
            say("   -> 🔴 the residual is ALREADY SATURATED at three classes. The")
            say("      registered prediction FAILS; the correction is a constant.")
            failures.append(f"G1 sub-clause {npass}/3 -- SATURATED, prediction refuted")
    winners = [t[4] for t in sub15]
    say(f"   nearest hypothesis per y_S level (N = 1500): {winners}")
    if winners and len(set(winners)) == 1 and winners[0] in ("H_K", "H_P"):
        say("   -> H_K / H_P. REGISTERED IN ADVANCE: these two are NOT separable by this")
        say("      design (1.25x at K = 4, 1.50x at K = 5, both below 1.6). The verdict")
        say("      is 'linear in K OR in pairs', and the block may not choose between them.")
        notinfo.append("G1 H_K vs H_P: not separable, registered in advance")
    say()

    # ------------------------------------------------------------------------- G2
    say("=== G2 -- the y_S axis at fixed K (five levels, 6.0x, vs the closed block's two) ===")
    for n in sorted(summ):
        for k in (2, 3):
            s = summ[n]
            pts = [(ys, s[f"K{k}_yS{ys:.2f}"]["R"]) for ys in kl.YS_LEVELS
                   if f"K{k}_yS{ys:.2f}" in s]
            if len(pts) < 3:
                continue
            rho, p = spearman([a for a, _ in pts], [b for _, b in pts])
            span = pts[-1][1] - pts[0][1]
            verdict = ("PASS" if (k == 3 and rho >= BAR_SPEARMAN and span >= BAR_RMS)
                       else ("FAIL" if k == 3 else "(reported)"))
            say(f"   N {n} K={k}: spearman(R, y_S) = {rho:+.4f} (p {p:.4f}) over "
                f"{len(pts)} levels; R(0.18) - R(0.03) = {span:+.2f} "
                f"(bar {BAR_RMS}) -> {verdict}")
            say("      " + "  ".join(f"y_S {a:.2f}: {b:+.2f}" for a, b in pts))
            if k == 3 and n == 1500 and verdict == "FAIL":
                failures.append(f"G2 K=3 N=1500 rho {rho:+.3f} span {span:+.2f}")
    say()

    # ------------------------------------------------------------------------- G3
    say("=== G3 -- the SPLIT control: K 3 -> 4 and P 3 -> 6 with the DISTRIBUTION held ===")
    for n in sorted(summ):
        s = summ[n]
        for ys in (0.05, 0.08):
            a, b = f"SPLIT_yS{ys:.2f}", f"K3_yS{ys:.2f}"
            if a not in s or b not in s:
                continue
            d = s[a]["R"] - s[b]["R"]
            se = math.hypot(s[a]["seR"], s[b]["seR"])
            dmod = (s[a]["models"]["ADDtab"] - s[b]["models"]["ADDtab"]) / SCALE
            say(f"   N {n} y_S {ys:.2f}: R(split) {s[a]['R']:+.2f} - R(K3) "
                f"{s[b]['R']:+.2f} = {d:+.2f} (se {se:.2f}, bar {BAR_RMS}); the "
                f"additive prediction itself moved {dmod:+.2f}")
            if n == 1500 and abs(d) >= BAR_RMS:
                failures.append(f"G3 y_S {ys}: split moved R by {d:+.2f}")
    say("   registered: |delta| < 2.02 on 2 of 2 means the residual is a FUNCTIONAL OF")
    say("   THE DISTRIBUTION; a jump comparable to G1's K 3 -> 4 increment means it is")
    say("   COUNTING CLASSES, i.e. an artefact of how a PSD is discretised.")
    say()

    # ------------------------------------------------------------------------- G4
    say("=== G4 -- N-doubling on the IDENTICAL cell list ===")
    if 3000 not in summ:
        say("   UNRUN: no N = 3000 records. Not passed, not failed.")
    else:
        say(f"   {'y_S':>6}{'D1500':>9}{'D3000':>9}{'|diff|':>9}{'2*se':>8}   verdict "
            "   (D = R(5) - R(3))")
        for ys in kl.YS_BALANCED:
            try:
                d15 = s15[f"K5_yS{ys:.2f}"]["R"] - s15[f"K3_yS{ys:.2f}"]["R"]
                s30 = summ[3000]
                d30 = s30[f"K5_yS{ys:.2f}"]["R"] - s30[f"K3_yS{ys:.2f}"]["R"]
                se = math.hypot(s15[f"K5_yS{ys:.2f}"]["seR"], s15[f"K3_yS{ys:.2f}"]["seR"],
                                s30[f"K5_yS{ys:.2f}"]["seR"], s30[f"K3_yS{ys:.2f}"]["seR"])
            except KeyError:
                continue
            okk = abs(d30 - d15) <= 2 * se
            say(f"   {ys:6.2f}{d15:+9.2f}{d30:+9.2f}{abs(d30 - d15):9.2f}{2 * se:8.2f}"
                f"   {'PASS' if okk else 'FAIL'}")
            if not okk:
                failures.append(f"G4 y_S {ys}: increment moved {d30 - d15:+.2f} > 2se")
        # the finite-size correlations, pooled over BOTH N, core arm only
        rows = [(v["R"], v["nmin"], v["dmoL"]) for n in summ for lbl, v in summ[n].items()
                if v["arm"] == "K"]
        rho_n, p_n = spearman([r[1] for r in rows], [r[0] for r in rows])
        rho_d, p_d = spearman([r[2] for r in rows], [r[0] for r in rows])
        say(f"   pooled over both N ({len(rows)} points): spearman(R, n_min) = "
            f"{rho_n:+.4f} (p {p_n:.4f}); spearman(R, d_max/L) = {rho_d:+.4f} "
            f"(p {p_d:.4f}); bar 0.6")
    say()

    # ------------------------------------------------------------------------- G5
    say("=== G5 -- the NA2 exception at six seeds ===")
    if "NA2_x6" in s15 and 3000 in summ and "NA2_x6" in summ[3000]:
        a = [float(c["phi"]) for c in data[1500]["NA2_x6"]]
        b = [float(c["phi"]) for c in data[3000]["NA2_x6"]]
        da = (np.mean(b) - np.mean(a)) / SCALE
        sed = math.hypot(np.std(a, ddof=1) / math.sqrt(len(a)),
                         np.std(b, ddof=1) / math.sqrt(len(b))) / SCALE
        t = stats.ttest_ind(b, a, equal_var=False)
        say(f"   phi(1500) {np.mean(a):.6f} (n {len(a)}, sd {np.std(a, ddof=1):.6f}), "
            f"phi(3000) {np.mean(b):.6f} (n {len(b)}, sd {np.std(b, ddof=1):.6f})")
        say(f"   shift {da:+.2f} scale, se OF THE DIFFERENCE {sed:.2f} -> "
            f"{abs(da) / sed if sed else float('nan'):.2f} se; Welch t = "
            f"{float(t.statistic):+.3f}, p = {float(t.pvalue):.4f}, df ~ "
            f"{float(t.df):.1f}")
        say("   the closed block measured +1.81 on 3 seeds against an se of 0.48.")
        say(f"   -> {'below 2 se: seed noise, as registered' if abs(da) < 2 * sed else 'SURVIVES at 6 seeds: a real finding about this composition'}")
    else:
        say("   UNRUN: needs both N. Not passed, not failed.")
    say()

    # ------------------------------------------------------------------------- G6
    say("=== G6 -- the SEQUENTIAL rule, convention SWEPT not chosen ===")
    core = [lbl for lbl, v in s15.items() if v["arm"] in ("K", "SPLIT", "NA2")
            and v["K"] > 1]
    say(f"   over {len(core)} admitted multi-class cells at N = 1500")
    say(f"   {'model':>10}{'rms':>8}{'max|res|':>10}{'mean':>8}{'above':>7}{'below':>7}")
    rmss = {}
    for nm in ("ADDtab", "ADDj", "SEQd10", "SEQd21", "SEQd32", "SEQd43", "IDEAL"):
        res = np.array([(s15[l]["phi"] - s15[l]["models"][nm]) / SCALE for l in core])
        rmss[nm] = float(np.sqrt(np.mean(res**2)))
        say(f"   {nm:>10}{rmss[nm]:8.2f}{float(np.max(np.abs(res))):10.2f}"
            f"{float(res.mean()):+8.2f}{int((res > 0).sum()):7d}"
            f"{int((res < 0).sum()):7d}")
    say("   'above' = the measurement sits ABOVE the model (the model under-fills).")
    seqs = {k: v for k, v in rmss.items() if k.startswith("SEQ")}
    best = min(seqs, key=lambda k: seqs[k])
    ratio = rmss["ADDtab"] / seqs[best] if seqs[best] else float("inf")
    say(f"   best SEQ convention: {best} (rms {seqs[best]:.2f}); ADDtab/{best} = "
        f"{ratio:.2f} (bar {SEP_BAR})")
    say(f"   convention spread among the four SEQ rms: "
        f"{min(seqs.values()):.2f} to {max(seqs.values()):.2f}, a factor "
        f"{max(seqs.values()) / max(min(seqs.values()), 1e-9):.2f}")
    say("   ADDtab vs SEQd10 was pre-declared NOT INFORMATIVE (separation 2.05 < 3.23).")
    resA = np.array([(s15[l]["phi"] - s15[l]["models"]["ADDtab"]) / SCALE for l in core])
    resS = np.array([(s15[l]["phi"] - s15[l]["models"]["SEQd43"]) / SCALE for l in core])
    fa_ = float((resA > 0).mean())
    fs_ = float((resS < 0).mean())
    say(f"   BRACKET: ADDtab below the measurement on {fa_:.0%} of cells, SEQd43 above "
        f"on {fs_:.0%} (registered: both >= 80%)")
    say("   The rms winner is NOT the rule: at K = 3 the winner carried the OPPOSITE")
    say("   bias, of comparable size, to the loser. The supportable claim is a bracket.")
    if not (fa_ >= 0.8 and fs_ >= 0.8):
        failures.append(f"G6 bracket {fa_:.0%}/{fs_:.0%}")
    say()

    # ------------------------------------------------------------------------- G7
    say("=== G7 -- does phi RISE with K at matched y_S? (the closure says it FALLS) ===")
    nrise = ntot = 0
    for ys in kl.YS_BALANCED:
        lab = [f"K{k}_yS{ys:.2f}" for k in (3, 4, 5)]
        if not all(v in s15 for v in lab):
            continue
        phis = [s15[v]["phi"] for v in lab]
        pred = [s15[v]["models"]["ADDtab"] for v in lab]
        ses = [s15[v]["se"] for v in lab]
        d = (phis[-1] - phis[0]) / SCALE
        se = math.hypot(ses[0], ses[-1]) / SCALE
        ntot += 1
        nrise += d > 0
        say(f"   y_S {ys:.2f}: measured phi K=3,4,5 = "
            f"{', '.join(f'{v:.4f}' for v in phis)}  (K5-K3 = {d:+.2f} scale, "
            f"se {se:.2f})")
        say(f"              ADDtab says   "
            f"{', '.join(f'{v:.4f}' for v in pred)}  "
            f"(K5-K3 = {(pred[-1] - pred[0]) / SCALE:+.2f})")
    say(f"   -> phi rises with K on {nrise} of {ntot} levels (registered: >= 2 of 3)")
    if ntot == 3 and nrise < 2:
        failures.append(f"G7 phi rises on only {nrise}/3")
    say()

    # ------------------------------------------------------------------------------
    say("=== SUMMARY ===")
    say(f"   FAILURES: {failures if failures else 'none'}")
    say(f"   NOT INFORMATIVE (registered in advance): {notinfo if notinfo else 'none'}")
    # One row per admitted cell, so a downstream script never has to re-fit the closure
    # (and so cannot fit it slightly differently). Post-hoc diagnostics read this.
    cpath = os.path.join(OUT, "kladder_cells.csv")
    with open(cpath, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["N", "label", "arm", "K", "y_S", "n_seeds", "phi", "sd", "se",
                    "ADDtab", "R", "seR", "SEQd32", "SEQd43", "IDEAL",
                    "n_min", "d_max_over_L", "beta_used"])
        for n in sorted(summ):
            for lbl in sorted(summ[n]):
                v = summ[n][lbl]
                w.writerow([n, lbl, v["arm"], v["K"], v["yS"], v["n"],
                            f"{v['phi']:.8f}", f"{v['sd']:.8f}", f"{v['se']:.8f}",
                            f"{v['models']['ADDtab']:.8f}", f"{v['R']:.6f}",
                            f"{v['seR']:.6f}", f"{v['models']['SEQd32']:.8f}",
                            f"{v['models']['SEQd43']:.8f}",
                            f"{v['models']['IDEAL']:.8f}",
                            v["nmin"], f"{v['dmoL']:.6f}", f"{v['beta_used']:.8f}"])
    say(f"   per-cell table written to {os.path.basename(cpath)}")
    with open(os.path.join(OUT, "grade_kladder.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
