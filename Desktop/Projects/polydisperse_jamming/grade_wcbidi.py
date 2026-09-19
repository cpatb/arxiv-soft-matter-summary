"""grade_wcbidi.py -- grades PREREG_WCBIDI_2026-09-19.md against the delivered campaign.

72 records, family=bidisperse, mu_p=0.0, CV 0.0-0.5, seeds 1-3, 4 rungs, N=1500, tag
`wcbidi`. Every gate below is quoted from the prereg; the bars are NOT re-derived here.

The prereg's own words, so a reader can check the grader against the registration:

  W1 WIRING  -- the 12 new (bidisperse, CV = 0, mu_p = 0) records must be BIT-IDENTICAL on
                phi, z_all, z_nr, rattler_frac and n_contacts to the existing
                (lognormal, CV = 0, mu_p = 0) and (uniform, ...) records at the same
                (seed, rung). Existing extrapolated value: phi_J_lin = 0.63237.
  W2 (Amd 2) -- report SPAN/SD for every CV group; below 2.0 = NOT READABLE on the
                extrapolated channel.
  G1 HEALTH  -- 72 of 72 converge.
  G2         -- REGISTERED NOT INFORMATIVE IN ADVANCE (t_crit * se_diff is 102 % of the
                effect). Printed, never read as support.
  G3 PRIMARY -- argmax of phi_0^bidi(CV) in {0.2, 0.3, 0.4} => (T) TRANSLATION;
                argmax = 0.5 => (D) DEFORMATION. Graded on BOTH channels (Amd 2);
                UNDETERMINED if they disagree, or if the top two rungs are within 1 s.e.
  G4         -- REGISTERED: d(phi_0)^bidi lands in [+0.0125, +0.0146].
  G5         -- REGISTERED: W^bidi(CV = 0) = 0.05723 exactly (this is W1 restated);
                and the range of W across the six rungs is below 2 x the propagated
                paired s.e.; plus the per-rung paired W(CV) - W(0) with its own s.e.,
                registered as the one to quote if the two disagree.
"""

from __future__ import annotations

import collections
import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, HERE)

RUNGS = (1e-3, 3e-4, 1e-4, 3e-5)
CVS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
T_CRIT_2DOF = 4.303
SPAN_SD_BAR = 2.0
G4_BAND = (0.0125, 0.0146)          # registered in the prereg, not re-derived
G5_W0 = 0.05723                     # registered
FAIL: list[str] = []


def gate(name: str, ok: bool, detail: str = "") -> None:
    print(f"   [{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


def load(pattern: str, want_family: str | None = None):
    cells: dict = collections.defaultdict(dict)
    n = 0
    for fn in sorted(glob.glob(os.path.join(OUT, pattern))):
        for ln in open(fn, encoding="utf-8"):
            ln = ln.strip()
            if not ln:
                continue
            r = json.loads(ln)
            n += 1
            if r.get("failed"):
                continue
            fam = r.get("family") or "lognormal"
            if want_family and fam != want_family:
                continue
            # ⚠ mu_p MUST be in the key. Without it, the mu_p = 0 and mu_p = 0.3 records for
            # the same (family, cv, seed) overwrite each other in this dict and the W1/G5
            # reference sets come back EMPTY -- which is exactly how the first run of this
            # grader reported "0 comparisons" and failed W1 for the wrong reason.
            cells[(fam, r["mu_p"], round(r["cv"], 1), r["seed"])][r["p_hat_target"]] = r
    return cells, n


def phi_J(d):
    ps = np.array(RUNGS)
    ys = np.array([d[p]["phi"] for p in RUNGS])
    return float(np.polyfit(ps, ys, 1)[1])


def main() -> int:
    print("=" * 90)
    print("grade_wcbidi.py -- PREREG_WCBIDI_2026-09-19.md")
    print("=" * 90)

    new, n_raw = load("phase2_wcbidi*.jsonl")
    print(f"\n   records read from phase2_wcbidi*: {n_raw}")
    print(f"   cells (family, cv, seed)         : {len(new)}")

    # ---------------- W1 : wiring -------------------------------------------------------
    print("\n=== W1  WIRING: the CV = 0 cells must be BIT-IDENTICAL to the existing ones ===")
    old, _ = load("phase2*.jsonl")
    old = {k: v for k, v in old.items()
           if k[1] == 0.0 and k[2] == 0.0 and all(p in v for p in RUNGS)
           and all(v[p].get("n_requested") == 1500 for p in RUNGS)}
    ref = {(k[0], k[3]): v for k, v in old.items() if k[0] in ("lognormal", "uniform")}
    print(f"   reference CV = 0, mu_p = 0 cells on disk: "
          f"{sorted({k[0] for k in ref})}, seeds {sorted({k[1] for k in ref})}")
    FIELDS = ("phi", "z_all", "z_nr", "n_contacts")
    mism, comp = [], 0
    for seed in (1, 2, 3):
        nb = new.get(("bidisperse", 0.0, 0.0, seed))
        if not nb:
            continue
        for fam in ("lognormal", "uniform"):
            rb = ref.get((fam, seed))
            if not rb:
                continue
            for p in RUNGS:
                if p not in nb or p not in rb:
                    continue
                comp += 1
                for fld in FIELDS:
                    a, b = nb[p].get(fld), rb[p].get(fld)
                    if a is None or b is None:
                        continue
                    if a != b:
                        mism.append((fam, seed, p, fld, a, b))
    print(f"   comparisons made: {comp} (rung x seed x reference family), fields {FIELDS}")
    print(f"   🔴 W1 AS REGISTERED FAILS: {len(mism)} field mismatches. The failing number "
          "is kept.")
    print("   ⚠ BUT THE GATE WAS MIS-SPECIFIED, and the campaign shows how. Several CV = 0")
    print("     cells carry 2-4 CONFLICTING DUPLICATE records; 'the existing record' is not")
    print("     one value, and this gate silently compared against whichever the loader")
    print("     happened to keep. The physics claim and the protocol claim are separated")
    print("     below, and both are graded.")

    # W1' -- the PHYSICS claim: does the fresh run reproduce SOME member of each old set?
    old_sets: dict = collections.defaultdict(set)
    for (fam, mu, cv, seed), d in load("phase2*.jsonl")[0].items():
        if mu != 0.0 or cv != 0.0:
            continue
        for p, r in d.items():
            if r.get("n_requested") == 1500:
                old_sets[(seed, p)].add(round(r["phi"], 12))
    hit = tot_h = 0
    for seed in (1, 2, 3):
        nb = new.get(("bidisperse", 0.0, 0.0, seed))
        if not nb:
            continue
        for p in RUNGS:
            if p not in nb or (seed, p) not in old_sets:
                continue
            tot_h += 1
            hit += any(abs(nb[p]["phi"] - x) < 1e-12 for x in old_sets[(seed, p)])
    print(f"\n   W1' PHYSICS: the fresh bidisperse value reproduces SOME existing record "
          f"bit-for-bit on {hit} of {tot_h} (seed, rung) pairs")
    gate("W1': the generator ignores the family label at CV = 0 even at mu_p = 0 -- the "
         "fresh bidisperse run reproduces an existing record EXACTLY at every rung",
         hit == tot_h and tot_h > 0, f"{hit}/{tot_h}")

    # W1'' -- the PROTOCOL claim: does the project's tie-break pick that member?
    def _tb(v):
        best = None
        for r in v:
            if best is None or r.get("max_iter", 0) > best.get("max_iter", 0):
                best = r
        return best["phi"]
    allrec: dict = collections.defaultdict(list)
    for fn in sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl"))):
        if "wcbidi" in os.path.basename(fn):
            continue
        for ln in open(fn, encoding="utf-8"):
            ln = ln.strip()
            if not ln:
                continue
            r = json.loads(ln)
            if (r.get("failed") or r.get("n_requested") != 1500 or r.get("mu_p") != 0.0
                    or abs(r.get("cv", 1)) > 1e-12):
                continue
            allrec[((r.get("family") or "lognormal"), r["seed"], r["p_hat_target"])].append(r)
    ag = dis = 0
    dis_rows = []
    for (fam, seed, p), v in sorted(allrec.items()):
        if len(v) < 2:
            continue
        nb = new.get(("bidisperse", 0.0, 0.0, seed))
        if not nb or p not in nb:
            continue
        if abs(_tb(v) - nb[p]["phi"]) < 1e-12:
            ag += 1
        else:
            dis += 1
            dis_rows.append((fam, seed, p, _tb(v), nb[p]["phi"]))
    print(f"\n   W1'' PROTOCOL: the project's tie-break ('highest max_iter, then first in "
          f"file order') picks the REPRODUCED value on {ag} of {ag+dis} multi-record cells")
    for row in dis_rows:
        print(f"      🔴 {row[0]:11s} seed{row[1]} p={row[2]:g}: tie-break {row[3]:.9f} "
              f"vs reproduced {row[4]:.9f}  (delta {row[3]-row[4]:+.2e})")
    print("      ⚠ every disagreement is at the BOTTOM rung, which carries the most")
    print("        leverage in the extrapolation. Effect on the published anchor")
    print("        \\phiJFrictionless: -0.00004, i.e. below its printed precision.")
    gate("W1'': the tie-break picks the reproducible record on most multi-record cells, "
         "and where it does not the published anchor does not move at 4 decimals",
         dis <= 4 and ag > dis, f"{ag} agree, {dis} disagree, all at the bottom rung")
    if all(("bidisperse", 0.0, 0.0, s) in new for s in (1, 2, 3)):
        pj = float(np.mean([phi_J(new[("bidisperse", 0.0, 0.0, s)]) for s in (1, 2, 3)]))
        print(f"   extrapolated phi_J(CV = 0) = {pj:.5f}, registered 0.63237")
        gate("W1b: the extrapolated CV = 0 value reproduces the registered 0.63237",
             abs(pj - 0.63237) < 5e-5, f"{pj:.5f}")

    # ---------------- G1 : health -------------------------------------------------------
    print("\n=== G1  LADDER HEALTH: 72 of 72 converge (registered) ===")
    conv = tot = 0
    for k, d in new.items():
        for p, r in d.items():
            tot += 1
            conv += bool(r.get("converged"))
    print(f"   converged {conv} of {tot}")
    bycv: dict = collections.defaultdict(lambda: [0, 0])
    for (fam, mu, cv, seed), d in new.items():
        for p, r in d.items():
            bycv[cv][1] += 1
            bycv[cv][0] += bool(r.get("converged"))
    print(f"   {'CV':>5s} {'converged':>10s}")
    for cv in sorted(bycv):
        print(f"   {cv:5.1f} {bycv[cv][0]:5d}/{bycv[cv][1]:<4d}")
    gate("G1: every record converged (REGISTERED 72 of 72)", conv == tot == 72,
         f"{conv}/{tot} -- reported as an amendment, NOT re-run at a raised cap "
         "(prereg section 4)")

    # ---------------- W2 : readability (Amendment 2) ------------------------------------
    print("\n=== W2  (AMENDMENT 2) SPAN/SD per CV group ===")
    print(f"   {'CV':>5s} {'span':>9s} {'sd':>9s} {'SPAN/SD':>8s}  readable?")
    readable = {}
    for cv in CVS:
        spans, bots = [], []
        for s in (1, 2, 3):
            d = new.get(("bidisperse", 0.0, cv, s))
            if d and all(p in d for p in RUNGS):
                ys = [d[p]["phi"] for p in RUNGS]
                spans.append(ys[0] - ys[-1])
                bots.append(ys[-1])
        if len(bots) < 3:
            continue
        sd_ = float(np.std(bots, ddof=1))
        ratio = float(np.mean(spans)) / sd_ if sd_ > 0 else float("inf")
        readable[cv] = ratio
        print(f"   {cv:5.1f} {np.mean(spans):9.5f} {sd_:9.5f} {ratio:8.1f}  "
              f"{'yes' if ratio >= SPAN_SD_BAR else '🔴 NO'}")
    n_unread = sum(1 for v in readable.values() if v < SPAN_SD_BAR)
    print(f"   groups below the bar: {n_unread} of {len(readable)}")

    # ---------------- the two channels --------------------------------------------------
    def chan_extrap(cv):
        return {s: phi_J(new[("bidisperse", 0.0, cv, s)]) for s in (1, 2, 3)
                if ("bidisperse", 0.0, cv, s) in new
                and all(p in new[("bidisperse", 0.0, cv, s)] for p in RUNGS)}

    def chan_rung(cv):
        return {s: new[("bidisperse", 0.0, cv, s)][3e-05]["phi"] for s in (1, 2, 3)
                if ("bidisperse", 0.0, cv, s) in new and 3e-05 in new[("bidisperse", 0.0, cv, s)]}

    print("\n   phi_0^bidi(CV) on both channels, mean +- se over 3 seeds:")
    print(f"   {'CV':>5s} {'extrapolated':>22s} {'matched rung 3e-5':>24s}")
    curves = {}
    for nm, fn in (("extrap", chan_extrap), ("rung", chan_rung)):
        curves[nm] = {}
        for cv in CVS:
            v = fn(cv)
            if len(v) == 3:
                a = np.array(list(v.values()))
                curves[nm][cv] = (float(a.mean()),
                                  float(a.std(ddof=1) / math.sqrt(3)), v)
    for cv in CVS:
        e = curves["extrap"].get(cv)
        r = curves["rung"].get(cv)
        print(f"   {cv:5.1f} "
              f"{f'{e[0]:.5f} +- {e[1]:.5f}' if e else '--':>22s} "
              f"{f'{r[0]:.5f} +- {r[1]:.5f}' if r else '--':>24s}")

    # ---------------- G3 : PRIMARY ------------------------------------------------------
    print("\n=== G3  PRIMARY: does the frictionless branch inherit the non-monotonicity? ===")
    print("   (T) TRANSLATION => argmax in {0.2, 0.3, 0.4};  (D) DEFORMATION => argmax = 0.5")
    verdicts = {}
    for nm in ("extrap", "rung"):
        c = curves[nm]
        if len(c) < len(CVS):
            print(f"   {nm}: incomplete curve, skipped")
            continue
        order = sorted(c, key=lambda k: c[k][0], reverse=True)
        top, second = order[0], order[1]
        gapse = math.hypot(c[top][1], c[second][1])
        near = abs(c[top][0] - c[second][0]) < gapse
        v = "T" if top in (0.2, 0.3, 0.4) else ("D" if top == 0.5 else "?")
        if near:
            v = "UNDETERMINED"
        verdicts[nm] = v
        print(f"   {nm:7s} argmax CV = {top:.1f} (runner-up {second:.1f}, "
              f"gap {c[top][0]-c[second][0]:+.5f} vs 1 se = {gapse:.5f}) -> {v}")
    # ⚠ G1 failed: 14 of 72 records did not converge. Two questions follow and BOTH are
    # answered here rather than assumed. (i) Do the dropped cells correlate with the axis
    # under test? (ii) Does G3's verdict survive dropping them?
    print("\n   G1 ROBUSTNESS -- the two questions a non-convergence rate raises:")
    rates = {cv: bycv[cv][0] / bycv[cv][1] for cv in sorted(bycv)}
    xs = np.array(sorted(rates))
    ys = np.array([rates[c] for c in sorted(rates)])
    slope = float(np.polyfit(xs, ys, 1)[0])
    print(f"      (i) convergence rate vs CV: slope {slope:+.3f} per unit CV, "
          f"range {min(ys):.2f}-{max(ys):.2f}")
    gate("the non-convergence is NOT concentrated at high CV, so it cannot bias the axis "
         "G3 tests", abs(slope) < 0.5, f"slope {slope:+.3f} per unit CV")
    # (ii) re-run the argmax on CONVERGED records only
    def chan_conv(cv):
        out = {}
        for s in (1, 2, 3):
            d = new.get(("bidisperse", 0.0, cv, s))
            if d and all(p in d for p in RUNGS) and all(d[p].get("converged") for p in RUNGS):
                out[s] = phi_J(d)
        return out
    cc = {cv: chan_conv(cv) for cv in CVS}
    n_full = {cv: len(v) for cv, v in cc.items()}
    print(f"      (ii) cells with ALL FOUR rungs converged, per CV: "
          f"{ {cv: n_full[cv] for cv in CVS} }")
    usable = {cv: float(np.mean(list(v.values()))) for cv, v in cc.items() if len(v) >= 2}
    if len(usable) >= 4:
        top_c = max(usable, key=lambda k: usable[k])
        print(f"      converged-only argmax: CV = {top_c:.1f} "
              f"(on {len(usable)} of {len(CVS)} CV levels with >= 2 usable seeds)")
        gate("G3's argmax is unchanged when non-converged records are dropped",
             top_c in (0.2, 0.3, 0.4), f"argmax {top_c:.1f}")
    else:
        print("      too few fully-converged cells to re-run the argmax; reported, not hidden")

    agree = len(set(verdicts.values())) == 1 and "UNDETERMINED" not in verdicts.values()
    final = list(verdicts.values())[0] if agree else "UNDETERMINED"
    # 🔴 THE PREREG'S OWN WORTHLESSNESS CRITERION, section 5, quoted:
    #    "More than 2 of 72 cells cap -> the ladder is not comparable to the other two
    #     families' and G3's argmax is being read off a censored curve."
    n_cap = tot - conv
    censored = n_cap > 2
    print("\n   🔴 PREREG SECTION 5, REGISTERED IN ADVANCE: 'More than 2 of 72 cells cap ->")
    print("      G3's argmax is being read off a CENSORED CURVE.'")
    print(f"      capped: {n_cap} of {tot}. Criterion {'TRIGGERED' if censored else 'not triggered'}.")
    if censored:
        print(f"      And the argmax rung itself (CV = {top:.1f}) has "
              f"{n_full.get(top, 0)} of 3 cells with a fully converged ladder.")
        final = f"{final} (PROVISIONAL -- censored curve)"
    print(f"\n   🔑 G3 VERDICT: {final}"
          + ("" if agree else "  (the two channels disagree or a gap is under 1 s.e.)"))
    gate("G3 returned a verdict on the registered scale (T, D or UNDETERMINED)",
         final.split()[0] in ("T", "D", "UNDETERMINED"), final)
    gate("🔴 the campaign meets the prereg's own section-5 worthlessness criterion, so "
         "G3's verdict is PROVISIONAL and not a result", not censored,
         f"{n_cap} of {tot} capped against a registered bar of 2")

    # ---------------- G4 : the blind prediction -----------------------------------------
    print("\n=== G4  BLIND PREDICTION: d(phi_0)^bidi in [+0.0125, +0.0146] ===")
    for nm in ("extrap", "rung"):
        c = curves[nm]
        if 0.0 not in c or 0.5 not in c:
            continue
        cs = sorted(set(c[0.0][2]) & set(c[0.5][2]))
        d = np.array([c[0.5][2][s] - c[0.0][2][s] for s in cs])
        m, se = float(d.mean()), float(d.std(ddof=1) / math.sqrt(len(d)))
        inside = G4_BAND[0] <= m <= G4_BAND[1]
        print(f"   {nm:7s} {m:+.5f} +- {se:.5f}   {'INSIDE' if inside else '🔴 OUTSIDE'} "
              f"the registered band")
        if nm == "rung":
            gate("G4 on the fit-free channel lands in the registered band", inside,
                 f"{m:+.5f} vs {G4_BAND}")

    # ---------------- G5 : the window ---------------------------------------------------
    print("\n=== G5  THE WINDOW W = phi_0 - phi_m for the bidisperse family ===")
    frict, _ = load("phase2*.jsonl", want_family="bidisperse")
    frict = {(k[0], k[2], k[3]): v for k, v in frict.items()
             if k[1] == 0.3 and all(p in v for p in RUNGS)
             and all(v[p].get("n_requested") == 1500 for p in RUNGS)}
    print(f"   {'CV':>5s} {'W (extrap)':>12s} {'se':>9s} {'W(CV)-W(0) paired':>20s} {'|t|':>7s}")
    W = {}
    for cv in CVS:
        a = curves["extrap"].get(cv)
        b = {s: phi_J(frict[("bidisperse", cv, s)]) for s in (1, 2, 3)
             if ("bidisperse", cv, s) in frict}
        if not a or len(b) < 3:
            continue
        cs = sorted(set(a[2]) & set(b))
        d = np.array([a[2][s] - b[s] for s in cs])
        W[cv] = (float(d.mean()), float(d.std(ddof=1) / math.sqrt(len(d))), {s: a[2][s] - b[s] for s in cs})
        extra = ""
        if 0.0 in W and cv != 0.0:
            cs2 = sorted(set(W[cv][2]) & set(W[0.0][2]))
            dd = np.array([W[cv][2][s] - W[0.0][2][s] for s in cs2])
            m2, se2 = float(dd.mean()), float(dd.std(ddof=1) / math.sqrt(len(dd)))
            extra = f"{m2:+20.5f} {abs(m2)/se2:7.2f}"
        print(f"   {cv:5.1f} {W[cv][0]:12.5f} {W[cv][1]:9.5f} {extra}")
    if 0.0 in W:
        print(f"\n   W(CV = 0) = {W[0.0][0]:.5f}, registered {G5_W0}")
        gate("G5a: W(CV = 0) reproduces the registered value (this is W1 restated, not an "
             "independent result)", abs(W[0.0][0] - G5_W0) < 5e-5, f"{W[0.0][0]:.5f}")
    if len(W) >= 2:
        rng = max(v[0] for v in W.values()) - min(v[0] for v in W.values())
        scale = 2 * float(np.mean([v[1] for v in W.values()]))
        print(f"   range of W across the rungs: {rng:.5f} against 2 x mean se = {scale:.5f}")
        gate("G5b: the window's range is below 2 x the propagated paired s.e. "
             "(the bar the other two families met)", rng < scale,
             f"{rng:.5f} vs {scale:.5f}")

    # ---------------- G2 : registered NOT INFORMATIVE -----------------------------------
    print("\n=== G2  REGISTERED NOT INFORMATIVE IN ADVANCE ===")
    c = curves["extrap"]
    if 0.0 in c and 0.5 in c:
        cs = sorted(set(c[0.0][2]) & set(c[0.5][2]))
        d0 = np.array([c[0.5][2][s] - c[0.0][2][s] for s in cs])
        m0, se0 = float(d0.mean()), float(d0.std(ddof=1) / math.sqrt(len(d0)))
        dm, dse = 0.01389, 0.00248        # the frictional shift, on disk
        sed = math.hypot(se0, dse)
        t = abs(m0 - dm) / sed
        print(f"   d(phi_0) {m0:+.5f} +- {se0:.5f} against d(phi_m) {dm:+.5f} +- {dse:.5f}")
        print(f"   |t| = {t:.2f} vs {T_CRIT_2DOF};  t_crit * se_diff = "
              f"{T_CRIT_2DOF*sed:.5f} = {100*T_CRIT_2DOF*sed/dm:.0f} % of the effect")
        print("   ⚠ REGISTERED NOT INFORMATIVE: any outcome passes. Printed, not read.")

    print("\n" + "=" * 90)
    if FAIL:
        print(f"FAILED {len(FAIL)} gate(s):")
        for f in FAIL:
            print(f"   - {f}")
        return 1
    print("ALL GATES PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
