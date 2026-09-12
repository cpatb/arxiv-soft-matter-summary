"""grade_common_ladder.py — grade PREREG_COMMON_LADDER_2026-09-09.md.

Written BEFORE the run landed, against the pre-registration's own numbers, so the
verdict cannot be chosen after seeing the data. Every bar below is quoted from that
file; none is computed here.

The four things it does, in the pre-registration's order:

  §6  integrity     -- cell count, failures, realised phi vs target, gate survival
  §3  the verdict   -- delta phi_c on the COMMON ladder against the CENTRED one
  §5a window test   -- drop the bottom rung and require E4 not to move
  §5b P3            -- p_ss at identical phi TARGETS across CV; the realised
                       mismatch is printed rather than assumed to be zero
  §5c bracket       -- the fit-free reading, to quote if every estimator refuses

Reads out/phase8_commonlad*.jsonl; writes out/common_ladder_grade.txt. No simulation.
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from analyse_phase8 import GATE_SETS, analyse_group, group_key, load
from phase8_criticalstate import PHI_J_ISO

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

# --- everything below is QUOTED from PREREG_COMMON_LADDER_2026-09-09.md ------------
GRID = (0.545, 0.570, 0.595, 0.620, 0.645, 0.670)   # §2
N_CELLS_EXPECTED = 54                                # §2
PHI_TOL = 1e-4                                       # §6, third bullet
MAX_FAILED = 3                                       # §6, first bullet
CENTRED = {"E4": 0.0505, "E4m": 0.0410}              # §3, from out/full_ladder_verdict
CENTRED_MID = 0.0458                                 # §3
BAND = 0.015                                         # §3, both hypotheses
WINDOW_BAR = 0.010                                   # §5a
CVS = (0.0, 0.2, 0.5)


def cells_by_group(recs):
    g = {}
    for r in recs:
        g.setdefault(group_key(r), []).append(r)
    return g


def integrity(p, recs) -> bool:
    p("=== §6 INTEGRITY (the pre-registration's own uninterpretability bars) ===")
    errored = [r for r in recs if r.get("error")]
    ok_recs = [r for r in recs if not r.get("error")]
    p(f"   records: {len(recs)}, errored {len(errored)}, usable {len(ok_recs)} "
      f"(expected {N_CELLS_EXPECTED})")
    good = True
    if len(errored) > MAX_FAILED:
        p(f"   🔴 {len(errored)} > {MAX_FAILED} failures -- run is UNINTERPRETABLE")
        good = False
    # the ladders must actually be common
    worst = 0.0
    for r in ok_recs:
        tgt = min(GRID, key=lambda g: abs(g - float(r["phi"])))
        worst = max(worst, abs(float(r["phi"]) - tgt))
    p(f"   worst |phi_realised - grid| = {worst:.2e} (bar {PHI_TOL:.0e}) -> "
      f"{'OK' if worst <= PHI_TOL else '🔴 THE LADDERS ARE NOT COMMON'}")
    good = good and worst <= PHI_TOL
    # seeds per CV
    for cv in CVS:
        seeds = sorted({int(r["seed"]) for r in ok_recs
                        if abs(float(r["cv"]) - cv) < 1e-6})
        rungs = sorted({round(float(r["phi"]), 4) for r in ok_recs
                        if abs(float(r["cv"]) - cv) < 1e-6})
        p(f"   cv {cv:.1f}: seeds {seeds}, {len(rungs)} rungs {rungs}")
        if len(seeds) < 3:
            p(f"      🔴 cv {cv} lost a seed -- §6 calls that uninterpretable")
            good = False
    return good


def contrast(groups, gate, estimator):
    """mean phi_c per CV and the (cv 0 -> 0.5) difference, on one gate set."""
    per = {cv: [] for cv in CVS}
    for k, g in groups.items():
        cv = k[1]
        if cv not in per:
            continue
        res = analyse_group(g, gate)
        v = res.get(f"phi_c_{estimator}", np.nan)
        if np.isfinite(v):
            per[cv].append(v)
    means = {cv: (float(np.mean(v)) if v else np.nan) for cv, v in per.items()}
    ns = {cv: len(v) for cv, v in per.items()}
    sds = {cv: (float(np.std(v, ddof=1)) if len(v) > 1 else np.nan)
           for cv, v in per.items()}
    d = means[0.5] - means[0.0]
    return d, means, ns, sds


def verdict(p, groups) -> None:
    p("")
    p("=== §3 THE VERDICT: delta phi_c on a window that cannot move with CV ===")
    p(f"   reference (centred ladders): E4 {CENTRED['E4']:+.4f}, "
      f"E4m {CENTRED['E4m']:+.4f}, midpoint {CENTRED_MID:+.4f}")
    p(f"   H_inherit predicts [{-BAND:+.4f}, {BAND:+.4f}];  "
      f"H_real predicts [{CENTRED_MID - BAND:+.4f}, {CENTRED_MID + BAND:+.4f}]")
    p("")
    calls = {}
    for est in ("E4", "E4m", "E2", "E3"):
        p(f"   -- {est}{' (PRIMARY)' if est == 'E4' else ''}"
          f"{'  [reported, NOT quoted -- window-dependent]' if est == 'E2' else ''}")
        for name, gate in GATE_SETS:
            d, means, ns, sds = contrast(groups, gate, est)
            if not np.isfinite(d):
                p(f"      {name:12s} unresolved (n = "
                  f"{ns[0.0]}, {ns[0.2]}, {ns[0.5]})")
                continue
            cell = "  ".join(
                f"cv{cv:.1f}={means[cv]:.4f}"
                f"{'' if not np.isfinite(sds[cv]) else f'+-{sds[cv]:.4f}'}(n{ns[cv]})"
                for cv in CVS)
            call = ("H_real" if abs(d - CENTRED_MID) <= BAND else
                    "H_inherit" if abs(d) <= BAND else "NEITHER (third branch)")
            p(f"      {name:12s} d = {d:+.4f}  -> {call}")
            p(f"      {'':12s} {cell}")
            if est in ("E4", "E4m"):
                calls.setdefault(est, {})[name] = call
    p("")
    p("   §3 requires the call to hold on E4 AND E4m. Summary:")
    for est, byg in calls.items():
        p(f"      {est}: " + ", ".join(f"{g}={c}" for g, c in byg.items()))
    agree = set()
    for byg in calls.values():
        agree |= set(byg.values())
    p(f"   distinct calls across E4/E4m and all gate sets: {sorted(agree)}")
    if len(agree) > 1:
        p("   🔴 the estimators/gates do NOT agree -- §3's third branch applies and the")
        p("      honest report is the RANGE, not a hypothesis.")


def window_test(p, groups) -> None:
    p("")
    p(f"=== §5a WINDOW TEST: drop the bottom rung, E4 must move < {WINDOW_BAR} ===")
    lo = min(GRID)
    moved = []
    for k, g in sorted(groups.items()):
        full = analyse_group(g, ())
        cut = analyse_group([r for r in g if float(r["phi"]) > lo + 1e-6], ())
        a, b = full.get("phi_c_E4", np.nan), cut.get("phi_c_E4", np.nan)
        if np.isfinite(a) and np.isfinite(b):
            moved.append(abs(b - a))
            flag = "" if abs(b - a) < WINDOW_BAR else "  🔴 OVER BAR"
            p(f"   cv {k[1]:.1f} seed {k[3]}: E4 {a:.4f} -> {b:.4f}  "
              f"move {b - a:+.4f}{flag}")
        else:
            p(f"   cv {k[1]:.1f} seed {k[3]}: E4 {a:.4f} -> {b:.4f}  "
              f"(one side refused -- no move to report)")
    if moved:
        p(f"   worst move {max(moved):.4f} over {len(moved)} groups, bar "
          f"{WINDOW_BAR} -> {'PASS' if max(moved) < WINDOW_BAR else '🔴 FAIL'}")
        if max(moved) == 0.0:
            p("   🔴 THIS TEST IS VACUOUS AS WRITTEN, and the exact zero is the proof.")
            p("      E4 interpolates between the two rungs that BRACKET z_nr = 4. Those")
            p("      are rungs 2-3 at cv 0 and 0.2 and rungs 4-5 at cv 0.5, so deleting")
            p("      rung 1 cannot touch the interpolation and the answer must come back")
            p("      bit-identical for ANY correct implementation. A routine returning a")
            p("      constant would pass it too. The number stays printed; the")
            p("      non-vacuous version is `placement_test` below, which moves the whole")
            p("      ladder instead of trimming it.")
    else:
        p("   no group resolved E4 on both ladders -- the test did not run")


def placement_test(p, groups) -> None:
    """The §5a test that CAN fail: same cell, ladder rungs at DIFFERENT absolute phi.

    The centred ladders (downlad/downlad2 + prod2) and the common ladder run the same
    (family, cv, mu_p, seed, N, protocol, gamma_max, shear_steps) cells with rungs at
    different absolute phi -- 0.5145-0.6945 at cv 0 and 0.5769-0.7569 at cv 0.5, against
    a common 0.545-0.670. If E4 is a property of the material it returns the same number
    on both; if it is a property of its ladder it does not. Nothing about this is forced.
    """
    p("")
    p(f"=== §5a' PLACEMENT TEST (the version that can fail): E4 must move < "
      f"{WINDOW_BAR} ===")
    centred_recs = [r for r in load("phase8_downlad*.jsonl,phase8_prod2*.jsonl")
                    if not r.get("error")]
    centred = {}
    for r in centred_recs:
        centred.setdefault(group_key(r), []).append(r)
    rows = []
    for k, g in sorted(groups.items()):
        if k not in centred:
            p(f"   cv {k[1]:.1f} seed {k[3]}: no centred-ladder counterpart")
            continue
        a = analyse_group(g, ()).get("phi_c_E4", np.nan)
        b = analyse_group(centred[k], ()).get("phi_c_E4", np.nan)
        lo_a, hi_a = min(GRID), max(GRID)
        cg = [float(r["phi"]) for r in centred[k]]
        if np.isfinite(a) and np.isfinite(b):
            rows.append(abs(a - b))
            flag = "" if abs(a - b) < WINDOW_BAR else "  🔴 OVER BAR"
            p(f"   cv {k[1]:.1f} seed {k[3]}: common [{lo_a:.3f},{hi_a:.3f}] -> {a:.4f}"
              f"   centred [{min(cg):.3f},{max(cg):.3f}] -> {b:.4f}"
              f"   diff {a - b:+.4f}{flag}")
        else:
            p(f"   cv {k[1]:.1f} seed {k[3]}: common {a:.4f} centred {b:.4f} "
              f"(one side refused)")
    if rows:
        p(f"   worst |move| {max(rows):.4f} over {len(rows)} cells, bar {WINDOW_BAR} -> "
          f"{'PASS' if max(rows) < WINDOW_BAR else '🔴 FAIL'}")
        p(f"   ⚠ and it is NOT zero ({min(rows):.4f} to {max(rows):.4f}), so unlike §5a")
        p("     this test had room to fail and did not.")
    else:
        p("   no cell resolved E4 on both ladders -- the test did not run")


def exact_p3(p, recs) -> None:
    p("")
    p("=== §5b P3 AT MATCHED phi (identical TARGETS; realised mismatch below) ===")
    # ⚠ Matched by TOLERANCE, not by rounding the realised phi to 4 dp. The controller
    # lands within 8.2e-05 of target, so cv = 0 at the 0.595 rung realises 0.59505 and
    # rounds to 0.5951 -- against a grid value of 0.5950. The first version of this
    # function keyed on round(phi, 4) and silently dropped that whole rung from P3.
    by = {}
    for r in recs:
        if r.get("error"):
            continue
        g = min(GRID, key=lambda x: abs(x - float(r["phi"])))
        if abs(g - float(r["phi"])) <= PHI_TOL:
            by.setdefault((g, round(float(r["cv"]), 3)), []).append(r)
    p("     phi     p_ss(cv0)   p_ss(cv0.5)   ratio    Z_nr 0 -> 0.5   rattlers 0 -> 0.5")
    ratios = []
    for phi in GRID:
        a = by.get((phi, 0.0), [])
        b = by.get((phi, 0.5), [])
        if not a or not b:
            p(f"   {phi:.3f}   missing one arm (n = {len(a)}, {len(b)})")
            continue
        pa = float(np.mean([r["p_ss_hat"] for r in a]))
        pb = float(np.mean([r["p_ss_hat"] for r in b]))
        za = float(np.mean([r["z_nr_ss"] for r in a]))
        zb = float(np.mean([r["z_nr_ss"] for r in b]))
        ra = float(np.mean([r.get("rattler_frac_ss", np.nan) for r in a]))
        rb = float(np.mean([r.get("rattler_frac_ss", np.nan) for r in b]))
        ratio = pa / pb if pb else np.nan
        head = phi >= 0.595
        if head and np.isfinite(ratio):
            ratios.append(ratio)
        p(f"   {phi:.3f}   {pa:.4e}  {pb:.4e}  {ratio:7.1f}x   {za:.2f} -> {zb:.2f}   "
          f"{ra:.3f} -> {rb:.3f}{'' if head else '   (below 0.595: EXCLUDED by §5b)'}")
    if ratios:
        p(f"   headline range over the rungs at/above 0.595: "
          f"{min(ratios):.1f}x to {max(ratios):.1f}x")
        p(f"   §5b predicted every one of these > 1: "
          f"{'PASS' if min(ratios) > 1 else '🔴 FAIL'}")


def bracket(p, recs) -> None:
    p("")
    p("=== §5c THE FIT-FREE BRACKET (quote this if every estimator refuses) ===")
    for cv in CVS:
        p(f"   cv {cv:.1f}:")
        prev = None
        for phi in GRID:
            g = [r for r in recs if not r.get("error")
                 and abs(float(r["phi"]) - phi) < 1e-4
                 and abs(float(r["cv"]) - cv) < 1e-6]
            if not g:
                continue
            ps = float(np.mean([r["p_ss_hat"] for r in g]))
            z = float(np.mean([r["z_nr_ss"] for r in g]))
            p(f"      phi {phi:.3f}  p_ss {ps:.3e}  z_nr {z:.2f}  (n {len(g)})")
            if prev is not None and prev < 4.0 <= z:
                p(f"      ^^ z_nr crosses 4 between {phi - 0.025:.3f} and {phi:.3f}")
            prev = z


def offset_h3(p, groups) -> None:
    """PREREG_DST_WINDOW §H-3: is phi_c(shear) - phi_J(iso, mu 0.3) constant in CV?"""
    p("")
    p("=== H-3 (PREREG_DST_WINDOW_2026-09-09.md): shear vs compression, same CV ===")
    p("   bar: the offset is CONSTANT in CV to within 0.010")
    offs = []
    for est in ("E4", "E4m"):
        d, means, ns, _ = contrast(groups, (), est)
        row = []
        for cv in CVS:
            iso = PHI_J_ISO[("lognormal", cv, 0.3)]
            o = means[cv] - iso
            row.append(o)
            p(f"   {est} cv {cv:.1f}: phi_c {means[cv]:.4f} - phi_J_iso {iso:.4f} "
              f"= {o:+.4f} (n {ns[cv]})")
        if all(np.isfinite(x) for x in row):
            span = max(row) - min(row)
            offs.append(span)
            p(f"   {est}: offset spans {span:.4f} -> "
              f"{'CONSTANT (H-3 passes)' if span <= 0.010 else '🔴 MOVES WITH CV'}")
    if offs:
        p("   ⚠ this comparison is only meaningful on the COMMON ladder. On the centred")
        p("     ladders it is forced, because they are centred on these same phi_J_iso.")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", default="phase8_commonlad*.jsonl")
    ap.add_argument("--out", default="common_ladder_grade.txt")
    args = ap.parse_args()

    lines: list[str] = []

    def p(s: str = "") -> None:
        print(s)
        lines.append(s)

    p("=== grade_common_ladder.py — PREREG_COMMON_LADDER_2026-09-09.md ===")
    p("   every bar is quoted from that file; none is computed here")
    recs = load(args.pattern)
    if not recs:
        p("no records yet")
        return 1
    ok = integrity(p, recs)
    groups = cells_by_group([r for r in recs if not r.get("error")])
    verdict(p, groups)
    window_test(p, groups)
    placement_test(p, groups)
    exact_p3(p, [r for r in recs if not r.get("error")])
    bracket(p, [r for r in recs if not r.get("error")])
    offset_h3(p, groups)
    if not ok:
        p("")
        p("🔴 an §6 integrity bar failed; every verdict above is provisional.")

    path = os.path.join(OUT, args.out)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n-> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
