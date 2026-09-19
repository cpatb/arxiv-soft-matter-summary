"""duplicate_tiebreak_robustness_2026-09-19.py

THE QUESTION
------------
`paper/main.tex` discloses in its Methods that the shard set contains keys whose records
DISAGREE on `phi` -- the same physical cell, re-run, returning a different answer -- and
resolves them with a tie-break rule ("highest `max_iter`, then first in file order"). It
quotes `\\conflictWorstDphi` = 0.0040 as the worst disagreement among the SELECTED cells.

What it does not do is ask the question a referee will: **does the choice of which
duplicate survives move the PUBLISHED NUMBERS?** The disclosure bounds the input; it does
not bound the output.

This module enumerates EVERY combination of duplicate choices among the selected cells and
recomputes the three published family shifts under each one. It is a direct, exhaustive
robustness test of the manuscript's headline and it costs nothing.

WHY IT MATTERS NOW
------------------
A referee pass on `wc_endpoints_move_together_2026-09-19.md` observed that the disagreeing
re-runs constitute a direct measurement of the protocol's own reproducibility, and that
this is a different -- and for some purposes more appropriate -- error model than the seed
scatter. Median |delta phi| on re-run is 0.00234 over the 14 disagreeing keys, against a
published seed-paired s.e. of 0.00080 for the uniform family. Those are different
quantities (a per-rung phi against a paired extrapolated shift), so the comparison cannot
be made directly -- which is exactly why the enumeration below is the right test instead of
an argument.

READ-ONLY. The live `out/phase2_wcbidi*.jsonl` campaign is excluded.
"""

from __future__ import annotations

import collections
import glob
import itertools
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(ROOT, "out")
sys.path.insert(0, ROOT)

RUNGS = (1e-3, 3e-4, 1e-4, 3e-5)
MU_REF = 0.3
FAMILIES = ("lognormal", "uniform", "bidisperse")
# the committed values, from paper/numbers.tex
COMMITTED = {"lognormal": 0.0621, "uniform": 0.0242, "bidisperse": 0.0139}
COMMITTED_RATIO = 4.47
FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"   [{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


def phi_J(cell: dict) -> float:
    """Linear extrapolation of phi(p_hat) to p_hat -> 0, the project's `phi_J_lin`."""
    ps = np.array(RUNGS)
    ys = np.array([cell[p]["phi"] for p in RUNGS])
    sl, ic = np.polyfit(ps, ys, 1)
    return float(ic)


def main() -> int:
    print("=" * 88)
    print("Does the duplicate tie-break move the published family shifts?")
    print("=" * 88)

    print("\n=== S0  LOADING, N printed at every filter ===")
    files = [f for f in sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl")))
             if "wcbidi" not in os.path.basename(f)]
    raw = []
    for fn in files:
        for ln in open(fn, encoding="utf-8"):
            ln = ln.strip()
            if ln:
                try:
                    raw.append(json.loads(ln))
                except json.JSONDecodeError:
                    continue
    print(f"   files (wcbidi excluded, mid-flight): {len(files)}")
    print(f"   N records                          : {len(raw)}")
    sel = [r for r in raw if not r.get("failed")]
    print(f"   N after dropping failed            : {len(sel)}")
    sel = [r for r in sel if r.get("n_requested") == 1500]
    print(f"   N after N = 1500                   : {len(sel)}")
    sel = [r for r in sel if r.get("mu_p") == MU_REF]
    print(f"   N after mu_p = {MU_REF}                : {len(sel)}")
    # ⚠ THE `converged` CUT IS LOAD-BEARING AND I LEFT IT OUT ON THE FIRST RUN.
    # `sweep_descriptive.load_shards` filters on it; without it this module saw records
    # like (uniform, CV 0.5, seed 1, p_hat = 3e-5) at phi = 0.478 against 0.603 -- a spread
    # of 0.125 -- and reported the paper's Methods claim ("every conflicting duplicate sits
    # at CV = 0") as FALSE. The claim is true on the paper's own selection; my filter was
    # wrong. Reproduce the selection before contradicting what it supports.
    n_before = len(sel)
    sel = [r for r in sel if r.get("converged")]
    print(f"   N after converged                  : {len(sel)}   "
          f"(dropped {n_before - len(sel)} non-converged)")

    recs: dict = collections.defaultdict(list)
    for r in sel:
        recs[((r.get("family") or "lognormal"), round(r["cv"], 1), r["seed"],
              r["p_hat_target"])].append(r)
    print(f"   N distinct (family, cv, seed, rung) cells: {len(recs)}")

    # ---- S1: which cells have DISAGREEING duplicates -----------------------------------
    print("\n=== S1  CELLS WITH DISAGREEING DUPLICATES, inside the selected set ===")
    multi = {k: v for k, v in recs.items()
             if len({round(x["phi"], 12) for x in v}) > 1}
    print(f"   cells with more than one record   : {sum(1 for v in recs.values() if len(v) > 1)}")
    print(f"   of those, DISAGREEING on phi      : {len(multi)}")
    for k, v in sorted(multi.items()):
        fam, cv, seed, p = k
        print(f"      {fam:11s} cv={cv:.1f} seed={seed} p_hat={p:g}  "
              f"phi = {[round(x['phi'], 6) for x in v]}  "
              f"spread {max(x['phi'] for x in v) - min(x['phi'] for x in v):.5f}")
    cvs = sorted({k[1] for k in multi})
    print(f"\n   CV levels affected: {cvs}")
    check("every disagreeing duplicate in the selected set sits at CV = 0, as main.tex's "
          "Methods states", cvs == [0.0] or not multi, f"{cvs}")
    check("there is something to enumerate -- the test is not vacuous",
          len(multi) > 0, f"{len(multi)} cells")

    # ---- S2: the exhaustive enumeration ------------------------------------------------
    n_combo = int(np.prod([len(v) for v in multi.values()])) if multi else 1
    # ⚠ The product over all disagreeing cells is 2^len(multi) and can be astronomical.
    # Enumerate exhaustively when it is small; otherwise draw a RECORDED, SEEDED random
    # sample and say so. A sampled bound is a bound on the sample, not on the space, and
    # that distinction is printed rather than buried.
    CAP = 4096
    SEED = 20260919
    exhaustive = n_combo <= CAP
    print(f"\n=== S2  DUPLICATE CHOICES: {n_combo} possible ===")
    if exhaustive:
        print(f"   ENUMERATING ALL {n_combo} exhaustively.")
        combos = itertools.product(*[range(len(multi[k])) for k in sorted(multi)])
    else:
        rng = np.random.default_rng(SEED)
        print(f"   🔴 {n_combo} exceeds the cap of {CAP}. Drawing {CAP} SEEDED random")
        print(f"      choices (numpy default_rng seed {SEED}) instead. THE RESULT IS A")
        print("      BOUND ON THE SAMPLE, NOT ON THE SPACE -- the true range can only be")
        print("      wider, so a range that is already negligible stays negligible, but a")
        print("      large one is not an upper bound.")
        ks = sorted(multi)
        combos = [tuple(int(rng.integers(len(multi[k]))) for k in ks) for _ in range(CAP)]
    base = {k: v[0] for k, v in recs.items()}
    keys = sorted(multi)
    out: dict = collections.defaultdict(list)
    for combo in combos:
        pick = dict(base)
        for k, i in zip(keys, combo):
            pick[k] = multi[k][i]
        for fam in FAMILIES:
            ds = []
            for seed in (1, 2, 3):
                try:
                    lo = {p: pick[(fam, 0.0, seed, p)] for p in RUNGS}
                    hi = {p: pick[(fam, 0.5, seed, p)] for p in RUNGS}
                except KeyError:
                    continue
                ds.append(phi_J(hi) - phi_J(lo))
            if len(ds) == 3:
                out[fam].append(float(np.mean(ds)))

    print(f"   {'family':11s} {'combos':>7s} {'min':>10s} {'max':>10s} {'range':>10s} "
          f"{'committed':>10s} {'range/committed':>16s}")
    ok_all = True
    for fam in FAMILIES:
        v = np.array(out[fam])
        rng = v.max() - v.min()
        frac = rng / COMMITTED[fam]
        print(f"   {fam:11s} {len(v):7d} {v.min():10.5f} {v.max():10.5f} {rng:10.5f} "
              f"{COMMITTED[fam]:10.4f} {100*frac:15.1f}%")
        # ⚠ The committed values are macros printed to FOUR DECIMALS, so they cannot be
        # compared against unrounded bounds: `\dPhiJUnif` = 0.0242 against a computed
        # 0.0242013 fails `min <= committed <= max` on the fifth decimal and says nothing.
        # Compare at the precision the manuscript actually prints, and state the tolerance.
        TOL = 5e-5          # half of the last printed digit
        if not (v.min() - TOL <= COMMITTED[fam] <= v.max() + TOL):
            ok_all = False
    check("the committed value of every family shift lies inside the range its duplicate "
          "choices span, at the 4-decimal precision the manuscript prints",
          ok_all, "committed values bracketed by the enumeration to within 5e-5")

    # ---- S3: the headline RATIO, which is what the paper quotes -------------------------
    print("\n=== S3  THE HEADLINE RATIO under every duplicate choice ===")
    lo_r = np.array(out["lognormal"]).min() / np.array(out["bidisperse"]).max()
    hi_r = np.array(out["lognormal"]).max() / np.array(out["bidisperse"]).min()
    print(f"   lognormal/bidisperse spans {lo_r:.2f}x to {hi_r:.2f}x "
          f"against the committed {COMMITTED_RATIO}x")
    print("   ⚠ that is a WORST-CASE pairing (best lognormal against worst bidisperse and")
    print("     vice versa), so it is an upper bound on the tie-break's influence, not an")
    print("     error bar: the two families share their CV = 0 cells, so the numerator and")
    print("     denominator move TOGETHER and the true span is narrower.")
    # the correlated version: same combo for both
    ratios = [a / b for a, b in zip(out["lognormal"], out["bidisperse"])]
    print(f"   taking the SAME combo for both families: {min(ratios):.2f}x to {max(ratios):.2f}x")
    check("the headline ratio's 'three to five' bracket survives every duplicate choice, "
          "on the correlated (correct) pairing",
          3.0 <= min(ratios) and max(ratios) <= 5.0,
          f"{min(ratios):.2f}x to {max(ratios):.2f}x against the paper's stated 'three to five'")

    # ---- S4: the protocol's own re-run scatter, for the record -------------------------
    print("\n=== S4  THE PROTOCOL'S OWN RE-RUN SCATTER (all cells, not just selected) ===")
    alldup: dict = collections.defaultdict(list)
    for r in raw:
        if r.get("failed") or "key" not in r:
            continue
        alldup[r["key"]].append(r)
    conf = {k: v for k, v in alldup.items()
            if len(v) > 1 and max(abs(a["phi"] - b["phi"]) for a in v for b in v) > 1e-12}
    d = [max(abs(a["phi"] - b["phi"]) for a in v for b in v) for v in conf.values()]
    print(f"   keys with >1 record: {sum(1 for v in alldup.values() if len(v) > 1)}; "
          f"DISAGREEING: {len(conf)}")
    print(f"   |delta phi| on re-run: median {np.median(d):.5f}, mean {np.mean(d):.5f}, "
          f"max {max(d):.5f}")
    print("   ⚠ This is a per-RUNG phi scatter. The published shifts are PAIRED differences")
    print("     of EXTRAPOLATED phi_J, so this number does not propagate to them directly --")
    print("     which is why S2's enumeration is the test and this is only context.")
    print(f"   S2 measured the actual propagation: at most "
          f"{max(np.array(out[f]).max() - np.array(out[f]).min() for f in FAMILIES):.5f} "
          "in any family shift.")

    print("\n" + "=" * 88)
    if FAILURES:
        print(f"FAILED {len(FAILURES)} check(s):")
        for f in FAILURES:
            print(f"   - {f}")
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
