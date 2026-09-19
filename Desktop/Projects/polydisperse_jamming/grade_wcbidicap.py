"""grade_wcbidicap.py -- PREREG_WCBIDI_2026-09-19.md, AMENDMENT 3 (the cap re-run).

The first `wcbidi` campaign delivered 58 of 72 converged, which triggers the prereg's own
section-5 worthlessness criterion and made G3's (T) verdict PROVISIONAL. This grades the
re-run at `PJ_MAX_ITER = 1000000`, tag `wcbidicap`.

The amendment's gates, quoted:

  C1 IDENTITY (PRIMARY) -- "the cells that ALREADY converged under the 200000 cap must come
     back BIT-IDENTICAL under the 1000000 one. A cell that stopped because it converged
     cannot be affected by raising a limit it never reached; if any such cell moves, the cap
     is not the only thing that changed and nothing from this re-run may be read."
     REGISTERED: C1 passes on all 58 previously-converged records.
  C2 -- the 14 previously-capped cells converge under the raised cap. REGISTERED: 14 of 14.
  C3 -- G3's argmax re-graded on CONVERGED-ONLY data. REGISTERED: it lands in
     {0.2, 0.3, 0.4} (the provisional (T) verdict survives); the alternative is 0.5.
  C4 -- d(phi_0) against G4's ORIGINAL registered band [+0.0125, +0.0146]. Not moved.

⚠ C1 IS GRADED FIRST AND IT GATES EVERYTHING ELSE. If it fails, C2-C4 are printed but
explicitly NOT READ, because the premise that only the cap changed is what licenses them.
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

RUNGS = (1e-3, 3e-4, 1e-4, 3e-5)
CVS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
G4_BAND = (0.0125, 0.0146)
OLD_CAP, NEW_CAP = 200_000, 1_000_000
FIELDS = ("phi", "z_all", "z_nr", "n_contacts")
FAIL: list[str] = []


def gate(name: str, ok: bool, detail: str = "") -> None:
    print(f"   [{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


def load(pat: str):
    cells: dict = collections.defaultdict(dict)
    n = 0
    for fn in sorted(glob.glob(os.path.join(OUT, pat))):
        for ln in open(fn, encoding="utf-8"):
            ln = ln.strip()
            if not ln:
                continue
            r = json.loads(ln)
            n += 1
            if r.get("failed"):
                continue
            cells[(round(r["cv"], 1), r["seed"])][r["p_hat_target"]] = r
    return cells, n


def phi_J(d):
    return float(np.polyfit(np.array(RUNGS), np.array([d[p]["phi"] for p in RUNGS]), 1)[1])


def main() -> int:
    print("=" * 92)
    print("grade_wcbidicap.py -- PREREG_WCBIDI_2026-09-19.md Amendment 3 (the cap re-run)")
    print("=" * 92)

    old, n_old = load("phase2_wcbidi_bidisperse_shard*.jsonl")
    new, n_new = load("phase2_wcbidicap*.jsonl")
    print(f"\n   original wcbidi records : {n_old}")
    print(f"   re-run wcbidicap records: {n_new}")
    if n_new == 0:
        print("\n   🔴 NO RE-RUN RECORDS ON DISK. The campaign has not delivered; nothing "
              "is graded.")
        return 1
    caps = {r.get("max_iter") for d in new.values() for r in d.values()}
    print(f"   caps present in the re-run: {sorted(c for c in caps if c is not None)}")
    gate("the re-run really used the raised cap", caps == {NEW_CAP}, f"{sorted(caps)}")

    # ---------------- C1 : IDENTITY, and it gates everything ----------------------------
    print("\n=== C1  IDENTITY (PRIMARY): previously-CONVERGED cells must be BIT-IDENTICAL ===")
    n_cmp = 0
    mism = []
    for key, od in old.items():
        for p, orec in od.items():
            if not orec.get("converged"):
                continue                      # C2's population, not C1's
            nrec = new.get(key, {}).get(p)
            if nrec is None:
                continue
            n_cmp += 1
            for fld in FIELDS:
                a, b = nrec.get(fld), orec.get(fld)
                if a is None or b is None:
                    continue
                if a != b:
                    mism.append((key, p, fld, a, b))
    print(f"   previously-converged records compared: {n_cmp}  (registered 58)")
    if mism:
        print(f"   🔴 {len(mism)} mismatches, first 8:")
        for m in mism[:8]:
            print(f"      cv={m[0][0]:.1f} seed={m[0][1]} p={m[1]:g} {m[2]}: "
                  f"new {m[3]!r} vs old {m[4]!r}")
    c1 = not mism and n_cmp > 0
    gate("C1: raising the cap did not move a single cell that had already converged",
         c1, f"{n_cmp} compared, {len(mism)} mismatches")
    if not c1:
        print("\n   🔴 C1 FAILED. Per the amendment, C2-C4 below are PRINTED BUT NOT READ:")
        print("      the premise that only the cap changed is what licenses them.")

    # ---------------- C2 : do the capped cells now converge? ----------------------------
    print("\n=== C2  the 14 previously-capped cells under the raised cap ===")
    was_capped = [(k, p) for k, d in old.items() for p, r in d.items()
                  if not r.get("converged")]
    print(f"   previously capped: {len(was_capped)}  (registered 14)")
    now_ok = sum(1 for k, p in was_capped
                 if new.get(k, {}).get(p, {}).get("converged"))
    still = [(k, p) for k, p in was_capped
             if not new.get(k, {}).get(p, {}).get("converged")]
    print(f"   of those, converged at {NEW_CAP}: {now_ok} of {len(was_capped)}")
    for k, p in still:
        it = new.get(k, {}).get(p, {}).get("iters")
        print(f"      🔴 STILL CAPPED: cv={k[0]:.1f} seed={k[1]} p={p:g}  iters={it}")
    gate("C2: every previously-capped cell converges under the raised cap",
         now_ok == len(was_capped) and len(was_capped) > 0,
         f"{now_ok}/{len(was_capped)}")

    tot = sum(len(d) for d in new.values())
    conv = sum(1 for d in new.values() for r in d.values() if r.get("converged"))
    print(f"\n   re-run overall: {conv} of {tot} converged")
    full = {cv: sum(1 for s in (1, 2, 3)
                    if (cv, s) in new and all(p in new[(cv, s)] for p in RUNGS)
                    and all(new[(cv, s)][p].get("converged") for p in RUNGS))
            for cv in CVS}
    print(f"   cells with ALL FOUR rungs converged, per CV: {full}")

    # ---------------- C3 : the point of the exercise ------------------------------------
    print("\n=== C3  G3's ARGMAX on CONVERGED-ONLY data (the registered test) ===")

    def curve(conv_only: bool, use_rung: bool):
        out = {}
        for cv in CVS:
            vals = {}
            for s in (1, 2, 3):
                d = new.get((cv, s))
                if not d or not all(p in d for p in RUNGS):
                    continue
                if conv_only and not all(d[p].get("converged") for p in RUNGS):
                    continue
                vals[s] = d[3e-05]["phi"] if use_rung else phi_J(d)
            if len(vals) >= 2:
                a = np.array(list(vals.values()))
                out[cv] = (float(a.mean()), float(a.std(ddof=1) / math.sqrt(len(a))), vals)
        return out

    verdicts = {}
    for label, (co, ur) in (("all, extrap", (False, False)),
                            ("all, rung", (False, True)),
                            ("CONVERGED, extrap", (True, False)),
                            ("CONVERGED, rung", (True, True))):
        c = curve(co, ur)
        if len(c) < 4:
            print(f"   {label:20s}: only {len(c)} CV levels usable -- not graded")
            continue
        order = sorted(c, key=lambda k: c[k][0], reverse=True)
        top, second = order[0], order[1]
        gapse = math.hypot(c[top][1], c[second][1])
        near = abs(c[top][0] - c[second][0]) < gapse
        v = "UNDETERMINED" if near else ("T" if top in (0.2, 0.3, 0.4) else
                                         ("D" if top == 0.5 else "?"))
        verdicts[label] = v
        print(f"   {label:20s}: {len(c)} CV levels, argmax {top:.1f} (runner-up {second:.1f},"
              f" gap {c[top][0]-c[second][0]:+.5f} vs 1 se {gapse:.5f}) -> {v}")
    conv_v = [v for k, v in verdicts.items() if k.startswith("CONVERGED")]
    if conv_v:
        agree = len(set(conv_v)) == 1
        print(f"\n   🔑 C3 VERDICT (converged-only): "
              f"{conv_v[0] if agree else 'UNDETERMINED (channels disagree)'}")
        gate("C3: the provisional (T) verdict SURVIVES on converged-only data",
             agree and conv_v[0] == "T", f"{conv_v}")
    else:
        gate("C3: enough converged cells to re-grade the argmax at all", False,
             "no CV level has >= 2 fully-converged cells across the ladder")

    # ---------------- C4 : the original blind band --------------------------------------
    print("\n=== C4  d(phi_0) against G4's ORIGINAL registered band ===")
    for label, (co, ur) in (("all, rung", (False, True)), ("CONVERGED, rung", (True, True))):
        c = curve(co, ur)
        if 0.0 not in c or 0.5 not in c:
            print(f"   {label}: endpoints unavailable")
            continue
        cs = sorted(set(c[0.0][2]) & set(c[0.5][2]))
        if len(cs) < 2:
            print(f"   {label}: fewer than 2 paired seeds")
            continue
        d = np.array([c[0.5][2][s] - c[0.0][2][s] for s in cs])
        m, se = float(d.mean()), float(d.std(ddof=1) / math.sqrt(len(d)))
        inside = G4_BAND[0] <= m <= G4_BAND[1]
        print(f"   {label:20s} {m:+.5f} +- {se:.5f} (n={len(cs)})  "
              f"{'INSIDE' if inside else '🔴 OUTSIDE'} {G4_BAND}")

    print("\n" + "=" * 92)
    if FAIL:
        print(f"FAILED {len(FAIL)} gate(s):")
        for f in FAIL:
            print(f"   - {f}")
        return 1
    print("ALL GATES PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
