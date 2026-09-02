"""audit_claims.py — re-check the report's headline claims against the current data.

Written as a self-audit after the confinement error was found. Every claim in
REPORT_2026-08-31.md §0 that rests on Phase 4 or Phase 5 data was computed BEFORE the
jamming gate existed, so each needs re-checking against the gated subset, or an explicit
statement of why the gate does not apply.

Prints, per claim: the number as reported, the number on gated data, and a verdict.
"""

from __future__ import annotations

import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
Z_RETENTION_MIN = 0.80


def load(pattern):
    by_key = {}
    for f in sorted(glob.glob(os.path.join(OUT, pattern))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if not r.get("failed") and r.get("key"):
                    by_key[r["key"]] = r
    return list(by_key.values())


def retention(r):
    zc = r.get("z_carved", r.get("z_initial"))
    zl = r.get("z_launch")
    if zc and zl is not None and zc > 0:
        return zl / zc
    return float("nan")


def main() -> None:
    recs = load("phase4*.jsonl") + load("phase5*.jsonl")
    for r in recs:
        r["ret"] = retention(r)
    print(f"loaded {len(recs)} phase-4/5 impact cells\n")

    print("=" * 78)
    print("CLAIM 5: capillary confinement is NOT irrelevant during impact")
    print("  reported as: capillary behaves like caged (spread ~1.0), not free (~3-4),")
    print("               at 13-45 sigma")
    print("=" * 78)
    print()
    print("Does the GATE apply here? NO, and the reason matters:")
    print("  the FREE aggregate disintegrating IS the finding, so gating on retention")
    print("  would exclude exactly the cells the claim is about. The gate is the wrong")
    print("  instrument for this comparison. What CAN be checked is whether the spread")
    print("  separation survives when the confined arms are restricted to jammed cells.")
    print()

    def arm(bc_pred, jammed_only):
        v = [r["spread"] for r in recs
             if bc_pred(r) and np.isfinite(r.get("spread", np.nan))
             and (not jammed_only or (np.isfinite(r["ret"])
                                      and r["ret"] >= Z_RETENTION_MIN))]
        return np.array(v, float)

    for jammed_only in (False, True):
        tag = "jammed cells only" if jammed_only else "all cells (as reported)"
        free = arm(lambda r: r["bc"] == "free", False)   # free is never gated
        cap = arm(lambda r: r["bc"] == "capillary", jammed_only)
        cag = arm(lambda r: r["bc"] == "caged", jammed_only)
        print(f"  --- {tag} ---")
        for name, v in (("free", free), ("capillary", cap), ("caged", cag)):
            if v.size:
                print(f"    {name:>10}: n={v.size:3d} spread {v.mean():.3f} "
                      f"+/- {v.std(ddof=1) if v.size > 1 else float('nan'):.3f} "
                      f"[{v.min():.3f}, {v.max():.3f}]")
        if free.size > 1 and cap.size > 1:
            sd = np.sqrt(free.var(ddof=1) / free.size + cap.var(ddof=1) / cap.size)
            print(f"    free - capillary = {free.mean() - cap.mean():+.3f}, "
                  f"{abs(free.mean() - cap.mean()) / sd:.1f} sigma")
        print()

    print("=" * 78)
    print("WHAT THE CLAIM SHOULD SAY")
    print("=" * 78)
    print("  The confinement in Phase 4/5 was too weak to keep the aggregate JAMMED,")
    print("  but it was still strong enough to keep the grains SPATIALLY TOGETHER.")
    print("  So the correct statement is about cohesion of the body, not about its")
    print("  jamming: any nonzero confinement holds the grains in place, while zero")
    print("  confinement lets them disperse. That is still a falsification of the")
    print("  We-ratio prediction, and for the same reason (zero tensile strength makes")
    print("  the free limit singular), but it is NOT evidence about a jammed drop.")


if __name__ == "__main__":
    main()
