"""diag_shape_cv0.py -- is the CV = 0 rung the SAME CELL in all three families?

At CV = 0 every PSD family is the same delta distribution, so `psd.make_psd` should
return identical radii for the same seed and the packing should be bit-identical. If so,
running CV = 0 once per family is running one cell three times -- a third of the shape
grid -- and the three families' D = phi_J(CV 0.5) - phi_J(CV 0) share a common reference
point, which makes their interactions CORRELATED rather than independent.

Both halves matter, and the second one is not a saving, it is a caveat. This script
settles the question from the records already on disk and from `psd.make_psd` directly,
and it prints the DISCRIMINATING check as well: at CV = 0.5 the families must DIFFER,
or an identity found at CV = 0 would just mean the family flag does nothing.
"""

from __future__ import annotations

import glob
import json
import os
from collections import defaultdict

import numpy as np

from psd import make_psd

HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> int:
    print("=== 1. the PSD itself, straight from make_psd (no Yade, no records) ===")
    for cv in (0.0, 0.5):
        ref = None
        print(f"  CV = {cv}")
        for fam in ("lognormal", "uniform", "bidisperse"):
            rng = np.random.default_rng(1)
            r = make_psd(fam, 1.0, cv).radii_for_yade(5400, rng, truncate=(0.2, 5.0))
            h = float(np.sum(r)), float(r.min()), float(r.max())
            if ref is None:
                ref = (r, h)
                print(f"    {fam:<11} sum={h[0]!r} min={h[1]!r} max={h[2]!r}")
            else:
                same = r.shape == ref[0].shape and bool(np.array_equal(r, ref[0]))
                print(f"    {fam:<11} sum={h[0]!r} min={h[1]!r} max={h[2]!r}"
                      f"   bit-identical to lognormal: {same}")

    print()
    print("=== 2. the RECORDS on disk ===")
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, "out", "shape",
                                           "phase2_shp_*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    r = json.loads(line)
                    if not r.get("failed"):
                        rows.append(r)
    by = defaultdict(dict)
    for r in rows:
        k = (round(float(r["cv"]), 3), round(float(r["alpha"]), 3), int(r["seed"]),
             float(r["p_hat_target"]))
        by[k][r["family"]] = r
    n_cmp = 0
    n_same = 0
    worst = 0.0
    print(f"  {'CV':>5}{'alpha':>7}{'seed':>5}{'p_hat':>9}  families        "
          f"max |dphi|   bit-identical")
    for k in sorted(by):
        fams = by[k]
        if len(fams) < 2:
            continue
        n_cmp += 1
        vals = {f: float(r["phi"]) for f, r in fams.items()}
        ref = list(vals.values())[0]
        d = max(abs(v - ref) for v in vals.values())
        bit = all(v == ref for v in vals.values())
        n_same += bool(bit)
        worst = max(worst, d)
        print(f"  {k[0]:>5.2f}{k[1]:>7.2f}{k[2]:>5}{k[3]:>9.0e}  "
              f"{','.join(sorted(vals)):<16}{d:>10.3e}   {bit}")
    print()
    print(f"  N: {n_cmp} cells present in more than one family; {n_same} bit-identical; "
          f"worst |dphi| = {worst:.3e}")
    print()
    print("READ THIS AS: an identity at CV = 0 is a fact about the MODEL (every family "
          "is the same delta distribution there), not about the solver. Two "
          "consequences, and only one of them is a saving:")
    print("  (a) CV = 0 need only be run ONCE, not once per family;")
    print("  (b) the three families' D = phi_J(CV_hi) - phi_J(CV 0) SHARE their "
          "reference point, so their interactions are CORRELATED and a cross-family "
          "comparison of D cannot treat them as three independent measurements.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
