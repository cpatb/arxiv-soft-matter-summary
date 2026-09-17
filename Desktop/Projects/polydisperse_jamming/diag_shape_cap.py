"""diag_shape_cap.py -- WHY the CV = 0.5 arm of the pressure ladder does not converge.

The shape campaign's registered ladder was shortened once already (amendment 5.6) because
`p_hat` = 3e-5 did not jam at CV = 0.5 with N = 5400. A referee pointed out that the
evidence covered only the 1e-3 rung there, and the production run has now shown the
**1e-4 rung fails too, 0 of 6**. So the shortened ladder has the same defect it was
introduced to remove, and the question is which of two very different things is
happening:

  (a) THE PACK NEVER JAMS -- `z_all` near zero, still a gas when the run stops. Then no
      iteration budget helps and the state is not a measurement of anything.
  (b) THE PACK JAMS BUT THE RUN IS CUT OFF -- `z_all` in the normal 3-5 range, the
      unbalanced-force criterion simply not met inside `MAX_ITER`. Then the cap is the
      binding constraint, and `MAX_ITER = 200000` -- calibrated when this project ran
      N = 1500 -- is an N-DEPENDENT criterion wearing a constant's clothes.

The two call for opposite responses, and `converged = False` cannot tell them apart.
`PackingState.tol_met` records whether the unbalanced bar was met IGNORING the cap, which
is exactly the discriminator.
"""

from __future__ import annotations

import glob
import json
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> int:
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, "out", "shape",
                                           "phase2_shp*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    r = json.loads(line)
                    if not r.get("failed"):
                        rows.append(r)
    print(f"N: {len(rows)} records that ran")

    bad = [r for r in rows if not r.get("converged")]
    print(f"N: -> {len(bad)} non-converged")
    print()
    print(f"{'fam':<11}{'a':>5}{'cv':>5}{'p_hat':>9}{'s':>3}{'iters':>8}"
          f"{'cap':>8}{'at cap':>8}{'unbal':>11}{'tol_met':>9}{'z_all':>8}"
          f"{'phi':>9}")
    for r in sorted(bad, key=lambda x: (x["cv"], x["alpha"], -x["p_hat_target"])):
        at_cap = int(r["iters"]) >= int(r.get("max_iter", 200_000))
        print(f"{r['family']:<11}{r['alpha']:>5.2f}{r['cv']:>5.2f}"
              f"{r['p_hat_target']:>9.0e}{r['seed']:>3}{r['iters']:>8}"
              f"{r.get('max_iter', 0):>8}{str(at_cap):>8}{r['unbalanced']:>11.3e}"
              f"{str(r.get('tol_met')):>9}{r['z_all']:>8.3f}{r['phi']:>9.4f}")

    print()
    n_cap = sum(1 for r in bad
                if int(r["iters"]) >= int(r.get("max_iter", 200_000)))
    n_tol = sum(1 for r in bad if r.get("tol_met"))
    z = np.array([float(r["z_all"]) for r in bad])
    n_gas = int(np.sum(z < 1.0))
    print(f"  of the {len(bad)} non-converged cells: {n_cap} stopped AT THE CAP, "
          f"{n_tol} met the unbalanced bar anyway, {n_gas} have z_all < 1 (a gas, "
          f"not a slow packing)")
    print()
    print("  z_all of the non-converged cells, split:")
    for lo, hi, label in ((0.0, 1.0, "z < 1      (never jammed)"),
                          (1.0, 3.0, "1 <= z < 3 (marginal)"),
                          (3.0, 99.0, "z >= 3     (a real packing, cut off)")):
        sel = [r for r in bad if lo <= float(r["z_all"]) < hi]
        if sel:
            zz = np.array([float(r["z_all"]) for r in sel])
            ps = sorted({float(r["p_hat_target"]) for r in sel})
            print(f"    {label}: {len(sel)} cells, z in [{zz.min():.3f}, {zz.max():.3f}]"
                  f", rungs {['%.0e' % p for p in ps]}")

    # the CONTROL: the same rungs at CV = 0, which do converge
    print()
    print("  CONTROL -- the same rungs at CV = 0.0, which is where they DO converge:")
    good = [r for r in rows if r.get("converged") and abs(float(r["cv"])) < 1e-12]
    by = defaultdict(list)
    for r in good:
        by[(round(float(r["alpha"]), 2), float(r["p_hat_target"]))].append(
            int(r["iters"]))
    for k in sorted(by):
        v = np.array(by[k])
        print(f"    alpha={k[0]:.2f} p_hat={k[1]:.0e}: n={v.size} iters "
              f"median={int(np.median(v))} max={v.max()} "
              f"(cap 200000, i.e. {100*v.max()/200000:.0f}% of it)")

    print()
    print("READ THIS AS: if the failures are AT THE CAP with a normal z_all, the cap is "
          "the binding constraint and it was calibrated at N = 1500; if z_all is near "
          "zero the pack never jammed and no budget fixes it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
