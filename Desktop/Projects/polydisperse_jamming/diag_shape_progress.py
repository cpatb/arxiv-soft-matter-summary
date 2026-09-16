"""diag_shape_progress.py -- read the shape campaign's shards as they land.

A PROGRESS reader, not a grader: it prints what is on disk, checks S0c's predicted
double count and the mass on every clump cell it can see, and says how much of the
registered grid is still missing. `grade_shape.py` is the only thing that issues a
verdict.

The S0c check here reproduces the grader's own arithmetic on the same records rather
than a private variant of it, so a disagreement between this and `grade_shape.py` is a
bug in one of them and not a matter of interpretation (lessons.md 2026-09-08: a
diagnostic must reproduce the number it audits).
"""

from __future__ import annotations

import glob
import json
import os
from collections import defaultdict

import numpy as np

import grade_shape
from shape import sub_sphere_diameter

HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> int:
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, "out", "shape",
                                           "phase2_shp_*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    rows.append(json.loads(line))
    ok = [r for r in rows if not r.get("failed")]
    bad = [r for r in rows if r.get("failed")]
    print(f"N: {len(rows)} production records read; {len(ok)} ran, {len(bad)} raised")
    for r in bad[:5]:
        print(f"  FAILED: {str(r.get('error'))[:200]}")

    if not ok:
        return 0

    print()
    print(f"{'family':<11}{'a':>5}{'cv':>5}{'p_hat':>9}{'s':>3}{'phi':>10}{'z_all':>8}"
          f"{'dup':>8}{'|ratio-pred|':>14}{'|mass-1|':>11}{'conv':>6}{'wall':>8}")
    worst_ratio = 0.0
    worst_mass = 0.0
    for r in sorted(ok, key=lambda x: (x["family"], x["alpha"], x["cv"],
                                       -x["p_hat_target"], x["seed"])):
        a = float(r["alpha"])
        # At alpha = 1 the grain is a SINGLE SPHERE, so there is no member overlap and
        # the predicted double count is 1, not 2*(ds/d_eq)^3 = 2. `grade_shape.g_s0c`
        # reads S0c at alpha > 1 only, for exactly this reason; the first version of
        # this reader applied the clump prediction at every alpha and printed a
        # |ratio - pred| of 1.00e+00 on sound sphere cells.
        pred = 2.0 * sub_sphere_diameter(1.0, a) ** 3 if a > 1.0 else 1.0
        ratio = float(r["phi_porosity"]) / float(r["phi_union"])
        dr = abs(ratio - pred)
        dm = abs(float(r["mass_over_vol_solid"]) - 1.0)
        if a > 1.0:
            worst_ratio = max(worst_ratio, dr)
        worst_mass = max(worst_mass, dm)
        print(f"{r['family']:<11}{a:>5.2f}{r['cv']:>5.1f}{r['p_hat_target']:>9.0e}"
              f"{r['seed']:>3}{r['phi']:>10.5f}{r['z_all']:>8.3f}"
              f"{r['dup_contact_frac']:>8.4f}{dr:>14.2e}{dm:>11.2e}"
              f"{str(r['converged']):>6}{r['wall_s']:>7.0f}s")

    print()
    n_clump = sum(1 for r in ok if float(r["alpha"]) > 1.0)
    print(f"S0c so far, over the {n_clump} cells at alpha > 1: "
          f"worst |ratio - 2(ds/d_eq)^3| = {worst_ratio:.3e} "
          f"(bar {grade_shape.S0C_RATIO_BAR:g}), worst |mass/vol - 1| = "
          f"{worst_mass:.3e} (bar {grade_shape.S0C_MASS_BAR:g})")
    if n_clump == 0:
        print("  -- NO CLUMP CELLS YET, so S0c is UNTESTED. The sphere rung's ratio is "
              "1.0 by construction and proves nothing about the clump construction.")

    # dup_contact_frac in a JAMMED pack, which is the number the loose probe could not
    # give. Reported per alpha because it is the confound's own size.
    by_a = defaultdict(list)
    for r in ok:
        by_a[round(float(r["alpha"]), 3)].append(float(r["dup_contact_frac"]))
    print()
    print("dup_contact_frac in the PRODUCTION packs (the loose probe measured 0.75 / "
          "0.21 / 0.08 / 0.13 / 0.05 at z ~ 0.2-1.0, which is not this):")
    for a in sorted(by_a):
        v = np.array(by_a[a])
        print(f"  alpha={a:.2f}  n={v.size}  mean={v.mean():.4f}  "
              f"min={v.min():.4f}  max={v.max():.4f}")

    # what is still missing from the registered wave-1 grid
    # The REGISTERED wave-1 grid: three rungs, not four (amendment 5.6). Read from the
    # grader's own constant so the two cannot drift.
    want = {(f, a, c, s, p)
            for f in ("lognormal", "uniform", "bidisperse")
            for a in (1.0, 2.0) for c in (0.0, 0.5) for s in (1, 2, 3)
            for p in grade_shape.P_HAT_LADDER_SHAPE}
    have = {(r["family"], round(float(r["alpha"]), 3), round(float(r["cv"]), 3),
             int(r["seed"]), float(r["p_hat_target"])) for r in ok}
    print()
    print(f"N: wave 1 is {len(have & want)} of {len(want)} cells "
          f"({100*len(have & want)/len(want):.0f}%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
