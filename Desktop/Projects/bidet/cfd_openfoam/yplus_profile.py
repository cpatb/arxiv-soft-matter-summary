"""yplus_profile.py -- where on the plate does y+ actually sit, per grid?

WHY THIS EXISTS
---------------
GRID_STUDY_2026-09-03.md explained the coarse family's non-monotonicity by y+ straddling 1
(medians 1.52 / 0.78 / 0.41 across refine 0.5 / 1 / 2), on the grounds that
`omegaWallFunction` blends the viscous and log layers across y+ ~ 1 while the quantity
being converged is the wall shear stress.

A physics referee objected that the MEDIAN hides the radial structure: if y+ is large only
near the stagnation point and small out at r/d = 15 -- which is where tau_w(15mm) lives and
where the non-monotonicity is -- then the mechanism predicts trouble in the wrong place.

That is a checkable objection and this script checks it, from the same raw sample
`fit_closure.py` reads. It reports y+ binned by r/d for every grid, so the mechanism can be
tested at the radius that matters rather than over the plate as a whole.

    python3 yplus_profile.py [case ...]
"""

from __future__ import annotations

import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CASES = ["case_Re10000_H6_r4", "case_Re10000_H6_r2_ext", "case_Re10000_H6",
                 "case_Re10000_H6_r0"]
BINS = [(2, 5), (5, 10), (10, 14), (14, 20), (20, 25)]


def read_info(case: str) -> dict:
    info = {}
    with open(os.path.join(HERE, case, "case_info.txt"), encoding="utf-8") as fh:
        for line in fh:
            if "=" in line:
                k, v = line.strip().split("=", 1)
                info[k] = v
    return info


def latest_raw(case: str, field: str):
    """IMPORT fit_closure's own selector so both read the SAME sample.

    This used to reimplement it, and got it wrong twice in one function: the primary glob
    needed three directory levels below `postProcessing` where the layout has two (so it
    matched nothing on every case and always fell through), and the fallback then took
    `sorted(...)[-1]`, a LEXICOGRAPHIC sort of paths containing the time directory. On this
    tree `sorted(['sampleWSS/11552/y', 'sampleWSS/3000/y'])[-1]` is **3000** -- so the
    moment a pre-restart sample sits beside a converged one, y+ would be read from the
    unconverged state while fit_closure reads the converged one, silently, with this
    module's docstring asserting they are the same sample. fit_closure sorts by
    `float(time_dir)`. Reuse it; do not copy it. (Code review, 2026-09-03.)
    """
    from fit_closure import latest_raw as _fc_latest_raw
    return _fc_latest_raw(os.path.join(HERE, case), field)


def main() -> int:
    cases = sys.argv[1:] or DEFAULT_CASES
    print("y+ on the plate, binned by r/d. The fit band is r/d in [2, 20]; the design")
    print("optimum evaluates the closure at r/d ~ 23-26, past the right edge.\n")
    hdr = "%-24s %8s" % ("case", "median")
    for lo, hi in BINS:
        hdr += " %11s" % ("%d-%d" % (lo, hi))
    print(hdr)
    for case in cases:
        path = latest_raw(case, "yPlus")
        if path is None:
            print("%-24s  (no yPlus sample)" % case)
            continue
        info = read_info(case)
        d = float(info.get("D", 1.0e-3)) if "D" in info else 1.0e-3
        arr = np.loadtxt(path)
        # raw surface format: x y z value...  radial coordinate is hypot(y, z)
        xyz, val = arr[:, :3], arr[:, 3]
        r_over_d = np.hypot(xyz[:, 1], xyz[:, 2]) / d
        row = "%-24s %8.3f" % (case, float(np.median(val)))
        for lo, hi in BINS:
            m = (r_over_d >= lo) & (r_over_d < hi)
            row += " %11s" % ("%.3f" % float(np.median(val[m])) if m.any() else "--")
        print(row)
    print("\nIf y+ is well below 1 in the 14-20 bin on EVERY grid, the wall-function")
    print("blending cannot be the mechanism for non-monotonicity at r/d = 15, and the")
    print("median is the wrong statistic to have argued from.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
