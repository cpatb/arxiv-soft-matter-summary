"""check_alpha_coverage.py -- how much of a flight actually leaves the coefficient table?

WHY
---
P1_AnalyticalModel/README.md states, as a known limitation, that "the tabulated curves
span -10..+30 deg; a realistic flight reaches ~33 deg in its final second, so ~7% of
samples are clamped". The manuscript's hammer discussion leaned on the same idea and got
it badly wrong (see wp4b_item6_retract.py), so the number is worth checking rather than
inheriting: the nominal backhand puts 17.6% of its samples outside the band, not 7%.

That discrepancy has an innocent explanation and a guilty one, and only a Monte Carlo
can tell them apart: "~7%" may be an ENSEMBLE average over the release distributions
(where many draws never reach +30 deg at all) rather than the NOMINAL throw. This script
reports both, per throw type, so the README can state which it means.

Seeded, so the numbers regenerate.

    python check_alpha_coverage.py [n_samples]
"""

from __future__ import annotations

import sys

import numpy as np

import disc_6dof as D
import trajectories as T

# Tabulated support of the Potts & Crowther coefficient curves.
ALPHA_LO_DEG, ALPHA_HI_DEG = -10.0, 30.0
SEED = 20260904
# The nominal constructors and the Monte Carlo distributions use DIFFERENT names for the
# same throw (right_handed_backhand vs right_backhand), so pair them explicitly.
NOMINAL = {"right_backhand": T.right_handed_backhand,
           "right_forehand": T.right_handed_forehand,
           "hammer": T.hammer}
HEADING_DEG = 0.0


def outside_fraction(res) -> tuple:
    n = res.landing_idx + 1
    a = np.degrees(np.array([s["alpha"] for s in res.diag[:n]]))
    out = float(np.mean((a < ALPHA_LO_DEG) | (a > ALPHA_HI_DEG)))
    return out, float(a.min()), float(a.max())


def main() -> int:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    rng = np.random.default_rng(SEED)

    print("Fraction of flight samples outside the tabulated alpha band "
          "[%.0f, %.0f] deg" % (ALPHA_LO_DEG, ALPHA_HI_DEG))
    print("seed = %d, n = %d per throw type\n" % (SEED, n))
    print("%-24s %10s %10s   %26s   %s"
          % ("throw", "NOMINAL", "ensemble", "ensemble alpha_max", "draws never leaving"))
    print("%-24s %10s %10s   %26s   %s"
          % ("", "outside", "mean out", "mean [min, max]", "the table"))

    for name, ctor in NOMINAL.items():
        nom_out, _, nom_max = outside_fraction(
            D.simulate(*ctor(), compute_diag=True))

        outs, maxes = [], []
        for _ in range(n):
            res = D.simulate(*T.sample_release(name, HEADING_DEG, rng=rng),
                             compute_diag=True)
            o, _, mx = outside_fraction(res)
            outs.append(o)
            maxes.append(mx)
        outs = np.array(outs)
        maxes = np.array(maxes)
        clean = float(np.mean(outs == 0.0))
        print("%-24s %9.1f%% %9.1f%%   %8.1f [%6.1f, %6.1f]   %11.0f%%"
              % (name, 100 * nom_out, 100 * outs.mean(),
                 maxes.mean(), maxes.min(), maxes.max(), 100 * clean))

    print()
    print("READING")
    print("-------")
    print("If the ensemble mean for the backhand lands near 7%% while its NOMINAL value")
    print("is ~18%%, the README's number is an ensemble figure and should say so.")
    print("If both are far from 7%%, the README's number does not reproduce and the")
    print("limitation needs restating from these numbers.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
