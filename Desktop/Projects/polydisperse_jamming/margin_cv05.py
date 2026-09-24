"""margin_cv05.py -- the closed-form CV claim of PREREG_CVMAP, restated in phi and re-graded.

At a fixed weighed phi, the CV-0.5 suspension can shear-jam iff phi >= phi_m(CV 0.5) (x >= chi means
phi/phi_0 >= phi_m/phi_0), and the CV-0 suspension flows at rest iff phi < phi_0(CV 0). A weighed phi at
which "CV 0 flows but CV 0.5 can jam" exists iff
    D = phi_m(CV 0.5, mu 0.3) - phi_0(CV 0, mu 0) < 0.
phi_0(CV 0.5) CANCELS, so its n = 1 converged seed was never the blocker (predict_cvmap.py propagated it
anyway; that 1.2x figure is superseded here). Seeds are PAIRED across CV and mu: packing.py draws radii
from default_rng(seed) and the lognormal reuses the same normal draws rescaled by sigma. So the error on D
comes from the per-seed difference, floored by the protocol's own re-run scatter, PHI_SCALE*sqrt(2/n)
(two measured phis, n seeds each). Writes out/margin_cv05.txt.
"""
import math
import os

import numpy as np
from scipy import stats

import jam_map as jm
import three_knobs as tk

HERE = os.path.dirname(os.path.abspath(__file__))
PHI_SCALE = 0.0021


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    best = tk.load()

    def seeds(cv, mu):
        return {s: r["phi"] for (m, c, s, p), r in best.items()
                if m == mu and c == cv and p == jm.RUNG and r.get("converged")}

    a = seeds(0.5, 0.3)
    b = seeds(0.0, 0.0)
    common = sorted(set(a) & set(b))
    d = np.array([a[s] - b[s] for s in common])
    n = len(d)
    se_pair = d.std(ddof=1) / math.sqrt(n)
    floor = PHI_SCALE * math.sqrt(2.0 / n)
    se = max(se_pair, floor)
    tcrit = stats.t.ppf(0.975, n - 1)
    say("margin_cv05 -- D = phi_m(CV 0.5, mu 0.3) - phi_0(CV 0, mu 0), rung p_hat 3e-5, converged seeds")
    say(f"  phi_m(0.5, mu .3) per seed: {[(s, round(a[s], 5)) for s in sorted(a)]}")
    say(f"  phi_0(0.0, mu 0)  per seed: {[(s, round(b[s], 5)) for s in sorted(b)]}")
    say(f"  paired seeds {common}: D per seed {[round(v, 5) for v in d]}")
    say(f"  mean D = {d.mean():+.5f};  se (paired) {se_pair:.5f};  protocol floor PHI_SCALE*sqrt(2/n) {floor:.5f}")
    say(f"  t = D/max(se, floor) = {d.mean() / se:+.2f}  against t_crit(df {n - 1}) = {tcrit:.3f}")
    say(f"  unpaired (hypot of per-arm se): t = {d.mean() / math.hypot(np.std(list(a.values()), ddof=1) / math.sqrt(len(a)), np.std(list(b.values()), ddof=1) / math.sqrt(len(b))):+.2f}")
    verdict = ("RESOLVED: no overlap (CV 0.5 cannot jam at any weighed phi where CV 0 flows)" if d.mean() / se > tcrit
               else "NOT RESOLVED at 95 %: the sign of D is the point estimate only")
    say(f"  => {verdict}")
    with open(os.path.join(HERE, "out", "margin_cv05.txt"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
