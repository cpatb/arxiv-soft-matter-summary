"""grade_cv05margin.py -- grades PREREG_CV05MARGIN_2026-09-24.md. Writes out/grade_cv05margin.txt.

Existing records (seeds 1-3) come from three_knobs.load() exactly as margin_cv05.py reads them; the new
seeds (4-9, plus the seed-1 wiring control) are read ONLY from out/cv05margin/phase2_cv05margin_*.jsonl
(the prereg keeps them out of out/ so no other analysis changes).
W1: the seed-1 re-runs reproduce the existing phi of both cells to |dphi| <= 1e-10.
Only converged records count; every non-converged new seed is printed and excluded.
Decision: t = mean D / max(se_paired, PHI_SCALE*sqrt(2/n)) vs t_crit(df n-1), D = phi_m(0.5, mu .3) - phi_0(0, mu 0).
"""
import glob
import json
import math
import os

import numpy as np
from scipy import stats

import jam_map as jm
import three_knobs as tk

HERE = os.path.dirname(os.path.abspath(__file__))
PHI_SCALE = 0.0021   # the protocol's own re-run scatter (mixture_model.py)
W1_TOL = 1e-10


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    best = tk.load()

    def old(cv, mu):
        return {s: r["phi"] for (m, c, s, p), r in best.items()
                if m == mu and c == cv and p == jm.RUNG and r.get("converged")}

    new = {}
    for f in sorted(glob.glob(os.path.join(HERE, "out", "cv05margin", "phase2_cv05margin_*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            new[(r["mu_p"], r["cv"], r["seed"])] = r
    say("grade_cv05margin -- PREREG_CV05MARGIN_2026-09-24.md")
    say(f"new records: {len(new)} (registered: 2 cells x seeds 1,4..9 = 14)")
    a_old, b_old = old(0.5, 0.3), old(0.0, 0.0)
    say(f"existing converged: phi_m(0.5, mu .3) seeds {sorted(a_old)}; phi_0(0, mu 0) seeds {sorted(b_old)}")
    # W1
    w1 = []
    for (mu, cv), ref in (((0.3, 0.5), a_old), ((0.0, 0.0), b_old)):
        r = new.get((mu, cv, 1))
        if r is None or 1 not in ref:
            say(f"  W1 mu {mu} cv {cv}: missing (new {r is not None}, old {1 in ref})")
            w1.append(False)
            continue
        dphi = r["phi"] - ref[1]
        w1.append(abs(dphi) <= W1_TOL)
        say(f"  W1 mu {mu} cv {cv}: new {r['phi']:.12f} old {ref[1]:.12f} dphi {dphi:+.2e} converged {r['converged']}")
    say(f"  W1 verdict: {'PASS' if all(w1) else 'FAIL'}")
    a, b = dict(a_old), dict(b_old)
    for (mu, cv, s), r in sorted(new.items()):
        if s == 1:
            continue
        tag = f"mu {mu} cv {cv} seed {s}: phi {r['phi']:.5f} converged {r['converged']} iters {r['iters']}"
        if not r["converged"]:
            say(f"  EXCLUDED (not converged) {tag}")
            continue
        say(f"  {tag}")
        (a if (mu, cv) == (0.3, 0.5) else b)[s] = r["phi"]
    common = sorted(set(a) & set(b))
    d = np.array([a[s] - b[s] for s in common])
    n = len(d)
    se_pair = d.std(ddof=1) / math.sqrt(n)
    floor = PHI_SCALE * math.sqrt(2.0 / n)
    se = max(se_pair, floor)
    t = d.mean() / se
    tcrit = stats.t.ppf(0.975, n - 1)
    say(f"\npaired seeds {common} (n {n}): D per seed {[round(v, 5) for v in d]}")
    say(f"mean D {d.mean():+.5f}; se paired {se_pair:.5f}; floor {floor:.5f}; t {t:+.2f} vs t_crit(df {n - 1}) {tcrit:.3f}")
    if t > tcrit:
        v = "RESOLVED-NO-OVERLAP: raising CV to 0.5 at fixed weighed phi removes shear-jamming for every suspension that flows at rest at CV 0 (mu_m 0.3)"
    elif t < -tcrit:
        v = "RESOLVED-OVERLAP: a band [phi_m(0.5), phi_0(0)) exists where both hold"
    else:
        v = (f"UNDETERMINED: |D| is below what this protocol resolves; CV 0.5 moves the frictional jamming point to "
             f"within +-{abs(d.mean()) + tcrit * se:.4f} of monodisperse frictionless jamming")
    say(f"VERDICT: {v}")
    with open(os.path.join(HERE, "out", "grade_cv05margin.txt"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
