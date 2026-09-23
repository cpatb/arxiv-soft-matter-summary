"""shearkmer_sjam.py -- does a doublet fraction LOWER the stress to jam across the whole window,
once the sheared activation softening is included? Writes out/shearkmer_sjam.txt. Reported, not
a registered gate (PREREG_SHEARKMER S2 FAILED: delta falls under shear, so this is the question
the failure raises). Fixed chemistry (kappa) and fixed primary size; sigma_jam/kappa =
1/Q_{1-f_c}(g), no Weibull assumption (derivations/jam_stress_2026-09-22.md eq. 8)."""
import glob
import json
import os
import sys

import numpy as np

import actshape as A
import grade_actshape as G

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FCS = (0.10, 0.25, 0.37, 0.50, 0.75, 0.90)


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    g = {}
    for f in glob.glob(os.path.join(OUT, "phase9_shearkmer_*.jsonl")):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            if "error" in r or not (r["steady"] and r["servo_held"] and r["quasistatic"]):
                continue
            d = np.load(os.path.join(OUT, "shearkmer_dump", r["key"] + ".npz"))
            gg = np.concatenate([A.g_values(s["fn"], s["r1"], s["r2"], s["p"], "derjaguin")
                                 for s in G.snapshots(d)])
            g.setdefault((r["mu_p"], r.get("doublet_w", 0.0)), {})[r["seed"]] = gg
    say("shearkmer_sjam -- sigma_jam(w2)/sigma_jam(0) at fixed kappa and primary, seed-paired")
    say("   (sigma = p; f_c is the frictional fraction the window position requires)")
    for mu in (0.3, 0.0):
        say(f"  mu_p {mu}:  f_c = " + "  ".join(f"{fc:.2f}" for fc in FCS))
        for w in (0.5, 1.0):
            a, b = g.get((mu, 0.0), {}), g.get((mu, w), {})
            s = sorted(set(a) & set(b))
            rows = []
            for fc in FCS:
                r = [A.sigma_jam_over_kappa(b[x], fc) / A.sigma_jam_over_kappa(a[x], fc) for x in s]
                rows.append(f"{np.mean(r):.3f}")
            say(f"    w2 {w:.1f} (n {len(s)}):  " + "  ".join(rows))
    say("")
    say("  < 1 means doublets lower the stress needed to jam at that position in the window.")
    with open(os.path.join(OUT, "shearkmer_sjam.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
