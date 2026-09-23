"""net_fixed_x.py -- the doublet effect on the stress to jam at FIXED weighed phi (fixed x),
the comparison referee pass 3 asked for. Writes out/net_fixed_x.txt.

At fixed x, singlets need f_s = (1-x)/(1-chi_s) and pure doublets f_d = (1-x)/(1-chi_d) < f_s.
Net ratio = sigma_jam,d(f_d) / sigma_jam,s(f_s) with sigma_jam/kappa = 1/Q_{1-f}(g) from the
SHEARED mu_p 0.3 per-contact dumps (shearkmer), seed-paired. Under stress-controlled activation.
"""
import glob
import json
import os
import sys

import numpy as np

import actshape as A
import grade_actshape as G

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
CHI = {"mu_p 0.3, p 3e-5 (clean)": (0.9083, 0.8631), "mu_p 0.3, p 1e-3 (clean)": (0.9284, 0.8978),
       "mu_m 3, p 3e-5 (K4-excluded)": (0.8856, 0.7874)}


def gs(w):
    out = {}
    for f in glob.glob(os.path.join(OUT, "phase9_shearkmer_*.jsonl")):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            if "error" in r or r["mu_p"] != 0.3 or r.get("doublet_w", 0.0) != w:
                continue
            d = np.load(os.path.join(OUT, "shearkmer_dump", r["key"] + ".npz"))
            out[r["seed"]] = np.concatenate([A.g_values(s["fn"], s["r1"], s["r2"], s["p"], "derjaguin")
                                             for s in G.snapshots(d)])
    return out


def main():
    s, d = gs(0.0), gs(1.0)
    seeds = sorted(set(s) & set(d))
    L = [f"net_fixed_x -- seeds {seeds}; stress-controlled activation assumed"]
    for lab, (cs, cd) in CHI.items():
        xs = np.linspace(cs + 1e-4, 0.995, 400)
        worst, rows = 0.0, []
        for x in xs:
            fs, fd = (1 - x) / (1 - cs), (1 - x) / (1 - cd)
            r = np.mean([A.sigma_jam_over_kappa(d[k], fd) / A.sigma_jam_over_kappa(s[k], fs)
                         for k in seeds])
            worst = max(worst, r)
            rows.append((x, fs, fd, r))
        L.append(f"{lab}: chi {cs} -> {cd}; f_d/f_s = {(1 - cs) / (1 - cd):.3f}; "
                 f"net ratio over x in [chi_s, 0.995]: worst (largest) {worst:.3f}")
        for x, fs, fd, r in rows[::80]:
            L.append(f"    x {x:.4f}: f_s {fs:.3f} -> f_d {fd:.3f}: net x{r:.3f}")
        L.append(f"    for x in [{cd}, {cs}) singlets cannot jam at all (f_s > 1); doublets can.")
    print("\n".join(L))
    with open(os.path.join(OUT, "net_fixed_x.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
