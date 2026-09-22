"""grade_kchain.py -- PREREG_KMER Amendment 1 (straight-chain trimers), with Amendment 2's K4
overlap gate. Writes out/grade_kchain.txt. Compares against the kmer triangle and dumbbell cells.
"""
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")


def load(pattern):
    out = {}
    for f in glob.glob(os.path.join(OUT, "shape", pattern)):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            if not r.get("failed"):
                out[r["key"]] = r
    return out


def k4_ok(key):
    p = os.path.join(OUT, "kmer_dump", key + ".npz")
    if not os.path.exists(p):
        return None
    fn = np.load(p)["fn"]
    return float((fn / 0.25 > 0.05).mean()) <= 1e-3


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    ch = load("phase2_kchain_*.jsonl")
    km = load("phase2_kmer_*.jsonl")
    say(f"grade_kchain -- {len(ch)} of 36 chain records")
    grp = defaultdict(dict)
    excl = 0
    for src, recs in (("chain", ch), ("kmer", km)):
        for k, r in recs.items():
            if not r.get("converged"):
                continue
            ok = k4_ok(k)
            if ok is False:
                excl += 1
                if src == "chain":
                    continue
                label = None
            if src == "chain":
                label = f"chain w3={r['triplet_w']}"
            elif r["doublet_w"] == 1.0 and r["triplet_w"] == 0:
                label = "dumbbell w2=1"
            elif r["triplet_w"] == 1.0:
                label = "triangle w3=1"
            elif r["doublet_w"] == 0 and r["triplet_w"] == 0:
                label = "singlet"
            else:
                continue
            if ok is False:
                label += " [K4-EXCLUDED]"
            grp[(label, r["mu_p"], r["p_hat_target"])][r["seed"]] = r
    say(f"K4-excluded chain cells: see below; total K4 exclusions touched: {excl}")
    for rung in (3e-5, 1e-3):
        say("")
        say(f"=== p_hat {rung:g}: chi = phi(mu m)/phi(mu 0), seed-paired; Z_nr ===")
        labels = sorted({k[0] for k in grp if k[2] == rung})
        for lab in labels:
            base = lab.replace(" [K4-EXCLUDED]", "")
            z0 = grp.get((lab, 0.0, rung), {}) or grp.get((base, 0.0, rung), {})
            line = f"  {lab:30s}"
            for mm in (0.3, 3.0):
                a = grp.get((lab, mm, rung), {})
                s = sorted(set(a) & set(z0))
                if s:
                    c = np.mean([a[x]["phi"] / z0[x]["phi"] for x in s])
                    line += f"  chi({mm:g}) {c:.4f} (n {len(s)})"
            for mm in (0.0, 3.0):
                a = grp.get((lab, mm, rung), {})
                if a:
                    line += f"  Z({mm:g}) {np.mean([x['z_nr'] for x in a.values()]):.3f}"
            say(line)
    say("")
    say("C1 (chi(chain) < 0.788 < 0.886 and chi(triangle) > chi(chain) at w3 = 1, p 3e-5, mu_m 3):")
    say("   read from the table; the mu_m 3 rows are K4-EXCLUDED wherever flagged, so C1 is graded")
    say("   only if the chain mu_m 3 cells pass K4 -- otherwise UNDETERMINED, with mu_m 0.3 printed.")
    with open(os.path.join(OUT, "grade_kchain.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
