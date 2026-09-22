"""kmer_overlap.py -- how much does inter-grain interpenetration inflate phi in the kmer packs?
Writes out/kmer_overlap.txt. A SENSITIVITY for PREREG_KMER Amendment 2 (gate K4), not a gate.

phi is the analytic union solid volume over the cell volume, which assumes grains do not
interpenetrate. Every contact overlap u between two primaries of radius r double-counts the lens
V_lens = pi (6r - u) u^2 / 12 (two equal spheres, centre distance 2r - u). With the linear
FrictPhys law between equal primaries u = fn/(E r) = fn/0.5 (E = 1, r = 1/2). So

    phi_corr = phi * (1 - sum V_lens / vol_solid),   vol_solid = N_primary (4/3) pi r^3.

Checked: at the small overlaps of a sphere pack the correction is ~1e-4 (the protocol's own
softness), so it measures only the anomalous deep overlaps of the frictional clump packs.
"""
import glob
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
R = 0.5


def lens(u):
    u = np.clip(u, 0.0, 2 * R)
    return math.pi * (6 * R - u) * u ** 2 / 12.0


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    rows = defaultdict(list)
    for f in glob.glob(os.path.join(OUT, "shape", "phase2_kmer_*.jsonl")):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            if r.get("failed") or not r.get("converged"):
                continue
            p = os.path.join(OUT, "kmer_dump", r["key"] + ".npz")
            if not os.path.exists(p):
                continue
            fn = np.load(p)["fn"]
            u = fn / (1.0 * R)                       # kn = E r
            dv = float(np.sum(lens(u)))
            corr = 1.0 - dv / r["vol_solid"]
            rows[(r["doublet_w"], r["triplet_w"], r["mu_p"], r["p_hat_target"])].append(
                (r["seed"], r["phi"], r["phi"] * corr, 1.0 - corr))
    say("kmer_overlap -- phi inflation by inter-grain interpenetration (lens double-count)")
    say("  w2    w3    mu   rung    n   phi      phi_corr  mean inflation (fraction of solid)")
    for k in sorted(rows):
        v = rows[k]
        say(f"  {k[0]:.2f}  {k[1]:.2f}  {k[2]:.1f}  {k[3]:.0e}  {len(v)}  "
            f"{np.mean([x[1] for x in v]):.4f}   {np.mean([x[2] for x in v]):.4f}    "
            f"{np.mean([x[3] for x in v]):.2e}")
    say("")
    say("chi = phi(mu 3)/phi(mu 0), seed-paired, RAW vs CORRECTED, p_hat 3e-5 and 1e-3:")
    for rung in (3e-5, 1e-3):
        for w in (0.0, 0.25, 0.5, 0.75, 1.0):
            a = {x[0]: x for x in rows.get((w, 0.0, 3.0, rung), [])}
            b = {x[0]: x for x in rows.get((w, 0.0, 0.0, rung), [])}
            s = sorted(set(a) & set(b))
            if not s:
                continue
            raw = np.mean([a[x][1] / b[x][1] for x in s])
            cor = np.mean([a[x][2] / b[x][2] for x in s])
            say(f"  p {rung:.0e} w2 {w:.2f}: chi raw {raw:.4f}  corrected {cor:.4f}  "
                f"(shift {cor - raw:+.4f}, n {len(s)})")
    with open(os.path.join(OUT, "kmer_overlap.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
