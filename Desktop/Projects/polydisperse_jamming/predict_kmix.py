"""predict_kmix.py -- blind predictions for PREREG_KMIX_2026-09-23.md (A2: mixed multiplet
populations). Writes out/predict_kmix.txt. Run BEFORE any kmix cell; do not regenerate after.

Generalised ideal specific-volume mixing from the PURE endpoints (kmer singlets and dumbbells,
kchain straight chains), per (mu_p, rung):  1/phi = w1/phi_s + w2/phi_d + w3/phi_c.
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
COMPS = ((0.3, 0.2), (0.2, 0.4), (0.5, 0.5))     # (w2 dumbbells, w3 chains); w1 = rest


def endpoints():
    e = {}
    for pat in ("phase2_kmer_*.jsonl", "phase2_kchain_*.jsonl"):
        for f in glob.glob(os.path.join(OUT, "shape", pat)):
            for line in open(f, encoding="utf-8"):
                r = json.loads(line)
                if r.get("failed") or not r.get("converged"):
                    continue
                k = (r["doublet_w"], r["triplet_w"], bool(r.get("triplet_chain", False)),
                     r["mu_p"], r["p_hat_target"])
                e.setdefault(k, []).append(r["phi"])
    return {k: float(np.mean(v)) for k, v in e.items()}


def main():
    E = endpoints()
    L = ["predict_kmix -- generalised ideal mixing from pure singlet / dumbbell / chain endpoints"]
    for p in (3e-5, 1e-3):
        for mu in (0.0, 0.3):
            s, d, c = (E[(0.0, 0.0, False, mu, p)], E[(1.0, 0.0, False, mu, p)],
                       E[(0.0, 1.0, True, mu, p)])
            L.append(f"p {p:g} mu {mu}: endpoints singlet {s:.4f} dumbbell {d:.4f} chain {c:.4f}")
            for w2, w3 in COMPS:
                w1 = 1 - w2 - w3
                phi = 1.0 / (w1 / s + w2 / d + w3 / c)
                L.append(f"   w2 {w2:.1f} w3 {w3:.1f}: phi_pred {phi:.4f}")
    L.append("")
    L.append("chi_pred(mu 0.3 / mu 0) at p 3e-5 per composition:")
    for w2, w3 in COMPS:
        w1 = 1 - w2 - w3
        ph = {}
        for mu in (0.0, 0.3):
            s, d, c = (E[(0.0, 0.0, False, mu, 3e-5)], E[(1.0, 0.0, False, mu, 3e-5)],
                       E[(0.0, 1.0, True, mu, 3e-5)])
            ph[mu] = 1.0 / (w1 / s + w2 / d + w3 / c)
        L.append(f"   w2 {w2:.1f} w3 {w3:.1f}: chi_pred {ph[0.3] / ph[0.0]:.4f}")
    print("\n".join(L))
    with open(os.path.join(OUT, "predict_kmix.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
