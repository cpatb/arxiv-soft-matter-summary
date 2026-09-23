"""grade_kmix.py -- PREREG_KMIX_2026-09-23.md (A1 kmu, A2 kmix). Writes out/grade_kmix.txt.
Written before the data. K4 (overlap gate of PREREG_KMER Amendment 2) applied throughout."""
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
PHI_SCALE = 0.0021


def load(pat):
    out = []
    for f in glob.glob(os.path.join(OUT, "shape", pat)):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            out.append(r)
    return out


def deep(key):
    p = os.path.join(OUT, "kmer_dump", key + ".npz")
    if not os.path.exists(p):
        return None
    fn = np.load(p)["fn"]
    return float((fn / 0.25 > 0.05).mean()) if fn.size else 0.0


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    def key3(r):
        return (r["doublet_w"], r["triplet_w"], bool(r.get("triplet_chain", False)))

    # ---- endpoints / mu 0 references from kmer + kchain ---------------------------------
    ref = defaultdict(dict)
    for r in load("phase2_kmer_*.jsonl") + load("phase2_kchain_*.jsonl"):
        if r.get("failed") or not r.get("converged"):
            continue
        ref[(key3(r), r["mu_p"], r["p_hat_target"])][r["seed"]] = r["phi"]

    # ---- A1 kmu --------------------------------------------------------------------------
    kmu = [r for r in load("phase2_kmu_*.jsonl")]
    say(f"=== A1 kmu: {len(kmu)} of 36 records ===")
    g = defaultdict(dict)
    for r in kmu:
        if r.get("failed") or not r.get("converged"):
            say(f"  UNHEALTHY {r['key']}: {r.get('error') or r.get('notes')}")
            continue
        g[(r["doublet_w"], r["mu_p"], r["p_hat_target"])][r["seed"]] = (r["phi"], deep(r["key"]))
    say("   mu    w2   rung    n   phi      deep-overlap frac (max over seeds)  K4")
    clean = {}
    for k in sorted(g):
        v = g[k]
        # a MISSING dump is not evidence of a clean pack: it makes the group unclean (stricter;
        # added before the data landed)
        miss = sum(x[1] is None for x in v.values())
        dmax = max((x[1] for x in v.values() if x[1] is not None), default=float("nan"))
        clean[k] = miss == 0 and dmax <= 1e-3
        if miss:
            say(f"  MISSING DUMP x{miss} in group {k} -> not clean")
        say(f"  {k[1]:.1f}  {k[0]:.2f}  {k[2]:.0e}  {len(v)}   {np.mean([x[0] for x in v.values()]):.4f}"
            f"   {dmax:.4f}   {'clean' if clean[k] else 'EXCLUDED'}")
    u1 = [clean.get((1.0, 1.0, p)) for p in (1e-3, 3e-5)]
    say(f"U1 (mu 1, w2 1): {'CLEAN' if all(u is True for u in u1) else 'DIRTY' if any(u is False for u in u1) else 'INCOMPLETE'}"
        f"  (per rung {u1})")
    for mm in (1.0, 0.6):
        for p in (3e-5, 1e-3):
            line = f"  chi(mu_m {mm}) p {p:.0e}:"
            chis = {}
            for w in (0.0, 0.5, 1.0):
                a = g.get((w, mm, p), {})
                b = ref.get(((w, 0.0, False), 0.0, p), {})
                s = sorted(set(a) & set(b))
                if s and clean.get((w, mm, p)):
                    chis[w] = float(np.mean([a[x][0] / b[x] for x in s]))
                    line += f"  w2 {w}: {chis[w]:.4f}"
                elif s:
                    line += f"  w2 {w}: [K4-excluded]"
            say(line)
            if mm == 1.0 and p == 3e-5 and 0.0 in chis and 1.0 in chis:
                r_ = chis[1.0] / chis[0.0]
                say(f"U2 chi_mu1(1)/chi_mu1(0) = {r_:.4f} (registered 0.889 +- 0.02) "
                    f"{'PASS' if abs(r_ - 0.889) <= 0.02 else 'FAIL'}")

    # ---- A2 kmix -------------------------------------------------------------------------
    kmix = load("phase2_kmix_*.jsonl")
    say("")
    say(f"=== A2 kmix: {len(kmix)} of 36 records ===")
    m = defaultdict(dict)
    for r in kmix:
        if r.get("failed") or not r.get("converged"):
            say(f"  UNHEALTHY {r['key']}")
            continue
        dd = deep(r["key"])
        if dd is None:
            say(f"  NO DUMP {r['key'][:60]} -> excluded (K4 cannot be applied)")
            continue
        if dd > 1e-3:
            say(f"  K4-EXCLUDED {r['key'][:60]}: {dd:.2%}")
            continue
        m[(r["doublet_w"], r["triplet_w"], r["mu_p"], r["p_hat_target"])][r["seed"]] = r["phi"]
    m1 = []
    for (w2, w3, mu, p), v in sorted(m.items()):
        s_, d_, c_ = (np.mean(list(ref[((0.0, 0.0, False), mu, p)].values())),
                      np.mean(list(ref[((1.0, 0.0, False), mu, p)].values())),
                      np.mean(list(ref[((0.0, 1.0, True), mu, p)].values())))
        pred = 1.0 / ((1 - w2 - w3) / s_ + w2 / d_ + w3 / c_)
        obs = float(np.mean(list(v.values())))
        m1.append(abs(obs - pred) < 2 * PHI_SCALE)
        say(f"  w2 {w2} w3 {w3} mu {mu} p {p:.0e}: phi {obs:.4f} vs ideal {pred:.4f} "
            f"({obs - pred:+.4f}, n {len(v)})")
    say(f"M1: {sum(m1)}/{len(m1)} within 2 PHI_SCALE" + ("  PASS" if len(m1) == 12 and all(m1)
                                                        else "  FAIL" if m1 and not all(m1)
                                                        else "  INCOMPLETE"))
    reg = {(0.3, 0.2): 0.8847, (0.2, 0.4): 0.8801, (0.5, 0.5): 0.8635}
    m2 = []
    for (w2, w3), want in reg.items():
        a, b = m.get((w2, w3, 0.3, 3e-5), {}), m.get((w2, w3, 0.0, 3e-5), {})
        s = sorted(set(a) & set(b))
        if not s:
            say(f"  M2 w2 {w2} w3 {w3}: no paired data")
            continue
        chi = float(np.mean([a[x] / b[x] for x in s]))
        m2.append(abs(chi - want) <= 0.01)
        say(f"  M2 w2 {w2} w3 {w3}: chi {chi:.4f} vs {want:.4f} {'ok' if abs(chi - want) <= 0.01 else 'MISS'}")
    say(f"M2: {sum(m2)}/{len(m2)}" + ("  PASS" if len(m2) == 3 and all(m2) else
                                    "  FAIL" if m2 and not all(m2) else "  INCOMPLETE"))
    with open(os.path.join(OUT, "grade_kmix.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
