"""grade_shearkmer.py -- PREREG_SHEARKMER_2026-09-22.md. Writes out/grade_shearkmer.txt.
Written before the data. Reuses actshape's estimators and grade_actshape's dump splitter."""
import glob
import json
import math
import os
import sys

import numpy as np

import actshape as A
import grade_actshape as G

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
T_CRIT = {1: 12.706, 2: 4.303}


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    recs = {}
    for f in glob.glob(os.path.join(OUT, "phase9_shearkmer_*.jsonl")):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            recs[r["key"]] = r
    ref = G.load_records(OUT, "phase9_actshape_c00a*.jsonl")
    say(f"grade_shearkmer -- {len(recs)} of 18 records")
    cells = {}
    for k, r in recs.items():
        if "error" in r:
            say(f"  ERROR {k}: {r['error'][:120]}")
            continue
        w = r.get("doublet_w", 0.0)
        p = os.path.join(OUT, "shearkmer_dump", k + ".npz")
        if not os.path.exists(p):
            say(f"  NO DUMP {k}")
            continue
        d = np.load(p)
        snaps = G.snapshots(d)
        fn_all = np.concatenate([s["fn"] for s in snaps])
        deep = float((fn_all / 0.25 > 0.05).mean())
        healthy = bool(r["steady"] and r["servo_held"] and r["quasistatic"])
        g = np.concatenate([A.g_values(s["fn"], s["r1"], s["r2"], s["p"], "derjaguin")
                            for s in snaps])
        cells[k] = dict(w=w, mu=r["mu_p"], seed=r["seed"], healthy=healthy, deep=deep,
                        sstar=A.sigma_star_over_kappa(g), d=A.delta_from_quantiles(g, .1, .9),
                        rec=r)
    # SK0
    sk0 = []
    for k, c in cells.items():
        if c["w"] == 0.0 and k in ref:
            a, b = c["rec"], ref[k]
            sk0.append(all(a[f] == b[f] for f in ("p_hat_series", "z_series", "force")))
    say(f"SK0 identity vs actshape (w2 = 0): {sum(sk0)}/{len(sk0)}"
        + ("  PASS" if sk0 and all(sk0) and len(sk0) == 6 else "  FAIL/INCOMPLETE"))
    bad = [k for k, c in cells.items() if not c["healthy"]]
    say(f"SK1 health: {len(cells) - len(bad)}/{len(cells)}")
    for k in bad:
        say(f"  EXCLUDED (health) {k}: {cells[k]['rec'].get('notes')}")
    k4 = [k for k, c in cells.items() if c["deep"] > 1e-3]
    say(f"SK2 overlap: {len(k4)} cells with > 0.1 % of contacts beyond u/r 0.05")
    for k in k4:
        say(f"  EXCLUDED (K4) {k}: {cells[k]['deep']:.2%}")
    ok = {k: c for k, c in cells.items() if k not in bad and k not in k4}
    say("")
    say("   mu   w2    n   sigma*/kappa    delta10/90")
    by = {}
    for c in ok.values():
        by.setdefault((c["mu"], c["w"]), {})[c["seed"]] = c
    for (mu, w), v in sorted(by.items()):
        say(f"  {mu:.1f}  {w:.2f}  {len(v)}   {np.mean([x['sstar'] for x in v.values()]):.4f}"
            f"        {np.mean([x['d'] for x in v.values()]):.3f}")
    say("")
    for mu in (0.3, 0.0):
        a, b = by.get((mu, 0.0), {}), by.get((mu, 1.0), {})
        s = sorted(set(a) & set(b))
        tag = "[GRADED]" if mu == 0.3 else "[reported]"
        if len(s) < 2:
            say(f"  mu {mu}: S1/S2 not gradable (n {len(s)}) {tag}")
            continue
        rat = float(np.mean([b[x]["sstar"] / a[x]["sstar"] for x in s]))
        dd = [b[x]["d"] - a[x]["d"] for x in s]
        se = float(np.std(dd, ddof=1) / math.sqrt(len(dd)))
        say(f"  mu {mu}: S1 sigma*(w2=1)/sigma*(0) = {rat:.3f} (bar [0.55, 0.95]) "
            f"{'PASS' if 0.55 <= rat <= 0.95 else 'FAIL'}   {tag}")
        say(f"  mu {mu}: S2 delta(w2=1) - delta(0) = {np.mean(dd):+.3f} +- {se:.3f} "
            f"(bar |.| < 0.05) {'PASS' if abs(np.mean(dd)) < 0.05 else 'FAIL'}   {tag}")
        m = by.get((mu, 0.5), {})
        if m:
            sm = np.mean([x["sstar"] for x in m.values()])
            lo, hi = sorted([np.mean([x["sstar"] for x in a.values()]),
                             np.mean([x["sstar"] for x in b.values()])])
            say(f"  mu {mu}: S3 sigma*(0.5) = {sm:.4f} between [{lo:.4f}, {hi:.4f}] "
                f"{'PASS' if lo <= sm <= hi else 'FAIL'}   {tag}")
    with open(os.path.join(OUT, "grade_shearkmer.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
