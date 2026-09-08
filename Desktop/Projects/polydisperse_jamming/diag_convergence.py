"""diag_convergence.py — is the non-converged set a BIASED sample?

20% of Phase-2 cells failed the unbalanced-force criterion and were dropped. If the
failures correlate with CV -- the axis the whole project measures -- then the surviving
cells are a biased sample and the phi_J(CV) trend is partly a survival artefact.

Prints the non-convergence rate resolved by every axis, so the answer is visible rather
than assumed.
"""
import collections
import glob
import json

import numpy as np


# System size is fixed here on purpose. out/ holds N = 1500 and N = 3000 cells; this
# diagnostic asks about the CONVERGENCE / LADDER structure at one size, and pooling two
# would make a system-size difference read as ladder scatter. See N_GROUPING_2026-09-08.md.
from analyse_phase2 import PRIMARY_N  # noqa: E402


def load():
    by_key = {}
    for f in glob.glob("out/phase2*.jsonl"):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if (r.get("key") and not r.get("failed")
                        and r.get("n_requested") == PRIMARY_N):
                    by_key[r["key"]] = r
    return list(by_key.values())


def rate_by(recs, keyfn, label):
    d = collections.defaultdict(lambda: [0, 0])
    for r in recs:
        k = keyfn(r)
        d[k][0] += 1
        if not r.get("converged"):
            d[k][1] += 1
    print(f"\n  non-convergence by {label}:")
    print(f"    {'value':>12} {'n':>5} {'failed':>7} {'rate':>7}")
    for k in sorted(d, key=lambda x: (str(type(x)), x)):
        n, bad = d[k]
        print(f"    {str(k):>12} {n:5d} {bad:7d} {bad/n:7.1%}")


def main():
    recs = load()
    conv = sum(1 for r in recs if r.get("converged"))
    print(f"N = {len(recs)} cells, {conv} converged, {len(recs)-conv} not "
          f"({(len(recs)-conv)/len(recs):.1%})")

    rate_by(recs, lambda r: r.get("family", "lognormal"), "family")
    rate_by(recs, lambda r: r["cv"], "CV")
    rate_by(recs, lambda r: r["mu_p"], "mu_p")
    rate_by(recs, lambda r: f"{r['p_hat_target']:.0e}", "p_hat target")

    # does dropping them move phi_J? compare converged-only vs all, per (cv, mu)
    print("\n  phi at p_hat=1e-4 with and without the non-converged cells:")
    print(f"    {'cv':>5} {'mu':>5} {'n_all':>6} {'n_conv':>7} "
          f"{'phi_all':>9} {'phi_conv':>9} {'shift':>9}")
    g = collections.defaultdict(list)
    for r in recs:
        if abs(np.log10(r["p_hat_target"]) - np.log10(1e-4)) < 0.01:
            g[(r["cv"], r["mu_p"])].append(r)
    for k in sorted(g):
        all_r = g[k]
        cv_r = [r for r in all_r if r.get("converged")]
        if not cv_r:
            continue
        pa = np.mean([r["phi"] for r in all_r])
        pc = np.mean([r["phi"] for r in cv_r])
        print(f"    {k[0]:5.2f} {k[1]:5.2f} {len(all_r):6d} {len(cv_r):7d} "
              f"{pa:9.4f} {pc:9.4f} {pc-pa:+9.4f}")

    # unbalanced force distribution of the failures
    bad = [r for r in recs if not r.get("converged")]
    if bad:
        u = np.array([r["unbalanced"] for r in bad])
        it = np.array([r["iters"] for r in bad])
        print(f"\n  failures: unbalanced median {np.median(u):.3f} "
              f"[{u.min():.3f}, {u.max():.3f}] against a 0.02 bar")
        print(f"            iters median {np.median(it):.0f}, "
              f"{(it >= 200000).sum()}/{len(it)} hit the 200k cap")


if __name__ == "__main__":
    main()
