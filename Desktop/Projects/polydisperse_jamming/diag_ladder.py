"""diag_ladder.py — print the raw phi(p_hat) ladders, to check the extrapolation model.

The free-exponent anchor came back b = 0.248 on the first complete ladder against a
predicted b = 1 for a linear contact law. Either the ladder is not in the critical
regime, or the extrapolation model is wrong, and either way phi_J is not trustworthy
until it is understood. This prints the numbers rather than the fit.
"""
import collections
import glob
import json

import numpy as np

LADDER = (1e-3, 3e-4, 1e-4, 3e-5)


def main() -> None:
    recs = {}
    for f in glob.glob("out/phase2*.jsonl"):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if r.get("key") and not r.get("failed"):
                    recs[r["key"]] = r

    g = collections.defaultdict(list)
    for r in recs.values():
        g[(r.get("family", "lognormal"), r["cv"], r["mu_p"], r["seed"])].append(r)

    hdr = f"{'fam':>9} {'cv':>5} {'mu':>5} {'s':>2} | "
    hdr += " ".join(f"{f'phi@{p:.0e}':>14}" for p in LADDER)
    print(hdr)
    print("-" * len(hdr))
    for k in sorted(g):
        fam, cv, mu, seed = k
        d = {round(np.log10(r["p_hat_target"]), 3): r for r in g[k]}
        if len(d) < 3:
            continue
        row = f"{fam[:9]:>9} {cv:5.2f} {mu:5.2f} {seed:2d} | "
        phis, ps = [], []
        for p in LADDER:
            r = d.get(round(np.log10(p), 3))
            if r:
                row += f"{r['phi']:14.4f} "
                phis.append(r["phi"])
                ps.append(r["p_hat"])
            else:
                row += f"{'-':>14} "
        # local slope in log-log of (phi - phi_min) against p, as a crude exponent check
        if len(phis) >= 3:
            phis = np.array(phis)
            ps = np.array(ps)
            dphi = phis - phis.min()
            ok = dphi > 0
            if ok.sum() >= 2:
                sl = np.polyfit(np.log(ps[ok]), np.log(dphi[ok]), 1)[0]
                row += f"  d(phi-phi_min)/dlog p slope ~ {sl:.3f}"
        print(row)

    print()
    print("Interpretation guide:")
    print("  A LINEAR (Hookean) contact law predicts phi - phi_J ~ p^1.")
    print("  If phi barely moves across 1.5 decades of p, the ladder is ALREADY at")
    print("  phi_J to within its own scatter, the free exponent is unidentifiable,")
    print("  and the linear extrapolation is harmless but the exponent check is")
    print("  meaningless -- it is measuring noise, not the contact law.")


if __name__ == "__main__":
    main()
