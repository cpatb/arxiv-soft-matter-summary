"""diag_shape_budget.py -- what a smaller CV endpoint would cost the shape campaign.

The shape grid's N is set by the periodic collider at the WORST cell, which is lognormal
at the top of the CV ladder with alpha = 2: N >= 5379 at CV 0.5, 3482 at CV 0.4, 1638 at
CV 0.3. Cost rises steeply with N, so moving the endpoint down is the only large lever
that does not touch a safety margin or the truncation.

What it costs is the SIGNAL: D = phi_J(CV_hi) - phi_J(CV_lo) shrinks as CV_hi falls,
while the seed scatter does not. This script measures both from the SPHERE records
already on disk -- the same records the project's 4.47x family ratio comes from -- so the
grid decision is made against a number rather than an intuition.

Reads the SHARED out/phase2*.jsonl glob, which is correct here: these are the existing
sphere cells and this script issues no shape verdict.
"""

from __future__ import annotations

import math
import os
from collections import defaultdict

import numpy as np

from analyse_phase2 import P_HAT_FIT_MAX, load_all
from phase2_phiJ import fit_phi_J

HERE = os.path.dirname(os.path.abspath(__file__))
MU = 0.3
N_WANT = 1500


def main() -> int:
    recs = load_all()
    print(f"N: {len(recs)} records from the shared sphere glob")
    recs = [r for r in recs if r.get("converged")]
    print(f"N: -> {len(recs)} converged")
    recs = [r for r in recs if abs(float(r.get("mu_p", -1)) - MU) < 1e-9]
    print(f"N: -> {len(recs)} at mu_p = {MU}")
    recs = [r for r in recs if int(r.get("n_requested", 0)) == N_WANT]
    print(f"N: -> {len(recs)} at N = {N_WANT}")
    recs = [r for r in recs if float(r.get("alpha", 1.0)) == 1.0]
    print(f"N: -> {len(recs)} at alpha = 1 (spheres; records written before the shape "
          f"axis existed carry no `alpha` and default to 1.0)")

    groups = defaultdict(list)
    for r in recs:
        groups[(r.get("family", "lognormal"), round(float(r["cv"]), 3),
                int(r["seed"]))].append(r)

    phij = defaultdict(dict)
    for (fam, cv, seed), rows in groups.items():
        rows = [r for r in rows if float(r["p_hat"]) <= P_HAT_FIT_MAX]
        if len(rows) < 3:
            continue
        p = np.array([float(r["p_hat"]) for r in rows])
        f = np.array([float(r["phi"]) for r in rows])
        v = fit_phi_J(p, f)["phi_J_lin"]
        if np.isfinite(v):
            phij[(fam, cv)][seed] = float(v)
    print(f"N: -> {sum(len(v) for v in phij.values())} (family, CV, seed) ladders fitted")

    fams = sorted({k[0] for k in phij})
    cvs = sorted({k[1] for k in phij})
    print()
    print(f"CV rungs available: {cvs}")
    print()
    print(f"{'family':<12}{'CV':>6}{'n':>4}{'mean phi_J':>13}{'sd':>11}{'se':>11}")
    sd_pool = []
    for fam in fams:
        for cv in cvs:
            d = phij.get((fam, cv), {})
            if len(d) < 2:
                continue
            v = np.array(list(d.values()))
            sd = float(v.std(ddof=1))
            sd_pool.append(sd)
            print(f"{fam:<12}{cv:>6.2f}{v.size:>4}{v.mean():>13.6f}{sd:>11.6f}"
                  f"{sd/math.sqrt(v.size):>11.6f}")

    sd_typ = float(np.median(sd_pool)) if sd_pool else float("nan")
    print()
    print(f"Median per-cell seed sd across every (family, CV): {sd_typ:.6f}")

    print()
    print("D = phi_J(CV_hi) - phi_J(CV_lo = 0.0), and what the endpoint choice costs:")
    print()
    print(f"{'family':<12}{'CV_hi':>7}{'D':>11}{'se(D)':>10}{'|D|/se':>9}"
          f"{'D/D(0.5)':>11}{'N needed':>10}")
    # N required at alpha = 2 for each CV_hi, from PREREG_SHAPE amendment 5.4
    N_NEEDED = {0.1: 256, 0.2: 686, 0.3: 1638, 0.4: 3482, 0.5: 5379}
    ref = {}
    for fam in fams:
        lo = phij.get((fam, 0.0), {})
        if len(lo) < 2:
            continue
        vlo = np.array(list(lo.values()))
        d05 = None
        for cv in cvs:
            if cv == 0.0:
                continue
            hi = phij.get((fam, cv), {})
            if len(hi) < 2:
                continue
            vhi = np.array(list(hi.values()))
            D = float(vhi.mean() - vlo.mean())
            se = float(math.sqrt(vhi.var(ddof=1) / vhi.size + vlo.var(ddof=1) / vlo.size))
            if abs(cv - 0.5) < 1e-9:
                d05 = D
            rel = D / d05 if d05 else float("nan")
            print(f"{fam:<12}{cv:>7.2f}{D:>11.6f}{se:>10.6f}"
                  f"{abs(D)/se if se else float('nan'):>9.2f}{rel:>11.3f}"
                  f"{N_NEEDED.get(round(cv, 1), 0):>10d}")
        ref[fam] = d05

    print()
    print("The cross-family ratio S3 asks about, at each CV endpoint:")
    print()
    print(f"{'CV_hi':>7}{'max D':>11}{'min D':>11}{'ratio':>9}{'N needed':>10}")
    for cv in cvs:
        if cv == 0.0:
            continue
        ds = []
        for fam in fams:
            lo, hi = phij.get((fam, 0.0), {}), phij.get((fam, cv), {})
            if len(lo) < 2 or len(hi) < 2:
                continue
            ds.append(float(np.mean(list(hi.values())) - np.mean(list(lo.values()))))
        if len(ds) < 2:
            continue
        mx, mn = max(ds), min(ds)
        print(f"{cv:>7.2f}{mx:>11.6f}{mn:>11.6f}"
              f"{mx/mn if mn else float('nan'):>9.2f}"
              f"{N_NEEDED.get(round(cv, 1), 0):>10d}")
    print()
    print("Read this as: moving the CV endpoint DOWN buys a large N reduction (cost "
          "falls faster than linearly in N) and sells signal, because D shrinks while "
          "the seed sd does not. The ratio column is the quantity S3's registered "
          "question is about.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
