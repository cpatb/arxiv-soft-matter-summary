"""diag_shape_iters.py -- HOW BIG a cap the CV = 0.5 arm needs, predicted before it is run.

Amendment 5.9 found every CV = 0.5 non-convergence sits AT the iteration cap. The cap
`MAX_ITER = 200000` counts ITERATIONS, and the timestep is set by the SMALLEST grain
(`PWaveTimeStep` reads the minimum radius), so at fixed physical dynamics the number of
iterations a cell needs scales as `1/dt`, i.e. with `d_mean/d_min` -- which is a function
of CV and of the family's truncation, not of the physics being measured.

Hypothesis tested here, on the cells already on disk:

    t_conv = iters * dt   is (approximately) CV-INDEPENDENT at fixed (alpha, p_hat)

If so, the cap a CV = 0.5 cell needs is `t_conv(CV = 0) / dt(CV = 0.5)`, computable from
records that exist, and the diagnostic can be checked against the one CV = 0.5 rung that
DID converge (alpha = 1, p_hat = 1e-3, 6/6) before it is believed for the rungs that did
not. That check is the "reproduce the number you audit" step; the prediction for the
unconverged rungs is only quoted after it passes.

Also printed: wall seconds per iteration, because the re-run has to fit in a session.
"""

from __future__ import annotations

import glob
import json
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PRED_TOL = 0.25   # a prediction of a converged rung's iteration count within 25 %
                  # is "reproduced"; the spread of iters across seeds at CV = 0 is < 2 %
                  # so 25 % is loose enough to be honest about family-to-family drift
                  # and tight enough to catch a wrong scaling (a 1/dt law off by the
                  # 5x d_mean/d_min ratio would miss by 400 %)


def load():
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, "out", "shape", "phase2_shp_*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    r = json.loads(line)
                    if not r.get("failed"):
                        rows.append(r)
    # dedup by key (identical duplicates exist across shard files)
    seen = {}
    for r in rows:
        seen.setdefault(r["key"], r)
    return list(seen.values())


def main() -> int:
    rows = load()
    print(f"N: {len(rows)} distinct production records")
    prod = [r for r in rows if r["p_hat_target"] >= 9e-5]
    print(f"N: -> {len(prod)} on the registered ladder")

    # 1. dt by (family, CV, alpha): the timestep is a property of the PSD's minimum
    print()
    print("1. timestep dt by (family, CV, alpha) -- set by the smallest grain")
    dt = defaultdict(list)
    for r in prod:
        dt[(r["family"], r["cv"], round(r["alpha"], 2))].append(float(r["dt"]))
    for k in sorted(dt):
        v = np.array(dt[k])
        print(f"   {k[0]:<11} CV={k[1]:.1f} alpha={k[2]:.2f}: dt median={np.median(v):.4e} "
              f"spread={v.max()/v.min()-1:.1e}  n={v.size}")
    dt_ref = {k: float(np.median(v)) for k, v in dt.items()}

    # 2. t_conv = iters*dt on CONVERGED cells, by (CV, alpha, p_hat)
    print()
    print("2. physical convergence time t_conv = iters*dt on CONVERGED cells")
    tc = defaultdict(list)
    for r in prod:
        if r["converged"]:
            tc[(r["family"], r["cv"], round(r["alpha"], 2), r["p_hat_target"])].append(
                int(r["iters"]) * float(r["dt"]))
    for k in sorted(tc, key=lambda k: (k[1], k[2], -k[3], k[0])):
        v = np.array(tc[k])
        print(f"   {k[0]:<11} CV={k[1]:.1f} alpha={k[2]:.2f} p_hat={k[3]:.0e}: "
              f"t_conv median={np.median(v):.3f} range=[{v.min():.3f},{v.max():.3f}] n={v.size}")

    # 3. THE CHECK: predict the iteration count of the CV = 0.5 rungs that converged from
    #    the CV = 0 t_conv and the CV = 0.5 dt, and compare with what they took.
    print()
    print("3. CHECK -- predict iters of the CONVERGED CV = 0.5 cells from CV = 0 t_conv "
          "and their own dt")
    ok_all = True
    n_checked = 0
    for r in sorted(prod, key=lambda r: (r["family"], r["alpha"], -r["p_hat_target"], r["seed"])):
        if r["cv"] == 0.0 or not r["converged"]:
            continue
        k0 = ("lognormal", 0.0, round(r["alpha"], 2), r["p_hat_target"])
        if k0 not in tc:
            continue
        t0 = float(np.median(tc[k0]))
        pred = t0 / float(r["dt"])
        got = int(r["iters"])
        err = pred / got - 1
        ok = abs(err) <= PRED_TOL
        ok_all &= ok
        n_checked += 1
        print(f"   {r['family']:<11} alpha={r['alpha']:.2f} p_hat={r['p_hat_target']:.0e} "
              f"s{r['seed']}: predicted {pred:9.0f}  actual {got:7d}  "
              f"({err:+.1%}) {'ok' if ok else 'MISS'}")
    print(f"   N: {n_checked} converged CV = 0.5 cells checked; "
          f"{'ALL within' if ok_all else 'NOT all within'} {PRED_TOL:.0%}")
    if n_checked == 0:
        print("   no converged CV = 0.5 cell to check against -- the prediction below is "
              "UNANCHORED")

    # 4. the prediction for the UNCONVERGED rungs, and the cap it implies
    print()
    print("4. PREDICTED iterations for every CV = 0.5 (family, alpha, p_hat), from CV = 0")
    print(f"   {'family':<11}{'alpha':>6}{'p_hat':>8}{'dt ratio':>10}{'t_conv(0)':>11}"
          f"{'pred iters':>12}{'cap 2e5':>9}")
    worst = 0.0
    preds = {}
    for fam in ("lognormal", "uniform", "bidisperse"):
        for a in (1.0, 2.0):
            k_dt = (fam, 0.5, a)
            if k_dt not in dt_ref:
                print(f"   {fam:<11}{a:>6.2f}   -- no CV = 0.5 record yet, dt unknown; "
                      f"scaling from the family's own d_min is not attempted here")
                continue
            for p in (1e-3, 3e-4, 1e-4):
                k0 = ("lognormal", 0.0, a, p)
                if k0 not in tc:
                    continue
                t0 = float(np.median(tc[k0]))
                ratio = dt_ref[("lognormal", 0.0, a)] / dt_ref[k_dt]
                pred = t0 / dt_ref[k_dt]
                preds[(fam, a, p)] = pred
                worst = max(worst, pred)
                print(f"   {fam:<11}{a:>6.2f}{p:>8.0e}{ratio:>10.2f}{t0:>11.3f}"
                      f"{pred:>12.0f}{'HITS' if pred > 200_000 else 'fits':>9}")
    print(f"   worst predicted: {worst:.0f} iterations")

    # 5. cost: wall seconds per iteration
    print()
    print("5. wall seconds per iteration (N = 5400), by CV -- a re-run must fit a session")
    spi = defaultdict(list)
    for r in prod:
        spi[(r["cv"], round(r["alpha"], 2))].append(float(r["wall_s"]) / int(r["iters"]))
    for k in sorted(spi):
        v = np.array(spi[k])
        print(f"   CV={k[0]:.1f} alpha={k[1]:.2f}: s/iter median={np.median(v):.4f} "
              f"max={v.max():.4f} n={v.size}")
    s_per_iter = float(np.median([x for k, v in spi.items() if k[0] == 0.5 for x in v]))
    print()
    for cap in (400_000, 600_000, 800_000, 1_000_000):
        print(f"   cap {cap:>9d}: worst-case wall {cap*s_per_iter/3600:.1f} h per cell "
              f"at {s_per_iter:.4f} s/iter; covers "
              f"{sum(1 for p in preds.values() if p*1.3 <= cap)}/{len(preds)} predicted "
              f"cells with 30 % margin")
    print()
    print("READ THIS AS: if section 3 reproduces the converged CV = 0.5 rung, the cap "
          "needed is a PSD-geometry number (t_conv / dt), and section 4 says what it is; "
          "the bidisperse family has no dt on disk yet and is not predicted.")
    return 0 if ok_all else 1


if __name__ == "__main__":
    raise SystemExit(main())
