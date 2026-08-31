"""Sweep phi over BOTH containment eta and airflow provisioning p, and save as it goes.

WHY
    `dc_thermal` models under-provisioning with `under = max(0, 1 - p)`, which makes B
    identical at p = 1 and p = 1.4: blowing extra cold air costs p^3 in fan power and
    buys exactly nothing, so the optimizer can only ever want less air. DESIGN.md §8
    calls this gap (1) and says the fix -- a two-sided phi(p) -- needs a CFD run to
    calibrate. This is that run.

    The p < 1 branch does not need CFD: mass balance already fixes it. If the CRAH
    supplies less air than the rack draws, the shortfall is made up from the hot aisle,
    so phi >= 1 - p exactly, and `dc_thermal` already carries that floor. What is
    missing, and what only CFD can give, is the p > 1 branch: how fast recirculation is
    suppressed when the cold aisle is flooded.

HOW p IS DEFINED, AND WHAT THAT DOES NOT CLAIM
    The cavity has no rack fan -- the rack is a solid partition -- so there is no
    airflow demand inside it to normalise against. p is therefore anchored on the cell's
    own buoyant circulation:

        Q_1(Ra) = max |psi| of the buoyancy-only solution at eta = 0
        p       = flux / Q_1

    p = 1 means "the CRAH pushes as much air through the aisle as the unconfined buoyant
    recirculation loop circulates on its own". That is a DEFINITION, fixed at eta = 0 so
    that the two axes stay independent, and it is not a measurement of any real hall's
    provisioning ratio. What transfers to the Tier-0 model is the SHAPE of phi(p) -- in
    particular whether it falls at all above p = 1, which is the question gap (1) is
    stuck on -- not the location of p = 1 on a real fan curve.

Run:  python run_provisioning.py          (resumable; skips points already in the JSON)
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

from aisle import AisleCell

OUT = "cfd_provisioning.json"
ANCHOR = "cfd_anchor.json"

ETAS = [0.0, 0.2, 0.4, 0.6, 0.8, 0.95]
PS = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
RA = 4.0e3
N = 40
INLET = (0.10, 0.40)    # floor grille under the cold aisle, entirely left of the rack
OUTLET = (0.60, 0.90)   # ceiling grille over the hot aisle


def anchor(Ra: float, n: int) -> float:
    """Q_1: the buoyancy-only circulation strength at eta = 0. Cached, because the
    whole p axis is defined against it and it must not drift between runs."""
    cache = json.load(open(ANCHOR)) if os.path.exists(ANCHOR) else {}
    key = f"Ra={Ra:.6g},N={n}"
    if key not in cache:
        a = AisleCell(N=n, eta=0.0, Ra=Ra)
        its = a.solve(max_iter=60000)
        cache[key] = {"Q1": float(np.abs(a.psi).max()), "iters": int(its),
                      "phi_eta0": a.recirc_fraction()}
        json.dump(cache, open(ANCHOR, "w"), indent=2)
    return float(cache[key]["Q1"])


def main() -> None:
    Q1 = anchor(RA, N)
    print(f"anchor Q_1(Ra={RA:.0e}, N={N}) = {Q1:.4f}")
    results = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(round(r["eta"], 4), round(r["p"], 4)) for r in results}

    for eta in ETAS:
        for p in PS:
            if (round(eta, 4), round(p, 4)) in done:
                continue
            t0 = time.time()
            a = AisleCell(N=N, H_room=1.0, eta=eta, Ra=RA,
                          flux=p * Q1, inlet=INLET, outlet=OUTLET)
            try:
                its = a.solve(max_iter=120000, tol=2e-5)
            except FloatingPointError:
                print(f"  eta={eta} p={p} DIVERGED -- recorded as such, not dropped")
                results.append({"eta": eta, "p": p, "Ra": RA, "N": N, "Q1": Q1,
                                "flux": p * Q1, "phi": None, "diverged": True})
                json.dump(results, open(OUT, "w"), indent=2)
                continue
            rec = {
                "eta": eta, "p": p, "Ra": RA, "N": N, "Q1": Q1, "flux": p * Q1,
                "phi": a.recirc_fraction(),
                "phi_tracer": a.recirc_fraction_tracer(),
                "iters": int(its), "diverged": False,
                "max_psi": float(np.abs(a.psi).max()),
                "seconds": round(time.time() - t0, 1),
            }
            results.append(rec)
            json.dump(results, open(OUT, "w"), indent=2)
            print(f"  eta={eta:<5} p={p:<5} phi={rec['phi']:.4f} "
                  f"tracer={rec['phi_tracer']:.4f} its={its} {rec['seconds']}s")


if __name__ == "__main__":
    main()
