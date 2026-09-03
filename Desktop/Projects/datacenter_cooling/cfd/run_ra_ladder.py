"""run_ra_ladder.py -- how far up in Rayleigh number can this solver actually go?

WHY THIS EXISTS
---------------
The Status card's item 1 is "re-measure the `eta` axis at hall-like Ra". `phi(eta)` was
calibrated at Ra = 4e3 and 8e3; a real hall is Ra ~ 1e9 to 1e12, and the 31 AUG Lewis
sweep disqualified the 7 % recirculation floor as thermal diffusion without replacing it,
so the ADVECTIVE `phi(eta)` at high Ra is the open number the containment ROI depends on.

This script does not pretend to reach a hall. It measures where this solver stops being
allowed to answer, which is a prerequisite for any statement about high Ra and is
currently written down nowhere.

THE ADMISSIBILITY CRITERION, fixed before the runs
--------------------------------------------------
In laminar Rayleigh-Benard-like convection the thermal boundary layer thins as

    delta / H  ~  Ra^(-1/4)

so resolving it needs at least `CELLS_IN_BL` cells across it:

    N  >=  CELLS_IN_BL * Ra^(1/4)

At the calibration's N = 40 with 3 cells demanded, that caps Ra at about 3e4 -- so the
existing points at 4e3 and 8e3 are inside it and anything much above 1e5 at N = 40 would
NOT have been, which is the thing worth knowing. Each rung here therefore picks its own N
from the criterion rather than inheriting 40, and every record carries `n_bl`, the number
of cells the criterion says are in the layer, so a reader can re-grade the whole ladder
against a different bar without re-running it.

Two further limits are recorded rather than hidden:

  * `TRANSITION_RA` -- above roughly 1e8-1e9 the flow is not laminar and a 2-D steady
    laminar solver is the wrong model no matter how fine the grid. Resolution cannot fix
    a missing physics. Rungs above it are marked inadmissible for that reason and are not
    run.
  * Wall clock. dt is set by the advective limit h/sqrt(Ra*Pr) once Ra is large, and the
    physical time to a steady state does not shrink as fast, so cost grows roughly as
    N^2 * N * sqrt(Ra). The ladder stops at `MAX_MINUTES` per cell and says so in the
    record instead of silently reporting an unconverged phi.

WHAT IS SWEPT
-------------
Only the two containment endpoints that the ROI actually needs -- eta = 0 (no containment)
and eta = 0.95 (near-full) -- because the question is whether the GAP between them survives
at high Ra, not the shape of the curve in between. A trend in the gap is the deliverable.

    python run_ra_ladder.py
"""
from __future__ import annotations

import json
import math
import os
import time

import numpy as np

from aisle import AisleCell

OUT = "cfd_ra_ladder.json"

RAS = [4e3, 8e3, 3e4, 1e5, 3e5, 1e6]
ETAS = [0.0, 0.95]
CELLS_IN_BL = 3.0          # cells demanded across the thermal boundary layer
N_MIN, N_MAX = 40, 160     # floor keeps the coarse rungs comparable to the calibration
TRANSITION_RA = 1e8        # above this a 2-D laminar model is the wrong physics, not a
                           # coarse one; recorded so the wall is attributed correctly
MAX_MINUTES = 12.0
PHI_TOL = 8e-4             # same convergence tolerance as run_calibration.py
CHECK_EVERY = 1000
MIN_STEPS = 6000


def n_required(ra: float) -> int:
    return int(math.ceil(CELLS_IN_BL * ra**0.25))


def main() -> None:
    results = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(r["eta"], r["Ra"]) for r in results}

    print(f"{'Ra':>8} {'N_req':>6} {'N_used':>7} {'eta':>5} "
          f"{'phi':>7} {'phi_trc':>8} {'n_bl':>5} {'steps':>8} {'min':>6}  status")
    for ra in RAS:
        n_req = n_required(ra)
        n_use = min(max(n_req, N_MIN), N_MAX)
        # cells actually across the layer at the N we can afford
        n_bl = n_use * ra**-0.25
        admissible = (n_bl >= CELLS_IN_BL) and (ra < TRANSITION_RA)
        for eta in ETAS:
            if (eta, ra) in done:
                continue
            if not admissible:
                why = ("laminar 2-D model invalid above Ra ~ %.0e" % TRANSITION_RA
                       if ra >= TRANSITION_RA
                       else "needs N >= %d, capped at %d -> %.2f cells in the BL"
                            % (n_req, N_MAX, n_bl))
                rec = dict(eta=eta, Ra=ra, N=n_use, n_bl=round(n_bl, 2),
                           admissible=False, reason=why, phi=None)
                results.append(rec)
                json.dump(results, open(OUT, "w"), indent=2)
                print(f"{ra:8.0e} {n_req:6d} {n_use:7d} {eta:5.2f} "
                      f"{'--':>7} {'--':>8} {n_bl:5.2f} {'--':>8} {'--':>6}  SKIP: {why}")
                continue

            a = AisleCell(N=n_use, H_room=1.0, eta=eta, Ra=ra)
            h = min(a.dx, a.dy)
            dt = 0.35 * min(h * h / (4 * max(a.Pr, 1.0)),
                            h / (np.sqrt(a.Ra * a.Pr) + 1.0))
            t0 = time.time()
            prev, phi, converged, n = 0.0, 0.0, False, 0
            while True:
                a.step(dt)
                n += 1
                if n % CHECK_EVERY == 0:
                    phi = a.recirc_fraction()
                    if not np.isfinite(phi):
                        break
                    if n >= MIN_STEPS and abs(phi - prev) < PHI_TOL:
                        converged = True
                        break
                    prev = phi
                    if (time.time() - t0) / 60.0 > MAX_MINUTES:
                        break
            mins = (time.time() - t0) / 60.0
            status = ("converged" if converged
                      else ("DIVERGED" if not np.isfinite(phi) else "TIMEOUT, phi is a"
                            " running value and must not be quoted"))
            rec = dict(eta=eta, Ra=ra, N=n_use, n_bl=round(n_bl, 2), admissible=True,
                       phi=None if not np.isfinite(phi) else round(float(phi), 4),
                       phi_tracer=None if not np.isfinite(phi)
                       else round(float(a.recirc_fraction_tracer()), 4),
                       steps=n, dt=dt, minutes=round(mins, 2), converged=converged,
                       status=status)
            results.append(rec)
            json.dump(results, open(OUT, "w"), indent=2)
            pt = rec["phi_tracer"]
            print(f"{ra:8.0e} {n_req:6d} {n_use:7d} {eta:5.2f} "
                  f"{(rec['phi'] if rec['phi'] is not None else float('nan')):7.4f} "
                  f"{(pt if pt is not None else float('nan')):8.4f} "
                  f"{n_bl:5.2f} {n:8d} {mins:6.2f}  {status}", flush=True)

    # the deliverable: does the containment GAP survive as Ra rises?
    print()
    print("phi(eta=0) - phi(eta=0.95), the quantity the containment ROI is built on:")
    by_ra = {}
    for r in results:
        if r.get("phi") is not None and r.get("converged"):
            by_ra.setdefault(r["Ra"], {})[r["eta"]] = r["phi"]
    for ra in sorted(by_ra):
        d = by_ra[ra]
        if 0.0 in d and 0.95 in d:
            print(f"  Ra = {ra:8.0e}   phi(0) = {d[0.0]:.4f}   phi(0.95) = {d[0.95]:.4f}"
                  f"   gap = {d[0.0] - d[0.95]:+.4f}")
    print("\nwrote ->", OUT)


if __name__ == "__main__":
    main()
