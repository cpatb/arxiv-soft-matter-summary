"""run_eta_curve.py -- the containment CURVE, not its two endpoints.

`run_ra_ladder.py` swept only eta = 0 and 0.95, and says so: "the question is whether the
GAP between them survives at high Ra, not the shape of the curve in between". But nobody
buys containment as a switch -- blanking panels, aisle doors and curtains go in
incrementally -- so what decides the spend is the SHAPE of phi(eta), and two points cannot
tell a straight line from a curve with a knee.

Gates and bars are fixed in ../PREREG_eta_curve_2026-09-08.md, written before this ran.
Every limit `run_ra_ladder.py` records still applies: flux = 0 throughout, Ra reaches 1e6
against 1e9-1e12 in a hall, and above Ra ~ 1e8 a 2-D steady laminar model is the wrong
model rather than a coarse one.

    python run_eta_curve.py
"""
from __future__ import annotations

import json
import math
import os
import time

import numpy as np

from aisle import AisleCell

OUT = "cfd_eta_curve.json"
LADDER = "cfd_ra_ladder.json"

# Nine eta values: the two endpoints the ladder already has (so GD1 has a control) plus
# seven interior points, denser near eta = 1 where the curvature must be if there is any.
ETAS = [0.0, 0.2, 0.4, 0.6, 0.75, 0.85, 0.9, 0.95, 0.98]
RAS = [4e3, 3e4, 1e5, 1e6]          # spans 250x, the same admissible range as the ladder

# Copied from run_ra_ladder.py so this script states what it uses. If they drift there,
# GD1 fails, which is the point of GD1.
CELLS_IN_BL = 3.0
N_MIN, N_MAX = 40, 160
TRANSITION_RA = 1e8
MAX_MINUTES = 12.0
PHI_TOL = 8e-4
CHECK_EVERY = 1000
MIN_STEPS = 6000

BAR_GD1 = 2 * PHI_TOL       # control: reproduce the ladder's endpoints
BAR_GD2_RATIO = 3.0         # proportional model rejected if its RMS is this much worse
BAR_GD3_A_DROP = 0.25       # the leak term may not fall by more than this over the range


def n_required(ra: float) -> int:
    return int(math.ceil(CELLS_IN_BL * ra ** 0.25))


def run_cell(ra: float, eta: float, n_use: int) -> dict:
    a = AisleCell(N=n_use, H_room=1.0, eta=eta, Ra=ra)
    dt = a.stable_dt()          # from the cell, never a copy of its formula
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
    ok = bool(np.isfinite(phi))
    return dict(eta=eta, Ra=ra, N=n_use,
                phi=round(float(phi), 5) if ok else None,
                phi_tracer=round(float(a.recirc_fraction_tracer()), 5) if ok else None,
                steps=n, dt=dt, minutes=round(mins, 2), converged=converged,
                status=("converged" if converged else
                        ("DIVERGED" if not ok else
                         "TIMEOUT, phi is a running value and must not be quoted")))


def main() -> None:
    results = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(r["eta"], r["Ra"]) for r in results}
    print(f"{'Ra':>8} {'N':>4} {'eta':>5} {'phi':>8} {'phi_trc':>9} {'min':>6}  status",
          flush=True)
    for ra in RAS:
        n_use = min(max(n_required(ra), N_MIN), N_MAX)
        if ra >= TRANSITION_RA:
            continue
        for eta in ETAS:
            if (eta, ra) in done:
                continue
            rec = run_cell(ra, eta, n_use)
            results.append(rec)
            json.dump(results, open(OUT, "w"), indent=2)
            print(f"{ra:8.0e} {n_use:4d} {eta:5.2f} "
                  f"{(rec['phi'] if rec['phi'] is not None else float('nan')):8.5f} "
                  f"{(rec['phi_tracer'] if rec['phi_tracer'] is not None else float('nan')):9.5f} "
                  f"{rec['minutes']:6.2f}  {rec['status']}", flush=True)
    print(f"\nwrote -> {OUT}  ({len(results)} cells). Grade with: python grade_eta_curve.py")


if __name__ == "__main__":
    main()
