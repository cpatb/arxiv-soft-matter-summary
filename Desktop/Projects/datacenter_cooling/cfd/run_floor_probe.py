"""Is the "~7% recirculation floor" recirculation, or is it conduction?

WHAT PROMPTED THIS
    The Tier-1 calibration's headline is that phi does NOT vanish at perfect
    containment -- `CFD_RECIRC` gives phi(eta=1) = 0.0715 -- and DESIGN.md calls the
    resulting flattening of the containment ROI "the real technical content" of the
    CFD step. The provisioning sweep run tonight put a number beside it that the
    eta-only sweep never had to face:

        eta = 0.95, p = 0:   phi (from T) = 0.0729   phi (from tracer) = 0.0002

    A factor of 360 between two channels that are supposed to measure the same thing.

TWO CHANNELS, AND WHY THEY CAN DISAGREE
    `recirc_fraction` reads the normalised temperature on the rack-inlet face.
    `recirc_fraction_tracer` reads a passive scalar released at the hot wall whose
    diffusivity is 0.02 -- one fiftieth of the thermal one. Heat therefore reaches the
    inlet by two routes, advection and conduction; the tracer effectively only by
    advection. Where the two agree the signal is advected air. Where they do not, the
    temperature channel is reading conduction across the cavity.

    That would matter, because the cell runs at Ra = 4e3. A real data hall is Ra ~ 1e9
    to 1e12. Conduction is over-represented here by orders of magnitude, and it is
    exactly the channel that survives when containment shuts the advective one off --
    which is to say, it is exactly the channel that would manufacture a floor.

THE TEST, pre-registered before the run
    Sweep Ra at fixed geometry and read both channels.

    H_cond  the floor is conduction. Then phi_T(eta = 0.95) FALLS as Ra rises, roughly
            as the conductive fraction does, while phi_tracer stays near zero and the
            two converge. The eta = 0 point, which is advection-dominated, moves much
            less.
    H_recirc the floor is a real advective short circuit over the top of the partition.
            Then phi_T(eta = 0.95) is roughly Ra-independent and phi_tracer must rise
            to meet it as the flow strengthens.

    Decision rule, fixed here: call it conduction if phi_T at eta = 0.95, p = 0 falls
    by more than 2x from Ra = 4e3 to Ra = 1e5 while phi_tracer stays below 0.01. Call
    it advection if phi_T changes by less than 30% over that range. Anything else is
    undecided and gets reported as undecided.

    A grid check runs at the same time: every case is solved at N = 40 and N = 60. A
    conduction signal is a boundary-layer quantity and the coarse grid does not resolve
    the wall layers, so a floor that moves strongly with N is another way of saying the
    same thing.

Run:  python run_floor_probe.py
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

from aisle import AisleCell
from run_provisioning import solve_to_phi

OUT = "cfd_floor_probe.json"

ETAS = [0.0, 0.4, 0.95]
RAS = [4.0e3, 1.0e4, 4.0e4, 1.0e5]
PS = [0.0, 1.0]
NS = [40, 60]
INLET = (0.10, 0.40)
OUTLET = (0.60, 0.90)


def anchor_for(Ra: float, n: int, cache: dict) -> float:
    """Q_1(Ra, N): the buoyancy-only circulation strength at eta = 0.

    Recomputed per Ra, because p is defined against it and the buoyant velocity scale
    grows like sqrt(Ra). Holding Q_1 fixed while Ra moved would sweep two things.
    """
    key = f"Ra={Ra:.6g},N={n}"
    if key not in cache:
        a = AisleCell(N=n, eta=0.0, Ra=Ra)
        solve_to_phi(a)
        cache[key] = float(np.abs(a.psi).max())
    return cache[key]


def main() -> None:
    results = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(r["eta"], r["Ra"], r["p"], r["N"]) for r in results}
    anchors: dict[str, float] = {}

    for n in NS:
        for Ra in RAS:
            Q1 = anchor_for(Ra, n, anchors)
            for eta in ETAS:
                for p in PS:
                    if (eta, Ra, p, n) in done:
                        continue
                    t0 = time.time()
                    a = AisleCell(N=n, eta=eta, Ra=Ra, flux=p * Q1,
                                  inlet=INLET, outlet=OUTLET)
                    try:
                        its = solve_to_phi(a)
                        rec = {"eta": eta, "Ra": Ra, "p": p, "N": n, "Q1": Q1,
                               "phi": a.recirc_fraction(),
                               "phi_tracer": a.recirc_fraction_tracer(),
                               "iters": int(its), "diverged": False,
                               "seconds": round(time.time() - t0, 1)}
                    except FloatingPointError:
                        rec = {"eta": eta, "Ra": Ra, "p": p, "N": n, "Q1": Q1,
                               "phi": None, "phi_tracer": None, "diverged": True,
                               "seconds": round(time.time() - t0, 1)}
                    results.append(rec)
                    json.dump(results, open(OUT, "w"), indent=2)
                    tag = "DIVERGED" if rec["diverged"] else (
                        f"phi={rec['phi']:.4f} tracer={rec['phi_tracer']:.4f}")
                    print(f"  N={n} Ra={Ra:.0e} eta={eta:<5} p={p:<4} {tag} "
                          f"{rec['seconds']}s")
    verdict(results)


def verdict(results: list[dict]) -> None:
    def get(eta, Ra, p, n):
        for r in results:
            if (r["eta"], r["Ra"], r["p"], r["N"]) == (eta, Ra, p, n) and not r["diverged"]:
                return r
        return None

    print("\nphi from T (phi from tracer in brackets), p = 0\n")
    for n in NS:
        print(f"  N = {n}")
        head = "    eta   " + "  ".join(f"Ra={Ra:.0e}      " for Ra in RAS)
        print(head)
        for eta in ETAS:
            cells = []
            for Ra in RAS:
                r = get(eta, Ra, 0.0, n)
                cells.append("   ---        " if r is None
                             else f"{r['phi']:.4f} [{r['phi_tracer']:.4f}]")
            print(f"    {eta:<5} " + "  ".join(cells))

    lo = get(0.95, RAS[0], 0.0, 40)
    hi = get(0.95, RAS[-1], 0.0, 40)
    zero_lo = get(0.0, RAS[0], 0.0, 40)
    zero_hi = get(0.0, RAS[-1], 0.0, 40)
    if lo is None or hi is None:
        print("\n  undecided: a case in the decision pair did not converge")
        return
    drop = lo["phi"] / hi["phi"] if hi["phi"] > 0 else float("inf")
    print(f"\n  eta = 0.95, p = 0, N = 40:  phi {lo['phi']:.4f} -> {hi['phi']:.4f} "
          f"over Ra {RAS[0]:.0e} -> {RAS[-1]:.0e}, a factor of {drop:.2f}")
    print(f"    tracer over the same range:  {lo['phi_tracer']:.4f} -> {hi['phi_tracer']:.4f}")
    if zero_lo and zero_hi:
        print(f"  eta = 0    , p = 0, N = 40:  phi {zero_lo['phi']:.4f} -> "
              f"{zero_hi['phi']:.4f}, a factor of "
              f"{zero_lo['phi'] / zero_hi['phi']:.2f}  (the advective control)")

    if drop > 2.0 and hi["phi_tracer"] < 0.01:
        print("\n  VERDICT: conduction. The floor falls with Ra while the tracer stays "
              "at zero,\n           so the temperature channel at eta -> 1 is reading "
              "heat that was\n           conducted, not air that was recirculated. A "
              "hall at Ra ~ 1e9 has no\n           such floor, and the flattening of "
              "the containment ROI that rests on\n           it does not transfer.")
    elif abs(hi["phi"] - lo["phi"]) / lo["phi"] < 0.30:
        print("\n  VERDICT: advection. The floor is Ra-independent, so it is a real "
              "short circuit\n           over the partition and the containment ROI "
              "flattening stands.")
    else:
        print("\n  VERDICT: undecided against the rule fixed before the run. Report as "
              "undecided;\n           do not pick whichever branch the numbers are "
              "nearer to.")


if __name__ == "__main__":
    main()
