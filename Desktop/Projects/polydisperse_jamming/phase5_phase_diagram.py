"""phase5_phase_diagram.py — sweep the impact number J and ask what actually happens.

The Phase-4 sweep ran ONE impact velocity (u0 = 0.02) against ONE carve pressure
(p_hat = 1e-4), i.e. a single J = rho u0^2 / p_carve = 4, and reported that nothing
rebounds. That was one point in a two-dimensional space, and the aggregate it used was
only marginally jammed. A pack carved just above its own jamming point is a marginal
solid; there is no reason to expect it to store and return impact energy.

This sweeps the controlling group over four decades:

    J = p_imp / p_carve = rho u0^2 / p

    J << 1   the pack is stiff compared with the impact  -> elastic response expected
    J >> 1   the impact overwhelms the confinement       -> fluidisation, spreading

crossed with CV and with the boundary condition, at fixed friction. The output is a
phase diagram in (J, CV) rather than a single verdict.

Axes
----
carve p_hat : 1e-4, 1e-3, 1e-2, 3e-2      how far above jamming the pack is carved
u0          : 0.02, 0.08                  impact velocity (see the note on U0 below)
CV          : 0.0, 0.3, 0.5
bc          : capillary (Lambda = 50), caged
seeds       : 2

The grid originally included u0 = 0.005 and was truncated mid-run on a measured cost
over-run; 21 of those cells completed and are kept. J is covered from 0.0025 to 64,
4.4 decades, which is the full range the sweep was designed to reach.

`free` is dropped from the production grid: Phase 4 already established at 13-45 sigma
that a free aggregate disintegrates regardless of anything else, so spending a third of
the cells re-establishing that buys nothing.
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

CARVE = (1e-4, 1e-3, 1e-2, 3e-2)
# u0 = 0.005 was TRUNCATED mid-sweep on a measured cost over-run. Flight time
# scales as 1/u0, so those cells cost ~16x a u0 = 0.08 cell and the arm alone
# was ~3.2 h of the remaining budget. It had already delivered 21 cells and
# the low-J corner: J was covered 0.0025 to 64 (4.4 decades) at truncation,
# which is the full designed range. Cells completed before truncation are
# kept and analysed; the grid is simply no longer generating more of them.
U0 = (0.02, 0.08)
CVS = (0.0, 0.3, 0.5)
BCS = (("capillary", 50.0), ("caged", 0.0))
SEEDS = (1, 2)
MU = 0.3


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--n", type=int, default=1500)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()

    from phase4_impact import run_impact

    os.makedirs(OUT_DIR, exist_ok=True)
    carve = (1e-4, 1e-2) if args.quick else CARVE
    u0s = (0.005, 0.08) if args.quick else U0
    cvs = (0.0,) if args.quick else CVS
    seeds = (1,) if args.quick else SEEDS

    suffix = "_quick" if args.quick else ""
    if args.nshards > 1:
        suffix += f"_shard{args.shard:02d}of{args.nshards:02d}"
    path = os.path.join(OUT_DIR, f"phase5{suffix}.jsonl")

    done = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["key"])
                except Exception:
                    continue

    combos = list(itertools.product(carve, u0s, cvs, BCS, seeds))
    print(f"phase5: {len(combos)} cells, N_parent={args.n}, shard {args.shard}/"
          f"{args.nshards} -> {path}", flush=True)

    t0 = time.time()
    for idx, (cph, u0, cv, (bc, lam), seed) in enumerate(combos):
        if idx % args.nshards != args.shard:
            continue
        key = (f"cv{cv:.3f}_mu{MU:.3f}_s{seed}_{bc}_lam{lam:.4g}"
               f"_cp{cph:.3e}_u{u0:.4g}_N{args.n}")
        if key in done:
            continue
        try:
            st = run_impact(cv, MU, seed, args.n, bc, lam, u0, carve_p_hat=cph)
            rec = json.loads(st.to_json())
        except Exception as exc:      # noqa: BLE001
            print(f"CELL FAILED {key}: {exc!r}", flush=True)
            rec = {"failed": True, "error": repr(exc), "cv": cv, "mu_p": MU,
                   "seed": seed, "bc": bc, "lam": lam, "carve_p_hat": cph,
                   "u0": u0, "n_parent": args.n}
        rec["key"] = key
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
        print(f"  [{idx+1}/{len(combos)}] elapsed {(time.time()-t0)/60:.1f} min",
              flush=True)
    print(f"phase5 done in {(time.time()-t0)/60:.1f} min", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
