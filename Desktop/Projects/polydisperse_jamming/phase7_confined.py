"""phase7_confined.py — the impact experiment, with the aggregate actually jammed.

WHAT WENT WRONG IN PHASES 4 AND 5
---------------------------------
`diag_carve_release.py` shows the aggregates were losing 50-100% of their coordination
during the 300-step settle, BEFORE being given any velocity. Worst case: capillary,
CV = 0.5, p_carve = 1e-3, Z falling 2.382 -> 0.100, a 95.8% loss with no impact
involved.

The cause is a pressure mismatch that was built into the design. The pack is carved out
of a periodic cell equilibrated at `p_carve`, but the confinement holding it afterwards
is `p_cap = p_imp / Lambda`. At the production settings that is 4e-4/50 = 8e-6 against a
carve pressure of 1e-3: the confinement is ~100x too weak, so the aggregate expands
until its contacts open.

So Phases 4 and 5 delivered their impacts to a loose cloud held roughly in place by
surface tension, not to a jammed solid. "Nothing rebounds" is true of what was simulated
and says little about a jammed drop, because nothing was jammed at the moment of impact.
It also explains why Phase 6 worked: bonding the contacts at step 1 freezes the network
before the settle can dismantle it.

THE FIX
-------
Hold the aggregate at a confining pressure equal to the pressure it was carved at:

    p_cap = p_carve      so that      Lambda = p_imp / p_cap = J

which collapses the two control parameters onto one. Sweeping J then genuinely means
"impact stress against confining stress" with the pack in mechanical equilibrium at
launch rather than expanding out of it.

THE GATE
--------
Every cell reports `z_retention = Z_launch / Z_carved`, and a cell is only admitted to
the analysis if it is above `Z_RETENTION_MIN`. That is the check Phases 4 and 5 lacked:
it makes "the aggregate was still jammed when it hit" a measured precondition rather
than an assumption. Cells that fail it are written out with the flag set, not silently
dropped.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import sys
import time

import numpy as np

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

POISSON = 0.30
DT_SAFETY = 0.4
DAMPING = 0.2
SHELL_FRACTION = 0.85
SETTLE_STEPS = 300
N_PROBES = 100

# A cell is "jammed at launch" if it kept this fraction of its as-carved coordination.
Z_RETENTION_MIN = 0.80

CARVE = (1e-4, 1e-3, 1e-2)
U0 = (0.02, 0.08)
CVS = (0.0, 0.25, 0.5)
SEEDS = (1, 2)
MU = 0.3


def run_confined(cv, mu_p, seed, n_parent, carve_p_hat, u0, young=1.0):
    """Carve at p_carve, confine at p_carve, then impact. Reports the jamming gate."""
    from yade import O, utils

    from phase4_impact import (_com_vel, _grain_graph, _ke, _lateral, carve_aggregate)

    t0 = time.time()
    pos, rad, phi_parent, mom = carve_aggregate(cv, mu_p, n_parent, seed, young,
                                                carve_p_hat)
    if pos.shape[0] < 50:
        return {"failed": True, "error": f"carve gave {pos.shape[0]} grains"}

    O.reset()
    O.periodic = False
    mat = O.materials.append(
        utils.FrictMat(young=young, poisson=POISSON, density=1.0,
                       frictionAngle=math.atan(mu_p), label="grain"))
    ids = [O.bodies.append(utils.sphere(tuple(p), float(r), material=mat))
           for p, r in zip(pos, rad)]
    r_agg = float(np.linalg.norm(pos, axis=1).max())
    O.bodies.append(utils.wall(-r_agg - 1.0, axis=2, sense=1))

    # THE FIX: confine at the pressure the pack was carved at.
    p_cap = carve_p_hat
    p_imp = u0**2 / young
    state = {"p": p_cap, "ids": ids}

    def confine():
        pc = state["p"]
        if pc <= 0:
            return
        i = state["ids"]
        P = np.array([[O.bodies[b].state.pos[j] for j in range(3)] for b in i])
        rel = P - P.mean(axis=0)
        rr = np.linalg.norm(rel, axis=1)
        rmax = rr.max()
        shell = rr > SHELL_FRACTION * rmax
        if not shell.any():
            return
        f = pc * 4.0 * math.pi * rmax**2 / shell.sum()
        for k, b in enumerate(i):
            if shell[k] and rr[k] > 0:
                O.forces.addF(b, tuple(-rel[k] / rr[k] * f))

    O.engines = [
        utils.ForceResetter(),
        utils.InsertionSortCollider([utils.Bo1_Sphere_Aabb(), utils.Bo1_Wall_Aabb()]),
        utils.InteractionLoop(
            [utils.Ig2_Sphere_Sphere_ScGeom(), utils.Ig2_Wall_Sphere_ScGeom()],
            [utils.Ip2_FrictMat_FrictMat_FrictPhys()],
            [utils.Law2_ScGeom_FrictPhys_CundallStrack()],
        ),
        utils.PyRunner(command="confine()", iterPeriod=1, label="conf"),
        utils.NewtonIntegrator(damping=DAMPING),
    ]
    O.dt = DT_SAFETY * utils.PWaveTimeStep()
    import __main__
    __main__.confine = confine

    O.run(1, wait=True)
    z_carved = 2.0 * (sum(_grain_graph(O, ids)[0].values()) // 2) / len(ids)
    O.run(SETTLE_STEPS, wait=True)
    z_launch = 2.0 * (sum(_grain_graph(O, ids)[0].values()) // 2) / len(ids)
    retention = z_launch / z_carved if z_carved > 0 else float("nan")

    for b in ids:
        O.bodies[b].state.vel = (0.0, 0.0, -u0)
    d0 = _lateral(O, ids)
    ke0 = _ke(O, ids)

    t_flight = (r_agg + 1.0) / u0
    n_steps = max(int(4.0 * t_flight / O.dt), 4000)
    every = max(1, n_steps // N_PROBES)
    z_peak, cpk = float("nan"), -1
    for _ in range(N_PROBES):
        O.run(every, wait=True)
        c, _n, _f = _grain_graph(O, ids)
        nreal = sum(c.values()) // 2
        if nreal > cpk:
            cpk, z_peak = nreal, 2.0 * nreal / len(ids)

    v = _com_vel(O, ids)
    out = {
        "cv": cv, "mu_p": mu_p, "seed": seed, "n_parent": n_parent,
        "n_grains": len(ids), "carve_p_hat": carve_p_hat, "p_cap": p_cap,
        "u0": u0, "p_imp": p_imp, "j_impact": p_imp / carve_p_hat,
        "phi_parent": phi_parent, "cv_realised": float(mom.get("cv", float("nan"))),
        "r_agg": r_agg, "z_carved": z_carved, "z_launch": z_launch,
        "z_retention": retention, "jammed_at_launch": bool(retention >= Z_RETENTION_MIN),
        "z_peak": z_peak, "contacts_peak": cpk,
        "cor": float(v[2] / u0), "cor_convention": "positive_is_rebound",
        "ke_ratio": float(_ke(O, ids) / ke0) if ke0 > 0 else float("nan"),
        "spread": float(_lateral(O, ids) / d0) if d0 > 0 else float("nan"),
        "wall_s": time.time() - t0,
    }
    print(f"  cv={cv:.2f} pc={carve_p_hat:.0e} u0={u0:.3f} J={out['j_impact']:.3g} -> "
          f"Zc={z_carved:.2f} Zl={z_launch:.2f} ret={retention:.2f} "
          f"{'JAMMED' if out['jammed_at_launch'] else 'NOT-JAMMED'} "
          f"Zpk={z_peak:.2f} COR={out['cor']:+.4f} spread={out['spread']:.3f} "
          f"{out['wall_s']:.0f}s", flush=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--n", type=int, default=1500)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    carve = (1e-3,) if args.quick else CARVE
    u0s = (0.02,) if args.quick else U0
    cvs = (0.0,) if args.quick else CVS
    seeds = (1,) if args.quick else SEEDS

    suffix = "_quick" if args.quick else ""
    if args.nshards > 1:
        suffix += f"_shard{args.shard:02d}of{args.nshards:02d}"
    path = os.path.join(OUT_DIR, f"phase7{suffix}.jsonl")

    done = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["key"])
                except Exception:
                    continue

    combos = list(itertools.product(carve, u0s, cvs, seeds))
    print(f"phase7: {len(combos)} cells, N_parent={args.n}, shard {args.shard}/"
          f"{args.nshards} -> {path}", flush=True)

    t0 = time.time()
    for idx, (cph, u0, cv, seed) in enumerate(combos):
        if idx % args.nshards != args.shard:
            continue
        key = f"cv{cv:.3f}_s{seed}_cp{cph:.3e}_u{u0:.4g}_N{args.n}"
        if key in done:
            continue
        try:
            rec = run_confined(cv, MU, seed, args.n, cph, u0)
        except Exception as exc:      # noqa: BLE001
            print(f"CELL FAILED {key}: {exc!r}", flush=True)
            rec = {"failed": True, "error": repr(exc), "cv": cv, "seed": seed,
                   "carve_p_hat": cph, "u0": u0, "n_parent": args.n}
        rec["key"] = key
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
    print(f"phase7 done in {(time.time()-t0)/60:.1f} min", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
