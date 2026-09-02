"""phase6_cohesion.py — how much tensile strength does a rebound actually need?

WHY THIS EXISTS
---------------
Phases 4 and 5 report that nothing rebounds: |COR| <= 0.035 everywhere. That result is
now anchored -- `test_cor_anchor.py` shows the same measurement path returns COR = 0.819
for a perfectly elastic sphere at the production damping, so a rebound was detectable
23x over. The negative is real.

But "a cohesionless frictional aggregate does not rebound" is only half an answer. The
useful question is the quantitative one: **how much tensile strength would it need?**
A dense suspension held together by a jammed contact network has no cohesion in this
model, yet a real shear-jammed drop is reported to bounce. Something supplies tension.

So this sweeps the cohesion, expressed as the ratio of the bond strength to the impact
stress,

    C = sigma_coh / p_imp = sigma_coh / (rho u0^2)

from 0 (which must reproduce Phase 4/5) upward, and reads off where COR lifts off the
floor. C = 0 is a DISCRIMINATING NEGATIVE: if the cohesive machinery reports a rebound
at zero cohesion, it is broken and nothing else here can be believed.

The output is a number the bench can use: the tensile strength, relative to the impact
stress, at which an aggregate of this kind starts to come back.

CONTACT MODEL
-------------
`CohFrictMat` with `Ip2_CohFrictMat_CohFrictMat_CohFrictPhys`, bonded by
`setCohesionNow` after carving. Every contact must produce `ScGeom6D`, because the cohesive law reads 6 degrees of
freedom. The first version used a `wall`, whose only functor emits a plain `ScGeom`, and
paired it with `Law2_ScGeom_CohFrictPhys_CohesionMoment` -- which turns out to be a
DEPRECATED ALIAS for the 6D law. Yade then read an `ScGeom` as an `ScGeom6D` and all six
workers died with `malloc(): invalid size (unsorted)`, i.e. heap corruption rather than a
clean error.

The floor is therefore a `box`, which has `Ig2_Box_Sphere_ScGeom6D`, so the dispatcher is
6D throughout and only one law is needed.
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
DAMPING = 0.2          # matches phase4/5, and test_cor_anchor measured its COR ceiling
SETTLE_STEPS = 300
N_PROBES = 100
CARVE_P_HAT = 1e-3     # a moderately jammed parent, mid-range of the Phase-5 ladder

# sigma_coh / p_imp. 0 is the discriminating negative; it must reproduce Phase 4/5.
#
# The ladder runs to 300 because a first probe at C = 10 still shattered the aggregate
# (spread 3.2, COR +0.016). A positive control has to BRACKET the transition: if the
# top of the ladder still does not rebound, the control has failed and says so, rather
# than leaving 'no rebound' resting on a range that was never wide enough.
COH_RATIOS = (0.0, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0)
# CV trimmed to the two ends: this sweep is about the cohesion axis, and 7 x 3 x 2 was
# more cells than the axis needs.
CVS = (0.0, 0.5)
SEEDS = (1, 2)
MU = 0.3
U0 = 0.02


def run_cohesive(cv, mu_p, seed, n_parent, coh_ratio, u0=U0, young=1.0):
    """Carve an aggregate, bond it at the given cohesion, drive it into a wall."""
    from yade import O, utils

    from phase4_impact import (_com_vel, _grain_graph, _ke, _lateral, _positions)

    t0 = time.time()
    # --- carve, using the SAME parent-packing machinery as Phase 4/5 --------
    from phase4_impact import carve_aggregate
    pos, rad, phi_parent, mom = carve_aggregate(cv, mu_p, n_parent, seed, young,
                                                CARVE_P_HAT)
    if pos.shape[0] < 50:
        return {"failed": True, "error": f"carve gave {pos.shape[0]} grains"}

    p_imp = u0**2 / young
    sigma_coh = coh_ratio * p_imp

    O.reset()
    O.periodic = False
    mat = O.materials.append(utils.CohFrictMat(
        young=young, poisson=POISSON, density=1.0,
        frictionAngle=math.atan(mu_p),
        normalCohesion=sigma_coh, shearCohesion=sigma_coh,
        momentRotationLaw=True, label="grain"))
    ids = [O.bodies.append(utils.sphere(tuple(p), float(r), material=mat))
           for p, r in zip(pos, rad)]

    r_agg = float(np.linalg.norm(pos, axis=1).max())
    wall_z = -r_agg - 1.0
    # A BOX, not a wall: only Box (and Facet) have a 6D geometry functor, and the
    # cohesive law needs one. Sized 20x the aggregate so it is effectively a half-space.
    half = 10.0 * r_agg
    floor = O.bodies.append(utils.box(center=(0, 0, wall_z - half),
                                      extents=(half, half, half),
                                      fixed=True, material=mat))

    O.engines = [
        utils.ForceResetter(),
        utils.InsertionSortCollider([utils.Bo1_Sphere_Aabb(), utils.Bo1_Box_Aabb()]),
        utils.InteractionLoop(
            [utils.Ig2_Sphere_Sphere_ScGeom6D(), utils.Ig2_Box_Sphere_ScGeom6D()],
            [utils.Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(
                setCohesionNow=True, setCohesionOnNewContacts=False)],
            [utils.Law2_ScGeom6D_CohFrictPhys_CohesionMoment()],
        ),
        utils.NewtonIntegrator(damping=DAMPING),
    ]
    O.dt = DT_SAFETY * utils.PWaveTimeStep()

    # one step to build contacts, which is when setCohesionNow bonds them
    O.run(1, wait=True)
    counts, neigh, _f = _grain_graph(O, ids)
    z_initial = 2.0 * (sum(counts.values()) // 2) / len(counts) if counts else float("nan")
    # bonding is a one-shot event; switch it off so new contacts during impact are NOT
    # cohesive. Otherwise the aggregate would re-bond to the wall and to itself on
    # every fresh contact, which is not a material property.
    O.engines[2].physDispatcher.functors[0].setCohesionNow = False

    O.run(SETTLE_STEPS, wait=True)
    z_launch = 2.0 * (sum(_grain_graph(O, ids)[0].values()) // 2) / len(ids)

    for i in ids:
        O.bodies[i].state.vel = (0.0, 0.0, -u0)
    d0_lat = _lateral(O, ids)
    ke0 = _ke(O, ids)

    t_flight = (r_agg + 1.0) / u0
    n_steps = max(int(4.0 * t_flight / O.dt), 4000)
    every = max(1, n_steps // N_PROBES)
    z_peak, contacts_peak = float("nan"), -1
    for _ in range(N_PROBES):
        O.run(every, wait=True)
        c, _n, _ff = _grain_graph(O, ids)
        nreal = sum(c.values()) // 2
        if nreal > contacts_peak:
            contacts_peak = nreal
            z_peak = 2.0 * nreal / len(ids)

    v = _com_vel(O, ids)
    out = {
        "cv": cv, "mu_p": mu_p, "seed": seed, "n_parent": n_parent,
        "n_grains": len(ids), "coh_ratio": coh_ratio, "sigma_coh": sigma_coh,
        "u0": u0, "p_imp": p_imp, "carve_p_hat": CARVE_P_HAT,
        "phi_parent": phi_parent, "cv_realised": float(mom.get("cv", float("nan"))),
        "r_agg": r_agg, "z_initial": z_initial, "z_launch": z_launch,
        "z_peak": z_peak, "contacts_peak": contacts_peak,
        "cor": float(v[2] / u0), "cor_convention": "positive_is_rebound",
        "ke_ratio": float(_ke(O, ids) / ke0) if ke0 > 0 else float("nan"),
        "spread": float(_lateral(O, ids) / d0_lat) if d0_lat > 0 else float("nan"),
        "wall_s": time.time() - t0,
    }
    _ = (_positions, floor)
    print(f"  cv={cv:.2f} C={coh_ratio:6.1f} N={len(ids)} -> Zi={z_initial:.2f} "
          f"Zpk={z_peak:.2f} COR={out['cor']:+.4f} spread={out['spread']:.3f} "
          f"KE={out['ke_ratio']:.4f} {out['wall_s']:.0f}s", flush=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--n", type=int, default=1500)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--ratios", default="",
                    help="comma-separated C values, overriding COH_RATIOS")
    ap.add_argument("--cvs", default="", help="comma-separated CV values")
    ap.add_argument("--seeds", default="", help="comma-separated seeds")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    ratios = (0.0, 10.0) if args.quick else COH_RATIOS
    cvs = (0.0,) if args.quick else CVS
    seeds = (1,) if args.quick else SEEDS
    if args.ratios:
        ratios = tuple(float(x) for x in args.ratios.split(","))
    if args.cvs:
        cvs = tuple(float(x) for x in args.cvs.split(","))
    if args.seeds:
        seeds = tuple(int(x) for x in args.seeds.split(","))

    suffix = "_quick" if args.quick else ""
    if args.tag:
        suffix += f"_{args.tag}"
    if args.nshards > 1:
        suffix += f"_shard{args.shard:02d}of{args.nshards:02d}"
    path = os.path.join(OUT_DIR, f"phase6{suffix}.jsonl")

    done = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["key"])
                except Exception:
                    continue

    combos = list(itertools.product(ratios, cvs, seeds))
    print(f"phase6: {len(combos)} cells, N_parent={args.n}, shard {args.shard}/"
          f"{args.nshards} -> {path}", flush=True)

    t0 = time.time()
    for idx, (ratio, cv, seed) in enumerate(combos):
        if idx % args.nshards != args.shard:
            continue
        key = f"cv{cv:.3f}_C{ratio:.4g}_s{seed}_N{args.n}"
        if key in done:
            continue
        try:
            rec = run_cohesive(cv, MU, seed, args.n, ratio)
        except Exception as exc:      # noqa: BLE001
            print(f"CELL FAILED {key}: {exc!r}", flush=True)
            rec = {"failed": True, "error": repr(exc), "cv": cv,
                   "coh_ratio": ratio, "seed": seed, "n_parent": args.n}
        rec["key"] = key
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
        print(f"  [{idx+1}/{len(combos)}] elapsed {(time.time()-t0)/60:.1f} min",
              flush=True)
    print(f"phase6 done in {(time.time()-t0)/60:.1f} min", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
