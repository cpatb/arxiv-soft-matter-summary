#!/bin/bash
# probe_seed3_crash.sh -- WHY does lognormal (alpha 2, CV 0.5, seed 3) abort in the collider?
#
#   wsl bash /home/cpb/pj/scripts/probe_seed3_crash.sh
#
# The cell dies with a C++ assertion in InsertionSortCollider::spatialOverlapPeri
# (`maxima - minima < .99*dim`), which is an abort, not a Python exception: the runner
# cannot catch it, no `failed` record is written, and the cell simply never lands. It did
# so at all three rungs in two separate launches (2026-09-16), so it is the INITIAL
# configuration. This probe builds that configuration through the shipped `insert_packing`
# and reports, before any step, the deepest member-member overlap in the pack and the
# smallest body involved -- then steps one iteration at a time printing max |v| until the
# abort, so the last line printed is the state that killed it. The alpha = 1 and seed 1-2
# controls run the same lines.
source "$HOME/pj/scripts/env.sh"
set -u
cd ~/pj || exit 1
export OMP_NUM_THREADS=1

python3 - "$@" <<'PYEOF'
import math, sys
import numpy as np
from yade import O, utils
import packing
from packing import PackingState, insert_packing, interaction_loop, DT_SAFETY, PHI_INSERT

def probe(family, cv, alpha, seed, n=5400, steps=60):
    O.reset()
    st = PackingState(cv=cv, mu_p=0.3, p_hat_target=1e-3, n_requested=n, seed=seed,
                      family=family, eta_roll=None, alpha=float(alpha))
    L, mat = insert_packing(st, truncate=(0.2, 5.0), young=1.0)
    bodies = [b for b in O.bodies if b is not None and not b.isClump]
    pos = np.array([b.state.pos for b in bodies])
    rad = np.array([b.shape.radius for b in bodies])
    ids = np.array([b.id for b in bodies])
    grain = np.array([b.clumpId if b.isClumpMember else b.id for b in bodies])
    # deepest overlap between members of DIFFERENT grains, periodic minimum image
    from scipy.spatial import cKDTree
    tree = cKDTree(np.mod(pos, L), boxsize=L)
    pairs = tree.query_pairs(2.0 * rad.max(), output_type="ndarray")
    d = pos[pairs[:, 0]] - pos[pairs[:, 1]]
    d -= L * np.round(d / L)
    dist = np.linalg.norm(d, axis=1)
    pen = rad[pairs[:, 0]] + rad[pairs[:, 1]] - dist
    other = grain[pairs[:, 0]] != grain[pairs[:, 1]]
    pen = np.where(other, pen, -np.inf)
    k = int(np.argmax(pen))
    i, j = pairs[k]
    rel = pen[k] / min(rad[i], rad[j])
    print(f"--- {family} alpha={alpha} CV={cv} seed={seed}: N_bodies={len(bodies)} L={L:.3f} "
          f"r_min={rad.min():.4f} r_max={rad.max():.4f}")
    print(f"    deepest cross-grain overlap: pen={pen[k]:.4f} between r={rad[i]:.4f} and "
          f"r={rad[j]:.4f} (pen / r_small = {rel:.3f}; > 2 means one body's centre lies "
          f"INSIDE the other)")
    n_inside = int(np.sum(pen > 2.0 * np.minimum(rad[pairs[:, 0]], rad[pairs[:, 1]])))
    print(f"    N: {n_inside} cross-grain pairs with a centre inside the other body; "
          f"{int(np.sum(pen > 0))} overlapping pairs of {len(pairs)} candidate pairs")

    O.engines = [
        utils.ForceResetter(),
        utils.InsertionSortCollider([utils.Bo1_Sphere_Aabb()], allowBiggerThanPeriod=False),
        interaction_loop(utils, None),
        utils.NewtonIntegrator(damping=0.4),
    ]
    O.dt = DT_SAFETY * utils.PWaveTimeStep()
    print(f"    dt={O.dt:.4e}", flush=True)
    for it in range(steps):
        O.run(1, True)
        v = np.array([np.linalg.norm(b.state.vel) for b in O.bodies if b is not None and not b.isClump])
        bmax = int(np.argmax(v))
        print(f"    iter {O.iter:3d}: max|v|={v.max():.3e} (body r={rad[bmax]:.4f}) "
              f"max|v|*dt/L={v.max()*O.dt/L:.3e}", flush=True)
        if v.max() * O.dt > 0.3 * L:
            print("    >>> a body would cross a third of the cell in one step; this is "
                  "the collider abort in the making", flush=True)
            break

args = sys.argv[1:]
if args:
    fam, cv, a, s = args[0], float(args[1]), float(args[2]), int(args[3])
    probe(fam, cv, a, s)
else:
    probe("lognormal", 0.5, 2.0, 3)
PYEOF
