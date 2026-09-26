"""ls_pack.py -- isotropic compression of CONVEX, FACETED-TENDING grains (level-set superellipsoids).

Why: the clump ladder (PREREG_ROLLCHI arm A) cannot say what a CONVEX faceted grain does, because a two-sphere
clump has a concave neck that neighbours nestle into. Cornstarch granules are convex and polyhedral-ish. A
superellipsoid with equal half-extents and exponents eps_e = eps_n = eps is a sphere at eps = 1 and a
progressively sharper rounded cube as eps -> 0: convex at every eps, no neck.

Protocol mirrors packing.build_and_compress: periodic cell, random loose insertion, PeriTriaxController to a
reduced pressure p_hat = p/E with the SAME tolerances (UNBALANCED_TOL, REL_STRESS_TOL), damping 0.4,
dt = DT_SAFETY * r_eq * sqrt(rho/E) (PWaveTimeStep's formula with the volume-equivalent radius, since
PWaveTimeStep only reads spheres). Contacts: Ig2_LevelSet_LevelSet_ScGeom (single contact point) with the
project's linear FrictMat / Cundall-Strack law. Monodisperse (CV 0) only.

Observables: phi = N V_grain / V_cell with V_grain from the level set's own volume; z_all and z_nr (iterative
k = 2 rattler removal, the project's default z_nr). Appends one JSON line per cell to out/lsdem/ls_<tag>.jsonl.

    python3 ls_pack.py --eps 1.0 --mu 0.3 --phat 1e-3 --n 200 --seed 1 --tag gate
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time

import numpy as np

import packing

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out", "lsdem")
R_EQ = 0.5          # volume-equivalent radius (d_eq = 1, as in packing.py)
RHO = 1.0
YOUNG = 1.0


def superellipsoid_volume(r: float, eps: float) -> float:
    """Volume of |x/r|^(2/eps) + |y/r|^(2/eps) + |z/r|^(2/eps) = 1 (eps_e = eps_n = eps):
    V = 2 r^3 eps^2 B(eps/2, eps/2) B(eps, eps/2 + 1)  (standard superquadric volume, a1 a2 a3 = r^3; source
    not verified here -- the ANCHORS are what certify it). Checked 2026-09-25: eps = 1 -> 4.18879 = 4/3 pi,
    eps = 1e-3 -> 8.000, and a 2e6-point Monte Carlo at eps = 0.5 -> 6.485 vs 6.482. A first version with a
    2/3 prefactor was low by exactly 3 at every eps and was caught by those anchors."""
    from scipy.special import beta
    return 2.0 * r ** 3 * eps * eps * beta(eps / 2, eps / 2) * beta(eps, eps / 2 + 1)


def z_nonrattler(neigh: dict[int, set[int]], k: int = 2) -> float:
    alive = {i for i, s in neigh.items() if len(s) >= k}
    changed = True
    while changed:
        changed = False
        for i in list(alive):
            if len(neigh[i] & alive) < k:
                alive.discard(i)
                changed = True
    if not alive:
        return float("nan")
    return sum(len(neigh[i] & alive) for i in alive) / len(alive)


def run(eps: float, mu: float, phat: float, n: int, seed: int, spacing_frac: float, nodes: int,
        max_iter: int) -> dict:
    from yade import O, Quaternion, Vector3, utils
    levelSetBody = utils.levelSetBody
    t0 = time.time()
    O.reset()
    # extents chosen so the grain's VOLUME equals the d_eq = 1 sphere at every eps
    r = R_EQ * ((4.0 / 3.0) * math.pi / superellipsoid_volume(1.0, eps)) ** (1.0 / 3.0)
    v_grain = superellipsoid_volume(r, eps)
    mat = O.materials.append(utils.FrictMat(young=YOUNG, poisson=packing.POISSON, density=RHO,
                                            frictionAngle=math.atan(mu)))
    phi0 = 0.10
    L = (n * v_grain / phi0) ** (1.0 / 3.0)
    O.periodic = True
    O.cell.hSize = np.diag([L, L, L])
    rng = np.random.default_rng(seed)
    q = rng.normal(size=(n, 4))
    q /= np.linalg.norm(q, axis=1)[:, None]
    for j in range(n):
        pos = rng.uniform(0, L, 3)
        ori = Quaternion(float(q[j, 0]), float(q[j, 1]), float(q[j, 2]), float(q[j, 3]))
        b = levelSetBody("superellipsoid", center=Vector3(*pos), extents=Vector3(r, r, r),
                         epsilons=(eps, eps), spacing=spacing_frac * 2 * r, nSurfNodes=nodes,
                         orientation=ori, material=mat)
        O.bodies.append(b)
    goal = -abs(phat) * YOUNG
    O.engines = [
        utils.ForceResetter(),
        utils.InsertionSortCollider([utils.Bo1_LevelSet_Aabb()], verletDist=0.0),
        utils.InteractionLoop([utils.Ig2_LevelSet_LevelSet_ScGeom()],
                              [utils.Ip2_FrictMat_FrictMat_FrictPhys()],
                              [utils.Law2_ScGeom_FrictPhys_CundallStrack()]),
        utils.PeriTriaxController(label="triax", goal=(goal, goal, goal), stressMask=7, dynCell=True,
                                  maxUnbalanced=packing.UNBALANCED_TOL, relStressTol=packing.REL_STRESS_TOL,
                                  growDamping=0.4, doneHook="O.pause()"),
        utils.NewtonIntegrator(damping=0.4),
    ]
    O.dt = packing.DT_SAFETY * R_EQ * math.sqrt(RHO / YOUNG)
    O.run(max_iter, wait=True)
    vol = O.cell.volume
    neigh: dict[int, set[int]] = {b.id: set() for b in O.bodies}
    nc = 0
    for i in O.interactions:
        if i.isReal:
            neigh[i.id1].add(i.id2)
            neigh[i.id2].add(i.id1)
            nc += 1
    unb = utils.unbalancedForce()
    s = utils.getStress()
    p = -float(s[0, 0] + s[1, 1] + s[2, 2]) / 3.0
    v_ls = float(np.mean([b.shape.volume() for b in O.bodies]))
    return dict(eps=eps, mu_p=mu, p_hat_target=phat, p_hat=p / YOUNG, n=n, seed=seed, r_extent=r,
                v_grain_analytic=v_grain, v_grain_levelset=v_ls, phi=n * v_grain / vol,
                phi_ls_volume=n * v_ls / vol, z_all=2.0 * nc / n, z_nr=z_nonrattler(neigh, 2),
                rattler_frac=1.0 - sum(1 for i in neigh if len(neigh[i]) >= 2) / n,
                unbalanced=float(unb), iters=int(O.iter), converged=bool(unb < packing.UNBALANCED_TOL and O.iter < max_iter),
                dt=float(O.dt), spacing_frac=spacing_frac, nodes=nodes, wall_s=time.time() - t0)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--eps", type=float, required=True)
    ap.add_argument("--mu", type=float, required=True)
    ap.add_argument("--phat", type=float, required=True)
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--spacing", type=float, default=0.05, help="grid spacing / extent diameter")
    ap.add_argument("--nodes", type=int, default=102)
    ap.add_argument("--max-iter", type=int, default=400000)
    ap.add_argument("--tag", default="probe")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    rec = run(a.eps, a.mu, a.phat, a.n, a.seed, a.spacing, a.nodes, a.max_iter)
    rec["tag"] = a.tag
    with open(os.path.join(OUT, f"ls_{a.tag}.jsonl"), "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec, sort_keys=True) + "\n")
    print(json.dumps(rec), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
