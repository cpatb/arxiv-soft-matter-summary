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
# Yade's Law2 default sphericalBodies=True computes torques from RADII, so a NORMAL force never torques a grain:
# right for spheres, wrong for any other shape (Yade's own examples/levelSet/levelSetBody.py: 'sphericalBodies=False
# is required in the general case'). Every level-set record before 2026-09-27 ran with the default, so frictionless
# grains could NOT ROTATE (scripts/ls_facepair.py: a 5-deg tilt stays 5.00 deg with True, relaxes to 0.02 with False).
SPHERICAL_BODIES = False
# ...except for the SPHERE itself (eps 1), where True is the exact physics (a normal force through the centre has no
# torque) and False is wrong: the contact point is the deepest surface NODE, laterally off the line of centres by up
# to ~half a node spacing, so False injects torque noise that spins frictionless and frictional spheres alike and the
# unbalanced force plateaus at 0.12-0.13 instead of reaching 0.02 (probe 2026-09-27, eps 1 mu 0.3 N 100 seed 1).
# Pass --spherical-bodies 1 for eps 1 (PREREG_FACET Am. 3).
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


def contact_census(O) -> dict:
    """Per-contact geometry, read AFTER the run (no effect on the trajectory).

    For each real contact: the contact normal n expressed in each grain's BODY frame, whose largest |component|
    (m1, m2) is 1 for a contact dead-centre on a cube face and ~0.58 on a corner diagonal; for a sphere the
    orientation is arbitrary, so the eps = 1 cells give the random-orientation NULL of the same statistic. Also
    the torque lever |(c - x_k) x n| of the normal force about each centre (0 for a sphere), with c the contact
    point wrapped to the minimum image of the (diagonal) periodic cell, and |F_n|."""
    L = np.array([O.cell.hSize[k, k] for k in range(3)])
    rows = []
    for i in O.interactions:
        if not i.isReal:
            continue
        n = i.geom.normal
        c = np.array(i.geom.contactPoint)
        rec = [i.id1, i.id2]
        for bid, sgn in ((i.id1, 1.0), (i.id2, -1.0)):
            st = O.bodies[bid].state
            nb = st.ori.conjugate() * (n * sgn)
            rec.append(max(abs(nb[0]), abs(nb[1]), abs(nb[2])))
            d = c - np.array(st.pos)
            d -= L * np.round(d / L)
            rec.append(float(np.linalg.norm(np.cross(d, np.array(n)))))
        rec.append(float(i.phys.normalForce.norm()))
        rows.append(rec)
    a = np.array(rows, dtype=float).reshape(-1, 7)
    return dict(id1=a[:, 0].astype(int), id2=a[:, 1].astype(int), m1=a[:, 2], arm1=a[:, 3], m2=a[:, 4],
                arm2=a[:, 5], fn=a[:, 6], cell=L)


def _state(O, n: int, v_grain: float, utils) -> dict:
    """phi, pressure, force imbalance and TORQUE imbalance: mean |torque| per body over (mean |F_n| * R_EQ)."""
    fns = [i.phys.normalForce.norm() for i in O.interactions if i.isReal]
    tq = [O.forces.t(b.id).norm() for b in O.bodies]
    s = utils.getStress()
    return dict(iter=int(O.iter), phi=n * v_grain / O.cell.volume, p_hat=-float(s[0, 0] + s[1, 1] + s[2, 2]) / 3.0 / YOUNG,
                unbalanced=float(utils.unbalancedForce()),
                torque_unb=float(np.mean(tq) / (np.mean(fns) * R_EQ)) if fns else float("nan"),
                z_all=2.0 * len(fns) / n)


def run(eps: float, mu: float, phat: float, n: int, seed: int, spacing_frac: float, nodes: int,
        max_iter: int, dump_path: str | None = None, post_iters: int = 0, sph: bool = SPHERICAL_BODIES) -> dict:
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
                              [utils.Law2_ScGeom_FrictPhys_CundallStrack(sphericalBodies=sph)]),
        utils.PeriTriaxController(label="triax", goal=(goal, goal, goal), stressMask=7, dynCell=True,
                                  maxUnbalanced=packing.UNBALANCED_TOL, relStressTol=packing.REL_STRESS_TOL,
                                  growDamping=0.4, doneHook="O.pause()"),
        utils.NewtonIntegrator(damping=0.4),
    ]
    if os.environ.get("LS_PROGRESS"):  # diagnostic only: prints state, does not touch the dynamics
        O.engines = O.engines + [utils.PyRunner(iterPeriod=int(os.environ["LS_PROGRESS"]), command=
            "from yade import O, utils; print('PROGRESS', O.iter, round(O.cell.volume,4), round(utils.unbalancedForce(),4), "
            "-sum(utils.getStress()[k,k] for k in range(3))/3, flush=True)")]

    O.dt = packing.DT_SAFETY * R_EQ * math.sqrt(RHO / YOUNG)
    O.run(max_iter, wait=True)
    post = []
    if post_iters > 0:
        # The stop rule above balances FORCES only; a non-sphere can stop with unrelaxed TORQUES. Keep the same
        # stress control running (doneHook disabled) and record phi, force and torque imbalance per chunk.
        post.append(_state(O, n, v_grain, utils))
        O.engines[3].doneHook = ""
        for _ in range(10):
            O.run(post_iters // 10, wait=True)
            post.append(_state(O, n, v_grain, utils))
    vol = O.cell.volume
    neigh: dict[int, set[int]] = {b.id: set() for b in O.bodies}
    nc = 0
    for i in O.interactions:
        if i.isReal:
            neigh[i.id1].add(i.id2)
            neigh[i.id2].add(i.id1)
            nc += 1
    census = contact_census(O) if dump_path else None
    if census is not None:
        np.savez_compressed(dump_path, **census)
    unb = utils.unbalancedForce()
    s = utils.getStress()
    p = -float(s[0, 0] + s[1, 1] + s[2, 2]) / 3.0
    v_ls = float(np.mean([b.shape.getVolume() for b in O.bodies]))
    return dict(eps=eps, mu_p=mu, p_hat_target=phat, p_hat=p / YOUNG, n=n, seed=seed, r_extent=r,
                v_grain_analytic=v_grain, v_grain_levelset=v_ls, phi=n * v_grain / vol,
                phi_ls_volume=n * v_ls / vol, z_all=2.0 * nc / n, z_nr=z_nonrattler(neigh, 2),
                rattler_frac=1.0 - sum(1 for i in neigh if len(neigh[i]) >= 2) / n,
                unbalanced=float(unb), iters=int(O.iter), converged=bool(unb < packing.UNBALANCED_TOL and O.iter < max_iter),
                dt=float(O.dt), spherical_bodies=sph, spacing_frac=spacing_frac, nodes=nodes, wall_s=time.time() - t0, post=post)


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
    ap.add_argument("--post-iters", type=int, default=0, help="extra stress-controlled iterations after convergence")
    ap.add_argument("--spherical-bodies", type=int, default=int(SPHERICAL_BODIES), help="Law2 flag; 1 only for eps 1")
    ap.add_argument("--census", action="store_true", help="write per-contact geometry to an npz beside the jsonl")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    dump = (os.path.join(OUT, f"contacts_{a.tag}_e{a.eps}_m{a.mu}_n{a.n}_s{a.seed}.npz") if a.census else None)
    if a.spherical_bodies and a.eps != 1.0:
        raise SystemExit("sphericalBodies=True freezes rotation of a non-sphere; refused")
    rec = run(a.eps, a.mu, a.phat, a.n, a.seed, a.spacing, a.nodes, a.max_iter, dump, a.post_iters, bool(a.spherical_bodies))
    rec["tag"] = a.tag
    rec["census"] = os.path.basename(dump) if dump else None
    with open(os.path.join(OUT, f"ls_{a.tag}.jsonl"), "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec, sort_keys=True) + "\n")
    print(json.dumps(rec), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
