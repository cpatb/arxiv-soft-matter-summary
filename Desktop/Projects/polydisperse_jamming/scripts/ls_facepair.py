"""ls_facepair.py -- does the SINGLE-POINT level-set contact hold a flush face-face contact?

Two volume-equivalent rounded cubes (superellipsoid eps), frictionless. A is fixed; B sits above A along z, tilted by
TILT0 about the y axis, and is pushed down by gravity acting on B alone (A's DOFs are blocked). Physically a
frictionless cube loaded through its centre onto a flat face is stable flat: the tilt must decay to ~0.

Compared: Ig2_LevelSet_LevelSet_ScGeom (one contact point = the deepest node; what ls_pack.py uses, because the
multi-point functor does not support periodic cells) against Ig2_LevelSet_LevelSet_MultiScGeom (a contact point at
every penetrating node). Printed: tilt (deg) of B's nearest body axis from z, lateral drift, and the lever arm of the
single contact about B's centre, every CHUNK steps.

    python3 scripts/ls_facepair.py --eps 0.5 --law single
"""
from __future__ import annotations

import argparse
import json
import math

import numpy as np


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--eps", type=float, default=0.5)
    ap.add_argument("--law", choices=["single", "multi"], required=True)
    ap.add_argument("--tilt0", type=float, default=5.0, help="initial tilt, degrees")
    ap.add_argument("--steps", type=int, default=60000)
    ap.add_argument("--nodes", type=int, default=1000)
    ap.add_argument("--load", type=float, default=1e-3, help="target overlap / d for the single law (sets g)")
    ap.add_argument("--spherical", type=int, default=0, help="Law2 sphericalBodies flag (Yade default 1; wrong for non-spheres)")
    a = ap.parse_args()

    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import ls_pack
    from yade import O, Quaternion, Vector3, utils

    eps = a.eps
    r = ls_pack.R_EQ * ((4.0 / 3.0) * math.pi / ls_pack.superellipsoid_volume(1.0, eps)) ** (1.0 / 3.0)
    mat = O.materials.append(utils.FrictMat(young=1.0, poisson=0.3, density=1.0, frictionAngle=0.0))
    half = r  # half-extent along a face normal
    ori_b = Quaternion(Vector3(0, 1, 0), math.radians(a.tilt0))
    A = utils.levelSetBody("superellipsoid", center=Vector3(0, 0, 0), extents=Vector3(r, r, r), epsilons=(eps, eps),
                           spacing=0.05 * 2 * r, nSurfNodes=a.nodes, material=mat)
    # lift B so its lowest corner just touches A's top face at the chosen tilt
    th = math.radians(a.tilt0)
    zlift = half + half * (math.cos(th) + math.sin(th))  # sharp-cube corner height; rounded corners leave a gap
    B = utils.levelSetBody("superellipsoid", center=Vector3(0, 0, zlift), extents=Vector3(r, r, r),
                           epsilons=(eps, eps), spacing=0.05 * 2 * r, nSurfNodes=a.nodes, orientation=ori_b,
                           material=mat)
    ia, ib = O.bodies.append(A), O.bodies.append(B)
    O.bodies[ia].state.blockedDOFs = "xyzXYZ"
    O.bodies[ib].state.blockedDOFs = "xy"  # B may fall and rotate but not slide off A (frictionless)
    kn_single = 2 * r * 1.0 * 0.5  # FrictPhys: kn = 2 E1R1 E2R2/(E1R1+E2R2) = E r for equal bodies
    m_b = O.bodies[ib].state.mass
    g = a.load * 1.0 * kn_single / m_b
    if a.law == "single":
        loop = utils.InteractionLoop([utils.Ig2_LevelSet_LevelSet_ScGeom()], [utils.Ip2_FrictMat_FrictMat_FrictPhys()],
                                     [utils.Law2_ScGeom_FrictPhys_CundallStrack(sphericalBodies=a.spherical)])
    else:
        # per-point stiffness: the face-on overlap engages ~ (nodes on a face); scale so the total is comparable
        loop = utils.InteractionLoop([utils.Ig2_LevelSet_LevelSet_MultiScGeom()],
                                     [utils.Ip2_FrictMat_FrictMat_MultiFrictPhys(kn=kn_single / 50, ks=0.0)],
                                     [utils.Law2_MultiScGeom_MultiFrictPhys_CundallStrack(sphericalBodies=a.spherical)])
    O.engines = [utils.ForceResetter(), utils.InsertionSortCollider([utils.Bo1_LevelSet_Aabb()], verletDist=0.0),
                 loop, utils.NewtonIntegrator(damping=0.4, gravity=(0, 0, -g))]
    O.dt = 0.5 * ls_pack.R_EQ * math.sqrt(1.0 / 1.0) * 0.2
    chunk = a.steps // 20
    out = []
    for _ in range(20):
        O.run(chunk, wait=True)
        st = O.bodies[ib].state
        axes = [st.ori * Vector3(*e) for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
        tilt = math.degrees(math.acos(min(1.0, max(abs(v[2]) for v in axes))))
        npts, arm = 0, float("nan")
        i = O.interactions[ia, ib] if O.interactions.has(ia, ib) else None
        if i is not None and i.isReal:
            if a.law == "single":
                npts = 1
                c = np.array(i.geom.contactPoint)
                arm = float(np.linalg.norm(c[:2] - np.array(st.pos)[:2]))
            else:
                npts = len(i.geom.contacts)
        out.append(dict(it=int(O.iter), tilt_deg=tilt, x=float(st.pos[0]), z=float(st.pos[2]), npts=npts,
                        arm_xy=arm))
    print(json.dumps(dict(eps=eps, law=a.law, spherical=a.spherical, tilt0=a.tilt0, g=g, traj=out)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
