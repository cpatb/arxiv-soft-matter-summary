"""ls_stiffness_probe.py -- what normal stiffness does a level-set contact get? (PREREG_LSPILOT P0 diagnosis)

Two grains pressed together along x by a small overlap, one step, then print the ScGeom radii and FrictPhys kn for:
FrictMat spheres (reference, kn = E * 2 r1 r2 / (r1 + r2) = E r for equal r), level-set spheres (eps 1), level-set
rounded cubes (eps 0.5) face-to-face and corner-first. Volume-equivalent to d = 1 throughout.
"""
import math

from yade import O, Quaternion, Vector3, utils

import ls_pack


def probe(kind, eps=1.0, rot=None):
    O.reset()
    mat = O.materials.append(utils.FrictMat(young=1.0, poisson=0.3, density=1.0, frictionAngle=0.3))
    if kind == "sphere":
        O.bodies.append([utils.sphere((0, 0, 0), 0.5, material=mat), utils.sphere((0.99, 0, 0), 0.5, material=mat)])
        ig, bo = utils.Ig2_Sphere_Sphere_ScGeom(), utils.Bo1_Sphere_Aabb()
        gap = 0.99
    else:
        r = ls_pack.R_EQ * ((4.0 / 3.0) * math.pi / ls_pack.superellipsoid_volume(1.0, eps)) ** (1.0 / 3.0)
        q = rot or Quaternion((1, 0, 0), 0)
        gap = 2 * r - 0.01 if rot is None else None
        if rot is not None:   # corner-first: rotate so a body diagonal points along x; place at contact by trial
            gap = 2 * r * 1.3
        for xc in (0.0, gap):
            O.bodies.append(utils.levelSetBody("superellipsoid", center=Vector3(xc, 0, 0), extents=Vector3(r, r, r),
                                               epsilons=(eps, eps), spacing=0.05 * 2 * r, orientation=q, material=mat))
        ig, bo = utils.Ig2_LevelSet_LevelSet_ScGeom(), utils.Bo1_LevelSet_Aabb()
    O.engines = [utils.ForceResetter(), utils.InsertionSortCollider([bo], verletDist=0.0),
                 utils.InteractionLoop([ig], [utils.Ip2_FrictMat_FrictMat_FrictPhys()],
                                       [utils.Law2_ScGeom_FrictPhys_CundallStrack()]), utils.NewtonIntegrator(damping=0.4)]
    O.dt = 1e-6
    O.step()
    for i in O.interactions:
        if i.isReal:
            g = i.geom
            print(f"{kind:8s} eps {eps}: refR1 {g.refR1:.4f} refR2 {g.refR2:.4f} penetration {g.penetrationDepth:.5f} "
                  f"kn {i.phys.kn:.4f}  (sphere reference kn = E r = 0.5)")
            return
    print(f"{kind} eps {eps}: no real interaction (gap {gap})")


probe("sphere")
probe("ls", 1.0)
probe("ls", 0.5)
probe("ls", 0.3)
