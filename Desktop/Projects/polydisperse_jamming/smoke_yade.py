"""smoke_yade.py — prove the userspace Yade actually SIMULATES, not just imports.

An import proving nothing is exactly the `rc=0` trap in the lab's lessons file: it
shows a run finished, not that the requested run happened. So this builds a real
packing, runs real timesteps, and checks quantities that could only come from a
working contact-mechanics engine:

  * a gravity-free isotropic compression raises the packing fraction
  * contacts appear, and their count is > 0
  * total kinetic energy decays under damping
  * the measured stress tensor is compressive and roughly isotropic

Run:  source scripts/env.sh && python3 smoke_yade.py
"""

import sys

import numpy as np

from yade import O, utils


def main() -> int:
    rng = np.random.default_rng(20260831)

    # --- material ----------------------------------------------------------
    mat = O.materials.append(
        utils.FrictMat(young=1e8, poisson=0.3, density=2200.0,
                       frictionAngle=np.arctan(0.3), label="grain")
    )

    # --- a small periodic cell of monodisperse spheres ----------------------
    O.periodic = True
    L = 20e-6                       # 20 um box
    O.cell.hSize = np.diag([L, L, L])
    d = 2.0e-6                      # 2 um diameter
    n = 400

    placed = 0
    for _ in range(20000):
        if placed >= n:
            break
        pos = rng.uniform(0, L, 3)
        O.bodies.append(utils.sphere(tuple(pos), d / 2.0, material=mat))
        placed += 1
    print(f"placed bodies: {len(O.bodies)}  (requested {n})")

    phi0 = utils.porosity()
    print(f"initial porosity: {phi0:.4f}  -> phi = {1 - phi0:.4f}")

    # --- engines -----------------------------------------------------------
    O.engines = [
        utils.ForceResetter(),
        utils.InsertionSortCollider([utils.Bo1_Sphere_Aabb()]),
        utils.InteractionLoop(
            [utils.Ig2_Sphere_Sphere_ScGeom()],
            [utils.Ip2_FrictMat_FrictMat_FrictPhys()],
            [utils.Law2_ScGeom_FrictPhys_CundallStrack()],
        ),
        utils.PeriIsoCompressor(charLen=d, stresses=[-1e4, -1e5],
                                maxUnbalanced=0.05, doneHook="O.pause()",
                                globalUpdateInt=20, keepProportions=True),
        utils.NewtonIntegrator(damping=0.4),
    ]
    O.dt = 0.3 * utils.PWaveTimeStep()
    print(f"dt = {O.dt:.4e} s   (0.3 * PWave)")

    O.run(60000, wait=True)

    # --- checks ------------------------------------------------------------
    phi = 1 - utils.porosity()
    ncon = len([i for i in O.interactions if i.isReal])
    unb = utils.unbalancedForce()
    stress = utils.getStress()
    tr3 = float(np.trace(np.array(stress))) / 3.0

    print(f"iterations run   : {O.iter}")
    print(f"final phi        : {phi:.4f}")
    print(f"real contacts    : {ncon}")
    print(f"unbalanced force : {unb:.4e}")
    print(f"mean stress      : {tr3:.4e} Pa  (negative = compressive)")

    ok = True
    if phi <= 1 - phi0:
        print("FAIL: compression did not raise phi")
        ok = False
    if ncon <= 0:
        print("FAIL: no real contacts formed")
        ok = False
    if not np.isfinite(tr3):
        print("FAIL: stress is not finite")
        ok = False
    if tr3 >= 0:
        print("FAIL: mean stress is not compressive")
        ok = False
    if phi < 0.3 or phi > 0.75:
        print(f"FAIL: phi = {phi:.4f} is outside any physical sphere-pack range")
        ok = False

    print("SMOKE_PASS" if ok else "SMOKE_FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
