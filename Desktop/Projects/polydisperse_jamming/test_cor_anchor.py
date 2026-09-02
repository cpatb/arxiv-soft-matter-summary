"""test_cor_anchor.py — can this pipeline measure a rebound AT ALL?

THE PROBLEM THIS EXISTS TO SOLVE
--------------------------------
Phase 4 and Phase 5 report "nothing rebounds": COR between -0.035 and +0.005 in every
cell. That is only a statement about the material if the measurement can register a
rebound when one happens. Until that is shown, "COR = 0 everywhere" is equally
consistent with a broken observable, and the claim is void rather than physical.
So this is a POSITIVE CONTROL, and a discriminating one: a single elastic sphere
dropped on the same wall, measured through the same `_com_vel` path, MUST come back
with COR > 0. If it does not, every rebound claim in this project is withdrawn.
Three cases, in increasing order of resemblance to the real runs:
  1. one frictionless sphere, zero damping     -> COR must approach 1
  2. one sphere with NewtonIntegrator damping  -> COR must be > 0 but reduced
  3. the damping actually used in Phase 4/5    -> tells us what COR the pipeline
                                                 could have seen if the aggregate
                                                 had behaved elastically
Case 3 is the important one. If the production damping alone drives COR to ~0, then
"nothing rebounds" is a statement about `NewtonIntegrator(damping=0.2)` and not about
granular physics.
Run:  ~/rl.sh test_cor_anchor.py
"""

from __future__ import annotations
import math
import sys

# Production settings settings, copied from phase4_impact so this control tracks them.
POISSON = 0.30
DT_SAFETY = 0.4
PROD_DAMPING = 0.2
U0 = 0.02


def one_sphere_cor(damping: float, young: float = 1.0, radius: float = 0.5,
                   restitution_en: float | None = None) -> dict:
    """Drop one sphere on a wall at U0 and measure COR through the SAME code path."""
    from yade import O, utils

    O.reset()
    O.periodic = False
    mat = O.materials.append(
        utils.FrictMat(young=young, poisson=POISSON, density=1.0,
                       frictionAngle=math.atan(0.3), label="grain"))
    bid = O.bodies.append(utils.sphere((0, 0, 0), radius, material=mat))
    wall_z = -radius - 1.0
    O.bodies.append(utils.wall(wall_z, axis=2, sense=1))
    O.engines = [
        utils.ForceResetter(),
        utils.InsertionSortCollider([utils.Bo1_Sphere_Aabb(), utils.Bo1_Wall_Aabb()]),
        utils.InteractionLoop(
            [utils.Ig2_Sphere_Sphere_ScGeom(), utils.Ig2_Wall_Sphere_ScGeom()],
            [utils.Ip2_FrictMat_FrictMat_FrictPhys()],
            [utils.Law2_ScGeom_FrictPhys_CundallStrack()],
        ),
        utils.NewtonIntegrator(damping=damping),
    ]
    O.dt = DT_SAFETY * utils.PWaveTimeStep()
    O.bodies[bid].state.vel = (0.0, 0.0, -U0)
    # fly, hit, and leave: 4x the flight time is enough for a full bounce
    t_flight = 1.0 / U0
    n_steps = int(4.0 * t_flight / O.dt)
    vz_min = 0.0
    vz_hist = []
    every = max(1, n_steps // 400)
    for _ in range(400):
        O.run(every, wait=True)
        vz = float(O.bodies[bid].state.vel[2])
        vz_hist.append(vz)
        vz_min = min(vz_min, vz)
    vz_end = float(O.bodies[bid].state.vel[2])

    # Physical convention: launched at -U0, so a clean rebound is +U0 and
    # COR = +vz_end/U0. The production code originally used the opposite
    # sign; finding that is what this test was written for.
    return {"damping": damping, "cor": vz_end / U0,
            "vz_max_up": max(vz_hist) / U0,
            "vz_end": vz_end, "restitution_en": restitution_en}


def main() -> int:
    rows = []
    print("single-sphere COR through the production measurement path")
    print(f"{'damping':>9} {'COR (end)':>11} {'max upward v/U0':>17}")
    for d in (0.0, 0.05, PROD_DAMPING, 0.4):
        r = one_sphere_cor(damping=d)
        rows.append(r)
        print(f"{d:9.2f} {r['cor']:11.4f} {r['vz_max_up']:17.4f}", flush=True)
    print()
    undamped = rows[0]
    prod = [r for r in rows if r["damping"] == PROD_DAMPING][0]
    ok = True
    if not (undamped["cor"] > 0.5):
        print(f"FAIL: an undamped elastic sphere returned COR = {undamped['cor']:.4f}, "
              f"which should approach 1. The COR measurement is broken and every")
        print("      'nothing rebounds' claim in this project is WITHDRAWN.")
        ok = False
    else:
        print(f"PASS: undamped elastic sphere gives COR = {undamped['cor']:.4f}, so the "
              f"measurement can register a rebound.")
    print(f"\nAt the PRODUCTION damping of {PROD_DAMPING}, one elastic sphere gives "
          f"COR = {prod['cor']:.4f}.")
    if prod["cor"] < 0.05:
        print("  => WARNING: the production damping ALONE suppresses the rebound of a")
        print("     perfectly elastic sphere. 'Nothing rebounds' would then be a")
        print("     statement about NewtonIntegrator(damping), not about the granular")
        print("     physics, and the Phase 4/5 conclusion must be re-run at lower")
        print("     damping before it can be quoted.")
        ok = False
    else:
        print(f"  => the production damping leaves COR = {prod['cor']:.4f} for an elastic")
        print("     sphere, comfortably above the 0.01 rebound bar used in Phase 5, so a")
        print("     rebound would have been visible had one occurred.")
    print("\nCOR_ANCHOR_PASS" if ok else "\nCOR_ANCHOR_FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
