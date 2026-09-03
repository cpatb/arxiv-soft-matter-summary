"""diag_cor_timing.py — is COR being read at a physical moment, or at an arbitrary one?

THE SUSPICION
-------------
`phase7_confined.py` measures COR from the centre-of-mass velocity after a fixed probe
loop: `n_steps = max(4 * t_flight / dt, 4000)`, sampled `N_PROBES` times, and COR is
whatever the last sample says. Nothing in that construction guarantees the aggregate has
finished interacting with the wall.

If it has not, COR is a snapshot of an ongoing collision rather than a restitution
coefficient, and its value depends on exactly where the loop happened to stop. That would
show up as large seed scatter -- and Phase 7 does report COR = +0.071 +/- 0.100 and
+0.168 +/- 0.236 at CV = 0.5, standard deviations larger than the means.

It would also masquerade as stiffness dependence: the E = 1 vs E = 10 "flip" that
test_scale_invariance.py reports at CV = 0.5 was measured on n = 2 seeds, against a seed
scatter of the same size as the effect.

WHAT THIS PRINTS
----------------
For a few cells, the full trajectory of centre-of-mass velocity, wall contact count and
kinetic energy over the probe loop, so the question can be answered by looking:

  * does wall contact actually end before the loop does?
  * has v_cm settled to a constant, or is it still changing at the last sample?
  * how much does COR move over the last 20% of the loop?

The last number is the diagnostic. If COR drifts by more than the E-difference under
test, then the "stiffness dependence" is a measurement artefact.
"""

from __future__ import annotations

import math
import sys

import numpy as np

POISSON = 0.30
DT_SAFETY = 0.4
DAMPING = 0.2
SHELL_FRACTION = 0.85
SETTLE_STEPS = 300
N_PROBES = 200          # finer than production, to see the shape


def trace(cv, mu_p, seed, n_parent, carve_p_hat, u0, young=1.0):
    from yade import O, utils

    from phase4_impact import _com_vel, _grain_graph, _ke, _lateral, carve_aggregate

    pos, rad, _phi, _mom = carve_aggregate(cv, mu_p, n_parent, seed, young, carve_p_hat)
    O.reset()
    O.periodic = False
    mat = O.materials.append(
        utils.FrictMat(young=young, poisson=POISSON, density=1.0,
                       frictionAngle=math.atan(mu_p), label="grain"))
    ids = [O.bodies.append(utils.sphere(tuple(p), float(r), material=mat))
           for p, r in zip(pos, rad)]
    r_agg = float(np.linalg.norm(pos, axis=1).max())
    wall_id = O.bodies.append(utils.wall(-r_agg - 1.0, axis=2, sense=1))

    p_cap = carve_p_hat * young
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
        # A confinement must exert ZERO NET FORCE: it squeezes, it does not push. The
        # first version applied an inward force to each shell grain and summed to
        # something non-zero as soon as the shell lost symmetry -- which it does the
        # moment the aggregate flattens against the wall. The result was a continuous
        # thrust: diag_cor_timing.py shows a CV = 0 aggregate ACCELERATING away from the
        # wall with no wall contact at all, COR climbing +0.031 -> +0.339 and ke_ratio
        # 0.0013 -> 0.115 over the probe loop. Every Phase 4/5/7 COR was contaminated by
        # it, and the "COR = 2.15, KE = 175" cells were just the visible extreme.
        #
        # Subtracting the mean makes the force set self-equilibrated, so it can compress
        # the aggregate but cannot translate it.
        dirs = -rel[shell] / rr[shell, None]
        forces = dirs * f
        forces -= forces.mean(axis=0)
        j = 0
        for k, b in enumerate(i):
            if shell[k] and rr[k] > 0:
                O.forces.addF(b, tuple(forces[j]))
                j += 1

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

    O.run(1 + SETTLE_STEPS, wait=True)
    for b in ids:
        O.bodies[b].state.vel = (0.0, 0.0, -u0)
    ke0 = _ke(O, ids)
    d0 = _lateral(O, ids)

    t_flight = (r_agg + 1.0) / u0
    n_steps = max(int(4.0 * t_flight / O.dt), 4000)
    every = max(1, n_steps // N_PROBES)
    idset = set(ids)
    hist = []
    for k in range(N_PROBES):
        O.run(every, wait=True)
        nwall = sum(1 for it in O.interactions
                    if it.isReal and (it.id1 == wall_id or it.id2 == wall_id)
                    and (it.id1 in idset or it.id2 in idset))
        v = _com_vel(O, ids)
        hist.append({"frac": (k + 1) / N_PROBES, "vz": float(v[2]),
                     "cor": float(v[2] / u0), "wall_contacts": nwall,
                     "ke_ratio": _ke(O, ids) / ke0 if ke0 else float("nan"),
                     "spread": _lateral(O, ids) / d0 if d0 else float("nan")})
    _ = _grain_graph
    return hist


def main() -> int:
    cases = [(0.0, 1.0), (0.5, 1.0), (0.5, 10.0)]
    print("Probe-loop trajectories. u0 scaled as sqrt(E) so p_imp/E is fixed.\n")
    for cv, young in cases:
        u0 = 0.02 * math.sqrt(young)
        h = trace(cv, 0.3, 1, 800, 1e-3, u0, young=young)
        last_contact = max((r["frac"] for r in h if r["wall_contacts"] > 0),
                           default=float("nan"))
        tail = [r["cor"] for r in h if r["frac"] > 0.8]
        drift = (max(tail) - min(tail)) if tail else float("nan")
        print(f"  --- CV = {cv:g}, E = {young:g} ---")
        print(f"    last wall contact at   {last_contact:.2f} of the loop")
        print(f"    final COR              {h[-1]['cor']:+.4f}")
        print(f"    COR range over last 20% of loop  {drift:.4f}")
        print(f"    final wall contacts    {h[-1]['wall_contacts']}")
        print("     frac  wall  COR       ke_ratio  spread")
        for r in h[::N_PROBES // 10]:
            print(f"    {r['frac']:5.2f} {r['wall_contacts']:5d} {r['cor']:+8.4f} "
                  f"{r['ke_ratio']:9.4f} {r['spread']:7.3f}")
        print()

    print("Reading: if wall contact persists to the end of the loop, or COR drifts over")
    print("the last 20% by more than the E-difference under test (~0.11), then COR is a")
    print("snapshot of an unfinished collision and the 'stiffness dependence' at CV = 0.5")
    print("is a measurement artefact rather than physics.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
