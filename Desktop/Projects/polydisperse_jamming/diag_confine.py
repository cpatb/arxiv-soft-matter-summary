"""diag_confine.py -- where does the inadmissible kinetic energy come from?

WHY THIS EXISTS
---------------
REPORT 7b/7d: `phase7_confined.py` returns `ke_ratio` up to 1.74 for a PASSIVE body
dropped on a rigid wall, which is impossible -- something is doing work on the aggregate.
The confining `PyRunner` is the only external agent besides the wall (there is no
gravity), so it is the suspect. 7d could only DETECT the artefact with an admissibility
gate; this script tries to LOCATE it, because the CV = 0.5 stiffness bracket cannot be
graded while an unphysical energy source is live.

TWO CANDIDATE MECHANISMS, and they make different predictions
-------------------------------------------------------------
(H1) GEOMETRIC PUMPING. The confining force is recomputed each step as
     f_tot = p_cap * 4*pi*rmax^2 with `rmax` the CURRENT aggregate radius. When the
     aggregate pancakes on the wall `rmax` grows, so the inward force grows, and the
     re-formed ball is squeezed harder than it was released. That is a non-conservative
     force loop: it can pump energy in over a cycle.
     PREDICTION: the cumulative confining work W is positive and comparable to the
     kinetic-energy excess; the excess is shared roughly in proportion to grain mass.

(H2) EQUAL-SHARE FLINGING. `f = p_cap*4*pi*rmax^2 / n_shell` is divided EQUALLY among
     shell grains, so every shell grain feels the same FORCE regardless of its size. A
     surface pressure should give each grain a force proportional to its share of the
     SURFACE, i.e. to r_i^2. At CV = 0 the two rules are identical (all radii equal) and
     the bug is invisible; at CV = 0.5 the smallest grains are a small fraction of the
     mean radius, so they get far too much force per unit area and far more per unit
     MASS. A small rattler in the shell has no contacts to resist it and is simply flung.
     PREDICTION: the kinetic-energy excess concentrates in the SMALLEST radius quartile,
     and the effect is absent at CV = 0. This would also explain why 7d's verdict flips
     happened only at CV = 0.5.

The two are not exclusive. The script measures both channels:
  * W_conf, the cumulative work of the confining force (sum of F.v dt over all grains)
  * KE split by radius quartile, normalised by that quartile's mass share
and prints, per cell, which of H1/H2 the numbers support.

Both weightings are run at both CVs, so the candidate fix ("area") is tested against its
own control (CV = 0, where the two rules must coincide to within seed scatter).

This is a diagnostic: it changes no production file and writes only out/diag_confine.jsonl.
"""

from __future__ import annotations

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
N_PARENT = 800
MU = 0.3
CARVE_P_HAT = 1e-3
BASE_U0 = 0.02


def run_instrumented(cv, seed, young, weight):
    """One confined impact with the confining force instrumented.

    `weight` selects how the total confining force is split over the shell grains:
      "equal" -- the production rule, f_i = f_tot / n_shell
      "area"  -- the surface-pressure rule, f_i proportional to r_i^2
    """
    from yade import O, utils

    from phase4_impact import _com_vel, _grain_graph, _ke, _lateral, carve_aggregate

    t0 = time.time()
    u0 = BASE_U0 * math.sqrt(young)
    pos, rad, phi_parent, mom = carve_aggregate(cv, MU, N_PARENT, seed, young,
                                               CARVE_P_HAT)
    if pos.shape[0] < 50:
        return {"failed": True, "error": "carve gave %d grains" % pos.shape[0]}

    O.reset()
    O.periodic = False
    mat = O.materials.append(
        utils.FrictMat(young=young, poisson=POISSON, density=1.0,
                       frictionAngle=math.atan(MU), label="grain"))
    ids = [O.bodies.append(utils.sphere(tuple(p), float(r), material=mat))
           for p, r in zip(pos, rad)]
    r_agg = float(np.linalg.norm(pos, axis=1).max())
    O.bodies.append(utils.wall(-r_agg - 1.0, axis=2, sense=1))

    radii = np.array([float(O.bodies[i].shape.radius) for i in ids])
    mass = np.array([float(O.bodies[i].state.mass) for i in ids])
    # radius quartiles, fixed once so the bookkeeping is the same at every probe
    qedge = np.percentile(radii, [25, 50, 75])
    qidx = np.digitize(radii, qedge)          # 0..3, 0 = smallest quarter

    p_cap = CARVE_P_HAT * young
    state = {"p": p_cap, "ids": ids, "W": 0.0}

    def confine():
        pc = state["p"]
        if pc <= 0:
            return
        i = state["ids"]
        P = np.array([[O.bodies[b].state.pos[j] for j in range(3)] for b in i])
        V = np.array([[O.bodies[b].state.vel[j] for j in range(3)] for b in i])
        rel = P - P.mean(axis=0)
        rr = np.linalg.norm(rel, axis=1)
        rmax = rr.max()
        shell = rr > SHELL_FRACTION * rmax
        if not shell.any():
            return
        f_tot = pc * 4.0 * math.pi * rmax**2
        if weight == "equal":
            share = np.where(shell, 1.0 / shell.sum(), 0.0)
        elif weight == "area":
            a = np.where(shell, radii**2, 0.0)
            share = a / a.sum()
        elif weight == "netzero":
            # The PRODUCTION rule as fixed on 2026-08-31: equal share, then the mean
            # force subtracted so the set exerts zero NET force. This arm exists to
            # answer the question the other two cannot: the fix was justified by the
            # net-force argument, but H1 is about WORK, and a self-equilibrated force
            # set can still do work on a deforming body. Measuring W_conf under the
            # production rule is the only way to know whether the fix addressed H1 too
            # or merely the translation.
            share = np.where(shell, 1.0 / shell.sum(), 0.0)
        else:
            raise ValueError(weight)
        W = 0.0
        fv = np.zeros((len(i), 3))
        for k in range(len(i)):
            if shell[k] and rr[k] > 0:
                fv[k] = -rel[k] / rr[k] * (f_tot * share[k])
        if weight == "netzero":
            m = shell & (rr > 0)
            if m.any():
                fv[m] -= fv[m].mean(axis=0)
        for k, b in enumerate(i):
            if not shell[k] or rr[k] <= 0:
                continue
            O.forces.addF(b, tuple(fv[k]))
            W += float(fv[k] @ V[k]) * O.dt
        state["W"] += W

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
    state["W"] = 0.0          # start the work budget at the moment of launch

    t_flight = (r_agg + 1.0) / u0
    n_steps = max(int(4.0 * t_flight / O.dt), 4000)
    every = max(1, n_steps // N_PROBES)
    trace = []
    z_peak, cpk = float("nan"), -1
    for _ in range(N_PROBES):
        O.run(every, wait=True)
        c, _n, _f = _grain_graph(O, ids)
        nreal = sum(c.values()) // 2
        if nreal > cpk:
            cpk, z_peak = nreal, 2.0 * nreal / len(ids)
        vel = np.array([[O.bodies[i].state.vel[j] for j in range(3)] for i in ids])
        ke_i = 0.5 * mass * (vel**2).sum(axis=1)
        trace.append((float(ke_i.sum()), state["W"]))

    vel = np.array([[O.bodies[i].state.vel[j] for j in range(3)] for i in ids])
    ke_i = 0.5 * mass * (vel**2).sum(axis=1)
    ke_q = np.array([ke_i[qidx == q].sum() for q in range(4)])
    m_q = np.array([mass[qidx == q].sum() for q in range(4)])
    # KE per unit mass in each radius quartile, relative to the whole body's mean.
    # H2 predicts this is strongly top-heavy toward the SMALLEST quartile (q = 0).
    if ke_i.sum() > 0:
        spec_ke = (ke_q / m_q) / (ke_i.sum() / mass.sum())
    else:
        spec_ke = ke_q * np.nan

    v = _com_vel(O, ids)
    out = {
        "cv": cv, "seed": seed, "young": young, "weight": weight, "u0": u0,
        "n_grains": len(ids), "r_agg": r_agg, "phi_parent": phi_parent,
        "cv_realised": float(mom.get("cv", float("nan"))),
        "z_carved": z_carved, "z_launch": z_launch, "z_retention": retention,
        "z_peak": z_peak,
        "cor": float(v[2] / u0), "ke0": ke0, "ke_end": float(ke_i.sum()),
        "ke_ratio": float(ke_i.sum() / ke0) if ke0 > 0 else float("nan"),
        "W_conf": state["W"],
        # H1's ruler: does the confining work account for the energy excess?
        "W_over_ke0": state["W"] / ke0 if ke0 > 0 else float("nan"),
        "excess_over_ke0": (float(ke_i.sum()) - ke0) / ke0 if ke0 > 0 else float("nan"),
        # H2's ruler: specific KE by radius quartile, smallest first
        "specific_ke_by_rquartile": [float(x) for x in spec_ke],
        "r_quartile_edges": [float(x) for x in qedge],
        "spread": float(_lateral(O, ids) / d0) if d0 > 0 else float("nan"),
        "ke_trace": [t[0] / ke0 for t in trace] if ke0 > 0 else [],
        "W_trace": [t[1] / ke0 for t in trace] if ke0 > 0 else [],
        "wall_s": time.time() - t0,
    }
    print("  cv=%.2f E=%-5g w=%-5s ret=%.2f COR=%+.4f ke_ratio=%.3f W/ke0=%+.3f "
          "excess=%+.3f specKE(q)=%s %.0fs"
          % (cv, young, weight, retention, out["cor"], out["ke_ratio"],
             out["W_over_ke0"], out["excess_over_ke0"],
             np.array2string(spec_ke, precision=2), out["wall_s"]), flush=True)
    return out


def main() -> int:
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, "diag_confine.jsonl")
    recs = []
    # CV = 0 is the control where H2 predicts NO effect (all radii equal).
    for cv in (0.5, 0.0):
        for young in (1.0,):
            for weight in ("equal", "netzero"):
                for seed in (1, 2):
                    try:
                        r = run_instrumented(cv, seed, young, weight)
                    except Exception as exc:      # noqa: BLE001
                        print("FAILED cv=%s E=%s w=%s s=%s: %r"
                              % (cv, young, weight, seed, exc), flush=True)
                        continue
                    recs.append(r)
                    with open(path, "a", encoding="utf-8") as fh:
                        fh.write(json.dumps(r, sort_keys=True) + "\n")
    print("\nwrote %d cells -> %s" % (len(recs), path), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
