"""diag_rmax_dynamic.py — does the rmax estimator inflate DURING the impact?

WHAT §7h LEFT OPEN
------------------
§7h screened the `rmax` hypothesis geometrically: on the CARVED aggregate at rest, the
force inflation `(rmax/r_p95)^2` is 1.07 at CV = 0 and 1.09 at CV = 0.5, a ratio of 1.02
against a pre-registered 1.5 bar. That clears the estimator AT LAUNCH.

But the confining force recomputes `rmax` EVERY STEP on an aggregate that pancakes against
the wall, and a deforming polydisperse aggregate could develop a larger `rmax/r_p95` than a
deforming monodisperse one. If it does, the estimator could inflate the WORK without
inflating the static ratio, and §7h's verdict would not survive.

So: log the radius estimators through the probe loop and compare.

WHAT IS MEASURED
----------------
At each probe, on the live aggregate:
    rmax        max |r_i - r_cm|           -- what the force actually uses
    r_p95       95th percentile of |r_i - r_cm|
    inflation   (rmax / r_p95)^2           -- the force-relevant quantity
plus the wall contact count, so the pancaked phase can be identified rather than assumed.

PRE-REGISTERED READING (fixed before the run)
---------------------------------------------
Let G_dyn = max-over-loop inflation(CV = 0.5) / max-over-loop inflation(CV = 0), using the
same 1.5 bar §7h used:

  G_dyn <= 1.5   the estimator is not differentially inflated during impact either.
                 §7h's verdict extends from launch to the whole trajectory.
  G_dyn >= 3.0   the estimator DOES inflate differentially once the aggregate deforms.
                 §7h must be withdrawn and the confinement re-run with a percentile radius.
  otherwise      partial; report the number and do not claim either way.

Reporting the WORST (max over the loop) rather than the mean is deliberate: the hypothesis
is about a transient geometric excursion during the pancake, which a mean would dilute.
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
N_PROBES = 120
MU = 0.3
CARVE_P_HAT = 1e-3
N_PARENT = 800
U0 = 0.02
BAR_SAME, BAR_DIFFER = 1.5, 3.0


def trace(cv, seed, young=1.0):
    from yade import O, utils

    from phase4_impact import _com_vel, _ke, _lateral, carve_aggregate

    t0 = time.time()
    pos, rad, _phi, _mom = carve_aggregate(cv, MU, N_PARENT, seed, young, CARVE_P_HAT)
    O.reset()
    O.periodic = False
    mat = O.materials.append(
        utils.FrictMat(young=young, poisson=POISSON, density=1.0,
                       frictionAngle=math.atan(MU), label="grain"))
    ids = [O.bodies.append(utils.sphere(tuple(p), float(r), material=mat))
           for p, r in zip(pos, rad)]
    r_agg = float(np.linalg.norm(pos, axis=1).max())
    wall_id = O.bodies.append(utils.wall(-r_agg - 1.0, axis=2, sense=1))

    p_cap = CARVE_P_HAT * young
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
        dirs = -rel[shell] / rr[shell, None]
        forces = dirs * f
        forces -= forces.mean(axis=0)      # the production net-zero rule
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
        O.bodies[b].state.vel = (0.0, 0.0, -U0)
    ke0 = _ke(O, ids)
    d0 = _lateral(O, ids)

    t_flight = (r_agg + 1.0) / U0
    n_steps = max(int(4.0 * t_flight / O.dt), 4000)
    every = max(1, n_steps // N_PROBES)
    idset = set(ids)
    hist = []
    for k in range(N_PROBES):
        O.run(every, wait=True)
        P = np.array([[O.bodies[b].state.pos[j] for j in range(3)] for b in ids])
        rel = P - P.mean(axis=0)
        rr = np.linalg.norm(rel, axis=1)
        rmax = float(rr.max())
        rp95 = float(np.percentile(rr, 95))
        nwall = sum(1 for it in O.interactions
                    if it.isReal and (it.id1 == wall_id or it.id2 == wall_id)
                    and (it.id1 in idset or it.id2 in idset))
        hist.append({"frac": (k + 1) / N_PROBES, "rmax": rmax, "r_p95": rp95,
                     "inflation": (rmax / rp95) ** 2 if rp95 > 0 else float("nan"),
                     "wall_contacts": nwall,
                     "spread": _lateral(O, ids) / d0 if d0 else float("nan")})
    v = _com_vel(O, ids)
    infl = np.array([h["inflation"] for h in hist], float)
    out = {"cv": cv, "seed": seed, "young": young, "n_grains": len(ids),
           "infl_launch": float(infl[0]), "infl_max": float(np.nanmax(infl)),
           "infl_mean": float(np.nanmean(infl)),
           "infl_at_peak_contact": float(
               infl[int(np.argmax([h["wall_contacts"] for h in hist]))]),
           "max_wall_contacts": int(max(h["wall_contacts"] for h in hist)),
           "cor": float(v[2] / U0),
           "ke_ratio": float(_ke(O, ids) / ke0) if ke0 > 0 else float("nan"),
           "wall_s": time.time() - t0, "hist": hist}
    print("  cv=%.2f s=%d  infl launch %.4f  max %.4f  at-peak-contact %.4f  (%.0fs)"
          % (cv, seed, out["infl_launch"], out["infl_max"],
             out["infl_at_peak_contact"], out["wall_s"]), flush=True)
    return out


def main() -> int:
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, "diag_rmax_dynamic.jsonl")
    recs = []
    for cv in (0.0, 0.5):
        for seed in (1, 2):
            try:
                r = trace(cv, seed)
            except Exception as exc:              # noqa: BLE001
                print("FAILED cv=%s seed=%s: %r" % (cv, seed, exc), flush=True)
                continue
            recs.append(r)
            with open(path, "a", encoding="utf-8") as fh:
                fh.write(json.dumps(r, sort_keys=True) + "\n")

    if not recs:
        print("no cells")
        return 1

    lines = ["=" * 78,
             "DOES THE rmax ESTIMATOR INFLATE DURING THE IMPACT?",
             "=" * 78, "",
             "inflation = (rmax / r_p95)^2, the force-relevant quantity, logged through the",
             "probe loop. 7h measured only the launch value.", ""]
    lines.append("  %5s %3s %11s %10s %10s %8s"
                 % ("CV", "n", "at launch", "max", "at contact", "COR"))
    per = {}
    for cv in (0.0, 0.5):
        g = [r for r in recs if r["cv"] == cv]
        if not g:
            continue
        per[cv] = {k: float(np.mean([r[k] for r in g]))
                   for k in ("infl_launch", "infl_max", "infl_at_peak_contact", "cor")}
        lines.append("  %5.2f %3d %11.4f %10.4f %10.4f %+8.4f"
                     % (cv, len(g), per[cv]["infl_launch"], per[cv]["infl_max"],
                        per[cv]["infl_at_peak_contact"], per[cv]["cor"]))
    lines.append("")
    tag = "INCONCLUSIVE"
    if 0.0 in per and 0.5 in per:
        g_dyn = per[0.5]["infl_max"] / per[0.0]["infl_max"]
        g_con = per[0.5]["infl_at_peak_contact"] / per[0.0]["infl_at_peak_contact"]
        lines.append("PRE-REGISTERED: G_dyn <= %.1f -> 7h extends; >= %.1f -> 7h withdrawn"
                     % (BAR_SAME, BAR_DIFFER))
        lines.append("  G_dyn (max over loop)      = %.3f" % g_dyn)
        lines.append("  G_dyn (at peak wall contact) = %.3f" % g_con)
        lines.append("")
        if g_dyn <= BAR_SAME and g_con <= BAR_SAME:
            tag = "EXTENDS"
            lines.append("  => 7h EXTENDS from launch to the whole trajectory. The estimator")
            lines.append("     is not differentially inflated once the aggregate deforms, so")
            lines.append("     the CV work gap remains attributable to the packing.")
        elif g_dyn >= BAR_DIFFER or g_con >= BAR_DIFFER:
            tag = "WITHDRAWN"
            lines.append("  => 7h MUST BE WITHDRAWN. The estimator inflates differentially")
            lines.append("     during the pancake, so it can raise the work at CV = 0.5")
            lines.append("     without showing up in the static screen. Re-run the")
            lines.append("     confinement with a percentile radius.")
        else:
            tag = "PARTIAL"
            lines.append("  => PARTIAL. Between the bars; report the number, claim neither.")
        lines.append("")
        lines.append("  RMAX_DYNAMIC_" + tag)
    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(OUT_DIR, "rmax_dynamic_verdict.txt"), "w",
              encoding="utf-8") as fh:
        fh.write(txt + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
