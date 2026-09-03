"""phase4_impact.py — a grain aggregate in a DROP-LIKE boundary, driven into a wall.

Attacks the field's own contested question, quoted in OPEN_PROBLEMS_MASTER.md §3:

    "intrinsic bulk transition against boundary-mediated dilatancy, unresolved FOR A
     FREE SURFACE, where no rigid wall supplies confinement"

Three boundary conditions are applied to the SAME carved aggregate, so the comparison
is within-configuration rather than between-preparations:

    free       nothing holds the aggregate together
    capillary  inward radial force on the outer shell = a Laplace pressure 4*gamma/D
    caged      a rigid spherical cage of blocked grains -- the rheometer-like limit

HOW THE AGGREGATE IS BUILT, AND WHY IT IS CARVED
------------------------------------------------
The first version inserted grains randomly into a sphere and compacted them with a
central pull applied from Python between `O.run()` calls. It exploded: COR came back
as 2.6e37. TWO bugs, both worth recording:

  1. `ForceResetter` zeroes O.forces at the start of EVERY timestep, so a force added
     from Python before `O.run(400)` acts for one step out of 400. Any body force in
     Yade must be applied from a `PyRunner` inside the engine list.
  2. Random insertion into a sphere leaves large overlaps, and with no controller to
     absorb them the contact law converts them straight into kinetic energy.

Both are avoided by CARVING: `packing.build_and_compress` already produces a validated,
jammed, non-overlapping periodic pack, so the aggregate is the subset of that pack whose
grain centres fall inside a sphere. This reuses the machinery the Phase-2 anchors test
rather than re-implementing compaction.

DIMENSIONLESS CONTROL
---------------------
In packing units (d10 = 1, rho = 1, E = 1):

    p_imp  = rho U0^2 / E         impact stress over contact stiffness
    p_cap  = P_L / E              capillary confinement over contact stiffness
    Lambda = p_imp / p_cap = We/4

The lab's drops sit at Lambda = 12.5 (U0 = 1 m/s) to 200 (U0 = 4 m/s), so surface
tension is one to two decades too weak to confine grains during impact. The sweep
brackets that and adds both limits.

A LIMITATION THAT MUST BE QUOTED WITH EVERY RESULT
--------------------------------------------------
The granular inertial number is I ~ gammadot*d/sqrt(p/rho_p); with gammadot ~ U0/D_agg
and p ~ rho U0^2 this collapses to I ~ d/D_agg, INDEPENDENT of U0. An aggregate of
N ~ 1000 grains has D_agg ~ 12 d, so I ~ 0.08, while a real 2.5 mm drop of 1 um grains
has d/D0 = 4e-4. **The model is ~200x more inertial per grain than the experiment, and
no affordable N fixes that.** This phase therefore compares BOUNDARY CONDITIONS and CV
at matched I. It does not predict the experiment's regime boundary and no number from
it may be quoted as one.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from dataclasses import asdict, dataclass, field

import numpy as np

from measure import fabric_tensor, force_moments, rattler_report

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

POISSON = 0.30
DT_SAFETY = 0.4
SHELL_FRACTION = 0.85     # grains beyond this fraction of r_max feel the capillary pull
# Reduced pressure of the parent periodic pack. Kept LOW: a pack carved at high
# p_hat carries stored elastic energy that blows the free aggregate apart before it
# ever reaches the wall (at 1e-3 the free case lost every contact during a 2000-step
# settle, Z_initial = 0.00). 1e-4 is just above jamming, so the aggregate is
# coherent but not pre-stressed.
CARVE_P_HAT = 1e-4

# The controlling group for whether an impact fluidises a pack or bounces off it is
# the impact stress against the pack's own confining stress:
#
#     J = p_imp / p_carve = (rho U0^2 / E) / (p / E) = rho U0^2 / p
#
# J << 1  the pack is stiff compared with the impact  -> expect elastic response
# J >> 1  the impact overwhelms the confinement       -> expect fluidisation, spread
#
# The first Phase-4 sweep ran ONE u0 (0.02) and ONE carve pressure (1e-4), i.e. a single
# J = 4, and found no rebound anywhere. That is one point, not a phase diagram.
N_PROBES = 100
SETTLE_STEPS = 300        # short: the free aggregate disintegrates if given longer


@dataclass
class ImpactState:
    cv: float
    mu_p: float
    seed: int
    n_parent: int
    bc: str
    lam: float
    u0: float
    carve_p_hat: float = float("nan")
    j_impact: float = float("nan")   # rho u0^2 / p_carve
    n_grains: int = 0
    cv_realised: float = float("nan")
    d32_realised: float = float("nan")
    r_agg: float = float("nan")
    i_inertial: float = float("nan")
    phi_parent: float = float("nan")
    z_initial: float = float("nan")   # as carved, before any settling
    z_launch: float = float("nan")    # after the brief settle, at launch
    z_peak: float = float("nan")
    z_final: float = float("nan")
    rattler_frac_peak: float = float("nan")
    largest_frac_peak: float = float("nan")
    contacts_peak: int = 0
    fabric_aniso_peak: float = float("nan")
    force_tail_peak: float = float("nan")
    fmax_peak: float = float("nan")
    t_peak: float = float("nan")
    cor: float = float("nan")          # +1 = perfect rebound; see cor_convention
    cor_convention: str = ""
    ke_ratio: float = float("nan")
    spread: float = float("nan")
    exploded: bool = False
    series: list = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    wall_s: float = float("nan")

    def to_json(self) -> str:
        return json.dumps(asdict(self), sort_keys=True)


def carve_aggregate(cv, mu_p, n_parent, seed, young=1.0, carve_p_hat=None):
    """Compress a periodic pack, then return the grains inside a sphere.

    Returns (positions, radii, phi_parent, moments). Uses the SAME validated
    compression as Phase 2 -- no bespoke compaction.
    """
    from yade import O

    from packing import build_and_compress
    from psd import LogNormalPSD

    st = build_and_compress(cv=cv, mu_p=mu_p,
                            p_hat=CARVE_P_HAT if carve_p_hat is None else carve_p_hat,
                            n=n_parent, seed=seed, young=young, verbose=False)
    pos, rad = [], []
    for b in O.bodies:
        if b is None:
            continue
        try:
            rad.append(float(b.shape.radius))
        except AttributeError:
            continue
        p = b.state.pos
        pos.append([p[0], p[1], p[2]])
    pos = np.array(pos)
    rad = np.array(rad)

    # carve the largest sphere that fits inside the periodic cell
    hs = np.array(O.cell.hSize)
    L = float(min(hs[0, 0], hs[1, 1], hs[2, 2]))
    centre = np.array([L / 2.0] * 3)
    r_cut = 0.5 * L * 0.98
    d = np.linalg.norm(pos - centre, axis=1)
    keep = d <= r_cut
    dist = LogNormalPSD(d10=1.0, cv=cv)
    mom = dist.realised_moments(2.0 * rad[keep]) if keep.sum() else {}
    return pos[keep] - centre, rad[keep], st.phi, mom


def run_impact(cv, mu_p, seed, n_parent, bc, lam, u0, young=1.0, carve_p_hat=None):
    """One carved aggregate, one boundary condition, one impact."""
    from yade import O, utils

    t_start = time.time()
    st = ImpactState(cv=cv, mu_p=mu_p, seed=seed, n_parent=n_parent, bc=bc,
                     lam=lam, u0=u0)

    cph = CARVE_P_HAT if carve_p_hat is None else carve_p_hat
    st.carve_p_hat = cph
    st.j_impact = (u0 ** 2 / young) / cph
    pos, rad, phi_parent, mom = carve_aggregate(cv, mu_p, n_parent, seed, young, cph)
    if pos.shape[0] < 50:
        st.notes.append(f"carve produced only {pos.shape[0]} grains")
        return st
    st.phi_parent = phi_parent
    st.n_grains = int(pos.shape[0])
    st.cv_realised = float(mom.get("cv", float("nan")))
    st.d32_realised = float(mom.get("d32", float("nan")))
    r_agg = float(np.linalg.norm(pos, axis=1).max())
    st.r_agg = r_agg
    st.i_inertial = 1.0 / (2.0 * r_agg)

    # --- rebuild as a NON-periodic scene ------------------------------------
    O.reset()
    O.periodic = False
    mat = O.materials.append(
        utils.FrictMat(young=young, poisson=POISSON, density=1.0,
                       frictionAngle=math.atan(mu_p), label="grain"))
    grain_ids = []
    for p, r in zip(pos, rad):
        grain_ids.append(O.bodies.append(
            utils.sphere(tuple(p), float(r), material=mat)))

    cage_ids: list[int] = []
    if bc == "caged":
        # Fibonacci lattice of blocked spheres: a rigid spherical wall, built without
        # depending on pack.regularHexa (which is slow and geometry-sensitive here).
        n_cage = 240
        r_cage_grain = 0.6
        mat_c = O.materials.append(
            utils.FrictMat(young=young * 100, poisson=0.3, density=1.0,
                           frictionAngle=0.5, label="cage"))
        gr = (1.0 + 5.0**0.5) / 2.0
        for k in range(n_cage):
            th = math.acos(1.0 - 2.0 * (k + 0.5) / n_cage)
            ph = 2.0 * math.pi * k / gr
            u = np.array([math.sin(th) * math.cos(ph),
                          math.sin(th) * math.sin(ph), math.cos(th)])
            bid = O.bodies.append(utils.sphere(
                tuple(u * (r_agg + r_cage_grain * 0.9)), r_cage_grain, material=mat_c))
            O.bodies[bid].state.blockedDOFs = "xyzXYZ"
            cage_ids.append(bid)

    wall_z = -r_agg - 1.0
    wall_id = O.bodies.append(utils.wall(wall_z, axis=2, sense=1))

    # --- capillary body force, applied EVERY STEP via PyRunner ---------------
    # p_imp is REDUCED (p/E). The confining force needs an ABSOLUTE stress, so p_cap is
    # multiplied by `young`. Identical at young = 1, which is all production used; see the
    # note in phase7_confined.py and test_scale_invariance.py.
    p_imp = u0**2 / young
    p_cap = ((p_imp / lam) * young
             if (bc == "capillary" and math.isfinite(lam) and lam > 0) else 0.0)

    state = {"p_cap": p_cap, "ids": grain_ids}

    def capillary():
        pc = state["p_cap"]
        if pc <= 0:
            return
        ids = state["ids"]
        P = np.array([[O.bodies[i].state.pos[j] for j in range(3)] for i in ids])
        com = P.mean(axis=0)
        rel = P - com
        rr = np.linalg.norm(rel, axis=1)
        rmax = rr.max()
        shell = rr > SHELL_FRACTION * rmax
        if not shell.any():
            return
        f_per = pc * 4.0 * math.pi * rmax**2 / shell.sum()
        # Zero net force: see the note in phase7_confined.py. An unbalanced shell force
        # translates the aggregate instead of only confining it.
        dirs = -rel[shell] / rr[shell, None]
        forces = dirs * f_per
        forces -= forces.mean(axis=0)
        j = 0
        for k, i in enumerate(ids):
            if shell[k] and rr[k] > 0:
                O.forces.addF(i, tuple(forces[j]))
                j += 1

    O.engines = [
        utils.ForceResetter(),
        utils.InsertionSortCollider([utils.Bo1_Sphere_Aabb(), utils.Bo1_Wall_Aabb()]),
        utils.InteractionLoop(
            [utils.Ig2_Sphere_Sphere_ScGeom(), utils.Ig2_Wall_Sphere_ScGeom()],
            [utils.Ip2_FrictMat_FrictMat_FrictPhys()],
            [utils.Law2_ScGeom_FrictPhys_CundallStrack()],
        ),
        # iterPeriod=1: the force must be re-added every step, since ForceResetter
        # clears it at the start of each one. This is the bug that broke v1.
        utils.PyRunner(command="capillary()", iterPeriod=1, label="cap"),
        utils.NewtonIntegrator(damping=0.2),
    ]
    O.dt = DT_SAFETY * utils.PWaveTimeStep()
    import __main__
    __main__.capillary = capillary

    # --- measure AS CARVED, then launch ---------------------------------------
    # Z is read here, before any settling: the free aggregate disintegrates within
    # a few hundred steps and a post-settle Z would report 0 for every BC and hide
    # the very difference this phase exists to measure.
    O.run(1, wait=True)                     # one step to build the interaction list
    st.z_initial = _z_all(O, grain_ids)
    O.run(SETTLE_STEPS, wait=True)          # brief settle under the active BC
    st.z_launch = _z_all(O, grain_ids)

    for i in grain_ids:
        O.bodies[i].state.vel = (0.0, 0.0, -u0)
    d0_lat = _lateral(O, grain_ids)
    ke0 = _ke(O, grain_ids)

    # --- run through the impact ----------------------------------------------
    t_flight = (r_agg + 1.0 - 0.0) / u0
    n_steps = max(int(4.0 * t_flight / O.dt), 4000)
    every = max(1, n_steps // N_PROBES)
    best = {"contacts": -1}
    series = []
    for _ in range(N_PROBES):
        O.run(every, wait=True)
        snap = _snapshot(O, grain_ids)
        snap["t"] = O.time
        series.append(snap)
        if snap["contacts"] > best["contacts"]:
            best = snap
        if not np.isfinite(snap.get("ke", 0.0)) or snap.get("ke", 0.0) > 1e6 * ke0:
            st.exploded = True
            st.notes.append("KE diverged -- run aborted")
            break

    st.series = series
    st.contacts_peak = int(best.get("contacts", 0))
    st.z_peak = float(best.get("z_all", float("nan")))
    st.rattler_frac_peak = float(best.get("rattler_frac_k2", float("nan")))
    st.largest_frac_peak = float(best.get("largest_frac_k2", float("nan")))
    st.fabric_aniso_peak = float(best.get("aniso", float("nan")))
    st.force_tail_peak = float(best.get("tail_slope", float("nan")))
    st.fmax_peak = float(best.get("fmax", float("nan")))
    st.t_peak = float(best.get("t", float("nan")))
    st.z_final = float(series[-1]["z_all"]) if series else float("nan")

    v = _com_vel(O, grain_ids)
    # COR sign convention. The aggregate is launched at v_z = -u0, so a clean rebound
    # leaves it at v_z = +u0 and the coefficient of restitution is +v_z/u0. The first
    # version used -v_z/u0, which reported a perfect bounce as -1. Caught by
    # test_cor_anchor.py: an undamped elastic sphere came back at -1.0009.
    #
    # `cor_convention` is stamped so that records written by runs already in flight
    # (which lack the field and carry the OLD sign) can be negated on read.
    st.cor = float(v[2] / u0)
    st.cor_convention = "positive_is_rebound"
    st.ke_ratio = float(_ke(O, grain_ids) / ke0) if ke0 > 0 else float("nan")
    lat = _lateral(O, grain_ids)
    st.spread = float(lat / d0_lat) if d0_lat > 0 else float("nan")
    st.wall_s = time.time() - t_start

    print(f"  cv={cv:.2f} bc={bc:9s} lam={lam:8.4g} N={st.n_grains} "
          f"J={st.j_impact:.3g} I={st.i_inertial:.4f} -> Zi={st.z_initial:.2f} Zl={st.z_launch:.2f} Zpk={st.z_peak:.2f} "
          f"rat={st.rattler_frac_peak:.3f} span={st.largest_frac_peak:.3f} "
          f"COR={st.cor:+.3f} spread={st.spread:.3f} "
          f"{'EXPLODED ' if st.exploded else ''}{st.wall_s:.0f}s", flush=True)
    _ = (wall_id, cage_ids)
    return st


# --- helpers ---------------------------------------------------------------------

def _grain_graph(O, ids):  # noqa: E741
    idset = set(ids)
    counts = {i: 0 for i in ids}
    neigh: dict[int, set[int]] = {i: set() for i in ids}
    forces = []
    for it in O.interactions:
        if not it.isReal:
            continue
        a, b = it.id1, it.id2
        if a in idset and b in idset:
            counts[a] += 1
            counts[b] += 1
            neigh[a].add(b)
            neigh[b].add(a)
            try:
                forces.append(float(it.phys.normalForce.norm()))
            except Exception:
                forces.append(float("nan"))
    return counts, neigh, np.asarray(forces, float)


def _z_all(O, ids) -> float:  # noqa: E741
    counts, _, _ = _grain_graph(O, ids)
    return 2.0 * (sum(counts.values()) // 2) / len(counts) if counts else float("nan")


def _positions(O, ids) -> np.ndarray:  # noqa: E741
    return np.array([[O.bodies[i].state.pos[j] for j in range(3)] for i in ids])


def _lateral(O, ids) -> float:  # noqa: E741
    p = _positions(O, ids)
    return float(2.0 * np.percentile(
        np.hypot(p[:, 0] - p[:, 0].mean(), p[:, 1] - p[:, 1].mean()), 95))


def _com_vel(O, ids) -> np.ndarray:  # noqa: E741
    return np.array([[O.bodies[i].state.vel[j] for j in range(3)]
                     for i in ids]).mean(axis=0)


def _ke(O, ids) -> float:  # noqa: E741
    tot = 0.0
    for i in ids:
        b = O.bodies[i]
        v = b.state.vel
        tot += 0.5 * b.state.mass * (v[0]**2 + v[1]**2 + v[2]**2)
    return tot


def _snapshot(O, ids) -> dict:  # noqa: E741
    counts, neigh, forces = _grain_graph(O, ids)
    n_all = len(counts)
    nreal = sum(counts.values()) // 2
    out = {"contacts": nreal,
           "z_all": 2.0 * nreal / n_all if n_all else float("nan"),
           "fmax": float(np.nanmax(forces)) if forces.size else float("nan"),
           "ke": _ke(O, ids)}
    out.update(rattler_report(counts, neigh, thresholds=(2,)))
    fm = force_moments(forces)
    out["tail_slope"] = fm["tail_slope"]
    idset = set(ids)
    normals = [[it.geom.normal[j] for j in range(3)] for it in O.interactions
               if it.isReal and it.id1 in idset and it.id2 in idset]
    out["aniso"] = fabric_tensor(np.array(normals))["aniso"] if normals else float("nan")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--n", type=int, default=1500)
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    cvs = (0.0, 0.3) if args.quick else (0.0, 0.2, 0.35, 0.5)
    seeds = (1,) if args.quick else (1, 2, 3)
    mus = (0.3,)
    bcs = ([("free", float("inf")), ("capillary", 50.0), ("caged", 0.0)]
           if args.quick else
           [("free", float("inf")), ("capillary", 200.0), ("capillary", 50.0),
            ("capillary", 12.5), ("caged", 0.0)])

    suffix = ("_quick" if args.quick else "")
    if args.nshards > 1:
        suffix += f"_shard{args.shard:02d}of{args.nshards:02d}"
    path = os.path.join(OUT_DIR, f"phase4{suffix}.jsonl")

    done = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["key"])
                except Exception:
                    continue

    u0 = 0.02
    combos = [(cv, mu, s, bc, lam)
              for cv in cvs for mu in mus for s in seeds for bc, lam in bcs]
    print(f"phase4: {len(combos)} cells, N_parent={args.n}, shard {args.shard}/"
          f"{args.nshards} -> {path}", flush=True)

    t0 = time.time()
    for idx, (cv, mu, seed, bc, lam) in enumerate(combos):
        if idx % args.nshards != args.shard:
            continue
        key = f"cv{cv:.3f}_mu{mu:.3f}_s{seed}_{bc}_lam{lam:.4g}_N{args.n}"
        if key in done:
            continue
        try:
            st = run_impact(cv, mu, seed, args.n, bc, lam, u0)
            rec = json.loads(st.to_json())
        except Exception as exc:      # noqa: BLE001
            print(f"CELL FAILED {key}: {exc!r}", flush=True)
            rec = {"failed": True, "error": repr(exc), "cv": cv, "mu_p": mu,
                   "seed": seed, "bc": bc, "lam": lam, "n_parent": args.n}
        rec["key"] = key
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
        print(f"  [{idx+1}/{len(combos)}] elapsed {(time.time()-t0)/60:.1f} min",
              flush=True)
    print(f"phase4 done in {(time.time()-t0)/60:.1f} min", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
