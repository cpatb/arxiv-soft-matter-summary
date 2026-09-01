"""packing.py — generate and compress polydisperse sphere packings in Yade.

Runs inside WSL under `scripts/env.sh`. `yade` is imported at call time so the
analysis and PSD layers stay importable on Windows without it.

WHY THIS WORKS IN DIMENSIONLESS UNITS
-------------------------------------
The first version ran at the physical scale (d10 = 2 um, E = 1e8 Pa, rho = 2200) and
took >3 minutes of 4-thread wall clock for a single N = 300 cell, because
`PWaveTimeStep` at micron size gives dt ~ 1.4e-9 s and equilibration then needs
~1e6 steps. That is a 150-cell sweep measured in days.

It is also unnecessary. For frictional hard spheres the packing observables this
project measures -- phi_J, Z, the rattler fraction, the fabric, the force
distribution -- are INVARIANT under a change of length, mass and time units. They
depend only on the dimensionless group

    p_hat = p / E        (reduced pressure; how soft the contacts are)

together with CV, mu_p and N. So the model is written with

    d10 = 1, rho = 1, E = 1

and `p_hat` is set directly as the compression target. Physical units are restored
only at the reporting boundary, and only for quantities that need them.

This is not a convenience: making `p_hat` the control parameter is what turns the
stiffness sweep demanded by SPEC §2.3 into a one-line sweep, and it removes the
temptation to quote a Young's modulus in Pa as though it were physical when the real
silica value (7e10 Pa) is unreachable by DEM anyway.

PROTOCOL, fixed here and not varied between cells
-------------------------------------------------
Frictional jamming points are preparation-dependent, so the preparation IS part of
the measurement:

  1. Random insertion of N spheres from a lognormal PSD at fixed NUMBER mean into a
     periodic cell sized for phi ~ 0.30. Insertion overlaps are relaxed away.
  2. Isotropic STRESS-controlled compression (PeriTriaxController, stressMask=7) to a
     target reduced pressure `p_hat`. Stress control is used rather than strain
     control because it converges to a well-defined mechanical state directly, and
     because p_hat is the physical control parameter.
  3. Interparticle friction is held at its production value throughout. No friction
     annealing -- annealing shifts phi_J by a few percent and the protocol would then
     be a hidden parameter.
  4. Read phi, pressure, Z, and the rattler set at the converged state.

phi_J is obtained in `phase2_phiJ.py` by running a LADDER of p_hat and extrapolating
phi(p_hat) -> p_hat = 0. It is never read off a single run.

CONVENTIONS
-----------
  * `d` is a DIAMETER everywhere; Yade wants radii. The conversion happens only in
    `psd.LogNormalPSD.radii_for_yade`.
  * Yade's `getStress()` is positive in TENSION. `pressure` here is -trace/3, so it
    is POSITIVE in compression, the granular convention. Getting this backwards
    silently inverts every phi_J fit.
  * `phi` INCLUDES rattlers (they occupy volume). `Z` EXCLUDES them (they carry no
    force). Mixing the two is the standard way to get Z != 6 at the frictionless
    isostatic point and then "explain" it.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field

import numpy as np

from measure import (
    RATTLER_MIN_CONTACTS,
    force_moments,
    rattler_report,
)
from psd import make_psd

# --- protocol constants, with their origin ---------------------------------------
UNBALANCED_TOL = 0.02     # Yade quasi-static criterion. Swept in test_packing.py;
                          # O'Hern-style studies use 1e-2..1e-3 on this ratio.
REL_STRESS_TOL = 2e-3     # PeriTriaxController's own stress-convergence bar
PHI_INSERT = 0.30
DT_SAFETY = 0.5
MAX_ITER = 200_000
POISSON = 0.30

# Periodic-collider geometry. Yade requires cell edge > 2.1 * maxBodySpan; 2.5 leaves
# margin. PHI_DENSE_BOUND is an upper bound on any phi this project reaches (frictionless
# monodisperse RCP is 0.64; nothing here should exceed 0.72), used to check the
# constraint at the DENSEST state the cell will pass through.
COLLIDER_CLEARANCE = 2.5
PHI_DENSE_BOUND = 0.72



@dataclass
class PackingState:
    """Everything measured from one relaxed packing. JSON-serialisable."""
    # inputs (dimensionless)
    n_requested: int
    cv: float
    mu_p: float
    p_hat_target: float
    seed: int
    family: str = "lognormal"
    # realised PSD -- NEVER assume the analytic moments survive truncation
    n_placed: int = 0
    d10_realised: float = float("nan")
    d32_realised: float = float("nan")
    d43_realised: float = float("nan")
    cv_realised: float = float("nan")
    d_span: float = float("nan")
    d_max_over_L: float = float("nan")   # finite-size diagnostic; >0.25 is a risk
    # state
    phi: float = float("nan")
    p_hat: float = float("nan")             # realised p/E, positive in compression
    stress_dev_hat: float = float("nan")    # von Mises deviator / E
    z_all: float = float("nan")
    z_nr: float = float("nan")              # over NON-RATTLERS only
    rattler_frac: float = float("nan")
    n_contacts: int = 0
    unbalanced: float = float("nan")
    converged: bool = False
    iters: int = 0
    wall_s: float = float("nan")
    rattler: dict = field(default_factory=dict)   # every threshold, see measure
    force: dict = field(default_factory=dict)     # contact-force moments
    protocol: str = "iso-stress-nondim-v2"
    notes: list[str] = field(default_factory=list)

    def to_json(self) -> str:
        return json.dumps(asdict(self), sort_keys=True)


def contact_graph(O):  # noqa: E741
    """(counts, neighbours, forces) over the REAL interactions of a Yade scene."""
    counts: dict[int, int] = {b.id: 0 for b in O.bodies if b is not None}
    neigh: dict[int, set[int]] = {b.id: set() for b in O.bodies if b is not None}
    forces: list[float] = []
    for i in O.interactions:
        if not i.isReal:
            continue
        counts[i.id1] += 1
        counts[i.id2] += 1
        neigh[i.id1].add(i.id2)
        neigh[i.id2].add(i.id1)
        try:
            fn = float(i.phys.normalForce.norm())
        except Exception:
            fn = float("nan")
        forces.append(fn)
    return counts, neigh, np.asarray(forces, dtype=float)


def measure(O, utils, st: PackingState, young: float) -> PackingState:  # noqa: E741
    """Read every observable off a relaxed Yade scene into `st`."""
    counts, neigh, forces = contact_graph(O)
    n_all = len(counts)
    nreal = sum(counts.values()) // 2
    st.n_contacts = nreal
    st.z_all = 2.0 * nreal / n_all if n_all else float("nan")

    # The rattler threshold is a definition that changes the answer, so all three are
    # recorded and none is privileged in the data. See measure.RATTLER_THRESHOLDS.
    st.rattler = rattler_report(counts, neigh)
    st.rattler_frac = st.rattler[f"rattler_frac_k{RATTLER_MIN_CONTACTS}"]
    st.z_nr = st.rattler[f"z_nr_k{RATTLER_MIN_CONTACTS}"]
    st.force = force_moments(forces)

    s = np.array(utils.getStress())
    st.p_hat = float(-np.trace(s) / 3.0 / young)
    dev = s - np.eye(3) * (np.trace(s) / 3.0)
    st.stress_dev_hat = float(np.sqrt(1.5 * np.sum(dev * dev)) / young)

    st.phi = float(1.0 - utils.porosity())
    st.unbalanced = float(utils.unbalancedForce())
    st.iters = int(O.iter)
    return st


def build_and_compress(cv: float, mu_p: float, p_hat: float, n: int, seed: int,
                       truncate: tuple[float, float] | None = (0.2, 5.0),
                       young: float = 1.0, verbose: bool = True,
                       family: str = "lognormal") -> PackingState:
    """Compress a packing isotropically to reduced pressure `p_hat` = p/E.

    Dimensionless throughout: d10 = 1, rho = 1, E = `young` (default 1).
    """
    import time

    from yade import O, utils

    t0 = time.time()
    O.reset()
    st = PackingState(n_requested=n, cv=cv, mu_p=mu_p, p_hat_target=p_hat,
                      seed=seed, family=family)

    rng = np.random.default_rng(seed)
    dist = make_psd(family, 1.0, cv)
    radii = dist.radii_for_yade(n, rng, truncate=truncate)
    diam = 2.0 * radii

    m = dist.realised_moments(diam)
    st.n_placed = int(m["n"])
    st.d10_realised = m["d10"]
    st.d32_realised = m["d32"]
    st.d43_realised = m["d43"]
    st.cv_realised = m["cv"]
    st.d_span = m["span"]

    vol_solid = float(np.sum(4.0 / 3.0 * math.pi * radii**3))
    d_max = float(diam.max())

    # Yade's periodic collider requires cell edge > 2.1 * maxBodySpan, and the cell
    # SHRINKS under compression, so the binding check is at the FINAL density, not the
    # initial one. Sizing only the insertion cell was not enough: N=200/CV=0.5 passed
    # at phi=0.30 and then died on the way to phi~0.64.
    #
    # Enlarging the cell is not available either -- the final cell edge is fixed by
    # vol_solid and the jammed phi -- so this is a constraint on N, and it is
    # reported as one rather than worked around by clipping the PSD (clipping would
    # silently move d32, the one moment this project turns on).
    L_final_min = (vol_solid / PHI_DENSE_BOUND) ** (1.0 / 3.0)
    if L_final_min < COLLIDER_CLEARANCE * d_max:
        v_mean = vol_solid / len(radii)
        n_min = math.ceil(PHI_DENSE_BOUND * (COLLIDER_CLEARANCE * d_max) ** 3 / v_mean)
        raise ValueError(
            f"cell too small for the periodic collider at cv={cv}, N={n}, "
            f"truncate={truncate}: largest grain d_max={d_max:.3f} needs a final "
            f"cell edge >= {COLLIDER_CLEARANCE * d_max:.3f} but phi={PHI_DENSE_BOUND} "
            f"gives {L_final_min:.3f}. Use N >= {n_min}, or tighten `truncate`."
        )

    L = (vol_solid / PHI_INSERT) ** (1.0 / 3.0)
    st.d_max_over_L = d_max / L_final_min      # worst case, at the densest state
    # A grain wider than ~1/4 of the cell interacts with its own periodic image through
    # very few intermediaries; flag it rather than trusting the number.
    if st.d_max_over_L > 0.25:
        st.notes.append(
            f"FINITE SIZE RISK: d_max/L_dense = {st.d_max_over_L:.3f} > 0.25"
        )

    mat = O.materials.append(
        utils.FrictMat(young=young, poisson=POISSON, density=1.0,
                       frictionAngle=math.atan(mu_p), label="grain")
    )
    O.periodic = True
    O.cell.hSize = np.diag([L, L, L])
    for r in radii:
        pos = rng.uniform(0, L, 3)
        O.bodies.append(utils.sphere(tuple(pos), float(r), material=mat))

    # Yade sign convention for the triax goal: compressive stress is NEGATIVE.
    goal = -abs(p_hat) * young
    O.engines = [
        utils.ForceResetter(),
        utils.InsertionSortCollider([utils.Bo1_Sphere_Aabb()],
                                    allowBiggerThanPeriod=False),
        utils.InteractionLoop(
            [utils.Ig2_Sphere_Sphere_ScGeom()],
            [utils.Ip2_FrictMat_FrictMat_FrictPhys()],
            [utils.Law2_ScGeom_FrictPhys_CundallStrack()],
        ),
        utils.PeriTriaxController(
            label="triax",
            goal=(goal, goal, goal),
            stressMask=7,            # all three goals are STRESSES
            dynCell=True,
            maxUnbalanced=UNBALANCED_TOL,
            relStressTol=REL_STRESS_TOL,
            growDamping=0.4,
            doneHook="O.pause()",
        ),
        utils.NewtonIntegrator(damping=0.4),
    ]
    O.dt = DT_SAFETY * utils.PWaveTimeStep()
    O.run(MAX_ITER, wait=True)

    st = measure(O, utils, st, young)
    st.wall_s = time.time() - t0
    st.converged = bool(st.unbalanced < UNBALANCED_TOL and O.iter < MAX_ITER)
    if not st.converged:
        st.notes.append(
            f"NOT CONVERGED: unbalanced={st.unbalanced:.3e} tol={UNBALANCED_TOL} "
            f"iters={st.iters}/{MAX_ITER}"
        )
    rel_p = abs(st.p_hat - p_hat) / p_hat if p_hat > 0 else float("nan")
    if np.isfinite(rel_p) and rel_p > 0.10:
        st.notes.append(f"pressure off target by {rel_p:.1%} "
                        f"(realised {st.p_hat:.3e} vs {p_hat:.3e})")
    if verbose:
        r = st.rattler
        print(f"[{family[:4]}] cv={cv:.2f} mu={mu_p:.2f} phat_t={p_hat:.1e} N={n} s={seed} -> "
              f"phi={st.phi:.4f} phat={st.p_hat:.3e} Zall={st.z_all:.3f} "
              f"rat(k2/3/4)={r['rattler_frac_k2']:.3f}/{r['rattler_frac_k3']:.3f}/"
              f"{r['rattler_frac_k4']:.3f} Znr(k2/3/4)={r['z_nr_k2']:.2f}/"
              f"{r['z_nr_k3']:.2f}/{r['z_nr_k4']:.2f} conv={st.converged} "
              f"it={st.iters} {st.wall_s:.1f}s", flush=True)
    return st
