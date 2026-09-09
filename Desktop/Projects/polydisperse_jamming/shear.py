"""shear.py — constant-volume (undrained) simple shear to the critical state.

WHY THIS EXISTS
---------------
Every packing number this project has measured so far comes from ISOTROPIC
compression: `phase2_phiJ.py` builds a ladder in reduced pressure and extrapolates
`phi(p_hat) -> 0` to get `phi_J`. That is the right instrument for "how dense can this
distribution get", and it is the wrong instrument for the question the mechanism now
being tested actually asks, which is:

    when this packing is SHEARED at fixed volume, does it try to expand
    (dilatant, so the pore fluid goes into suction) or to contract
    (so the pore fluid goes into overpressure)?

That sign is set by the CRITICAL-STATE packing fraction `phi_c`: the volume fraction a
granular material approaches under sustained shear, independent of its initial state.
Above `phi_c` a pack must dilate to keep flowing, below it a pack contracts. `phi_c` is
not `phi_J`, it is not measurable from an isotropic ladder, and it is the quantity that
appears in the denominator of the inertial-dilatant number of

  [M26] Mobaseri, Gordillo, Burton, Yoon, Lee, Kumar, Driscoll & Cheng,
        PRL 136, 148201 (2026)  --  ID = rho kappa U0 / (A mu_w (phi - phi_c) D0)

THE PROTOCOL, AND WHY CONSTANT VOLUME
-------------------------------------
Constant-volume simple shear IS the undrained test. Hold the total volume of a
saturated packing fixed and shear it: whatever the grain skeleton wants to do to its
volume, it cannot, so the wish is expressed as a change in the EFFECTIVE (grain-borne)
pressure instead, and the pore fluid takes the complement. A pack above `phi_c` builds
grain pressure (the fluid goes into suction -- the solid-like, jammed branch); a pack
below `phi_c` loses it (the fluid goes into overpressure -- the liquefied branch). So
the DEM observable `p_ss(phi)` at constant volume is, up to the pore-fluid bookkeeping,
the pore-pressure sign the experiment sees.

    stage A  insert a random loose packing              (packing.insert_packing)
    stage B  compress isotropically at constant strain rate to phi_target
    stage C  relax at constant volume
    stage D  simple shear at CONSTANT VOLUME to gamma_max, sampling throughout
    stage E  average over the last (1 - STEADY_FRAC) of the strain, with a
             steadiness gate that can FAIL

`velGrad` for stage D is Matrix3(0, gdot, 0,  0, 0, 0,  0, 0, 0). Its trace is exactly
zero, so the cell volume is conserved to machine precision rather than to a controller
tolerance -- `phi` during the shear is an INPUT, known exactly, not a measurement.

WHAT CAN GO WRONG, AND WHAT IS GATED
------------------------------------
  * NOT QUASI-STATIC. The inertial number I = gdot d50 sqrt(rho/p) grows without bound
    as p -> 0, which is exactly the interesting end of the ladder. `I_ss` is recorded on
    every cell and `analyse_phase8.py` excludes cells above `I_MAX`. A cell that is
    excluded is PRINTED, never dropped silently.
  * NOT STEADY. Reaching the critical state takes strain, and how much depends on how
    far the initial state is from it. The steadiness gate compares the two halves of the
    averaging window and flags disagreement; it is a gate that can fail, and it does.
  * SIGN. Yade's `getStress()` is positive in TENSION. `p_hat` here is -trace/3/E, so it
    is POSITIVE IN COMPRESSION, matching `packing.measure`. Getting this backwards
    inverts the whole result.

UNITS: dimensionless throughout, exactly as `packing.py` -- d10 = 1, rho = 1, E = young.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field

import numpy as np

from measure import RATTLER_MIN_CONTACTS, force_moments, rattler_report
from packing import (
    DT_SAFETY,
    UNBALANCED_TOL,
    PackingState,
    contact_graph,
    insert_packing,
)

# --- protocol constants, with their origin ---------------------------------------
# Compression rate for stage B, as a strain rate in units of 1/(time). Chosen so the
# whole compression from PHI_INSERT to the densest phi on the ladder takes O(1e4)
# steps: fast enough not to dominate the cell cost, slow enough that stage C can
# relax what it leaves behind. The shear (stage D) erases the preparation anyway --
# that is the defining property of a critical state, and STEADY_FRAC is what checks it.
COMPRESS_EDOT = 2.0e-4
COMPRESS_MAX_STEPS = 400_000
RELAX_STEPS = 20_000          # stage C cap; exits early on unbalancedForce < tol
# How close stage B must land on the requested phi. One refinement step at
# COMPRESS_EDOT/100 moves phi by ~1e-6, so 1e-5 is achievable and 1e-5 is 1/1500 of the
# ladder spacing and ~1/200 of the seed scatter in phi_J. A cell that misses is FLAGGED,
# and `analyse_phase8.py` fits against the REALISED phi in either case.
PHI_TARGET_TOL = 1.0e-5

STEADY_FRAC = 0.5             # average over the last half of the strain
STEADY_SPLIT_TOL = 0.15       # the two halves of that window must agree to 15 %
I_MAX = 1.0e-2                # quasi-static bar on the inertial number

# Sampling. 200 samples over the shear is enough to see a plateau and cheap.
N_SAMPLES = 200


@dataclass
class ShearState:
    """Everything measured from one constant-volume shear run. JSON-serialisable."""

    # inputs
    n_requested: int
    cv: float
    mu_p: float
    phi_target: float
    seed: int
    family: str = "lognormal"
    young: float = 1.0
    gamma_max: float = 1.0
    shear_steps: int = 0
    # realised PSD (filled by insert_packing via a PackingState)
    n_placed: int = 0
    d10_realised: float = float("nan")
    d32_realised: float = float("nan")
    d43_realised: float = float("nan")
    cv_realised: float = float("nan")
    d50_realised: float = float("nan")
    d_max_over_L: float = float("nan")
    # state reached
    phi: float = float("nan")          # realised, should equal phi_target to ~1e-12
    phi_drift: float = float("nan")    # |phi_end - phi_start|, the constant-volume check
    gdot: float = float("nan")
    dt: float = float("nan")
    # steady-state readings (averages over the last STEADY_FRAC of the strain)
    p_ss_hat: float = float("nan")     # p/E, POSITIVE IN COMPRESSION
    tau_ss_hat: float = float("nan")   # -sigma_xy / E
    mu_macro: float = float("nan")     # tau/p, the macroscopic friction
    z_ss: float = float("nan")
    z_nr_ss: float = float("nan")
    rattler_frac_ss: float = float("nan")
    i_ss: float = float("nan")         # inertial number at the steady state
    # gates
    steady: bool = False
    steady_split: float = float("nan")
    quasistatic: bool = False
    compressed: bool = False           # stage B reached phi_target
    # series, for plotting and for re-reading a gate without re-running the cell
    gamma: list[float] = field(default_factory=list)
    p_hat_series: list[float] = field(default_factory=list)
    tau_hat_series: list[float] = field(default_factory=list)
    z_series: list[float] = field(default_factory=list)
    force: dict = field(default_factory=dict)
    wall_s: float = float("nan")
    protocol: str = "cv-shear-v2"
    notes: list[str] = field(default_factory=list)

    def to_json(self) -> str:
        return json.dumps(asdict(self), sort_keys=True)


def _stress_hat(utils, young: float) -> tuple[float, float]:
    """(p_hat, tau_hat) from the volume-averaged Cauchy stress.

    Yade's getStress() is positive in TENSION, so p = -trace/3 is positive in
    compression.

    The shear component does NOT flip with the pressure. With the imposed velocity
    gradient L_xy = +gdot (v_x = gdot*y), a material that RESISTS the imposed motion
    carries Cauchy sigma_xy > 0 -- the same sign a Newtonian fluid would give,
    sigma_xy = mu * dv_x/dy. Measured on the 2026-09-08 smoke cell: sigma_xy =
    +4.38e-3 against p = +1.21e-2, i.e. mu_macro = +0.362, which is the physically
    admissible sign and a plausible magnitude for mu_p = 0.3. The first version of
    this function flipped it and returned mu_macro = -0.362; a negative macroscopic
    friction is impossible, which is how the error was caught in one run.
    """
    s = np.array(utils.getStress())
    p_hat = float(-np.trace(s) / 3.0 / young)
    tau_hat = float(s[0, 1] / young)
    return p_hat, tau_hat


def _measure_network(O):  # noqa: E741
    counts, neigh, forces = contact_graph(O)
    n_all = len(counts)
    nreal = sum(counts.values()) // 2
    z_all = 2.0 * nreal / n_all if n_all else float("nan")
    rep = rattler_report(counts, neigh)
    return (z_all, rep[f"z_nr_k{RATTLER_MIN_CONTACTS}"],
            rep[f"rattler_frac_k{RATTLER_MIN_CONTACTS}"], forces)


def shear_at_constant_volume(cv: float, mu_p: float, phi_target: float, n: int,
                             seed: int, *,
                             truncate: tuple[float, float] | None = (0.2, 5.0),
                             young: float = 1.0, family: str = "lognormal",
                             gamma_max: float = 1.0, shear_steps: int = 300_000,
                             verbose: bool = True) -> ShearState:
    """Insert, compress to `phi_target`, then shear at constant volume."""
    import time

    from yade import O, Matrix3, utils

    t0 = time.time()
    O.reset()

    st = ShearState(n_requested=n, cv=cv, mu_p=mu_p, phi_target=phi_target,
                    seed=seed, family=family, young=young, gamma_max=gamma_max,
                    shear_steps=shear_steps)

    # --- stage A: insert -----------------------------------------------------------
    pk = PackingState(n_requested=n, cv=cv, mu_p=mu_p, p_hat_target=float("nan"),
                      seed=seed, family=family)
    insert_packing(pk, truncate=truncate, young=young)
    st.n_placed = pk.n_placed
    st.d10_realised = pk.d10_realised
    st.d32_realised = pk.d32_realised
    st.d43_realised = pk.d43_realised
    st.cv_realised = pk.cv_realised
    st.d_max_over_L = pk.d_max_over_L
    st.notes.extend(pk.notes)
    radii = np.array([b.shape.radius for b in O.bodies if b is not None])
    st.d50_realised = float(np.median(2.0 * radii))

    O.engines = [
        utils.ForceResetter(),
        utils.InsertionSortCollider([utils.Bo1_Sphere_Aabb()],
                                    allowBiggerThanPeriod=False),
        utils.InteractionLoop(
            [utils.Ig2_Sphere_Sphere_ScGeom()],
            [utils.Ip2_FrictMat_FrictMat_FrictPhys()],
            [utils.Law2_ScGeom_FrictPhys_CundallStrack()],
        ),
        utils.NewtonIntegrator(damping=0.4),
    ]
    O.dt = DT_SAFETY * utils.PWaveTimeStep()
    st.dt = float(O.dt)

    # --- stage B: isotropic compression to phi_target, to a DECLARED tolerance -------
    # velGrad = -edot * I contracts the cell homogeneously, so phi is a pure function
    # of the number of steps taken and can be driven onto a target exactly.
    #
    # 🔴 The first version of this ran in fixed 500-step chunks and stopped on the
    # first chunk that crossed the target. One chunk is a linear strain of
    # edot*dt*500 = 2.5e-2, i.e. a jump of ~0.04 in phi -- nearly THREE TIMES the
    # 0.015 spacing of the ladder this function exists to measure. Two rungs 0.015
    # apart then landed on the SAME state and returned bit-identical p_hat, Z and
    # rattler fractions (measured 2026-09-08, cv=0: phi_target 0.5595 and 0.5745 both
    # gave p_hat = 1.6786e-04). The stop condition was quantized by its own chunk.
    #
    # The fix computes the number of steps to the target from the CURRENT volume every
    # iteration and FLOORS it, so the state is always approached from below and never
    # overshoots, then refines with a 100x slower rate. `phi` is measured afterwards
    # and reported; `phi_target` is never used downstream as if it were achieved.
    # The rate magnitude is stepped down 20x twice. A coarse pass gets most of the way
    # in few steps; each finer pass resolves what the previous one could not, because
    # `floor` can never place a step smaller than one step's worth of strain. It is
    # SIGNED -- an overshoot is corrected by expanding rather than accepted -- which the
    # first version was not: Yade advances the cell by exp(velGrad*dt) while the step
    # count here assumes (1 + velGrad*dt), so every cell landed ~8.4e-5 ABOVE target,
    # eight times the tolerance, and the `compressed` flag would have fired on all of
    # them (a flag that always fires reports nothing).
    budget = COMPRESS_MAX_STEPS
    for mag in (COMPRESS_EDOT, COMPRESS_EDOT / 20.0, COMPRESS_EDOT / 400.0):
        for _ in range(400):
            phi_now = 1.0 - utils.porosity()
            if abs(phi_now - phi_target) <= PHI_TARGET_TOL or budget <= 0:
                break
            # linear scale factor required: <1 compresses, >1 expands
            want = (phi_now / phi_target) ** (1.0 / 3.0)
            edot = -mag if want < 1.0 else +mag
            per_step_scale = 1.0 + edot * O.dt
            nsteps = int(math.floor(math.log(want) / math.log(per_step_scale)))
            if nsteps < 1:
                break
            nsteps = min(nsteps, budget)
            O.cell.velGrad = Matrix3(edot, 0, 0, 0, edot, 0, 0, 0, edot)
            O.run(nsteps, wait=True)
            budget -= nsteps
    O.cell.velGrad = Matrix3.Zero
    phi_reached = 1.0 - utils.porosity()
    st.compressed = bool(abs(phi_reached - phi_target) <= PHI_TARGET_TOL)
    if not st.compressed:
        st.notes.append(
            f"STAGE B MISSED TARGET: phi={phi_reached:.6f} vs {phi_target:.6f} "
            f"(tol {PHI_TARGET_TOL:.0e}, {COMPRESS_MAX_STEPS - budget} steps used)"
        )

    # --- stage C: relax at constant volume -----------------------------------------
    chunk = 500
    relaxed = 0
    while relaxed < RELAX_STEPS:
        O.run(chunk, wait=True)
        relaxed += chunk
        if float(utils.unbalancedForce()) < UNBALANCED_TOL:
            break

    phi_start = float(1.0 - utils.porosity())

    # --- stage D: constant-volume simple shear -------------------------------------
    # trace(velGrad) = 0 exactly, so the volume is conserved to machine precision.
    gdot = gamma_max / (shear_steps * O.dt)
    st.gdot = float(gdot)

    # 🔴 Let Yade flip the periodic cell. A cell sheared to gamma ~ 1 becomes skewed
    # enough that the shortest distance between periodic images falls below twice the
    # largest grain, and `InsertionSortCollider` then aborts the run with
    #   "Body larger than half of the cell size encountered."
    # Measured 2026-09-08: EVERY rung of cv = 0.2 seed 1 died this way -- five cells,
    # deterministically, because that seed's largest grain is near the truncation limit
    # -- while cv = 0 (where d_max = d10) ran fine. Left unfixed it would have taken out
    # the cv = 0.5 arm, i.e. the arm the whole measurement exists to compare against.
    #
    # Flipping replaces the cell by an equivalent, more compact one. It is a relabelling
    # of the periodic lattice: the volume, the stress and the contact network are
    # unchanged, so it is not a physical knob. It DOES perturb the trajectory at
    # floating-point level, so it must be on for every cell in a comparison or none --
    # which is why the protocol string is bumped to v2 and cv-shear-v1 records are not
    # pooled with these.
    O.cell.flipFlippable = True
    O.cell.velGrad = Matrix3(0, gdot, 0,
                             0, 0, 0,
                             0, 0, 0)

    every = max(1, shear_steps // N_SAMPLES)
    steps_done = 0
    while steps_done < shear_steps:
        O.run(min(every, shear_steps - steps_done), wait=True)
        steps_done += min(every, shear_steps - steps_done)
        p_hat, tau_hat = _stress_hat(utils, young)
        z_all, _, _, _ = _measure_network(O)
        st.gamma.append(gdot * O.dt * steps_done)
        st.p_hat_series.append(p_hat)
        st.tau_hat_series.append(tau_hat)
        st.z_series.append(z_all)
    O.cell.velGrad = Matrix3.Zero

    st.phi = float(1.0 - utils.porosity())
    st.phi_drift = abs(st.phi - phi_start)

    # --- stage E: steady-state average and its gate ---------------------------------
    g = np.asarray(st.gamma)
    p = np.asarray(st.p_hat_series)
    tau = np.asarray(st.tau_hat_series)
    win = g >= (1.0 - STEADY_FRAC) * gamma_max
    st.p_ss_hat = float(np.mean(p[win]))
    st.tau_ss_hat = float(np.mean(tau[win]))
    st.mu_macro = float(st.tau_ss_hat / st.p_ss_hat) if st.p_ss_hat > 0 else float("nan")

    idx = np.flatnonzero(win)
    half = len(idx) // 2
    a, b = float(np.mean(p[idx[:half]])), float(np.mean(p[idx[half:]]))
    scale = max(abs(a), abs(b))
    st.steady_split = float(abs(a - b) / scale) if scale > 0 else 0.0
    st.steady = bool(st.steady_split <= STEADY_SPLIT_TOL)
    if not st.steady:
        st.notes.append(
            f"NOT STEADY: the two halves of the averaging window differ by "
            f"{st.steady_split:.1%} (bar {STEADY_SPLIT_TOL:.0%})"
        )

    z_all, z_nr, rat, forces = _measure_network(O)
    st.z_ss, st.z_nr_ss, st.rattler_frac_ss = z_all, z_nr, rat
    st.force = force_moments(forces)

    # I = gdot d50 sqrt(rho/p); rho = 1 in these units and p is the DIMENSIONAL
    # pressure p_hat * young.
    p_dim = st.p_ss_hat * young
    st.i_ss = (float(gdot * st.d50_realised / math.sqrt(p_dim))
               if p_dim > 0 else float("inf"))
    st.quasistatic = bool(st.i_ss <= I_MAX)
    if not st.quasistatic:
        st.notes.append(f"NOT QUASI-STATIC: I = {st.i_ss:.3e} > {I_MAX:.0e}")

    st.wall_s = time.time() - t0
    if verbose:
        print(f"[shear {family[:4]}] cv={cv:.2f} mu={mu_p:.2f} "
              f"phi={st.phi:.6f}(tgt {phi_target:.6f}) "
              f"N={n} s={seed} -> p_hat={st.p_ss_hat:.4e} tau_hat={st.tau_ss_hat:.4e} "
              f"mu_M={st.mu_macro:.3f} Z={st.z_ss:.2f} rat={st.rattler_frac_ss:.3f} "
              f"I={st.i_ss:.2e} steady={st.steady}({st.steady_split:.2f}) "
              f"drift={st.phi_drift:.1e} {st.wall_s:.1f}s", flush=True)
    return st
