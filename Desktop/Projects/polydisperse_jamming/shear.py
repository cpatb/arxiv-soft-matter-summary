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

# --- constant-PRESSURE shear (added 2026-09-12) -----------------------------------
# Why it exists, in one line: at fixed pressure the steady-state solid fraction IS
# phi_c(p) -- no ladder, no root, no fitting window, and no isostatic target. Every
# estimator this project has used so far (E2, E4, E4m) needs at least one of those, and
# at mu_p = 0.3 all of them end up reading the isotropic phi_J the shear ladders are
# centred on (PREREG_ROOT_PINNING_2026-09-09.md, PREREG_ZISO_2026-09-12.md §8).
#
# The servo is strain-clocked, not time-clocked: the volumetric strain rate is set as a
# multiple of the SHEAR rate, so `gamma` remains the only clock and the gain is
# dimensionless. Sign: Yade's velGrad diagonal is positive in EXPANSION and
# `_stress_hat` returns p positive in COMPRESSION, so an over-pressured cell
# (p > p_target) must expand, i.e. edot_v > 0. Getting this backwards runs away in one
# sample, which is the cheap way it is caught.
SERVO_GAIN = 0.5          # edot_v = GAIN * gdot * (p/p_target - 1)
SERVO_EDOT_CLAMP = 0.5    # |edot_v| <= CLAMP * gdot, so one sample cannot blow the cell
SERVO_P_TOL = 0.05        # |<p>_ss/p_target - 1| for the servo to count as HELD
# The dilatancy fit uses the FIRST part of the strain, where phi is still approaching
# phi_c; the steady window is the last STEADY_FRAC. These must not overlap.
ALPHA_FIT_FRAC = 0.4      # fit d ln|phi - phi_ss| / d gamma over the first 40 % of strain
ALPHA_MIN_DECADES = 0.3   # below this, |phi - phi_ss| never moved enough to fit a rate


def servo_edot(p_hat: float, p_hat_target: float, gdot: float,
               gain: float = SERVO_GAIN, clamp: float = SERVO_EDOT_CLAMP) -> float:
    """Volumetric strain rate the pressure servo asks for. Positive = EXPAND.

    Pure arithmetic, deliberately separate from the Yade loop so that the sign and the
    clamp can be tested without a DEM engine (Yade does not import on the Windows side
    where the test suite runs).

    An unjammed cell reads p_hat <= 0; the relative error is then floored at -1, which
    asks for the maximum CONTRACTION, which is the physically right response.
    """
    if not (p_hat_target > 0):
        raise ValueError(f"p_hat_target must be positive, got {p_hat_target!r}")
    if not math.isfinite(p_hat):
        return 0.0
    rel = max(p_hat, 0.0) / p_hat_target - 1.0
    edot = gain * gdot * rel
    lim = clamp * abs(gdot)
    return float(max(-lim, min(lim, edot)))


def dilatancy_rate(gamma: list[float] | np.ndarray, phi: list[float] | np.ndarray,
                   phi_ss: float, *, fit_frac: float = ALPHA_FIT_FRAC,
                   min_decades: float = ALPHA_MIN_DECADES) -> dict:
    """The Reynolds dilatancy rate from the APPROACH of phi to its steady value.

    The constitutive law this measures is

        edot_v = alpha * gdot * (phi - phi_c)                                    (1)

    and at constant pressure the steady state is phi = phi_c, so writing the volumetric
    strain rate in terms of phi -- with V phi = const (the grains are incompressible),
    so edot_v = d ln V / dt = -d ln phi / dt -- gives

        -(1/phi) dphi/dt = alpha gdot (phi - phi_c)
        =>  dphi/dgamma  = -alpha phi (phi - phi_c)                             (2)

    Linearised about phi_c, (2) is an exponential in the STRAIN with rate
    alpha*phi_c, so

        alpha = -[d ln|phi - phi_c| / d gamma] / phi_c                          (3)

    and the sign structure is a check rather than an assumption: a pack denser than
    phi_c must DILATE (phi falls) and a looser one must CONTRACT (phi rises), so the
    fitted rate must be NEGATIVE in |phi - phi_c| from either side, and the two sides
    must return the same alpha.

    Returns the fitted rate, alpha, and the refusal reason if any. `alpha` is NaN
    unless the fit is admissible -- a rate fitted to a signal that never moved is the
    2026-09-01 lesson (gate on the SHAPE the functional form requires, not only on
    signal-to-noise).
    """
    g = np.asarray(gamma, float)
    ph = np.asarray(phi, float)
    out = {"alpha": float("nan"), "alpha_rate": float("nan"), "alpha_decades": 0.0,
           "alpha_n": 0, "alpha_note": "", "alpha_sign_ok": False}
    if g.size != ph.size or g.size < 4 or not (phi_ss > 0):
        out["alpha_note"] = "too few samples, or phi_ss is not positive"
        return out
    win = g <= fit_frac * float(g.max())
    d = np.abs(ph - phi_ss)
    ok = win & np.isfinite(d) & (d > 0)
    if ok.sum() < 4:
        out["alpha_note"] = f"only {int(ok.sum())} usable samples in the fit window"
        return out
    y = np.log(d[ok])
    decades = float((y.max() - y.min()) / math.log(10.0))
    out["alpha_decades"] = decades
    out["alpha_n"] = int(ok.sum())
    if decades < min_decades:
        out["alpha_note"] = (f"|phi - phi_ss| moved {decades:.2f} decades, below "
                             f"{min_decades} -- there is no rate to fit")
        return out
    rate = float(np.polyfit(g[ok], y, 1)[0])
    out["alpha_rate"] = rate
    out["alpha_sign_ok"] = rate < 0
    if rate >= 0:
        out["alpha_note"] = ("|phi - phi_ss| GREW over the fit window: the state is "
                             "moving away from its own steady value")
        return out
    out["alpha"] = float(-rate / phi_ss)
    return out


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
    doublet_w: float = 0.0            # PREREG_SHEARKMER: tangent doublets of ONE primary
    young: float = 1.0
    gamma_max: float = 1.0
    shear_steps: int = 0
    # realised PSD (filled by insert_packing via a PackingState)
    n_placed: int = 0
    n_doublets: int = 0               # realised count (PREREG_SHEARKMER)
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
    # --- constant-PRESSURE fields; all NaN/empty for a constant-volume run ----------
    p_hat_target: float = float("nan")   # the servo's set point (constant-p runs only)
    servo_gain: float = float("nan")
    phi_series: list[float] = field(default_factory=list)
    edot_v_series: list[float] = field(default_factory=list)
    # 🔴 Added 2026-09-12 after the production ladder. `z_series` records z_ALL, which
    # includes rattlers, and at CV 0.5 a pack with a rattler fraction ~0.23 sits at
    # z_all ~ 2 BY CONSTRUCTION (the project's own M3 result). An admission rule with an
    # absolute floor on z_all therefore refuses the high-CV arm for being polydisperse
    # rather than for unjamming -- the genuine unjamming signature in the shakeout was
    # z_all -> 0.02 WITH p -> 1e-08, two orders lower. The non-rattler count is recorded
    # beside it so the next run's rule can use the quantity it meant.
    z_nr_series: list[float] = field(default_factory=list)
    phi_ss: float = float("nan")         # THE MEASUREMENT: phi_c(p) at steady state
    phi_ss_split: float = float("nan")   # the two halves of the steady window
    p_servo_err: float = float("nan")    # |<p>_ss/p_target - 1|
    servo_held: bool = False
    alpha: float = float("nan")          # Reynolds dilatancy coefficient
    alpha_rate: float = float("nan")     # d ln|phi - phi_ss| / d gamma
    alpha_decades: float = float("nan")
    alpha_n: int = 0
    alpha_sign_ok: bool = False
    alpha_note: str = ""
    # --- strain-clock fields (added 2026-09-12, PREREG_STRAINCLOCK_2026-09-13.md) ------
    # `mu_prep` is the contact friction during stages A-C. It equals `mu_p` on every
    # record before this date (a frictional pack prepared frictionally). The strain-clock
    # campaign prepares FRICTIONLESS (mu_prep = 0), which is the OCR -> infinity limit of
    # the unloading test cv_unloading.md §7 asks for: below phi_0(CV) such a pack carries
    # p ~ 0, so the whole pressure the shear builds is the clock. Friction is switched on
    # AFTER stage C and BEFORE stage D, and the switch is read back from the engine state
    # (`friction_switch_ok`) rather than trusted.
    mu_prep: float = float("nan")
    friction_switch_ok: bool = False
    p_hat_prep: float = float("nan")       # p/E after stage C, before any shear
    z_nr_prep: float = float("nan")
    rattler_frac_prep: float = float("nan")
    n_samples: int = 0                     # samples over the shear (N_SAMPLES if 0)
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


# --- per-contact dump (added 2026-09-22, PREREG_ACTSHAPE_2026-09-22.md) -------------
# `force_moments` keeps only a summary of f/<f>. The Wyart-Cates activation threshold is
# a force F_rep that, for every surface force this project would map onto (screened
# electrostatics, steric brushes, DLVO in the Derjaguin approximation), scales with the
# REDUCED RADIUS of the pair, R* = R1 R2/(R1 + R2). The variable that decides whether a
# contact is frictional is then f/R*, not f/<f>, and no summary on disk can form it.
# So the dump records the raw per-contact arrays, and every normalisation is chosen in
# the analysis (actshape.py), not here.
# Read-only: it iterates O.interactions between O.run calls and cannot drive the solver.
N_DUMP = 10                    # snapshots, evenly spaced over the steady window


def _contact_arrays(O):  # noqa: E741
    """(fn, ft, r1, r2) over the real sphere-sphere interactions, float64 arrays."""
    fn, ft, r1, r2 = [], [], [], []
    for i in O.interactions:
        if not i.isReal:
            continue
        b1, b2 = O.bodies[i.id1], O.bodies[i.id2]
        if b1 is None or b2 is None:
            continue
        fn.append(float(i.phys.normalForce.norm()))
        ft.append(float(i.phys.shearForce.norm()))
        r1.append(float(b1.shape.radius))
        r2.append(float(b2.shape.radius))
    return (np.asarray(fn), np.asarray(ft), np.asarray(r1), np.asarray(r2))


def dump_sample_indices(n_samples: int, steady_frac: float = STEADY_FRAC,
                        n_dump: int = N_DUMP) -> list[int]:
    """1-based sample indices (of `n_samples`) at which the dump is taken.

    Evenly spaced over the steady window and ending on the LAST sample, which is the
    state `force_moments` summarises -- so the final snapshot of every dump can be
    checked against the record's own `force` block.
    """
    first = int(math.ceil((1.0 - steady_frac) * n_samples))
    idx = np.linspace(first, n_samples, n_dump + 1)[1:]
    return sorted({int(round(v)) for v in idx})


def _prepare(st: ShearState, O, Matrix3, utils, *,  # noqa: E741
             truncate: tuple[float, float] | None) -> float:
    """Stages A-C: insert, compress to `st.phi_target`, relax. Returns the phi reached.

    ⚠ EXTRACTED VERBATIM from `shear_at_constant_volume` on 2026-09-12 so that the
    constant-PRESSURE driver runs the identical preparation. Not a rewrite: the block
    was moved, and the move is verified by re-running a committed cell and requiring the
    record to come back bit-identical (scripts/check_prepare_refactor.sh). If that check
    is not green, nothing measured through either driver is comparable with the records
    already on disk.
    """
    n, cv = st.n_requested, st.cv
    # The material is built with the PREPARATION friction. With `mu_prep` unset (every
    # record before 2026-09-12) this is `st.mu_p` and the code path is unchanged.
    mu_p = st.mu_p if math.isnan(st.mu_prep) else st.mu_prep
    seed, family, young = st.seed, st.family, st.young
    phi_target = st.phi_target
    # --- stage A: insert -----------------------------------------------------------
    # doublet_w must be carried INTO this PackingState: insert_packing reads it from here,
    # not from `st` (the first shearkmer launch built singlets for every w2 > 0 because of this)
    pk = PackingState(n_requested=n, cv=cv, mu_p=mu_p, p_hat_target=float("nan"),
                      seed=seed, family=family,
                      doublet_w=float(getattr(st, "doublet_w", 0.0)))
    insert_packing(pk, truncate=truncate, young=young)
    st.n_doublets = pk.n_doublets
    st.n_placed = pk.n_placed
    st.d10_realised = pk.d10_realised
    st.d32_realised = pk.d32_realised
    st.d43_realised = pk.d43_realised
    st.cv_realised = pk.cv_realised
    st.d_max_over_L = pk.d_max_over_L
    st.notes.extend(pk.notes)
    # clump BODIES carry no shape radius; read the member spheres (identical to before for a
    # sphere pack, which has no clumps)
    radii = np.array([b.shape.radius for b in O.bodies if b is not None and not b.isClump])
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

    return float(1.0 - utils.porosity())


def _switch_friction(st: ShearState, O, utils) -> None:  # noqa: E741
    """Set the contact friction to `st.mu_p` on every material AND every existing
    interaction, then READ IT BACK.

    `utils.setContactFriction(angle)` is compiled (`_utils.so`); the Yade docs describe
    it as writing `frictionAngle` into the FrictMat of every dynamic body and
    `tangensOfFrictionAngle` into every existing FrictPhys. That description is a
    claim about a library, so each call verifies it on the state it just produced:
    one body's material and every real interaction must carry tan(atan(mu_p)) to
    machine precision, or the record is flagged and the cell is not a measurement.
    Interactions created AFTER the switch get their friction from the (now updated)
    materials through Ip2_FrictMat_FrictMat_FrictPhys, so the read-back on the
    existing set is the only part that can silently fail.
    """
    tan_target = math.tan(math.atan(st.mu_p))
    utils.setContactFriction(math.atan(st.mu_p))
    mats_ok = all(
        math.isclose(math.tan(b.mat.frictionAngle), tan_target, rel_tol=1e-12, abs_tol=1e-15)
        for b in O.bodies if b is not None
    )
    n_real = 0
    ints_ok = True
    for i in O.interactions:
        if not i.isReal:
            continue
        n_real += 1
        if not math.isclose(i.phys.tangensOfFrictionAngle, tan_target,
                            rel_tol=1e-12, abs_tol=1e-15):
            ints_ok = False
    st.friction_switch_ok = bool(mats_ok and ints_ok)
    st.notes.append(f"friction switch: mu {st.mu_prep:g} -> {st.mu_p:g}, "
                    f"materials_ok={mats_ok} interactions_ok={ints_ok} on {n_real} real")
    if not st.friction_switch_ok:
        st.notes.append("FRICTION SWITCH FAILED READ-BACK: not a measurement")


def shear_at_constant_volume(cv: float, mu_p: float, phi_target: float, n: int,
                             seed: int, *,
                             truncate: tuple[float, float] | None = (0.2, 5.0),
                             young: float = 1.0, family: str = "lognormal",
                             gamma_max: float = 1.0, shear_steps: int = 300_000,
                             mu_prep: float | None = None,
                             n_samples: int = N_SAMPLES,
                             verbose: bool = True) -> ShearState:
    """Insert, compress to `phi_target`, then shear at constant volume.

    `mu_prep` (default None -> `mu_p`, the historical path) is the friction during
    stages A-C; when it differs from `mu_p` the friction is switched after stage C and
    the protocol string becomes "cv-shear-v3-fprep" so these records are never pooled
    with cv-shear-v2 ones. `n_samples` sets the sampling cadence of the series; the
    strain-clock campaign needs ~600 over gamma_max = 3 (0.005 per sample) where the
    default 200 gives 0.015, coarser than the half-strains it is reading.
    """
    import time

    from yade import O, Matrix3, utils

    t0 = time.time()
    O.reset()

    st = ShearState(n_requested=n, cv=cv, mu_p=mu_p, phi_target=phi_target,
                    seed=seed, family=family, young=young, gamma_max=gamma_max,
                    shear_steps=shear_steps)
    switch = False
    if mu_prep is not None and mu_prep != mu_p:
        switch = True
        st.mu_prep = float(mu_prep)
        st.protocol = "cv-shear-v3-fprep"
    st.n_samples = int(n_samples)

    phi_start = _prepare(st, O, Matrix3, utils, truncate=truncate)

    if switch:
        # The prepared state, BEFORE friction and before shear: the pressure the
        # frictionless pack carries at this phi is what the strain clock starts from.
        st.p_hat_prep, _ = _stress_hat(utils, young)
        _, st.z_nr_prep, st.rattler_frac_prep, _ = _measure_network(O)
        _switch_friction(st, O, utils)

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

    every = max(1, shear_steps // st.n_samples)
    steps_done = 0
    while steps_done < shear_steps:
        O.run(min(every, shear_steps - steps_done), wait=True)
        steps_done += min(every, shear_steps - steps_done)
        p_hat, tau_hat = _stress_hat(utils, young)
        z_all, z_nr_now, _, _ = _measure_network(O)
        st.gamma.append(gdot * O.dt * steps_done)
        st.p_hat_series.append(p_hat)
        st.tau_hat_series.append(tau_hat)
        st.z_series.append(z_all)
        st.z_nr_series.append(z_nr_now)
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


def shear_at_constant_pressure(cv: float, mu_p: float, p_hat_target: float, n: int,
                               seed: int, *, phi_start_target: float,
                               truncate: tuple[float, float] | None = (0.2, 5.0),
                               young: float = 1.0, family: str = "lognormal",
                               gamma_max: float = 3.0, shear_steps: int = 900_000,
                               servo_gain: float = SERVO_GAIN,
                               dump_path: str | None = None,
                               doublet_w: float = 0.0,
                               verbose: bool = True) -> ShearState:
    """Insert, compress to `phi_start_target`, then shear at CONSTANT PRESSURE.

    THE MEASUREMENT: `st.phi_ss` is `phi_c(p_hat_target)` directly. There is no ladder,
    no power-law root, no fitting window and no isostatic contact target -- which is
    what makes this the first `phi_c` estimator in the project that does not read the
    isotropic `phi_J` the shear ladders are centred on.

    `phi_start_target` is the state the pack is prepared in, NOT an answer. The critical
    state is by definition independent of it, so two cells prepared on opposite sides of
    `phi_c` must converge to the same `phi_ss` -- the strongest internal check available
    here, and one that needs nothing from outside the run. It also gives `alpha` from
    both the dilating and the contracting side, where the magnitudes must agree.

    ⚠ THE SERVO IS A KNOB AND ITS ANSWER MUST NOT DEPEND ON IT. `servo_gain` is recorded
    in every record; a gain sweep at fixed everything else is the validity test, and
    `phi_ss` moving with the gain means the cell is tracking the controller rather than
    the material. The answer is not quotable until that sweep is flat.
    """
    import time

    from yade import O, Matrix3, utils

    t0 = time.time()
    O.reset()

    st = ShearState(n_requested=n, cv=cv, mu_p=mu_p, phi_target=phi_start_target,
                    seed=seed, family=family, young=young, gamma_max=gamma_max,
                    shear_steps=shear_steps)
    st.doublet_w = float(doublet_w)
    st.protocol = "cp-shear-v1"
    st.p_hat_target = float(p_hat_target)
    st.servo_gain = float(servo_gain)

    phi_start = _prepare(st, O, Matrix3, utils, truncate=truncate)
    st.dt = float(O.dt)

    # --- stage D': simple shear with a pressure servo on the volumetric part --------
    # velGrad = [[e/3, gdot, 0], [0, e/3, 0], [0, 0, e/3]] -- trace = e, so the cell
    # volume follows dV/V = e dt while the shear component is untouched. `gamma` is
    # still gdot*t, so the strain clock is unchanged from the constant-volume protocol.
    gdot = gamma_max / (shear_steps * O.dt)
    st.gdot = float(gdot)
    O.cell.flipFlippable = True     # same reason as the constant-volume driver

    every = max(1, shear_steps // N_SAMPLES)
    steps_done = 0
    edot = 0.0
    n_taken = 0
    dump_at = set(dump_sample_indices(-(-shear_steps // every))) if dump_path else set()
    dump: dict[str, list] = {"gamma": [], "p_hat": [], "fn": [], "ft": [],
                             "r1": [], "r2": [], "n": []}
    while steps_done < shear_steps:
        third = edot / 3.0
        O.cell.velGrad = Matrix3(third, gdot, 0,
                                 0, third, 0,
                                 0, 0, third)
        take = min(every, shear_steps - steps_done)
        O.run(take, wait=True)
        steps_done += take
        p_hat, tau_hat = _stress_hat(utils, young)
        z_all, z_nr_now, _, _ = _measure_network(O)
        phi_now = float(1.0 - utils.porosity())
        st.gamma.append(gdot * O.dt * steps_done)
        st.p_hat_series.append(p_hat)
        st.tau_hat_series.append(tau_hat)
        st.z_series.append(z_all)
        st.z_nr_series.append(z_nr_now)
        st.phi_series.append(phi_now)
        st.edot_v_series.append(float(edot))
        n_taken += 1
        if n_taken in dump_at:
            fn_, ft_, r1_, r2_ = _contact_arrays(O)
            dump["gamma"].append(st.gamma[-1])
            dump["p_hat"].append(p_hat)
            dump["n"].append(len(fn_))
            for k, v in (("fn", fn_), ("ft", ft_), ("r1", r1_), ("r2", r2_)):
                dump[k].append(v)
        edot = servo_edot(p_hat, p_hat_target, gdot, gain=servo_gain)
    O.cell.velGrad = Matrix3.Zero
    if dump_path:
        np.savez_compressed(
            dump_path, gamma=np.asarray(dump["gamma"]), p_hat=np.asarray(dump["p_hat"]),
            n=np.asarray(dump["n"], dtype=np.int64),
            **{k: np.concatenate(dump[k]) if dump[k] else np.zeros(0)
               for k in ("fn", "ft", "r1", "r2")},
            young=young, mu_p=mu_p, cv=cv, seed=seed, p_hat_target=p_hat_target)

    st.phi = float(1.0 - utils.porosity())
    st.phi_drift = abs(st.phi - phi_start)      # here a MEASUREMENT, not a gate: the
                                                # cell is supposed to change volume

    # --- stage E': the steady state -------------------------------------------------
    g = np.asarray(st.gamma)
    p = np.asarray(st.p_hat_series)
    tau = np.asarray(st.tau_hat_series)
    ph = np.asarray(st.phi_series)
    win = g >= (1.0 - STEADY_FRAC) * gamma_max
    st.p_ss_hat = float(np.mean(p[win]))
    st.tau_ss_hat = float(np.mean(tau[win]))
    st.mu_macro = float(st.tau_ss_hat / st.p_ss_hat) if st.p_ss_hat > 0 else float("nan")
    st.phi_ss = float(np.mean(ph[win]))

    idx = np.flatnonzero(win)
    half = len(idx) // 2
    a, b = float(np.mean(p[idx[:half]])), float(np.mean(p[idx[half:]]))
    scale = max(abs(a), abs(b))
    st.steady_split = float(abs(a - b) / scale) if scale > 0 else 0.0
    pa, pb = float(np.mean(ph[idx[:half]])), float(np.mean(ph[idx[half:]]))
    st.phi_ss_split = float(abs(pa - pb) / max(abs(pa), abs(pb)))
    # BOTH channels have to be steady. p being held is the servo's job and phi being
    # settled is the physics; a cell where the servo is still hunting can hold <p> on
    # average while phi drifts, and that cell's phi_ss is not phi_c.
    st.steady = bool(st.steady_split <= STEADY_SPLIT_TOL
                     and st.phi_ss_split <= STEADY_SPLIT_TOL)
    if not st.steady:
        st.notes.append(
            f"NOT STEADY: p halves differ by {st.steady_split:.1%}, phi halves by "
            f"{st.phi_ss_split:.1%} (bar {STEADY_SPLIT_TOL:.0%})")

    st.p_servo_err = abs(st.p_ss_hat / p_hat_target - 1.0)
    st.servo_held = bool(st.p_servo_err <= SERVO_P_TOL)
    if not st.servo_held:
        st.notes.append(
            f"SERVO DID NOT HOLD: <p>_ss/p_target - 1 = "
            f"{st.p_ss_hat / p_hat_target - 1:+.1%} (bar {SERVO_P_TOL:.0%})")

    # --- the dilatancy coefficient, from the approach rather than the steady state ---
    al = dilatancy_rate(st.gamma, st.phi_series, st.phi_ss)
    st.alpha = al["alpha"]
    st.alpha_rate = al["alpha_rate"]
    st.alpha_decades = al["alpha_decades"]
    st.alpha_n = al["alpha_n"]
    st.alpha_sign_ok = al["alpha_sign_ok"]
    st.alpha_note = al["alpha_note"]
    if al["alpha_note"]:
        st.notes.append("ALPHA: " + al["alpha_note"])

    z_all, z_nr, rat, forces = _measure_network(O)
    st.z_ss, st.z_nr_ss, st.rattler_frac_ss = z_all, z_nr, rat
    st.force = force_moments(forces)

    p_dim = st.p_ss_hat * young
    st.i_ss = (float(gdot * st.d50_realised / math.sqrt(p_dim))
               if p_dim > 0 else float("inf"))
    st.quasistatic = bool(st.i_ss <= I_MAX)
    if not st.quasistatic:
        st.notes.append(f"NOT QUASI-STATIC: I = {st.i_ss:.3e} > {I_MAX:.0e}")

    st.wall_s = time.time() - t0
    if verbose:
        print(f"[shearP {family[:4]}] cv={cv:.2f} mu={mu_p:.2f} "
              f"p*={p_hat_target:.3e} K={servo_gain:.2f} phi0={phi_start:.6f} "
              f"N={n} s={seed} -> phi_ss={st.phi_ss:.6f} "
              f"p_err={st.p_servo_err:+.1%} held={st.servo_held} "
              f"alpha={st.alpha:.3g}({st.alpha_decades:.2f}dec) "
              f"mu_M={st.mu_macro:.3f} Z={st.z_ss:.2f} I={st.i_ss:.2e} "
              f"steady={st.steady}({st.steady_split:.2f}/{st.phi_ss_split:.2f}) "
              f"{st.wall_s:.1f}s", flush=True)
    return st
