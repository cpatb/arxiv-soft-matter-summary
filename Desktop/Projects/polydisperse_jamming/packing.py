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
import os
from dataclasses import asdict, dataclass, field

import numpy as np

from measure import (
    RATTLER_MIN_CONTACTS,
    RATTLER_THRESHOLDS,
    fabric_tensor,
    force_moments,
    iterated_rattlers,
    rattler_report,
)
from psd import make_psd
from shape import (
    clump_centres,
    corrugation,
    random_axes,
    span_factor,
    sub_sphere_diameter,
)

# --- protocol constants, with their origin ---------------------------------------
UNBALANCED_TOL = 0.02     # Yade quasi-static criterion. Swept in test_packing.py;
                          # O'Hern-style studies use 1e-2..1e-3 on this ratio.
REL_STRESS_TOL = 2e-3     # PeriTriaxController's own stress-convergence bar
PHI_INSERT = 0.30
DT_SAFETY = 0.5
MAX_ITER = 200_000
# ⚠ OVERRIDABLE, added 2026-09-09. 82 of the isotropic cells on disk are marked
# non-converged and EVERY ONE of them stopped at exactly MAX_ITER -- they are slow,
# not wrong. The drop rate is not flat in the axis this project measures: for the
# lognormal FRICTIONLESS arm it runs 0.92 (CV 0) to 0.33 (CV 0.5), against 1.00 -> 0.75
# at mu_p = 0.3, so the two arms of any frictionless-vs-frictional difference are
# censored at different, CV-dependent rates. `PJ_MAX_ITER` raises the cap so those
# cells can be re-run; the value actually used is recorded in every record as
# `max_iter`, so a run at a non-default cap can never be pooled with one at the
# default without it being visible.
def max_iter() -> int:
    """The iteration cap, read at CALL time so a launcher can set it per run."""
    v = os.getenv("PJ_MAX_ITER")
    return MAX_ITER if not v else int(v)
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
    # --- the SHAPE axis, added 2026-09-18. alpha = 1.0 on every record written before
    # that date and on every sphere run since, so a reader can filter on it safely.
    #
    # alpha = 1.0 builds a SINGLE SPHERE, not a clump of two coincident spheres, and
    # that is a correction to PREREG_SHAPE_2026-09-18.md forced by a measurement:
    # scripts/probe_clump_s0.sh Q2 found that two clumps of coincident members produce
    # FOUR interactions per grain pair (25 distinct pairs behind 100 interactions), so
    # the naive alpha = 1 clump is not the sphere this project has run -- it is a body
    # with every contact quadrupled. The union of two coincident spheres IS one sphere,
    # so the single sphere is the correct construction at u = 0, and the alpha ladder
    # has a SEAM there rather than a limit: dup_contact_frac jumps from 0 to 0.75 the
    # moment two members exist. That seam is registered, not smoothed over.
    alpha: float = 1.0
    # --- the MULTIPLET axis, added 2026-09-22 (PREREG_KMER_2026-09-22.md). ONE primary
    # sphere (CV must be 0, alpha must be 1); a grain is a rigid cluster of k TANGENT
    # primaries: k = 2 a dumbbell, k = 3 an equilateral triangle. `doublet_w` and
    # `triplet_w` are SOLID-VOLUME fractions; the rest is singlets. 0/0 is the historical
    # path and draws nothing from any RNG, so every earlier record reproduces bit-for-bit.
    doublet_w: float = 0.0
    triplet_w: float = 0.0
    n_doublets: int = 0
    n_triplets: int = 0
    eps_corrugation: float = 0.0      # shape.corrugation; the INTERLOCKING half
    span_over_deq: float = 1.0        # shape.span_factor; what S2's clearance uses
    density_scale: float = 1.0        # 1/(2*(ds/d_eq)^3) for clumps, 1.0 for spheres
    # Total body mass read back from YADE divided by the analytic union volume. With the
    # density correction applied this must be 1.0 for both spheres and clumps, and it is
    # the only channel on which a wrong `ds`, a wrong `g(u)` or a wrong density is
    # visible without another simulation -- gate S0c reads it. Read from the scene, not
    # recomputed from the same formula that set it (which would be a tautology): the
    # density goes IN as 1/(2*(ds/d_eq)^3) and the mass comes OUT of Yade's own clump
    # assembly, so agreement tests that assembly.
    mass_over_vol_solid: float = float("nan")
    dt: float = float("nan")          # DT_SAFETY*PWaveTimeStep(); moves with density
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
    # phi is ALWAYS the analytic union volume over the cell volume -- one estimator
    # across the whole (CV, alpha) grid, so there is no seam at alpha = 1 to explain.
    # `utils.porosity()` is recorded beside it and NEVER used as phi: probe_clump.sh
    # measured it DOUBLE-COUNTING the member overlap by 2*(ds/d_eq)^3, which runs 2.000
    # at alpha = 1 to 1.000 at alpha = 2 -- an uncorrected campaign would have reported
    # phi falling by a factor of two along the aspect ladder out of pure bookkeeping.
    phi_union: float = float("nan")      # analytic; == phi, kept so the pair is visible
    phi_porosity: float = float("nan")   # 1 - utils.porosity(); the DEFECTIVE channel
    vol_solid: float = float("nan")      # sum (pi/6) d_eq^3 over grains
    # For SPHERES the two agree to the last ulp and not exactly: probe_clump_s0.sh Q4
    # measured |analytic - utils.porosity()| = 5.55e-17 at phi = 0.171 (1 ulp), because
    # Yade sums body volumes in its own order. That is why gate S0a reads bit identity
    # on `phi_porosity` and not on `phi`.
    p_hat: float = float("nan")             # realised p/E, positive in compression
    stress_dev_hat: float = float("nan")    # von Mises deviator / E
    z_all: float = float("nan")
    z_nr: float = float("nan")              # over NON-RATTLERS only
    rattler_frac: float = float("nan")
    n_contacts: int = 0
    # THE CONTACT-MULTIPLICITY CONFOUND, found by scripts/probe_clump_s0.sh and absent
    # from PREREG_SHAPE_2026-09-18.md, which asked only whether clump members contact
    # EACH OTHER (they do not). The mirror question is whether an EXTERNAL neighbour
    # contacts BOTH members of one clump, and it does: two clumps in contact can raise
    # up to FOUR member-member interactions. `n_contacts` is therefore the count of
    # DISTINCT GRAIN PAIRS -- which is the physical contact number and is what z_all,
    # the rattler set and the backbone are built from -- and the raw interaction count
    # is recorded beside it so the artefact has a number rather than a caveat.
    #
    # It cannot be removed: two overlapping members in contact with the same neighbour
    # are two springs in parallel, so the MECHANICS carry it even though the counting
    # does not. It is a third member of this axis's confound list, beside elongation
    # and corrugation, and S4's label says so.
    n_contacts_raw: int = 0           # real interactions, before grain-pair dedup
    n_contacts_internal: int = 0      # member-member within one clump; measured 0
    dup_contact_frac: float = 0.0     # (raw - internal - distinct)/(raw - internal)
    unbalanced: float = float("nan")
    converged: bool = False
    iters: int = 0
    eta_roll: float | None = None  # None = FrictMat (every record before 2026-09-13)
    max_iter: int = MAX_ITER      # the cap this cell actually ran under
    tol_met: bool = False         # unbalanced < UNBALANCED_TOL, ignoring the cap
    wall_s: float = float("nan")
    rattler: dict = field(default_factory=dict)   # every threshold, see measure
    force: dict = field(default_factory=dict)     # contact-force moments
    fabric: dict = field(default_factory=dict)    # contact-orientation anisotropy.
    #                                               EMPTY on every record written before
    #                                               2026-09-18, so read it with .get()
    fabric_nr: dict = field(default_factory=dict)  # the same, over the NON-RATTLER
    #                                                BACKBONE only. Separate because the
    #                                                rattler fraction runs 0.040 -> 0.234
    #                                                across CV, so an all-contact fabric
    #                                                mixes a change in the rattler
    #                                                POPULATION with a change in the
    #                                                backbone's STRUCTURE. EMPTY on every
    #                                                record written before 2026-09-18.
    protocol: str = "iso-stress-nondim-v2"
    notes: list[str] = field(default_factory=list)

    def to_json(self) -> str:
        return json.dumps(asdict(self), sort_keys=True)


def contact_graph(O, grain_of=None, stats=None):  # noqa: E741
    """(counts, neighbours, forces) over the REAL interactions of a Yade scene.

    `grain_of` maps a BODY id to the id of the GRAIN it belongs to. For a sphere pack it
    is the identity and this function does exactly what it did before the shape axis
    existed; for a clump pack the bodies are the two member spheres and the grain is the
    clump, so the graph has to be built at the grain level or `z_all`'s denominator
    counts clumps AND members and is wrong by more than a factor of two.

    The graph is built over DISTINCT GRAIN PAIRS. For spheres that changes nothing --
    two spheres raise exactly one interaction, so the deduplicated counts equal the old
    per-interaction counts body for body, which is what gate S0a checks rather than
    assumes. For clumps it is the whole point: probe_clump_s0.sh Q2 measured 100 real
    interactions behind 25 distinct clump pairs at alpha = 1, i.e. every contact
    appearing four times, and a `z_all` built from the raw count would have been that
    factor wrong.

    `stats`, if given, is a dict filled IN PLACE with the raw and internal interaction
    counts and the duplicate fraction. It is an out-parameter rather than a fourth
    return value because six call sites unpack this as a 3-tuple and widening it would
    break every one of them for a diagnostic none of them wants.
    """
    if grain_of is None:
        grain_of = {b.id: b.id for b in O.bodies if b is not None}
    grains = sorted(set(grain_of.values()))
    counts: dict[int, int] = {g: 0 for g in grains}
    neigh: dict[int, set[int]] = {g: set() for g in grains}
    forces: list[float] = []
    n_raw = 0
    n_internal = 0
    seen: set[tuple[int, int]] = set()
    for i in O.interactions:
        if not i.isReal:
            continue
        g1 = grain_of.get(i.id1)
        g2 = grain_of.get(i.id2)
        if g1 is None or g2 is None:
            continue
        if g1 == g2:
            n_internal += 1
            continue
        n_raw += 1
        try:
            fn = float(i.phys.normalForce.norm())
        except Exception:
            fn = float("nan")
        forces.append(fn)
        key = (g1, g2) if g1 < g2 else (g2, g1)
        if key in seen:
            continue
        seen.add(key)
        counts[g1] += 1
        counts[g2] += 1
        neigh[g1].add(g2)
        neigh[g2].add(g1)
    if stats is not None:
        stats["n_raw"] = n_raw
        stats["n_internal"] = n_internal
        stats["n_distinct"] = len(seen)
        stats["dup_frac"] = (n_raw - len(seen)) / n_raw if n_raw else 0.0
    return counts, neigh, np.asarray(forces, dtype=float)


def contact_normals(O, with_ids: bool = False):  # noqa: E741
    """Unit contact normals over the REAL interactions, as an (n, 3) array.

    A SEPARATE function rather than a fourth return from `contact_graph`, because six
    call sites unpack that as a 3-tuple and widening it would break every one of them
    for an observable none of them wants. The extra pass over `O.interactions` is free
    beside the relaxation that produced the scene.

    Read AFTER relaxation and used only by `fabric_tensor`, so it observes the solver
    and cannot drive it -- the failure mode lessons.md records for MEASURE_NSAMP, where
    an output probe became the binding constraint on dt.

    With `with_ids` it also returns the (n, 2) array of the two body ids behind each
    normal, which is what lets `measure` read the fabric over the NON-RATTLER BACKBONE
    as well as over every contact. That split is the whole point of recording this: the
    rattler fraction runs 0.040 -> 0.234 across CV, so an all-contact fabric mixes a
    change in the rattler POPULATION with a change in the backbone's STRUCTURE, and only
    the second one answers the question the rattler result left open. The two lists are
    appended in the SAME iteration behind the SAME guard, never filtered independently
    (lessons.md 2026-09-05: pairing two columns by filtering each one separately
    misaligns them the moment one has a gap).
    """
    out = []
    ids = []
    for i in O.interactions:
        if not i.isReal:
            continue
        try:
            n = i.geom.normal
            row = [float(n[0]), float(n[1]), float(n[2])]
        except Exception:
            continue
        out.append(row)
        ids.append([int(i.id1), int(i.id2)])
    arr = np.asarray(out, dtype=float).reshape(-1, 3)
    if with_ids:
        return arr, np.asarray(ids, dtype=int).reshape(-1, 2)
    return arr


def backbone_mask(ids, backbone) -> np.ndarray:
    """Boolean mask over contacts whose BOTH ends are in `backbone`.

    Pure and id-based, so it is testable without a Yade scene.
    """
    ids = np.asarray(ids, dtype=int).reshape(-1, 2)
    if ids.shape[0] == 0:
        return np.zeros(0, dtype=bool)
    bb = set(int(b) for b in backbone)
    return np.array([(int(a) in bb) and (int(b) in bb) for a, b in ids], dtype=bool)


def measure(O, utils, st: PackingState, young: float,  # noqa: E741
            grain_of=None) -> PackingState:
    """Read every observable off a relaxed Yade scene into `st`.

    `grain_of` is None for a sphere pack (the identity map) and maps member body id ->
    clump id for a clump pack. Everything downstream -- z, the rattler set, the
    backbone, the fabric -- is then built over GRAINS, which is the only level at which
    those observables mean what their names say.
    """
    cstats: dict[str, float] = {}
    counts, neigh, forces = contact_graph(O, grain_of=grain_of, stats=cstats)
    st.n_contacts_raw = int(cstats["n_raw"])
    st.n_contacts_internal = int(cstats["n_internal"])
    st.dup_contact_frac = float(cstats["dup_frac"])
    n_all = len(counts)
    nreal = sum(counts.values()) // 2
    st.n_contacts = nreal
    st.z_all = 2.0 * nreal / n_all if n_all else float("nan")

    # The rattler threshold is a definition that changes the answer, so all three are
    # recorded and none is privileged in the data. See measure.RATTLER_THRESHOLDS.
    # Per-GRAIN volume, for the rattler VOLUME fraction (added 2026-09-19). Only for a
    # sphere pack: for a clump pack the grain volume is the analytic union
    # (`shape.clump_volume`), and `utils.porosity()` double-counts the member overlap by
    # 2(ds/d_eq)^3 -- a trap this project has already been caught by once. Rather than
    # risk repeating it, clump packs get NO volume field at all and say so by its absence.
    _vols = None
    if grain_of is None:
        _vols = {}
        for b in O.bodies:
            if b is None:
                continue
            r = getattr(getattr(b, "shape", None), "radius", None)
            if r is not None:
                _vols[b.id] = 4.0 / 3.0 * math.pi * float(r) ** 3
    st.rattler = rattler_report(counts, neigh, volumes=_vols)
    st.rattler_frac = st.rattler[f"rattler_frac_k{RATTLER_MIN_CONTACTS}"]
    st.z_nr = st.rattler[f"z_nr_k{RATTLER_MIN_CONTACTS}"]
    st.force = force_moments(forces)
    # Contact-orientation anisotropy, over ALL real contacts and over the NON-RATTLER
    # BACKBONE. Recorded from 2026-09-18; absent on every earlier record.
    #
    # Two things about how to read it, both registered in PREREG_FABRIC_2026-09-18.md.
    # (1) `aniso` is the von Mises norm of the deviator, and under isotropic compression
    #     the ENSEMBLE-mean deviator is zero by symmetry at every CV -- so the deviator
    #     itself carries no CV signal that symmetry permits. What does carry one is its
    #     MAGNITUDE relative to the finite-sample floor: for n INDEPENDENT unit normals
    #     E[aniso^2] = 1/n exactly, so T = n * aniso^2 has E[T] = 1 and is free of the
    #     contact-count confound. T > 1 is orientational CORRELATION in the network.
    # (2) n_contacts falls 3682 -> 2050 across CV here, so a raw `aniso` vs CV trend is
    #     manufactured by the floor alone (0.0157 -> 0.0210). Never read `aniso` against
    #     CV without dividing the floor out.
    # (3) The backbone fabric is recorded at EVERY rattler threshold, not only the
    #     primary one, because `rattler_report` already establishes that convention for
    #     this project: a conclusion that holds at one k is a statement about the
    #     definition. It matters more here than it looks. At k = 2 a rattler is a grain
    #     with <= 1 contact, so removing 404 of 1500 grains removed TWO of 2762 contacts
    #     in the gate cell -- the k = 2 backbone fabric is all but identical to the
    #     all-contact one, and the conflation this split exists to avoid is negligible
    #     THERE and real at k = 4.
    normals, cids = contact_normals(O, with_ids=True)
    if grain_of is not None:
        # The backbone mask is keyed on GRAIN ids, so the contact ids have to be lifted
        # from members to clumps before the mask is applied. Mapping in place keeps the
        # normals and the ids paired row for row (lessons.md 2026-09-05: filtering two
        # columns independently misaligns them the moment one has a gap).
        cids = np.array([[grain_of.get(int(a), int(a)), grain_of.get(int(b), int(b))]
                         for a, b in cids], dtype=int).reshape(-1, 2)
    st.fabric = fabric_tensor(normals)
    fnr: dict[str, float] = {}
    for k in RATTLER_THRESHOLDS:
        bb_k = set(counts) - iterated_rattlers(counts, neigh, min_contacts=k)
        f_k = fabric_tensor(normals[backbone_mask(cids, bb_k)])
        fnr[f"aniso_k{k}"] = f_k["aniso"]
        fnr[f"n_k{k}"] = f_k["n"]
    # The primary threshold is aliased to the bare keys so a reader who does not know
    # about the threshold ladder still gets the PRIMARY definition and not an arbitrary
    # one -- the same choice `rattler_frac` and `z_nr` already make.
    fnr["aniso"] = fnr[f"aniso_k{RATTLER_MIN_CONTACTS}"]
    fnr["n"] = fnr[f"n_k{RATTLER_MIN_CONTACTS}"]
    st.fabric_nr = fnr

    s = np.array(utils.getStress())
    st.p_hat = float(-np.trace(s) / 3.0 / young)
    dev = s - np.eye(3) * (np.trace(s) / 3.0)
    st.stress_dev_hat = float(np.sqrt(1.5 * np.sum(dev * dev)) / young)

    # phi from the ANALYTIC UNION VOLUME, always. `utils.porosity()` is recorded beside
    # it and never used: it double-counts the member overlap by 2*(ds/d_eq)^3, which is
    # alpha-dependent (2.000 -> 1.000 along this ladder). The ratio is a PREDICTION with
    # a value, so gate S0c can fail on it rather than merely observing it.
    st.phi_union = float(st.vol_solid / O.cell.volume)
    st.phi_porosity = float(1.0 - utils.porosity())
    st.phi = st.phi_union
    st.unbalanced = float(utils.unbalancedForce())
    st.iters = int(O.iter)
    return st


# --- the rolling-resistance path -----------------------------------------------------
# ALPHA_KR is Yade's dimensionless rolling-stiffness coefficient (k_r = alphaKr * k_n *
# r1 * r2). 2.0 is the library default and is kept rather than tuned: this arm exists to
# move phi_J, not to fit anything, and a tuned stiffness would be a second knob.
ALPHA_KR = 2.0


def interaction_loop(utils, eta_roll):
    """The InteractionLoop for this contact model.

    `eta_roll is None` returns EXACTLY the engine list every record on disk was taken
    with. Anything else returns the 6-degree-of-freedom cohesive-frictional stack with
    cohesion off and the moment law forced on.

    🔴 THE TRAP, verified against this Yade build by probe_roll_api.py:
    `Law2_ScGeom6D_CohFrictPhys_CohesionMoment.always_use_moment_law` DEFAULTS TO FALSE,
    and with `isCohesive=False` nothing is ever bonded -- so without this flag the
    rolling-resistance knob would be silently INERT and the whole arm would measure
    nothing while returning plausible numbers.
    """
    if eta_roll is None:
        return utils.InteractionLoop(
            [utils.Ig2_Sphere_Sphere_ScGeom()],
            [utils.Ip2_FrictMat_FrictMat_FrictPhys()],
            [utils.Law2_ScGeom_FrictPhys_CundallStrack()],
        )
    return utils.InteractionLoop(
        [utils.Ig2_Sphere_Sphere_ScGeom6D()],
        [utils.Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(
            setCohesionNow=False, setCohesionOnNewContacts=False)],
        [utils.Law2_ScGeom6D_CohFrictPhys_CohesionMoment(
            always_use_moment_law=True, useIncrementalForm=True)],
    )


def insert_packing(st: "PackingState", *, truncate: tuple[float, float] | None,
                   young: float) -> tuple[float, int]:
    """Insert one random loose packing into a fresh periodic Yade cell.

    EXTRACTED 2026-09-08 from `build_and_compress` so the constant-volume shear
    protocol (`shear.py`) can reuse the *identical* preparation. Nothing about the
    random-number order changed in the extraction: `rng` is seeded, radii are drawn
    first, then one uniform position per body, exactly as before. That is the whole
    reproducibility surface, and `test_shear.py::test_insertion_rng_order_unchanged`
    pins it.

    Fills the realised-PSD fields of `st` in place and returns (cell edge L,
    material id).
    """
    from yade import O, utils

    rng = np.random.default_rng(st.seed)
    dist = make_psd(st.family, 1.0, st.cv)
    radii = dist.radii_for_yade(st.n_requested, rng, truncate=truncate)
    diam = 2.0 * radii

    m = dist.realised_moments(diam)
    st.n_placed = int(m["n"])
    st.d10_realised = m["d10"]
    st.d32_realised = m["d32"]
    st.d43_realised = m["d43"]
    st.cv_realised = m["cv"]
    st.d_span = m["span"]

    vol_solid = float(np.sum(4.0 / 3.0 * math.pi * radii**3))
    st.vol_solid = vol_solid
    d_max = float(diam.max())

    # --- the SHAPE axis ---------------------------------------------------------------
    # `d_eq` is the volume-equivalent diameter and is drawn from the SAME psd.make_psd at
    # every alpha, so `vol_solid` above IS the union volume of the clumps and every PSD
    # descriptor is invariant along this axis by construction (gate S1 checks it anyway).
    # What is NOT invariant is the LONGEST EXTENT: an elongated grain of the same volume
    # is longer than the sphere by span_factor(alpha), and the periodic collider's
    # clearance is set by the longest extent. So `d_span_max` -- not `d_max` -- is what
    # the clearance check below must use.
    alpha = float(st.alpha)
    st.span_over_deq = span_factor(alpha)
    st.eps_corrugation = corrugation(1.0, alpha)
    ds_over_deq = sub_sphere_diameter(1.0, alpha)
    is_clump = alpha > 1.0
    # The mass correction. A clump's mass is the SUM of its member sphere masses, and at
    # alpha < 2 the members overlap, so that sum exceeds the union mass by exactly
    # 2*(ds/d_eq)^3 (probe_clump.sh measured 0.62056151 against 0.52359878 at alpha=1.5,
    # ratio 1.185185). The ratio is independent of grain size, so ONE density per alpha
    # fixes it for the whole polydisperse pack.
    st.density_scale = 1.0 / (2.0 * ds_over_deq ** 3) if is_clump else 1.0
    d_span_max = st.span_over_deq * d_max

    # --- MULTIPLETS of one primary (see PackingState.doublet_w). Tangent members do not
    # overlap, so the union volume is k primary volumes EXACTLY, the member density needs
    # no correction, and utils.porosity() is not double-counting (the overlap is zero).
    kmer = {}
    if st.doublet_w > 0.0 or st.triplet_w > 0.0:
        if st.cv != 0.0 or alpha != 1.0:
            raise ValueError("multiplets are built from ONE primary: need cv = 0, alpha = 1")
        w2, w3 = float(st.doublet_w), float(st.triplet_w)
        w1 = 1.0 - w2 - w3
        if w1 < -1e-12:
            raise ValueError(f"doublet_w + triplet_w = {w2 + w3} > 1")
        # number fractions: a k-mer carries k primary volumes
        z = w1 + w2 / 2.0 + w3 / 3.0
        n = len(radii)
        n2 = int(round(n * (w2 / 2.0) / z))
        n3 = int(round(n * (w3 / 3.0) / z))
        # a SEPARATE stream, drawn only on this path, so the sphere path is untouched
        pick = np.random.default_rng(st.seed + 2_000_003).permutation(n)
        kmer = {int(j): 2 for j in pick[:n2]}
        kmer.update({int(j): 3 for j in pick[n2:n2 + n3]})
        st.n_doublets, st.n_triplets = n2, n3
        v1 = 4.0 / 3.0 * math.pi * radii ** 3
        vol_solid = float(np.sum(v1) + sum((k - 1) * v1[j] for j, k in kmer.items()))
        st.vol_solid = vol_solid
        # longest extent: a dumbbell spans 2 primary diameters; a tangent triangle's
        # enclosing length is also 2 (two members side by side along one edge)
        d_span_max = 2.0 * d_max

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
    if L_final_min < COLLIDER_CLEARANCE * d_span_max:
        v_mean = vol_solid / len(radii)
        n_min = math.ceil(PHI_DENSE_BOUND * (COLLIDER_CLEARANCE * d_span_max) ** 3
                          / v_mean)
        raise ValueError(
            f"cell too small for the periodic collider at cv={st.cv}, N={st.n_requested}, "
            f"alpha={alpha}, truncate={truncate}: largest grain d_max={d_max:.3f} spans "
            f"{d_span_max:.3f} at this alpha and needs a final cell edge >= "
            f"{COLLIDER_CLEARANCE * d_span_max:.3f} but phi={PHI_DENSE_BOUND} "
            f"gives {L_final_min:.3f}. Use N >= {n_min}, or tighten `truncate`."
        )

    L = (vol_solid / PHI_INSERT) ** (1.0 / 3.0)
    # The finite-size diagnostic is read on the SPAN, not the equivalent diameter: what
    # interacts with its own periodic image is the longest extent. At alpha = 1 the two
    # coincide exactly, so every sphere record's flag is unchanged.
    st.d_max_over_L = d_span_max / L_final_min      # worst case, at the densest state
    # A grain wider than ~1/4 of the cell interacts with its own periodic image through
    # very few intermediaries; flag it rather than trusting the number.
    if st.d_max_over_L > 0.25:
        st.notes.append(
            f"FINITE SIZE RISK: d_max/L_dense = {st.d_max_over_L:.3f} > 0.25"
        )

    if st.eta_roll is None:
        mat = O.materials.append(
            utils.FrictMat(young=young, poisson=POISSON, density=st.density_scale,
                           frictionAngle=math.atan(st.mu_p), label="grain")
        )
    else:
        # ROLLING-RESISTANCE PATH (transfer_law.py's lever; see module docstring).
        # `isCohesive=False` makes this a COHESIONLESS pack that merely resists
        # rolling -- not a bonded aggregate. `momentRotationLaw=True` is what turns
        # the moment branch on at all; it defaults to False.
        # `etaRoll` is the plastic moment limit in units of (normal force x radius):
        # 0 means the moment is always at its limit, i.e. NO rolling resistance, which
        # is this path's own control against the FrictMat default.
        mat = O.materials.append(
            utils.CohFrictMat(young=young, poisson=POISSON, density=1.0,
                              frictionAngle=math.atan(st.mu_p),
                              isCohesive=False,
                              normalCohesion=0.0, shearCohesion=0.0,
                              momentRotationLaw=True,
                              etaRoll=float(st.eta_roll),
                              etaTwist=float(st.eta_roll),
                              alphaKr=ALPHA_KR, alphaKtw=ALPHA_KR,
                              label="grain")
        )
    O.periodic = True
    O.cell.hSize = np.diag([L, L, L])
    # Orientations come from a SEPARATE, independently seeded stream, so drawing them
    # cannot perturb the (radii, then one position per grain) order that every record on
    # disk was written with. That is not a stylistic choice: gate S0a is a bit-identity
    # reproduction of three committed sphere cells, and it can only fail on the RNG
    # order if the RNG order is capable of moving. Drawn UNCONDITIONALLY, including at
    # alpha = 1 where they are unused, so the alpha = 1 path executes the same lines.
    axes = random_axes(len(radii), st.seed + 1_000_003)
    for j, r in enumerate(radii):
        pos = rng.uniform(0, L, 3)
        k = kmer.get(j, 1)
        if k > 1:
            # centres about the centroid: dumbbell along the random axis; triangle in the
            # plane spanned by the axis and a perpendicular, side = one primary diameter
            a = np.asarray(axes[j], dtype=float)
            a = a / np.linalg.norm(a)
            if k == 2:
                off = [a * float(r), -a * float(r)]
            else:
                t = np.cross(a, [1.0, 0.0, 0.0])
                if np.linalg.norm(t) < 1e-6:
                    t = np.cross(a, [0.0, 1.0, 0.0])
                t = t / np.linalg.norm(t)
                side = 2.0 * float(r)
                verts = [np.zeros(3), a * side,
                         a * side / 2.0 + t * side * math.sqrt(3.0) / 2.0]
                cen = sum(verts) / 3.0
                off = [v - cen for v in verts]
            ids = [O.bodies.append(utils.sphere(tuple(pos + o), float(r), material=mat))
                   for o in off]
            O.bodies.clump(ids)
            continue
        if not is_clump:
            O.bodies.append(utils.sphere(tuple(pos), float(r), material=mat))
            continue
        # A clump of two sub-spheres of diameter `ds`, centres `u*ds` apart, whose UNION
        # volume is exactly (pi/6)*d_eq^3. `clump_centres` returns the offsets about the
        # centroid, so `pos` is the centroid and the grain occupies the same mean
        # position the sphere would have.
        d_eq = 2.0 * float(r)
        ds = sub_sphere_diameter(d_eq, alpha)
        off = clump_centres(d_eq, alpha, axes[j])
        ids = [O.bodies.append(utils.sphere(tuple(pos + off[k]), ds / 2.0,
                                            material=mat)) for k in (0, 1)]
        O.bodies.clump(ids)
    return L, mat


def kmer_contact_dump(O):  # noqa: E741
    """Per-contact (fn, k1, k2) over real interactions between DIFFERENT grains, where k is
    the number of primaries in the grain each member belongs to (1 for a free sphere).
    For a one-primary multiplet pack every contact is primary-primary, so R* = a/2 is the
    same for every contact and fn alone is the activation variable (PREREG_KMER)."""
    ksize, gid = {}, {}
    for b in O.bodies:
        if b is None or b.isClump:
            continue
        if b.isClumpMember:
            gid[b.id] = b.clumpId
            ksize[b.id] = len(O.bodies[b.clumpId].shape.members)
        else:
            gid[b.id], ksize[b.id] = b.id, 1
    fn, k1, k2 = [], [], []
    for i in O.interactions:
        if not i.isReal or i.id1 not in gid or i.id2 not in gid:
            continue
        if gid[i.id1] == gid[i.id2]:
            continue
        fn.append(float(i.phys.normalForce.norm()))
        k1.append(ksize[i.id1])
        k2.append(ksize[i.id2])
    return (np.asarray(fn), np.asarray(k1, dtype=np.int64), np.asarray(k2, dtype=np.int64))


def grain_map(O):  # noqa: E741
    """body id -> GRAIN id, for a scene that may contain clumps.

    A clump member's grain is its clump; a standalone sphere is its own grain; the
    clump body itself carries no shape and is NOT a graph node (including it would put
    a zero-contact entry into the rattler count for every grain in the pack).

    Returns None for a pure sphere scene, which is what `contact_graph` and `measure`
    take to mean "the identity map" -- so the sphere path runs the code it always ran.
    """
    gm = {}
    any_clump = False
    for b in O.bodies:
        if b is None:
            continue
        if b.isClump:
            any_clump = True
            continue
        gm[b.id] = b.clumpId if b.isClumpMember else b.id
    return gm if any_clump else None


def build_and_compress(cv: float, mu_p: float, p_hat: float, n: int, seed: int,
                       truncate: tuple[float, float] | None = (0.2, 5.0),
                       young: float = 1.0, verbose: bool = True,
                       family: str = "lognormal",
                       eta_roll: float | None = None,
                       alpha: float = 1.0, doublet_w: float = 0.0,
                       triplet_w: float = 0.0) -> PackingState:
    """Compress a packing isotropically to reduced pressure `p_hat` = p/E.

    Dimensionless throughout: d10 = 1, rho = 1, E = `young` (default 1).
    """
    import time

    from yade import O, utils

    t0 = time.time()
    O.reset()
    st = PackingState(n_requested=n, cv=cv, mu_p=mu_p, p_hat_target=p_hat,
                      seed=seed, family=family, eta_roll=eta_roll, alpha=float(alpha),
                      doublet_w=float(doublet_w), triplet_w=float(triplet_w))
    insert_packing(st, truncate=truncate, young=young)
    grain_of = grain_map(O)
    # Read the mass back off the SCENE, before anything moves. Clump members carry
    # `state.mass` = 0 once clumped (the clump body carries it), so summing over every
    # body with a finite mass counts each grain exactly once either way.
    st.mass_over_vol_solid = float(
        sum(b.state.mass for b in O.bodies if b is not None and not b.isClumpMember)
        / st.vol_solid)

    # Yade sign convention for the triax goal: compressive stress is NEGATIVE.
    goal = -abs(p_hat) * young
    O.engines = [
        utils.ForceResetter(),
        utils.InsertionSortCollider([utils.Bo1_Sphere_Aabb()],
                                    allowBiggerThanPeriod=False),
        interaction_loop(utils, st.eta_roll),
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
    # PWaveTimeStep reads the MEMBER radius and the MEMBER density, and the mass
    # correction changes the density -- so the timestep moves along the alpha ladder.
    # Measured (probe_clump_s0.sh Q3), as a ratio to the equal-volume sphere's:
    # 0.7071 / 0.7449 / 0.7715 / 0.7879 / 0.7937 at alpha = 1 / 1.25 / 1.5 / 1.75 / 2.
    # It moves in the SAFE direction (smaller dt) at every rung, and it is recorded on
    # every cell rather than inferred. It is also why gate S0 as originally written --
    # a bit-identity reproduction of the sphere through a CLUMP built at alpha = 1 --
    # cannot be met by any correct implementation: the density that makes the mass right
    # makes the timestep 0.7071x, which is a different trajectory.
    O.dt = DT_SAFETY * utils.PWaveTimeStep()
    st.dt = float(O.dt)
    cap = max_iter()
    st.max_iter = cap
    O.run(cap, wait=True)

    st = measure(O, utils, st, young, grain_of=grain_of)
    st.wall_s = time.time() - t0
    st.converged = bool(st.unbalanced < UNBALANCED_TOL and O.iter < cap)
    # Recorded separately: `unbalanced` alone met the physical bar. 18 of the 82
    # non-converged cells on disk are in this state -- rejected ONLY on the iteration
    # count -- and 12 of those 18 sit in the starved frictionless arm.
    st.tol_met = bool(st.unbalanced < UNBALANCED_TOL)
    if not st.converged:
        st.notes.append(
            f"NOT CONVERGED: unbalanced={st.unbalanced:.3e} tol={UNBALANCED_TOL} "
            f"iters={st.iters}/{cap} tol_met={st.tol_met}"
        )
    rel_p = abs(st.p_hat - p_hat) / p_hat if p_hat > 0 else float("nan")
    if np.isfinite(rel_p) and rel_p > 0.10:
        st.notes.append(f"pressure off target by {rel_p:.1%} "
                        f"(realised {st.p_hat:.3e} vs {p_hat:.3e})")
    if verbose:
        r = st.rattler
        print(f"[{family[:4]}] a={st.alpha:.2f} eps={st.eps_corrugation:.3f} "
              f"dup={st.dup_contact_frac:.3f} nraw={st.n_contacts_raw} "
              f"phi_por={st.phi_porosity:.6f} "
              f"cv={cv:.2f} mu={mu_p:.2f} phat_t={p_hat:.1e} N={n} s={seed} -> "
              f"phi={st.phi:.4f} phat={st.p_hat:.3e} Zall={st.z_all:.3f} "
              f"rat(k2/3/4)={r['rattler_frac_k2']:.3f}/{r['rattler_frac_k3']:.3f}/"
              f"{r['rattler_frac_k4']:.3f} Znr(k2/3/4)={r['z_nr_k2']:.2f}/"
              f"{r['z_nr_k3']:.2f}/{r['z_nr_k4']:.2f} conv={st.converged} "
              f"it={st.iters} {st.wall_s:.1f}s", flush=True)
    return st
