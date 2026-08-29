"""
Monte Carlo trajectory sweep over launch positions, headings, throw types,
and release-parameter scatter.

For each (launch_pos, heading_wind, throw_type, throw_bearing) cell, draws
n_throws release realizations, integrates each through the 6-DOF model with
CFD-mean + gust wind, and aggregates landing statistics including the
"delivered-to-target" success rate that feeds the P5 throwability metric.

Reproducibility
---------------
Every cell derives its RNG from a single `--seed` via `np.random.SeedSequence`
spawning, so results are bit-reproducible and independent of how the work is
distributed across processes. The seed is recorded in the output file.

Parallelism
-----------
Cells are independent, so the sweep is embarrassingly parallel. Workers build
their own `CFDField` from the npz path (interpolator objects are expensive to
pickle) and cache it per process.

Common random numbers
---------------------
The RNG spawn key deliberately EXCLUDES the heading index and the throw-type
index is retained only because the two throw types consume different release
distributions. Heading is a non-random parameter of the cell, so cells that
differ only in heading share release draws and gust phases exactly. The
between-heading contrast is the comparison this study exists to make, and it is
the one that needs the variance reduction most: ~98% of per-throw variance is
the shared release draw. See `test_coupling.py::test_headings_share_random_draws`.

REFERENCE FRAMES (read before using any saved array)
----------------------------------------------------
Two frames appear in the output and they are NOT interchangeable:

  world      : the CFD site frame. x east, y north, z up, metres. Absolute.
  throw      : origin at the launch point, +downrange along `bearing_deg`,
               +crossrange 90 deg counterclockwise from it (i.e. to the
               thrower's left). Relative.

  `landing_mean`  is WORLD-frame absolute (x, y, z).
  `landing_bias`  is THROW-frame, relative to the intended target
                  (downrange - target_range_m, crossrange).
  `target_xy`     is WORLD-frame, and is now saved so that no consumer has to
                  reconstruct it. It is bearing-rotated:
                      target = launch_xy + R * (cos b, sin b).
  `landing_cov`, `ellipse_axes` are translation-invariant and therefore frame-
                  free, except that the covariance axes are world-aligned.

Reconstructing the target as a FIXED point at bearing 0 -- which the old
`P4_Coupling/README.md` invited by describing `landing_mean` as "(dx, dy) from
target" -- produces spurious 25-50 m "biases" for the five non-zero bearings.
The saved `landing_bias` and `target_xy` exist to make that mistake impossible.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "P1_AnalyticalModel"))
import disc_6dof as d6                              # noqa: E402
from coefficients import DEFAULT_COEFFS, ULTRASTAR  # noqa: E402
from trajectories import sample_release             # noqa: E402

from cfd_field import StubField, CFDField           # noqa: E402
import synthetic_turbulence as st                 # noqa: E402
from synthetic_turbulence import (CorrelatedGusts, wind_with_gusts,  # noqa: E402
                                  wind_mean_only)

# Canonical scoring target, per the manuscript: a 25 m throw is "delivered" if
# it lands within 3 m of the intended target point.
TARGET_RANGE_M = 25.0
TARGET_RADIUS_M = 3.0

# Monte Carlo throws may be integrated at a looser tolerance than the
# validation runs; the convergence study in P1 shows the landing point moves by
# ~0.02 mm, six orders of magnitude below TARGET_RADIUS_M.
MC_RTOL = 1e-4
MC_ATOL = 1e-6

# A throw is only a landing if the terminal ground event actually fired.
# `solve_ivp` returns the state at `t_max` when it does not, and that state is
# finite, so a finiteness test alone would score a disc still in the air as
# having "landed" at its current altitude.
GROUND_TOL_M = 1e-3

# Integration horizon. A throw still airborne at T_MAX is not scored a landing.
T_MAX = 8.0

# Saved into every artifact: which reference frame each field lives in, and
# which population it was computed over. Written because a reader who assumed
# `landing_mean` was target-relative measured 25-50 m of nonexistent bias.
SCHEMA_FRAMES = (
    "landing_mean:   world frame (x,y,z), absolute site coordinates, m; "
    "population = delivered (landed, not blocked)\n"
    "landing_bias:   throw frame (downrange-minus-target, crossrange), m, "
    "relative to target_xy; crossrange positive to thrower's left; "
    "population = delivered\n"
    "target_xy:      world frame, = launch_xy + target_range_m * "
    "(cos bearing, sin bearing)\n"
    "landing_cov:    world-aligned axes, translation-invariant; "
    "population = delivered\n"
    "ellipse_axes:   90% containment semi-axes from landing_cov, m; "
    "population = delivered\n"
    "success_rate:   fraction of VALID throws landing within target_radius_m "
    "of target_xy; blocked throws are counted in the denominator as failures\n"
    "n_valid:        throws that reached the ground inside t_max\n"
    "n_delivered:    n_valid minus those landing inside a building footprint\n"
    "n_airborne:     throws still in flight at t_max; excluded everywhere"
)


def sweep_one_launch_point(
    launch, heading_deg, throw_type, n_throws, field,
    bearing_deg: float = 0.0,
    seed=0,
    t_max: float = 8.0,
    target_range_m: float = TARGET_RANGE_M,
    target_radius_m: float = TARGET_RADIUS_M,
    use_gusts: bool = True,
    gust_model: str = "sweeping",
    sigma_ratios=None,
):
    """Aggregated landing statistics for one (launch, heading, throw, bearing).

    Parameters
    ----------
    seed : int or np.random.SeedSequence
        Anything `np.random.default_rng` accepts.

    Returns
    -------
    dict with keys: landings, mean, cov, ellipse_axes, success_rate,
    alpha_landing_mean, target_xy, n_valid.
    """
    rng = np.random.default_rng(seed)
    launch = np.asarray(launch, dtype=float)
    if sigma_ratios is None:
        sigma_ratios = st.SIGMA_RATIOS

    landings = np.full((n_throws, 3), np.nan)
    alpha_lands = np.full(n_throws, np.nan)
    scales = np.full(n_throws, np.nan)
    n_clipped = 0

    for i in range(n_throws):
        # A fresh gust realization per throw, drawn up front so that the wind
        # is a deterministic function of (r, t) during the integration. The
        # integral SCALE comes from the field at the launch point; which scale
        # (a length for 'sweeping', a time for 'temporal') is the modelling
        # choice documented in synthetic_turbulence.
        gusts = CorrelatedGusts(field, rng, t_max=t_max, ref_pos=launch,
                                model=gust_model, sigma_ratios=sigma_ratios)
        scales[i] = gusts.scale
        n_clipped += int(gusts.clipped)
        # NB: the gust realization is drawn either way so the two arms of the
        # paired experiment consume the same random numbers (see wind_mean_only).
        wind = (wind_with_gusts if use_gusts else wind_mean_only)(field, gusts)

        r0, v0, q0, w0 = sample_release(
            throw_type=throw_type,
            heading_deg=bearing_deg,
            launch_pos=launch,
            rng=rng,
        )
        res = d6.simulate(r0, v0, q0, w0,
                          coeffs=DEFAULT_COEFFS, geom=ULTRASTAR,
                          wind=wind, t_max=t_max,
                          rtol=MC_RTOL, atol=MC_ATOL,
                          compute_diag=False)
        landings[i] = res.r[-1]
        alpha_lands[i] = res.diag[-1].get("alpha_deg", np.nan)

    # A throw that never reached the ground inside t_max is not a landing.
    # Finiteness is not sufficient: check the disc is actually ON the ground.
    valid = np.isfinite(landings[:, 0]) & (np.abs(landings[:, 2]) <= GROUND_TOL_M)
    n_airborne = int((np.isfinite(landings[:, 0])
                      & ~(np.abs(landings[:, 2]) <= GROUND_TOL_M)).sum())

    # A throw whose landing point lies inside a building was never delivered.
    # The 6-DOF model has no collision detection, so without this such throws
    # would be scored on wind interpolated from the nearest open cell -- i.e.
    # credited with flying through the building. Campus Towers sits 9 m off the
    # field's eastern edge, so this is not a corner case for easterly bearings.
    blocked = np.zeros(n_throws, dtype=bool)
    if hasattr(field, "is_solid"):
        for i in range(n_throws):
            if valid[i]:
                blocked[i] = field.is_solid(landings[i, 0], landings[i, 1])
    n_blocked = int(blocked.sum())
    n_valid = int(valid.sum())

    # Two populations, and they must not be confused:
    #   `valid`     -- reached the ground. Used as the DELIVERY denominator,
    #                  because a throw into a building is a failed delivery,
    #                  not a missing observation.
    #   `delivered` -- reached the ground in open air. Used for all DISPERSION
    #                  statistics, because the landing point of a blocked throw
    #                  is fictitious: the 6-DOF model has no collision
    #                  detection, so the disc was flown through the building
    #                  and "landed" somewhere on the far side. Including those
    #                  points inflates the covariance with coordinates that
    #                  describe no physical trajectory.
    # Before this was separated, `mean`/`cov`/`ellipse_axes` used `valid` while
    # `success_rate` used `valid & ~blocked` -- two different populations
    # reported side by side as if they described the same throws.
    delivered = valid & ~blocked
    n_delivered = int(delivered.sum())
    good = landings[delivered]

    # Intended target: target_range_m downrange along the throw bearing,
    # in the WORLD frame. `success_rate` and `landing_bias` are both measured
    # against THIS point; see the module docstring on reference frames.
    bearing = np.radians(bearing_deg)
    e_down = np.array([np.cos(bearing), np.sin(bearing)])
    e_cross = np.array([-np.sin(bearing), np.cos(bearing)])
    target_xy = launch[:2] + target_range_m * e_down

    if n_delivered >= 2:
        mean = good.mean(axis=0)
        cov = np.cov(good[:, :2].T)
        eigvals = np.linalg.eigvalsh(cov)
        # 90% containment for an isotropic-ish 2D Gaussian: r^2/sigma^2 =
        # -2 ln(0.1). eigvals can go slightly negative from round-off.
        chi2_90 = -2.0 * np.log(0.1)
        ellipse_axes = np.sqrt(np.maximum(eigvals, 0.0) * chi2_90)
        # Throw-frame bias of the landing centroid relative to the target.
        d_world = mean[:2] - target_xy
        landing_bias = np.array([d_world @ e_down, d_world @ e_cross])
    else:
        mean = np.full(3, np.nan)
        cov = np.full((2, 2), np.nan)
        ellipse_axes = np.full(2, np.nan)
        landing_bias = np.full(2, np.nan)

    if n_valid >= 2:
        d_target = np.linalg.norm(landings[valid][:, :2] - target_xy, axis=1)
        hit = (d_target <= target_radius_m) & (~blocked[valid])
        success_rate = float(np.mean(hit))
    else:
        success_rate = np.nan

    return dict(
        landings=landings,
        mean=mean, cov=cov,
        ellipse_axes=ellipse_axes,
        landing_bias=landing_bias,
        success_rate=success_rate,
        alpha_landing_mean=float(np.nanmean(alpha_lands)) if n_valid else np.nan,
        target_xy=target_xy,
        n_valid=n_valid,
        n_delivered=n_delivered,
        n_airborne=n_airborne,
        n_blocked=n_blocked,
        blocked_fraction=n_blocked / max(n_throws, 1),
        # What the gust model ACTUALLY did in this cell. None of this was
        # recorded before, so silent clipping at T_L_MAX was undetectable and
        # the realised correlation scale could not be reconstructed at all.
        gust_scale=float(np.mean(scales)) if n_throws else np.nan,
        gust_clipped_fraction=n_clipped / max(n_throws, 1),
    )


# ---------------------------------------------------------------------------
# Parallel driver
# ---------------------------------------------------------------------------

_FIELD_CACHE: dict = {}


def _get_field(field_path, heading):
    """Build (and per-process cache) the wind field for one heading."""
    key = (field_path, heading)
    if key not in _FIELD_CACHE:
        if field_path is None:
            _FIELD_CACHE[key] = StubField(np.zeros(3))
        elif str(field_path).startswith("stub:"):
            # "stub:vx,vy,vz" -- a uniform analytic field for smoke testing.
            vec = np.array([float(v) for v in str(field_path)[5:].split(",")])
            _FIELD_CACHE[key] = StubField(vec)
        else:
            _FIELD_CACHE[key] = CFDField(field_path, heading=heading)
    return _FIELD_CACHE[key]


def crn_spawn_key(idx):
    """Map a cell index to its common-random-number spawn key.

    `idx` is (ix, iy, i_h, it, ib). The HEADING index is dropped: heading is a
    non-random property of the cell, so the three headings at one (launch,
    throw, bearing) share release draws and gust phases exactly, and the
    between-heading difference is estimated on paired throws. Keeping i_h in
    the key gave the study's central comparison independent samples and hence
    no variance reduction at all.

    `it` (throw type) is retained because the throw types draw from different
    release distributions, so pairing across them is not meaningful anyway.
    """
    ix, iy, _i_h, it, ib = idx
    return (ix, iy, it, ib)


def _run_cell(task):
    """Worker entry point: run one cell and return its summary scalars."""
    (idx, launch, heading, throw, bearing, n_throws, field_path,
     seed_entropy, t_max, use_gusts, gust_model, sigma_ratios) = task
    field = _get_field(field_path, heading)
    seed = np.random.SeedSequence(entropy=seed_entropy,
                                  spawn_key=crn_spawn_key(idx))
    stats = sweep_one_launch_point(
        launch=launch, heading_deg=heading, throw_type=throw,
        bearing_deg=bearing, n_throws=n_throws, field=field,
        seed=seed, t_max=t_max, use_gusts=use_gusts,
        gust_model=gust_model, sigma_ratios=sigma_ratios,
    )
    return idx, (stats["mean"], stats["cov"], stats["ellipse_axes"],
                 stats["landing_bias"], stats["target_xy"],
                 stats["success_rate"], stats["alpha_landing_mean"],
                 stats["n_valid"], stats["n_delivered"], stats["n_airborne"],
                 stats["blocked_fraction"], stats["gust_scale"],
                 stats["gust_clipped_fraction"])


def _git_sha() -> str:
    """Short SHA of the tree that produced this artifact, '+dirty' if unclean."""
    import subprocess
    repo = Path(__file__).resolve().parent
    try:
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             cwd=repo, capture_output=True, text=True,
                             timeout=10).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"],
                               cwd=repo, capture_output=True, text=True,
                               timeout=10).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    if not sha:
        return "unknown"
    return sha + ("+dirty" if dirty else "")


def _file_hash(path) -> str:
    """SHA-256 of an input artifact, so downstream stages can bind to it.

    `field_source` was previously recorded as a bare path string, which is not
    a binding: re-running the CFD to the same path and re-running only P5 mixes
    fresh w_f/t_f with stale success rates and nothing errors or warns.
    """
    import hashlib
    if path is None or str(path).startswith("stub:"):
        return "none"
    p = Path(path)
    if not p.is_file():
        return "missing"
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _provenance(args) -> dict:
    """Everything needed to reconstruct what this run actually did.

    None of this was recorded before, so the realised gust parameters -- the
    single most consequential modelling choice in the pipeline -- could not be
    reconstructed from the artifact at all.
    """
    import platform
    import scipy
    import synthetic_turbulence as st
    return dict(
        prov_git_sha=_git_sha(),
        prov_python=platform.python_version(),
        prov_numpy=np.__version__,
        prov_scipy=scipy.__version__,
        prov_field_sha256=_file_hash(args.field),
        prov_mc_rtol=MC_RTOL,
        prov_mc_atol=MC_ATOL,
        prov_t_max=T_MAX,
        prov_ground_tol_m=GROUND_TOL_M,
        prov_n_modes=st.CorrelatedGusts.N_MODES,
        prov_gust_dt=st.GUST_DT,
        prov_T_L_min=st.T_L_MIN,
        prov_T_L_max=st.T_L_MAX,
    )


def parse_grid(args_grid):
    xa, xb, dx, ya, yb, dy = args_grid
    xs = np.arange(xa, xb + dx / 2, dx)
    ys = np.arange(ya, yb + dy / 2, dy)
    return xs, ys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--field", default=None,
                    help="P3 flow_fields.npz, or 'stub:vx,vy,vz' for a "
                         "uniform analytic field, or omitted for still air")
    ap.add_argument("--launch-grid", nargs=6, type=float,
                    metavar=("XA", "XB", "DX", "YA", "YB", "DY"),
                    default=[-25, 25, 10, -25, 25, 10])
    ap.add_argument("--launch-z", type=float, default=1.5)
    ap.add_argument("--headings", nargs="+", type=int, required=True)
    ap.add_argument("--throws", nargs="+",
                    default=["right_backhand", "right_forehand"])
    ap.add_argument("--throw-bearings", nargs="+", type=float, default=[0])
    ap.add_argument("--n-throws", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20260725)
    ap.add_argument("--jobs", type=int, default=0,
                    help="worker processes (0 = cpu_count-2)")
    ap.add_argument("--no-gusts", action="store_true",
                    help="use the RANS mean flow only, suppressing the "
                         "synthetic gust field. The counterfactual arm of a "
                         "paired experiment: release draws are unchanged, so "
                         "the scatter difference isolates turbulence.")
    ap.add_argument("--gust-model", choices=["sweeping", "temporal"],
                    default="sweeping",
                    help="'sweeping' (default) decorrelates the gust over path "
                         "arclength with integral LENGTH L_int -- the frozen-"
                         "field limit appropriate to a body that supplies its "
                         "own advection. 'temporal' decorrelates over "
                         "wall-clock time with T_L = k/eps at the launch "
                         "point, the previous behaviour; at this site that is "
                         "10-18 s against a 2.5 s flight, i.e. a randomly "
                         "drawn constant wind. Retained for comparison.")
    ap.add_argument("--isotropic-gusts", action="store_true",
                    help="force sigma_v = sigma_w = sigma_u instead of the "
                         "surface-layer ratios (0.75, 0.5). Use to isolate the "
                         "correlation-model change from the anisotropy change.")
    ap.add_argument("--out", default="landing_stats.npz")
    args = ap.parse_args()
    sigma_ratios = (st.ISOTROPIC_RATIOS if args.isotropic_gusts
                    else st.SIGMA_RATIOS)

    xs, ys = parse_grid(args.launch_grid)
    headings, throws, bearings = args.headings, args.throws, args.throw_bearings
    n_lx, n_ly = len(xs), len(ys)
    n_h, n_t, n_b = len(headings), len(throws), len(bearings)

    shape = (n_lx, n_ly, n_h, n_t, n_b)
    mean_arr = np.full(shape + (3,), np.nan)
    cov_arr = np.full(shape + (2, 2), np.nan)
    ax_arr = np.full(shape + (2,), np.nan)
    bias_arr = np.full(shape + (2,), np.nan)
    target_arr = np.full(shape + (2,), np.nan)
    succ_arr = np.full(shape, np.nan)
    alpha_arr = np.full(shape, np.nan)
    nvalid_arr = np.zeros(shape, dtype=int)
    ndeliv_arr = np.zeros(shape, dtype=int)
    nair_arr = np.zeros(shape, dtype=int)
    blocked_arr = np.full(shape, np.nan)
    gscale_arr = np.full(shape, np.nan)
    gclip_arr = np.full(shape, np.nan)

    tasks = []
    for ix, x_l in enumerate(xs):
        for iy, y_l in enumerate(ys):
            for i_h, h in enumerate(headings):
                for it, throw in enumerate(throws):
                    for ib, bear in enumerate(bearings):
                        tasks.append((
                            (ix, iy, i_h, it, ib),
                            (float(x_l), float(y_l), args.launch_z),
                            int(h), throw, float(bear), args.n_throws,
                            args.field, args.seed, T_MAX,
                            not args.no_gusts, args.gust_model,
                            sigma_ratios,
                        ))

    n_jobs = args.jobs or max(1, (os.cpu_count() or 2) - 2)
    total_throws = len(tasks) * args.n_throws
    print(f"{len(tasks)} cells x {args.n_throws} throws "
          f"= {total_throws:,} trajectories on {n_jobs} workers")

    t0 = time.perf_counter()
    if n_jobs == 1:
        results = map(_run_cell, tasks)
        results = list(results)
    else:
        import multiprocessing as mp
        with mp.Pool(n_jobs) as pool:
            results = []
            for n, out in enumerate(pool.imap_unordered(_run_cell, tasks,
                                                        chunksize=1), 1):
                results.append(out)
                if n % max(1, len(tasks) // 20) == 0 or n == len(tasks):
                    el = time.perf_counter() - t0
                    rate = n / el
                    eta = (len(tasks) - n) / rate if rate > 0 else float("nan")
                    print(f"  {n}/{len(tasks)} cells  "
                          f"[{el/60:.1f} min elapsed, ETA {eta/60:.1f} min]",
                          flush=True)

    for (idx, (mean, cov, axes, bias, tgt, succ, alpha,
               nvalid, ndeliv, nair, blk, gsc, gcl)) in results:
        mean_arr[idx] = mean
        cov_arr[idx] = cov
        ax_arr[idx] = axes
        bias_arr[idx] = bias
        target_arr[idx] = tgt
        succ_arr[idx] = succ
        alpha_arr[idx] = alpha
        nvalid_arr[idx] = nvalid
        ndeliv_arr[idx] = ndeliv
        nair_arr[idx] = nair
        blocked_arr[idx] = blk
        gscale_arr[idx] = gsc
        gclip_arr[idx] = gcl

    elapsed = time.perf_counter() - t0
    print(f"swept {total_throws:,} trajectories in {elapsed/60:.1f} min "
          f"({total_throws/max(elapsed,1e-9):.0f} throws/s)")

    np.savez_compressed(
        args.out,
        launch_x=xs, launch_y=ys, launch_z=args.launch_z,
        headings=np.array(headings),
        throws=np.array(throws),
        bearings=np.array(bearings),
        landing_mean=mean_arr,
        landing_cov=cov_arr,
        ellipse_axes=ax_arr,
        landing_bias=bias_arr,
        target_xy=target_arr,
        success_rate=succ_arr,
        alpha_landing=alpha_arr,
        n_valid=nvalid_arr,
        n_delivered=ndeliv_arr,
        n_airborne=nair_arr,
        blocked_fraction=blocked_arr,
        gust_scale=gscale_arr,
        gust_clipped_fraction=gclip_arr,
        gust_model=args.gust_model,
        sigma_ratios=np.array(sigma_ratios),
        n_throws=args.n_throws,
        seed=args.seed,
        target_range_m=TARGET_RANGE_M,
        target_radius_m=TARGET_RADIUS_M,
        field_source=str(args.field),
        gusts_enabled=bool(not args.no_gusts),
        coeff_set=DEFAULT_COEFFS.name,
        schema_frames=SCHEMA_FRAMES,
        **_provenance(args),
    )
    print(f"wrote {args.out}")
    with np.errstate(invalid="ignore"):
        print(f"  success_rate: min={np.nanmin(succ_arr):.3f} "
              f"mean={np.nanmean(succ_arr):.3f} max={np.nanmax(succ_arr):.3f}")
        print(f"  blocked by buildings: mean "
              f"{100*np.nanmean(blocked_arr):.1f}% of throws")
    n_air = int(nair_arr.sum())
    if n_air:
        print(f"  WARNING: {n_air} throws were still airborne at "
              f"t_max={T_MAX} s and were excluded from all statistics")


if __name__ == "__main__":
    main()
