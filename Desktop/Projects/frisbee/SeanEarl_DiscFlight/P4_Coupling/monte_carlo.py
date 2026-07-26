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
from synthetic_turbulence import CorrelatedGusts, wind_with_gusts  # noqa: E402

# Canonical scoring target, per the manuscript: a 25 m throw is "delivered" if
# it lands within 3 m of the intended target point.
TARGET_RANGE_M = 25.0
TARGET_RADIUS_M = 3.0

# Monte Carlo throws may be integrated at a looser tolerance than the
# validation runs; the convergence study in P1 shows the landing point moves by
# ~0.02 mm, six orders of magnitude below TARGET_RADIUS_M.
MC_RTOL = 1e-4
MC_ATOL = 1e-6


def sweep_one_launch_point(
    launch, heading_deg, throw_type, n_throws, field,
    bearing_deg: float = 0.0,
    seed=0,
    t_max: float = 8.0,
    target_range_m: float = TARGET_RANGE_M,
    target_radius_m: float = TARGET_RADIUS_M,
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

    landings = np.full((n_throws, 3), np.nan)
    alpha_lands = np.full(n_throws, np.nan)

    for i in range(n_throws):
        # A fresh gust realization per throw, drawn up front so that the wind
        # is a deterministic function of (r, t) during the integration. The
        # integral time scale comes from k/eps at the launch point.
        gusts = CorrelatedGusts(field, rng, t_max=t_max, ref_pos=launch)
        wind = wind_with_gusts(field, gusts)

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
    valid = np.isfinite(landings[:, 0])

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
    good = landings[valid]

    # Intended target: target_range_m downrange along the throw bearing.
    bearing = np.radians(bearing_deg)
    target_xy = launch[:2] + target_range_m * np.array(
        [np.cos(bearing), np.sin(bearing)])

    if n_valid >= 2:
        mean = good.mean(axis=0)
        cov = np.cov(good[:, :2].T)
        eigvals = np.linalg.eigvalsh(cov)
        # 90% containment for an isotropic-ish 2D Gaussian: r^2/sigma^2 =
        # -2 ln(0.1). eigvals can go slightly negative from round-off.
        chi2_90 = -2.0 * np.log(0.1)
        ellipse_axes = np.sqrt(np.maximum(eigvals, 0.0) * chi2_90)
        d_target = np.linalg.norm(good[:, :2] - target_xy, axis=1)
        hit = (d_target <= target_radius_m) & (~blocked[valid])
        success_rate = float(np.mean(hit))
    else:
        mean = np.full(3, np.nan)
        cov = np.full((2, 2), np.nan)
        ellipse_axes = np.full(2, np.nan)
        success_rate = np.nan

    return dict(
        landings=landings,
        mean=mean, cov=cov,
        ellipse_axes=ellipse_axes,
        success_rate=success_rate,
        alpha_landing_mean=float(np.nanmean(alpha_lands)) if n_valid else np.nan,
        target_xy=target_xy,
        n_valid=n_valid,
        n_blocked=n_blocked,
        blocked_fraction=n_blocked / max(n_throws, 1),
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


def _run_cell(task):
    """Worker entry point: run one cell and return its summary scalars."""
    (idx, launch, heading, throw, bearing, n_throws, field_path,
     seed_entropy, t_max) = task
    field = _get_field(field_path, heading)
    seed = np.random.SeedSequence(entropy=seed_entropy, spawn_key=tuple(idx))
    stats = sweep_one_launch_point(
        launch=launch, heading_deg=heading, throw_type=throw,
        bearing_deg=bearing, n_throws=n_throws, field=field,
        seed=seed, t_max=t_max,
    )
    return idx, (stats["mean"], stats["cov"], stats["ellipse_axes"],
                 stats["success_rate"], stats["alpha_landing_mean"],
                 stats["n_valid"], stats["blocked_fraction"])


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
    ap.add_argument("--out", default="landing_stats.npz")
    args = ap.parse_args()

    xs, ys = parse_grid(args.launch_grid)
    headings, throws, bearings = args.headings, args.throws, args.throw_bearings
    n_lx, n_ly = len(xs), len(ys)
    n_h, n_t, n_b = len(headings), len(throws), len(bearings)

    shape = (n_lx, n_ly, n_h, n_t, n_b)
    mean_arr = np.full(shape + (3,), np.nan)
    cov_arr = np.full(shape + (2, 2), np.nan)
    ax_arr = np.full(shape + (2,), np.nan)
    succ_arr = np.full(shape, np.nan)
    alpha_arr = np.full(shape, np.nan)
    nvalid_arr = np.zeros(shape, dtype=int)
    blocked_arr = np.full(shape, np.nan)

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
                            args.field, args.seed, 8.0,
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

    for idx, (mean, cov, axes, succ, alpha, nvalid, blk) in results:
        mean_arr[idx] = mean
        cov_arr[idx] = cov
        ax_arr[idx] = axes
        succ_arr[idx] = succ
        alpha_arr[idx] = alpha
        nvalid_arr[idx] = nvalid
        blocked_arr[idx] = blk

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
        success_rate=succ_arr,
        alpha_landing=alpha_arr,
        n_valid=nvalid_arr,
        blocked_fraction=blocked_arr,
        n_throws=args.n_throws,
        seed=args.seed,
        target_range_m=TARGET_RANGE_M,
        target_radius_m=TARGET_RADIUS_M,
        field_source=str(args.field),
        coeff_set=DEFAULT_COEFFS.name,
    )
    print(f"wrote {args.out}")
    with np.errstate(invalid="ignore"):
        print(f"  success_rate: min={np.nanmin(succ_arr):.3f} "
              f"mean={np.nanmean(succ_arr):.3f} max={np.nanmax(succ_arr):.3f}")
        print(f"  blocked by buildings: mean "
              f"{100*np.nanmean(blocked_arr):.1f}% of throws")


if __name__ == "__main__":
    main()
