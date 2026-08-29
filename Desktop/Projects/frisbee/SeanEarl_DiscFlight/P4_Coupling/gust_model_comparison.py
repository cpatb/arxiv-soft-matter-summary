"""
Measure turbulence-induced landing dispersion, isolated by common random
numbers, for each gust correlation model. This is the GATE G1 instrument.

WHAT IS BEING MEASURED
----------------------
Total landing scatter is dominated by release-parameter scatter: ~98% of the
per-throw variance is the thrower, not the wind. Comparing total sigma_land
between two gust models therefore measures almost nothing. Instead every throw
is run TWICE from an identical release draw -- once with gusts, once with the
RANS mean flow alone -- and the paired displacement

    Delta_i = landing_i(mean + gust)  -  landing_i(mean only)

is the landing-point displacement caused by turbulence for that specific throw.
Its RMS, per horizontal axis,

    sigma_turb = sqrt( mean_i |Delta_i|^2 / 2 )

is the turbulence-induced dispersion. Because the arms share release draws
exactly (`CorrelatedGusts` consumes a fixed 3*N_MODES uniforms whatever the
model or field, so the RNG stream stays in lockstep), this estimator has a tiny
fraction of the variance of the unpaired difference and needs no large sample.

The G1 question is the ratio sigma_turb(temporal) / sigma_turb(sweeping):
how much did the Eulerian-decay-time gust model overpredict?

Usage
-----
    python3 gust_model_comparison.py --field ../P3_CFD/flow_fields.npz \
        --n-throws 40 --jobs 20 --out gust_model_comparison.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "P1_AnalyticalModel"))
import disc_6dof as d6                              # noqa: E402
from coefficients import DEFAULT_COEFFS, ULTRASTAR  # noqa: E402
from trajectories import sample_release             # noqa: E402

import synthetic_turbulence as st                   # noqa: E402
from synthetic_turbulence import (CorrelatedGusts, wind_with_gusts,  # noqa: E402
                                  wind_mean_only)
from monte_carlo import (crn_spawn_key, _get_field, MC_RTOL, MC_ATOL,  # noqa: E402
                         T_MAX, GROUND_TOL_M)


def paired_cell(launch, bearing_deg, throw_type, field, seed,
                n_throws, gust_model, sigma_ratios):
    """Run one cell twice from identical release draws.

    Returns (landings_gust, landings_mean, realised_scale, n_clipped), each
    landing array (n_throws, 3) with NaN where the throw did not reach ground.
    """
    launch = np.asarray(launch, dtype=float)
    out = {}
    scale = np.nan
    n_clipped = 0

    for arm, use_gusts in (("gust", True), ("mean", False)):
        # A FRESH generator seeded identically for each arm. This is what makes
        # the pairing exact: arm 2 replays arm 1's release draws.
        rng = np.random.default_rng(seed)
        land = np.full((n_throws, 3), np.nan)
        for i in range(n_throws):
            g = CorrelatedGusts(field, rng, t_max=T_MAX, ref_pos=launch,
                                model=gust_model, sigma_ratios=sigma_ratios)
            if arm == "gust":
                scale = g.scale
                n_clipped += int(g.clipped)
            wind = (wind_with_gusts if use_gusts else wind_mean_only)(field, g)
            r0, v0, q0, w0 = sample_release(
                throw_type=throw_type, heading_deg=bearing_deg,
                launch_pos=launch, rng=rng)
            res = d6.simulate(r0, v0, q0, w0,
                              coeffs=DEFAULT_COEFFS, geom=ULTRASTAR,
                              wind=wind, t_max=T_MAX,
                              rtol=MC_RTOL, atol=MC_ATOL, compute_diag=False)
            land[i] = res.r[-1]
        out[arm] = land

    return out["gust"], out["mean"], scale, n_clipped


def _run(task):
    (idx, launch, heading, bearing, throw, n_throws, field_path,
     entropy, gust_model, sigma_ratios) = task
    field = _get_field(field_path, heading)
    seed = np.random.SeedSequence(entropy=entropy, spawn_key=crn_spawn_key(idx))
    lg, lm, scale, nclip = paired_cell(
        launch, bearing, throw, field, seed, n_throws, gust_model,
        sigma_ratios)

    # A pair is usable only if BOTH arms landed. Airborne-at-t_max throws are
    # dropped, not scored at altitude (see monte_carlo.GROUND_TOL_M).
    ok = (np.isfinite(lg[:, 0]) & (np.abs(lg[:, 2]) <= GROUND_TOL_M)
          & np.isfinite(lm[:, 0]) & (np.abs(lm[:, 2]) <= GROUND_TOL_M))
    d = lg[ok, :2] - lm[ok, :2]
    return dict(idx=list(idx), heading=int(heading), bearing=float(bearing),
                throw=throw, n_pairs=int(ok.sum()),
                n_dropped=int((~ok).sum()),
                sum_sq=float((d ** 2).sum()),
                sum_d=[float(d[:, 0].sum()), float(d[:, 1].sum())],
                scale=float(scale), n_clipped=int(nclip))


def summarise(rows):
    """Pooled sigma_turb over all cells, per horizontal axis."""
    n = sum(r["n_pairs"] for r in rows)
    ss = sum(r["sum_sq"] for r in rows)
    if n == 0:
        return dict(sigma_turb=float("nan"), n_pairs=0)
    return dict(
        # per-axis RMS: |Delta|^2 pools both axes, hence the factor 2
        sigma_turb=float(np.sqrt(ss / (2 * n))),
        sigma_turb_radial=float(np.sqrt(ss / n)),
        n_pairs=int(n),
        n_dropped=int(sum(r["n_dropped"] for r in rows)),
        scale_mean=float(np.mean([r["scale"] for r in rows])),
        scale_min=float(np.min([r["scale"] for r in rows])),
        scale_max=float(np.max([r["scale"] for r in rows])),
        clipped_cells=int(sum(1 for r in rows if r["n_clipped"] > 0)),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--field", default="../P3_CFD/flow_fields.npz")
    ap.add_argument("--launch-grid", nargs=6, type=float,
                    metavar=("XA", "XB", "DX", "YA", "YB", "DY"),
                    default=[-25, 40, 13, -5, 5, 10])
    ap.add_argument("--launch-z", type=float, default=1.5)
    ap.add_argument("--headings", nargs="+", type=int, default=[90, 210, 270])
    ap.add_argument("--throws", nargs="+", default=["right_backhand"])
    ap.add_argument("--throw-bearings", nargs="+", type=float,
                    default=[0, 120, 240])
    ap.add_argument("--n-throws", type=int, default=40)
    ap.add_argument("--seed", type=int, default=20260726)
    ap.add_argument("--jobs", type=int, default=0)
    ap.add_argument("--models", nargs="+", default=["temporal", "sweeping"])
    ap.add_argument("--isotropic", action="store_true",
                    help="isotropic sigma in BOTH models, so the comparison "
                         "isolates the correlation-model change alone")
    ap.add_argument("--out", default="gust_model_comparison.json")
    args = ap.parse_args()

    xa, xb, dx, ya, yb, dy = args.launch_grid
    xs = np.arange(xa, xb + dx / 2, dx)
    ys = np.arange(ya, yb + dy / 2, dy)
    ratios = st.ISOTROPIC_RATIOS if args.isotropic else st.SIGMA_RATIOS

    results = {}
    n_jobs = args.jobs or max(1, (os.cpu_count() or 2) - 2)

    for model in args.models:
        tasks = []
        for ix, x_l in enumerate(xs):
            for iy, y_l in enumerate(ys):
                for i_h, h in enumerate(args.headings):
                    for it, throw in enumerate(args.throws):
                        for ib, bear in enumerate(args.throw_bearings):
                            tasks.append((
                                (ix, iy, i_h, it, ib),
                                (float(x_l), float(y_l), args.launch_z),
                                int(h), float(bear), throw, args.n_throws,
                                args.field, args.seed, model, ratios))

        print(f"[{model}] {len(tasks)} cells x {args.n_throws} paired throws "
              f"= {2*len(tasks)*args.n_throws:,} trajectories on {n_jobs} workers",
              flush=True)
        t0 = time.perf_counter()
        if n_jobs == 1:
            rows = [_run(t) for t in tasks]
        else:
            import multiprocessing as mp
            with mp.Pool(n_jobs) as pool:
                rows = list(pool.imap_unordered(_run, tasks, chunksize=1))
        el = time.perf_counter() - t0
        s = summarise(rows)
        s["elapsed_s"] = el
        s["model"] = model
        results[model] = dict(summary=s, cells=rows)
        print(f"[{model}] sigma_turb = {s['sigma_turb']:.4f} m per axis "
              f"({s['sigma_turb_radial']:.4f} m radial) from {s['n_pairs']:,} "
              f"pairs, scale {s['scale_min']:.2f}-{s['scale_max']:.2f} "
              f"(mean {s['scale_mean']:.2f}), {el/60:.1f} min", flush=True)
        if s["n_dropped"]:
            print(f"[{model}] dropped {s['n_dropped']} pairs "
                  "(a throw in one arm did not reach the ground)", flush=True)
        if s["clipped_cells"]:
            print(f"[{model}] WARNING: scale clipped in "
                  f"{s['clipped_cells']} cells", flush=True)

    print()
    print("=" * 68)
    if len(args.models) == 2:
        a, b = args.models
        ra = results[a]["summary"]["sigma_turb"]
        rb = results[b]["summary"]["sigma_turb"]
        ratio = ra / rb if rb > 0 else float("inf")
        print(f"sigma_turb({a})   = {ra:.4f} m per axis")
        print(f"sigma_turb({b})   = {rb:.4f} m per axis")
        print(f"RATIO {a}/{b} = {ratio:.3f}")
        print()
        print("GATE G1 acceptance: dispersion falls by a factor of 2.0-2.6")
        if ratio < 1.2:
            verdict = ("FAIL-NULL: change < 20%. The defect was harmless and "
                       "the panels' central prediction did not reproduce.")
        elif ratio < 2.0:
            verdict = f"BELOW BAND: fell {ratio:.2f}x, less than the predicted 2.0-2.6x"
        elif ratio <= 2.6:
            verdict = f"PASS: fell {ratio:.2f}x, inside the predicted 2.0-2.6x"
        elif ratio <= 3.0:
            verdict = f"ABOVE BAND: fell {ratio:.2f}x, above the predicted 2.0-2.6x"
        else:
            verdict = (f"FAIL-HIGH: fell {ratio:.2f}x (>3x). The disc is not "
                       "behaving as a linear filter.")
        print(f"VERDICT: {verdict}")
        results["gate_G1"] = dict(ratio=ratio, verdict=verdict,
                                  sigma_turb=dict(zip(args.models, (ra, rb))))
    print("=" * 68)

    results["config"] = dict(
        n_throws=args.n_throws, seed=args.seed, field=str(args.field),
        launch_x=list(map(float, xs)), launch_y=list(map(float, ys)),
        headings=args.headings, throws=args.throws,
        bearings=args.throw_bearings, sigma_ratios=list(ratios),
        isotropic=bool(args.isotropic),
    )
    Path(args.out).write_text(json.dumps(results, indent=1))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
