"""Does the production field's z = 6 m ceiling change any landing?

THE DEFECT
----------
`P3_CFD/flow_fields.npz` -- the field the manuscript's 97,200-throw sweep was
flown through -- is sampled at z = 1, 2, 3, 4, 5, 6 m only. `CFDField` CLAMPS
outside its box rather than extrapolating, which is the right choice, but it
means a disc above 6 m flies through wind frozen at its 6 m value. The
behind-the-thrower render's new z-coverage guard caught this on its first run:
a 30-throw backhand ensemble peaks at 6.6 m, with 1.4% of trajectory samples
above the ceiling.

THE TEST
--------
Not "production field vs vertical field": those differ in HORIZONTAL resolution
too (1 m vs 2 m), so any difference would confound the ceiling with the
sampling. Instead this truncates `vertical_fields.npz` at 6 m and compares it
against the same file untruncated -- identical solution, identical horizontal
grid, identical release draws and gust phases, one difference.

Paired, so the answer is a per-throw displacement rather than a difference of
two scatters.

WHAT A NULL AND A NON-NULL MEAN
-------------------------------
A null result closes the question: the ceiling is harmless for throws of this
length and the production numbers stand. A non-null result does NOT
automatically invalidate the production sweep -- it bounds the error the ceiling
introduces, and the sweep also covers launch points and throw types not sampled
here. Either way the number is what the paper needs to be able to quote.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "P1_AnalyticalModel"))
sys.path.insert(0, str(ROOT / "P4_Coupling"))
sys.path.insert(0, str(ROOT / "P5_Synthesis"))


def truncate_field(src, dst, z_max=6.0):
    """Write a copy of `src` keeping only levels with z <= z_max."""
    F = np.load(src, allow_pickle=False)
    z = F["z"]
    keep = z <= z_max + 1e-9
    if keep.all():
        raise SystemExit(f"{src} already stops at or below z = {z_max}")
    out = {}
    for name in F.files:
        a = F[name]
        if name == "z":
            out[name] = a[keep]
        elif name.endswith("__U") or name.endswith("__k") or name.endswith("__eps"):
            out[name] = a[:, :, keep, ...]
        elif name.endswith("__valid3d"):
            out[name] = a[:, :, keep]
        else:
            out[name] = a
    np.savez(dst, **out)
    return float(z[keep].max()), int(keep.sum()), int(len(z))


def fly(field_path, heading, throw, launch, bearing, n_throws, seed=20260830,
        t_max=15.0):
    import disc_6dof as d6
    import monte_carlo as mc
    import synthetic_turbulence as st
    from cfd_field import CFDField
    from throw_catalogue import build_catalogue

    dist = build_catalogue()[throw]
    field = CFDField(str(field_path), heading=heading)
    launch = np.asarray(launch, dtype=float)
    lands, apex, above = [], [], []
    for i in range(n_throws):
        rng = np.random.default_rng(
            np.random.SeedSequence(entropy=seed, spawn_key=(int(i),)))
        r0, v0, q0, w0 = d6.build_initial_conditions(
            launch_pos=launch,
            speed=rng.normal(dist["speed_mean"], dist["speed_std"]),
            heading_deg=bearing,
            elevation_deg=rng.normal(dist["elevation_mean"], dist["elevation_std"]),
            release_alpha_deg=rng.normal(dist["release_alpha_mean"],
                                         dist["release_alpha_std"]),
            roll_deg=rng.normal(dist["roll_mean"], dist["roll_std"]),
            spin_rad_s=max(20.0, rng.normal(dist["spin_mean"], dist["spin_std"])),
            handedness=dist["handedness"])
        grng = np.random.default_rng(
            np.random.SeedSequence(entropy=seed + 1, spawn_key=(int(i),)))
        g = st.CorrelatedGusts(field, grng, t_max=t_max, ref_pos=launch,
                               model="sweeping", sigma_ratios=st.SIGMA_RATIOS)
        res = d6.simulate(r0, v0, q0, w0, coeffs=mc.DEFAULT_COEFFS,
                          geom=mc.ULTRASTAR, wind=st.wind_with_gusts(field, g),
                          t_max=t_max, rtol=mc.MC_RTOL, atol=mc.MC_ATOL,
                          compute_diag=False)
        r = np.asarray(res.r)
        lands.append(r[-1, :2])
        apex.append(float(r[:, 2].max()))
        above.append(float((r[:, 2] > 6.0).mean()))
    return np.asarray(lands), np.asarray(apex), np.asarray(above)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--field", default=str(ROOT / "P3_CFD" / "vertical_fields.npz"))
    ap.add_argument("--headings", nargs="+", type=int, default=[90, 210, 270])
    ap.add_argument("--throws", nargs="+", default=["right_backhand",
                                                    "right_forehand"])
    ap.add_argument("--n-throws", type=int, default=64)
    ap.add_argument("--out", default=str(HERE / "out" / "z_ceiling.json"))
    args = ap.parse_args()

    import site_view as sv
    trunc = HERE / "out" / "_truncated_z6.npz"
    trunc.parent.mkdir(parents=True, exist_ok=True)
    zmax, nk, nz = truncate_field(args.field, trunc)
    print(f"truncated {Path(args.field).name}: {nz} levels -> {nk}, "
          f"top z = {zmax} m")

    rows = []
    for h in args.headings:
        for th in args.throws:
            full_l, apex, above = fly(args.field, h, th, sv.LAUNCH, sv.EAST,
                                      args.n_throws)
            cut_l, _, _ = fly(trunc, h, th, sv.LAUNCH, sv.EAST, args.n_throws)
            d = cut_l - full_l
            dist = np.hypot(d[:, 0], d[:, 1])
            row = dict(heading=h, throw=th, n=args.n_throws,
                       apex_mean=float(apex.mean()), apex_max=float(apex.max()),
                       frac_throws_over_6m=float((apex > 6.0).mean()),
                       frac_samples_over_6m=float(above.mean()),
                       shift_mean_m=float(dist.mean()),
                       shift_max_m=float(dist.max()),
                       shift_p95_m=float(np.percentile(dist, 95)),
                       group_centroid_shift_m=float(
                           np.hypot(*(cut_l.mean(axis=0) - full_l.mean(axis=0)))),
                       range_mean_m=float(np.hypot(
                           *(full_l.mean(axis=0) - sv.LAUNCH[:2]))))
            rows.append(row)
            print(f"  heading {h:>3} {th:<16} apex {row['apex_mean']:.2f} m "
                  f"(max {row['apex_max']:.2f}), {100*row['frac_throws_over_6m']:.0f}% "
                  f"of throws break 6 m")
            print(f"      per-throw landing shift {row['shift_mean_m']:.3f} m "
                  f"mean, {row['shift_max_m']:.3f} m max, "
                  f"{row['shift_p95_m']:.3f} m p95; group centroid moves "
                  f"{row['group_centroid_shift_m']:.3f} m "
                  f"on a {row['range_mean_m']:.1f} m throw")

    worst = max(r["shift_mean_m"] for r in rows)
    rng = min(r["range_mean_m"] for r in rows)
    print(f"\nWorst mean per-throw shift {worst:.3f} m = "
          f"{100*worst/rng:.2f}% of the shortest mean range.")
    Path(args.out).write_text(json.dumps(dict(
        field=args.field, truncated_to_z=zmax, rows=rows,
        worst_mean_shift_m=worst), indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
