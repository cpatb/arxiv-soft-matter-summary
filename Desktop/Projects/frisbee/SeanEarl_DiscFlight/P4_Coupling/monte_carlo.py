"""
Monte Carlo trajectory sweep over launch positions, headings, throw types,
and release-parameter scatter.

For each (launch_pos, heading_wind, throw_type, throw_bearing) tuple,
draws n_throws samples from the throw-type distribution, runs each
through the 6-DOF integrator with CFD-mean + gust wind, and aggregates
landing statistics.
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "P1_AnalyticalModel"))
import disc_6dof as d6                              # noqa: E402
from coefficients import HUMMEL_SHORT, ULTRASTAR    # noqa: E402
from trajectories import sample_release             # noqa: E402

from cfd_field import StubField, CFDField           # noqa: E402
from synthetic_turbulence import WhiteNoiseGusts, wind_with_gusts  # noqa: E402


def sweep_one_launch_point(
    launch, heading_deg, throw_type, n_throws, field,
    bearing_deg: float = 0.0,
    seed: int = 0,
    t_max: float = 8.0,
):
    """Return aggregated landing statistics for one (launch, heading, throw) cell."""
    rng = np.random.default_rng(seed)
    landings = np.zeros((n_throws, 3))
    alpha_lands = np.zeros(n_throws)

    gusts = WhiteNoiseGusts(field, rng)
    wind = wind_with_gusts(field, gusts)

    for i in range(n_throws):
        r0, v0, q0, w0 = sample_release(
            throw_type=throw_type,
            heading_deg=bearing_deg,
            launch_pos=launch,
            rng=rng,
        )
        res = d6.simulate(r0, v0, q0, w0,
                          coeffs=HUMMEL_SHORT, geom=ULTRASTAR,
                          wind=wind, t_max=t_max)
        landings[i] = res.r[-1]
        alpha_lands[i] = res.diag[-1].get("alpha_deg", np.nan)

    mean = landings.mean(axis=0)
    cov = np.cov(landings[:, :2].T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    # 90% containment for 2D Gaussian: r^2/sigma^2 = -2 ln(0.1)
    chi2_90 = -2.0 * np.log(0.1)
    ellipse_axes = np.sqrt(eigvals * chi2_90)

    return dict(
        landings=landings,
        mean=mean, cov=cov,
        ellipse_axes=ellipse_axes,
        alpha_landing_mean=float(np.nanmean(alpha_lands)),
    )


def parse_grid(args_grid):
    xa, xb, dx, ya, yb, dy = args_grid
    xs = np.arange(xa, xb + dx / 2, dx)
    ys = np.arange(ya, yb + dy / 2, dy)
    return xs, ys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--field", default=None,
                    help="P3 flow_fields.npz; if omitted, use StubField with no wind")
    ap.add_argument("--launch-grid", nargs=6, type=float,
                    metavar=("XA","XB","DX","YA","YB","DY"),
                    default=[-25, 25, 10, -25, 25, 10])
    ap.add_argument("--launch-z", type=float, default=1.5)
    ap.add_argument("--headings", nargs="+", type=int, required=True)
    ap.add_argument("--throws", nargs="+",
                    default=["right_backhand", "right_forehand"])
    ap.add_argument("--throw-bearings", nargs="+", type=float,
                    default=[0])
    ap.add_argument("--n-throws", type=int, default=200)
    ap.add_argument("--out", default="landing_stats.npz")
    args = ap.parse_args()

    xs, ys = parse_grid(args.launch_grid)
    headings = args.headings
    throws = args.throws
    bearings = args.throw_bearings

    n_lx, n_ly = len(xs), len(ys)
    n_h, n_t, n_b = len(headings), len(throws), len(bearings)

    mean_arr = np.zeros((n_lx, n_ly, n_h, n_t, n_b, 3))
    cov_arr = np.zeros((n_lx, n_ly, n_h, n_t, n_b, 2, 2))
    ax_arr = np.zeros((n_lx, n_ly, n_h, n_t, n_b, 2))
    alpha_arr = np.zeros((n_lx, n_ly, n_h, n_t, n_b))

    for i_h, h in enumerate(headings):
        field = (StubField(np.zeros(3)) if args.field is None
                 else CFDField(args.field, heading=h))
        for ix, x_l in enumerate(xs):
            for iy, y_l in enumerate(ys):
                for it, throw in enumerate(throws):
                    for ib, bear in enumerate(bearings):
                        seed = hash((h, ix, iy, it, ib)) & 0xFFFFFFFF
                        stats = sweep_one_launch_point(
                            launch=(x_l, y_l, args.launch_z),
                            heading_deg=h, throw_type=throw,
                            bearing_deg=bear,
                            n_throws=args.n_throws,
                            field=field, seed=seed,
                        )
                        mean_arr[ix, iy, i_h, it, ib] = stats["mean"]
                        cov_arr[ix, iy, i_h, it, ib] = stats["cov"]
                        ax_arr[ix, iy, i_h, it, ib] = stats["ellipse_axes"]
                        alpha_arr[ix, iy, i_h, it, ib] = stats["alpha_landing_mean"]
        print(f"finished heading {h}°")

    np.savez_compressed(
        args.out,
        launch_x=xs, launch_y=ys, launch_z=args.launch_z,
        headings=np.array(headings),
        throws=np.array(throws),
        bearings=np.array(bearings),
        landing_mean=mean_arr,
        landing_cov=cov_arr,
        ellipse_axes=ax_arr,
        alpha_landing=alpha_arr,
    )
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
