"""
Compute Q(x, y, theta_wind) from P3's CFD output and P4's landing stats.

Q in [0, 1] per the manuscript definition. Saves a single .npz keyed by
launch grid + heading.
"""

from __future__ import annotations
import argparse
import numpy as np


def wind_factor(U_2m, U_ref):
    mag = np.linalg.norm(U_2m, axis=-1)
    return np.clip(1.0 - mag / U_ref, 0.0, 1.0)


def turbulence_factor(k_2m, U_2m):
    sigma = np.sqrt(2.0 * np.maximum(k_2m, 0.0) / 3.0)
    mag = np.linalg.norm(U_2m, axis=-1)
    I = np.where(mag > 1e-3, sigma / np.maximum(mag, 1e-3), 0.0)
    return np.clip(1.0 - I, 0.0, 1.0)


def success_rate_from_landings(landing_mean, ellipse_axes, target_radius=3.0):
    """Crude success rate proxy: P(landing inside disc of radius target_radius
    centred at mean) assuming bivariate Gaussian with given ellipse semiaxes.
    """
    # 2D-Gaussian probability of falling inside a disc: approximate as
    # min(1, (target_radius / mean_axis)^2)
    a, b = ellipse_axes[..., 0], ellipse_axes[..., 1]
    mean_axis = 0.5 * (a + b)
    return np.clip((target_radius / np.maximum(mean_axis, 0.1)) ** 2, 0.0, 1.0)


def compute_Q(fields_npz, landings_npz, U_ref):
    F = np.load(fields_npz, allow_pickle=False)
    L = np.load(landings_npz, allow_pickle=False)

    headings = L["headings"]
    launch_x = L["launch_x"]
    launch_y = L["launch_y"]
    Q = np.zeros((len(launch_x), len(launch_y), len(headings)))

    for ih, h in enumerate(headings):
        key = f"dir_{int(h):03d}"
        # 2 m AGL slice from CFD field
        zk = int(np.argmin(np.abs(F["z"] - 2.0)))
        U_2m = F[f"{key}__U"][:, :, zk, :]
        k_2m = F[f"{key}__k"][:, :, zk]

        # Interpolate to launch grid
        from scipy.interpolate import RegularGridInterpolator
        Ui = [RegularGridInterpolator((F["x"], F["y"]), U_2m[..., i],
                                      bounds_error=False, fill_value=0.0)
              for i in range(3)]
        ki = RegularGridInterpolator((F["x"], F["y"]), k_2m,
                                     bounds_error=False, fill_value=0.0)

        XX, YY = np.meshgrid(launch_x, launch_y, indexing="ij")
        pts = np.column_stack([XX.ravel(), YY.ravel()])
        Ulaunch = np.stack([ui(pts) for ui in Ui], axis=-1).reshape(
            len(launch_x), len(launch_y), 3)
        klaunch = ki(pts).reshape(len(launch_x), len(launch_y))

        wf = wind_factor(Ulaunch, U_ref)
        tf = turbulence_factor(klaunch, Ulaunch)

        # success rate averaged across throw types and bearings
        ell = L["ellipse_axes"][:, :, ih, ...]   # (nx, ny, nt, nb, 2)
        mean = L["landing_mean"][:, :, ih, ...]  # (nx, ny, nt, nb, 3)
        sr_per = success_rate_from_landings(mean, ell)
        # average over throw types and bearings
        sr = sr_per.mean(axis=(-1, -2))

        Q[:, :, ih] = wf * tf * sr

    return Q, launch_x, launch_y, headings


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fields", required=True)
    ap.add_argument("--landings", required=True)
    ap.add_argument("--u-ref", type=float, default=5.0)
    ap.add_argument("--out", default="Q.npz")
    args = ap.parse_args()

    Q, x, y, h = compute_Q(args.fields, args.landings, args.u_ref)
    np.savez_compressed(args.out, Q=Q, x=x, y=y, headings=h)
    print(f"wrote {args.out}: Q shape {Q.shape}, range [{Q.min():.3f}, {Q.max():.3f}]")


if __name__ == "__main__":
    main()
