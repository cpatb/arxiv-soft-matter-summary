"""
Generate a synthetic `flow_fields.npz` matching P3's output schema.

Purpose: exercise and pin down the P3 -> P4 -> P5 data contract before any CFD
exists, so that when the real OpenFOAM output arrives only the file changes.
This is a TEST DOUBLE. It is not a physical model of the site and must never be
used to produce results for the paper -- `main()` stamps `is_stub=True` into
the archive, and `throwability.py` refuses to treat a stub as real data.

The synthetic field is deliberately non-trivial so that downstream maps have
structure worth plotting:
  * a log-law vertical profile, U(z) = (u*/kappa) ln(z/z0), rotated to the
    requested inflow heading;
  * a sheltered wake behind a rectangular obstacle west of the field, with a
    velocity deficit that recovers downstream;
  * a turbulent kinetic energy field elevated in the wake shear layer.

Schema (consumed by P4 `CFDField` and P5 `throwability.py`):
    x, y, z            : 1-D coordinate axes, m
    dir_<hhh>__U       : (nx, ny, nz, 3) mean velocity, m/s
    dir_<hhh>__k       : (nx, ny, nz)    turbulent kinetic energy, m^2/s^2
    dir_<hhh>__eps     : (nx, ny, nz)    dissipation rate, m^2/s^3
where <hhh> is the inflow heading in degrees, zero-padded to 3 digits.
"""

from __future__ import annotations

import argparse

import numpy as np

KAPPA = 0.41           # von Karman constant
Z0_URBAN = 0.5         # m, roughness length used for the stub profile
U_REF = 5.0            # m/s at Z_REF
Z_REF = 10.0           # m

# Rectangular obstacle standing in for Mertz Hall, in local field coords (m).
OBSTACLE = dict(x0=-80.0, x1=-50.0, y0=-30.0, y1=30.0, height=22.0)


def log_profile(z, u_ref=U_REF, z_ref=Z_REF, z0=Z0_URBAN):
    """Neutral log-law wind speed. Clamped below z0 to avoid log(0)."""
    z = np.maximum(z, z0 * 1.001)
    u_star = u_ref * KAPPA / np.log(z_ref / z0)
    return (u_star / KAPPA) * np.log(z / z0)


def build_field(x, y, z, heading_deg):
    """Mean velocity, k and eps on the (x, y, z) grid for one inflow heading.

    `heading_deg` is the direction the wind blows TOWARD, measured
    counter-clockwise from +x, consistent with the trajectory model's heading
    convention.
    """
    nx, ny, nz = len(x), len(y), len(z)
    XX, YY, ZZ = np.meshgrid(x, y, z, indexing="ij")

    speed = log_profile(ZZ)
    th = np.radians(heading_deg)
    ehat = np.array([np.cos(th), np.sin(th), 0.0])

    # Streamwise coordinate, so "downstream of the obstacle" is well defined
    # for any heading.
    s = XX * ehat[0] + YY * ehat[1]

    # Obstacle centre and its streamwise station.
    cx = 0.5 * (OBSTACLE["x0"] + OBSTACLE["x1"])
    cy = 0.5 * (OBSTACLE["y0"] + OBSTACLE["y1"])
    half_w = 0.5 * max(OBSTACLE["x1"] - OBSTACLE["x0"],
                       OBSTACLE["y1"] - OBSTACLE["y0"])
    s_obs = cx * ehat[0] + cy * ehat[1]
    # Cross-stream distance from the obstacle centreline.
    n_coord = -(XX - cx) * ehat[1] + (YY - cy) * ehat[0]

    ds = s - s_obs                       # >0 downstream of the obstacle
    H = OBSTACLE["height"]

    # Wake deficit: only downstream, decaying over ~10 H, confined laterally to
    # the obstacle width plus spreading, and to heights below ~1.5 H.
    downstream = ds > 0
    decay = np.exp(-np.maximum(ds, 0.0) / (10.0 * H))
    lateral = np.exp(-(n_coord / (half_w + 0.15 * np.maximum(ds, 0.0))) ** 2)
    vertical = np.exp(-(ZZ / (1.5 * H)) ** 2)
    deficit = 0.65 * decay * lateral * vertical * downstream

    speed_wake = speed * (1.0 - deficit)

    U = np.zeros((nx, ny, nz, 3))
    U[..., 0] = speed_wake * ehat[0]
    U[..., 1] = speed_wake * ehat[1]
    # A weak vertical component in the recovering shear layer.
    U[..., 2] = 0.05 * speed * deficit

    # Background TKE from the log law, boosted in the wake shear layer.
    u_star = U_REF * KAPPA / np.log(Z_REF / Z0_URBAN)
    k_bg = np.full_like(speed, u_star ** 2 / np.sqrt(0.09))
    k = k_bg * (1.0 + 4.0 * deficit)
    eps = np.maximum(u_star ** 3 / (KAPPA * np.maximum(ZZ, Z0_URBAN)), 1e-6)
    eps = eps * (1.0 + 4.0 * deficit)

    return U, k, eps


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--headings", nargs="+", type=int,
                    default=[0, 90, 180])
    ap.add_argument("--extent", nargs=4, type=float,
                    metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
                    default=[-120.0, 120.0, -120.0, 120.0])
    ap.add_argument("--dxy", type=float, default=2.0)
    ap.add_argument("--z", nargs="+", type=float,
                    default=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    ap.add_argument("--out", default="stub_flow_fields.npz")
    args = ap.parse_args()

    xmin, xmax, ymin, ymax = args.extent
    x = np.arange(xmin, xmax + args.dxy / 2, args.dxy)
    y = np.arange(ymin, ymax + args.dxy / 2, args.dxy)
    z = np.array(args.z, dtype=float)

    out = dict(x=x, y=y, z=z, is_stub=True,
               headings=np.array(args.headings))
    for h in args.headings:
        U, k, eps = build_field(x, y, z, h)
        out[f"dir_{h:03d}__U"] = U
        out[f"dir_{h:03d}__k"] = k
        out[f"dir_{h:03d}__eps"] = eps

    np.savez_compressed(args.out, **out)
    print(f"wrote {args.out}")
    print(f"  grid {len(x)} x {len(y)} x {len(z)}, headings {args.headings}")
    for h in args.headings:
        mag = np.linalg.norm(out[f"dir_{h:03d}__U"], axis=-1)
        print(f"  dir_{h:03d}: |U| in [{mag.min():.2f}, {mag.max():.2f}] m/s, "
              f"k in [{out[f'dir_{h:03d}__k'].min():.3f}, "
              f"{out[f'dir_{h:03d}__k'].max():.3f}]")
    print("\n  NOTE: this is a STUB field for pipeline testing, not CFD.")


if __name__ == "__main__":
    main()
