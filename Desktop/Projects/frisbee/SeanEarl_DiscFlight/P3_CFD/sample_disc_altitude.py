"""
Sample converged CFD fields onto the regular disc-altitude grid P4 consumes,
and write `flow_fields.npz` in the schema P4/P5 expect.

METHOD
------
Uses OpenFOAM's own `sets` sampling with a `cloud` set: we hand it the exact
grid points we want and it returns interpolated values there. This avoids the
`foamToVTK` + pyvista route (a heavyweight extra dependency and a lossy
round-trip) and any hand-rolled unstructured interpolation.

Points that fall INSIDE a building have no owning cell. OpenFOAM drops them
from the output rather than returning a sentinel, so returned coordinates are
matched back onto grid indices rather than assuming the output is complete and
ordered. Dropped points become NaN and are then filled by nearest neighbour,
because P4's integrator must be able to evaluate the wind anywhere in the box.
The filled fraction is reported and stored: a large value would mean the launch
grid overlaps buildings and the throwability map is partly meaningless there.

NOTE ON THE ROTATION. Each case is solved with the geometry rotated so the wind
blows along +x. Sampled velocities are in the ROTATED frame and are rotated BACK
into site coordinates here, so P4/P5 always work in the site frame (x East,
y North) whatever the wind direction.

Output schema (consumed by P4 `CFDField` and P5 `throwability.py`):
    x, y, z            1-D axes, m (z is AGL)
    dir_<hhh>__U       (nx, ny, nz, 3) mean velocity, m/s
    dir_<hhh>__k       (nx, ny, nz)    TKE, m^2/s^2
    dir_<hhh>__eps     (nx, ny, nz)    dissipation, m^2/s^3
    headings           integer headings present
    is_stub            False
    meta_json          per-direction case metadata
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np

# epsilon = Cmu * k * omega  (standard k-omega <-> k-epsilon relation)
CMU = 0.09
FOAM_BASHRC = "/usr/lib/openfoam/openfoam2606/etc/bashrc"

SAMPLE_DICT = """/*--------------------------------*- C++ -*----------------------------------*/
FoamFile {{ version 2.0; format ascii; class dictionary; object sampleDisc; }}

type            sets;
libs            (sampling);
setFormat       raw;
interpolationScheme cellPoint;
fields          (U k omega);

sets
{{
    discGrid
    {{
        type    cloud;
        axis    xyz;
        points  (
{points}
        );
    }}
}}
"""


def write_sample_dict(case: Path, pts: np.ndarray):
    body = "\n".join(f"            ({p[0]:.4f} {p[1]:.4f} {p[2]:.4f})"
                     for p in pts)
    (case / "system" / "sampleDisc").write_text(SAMPLE_DICT.format(points=body))


def run_sampling(case: Path):
    cmd = (f"source {FOAM_BASHRC} && cd {case} && "
           "postProcess -func sampleDisc -latestTime")
    r = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"sampling failed in {case}:\n{r.stdout[-3000:]}\n"
                         f"{r.stderr[-2000:]}")
    return r.stdout


def read_samples(case: Path):
    """Return (coords (n,3), U (n,3), k (n,), omega (n,))."""
    base = case / "postProcessing" / "sampleDisc"
    times = sorted(base.glob("*"), key=lambda p: float(p.name))
    if not times:
        raise SystemExit(f"no sampled output under {base}")
    tdir = times[-1]

    def load(name):
        cands = (sorted(tdir.glob(f"discGrid_{name}.xy"))
                 or sorted(tdir.glob(f"discGrid*_{name}*"))
                 or sorted(tdir.glob(f"*{name}*")))
        if not cands:
            raise SystemExit(f"missing sampled field '{name}' in {tdir}: "
                             f"{[p.name for p in tdir.iterdir()]}")
        return np.loadtxt(cands[0])

    aU = load("U")
    return aU[:, :3], aU[:, 3:6], load("k")[:, 3], load("omega")[:, 3]


def rotate_xy(a, phi):
    """Rotate the x,y components of points or vectors by phi radians (CCW)."""
    c, s = np.cos(phi), np.sin(phi)
    out = np.array(a, dtype=float, copy=True)
    x, y = out[..., 0].copy(), out[..., 1].copy()
    out[..., 0] = c * x - s * y
    out[..., 1] = s * x + c * y
    return out


def fill_nearest(arr, valid):
    """Fill invalid entries from the nearest valid grid point."""
    from scipy.ndimage import distance_transform_edt
    if valid.all():
        return arr
    idx = distance_transform_edt(~valid, return_distances=False,
                                 return_indices=True)
    return arr[tuple(idx)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("cases", nargs="+", help="case directories")
    ap.add_argument("--out", default="flow_fields.npz")
    ap.add_argument("--extent", nargs=4, type=float,
                    metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
                    default=[-90.0, 110.0, -70.0, 70.0])
    ap.add_argument("--dxy", type=float, default=1.0)
    ap.add_argument("--z", nargs="+", type=float,
                    default=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    args = ap.parse_args()

    xmin, xmax, ymin, ymax = args.extent
    x = np.arange(xmin, xmax + args.dxy / 2, args.dxy)
    y = np.arange(ymin, ymax + args.dxy / 2, args.dxy)
    z = np.array(args.z, dtype=float)
    nx, ny, nz = len(x), len(y), len(z)

    XX, YY, ZZ = np.meshgrid(x, y, z, indexing="ij")
    site_pts = np.column_stack([XX.ravel(), YY.ravel(), ZZ.ravel()])
    print(f"grid {nx} x {ny} x {nz} = {len(site_pts):,} points "
          f"({args.dxy} m horizontal, z = {list(z)})")

    out = dict(x=x, y=y, z=z, is_stub=False)
    headings, meta = [], {}

    for cdir in args.cases:
        case = Path(cdir)
        info = json.loads((case / "case_info.json").read_text())
        theta = int(info["theta_wind_deg"])
        phi = np.radians(info["geometry_rotation_deg"])
        print(f"\n[dir_{theta:03d}] geometry rotation {np.degrees(phi):+.1f} deg")

        solve_pts = rotate_xy(site_pts, phi)
        write_sample_dict(case, solve_pts)
        run_sampling(case)
        coords, U_rot, kk, om = read_samples(case)
        print(f"  sampled {len(coords):,} / {len(site_pts):,} points "
              f"({100*len(coords)/len(site_pts):.1f}%)")

        site_back = rotate_xy(coords, -phi)
        ix = np.rint((site_back[:, 0] - xmin) / args.dxy).astype(int)
        iy = np.rint((site_back[:, 1] - ymin) / args.dxy).astype(int)
        iz = np.abs(site_back[:, 2][:, None] - z[None, :]).argmin(axis=1)
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        ix, iy, iz = ix[ok], iy[ok], iz[ok]

        Ug = np.full((nx, ny, nz, 3), np.nan)
        kg = np.full((nx, ny, nz), np.nan)
        eg = np.full((nx, ny, nz), np.nan)

        Ug[ix, iy, iz] = rotate_xy(U_rot[ok], -phi)
        kg[ix, iy, iz] = kk[ok]
        eg[ix, iy, iz] = CMU * kk[ok] * om[ok]

        valid = np.isfinite(kg)
        frac = 1.0 - float(valid.mean())
        print(f"  {100*frac:.2f}% of grid points unresolved (inside "
              f"buildings) -> nearest-neighbour filled")
        if frac > 0.25:
            print("  ** WARNING: over a quarter of the sampling grid lies "
                  "inside buildings **")

        for c in range(3):
            Ug[..., c] = fill_nearest(Ug[..., c], valid)
        kg = fill_nearest(kg, valid)
        eg = np.maximum(fill_nearest(eg, valid), 1e-6)

        out[f"dir_{theta:03d}__U"] = Ug
        out[f"dir_{theta:03d}__k"] = kg
        out[f"dir_{theta:03d}__eps"] = eg
        headings.append(theta)
        meta[str(theta)] = dict(info, unresolved_fraction=frac)

        mag = np.linalg.norm(Ug, axis=-1)
        print(f"  |U| at 2 m: {mag[:, :, 1].min():.2f} .. "
              f"{mag[:, :, 1].max():.2f} m/s (mean {mag[:, :, 1].mean():.2f})")
        print(f"  k   at 2 m: {kg[:, :, 1].min():.3f} .. "
              f"{kg[:, :, 1].max():.3f} m2/s2")

    out["headings"] = np.array(sorted(headings))
    out["meta_json"] = json.dumps(meta)
    np.savez_compressed(args.out, **out)
    print(f"\nwrote {args.out}: {len(headings)} direction(s) "
          f"{sorted(headings)}")


if __name__ == "__main__":
    main()
