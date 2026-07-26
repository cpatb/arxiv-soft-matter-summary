"""
Post-process: extract U_mean, k, epsilon (omega) on a 0.5 m horizontal
grid at 1, 2, 3 m AGL from each completed OpenFOAM case.

Writes everything into a single .npz with keys per direction:
    fields['dir_NNN']  -> dict(U=(nx,ny,nz,3), k=(nx,ny,nz), eps=(nx,ny,nz))
    fields['x']        -> (nx,)  local field x in meters
    fields['y']        -> (ny,)
    fields['z']        -> (nz,)  AGL heights, e.g. [1.0, 2.0, 3.0]
    fields['headings'] -> sorted list of integer headings

Uses PyFoam (or foamToVTK + pyvista as fallback).
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np


def sample_one_case(case_dir: Path, x_grid, y_grid, z_grid):
    """Return dict of U, k, eps interpolated onto the grid for one case."""
    try:
        import pyvista as pv
    except ImportError:
        raise SystemExit("install pyvista, then `foamToVTK` on each case")

    # Run foamToVTK first if VTK directory missing
    vtk_dir = case_dir / "VTK"
    if not vtk_dir.exists():
        raise FileNotFoundError(
            f"{vtk_dir} missing -- run `foamToVTK` in {case_dir} first"
        )

    # Latest time step VTK
    vtk_files = sorted(vtk_dir.glob(f"{case_dir.name}_*.vtm"))
    if not vtk_files:
        raise FileNotFoundError(f"no .vtm files in {vtk_dir}")
    mb = pv.read(vtk_files[-1])

    # Build a structured probe grid
    X, Y, Z = np.meshgrid(x_grid, y_grid, z_grid, indexing="ij")
    probe = pv.PolyData(np.column_stack([X.ravel(), Y.ravel(), Z.ravel()]))
    sampled = probe.sample(mb)

    U   = sampled["U"].reshape(X.shape + (3,))
    k   = sampled["k"].reshape(X.shape)
    # OpenFOAM kOmegaSST stores omega; epsilon = beta* k omega with beta*=0.09
    om  = sampled["omega"].reshape(X.shape)
    eps = 0.09 * k * om

    return dict(U=U, k=k, eps=eps)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="+", help="case directories")
    ap.add_argument("--out", default="flow_fields.npz")
    ap.add_argument("--xrange", nargs=2, type=float, default=[-50, 50])
    ap.add_argument("--yrange", nargs=2, type=float, default=[-50, 50])
    ap.add_argument("--dx", type=float, default=0.5)
    ap.add_argument("--z-agl", nargs="+", type=float, default=[1.0, 2.0, 3.0])
    args = ap.parse_args()

    x = np.arange(args.xrange[0], args.xrange[1] + args.dx / 2, args.dx)
    y = np.arange(args.yrange[0], args.yrange[1] + args.dx / 2, args.dx)
    z = np.array(args.z_agl)

    out: dict = {"x": x, "y": y, "z": z}
    headings = []
    for case_str in args.cases:
        case = Path(case_str)
        name = case.name                # e.g. dir_045
        heading = int(name.split("_")[1])
        headings.append(heading)
        print(f"sampling {case} (heading {heading})...")
        fields = sample_one_case(case, x, y, z)
        for k, v in fields.items():
            out[f"{name}__{k}"] = v
    out["headings"] = np.array(sorted(headings))

    np.savez_compressed(args.out, **out)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
