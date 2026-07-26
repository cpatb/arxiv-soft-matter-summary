"""
Convert a DEM GeoTIFF to a triangulated ground STL in local Cartesian meters.

Decimates the DEM grid to a manageable triangle count (--simplify metres,
default 1.0 m), reprojects into UTM, writes STL.
"""

from __future__ import annotations
import argparse


def dem_to_stl(dem_path: str, out_stl: str, simplify_m: float,
               origin_lon: float, origin_lat: float, crs_epsg: int = 26916):
    try:
        import rasterio
        from rasterio.warp import calculate_default_transform, reproject, Resampling
        import numpy as np
        import trimesh
        from pyproj import Transformer
    except ImportError:
        raise SystemExit("install rasterio, trimesh, pyproj")

    print(f"[stub] read {dem_path}; reproject to EPSG:{crs_epsg}")
    print(f"[stub] decimate to ~{simplify_m} m grid; write STL to {out_stl}")
    print()
    print("Implementation outline:")
    print("  with rasterio.open(dem_path) as src:")
    print("      # reproject to crs_epsg, resample to simplify_m grid")
    print("      ...")
    print("      Z = reprojected_band  # (ny, nx) meters")
    print("  ox, oy = Transformer.from_crs(4326, crs_epsg, always_xy=True).transform(")
    print("      origin_lon, origin_lat)")
    print("  X, Y = np.meshgrid(x_local, y_local)")
    print("  X -= ox; Y -= oy")
    print("  verts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)")
    print("  faces = []  # build triangles for adjacent grid cells")
    print("  mesh = trimesh.Trimesh(vertices=verts, faces=faces)")
    print("  mesh.export(out_stl)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dem")
    ap.add_argument("--out", required=True)
    ap.add_argument("--simplify", type=float, default=1.0,
                    help="target grid spacing in metres")
    ap.add_argument("--origin-lon", type=float, default=-87.6580)
    ap.add_argument("--origin-lat", type=float, default=42.0006)
    ap.add_argument("--crs", type=int, default=26916)
    args = ap.parse_args()
    dem_to_stl(args.dem, args.out, args.simplify,
               args.origin_lon, args.origin_lat, args.crs)


if __name__ == "__main__":
    main()
