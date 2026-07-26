"""
Extrude a building-polygon GeoPackage to a triangulated STL.

Each polygon in the input is extruded vertically from its base height
(read from the DEM at the polygon centroid, or assumed 0) to base+height_m.
Resulting solid is triangulated and merged into one STL.

Inputs:
    --polygons  : GeoPackage with 'geometry' and 'height_m' columns
    --dem       : (optional) GeoTIFF DEM to sample base elevation
    --crs       : EPSG code for local Cartesian projection (default 26916 = UTM 16N)

Output:
    --out       : STL written in local Cartesian meters, origin = field center
"""

from __future__ import annotations
import argparse


def extrude_polygons(polygons_path: str, dem_path: str | None,
                     crs_epsg: int, out_stl: str,
                     origin_lon: float, origin_lat: float):
    try:
        import geopandas as gpd
        import numpy as np
        import shapely.geometry as sg
        import trimesh
        from pyproj import Transformer
    except ImportError:
        raise SystemExit("install geopandas, shapely, trimesh, pyproj")

    print(f"[stub] read {polygons_path}, extrude to STL at {out_stl}")
    print(f"[stub] origin (lon, lat) = ({origin_lon}, {origin_lat}), CRS EPSG:{crs_epsg}")
    print()
    print("Implementation outline:")
    print("  gdf = gpd.read_file(polygons_path).to_crs(epsg=crs_epsg)")
    print("  T = Transformer.from_crs(4326, crs_epsg, always_xy=True)")
    print("  ox, oy = T.transform(origin_lon, origin_lat)")
    print("  meshes = []")
    print("  for _, row in gdf.iterrows():")
    print("      poly = row.geometry  # in EPSG:crs_epsg")
    print("      coords = np.array(poly.exterior.coords)")
    print("      coords[:,0] -= ox; coords[:,1] -= oy   # to local origin")
    print("      base_z = sample_dem(dem_path, poly.centroid)  # or 0")
    print("      h = row['height_m']")
    print("      m = trimesh.creation.extrude_polygon(")
    print("          shapely.Polygon(coords[:,:2]), height=h)")
    print("      m.apply_translation([0, 0, base_z])")
    print("      meshes.append(m)")
    print("  combined = trimesh.util.concatenate(meshes)")
    print("  combined.export(out_stl)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("polygons")
    ap.add_argument("--dem", default=None)
    ap.add_argument("--crs", type=int, default=26916, help="EPSG (UTM 16N default)")
    ap.add_argument("--origin-lon", type=float, default=-87.6580)
    ap.add_argument("--origin-lat", type=float, default=42.0006)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    extrude_polygons(args.polygons, args.dem, args.crs, args.out,
                     args.origin_lon, args.origin_lat)


if __name__ == "__main__":
    main()
