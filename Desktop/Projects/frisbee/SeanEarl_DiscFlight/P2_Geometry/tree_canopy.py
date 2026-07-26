"""
Hand-traced tree-canopy polygons over the Sean Earl Field lakefront strip,
extruded to a porous-zone STL.

Coordinates below are placeholders in local Cartesian (origin = field
center). Replace with QGIS-traced polygons from a recent aerial image
(NAIP 2024 or Cook County orthophotos).

The OpenFOAM porous-zone treatment (P3) uses leaf area density
alpha_LAD ~ 0.5 m^2/m^3 (Gromke & Ruck 2008) inside the volume swept by
this STL.
"""

from __future__ import annotations
import argparse


# Placeholder polygons in local meters, origin = Sean Earl Field center.
# +x East, +y North.
TREE_STRIPS = [
    {
        "name": "lakefront_north",
        "polygon_xy": [(50, 80), (180, 80), (180, 110), (50, 110)],
        "canopy_base_m": 4.0,
        "canopy_top_m":  12.0,
    },
    {
        "name": "lakefront_south",
        "polygon_xy": [(50, -110), (180, -110), (180, -80), (50, -80)],
        "canopy_base_m": 4.0,
        "canopy_top_m":  10.0,
    },
]


def build_canopy_stl(out_stl: str):
    try:
        import numpy as np
        import trimesh
    except ImportError:
        raise SystemExit("install trimesh")

    print(f"[stub] build canopy STL with {len(TREE_STRIPS)} strips -> {out_stl}")
    print()
    print("Implementation outline:")
    print("  meshes = []")
    print("  for strip in TREE_STRIPS:")
    print("      base = strip['canopy_base_m']")
    print("      top  = strip['canopy_top_m']")
    print("      poly2d = shapely.Polygon(strip['polygon_xy'])")
    print("      m = trimesh.creation.extrude_polygon(poly2d, height=top - base)")
    print("      m.apply_translation([0, 0, base])")
    print("      meshes.append(m)")
    print("  trimesh.util.concatenate(meshes).export(out_stl)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    build_canopy_stl(args.out)


if __name__ == "__main__":
    main()
