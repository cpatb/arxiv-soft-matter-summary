"""
Compose the final `seanEarl.stl` for snappyHexMesh and report its quality.

WHAT IS AND IS NOT IN THE STL
-----------------------------
Buildings only.

  * GROUND is NOT an STL surface. USGS 3DEP point queries give +-1 m of relief
    across the field and the near-field buildings (campus plateau 179.1-181.4 m
    MSL), so the ground is represented as the flat bottom patch of the
    blockMesh background domain with a rough-wall function. That is both more
    robust than snapping to a near-planar triangulation and avoids the sliver
    cells a decimated DEM produces. The one real topographic feature, a ~4 m
    step down to the lake at the shoreline (east of about x = +150 m), is
    represented where it actually matters -- as the inlet roughness transition
    that generates the internal boundary layer -- rather than as topography.

  * TREE CANOPY is NOT included. The manuscript proposes a lakefront canopy
    porous zone with leaf area density ~0.5 m^2/m^3 (Gromke & Ruck 2008). OSM
    maps 453 individual trees in the box but no contiguous canopy polygon east
    of the field; the two mapped woodland polygons lie ~200 m SOUTHWEST, and
    the lakefront parks are to the NORTHEAST. Rather than hand-trace a canopy
    we could not validate, the porous zone is omitted and recorded as scoped
    down. This is defensible here because the dominant obstruction between the
    field and the lake is not vegetation but Campus Towers, a 45 m apartment
    block whose footprint starts 9 m from the field's eastern edge.

Both omissions are stated in the manuscript.
"""

from __future__ import annotations

import argparse
import shutil

import numpy as np


def report(stl_path: str):
    """Load the STL and report the checks snappyHexMesh cares about."""
    import trimesh

    scene = trimesh.load(stl_path, process=False)
    if hasattr(scene, "geometry"):
        parts = dict(scene.geometry)
    else:
        parts = {"all": scene}

    n_solids = len(parts)
    n_tri = sum(len(m.faces) for m in parts.values())
    watertight = sum(1 for m in parts.values() if m.is_watertight)
    inverted = sum(1 for m in parts.values() if m.volume < 0)

    allv = np.vstack([m.vertices for m in parts.values()])
    heights = np.array([m.vertices[:, 2].max() for m in parts.values()])

    print(f"  solids      : {n_solids}")
    print(f"  triangles   : {n_tri}")
    print(f"  watertight  : {watertight}/{n_solids}"
          f"{'  ** LEAKS **' if watertight != n_solids else ''}")
    print(f"  inverted    : {inverted}"
          f"{'  ** NORMALS **' if inverted else ''}")
    print(f"  bounds x    : [{allv[:,0].min():8.1f}, {allv[:,0].max():8.1f}] m")
    print(f"  bounds y    : [{allv[:,1].min():8.1f}, {allv[:,1].max():8.1f}] m")
    print(f"  bounds z    : [{allv[:,2].min():8.1f}, {allv[:,2].max():8.1f}] m")
    print(f"  tallest H   : {heights.max():.1f} m")

    # Smallest triangle edge drives the snappyHexMesh refinement level needed
    # to resolve the surface without producing slivers.
    edges = []
    for m in parts.values():
        v = m.vertices[m.faces]
        for a, b in ((0, 1), (1, 2), (2, 0)):
            edges.append(np.linalg.norm(v[:, a] - v[:, b], axis=1))
    edges = np.concatenate(edges)
    edges = edges[edges > 0]
    print(f"  edge length : min {edges.min():.3f} m, "
          f"median {np.median(edges):.2f} m")
    if edges.min() < 0.1:
        print("  ** sub-0.1 m edges present; snappyHexMesh may need a high "
              "refinement level or these slivers will be lost **")
    return heights.max()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--buildings", default="buildings.stl")
    ap.add_argument("--trees", default=None,
                    help="optional canopy STL; omitted by default (see module "
                         "docstring)")
    ap.add_argument("--ground", default=None,
                    help="optional ground STL; omitted by default, the ground "
                         "is the blockMesh bottom patch")
    ap.add_argument("--out", default="seanEarl.stl")
    args = ap.parse_args()

    if args.ground or args.trees:
        # Multi-solid ASCII STL concatenation preserves the solid names that
        # snappyHexMesh uses as region handles, so a plain concatenation is
        # the correct merge here.
        with open(args.out, "w") as out:
            for src in filter(None, (args.buildings, args.ground, args.trees)):
                with open(src) as fh:
                    shutil.copyfileobj(fh, out)
    else:
        shutil.copyfile(args.buildings, args.out)

    print(f"wrote {args.out}")
    H = report(args.out)

    print("\ndomain guidance (COST 732 / Tominaga 2008) for H = "
          f"{H:.0f} m:")
    print(f"  >= 5H  above tallest building : {5*H:6.0f} m")
    print(f"  >= 5H  upstream and lateral   : {5*H:6.0f} m")
    print(f"  >= 15H downstream             : {15*H:6.0f} m")
    print("  (P3 records which of these are met and which are relaxed for "
          "compute)")


if __name__ == "__main__":
    main()
