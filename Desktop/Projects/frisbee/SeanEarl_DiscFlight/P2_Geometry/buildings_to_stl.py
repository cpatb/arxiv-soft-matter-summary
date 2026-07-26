"""
Extrude OSM building footprints to a watertight, multi-region STL for
snappyHexMesh.

Input is the JSON written by `fetch_buildings.py` (polygons already in local
metres, x East / y North, origin at the Sean Earl Field centroid).

GROUND MODEL
------------
Buildings are extruded from a FLAT ground plane at z = 0, which is the field
elevation (~179.5 m MSL). This is justified by USGS 3DEP point queries over a
600 m box around the field:

    campus plateau   179.1 - 181.4 m     (+-1 m across the field itself)
    lake / beach     175.2 - 175.6 m     (east of about x = +150 m)

so relief over the playing surface and the near-field buildings is about +-1 m,
an order of magnitude below the building heights that set the flow. The one
real topographic feature is a ~4 m step down to the lake at the shoreline; that
is NOT represented as topography. It is instead represented where it matters,
as the inlet roughness transition (z0 = 0.001 m over water vs 0.5 m urban) that
generates the internal boundary layer the paper is about. Neglecting the step
itself displaces the IBL origin by a few metres of vertical offset over 150+ m
of fetch; the roughness change dominates. This is stated in the paper.

TRIANGULATION
-------------
Footprints are not necessarily convex, so caps are triangulated with earcut
rather than a naive fan, which would produce self-intersecting caps on L- and
U-shaped buildings and leak in snappyHexMesh.

Output is an ASCII multi-solid STL: one `solid` block per building, named
`bldg_<osm_id>`, so snappyHexMesh can refine named buildings individually.
"""

from __future__ import annotations

import argparse
import json

import numpy as np


def _signed_area(ring: np.ndarray) -> float:
    x, y = ring[:, 0], ring[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def triangulate(ring: np.ndarray) -> np.ndarray:
    """Triangulate a simple polygon; returns (n_tri, 3) vertex indices."""
    import mapbox_earcut as earcut
    idx = earcut.triangulate_float64(ring, np.array([len(ring)]))
    return np.asarray(idx, dtype=np.int64).reshape(-1, 3)


def extrude(ring: np.ndarray, z0: float, z1: float):
    """Return (vertices, faces) for a closed extruded prism.

    The ring is forced counter-clockwise so that the outward normal
    convention is consistent: side faces wind outward, the roof cap points up
    and the floor cap points down. snappyHexMesh does not strictly require
    consistent normals, but an inconsistent surface makes `surfaceCheck`
    report the region as non-orientable and makes inside/outside tests for
    refinement regions unreliable.
    """
    if _signed_area(ring) < 0:
        ring = ring[::-1]
    n = len(ring)

    bottom = np.column_stack([ring, np.full(n, z0)])
    top = np.column_stack([ring, np.full(n, z1)])
    verts = np.vstack([bottom, top])

    faces = []
    # Side walls: quad (i, i+1, i+1+n, i+n) -> two triangles, outward normals.
    for i in range(n):
        j = (i + 1) % n
        faces.append([i, j, j + n])
        faces.append([i, j + n, i + n])

    cap = triangulate(ring)
    for tri in cap:                       # roof, +z outward
        faces.append([tri[0] + n, tri[1] + n, tri[2] + n])
    for tri in cap:                       # floor, -z outward (reversed)
        faces.append([tri[2], tri[1], tri[0]])

    return verts, np.array(faces, dtype=np.int64)


def write_multisolid_stl(path, solids):
    """Write an ASCII STL with one named solid per entry.

    solids: list of (name, vertices (nv,3), faces (nf,3))
    """
    with open(path, "w") as fh:
        for name, verts, faces in solids:
            fh.write(f"solid {name}\n")
            v = verts[faces]                                # (nf, 3, 3)
            n = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
            ln = np.linalg.norm(n, axis=1, keepdims=True)
            n = np.divide(n, ln, out=np.zeros_like(n), where=ln > 0)
            for tri, nrm in zip(v, n):
                fh.write(f"  facet normal {nrm[0]:.6e} {nrm[1]:.6e} "
                         f"{nrm[2]:.6e}\n    outer loop\n")
                for p in tri:
                    fh.write(f"      vertex {p[0]:.6e} {p[1]:.6e} "
                             f"{p[2]:.6e}\n")
                fh.write("    endloop\n  endfacet\n")
            fh.write(f"endsolid {name}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("geometry_json")
    ap.add_argument("--out", default="buildings.stl")
    ap.add_argument("--half-size", type=float, default=None,
                    help="clip to a square half-width (m) about the origin")
    ap.add_argument("--min-height", type=float, default=2.0)
    ap.add_argument("--base-z", type=float, default=0.0,
                    help="flat ground elevation in local coords")
    args = ap.parse_args()

    with open(args.geometry_json) as fh:
        G = json.load(fh)

    solids = []
    skipped_far = skipped_small = 0
    degenerate = 0
    for b in G["buildings"]:
        ring = np.asarray(b["ring"], dtype=float)
        if len(ring) < 3:
            degenerate += 1
            continue

        cx, cy = ring[:, 0].mean(), ring[:, 1].mean()
        if args.half_size is not None and (abs(cx) > args.half_size
                                           or abs(cy) > args.half_size):
            skipped_far += 1
            continue

        h = float(b["height_m"])
        if h < args.min_height:
            skipped_small += 1
            continue

        # Drop consecutive duplicate vertices, which earcut cannot handle.
        keep = np.ones(len(ring), dtype=bool)
        d = np.linalg.norm(ring - np.roll(ring, -1, axis=0), axis=1)
        keep[d < 1e-6] = False
        ring = ring[keep]
        if len(ring) < 3 or abs(_signed_area(ring)) < 1.0:
            degenerate += 1
            continue

        try:
            verts, faces = extrude(ring, args.base_z, args.base_z + h)
        except Exception as exc:                       # noqa: BLE001
            print(f"  skipping {b.get('osm_id')}: {exc}")
            degenerate += 1
            continue
        solids.append((f"bldg_{b['osm_id']}", verts, faces))

    write_multisolid_stl(args.out, solids)

    n_tri = sum(len(f) for _, _, f in solids)
    print(f"wrote {args.out}")
    print(f"  {len(solids)} buildings, {n_tri} triangles")
    print(f"  skipped: {skipped_far} outside clip, {skipped_small} below "
          f"{args.min_height} m, {degenerate} degenerate")
    if solids:
        allv = np.vstack([v for _, v, _ in solids])
        print(f"  bounds x [{allv[:,0].min():.1f}, {allv[:,0].max():.1f}] "
              f"y [{allv[:,1].min():.1f}, {allv[:,1].max():.1f}] "
              f"z [{allv[:,2].min():.1f}, {allv[:,2].max():.1f}] m")


if __name__ == "__main__":
    main()
