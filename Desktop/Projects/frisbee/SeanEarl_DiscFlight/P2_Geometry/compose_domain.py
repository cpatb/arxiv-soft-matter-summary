"""
Concatenate ground, buildings, and trees into a single STL for snappyHexMesh.

Each input STL is loaded; named regions are added to the final STL so that
snappyHexMesh can target them individually in `geometry` blocks.
"""

from __future__ import annotations
import argparse


def compose(buildings_stl, ground_stl, trees_stl, out_stl):
    try:
        import trimesh
    except ImportError:
        raise SystemExit("install trimesh")

    parts = {}
    if ground_stl:
        parts["ground"] = trimesh.load(ground_stl)
    if buildings_stl:
        parts["buildings"] = trimesh.load(buildings_stl)
    if trees_stl:
        parts["trees"] = trimesh.load(trees_stl)

    print(f"composing: {list(parts.keys())} -> {out_stl}")
    # trimesh's Scene-export-to-STL writes a single mesh with named regions
    scene = trimesh.Scene()
    for name, mesh in parts.items():
        scene.add_geometry(mesh, node_name=name)
    scene.export(out_stl)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--buildings", default=None)
    ap.add_argument("--ground", default=None)
    ap.add_argument("--trees", default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    compose(args.buildings, args.ground, args.trees, args.out)


if __name__ == "__main__":
    main()
