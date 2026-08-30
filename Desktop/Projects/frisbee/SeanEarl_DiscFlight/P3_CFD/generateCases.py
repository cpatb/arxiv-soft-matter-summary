"""
Generate one OpenFOAM case per inflow direction.

DIRECTION CONVENTION
--------------------
`theta_wind` is METEOROLOGICAL: the direction the wind blows FROM, in degrees
clockwise from north. theta = 90 is an easterly, i.e. wind off the lake. This
matches the manuscript ("easterly over-lake ~90 deg").

HOW DIRECTIONS ARE REALISED
---------------------------
Rather than move the inlet from patch to patch (which only works for the four
axis-aligned directions and gives a different mesh for each), the SITE GEOMETRY
IS ROTATED so the wind always blows along +x. The blockMesh, the boundary
conditions and the inlet patch are then identical for every direction, and each
direction differs only in the rotated STL. Mesh-induced differences between
directions are therefore minimised, which matters because the paper compares
directions against each other.

The building set is clipped to a CIRCLE (r <= 300 m) in P2 precisely so it is
invariant under this rotation: every direction sees the same buildings, and
none rotate in or out of the modelled set.

Rotation angle: the wind blows toward d = (-sin(theta), -cos(theta)) in
(East, North). We need d to map onto +x, so the geometry is rotated by
phi = -atan2(d_y, d_x).

INLET ROUGHNESS
---------------
Per the manuscript, z0 = 0.001 m (open water) for theta in [60, 120] and
z0 = 0.5 m (urban) otherwise, blended linearly over 30 deg transition bands so
neighbouring directions in the sweep are not discontinuous.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path

# --- Site / domain constants ------------------------------------------------
# H is the tallest building inside the modelled r <= 300 m set (Mertz Hall).
H_TALLEST = 60.8
BUILDING_RADIUS = 300.0

# Ground roughness BETWEEN the resolved buildings: grass and pavement.
# Deliberately NOT the 0.5 m urban value, which is an equivalent roughness
# standing in for buildings that we resolve explicitly.
Z0_GROUND = 0.03

Z0_WATER = 0.001
Z0_URBAN = 0.5
KAPPA = 0.41
CMU = 0.09


def inlet_z0(theta_deg: float) -> float:
    """Blended inlet roughness: open water to the east, urban elsewhere."""
    t = theta_deg % 360.0
    if 60.0 <= t <= 120.0:
        return Z0_WATER
    if 30.0 <= t < 60.0:
        f = (t - 30.0) / 30.0
        return Z0_URBAN * (1 - f) + Z0_WATER * f
    if 120.0 < t <= 150.0:
        f = (t - 120.0) / 30.0
        return Z0_WATER * (1 - f) + Z0_URBAN * f
    return Z0_URBAN


def rotation_for(theta_deg: float) -> float:
    """Geometry rotation (radians, CCW) that puts the wind along +x."""
    th = math.radians(theta_deg)
    dx, dy = -math.sin(th), -math.cos(th)
    return -math.atan2(dy, dx)


def rotate_stl(src: Path, dst: Path, phi: float):
    """Rotate an ASCII STL about the z axis by phi radians."""
    c, s = math.cos(phi), math.sin(phi)
    out = []
    with open(src) as fh:
        for line in fh:
            t = line.split()
            if t and t[0] == "vertex":
                x, y, z = float(t[1]), float(t[2]), float(t[3])
                out.append(f"      vertex {c*x - s*y:.6e} "
                           f"{s*x + c*y:.6e} {z:.6e}\n")
            elif t and t[0] == "facet":
                x, y, z = float(t[2]), float(t[3]), float(t[4])
                out.append(f"  facet normal {c*x - s*y:.6e} "
                           f"{s*x + c*y:.6e} {z:.6e}\n")
            else:
                out.append(line)
    dst.write_text("".join(out))


def abl_reference(u_ref, z_ref, z0):
    """u*, and the k / omega the Richards-Hoxey profile implies at z_ref."""
    u_star = u_ref * KAPPA / math.log((z_ref + z0) / z0)
    k = u_star ** 2 / math.sqrt(CMU)
    omega = u_star / (KAPPA * math.sqrt(CMU) * (z_ref + z0))
    return u_star, k, omega


def solve_grading(zmax, dz1):
    """Pick (nz, ratio) for a geometric vertical distribution."""
    for n in range(20, 121):
        lo, hi = 1.0001, 1.30
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            tot = dz1 * (mid ** n - 1.0) / (mid - 1.0)
            if tot < zmax:
                lo = mid
            else:
                hi = mid
        r = 0.5 * (lo + hi)
        if 1.02 <= r <= 1.12:
            return n, r
    return 45, 1.06


# Resolution ladder. "coarse" and "production" differ ONLY in the fieldBox
# refinement level, so comparing them is a clean factor-2 grid-independence
# test in exactly the region the paper reports, with everything else fixed.
#
# Do not raise `surf` and `feat` together with a large `dxy` reduction: that
# combination (dxy=12, surf=(2,3), feat=3) makes hexRef8 fail with
# "nAnchors:2 facei:..." on this geometry, a refinement-consistency error from
# too abrupt a level jump between the graded background cells and the feature
# refinement. The ladder below keeps the surface/feature levels fixed and moves
# only the volume refinement.
RESOLUTIONS = {
    "coarse":     dict(dxy=16.0, dz1=4.0, surf=(1, 2), near=1, field=2, feat=2,
                       maxcells=4_000_000),
    "production": dict(dxy=16.0, dz1=4.0, surf=(1, 2), near=1, field=3, feat=2,
                       maxcells=8_000_000),
    "fine":       dict(dxy=12.0, dz1=3.0, surf=(1, 2), near=2, field=3, feat=2,
                       maxcells=16_000_000),
    # --- Grid-convergence family, added 2026-08-29 -----------------------------
    # `coarse` and `production` are a clean factor-2 pair, but TWO grids give a
    # DIFFERENCE, not an ERROR: they cannot yield an observed order or a
    # Richardson extrapolation. These two extend that pair into a four-member
    # family that moves ONLY the fieldBox level, exactly as the note above
    # requires, so the meshes stay geometrically similar and GCI is defined:
    #
    #     gridL1  field=1   8.0 x 2.0 m      gridL4  field=4   1.0 x 0.25 m
    #     coarse  field=2   4.0 x 1.0 m      production field=3  2.0 x 0.50 m
    #
    # `fine` is NOT a member: it moves dxy and near as well, so it is a
    # different mesh family and must not be used in a Richardson fit.
    "gridL1":     dict(dxy=16.0, dz1=4.0, surf=(1, 2), near=1, field=1, feat=2,
                       maxcells=4_000_000),
    "gridL4":     dict(dxy=16.0, dz1=4.0, surf=(1, 2), near=1, field=4, feat=2,
                       maxcells=16_000_000),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stl", default="../P2_Geometry/seanEarl.stl")
    ap.add_argument("--out", default="cases")
    ap.add_argument("--headings", nargs="+", type=int, required=True)
    ap.add_argument("--u-ref", type=float, default=5.0)
    ap.add_argument("--z-ref", type=float, default=10.0)
    ap.add_argument("--nproc", type=int, default=20)
    ap.add_argument("--end-time", type=int, default=2500)
    ap.add_argument("--resolution", choices=list(RESOLUTIONS),
                    default="medium")
    ap.add_argument("--tag", default="", help="suffix for case dir names")
    ap.add_argument("--z0-ground", type=float, default=None,
                    help="override the ground roughness. Set it equal to the "
                         "inlet z0 for a Hargreaves & Wright homogeneity test: "
                         "with mismatched roughness the profile is SUPPOSED to "
                         "evolve, so the test would measure the intended "
                         "roughness transition rather than numerical error.")
    ap.add_argument("--empty", action="store_true",
                    help="omit the buildings entirely: an empty-fetch case for "
                         "the Hargreaves & Wright (2007) ABL-recovery test")
    args = ap.parse_args()

    stl = Path(args.stl).resolve()
    if not stl.exists():
        raise SystemExit(f"missing STL: {stl}  (run P2 first)")

    template = Path(__file__).parent / "caseTemplate"
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    pad = 5.0 * H_TALLEST
    xmin = -(BUILDING_RADIUS + pad)
    xmax = BUILDING_RADIUS + 2.0 * pad
    ymin = -(BUILDING_RADIUS + pad)
    ymax = BUILDING_RADIUS + pad
    zmin, zmax = 0.0, 5.0 * H_TALLEST

    res = RESOLUTIONS[args.resolution]
    nx = int(round((xmax - xmin) / res["dxy"]))
    ny = int(round((ymax - ymin) / res["dxy"]))
    nz, ratio = solve_grading(zmax, res["dz1"])
    z_last_over_first = ratio ** (nz - 1)

    print(f"domain  x [{xmin:.0f}, {xmax:.0f}]  y [{ymin:.0f}, {ymax:.0f}]  "
          f"z [0, {zmax:.0f}] m")
    print(f"background {nx} x {ny} x {nz} = {nx*ny*nz:,} cells "
          f"({res['dxy']} m horizontal, {res['dz1']} m first vertical cell, "
          f"grading {z_last_over_first:.1f})")
    print(f"resolution '{args.resolution}': field box level {res['field']} -> "
          f"{res['dxy']/2**res['field']:.2f} m horizontal, "
          f"{res['dz1']/2**res['field']:.2f} m vertical near the ground")

    for theta in args.headings:
        name = f"dir_{theta:03d}{args.tag}"
        case = outdir / name
        if case.exists():
            shutil.rmtree(case)
        shutil.copytree(template, case)
        (case / "constant" / "triSurface").mkdir(parents=True, exist_ok=True)

        phi = rotation_for(theta)
        rotate_stl(stl, case / "constant" / "triSurface" / "seanEarl.stl", phi)

        z0i = inlet_z0(theta)
        z0g = Z0_GROUND if args.z0_ground is None else args.z0_ground
        us, kk, om = abl_reference(args.u_ref, args.z_ref, z0i)

        subs = {
            "@@XMIN@@": f"{xmin:.1f}", "@@XMAX@@": f"{xmax:.1f}",
            "@@YMIN@@": f"{ymin:.1f}", "@@YMAX@@": f"{ymax:.1f}",
            "@@ZMIN@@": f"{zmin:.1f}", "@@ZMAX@@": f"{zmax:.1f}",
            "@@NX@@": str(nx), "@@NY@@": str(ny), "@@NZ@@": str(nz),
            "@@ZRATIO@@": f"{z_last_over_first:.4f}",
            "@@ZFIRST@@": f"{res['dz1']:.2f}",
            "@@UREF@@": f"{args.u_ref:g}", "@@ZREF@@": f"{args.z_ref:g}",
            "@@Z0_INLET@@": f"{z0i:g}", "@@Z0_GROUND@@": f"{z0g:g}",
            "@@KINIT@@": f"{kk:.5g}", "@@OMEGAINIT@@": f"{om:.5g}",
            "@@ENDTIME@@": str(args.end_time),
            "@@WRITEINTERVAL@@": str(max(args.end_time // 2, 100)),
            "@@NPROC@@": str(args.nproc),
            "@@MAXCELLS@@": str(res["maxcells"]),
            "@@SURFMIN@@": str(res["surf"][0]),
            "@@SURFMAX@@": str(res["surf"][1]),
            "@@NEARLEVEL@@": str(res["near"]),
            "@@FIELDLEVEL@@": str(res["field"]),
            "@@FEATLEVEL@@": str(res["feat"]),
            "@@FBXMIN@@": "-90", "@@FBXMAX@@": "110",
            "@@FBYMIN@@": "-70", "@@FBYMAX@@": "70",
            "@@FBZMAX@@": "20",
            "@@NBXMIN@@": "-260", "@@NBXMAX@@": "260",
            "@@NBYMIN@@": "-260", "@@NBYMAX@@": "260",
            "@@NBZMAX@@": "90",
            "@@LOCX@@": "0", "@@LOCY@@": "0", "@@LOCZ@@": "120",
        }

        for f in list((case / "system").iterdir()) + \
                 list((case / "0.orig").iterdir()):
            if f.is_file():
                txt = f.read_text()
                for kx, v in subs.items():
                    txt = txt.replace(kx, v)
                f.write_text(txt)

        if args.empty:
            # Empty-fetch validation case: no buildings, so snappyHexMesh is
            # skipped entirely and the ABL profile must survive advection
            # across the bare domain (Hargreaves & Wright 2007).
            (case / "EMPTY_FETCH").write_text(
                "no buildings; run blockMesh + simpleFoam only\n")
            for fld in ("U", "k", "omega", "nut", "p"):
                p = case / "0.orig" / fld
                p.write_text(p.read_text().replace(
                    '"(ground|buildings)"', "ground").replace(
                    "    buildings\n    {\n", "    _unused_buildings\n    {\n"))

        info = dict(
            theta_wind_deg=theta,
            convention="meteorological: direction the wind blows FROM, "
                       "clockwise from north",
            geometry_rotation_deg=math.degrees(phi),
            empty_fetch=bool(args.empty),
            u_ref=args.u_ref, z_ref=args.z_ref,
            z0_inlet=z0i, z0_ground=z0g,
            u_star=us, k_inlet_ref=kk, omega_inlet_ref=om,
            domain=dict(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                        zmin=zmin, zmax=zmax),
            background_cells=[nx, ny, nz],
            background_dxy=res["dxy"], first_dz=res["dz1"],
            z_grading_last_over_first=z_last_over_first,
            resolution=args.resolution,
            refinement=dict(surface=list(res["surf"]),
                            nearFieldBox=res["near"], fieldBox=res["field"]),
            effective_field_cell_m=[res["dxy"] / 2 ** res["field"],
                                    res["dz1"] / 2 ** res["field"]],
            building_set_radius_m=BUILDING_RADIUS,
            H_tallest_m=H_TALLEST,
            cost732=dict(
                upstream_H=(abs(xmin) - BUILDING_RADIUS) / H_TALLEST,
                lateral_H=(abs(ymin) - BUILDING_RADIUS) / H_TALLEST,
                top_H=zmax / H_TALLEST,
                downstream_H=(xmax - BUILDING_RADIUS) / H_TALLEST,
                downstream_guideline_H=15.0,
                note="upstream/lateral/top meet the 5H guideline; downstream "
                     "is relaxed from 15H for compute and recorded here",
            ),
        )
        (case / "case_info.json").write_text(json.dumps(info, indent=2))
        print(f"  {name}: rotate {math.degrees(phi):+7.1f} deg, "
              f"z0_inlet {z0i:.4f} m, u* {us:.3f} m/s")

    print(f"\nwrote {len(args.headings)} case(s) under {outdir}/")


if __name__ == "__main__":
    main()
