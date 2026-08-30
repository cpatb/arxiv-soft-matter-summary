"""Add a porous windbreak to a solved case, as a design study.

THE PROBLEM THE BARRIER IS FOR
------------------------------
Under a westerly (`theta_wind` = 270) the disc-band turbulence over Sean Earl
Field is not uniform: `sigma_u` at z = 2 m rises monotonically west to east,
0.39 m/s at x = -30 to 0.774 m/s at x = +40, a factor **2.03**. Campus Towers'
footprint starts at x = +56 on the pitch's y-band, so this is the windward
stagnation and horseshoe-vortex system of a 45 m block spilling roughly 16 m
back onto the playing surface. It is confined to the eastern third: the middle
of the pitch is the quietest part of it.

The same map for the other two headings says this is a westerly phenomenon --
east/west `sigma_u` ratio 1.80 at 270 against 1.17 at 90 and 1.03 at 210 -- so
a barrier designed for it must be checked for what it does to the other two,
not only for what it does to the one it was designed for.

WHY POROUS AND NOT SOLID
------------------------
The agricultural shelterbelt literature is unambiguous that a *porous* barrier
gives a longer and quieter shelter than a solid one of the same height: a solid
fence produces a deeper velocity deficit but sheds a strong shear layer that
reattaches with high turbulence downstream. This project has already made the
same point in its own terms -- `FRISBEE_GENERALIZATION.md`'s headline is that
shelter is not protection, because buildings remove mean wind and ADD
turbulence, and a disc cares about the fluctuation rather than the mean. A
solid fence here would risk trading a stagnation zone for a shear layer sitting
in exactly the 1-6 m band the disc flies in.

So the barrier is modelled as a DarcyForchheimer porous cell zone, which is
the standard way shelterbelts are represented in CFD, and its resistance is a
swept parameter rather than a fixed choice.

    pressure loss over the barrier:   dp = 0.5 * rho * k_r * U^2
    DarcyForchheimer inertial term:   dp/L = 0.5 * rho * f * |U| * U
    hence                             f = k_r / t          (t = thickness)

`k_r` is the barrier's pressure-loss coefficient. For a fence of optical
porosity beta, k_r is roughly (1 - beta^2)/beta^2 for a thin screen, so
beta = 0.5 gives k_r = 3.0 and beta = 0.6 gives k_r = 1.78. Those are the
values the sweep uses; the mapping is a standard orifice-plate estimate and is
quoted as such, not measured here.

THE FRAME TRAP
--------------
Each case is solved with the SITE geometry rotated so the wind blows along +x
(`generateCases.py`). A structure fixed to the site therefore rotates with it,
and a barrier specified in site coordinates must be rotated by the case's own
`geometry_rotation_deg` before it becomes a cell zone. `boxToCell` is
axis-aligned and cannot express that, so `rotatedBoxToCell` is used throughout
-- including for the 270 case, whose rotation happens to be zero. Writing the
zero case with the general path is deliberate: an axis-aligned shortcut that
works for the one heading the barrier was designed for and silently mislocates
it for the other two is exactly the failure this study would not survive.

MESH
----
The production mesh's effective cell over the pitch is 2.0 m horizontal by
0.5 m vertical, so the barrier is given a streamwise thickness of 4 m -- two
cells -- rather than a physically thin screen. The resistance is set from the
thickness, so a thicker zone with a proportionally smaller `f` is the same
barrier; what it cannot represent is the near-field detail within a fence
thinner than the mesh. That is a stated limitation, not a hidden one.

Usage:
    python3 makeWindbreakCase.py cases/dir_270 W1 --x 21 --height 5 --kr 3.0
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path

# The playing surface, from site_view.py / the manuscript.
PITCH_X = (-31.0, 47.0)
PITCH_Y = (-10.0, 15.0)

TOPOSET = """/*--------------------------------*- C++ -*----------------------------------*/
FoamFile {{ version 2.0; format ascii; class dictionary; object topoSetDict; }}

// The barrier, expressed as a ROTATED box so the same site-frame geometry
// lands correctly in every direction's solve frame.
actions
(
    {{
        name    windbreakCells;
        type    cellSet;
        action  new;
        source  rotatedBoxToCell;
        origin  ({ox:.6f} {oy:.6f} {oz:.6f});
        i       ({ix:.6f} {iy:.6f} {iz:.6f});
        j       ({jx:.6f} {jy:.6f} {jz:.6f});
        k       ({kx:.6f} {ky:.6f} {kz:.6f});
    }}
    {{
        name    windbreak;
        type    cellZoneSet;
        action  new;
        source  setToCellZone;
        set     windbreakCells;
    }}
);
"""

FVOPTIONS = """/*--------------------------------*- C++ -*----------------------------------*/
FoamFile {{ version 2.0; format ascii; class dictionary; object fvOptions; }}

// Porous windbreak. DarcyForchheimer with the viscous term set to zero: at
// these Reynolds numbers a fence's loss is entirely form drag, and a nonzero
// Darcy term would make the barrier's resistance depend on wind speed, which
// would break the Reynolds-independence the wind-speed theory rests on.
//
//   f = k_r / t  with k_r = {kr:g} and t = {thickness:g} m
//
// The local axes are the SOLVE frame's, in which the wind blows along +x, so
// the streamwise resistance is the x component.
windbreak
{{
    type            explicitPorositySource;
    active          true;

    explicitPorositySourceCoeffs
    {{
        type            DarcyForchheimer;
        selectionMode   cellZone;
        cellZone        windbreak;

        d               (0 0 0);
        f               ({f:.6f} 0 0);

        coordinateSystem
        {{
            origin      (0 0 0);
            rotation
            {{
                type    axes;
                e1      (1 0 0);
                e2      (0 1 0);
            }}
        }}
    }}
}}
"""


def rot(v, phi):
    c, s = math.cos(phi), math.sin(phi)
    return (c * v[0] - s * v[1], s * v[0] + c * v[1], v[2])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src", help="a case directory to clone, e.g. cases/dir_270")
    ap.add_argument("tag", help="short name for the design, e.g. W1")
    ap.add_argument("--x", type=float, required=True,
                    help="barrier centre, site-frame x (m)")
    ap.add_argument("--height", type=float, required=True, help="metres")
    ap.add_argument("--kr", type=float, required=True,
                    help="pressure-loss coefficient, dp = 0.5 rho k_r U^2")
    ap.add_argument("--thickness", type=float, default=4.0,
                    help="streamwise thickness, m (2 production cells)")
    ap.add_argument("--y-span", nargs=2, type=float, default=None,
                    help="site-frame y extent; default is the pitch width "
                         "plus 10 m of overhang each side")
    ap.add_argument("--z-bottom", type=float, default=0.0)
    args = ap.parse_args()

    src = Path(args.src)
    info = json.loads((src / "case_info.json").read_text())
    phi = math.radians(info["geometry_rotation_deg"])
    dst = src.parent / f"{src.name}_{args.tag}"
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)

    for sub in ("0.orig", "constant", "system"):
        shutil.copytree(src / sub, dst / sub)
    shutil.rmtree(dst / "constant" / "polyMesh", ignore_errors=True)
    shutil.copy2(src / "case_info.json", dst / "case_info.json")

    ya, yb = args.y_span if args.y_span else (PITCH_Y[0] - 10.0,
                                              PITCH_Y[1] + 10.0)
    # Site-frame box: origin at one corner, then the three edge vectors.
    origin = (args.x - args.thickness / 2.0, ya, args.z_bottom)
    ivec = (args.thickness, 0.0, 0.0)
    jvec = (0.0, yb - ya, 0.0)

    o = rot(origin, phi)
    i = rot(ivec, phi)
    j = rot(jvec, phi)
    k = (0.0, 0.0, args.height)          # z is unaffected by the rotation

    (dst / "system" / "topoSetDict").write_text(TOPOSET.format(
        ox=o[0], oy=o[1], oz=o[2], ix=i[0], iy=i[1], iz=i[2],
        jx=j[0], jy=j[1], jz=j[2], kx=k[0], ky=k[1], kz=k[2]))

    f = args.kr / args.thickness
    (dst / "constant" / "fvOptions").write_text(FVOPTIONS.format(
        kr=args.kr, thickness=args.thickness, f=f))

    info["windbreak"] = dict(tag=args.tag, x_site_m=args.x,
                             height_m=args.height, thickness_m=args.thickness,
                             k_r=args.kr, forchheimer_f=f,
                             y_span_site_m=[ya, yb],
                             z_bottom_m=args.z_bottom,
                             geometry_rotation_deg=info["geometry_rotation_deg"],
                             model="DarcyForchheimer porous cellZone, "
                                   "viscous term zero",
                             porosity_note="k_r ~ (1-beta^2)/beta^2 for a thin "
                                           "screen: beta 0.5 -> 3.0, "
                                           "0.6 -> 1.78")
    (dst / "case_info.json").write_text(json.dumps(info, indent=2))

    print(f"wrote {dst}")
    print(f"  barrier at site x = {args.x:+.1f} m, height {args.height:.1f} m, "
          f"thickness {args.thickness:.1f} m, k_r = {args.kr:g}, "
          f"f = {f:.4f} 1/m")
    print(f"  y span {ya:+.1f} to {yb:+.1f} m, geometry rotation "
          f"{info['geometry_rotation_deg']:+.1f} deg")
    print(f"  solve-frame box origin ({o[0]:+.2f} {o[1]:+.2f} {o[2]:+.2f}), "
          f"i ({i[0]:+.2f} {i[1]:+.2f} {i[2]:+.2f}), "
          f"j ({j[0]:+.2f} {j[1]:+.2f} {j[2]:+.2f})")


if __name__ == "__main__":
    main()
