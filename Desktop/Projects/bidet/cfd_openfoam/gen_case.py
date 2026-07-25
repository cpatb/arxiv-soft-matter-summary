#!/usr/bin/env python3
"""
gen_case.py -- generate an axisymmetric-wedge OpenFOAM v2606 case for a round
water jet (nozzle diameter d) impinging normally on a flat plate at standoff H,
solved with simpleFoam (steady, incompressible) + k-omega SST.

Geometry (wedge about the x-axis; x = axial jet direction, r = radial):
    x = 0    : nozzle exit plane.  r in [0, d/2] = inlet (uniform V0 in +x).
                                   r in [d/2, R] = entrainment (top boundary).
    x = H    : the plate (no-slip wall) -> tau_w(r) sampled here.
    r = R    : outer radial boundary = outlet.
    r = 0    : axis (collapsed wedge apex -- shared vertices, no patch/faces).
    +/- wedge faces : wedge patches.

Assumption (stated in the paper): single liquid phase only; the free surface of
the spreading film is not modelled. This is valid for the wall shear under the
impinging core / developed wall jet, which is what the cleaning criterion needs.

Wall treatment: standard OpenFOAM blended wall functions (nutkWallFunction,
kqRWallFunction, omegaWallFunction). omegaWallFunction blends viscous/log so the
case tolerates a range of y+; actual y+ on the plate is reported by fit_closure.

Usage:  python3 gen_case.py <Re_d> <H_over_d> [refine]
        refine = mesh multiplier (1 = base, 2 = refined for grid-independence).
Writes case dir  case_Re{Re}_H{Hd}[_r{refine}]  and does NOT run it (run_case.sh does).
"""
import sys, os, math

# ---- physical constants (match model/bidet_model.py) ----
NU  = 8.9e-7      # kinematic viscosity of water [m^2/s]
RHO = 1000.0      # density [kg/m^3]  (post-processing multiplies kinematic tau by this)
D   = 1.0e-3      # nozzle diameter [m]

WEDGE_DEG = 2.5                       # wedge half-angle
WA = math.radians(WEDGE_DEG)
COS, SIN = math.cos(WA), math.sin(WA)


def V0_of_Re(Re):
    return Re * NU / D


def geom(H_over_d):
    """Return domain sizes in metres."""
    H = H_over_d * D
    Rjet = D / 2.0
    Rdom = 25.0e-3            # radial extent >= 20 mm (cleaning radius r_max = 15 mm)
    return H, Rjet, Rdom


def head(cls, obj):
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""


def write(path, txt):
    with open(path, "w") as f:
        f.write(txt)


def blockmesh(H, Rjet, Rdom, refine):
    """Two radial blocks (jet, outer), 1 wedge cell thick, graded to the plate."""
    # cell counts (base), scaled by refine
    nx   = int(120 * refine)     # axial (0 -> H), graded fine at the plate
    nrj  = int(16  * refine)     # radial across the jet   [0, Rjet]
    nro  = int(240 * refine)     # radial outer            [Rjet, Rdom]
    # grading (expansionRatio = last/first cell width)
    gx   = 0.02                  # cluster cells at the plate (x = H, the "end")
    gro  = 60.0                  # radial cells grow outward (fine near impingement)

    def v(x, r):  # two wedge vertices (-z, +z) for axial x, radius r
        y = r * COS
        z = r * SIN
        return (x, y, -z), (x, y, z)

    # radius stations: 0, Rjet, Rdom ; axial: 0, H
    pts = []
    idx = {}
    def add(x, r):
        if r == 0.0:
            # axis: minus-z and plus-z collapse to one shared vertex -> prism cells,
            # axis becomes an edge (no patch), wedge apex handled by front/back faces.
            i = len(pts); pts.append((x, 0.0, 0.0))
            idx[(x, r)] = (i, i)
            return
        m, p = v(x, r)
        i0 = len(pts); pts.append(m)
        i1 = len(pts); pts.append(p)
        idx[(x, r)] = (i0, i1)   # (minus-z index, plus-z index)
    for r in (0.0, Rjet, Rdom):
        for x in (0.0, H):
            add(x, r)

    def blk(r0, r1, ncells_r, grade_r):
        (a0m, a0p) = idx[(0.0, r0)]; (a1m, a1p) = idx[(H, r0)]
        (b0m, b0p) = idx[(0.0, r1)]; (b1m, b1p) = idx[(H, r1)]
        # hex order: 4 minus-z verts then 4 plus-z verts
        # face (x0,r0)-(xH,r0)-(xH,r1)-(x0,r1)
        hexv = f"({a0m} {a1m} {b1m} {b0m} {a0p} {a1p} {b1p} {b0p})"
        cells = f"({nx} {ncells_r} 1)"
        grad  = f"simpleGrading ({gx} {grade_r} 1)"
        return f"    hex {hexv} {cells} {grad}\n"

    vstr = "".join(f"    ({p[0]:.8g} {p[1]:.8g} {p[2]:.8g})\n" for p in pts)
    blocks = blk(0.0, Rjet, nrj, 1.0) + blk(Rjet, Rdom, nro, gro)

    # boundary faces (reference vertices by station)
    (aj_m, aj_p) = idx[(0.0, Rjet)]     # x0, Rjet
    (o0_m, o0_p) = idx[(0.0, 0.0)]      # x0, axis
    (oj_m, oj_p) = idx[(H, Rjet)]       # xH, Rjet (not directly used)
    (ax0m, ax0p) = idx[(0.0, 0.0)]
    (axHm, axHp) = idx[(H, 0.0)]
    (jx0m, jx0p) = idx[(0.0, Rjet)]
    (jxHm, jxHp) = idx[(H, Rjet)]
    (Rx0m, Rx0p) = idx[(0.0, Rdom)]
    (RxHm, RxHp) = idx[(H, Rdom)]

    def face(a, b, c, d):
        return f"            ({a} {b} {c} {d})\n"

    # inlet: x=0, r in [0,Rjet] ; quad on minus-z? we give one representative face
    # (blockMesh matches by patch face lists; use the block boundary faces)
    inlet = face(ax0m, jx0m, jx0p, ax0p)
    entr  = face(jx0m, Rx0m, Rx0p, jx0p)
    plate = face(axHm, jxHm, jxHp, axHp) + face(jxHm, RxHm, RxHp, jxHp)
    outlet= face(Rx0m, RxHm, RxHp, Rx0p)
    # axis edge (r=0): shared vertices collapse the face to a line -> no patch needed.
    # wedge faces: minus-z side (all -z verts of both blocks), plus-z side
    (j0m,_)=idx[(0.0,Rjet)]; (jHm,_)=idx[(H,Rjet)]
    front = ("            (%d %d %d %d)\n" % (ax0p, axHp, jxHp, jx0p) +
             "            (%d %d %d %d)\n" % (jx0p, jxHp, RxHp, Rx0p))
    back  = ("            (%d %d %d %d)\n" % (ax0m, jx0m, jxHm, axHm) +
             "            (%d %d %d %d)\n" % (jx0m, Rx0m, RxHm, jxHm))

    boundary = f"""
boundary
(
    inlet
    {{
        type patch;
        faces
        (
{inlet}        );
    }}
    entrainment
    {{
        type patch;
        faces
        (
{entr}        );
    }}
    plate
    {{
        type wall;
        faces
        (
{plate}        );
    }}
    outlet
    {{
        type patch;
        faces
        (
{outlet}        );
    }}
    front
    {{
        type wedge;
        faces
        (
{front}        );
    }}
    back
    {{
        type wedge;
        faces
        (
{back}        );
    }}
);
"""
    body = head("dictionary", "blockMeshDict")
    body += "\nscale 1;\n\nvertices\n(\n" + vstr + ");\n\nblocks\n(\n" + blocks + ");\n\nedges\n(\n);\n" + boundary + "\n// ************************************************************************* //\n"
    return body


def turb_inlet(V0):
    I = 0.05
    k = 1.5 * (I * V0) ** 2
    L = 0.07 * D
    Cmu = 0.09
    omega = math.sqrt(k) / (Cmu ** 0.25 * L)
    return k, omega


def field_U(V0):
    b = head("volVectorField", "U")
    b += f"""
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);

boundaryField
{{
    inlet        {{ type fixedValue; value uniform ({V0:.6g} 0 0); }}
    entrainment  {{ type pressureInletOutletVelocity; value uniform (0 0 0); }}
    outlet       {{ type pressureInletOutletVelocity; value uniform (0 0 0); }}
    plate        {{ type noSlip; }}
    front        {{ type wedge; }}
    back         {{ type wedge; }}
}}
// ************************************************************************* //
"""
    return b


def field_p():
    b = head("volScalarField", "p")
    b += """
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;

boundaryField
{
    inlet        { type zeroGradient; }
    entrainment  { type totalPressure; p0 uniform 0; value uniform 0; }
    outlet       { type fixedValue; value uniform 0; }
    plate        { type zeroGradient; }
    front        { type wedge; }
    back         { type wedge; }
}
// ************************************************************************* //
"""
    return b


def field_k(k):
    b = head("volScalarField", "k")
    b += f"""
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform {k:.6g};

boundaryField
{{
    inlet        {{ type fixedValue; value uniform {k:.6g}; }}
    entrainment  {{ type inletOutlet; inletValue uniform {k:.6g}; value uniform {k:.6g}; }}
    outlet       {{ type inletOutlet; inletValue uniform {k:.6g}; value uniform {k:.6g}; }}
    plate        {{ type kqRWallFunction; value uniform {k:.6g}; }}
    front        {{ type wedge; }}
    back         {{ type wedge; }}
}}
// ************************************************************************* //
"""
    return b


def field_omega(omega):
    b = head("volScalarField", "omega")
    b += f"""
dimensions      [0 0 -1 0 0 0 0];
internalField   uniform {omega:.6g};

boundaryField
{{
    inlet        {{ type fixedValue; value uniform {omega:.6g}; }}
    entrainment  {{ type inletOutlet; inletValue uniform {omega:.6g}; value uniform {omega:.6g}; }}
    outlet       {{ type inletOutlet; inletValue uniform {omega:.6g}; value uniform {omega:.6g}; }}
    plate        {{ type omegaWallFunction; value uniform {omega:.6g}; }}
    front        {{ type wedge; }}
    back         {{ type wedge; }}
}}
// ************************************************************************* //
"""
    return b


def field_nut():
    b = head("volScalarField", "nut")
    b += """
dimensions      [0 2 -1 0 0 0 0];
internalField   uniform 0;

boundaryField
{
    inlet        { type calculated; value uniform 0; }
    entrainment  { type calculated; value uniform 0; }
    outlet       { type calculated; value uniform 0; }
    plate        { type nutkWallFunction; value uniform 0; }
    front        { type wedge; }
    back         { type wedge; }
}
// ************************************************************************* //
"""
    return b


def transportProperties():
    b = head("dictionary", "transportProperties")
    b += f"""
transportModel  Newtonian;
nu              {NU:.6g};
// ************************************************************************* //
"""
    return b


def turbulenceProperties():
    b = head("dictionary", "turbulenceProperties")
    b += """
simulationType  RAS;

RAS
{
    RASModel        kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}
// ************************************************************************* //
"""
    return b


def controlDict(endTime):
    b = head("dictionary", "controlDict")
    b += f"""
application     simpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {endTime};
deltaT          1;
writeControl    timeStep;
writeInterval   {endTime};
purgeWrite      1;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
runTimeModifiable true;

functions
{{
    wallShearStress
    {{
        type            wallShearStress;
        libs            (fieldFunctionObjects);
        patches         (plate);
        writeControl    onEnd;
    }}
    yPlus
    {{
        type            yPlus;
        libs            (fieldFunctionObjects);
        writeControl    onEnd;
    }}
    sampleWSS
    {{
        type            surfaces;
        libs            (sampling);
        writeControl    onEnd;
        surfaceFormat   raw;
        fields          (wallShearStress yPlus);
        interpolate     false;
        surfaces
        {{
            plateSurf
            {{
                type    patch;
                patches (plate);
            }}
        }}
    }}
}}
// ************************************************************************* //
"""
    return b


def fvSchemes():
    b = head("dictionary", "fvSchemes")
    b += """
ddtSchemes      { default steadyState; }
gradSchemes     { default Gauss linear; }
divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,k)      bounded Gauss upwind;
    div(phi,omega)  bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
wallDist        { method meshWave; }
// ************************************************************************* //
"""
    return b


def fvSolution():
    b = head("dictionary", "fvSolution")
    b += """
solvers
{
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          0.01;
    }
    "(U|k|omega)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0.1;
    }
}
SIMPLE
{
    nNonOrthogonalCorrectors 1;
    consistent      yes;
    residualControl
    {
        p           1e-5;
        U           1e-5;
        "(k|omega)" 1e-5;
    }
}
relaxationFactors
{
    equations
    {
        U           0.9;
        "(k|omega)" 0.7;
    }
}
// ************************************************************************* //
"""
    return b


def main():
    Re = float(sys.argv[1])
    Hd = float(sys.argv[2])
    refine = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    V0 = V0_of_Re(Re)
    H, Rjet, Rdom = geom(Hd)
    k, omega = turb_inlet(V0)
    endTime = 3000

    tag = f"case_Re{int(Re)}_H{int(Hd)}" + ("" if refine == 1.0 else f"_r{int(refine)}")
    case = os.path.join(os.path.dirname(os.path.abspath(__file__)), tag)
    for sub in ("0", "constant", "system"):
        os.makedirs(os.path.join(case, sub), exist_ok=True)

    write(os.path.join(case, "system", "blockMeshDict"), blockmesh(H, Rjet, Rdom, refine))
    write(os.path.join(case, "system", "controlDict"), controlDict(endTime))
    write(os.path.join(case, "system", "fvSchemes"), fvSchemes())
    write(os.path.join(case, "system", "fvSolution"), fvSolution())
    write(os.path.join(case, "constant", "transportProperties"), transportProperties())
    write(os.path.join(case, "constant", "turbulenceProperties"), turbulenceProperties())
    write(os.path.join(case, "0", "U"), field_U(V0))
    write(os.path.join(case, "0", "p"), field_p())
    write(os.path.join(case, "0", "k"), field_k(k))
    write(os.path.join(case, "0", "omega"), field_omega(omega))
    write(os.path.join(case, "0", "nut"), field_nut())

    with open(os.path.join(case, "case_info.txt"), "w") as f:
        f.write(f"Re_d={Re}\nH_over_d={Hd}\nrefine={refine}\nV0={V0}\nH={H}\n"
                f"Rjet={Rjet}\nRdom={Rdom}\nk_in={k}\nomega_in={omega}\nNU={NU}\nRHO={RHO}\nD={D}\n")
    print(f"wrote {case}  (Re={Re:g}, H/d={Hd:g}, V0={V0:.4g} m/s, refine={refine:g})")


if __name__ == "__main__":
    main()
