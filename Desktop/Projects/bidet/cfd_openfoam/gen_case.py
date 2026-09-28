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


RDOM_DEFAULT = 25.0e-3        # radial extent >= 20 mm (cleaning radius r_max = 15 mm)


def geom(H_over_d, Rdom=RDOM_DEFAULT):
    """Return domain sizes in metres.

    Rdom was fixed at 25 mm for the 2026-07-25 sweep. That is only ~1.7 r_max, and
    tau_w has fallen just to 1.957 Pa at the outer boundary, so the tau_c = 1 Pa
    cleaning contour lies OUTSIDE the domain and the paper can only bound it
    (r > 25 mm) or extrapolate the fitted closure to it. Rdom is a parameter so
    that case can be re-run on a domain wide enough to compute it directly.
    """
    H = H_over_d * D
    Rjet = D / 2.0
    return H, Rjet, float(Rdom)


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


def _first_cell(length, n, expansion):
    """Width of the first cell of a blockMesh simpleGrading block.

    blockMesh's expansionRatio is last/first cell width, so the per-cell ratio is
    q = expansion**(1/(n-1)) and the block is a geometric series of n terms.
    """
    if n <= 1:
        return float(length)
    q = expansion ** (1.0 / (n - 1))
    return float(length) * (q - 1.0) / (q ** n - 1.0) if abs(q - 1.0) > 1e-15 \
        else float(length) / n


def _grading_for_first_cell(length, n, first, tol=1e-14, iters=200):
    """Inverse of _first_cell: the expansionRatio giving `first` over `length`.

    Monotone in q, so bisect on the per-cell ratio rather than trusting a solver.
    Bracketed well above any grading this mesh uses.
    """
    if n <= 1:
        return 1.0
    lo, hi = 1.0 + 1e-12, 2.0
    while _first_cell(length, n, hi ** (n - 1)) > first and hi < 1e6:
        hi *= 2.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if _first_cell(length, n, mid ** (n - 1)) > first:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return float((0.5 * (lo + hi)) ** (n - 1))


def blockmesh(H, Rjet, Rdom, refine):
    """Two radial blocks (jet, outer), 1 wedge cell thick, graded to the plate."""
    # cell counts (base), scaled by refine
    nx   = int(120 * refine)     # axial (0 -> H), graded fine at the plate
    nrj  = int(16  * refine)     # radial across the jet   [0, Rjet]
    # Radial outer [Rjet, Rdom]. The count scales with the span so that widening the
    # domain does not silently coarsen the wall jet: at the default 25 mm this is the
    # original 240, and a 50 mm domain gets 485.
    span_scale = (Rdom - Rjet) / (RDOM_DEFAULT - D / 2.0)
    nro_default = int(240 * refine)
    nro  = int(round(nro_default * span_scale))
    # grading (expansionRatio = last/first cell width)
    gx   = 0.02                  # cluster cells at the plate (x = H, the "end")
    # Widening the domain must not move the FIRST radial cell, which is what resolves
    # the near-impingement peak where tau_w is 390 Pa. So hold the first-cell width at
    # whatever the 25 mm mesh produces and re-solve the grading for the longer block.
    # Preserving the per-cell growth rate instead would compound over the extra cells
    # and drive the first cell to 0.2 um, which buys nothing and wrecks the aspect
    # ratio. At the default Rdom this returns exactly 60.0, so the validated 2026-07-25
    # mesh is reproduced bit for bit.
    # Rounded to 10 significant figures so the default case writes "60.0" rather than
    # the solver's 59.9999999999796 and its blockMeshDict stays byte-identical to the
    # validated 2026-07-25 one.
    gro  = float(f"{_grading_for_first_cell(Rdom - Rjet, nro, _first_cell(RDOM_DEFAULT - D / 2.0, nro_default, 60.0)):.10g}")

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


def fvSchemes(turb_scheme="upwind"):
    # turb_scheme: "upwind" (first order, the shipped 2026-07-25 default) or
    # "limitedLinear 1" (second-order bounded TVD; PREREG_so_triplet_2026-09-27.md).
    b = head("dictionary", "fvSchemes")
    b += f"""
ddtSchemes      {{ default steadyState; }}
gradSchemes     {{ default Gauss linear; }}
divSchemes
{{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,k)      bounded Gauss {turb_scheme};
    div(phi,omega)  bounded Gauss {turb_scheme};
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}}
laplacianSchemes {{ default Gauss linear corrected; }}
interpolationSchemes {{ default linear; }}
snGradSchemes   {{ default corrected; }}
wallDist        {{ method meshWave; }}
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
    # Flags (2026-09-27, PREREG_so_triplet): --scheme=fo|so selects the k/omega advection
    # scheme AND switches to a non-truncating tag `case_Re.._H.._{fo|so}_r{refine}` (dots
    # as 'p'), so these cases can never collide with the legacy `_r{int(refine)}` names.
    # --end=N sets endTime. Without --scheme the legacy behaviour is byte-identical.
    flags = {a.split("=", 1)[0]: (a.split("=", 1)[1] if "=" in a else "")
             for a in sys.argv[1:] if a.startswith("--")}
    pos = [a for a in sys.argv[1:] if not a.startswith("--")]
    Re = float(pos[0])
    Hd = float(pos[1])
    refine = float(pos[2]) if len(pos) > 2 else 1.0
    Rdom_mm = float(pos[3]) if len(pos) > 3 else RDOM_DEFAULT * 1e3
    V0 = V0_of_Re(Re)
    H, Rjet, Rdom = geom(Hd, Rdom_mm * 1e-3)
    k, omega = turb_inlet(V0)
    endTime = int(flags.get("--end") or 3000)
    scheme = flags.get("--scheme")
    if scheme not in (None, "fo", "so"):
        print(f"unknown --scheme={scheme}; use fo or so", file=sys.stderr)
        return 1
    turb_scheme = "limitedLinear 1" if scheme == "so" else "upwind"

    # A widened domain gets its own tag, so it can never overwrite the 2026-07-25
    # sweep the manuscript's numbers are read from.
    if scheme is None:
        tag = f"case_Re{int(Re)}_H{int(Hd)}" + ("" if refine == 1.0 else f"_r{int(refine)}")
    else:
        tag = f"case_Re{int(Re)}_H{int(Hd)}_{scheme}_r{refine:g}".replace(".", "p")
    if abs(Rdom - RDOM_DEFAULT) > 1e-12:
        tag += f"_R{int(round(Rdom_mm))}"
    case = os.path.join(os.path.dirname(os.path.abspath(__file__)), tag)

    # REFUSE TO OVERWRITE A CASE THAT WAS BUILT AT A DIFFERENT REFINEMENT.
    #
    # The tag above truncates: int(1.414) and int(1.5) are both 1, so every refine in
    # [1, 2) collides on `_r1`, and 0.5 lands on `_r0`. On 2026-09-03 a refine = 1.5 run
    # silently overwrote the refine = 1.414 case's blockMeshDict and case_info.txt while
    # leaving its log.blockMesh (61 009 cells) in place -- an internally inconsistent case
    # that would have been read by gci.py as if it were sound. Recovered from git.
    #
    # The naming scheme is left alone deliberately: renaming would orphan the case
    # directories gci.py's families already reference. This guard is the cheap half and it
    # is the half that stops the damage. Pass --force to override, which is then a
    # deliberate act visible in a shell history.
    info = os.path.join(case, "case_info.txt")
    if os.path.exists(info) and "--force" not in sys.argv:
        prev = {}
        with open(info, encoding="utf-8") as fh:
            for line in fh:
                if "=" in line:
                    k_, v_ = line.strip().split("=", 1)
                    prev[k_] = v_
        was = float(prev.get("refine", "nan"))
        if not (was == refine or (was != was and refine != refine)):
            print(f"REFUSING to overwrite {tag}: it was built at refine={was:g} and you "
                  f"asked for refine={refine:g}. The tag truncates, so these collide. "
                  f"Move or delete the existing case, or pass --force.", file=sys.stderr)
            return 1

    for sub in ("0", "constant", "system"):
        os.makedirs(os.path.join(case, sub), exist_ok=True)

    write(os.path.join(case, "system", "blockMeshDict"), blockmesh(H, Rjet, Rdom, refine))
    write(os.path.join(case, "system", "controlDict"), controlDict(endTime))
    write(os.path.join(case, "system", "fvSchemes"), fvSchemes(turb_scheme))
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
        if scheme is not None:
            f.write(f"turb_scheme={turb_scheme}\nendTime={endTime}\n")
    print(f"wrote {case}  (Re={Re:g}, H/d={Hd:g}, V0={V0:.4g} m/s, refine={refine:g})")
    return 0


if __name__ == "__main__":
    # main() returns 1 when it refuses to overwrite a colliding case; that has to
    # reach the shell, or `set -e` in run_case*.sh will not stop the run.
    sys.exit(main() or 0)
