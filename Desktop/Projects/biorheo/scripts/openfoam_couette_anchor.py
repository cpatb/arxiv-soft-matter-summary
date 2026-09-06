"""External anchor for the L3 OpenFOAM bridge: run OpenFOAM itself on an exported record.

Why this exists. `tests/test_openfoam_export.py` round-trips a record through
`biorheo.solver.openfoam.nu_openfoam`, which is THIS REPO's transcription of OpenFOAM's
`calcNu`. That certifies the dictionary and the parameter map against a transcription,
not against OpenFOAM. This script closes the loop: it writes the exported
`transportProperties` into a real case, lets `simpleFoam` solve, and reads back the
viscosity field OpenFOAM's own `viscosityModel` computed. Nothing in the comparison
passes through `nu_openfoam`.

The case is 1-D plane Couette: gap H, lower wall at rest, upper wall at U0, cyclic in x,
empty in z, laminar, steady. In steady plane Couette the momentum equation is
d(tau)/dy = 0, so tau is uniform; where the STRESS tau(gd) = eta(gd)*gd is strictly
increasing that forces gd uniform, hence gd = U0/H EXACTLY and independent of the
constitutive law. The condition is on tau, not on eta -- every law here is
shear-thinning, so eta is DECREASING and a condition on it would be the wrong one. It is
derived, with the viscosity clamp included and with numeric spot checks at the shipped
record parameters, in `derivations/couette_uniform_shear.py`. The shear rate is
therefore imposed by the boundary condition, not measured, and there is no
closure-dependent profile to converge; A1 below is the check that the solver found this
branch and not a banded one. The discrete problem is linear once nu is uniform, so the
SIMPLE iteration should drive residuals to round-off, not to a tolerance.

PRE-REGISTERED GATES (written before the first run, 2026-09-05; do not loosen):

  A1  U-profile linearity: the cell-centre shear rate (dU/dy by central difference on the
      uniform mesh) is uniform across the gap to <= 1e-6 relative spread. If this fails
      the comparison rate gd = U0/H is not well posed and NOTHING below may be quoted.

  A1 IS NOT SUFFICIENT AND IS KEPT ONLY AS A REPORTED NUMBER (raised 2026-09-05 by the
  code review, before the shipped run). A central difference has the odd-even
  (checkerboard) mode in its null space -- exactly the mode a colocated FV scheme
  admits -- so A1 returns a BIT-IDENTICAL 1.95e-15 on a linear profile and on the same
  profile plus a 0.3-cell checkerboard whose true cell-to-cell shear rate varies by 1.6x.
  A1 also tests uniformity only: nothing in it says the shear rate is U0/H, which is what
  every other gate compares against. Replaced by a strictly stronger test, added before
  the shipped run and never loosened:

  A1' Profile residual: max |Ux_i - U0 * y_i / H| / U0 <= 1e-9 over all cells, against
      the analytic line derived in `derivations/couette_uniform_shear.py`. This sees the
      checkerboard (it is a pointwise residual, not a difference operator) and it sees a
      wrong slope, which A1 cannot.
  A2  Viscosity: OpenFOAM's own `nu` field equals record.eta_at(U0/H)/rho to <= 1e-6
      relative, at every cell, for every (record, U0). Expectation, stated in advance:
      much tighter than the bar, because the discrete solution is exact -- if the measured
      value lands near 1e-6 rather than near round-off, that is a finding to chase, not a
      pass to bank.
  A3  Wall shear stress: the `wallShearStress` function object at the moving wall equals
      eta_at(gd)*gd/rho to <= 1e-6 relative. This is a second, independent channel (it
      goes through the momentum assembly and the FO's own stress convention, A2 does not).
  A4  Newtonian control: a hand-made Newtonian record with nu = 1e-4 m^2/s at U0/H = 10
      must give nu = 1e-4 and tau/rho = 1e-3 to the same bar. This is what shows A2/A3
      would have caught a sign or factor-of-2 convention error in the FO, rather than the
      whole comparison being self-consistently wrong.
  A5  Convergence, not the iteration cap: every case must stop on residualControl before
      endTime. A case that hits the cap is reported FAILED, not quoted.

  A5 FAILED ON THE FIRST RUN (2026-09-05) AND WAS THE WRONG TEST. Recorded, not erased:
  with residualControl 1e-12 and a 2000-iteration cap, 6 of 7 cases hit the cap; the one
  that did not was the Newtonian control (1938 its). Diagnosis (`openfoam_couette_diagnose.py`,
  table in the .md): the outer SIMPLE iteration's convergence RATE is set by the local
  Bingham number B = tau_y/(k gd^n) -- at B = 2.0 it meets 1e-12 in 8123 iterations, at
  B = 7.6 and above it does not meet it in 20000, and at B = 28.9 the normalised residual
  stalls near 7e-5 while the SOLUTION is still improving (relative nu error 1.7 at 2000
  its, 1.4e-5 at 20000, 9.6e-13 at 20000 with relaxation 0.95). So the residual is a bad
  proxy for the answer here, and a bar on it measures the normalisation.

  The replacement is STRICTER, not looser, and is applied to every case uniformly:

  A5' Cauchy convergence of the ANSWER: residualControl is removed, every case runs the
      full 20000 iterations, and the cell shear-rate field at iteration 10000 must equal
      the field at 20000 to <= 1e-9 relative, everywhere. A residual can be small while
      the solution drifts; this cannot. The final Ux initial residual is printed beside
      it as information, never as the criterion.

  A5' ALONE IS ALSO NOT SUFFICIENT (same review). A Cauchy increment is not an error
  bound: for an outer iteration contracting by rho per window the increment is about
  err*(1-rho), and rho -> 1 is precisely the high-Bingham regime under study, so a
  stalled solve sitting on an arbitrarily wrong answer can pass. Five of seven rows
  report exactly 0.0. Two things are done about it, both added before the shipped run:
  the increment is measured over THREE windows (5000, 10000, 20000) and the contraction
  ratio is reported, so a stall is visible as a ratio near 1 rather than as a small
  number; and the honest statement is recorded that A2/A3/A1' are the binding gates here
  because the answer is ANALYTIC -- A5' is a convergence diagnostic, not the correctness
  test, and it is not "stricter than A5", it tests a different thing.

  Two numerics settings changed with A5' and are not bars: relaxation on U 0.7 -> 0.95
  (the sweep says stronger relaxation is WORSE here -- 0.3 gives a relative error of 0.29
  where 0.9 gives 4.5e-12 at the same 20000 iterations) and the cap 2000 -> 20000. There
  are TWO regimes and not one mechanism, which an earlier draft of this docstring got
  wrong: below the shipped setting the outer loop is a slowly damped relaxation, so
  under-relaxing harder only costs iterations; at 0.99 it is genuinely UNSTABLE and
  simpleFoam dies with rc = 136 (SIGFPE) after the continuity error blows up past 1e73.
  The shipped 0.95 therefore sits between a slow regime and a divergent one.

Run:  python scripts/openfoam_couette_anchor.py
Needs WSL with OpenFOAM v2606 at /usr/lib/openfoam/openfoam2606. Writes
`data/openfoam/couette_anchor.md` and `couette_anchor.json`; cases go to
`data/openfoam/runs/` (gitignored). Pinned in `tests/test_openfoam_couette.py`.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np

from biorheo.registry.loader import load_all
from biorheo.solver.openfoam import (coefficients, export_transport_properties,
                                    nu_openfoam, parse_transport_properties)

ROOT = Path(__file__).resolve().parent.parent
# BIORHEO_OF_OUT lets a test run the whole anchor without overwriting the COMMITTED
# artifacts it is validating (raised by the code review 2026-09-05).
OUT = Path(os.environ.get("BIORHEO_OF_OUT") or (ROOT / "data" / "openfoam"))
RUNS = OUT / "runs"
OF_BASHRC = "/usr/lib/openfoam/openfoam2606/etc/bashrc"

RHO = 1060.0        # kg/m^3, the test density used throughout the L3 bridge
H = 1.0e-3          # m, gap
NY = 20             # cells across the gap
LX = LZ = 1.0e-3    # m, one cell each (cyclic / empty)
N_RATES = 8         # wall speeds per record, geometric over the whole measured envelope
END_TIME = 20000    # fixed iteration count; residualControl is NOT used (see A5')
RELAX_U = 0.95      # 0.99 dies with SIGFPE; 0.3 converges far more slowly than 0.9

# pre-registered bars
BAR_A1 = 1e-6
BAR_A1P = 1e-9      # A1': max |Ux - U0 y/H| / U0, the test A1 should have been
BAR_A2 = 1e-6
BAR_A3 = 1e-6
BAR_A5 = 1e-9       # A5': relative change of the shear-rate field over the last half run
WINDOWS = (4, 2, 1)  # END_TIME // w: the three writes A5' measures contraction over
ROUNDOFF = 1e-13     # below this a relative increment is float noise, not a signal

RECORD_IDS = ["carbopol_etd2020_0p1wt_30uL", "blood_horner2020_donor_A",
              # Added 2026-09-06 (NEXT_STEPS item 6): one record of each family that has
              # NO yield stress, so B = tau_y/(k gd^n) is identically zero. They are the
              # control for the finding that simpleFoam's convergence failure is about
              # YIELD and not about "hard fluid".
              #
              # Neither is an arbitrary pick. Both are the record of their family with
              # every parameter graded `measurable` (so no allow_unidentified override is
              # needed) that maximises the separation between the shipped map and the
              # candidate the gate has to reject: for `_cross` that is the name-swapped
              # map (G3), and for `_carreau_yasuda` it is `a` forced to 2 (G4, 1.01e-1
              # relative -- 1e5 times BAR_A2, so A2 can decide whether OpenFOAM honours
              # `a` at all). Both numbers are printed by the script.
              "blood_horner2020_donor_N_cross",
              "blood_horner2020_donor_S_carreau_yasuda"]
# The families whose `nu0` is a FITTED PARAMETER rather than a clamp; `eta_max` is None
# for them and `clamped_frac` is meaningless, so it is not read.
NO_CLAMP_MODELS = ("cross", "carreau_yasuda")


# --------------------------------------------------------------------------- case files

_FOAM = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
"""


def _block_mesh_dict(h: float, ny: int) -> str:
    return _FOAM.format(cls="dictionary", obj="blockMeshDict") + f"""
scale   1;
vertices
(
    (0 0 0) ({LX:.17g} 0 0) ({LX:.17g} {h:.17g} 0) (0 {h:.17g} 0)
    (0 0 {LZ:.17g}) ({LX:.17g} 0 {LZ:.17g}) ({LX:.17g} {h:.17g} {LZ:.17g}) (0 {h:.17g} {LZ:.17g})
);
blocks ( hex (0 1 2 3 4 5 6 7) (1 {ny} 1) simpleGrading (1 1 1) );
edges ();
boundary
(
    lowerWall {{ type wall; faces ((0 1 5 4)); }}
    upperWall {{ type wall; faces ((3 7 6 2)); }}
    inlet  {{ type cyclic; neighbourPatch outlet; faces ((0 4 7 3)); }}
    outlet {{ type cyclic; neighbourPatch inlet;  faces ((1 2 6 5)); }}
    frontAndBack {{ type empty; faces ((0 3 2 1) (4 5 6 7)); }}
);
mergePatchPairs ();
"""


def _u_file(u0: float) -> str:
    return _FOAM.format(cls="volVectorField", obj="U") + f"""
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);
boundaryField
{{
    lowerWall    {{ type fixedValue; value uniform (0 0 0); }}
    upperWall    {{ type fixedValue; value uniform ({u0:.17g} 0 0); }}
    inlet        {{ type cyclic; }}
    outlet       {{ type cyclic; }}
    frontAndBack {{ type empty; }}
}}
"""


_P_FILE = _FOAM.format(cls="volScalarField", obj="p") + """
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{
    lowerWall    { type zeroGradient; }
    upperWall    { type zeroGradient; }
    inlet        { type cyclic; }
    outlet       { type cyclic; }
    frontAndBack { type empty; }
}
"""

_TURB = _FOAM.format(cls="dictionary", obj="turbulenceProperties") + """
simulationType  laminar;
"""

_FVSCHEMES = _FOAM.format(cls="dictionary", obj="fvSchemes") + """
ddtSchemes      { default steadyState; }
gradSchemes     { default Gauss linear; }
divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linear;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
"""

_FVSOLUTION = _FOAM.format(cls="dictionary", obj="fvSolution") + f"""
solvers
{{
    p {{ solver PCG; preconditioner DIC; tolerance 1e-14; relTol 0; }}
    U {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-14; relTol 0; }}
}}
SIMPLE
{{
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       0;
    // no residualControl: every case runs the full endTime and is judged on A5'
}}
relaxationFactors {{ fields {{ p 0.3; }} equations {{ U {RELAX_U:g}; }} }}
"""


def _control_dict() -> str:
    return _FOAM.format(cls="dictionary", obj="controlDict") + f"""
application     simpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {END_TIME};
deltaT          1;
writeControl    timeStep;
writeInterval   {END_TIME // max(WINDOWS)};
purgeWrite      0;
writeFormat     ascii;
writePrecision  17;
writeCompression off;
timeFormat      general;
runTimeModifiable false;
functions
{{
    writeNu
    {{
        type            writeObjects;
        libs            (utilityFunctionObjects);
        objects         (nu);
        writeOption     anyWrite;
        writeControl    writeTime;   // every write, not just the last: the clamp
                                     // fraction is read at the first window too
    }}
    wallShearStress
    {{
        type            wallShearStress;
        libs            (fieldFunctionObjects);
        patches         (upperWall lowerWall);
        writeControl    onEnd;
    }}
}}
"""


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write(text)


def build_case(case: Path, transport_text: str, u0: float) -> None:
    if case.exists():
        shutil.rmtree(case)
    _write(case / "system" / "blockMeshDict", _block_mesh_dict(H, NY))
    _write(case / "system" / "controlDict", _control_dict())
    _write(case / "system" / "fvSchemes", _FVSCHEMES)
    _write(case / "system" / "fvSolution", _FVSOLUTION)
    _write(case / "constant" / "transportProperties", transport_text)
    _write(case / "constant" / "turbulenceProperties", _TURB)
    _write(case / "0" / "U", _u_file(u0))
    _write(case / "0" / "p", _P_FILE)


# ------------------------------------------------------------------------------- runner

def wsl_path(p: Path) -> str:
    p = Path(p).resolve()
    return f"/mnt/{p.drive[0].lower()}{p.as_posix()[2:]}"


def has_openfoam() -> bool:
    if shutil.which("wsl.exe") is None:
        return False
    r = subprocess.run(["wsl.exe", "-e", "test", "-f", OF_BASHRC], capture_output=True)
    return r.returncode == 0


def run_case(case: Path) -> str:
    """blockMesh then simpleFoam, launched from a script FILE in the case (never an
    inline wsl.exe expansion). Returns the simpleFoam log."""
    sh = case / "run.sh"
    _write(sh, "#!/bin/bash\n"
               f"source {OF_BASHRC} > /dev/null 2>&1\n"
               "set -e\n"
               f"cd {wsl_path(case)}\n"
               "blockMesh > log.blockMesh 2>&1\n"
               "simpleFoam > log.simpleFoam 2>&1\n")
    r = subprocess.run(["wsl.exe", "-e", "bash", wsl_path(sh)],
                       capture_output=True, text=True, timeout=900)
    log = (case / "log.simpleFoam")
    if r.returncode != 0:
        tail = log.read_text(errors="replace")[-2000:] if log.exists() else r.stderr
        raise RuntimeError(f"{case.name} failed (rc={r.returncode}):\n{tail}")
    return log.read_text(errors="replace")


def read_internal_field(path: Path) -> np.ndarray:
    """Cell values of an ascii volScalarField (uniform or nonuniform List<scalar>)."""
    txt = path.read_text()
    i = txt.index("internalField")
    body = txt[i:]
    if body[: body.index(";")].strip().split()[1] == "uniform":
        return np.array([float(body[: body.index(";")].split()[-1])])
    lo = body.index("(")
    hi = body.index(")", lo)
    return np.array([float(v) for v in body[lo + 1: hi].split()])


_PAREN_GROUP = re.compile(r"\(([^()]*)\)")


def _vectors(body: str) -> np.ndarray:
    """Vectors from `uniform (a b c)`, a multi-line `nonuniform List<vector> N (...)`, or
    the ONE-LINE list form OpenFOAM writes for short lists (`N((a b c)(d e f))`). The
    one-line form is not reachable from the cases here -- a 1-face patch is written
    `uniform` -- but it was a hard failure rather than a wrong number if it ever were."""
    lo = body.index("(")
    head = body[:lo]
    if "uniform" in head and "nonuniform" not in head:
        return np.array([[float(x) for x in body[lo + 1: body.index(")", lo)].split()]])
    if "\n)" not in body[lo:]:
        depth, end = 0, len(body)
        for i in range(lo, len(body)):
            depth += (body[i] == "(") - (body[i] == ")")
            if depth == 0:
                end = i
                break
        return np.array([[float(x) for x in g.split()]
                         for g in _PAREN_GROUP.findall(body[lo + 1: end])])
    hi = body.index("\n)", lo)
    rows = [ln.strip().strip("()").split() for ln in body[lo + 1: hi].splitlines() if "(" in ln]
    return np.array([[float(x) for x in r] for r in rows])


def read_vector_field_x(path: Path) -> np.ndarray:
    txt = path.read_text()
    return _vectors(txt[txt.index("internalField"):])[:, 0]


def read_wall_shear_x(case: Path, time: str, patch: str) -> float:
    """SIGNED x-component of the boundary wallShearStress (tau/rho, m^2/s^2) on `patch`.

    Signed, not a magnitude: a magnitude cannot detect a sign error, and the Newtonian
    control's whole job is to catch a sign or factor-of-2 in OpenFOAM's stress
    convention. It also excludes the y and z components, which are ~1e-20 here but are
    not negligible in a badly converged case, where a norm would inflate the value.
    """
    txt = (case / time / "wallShearStress").read_text()
    i = txt.index("\n    " + patch + "\n")
    seg = txt[i: txt.index("}", i)]
    v = _vectors(seg[seg.index("value"):])
    return float(v[:, 0].mean())


def last_time_dir(case: Path) -> str:
    times = [d.name for d in case.iterdir() if d.is_dir() and d.name.replace(".", "").isdigit()]
    return max(times, key=float)


def final_ux_residual(log: str) -> float:
    """The last `Solving for Ux, Initial residual = ...` in the log. INFORMATION only:
    the first run showed this quantity stalling near 7e-5 while the solution was still
    converging, which is why A5' is a Cauchy test on the answer instead."""
    vals = [ln.split("Initial residual = ")[1].split(",")[0]
            for ln in log.splitlines() if "Solving for Ux" in ln]
    return float(vals[-1]) if vals else float("nan")


def shear_rate_field(case: Path, time: str) -> np.ndarray:
    """Cell-centre dU/dy on the uniform mesh (central differences)."""
    return np.gradient(read_vector_field_x(case / time / "U"), H / NY)


# -------------------------------------------------------------------------------- cases

def one_case(name: str, transport_text: str, gd: float, eta_ref: float,
             has_clamp: bool = True) -> dict:
    """Run at wall speed gd*H; compare OpenFOAM's nu and tau_w with eta_ref (Pa s)."""
    case = RUNS / name
    build_case(case, transport_text, u0=gd * H)
    log = run_case(case)
    t = last_time_dir(case)
    gd_cells = shear_rate_field(case, t)
    ux = read_vector_field_x(case / t / "U")
    nu = read_internal_field(case / t / "nu")
    assert nu.size in (1, NY), f"nu has {nu.size} cells, expected 1 (uniform) or {NY}"
    assert ux.size == NY, f"U has {ux.size} cells, expected {NY}"
    nu_ref = eta_ref / RHO
    err_nu = float(np.abs(nu - nu_ref).max() / nu_ref)

    # A1': pointwise residual against the analytic line. Unlike A1 this sees a
    # checkerboard mode (it is not a difference operator) and it sees a wrong slope.
    y_c = (np.arange(NY) + 0.5) * (H / NY)
    err_profile = float(np.abs(ux - gd * y_c).max() / (gd * H))

    # nu0 is the clamp the dictionary actually carries. The Newtonian control has no
    # Coeffs block at all and therefore no clamp; that is a legitimate None, not a
    # parse failure, so only that one case is allowed to raise here.
    # The clamp key is `nu0` for HerschelBulkley and `nuMax` for Casson. CrossPowerLaw
    # and BirdCarreau also carry a `nu0`, but there it is the record's own fitted
    # zero-shear plateau and NOT a clamp, so `has_clamp=False` keeps this from reporting
    # "the cap binds on x % of cells" about a number that is not a cap.
    try:
        c = parse_transport_properties(
            (case / "constant" / "transportProperties").read_text())[1]
        nu0 = c.get("nu0", c.get("nuMax")) if has_clamp else None
    except ValueError:
        nu0 = None

    def _clamped(field: np.ndarray) -> float | None:
        return (float(np.mean(np.isclose(field, nu0, rtol=1e-9))) if nu0 is not None
                else None)

    # Setting eta_max well above the CONVERGED viscosity guarantees nothing about the
    # transient, and the transient is where the low-rate case spends its whole budget.
    # Reported at the final time (where err_nu is read) and at the first window.
    clamped_final = _clamped(nu)
    windows = [case / str(END_TIME // w) for w in WINDOWS]
    if all(d.is_dir() for d in windows) and t == windows[-1].name:
        fields = [shear_rate_field(case, d.name) for d in windows]
        scale = abs(fields[-1]).mean()
        incs = [float(np.abs(fields[i + 1] - fields[i]).max() / scale) for i in range(2)]
        cauchy = incs[-1]                                    # the last window: A5'
        # A contraction ratio near 1 is a STALL; near 0 is convergence -- but ONLY where
        # the increments are above round-off. Once a case has converged both increments
        # are ~1e-15 and their ratio is noise that happens to sit near 1, which is the
        # same value a stall gives. Report NaN there rather than a number that reads as
        # a stall: the honest statement is that the ratio is undefined at round-off.
        contraction = (float(incs[1] / incs[0]) if incs[0] > ROUNDOFF else float("nan"))
        clamped_early = _clamped(read_internal_field(windows[0] / "nu"))
    else:
        # only reachable from the diagnosis script, which restores residualControl and so
        # can stop before the writes. NaN fails the A5' comparison, never passes it.
        cauchy, contraction, clamped_early = float("nan"), float("nan"), None
    interior = gd_cells[1:-1]                            # ends are one-sided by construction
    spread = float((interior.max() - interior.min()) / abs(interior.mean()))

    tau_ref = eta_ref * gd / RHO
    # SIGNED x-traction on each wall. A magnitude cannot see a sign error, which is what
    # the Newtonian control is supposed to be able to catch; the two walls must carry
    # equal and opposite traction because d(tau)/dy = 0.
    tau_up = read_wall_shear_x(case, t, "upperWall")
    tau_lo = read_wall_shear_x(case, t, "lowerWall")
    err_tau = float(abs(tau_up - tau_ref) / tau_ref)
    imbalance = float(abs(tau_up + tau_lo) / (0.5 * (abs(tau_up) + abs(tau_lo))))
    return {"case": name, "gammadot": gd, "U0": gd * H, "iterations": int(float(t)),
            "cauchy": cauchy, "contraction": contraction,
            "final_ux_residual": final_ux_residual(log),
            "err_profile": err_profile, "wall_imbalance": imbalance,
            "cap_bound_fraction_early": clamped_early,
            "cap_bound_fraction_final": clamped_final,
            "shear_spread": spread, "nu_openfoam_min": float(nu.min()),
            "nu_openfoam_max": float(nu.max()), "nu_expected": nu_ref,
            "err_nu": err_nu, "tau_openfoam": tau_up, "tau_lower_wall": tau_lo,
            "tau_expected": tau_ref, "err_tau": err_tau, "n_cells_nu": int(nu.size)}


def newtonian_control() -> dict:
    """A4: nu = 1e-4 m^2/s at gd = 10 1/s -> nu = 1e-4, tau/rho = 1e-3, by hand."""
    nu0 = 1.0e-4
    txt = ("FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
           "    class       dictionary;\n    object      transportProperties;\n}\n"
           "// biorheo Newtonian control, hand-written, not an export\n\n"
           "transportModel  Newtonian;\n\n"
           f"nu              [0 2 -1 0 0 0 0] {nu0:.17g};\n")
    r = one_case("newtonian_control", txt, gd=10.0, eta_ref=nu0 * RHO)
    r["hand_nu"] = nu0
    r["hand_tau"] = nu0 * 10.0
    return r


def diagnosis_section() -> list[str]:
    """The 'why A5 failed' table, rendered from `couette_diagnosis.json` (written by
    `openfoam_couette_diagnose.py`). Every number comes from that file; if it is missing
    the report says so instead of quoting remembered values."""
    path = ROOT / "data" / "openfoam" / "couette_diagnosis.json"
    head = ["## The finding the first run produced (kept, because it is the product point)",
            "",
            "An identified record exported to OpenFOAM is not automatically usable in "
            "OpenFOAM. On this trivial 20-cell 1-D problem the convergence rate of "
            "`simpleFoam`'s outer iteration is set by the local Bingham number "
            "`B = tau_y/(k gd^n)`, and the failure is SILENT: `rc = 0`, no warning, a "
            "written result."]
    if not path.exists():
        return head + ["", "*(No `couette_diagnosis.json` on disk. Run "
                           "`python scripts/openfoam_couette_diagnose.py` and re-run this "
                           "script; the table is not reproduced from memory.)*"]
    d = json.loads(path.read_text())
    p = d["params"]
    out = head + ["", f"Measured on `{d['record']}` (tau_y = {p['tau_y']:g} Pa, "
                      f"k = {p['k']:g} Pa s^n, n = {p['n']:g}) by "
                      "`scripts/openfoam_couette_diagnose.py`, with `residualControl "
                      f"{d['residual_control']:g}` restored so 'did it meet the "
                      "criterion' is measurable.", "",
                  f"### Numerics knobs at gd = {d['knob_sweep'][0]['gammadot']:.5g} 1/s "
                  f"(B = {d['knob_sweep'][0]['B']:.4g}), the hardest corner", "",
                  "| relax U | cap | its | met 1e-12 | clamp saturated | shear spread | "
                  "nu rel err |", "|---|---|---|---|---|---|---|"]
    for r in d["knob_sweep"]:
        f = r.get("cap_bound_fraction_final")
        out.append(f"| {r['relax_U']:g} | {r['end_time']} | {r['iterations']} | "
                   f"{'yes' if r['residual_control_met'] else 'no'} | "
                   f"{'n/a' if f is None else f'{f:.2f}'} | "
                   f"{r['shear_spread']:.3e} | {r['err_nu']:.3e} |")
    for c in d["crashes"]:
        out.append(f"| {c['relax_U']:g} | {c['end_time']} | CRASH | - | - | - | "
                   f"simpleFoam died, {c['error'].splitlines()[0][-40:]} |")
    out += ["", f"### Shear rate at relaxation {d['relax_for_rate_sweep']:g}, "
                "cap 20000: where it starts to fail", "",
            "| gd (1/s) | B | its | met 1e-12 | shear spread | nu rel err |",
            "|---|---|---|---|---|---|"]
    for r in d["rate_sweep"]:
        out.append(f"| {r['gammadot']:.5g} | {r['B']:.4g} | {r['iterations']} | "
                   f"{'yes' if r['residual_control_met'] else 'no'} | "
                   f"{r['shear_spread']:.3e} | {r['err_nu']:.3e} |")
    # The headline must come from a row where the clamp did NOT saturate. err_nu = 9.0
    # exactly on the relax-0.3 / 2000 row is nu0/nu_ref - 1 = 10 - 1, i.e. a readout of
    # the caller's own eta_max = 10 * eta_at(gmin), not a measure of the solver: with
    # eta_max = 100x the same run would print "9900 % wrong". Flagged by the code review
    # 2026-09-05. Saturated rows stay in the table, marked, and are never the headline.
    SAT = 0.5   # a row is clamp-saturated if this fraction of cells sat on the cap
    unsat = [r for r in d["knob_sweep"]
             if (r.get("cap_bound_fraction_final") or 0.0) < SAT]
    worst = max(unsat or d["knob_sweep"], key=lambda r: r["err_nu"])
    sat = [r for r in d["knob_sweep"] if r not in unsat]
    # "at the same cap" must be true BY CONSTRUCTION, not by coincidence in this data
    cap = max(r["end_time"] for r in d["knob_sweep"])
    at_cap = [r for r in d["knob_sweep"] if r["end_time"] == cap]
    soft = min((r for r in at_cap if r in unsat and r["relax_U"] < 0.5),
               key=lambda r: r["err_nu"], default=None)
    best = min(at_cap, key=lambda r: r["err_nu"])
    ok_rows = [r for r in d["rate_sweep"] if r["residual_control_met"]]
    out += ["",
            f"Read three things off it. (1) At B = {worst['B']:.3g} with relaxation "
            f"{worst['relax_U']:g} and a {worst['end_time']}-iteration budget the solver "
            f"returns a viscosity {100 * worst['err_nu']:.0f} % wrong with no error. "
            + (f"({len(sat)} row(s) with a LARGER error are excluded from that headline "
               "because the clamp saturated there, which makes the error a readout of "
               "the chosen `eta_max` rather than of the solve; they are marked in the "
               "table.) " if sat else "")
            + (f"(2) Under-relaxing HARDER makes it worse: relaxation {soft['relax_U']:g} "
               f"is {soft['err_nu'] / best['err_nu']:.3g}x worse than "
               f"{best['relax_U']:g} at the same {cap}-iteration cap -- the usual reflex "
               "is the wrong one here. " if soft is not None and best["err_nu"] > 0 else "")
            + (f"But there are TWO regimes, not one mechanism: below the shipped setting "
               "the outer loop is a slowly damped relaxation (small factor, many "
               f"iterations), while at relaxation {d['crashes'][0]['relax_U']:g} it is "
               "genuinely unstable -- simpleFoam dies with SIGFPE after the continuity "
               "error blows up. The working 0.95 sits between a slow regime and an "
               "unstable one, which is a narrow place to be and is stated rather than "
               "hidden. " if d["crashes"] else "")
            + "(3) The normalised residual is not a usable stopping criterion here: "
            + (f"only the B = {min(r['B'] for r in ok_rows):.3g} case ever meets 1e-12, "
               if ok_rows else "no case in the rate sweep meets 1e-12, ")
            + "while the answer improves by many decades past the point where the "
              "residual stops moving. That is why A5' is a Cauchy test on the field.",
            "",
            "The practical consequence for the bridge: the low-rate end of a yield-stress "
            "record is where a solver hand-off goes wrong, and it is exactly the end a "
            "rheometer measures worst. `export_report`'s `cap_binds` flag does not catch "
            "it -- the clamp was 10x clear of the envelope in every row above."]
    # Is there a detector a user WITHOUT the analytic answer could have run? The wall
    # traction imbalance is the only quantity here that needs no reference value.
    det = [r for r in d["knob_sweep"] + d["rate_sweep"] if "wall_imbalance" in r]
    if det:
        err = np.array([r["err_nu"] for r in det])
        imb = np.array([r["wall_imbalance"] for r in det])
        rho_s = float(np.corrcoef(np.argsort(np.argsort(err)),
                                  np.argsort(np.argsort(imb)))[0, 1])
        ratio = imb / err
        out += ["", "### A detector that needs no reference answer", "",
                "Every gate above compares against `record.eta_at`, which a real case "
                "does not have. The one quantity here that does not is the WALL TRACTION "
                "IMBALANCE |tau_upper + tau_lower| / mean|tau|: the momentum balance that "
                "makes tau uniform requires the two walls to carry equal and opposite "
                "traction, so any departure is self-evident error. Over the "
                f"{len(det)} diagnosis cases, spanning "
                f"{err.max() / err.min():.3g}x in the true error, its Spearman rank "
                f"correlation with the true `nu` error is {rho_s:.4g}, and it runs "
                f"{ratio.min():.3g} to {ratio.max():.3g} times that error -- i.e. it is a "
                "conservative flag, one to two decades BELOW the error it is signalling, "
                "not a calibrated estimate of it. Used as `imbalance >> round-off`, it "
                "would have caught every bad case in this study.", "",
                "| relax U | cap | gd (1/s) | true nu rel err | wall imbalance |",
                "|---|---|---|---|---|"]
        for r in det:
            out.append(f"| {r['relax_U']:g} | {r['end_time']} | {r['gammadot']:.5g} | "
                       f"{r['err_nu']:.3e} | {r['wall_imbalance']:.3e} |")
    g = d.get("discriminator")
    if g:
        out += ["", "### Is it really B, or just the viscosity? (the confound, tested)", "",
                "The rate sweep above moves B, the viscosity eta and the wall speed U0 "
                "together, so on its own it cannot attribute anything to B. This 2x2 "
                f"fixes gd = {g['gammadot']:g} 1/s (hence U0, the mesh and the boundary "
                "conditions) and picks synthetic Herschel-Bulkley parameters that hit a "
                "target (B, eta) pair: at fixed gd, `k gd^n = eta gd/(1+B)` and "
                "`tau_y = B eta gd/(1+B)`, so scaling tau_y and k by a common factor "
                "moves eta at constant B. The verdict rule "
                f"(a {g['separation_factor']:g}x separation) was fixed before the run.",
                "", "| B | eta (Pa s) | tau_y (Pa) | k (Pa s^n) | nu rel err |",
                "|---|---|---|---|---|"]
        for c in g["cells"]:
            out.append(f"| {c['B']:g} | {c['eta']:g} | {c['tau_y']:.6g} | {c['k']:.6g} | "
                       f"{c['err_nu']:.4e} |")
        rb = ", ".join(f"{v:.4g}x at eta = {kk}" for kk, v in g["R_B_at_fixed_eta"].items())
        re_ = ", ".join(f"{v:.6g}x at B = {kk}" for kk, v in g["R_eta_at_fixed_B"].items())
        out += ["", f"Across B at fixed eta the error moves {rb}; across a 100x change in "
                    f"eta at fixed B it moves {re_}. **Verdict: {g['verdict']}.** The two "
                    "eta rows agree to the digits printed, which is what a genuine "
                    "invariance looks like -- the solve is unchanged by a common "
                    "rescaling of (tau_y, k), so it is the yield-to-viscous RATIO that "
                    "sets the convergence, not the size of the viscosity, the wall speed "
                    "or the matrix conditioning."]
    return out


def item6_section(item6: dict, nc_rows: list[dict]) -> list[str]:
    """The families with no yield stress: G3/G4/G7 of
    `data/openfoam/nonyield_bridge_prereg.md`, written 2026-09-06 before the code."""
    if not nc_rows:
        return ["## The families with no yield stress", "",
                "No `cross` or `carreau_yasuda` record was run; G3/G4/G7 not evaluated."]
    seps = item6["candidate_separation"]
    out = ["## The families with no yield stress (Cross, Carreau-Yasuda)", "",
           f"{item6['n_cases']} cases over {len(item6['records'])} records: "
           + ", ".join(f"`{r}`" for r in item6["records"]) + ".", "",
           "| gate | what | number | bar | verdict |", "|---|---|---|---|---|",
           f"| G7 | worst OpenFOAM `nu` vs `eta_at/rho` on these records | "
           f"{item6['G7_worst_err_nu']:.3e} | {item6['G7_bar']:g} | "
           f"{'PASS' if item6['G7_pass'] else 'FAILED'} |",
           f"| G3/G4 | smallest separation between the shipped map and the rejected "
           f"candidate | {item6['G3_G4_min_separation']:.3e} | >= "
           f"{item6['G3_G4_bar']:g} | "
           f"{'informative' if item6['G3_G4_informative'] else 'UNINFORMATIVE'} |", ""]
    for rid, s in seps.items():
        what = ("the name-swapped Cross map (biorheo `m` into OpenFOAM's `m`, `lam` into "
                "`n`)" if recs_model(rid) == "cross" else
                "`a` forced to 2, i.e. what OpenFOAM would compute if it ignored the entry")
        out.append(f"- `{rid}`: {what} differs from the record by {s:.3e} relative "
                   f"over the {len(nc_rows) // max(len(seps), 1)} anchor rates, "
                   f"{s / BAR_A2:.1e} times A2's bar.")
    out += ["",
            "**G7 is the control for the Bingham-number finding above.** Neither family "
            "has a yield stress, so `B = tau_y/(k gd^n)` is identically zero on every one "
            f"of these cases, and the worst viscosity error is {item6['G7_worst_err_nu']:.3e} "
            "against the 1.72 (172 %) the Herschel-Bulkley record returns at B = 29 with "
            "`rc = 0`. The failure is about YIELD, not about a large viscosity: the "
            "Cross record's zero-shear plateau is itself a large viscosity and it "
            "converges. The worst wall-traction imbalance on these cases is "
            f"{item6['worst_wall_imbalance']:.3e}.", "",
            "**G4 decides a premise of NEXT_STEPS item 6 by running the solver.** That "
            "item was written believing OpenFOAM v2606's `BirdCarreau` has no `a` "
            "exponent, so that the export would have to REFUSE a Carreau-Yasuda record "
            "with `a != 2`. It has one -- `a_` is "
            "`getOrDefault(\"a\", dimensionedScalar(\"a\", dimless, 2))` "
            "(BirdCarreau.C:81-88) and `calcNu` uses it. The Carreau-Yasuda record run "
            "here has a fitted `a` far from 2, so if OpenFOAM ignored the entry A2 would "
            f"fail by about {seps.get(item6['records'][-1], float('nan')):.1e}; it passes "
            f"at {max(r['err_nu'] for r in nc_rows):.1e}. **The premise is refuted and no "
            "`a != 2` refusal is shipped.** What is shipped instead is the mirror hazard: "
            "because `a` is optional with a silent default, the export always WRITES it "
            "and the parser fills a missing one with 2 rather than reporting it missing.",
            "",
            "**G3 is the discriminator for the Cross map.** biorheo's `lam` (a time) is "
            "OpenFOAM's `m` and biorheo's `m` (an exponent) is OpenFOAM's `n` -- the names "
            "collide with the meanings swapped, and a map written by matching names would "
            "raise nothing, because both slots take a plain scalar. Without this number "
            "the round trip passing at 1e-10 would not show the round trip can tell the "
            "two maps apart.", "",
            "**These records take no `eta_max`.** Their zero-shear plateau `eta_0` is a "
            "fitted parameter of the record, so the exported `nu0` IS the record and there "
            "is nothing to clamp; `coefficients` refuses an `eta_max` for them rather than "
            "accepting and ignoring one."]
    return out


_MODEL_OF_ID: dict[str, str] = {}


def recs_model(rid: str) -> str:
    return _MODEL_OF_ID.get(rid, "")


def main() -> int:
    if not has_openfoam():
        print("SKIP: no WSL OpenFOAM v2606")
        return 0
    OUT.mkdir(parents=True, exist_ok=True)
    RUNS.mkdir(parents=True, exist_ok=True)
    recs = load_all()
    _MODEL_OF_ID.update({rid: recs[rid].model.name for rid in RECORD_IDS})
    rows = [newtonian_control()]
    for rid in RECORD_IDS:
        rec = recs[rid]
        gmin, gmax = rec.envelope.gammadot_min, rec.envelope.gammadot_max
        has_clamp = rec.model.name not in NO_CLAMP_MODELS
        # cap far above the envelope for the yield families; the other two take none
        eta_max = 10.0 * float(rec.eta_at(gmin)) if has_clamp else None
        text = export_transport_properties(rec, RHO, eta_max)
        for k, gd in enumerate(np.geomspace(gmin, gmax, N_RATES)):
            eta = float(rec.eta_at(np.array([gd]))[0])
            rows.append(one_case(f"{rid}_g{k}", text, float(gd), eta, has_clamp=has_clamp)
                        | {"record": rid, "model": rec.model.name})

    a1 = max(r["shear_spread"] for r in rows)
    a1p = max(r["err_profile"] for r in rows)
    a2 = max(r["err_nu"] for r in rows)
    a3 = max(r["err_tau"] for r in rows)
    ctl = rows[0]
    a4 = max(abs(ctl["nu_openfoam_max"] - ctl["hand_nu"]) / ctl["hand_nu"],
             abs(ctl["nu_openfoam_min"] - ctl["hand_nu"]) / ctl["hand_nu"],
             abs(ctl["tau_openfoam"] - ctl["hand_tau"]) / ctl["hand_tau"])
    a5 = float(np.max([r["cauchy"] for r in rows]))  # np.max propagates NaN; max() may not
    # a sign error on either wall is a failure, not a small number: the moving wall must
    # pull in +x and the stationary one must react in -x
    signs_ok = all(r["tau_openfoam"] > 0 > r["tau_lower_wall"] for r in rows)

    # ---- item 6: the families with no yield stress (G3, G4, G7 of the pre-registration)
    # A separation is measured between the SHIPPED law -- taken as `record.eta_at`, which
    # is not nu_openfoam -- and the candidate the gate has to reject, which is evaluated
    # through nu_openfoam. Only the rejected candidate goes through this repo's
    # transcription, so a bug there could understate a separation but could not manufacture
    # a passing verdict.
    nc_rows = [r for r in rows if r.get("model") in NO_CLAMP_MODELS]
    seps: dict[str, float] = {}
    for rid in RECORD_IDS:
        rec = recs[rid]
        if rec.model.name not in NO_CLAMP_MODELS:
            continue
        gds = np.geomspace(rec.envelope.gammadot_min, rec.envelope.gammadot_max, N_RATES)
        nu_ref = rec.eta_at(gds) / RHO
        c = coefficients(rec, RHO, None)
        if rec.model.name == "cross":
            # G3: the map a name-matcher would write -- biorheo's exponent `m` into
            # OpenFOAM's `m` (a time) and biorheo's `lam` into `n` (an exponent).
            with np.errstate(over="ignore"):
                bad = nu_openfoam("CrossPowerLaw",
                                  {**c, "m": c["n"], "n": c["m"]}, gds)
        else:
            # G4: `a` forced to the value it would take if OpenFOAM ignored the entry.
            bad = nu_openfoam("BirdCarreau", {**c, "a": 2.0}, gds)
        seps[rid] = float(np.max(np.abs(bad / nu_ref - 1.0)))
    g7_worst = max((r["err_nu"] for r in nc_rows), default=float("nan"))
    item6 = {
        "records": [r for r in RECORD_IDS if recs[r].model.name in NO_CLAMP_MODELS],
        "n_cases": len(nc_rows),
        "candidate_separation": seps,
        "G3_G4_min_separation": min(seps.values()) if seps else float("nan"),
        "G3_G4_bar": 1e-2,
        "G3_G4_informative": bool(seps) and min(seps.values()) >= 1e-2,
        "G7_worst_err_nu": g7_worst,
        "G7_bar": BAR_A2,
        "G7_pass": bool(g7_worst <= BAR_A2),
        "worst_wall_imbalance": max((r["wall_imbalance"] for r in nc_rows),
                                    default=float("nan")),
        "note": ("These records have no yield stress, so the local Bingham number "
                 "B = tau_y/(k gd^n) that controls simpleFoam's outer iteration for the "
                 "Herschel-Bulkley record is identically zero here. G7 is the control for "
                 "the reading that the convergence failure is about YIELD, not about a "
                 "large viscosity; the separations say the A2 verdict on these records "
                 "can discriminate the map the gate rejects."),
    }

    verdict = {"A1_shear_uniformity": {"worst": a1, "bar": BAR_A1, "pass": a1 <= BAR_A1,
                                       "note": "reported only; A1' is the real test"},
               "A1prime_profile": {"worst": a1p, "bar": BAR_A1P, "pass": a1p <= BAR_A1P},
               "A2_nu": {"worst": a2, "bar": BAR_A2, "pass": a2 <= BAR_A2},
               "A3_tau_wall": {"worst": a3, "bar": BAR_A3, "pass": a3 <= BAR_A3,
                               "signs_opposite_on_the_two_walls": signs_ok,
                               "worst_wall_imbalance":
                                   max(r["wall_imbalance"] for r in rows)},
               "A4_newtonian_control": {"worst": a4, "bar": BAR_A2, "pass": a4 <= BAR_A2},
               "A5prime_cauchy": {"worst": a5, "bar": BAR_A5,
                                  "pass": bool(a5 <= BAR_A5 and signs_ok),
                                  "worst_contraction_ratio":
                                      float(np.nanmax([r["contraction"] for r in rows])),
                                  "worst_final_ux_residual":
                                      max(r["final_ux_residual"] for r in rows)}}
    payload = {"rho": RHO, "H": H, "NY": NY, "n_rates": N_RATES, "rows": rows,
               "relax_U": RELAX_U, "end_time": END_TIME, "windows": list(WINDOWS),
               "verdict": verdict, "openfoam": "v2606"}
    (OUT / "couette_anchor.json").write_text(json.dumps(payload, indent=2), encoding="ascii")

    lines = ["# OpenFOAM plane-Couette anchor for the L3 bridge", "",
             f"Generated by `scripts/openfoam_couette_anchor.py`. OpenFOAM v2606 in WSL, "
             f"`simpleFoam`, laminar, gap H = {H:g} m, {NY} cells, rho = {RHO:g} kg/m^3, "
             f"{N_RATES} wall speeds per record.", "",
             "The comparison uses OpenFOAM's own `nu` field and its `wallShearStress` "
             "function object against `record.eta_at`. It does not touch "
             "`biorheo.solver.openfoam.nu_openfoam`, so it is an anchor on the "
             "transcription of `calcNu`, not a restatement of it.", "",
             "## Pre-registered gates (see the script docstring)", "",
             "| gate | what | worst | bar | verdict |", "|---|---|---|---|---|",
             f"| A1 | shear rate uniform across the gap (reported only, see below) | "
             f"{a1:.3e} | {BAR_A1:g} | {'PASS' if a1 <= BAR_A1 else 'FAILED'} |",
             f"| A1' | max \\|Ux - U0 y/H\\| / U0 against the analytic line | {a1p:.3e} | "
             f"{BAR_A1P:g} | {'PASS' if a1p <= BAR_A1P else 'FAILED'} |",
             f"| A2 | OpenFOAM `nu` vs `eta_at/rho` | {a2:.3e} | {BAR_A2:g} | "
             f"{'PASS' if a2 <= BAR_A2 else 'FAILED'} |",
             f"| A3 | wall shear stress vs `eta_at*gd/rho` | {a3:.3e} | {BAR_A3:g} | "
             f"{'PASS' if a3 <= BAR_A3 else 'FAILED'} |",
             f"| A4 | Newtonian control vs hand value | {a4:.3e} | {BAR_A2:g} | "
             f"{'PASS' if a4 <= BAR_A2 else 'FAILED'} |",
             f"| A5' | shear field at {END_TIME // 2} its vs {END_TIME} its (Cauchy) | "
             f"{a5:.3e} | {BAR_A5:g} | {'PASS' if a5 <= BAR_A5 else 'FAILED'} |", "",
             "**A5 as originally written FAILED and was replaced by A5'; A1 was shown by "
             "the code review to be blind to a checkerboard mode and was replaced by A1'. "
             "Both originals are kept, reported, and explained in the script docstring; "
             "no bar was moved.** Read A1', A2, A3 as the correctness gates -- the answer "
             "here is analytic, so they measure the error itself. A5' is a convergence "
             "diagnostic and can pass on a stall, which is why the contraction ratio over "
             "three windows is reported beside it: a ratio near 1 means the iteration is "
             "not contracting. It is reported as `nan` where the first increment is "
             f"already below {ROUNDOFF:g} -- there the ratio is float noise that happens "
             "to land near 1, and quoting it as a stall would be reading a number the "
             f"data cannot carry. Worst defined value here: "
             f"{verdict['A5prime_cauchy']['worst_contraction_ratio']:.3g}. The residual "
             "is not a criterion at all: its worst final value is "
             f"{max(r['final_ux_residual'] for r in rows):.3e}.",
             "", "## Every case", "",
             "`wall imbal` is |tau_upper + tau_lower| / mean|tau| on the SIGNED "
             "x-tractions: in steady plane Couette the two walls carry equal and opposite "
             "traction, so it is a consistency check needing NO reference value -- the "
             "one diagnostic here that a user with no analytic answer could still run. "
             "`clamp early/final` is the fraction of cells sitting on the viscosity cap "
             "at the first window and at the end: setting `eta_max` 10x clear of the "
             "converged viscosity says nothing about the transient, and the transient is "
             "where the low-rate case spends its whole budget.", "",
             "| case | gd (1/s) | U0 (m/s) | Cauchy | contraction | final Ux resid | "
             "wall imbal | clamp early | clamp final | shear spread | nu_OF min | "
             "nu_OF max | nu_expected | rel err | tau_up/rho | tau_lo/rho | rel err |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        ce, cf = r["cap_bound_fraction_early"], r["cap_bound_fraction_final"]
        lines.append(f"| {r['case']} | {r['gammadot']:.4g} | {r['U0']:.4g} | "
                     f"{r['cauchy']:.2e} | {r['contraction']:.2e} | "
                     f"{r['final_ux_residual']:.2e} | {r['wall_imbalance']:.2e} | "
                     f"{'n/a' if ce is None else f'{ce:.2f}'} | "
                     f"{'n/a' if cf is None else f'{cf:.2f}'} | "
                     f"{r['shear_spread']:.2e} | {r['nu_openfoam_min']:.10g} | "
                     f"{r['nu_openfoam_max']:.10g} | "
                     f"{r['nu_expected']:.10g} | {r['err_nu']:.2e} | "
                     f"{r['tau_openfoam']:.10g} | {r['tau_lower_wall']:.10g} | "
                     f"{r['err_tau']:.2e} |")
    lines += ["",
              "The Newtonian control is hand-written, not an export: it fixes the "
              "`wallShearStress` convention (a factor of 2 or a sign there would make A3 "
              "self-consistently wrong on every non-Newtonian case and A4 is the only "
              "row where the answer is known without OpenFOAM).",
              "",
              ] + diagnosis_section() + [""] + item6_section(item6, nc_rows) + [""]
    (OUT / "couette_anchor.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines[8:]))
    return 0 if all(v.get("pass") for v in verdict.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
