"""
read_openfoam.py — read the completed v1 buoyantSimpleFoam run and extract the
numbers the degenerate `coldFaceT` probe could not.

WHY THIS EXISTS
    `postProcessing/coldFaceT/0/surfaceFieldValue.dat` area-averages T on the
    *coldWall* patch, which carries a fixedValue (Dirichlet) BC at 291 K. It
    therefore returns 2.91e+02 at every one of its 30 samples, by construction
    (README.md caveat 1). The information is not lost: the full ASCII volume
    fields were written to `2500/` and `3000/`. Nothing in this project read
    them until this script.

METHOD
    Cell centres are computed from `constant/polyMesh` (points + faces + owner
    + neighbour): each cell's vertex set is accumulated from its faces and
    averaged, which is exact for the planar-faced hexes blockMesh produced
    here. Nothing about the cell ordering is assumed -- the centres are then
    lexsorted onto the (ny, nx) grid, and the reconstruction is asserted
    against the uniform grid blockMeshDict specifies before any number is
    reported. Means over a uniform mesh are plain numpy reductions, which is
    the correct area average here.

WHAT IT FOUND (2026-08-09), and why it does not rescue the run
    The v2 probe fix specced in README.md -- "sample T at the cold-aisle
    inlet/interior, not on the wall patch" -- cannot be executed on this mesh.
    Two independent reasons, both structural:

    1. v1 geometry is a CLOSED differentially-heated cavity: four walls, no
       inlet, no outlet, no rack, no containment baffle. There is no
       cold-aisle inlet to sample and no eta to vary, so no phi(eta) exists
       here even in principle.
    2. Sampling the interior instead does not help, because the 80x100 uniform
       mesh does not resolve the wall thermal boundary layers. The first cell
       centre off each wall (15 mm) already sits at the well-mixed core
       temperature: 303.6 K next to the 291 K cold wall and 304.0 K next to
       the 313 K hot wall -- 0.44 K apart across a 22 K driving difference.
       An interior probe on this grid reads the core, not the aisle.

    So this is a toolchain validation, as its own README says, and the number
    it yields is a cavity-mixing diagnostic, NOT the recirculation fraction
    phi(eta) that `../cfd/` measures. It must not be fed into `dc_thermal.py`.
    See the handoff note at the bottom of ../DESIGN.md for what v2 requires.

Usage:  python read_openfoam.py
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

CASE = Path(__file__).resolve().parent
MESH = CASE / "constant" / "polyMesh"

# --- expected geometry, transcribed from system/blockMeshDict (cross-checked below) ---
NX, NY, NZ = 80, 100, 1
X_MAX, Y_MAX = 2.4, 3.0
T_COLD, T_HOT = 291.0, 313.0            # coldWall / hotWall fixedValue [K] (0/T)


# --------------------------------------------------------------------------- parsing
def _top_level_list(path: Path) -> str:
    """Return the body of the first top-level `N ( ... )` block in an OpenFOAM file."""
    txt = path.read_text()
    m = re.search(r"\n(\d+)\s*\n\(", txt)
    if m is None:
        raise ValueError(f"no top-level list found in {path}")
    start, depth, i = m.end(), 1, m.end()
    while depth:
        if txt[i] == "(":
            depth += 1
        elif txt[i] == ")":
            depth -= 1
        i += 1
    return txt[start:i - 1]


def read_internal_field(path: Path) -> np.ndarray:
    """Parse `internalField` from an ASCII volScalarField/volVectorField -> (n,) or (n, 3)."""
    txt = path.read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s*\n?\s*(\d+)\s*\(", txt)
    if m is None:
        raise ValueError(f"{path}: internalField is uniform or unparsable; nothing to sample")
    kind, n = m.group(1), int(m.group(2))
    start, depth, i = m.end(), 1, m.end()
    while depth:
        if txt[i] == "(":
            depth += 1
        elif txt[i] == ")":
            depth -= 1
        i += 1
    body = txt[start:i - 1]

    if kind == "scalar":
        arr = np.fromstring(body, sep=" ")
        assert arr.shape == (n,), f"{path}: parsed {arr.shape}, header says {n}"
    else:
        arr = np.fromstring(body.replace("(", " ").replace(")", " "), sep=" ")
        assert arr.size == 3 * n, f"{path}: parsed {arr.size} components, header says 3*{n}"
        arr = arr.reshape(n, 3)
    assert np.isfinite(arr).all(), f"{path}: non-finite values in internalField"
    return arr


def cell_centres() -> np.ndarray:
    """Cell centroids from polyMesh, as the mean of each cell's unique vertices.

    Exact for planar-faced hexahedra (what blockMesh wrote here). Makes no
    assumption about how blockMesh ordered the cells.
    """
    pts_body = _top_level_list(MESH / "points")
    pts = np.fromstring(pts_body.replace("(", " ").replace(")", " "), sep=" ").reshape(-1, 3)

    faces = [np.fromstring(g, sep=" ", dtype=int)
             for g in re.findall(r"\d+\(([\d\s]+)\)", _top_level_list(MESH / "faces"))]
    owner = np.fromstring(_top_level_list(MESH / "owner"), sep=" ", dtype=int)
    neigh = np.fromstring(_top_level_list(MESH / "neighbour"), sep=" ", dtype=int)
    assert len(faces) == owner.size, "faces/owner length mismatch"

    n_cells = int(owner.max()) + 1
    verts: list[set[int]] = [set() for _ in range(n_cells)]
    for fi, fp in enumerate(faces):
        verts[owner[fi]].update(fp.tolist())
        if fi < neigh.size:
            verts[neigh[fi]].update(fp.tolist())
    cc = np.array([pts[sorted(s)].mean(axis=0) for s in verts])
    assert cc.shape == (n_cells, 3)
    return cc


def build_grid_index(cc: np.ndarray) -> np.ndarray:
    """Return an index array of shape (ny, nx) mapping grid position -> cell id.

    Verifies against blockMeshDict's uniform (NX, NY, NZ) block; asserts on mismatch.
    """
    xs, ys = np.unique(np.round(cc[:, 0], 9)), np.unique(np.round(cc[:, 1], 9))
    assert xs.size == NX and ys.size == NY, (
        f"mesh is not the {NX}x{NY} block blockMeshDict specifies "
        f"(found {xs.size} x-planes, {ys.size} y-planes)"
    )
    np.testing.assert_allclose(xs, (np.arange(NX) + 0.5) * X_MAX / NX, rtol=0, atol=1e-9,
                               err_msg="x cell centres are not the uniform grid expected")
    np.testing.assert_allclose(ys, (np.arange(NY) + 0.5) * Y_MAX / NY, rtol=0, atol=1e-9,
                               err_msg="y cell centres are not the uniform grid expected")
    order = np.lexsort((np.round(cc[:, 0], 9), np.round(cc[:, 1], 9)))   # y major, x minor
    return order.reshape(NY, NX)


# ----------------------------------------------------------------------------- report
def main() -> None:
    cc = cell_centres()
    gid = build_grid_index(cc)
    x = (np.arange(NX) + 0.5) * X_MAX / NX

    times = sorted(int(p.name) for p in CASE.iterdir()
                   if p.is_dir() and p.name.isdigit() and p.name != "0" and (p / "T").exists())
    assert times, "no written time directories with a T field"
    fields = {t: (read_internal_field(CASE / str(t) / "T")[gid],
                  read_internal_field(CASE / str(t) / "U")[gid]) for t in times}
    T, U = fields[times[-1]]

    assert T.min() > 250.0 and T.max() < 350.0, "T outside a physical range [K]"
    assert T_COLD - 1.0 <= T.min() and T.max() <= T_HOT + 1.0, (
        "interior T escaped the wall BCs — a maximum principle violation, not a probe issue"
    )

    print("=" * 74)
    print("OpenFOAM v1 run — volume fields read from disk (the coldFaceT probe could not)")
    print("=" * 74)
    print(f"mesh {NX} x {NY} x {NZ} = {NX*NY*NZ} cells over {X_MAX} x {Y_MAX} m; "
          f"cell size {X_MAX/NX*1000:.0f} mm; times on disk {times}")
    print("cell centres taken from constant/polyMesh and verified uniform: PASS")

    print(f"\n-- interior temperature, t = {times[-1]} --")
    print(f"  min {T.min():.3f}  mean {T.mean():.3f}  max {T.max():.3f} K "
          f"(wall BCs {T_COLD} / {T_HOT})")
    print(f"  |U| max {np.linalg.norm(U, axis=-1).max():.4f} m/s, "
          f"mean {np.linalg.norm(U, axis=-1).mean():.4f} m/s")
    print(f"  vertical stratification: bottom row {T[0].mean():.3f} K, "
          f"top row {T[-1].mean():.3f} K, delta {T[-1].mean()-T[0].mean():+.3f} K")

    print("\n-- horizontal profile at mid-height (y = 1.5 m), the 'aisle' direction --")
    j = NY // 2
    for i in (0, 1, 2, 5, 10, NX // 2, NX - 11, NX - 6, NX - 3, NX - 2, NX - 1):
        print(f"   x = {x[i]:5.3f} m   T = {T[j, i]:8.3f} K   "
              f"theta = {(T[j, i]-T_COLD)/(T_HOT-T_COLD):6.4f}")

    print("\n-- boundary-layer resolution check (why an interior probe cannot help) --")
    print(f"  first cell off coldWall (x = {x[0]*1000:.0f} mm): {T[:, 0].mean():.3f} K "
          f"vs the {T_COLD} K wall  -> {T[:, 0].mean()-T_COLD:+.2f} K")
    print(f"  first cell off hotWall  (x = {(X_MAX-x[-1])*1000:.0f} mm): {T[:, -1].mean():.3f} K "
          f"vs the {T_HOT} K wall  -> {T[:, -1].mean()-T_HOT:+.2f} K")
    print(f"  the two near-wall columns differ by only {abs(T[:, -1].mean()-T[:, 0].mean()):.3f} K "
          f"across a {T_HOT-T_COLD:.0f} K drive")
    print("  => the thermal boundary layers are unresolved; every interior cell reads the")
    print("     well-mixed core. Moving the probe off the wall does NOT recover a gradient.")

    # ---- flow regime: is a laminar steady solver even admissible here? ----
    # Boussinesq-style scales for a differentially-heated cavity of height H.
    # Properties from constant/thermophysicalProperties (const transport, perfectGas)
    # evaluated at the measured core temperature.
    R_SPECIFIC = 8314.462 / 28.96       # J/(kg K), molWeight from thermophysicalProperties
    MU, PR, P_REF = 1.831e-05, 0.705, 101325.0
    T_core = float(T.mean())
    rho = P_REF / (R_SPECIFIC * T_core)
    nu = MU / rho
    alpha = nu / PR
    beta = 1.0 / T_core                                   # ideal gas
    dT_wall = T_HOT - T_COLD
    Ra = 9.81 * beta * dT_wall * Y_MAX ** 3 / (nu * alpha)
    U_buoy = np.sqrt(9.81 * beta * dT_wall * Y_MAX)       # free-fall velocity scale
    U_max = float(np.linalg.norm(U, axis=-1).max())
    U_mean = float(np.linalg.norm(U, axis=-1).mean())

    print("\n-- flow regime (properties at the measured core T) --")
    print(f"  rho {rho:.4f} kg/m^3   nu {nu:.3e}   alpha {alpha:.3e} m^2/s")
    print(f"  Ra_H = g*beta*dT*H^3/(nu*alpha) = {Ra:.3e}   (transition to turbulence ~1e9)")
    print("  constant/turbulenceProperties: simulationType laminar")
    print(f"  free-fall velocity scale sqrt(g*beta*dT*H) = {U_buoy:.3f} m/s")
    print(f"  computed |U|: max {U_max:.3f} = {U_max/U_buoy:.1f}x that scale, "
          f"mean {U_mean:.3f} = {U_mean/U_buoy:.1f}x")
    print("  a laminar model at Ra ~ 1e10 has no steady solution to converge to, and the")
    print("  computed speeds exceed the only velocity scale the problem has. Relaxation")
    print("  tuning and more iterations cannot fix either — see the handoff in ../DESIGN.md.")

    if len(times) >= 2:
        Ta, Ua = fields[times[0]]
        dT, dU = T - Ta, np.linalg.norm(U - Ua, axis=-1)
        print(f"\n-- steadiness: change between t = {times[0]} and t = {times[-1]} --")
        print(f"  max |dT| {np.abs(dT).max():.4f} K   rms |dT| {np.sqrt((dT**2).mean()):.4f} K")
        print(f"  max |dU| {dU.max():.5f} m/s   rms |dU| {np.sqrt((dU**2).mean()):.5f} m/s")
        verdict = "NOT steady" if np.abs(dT).max() > 0.01 else "changes are small"
        print(f"  -> {verdict} over the last {times[-1]-times[0]} SIMPLE iterations")

    log = CASE / "log.buoyantSimpleFoam"
    if log.exists():
        pat = re.compile(r"Solving for (\w+), Initial residual = ([\d.eE+-]+)")
        hist: dict[str, list[float]] = {}
        for line in log.read_text(errors="ignore").splitlines():
            m = pat.search(line)
            if m:
                hist.setdefault(m.group(1), []).append(float(m.group(2)))
        print("\n-- initial residuals from log.buoyantSimpleFoam (last 500 iterations) --")
        for k, v in hist.items():
            tail = np.array(v[-500:])
            print(f"  {k:8s} final {v[-1]:.3e}   last-500 min {tail.min():.3e} "
                  f"max {tail.max():.3e}   spread x{tail.max()/max(tail.min(), 1e-30):.1f}")
        print("  a residual cycling over a fixed band instead of falling is a limit cycle,")
        print("  not convergence — matches README.md caveat 2.")


if __name__ == "__main__":
    main()
