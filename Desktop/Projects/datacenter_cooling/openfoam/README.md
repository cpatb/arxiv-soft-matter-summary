# Tier 2 (higher fidelity) — OpenFOAM aisle case

**v1 status: TOOLCHAIN VALIDATED on OpenFOAM v2606 (2026-07-25).** A closed
differentially-heated 2D aisle box (cold-aisle wall 18 C, hot-aisle wall 40 C) that
exercises `buoyantSimpleFoam` end to end. `blockMesh` -> `buoyantSimpleFoam` (3000 SIMPLE
iters, ~70 s) -> `surfaceFieldValue` post-processing all run. This proves the v2606
toolchain works on the machine; v2 turns it into a physical result.

## Two dict fixes needed on first blind run (v2606, recorded)
1. **`system/fvSolution` SIMPLE:** the box is fully enclosed (all walls, no pressure-fixed
   outlet), so pressure floats -> `FATAL: Unable to set reference cell for field p`. Fixed
   by adding `pRefCell 0; pRefValue 0;`.
2. **`system/controlDict` `coldFaceT`:** v2606's `surfaceFieldValue` requires an explicit
   `writeFields` entry -> `FATAL: Entry 'writeFields' not found`. Fixed by adding
   `writeFields false;`.

## Two honest caveats before trusting the physics (v2 to-do, NOT toolchain blockers)
1. **The proxy is degenerate.** It reads 291.0 K = exactly the 18 C cold-wall BC, because
   area-averaging T on a fixed-temperature (Dirichlet) wall just returns that wall's set
   value. To measure recirculation, sample T at the cold-aisle inlet/interior, not on the
   wall patch.
2. **Not steady-converged.** The run hit the 3000-iteration cap with `p_rgh` limit-cycling
   at ~4e-2 initial residual (U, h ~1e-4; 0 convergence flags) — the signature of a mildly
   unsteady buoyant cavity. Needs relaxation tuning / a transient check / more iters before
   any number here is a physical result.

The validated Tier-1 recirculation calibration remains the Python `cfd/` solver (de Vahl
Davis, φ(η)); this OpenFOAM case is the 3D-upgrade scaffold, per DESIGN.md.

## Run
```
source /usr/lib/openfoam/openfoam2606/etc/bashrc   # load OpenFOAM env
cd openfoam
./Allrun
```
Watch residuals fall in `log.buoyantSimpleFoam`; the cold-wall T history is written to
`postProcessing/coldFaceT/`. Visualize with `paraFoam` (or `touch case.foam && paraview case.foam`).

> First blind run note: OpenFOAM cases rarely run first-try when hand-authored. If `./Allrun`
> errors, paste the tail of `log.blockMesh` or `log.buoyantSimpleFoam` and we'll fix it fast
> (usual suspects: a dict keyword renamed in v2606, or a solver/scheme name).

---

# Tier 2 (higher fidelity) — OpenFOAM aisle/hall scaffold

**Status: scaffold, not yet run.** The Python CFD in `../cfd/` is validated (de Vahl Davis)
and already calibrates the dominant (self-)recirculation. This OpenFOAM case is the upgrade
that adds what 2D cannot: **3D rack-to-rack coupling** (the off-diagonal of the influence
matrix `B`) and a real turbulence model at hall Reynolds numbers.

## Why OpenFOAM here
- Solver: **`buoyantSimpleFoam`** (steady, buoyant, incompressible-ish via Boussinesq or full
  `rhoThermo`). Standard for room-scale HVAC / data-hall CFD.
- Turbulence: `kOmegaSST` (robust for wall-bounded recirculating flow at Re ~ 10^4–10^5).
- This is the open-source stand-in for the commercial DC-CFD tools (6SigmaDCX, Future Facilities).

## Case setup (target)
- Geometry: one cold aisle + a row of racks + one hot aisle; symmetry planes on the row ends
  for a periodic slice, or the full room for a specific layout (Tier 2 with real plans).
- Racks: momentum sources (fans) + volumetric heat sources sized to each rack's kW; inlet on
  the cold-aisle face, outlet on the hot-aisle face.
- CRAH: fixed-flux cold supply inlet (raised-floor tiles) + return outlet at the ceiling.
- Containment: a baffle patch of variable height over the cold aisle; sweep its height to get
  φ(η) in 3D, directly comparable to `../cfd/cfd_recirc.json`.

## Files here
- `system/blockMeshDict` — a parametric 2D-slice mesh matching the Python cell (a concrete
  starting mesh; extend to 3D by adding rack depth in `z`).
- `constant/g` — gravity.
- `Allrun` — mesh → run → post skeleton.

## What still needs adding before it runs (standard dicts)
`system/{controlDict,fvSchemes,fvSolution}`, `constant/{thermophysicalProperties,
turbulenceProperties}`, and `0/{U,T,p,p_rgh,k,omega,nut,alphat}`. Copy from the OpenFOAM
tutorial `heatTransfer/buoyantSimpleFoam/hotRoom` and adapt the rack heat sources
(`fvOptions`: `scalarSemiImplicitSource` for heat, `meanVelocityForce`/`fanPressureJump` for
rack fans).

## Recirculation extraction (the φ that feeds the reduced-order model)
Add a `functionObjects` block to `controlDict`:
- `surfaceFieldValue` on the rack-inlet patch → area-weighted mean `T` → normalize by
  (T_exhaust − T_supply) to get φ, exactly as in the 2D calc.
- Do this per rack to populate the full influence matrix `B_ij` (release a passive scalar in
  rack j's exhaust, integrate its capture at each rack i's inlet).

## Validation before trusting it
1. Reproduce the `../cfd/` 2D φ(η) trend on a 2D slice (sign + shape must match).
2. Grid-independence: refine once, φ must move < ~5%.
3. Then extend to 3D and report the off-diagonal coupling the 2D model cannot see.
