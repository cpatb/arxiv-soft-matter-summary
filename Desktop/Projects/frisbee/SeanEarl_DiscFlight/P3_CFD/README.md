# P3 — OpenFOAM CFD Sweep

13 steady-state RANS runs of atmospheric flow over Sean Earl Field, one
per inflow direction (0°, 15°, ..., 180° at 15° increments).

## Solver stack

- **OpenFOAM** v2406 (or 11) — `simpleFoam` steady incompressible.
- **Turbulence**: k-ω SST.
- **Inlet**: Richards–Hoxey ABL with Yang et al. (2009) corrected wall functions.
- **Mesh**: `snappyHexMesh` on a hex-dominant background.
- **Validation**: AIJ-B cube benchmark (Tominaga 2008) using identical
  mesh and BC stack, on a separate test mesh.

## Files

- `caseTemplate/` — single OpenFOAM case directory, fully wired
  for one inflow direction. Used by `generateCases.py` as the source.
- `generateCases.py` — clones the template for each of 13 directions,
  patches the inlet velocity vector, writes to `cases/dir_NNN/`.
- `runSweep.sh` — drives `Allmesh && Allrun` for each case sequentially
  (or in parallel via `xargs -P`).
- `sample_disc_altitude.py` — post-process: extract `U_mean`, `k`,
  `epsilon` on a 0.5-m grid at 1, 2, 3 m AGL across the field.
- `validate_AIJ.md` — checklist for the cube benchmark.

## Workflow

```bash
# from this directory
foamSource    # source your OpenFOAM env

# (1) put seanEarl.stl from P2 into caseTemplate/constant/triSurface/
cp ../P2_Geometry/seanEarl.stl caseTemplate/constant/triSurface/

# (2) generate the 13 case directories
python generateCases.py --out cases/ \
    --headings 0 15 30 45 60 75 90 105 120 135 150 165 180 \
    --u-ref 5.0 --z-ref 10.0

# (3) mesh and solve each case
bash runSweep.sh cases/

# (4) post-process to extract the disc-altitude flow field
python sample_disc_altitude.py cases/dir_*/ --out flow_fields.npz
```

The final `flow_fields.npz` is the input to P4 (one entry per direction
plus a metadata header).

## Domain and mesh notes

- Domain: 500 m × 500 m × 80 m (per Tominaga 2008 guidelines: 5H upstream/lateral/top, 15H downstream of H = 22 m tallest building, blockage <3%).
- Background cell: 4 m. Three refinement levels in a 250-m-centered box.
- y+: 30–100 (wall functions).
- Total cell count target: 5–8 M.

## Inlet

z0 depends on inflow heading θ_wind. For θ ∈ [60°, 120°] (easterly,
over lake): z0 = 0.001 m. Else: z0 = 0.5 m (urban). Blended in the
transitional bands.

Reference: u_ref = 5 m/s at z_ref = 10 m. Increase / decrease this
uniformly to study wind-speed sensitivity; the disc Reynolds number
is essentially Re-independent in operational range so a single u_ref
sweep gives the dominant signal.
