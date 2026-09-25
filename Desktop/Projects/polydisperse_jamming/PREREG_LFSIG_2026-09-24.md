# PREREG — `lfsig`: σ\* and the jamming window from a LUBRICATED simulation, at two CVs

Written 2026-09-24 ~21:50 CDT, before any lfsig cell ran (commit time is the record). Code:
`lfdem/make_conf.py`, `lfdem/params_base.txt`, `lfdem/run_cell.sh`, `lfdem/launch_lfsig.sh`,
`lfdem/queue_lfsig.txt`. Grader `lfdem/grade_lfsig.py` is written before the first cell is read.

## Why
Every σ\* number in this project comes from the force tails of rate-controlled, hydrodynamics-free
Yade packs (`actshape`, `predict_cvmap.txt`): σ\*(CV)/σ\*(0) = 0.957 (CV 0.2), 0.735 (CV 0.5) at fixed
number-mean size, Derjaguin threshold `F_rep = κR*`. σ\* is the one axis the programme names as
unmeasured (`PROGRAMME_2026-09-22.md` §2; critique of 2026-09-24 §2). LF_DEM (Seto, Mari, Morris;
github.com/ryseto/LF_DEM commit 00fa7cc0, built 2026-09-24 against conda-forge SuiteSparse 5.3 CHOLMOD)
has lubrication, a Derjaguin-scaled repulsion `F = F_R·a₀a₁/(a₀+a₁)·exp(−h/λ)` (`RepulsiveForce.cpp`)
and Coulomb friction once surfaces touch — i.e. it produces a genuine rheological σ\*, with the SAME
size convention as the Yade prediction. It is therefore an independent test of the one σ\* prediction
the project has.

## Design
- N = 500, 3D periodic Lees–Edwards, lognormal radii from `psd.LogNormalPSD` (project sampler,
  truncation (0.2, 5)·d₁₀), number-mean radius 1 = LF_DEM length unit. Stress unit F_R/a², a = number
  mean radius (the d₁₀ convention). Realised CV is recorded per config (CV 0.2 draws 0.184 at N 500,
  seed 1 — truncation and finite N; the realised value is the one quoted).
- CV ∈ {0.2, 0.5}. CV 0 is excluded: monodisperse spheres order under shear at these φ, which is a
  different state, not a colder version of the same one.
- φ at x_m = φ/φ_m^Yade(CV) ∈ {0.90, 0.95}, with φ_m^Yade(μ 0.3, p̂ 3e-5) = 0.58755 (CV 0.2), 0.63667
  (CV 0.5) from `out/predict_cvmap.txt`: φ = 0.52880, 0.55817 (CV 0.2); 0.57300, 0.60484 (CV 0.5).
- σ ∈ {0.1, 0.3, 1, 3, 10, 30, 100} F_R/a², stress-controlled. μ = 0.3 (the project's μ_p), λ = 0.02 a,
  lub_reduce 1e-3, kn = max(1e4, 1000σ) F_R/a (overlap ≲ 1 %), kt = kn/2, strain 6 per cell, one
  configuration (seed 1) per (CV, φ) shared across σ (paired).
- 28 cells, ≤ 10 concurrent (core budget agreed with the concurrent basilisk session). Wall cap 4 h per
  cell: a capped cell (rc 124) is CENSORED and printed, never imputed.

## Observables (per cell, strain window [2, 6] unless it is shorter; block-averaged, 4 blocks)
- η_r (column 4 of `data_*`), contact number Z_c (col 19), min gap, max overlap.
- Health: rc 0, strain reached ≥ 4, max overlap (−min gap) < 0.05 a. Unhealthy cells are excluded.

## Readouts
- **σ\*_η**: at each (CV, φ), the stress at which ln η_r crosses the midpoint between its σ = 0.1 and
  σ = 100 values, by linear interpolation in (ln σ, ln η). Defined only if ln η_r is monotone
  non-decreasing across the ladder within ±2 block s.e.; otherwise printed as UNDEFINED.
- **WC inversion**: φ_J(σ) from each φ by η_r = (1 − φ/φ_J)^−2; φ₀ ≈ φ_J(0.1), φ_m ≈ φ_J(100) (the
  plateau readings); f(σ) = (φ₀ − φ_J)/(φ₀ − φ_m); χ_LF = φ_m/φ₀. The two φ must give φ₀ and φ_m
  agreeing to 0.01 for the inversion to be used (the prefactor-1 WC form is assumed, and this checks it).

## Predictions and gates
- **P1 (the Yade prediction, primary).** R = σ\*_η(CV 0.5)/σ\*_η(CV 0.2) at matched x_m, both x_m.
  Yade (d₁₀, Derjaguin): R = 0.735/0.957 = **0.768**. PASS if R < 1 at BOTH x_m and |ln R − ln 0.768| <
  ln 1.3 at both. The rival (the d₃₂ convention, σ\* ∝ 1/d₃₂²-like, `predict_cvmap` H2 bracket 1.145 for
  CV 0.5 vs 0) predicts R > 1. Separation of the two candidates: 1.145/0.768 ≈ 1.49 (below the 1.6 bar
  for a sharp discriminator), so a P1 verdict decides the SIDE of 1, and the magnitude only within ×1.3.
- **P2 (window).** χ_LF(CV 0.5) − χ_LF(CV 0.2) has the Yade sign (+0.005) within ±0.02 — i.e. CV does
  not open the window. Registered as |Δχ_LF| < 0.02.
- **P3 (anchor that does not pass through this code).** The low-stress branch must be frictionless-like:
  Z_c(σ 0.1) < 1 at x_m 0.90 (repulsion keeps surfaces apart), and η_r must RISE from σ 0.1 to 100 at
  every (CV, φ) (shear thickening exists). A cell set that does not thicken is a wiring failure, not a
  result.
- **P4 (the geometric closed form, `derivations/sigma_star_geometry_2026-09-24.md`)** — registered in
  that file before any lfsig cell is read; graded against the same R.

## What would change the plan
If > 1/3 of cells are censored by the wall cap, the strain window is shortened to what landed and the
change is written here as an amendment before grading.

## Amendment 1 (2026-09-24, ~22:10 CDT, before any cell finished or was graded)
Cost at launch: σ = 100 cells advance 0.3–0.85 strain per 11 min with the machine at load 24/22
(shared with the basilisk session), so strain 6 would reach the 4 h wall cap and the 28-cell queue
would need > 10 h. Seen before this amendment: the progress lines (strain, η, Z_c, min gap) of the
first 5 cells at strain ≤ 0.85 — transient values, not graded observables.
- Cells that have not started get STRAIN = 4 (override file read at cell start); the 5 running
  σ ≥ 30 cells keep time_end 6 and are graded on the SAME window as everything else.
- Averaging window for EVERY cell: strain [1.5, 4]. Health: strain reached ≥ 3 (was ≥ 4).
- Nothing else changes (grids, predictions, P1–P4 bars).

## Erratum to P4 (2026-09-24, before any cell finished)
`derivations/sigma_star_geometry_2026-09-24.md` §6 wrote the Yade C/q ratio for CV 0.5 vs 0.2 as "1.05",
with "recomputed exactly in the grader". The exact value from `out/sigma_star_geometry.txt`
(μ 0.3, p̂ 3e-4: (1.0701·1.0700)/(1.0157·1.0055)) is **1.1211**; the grader (`lfdem/grade_lfsig.py`,
committed now) parses it from that file. The ×1.2 bar is unchanged.
