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

## Amendment 2 (2026-09-24, before any cell finished) — stress unit, and an external anchor
- **Unit correction.** The `-s` input is in σ₀ = F_R/(6πa²), not F_R/a² (`Simulation.cpp:446`,
  `SimulationInit.cpp:289`: `target_stress = input/6/M_PI`). This is the Mari et al. (2014) unit. Ratios
  (P1, P4) are unaffected.
- **P5 (external anchor, graded descriptively).** Sources resolved by a citation agent (full text read):
  Mari, Seto, Morris & Denn, J. Rheol. 58, 1693 (2014), doi:10.1122/1.4890747 — onset of thickening
  "σ_ST/(η₀γ̇₀) ≈ 0.3 … roughly independent of the volume fraction" (bidisperse 1:1.4, equal volumes);
  Singh, Mari, Denn & Morris, J. Rheol. 62, 457 (2018), doi:10.1122/1.4999237 — Table I φ_J⁰ = 0.646,
  φ_J^∞ = 0.562; f = exp(−σ̃*/σ̃) with σ̃* = 1.45. The closest of my arms is CV 0.2 (realised 0.184;
  1:1.4 at equal volumes has number CV ≈ 0.17). Anchor, printed not gated: (a) the CV-0.2 onset (first σ
  with η_r > 1.1 × η_r(0.1)) lies in [0.1, 1] σ₀; (b) the CV-0.2 low-stress inversion φ₀ lies within
  0.02 of 0.646; (c) φ_m(μ 0.3) lies between 0.562 (μ = ∞) and 0.646. My λ = 0.02 a; the citation agent
  reports Mari's ERM used λ = 0.05 a — onset stress depends on λ, so (a) is a loose check.

## Amendment 3 (2026-09-24 ~22:05 CDT, before any cell finished)
The five σ ≥ 30 cells launched before Amendment 1 run to strain 6; at the observed rate
(c20_x95_s100: strain 0.8 after 38 min) at least one will hit the 4 h wall cap near strain 5 — after the
graded window [1.5, 4] is complete. Health rule, for every cell: **the window is complete (strain ≥ 4)
and max overlap < 0.05; rc 0 OR rc 124 (wall cap after the window).** A cell capped before strain 4 is
censored. (Replaces "rc 0, strain ≥ 3"; the registered-vs-post-hoc split in the grader becomes this
single rule. Nothing else changes.)

## Event (2026-09-24 ~23:05 CDT) — queue driver killed by the harness (low Windows memory)
Claude Code stopped the backgrounded `wsl.exe` launcher because Windows was critically low on memory
(0.9 GB free of 31; the WSL VM held 14.8 GB, mostly page cache). Inside WSL the nine running LF_DEM cells
SURVIVED (vmIdleTimeout keeps the VM up) and will finish; the xargs driver died, so the unstarted cells
(σ ≤ 3 for most arms) will not start until relaunched. The harness instruction is not to restart it
without the author's OK; `lfdem/relaunch_lfsig.sh` relaunches only cells whose run directory does not
exist. Nothing about the design, gates or data changes; cells that did run are unaffected.

## Amendment 4 (2026-09-24 ~23:50 CDT, before grading) — reading cells whose runner was killed
At ~23:31 every remaining LF_DEM process stopped (all live data files' last write is 23:31; the WSL VM
was not restarted; cause not identified from inside the session). Because the `run_cell.sh` shells died
too, cells that had finished never got `rc=` written to META. Rule, fixed before any cell is graded:
- a cell whose `run.log` ends with LF_DEM's own "Job done ok" is treated as rc 0;
- a cell killed mid-run is treated as the wall-cap case (rc 124): healthy iff its window is complete
  (strain ≥ 4, Amendment 3) and overlap < 0.05. That admits c50_x95_s100 (killed at strain 5.4) and
  censors c20_x95_s10 (2.9), c20_x95_s100 (3.05), c50_x90_s3 (2.45), c50_x95_s10 (2.2).
Ran to completion: c20_x90_s10/s30/s100, c20_x95_s30, c50_x90_s10/s30/s100, c50_x95_s30. Never started
(driver killed at ~23:05): every σ ≤ 3 cell except c50_x90_s3, and c20_x90_s3 … So P1–P5 are INCOMPLETE
until the queue is relaunched with the author's OK.

## Event 2 (2026-09-25 ~07:00 CDT) — relaunch
The author's standing instruction for this session was all-night work (recorded as feedback in the other session:
an all-night request is itself the OK to relaunch after a harness reap). Windows memory had recovered (13.5 GB free,
WSL page cache ~1 GB). `lfdem/relaunch_lfsig.sh` relaunched at P = 8: it re-runs the four cells killed before
strain 4 (their partial directories moved aside to `<label>.killed_<time>`) and every never-started cell, and skips
the nine finished/complete-window cells. Nothing in the design, gates or predictions changes.
