# PREREG_ROLLCHI_2026-09-25 — is non-sphericity "blocked rolling"? and where does the window open along the aspect ratio?

Written 2026-09-25 ~19:15 CDT, **before** any graded cell ran. Runner: `phase2_shape.py`
(new `--eta-rolls` axis, key unchanged when absent). Launcher: `scripts/launch_rollchi.sh`.
Grader: `grade_rollchi.py` (written after this file, before the data is read).

## 0. Why, and what the author asked

The author asked (2026-09-25) for **a geometric argument that accounts for shape as well as size**:
why polydisperse silica spheres fail to jam solidly where monodisperse ones do, and whether
**filtered, size-matched cornstarch** could show ball-like jamming or whether the **sphericity** of
Stöber silica matters. The project already has: CV is an x-knob (it raises φ₀ and φ_m together,
χ ≈ CV-invariant 0.911/0.912/0.917); friction and "shape" (tangent doublets, α = 2 clumps) are
χ-knobs acting almost only on the FRICTIONAL branch (shapechi: φ₀ −0.9 %, φ_m −11.7 %). What is
missing is (i) WHY shape acts on the frictional branch, stated geometrically and tested; (ii) any
measurement between α = 1 and α = 2, which is where cornstarch-like near-spheres live (only one
α = 1.5 timing probe exists: φ 0.6289 at μ_p 0.3, p̂ 1e-3 — DENSER than spheres, 0.5957).

## 1. The geometric hypothesis under test (H_roll)

Isostatic counting `Z_iso = 2 n_f / n_c` (derivations/geom_shape E1). A frictional sphere contact
constrains 3 of the 6 relative DOF (normal + 2 sliding); ROLLING (2) and TWISTING (1) remain free,
and rolling without sliding is the mechanism by which frictional spheres relax to φ_m ≈ 0.56–0.57
rather than lower. A non-spherical grain cannot roll for free: its contact normals miss its
centroid by an eccentricity `e`, so rotation moves its centroid against its neighbours. **H_roll:
the frictional-branch effect of shape is the blocking of the rolling channel.** It predicts:
  (a) spheres given an explicit contact moment limit (Yade `etaRoll`, rolling AND twisting,
      `M ≤ η r F_n`) reproduce the doublet's frictional φ at realistic friction;
  (b) the effect needs sliding to be blocked first, so it is ZERO at μ_p = 0 and grows with μ_p;
  (c) it saturates once η exceeds the scale at which rolling is no longer the weak link.
Geometric map (heuristic, not derived): a two-sphere clump's normals pass through member centres
at ±(α−1)d_s/2 from the centroid, so the mean lever arm over isotropic normals is
`⟨e⟩/r_eq = (π/4)(α − 1)(d_s/d_eq)` → η_eff = 0.076 / 0.15 / 0.22 / 0.34 / 0.50 / 0.62 at
α = 1.1 / 1.2 / 1.3 / 1.5 / 1.75 / 2.0.

The competing geometric effect of shape is **orientational filling** (Donev 2004,
doi:10.1126/science.1093010, frictionless ellipsoids pack DENSER than spheres, max near α ≈ 1.5):
it acts on both branches. The α = 1.5 probe suggests it is strong at μ_p 0.3.

## 2. ⚠ Disclosure — an unregistered peek

`scripts/rollchi_smoke.sh` (shakeout, OFF-GRID seed 99, `out/rollchi_smoke/`) ran BEFORE this file:
μ 0 η 0 vs η 1: bit-identical φ 0.6432, it 18281 (N 600, p̂ 1e-3). μ 0.3 η 0.1: φ 0.5741 (1e-3),
0.5408 (3e-5); μ 3 η 1: φ 0.4725 (1e-3); z_nr 4.55 / 4.03 / 3.18. So R0 and R3 below were
**seen at seed 99** and are CONFIRMATORY, not blind; R4–R7 and the whole aspect arm are blind
(the seed-99 point touches R5 at four cells, disclosed). Seed 99 is excluded from every gate.

## 3. Grid (N = 1700, lognormal, PJ_MAX_ITER = 1e6, rungs p̂ ∈ {1e-3, 3e-5}, seeds 1, 2, 3)

**Arm R — rolling-resistant spheres (α = 1, CohFrictMat path).**
η ∈ {0, 0.01, 0.03, 0.1, 0.3, 1.0} × μ_p ∈ {0.1, 0.3, 1.0, 3.0} × CV ∈ {0, 0.2, 0.5}; plus μ_p = 0
at η ∈ {0, 1.0} × CV ∈ {0, 0.2, 0.5}. 432 + 36 = 468 cells. tag `rchi`.
**Arm A — two-sphere clumps (FrictMat path, the shape campaign's builder).**
α ∈ {1.1, 1.2, 1.3, 1.5, 1.75} × μ_p ∈ {0, 0.3, 3.0} × CV 0, plus α × μ_p ∈ {0, 0.3} × CV 0.5.
90 + 60 = 150 cells, contact dumps written for gate K4. tag `achi`. α = 1 and α = 2 at CV 0 are
taken from `shapechi` (same N, seeds, rungs, cap); one wiring cell re-runs α 2, μ 0.3, p̂ 1e-3, s 1.

## 4. Gates (health, not physics)

- **G0 wiring**: the achi re-run of shapechi α 2 μ .3 p̂ 1e-3 s1 is bit-identical in (φ, iters, n_contacts).
- **G1 health**: converged at 1e6; a non-converged cell is EXCLUDED and counted, never failed.
- **K4 (arm A)**: a clump cell is excluded if > 0.1 % of inter-grain contacts carry overlap
  u/r = f_n/(E r²) > 0.05 (r the member radius). Same bar as PREREG_KMER Am. 2. μ_p 3 clump cells
  are expected to fail it; they are secondary.

## 5. Registered tests

PHI_SCALE = 2.1e-3 (the protocol's own re-run scatter). Unpaired comparisons are Welch t at the
two-sided 95 % t_crit; seed-paired ones use the paired t. A bar and a t are both printed.

- **R0 (control, two-sided)**: at μ_p = 0, η = 1.0 vs η = 0 is BIT-IDENTICAL in φ and iters on all
  18 (CV, rung, seed) pairs; AND at μ_p = 0.3 η 1.0 vs η 0 moves φ by more than 2·PHI_SCALE at every
  (CV, rung). Either half failing fails R0.
- **R1 (path control)**: CohFrictMat η = 0 vs the FrictMat records (kmer/shapechi, CV 0, μ 0 / 0.3 / 3,
  both rungs) agree within 2·PHI_SCALE (seed means).
- **R2 (saturation at realistic friction)**: μ_p 0.3, CV 0: |φ(η) − φ(1.0)| < PHI_SCALE for η = 0.1
  and 0.3 at both rungs.
- **R3 (H_roll at realistic friction; CONFIRMATORY after the peek)**: μ_p 0.3, CV 0: the saturated
  rolling sphere and the pure doublet (kmer w₂ = 1) agree in φ within 2·PHI_SCALE at BOTH rungs.
  Expected from the peek: pass at 3e-5, marginal at 1e-3.
- **R4 (friction dependence, blind at μ 0.1, 1.0)**: Δ_roll(μ) = φ(η 1.0) − φ(η 0) is monotone
  non-increasing (more negative) in μ ∈ {0, 0.1, 0.3, 1.0, 3.0} at both rungs, CV 0; with Δ_roll(0) = 0.
- **R5 (a geometric Z → φ map, blind except four seed-99 points)**: the mean-field Voronoi
  relation for monodisperse spheres φ = Z/(Z + 2√3) [attributed to Song, Wang & Makse, Nature 2008 —
  UNVERIFIED until lit-verify] with Z = the measured non-rattler coordination z_nr predicts φ within
  0.015 on ≥ 90 % of converged CV-0 sphere cells of arm R (all μ, η, rungs). Printed: the residual
  distribution and whether it trends with μ or η.
- **R6 (size ⟂ constraint separability, blind)**: per CV, one filling factor w(CV) fitted on the μ_p
  = 0 cells only via φ = Z/(Z + 2√3·w) then predicts every frictional/rolling cell at that CV within
  0.015 on ≥ 90 % of cells. Pass = the size distribution acts on the Z → φ map as one number,
  i.e. size and constraints separate.
- **R7 (the size × shape interaction, blind; the discriminator)**: at μ_p 0.3, D(η) = φ(CV 0.5) − φ(CV 0)
  per rung. Ratio ρ = D(1.0)/D(0). Candidates: **clump-like** ρ ≈ 0.62 (shape campaign, lognormal,
  α 2 vs 1, φ_J) vs **separable** ρ ≈ φ_m(η 1)/φ_m(η 0) ≈ 0.94–0.95 (or additive, ρ ≈ 1). Verdict per rung:
  clump-like if the 95 % CI of ρ lies in [0.45, 0.80]; not clump-like if it lies above 0.80;
  otherwise UNDETERMINED. Separation 0.62 vs 0.95 = 0.33; informative only if the CI half-width < 0.16
  (printed).
- **A1 (orientational filling, blind)**: φ₀(α) (μ_p 0) at p̂ 3e-5 exceeds the sphere value at some α
  in {1.1 … 1.75} by more than 2·PHI_SCALE (Donev-like rise). Two-sided report: the α of the max.
- **A2 (where the window opens, blind; the cornstarch question)**: χ(α) = φ(μ .3)/φ(μ 0), seed-matched
  means, at p̂ 3e-5. H_roll + saturation predicts the window opens EARLY: χ(1.2) ≤ χ(1) − ½[χ(1) − χ(2)].
  Rival: late or non-monotone opening (χ(1.2) above that line). Report the full χ(α) curve both rungs.
- **A3 (the eccentricity map, blind)**: χ_clump(α) − χ_roll(η_eff(α)) within ±0.015 at α ≥ 1.2, p̂ 3e-5,
  where χ_roll(η) = φ(μ .3, η)/φ(μ 0, η 0) is log-interpolated in η on arm R (CV 0). Heuristic map;
  a FAIL says the lever-arm picture does not transfer, not that H_roll is false.
- **A4 (size × aspect)**: D(α)/D(1) at μ_p 0.3 and 0 per rung (the shape campaign's 0.62 generalised). Report only.

Every verdict is printed beside the numbers it rests on. A test whose candidates are closer than
their resolving power is reported NOT INFORMATIVE, whatever its nominal outcome.
