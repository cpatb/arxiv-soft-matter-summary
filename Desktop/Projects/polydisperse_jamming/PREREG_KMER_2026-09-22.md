# PREREG — `kmer`: singlet / doublet / triplet suspensions of ONE primary particle

Written 2026-09-22 (fourth block), **before any mixture cell ran**. Predictions:
`predict_kmer.py` → `out/predict_kmer.txt` (committed with this file; must not be regenerated).
Theory: `derivations/multiplet_knobs_2026-09-22.md`.

## 0. Why, from the question down

A real dense colloid is one primary particle plus whatever it has aggregated into. The
existing multiplet model (`derivations/multiplet_psd_2026-09-22.md`) maps each k-mer onto an
**effective sphere** (its A1), so aggregation can only act through `x` (packing position) and
the trapped-fluid factor Λ. But a tangent doublet is a **dumbbell**, and this project has
since MEASURED that dumbbells move the window-knob χ (`shapechi`: 0.886 → 0.788 at μ_m = 3,
window ×1.85). And every contact in a one-primary suspension is primary–primary, so the
Derjaguin threshold `F_rep = κR*` is the SAME for every contact while the force per
contact is not. So aggregation of one primary acts on **χ and σ\***, not only on `x`, and
the effective-sphere model cannot see either.

## 1. Construction and grid

`packing.PackingState.doublet_w`, `.triplet_w` (added 2026-09-22): one primary sphere
(diameter 1), CV = 0; a doublet is two TANGENT primaries, a triplet an equilateral triangle
of three. `w` is a SOLID-VOLUME fraction. Isotropic ladder, `PJ_MAX_ITER` = 1e6, N = 1700
GRAINS, tag `kmer`, per-contact dumps `(f_n, k₁, k₂)` on every cell.

&nbsp;&nbsp;&nbsp;&nbsp;(w₂, w₃) ∈ {(0,0), (0.25,0), (0.5,0), (0.75,0), (1,0), (0,0.5), (0,1)} × μ_p ∈ {0, 0.3, 3} × seeds {1,2,3} × p̂ ∈ {1e-3, 3e-5} = **126 cells**.

## 2. Gates

- **K0 — IDENTITY (passed before launch).** The historical path through the new code
  reproduces `shapechi` α = 1 and α = 2 cells BIT-IDENTICALLY on `phi`, `iters`,
  `n_contacts`, `z_nr`, `force` (2/2, `scripts/gate_kmer.sh`).
- **K0′ — the (0,0) composition** is the historical singlet path: its 18 cells must equal the
  `shapechi` α = 1 records on `phi` and `iters` exactly.
- **K2 — THE SCALED COPY.** With the same seed, a w₂ = 1 pack draws the same radii, the same
  positions (uniform in a cell whose edge scales with the solid volume) and the same axes as
  `shapechi` α = 2, scaled by 2^(1/3). Packing is scale-free, so **`phi` must agree to 1e-4**
  per cell. Worse than that means the two constructions are not the same body.
- **K3 — HEALTH.** converged; mass/volume = 1 to 1e-9; zero internal contacts.

## 3. Blind predictions (`out/predict_kmer.txt`)

- **H_mix — ideal specific-volume mixing**, `1/φ_J(w) = (1−w)/φ_J(0) + w/φ_J(1)`, separately for
  φ₀ (μ 0), φ_m (μ 3) and the μ 0.3 branch. Registered bar: at every (μ, rung, w₂ ∈ {0.25,
  0.5, 0.75}) the seed-mean `φ` is within **2·PHI_SCALE = 0.0042** of H_mix. The rival, also
  registered: a **positive excess** (mixed shapes pack denser than the harmonic mean).
- **H_χ — the window opens linearly-ish in w₂**: χ(w₂) at p̂ 3e-5, μ_m 3 = 0.886 / 0.859 /
  0.834 / 0.811 / 0.788 at w₂ 0 / .25 / .5 / .75 / 1 (±0.01), and so the doublet fraction
  that opens the window at `x` is `w*(x)` = 0.87 / 0.54 / 0.34 / 0.15 at x = 0.80 / 0.83 /
  0.85 / 0.87. At μ_m = 0.3 the window widens only 0.092 → 0.138.
- **H_Z — isostatic counting.** μ_p = 3: `Z_nr` = 4 ± 0.25 at **every** composition, doublets
  and triplets included (any rigid body: Z = 2·6/3). μ_p = 0: linear in the doublet NUMBER
  fraction between 5.97 and 7.92 (± 0.30); hypostatic (isostatic would be 6 + 4n₂).
- **H_σ\* — the onset moves at fixed primary and fixed potential.** Every contact is
  primary–primary, so `F_rep` is uniform and `σ*/κ = 1/Q_{1−1/e}(f/(R* p))`, R* = 1/4. With
  fewer external contacts per primary (Z_body/k) the mean force per contact rises:
  **σ\*(w₂ = 1)/σ\*(0) at μ_p 0.3 lies in [0.55, 0.95]** (Love–Weber bracket [0.61, 0.86]).
  ⚠ Isotropic packs, not sheared — the same caveat the CV lead carried; a sheared re-run is
  the follow-up if this passes.
- **H_δ — sharpness.** δ₁₀/₉₀ of `f` at μ_p 0.3: **δ(w₂ = 1) < δ(0)** (the frictional doublet
  records already show force CV 1.3–8.8 against 0.7 for spheres). Reported per composition;
  no curvature predicted.
- **Triplets (w₃)** carry no numerical prediction: their endpoint is unmeasured. Reported.

## 4. Not claimed

Breakage (rigid clusters only), fractal k ≥ 4 aggregates (the Λ lever), and anything
sheared. The CV 0.5 / N = 1000 finite-size lesson applies: every cell's `d_max/L` is printed.
