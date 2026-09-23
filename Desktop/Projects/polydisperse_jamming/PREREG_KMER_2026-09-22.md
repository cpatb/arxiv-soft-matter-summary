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

---

## Amendment 1 — the straight-chain trimer, written 2026-09-22 ~16:40 CDT, before its campaign

`packing.PackingState.triplet_chain` (added after launch; the running campaign uses the
triangle): k = 3 as three collinear tangent primaries, span 3 diameters. It asks whether
GEOMETRY at fixed k matters as much as k does. One smoke cell (w₃ 0.5, μ_p 0.3, p̂ 1e-3,
seed 1) ran to check the construction: converged, mass/volume 1 to 2e-15, 0 internal
contacts, φ = 0.5547 (the triangle's gate cell: 0.5773). Seen, one friction, one seed; it
decides none of the predictions below, which need χ (two frictions).

Grid (tag `kchain`): (w₂, w₃) ∈ {(0, 0.5), (0, 1)} as CHAINS × μ_p {0, 0.3, 3} × seeds
{1,2,3} × p̂ {1e-3, 3e-5} = 36 cells.

**Predictions (ordinal; elongation is the hypothesis):**
- **C1** at w₃ = 1, p̂ 3e-5, μ_m 3: `χ(chain) < χ(dumbbell, pure) = 0.788 < χ(singlet) = 0.886`,
  and `χ(triangle) > χ(chain)`. I.e. the window-knob tracks ELONGATION, not k.
- **C2** rough (μ_p 3) `Z_nr` = 4 ± 0.25 for pure chains (isostatic counting; any rigid body).
- **C3** frictionless chains are hypostatic: `Z_nr` < 10 (a chain, like a dumbbell, has an
  unconstrained spin about its axis: 5 dofs).

## Amendment 2 — written 2026-09-22 ~18:45 CDT AFTER 109 of 126 cells had landed (post-hoc; stricter)

**Found in the per-contact dumps:** in FRICTIONAL clump packs a fraction of contacts carry forces
implying overlaps far beyond the protocol's softness — with the linear law between equal
primaries, `u/r = fn/(E r²)` — up to 0.7–1.15 of the radius (interpenetrated primaries) and
7–14 % of contacts above u/r = 0.05 in pure frictional triplets. Sphere packs: max u/r 1e-3.
Such contacts are not a rigid-body packing, and they carry unphysical forces.

- **K4 (new, post-hoc, stricter):** a cell is excluded if > 0.1 % of its contacts exceed
  u/r = 0.05. It removes **every μ_p = 3 cell with any clusters** and the μ_p = 0.3 triplet
  cells; μ_p = 0 and 0.3 doublet cells pass. So under K4, **H_χ at μ_m = 3 cannot be graded
  (UNDETERMINED)**; H_σ\* and H_δ are graded at μ_p 0.3 and are unaffected.
- **Sensitivity (`kmer_overlap.py` → `out/kmer_overlap.txt`), reported beside it, not a gate:**
  subtracting the double-counted lens volume moves φ by ≤ 0.2 % of the solid for doublets and
  moves χ(w₂) at μ_m 3 by at most −0.0017 — 60× smaller than the doublet effect (≈ −0.10) and in
  the CONSERVATIVE direction (a smaller χ, i.e. a wider window). The registered H_χ numbers
  therefore stand under the correction; they do not stand under K4, and both are printed.
- **What is clean either way:** the realistic-friction window, μ_p 0.3: χ 0.908 → ~0.864 for
  singlets → pure doublets, all cells passing K4.
- ⚠ The same artefact is expected in `shapechi`'s α = 2, μ_p 3 cells (no dumps there). Its χ
  0.788 headline inherits this caveat; the kmer w₂ = 1 cells (near-copies) bound the effect.

**Amendment 2 — diagnosis (added the same evening, `scripts/diag_dt_overlap.sh`,
`scripts/diag_penetration.py`).** One bad cell (μ_p 3, w₂ = 1, p̂ 1e-3, seed 1) re-run:
- **Not a time-step artefact:** DT_SAFETY 0.5 / 0.25 / 0.1 gives 3.66 / 3.54 / 3.43 % of contacts
  beyond u/r = 0.05, max u/r 0.735 / 0.740 / 0.737 (φ 0.4952 / 0.4970 / 0.4990). The default run
  reproduces the campaign cell's 3.66 %.
- **The overlap is geometric, not inferred:** Yade's `penetrationDepth` equals fn/kn to 5.6e-17;
  worst contact 0.735 of a radius; tangential/normal 0.39–0.81, far below μ = 3 (locked, not sliding).
- **It grows as the pressure falls** (5.3 % at p̂ 3e-5 vs 3.7 % at 1e-3), which a softness artefact
  at the final state would not do. Reading (hypothesis): frictional rigid clumps interlock during
  stress-controlled compression and LOCK IN deep penetrations that μ = 3 prevents from relaxing —
  a preparation-history artefact of rigid, soft-contact clumps at very high friction.
- **Consequence:** clump χ is quoted at **μ_m = 0.3** (clean, and the realistic friction);
  μ_m = 3 clump numbers are reported only with the overlap-volume sensitivity beside them.

---

## RESULT (graded 2026-09-22 ~20:00 CDT, `out/grade_kmer.txt`, 126/126 records)

- K3 health 126/126; **K0′ 18/18 bit-identical**; **K2 FAIL** (|Δφ| up to 6.0e-3: a geometric scaled
  copy is not a dynamic one — the time step scales with the member radius); **K4 (post-hoc)**
  excludes 46 cells: every μ_p 3 cluster cell and the μ_p 0.3 triplet cells.
- **H_mix: 12/12 K4-clean cells within 2·PHI_SCALE** (worst −0.0041 at μ 0.3, p̂ 3e-5, w₂ 0.75);
  INCOMPLETE as registered (the μ 3 half is excluded). No systematic excess of either sign.
- **H_χ at μ_m 3: UNDETERMINED** under K4 (the overlap-volume sensitivity moves χ ≤ 0.0017, so
  the registered values would pass). **Clean, realistic friction (μ_p 0.3, p̂ 3e-5):**
  χ = 0.908 / 0.898 / 0.881 / 0.868 / 0.863 at w₂ = 0 → 1 — the window widens ×1.49
  (0.092 → 0.137), on the ideal-mixing curve (prediction 0.908 / 0.896 / 0.885 / 0.873 / 0.862).
- **H_Z PASS** (rough Z = 4 ± 0.25; frictionless Z linear in doublet number fraction, 3/3).
- **H_σ\* PASS:** σ\*(w₂ = 1)/σ\*(0) = **0.783** at μ_p 0.3, p̂ 3e-5 (registered [0.55, 0.95]).
- **H_δ FAIL:** δ(w₂ = 1) − δ(0) = +0.004, t = 0.05 — doublets do not change the activation
  sharpness (the mechanism behind the prediction had already been withdrawn).
- Triplets: μ_p 0 cells clean; μ_p 0.3 and 3 triangle cells K4-excluded — no clean frictional
  triplet number exists yet. Chains: `grade_kchain.py` (Amendment 1).

## RESULT — chains (Amendment 1), graded 2026-09-22 ~20:30 CDT, `out/grade_kchain.txt`, 36/36

- **C1 UNDETERMINED:** every μ_p 3 chain cell fails K4 (2–15 % of contacts beyond u/r 0.05).
- **C2 UNDETERMINED** (same exclusion). **C3 PASS:** frictionless pure chains Z_nr = 8.06 < 10.
- **Clean, realistic friction (μ_p 0.3; every chain cell passes K4, zero deep overlaps):**
  pure chains χ = **0.8638** vs pure dumbbells **0.8631** at p̂ 3e-5 — indistinguishable; at
  p̂ 1e-3 chains are wider (0.8805 vs 0.8978). So at the quoted rung a third collinear primary
  does NOT widen the window beyond a dumbbell's; the elongation ordering is rung-dependent and
  unresolved. What is robust: **any** cluster (dumbbell or chain) at w = 1 widens the μ 0.3 window
  from 0.092 to ~0.137.
- Note: frictional TRIANGLES interpenetrate (K4-excluded even at μ_p 0.3) while chains do not —
  the compact triangle's pockets are what lock in overlap, consistent with the history-locking
  reading of Amendment 2.
