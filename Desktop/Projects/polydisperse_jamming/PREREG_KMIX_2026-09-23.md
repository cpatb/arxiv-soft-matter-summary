# PREREG — `kmu` (A1) and `kmix` (A2): a clean saturated-friction window, and mixed populations

Written 2026-09-22, committed 21:34 CDT (git commit time), before either campaign ran. Plan: `OVERNIGHT_PLAN_2026-09-23.md`.
Blind numbers: `predict_kmix.py` → `out/predict_kmix.txt` (committed with this file).

## A1 `kmu` — does frictional clump interpenetration go away below μ = 3?
Isotropic ladder, N = 1700 grains, `PJ_MAX_ITER` 1e6, dumps on: w₂ ∈ {0, 0.5, 1} (tangent doublets)
× μ_p ∈ {0.6, 1.0} × seeds {1,2,3} × p̂ {1e-3, 3e-5} = 36 cells, tag `kmu`.
- **U1 (the question):** the K4 deep-overlap fraction (contacts beyond u/r 0.05) at μ_p 1.0, w₂ = 1.
  Two registered outcomes: CLEAN (≤ 0.1 %, K4 passes) → report χ(w₂) at μ_m = 1 as the clean
  saturated-friction window; DIRTY → interpenetration is a property of frictional clumps at any
  high μ, and the clean window stays at μ_m 0.3. No direction is predicted.
- **U2 (if clean):** χ_{μ1}(w₂ = 1)/χ_{μ1}(0) = **0.889 ± 0.02** (the μ_m 3 ratio 0.788/0.886, on
  the premise that φ_m saturates by μ = 1).

## A2 `kmix` — do mixed singlet/dumbbell/chain populations follow ideal mixing of the pure ends?
Compositions (w₂ dumbbells, w₃ straight chains; w₁ = rest) ∈ {(0.3, 0.2), (0.2, 0.4), (0.5, 0.5)} ×
μ_p {0, 0.3} × seeds {1,2,3} × p̂ {1e-3, 3e-5} = 36 cells, tag `kmix`, dumps on, K4 applied.
- **M1:** φ within **2·PHI_SCALE = 0.0042** of `1/φ = Σ w_k/φ_k` (endpoints: kmer singlets and
  dumbbells, kchain chains, same μ and rung) at every (composition, μ, rung) — values in
  `out/predict_kmix.txt` (e.g. p̂ 3e-5, μ 0.3: 0.5448 / 0.5296 / 0.5124).
- **M2:** χ(μ 0.3/μ 0) at p̂ 3e-5 within ±0.01 of 0.8847 / 0.8801 / 0.8635.
- The rival for both: a non-additive excess (shape mixing). A systematic sign is reported either way.

## RESULT A1 `kmu` (graded 2026-09-22, committed 23:06 CDT (git commit time), 36/36, `grade_kmix.py` → `out/grade_kmix.txt`)

**U1: DIRTY**, as registered. At μ_p 1.0, w₂ 1 the deep-overlap fraction is 1.84 % (p̂ 3e-5) and
1.28 % (1e-3), against the 0.1 % K4 gate. The interpenetration also starts EARLIER than the
question assumed: μ 0.6 pure dumbbells fail too (0.70 % / 0.51 %), and so does μ 1.0 at w₂ 0.5
(0.42 % / 0.32 %). Singlets are clean at both frictions. Clean cells: singlets at μ 0.6 and 1.0, and
μ 0.6 w₂ 0.5. So frictional tangent clumps interpenetrate from μ_p ~0.6 upward, growing with
friction and with clump content, and **the clean window stays at μ_m 0.3**, as registered for
this outcome. U2 is not graded (conditional on CLEAN).

The one clean frictional mixture (reported, not graded): χ(μ_m 0.6, w₂ 0.5) = 0.8558 vs 0.8930 for
singlets at p̂ 3e-5 (0.8732 vs 0.9021 at 1e-3), i.e. the window widens 0.107 → 0.144. At μ_m 0.3 the
same composition goes 0.092 → 0.119. Whether higher friction amplifies the doublet effect is not
resolved at n = 3 and is not claimed.

## RESULT A2 `kmix` (graded 2026-09-22, committed 23:35 CDT (git commit time), 36/36, no cell K4-excluded)

**M1 FAIL: 5/12 within 2·PHI_SCALE.** The deviation from ideal specific-volume mixing is
**positive in 12/12** (+0.0016 to +0.0069, median +0.0048): mixed singlet/dumbbell/chain populations
pack DENSER than their pure ends predict. This is the registered rival (a non-additive excess),
with a systematic positive sign.

**Where the excess comes from.** The 50:50 singlet:chain binary already on disk from `kchain`
shows the same excess, +0.0050 to +0.0082 at all four (μ, rung), seed-paired sd ≤ 0.002, K4-clean.
The singlet:dumbbell binaries of `kmer` were ideal within the bar. So the excess enters with the
straight CHAINS (the endpoint that packs worst, φ 0.486 against 0.575 for singlets at μ 0.3 p̂ 3e-5),
consistent with compact grains filling the voids chains leave. ⚠ That binary existed before M1 was
registered and was not checked; had it been, M1 would have been registered with the excess.
Protocol confounds checked and excluded: every endpoint and mixture cell ran at the same cap (1e6)
and converged.

**M2 PASS: 3/3.** χ(μ 0.3)/χ(μ 0) = 0.8849 / 0.8729 / 0.8662 against 0.8847 / 0.8801 / 0.8635
(±0.01). The excess is present on both friction branches and largely cancels in the ratio, so
**the window prediction from pure endpoints holds even where the packing prediction fails.** For
the question that is the quantity that matters: χ, not φ, decides whether a drop can jam.
