# Readiness — bio_fluids_bridge

**Assessed:** 2026-07-22  ·  **Canonical artifact:** code — `blood_rheology.py`, `run_blood_bridge.py`
**One-line status:** A working, validated demonstration that the suspension closure ports to blood — but it demonstrates transferability, not a new result, and it has no manuscript.

## Scores (0–5)

| Stage | Score | Primary blocker to next point |
|---|---|---|
| Preprint / arXiv | 1 | No manuscript, and the standalone claim ("blood is shear-thinning and a structural closure captures it") is well-established textbook rheology. |
| Experimental validation | 3 | Needs calibration to real blood viscometry — but that is literature digitization, not a new experiment, which caps this bar. |
| Journal submission | 1 | No manuscript; nothing here would survive a rheology referee as a standalone contribution. |

**Recommended next stage:** Experimental validation (calibration against published viscometry)  ·  **Why:** it is the only bar with a concrete, cheap, currently-doable action — but see the honest note below.

**Honest framing:** this project's real value is as *one leg* of the three-fluid transferability story, not as a standalone unit. Scored alone it is weak; scored as a component of the combined closure paper it is essential. The master report recommends the combined route.

## What exists (evidence)

- **Re-verified tonight: `python run_blood_bridge.py` → VALIDATION PASS.** All five gates green: blood viscosity decreasing with shear rate; η₀/η_∞ = 112× (physiological ballpark is tens); aggregation parameter λ ∈ [0,1] and decreasing; silica thickening as a contrast; thixotropic loop area > 0 (genuine hysteresis).
- Two current figures: `fig1_blood_vs_silica.png` (blood thins ~500 → ~4 mPa·s while silica shows a DST jump, from one closure skeleton), `fig2_thixotropic_loop.png` (up/down shear-ramp hysteresis).
- Model is of the standard Moore–Coussot structural/thixotropic lineage: `λ_ss(γ̇)=1/(1+γ̇/γ_c)`, `η(λ)=η_∞(1+(η₀/η_∞−1)λ)`, Casson-style yield offset, transient `dλ/dt = k_build(1−λ) − k_break·γ̇·λ`.
- `requirements.txt` present. README states caveats without hedging.

## Blockers to the recommended next stage (ranked)

1. **[HIGH]** **Parameters are physiological-ballpark, not fitted.** The README says so explicitly ("illustrative not fitted"). No comparison to a measured blood flow curve exists.
2. **[HIGH]** **Novelty problem, not a completeness problem.** Structural thixotropic models of blood are a mature literature (Moore; Coussot; Apostolidis & Beris; the Casson and Quemada laws). Reproducing shear thinning from one is a sanity check, not a finding. The transferable-framework claim only becomes interesting alongside the shear-*thickening* and thrombosis legs.
3. **[MEDIUM]** 0-D constitutive toy by construction: homogeneous shear, single relaxation time, no Fåhræus–Lindqvist confinement, no margination, no cell deformation. All self-flagged, all real limits on what can be claimed.
4. **[LOW]** No tests. Unlike `tier0_wyart_cates_closure` (6 invariant tests) and `bloodstain_bpa` (6 correctness tests), this project has only the run-time gate report — there is no `test_*.py`. For a portfolio piece that is a visible gap.

## Carried-forward flags from prior audits (status)

- No prior formal audit exists (built 2026-07-16). This is the first assessment.
- `OVERNIGHT_BUILDS_2026-07-16.md` claim "validation PASS" — **verified tonight, still true.**
- README next steps (calibrate to real blood viscometry; margination; clotting as jamming; Fåhræus–Lindqvist) — **all still open.** The "clotting as jamming" item was subsequently built out as `thrombosis_jamming/`.

## Concrete next actions (shortest path to experimental validation)

1. **The one measurement (literature route, free):** digitize a published whole-blood flow curve at known haematocrit — Chien et al. is the classical source, with modern LAOS/thixotropy data from the Beris group — and fit `(η₀, η_∞, γ̇_c, τ_y)` plus the kinetics `(k_build, k_break)`.
   - **Observable:** apparent viscosity vs shear rate over γ̇ ≈ 0.01–1000 s⁻¹, plus one up/down ramp for the hysteresis loop.
   - **Expected signal:** a single parameter set should fit both the steady curve and the loop area. If the steady fit cannot also reproduce the measured hysteresis, the single-relaxation-time assumption is the thing that breaks — which is itself worth reporting.
   - **Control that would falsify it:** fit at two haematocrits. The structural closure predicts η₀ and the loop area both rise steeply with Hct while γ̇_c moves weakly. If γ̇_c has to be refitted freely per Hct, the model is descriptive curve-fitting, not a closure.
   - **Wet-lab version (if desired):** cone-and-plate rheometry on fresh anticoagulated blood at 37 °C, with a solvent trap — but the published-data route settles the same question at zero cost and should come first.
2. Add a `test_bio_fluids.py` with invariants matching the tier0 suite (λ bounds, monotonicity, η limits at γ̇→0 and →∞, loop-area sign). It is a half-hour of work and it is what makes the pair look like one library.
3. Fold into the combined closure paper rather than writing this up alone.
