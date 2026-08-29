# Readiness — tier0_wyart_cates_closure

**Assessed:** 2026-07-22  ·  **Canonical artifact:** code — `wyart_cates.py`, `run_sanity_gate.py`, `test_wyart_cates.py`, `diff_calibrate.py`
**One-line status:** Clean, tested, reproducing closure prototype with no manuscript and illustrative-only parameters — one calibration run from being a real result.

## Scores (0–5)

| Stage | Score | Primary blocker to next point |
|---|---|---|
| Preprint / arXiv | 2 | No manuscript exists at all — only a README. The physics is textbook Wyart–Cates, so a write-up needs a claim of its own. |
| Experimental validation | 4 | All parameters are explicitly illustrative; needs one real flow-curve dataset run through the existing autodiff calibrator. |
| Journal submission | 1 | No manuscript, no calibrated numbers, no novelty claim positioned against Wyart & Cates (2014) and the DST literature. |

**Recommended next stage:** Experimental validation  ·  **Why:** the calibration machinery already works on synthetic data with a passing ground-truth recovery; pointing it at measured corn-starch and silica flow curves is a single bounded run.

## What exists (evidence)

- **Re-verified tonight, all green:**
  - `python run_sanity_gate.py` → **GATE 1 (rheology contrast) PASS, GATE 2 (reversibility) PASS, NUMERICS (ODE consistency) PASS.** Silica-like back-bending flow curve with a 7358.7× viscosity span; corn-starch-like continuous CST with 174.5× span and no back-bending. Reversibility: silica retains f = 0.598 after the pulse, corn starch relaxes to f = 0.006. `df/dt` converges to `f_ss` to 4 d.p. (0.5134 vs 0.5134).
  - `python test_wyart_cates.py` → **6/6 invariant tests pass** (flow-curve DST vs CST, f_ss bounds/monotonicity, jamming-fraction endpoints, viscosity divergence at φ_J, retention vs recovery, state convergence).
- **Differentiable calibration** (`diff_calibrate.py`): PyTorch-autodiff twin recovers ground-truth parameters from noisy synthetic viscometry — φ = 0.520 (true 0.520), φ_m = 0.559 (0.560), σ* = 8.13 (8.0), loss gradient ~1e-10 at the optimum. Not re-run tonight (torch startup is not a fast check); graded from the recorded result and the present `diff_calibration.png`.
- Three current figures: `gate1_flow_curves.png`, `gate2_reversibility.png`, `diff_calibration.png`. `requirements.txt` present.
- Model documented in the README and in `wyart_cates.py`'s docstring: `f_ss(σ)=exp(−σ*/σ)`, `φ_J(f)=(1−f)φ₀+fφ_m`, `η_r=(1−φ/φ_J)^(−2)`, first-order memory with an irreversible latch floor.

## Blockers to the recommended next stage (ranked)

1. **[HIGH]** **Every parameter is illustrative.** The README states this plainly: `(φ, φ_m, σ*, σ_yield, f_res, τ)` were chosen to reproduce a qualitative contrast, not fitted. Until they come from measured flow curves, the gates demonstrate that the *model can* produce DST-vs-CST, not that *these materials do*.
2. **[MEDIUM]** No real flow-curve data is present in this folder. (Third-party corn-starch and 800 nm silica rheometry is reported to exist on an external drive; **unverified from within this tree** — locate and vendor a copy, with provenance, before calibrating.)
3. **[MEDIUM]** The reversibility latch (`f_res`) is the load-bearing physical claim for the drop-impact story and has **no independent experimental anchor** — it is currently a modelling choice that produces the desired retention/recovery split by construction.
4. **[LOW]** Open roadmap items from the README are untouched: coupling into the Herschel–Bulkley Basilisk drop sim, and the κ ∝ a² onset-timescale lever as a second discriminator.

## Carried-forward flags from prior audits (status)

- No prior formal audit exists for this project (built 2026-07-16; `OVERNIGHT_BUILDS_2026-07-16.md` records the build). This is the first assessment.
- `OVERNIGHT_BUILDS_2026-07-16.md` claim "All gates PASS" — **verified tonight, still true.**
- README next-step checklist: items 1–2 checked (both gates) — **confirmed**; items 3–5 (calibrate to real flow curves; couple into Basilisk; add κ ∝ a²) — **still open**.

## Concrete next actions (shortest path to experimental validation)

1. **The one measurement:**
   - **Observable:** steady-state stress-controlled flow curve, η(σ) or η(γ̇), for *your own* corn-starch and silica preps at 2–3 volume fractions each spanning the suspected φ_m.
   - **Apparatus:** rotational rheometer, roughened parallel plate or cone-and-plate (roughened surfaces are essential — smooth tools wall-slip badly in dense suspensions and will fake a lower viscosity plateau).
   - **Parameter regime:** σ from below to well above the expected σ* (the prototype's illustrative value is 8 Pa; sweep ~0.1–1000 Pa); φ chosen to bracket φ_m ≈ 0.56.
   - **Expected signal:** silica in the shear-jamming window shows S-shaped/back-bending η(σ); corn starch below φ_m shows continuous thickening with a yield offset. The calibrator should then return φ_m and σ* consistently across φ for one material.
   - **Control that would falsify it:** fit each φ *independently*. Wyart–Cates requires a **single** (φ_m, σ*) pair to describe the whole φ series for one material. If the recovered σ* drifts systematically with φ, the closure is wrong for this system, not merely mis-tuned.
2. Feed those curves to `diff_calibrate.py` (it already does exactly this on synthetic data) and report recovered parameters **with uncertainties**, not point estimates.
3. Add a provenance block (instrument, geometry, gap, temperature, date, operator) to whatever data is vendored in — without it the calibration is not reproducible.
4. *For a preprint, later:* this project alone is thin — see the master report's cross-cutting note recommending it be written up **together** with `bio_fluids_bridge` and `thrombosis_jamming` as one transferability paper rather than three notes.
