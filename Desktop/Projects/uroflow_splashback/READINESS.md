# Readiness — uroflow_splashback

**Assessed:** 2026-07-22  ·  **Canonical artifact:** code — `uroflow.py`, `run_uroflow.py`; figures `fig1`–`fig3`
**One-line status:** Correct physics reaching conclusions that were already published experimentally — its value is as a calibratable model, not as a finding.

## Scores (0–5)

| Stage | Score | Primary blocker to next point |
|---|---|---|
| Preprint / arXiv | 2 | No manuscript, and all three headline recommendations reproduce the Splash Lab's published experimental advice (APS-DFD 2013). |
| Experimental validation | 4 | The calibration measurement is crisp and needs only a phone camera — it simply has not been done. |
| Journal submission | 1 | No manuscript; the splash index is an uncalibrated composite, not a measured quantity. |

**Recommended next stage:** Experimental validation  ·  **Why:** two cheap measurements would convert an uncalibrated ranking model into a quantitative one, and both are specified in the README already.

## What exists (evidence)

- **Re-verified tonight: `python run_uroflow.py` runs clean.** Exit conditions Q = 21 mL/s, d = 4 mm → v = 1.67 m/s, We = 173, Re = 6818, Oh = 0.0019. Ideal breakup length L_b = 63 cm. Driest configuration found: height 15 cm, distance 12 cm, aim +0°, splash index 0.9. Wettest: height 75 cm, distance 39 cm, aim −38°, path 85 cm (stream broken into drops), splash index 28.4 — a ~30× span.
- **The dispersion relation is genuinely validated:** Weber's viscous long-wave result reduces to Rayleigh's inviscid ω_max = 0.354√(σ/ρr³) at x = 1/√2, matching the classical 0.343 constant **to 2.5%**. That is a real analytic-limit check, not a self-consistency check.
- Dimensional spot-check passes: We = ρv²d/σ with the quoted inputs gives ~170 ✓; Oh = µ/√(ρσd) ≈ 0.002 for water at d = 4 mm ✓; L_b = (v/ω_max)·ln(r/δ₀) has units of length ✓.
- Three current figures. Physics chain is correctly staged: exit conditions → Plateau–Rayleigh breakup → coherence-vs-drops splash regime → projectile trajectory.

## Blockers to the recommended next stage (ranked)

1. **[HIGH]** **The splash index is a physically-motivated composite, not a calibrated observable.** `(coherence penalty) × √We_impact × sin(impact angle)`. The README is honest that it is "monotonic in the right variables and good for ranking, but the absolute number is not a droplet count." Every quantitative statement — including the headline "~30×" — inherits this. Ranking claims survive; magnitude claims do not.
2. **[HIGH]** **L_b from ideal theory over-predicts by ~3–4×.** The model gives 63 cm; the Splash Lab observed breakup at 15–20 cm. The README flags this and marks it on fig1. Since the entire "stand close" recommendation hinges on whether impact happens before breakup, a 3–4× error in L_b is the difference between a quantitative and a qualitative claim. `ln(r/δ₀) ≈ 12` is the free parameter absorbing it.
3. **[MEDIUM]** **Novelty problem.** All three recommendations (stand close, aim shallow, sit) reproduce the Splash Lab's published experimental findings. Agreement is a validation success and a publication liability at the same time.
4. **[MEDIUM]** **The real design lever is not modelled.** Standing water in the target — pool depth vs splash — is, as the README says, "the real urinal-design lever," and it is absent. A model of splashback that omits the receiving pool is missing the dominant term for the actual application.
5. **[LOW]** The meatus is a slit, not a circle; jet asymmetry and stream splitting are lumped into one effective diameter. No test suite.

## Carried-forward flags from prior audits (status)

- No prior formal audit exists (built 2026-07-08). `OVERNIGHT_fluids_2026-07-08.md` is the build record.
- That note's claims — validation to 2.5% against inviscid Rayleigh; ~30× span between driest and wettest configs — **verified tonight, both still true.**
- Its honest caveats (ideal-theory over-prediction of L_b; composite splash index) — **still open, unchanged.**
- Its four named next steps (standing-water target; two-stream/fan asymmetry; calibrate `ln(r/δ₀)` and the splash index against a phone high-speed clip; reuse the bidet Watson wall-jet machinery for the post-impact film) — **all still open.**

## Concrete next actions (shortest path to experimental validation)

1. **The one measurement (phone camera, ~30 minutes):**
   - **Observable:** the actual intact-stream length L_b, measured as the distance from the orifice to the first pinch-off, for a controlled water jet.
   - **Apparatus:** a syringe pump or gravity-fed reservoir through a 3–6 mm orifice; a phone at 240 fps against a ruled backdrop. No specialised equipment.
   - **Parameter regime:** Q = 15–25 mL/s, matching the modelled peak-uroflow range, at 2–3 flow rates.
   - **Expected signal:** L_b well below the ideal 63 cm — the Splash Lab's 15–20 cm is the prior. Fitting the measured L_b back through `L_b = (v/ω_max)·ln(r/δ₀)` calibrates the single free parameter `ln(r/δ₀)`, which is currently set to Weber's textbook 12.
   - **Control that would falsify the model:** measure L_b at two flow rates. The dispersion-relation model requires L_b ∝ v at fixed r once `ln(r/δ₀)` is fixed. If measured L_b is flat or non-monotonic in v, the jet is turbulent-disturbance-dominated rather than Rayleigh-mode-dominated, and the whole coherence framing is the wrong physics for this regime.
2. **Second measurement, same session:** weigh dry-vs-wet absorbent paper placed around the target at each configuration to get a real splash *mass*. That is what calibrates the splash index from a ranking score into a quantity, and it is the difference between scores 2 and 4 on the preprint bar.
3. Add the standing-water target (pool depth as a parameter). This is the modelling gap that matters most for the application.
4. Only after 1–3 consider writing. With calibrated L_b and measured splash mass, the contribution becomes "a predictive, calibrated splash model" rather than "a model that agrees with known advice" — the latter is not publishable, the former is.
