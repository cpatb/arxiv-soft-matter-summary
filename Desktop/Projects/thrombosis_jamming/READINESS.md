# Readiness — thrombosis_jamming

**Assessed:** 2026-07-22  ·  **Canonical artifact:** code — `thrombosis_jamming.py`, `run_thrombosis.py`
**One-line status:** The sharpest falsifiable prediction in the entire folder, testable against already-published data for free — and it has no manuscript.

## Scores (0–5)

| Stage | Score | Primary blocker to next point |
|---|---|---|
| Preprint / arXiv | 2 | No manuscript. The result is one figure and a signature check; the prior-art positioning that makes it a *claim* has not been written or verified. |
| Experimental validation | 5 | Nothing blocks it. The discriminating test is named in the README, the data source is named, and it requires no wet lab. |
| Journal submission | 1 | No manuscript; parameters illustrative; the prior-art risk flagged in the README is unresolved. |

**Recommended next stage:** Experimental validation  ·  **Why:** "here is the one measurement that would confirm or kill this" is already essentially written, and the zero-cost version can be run this week from published figures.

## What exists (evidence)

- **Re-verified tonight: `python run_thrombosis.py` → SIGNATURE VALIDATION PASS.** VWF activation shear rate γ̇_c = 5000 s⁻¹; predicted occlusion threshold γ̇_crit = 1475 s⁻¹; occlusion occurs only in 1477–39811 s⁻¹; t_occ spans 0.37–103.7 s (physiologically plausible). All three signature checks green: no occlusion below threshold, t_occ monotone decreasing above it, t_occ diverging near it.
- Figure `thrombosis_occlusion.png` current.
- Model is the Wyart–Cates skeleton with microphysics relabelled: `f(γ̇)=1−exp(−(γ̇/γ̇_c)^β)` for VWF-bridged adhesive contacts, migrating jamming point `φ_J(f)=φ_m f+φ₀(1−f)`, divergence `η_eff=η_plasma(1−φ/φ_J)^(−2)`, 1-D throat densification → occlusion time.
- **The distinguishing claim is stated crisply and correctly** in the README: incumbent continuum thrombosis models (Fogelson; Ku/Mehrabadi two-phase; Antaki VWF-field) produce arrest by *deposition kinetics*, giving smooth monotonic speed-up with shear. A jamming mechanism instead predicts a **critical point** — no occlusion below γ̇_crit, and `t_occ ∝ (γ̇ − γ̇_crit)^(−α)` above it. That is a genuine, falsifiable, mechanism-discriminating signature.
- `requirements.txt` present.

## Blockers to the recommended next stage (ranked)

1. **[HIGH — do this before writing any equations]** **Scout the data before committing.** The README flags this itself: published occlusion-time-vs-shear data may be too sparse to fit an exponent. Pull the actual figures from the Ku group (*Sci. Rep.* 2020; *Biomicrofluidics* 2019) and count the points *before* building a paper on a power-law collapse. If there are six points spanning one decade, the exponent claim is not supportable and the paper must be reframed around the threshold alone.
2. **[HIGH]** **Prior-art risk, unresolved.** The README names it: "prior art closer than it looks (must show fewer params + a real critical point, not relabeling)." It also cites a 2026 Wyart–Cates closure for epithelial *tissue* (arXiv 2602.20886) as genre precedent, marked "verify" — **still unverified**. This must be checked before any writing; if a jamming-thrombosis paper already exists the whole project changes shape.
3. **[MEDIUM]** **γ̇_crit is an output of illustrative parameters, not a prediction yet.** For the test to be a genuine prediction rather than a fit, γ̇_c (≈5000 s⁻¹, attributed to Schneider et al. 2007 — marked "verify" in the README, **still unverified**) and the modulus scale must be fixed from *independent* literature, then γ̇_crit computed and compared. If those are tuned to match the occlusion data, the test proves nothing.
4. **[MEDIUM]** No test suite (`test_*.py` absent) — only the run-time signature gate.
5. **[LOW]** Deliberately minimal physics: no shear-*gradient* effects, no platelet contraction, no embolization. Self-flagged and defensible for a Tier-0 toy, but referees will raise all three.

## Carried-forward flags from prior audits (status)

- No prior formal audit exists (built 2026-07-16). This is the first assessment.
- `OVERNIGHT_BUILDS_2026-07-16.md` claim "signature validation PASS" — **verified tonight, still true.**
- README's own two "verify" markers (Schneider 2007 γ̇_c ≈ 5000 s⁻¹; arXiv 2602.20886 tissue-closure precedent) — **both still unverified.** These are the two citations the argument leans on hardest.
- README's three stated risks (prior art; data sparsity; missing physics) — **all still open.**

## Concrete next actions (shortest path to experimental validation)

1. **The one measurement — do the free version first.**
   - **Observable:** occlusion time t_occ vs wall shear rate γ̇_wall in a stenotic channel.
   - **Apparatus:** none needed for the first pass — digitize published t_occ-vs-γ̇ figures (Ku group, as named in the README). The wet-lab version is a microfluidic stenosis with whole blood at controlled wall shear, occlusion detected by flow cessation or pressure rise.
   - **Parameter regime:** γ̇_wall ≈ 500–20,000 s⁻¹, spanning below and well above the predicted γ̇_crit ≈ 1500 s⁻¹.
   - **Expected signal:** replot t_occ against (γ̇ − γ̇_crit) on log–log axes. A jamming mechanism collapses to a straight line with slope −α and shows **no occlusion at all** below γ̇_crit. Deposition kinetics does neither — it predicts occlusion at every shear rate, just slower, with no collapse and no threshold.
   - **The control that would falsify it:** the existence of *any* reproducible occlusion below the predicted threshold kills the jamming picture outright. This is the cleanest falsifier in the folder — an unambiguous, pre-registered "if we see X, the model is wrong."
2. Verify the two flagged citations before anything else is written.
3. Fix γ̇_c and the modulus scale from independent sources, then *predict* γ̇_crit and compare — do not fit it.
4. Add invariant tests mirroring the tier0 suite.
5. Write it up as the falsifiable-prediction leg of the combined closure paper (see master report) — this is the leg that gives the other two their scientific stakes.
