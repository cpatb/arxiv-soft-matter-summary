# Readiness — chicago_run_bike_monitor.TEX (DIY hyperlocal air-quality sensing)

**Assessed:** 2026-07-22  ·  **Canonical artifact:** `chicago_run_bike_monitor.TEX` (20 pp, PDF current 2026-06-30)
**One-line status:** A well-built citizen-science project proposal whose central physics hypothesis is still analysed with the wrong dispersion coordinate.

## Scores (0–5)

| Stage | Score | Primary blocker to next point |
|---|---|---|
| Preprint / arXiv | 1 | A project proposal with hardware BOM and build instructions, not a scientific manuscript. No results. |
| Experimental validation | 3 | The measurement campaign is fully specified, but H2's analysis plan fits the wrong Gaussian coordinate and would produce an uninterpretable σ. |
| Journal submission | 1 | No data, no manuscript, no target journal. |

**Recommended next stage:** Experimental validation  ·  **Why:** the instrument, pipeline, and campaign design are complete enough to execute — but fix H2 first, or one of the four hypotheses generates a meaningless number.

## What exists (evidence)

- **Compiles cleanly** (20 pp, no undefined refs/citations) after two fatal errors were repaired in the 2026-06-30 audit. Current PDF present.
- **Complete build package:** Raspberry Pi + Sensirion SEN54 + GPS hardware BOM, build and software instructions, data pipeline, and a statistical-power table.
- **Four physics hypotheses:** H1 street-canyon recirculation, H2 Gaussian-plume dispersion from highway line sources, H3 Stokes-number particle separation, H4 lake-breeze pollutant dilution. Plus an open-dataset deliverable.
- **A real motivating gap:** ~5 EPA monitoring stations for 2.7 M people, versus ~10 m resolution from a wearable — the arithmetic (2.7M/5 = 540k people per station) checks out even though the station count itself is unverified.
- **Physics spot-checks that passed the prior audit:** `τ_p = ρ_p d_p²/18µ` dimensionally seconds ✓; `v_s = ρ_p d_p² g/18µ` → ≈ 4.5 mm/s for a 10 µm, 1500 kg/m³ particle, which is reasonable ✓; `St = τ_p U/L` dimensionless ✓. The line-source Gaussian-plume formula is dimensionally consistent if Q is a per-unit-length line strength ✓.
- SEN54 default I²C address 0x69 confirmed correct.

## Blockers to the recommended next stage (ranked)

1. **[HIGH — verified still open tonight]** **H2 fits σ_y where σ_z governs.** Lines 190–207 and 713–716 model a highway as a continuous **line source**, then fit a **crosswind** Gaussian `C(y) = A·exp(−y²/2σ_y²)` and interpret the fitted spread as σ_y, extracting σ_y(x) and comparing its power-law growth to Pasquill–Gifford classes. For an effectively infinite crosswind line source, **concentration does not fall off in the crosswind direction at all** — the fall-off along a perpendicular transect with distance from the road is governed by **vertical** dispersion σ_z(x). As written, the pipeline (`scipy.optimize.curve_fit` on the crosswind form, line 715) will return a number, and that number will not be σ_y or σ_z or anything physically interpretable. Either switch to the line-source Gaussian in σ_z(x), or treat the road as a finite-length source and integrate along it explicitly.
2. **[MEDIUM]** **H1, H3, H4 are "measure and see," not falsifiable predictions.** Street-canyon recirculation, Stokes separation, and lake-breeze dilution are each posed as effects to look for rather than as quantitative predictions with a stated magnitude and a way to be wrong. H2 is the only one framed quantitatively — which is why its analysis error matters disproportionately.
3. **[MEDIUM]** **Unverified factual and specification claims.** The "five EPA monitoring stations" figure is uncorroborated, and "PM2.5 R² > 0.95 vs FEM" for the SEN5x is optimistic and strongly condition-dependent. The document already flags high-RH overcounting elsewhere, which is the right instinct — extend it to the headline accuracy claim and cite the specific evaluations.
4. **[MEDIUM — practical]** No co-location calibration step is specified. A low-cost PM sensor that has never been run beside a reference-grade instrument produces a dataset nobody can use quantitatively. This is the single biggest determinant of whether the open dataset — the proposal's main deliverable — has value.
5. **[LOW — verified still open tonight]** Line 6 header comment reads `Compile with: pdflatex chicago_airquality_project.tex`; the file is `chicago_run_bike_monitor.TEX`. Stale filename in a comment only.

## Carried-forward flags from prior audits (status)

All `AUDIT_papers.md` (2026-06-30) items for Paper 3 re-checked tonight:

- **#1 [CRITICAL] `\titlerule[0.8pt]` inside titlesec's optional argument → `! Missing \begin{document}`** — **resolved** (line 45 now `{\titlerule[0.8pt]}`, verified).
- **#2 [CRITICAL] `\rowcolor` undefined (xcolor loaded without `table`)** — **resolved**.
- **#3 [Medium — physics] σ_y vs σ_z in the highway plume model** — **still open**, verified at lines 190–207 and 713–716, unchanged.
- **#4 [Low] BOM total \$187 vs \$177 sum** — **resolved** in the audit.
- **#5 [Low] Stale compile-instruction filename** — **still open**, verified at line 6.
- **#6 [Author-judgment] "Five EPA stations" and SEN5x "R² > 0.95" unverified** — **still open.**

## Concrete next actions (shortest path to experimental validation)

1. **Reframe H2 before collecting any transect data.** Replace the crosswind fit with the line-source Gaussian in vertical dispersion: concentration along a perpendicular transect falls off through σ_z(x), so fit C(x) with σ_z(x) = a·x^b and compare *b* to Pasquill–Gifford σ_z curves — not σ_y. If a crosswind treatment is genuinely wanted, model the highway as a finite line and integrate along it explicitly. Fixing this after the campaign means re-running transects in matched meteorology, which is far more expensive than fixing it now.
2. **Add a co-location calibration protocol** — park the sensor beside one of the EPA reference stations for 24–48 h before and after the campaign, and report the regression. Without this, H2's fitted emission rate Q cannot be compared to NEI emission factors in any meaningful way.
3. **The one measurement (after 1 and 2):**
   - **Observable:** PM2.5 concentration along perpendicular transects from a highway edge out to ~300 m.
   - **Apparatus:** the specified Pi + SEN54 + GPS package, on foot or by bike, with simultaneous wind speed and direction (from the nearest ASOS station at minimum).
   - **Parameter regime:** repeat transects under neutral (Pasquill class D) conditions, at multiple times of day to vary traffic volume; hold wind direction roughly perpendicular to the road.
   - **Expected signal:** an exponential-family fall-off with distance whose fitted σ_z(x) grows as a power law with an exponent matching the PG class-D curve.
   - **The control that would falsify it:** run the identical transect protocol at a location with no line source under the same meteorology. A flat profile there, and a fitted σ_z(x) that matches PG at the highway, together make the case. A fall-off at *both* sites means the signal is background gradient or instrument drift, not plume dispersion.
4. Verify the EPA station count; temper the SEN5x accuracy claim with specific citations and condition caveats.
5. Update the stale filename comment on line 6.
