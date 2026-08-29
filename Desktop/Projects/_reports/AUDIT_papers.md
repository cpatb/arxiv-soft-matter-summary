# Audit of Three Standalone LaTeX Physics Papers

Date: 2026-06-30. Auditor: automated LaTeX + physics review.
Scope: compile-check, safe fixes only, physics sanity-check. Scientific claims were **flagged, not altered**.

Compiler: `pdflatex` (MiKTeX). No `\bibliography`/BibTeX is used by any file; all three use inline `thebibliography`, so the `.bib` files present in the tree (including `friday_experimentsNotes.bib`) are **not** consumed at compile time.

---

## Paper 1 — Drop Impact of Dense TPM vs. Cornstarch Suspensions
Files: `paper1/paper1_research_approach.tex`, `paper1/paper1_referee_review.tex`, `paper1/paper1_referee_review.sty`

### (a) Summary
A detailed research-approach / project-plan document (not yet a data paper) proposing a controlled comparison of drop impact between non-porous TPM colloidal spheres (~830 nm) and porous cornstarch at *matched effective volume fraction* and matched bulk rheology. Central claim: jamming/unjamming dynamics are governed by particle identity (size, porosity, surface chemistry) and cannot be predicted from bulk shear-thickening rheology alone. It introduces an "unjamming Péclet number" Pe_unjam = τ_B/τ_KWW, proposes six figures, a rheology + high-speed-imaging protocol, a NaCl ionic-strength control, and anticipates reviewer objections. A companion referee-style review is included.

### (b) Compile status
- `paper1_research_approach.tex`: **compiles cleanly** (14 pp), no undefined refs/citations, no fatal errors, no missing figures (document contains no `\includegraphics`; all figures are described, not embedded).
- `paper1_referee_review.tex`: **compiles cleanly** (7 pp).
- Note: `paper1_referee_review.sty` is **not a style package**: it is a byte-for-byte duplicate of the full `paper1_referee_review.tex` article (with `\documentclass`/`\begin{document}`). It is never `\input`/`\usepackage`d by anything and is inert. Harmless but misleading; consider deleting or converting to a real `.sty`.

### (c) Issues (severity-ranked)
1. **[Low/None — compile]** No LaTeX errors. Uncited bibitem `Brown2014` sits in the reference list but is never `\cite`d (cosmetic).
2. **[Medium — internal consistency]** Threshold ε_jam inconsistency introduced by revision history. The body now sets the **primary** jamming threshold ε_jam = 0.07 and declares values below it (0.05, and the sensitivity band [0.05, 0.10]) as "at/below the measurement floor." But **Objection 5's response** (§ "Anticipated Reviewer Objections") still says "We test three threshold definitions (ε = 0.02, 0.03, 0.05)" — two of which the paper itself calls unresolvable. Residual from the ε 0.03→0.07 bump. Flagged, not fixed (methodological claim).
3. **[Medium — consistency between the two documents]** The referee review critiques an **older draft**: it discusses "colloidal silica" (the approach doc is v4 "TPM replaces silica") and ε_jam = 0.03 (now 0.07). Several review points (salt screening "optional," Pe framing, silica-mechanism caution) have already been addressed in the current approach doc. Expected given the version history, but the two files are out of sync.
4. **[Author-judgment — overstated claims, already self-flagged]** "Orders of magnitude more slowly" for TPM relaxation: the paper's own predicted ranges give τ_KWW^TPM = 300–600 ms vs cornstarch 50–200 ms, i.e. ~2–10×, not "orders of magnitude." The referee review flags exactly this. Likewise "order of magnitude earlier" t_jam. Not altered.

### (d) Fixes applied
None required. The document compiles cleanly and its key equations are dimensionally and numerically sound (verified below).

**Physics spot-checks (all consistent):** Pe_impact (TPM ≈ 220, cornstarch ≈ 3.8×10⁵), D_B via Stokes–Einstein (TPM 5.3×10⁻¹³, cornstarch 4.4×10⁻¹⁴ m²/s), τ_B (325 ms / 570 s), Pe_unjam (~0.6–0.8 / ~3000–10000), porosity correction φ_eff ≈ φ_nom/(1−0.31) ≈ 1.45 φ_nom, Debye lengths (~9.6 nm at 1 mM, ~3 nm at 10 mM), and τ_ic = √(ρD₀³/8σ) (dimensionally seconds) all check out.

### (e) Recommendations NOT applied (need author judgment)
- Reconcile Objection-5 ε list (0.02/0.03/0.05) with the primary ε = 0.07 and sensitivity band [0.05, 0.10].
- Regenerate the referee review against the current TPM (v4) draft, or add a note that it reviews an earlier silica version.
- Soften "orders of magnitude"/"order of magnitude earlier" to match the quoted prediction ranges (as the referee recommends).
- Replace the placeholder reference `Sacanna2013` with the actual TPM synthesis citation before submission.
- Either delete `paper1_referee_review.sty` or make it a genuine shared style file (currently a stray duplicate).

---

## Paper 2 — "Friday Experiments" on Shear-Thickening Free-Surface Flows
File: `friday_experiments.tex` (RevTeX 4-2; `friday_experimentsNotes.bib` present but unused)

### (a) Summary
Proposes five preliminary experiments on dense cornstarch ("oobleck") free-surface flows, motivated by Darbois Texier et al. (inertia-free oobleck-wave instability) and Bougouin et al. (the constant-speed "liquid dam"). The five: (1) 3D radial dam break, (2) elastic-membrane upward "reverse impact," (3) controlled-perturbation mapping of the oobleck-wave dispersion relation, (4) stairwave wavelength selection vs. tilt rate, (5) liquid-dam interaction with a step discontinuity (+ indirect measurement of overflow thickness δ). Derives scalings and measurement strategy for each. A top-of-file changelog documents self-corrections in a prior revision.

### (b) Compile status
**Compiles cleanly** (19 pp), no undefined refs/citations, no fatal errors, no missing figures (none embedded). BibTeX not needed (inline `thebibliography`).

### (c) Issues (severity-ranked)
1. **[HIGH — physics, flagged not fixed]** **Axisymmetric Newtonian spreading exponent looks wrong.** Experiment 1 repeatedly states the 3D (axisymmetric) constant-volume *viscous* gravity current follows **R ∝ t^{1/2}** (lines ~190, 214, 255). For a constant-**volume** viscous (lubrication) current, integrating the paper's *own* front-speed relation Ṙ ~ ρg H³/(η₀R) with H = V/(πR²) gives R⁸ ∝ t, i.e. **R ∝ t^{1/8}** (Huppert 1982). t^{1/2} is the constant-**flux** axisymmetric result, not the dam-break (fixed-volume) case. This underpins the "transition to Newtonian spreading" argument, so it matters. Not auto-fixed because it is a scientific claim, but the correction is unambiguous and follows from the author's own equation.
2. **[MEDIUM — physics, flagged not fixed]** **2D height-decay exponent.** Lines ~237 and ~816 state the 2D Newtonian case decays as **H ∝ t^{−2/5}**. For a 2D constant-volume viscous current, X ∝ t^{1/5} (which the paper correctly states at line 191) and area conservation gives **H ∝ t^{−1/5}**, not t^{−2/5}. Inconsistent with the paper's own X ∝ t^{1/5}.
3. **[Low]** Uncited bibitem `Rocha2023JFM` present in the reference list but never `\cite`d (cosmetic).
4. **[Low — minor numeric]** Exp 2 quotes v_c ~ 10–100 mm/s; using the stated ranges (τ_c 2–4 Pa, η₀ 0.2–1.5 Pa·s, h₀ = 1 cm) the extremes are ~13–200 mm/s. Order-of-magnitude fine; upper bound slightly understated.

### (d) Fixes applied
None. The two exponent problems (#1, #2) are scientific claims and were flagged rather than edited per the conservative-fix policy. No safe typo/units/delimiter fixes were needed.

**Physics spot-checks that PASS:** U_LD = 2k τ_c²/(ρgη₀) is dimensionally a velocity and numerically ~30 mm/s (in the quoted 5–50 mm/s band); the internal geometry is self-consistent (ℓ/H = k/k₀ = 10/2.5 ≈ 4, matching the quoted ℓ/H ≈ 3–5); overflow thickness δ ~ (2ℓτ_c/ρg)^{1/2} is dimensionally a length and gives ≈ 8 mm for the stated inputs, matching Eq. (delta_num); v_c = τ_c h₀/η₀ is a velocity. The 1/2-power δ scaling (corrected in the changelog from an earlier 1/3-power) is dimensionally correct.

### (e) Recommendations NOT applied (need author judgment)
- Correct the axisymmetric viscous spreading law to R ∝ t^{1/8} (and the derived H ∝ t^{−4/8}=t^{−1/2}) or clearly state a constant-flux assumption if t^{1/2} is intended.
- Correct the 2D height decay to H ∝ t^{−1/5} for consistency with the stated X ∝ t^{1/5}.
- Re-examine the "H ∝ t^{−2} in 3D liquid dam vs. Newtonian" comparison once the Newtonian exponents are fixed: the *contrast* survives (the liquid-dam H ∝ t^{−2} follows correctly from R ∝ U_LD t), but the specific Newtonian reference exponents it is compared against are wrong.
- Delete or cite `Rocha2023JFM`.

---

## Paper 3 — Hyperlocal Air Quality Mapping in Chicago
File: `chicago_run_bike_monitor.TEX`

### (a) Summary
A solo citizen-science PhD project proposal: build a wearable Raspberry-Pi + Sensirion SEN54 + GPS sensor package, log street-level PM/VOC/T/RH data on running and biking routes, and test four physics hypotheses — (H1) street-canyon recirculation, (H2) Gaussian-plume dispersion from highway line sources, (H3) Stokes-number particle separation, (H4) lake-breeze pollutant dilution — plus publish an open dataset. Includes hardware BOM, build/software instructions, data pipeline, statistical-power table, and a publication/dissemination plan.

### (b) Compile status
**Initially FAILED (fatal, no PDF).** Two independent fatal errors, now fixed; the file **now compiles cleanly** (20 pp), no undefined refs/citations, no missing figures (none embedded).

### (c) Issues (severity-ranked)
1. **[CRITICAL — fatal, FIXED]** `\titleformat{\section}` after-code (line 45) contained `\titlerule[0.8pt]`; the inner `]` prematurely closed titlesec's optional-argument bracket, leaking tokens into the body → `! Missing \begin{document}`. Reproduced and isolated in a minimal example.
2. **[CRITICAL — fatal, FIXED]** `\rowcolor` (first used line ~275) was undefined: `xcolor` was loaded **without** the `table` option (and `colortbl` was not loaded). → `! Undefined control sequence`. This was masked by error #1 until it was fixed.
3. **[Medium — physics, flagged not fixed]** **Highway plume model — σ_y vs σ_z.** H2 models a highway as a *line source* but fits a **crosswind** Gaussian C(y) = A·exp(−y²/2σ_y²) and interprets the fitted "spread" as σ_y. For an (effectively) infinite crosswind line source, concentration does not fall off in the crosswind direction; the perpendicular-transect fall-off with distance from the road is governed by **vertical** dispersion σ_z(x), not σ_y. The functional form and the σ_y interpretation should be reconsidered (or the highway treated as a finite-length source with an explicit along-road integration). Dimensionally the fitted formula is consistent if Q is a per-unit-length line strength.
4. **[Low — arithmetic, FIXED]** Hardware BOM "Total" row said **\$187**, but the 13 listed component costs sum to **\$177** (verified), and the intro states \$165–\$185. Corrected 187 → 177 (now internally consistent and within the stated range).
5. **[Low — cosmetic, not changed]** Header comment says "Compile with: pdflatex chicago_airquality_project.tex", but the file is `chicago_run_bike_monitor.TEX`. Stale filename in a comment only.
6. **[Author-judgment]** Unverifiable factual/spec claims: "five EPA monitoring stations" for Chicago (2.7M/5 = 540k arithmetic checks out, but the station count itself is unverified); "PM2.5 R² > 0.95 vs FEM" for the SEN5x is optimistic and condition-dependent. SEN54 default I²C address 0x69 is correct.

### (d) Fixes applied
- Line 45: wrapped the inner rule as `{\titlerule[0.8pt]}` so its `]` no longer closes titlesec's bracket. (Fatal → resolved.)
- Line ~22: `\usepackage{xcolor}` → `\usepackage[table]{xcolor}` to define `\rowcolor`. (Fatal → resolved.)
- Line 327: BOM total `\$187` → `\$177` (matches sum of listed rows; verified = 177).

**Physics spot-checks that PASS:** τ_p = ρ_p d_p²/18μ (dimensionally seconds), v_s = ρ_p d_p² g/18μ (m/s; ≈ 4.5 mm/s for a 10 µm, 1500 kg/m³ particle, reasonable), St = τ_p U/L (dimensionless). Gaussian-plume formula is dimensionally consistent for a line-source Q.

### (e) Recommendations NOT applied (need author judgment)
- Reframe H2: for a highway (line source) with perpendicular transects, model concentration vs downwind distance via vertical dispersion σ_z(x) (line-source Gaussian), or explicitly treat the road as a finite line and integrate; do not fit/interpret a crosswind σ_y from perpendicular transects.
- Verify the "five EPA stations" figure and temper the SEN5x "R² > 0.95" accuracy claim (cite the specific evaluations, note RH/condition dependence; the doc already flags high-RH overcounting elsewhere, which is good).
- Update the stale compile-instruction filename comment.

---

## Summary of edits made to source files
| File | Line | Change | Reason |
|------|------|--------|--------|
| `chicago_run_bike_monitor.TEX` | 45 | `\titlerule[0.8pt]` → `{\titlerule[0.8pt]}` | Fatal: inner `]` closed titlesec bracket |
| `chicago_run_bike_monitor.TEX` | ~22 | `\usepackage{xcolor}` → `\usepackage[table]{xcolor}` | Fatal: `\rowcolor` undefined |
| `chicago_run_bike_monitor.TEX` | 327 | Total `\$187` → `\$177` | Arithmetic: sum of BOM rows = 177 |

Paper 1 and Paper 2 required no source edits (both compiled cleanly). All substantive physics concerns were **flagged for author judgment**, not silently changed.
