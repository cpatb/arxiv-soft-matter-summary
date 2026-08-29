# Readiness — hydrogels (PEGDA/PEGMEA drying-fracture program)

**Assessed:** 2026-07-22  ·  **Canonical artifacts:** `Documents/hydrogel_drying_roadmap.tex` (program plan); `drying_analysis/` (validated measurement package, 36/36); `coumarin_penetration/coumarin_penetration_analysis.tex` (collaborator analysis note); legacy `align/analyze/annotate/visualize.py` (wrinkle toolkit)
**One-line status:** A fully-built, validated measurement layer with no data in it — the entire program is blocked on acquiring footage and three cheap bench measurements, not on code or theory.

## Scores (0–5)

| Stage | Score | Primary blocker to next point |
|---|---|---|
| Preprint / arXiv | 1 | No results manuscript exists. The roadmap is a plan; the coumarin document is a working note. Nothing here reports a measured result. |
| Experimental validation | 4 | The measurements are specified down to the instrument; they simply have not been run. One is a ~10-minute bench task. |
| Journal submission | 1 | No data, no manuscript, no target journal named. |

**Recommended next stage:** Experimental validation  ·  **Why:** every observable in the roadmap now has validated measurement code waiting for it, and the coumarin question turns on a single 10-minute transmission measurement.

## What exists (evidence)

- **`drying_analysis/` — validated, re-verified tonight: `python -m drying_analysis.validate` → 36/36 checks pass.** Seven modules covering the roadmap observables: `calibration` (px/mm from fiducial), `swelling_kinetics` (D, k_m, Bi), `indentation_modulus` (E, Hertz + Dimitriadis finite-thickness), `curl_radius` (ρ, Kåsa fit), `crack_metrics` (A_p, Λ), `toughness` (Γ pure-shear, l_fc, ℓ), `kymograph` (crack-vs-peel ordering). The suite includes convention cross-checks, a hand-computed reference, and discriminating cases (crack-only vs crack→peel; surface-at-top vs -bottom) so a wrong constant or a cancelling pair of bugs cannot pass silently. A science-code-review pass caught and fixed 5 silent-correctness bugs during the build.
- **`coumarin_penetration/`** — 25 kB LaTeX analysis + PDF + `penetration.py`. Quantitative Beer–Lambert / photon-budget diagnosis of a collaborator's 1:9:90 gel, concluding **dose starvation is not the cause** (incident dose exceeds the photon-budget bound by ~20×), with the single dominant unknown named (ε₃₆₅ of CoumAc) and three cheap experiments specified.
- **`Documents/hydrogel_drying_roadmap.tex`** — authored (Charlie Burton, Driscoll Lab), the camera-first crack/peel/curl program that defines the headline deliverable (a drying-mode map).
- **Legacy wrinkle toolkit** — coherent 4-script napari pipeline, all imports resolve, no hardcoded scale, no hardcoded Shamah numbers (quarantine respected).

## Blockers to the recommended next stage (ranked)

1. **[HIGH]** **No real data anywhere in the tree.** `drying_analysis` takes crack masks, traced arcs, and frame stacks; none exist. Everything downstream — the mode map, t_c, Λ, ρ — waits on footage.
2. **[HIGH]** **The segmentation/tracing step is unautomated and unbuilt.** Producing crack masks from raw footage (threshold, edge-detect, or hand-trace via `annotate.py`) is explicitly the one interactive step `drying_analysis/README.md` says is not covered. Without it there is no path from camera to metric.
3. **[MEDIUM]** ε₃₆₅ for CoumAc is unmeasured and spans two orders of magnitude in the note's own Table (`tab:delta`); its value alone decides whether optics matter at all (δ₁/ₑ ranges 3.4 mm → 17 µm across the plausible band).
4. **[MEDIUM]** **The mode classifier does not exist.** `drying_analysis` produces the ingredients (Λ, ρ, ordering) but the crack/peel/curl label per run is a downstream rule that needs real morphologies to define.
5. **[LOW]** Legacy-toolkit defects S2/S3/S4 from `AUDIT_camera.md` remain open: `count_forks` over-counts junctions (a synthetic "+" reports 5), the `condition/sample/face/…min` folder contract is inferred from path depth and fails silently, and annotate/analyze disagree on 4- vs 8-connectivity. Harmless today; a trap if the wrinkle outputs are ever used quantitatively.

## Carried-forward flags from prior audits (status)

- `AUDIT_camera.md` **S1 (CRITICAL** — the pipeline measured wrinkle morphology, not the roadmap's drying-fracture observables) — **resolved 2026-07-06** via the separate `drying_analysis/` package (option (a): legacy scripts left intact). Verified tonight, 36/36.
- `AUDIT_camera.md` **Recommendation 2** (per-run px/mm calibration) — **resolved**: `drying_analysis.calibration` implements it and every physical output routes through it.
- `AUDIT_camera.md` **Recommendation 3** (side/kymograph path) — **resolved**: `drying_analysis.kymograph`.
- `AUDIT_camera.md` **S2, S3, S4, S5** (legacy scripts) — **still open**, unchanged.
- `NOTES_moc_flags.md` **Flag 1** (MOC §5 prose has `ℓ ≈ l_fc for brittle low-stretch` inverted; also writes the uniaxial neo-Hookean *f* where the pipeline is pure shear) — **still open in the MOC prose**. The code is correct and asserts the right direction; the vault text needs the author's edit.
- `NOTES_moc_flags.md` **Flag 2** (Bi length convention: §5A full thickness vs §6.3 half-thickness, a 2× discrepancy) — **still open**. Pick one convention before any Bi_c comparison is reported.
- Shamah thesis quarantine — **respected**; no Shamah-derived constants anywhere in the code.

## Concrete next actions (shortest path to experimental validation)

1. **Run the one measurement (tonight, ~10 min).** Coumarin note §"Three cheap experiments", Experiment 1:
   - **Observable:** direct optical attenuation, α = −ln(I/I₀)/L, giving δ₁/ₑ **with no assumption about ε₃₆₅**.
   - **Apparatus:** ThorLabs power sensor; 365 nm source at the working intensity; a 2 mm 1:9:90 gel slab and a blank.
   - **Parameter regime:** L = 2 mm, c_CoumAc = 0.25 M, I₀ ≈ 0.1 W/cm².
   - **Expected signal:** if ε₃₆₅ ≲ 20 M⁻¹cm⁻¹, T_back ≳ 10% and optics are irrelevant — the "won't cure through" report must then be a conversion ceiling or a readout artifact. If ε₃₆₅ ≳ 50, T_back ≲ 0.3% and the gel is genuinely surface-limited.
   - **Control that would falsify it:** repeat the reading *during* irradiation. A photobleaching absorber must show transmission climbing as the front self-clears; flat transmission means a non-bleaching absorber floor, which is a different physical diagnosis.
2. Run the note's Experiment 2 (cut a cross-section, image the 2 mm edge under blue/UV — monomer fluoresces, dimer does not). This is the depth probe the project currently lacks and it directly settles bulk-vs-surface.
3. Shoot one drying run (top + side view, fiducial in frame) and push it through `annotate.py` → `drying_analysis`. One real dataset converts a validated-on-synthetic package into a measured result and is the gate for everything else.
4. Fix the two `NOTES_moc_flags.md` prose discrepancies in the vault MOC so Bi and ℓ/l_fc comparisons are unambiguous before any number is quoted.
5. *Only after 1–3:* decide whether the first paper is the drying mode map or the programmed-stiffness-gradient direction. Both are viable; neither is writable yet.
