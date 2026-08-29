# Readiness — friday_experiments.tex (five oobleck free-surface experiments)

**Assessed:** 2026-07-22  ·  **Canonical artifact:** `friday_experiments.tex` (RevTeX 4-2, 19 pp, PDF current 2026-07-06)
**One-line status:** A fully-specified five-experiment program whose two physics errors have been corrected since the last audit — the best-prepared, ready-to-run experimental plan in this folder.

## Scores (0–5)

| Stage | Score | Primary blocker to next point |
|---|---|---|
| Preprint / arXiv | 2 | It is an experiment *proposal*, not a results manuscript — posting it advertises unrun experiments with no data behind them. |
| Experimental validation | 5 | Nothing blocks it. Five experiments with derived scalings, named observables, measurement strategies, and falsifiable predictions. Correctness issues resolved. |
| Journal submission | 1 | No data, no results, no target journal. |

**Recommended next stage:** Experimental validation  ·  **Why:** this document is precisely what bar 2 asks for — "here is the one measurement that would confirm or kill this," five times over — and the exponent errors that would have poisoned Experiment 1 are now fixed.

## What exists (evidence)

- **Compiles cleanly** (19 pp, no undefined refs or citations; inline `thebibliography`, no BibTeX needed). Current PDF 2026-07-06, matching the `.tex` timestamp.
- **Properly attributed:** `\author{C.~Burton}`, `\affiliation{Department of Mechanical Engineering, Northwestern University}`. No placeholders.
- **Five experiments, each with derived scalings and a measurement strategy:** (1) 3D radial dam break; (2) elastic-membrane upward "reverse impact"; (3) controlled-perturbation mapping of the oobleck-wave dispersion relation; (4) stairwave wavelength selection vs tilt rate; (5) liquid-dam interaction with a step discontinuity, with an indirect measurement of overflow thickness δ.
- **Physics spot-checks that passed the prior audit:** `U_LD = 2kτ_c²/(ρgη₀)` is dimensionally a velocity and numerically ~30 mm/s (inside the quoted 5–50 mm/s band); internal geometry self-consistent (ℓ/H = k/k₀ = 10/2.5 ≈ 4, matching the quoted ℓ/H ≈ 3–5); `δ ~ (2ℓτ_c/ρg)^{1/2}` is dimensionally a length and gives ≈ 8 mm for the stated inputs; `v_c = τ_c h₀/η₀` is a velocity. The ½-power δ scaling (self-corrected in an earlier revision from a ⅓ power) is dimensionally correct.
- Grounded in real recent literature: Darbois Texier et al. (inertia-free oobleck-wave instability) and Bougouin et al. (the constant-speed "liquid dam").

## Blockers to the recommended next stage (ranked)

**There are no substantive blockers to running these experiments.** The items below are minor.

1. **[LOW]** `\bibitem{Rocha2023JFM}` (line 899) remains in the reference list and is never cited. Cosmetic; cite or drop.
2. **[LOW]** Experiment 2 quotes v_c ~ 10–100 mm/s, but the stated input ranges (τ_c 2–4 Pa, η₀ 0.2–1.5 Pa·s, h₀ = 1 cm) give extremes of ~13–200 mm/s. Order-of-magnitude fine; the upper bound is slightly understated.
3. **[LOW — practical, not documentary]** Cornstarch suspensions are notoriously batch- and humidity-sensitive. The document specifies τ_c and η₀ ranges but no protocol for characterising the specific batch used. Measure the flow curve of the actual suspension on the day, or the quantitative predictions cannot be checked against anything.
4. **[MEDIUM — if a preprint is wanted]** No results exist. To become postable this needs to become a results paper on one or two of the five experiments, not a proposal.

## Carried-forward flags from prior audits (status)

All `AUDIT_papers.md` (2026-06-30) items for Paper 2 re-checked tonight:

- **#1 [HIGH] Axisymmetric constant-volume viscous spreading written as R ∝ t^{1/2} instead of t^{1/8}** — **RESOLVED 2026-07-06.** The file now carries an explicit correction note at lines 190–192 ("fixed-VOLUME axisymmetric current is R ~ t^{1/8} (Huppert 1982); t^{1/2} is the constant-FLUX case"), and the corrected exponent appears at lines 193, 217, 258 with `\cite{Huppert1982}`. The `Huppert1982` bibitem was added (line 924). This was the audit's most serious finding and it is fully fixed at every flagged site.
- **#2 [MEDIUM] 2D height decay written as H ∝ t^{−2/5} instead of t^{−1/5}** — **RESOLVED 2026-07-06.** Now reads t^{−1/5} at lines 241 and 819, consistent with the paper's own X ∝ t^{1/5}.
- **#3 [Low] Uncited `Rocha2023JFM`** — **still open.**
- **#4 [Low] Exp 2 v_c upper bound slightly understated** — **still open** (not re-derived tonight; carried forward from the audit).
- Audit recommendation to "re-examine the H ∝ t^{−2} liquid-dam vs Newtonian comparison once the exponents are fixed" — **addressed**: line 819 now contrasts H ~ t^{−1/5} (Newtonian 2D) against the liquid dam's H ~ t^{−2}, which is the corrected comparison.

**This is the only project in the sweep where a previously-flagged physics error has been fully corrected.** The correction is documented in-file with a dated comment and the correct citation, which is exactly the right way to do it.

## Concrete next actions (shortest path to experimental validation)

1. **Run Experiment 1 first — it is the cleanest falsifier and the correction makes it sharp.**
   - **Observable:** front radius R(t) of a fixed-volume radial dam break, on log–log axes.
   - **Apparatus:** a cylindrical gate release onto a level plate, overhead camera, ruled backdrop; plus a rheometer flow curve of the same batch.
   - **Parameter regime:** the τ_c 2–4 Pa, η₀ 0.2–1.5 Pa·s range the document specifies; fixed release volume; a Newtonian control liquid (glycerol–water) at matched η₀.
   - **Expected signal:** the shear-thickening suspension should depart from the Newtonian fixed-volume law; the Newtonian control must reproduce **R ∝ t^{1/8}** (Huppert 1982) — and *this is the check that the corrected exponent bought you*. Against the old t^{1/2} it would have looked like a spectacular anomaly in the control.
   - **The control that would falsify the framing:** if the Newtonian control does not give t^{1/8}, the current is not in the fixed-volume viscous-gravity regime (likely inertial or surface-tension-affected at early times), and the comparison baseline is invalid before any oobleck result can be interpreted.
2. Measure and record the flow curve of the specific cornstarch batch on the day of each run.
3. Cite or delete `Rocha2023JFM`; nudge the Exp-2 v_c upper bound to ~200 mm/s.
4. **After running:** Experiments 1 and 5 together (radial dam break + liquid dam meeting a step, with the indirect δ measurement) are the natural first results paper — that is the route to a real preprint and a real journal submission, rather than posting the proposal.
