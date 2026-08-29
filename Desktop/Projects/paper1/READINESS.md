# Readiness — paper1 (Drop impact: dense TPM vs cornstarch suspensions)

**Assessed:** 2026-07-22  ·  **Canonical artifact:** `paper1_research_approach.tex` (14 pp, compiles clean); companion `paper1_referee_review.tex` (7 pp)
**One-line status:** A rigorous, fully-specified experimental program — not a paper — carrying one live methodological defect in the very threshold that defines its headline measurement.

## Scores (0–5)

| Stage | Score | Primary blocker to next point |
|---|---|---|
| Preprint / arXiv | 1 | It is explicitly titled "Research Approach" and contains no data. Posting it would publish the experimental design before the experiments, for no gain. |
| Experimental validation | 4 | The ε_jam threshold inconsistency must be reconciled before the protocol is executed — it is a bounded, named, one-evening fix. |
| Journal submission | 1 | No data. The document itself targets PRF with an August 2026 submission, which requires the experiments to have been run. |

**Recommended next stage:** Experimental validation  ·  **Why:** this document *is* the experimental design — matched φ_eff, matched bulk rheology, six defined figures, a full rheology + high-speed-imaging protocol, an ionic-strength control, and pre-empted referee objections.

## What exists (evidence)

- **Compiles cleanly** (14 pp, no undefined refs/citations, no missing figures — all six figures are specified in prose, none embedded). Companion referee review compiles cleanly (7 pp).
- **Physics independently verified** in `AUDIT_papers.md` (2026-06-30), all consistent: Pe_impact (TPM ≈ 220, cornstarch ≈ 3.8×10⁵); D_B via Stokes–Einstein (5.3×10⁻¹³ and 4.4×10⁻¹⁴ m²/s); τ_B (325 ms / 570 s); Pe_unjam (~0.6–0.8 / ~3000–10000); porosity correction φ_eff ≈ φ_nom/(1−0.31) ≈ 1.45 φ_nom; Debye lengths ~9.6 nm at 1 mM and ~3 nm at 10 mM; τ_ic = √(ρD₀³/8σ) dimensionally seconds. Nothing in the numerics is wrong.
- **The design is genuinely complete:** central claim (jamming/unjamming controlled by particle identity, not predictable from bulk rheology), a new dimensionless group (Pe_unjam = τ_B/τ_KWW), matched-φ_eff comparison, NaCl ionic-strength control, and a self-authored adversarial review.
- Target stated in the document: *Physical Review Fluids* (primary), *Soft Matter* (secondary); submission target August 2026.

## Blockers to the recommended next stage (ranked)

1. **[HIGH — verified still open tonight]** **ε_jam inconsistency.** The body sets the primary threshold ε_jam = **0.07** (line 595), justifies it as ≈2 pixels at 50 µm/pixel (line 600), declares values below 0.07 to fall below the measurement floor (line 601), and gives a sensitivity band [0.05, 0.10] (line 624). But Objection 5's response (line 928) still reads "We test three threshold definitions (ε = 0.02, 0.03, 0.05)" — **two of which the paper itself calls unresolvable.** This is a residual from the ε 0.03→0.07 revision. It is not cosmetic: ε_jam defines t_jam, which is the headline observable. Reconcile before running the experiments, not after.
2. **[HIGH — verified still open tonight]** **Overstated magnitude language.** Line 191 "an order of magnitude *sooner*", line 198 "orders of magnitude *more slowly*". The paper's own predicted ranges give τ_KWW^TPM = 300–600 ms vs cornstarch 50–200 ms — a factor of **~2–10×**, not orders of magnitude. The companion referee review flags exactly this. Left unfixed, it sets a target the experiment cannot hit and a referee will catch immediately.
3. **[HIGH — verified still open tonight]** **Placeholder reference.** `\bibitem{Sacanna2013}` (line 1023) reads verbatim: *"[PLACEHOLDER --- replace with the specific TPM synthesis protocol used in your lab. Candidates: Sacanna & Pine (2010), van der Wel et al. (2017), or your own adapted procedure. Verify DOI before submission.]"* The synthesis protocol is the basis of the entire TPM-vs-cornstarch comparison.
4. **[MEDIUM]** **The companion referee review reviews a superseded draft.** It discusses "colloidal silica" (the approach doc is v4, TPM replaces silica) and ε_jam = 0.03 (now 0.07). Several of its points are already addressed. Either regenerate it against the current draft or add a header stating which version it reviews.
5. **[MEDIUM — external dependency, unverified from this tree]** The drop-impact **data** pipeline that would populate this paper lives outside `Desktop/Projects`. A normalization correction reported on 2026-07-22 changes shipped under-view d/d₀ values by a factor of ~2.14, with reprocessing performed. This is **reported, not verified from within this folder** — but if accurate, every previously-quoted spreading number from that pipeline must be re-checked before it enters this manuscript. Confirm before drafting results.
6. **[LOW]** `paper1_referee_review.sty` is not a style package — it is a byte-for-byte duplicate of `paper1_referee_review.tex` (both 13,349 bytes), never `\input` or `\usepackage`d by anything. Inert but actively misleading. Also `\bibitem{Brown2014}` is uncited, and a zero-byte file `paper1_referee_review` (no extension) sits in the folder.

## Carried-forward flags from prior audits (status)

All four `AUDIT_papers.md` (2026-06-30) items for Paper 1 re-checked tonight:

- **ε_jam inconsistency (Medium)** — **still open**, verified at line 928.
- **Referee review out of sync with the v4 TPM draft (Medium)** — **still open**; neither file has been modified since 2026-04-29 / 2026-05-13.
- **"Orders of magnitude" overstatement (author judgment)** — **still open**, verified at lines 191 and 198.
- **`Sacanna2013` placeholder reference** — **still open**, verified at line 1023.
- **`paper1_referee_review.sty` stray duplicate; uncited `Brown2014`** — **both still open.**

Nothing flagged in the June audit has been addressed. The files predate the audit and have not been touched since.

## Concrete next actions (shortest path to experimental validation)

1. **Reconcile ε_jam (one evening, highest scientific value per hour in this folder).** Decide the threshold set once — most likely ε_jam = 0.07 primary with the [0.05, 0.10] sensitivity band — and rewrite Objection 5's response to test *that* set. Running the August experiments against an internally contradictory threshold definition would compromise the headline measurement.
2. Soften "order(s) of magnitude" at lines 191 and 198 to the ~2–10× the paper's own numbers support. This also makes the prediction falsifiable rather than unreachable.
3. Replace the `Sacanna2013` placeholder with the actual synthesis protocol and verify the DOI.
4. Confirm the status of the external under-view normalization correction and whether any number destined for this manuscript is affected.
5. Delete or convert `paper1_referee_review.sty`; remove the zero-byte stray file; cite or drop `Brown2014`.
6. **Then execute.** The one measurement this document is built around:
   - **Observable:** t_jam — the time after impact at which the measured spreading radius r(t) departs from the Newtonian reference r_N(t) by more than ε_jam·R₀ — measured for TPM and cornstarch at **matched φ_eff and matched bulk rheology**.
   - **Apparatus:** high-speed imaging at 50 µm/pixel (the resolution that sets the 0.07 floor), plus rotational rheometry to establish the matched-rheology condition.
   - **Parameter regime:** matched φ_eff via the porosity correction φ_eff ≈ 1.45 φ_nom for cornstarch; NaCl at 1 mM and 10 mM for the ionic-strength control (Debye lengths ~9.6 nm and ~3 nm).
   - **Expected signal:** TPM jams earlier after impact and relaxes ~2–10× more slowly, at matched bulk rheology — i.e. particle identity, not the flow curve, sets the dynamics.
   - **The control that would falsify the central claim:** if TPM and cornstarch at matched φ_eff *and* matched bulk rheology show the same t_jam and the same relaxation within error, the claim that bulk rheology is insufficient is dead. The NaCl series is the second falsifier: screening the electrostatics should shift the TPM result if surface chemistry is genuinely part of "particle identity," and should barely move cornstarch.
