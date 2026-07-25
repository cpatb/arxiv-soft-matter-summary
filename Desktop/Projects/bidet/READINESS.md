# Readiness — bidet (nozzle optimization + adoption policy, matched pair)

**Assessed:** 2026-07-22  ·  **Canonical artifacts:** `bidet_optimization_1.tex` (revtex4-2, REVISION v3) and `bidet_paper_1.tex` (revtex4-2, REVISION 30 May 2026). Superseded: `bidet_optimization.tex`, `bidet_paper.tex`, both `*_pre-review.tex`, both `_test_*.tex` (article-class compile-portability copies).
**One-line status:** Two complete, compiling, numerically-verified manuscripts blocked by an author placeholder and three number reconciliations — the closest thing in this folder to "post this week."

## Scores (0–5)

| Stage | Score | Primary blocker to next point |
|---|---|---|
| Preprint / arXiv | 4 | `bidet_optimization_1.tex` still carries `[Author Name]` / `[email@institution.edu]` / `[Institution]`; a 164 pt table overflow bleeds into the margin in the policy paper. |
| Experimental validation | 3 | The key measurement (fecal detachment shear τ_c) is correctly identified as *the* unknown but has no protocol, apparatus, or surrogate specified. |
| Journal submission | 3 (optimization) / 2 (policy) | Novelty not positioned; τ_c unmeasured is the referee-bait; no target journal named for either paper. |

**Recommended next stage:** Preprint / arXiv  ·  **Why:** both manuscripts compile, their headline numbers independently reproduce, and every remaining blocker is named, bounded, and fixable in one evening.

## What exists (evidence)

- **Both canonical files compile cleanly** (`AUDIT.md`, 2026-06-30, two passes each, exit 0). No undefined refs, no undefined cites, no multiply-defined labels. Current PDFs present (440 kB and 416 kB, both 2026-06-30).
- **Numbers independently recomputed and confirmed** in the prior audit: optimization V₀ ≈ 11.48 m/s, d* ≈ 0.80 mm, Q* ≈ 0.347 L/min; policy r_clean ≈ 26 mm (from τ_c = 5 Pa, V₀ = 11.2 m/s, a = 0.65 mm, ν = 1e-6 → 0.0258 m ✓).
- **Dimensional consistency verified** for the core results: `Q* = π(τ_c r_max⁵ ν / c₂ρ)^{1/3}` → m³/s ✓; `a* = (τ_c r_max⁵ ν / c₂ ρ V₀³)^{1/6}` → m ✓; hydraulic-jump scaling `r_h ~ (Q⁵/ν³g)^{1/8}` → m ✓ (this correctly replaced a previously dimensionally-broken Watson gravity formula).
- **Optimization structure is clean and stated:** P1 (min Q at fixed ΔP), P2 (max r_clean at fixed Q), constraints C1 turbulent impingement / C2 jet coherence / C3 pressure ceiling.
- **Honest self-flagging already present:** both papers mark the laminar Region-II asymptotic as a *lower bound* on wall shear and openly note the tension that the standoff-corrected 1.6 L/min laminar prediction exceeds commercial flow rates. No fabricated data was found; environmental and economic claims are cited.
- The two papers are a matched pair and cite each other.

## Blockers to the recommended next stage (ranked)

1. **[HIGH — verified still open tonight]** **Author placeholders.** `bidet_optimization_1.tex` lines 69–71: `\author{[Author Name]}`, `\email{[email@institution.edu]}`, `\affiliation{Department of Chemical and Biomolecular Engineering, [Institution]}`. The policy paper is fully attributed (Charlie Burton, Northwestern); the optimization paper is not. **Note the departmental mismatch** — the placeholder says Chemical and Biomolecular Engineering while `friday_experiments.tex` uses Mechanical Engineering. Pick the right one.
2. **[HIGH — verified still open tonight]** **164 pt table overflow.** `bidet_paper_1.tex` `tab:stresscomparison` (lines ~759–770) overflows the two-column measure by 164 pt — in revtex `reprint` that bleeds well past the column into the gutter. Six overfull boxes remain in the current log. Fix: promote to `table*` full-width, or `\resizebox{\columnwidth}{!}{…}`, or shorten the "Ratio (TP/Bidet)" column. `tab:commercial` in the optimization paper overflows by 51.7 pt — same fix, milder.
3. **[MEDIUM — verified still open tonight]** **Retracted 16 mm value still feeding a calculation.** `bidet_paper_1.tex` lines 685 and 689 still use `L ≈ r_clean ≈ 0.016 m` in the Péclet estimate (`Pe ≈ 3.2×10⁸`), even though the file's own revision header (line 14) records that 16 mm was a numerical slip corrected to 26 mm. The stated conclusion `Pe ~ 10⁸–10⁹` is unaffected by whether L is 16, 26, or 6.5 mm — but a reader who checks will find a retracted number still in use.
4. **[MEDIUM — verified still open tonight]** **Internal number mismatch.** Conclusions item 4 (line 1042) says debris is removed in `~8 ms`; the body (line 695) gives `≈ 3.4 ms`. Both are "~ms", but they should trace to one parameter set.
5. **[LOW — verified still open tonight]** **Two orphan references.** `BohrEtAl1993` and `BushAristoff2003` appear in `bidet_optimization_1.tex`'s bibliography (lines 1120, 1125) but are never cited in its body — that paper has no hydraulic-jump section. They are cited in the policy paper. Cite or drop.
6. **[LOW]** Version clutter: 8 `.tex` files for 2 papers, plus four empty (0-byte) `*Notes.bib` files that nothing consumes (both papers use inline `thebibliography`; the policy paper even sets `% !BIB program = none`).
7. **[MEDIUM — for the policy paper specifically]** It is an advocacy/review document making health and economic claims (hemorrhoid prevalence, \$770M healthcare costs, 15–31M trees). That is a different publication object from the optimization paper and a different decision. Posting it to arXiv alongside a fluid-mechanics paper is a judgment call about how the pair reads.

## Carried-forward flags from prior audits (status)

All `AUDIT.md` (2026-06-30) items re-checked tonight:

- **S1 — stray `$` in `bidet_optimization_1.tex` (BLOCKER, was fatal / no PDF)** — **resolved** in the audit.
- **S2 — stale retracted `r_clean ≈ 16 mm` in the policy conclusions** — **resolved** in the audit.
- **S3 — 164 pt / 51.7 pt table overflows** — **still open**, 6 overfull boxes in the current log.
- **S4 — other cosmetic overfull boxes** — **still open** (not layout-breaking).
- **S5 — residual 16 mm Péclet inputs** — **still open**, verified at lines 685/689.
- **S6 — 8 ms vs 3.4 ms mismatch** — **still open**, verified at lines 1042 and 695.
- **S7 — two uncited hydraulic-jump references** — **still open**, verified at lines 1120/1125.
- **Author/affiliation placeholders** — **still open**, verified at lines 69–71.
- **Empty `*Notes.bib` files** — **still open** (inert).

## Concrete next actions (shortest path to preprint)

1. Fill the author block in `bidet_optimization_1.tex` (and settle the department).
2. Fix `tab:stresscomparison` — promote to `table*` full-width is the least disruptive of the three options. Recompile and confirm the overfull count drops.
3. Reconcile the three numbers so everything traces to one parameter set: the Péclet L (lines 685/689), the 8 ms vs 3.4 ms removal time, and confirm r_clean = 26 mm is used consistently throughout.
4. Cite or delete `BohrEtAl1993` and `BushAristoff2003`; delete the four empty `.bib` files.
5. Move the six superseded/test `.tex` files into an `archive/` subfolder and add a one-line `% CANONICAL` comment at the top of each of the two live files. Eight `.tex` files for two papers is how the wrong one gets edited.
6. Post the optimization paper (physics.flu-dyn). Decide separately about the policy paper.
7. *For experimental validation, later:* the one measurement is already correctly identified by the papers themselves — **the critical detachment shear stress τ_c for adherent fecal material**, which sets everything since Q* ∝ τ_c^{1/3}. **Observable:** the wall shear stress at which a standardised soil surrogate detaches from a skin-mimicking substrate. **Apparatus:** an impinging water jet at controlled Q and standoff onto a coated substrate, with detachment scored optically; wall shear from the Watson solution the paper already uses. **Regime:** τ from ~0.5 to ~50 Pa, bracketing the assumed 5 Pa. **Expected signal:** a detachment threshold; the paper's entire optimum shifts as τ_c^{1/3}, so even a factor-of-3 error in τ_c moves Q* by only ~1.4× — worth stating, since it means the design conclusion is robust to a fairly wide τ_c uncertainty. **The control that would falsify the framework:** if detachment depends strongly on impingement *time* or on normal pressure rather than on wall shear, then a shear-stress criterion is the wrong model and the whole optimization is built on the wrong objective.
