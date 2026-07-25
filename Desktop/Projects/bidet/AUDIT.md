# Bidet Paper Audit

Date: 2026-06-30
Toolchain: MiKTeX-pdfTeX 4.21 (MiKTeX 25.4), `pdflatex` on PATH (working).
Scope: conservative audit, compile-check, and minimal safe fixes of the LaTeX
manuscripts in `C:\Users\compl\Desktop\Projects\bidet`.

---

## (a) Version lineage — which file is canonical

There are two independent manuscripts ("families"), each with several versions.
Nothing was deleted; lineage is documented below.

### Family 1 — "Optimization" paper
*Title:* "Towards the Optimal Bidet: A Fluid Mechanical Framework for Nozzle
Design and Flow Parameter Optimization."

| File | Role | Notes |
|---|---|---|
| `bidet_optimization.tex` (Mar 22) | Original draft | Superseded. |
| `bidet_optimization_1_pre-review.tex` (May 30 12:23) | Pre-review snapshot | Superseded; lacks the v3 "Verification Roadmap" section and later fixes. |
| **`bidet_optimization_1.tex` (May 30 12:32)** | **CANONICAL** | revtex4-2, REVISION v3. Most complete. |
| `_test_bidet_optimization_1.tex` (May 30 12:37) | Test artifact | Same v3 text but `\documentclass{article}`+natbib instead of revtex: a compile-portability test copy, not the real paper. |

### Family 2 — "Policy/Case" paper
*Title:* "On the Superiority of Hydrodynamic Perianal Cleansing … Bidet Adoption
Policy in the United States."

| File | Role | Notes |
|---|---|---|
| `bidet_paper.tex` (Mar 22) | Original draft | **Byte-identical** to `bidet_paper_pre-review.tex`. Superseded. |
| `bidet_paper_pre-review.tex` (May 30 12:23) | Pre-review snapshot | Identical to the original; superseded. |
| **`bidet_paper_1.tex` (May 30 12:36)** | **CANONICAL** | revtex4-2, REVISION (30 May 2026). Most complete. |
| `_test_bidet_paper_1.tex` (May 30 12:37) | Test artifact | Same content, `article`+natbib class. Compile-portability test copy. |

The two canonical papers are a matched pair and cite each other
(`\citep{BurtonOptPaper}` / `\citep{BurtonCasePaper}`).

`.bib` files (`*Notes.bib`) are all empty (0 bytes); both papers use an inline
`thebibliography`, so **no bibtex run is required** (the policy paper even sets
`% !BIB program = none`). `.aux/.log/.out/.fls/.fdb_latexmk/.synctex.gz/.pdf`
are all regenerable build artifacts.

---

## (b) Compile status of each canonical file

Command used (twice each, from the folder):
`pdflatex -interaction=nonstopmode -halt-on-error <file>.tex`

| File | Status | PDF |
|---|---|---|
| `bidet_optimization_1.tex` | **FAILED as received**, now **compiles cleanly** after 1 fix (see (d)). Exit 0 both passes. | produced |
| `bidet_paper_1.tex` | **Compiles cleanly** (exit 0 both passes). | produced |

No undefined `\ref`, no undefined `\cite`, no multiply-defined labels in either
file after fixes.

---

## (c) Issues found (severity-ranked)

**S1 — BLOCKER (fixed): unbalanced math delimiter in `bidet_optimization_1.tex`.**
Line 1104 contained a stray `$`: ``` ``$\sim$\SI{1.6}{L/min}$ exceeds commercial'' ```.
This threw `! Missing $ inserted` and produced **no PDF at all**. Fixed (see (d)).

**S2 — MEDIUM (fixed): stale retracted number in `bidet_paper_1.tex` conclusion.**
The revision header of that same file states the cleaning radius "16 mm was a
numerical slip; the correct value … is approximately 26 mm," and the body
Eq.~(rcleanval) derives ≈26 mm. But Conclusions item 2 still read
`r_clean ≈ 16 mm`. Corrected to 26 mm with an eqref cross-reference (see (d)).
Independently arithmetic-verified: for the stated parameters
(τ_crit=5 Pa, V0=11.2 m/s, a=0.65 mm, ν=1e-6) Eq.(rclean) gives ≈0.0258 m = 26 mm. ✓

**S3 — MEDIUM (NOT fixed, layout): wide tables overflow the two-column measure.**
- `bidet_paper_1.tex`, `tab:stresscomparison` (lines 759–770): Overfull \hbox
  **164 pt too wide**: in a revtex `reprint` (two-column) this table bleeds
  well past the column into the gutter/margin. This one genuinely hurts layout.
- `bidet_optimization_1.tex`, `tab:commercial` (lines 848–862): Overfull \hbox
  **51.7 pt**, same problem, milder.
  Recommended remedy left to the author (see (e)).

**S4 — LOW: other overfull boxes (cosmetic, do not break layout).**
Optimization: 18.9 pt at eq. line 576; 2.8 pt at line 737.
Policy: 12.2 pt (l.334), 3.4 pt (l.380), 19.1 pt (l.639), 32.3 pt (l.745),
29.2 pt (l.891). These are ordinary tight-line overflows, typical of revtex
two-column math, and are not layout-breaking.

**S5 — LOW: residual `16 mm` values in the policy paper's Péclet estimate
(NOT fixed).** Lines 685/689 still use `L ≈ r_clean ≈ 0.016 m` (the retracted
16 mm) in the order-of-magnitude Péclet calculation (`Pe ≈ 3.2×10^8`). The
stated conclusion (`Pe ~ 10^8–10^9`) is unaffected by whether L is 16, 26, or
6.5 mm, so this was flagged rather than changed (the "correct" L here is
genuinely ambiguous between the ideal 26 mm and realistic ~6.5 mm cases).

**S6 — LOW: minor internal number mismatch (NOT fixed).** Policy paper
Conclusions item 4 says debris is removed in "~8 ms"; the body Eq.(peval)
neighbourhood gives `t ≈ 3.4 ms`. Both are "~ms"; the discrepancy is cosmetic.
Author should reconcile.

**S7 — LOW: two uncited references in `bidet_optimization_1.tex`.**
`BohrEtAl1993` and `BushAristoff2003` appear in that paper's bibliography (the
v3 note says they were "added as modern hydraulic-jump references") but are
**never `\cite`d in its body**: the optimization paper has no hydraulic-jump
section. (They *are* cited in the policy paper.) Harmless with a manual
`thebibliography` (no error), but they are dead entries. Author may cite or drop.

---

## (d) Fixes applied (minimal, reversible)

1. **`bidet_optimization_1.tex`, line 1104** — removed the stray `$`:
   - before: `` ``$\sim$\SI{1.6}{L/min}$ exceeds commercial'' ``
   - after:  `` ``$\sim$\SI{1.6}{L/min} exceeds commercial'' ``
   Result: paper now compiles (was previously a fatal error / no PDF).

2. **`bidet_paper_1.tex`, Conclusions item 2 (line ~1032)** — replaced the
   author-retracted `r_clean ≈ 16 mm` with `≈ 26 mm` and added a cross-reference
   to Eq.(rcleanval), matching the file's own revision note and derived value.

No other content was altered. Both edits are one-line and trivially reversible.

---

## (e) Recommendations NOT applied (need author judgment)

1. **Fix the 164 pt table overflow (S3).** `tab:stresscomparison` should either
   be promoted to a full-width float (`table*` + `\begin{tabular*}` or
   `\resizebox{\columnwidth}{!}{…}`) or have its "Ratio (TP/Bidet)" column
   shortened. Same, more mildly, for `tab:commercial`. This changes float
   placement/pagination, so it is left to the author.
2. **Reconcile the residual 16 mm Péclet inputs (S5)** and the 8 ms vs 3.4 ms
   removal time (S6) so all numbers trace to one parameter set.
3. **Cite or remove the two orphan hydraulic-jump references (S7)** in the
   optimization paper.
4. **Author/affiliation placeholders.** `bidet_optimization_1.tex` still has
   `[Author Name]`, `[email@institution.edu]`, `[Institution]`. The policy
   paper is fully attributed (Charlie Burton, Northwestern). Fill in before
   submission.
5. **Empty `*Notes.bib` files** can be deleted or populated; currently inert.

### Physics / optimization sanity check (spot audit — no changes needed)
- **Dimensional consistency verified** for the core results:
  `Q* = π(τc r_max^5 ν / c2 ρ)^{1/3}` → m³/s ✓;
  `a* = (τc r_max^5 ν / c2 ρ V0³)^{1/6}` → m ✓;
  new hydraulic-jump scaling `r_h ~ (Q^5/ν^3 g)^{1/8}` → m ✓ (this correctly
  replaces the previously dimensionally-broken Watson gravity formula).
- **Numerical results reproduce** from the stated nominal parameters:
  optimization `V0≈11.48 m/s`, `d*≈0.80 mm`, `Q*≈0.347 L/min` all recomputed and
  agree; policy `r_clean≈26 mm` recomputed and agrees.
- **Objective/constraints clearly stated** (P1: min Q at fixed ΔP; P2: max
  r_clean at fixed Q; constraints C1 turbulent impingement, C2 jet coherence,
  C3 pressure ceiling). Clear and self-consistent.
- **Honestly-sourced caveats:** both papers now flag the laminar Region-II
  asymptotic as a *lower bound* on wall shear and openly note the tension that
  the standoff-corrected 1.6 L/min laminar prediction exceeds commercial flow
  rates (turbulent enhancement invoked to reconcile). The most load-bearing
  empirical unknown, the fecal detachment shear τ_c (Q* ∝ τ_c^{1/3}), is
  explicitly identified as the key measurement. No fabricated data detected;
  quantitative environmental/economic claims are cited to NRDC/Peery/etc.
