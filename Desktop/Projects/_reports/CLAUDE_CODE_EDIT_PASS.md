# Overnight Edit Pass — Claude Code Delegation Brief (Stage 2)

**Run this AFTER the readiness sweep finishes.** It depends on the `READINESS.md` files and `READINESS_SWEEP_<date>.md` that the sweep produces.

**Paste everything below the line into Claude Code, with `C:\Users\compl\Desktop\Projects` as the working directory.**
Hands-off launch: `claude -p "$(cat CLAUDE_CODE_EDIT_PASS.md)" --dangerously-skip-permissions`

---

## ROLE & MODE

You are applying fixes to a set of soft-matter / fluids manuscripts and code, on behalf of PI Charlie Burton. The readiness sweep has already run and graded every project. Your job now is to **actually make the fixes** on the projects that are close to a publishable stage — but under a strict safety model:

> **NEVER modify an original file.** Every edit is written to a **copy that sits next to the original**, with `.EDITED` inserted before the extension (e.g. `bidet_optimization_1.tex` → `bidet_optimization_1.EDITED.tex`; `friday_experiments.tex` → `friday_experiments.EDITED.tex`; `physarum.py` → `physarum.EDITED.py`). The original is read-only ground truth. This folder is **not** under version control, so the untouched original IS the backup.

If an original with a `.EDITED` twin already exists from a prior pass, don't clobber it — make `.EDITED2` and note it.

## PICK TARGETS FROM THE SWEEP

1. Read `READINESS_SWEEP_<date>.md` (newest) and each project's `READINESS.md`.
2. Edit **only** projects whose recommended stage is preprint or journal AND whose relevant score is **≥3** — i.e. things close enough that fixes materially move them. Skip 0–2 projects, reference libraries (`LDV Paper/`), and archives (`ARCHIVE/`). List which projects you're editing and which you're skipping, with the score that decided it.
3. Within each target, work only the blockers the `READINESS.md` already identified plus the carried-forward flags from the prior `AUDIT*.md` files. Don't hunt for new issues to "improve" — fix the known, documented ones.

## TWO TIERS OF EDIT — BOTH ALLOWED, BOTH LOGGED, BOTH ANNOTATED

**Tier A — mechanical / safe.** Apply confidently:
- Fill author/affiliation/email placeholders (author is **Charlie Burton, Northwestern University**; use that unless the file clearly indicates otherwise). Do not invent an email — insert `[email]` only if none exists anywhere in the tree, and note it.
- Fatal LaTeX compile errors (stray `$`, bracket leaks, undefined `\rowcolor` → `\usepackage[table]{xcolor}`, etc.).
- Overfull tables / boxes that break layout (e.g. bidet `tab:stresscomparison`, 164 pt): promote to `table*` full-width or wrap in `\resizebox{\columnwidth}{!}{…}`.
- Numbers the author's own revision notes already retracted (e.g. bidet policy paper `r_clean ≈ 16 mm` → `26 mm`, and the residual 16 mm Péclet inputs the AUDIT flagged as stale).
- Cross-reference / consistency sync within a doc (e.g. paper1 Objection-5 ε list vs the primary ε_jam = 0.07 and band [0.05, 0.10]).
- Uncited-bibitem cleanup, stale filename comments, empty `.bib` housekeeping.

**Tier B — physics / scientific claims.** Apply, but flag loudly and show your work:
- `friday_experiments`: axisymmetric constant-volume viscous spreading **R ∝ t^{1/2} → t^{1/8}** (Huppert 1982), and the derived 2D height decay **H ∝ t^{−2/5} → t^{−1/5}** (consistent with the paper's own X ∝ t^{1/5}). Re-examine the downstream liquid-dam-vs-Newtonian comparison once corrected.
- `chicago_run_bike_monitor`: H2 highway plume — reframe the crosswind σ_y fit to a **line-source model governed by vertical dispersion σ_z(x)**, or treat the road as a finite line with along-road integration.
- `paper1`: soften "orders of magnitude" / "order of magnitude earlier" to match the paper's own quoted ranges (τ_KWW 300–600 ms vs 50–200 ms ≈ 2–10×), as the paper's own referee review recommends.
- Any other physics blocker a `READINESS.md` raised **with an unambiguous correction**. If the correct value is genuinely ambiguous or you're not confident, DO NOT edit it — leave it flagged in the changelog as "needs author decision" and explain the options.

For **every Tier-B edit**, insert an inline comment right at the change in the `.EDITED` file:
`% EDIT [physics, 2026-..]: R∝t^{1/2} → t^{1/8}. Reason: constant-volume lubrication current, Huppert 1982; follows from the paper's own Ṙ relation. Confidence: high.`
(Use `#` for Python.) This makes every scientific change greppable and reversible by eye.

## PER-PROJECT DELIVERABLES

For each edited project, alongside the `.EDITED` files write **`EDITS_<date>.md`** with:

```markdown
# Edits — <project> (<date>)
Targets chosen because: <readiness score / stage>
Files edited (original → copy): <list>

## Tier A (mechanical) — applied
| File | Location | Before → After | Why |

## Tier B (physics) — applied, needs your eyes
| File | Location | Before → After | Basis (citation/derivation) | Confidence |

## NOT edited — needs author decision
- <flag> — <why left alone, and the options>

## Compile / run check
- <.EDITED.tex compiles? pp count? or code runs / tests pass?>
```

**Verify before you call it done:** compile each `.EDITED.tex` (pdflatex, in place so assets/`.sty` resolve) and confirm it still builds and the page count is sane; for any `.EDITED.py`, run its tests / a fast sanity run if safe. Record pass/fail. If a Tier-B physics edit breaks the build or a test, revert that single edit in the copy and downgrade it to "needs author decision."

## MASTER SUMMARY

Write `EDIT_PASS_<date>.md` at the top level:
1. Table: project | files edited | # Tier-A | # Tier-B | # left for author | compiles?
2. The **highest-leverage Tier-B calls** across all projects, in one place, so the PI can accept/reject them quickly.
3. Anything you deliberately did NOT touch and why.
4. For each edited project, the remaining shortest path to the recommended stage after these edits.

## GUARDRAILS

- Originals are immutable. If you ever catch yourself about to write to a non-`.EDITED` source file, stop.
- No data files touched, ever. No new dependencies installed, no long jobs.
- No fabricated citations or numbers. A placeholder reference (e.g. `Sacanna2013` TPM synthesis) stays a flagged placeholder unless you can cite the real source with confidence — note it in "needs author decision."
- Physics you're not sure about is flagged, not edited. Better to leave a clean flag than a confident wrong "fix."
- End with a 10-line stdout summary: projects edited, count of Tier-A vs Tier-B, and the top 3 physics calls the PI should review first.
