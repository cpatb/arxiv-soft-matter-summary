# Overnight Publication-Readiness Sweep — Claude Code Delegation Brief

**Paste everything below the line into Claude Code, running with `C:\Users\compl\Desktop\Projects` as the working directory.**
Recommended launch: `claude` in that folder, then paste. Or run headless:
`claude -p "$(cat CLAUDE_CODE_READINESS_SWEEP.md)" --dangerously-skip-permissions` (only if you trust an unattended run to write files into this tree).

This is a **read-and-assess** job, not a fix job. You are grading, not editing. The only files you create are the ones this brief tells you to create (`READINESS.md` per project + one master report). **Do not modify any manuscript, code, or data file.**

---

## ROLE

You are a soft-matter / fluids research assistant auditing a working-area folder of ~12 projects to decide **which are ready to move to a next publication stage**. The PI is Charlie Burton (soft matter & fluids; Driscoll Lab / Northwestern context). You have deep-dived this folder before — respect the existing `AUDIT*.md`, `README.md`, `REVIEW.md`, and `OVERNIGHT_*.md` files as prior ground truth and build on them rather than redoing them.

## SCOPE

Folder: `C:\Users\compl\Desktop\Projects`

Assess every project directory and every standalone paper. Known inventory (verify against the actual tree — folders may have changed since this brief was written):

**Code / computational projects**
- `slime mold/` — Physarum CTA-rail growth solver + `bike/` sub-project on Chicago protected-bike-lane gaps. Mature, 5× reviewed, headline reproduces.
- `hydrogels/` — PEGDA/PEGMEA drying-fracture: napari wrinkle-morphology pipeline + drying roadmap LaTeX + (live physics solver lives in the Obsidian vault, not here).
- `tier0_wyart_cates_closure/` — reduced-order Wyart–Cates suspension closure; tests + autodiff calibration; "portfolio piece #1."
- `bio_fluids_bridge/` — same closure ported to blood shear-thinning; validated.
- `thrombosis_jamming/` — SIPA-as-jamming toy with a falsifiable critical-occlusion prediction.
- `glymphatic_transport/` — Péclet regime map (dimensional analysis, not a sim).
- `bloodstain_bpa/` — ballistic region-of-origin reconstruction; quantifies straight-line-method bias.
- `uroflow_splashback/` — Rayleigh-breakup + stream-tracking splash-minimization; bidet sibling.

**LaTeX manuscripts / papers**
- `paper1/` — drop impact TPM vs cornstarch: a research *plan*, not a data paper (target PRF).
- `bidet/` — two matched revtex papers (optimization + policy), 4 versions each; canonical = `bidet_optimization_1.tex`, `bidet_paper_1.tex`.
- `friday_experiments.tex` (top level) — five oobleck free-surface experiments (proposal).
- `chicago_run_bike_monitor.TEX` (top level) — DIY air-quality sensing proposal.

**Reference / archived (assess but expect LOW readiness)**
- `LDV Paper/` — reference PDFs + slide decks, **no manuscript source**.
- `ARCHIVE/SCM/` — old C. elegans work, correctly archived. Note and skip deep analysis.

Ignore build artifacts (`.aux .log .fls .fdb_latexmk .synctex.gz .out .toc .bbl .blg`), `__pycache__/`, `.ruff_cache/`, `.pytest_cache/`, `.git/`.

## WHAT "READY" MEANS — SCORE AGAINST THREE BARS

For **each** project, score readiness 0–5 (rubric below) against **all three** target stages, because different projects are ripe for different next steps:

1. **PREPRINT / arXiv** — could go up now (or with ≤1 day polish) as a preprint. Needs: a coherent manuscript or a write-up that is ≥80% a manuscript, defensible core claims, figures present, no fatal internal contradictions, author/affiliation fillable. Novelty need not be fully defended; correctness must be.
2. **EXPERIMENTAL VALIDATION** — a theory/model/sim project whose next step is a concrete, hand-to-a-lab experiment. Needs: a *specific, falsifiable prediction*, an identified measurable observable, a realistic apparatus/parameter regime, and named controls. Score high when "here is the one measurement that would confirm or kill this" is already essentially written.
3. **JOURNAL SUBMISSION** — referee-ready manuscript for a named journal. Needs everything for preprint PLUS: novelty positioned against the literature, methods complete and reproducible, all self-flagged physics/consistency issues resolved (not just flagged), figures publication-grade, references complete (no placeholders), a plausible target journal.

### 0–5 scale (apply per bar)
- **5 — Ready now.** Submit/post/hand-off with only cosmetic touch-up. No open blockers.
- **4 — Ready after a named, bounded task** (≤1–2 days): e.g. fill author block, resolve one flagged inconsistency, regenerate one figure.
- **3 — Close; a real but well-defined gap remains** (a section, a calibration run, a novelty positioning pass).
- **2 — Substantive work remaining** but the spine exists (core result or working code present).
- **1 — Early.** Idea + partial artifact; not near this bar.
- **0 — N/A or not applicable to this bar** (e.g. a reference library scored for journal submission; a pure-sim toy scored as-is for experimental validation when no prediction is articulated).

Give each project a **primary recommended next stage** = the bar with the highest score, and a one-line justification.

## HOUSE RULES (match the PI's existing convention — this is important)

- **Conservative-fix / flag-don't-edit policy.** This is an *assessment*. Do not alter scientific claims, code, or manuscripts. Where you spot a physics error, inconsistency, placeholder, or overstated claim, **record it as a blocker/flag in the READINESS output** — never fix it in place. (The existing `AUDIT_papers.md`, `bidet/AUDIT.md` set this tone: fatal LaTeX-compile fixes were the only edits, and every physics concern was flagged for author judgment.)
- **Respect prior audits.** Where an `AUDIT*.md` already lists issues (e.g. paper1's ε_jam inconsistency and "orders of magnitude" overstatement; friday_experiments' R∝t^{1/8} vs t^{1/2} exponent; chicago monitor's σ_y-vs-σ_z plume model; bidet's 164 pt table overflow, residual 16 mm Péclet inputs, and author placeholders), **carry those forward as known blockers** and check whether they've since been resolved. Do not re-derive from scratch what's already documented; do verify status.
- **Reproducibility is evidence, not a requirement to run.** You may attempt a quick `python run_*.py` / compile check **only if it is fast and safe** and helps grade "does the headline still reproduce." If a run is slow, needs GUI (napari), or needs data you don't have, **do not run it** — grade from code + prior REVIEW/AUDIT logs and say so. Never install heavy dependencies or start long jobs.
- **No fabricated data or results.** If you can't verify a number, say "unverified" and cite where it would come from.
- **Physics sanity, lightly.** Spot-check dimensional consistency and order-of-magnitude on each project's headline claim; flag anything that fails. Don't launch a full re-derivation.
- **LaTeX manuscripts:** identify the **canonical file** among versions (e.g. bidet's `_1.tex` files; ignore `_test_*`, `_pre-review`, byte-identical duplicates). Note compile status if a `.pdf` exists and is recent; do not run pdflatex unless trivial.

## DELIVERABLE 1 — per-project `READINESS.md`

Write one `READINESS.md` **inside each assessed project folder** (for the two top-level standalone papers, write `friday_experiments_READINESS.md` and `chicago_run_bike_monitor_READINESS.md` at the top level). Use exactly this template:

```markdown
# Readiness — <project name>

**Assessed:** <date>  ·  **Canonical artifact:** <file(s) or "code: run_*.py">
**One-line status:** <e.g. "Complete, self-validated model; needs one real measurement to submit.">

## Scores (0–5)
| Stage | Score | Primary blocker to next point |
|---|---|---|
| Preprint / arXiv | X | ... |
| Experimental validation | X | ... |
| Journal submission | X | ... |

**Recommended next stage:** <one of the three>  ·  **Why:** <one line>

## What exists (evidence)
- <artifacts, tests, figures, prior review/audit status, reproduction result if run>

## Blockers to the recommended next stage (ranked)
1. **[severity]** <specific, actionable — e.g. "Fill author/affiliation placeholders in bidet_optimization_1.tex">
2. ...

## Carried-forward flags from prior audits (status)
- <flag> — <resolved / still open / not checked>

## Concrete next actions (the shortest path to the recommended stage)
1. <the single most valuable next step>
2. ...
- If experimental validation is the recommended stage, state **the one measurement**: observable, apparatus, parameter regime, expected signal, and the control that would falsify it.
```

Keep each READINESS.md tight (roughly one page). Rank blockers by what stands between the project and its recommended stage.

## DELIVERABLE 2 — master report `READINESS_SWEEP_<date>.md` (top level)

One file at `C:\Users\compl\Desktop\Projects\READINESS_SWEEP_<date>.md` containing:

1. **Executive ranking table** — all projects, sorted by "closest to a publishable next step first":

   | Rank | Project | Best stage | Score | One-line: what it needs |

2. **Three shortlists:**
   - **Post as preprint now / this week** (preprint score ≥4).
   - **Ready to design the validating experiment** (experimental-validation score ≥3, with the one measurement named).
   - **Push toward journal submission** (journal score ≥3), with the named target journal and the gap.

3. **Per-project one-paragraph rationale** (2–4 sentences each) — the reasoning behind the scores, so the PI can argue with it.

4. **Cross-cutting observations** — e.g. the "one closure, many fluids" cluster (`tier0` → `bio_fluids_bridge` → `thrombosis_jamming`) as a possible *single* methods paper vs. three notes; version-hygiene issues; anything blocking multiple projects.

5. **A "do this first" recommendation** — if the PI has one evening, which project + which action yields the most publication progress.

## EXECUTION ORDER

1. `ls`/tree the folder; reconcile against the inventory above; note additions/removals.
2. For each project: read its README/AUDIT/REVIEW + canonical artifact(s) + skim code; (optional, fast/safe only) reproduce headline; score against the three bars; write its `READINESS.md`.
3. Write the master `READINESS_SWEEP_<date>.md`.
4. End with a 10-line summary to stdout: the ranking table + the "do this first" line.

Work through **every** project before writing the master report. Be decisive with scores and blunt about blockers — the value is in a defensible ranking, not diplomacy. If something is genuinely not ready, say so and say exactly why.
