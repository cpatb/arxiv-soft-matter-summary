# Readiness — slime mold (Physarum CTA-rail growth + `bike/` gap analysis)

**Assessed:** 2026-07-22  ·  **Canonical artifact:** `doc/chicago_slime_mold.tex` (9 pp) + `code/run_experiments.py`; sub-project `bike/doc/chicago_bike_gap_brief.tex` (6 pp) + `bike/code/bike_run.py`
**One-line status:** The most finished thing in this folder — a complete, 5×-reviewed, exactly-reproducing manuscript held back from arXiv by a placeholder author line.

## Scores (0–5) — main rail project

| Stage | Score | Primary blocker to next point |
|---|---|---|
| Preprint / arXiv | 4 | `\author{Prepared for Charlie}` placeholder; Table 1 Steiner value no longer matches current networkx |
| Experimental validation | 0 | N/A — an urban-network optimization model; there is no lab experiment to hand off |
| Journal submission | 3 | Novelty not positioned against a crowded post-Tero transit-network literature; demand is normalized, not ridership-calibrated |

## Scores (0–5) — `bike/` sub-project

| Stage | Score | Primary blocker to next point |
|---|---|---|
| Preprint / arXiv | 3 | It is a policy *brief*, not a paper: placeholder byline, no abstract/intro/related-work in paper form. **Also: it does not compile from its own directory** (see blocker 5) |
| Experimental validation | 0 | N/A |
| Journal submission | 3 | Fits *Transport Findings* (short format) almost as-is; needs a real byline and paper-form framing |

> **Correction (made during the Stage-2 edit pass):** an earlier draft of this file listed
> "restate modeled-vs-measured in the body, not just in `REVIEW_bike.md`" as a blocker.
> That was wrong — the brief already states it plainly in its own body at lines 179–180
> ("The mileage and crash figures are measured. The *ideal network* and the gap ranking
> are modeled…"). No such blocker exists.

**Recommended next stage:** Preprint / arXiv (physics.soc-ph or cs.SI)  ·  **Why:** the manuscript is complete, compiles, and its headline reproduces to 2 d.p.; the only thing between it and a posting is an author block.

## What exists (evidence)

- **Manuscript:** `doc/chicago_slime_mold.tex` + current PDF (9 pp). Compiles with zero undefined refs/citations (per `REVIEW.md` Pass 5).
- **Reproduction:** headline verified **exactly** on 2026-06-30 (`AUDIT_2026-06-30.md`): Physarum γ=1.3 → TL 247.49±11.10 km, MD 18.69±0.34 km, FT 0.884±0.019, matching README and paper. Not re-run tonight (multi-seed sweep is not a fast/safe check); graded from that audit.
- **Verification depth:** `REVIEW.md` records 5 independent passes — model fidelity to Tero (2010), Kirchhoff residual = 0.000, geography/lake-constraint checks, figures-vs-numbers, reproducibility. Three real issues were found and fixed in those passes.
- **Figures:** fig1–fig6 present and regenerated 2026-06-30.
- **`bike/`:** real City-of-Chicago Socrata data with every measured number traced to its dataset ID in `REVIEW_bike.md` (68.74 mi protected = 15.4% of 445.9 mi; 17,677 cyclist crashes). figB1–figB4 regenerate. 6-page brief PDF present.

## Blockers to the recommended next stage (ranked)

1. **[HIGH]** `doc/chicago_slime_mold.tex` line 16 reads `\author{Prepared for Charlie}`. Replace with a real author/affiliation block before any posting.
2. **[MEDIUM]** **Known numeric drift, unresolved.** Table 1 / README report Steiner = 145.2 km / MD 29.29 km, but under networkx 3.6.1 the Mehlhorn Steiner tree now coincides exactly with the MST (147.0 / 26.9). Either pin the networkx version in a reproducibility note or add one line to Table 1 saying the two tree baselines coincide under current networkx. Conclusions are unaffected — both are FT = 0 trees the Physarum family dominates — but a reader who re-runs will hit the mismatch.
3. **[LOW]** `\documentclass[11pt]{article}` — fine for arXiv, but a journal will want its own class.
4. **[LOW]** `bike/doc/chicago_bike_gap_brief.tex` has no `\author`/`\title` block — its byline reads "Prepared for Charlie", the same placeholder as the rail paper.
5. **[HIGH — found during the Stage-2 edit pass, not previously flagged in any audit]** **The bike brief does not compile from its own directory.** Its four figures are `\includegraphics`'d by bare filename but live in `../figures/`, so `pdflatex` aborts with *"File `figB1_existing_protected.png' not found … Fatal error occurred, no output PDF file produced"*. Verified against the unmodified original. One-line fix: `\graphicspath{{../figures/}{./}}` (applied in the `.EDITED` copy).

## Carried-forward flags from prior audits (status)

- `bike_run.py` hardcoded Unix `/tmp/` output paths (`AUDIT_2026-06-30.md`) — **resolved**; outputs now land in `bike/figures/` and `bike/data/`, both populated.
- Watanabe citation year 2014→2011 (`REVIEW.md` Pass 1) — **resolved**.
- Steiner over-attribution to Bonifaci et al. (`REVIEW.md` Pass 1) — **resolved** (softened to "empirically observed").
- Steiner-baseline networkx drift (`AUDIT_2026-06-30.md`) — **still open**, documented but not addressed in the paper.
- Residual limitations in paper §7 (schematic line geometry, normalized demand, structural rather than operational fault tolerance) — **still open by design**, disclosed in the manuscript. These are what cap the journal score.

## Concrete next actions (shortest path to preprint)

1. Fill the author/affiliation block in `doc/chicago_slime_mold.tex`; recompile.
2. Add a one-sentence footnote to Table 1 recording the networkx Steiner/MST coincidence, plus a pinned `requirements.txt` with the version used for the reported numbers.
3. Post to arXiv (physics.soc-ph). Include the `bike/` study as a companion or fold it in as a §8 application.
4. *For journal, later:* the real gap is positioning. Physarum-on-real-transit-networks has been done repeatedly since Tero (2010); the defensible novelty here is the **combination** of a hard no-build constraint (Lake Michigan), a construction-cost field folded into effective tube length, and transit-desert-weighted demand — say that explicitly in the introduction, with citations to the closest prior applications.
