# Multi-Tier Preprint Review Pipeline — Claude Code Delegation Brief

A staged, gated review harness for any manuscript you intend to post as a preprint. It takes one manuscript, runs it through four independent tiers, and produces a single dossier telling you exactly what to fix before you post — and, at the end, a humanized final draft.

**Current preprint candidates (from `READINESS_SWEEP_2026-07-22.md`):**
- `bidet/bidet_optimization_1.tex` (Preprint 4 — the arXiv-shaped object)
- `slime mold/doc/chicago_slime_mold.tex` (Preprint 4)
- `slime mold/bike/doc/chicago_bike_gap_brief.tex` (Preprint 3 / *Transport Findings*) — optional

**Paste below the line into Claude Code with `C:\Users\compl\Desktop\Projects` as the working directory.** Run it once per manuscript. Set `TARGET` to the manuscript on the first line of your instruction, e.g.:
`TARGET = bidet/bidet_optimization_1.tex. Read CLAUDE_CODE_PREPRINT_REVIEW_PIPELINE.md and follow it exactly.`

---

## SAFETY MODEL (same as the edit pass)

Never modify the original. All tier outputs and the final humanized draft go to **new files** beside the source: `<stem>.REVIEWED.tex` (the working copy the pipeline edits) and `<stem>_REVIEW_DOSSIER_<date>.md` (the report). The original `.tex` is immutable ground truth. Not under git → the untouched original is the backup. No data/code files touched. No dependencies installed. No long jobs.

## PIPELINE SHAPE

Four tiers, run **in order**, each a subagent with a fresh context and a single mandate. Tiers 1–3 are **read-and-report** (they do not rewrite the manuscript; they produce findings). After Tier 3 you get a consolidated fix list; the pipeline applies the accepted fixes to the `.REVIEWED` copy. Tier 4 (humanizer) runs **last, on the near-final text**, because reviewer fixes change content and you humanize once the science is settled.

```
Tier 0  Intake & compile-baseline  (must pass to proceed)
Tier 1  Validation / technical correctness   → gate: no OPEN blockers
Tier 2  Post-doc content screen (rigor)       → findings
Tier 3  Full-professor content screen (significance/framing/venue) → findings
        ── consolidate → apply accepted fixes to .REVIEWED ──
Tier 4  Humanizer voice pass (near-final)      → .REVIEWED final + diff notes
Tier 5  Final proof (recompile, refs, consistency)  → PASS/FAIL
```

Use the Task tool to spawn each tier as its own subagent so contexts stay clean and the personas don't bleed. Give each subagent ONLY the manuscript text + the prior tier's dossier section, not this whole brief.

---

## TIER 0 — INTAKE & BASELINE

- Identify the canonical source (ignore `_test_`, `_pre-review`, `.EDITED` unless the `.EDITED` is the intended input — if a `.EDITED` copy exists for this TARGET, **start from it**, since the edit pass already applied fixes; note that you did).
- Compile it (pdflatex, in a scratch dir so no original build artifact is destroyed — learn from the edit pass's one mistake). Record page count, undefined refs/cites. If it does not compile, STOP and report; a broken source can't be reviewed.
- Extract: title, abstract, claimed contributions, target venue if stated, figure/table inventory, reference count + any `[PLACEHOLDER]` strings.

## TIER 1 — VALIDATION (technical correctness)  ·  *gate*

Subagent mandate: *"You are a meticulous technical referee. Find everything that is wrong, unsupported, or unreproducible. You are not here to judge importance — only correctness."* Check:
- **Dimensional consistency** of every displayed equation; flag any that don't balance.
- **Numerical claims** recomputed from stated inputs where possible (the sweep already did this for several — carry forward and spot-check, don't blindly trust).
- **Internal consistency**: every number that appears twice agrees; every threshold/definition used consistently (this is where `paper1`'s ε_jam class of bug lives).
- **Citations**: no `[PLACEHOLDER]`; no uncited bibitems; claims attributed to real, correctly-described sources; no citation of a result the source doesn't contain.
- **Figures**: referenced figures exist; axes/units sane; no figure promised in text but absent.
- **Reproducibility**: if the paper leans on code in this repo, confirm the headline still reproduces (cheap runs only).
- **Overclaiming**: quantitative words ("orders of magnitude", "dramatically") checked against the paper's own numbers.

Output: a severity-ranked table (BLOCKER / MAJOR / MINOR), each with location, the problem, and the specific fix. **Gate: any BLOCKER stops the pipeline** — report and halt; those are author decisions (e.g. a real citation only you can supply). MAJOR/MINOR flow downstream.

## TIER 2 — POST-DOC CONTENT SCREEN (rigor & methods)

Subagent persona: *"You are a sharp senior postdoc in soft matter / fluids, one desk over from the author. You will be the one who has to defend this at group meeting. Read for whether it would survive a careful, hostile referee."* Focus:
- Are the **methods complete enough to reproduce**? What's missing (a parameter, a protocol step, an error bar, a control)?
- Are the **assumptions stated and defended**? Which load-bearing assumption is most likely to draw fire, and is it addressed?
- Is the **error/uncertainty treatment** honest and sufficient? Any result quoted without a stated uncertainty or sensitivity?
- Are there **obvious missing controls or alternative explanations** the paper doesn't rule out?
- **Robustness**: does any conclusion hinge on an un-measured parameter? (If so, is the sensitivity shown — e.g. bidet's Q* ∝ τ_c^{1/3} robustness argument.)

Output: findings ranked by "how likely to sink it in review," each with a concrete strengthening action. Not prose edits — substance.

## TIER 3 — FULL-PROFESSOR CONTENT SCREEN (significance, framing, venue)

Subagent persona: *"You are a senior PI and frequent journal referee deciding whether this is worth the field's attention and where it belongs. Be direct about the 'so what.'"* Focus:
- **The 'so what' test.** In one sentence: what does the field know after this that it didn't before? If that sentence is weak, say so plainly.
- **Novelty positioning.** What's the nearest prior work, and is the delta argued? (For bidet: is the framework's contribution positioned? For slime mold: is the Lake-Michigan-constraint + cost-field + transit-desert-demand *combination* argued as the novelty, vs. yet-another-Physarum-on-transit?)
- **Framing & narrative.** Does the abstract/intro promise what the body delivers? Is the strongest result foregrounded?
- **Venue fit.** Name 1–3 realistic venues (incl. preprint category, e.g. physics.flu-dyn / physics.soc-ph) and what each would demand. Is this better as a short-format (*Transport Findings*, *Letters*) or full paper?
- **Referee objections.** List the top 3 objections a real referee will raise, and whether the paper pre-empts them.
- **Post-vs-wait call.** Is this ready to post as a preprint now, or does one gap make posting premature? (Preprint tolerance is high; the bar is "correct and coherent," not "certain to be accepted.")

Output: a go/refine recommendation with the 3 highest-leverage improvements, and a named venue.

## CONSOLIDATE → APPLY

Merge Tiers 1–3 into a single ranked fix list, tagged: **[mechanical]** (apply now to `.REVIEWED`), **[author decision]** (leave flagged — placeholders, novelty claims, ambiguous physics), **[content — apply if unambiguous]**. Apply the mechanical + unambiguous-content fixes to `<stem>.REVIEWED.tex` using the edit pass's inline-comment convention (`% REVIEW [tier N]: …`). Everything else stays flagged in the dossier for you.

## TIER 4 — HUMANIZER PASS (voice & authenticity, near-final)

Runs on the near-final `.REVIEWED` text. **This is a voice/quality pass, not a deception pass** — the goal is prose that reads as a human physicist wrote it, which also removes the tells that trip false-positive AI-detectors. Mandate:
- Kill LLM tells: uniform paragraph lengths, "it's worth noting," "moreover/furthermore" scaffolding, tricolon overuse, hedge-everything qualifiers, the "not only X but also Y" reflex, section-ending summary sentences that restate the section.
- Restore human cadence: variable sentence length, occasional short declaratives, field-idiomatic phrasing, a real authorial voice consistent with the author's other papers **in this repo** (sample `friday_experiments.tex` and the slime mold write-up for the author's actual register — match it).
- **Preserve every number, symbol, equation, citation, and claim exactly.** Voice only. Any change that touches meaning is out of scope and goes back to the content tiers.
- Do not inflate. The author's register is precise and understated; keep it.

Output: the humanized `.REVIEWED.tex` + a short note listing the *kinds* of changes made (not a full diff), and a flag if any passage resisted humanizing because it was doing real technical work.

## TIER 5 — FINAL PROOF

Recompile the humanized `.REVIEWED.tex` in a scratch dir. Confirm: compiles, page count sane, zero undefined refs/cites, no `[PLACEHOLDER]` survived, every author-decision flag from the dossier is either resolved or explicitly listed as outstanding. Produce the final line: **"POST-READY" or "N author decisions outstanding."**

---

## THE DOSSIER (single output file per manuscript)

`<stem>_REVIEW_DOSSIER_<date>.md`:
1. **Verdict line** — POST-READY / N decisions outstanding, and the recommended venue.
2. **Author-decision queue** — the short list only you can resolve, each with options and the pipeline's recommendation.
3. **Tier 1 validation table** (severity-ranked, with which were auto-applied).
4. **Tier 2 postdoc findings** (rigor).
5. **Tier 3 professor findings** (significance/framing/venue + top 3 referee objections).
6. **Humanizer notes.**
7. **Final proof result.**
8. **Diff summary** — what changed between original and `.REVIEWED`, mechanical vs content vs voice.

Keep it skimmable. The author-decision queue at the top is the part you act on; everything else is evidence.

## GUARDRAILS

- Personas are lenses, not authorities — a subagent that can't support a finding drops it. Better a short honest dossier than a padded one.
- No fabricated citations, data, or numbers, ever. Placeholders stay flagged.
- Preprint bar = correct + coherent + honestly scoped, NOT guaranteed-to-pass-review. Don't let the professor tier talk you out of posting sound work; its job is to sharpen, and to catch the rare "posting this now would embarrass you" case.
- End with a 10-line stdout summary: verdict, venue, and the author-decision queue.
