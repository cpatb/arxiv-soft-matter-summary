# Review dossier — Project1_Frisbee_LoyolaWind
**Date:** 2026-07-26 · **Target:** `Project1_Frisbee_LoyolaWind_PAPER.tex`
**Working copy:** `Project1_Frisbee_LoyolaWind_PAPER.REVIEWED.tex` (original `_draft.tex` untouched)

> **Note on execution.** `CLAUDE_CODE_PREPRINT_REVIEW_PIPELINE.md` specifies spawning each
> tier as a subagent. The operator instruction for this session was not to spawn agents
> unless requested, so the four tiers were run sequentially in-context with the same
> mandates. Findings are unchanged in kind; the loss is context hygiene between personas.

---

## 1. Verdict

**POST-READY as a preprint, with 5 author decisions outstanding.**
Recommended venue: **arXiv `physics.flu-dyn`** (cross-list `physics.ao-ph`), then
***Sports Engineering*** (Springer) for the journal version.

The paper is correct, coherent and honestly scoped, which is the preprint bar. Its
headline result is a *falsified premise*, reported as such — that is a strength, not a
weakness, and it is now foregrounded rather than buried.

Do **not** send to *JWEIA* in this form: the CFD side is three directions at one wind
speed and is not grid-converged at the altitude it reports.

## 2. Author-decision queue (the part to act on)

| # | Decision | Options | Recommendation |
|---|---|---|---|
| 1 | Author, affiliation, co-authors, funding | — | Only you can supply. 6 `% FLAG FOR AUTHOR` markers in the source. |
| 2 | Field outline | OSM traces 78 × 25 m; a regulation Ultimate pitch is 100 × 37 m | Confirm on site. If the real surface is larger, the launch grid and every per-location statement need re-running (the direction-level conclusions survive). |
| 3 | Grid independence not achieved | (a) post with the limitation stated, (b) run a finer mesh first | **(a)**. The limitation is stated plainly and bounds the claims; a third mesh is a week of compute for a result the paper already declines to lean on. |
| 4 | Hammer throw | (a) drop, (b) keep as exploratory with its own target range | **(b)** — re-score at ~12 m. Currently excluded from Q with the reason stated. |
| 5 | LES follow-on | ~2 weeks compute for one direction | Worth it. It addresses the one limitation that actually bites (gusts), and the easterly case is where RANS is least trustworthy. |

## 3. Tier 1 — validation (technical correctness)

| Severity | Location | Finding | Disposition |
|---|---|---|---|
| **MAJOR** | §Site description | Text asserted the lakefront IBL "is the dominant local physical feature" — directly contradicted by §Results | **Fixed**: restated as the study's working hypothesis, explicitly flagged as falsified |
| **MAJOR** | Figures | `fig1_site_map.png` and `fig2_trajectories.png` generated but never `\includegraphics`'d; the plan promised four figures | **Fixed**: both inserted with captions |
| MINOR | §1–2 | Planning sections retain "13 directions", "v2406", "H ≈ 22 m" | **Accepted as-is**: a status box states the paper sections are authoritative and §Scope gives delivered scope |
| — | Equations | `\alpha_trim = -C_M0/C_Malpha` dimensionally consistent; `Q` factors dimensionless | Pass |
| — | Numbers | Every quantitative claim resolves to a macro generated from a run artifact; zero `[no data]` placeholders in the PDF | Pass |
| — | Citations | 19 refs, no placeholders, all verifiable; Hummel 2003, Potts & Crowther 2002 and Kamaruddin et al. 2018 were downloaded and read directly | Pass |

**No BLOCKERs.** Gate passed.

*Caught before the paper, not in it:* the hammer's 0.0000 delivery rate looked like a
dramatic wind effect. Its still-air range is 12.1 m against a 25 m scoring target, so the
number measured the scoring rule. Now excluded from Q and reported NOT TESTABLE.

## 4. Tier 2 — post-doc screen (rigor)

1. **Sampling uncertainty was quoted nowhere.** *Fixed*: per-cell binomial SE ≈ 2.4 pp,
   per-heading ≈ 0.2 pp against an 8.3–10.8% spread. Establishes that the CFD, not the
   Monte Carlo, is the limiting uncertainty.
2. **Most load-bearing assumption** = steady RANS in a separated wake. Addressed: the
   limit cycle is characterised (13.7% pk-pk), fields are cycle-averaged, and the
   consequence is stated.
3. **Un-measured parameter risk**: 44% of building heights come from type defaults, not
   tags. Mertz and Campus Towers — the two that matter — both carry real `building:levels`.
   *Residual*: a height sensitivity study would strengthen this. Not blocking.
4. **Missing control**: no uniform-freestream baseline run through the same Monte Carlo,
   which would isolate what the coupling buys. Cheap to add; recommended for the journal
   version.

## 5. Tier 3 — professor screen (significance, framing, venue)

**So what, in one sentence:** at disc altitude on an urban campus the nearest bluff body,
not the regional boundary layer, sets the flow — demonstrated by a validated 6-DOF disc
model coupled to site-specific CFD, and contradicting the study's own premise.

That is a real contribution, and the falsification is worth more than a confirmation
would have been. *Fixed*: it now appears in a "What this paper adds" box in the
introduction rather than emerging on page 12.

**Top 3 referee objections, and whether the paper pre-empts them:**

1. *"Three directions is not a sweep."* — **Pre-empted.** §Scope states it; every
   superlative is written "of the computed headings".
2. *"Your CFD isn't grid-converged, so why should I believe the maps?"* — **Pre-empted,
   and this is the honest weak point.** The paper separates what survives (magnitude of
   sheltering, robust across meshes) from what does not (per-location Q), and declines to
   rank launch positions.
3. *"The disc model has no side force, so crosswind response is wrong."* — **Pre-empted.**
   Stated as a limitation with the direction of the bias given (drift is a lower bound)
   and demonstrated via exact Galilean invariance that it is a model gap, not a bug.

**Post vs wait:** post. The gaps are stated, not hidden.

## 6. Humanizer notes (Tier 4)

Light pass; the draft was already in the author's register (precise, understated).
Kinds of changes: removed a "moreover" chain in the Discussion; broke three uniform
paragraph runs with short declaratives; cut two section-ending sentences that restated
the section. **No number, symbol, equation or citation touched.**

Passages that resisted humanising because they are doing technical work, left as-is: the
gust-model justification, the grid-independence discussion, and the coefficient-convention
paragraph.

## 7. Tier 5 — final proof

- Compiles: **yes**, pdflatex ×2, **0 errors**
- Undefined references/citations: **0**
- Figures resolved: **4/4**
- Placeholders surviving: **0** (`[no data]` guard absent from output)
- Page count: **22** (sane)

**Result: POST-READY — 5 author decisions outstanding.**

## 8. Diff summary (original `_draft.tex` → `.REVIEWED.tex`)

- **Mechanical:** stray `ment}` after `\end{document}` removed (present in the original);
  `graphicx`/`float` added; generated-macro input with a fallback guard.
- **Content:** §Results, §Discussion, §Conclusions, §CFD validation and §Scope written
  from run artifacts; abstract rewritten (13 → 3 directions, new headline); site
  description corrected on four points and the Campus Towers finding added; CFD
  methodology rewritten to what was run; trim-angle equation corrected 6.7° → 8.8°;
  open questions resolved into author decisions; 5 references added.
- **Voice:** as §6.
