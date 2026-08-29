# Review Dossier — `chicago_slime_mold.tex` (Physarum CTA-rail expansion)

**Date:** 2026-07-23 · **Pipeline:** validation → postdoc rigor → full-professor significance → reviser (+ light humanizer) → proof
**Input:** `doc/chicago_slime_mold.EDITED.tex` · **Output:** `doc/chicago_slime_mold.REVIEWED.tex` (originals untouched)
**Proof:** braces balanced (266/266), environments balanced (19/19), all `\cite`/`\ref` resolve, one intended `[INSERT REPOSITORY URL/DOI]` placeholder remains.

## Verdict

**POST-READY for arXiv (physics.soc-ph, cross-list cs.SI) after you fill two fields** — the repository URL and your department. All correctness, consistency, framing, and prior-art fixes are applied. The paper is correct, coherent, and honestly scoped; the three reviewers agree it clears the preprint bar. One genuine methodological strengthening (a cost-matched baseline) is a *journal* item, and the manuscript now says so explicitly instead of overclaiming.

## Author-decision queue (what only you can resolve)

1. **Repository URL/DOI** — the manuscript now has a "Code and data availability" paragraph with a `[INSERT REPOSITORY URL/DOI]` placeholder. Fill it (a GitHub/Zenodo link) before posting; reproducibility is this paper's main asset.
2. **Department** — the affiliation reads "Northwestern University" with the department deliberately left off (it's written three ways across your papers). You confirmed **Physics and Astronomy** for the bidet paper; apply the same here unless this project belongs to a different department.
3. **Confirm the track-accounting sentence** — I added a sentence stating Physarum `TL` counts *new track only* (existing corridors reused at reduced cost), which is what reconciles the γ-sweep minimal tree (117 km) sitting below the MST (147 km). Confirm that matches what the code actually does. It's flagged inline.
4. **Post now, or add the cost-matched baseline first?** See below — this is the one real call.

## The one real call: the cost-matched baseline

All three reviewers converged on the same soft spot. The paper benchmarks a redundant 247 km Physarum network against minimum-cost *trees* (~145 km), and "raises fault tolerance from 0 to 0.89" is partly **definitional** — any tree has FT=0. The fair control is a **cost-matched redundant baseline** (MST augmented with the cheapest cycles up to equal TL), which is currently absent.

- **The postdoc** calls this a must-fix before arXiv.
- **The professor** calls it a journal item and says GO for arXiv now, because the framing is honest.
- **What I did:** softened the untested "dominates the trees at a given cost" claim to "occupies a region the baselines never reach," and added an explicit note that the cost-matched baseline is the key future comparison. That makes posting *now* honest rather than overclaiming.

**So:** you can post as-is (correct, honestly scoped), **or** spend an afternoon running one cost-matched baseline (it runs on your existing candidate graph — no new data) and an edge-count-independent robustness metric, which would materially harden the central claim and move it toward journal quality. My recommendation: **post now** to establish the preprint, then add the baseline for the journal version — unless you'd rather hold and strengthen first.

## What was applied to the .REVIEWED copy (12 edits, all commented `% REVIEW`)

- **Abstract headline reconciled** to Table 1: `1.74×`→`≈1.7×`, quoting the 8-seed mean (247.5±11.1 km) and the single-run 256 km; removed the cherry-picked "30% vs MST"; FT stated as 0→0.88±0.02.
- **Equity foregrounded** — one abstract clause + contribution (i) reworded to name transit-desert/equity weighting as the distinguishing element (framed as a combination, not a priority claim).
- **Prior-art added** — Kay et al. *Sci. Rep.* 12:1322 (2022) [the Toronto TTC slime-mould paper] and Zhang et al. *Sci. Rep.* 5:10794 (2015), with one sentence distinguishing this work. Real, verified references.
- **Reproducibility** — stated the missing decay rate **μ=1** (confirmed from `physarum.py`), plus batch 12 and seed 3.
- **Defined TC** in the metrics list (was used in Table 1 but undefined).
- **Track accounting** sentence added (resolves 117<147 km).
- **Pareto overclaim softened** + cost-matched-baseline future-work note.
- **Conclusion overclaim fixed** — hub-and-spoke matches Physarum on trip distance (~3%), so the claim is now robustness-focused, not "markedly more efficient."
- **Code/data availability** paragraph added (with the repo placeholder).
- **Humanizer:** no changes needed — the prose is already your understated voice.

## Deferred to the journal version (flagged in-file, not applied — need a code re-run or your decision)

(a) cost-matched redundant baseline; (b) an E-independent robustness metric alongside FT (FT currently scales with edge count); (c) demand-weight sensitivity over (α,β,κ); (d) quantify the γ=1.3 knee and reconcile the "γ>1 ⇒ tree" linearization with the observed loopy regime; (e) ridership (LODES/ACS) calibration of demand.

## Venue guidance (from the professor screen)

arXiv **physics.soc-ph** + cross-list **cs.SI** now. Best *published* home is **Transport Findings** (short format — the honest scope is an asset there) after the deferred items; a full journal (*PLOS ONE*, *Sci. Rep.*, *Environment and Planning B*) would want the cost-matched baseline and weight sensitivity.

## Sources
- [Kay et al., Stepwise slime mould growth as a template for urban design, Sci. Rep. 12:1322 (2022)](https://www.nature.com/articles/s41598-022-05439-w)
- [Zhang et al., A Biologically Inspired Network Design Model, Sci. Rep. 5:10794 (2015)](https://www.nature.com/articles/srep10794)
