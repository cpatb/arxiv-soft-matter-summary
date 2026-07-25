# Review Dossier — `bidet_optimization_1.tex` (Optimal-Bidet nozzle framework)

**Date:** 2026-07-23 · **Pipeline:** validation → postdoc rigor → full-professor significance → consolidation/revise
**Input:** the solo-pass `bidet_optimization_1.REVIEWED.tex` (author block + abstract oblique fix already in)
**Output:** updated `bidet_optimization_1.REVIEWED.tex` (safe fixes applied; blockers flagged inline, originals untouched)
**Proof:** braces balanced (852/852), environments balanced (46/46).

## Verdict — this corrects the overnight sweep

**NOT post-ready.** The multi-agent review found problems the overnight sweep (which ranked bidet "do this first, post tonight") and my earlier solo pass did **not** catch. The paper's algebra is *not* clean, and its central novelty claim is close to false. This is exactly the case the reviewer pipeline exists for: three independent expert lenses caught what one pass missed.

The framework and the optimization *idea* are fine. But several headline results sit on a laminar law the paper itself disavows, one whole section (P2) has no valid optimum, a displayed equation contains a literal `?`, and the "no prior work" claim ignores an established impinging-jet cleaning literature that does essentially this. None is fatal to the project; all are fatal to posting *tonight*.

## Author work/decision queue (ranked — these block posting)

1. **[B1 · blocker, trivial] A literal `a^{?}` appears in a displayed equation** (the standoff derivation, after Eq. tauIIH). Never post an equation with a `?` in it. But don't just patch it to `a^5` — it lives inside B3, which is broken anyway.
2. **[B4 · the big one] The laminar/turbulent incoherence.** Every headline number uses the laminar Region-II law τ_w ∝ r⁻⁵, but the real wall jet is turbulent (r⁻¹·² to r⁻², per the Phares paper you already cite). At r_max these differ by ~10²–10³×. The "conservative lower bound" (1.6 L/min) sits *above* commercial flow rates that demonstrably work — i.e. the model misfits reality at the operating point, and is rescued by an unquantified "turbulence helps" hand-wave. **Key action (one afternoon, no new hardware, uses a reference already in your bib): plug the Phares r⁻¹·² correlation into the cleaning criterion and recompute d\*, Q\*, and every lever.** The postdoc's estimate is that several headline levers (oscillation 55%, standoff 4.6×, the Pareto 5/3 exponent) do **not** survive this. You need to know before posting.
3. **[B3 · blocker] The standoff derivation is non-reproducible and possibly wrong-sign.** The a→σⱼ, V₀→Vᶜ substitution isn't shown term-by-term, may make τ_w *increase* with H, and the σⱼ≈7 mm geometry pushes cleaning into Region I — contradicting the paper's own "constraint lives in Region II" premise. Your #1 design conclusion ("standoff is the primary lever, 4.6×") rests entirely on this. Re-derive or cut.
4. **[B2 · blocker] The entire P2 section has no valid optimum.** Under your Eq. (tauII), τ_w at fixed Q depends on Q alone — not on the a/V₀ split — so there is no interior optimum to find; the pressure-ceiling optimization (Eqs. aP2, rstar) is spurious, and (rstar) is dimensionally not a length (spurious `(V₀*)` and `a^{1/5}` factors; rcleanQ exponent `a^{-4}` should be `a^{-5}`). Decide: delete P2, or recast it as the Pareto relation you already have.
5. **[B5 · decision] C_d dropped from Q but C_v kept.** Your headline Q\*=0.35 L/min omits the 0.62 discharge coefficient (physical ≈0.22 L/min) yet is compared directly to commercial 0.42 L/min. Reconcile the coefficient bookkeeping or flag the comparison as un-normalized.
6. **[B6 · decision] Oscillation κ arithmetic doesn't close.** 0.6 × 73% = 44%, not the stated 55%; and "~80% for the ideal stadium" contradicts the 73% you compute for that same κ=1 case. Pick a κ and propagate it to abstract and conclusion.
7. **[false-claim, FIXED] The Region I/II "agree to within 30% / internally consistent" statement** was wrong by ~200× — I removed it and replaced it with the honest order-of-magnitude-matching statement (see below).
8. **[novelty, softened + needs citations] The "no prior work" claim.** Softened in the manuscript; you must add the cleaning-jet citations (list below) and the department is set to Physics & Astronomy per your earlier call.

## What was safely fixed in the .REVIEWED copy (with `% REVIEW` tags)

- **Removed the false "within 30% / internally consistent" matching claim**; replaced with the true O(10²) mismatch and a pointer that any *absolute* Q\* needs the Tier-A numerical solve. (Scaling exponents don't depend on it.)
- **Softened the novelty claim** in the abstract and intro to "no treatment frames bidet nozzle design as explicit water-minimization," acknowledging the impinging-jet cleaning literature; inserted a `% FLAG` + `[CITE cleaning-jet literature]` marker where the citations go.
- **τ_c "geometric mean of 1–10 Pa" → "mid-range"** (geometric mean is 3.16, not 5).
- **Table typos:** standoff H/d=40 penalty 3.4→3.5 (1.20→1.22); pressure ΔP=0.05 d\* 0.89→0.87.
- **Added the missing `Table~\ref{tab:oblique}`** callout; marked the ν=10⁻⁶ estimate as order-of-magnitude.
- **Added a prominent top-of-file `% FLAG` block** listing B1–B6 so they can't be missed.
- *(Not touched: the blockers above — deliberately, since they need your derivation/decision.)*

## Prior-art to add (the professor's P0 item; add real \bibitems)

The Watson-wall-jet-plus-wall-shear-detachment method is the standard toolkit of the Cambridge impinging-jet cleaning group — citing it is what stops the paper looking naive:

- Wilson, D.I. and co-workers — cleaning of soft-solid soil layers by stationary coherent impinging liquid jets, *Chemical Engineering Science* (2014). https://www.sciencedirect.com/science/article/abs/pii/S000925091400044X
- Bhagat, R.K. & Wilson, D.I. — flow in the thin film created by a coherent turbulent water jet impinging on a wall, *Chemical Engineering Science* 152 (2016). https://www.sciencedirect.com/science/article/abs/pii/S0009250916303086
- Impinging-jet cleaning of tank walls: jet length, wall curvature, *Food and Bioproducts Processing* (2018). https://www.sciencedirect.com/science/article/abs/pii/S0960308518303973
- CFD study of a bidet water-injection nozzle (2007) — the direct "someone has CFD'd a bidet nozzle" reference. https://www.koreascience.or.kr/article/JAKO200735836575826.do

(I couldn't auto-fetch the full author lists — ScienceDirect blocks scraping — so I did not fabricate \bibitems. Pull the author lists from the DOIs above.)

## Condensed tier findings

**Validation** — 3 BLOCKERS (the `a^{?}`, the dimensionally-broken rstar, the a⁻⁴→a⁻⁵ typo), 3 MAJORS (P2 has no optimum; the 30%-matching falsehood; κ arithmetic), plus ~6 MINOR numeric/label/ref issues. Verdict: clear the blockers + majors before uploading.

**Postdoc rigor** — the quantitative spine (r⁻⁵) is used outside its regime and disavowed only where it embarrasses; the "levers survive because ratios cancel" claim is false for r-dependent corrections; standoff derivation wrong-sign; τ_c sourced to bulk fecal rheology (Woolley) rather than an adhesion/detachment measurement. Single most valuable addition: the Phares cross-check (free).

**Professor significance** — competent applied/pedagogical exercise, *not* new fluid mechanics; the method is the established cleaning-jet framework. "So what" is thin. Honest but over-caveated to the point of dissolving its own headline. **Venue: European Journal of Physics** (pedagogical) is the right journal home, not PRF/JFM (which would reject on significance + the failed reality check). arXiv `physics.flu-dyn` (cross-list `physics.ed-ph`) is fine *after* the novelty fix — but the blockers above make even a preprint premature right now.

## Bottom line

Do **B4 first** (the Phares recompute — one afternoon). It tells you which headline numbers are real. Then resolve B1–B3 (the standoff/P2 derivations), decide B5/B6, add the citations, and *then* post. This is a real project worth finishing — it's just an afternoon or two of physics from posting, not a compile-and-upload. If you want, I can do the Phares cross-check computation now and hand you the corrected d\*/Q\*/lever numbers so you can see how much of the paper survives.
