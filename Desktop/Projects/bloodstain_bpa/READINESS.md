# Readiness — bloodstain_bpa

**Assessed:** 2026-07-22  ·  **Canonical artifact:** code — `bpa.py`, `run_experiments.py`, `test_bpa.py`; figures `fig1`–`fig4`
**One-line status:** The strongest quantified *result* in the folder that has no manuscript — and the one whose prior-art risk is most likely to be fatal if not checked first.

## Scores (0–5)

| Stage | Score | Primary blocker to next point |
|---|---|---|
| Preprint / arXiv | 3 | A real quantified result with figures exists, but no prose at all, and the headline bias magnitude is admitted to be geometry-dependent. |
| Experimental validation | 3 | The validating test (a documented spatter pattern with known origin) exists in the forensic literature but has not been located or run. |
| Journal submission | 2 | Prior art unchecked — Attinger's group has published extensively on exactly this problem, including drag and gravity corrections. |

**Recommended next stage:** Preprint / arXiv  ·  **Why:** it has the one thing every other code project here lacks — a specific, quantified, non-obvious finding with a bias magnitude, a mechanism, and a noise study — but see blocker 1 before writing a word.

## What exists (evidence)

- **Re-verified tonight: `python test_bpa.py` → 6/6 passed.** Free-fall vs analytic, forward/backward reversibility, Balthazard angle, spread-model monotonicity and invertibility, round-trip origin recovery, spread increases with speed.
- **The finding:** the classical forensic straight-line ("stringing") method **overestimates height of origin by ~1.3–2.2×** for impact-spatter-onto-floor geometry, and the bias **grows with impact speed** (flatter trajectories accumulate more gravitational curvature). Mean bias +157% straight-line vs +0.7% ballistic — a ~230× reduction. Under realistic measurement noise (5% ellipse, 3° azimuth, 10% drop size) the ballistic estimate keeps an honest −18 ± 5.5 cm bias, driven almost entirely by drop-size uncertainty, against the straight-line method's +212 cm.
- Four current figures: scene, bias sweep, drop-impact spread bridge, noise study. `run_output.txt` records a run.
- Method anchors are real and correctly attributed: Balthazard (1939) angle relation, Pasandideh-Fard et al. (1996) max-spread correlation, Schiller–Naumann (1935) drag, Attinger et al. (2013) BPA review.
- Library-first and reproducible per lab conventions: `scipy.solve_ivp`, `brentq`, seeded `default_rng`.

## Blockers to the recommended next stage (ranked)

1. **[CRITICAL — check this before writing anything]** **Prior art is unchecked and the risk is high.** Daniel Attinger's group has published repeatedly on the fluid dynamics of bloodstain pattern analysis, explicitly including the effect of drag and gravity on trajectory reconstruction and the resulting bias in region-of-origin estimates. The 2013 review is already cited *in this project's own README*. If the quantified bias is already in that literature, this is a replication, not a finding. **One hour of literature search decides whether this project is a paper or a well-built exercise.**
2. **[HIGH]** **No manuscript exists.** Not a draft, not an outline — the README is the only prose.
3. **[HIGH]** **The headline number is geometry-dependent and the project says so.** "The exact bias magnitude is geometry-dependent (emission cone, floor vs. wall, source height). The direction and mechanism are robust; the percentage is illustrative of this scenario, not a universal constant." A referee will ask for the bias mapped across geometries, not one scenario. That sweep is cheap to run and would convert the weakest part of the claim into the strongest.
4. **[MEDIUM]** **Parent drop diameter D₀ is the load-bearing unknown.** It is known by construction in the synthetic setup but uncertain in real casework, and the noise study shows drop-size uncertainty dominates the residual bias. The included spine-count cross-check (`n ≈ 1.14·We^{1/2}`) is acknowledged as noisy. Any real-world claim stands or falls on this.
5. **[MEDIUM]** Never run against real data — validation is round-trip and analytic-limit only, by the project's own statement.
6. **[LOW]** Spread model is the smooth-surface correlation; real substrates (fabric, textured drywall) deviate.

## Carried-forward flags from prior audits (status)

- No prior formal audit exists (built 2026-07-08). `OVERNIGHT_fluids_2026-07-08.md` is the build record.
- `OVERNIGHT_fluids_2026-07-08.md` claim "6/6 correctness tests pass" — **verified tonight, still true.**
- That note's closing statement — "Nothing here has been run against real data yet. The validation is internal (round-trip / analytic-limit) and the next step for both is a single real measurement to calibrate" — **still true, unchanged.**
- Its four named next steps (wall+floor combined patterns; Monte-Carlo posterior over the origin region; validate against a public BPA benchmark set; replace the analytic spread law with one fit to your own drop-impact simulation data) — **all still open.**

## Concrete next actions (shortest path to preprint)

1. **Literature check first (1 hour, gates everything else).** Search Attinger, Comiskey, Yarin on BPA trajectory reconstruction with drag/gravity. Determine whether the quantified straight-line bias is already published. Do not write until this is answered — the answer determines whether the paper is "we quantify a known bias" (weak), "we quantify an unquantified bias" (good), or "this is already done" (stop).
2. **Run the geometry sweep.** Map the bias across emission cone angle, source height, and floor-vs-wall target. This is a loop over existing code and it converts the project's most-criticised claim into a defensible map. Do this before writing regardless of what step 1 finds.
3. **Replace the point estimate with a Monte-Carlo posterior** over the region of origin (README's own item 2). Forensic reviewers care about admissibility and uncertainty far more than about a better point estimate; a posterior is the contribution that community would actually value.
4. Write the manuscript. Target *Forensic Science International* or *Journal of Forensic Sciences*; the drop-impact-spread bridge gives it a legitimate second home in the fluids literature.
5. *If pursuing experimental validation instead:* the one measurement is a **documented experimental spatter pattern with an independently known origin** — controlled-height blood-drop or impact-spatter sets from a forensic teaching lab, with the source position measured directly. **Observable:** reconstructed height of origin vs the surveyed true height. **Expected signal:** straight-line reconstruction high by 30–120%; ballistic within noise. **Control that would falsify it:** if the straight-line method is *not* biased high on real patterns, the synthetic emission-cone geometry — not the physics — is producing the result.
