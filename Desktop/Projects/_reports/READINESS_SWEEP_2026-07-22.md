# Publication-Readiness Sweep — `Desktop/Projects`

**Date:** 2026-07-22 · **Scope:** every project directory and standalone paper in `C:\Users\compl\Desktop\Projects`
**Mode:** read-and-assess. **No manuscript, code, or data file was modified.** The only files written are the 15 per-project `READINESS.md` files and this report.
**Basis:** prior `AUDIT*.md` / `README.md` / `REVIEW.md` / `OVERNIGHT_*.md` treated as ground truth and re-checked, plus fresh reproduction runs where fast and safe.

## Inventory reconciliation vs the brief

- **New since the brief was written:** `axidrop/` — a Basilisk drop-impact solver harness created within the last day (files dated 2026-07-23 by the local clock). Assessed; it is a scaffold with no solver code.
- **Present but not itemised in the brief:** `hydrogels/drying_analysis/` (validated measurement package, 2026-07-06) and `hydrogels/coumarin_penetration/` (collaborator analysis note, 2026-07-08). Both are material to the hydrogels assessment and are covered in its `READINESS.md`.
- **Nothing in the brief's inventory is missing.** All listed projects exist.
- `slime mold/bike/` is scored separately from the rail project below, because it is a different deliverable with a different natural venue.

## Reproductions run tonight (all PASS)

| Project | Command | Result |
|---|---|---|
| `tier0_wyart_cates_closure` | `run_sanity_gate.py`; `test_wyart_cates.py` | Gate 1, Gate 2, Numerics **PASS**; **6/6** invariant tests |
| `bio_fluids_bridge` | `run_blood_bridge.py` | **PASS** (5/5 gates) |
| `thrombosis_jamming` | `run_thrombosis.py` | Signature **PASS** (threshold 1475 s⁻¹, t_occ 0.37–103.7 s) |
| `glymphatic_transport` | `glymphatic_peclet.py` | **PASS**, matches field synthesis |
| `bloodstain_bpa` | `test_bpa.py` | **6/6 passed** |
| `hydrogels` | `python -m drying_analysis.validate` | **36/36 checks passed** |
| `uroflow_splashback` | `run_uroflow.py` | Runs clean; ~30× splash-index span reproduces |

Not re-run (graded from prior audit + code, as the brief permits): `slime mold` (multi-seed sweep, reproduced exactly on 2026-06-30), `diff_calibrate.py` (torch startup), the napari GUI scripts (no footage present), and all LaTeX compiles except where a log was inspected.

---

## 1. Executive ranking table

Sorted by closest to a publishable next step first.

| Rank | Project | Best stage | Score | One-line: what it needs |
|---|---|---|---|---|
| 1 | **bidet/** (optimization + policy) | Preprint | **4** | Fill author block; fix the 164 pt table; reconcile three stale numbers. One evening. |
| 2 | **slime mold/** (rail) | Preprint | **4** | Replace `\author{Prepared for Charlie}`; footnote the networkx Steiner drift. One hour. |
| 3 | **friday_experiments.tex** | Exp. validation | **5** | Nothing blocking — run Experiment 1. Both physics errors already corrected. |
| 4 | **thrombosis_jamming/** | Exp. validation | **5** | Verify two citations and prior art, then digitize published occlusion-vs-shear data. Free. |
| 5 | **paper1/** | Exp. validation | **4** | Reconcile the ε_jam threshold contradiction before running the August experiments. |
| 6 | **hydrogels/** | Exp. validation | **4** | One 10-minute transmission measurement settles the coumarin question; then shoot one drying run. |
| 7 | **tier0_wyart_cates_closure/** | Exp. validation | **4** | Measure real corn-starch and silica flow curves; the autodiff calibrator is already built. |
| 8 | **uroflow_splashback/** | Exp. validation | **4** | One phone-camera breakup-length measurement calibrates the single free parameter. |
| 9 | **slime mold/bike/** | Preprint / *Findings* | **3** | Add an author block; restate modeled-vs-measured in the body, not just in `REVIEW_bike.md`. |
| 10 | **bloodstain_bpa/** | Preprint | **3** | **Check Attinger prior art first** — that one hour decides whether this is a paper. |
| 11 | **chicago_run_bike_monitor.TEX** | Exp. validation | **3** | Reframe H2 from crosswind σ_y to line-source σ_z *before* collecting transects. |
| 12 | **bio_fluids_bridge/** | Exp. validation | **3** | Calibrate against published blood viscometry. Weak alone; essential to the combined paper. |
| 13 | **glymphatic_transport/** | (Exp. validation) | **2** | Not a publication unit. A field-entry decision, not a readiness gap. |
| 14 | **LDV Paper/** | — | **1** | No manuscript exists. Starting one is a writing project, not a fix. |
| 15 | **axidrop/** | — | **0** | Scaffold only; run Rung 0. Correctly months from any paper. |
| — | **ARCHIVE/SCM/** | — | **0** | Leave archived. |

---

## 2. Three shortlists

### A. Post as preprint now / this week (preprint score ≥4)

**`bidet/`** — two complete, compiling, numerically-verified revtex manuscripts. Canonical files are `bidet_optimization_1.tex` and `bidet_paper_1.tex`. Blockers are entirely mechanical: author placeholders at lines 69–71 of the optimization paper, a 164 pt overfull table in the policy paper, a retracted 16 mm value still feeding the Péclet estimate at lines 685/689, and an 8 ms / 3.4 ms mismatch. *Caveat:* the optimization paper is the arXiv-shaped object. The policy paper is an advocacy/review document making health and economic claims — posting it is a separate judgment call, not a technical one.

**`slime mold/` (rail)** — a 9-page, 5×-reviewed manuscript whose headline reproduces to 2 d.p. Blocked by a literal placeholder author line and one documented numeric drift (under networkx 3.6.1 the Steiner baseline now coincides exactly with the MST, against the paper's reported 145.2 km). Conclusions are unaffected, but a reader who re-runs will hit it.

### B. Ready to design the validating experiment (exp-validation score ≥3, with the one measurement named)

| Project | The one measurement | Falsifier |
|---|---|---|
| **friday_experiments** (5) | R(t) for a fixed-volume radial dam break, log–log, with a Newtonian control at matched η₀ | Newtonian control must give R ∝ t^{1/8}; if not, the current isn't in the viscous-gravity regime and the baseline is invalid |
| **thrombosis_jamming** (5) | t_occ vs wall shear rate, 500–20,000 s⁻¹ — digitized from published Ku-group figures first, free | Any reproducible occlusion **below** γ̇_crit ≈ 1500 s⁻¹ kills the jamming picture outright |
| **paper1** (4) | t_jam for TPM vs cornstarch at **matched φ_eff and matched bulk rheology**, high-speed at 50 µm/px | Same t_jam and same relaxation within error ⇒ the "bulk rheology is insufficient" claim is dead |
| **hydrogels** (4) | Direct 365 nm transmission through a 2 mm gel on a ThorLabs power sensor, ~10 min — gives δ₁/ₑ with no assumption about ε₃₆₅ | Transmission flat during irradiation ⇒ non-bleaching absorber floor, a different diagnosis than the model assumes |
| **tier0_wyart_cates_closure** (4) | Steady flow curves for your own corn starch and silica at 2–3 φ each, roughened tools | σ* drifting systematically with φ ⇒ the closure is wrong for this system, not merely mis-tuned |
| **uroflow_splashback** (4) | Intact-stream length L_b at 240 fps against a ruled backdrop, ~30 min | L_b flat or non-monotonic in v ⇒ turbulent-disturbance-dominated, not Rayleigh; the coherence framing is wrong physics |
| **bloodstain_bpa** (3) | A documented spatter pattern with an independently surveyed origin | Straight-line method *not* biased high on real patterns ⇒ the synthetic emission geometry produced the result, not the physics |
| **chicago_run_bike_monitor** (3) | PM2.5 perpendicular transects 0–300 m from a highway, neutral conditions, with co-location calibration | Identical fall-off at a no-line-source control site ⇒ background gradient or drift, not plume dispersion |
| **bio_fluids_bridge** (3) | Digitized whole-blood flow curve + one up/down ramp; fit steady curve and hysteresis with one parameter set | γ̇_c needing free refitting per haematocrit ⇒ descriptive curve-fitting, not a closure |

### C. Push toward journal submission (journal score ≥3)

Only three units clear the bar, and all three clear it barely.

- **`bidet/` optimization paper — 3.** Target: *Physics of Fluids* (or *European Journal of Physics* for the pedagogical framing). **The gap:** novelty is not positioned, and the load-bearing empirical unknown — the fecal detachment shear τ_c — is unmeasured. Worth noting in the paper's defence: Q* ∝ τ_c^{1/3}, so even a factor-of-3 error in τ_c moves the optimum by only ~1.4×. Say that explicitly; it converts the biggest referee objection into a robustness result.
- **`slime mold/` rail — 3.** Target: *Physica A*, *EPJ Data Science*, or *Journal of Transport Geography*. **The gap:** Physarum-applied-to-real-transit-networks is a crowded post-Tero literature. The defensible novelty is the *combination* of a hard no-build constraint (Lake Michigan), a construction-cost field folded into effective tube length, and transit-desert-weighted demand — none of which is currently argued as novel in the introduction. Plus demand is normalized, not ridership-calibrated.
- **`slime mold/bike/` — 3.** Target: *Transport Findings* (short format), which this brief nearly fits as-is. **The gap:** author block, and the modeled-vs-measured distinction currently lives in `REVIEW_bike.md` rather than in the brief's own body.

**Nothing in this folder is journal-ready at 4 or 5.** The reason is uniform and is stated in §4.

---

## 3. Per-project rationale

**bidet/ — Preprint 4, Exp 3, Journal 3.** Two matched revtex manuscripts that compile cleanly and whose numbers independently recompute (V₀ ≈ 11.48 m/s, d* ≈ 0.80 mm, Q* ≈ 0.347 L/min, r_clean ≈ 26 mm), with core dimensional consistency verified including a hydraulic-jump scaling that correctly replaced an earlier dimensionally-broken Watson gravity formula. The papers self-flag their weakest points honestly — the laminar Region-II asymptotic as a lower bound, and the awkwardness that the standoff-corrected 1.6 L/min prediction exceeds commercial flow rates. Every remaining blocker is named, bounded, and mechanical, which is exactly what a 4 means.

**slime mold/ (rail) — Preprint 4, Exp 0 (N/A), Journal 3.** The most thoroughly verified project here: five independent review passes covering literature fidelity, numerics (Kirchhoff residual exactly 0.000), geography constraints, figures-vs-numbers, and reproducibility, with three real issues found and fixed. The headline reproduced exactly on re-run. It is held back from arXiv by a placeholder author line, which is close to the purest possible example of a score-4 blocker. Journal is capped at 3 by novelty positioning in a crowded literature, not by any defect.

**friday_experiments.tex — Preprint 2, Exp 5, Journal 1.** The standout finding of this sweep: **the two physics errors flagged in the June audit have been fully corrected.** The axisymmetric constant-volume spreading law now reads R ∝ t^{1/8} with a dated in-file correction comment and a Huppert 1982 citation at every one of the three flagged sites, and the 2D height decay now reads t^{−1/5}, consistent with the paper's own X ∝ t^{1/5}. The downstream liquid-dam comparison was updated too. This is the only place in the folder where a flagged physics error was closed, and it was closed properly. Preprint stays at 2 only because the document is a proposal, not a results paper.

**thrombosis_jamming/ — Preprint 2, Exp 5, Journal 1.** It has what nothing else here has: a genuine mechanism-discriminating prediction. Incumbent thrombosis models produce arrest by deposition kinetics — smooth, monotonic, no threshold. This predicts a critical point: no occlusion below γ̇_crit, and t_occ diverging above it. That is falsifiable, cheap to test against published data, and unambiguous about what would kill it. The two risks are equally real and both self-flagged: the prior-art check is unresolved, and the data may be too sparse to fit an exponent. Scout the figures before writing equations.

**paper1/ — Preprint 1, Exp 4, Journal 1.** A rigorous, complete experimental design whose physics was independently verified in June (Pe_impact, Stokes–Einstein D_B, τ_B, Pe_unjam, porosity correction, Debye lengths, τ_ic all check out). But **none of the four June flags has been addressed** — the files have not been touched since May. The ε_jam contradiction is the one that matters: the body sets ε_jam = 0.07 and calls anything below it unresolvable, while Objection 5 still promises to test ε = 0.02, 0.03, 0.05. Since ε_jam defines t_jam, which is the headline observable, this must be settled before the August experiments run, not after.

**hydrogels/ — Preprint 1, Exp 4, Journal 1.** The June audit's CRITICAL S1 finding (the pipeline measured wrinkle morphology, not the drying-fracture observables the roadmap needs) was **properly resolved** by building a separate `drying_analysis/` package rather than by contorting the legacy scripts. That package re-verified at 36/36 tonight, with convention cross-checks and discriminating cases designed so a cancelling pair of bugs cannot pass silently. The entire program is now blocked on data, not code — and the coumarin sub-project has reduced its central question to a single ~10-minute measurement.

**tier0_wyart_cates_closure/ — Preprint 2, Exp 4, Journal 1.** Clean, tested, reproducing, with a working autodiff calibrator that recovers ground-truth parameters from noisy synthetic viscometry. The honest limit is stated in its own README: every parameter is illustrative. The gates prove the model *can* produce a DST-vs-CST contrast, not that these materials do. One flow-curve dataset changes that, and the machinery to consume it already exists.

**uroflow_splashback/ — Preprint 2, Exp 4, Journal 1.** Genuinely validated physics (the viscous dispersion relation reduces to Rayleigh's inviscid result within 2.5%), reaching conclusions that the Splash Lab already published experimentally. Agreement is simultaneously a validation success and a novelty liability. Two cheap measurements — real breakup length, real splash mass — would convert a ranking model into a calibrated predictive one, which is the difference between not publishable and publishable.

**slime mold/bike/ — Preprint 3, Exp 0 (N/A), Journal 3.** Every measured number traces to a named City of Chicago dataset ID with recomputation checks, and the modeled-vs-measured boundary is drawn explicitly. Its handicap is format, not quality: it is a policy brief without an author block or paper structure. *Transport Findings* is a real, plausible, short-format home.

**bloodstain_bpa/ — Preprint 3, Exp 3, Journal 2.** The most publication-shaped *result* among the code projects: a quantified bias (straight-line stringing overestimates origin height by 1.3–2.2×, growing with impact speed), a mechanism, a 230× improvement, and an honest noise study. Two things hold it back, one fatal if unlucky. The prior-art risk is real — Attinger's group has published extensively on BPA trajectory reconstruction including drag and gravity, and their review is cited in this project's own README. One hour of searching decides whether this is a contribution or a replication.

**chicago_run_bike_monitor.TEX — Preprint 1, Exp 3, Journal 1.** Both fatal LaTeX errors from June are fixed and it compiles. The σ_y-vs-σ_z problem is not fixed, and it is worse than a documentation error: the analysis plan pipes a crosswind Gaussian into `curve_fit` for a perpendicular transect from a line source, which will return a number that is not σ_y, not σ_z, and not interpretable. Fixing it after a field campaign means re-running transects in matched meteorology.

**bio_fluids_bridge/ — Preprint 1, Exp 3, Journal 1.** Validated and correct, but reproducing shear thinning from a Moore–Coussot structural model is a sanity check against a mature literature, not a finding. Its value is as one leg of the transferability story. Also the only one of the three closure projects with no test suite.

**glymphatic_transport/ — Preprint 1, Exp 2, Journal 0.** Correct, well-made, and validated by *agreeing with the consensus* — which is precisely why it cannot yield a testable new claim. The genuinely contested physics (oscillatory Taylor–Aris dispersion in eccentric perivascular annuli) is explicitly excluded as "an uncertainty band, not modeled precisely here." That exclusion is where a contribution would live.

**LDV Paper/ — 1/0/0.** A reference shelf and two talks, static since June 2025. The 29MAY25 deck contains a coherent tutorial arc, so the skeleton of a methods tutorial exists — but writing it is a project, not a fix.

**axidrop/ — 0/0/0.** A scaffold created in the last day, with an all-TBD validation spec and an empty `src/`. Judged against publication bars it scores zero; judged as infrastructure its "validated" definition (named source + reproducing case + stated tolerance + one command, or else say "implemented") is the most rigorous standard in this folder and worth propagating to the other projects.

**ARCHIVE/SCM/ — 0/0/0.** Correctly archived, untouched since October 2025. Noted and skipped per the brief.

---

## 4. Cross-cutting observations

### 4.1 The "one closure, many fluids" cluster — write ONE paper, not three

`tier0_wyart_cates_closure` → `bio_fluids_bridge` → `thrombosis_jamming` are three validated, README-only projects sharing one state-variable structural closure with the microphysics relabelled: frictional contacts (silica, thickening), rouleaux aggregation (blood, thinning), VWF-bridged adhesive contacts (platelets, jamming).

**Recommendation: one methods paper, not three notes.** Individually, each is thin enough that a referee will call it a relabeling — and for `bio_fluids_bridge` that criticism is basically correct, since structural thixotropic blood models are mature. Together, the *transferability* is the claim, and `thrombosis_jamming` supplies the falsifiable prediction that gives the other two their stakes. Plausible venues: *Journal of Rheology*, *Soft Matter*, or *PRE*.

**The blocker is common to all three, which is what makes this efficient:** every parameter in all three projects is illustrative rather than calibrated. **One calibration pass unblocks all three at once** — tier0 against real flow curves, blood against digitized published viscometry, thrombosis against digitized Ku-group occlusion data. This is the single highest-leverage action in the folder measured by projects-moved-per-hour, and it is also on `axidrop`'s critical path (Rung 3 couples the tier0 closure).

### 4.2 The bottleneck is writing, not science

Eight of fifteen units have validated, reproducing results and **zero manuscript prose**. Every code project's preprint score is capped at 2–3 for one reason only: no manuscript exists. There is no science blocker on `tier0`, `bio_fluids_bridge`, `thrombosis_jamming`, `bloodstain_bpa`, `uroflow_splashback`, or `glymphatic_transport` — there is a writing blocker. Conversely, the two projects nearest to posting (`bidet`, `slime mold`) are blocked by an author block and a table width. **The folder's constraint is prose and polish, not physics.**

### 4.3 Nothing here has touched real data

Every computational project self-reports internal validation only — round-trip, analytic-limit, or synthetic ground truth. That is honest and consistently disclosed, and the validation is genuinely good (the Rayleigh 2.5% check, the Kirchhoff residual, the 36/36 discriminating-case suite). But it means **no journal score in this folder exceeds 3, and none can until a real measurement enters somewhere.** Six of the nine bar-2 entries in shortlist B need only a phone camera, a power meter, a rheometer session, or a literature digitization. This is the cheapest structural upgrade available.

### 4.4 Version hygiene invites editing the wrong file

- `bidet/`: 8 `.tex` files for 2 papers (originals, `_pre-review`, `_test_` article-class copies), plus 4 empty 0-byte `*Notes.bib` files that nothing consumes.
- `paper1/`: `paper1_referee_review.sty` is a **byte-identical duplicate** of `paper1_referee_review.tex` (both 13,349 bytes), never `\input` or `\usepackage`d — plus a zero-byte extensionless `paper1_referee_review` file.
- Top level: ~20 loose LaTeX build artifacts.

None is dangerous; all are ways to lose an evening editing a superseded file. Recommend an `archive/` subfolder per project and a `% CANONICAL` comment at the top of each live source. **Stage 2's `.EDITED` convention makes this more urgent, not less** — adding a third naming axis to a folder that already has four versions per paper needs the superseded ones moved aside.

### 4.5 Prior audits: a mixed but improving record

Of the flags carried forward, one project closed its physics flags properly (`friday_experiments`, both exponents corrected in-file with a citation and a dated comment), one closed its CRITICAL architectural flag by building the right thing instead of patching the wrong thing (`hydrogels` S1 → `drying_analysis`), and one closed nothing at all (`paper1`, all four flags open, files untouched since May). `bidet` closed both fatal/blocking items and left all five cosmetic-to-medium ones open. **The pattern: fatal and structural issues get fixed; medium-severity consistency issues accumulate.** Those accumulated medium issues are precisely what now blocks the two highest-ranked projects from posting.

### 4.6 The Chicago cluster has no home

`slime mold/` (rail), `slime mold/bike/`, and `chicago_run_bike_monitor` form a coherent urban-systems track with real data and real deliverables — and no connection to the soft-matter/fluids identity that every other project serves. This is a portfolio-strategy question, not a readiness question, but it is worth deciding deliberately: either commit to it as a genuine second track, or post the two finished pieces and stop.

### 4.7 External dependency worth confirming (unverified from this tree)

The drop-impact **data** pipeline that would populate `paper1/` lives outside `Desktop/Projects`. A normalization correction reported on 2026-07-22 changes shipped under-view d/d₀ values by a factor of ~2.14, with reprocessing performed. **This could not be verified from within this folder** and is recorded only so it is not forgotten: if accurate, any spreading number from that pipeline needs re-checking before it enters the paper1 manuscript.

---

## 5. Do this first

> **If you have one evening: `bidet/`. Fill the author block, fix the 164 pt table, reconcile the three stale numbers, and post the optimization paper to arXiv.**

**Why this and not something else.** It is the only action that converts a finished artifact into a public one in a single sitting. Two manuscripts already compile, their numbers already reproduce, their dimensional consistency is already verified, and every remaining blocker is named and mechanical — no physics decisions, no new analysis, no data. Every other candidate needs either a measurement (which needs daylight and a lab) or net-new prose (which needs more than an evening). Post the optimization paper; treat the policy paper as a separate decision.

**Runner-up, and it is close: `paper1`'s ε_jam reconciliation.** It yields no publication, so it loses on the stated criterion — but it has the highest *scientific* value per hour in the folder. The August experiments are designed around a threshold definition that currently contradicts itself in two places, and ε_jam defines t_jam, the headline observable. An evening spent settling it now prevents executing a compromised protocol later. **If the August timeline is real, do this one first and bidet second.**

**Third: the one-hour Attinger prior-art check on `bloodstain_bpa`.** Cheapest decisive action available. It either promotes a well-built project to "write it up" or retires it, and right now it is consuming shelf space in an undecided state.
