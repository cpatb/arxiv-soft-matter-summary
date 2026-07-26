# FINISH FRISBEE — drive the Sean Earl disc-flight project to a finished paper

**Paste into Claude Code with `C:\Users\compl\Desktop\Projects\frisbee` as the working directory.**
Environment: Windows + WSL2 with **OpenFOAM v2606** installed (`source /usr/lib/openfoam/openfoam2606/etc/bashrc`), Python 3.10+ (numpy/scipy/matplotlib; `pip install pyvista trimesh rasterio requests` as needed), MiKTeX (`pdflatex`), internet access for the public geospatial data. Use your engineering/testing/code-review skills and git.

## What this project is (read first)
A first-of-its-kind coupling of **environmental CFD** (wind over Loyola's lakefront Sean Earl Field) to a validated **6-DOF flying-disc flight model**, producing a spatial "throwability" map. The manuscript (`Project1_Frisbee_LoyolaWind_draft.{tex,md,pdf}`) is a complete methodology skeleton — abstract, intro, physics, CFD method, coupling, metric all drafted — with **§3.7 Results, §3.8 Discussion, §3.9 Conclusions as placeholders "to be filled by simulation."** The code to fill them is the 5-stage pipeline in `SeanEarl_DiscFlight/` (P1→P5). `Project1_AnalyticalTreatment.{tex,md,pdf}` is the companion full derivation of the 6-DOF model. **Your job: run the pipeline to generate the data, fill the placeholder sections, and finish the paper.**

## The key strategic order (de-risk before the expensive CFD)
The one heavy step is P3 (OpenFOAM urban CFD sweep). Do NOT start there. **Validate the entire downstream pipeline first using the built-in stub wind field**, so that when real CFD data arrives you only swap the field source. Order: P1 → P4/P5 on a stub field (end-to-end smoke test) → P2 (geometry) → P3 (CFD, start with ONE direction) → swap real field into P4/P5 → scale P3 to the full sweep → fill paper.

---

## STAGE 1 — P1: finalize and validate the 6-DOF disc model
`cd SeanEarl_DiscFlight/P1_AnalyticalModel`
1. `python -m pip install numpy scipy matplotlib pytest`; run `python validation.py` and `pytest test_disc_6dof.py -v`. Confirm the 4 sanity tests pass (range, headwind, crosswind drift, spin decay) and the trim-glide is physical.
2. Address the P1 README's "what still needs doing": (a) confirm the linear CL/CD/CM model is used only within its stated α range; (b) if you can find/digitize the Potts–Crowther (2002) or Hummel (2003) tabulated coefficient curves, add them as an option and show the trajectory sensitivity; otherwise keep linear and state the α restriction in the paper. (c) Add at least one **quantitative** validation against a published throw reconstruction (Hummel 2003 Fig. 5-1) — even a digitized-points RMS comparison — and record it; this is what a referee will demand.
3. Commit. Output: a validated 6-DOF integrator with a documented validation figure.

## STAGE 2 — P4+P5 end-to-end on the STUB field (smoke test, no CFD yet)
Prove the coupling → Monte Carlo → throwability → figures chain runs before spending CFD time.
1. `cd ../P4_Coupling`; run the "Quick test (no CFD needed)" snippet in its README (StubField uniform wind). Then run a small `monte_carlo.py` sweep against a synthetic uniform/analytic field (write a tiny `StubField` npz if needed) to produce a `landing_stats.npz` with the documented schema.
2. `cd ../P5_Synthesis`; run `throwability.py` + `plot_throwability_map.py` + `paper_figures.py` on the stub outputs. Confirm Q(x,θ) maps and the 4 paper figures generate with correct axes/labels.
3. Fix any bugs in the chain now (interpolator, MC driver, metric normalization). Commit. **Gate: the whole downstream pipeline runs start-to-finish on fake data before P3.**

## STAGE 3 — P2: build the site geometry
`cd ../P2_Geometry`
1. Follow the P2 README to fetch public data (Microsoft GlobalMLBuildingFootprints, USGS 3DEP LiDAR DEM, OSM) for a ~500 m box around Sean Earl Field, extrude building heights, add the lakefront tree canopy as a porous zone, and export `seanEarl.stl` (+ ground DEM). If a data source is unreachable, substitute the nearest open equivalent and note it.
2. Sanity-check the STL against the Loyola campus map (named buildings: Mertz, Norville, Damen). Commit. Output: `seanEarl.stl` for P3.

## STAGE 4 — P3: the OpenFOAM wind CFD (the simulation data)
`cd ../P3_CFD`; `source /usr/lib/openfoam/openfoam2606/etc/bashrc`
1. **Version-fix the case first.** The `caseTemplate/` targets OpenFOAM v2406/11; you have **v2606** (ESI). `simpleFoam`, `kOmegaSST`, `snappyHexMesh`, and the Richards–Hoxey/`atmBoundaryLayerInletVelocity` ABL BCs all exist in v2606, but a dict keyword or BC name may have changed. Run one case, fix any dict errors (record them), and get `blockMesh`/`snappyHexMesh`/`simpleFoam` converging on a SINGLE representative direction first (suggest θ=90°, easterly over-lake — the paper's headline internal-boundary-layer case).
2. **Validate before trusting** (the paper promises this): (a) recover the ABL inlet profile in an empty upstream fetch (Hargreaves & Wright 2007 protocol); (b) run the Tominaga (2008) AIJ-B cube benchmark on the same mesh/BC stack and compare; (c) grid-independence: refine once, mean velocity at disc height within ~10%. Do (a) and (c) at minimum; (b) is the credibility anchor if time allows.
3. **Scope the sweep to compute reality.** The full 13-direction × 5–8 M-cell sweep is ~1–4 days on a 32-core box. Check the machine's cores (`nproc`); if it's a laptop, DON'T run all 13 at full resolution. Instead run a **reduced, decision-useful set** first: the three regimes the paper hypothesizes — easterly over-lake (~90°), westerly urban (~270° if you extend past 180°, else the 0°/180° endpoints), and a prevailing SSW-ish case — at a coarser mesh (2–4 M cells). `generateCases.py` + `runSweep.sh` drive it. Log exactly which directions/resolution you ran; the paper must state the subset honestly.
4. `python sample_disc_altitude.py cases/dir_*/ --out flow_fields.npz` → the disc-altitude flow field. Save a CFD mean-flow figure (paper Fig. 3). Commit.

## STAGE 5 — swap real CFD into P4, run the real Monte Carlo
`cd ../P4_Coupling`
1. `python monte_carlo.py --field ../P3_CFD/flow_fields.npz ...` (see README) over the launch grid × the CFD headings you actually ran × the three throw types. Start with `--n-throws 1000` to shake it out, then raise if time allows. Output: real `landing_stats.npz`.
2. Note the P4 README's honest caveats: the synthetic-turbulence gust generator is a white-noise first pass (fine for v1); state this limitation in the paper.

## STAGE 6 — P5: throwability maps + paper figures on real data
`cd ../P5_Synthesis`
1. `throwability.py` → `Q.npz`; `plot_throwability_map.py` and `paper_figures.py` → the four manuscript figures on REAL data (site map, trim-glide trajectories, CFD mean-flow, throwability heatmap best-vs-worst heading).
2. Check the results against the draft §3.7 hypotheses (easterly Q-gradient across the field, sheltered westerly throw-lanes, best-condition ~light SSW). Report where they hold and where the CFD says otherwise — divergence from the hypothesis is a finding, not a failure.

## STAGE 7 — fill and finish the manuscript
Edit `Project1_Frisbee_LoyolaWind_draft.tex` (keep a copy of the original):
1. **Write §3.7 Results** from the P5 outputs — the throwability maps, the regime comparison, quantitative Q values, which throws/positions win under which wind. Insert the four real figures.
2. **Write §3.8 Discussion** (generalizability, steady-RANS limits — no peak gusts, the white-noise turbulence caveat, coaching implications) and **§3.9 Conclusions**.
3. **Resolve the §5 open questions** with sensible defaults, flagging each for the author: site = Sean Earl (assumed); LES = out of scope, noted as follow-on; anemometer validation = state as future field-campaign unless real data exists; hammer throw = keep as "exploratory" (least-validated coefficients); fill author/affiliation/co-authors (`% FLAG FOR AUTHOR`).
4. State the **honest scope of what was actually simulated** (the CFD direction subset and mesh resolution) — do not imply the full 13-direction sweep if you ran fewer.
5. Verify the §4 references resolve; keep the paper's own "verified references" list.

## STAGE 8 — review, compile, finalize
1. Run the multi-tier review harness: `Read ../CLAUDE_CODE_PREPRINT_REVIEW_PIPELINE.md and follow it, TARGET = frisbee/Project1_Frisbee_LoyolaWind_draft.tex`. Resolve blockers.
2. `pdflatex` twice; confirm it compiles, figures resolve, no undefined refs/placeholders. Sane page count.
3. Write `FRISBEE_STATUS.md` at the frisbee/ top level: what ran vs what was scoped down, the headline throwability findings, which draft sections are now complete, and the remaining author decisions + any second-pass compute (full 13-direction sweep, LES spot-check) worth doing for journal submission.
4. Commit throughout so the pipeline runs, the data, and the manuscript edits are separable.

## Venue (from the draft)
*Sports Engineering* (Springer) for the coupled paper; or split into a disc-trajectory paper (AIAA Applied Aero) + an urban-CFD paper (*JWEIA*). Decide at the end based on the strength of the result.

## Guardrails
- Don't fabricate results — if you scoped the CFD down, the paper says so, and §3.7 reports only what was computed.
- Validate the CFD (ABL recovery + grid independence at minimum) before trusting any flow field.
- Preserve original drafts and the pipeline; work on copies where you edit manuscripts.
- Every headline number traces to a converged, validated run or a passing test.
- Flag, don't invent: author identity, co-authors, and any real field measurements are the author's to supply.
