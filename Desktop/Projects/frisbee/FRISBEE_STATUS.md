# FRISBEE_STATUS — what ran, what was scoped down, what remains

Generated at the end of the `FINISH_FRISBEE.md` execution.
**Results sections are filled from the pipeline artifacts; see "Provenance".**

---

## 1. Headline: what this project now is

A coupled disc-flight / urban-CFD study of Sean Earl Field with a **validated
6-DOF disc model**, a **real site geometry built from public data**, and a
**deliberately scoped three-direction RANS sweep**. The disc model is the
strongest part and is quantitatively validated against published flights. The
CFD is honest but limited, and the limitations are stated in the paper rather
than papered over.

The single most important scientific outcome is a **contradiction of the
draft's central premise** — see §4.

---

## 2. What actually ran

### Stage 1 — P1, 6-DOF disc model ✅ COMPLETE

- `validation.py`, `validation_hummel.py`, `pytest test_disc_6dof.py` — **26 tests passing**.
- Validated against two independently published, fully specified results in
  Hummel (2003): the §2.4.3 / Fig. 2-9 "common throw" and the §2.4.2 analytic
  trim glide.

| Quantity | Hummel (2003) | This model | Error |
|---|---|---|---|
| Range | ~25 m | 26.99 m | +8% |
| Apex height | 2.55 m | 2.45 m | −4% |
| Flight time | 3.5 s | 3.48 s | −1% |
| Terminal speed | ~V₀/3 | 4.83 m/s | +4% |
| Lateral reversal | ~3 s | 3.34 s | ✓ |
| Trim α (C_M=0) | 9° | 9.00° | <0.1% |
| Glide angle γ | −19.3° | −19.32° | <0.2% |
| C_L at trim | 0.67 | 0.672 | +0.3% |
| Trim speed | 9.1 m/s | 8.28 m/s | −9% |

**Blocker fixed first:** `disc_6dof.py` contained a duplicated fragment after
its final `return`, making the module a syntax error — nothing in P1 imported.

**Five substantive model errors found and corrected** against Hummel's own
MATLAB (thesis Appendix B), which is the code the Table 3-2 coefficients were
*fitted with*:

1. `CL` and `CM` used `(α − α₀)`; Hummel uses `α`. α₀ belongs only to the drag
   parabola. This moved the zero-pitching-moment angle from 8.8° to 4.8° — the
   number that sets where turn flips to fade.
2. α₀ was hardcoded to −4°; it is derived, `−CL0/CLα`, i.e. **−9.9°** for the
   Hummel short-flight set.
3. Rate-damping terms were non-dimensionalised by `d/(2V)`; Hummel's `CMq`,
   `CRp` multiply raw rad/s and `CRr` carries an explicit `√(d/g)` scale. The
   normalisation under-damped pitch and roll by ~2 orders of magnitude.
4. `CRr`: Table 3-2 lists `1.4e-2`, but Appendix A carries `0.00171` *"for
   nondimensionalization = sqrt(d/g)"*. Only the latter matches the form used.
5. Moments of inertia used the thin-uniform-disc idealisation, **33% and 35%
   below** Hummel's measured Table 3-1 values. Precession goes as
   `M/(I_s ω_z)`, so this overstated turn and fade by roughly 50%.

Also: the tabulated Potts & Crowther CL/CD/CM curves and the
`CRr(α, advance ratio)` surface were transcribed from the thesis appendix and
are now the default, because only the tabulated set reproduces the published
flights — the linear fits miss range/apex/time by 30–57% and never fade back
left, since a realistic throw reaches α ≈ 33° in its final second.

### Stage 2 — P4 + P5 on the stub field ✅ COMPLETE (the de-risking gate)

18 tests passing. The gate caught three defects that would each have
invalidated or blocked the real sweep:

1. **The white-noise gust generator was unusable.** It drew fresh randoms
   inside every RHS evaluation, so the ODE was non-deterministic; the adaptive
   solver's error estimate measured noise rather than truncation error and
   shrank the step indefinitely. **>100 s for one throw** vs ~0.1 s without.
   A 10⁴-throw sweep would never have finished, and the landing scatter would
   have been a step-size artefact rather than a solution. Replaced with a
   Shinozuka–Deodatis spectral representation (autocorrelation matches
   `exp(−Δt/T_L)` to within 0.06); gusts now cost essentially nothing.
2. **`success_rate` was in the documented schema but never computed**, so P5
   fell back to a proxy `min(1, (r/ā)²)` that is not a probability — it ignores
   the offset between landing centroid and target, scoring a tight but badly
   aimed group as a success. Now measured directly.
3. **Field interpolation dominated the sweep** (5 scipy interpolator calls per
   RHS evaluation). Replaced with hand-rolled trilinear, asserted bit-equal.

P1 speed work needed to make the sweep feasible, results unchanged: DOP853 at
rtol=1e-5 after a convergence study against a rtol=1e-9 reference (landing
converged to ~5e-6 m), pre-sampled coefficient splines, optional diagnostics.
**866 → ~140 ms per throw** with an identical landing point.

### Stage 3 — P2, site geometry ✅ COMPLETE

The P2 scripts were pure stubs that printed `[stub] would fetch …`. Implemented
against real public data.

- **`seanEarl.stl`**: 116 buildings, 6 324 triangles, all closed (zero open
  edges), tallest H = 60.8 m (Mertz Hall).
- Source: **OpenStreetMap via Overpass** (not Microsoft GlobalMLBuildingFootprints
  — OSM has the *named* Loyola buildings needed for the campus-map cross-check
  and carries height tags). 171/303 buildings (56%) have a `building:levels`
  tag; the rest use documented type defaults, each recording which rule applied.
- Clipped to a **circle** (r ≤ 300 m) so the modelled building set is invariant
  under the rotation used to realise each wind direction.
- Terrain from **USGS 3DEP** point queries: 179.1–181.4 m MSL across the field
  (±1 m) against 12–61 m buildings, so ground is modelled flat.

### Stage 4 — P3, OpenFOAM CFD ✅ COMPLETE (scoped — see §3)

Ported to **OpenFOAM v2606**. Dictionary/setup errors found and fixed:

1. `0/omega` used `atmBoundaryLayerInletEpsilon` — that BC returns epsilon,
   dimensions `[0 2 -3]`, on a field whose dimensions are `[0 0 -1]`.
2. The proposed 500×500×**80** m domain has its top **below** the tallest
   building. Domain re-sized from the real H to 1512×1208×304 m.
3. Ground roughness was to equal the inlet's 0.5 m urban z₀ — that
   double-counts roughness, since 0.5 m is an *equivalent* roughness standing
   in for buildings we resolve explicitly, and it is numerically incompatible
   with sub-metre near-ground cells. Ground now uses z₀ = 0.03 m.
4. Meshing made serial: the `decomposePar → parallel snappy → solve` pipeline
   leaves each processor's `0/` stale after `snappy -overwrite`, and re-seeding
   from `0.orig` then fails for want of processor-boundary patch entries.

Directions are realised by **rotating the site geometry** so the wind always
blows along +x, which gives every direction the same mesh and BCs — important
because the paper's claims are comparisons *between* directions.

---

### Stage 5 — P4, real Monte Carlo ✅ COMPLETE

**97,200 trajectories in 59 min** (27/s on 20 workers), through the time-averaged
CFD fields: 12 launch points × 3 headings × 3 throw types × 6 bearings × 150
release realisations. 10.2% of throws were scored as blocked by a building
rather than credited with flying through it.

### Stage 6 — P5, throwability ✅ COMPLETE

| θ_wind | Q̄ | Q_min | Q_max | wind | turb | delivery |
|---|---|---|---|---|---|---|
| 90° (easterly, over-lake) | 0.0443 | 0.0303 | 0.0563 | 0.727 | 0.744 | 8.3% |
| 210° (SSW, prevailing) | **0.0483** | 0.0261 | 0.0643 | 0.715 | 0.634 | 10.8% |
| 270° (westerly, urban) | 0.0356 | 0.0131 | 0.0471 | 0.723 | 0.563 | 9.2% |

**Best of the three computed = 210° (SSW). Worst = 270° (westerly).**

Mean disc-altitude (2 m) wind speed over the playing surface: **1.09, 1.13 and
1.01 m/s** for 90°/210°/270°, against a 5 m/s reference at 10 m — **roughly one
fifth, in every direction**.

### Hypothesis verdicts (Stage 6, item 2)

| # | Draft prediction | Verdict |
|---|---|---|
| H1a | Easterly Q gradient across the field | **SUPPORTED** (Q 0.0475 west vs 0.0410 east) — *but by the wrong mechanism, see §4* |
| H1b | Hammer systematically degraded | **NOT TESTABLE** — scoring artefact, see below |
| H2a | Westerly Q patchier than others | **CONTRADICTED** (CV 0.189 / 0.252 / 0.252 for 90/210/270 — easterly is smoothest) |
| H2b | Forehand benefits most under westerlies | **CONTRADICTED** (forehand advantage largest under the *easterly*, +0.032 vs +0.021) |
| H3 | Best conditions at a light SSW (~210°) | **SUPPORTED** (of the three computed) |

**A scoring artefact caught before it reached the paper.** The hammer's delivery
rate came out as exactly 0.0000 at every heading, which reads like a dramatic
wind effect. It is not: the hammer's still-air range in this model is **12.1 m**
and the scoring target is **25 m** downrange, so its success rate is
structurally zero whatever the wind does. The number measured the scoring rule,
not the physics. Hammer is now excluded from the Q success factor (it was
depressing every heading by about a third) and the verdict is NOT TESTABLE.

### Stage 7–8 — manuscript, review, PDF ✅ COMPLETE

- `Project1_Frisbee_LoyolaWind_PAPER.tex` — working copy (original draft untouched)
- `Project1_Frisbee_LoyolaWind_PAPER.REVIEWED.tex` → **22-page PDF, 0 errors, 0 undefined refs, 4/4 figures**
- `Project1_Frisbee_LoyolaWind_REVIEW_DOSSIER_2026-07-26.md` — four-tier review
- **Verdict: POST-READY, 5 author decisions outstanding.** Venue: arXiv
  `physics.flu-dyn`, then *Sports Engineering*.

### CFD validation results

| Test | Result | Verdict |
|---|---|---|
| ABL horizontal homogeneity (Hargreaves & Wright) | 0.7% worst streamwise drift over 1400 m; 1.3% vs analytic Richards–Hoxey | **PASS** |
| Limit-cycle amplitude at disc altitude | 13.7% (fine) / 7.4% (coarse) peak-to-peak | characterised; fields are cycle-averaged |
| Grid independence at disc altitude | 57.1% of local speed / 10.3% of u_ref | **FAIL** (10% criterion) |
| Tominaga AIJ-B cube benchmark | not run | — |

**The grid-independence failure is reported as a failure.** The criterion was
not relaxed after seeing the result. What survives it: the *magnitude* of the
sheltering (both meshes agree, every direction). What does not: the spatial
detail of Q within the field — so **no ranking of individual launch positions is
claimed**.

## 3. What was scoped down (stated in the paper, §"Scope and limitations")

| Planned | Delivered | Why |
|---|---|---|
| 13 directions, 15° increments | **3 directions: 90°, 210°, 270°** | compute |
| 5–8 M cells | **~4.5 × 10⁵ cells** | 15 GB RAM / 11 physical cores |
| 15H downstream (COST 732) | **10H** (5H upstream/lateral/top met) | compute |
| Tominaga AIJ-B cube benchmark | **not run** | time |
| Anemometer field validation | **not run** | no field data exists |
| Tree-canopy porous zone | **not included** | OSM maps 453 individual trees but no contiguous canopy east of the field; hand-tracing one was rejected as unvalidatable |
| DEM ground triangulation | **flat ground** | 3DEP shows ±1 m relief across the field |
| 10⁴ throws per direction | **150 per cell, 97 200 total** | compute |

**Machine reality:** 11 physical cores (22 threads), 15 GB RAM — against the
plan's assumed 32 cores / 128 GB.

Any statement about a "best" or "worst" heading in the paper means **best or
worst of the three computed**, and is written that way throughout.

---

## 4. The finding that contradicts the draft

The draft's §3.3 asserts that easterly inflow crosses ~80 km of unbroken Lake
Michigan fetch and that the resulting internal boundary layer is *"the dominant
local physical feature."*

**It is not.** OSM places **Campus Towers**, a 14-storey (44.8 m) apartment
block, at x ∈ [55.8, 79.3] m, y ∈ [−15.9, 21.5] m — its footprint begins **9 m
east of the field's eastern edge** and spans the full width of the field. Any
easterly flow reaching the playing surface has just passed over a 45 m bluff
body. The disc-altitude flow is a building wake, not a marine IBL.

Three further site-description corrections:

- The assumed origin was ~75 m south of the true field centroid.
- **Mertz Hall is 19 storeys (~61 m), not ~22 m** — corroborated by Loyola's
  own residence-life page. This invalidated the draft's COST 732 domain sizing.
- The field's long axis runs **E–W (78 × 25 m), not N–S**. It is also smaller
  than a regulation Ultimate pitch, so the OSM outline is likely a partial
  trace — **flagged for the author**.
- Cudahy Library and the Information Commons are **ESE** of the field, not west.

---

## 5. Provenance of every number

No number in the manuscript is transcribed by hand.
`P5_Synthesis/make_results_tex.py` generates `results_generated.tex` — LaTeX
macros and tables computed directly from:

- `P3_CFD/flow_fields.npz` — sampled disc-altitude fields (time-averaged)
- `P3_CFD/cfd_validation.json` — ABL, grid-independence, convergence verdicts
- `P4_Coupling/landing_stats.npz` — Monte Carlo landing statistics
- `P5_Synthesis/Q.npz` — throwability metric
- `P5_Synthesis/hypothesis_check.json` — hypothesis verdicts

Re-running `SeanEarl_DiscFlight/runDownstream.sh` regenerates all of them.

---

## 6. Reproducing

```bash
cd SeanEarl_DiscFlight
# P1 + P4 test suites
(cd P1_AnalyticalModel && python3 -m pytest test_disc_6dof.py -q)
(cd P4_Coupling      && python3 -m pytest test_coupling.py -q)

# geometry
(cd P2_Geometry && python3 fetch_buildings.py --half-size 350 --out site_geometry.json \
   && python3 buildings_to_stl.py site_geometry.json --out buildings_r300.stl \
        --radius 300 --single-solid buildings --bury 5 \
   && cp buildings_r300.stl seanEarl.stl && python3 verify_site.py)

# CFD (needs OpenFOAM v2606)
cd P3_CFD
python3 generateCases.py --headings 90 210 270 --resolution production --nproc 11
bash runSweep.sh 11 cases/dir_090 cases/dir_210 cases/dir_270
bash averagePass.sh 11 800 cases/dir_090 cases/dir_210 cases/dir_270

# everything downstream
cd .. && bash runDownstream.sh 20 150
```

---

## 7. Remaining author decisions

| # | Decision | Recommendation |
|---|---|---|
| 1 | Author, affiliation, co-authors, funding | Only you can supply — 6 `% FLAG FOR AUTHOR` markers |
| 2 | Field outline (OSM 78×25 m vs regulation 100×37 m) | Confirm on site; direction-level conclusions survive either way |
| 3 | Grid independence not achieved | Post with the limitation stated |
| 4 | Hammer throw scoring | Re-score at ~12 m, its own characteristic range |
| 5 | Venue | arXiv `physics.flu-dyn` now; *Sports Engineering* for the journal version |

## 8. Second-pass compute worth doing, in priority order

1. **A third, finer mesh** (~4 h/direction). Grid independence at disc altitude
   is the one failed validation and the only thing gating per-location claims.
2. **LES of the easterly case** (~2 weeks). Steady RANS cannot give peak gusts
   or spatially correlated gust structure, which is what a thrower actually
   experiences. Highest scientific value.
3. **The full 13-direction sweep** (~13 × 25 min meshing + solve ≈ 1 day here).
   Cheap now that the harness works; would turn "best of three" into a real
   direction-dependence curve.
4. **A uniform-freestream control** through the same Monte Carlo (~1 h). Isolates
   what the site-specific coupling actually buys. Cheapest item on this list and
   a referee will ask for it.
5. **A second reference wind speed** (~2 h). The Q maps rescale approximately
   with `u_ref`, but the trajectory model is not linear in wind speed, so this
   is an assumption currently untested.
6. **Building-height sensitivity** (~1 h). 44% of heights are type defaults;
   Mertz and Campus Towers are tagged, but a ±20% perturbation study would
   bound the effect.

## 9. Honest bottom line

The disc model is the strong part: quantitatively validated, five real errors
found and fixed against the primary source, 26 tests. The coupling machinery is
sound and now well-tested (18 tests). The CFD is real, validated for ABL
homogeneity, **and not grid-converged at the altitude it reports** — which is
stated in the paper, bounds what is claimed, and is the first thing a second
pass should fix.

The scientific result — that a 45 m building nine metres off the sideline, not
the lakefront boundary layer, controls disc-altitude flow — is robust to that
limitation, because it is a statement about magnitude that both meshes agree on.
