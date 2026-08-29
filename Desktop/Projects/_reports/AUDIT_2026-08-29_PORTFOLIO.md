# Portfolio audit, 2026-08-29

A sweep of `Desktop\Projects`: recompute the headline numbers of the near-finished
projects from their own raw outputs, run every test suite, fix what is wrong, and clear
the build litter. Everything below was recomputed in this session, not read off a status
card.

---

## 1. Test suites: 509 passing across 13 suites, one suite rescued

| project | result |
|---|---|
| `dst_fit` | 205 passed |
| `frisbee/P4_Coupling` | 27 passed (2 min) |
| `frisbee/P1_AnalyticalModel` | 26 passed |
| `frisbee/P5_Synthesis` | 13 passed |
| `impact_stress` | 86 passed |
| `uroflow_splashback` | 44 passed |
| `bloodstain_bpa` | 23 passed |
| `silver_to_color` | 23 passed (was 10 passed, 13 errors) |
| `beach_sand_dilatancy` | 22 passed |
| `datacenter_cooling` | 13 passed (was 10) |
| `bidet/model` | 13 passed |
| `catfish school` | 8 passed |
| `tier0_wyart_cates_closure` | 6 passed |

`silver_to_color/tools/test_pipeline.py` was written for its own `main()` runner, which
passes each test a scratch directory as a positional argument named `tmp`. Under pytest
that reads as a fixture request, no such fixture exists, and 13 of the 23 tests errored at
*setup*. They had never run under pytest at all, in either state. Added
`tools/conftest.py` supplying a `tmp` fixture that aliases `tmp_path`. Both runners now
agree; the standalone runner still reports 17/17.

`dst_fit`'s twelve gates all reproduce the exit codes recorded in `NEXT.md`
(`n1 n2 n3 n7 g5 g10 g12` → 0; `n4 n5 n6 g4 g11` → 1, each a documented
fail-by-design). Nothing rotted. `README.md` still advertised 149 tests; corrected to 205.

---

## 2. Bidet: verified end to end, two small corrections

**The mesh-truncation rewrite is correct.** Every number the uncommitted edit added was
re-derived from `case_Re10000_H6/tau_w_r.csv` without reusing `fit_closure.py`:

| claim | paper | recomputed |
|---|---|---|
| outermost cell centre | 24.77 mm | 24.7655 mm |
| τ_w there | 1.96 Pa | 1.9570 Pa |
| cleaning radius, τ_c = 3, 5, 10 Pa | 21, 16, 12 mm | 20.63, 16.18, 11.63 |
| ratio to r_max = 15 mm | 1.4, 1.1, 0.8 | 1.38, 1.08, 0.78 |
| full-band fit | n_t = 1.90, C_P = 1.35, r² = 0.997 | 1.902, 1.352, 0.99678 |
| τ_c = 1 Pa root, full closure | 39 mm | 39.14 mm |
| outer-half refit, r/d ∈ [10, 24.8] | n_t = 2.10 | 2.103 |
| τ_c = 1 Pa root, outer half | 35 mm | 34.75 mm |
| over-provisioning | 2.3-2.6× | 2.32-2.61× |

Two things fixed:

- **`r² = 0.9998` → `0.9997`** for the outer-half refit. The computed value is 0.99970.
- **`ΔP ≈ 40 kPa` → `42 kPa`.** The abstract's 8.6 Pa at ΔP = 0.07 MPa *is* right, but only
  from the stated V₀ = 8.9 m/s, which needs ΔP = 42.1 kPa at the model's C_v = 0.97. A
  reader who took the rounded 40 kPa literally reproduced 8.9 Pa and would have flagged the
  paper. Now self-checkable.

Root-finding was audited too: the outer branch has one non-monotone step (a 0.08 Pa wiggle
at r = 23.9 mm, at the outflow boundary), so `np.interp` on a reversed axis could in
principle return a wrong root. It does not here, each of τ_c = 3, 5, 10 Pa has exactly one
crossing, all far from 23.9 mm.

Cross-checks that passed with nothing to fix: oblique-jet penalty 1.95× ("~95%"),
oscillation saving 0.518 ("~52%"), τ_w(r_max) = 5.866 Pa ("≈5.9"), C_P/n_t/r² in
`bidet_results.json` identical to the CFD sweep.

**Builds clean: 12 pp, 0 errors, 0 undefined references, no overfull box over 15 pt.**

---

## 3. Datacenter cooling: a conservation-law bug in the headline number

`provisioning` *p* is defined as supplied ÷ demanded airflow, and the row sum of `B` as the
recirculated fraction of inlet air. If the CRAHs deliver only *p* of what the servers
ingest, then at least (1 − *p*) of every inlet stream has to be recirculated exhaust ,
there is no third air source in a sealed hall. So `rowsum(B) ≥ 1 − p`.

The model charged a mild `φ(η)·(1 + (1 − p))` bump instead:

| p | rowsum(B) | mass-balance floor | |
|---|---|---|---|
| 1.0 | 0.163 | 0.000 | ok |
| 0.9 | 0.179 | 0.100 | ok |
| 0.8 | 0.196 | 0.200 | **violated** |
| 0.7 | 0.212 | 0.300 | **violated** |
| 0.6 | 0.228 | 0.400 | **violated** |

Fan power goes as *p*³, so under-blowing was nearly free and **the optimum sat exactly on
the `p_min` = 0.6 bound**, the deepest point of the violated region. The headline was
being set by the part of the model that broke conservation of mass.

**Fixed** (`φ = max(φ, 1 − p)` in `influence_matrix`, two new tests). Consequences:

- The optimum becomes *interior*: *p* = 0.684, T_s = 21.0 °C. Verified against a
  brute-force 801 × 401 grid (agrees to 2 × 10⁻⁵) and (I − B) stays invertible over
  η ∈ [0,1] × p ∈ [0.05, 1.4], max spectral radius 0.95.
- **Headline PUE 1.215 → 1.236; cooling-energy cut 33% → 26%.**
- The containment sweep under the free optimum goes **flat**: once the floor binds it
  exceeds the geometric term at *every* η, so η leaves `B` entirely. The containment ROI is
  now reported at fixed *p* = 1, **PUE 1.294 → 1.281, ΔPUE 0.0126**, which supersedes the
  0.022 in the README.

**Two further gaps found and documented (`DESIGN.md` §8), not fixed**, they need the 3D
run or a modelling decision:

1. **Over-provisioning carries no thermal benefit.** `under = max(0, 1−p)` makes `B`
   identical at *p* = 1 and *p* = 1.4, so extra airflow costs *p*³ and buys nothing. The
   airflow axis is one-sided, which is why the optimizer can only ever want less of it.
2. **Server-fan power is not charged.** Cutting CRAH airflow raises inlet temperatures,
   server fans ramp, and that power lands *inside* the IT denominator, PUE improves while
   hall energy worsens. This is the textbook way a cooling model games PUE.

Gap (2) was then **closed the same day** (see §9). Gap (1) stands, so quote the fixed-*p*
containment ROI, not the free optimum. That is what `optimize_cooling.py` prints and what
`fig2` draws.

`fig2_containment_sweep.png` was rebuilt to your standing figure rules: two panels, no
titles, short legend, PUE anchored at its physical floor of 1.0 rather than the data range,
temperature anchored at 0. The explanation moved to a sibling `fig2_containment_sweep.md`.
The flat free-airflow curve is drawn rather than hidden, because its flatness is the caveat.

---

## 4. Slime mold: Table 1's MD column was not comparing like with like

`network_metrics()` averages mean inter-terminal distance over the terminals a network
*reaches*. That is not the same set for every row:

- MST, Steiner, Physarum reach all **38** terminals.
- Hub-and-spoke reaches **28**, because `baseline_hub_spoke()` routes only to the 22
  `DEMAND_NODES`.

The 10 it misses are exactly the line ends, O'Hare, Midway, 95th/Dan Ryan, Howard, Linden,
Dempster-Skokie, Forest Park, Harlem/Lake, Cermak-54th, Clark/Lake. Mean inter-terminal
distance is precisely the quantity that shrinks when the most peripheral stops are dropped,
so the hub-and-spoke row was flattered and the column read as a fair comparison.

Rescored on the common 28 terminals (headline network rebuilt from `_headline_D.npy`, which
reproduces `_metrics.json` exactly):

| network | TL | MD tabulated | n | MD common | n |
|---|---|---|---|---|---|
| MST | 147.03 | 26.95 | 38 | **22.28** | 28 |
| Steiner | 147.03 | 26.95 | 38 | **22.28** | 28 |
| Hub-and-spoke | 132.90 | 19.23 | 28 | 19.23 | 28 |
| Physarum γ=1.3 | 256.16 | 18.81 | 38 | **15.47** | 28 |

The correction strengthens the paper. The Physarum advantage over hub-and-spoke goes
from a marginal **0.42 km to 3.76 km**. The MST-to-hub gap, meanwhile, shrinks from a
tabulated 7.71 km to a real 3.04 km.

A second, plain error: §"Comparison with classical baselines" opened *"The MST and Steiner
trees are cheapest (∼145 km)"*, contradicted by the paper's own Table 1, where
hub-and-spoke is 132.9 km TL and 98.9 TC, cheaper on both. Root cause is the same terminal
mismatch.

Fixed: the paragraph rewritten, a `‡` footnote added to Table 1's MD header, and
`code/common_terminal_md.py` written so the numbers are reproducible rather than quoted.
That script asserts the rebuild matches `_metrics.json` before it reports anything.

A LaTeX bug came with it and is fixed: `\\` inside the `\multicolumn{6}{p{...}}` footnote
cell ends the *tabular row*, not the line. Uses `\par\vspace{2pt}`.

**Builds clean: 10 pp, 0 errors, 0 undefined references.** The two overfull boxes (26.8 pt,
24.5 pt) are pre-existing, identical values in the older `.EDITED.log`.

---

## 5. Frisbee: the uncommitted work checks out

The 21 modified files carry two real defect fixes, both correctly reasoned:

- `landing_stats.py` was averaging `landing_mean`, an **absolute world-frame position** ,
  over a 65 m launch grid and six bearings, and printing the result as "mean drift". The
  new code uses `landing_bias` and refuses to run on artifacts that predate it.
- `throwability.py` now records which throws fed the success factor. The shipped `Q.npz`
  had been built with `--exclude-throws hammer` while `runDownstream.sh` passed no such
  flag, so **the committed script did not reproduce the committed artifact**, a factor of
  ~1.5 in Q.

Verified independently: `q_significance.json` matches `CORRECTION_APPLIED_2026-08-09.md` to
the digit (90 vs 210 p = 0.206, CI [−0.0095, +0.0018], favoured 4/12; 90 vs 270
p_holm = 0.031, 11/12; 210 vs 270 p_holm = 0.026, 10/12; spatial-to-MC ratios 8.8-13.6).
The 24-28% sheltering fraction recomputes exactly from 1.21/1.42/1.40 against 5 m/s.

**Builds clean: 26 pp, 0 errors, 0 undefined references, 0 `\MISSINGNUM` renders.**

---

## 5b. Frisbee: GATE G1 was run for the first time, and it fails

WP1's instrument had been written and never executed; there was no
`gust_model_comparison.json` on disk. Full write-up in
`frisbee/GATE_G1_RESULT_2026-08-29.md`.

`sigma_turb(temporal)/sigma_turb(sweeping)` = **1.735** against a pre-registered acceptance
band of 2.0-2.6×. A cluster bootstrap over the 108 cells (20,000 resamples, same cells drawn
for both arms so the paired design survives) gives a 95% CI of **[1.638, 1.838]** and
**P(ratio ≥ 2.0) = 0.0000**. The band is excluded, not missed.

Both documented failure branches are excluded too: it is not the "< 20%, the defect was
harmless" case (P = 0.0000), and it is not the "≫ 3×, not a linear filter" case. Inverting
the brief's own `√N_eddy` argument: **N_eddy = 3.01, 95% CI [2.68, 3.38]**, against the
predicted N = 4-6. The disc integrates about three independent eddies, not four to six.

A robustness re-run with `--isotropic` (isotropic σ in both arms) gives 1.868, 95% CI
[1.789, 1.955], P(ratio ≥ 2.0) = 0.0018, N_eddy = 3.49. **The verdict does not depend on
the σ choice**, and anisotropy moves the ratio the wrong way (1.87 → 1.74) to explain the
gap.

So the correction is real and large, the old temporal model overstated turbulence-induced
landing dispersion by 74%, but the theory that predicted its *size* is wrong. Per the
brief, WP2 should not start on the 4-6 eddy count.

**The run exposed a separate blocker.** `landing_stats.npz` (2026-07-26 03:19) predates
`monte_carlo.py` (14:19 the same day) by eleven hours and has **none** of the schema the
current `P4_Coupling/README.md` documents, no `landing_bias`, no `target_xy`, no
`schema_frames`, no `prov_*` block. In particular **it does not record which gust model
built it**, so it cannot be established from the artifact whether the manuscript's Q rests
on `temporal` or `sweeping`. Since the file predates the rewrite, `temporal` is the likely
answer, and that is the model G1 just measured as a 74% overstatement. `Q.npz` derives from
this file, and every manuscript macro derives from `Q`.

The fix is `bash SeanEarl_DiscFlight/runDownstream.sh 20 150`, 97,200 trajectories, about
1.5 h at the 19.5 traj/s measured here. Not run in this session because it moves every
downstream number and the paired-significance verdict (`worst_supported: 270`) could move
with them.

---

## 6. Cleanup

- **115 LaTeX build artifacts (3.07 MB) removed**, plus 6 orphans in the Projects root.
  Selected by matching each artifact against a `.tex` of the same basename *in the same
  directory*, so run logs could not be caught by an extension glob ,
  `catfish school/logs`, `dst_fit/logs`, `cfd_openfoam/sweep_console.log` and
  `slime mold/code/run.log` were all untouched. No tracked file was deleted (git tracked
  zero build artifacts).
- **61 `__pycache__` / `.pytest_cache` / `.ruff_cache` directories (3.0 MB) removed.**
- `paper1_referee_review.pdf` and `paper1_research_approach.pdf` kept in the Projects root:
  they are **newer** than the copies in `paper1/` (2026-06-30 vs 2026-04-29 / 2026-05-13),
  so the root holds the current builds while the sources live in `paper1/`. Worth
  reconciling by hand.

### The repository had no active ignore rules at all

The repo is rooted at `C:\Users\compl`, your home directory, and its ignore file is named
**`gitignore`, without the leading dot**, so git has never read it. It even has an
`# API Keys / .env` section that was doing nothing. `git status` was listing `.ssh/`,
`.claude.json`, `AppData/`, `Calibre Library/` and the rest of `$HOME` as untracked, and a
bare `git add -A` would have staged them. (`git add -An` did not even finish inside two
minutes.)

Added a `.gitignore` that ignores `$HOME` by default and re-includes only the paths the
repository actually tracks. Verified: `.ssh/id_rsa`, `.claude.json`, `AppData/**`,
`Documents/**`, `Downloads/**` are ignored; every file under `Desktop/Projects` stays
trackable; **no already-tracked file became ignored**; `git status` went from timing out to
0.11 s and 114 entries. The old `gitignore` was left in place, and its rules folded in.

---

---

## 8. Online data that would make the story better

Ordered by how much each buys per hour spent. Everything here is public and free.

### Bidet: the cheapest win needs no new data at all

The τ_c = 1 Pa contour is the one the domain truncates (τ_w has only fallen to 1.957 Pa at
R = 25 mm). It is **also exactly where the paper's own best-supported literature value
sits**: `RESEARCH_NOTES_citations_tauc.md` settles on Nejadnik et al. (2008), a genuine
wall-detachment measurement, at 1.0 ± 0.2 Pa for *S. aureus* and 1.1 ± 0.2 Pa for
*P. aeruginosa*. So the most defensible cleaning threshold in the paper is the only one the
computation cannot resolve. **Re-running `case_Re10000_H6` at R ≈ 50 mm is free compute and
converts the paper's weakest sentence, a bound, r > 25 mm, into a computed number.** Do
this before submitting; it is the one open item the paper flags against itself.

Data that would still help, in order:

1. **Manufacturer spec sheets** (Toto Washlet, Bio Bidet, Brondell, Tushy): flow rate in
   L/min, nozzle diameter, operating pressure. The paper currently calibrates against a
   single "commercial datum". A dozen products turns *"mainstream designs sit within a
   factor of ~2 of the cleaning constraint"* from an n = 1 claim into a distribution, which
   is a much harder statement to referee away.
2. **Residential supply-pressure distributions** from plumbing-code and utility sources, to
   justify the ΔP = 0.07 MPa operating point rather than asserting it.
3. **A non-bacterial soil-removal shear threshold.** The current bracket is bacterial
   adhesion at the low end and *Listeria* biofilm (24-144 Pa) at the high end. Cleaning-in-
   place literature for soft organic soils would put a defensible number in the middle of a
   two-decade-wide band.

### Datacenter: both documented gaps have published closures

The two `DESIGN.md` §8 gaps are not open research questions; the numbers exist.

- **Server power vs inlet temperature (gap 2, the PUE-gaming one).** Published measurements
  give a room-wide sensitivity of **0.35-0.5 %/°C** across the 23-30 °C ASHRAE band, a
  **7-20% total server-power increase** going 15 → 35 °C, and required air volume rising to
  **1.2-2.5×** over the same span, with server fans typically beginning to ramp at about
  **25 °C**. That last number matters: the model's optimum places the hottest rack at
  exactly 27 °C, *inside* the ramping band, so the omission is load-bearing rather than
  academic. Sources: Schneider Electric's *Energy Impact of Increased Server Inlet
  Temperature*, the DOE/LBNL wider-operating-range white paper, and the measurement-based
  server power model in *Energy* (2022).
- **A warning, because it is an easy wrong turn:** ASHRAE's **X-factor is not this**. It is
  a dimensionless *hardware-failure-rate* multiplier versus inlet temperature, baselined at
  20 °C. Useful for a reliability caveat, useless for the fan-power term.
- **Real halls, for the outreach step.** ENERGY STAR / LBNL data-centre benchmarking, the
  Uptime Institute annual PUE survey, and Open Compute Project reference designs give
  published layouts and PUE distributions, enough to say where the demonstrator's 1.24-1.29
  band sits against the field, which is the first question any company will ask.

### Slime mold: replace expert scores, and get a validation test for free

The paper's limitation (ii) already concedes the demand sub-scores are "normalized expert
estimates". They are: `data_chicago.py` assigns `pop`/`activity`/`underserved` by hand on a
0-1 scale. Two levels of fix:

1. **Better inputs.** ACS block-group population → ρ_i; LEHD **LODES** workplace-area
   characteristics (`lehd.ces.census.gov`) → a_i. The Chicago Data Portal also has
   `CTA-L-Rail-Stations-Shapefile` (`vmyy-m9qj`) for real station geometry in place of the
   schematic anchors.
2. **The better idea, and the biggest single upgrade available to this paper.** Ridership
   data cannot calibrate demand at *un-served* nodes, but it enables a **validation the
   paper currently has none of**. `CTA - Ridership - 'L' Station Entries - Daily Totals`
   (`5neh-572f`, 1.31 M rows, daily entries per station back to 2001, still updated) lets
   you run the same Physarum solver on the **existing** L topology and test whether the
   converged tube conductivities correlate with measured station entries. That moves the
   model from "plausible and seed-stable" to "reproduces a network we can check", which is
   the difference between a preprint and a paper, and it costs no new physics, only a
   download and a rank correlation.

### Frisbee: wind climatology turns a ranking into an expected cost

The study weights 90°, 210° and 270° equally, and after the August correction its surviving
claim is "270° is worst". Real wind statistics change what that claim is worth. NOAA's
**Integrated Surface Database** (hourly wind direction and speed, 20,000+ stations, also
mirrored on AWS Open Data) plus NOAA **GLERL's Chicago lakefront** met station, the
relevant one for a field this close to the water, would let you weight the three headings
by measured frequency, restricted to the hours and season the sport is actually played.
"270° is worst" then becomes "the field is degraded for X% of playing hours", which is a
statement a reader can use. One data pull, no new simulation.

### dst_fit: more reversal transients, not more flow curves

The strain clock is the result that a flow curve provably cannot contain, and it currently
rests on Lin et al.'s reversal transients: two material systems. Flow-curve deposits will
not extend it, reversal or step-strain transients will. Worth sweeping Zenodo, journal
supplementary material, and further Edinburgh DataShare deposits specifically for
**shear-reversal viscosity transients**. Measuring γ₀ on a third and fourth system is what
turns "0.5-1.6 on these two deposits" into a claim about dense suspensions.

## 7. What still needs you

1. **Nothing is committed.** 22 tracked files are modified, including the whole verified
   bidet and frisbee body of work. Say the word and it goes in as separate commits.
2. **`datacenter_cooling` and `slime mold` are not in git at all**, they are untracked
   directories. The two projects where this session found the most significant defects are
   the two with no version history.
3. **Bidet: venue.** It is the most submittable thing here and nothing technical is
   blocking it.
4. **Slime mold: two placeholders**, the department in the author block (the tree gives
   three different ones for you) and `[INSERT REPOSITORY URL/DOI]`. Plus four items the
   reviser deliberately deferred (§(a)-(d) in the file header).
5. **Datacenter: the airflow modelling decision** in §3 above.

---

## 9. Implementation: what was built after the audit

### Datacenter: server-fan power is charged, and PUE turns out to be conditionally safe

Gap (2) above is closed. `dc_thermal` now adds a server-fan ramp of 0.4% of rack IT power
per K above a 25 °C inlet knee, applied per rack against its own inlet temperature, and
charges the chiller for it. `pue()` puts the ramped IT in both numerator and denominator;
`total_kw()` sits beside it; `optimize(..., objective="total")` ranks by energy.

Constants are published measurements: 0.35-0.5 %/K room-wide sensitivity across the
23-30 °C ASHRAE band, 7-20% total server-power increase over 15 to 35 °C, air volume
1.2-2.5x over the same span, ramp knee near 25 °C. **ASHRAE's x-factor is not this** and is
the wrong thing to reach for: it is a hardware failure-rate multiplier versus inlet
temperature, with no power content.

The finding is more useful than the term. Extra hall power from optimizing PUE rather than
energy:

| ramp %/K | T_allow 27 | 32 | 35 | 40 |
|---|---|---|---|---|
| 0.35 | 0.00 kW | 0.00 | 0.07 | 0.45 |
| **0.40** | **0.00 kW** | 0.15 | 0.54 | 1.34 |
| 0.50 | 0.25 kW | 1.48 | 2.56 | 4.19 |
| 0.75 | 1.44 kW | 5.57 | 8.34 | 12.03 |
| 1.00 | 2.62 kW | 9.64 | 14.12 | 19.87 |

At the mid-band coefficient and ASHRAE's 27 °C **recommended** limit the two objectives
pick the same design: **the inlet constraint, not the objective, is what prevents the
gaming.** The gap opens on the **allowable** envelopes and at the steep end of the ramp,
where the lowest-PUE design in the whole study, **PUE 1.197** at 1.0 %/K and T_allow 35,
buys that number by moving **20 kW, 10% of nominal IT**, into the denominator and burns
14 kW more than the energy optimum.

So PUE is a safe objective for this hall only while the 27 °C recommended limit holds.
That also names the first constant to ask a real company for: their ramp coefficient.

Headline moves 1.236 → 1.235, cooling cut 26.4 → 25.9%. Tests 13 → 16, including one that
would catch the term being neutered.

### Bidet: the 1 Pa contour is computed, not extrapolated

`gen_case.py` had `Rdom` hardcoded at 25 mm. It is now a parameter, and the radial cell
count scales with the span so widening cannot silently coarsen the wall jet. The grading is
re-solved to hold the **first cell width** fixed at the 7.0512 µm the validated mesh
produces, rather than holding the growth rate, which would have compounded over the extra
cells and driven the first cell to 0.2 µm. At the default radius the generator reproduces
the 2026-07-25 `blockMeshDict` **byte for byte**, which is the regression check that makes
the widened run comparable.

`compare_domain_extent.py` does the two-part test in that order: the wide case must agree
with the narrow one on the overlap, and only then is the root read off the computed field.

**Result. τ_c = 1 Pa is crossed at r = 34.79 mm**, inside the domain, with τ_w reaching
0.455 Pa at the outer boundary. Validation passed: median 0.7% and max 1.8% agreement in
τ_w for r ≤ 22.5 mm, diverging only in the outer 10% of the narrow domain (5.7% at its last
cell), the signature of an outflow boundary reaching back into the profile it was meant to
terminate. Converged in 2689 SIMPLE iterations, y⁺ 0.047-1.39, max aspect ratio 106.6
against the narrow case's 105.9.

The two extrapolations can now be scored, and the answer is transferable:

| method | r at τ_c = 1 Pa | error |
|---|---|---|
| **computed, R = 50 mm** | **34.79 mm** |, |
| outer-half refit, r/d ∈ [10, 24.8] | 34.75 mm | **−0.1%** |
| full-band closure pushed past r/d = 20 | 39.14 mm | **+12.5%** |

The far field is steeper than the whole-band fit, so extrapolating a single power law
through a profile whose local slope is still steepening overestimates the reach; a local
refit does not.

Manuscript updated: all four cleaning radii now read off the computed field, over-provisioning
at 1 Pa is a computed 2.3× rather than an extrapolated 2.3-2.6×, τ_w(r_max) 5.9 → 6.0 Pa and
the 0.07 MPa value 8.6 → 8.7 Pa, and the "Radial domain truncation" limitation becomes
"Radial domain extent". The figure is redrawn to your figure rules with a sibling
`fig_cfd_tau_w.md`. **12 pp, 0 errors, 0 undefined references.**

### Frisbee: the chain is regenerated, and the verdict survives it

`rerunMonteCarloOnward.sh` re-ran `runDownstream.sh` steps 3-8 and deliberately skipped
steps 1-2, so `flow_fields.npz` stayed byte-identical to the field both the manuscript and
GATE G1 were computed on. The Monte Carlo was then the only variable that moved. 97,200
trajectories, 648 cells × 150 throws, ~110 min.

The provenance gap is closed: `landing_stats.npz` has 42 keys instead of 19, and
**`gust_model = sweeping`** is now recorded in the artifact rather than inferred from file
timestamps.

**The conclusions did not move, and that is the result.** 17 of 58 results macros and 19 of
35 significance macros changed value; none were added or removed, so the `\input` wiring and
the `\MISSINGNUM` guard are intact.

| comparison | before | after | verdict |
|---|---|---|---|
| 90 vs 210 | p = 0.206, 4/12 | p = 0.203, 4/12 | **not resolved**, unchanged |
| 90 vs 270 | p = 0.031, 11/12 | **p = 0.022**, 10/12 | resolved, unchanged |
| 210 vs 270 | p = 0.026, 10/12 | **p = 0.0038**, 11/12 | resolved, unchanged |

Mean Q at 270° fell 0.0356 → 0.0324, which is why both comparisons against it tightened.
That the ordering survives a full regeneration under a corrected gust model, a corrected
reference frame and a corrected throw-exclusion flag is the strongest evidence yet that
"270° is the worst of the three computed" is a property of the field, not of the pipeline.

It also caught one more typed number: the manuscript hardcoded "9 to 14 times the Monte
Carlo standard error" in a project whose whole convention is that manuscript numbers are
generated. The regenerated values are 12, 13 and 13, and the sentence now reads them from
`\QsdOverSeH*`. **26 pp, 0 errors, 0 undefined references, 0 `\MISSINGNUM`.**

GATE G1 is still failed and untouched by this. Do not start WP2 on the 4-6 eddy count.
