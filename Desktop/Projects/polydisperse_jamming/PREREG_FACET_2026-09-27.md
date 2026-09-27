# PREREG_FACET_2026-09-27 — do convex facets narrow the shear-jamming window? (N 200 replication + near-sphere ladder + contact census)

Written 2026-09-27 ~10:55 CDT, **before any of these cells ran.** Follows `PREREG_LSPILOT_2026-09-26` (pilot, N 100),
whose registered REVERSAL test ended INCOMPLETE and whose Amendment 5 left the window change UNRESOLVED
(+0.022 ± 0.018, narrower on 2/4 seeds) while resolving a frictionless-edge drop φ₀(1.0) − φ₀(0.5) = 0.026 (t 17.5).
Machine idle at launch (load 0.00, 22 threads); 20 workers.

## Cells (script `ls_pack.py`, launcher `scripts/ls_facet.sh`)
Level-set superellipsoids, eps_e = eps_n = eps, volume-equivalent to d = 1, 1000 surface nodes, grid spacing 0.05,
**N 200**, p̂ 1e-3, max 300000 iterations, `--census` on (post-run contact dump; verified 2026-09-27 on a smoke cell
that the census leaves φ and Z bit-identical: 0.6138650806557239 with and without). Tag `facet`. Queue order =
priority order:
- **A** eps {1.0, 0.5} × μ_p {0, 0.3} × seeds 2–5 (16 cells). With the existing N 200 seed-1 cells (`ls_n200`) this is
  5 seeds.
- **B** eps {0.9, 0.8} × μ_p 0 × seeds 1–3 (6); **C** the same at μ_p 0.3 (6); **D** eps 0.7 × μ_p {0, 0.3} × seeds 1–3 (6).
- **E** seed 1 of block A's four (eps, μ) re-run with `--census` (4): expected BIT-IDENTICAL φ to `ls_n200` (wiring gate W0).

Pairing: the RNG seed sets both the insertion positions and the orientations regardless of eps and μ_p, so
cells at one seed are paired across eps and friction; every difference is taken seed-paired.
Health: converged (unbalanced < 0.02 before the cap). Unconverged cells are excluded and reported, never graded;
a test missing a member is INCOMPLETE, never PASS/FAIL. Graded only when every member cell has ended.

Error scale: se_eff = max(seed-paired standard error, √2·PHI_SCALE/√n) on φ differences, PHI_SCALE = 2.1e-3 (the
project's declared run-to-run reproducibility, `mixture_model.py`); for W the φ scale is divided by φ₀.

## Registered tests
- **W0 (wiring):** block E reproduces `ls_n200` φ exactly on all 4. Failure is reported and block E's census used anyway.
- **F1 (the question — window):** D_s = W_s(1.0) − W_s(0.5), W = 1 − φ_m/φ₀, over seeds 1–5 at N 200.
  **NARROWS** if mean D > 3 se_eff and D_s > 0 on ≥ 4/5 seeds; **WIDENS** if mean D < −3 se_eff and D_s < 0 on ≥ 4/5;
  otherwise **UNRESOLVED**.
- **F2 (frictionless edge, replication):** φ₀(1.0) − φ₀(0.5) > 3 se_eff and positive on ≥ 4/5 seeds → REPLICATES.
- **F3 (branch attribution, only if F1 is not UNRESOLVED):** the frictionless branch carries it if
  |mean Δφ₀| > |mean Δφ_m| (the pilot's mechanism check).
- **F4 (near-sphere shape of φ₀(eps) — two competing readings registered):**
  (i) the pilot's reading: facets lower φ₀ monotonically from the sphere, so φ₀(0.9) ≤ φ₀(1.0);
  (ii) the generic non-spherical-grain reading (ellipsoids, dimers, and — as recalled, **[UNVERIFIED]** until the
  literature lookup running at registration returns — superballs): small deviations from the sphere RAISE the
  frictionless jammed φ, so φ₀(0.9) > φ₀(1.0).
  Δ = φ₀(0.9) − φ₀(1.0) seed-paired over seeds 1–3: **RISE** if Δ > 3 se_eff, **FALL** if Δ < −3 se_eff, else
  UNRESOLVED. Same statistic reported at 0.8 and 0.7. A RISE at 0.9 followed by a fall by 0.5 would make the pilot's
  "facets lower φ₀" a statement about sharp facets only.
- **F5 (contact census — are faces doing it?):** m = max|nᵢ| of the contact normal in a grain's body frame (1 at a face
  centre, 1/√3 at a corner). Null for a given shape: the distribution of m under UNIFORM-AREA sampling of that
  grain's surface (computed by Monte Carlo from the analytic superellipsoid, script `facet_census.py`), and for two
  independent contact partners the face–face fraction f_ff⁰ = P_area(m > 0.95)². Measured f_ff = fraction of contacts
  with m > 0.95 on both grains. Face-contact reading: **E_ff = f_ff/f_ff⁰ > 1.25 at eps 0.5, μ_p 0** (pooled seeds,
  reported with the per-seed range). The eps 1.0 cells are the sanity null (E_ff must be 1 within scatter there,
  since a sphere has no body frame to align). Reported, not graded: E_ff at μ_p 0.3; the torque-lever distribution
  |(c − x) × n| (on spheres it measures the level-set's discretisation, 0.03–0.09 on a smoke cell, so only its excess
  over eps 1.0 is read); Z₀ against the non-spherical isostatic count 12 (a pack well below 12 is hypostatic).

## What each outcome would mean for the author's question 2
F1 NARROWS with F3 on the frictionless branch → convex facets act against shear jamming, opposite to
elongation/concavity, and a faceted grain like cornstarch sits on the narrow side unless its irregularity dominates.
F1 UNRESOLVED at N 200 × 5 seeds → at this protocol the window change of convex faceting is below ~3 se_eff, i.e.
faceting moves the window's POSITION (φ₀) but not measurably its WIDTH; the rheometer's χ remains the decider.

## Amendment 1 (2026-09-27 ~11:35, BEFORE any registered cell ran) — a method bug voids the pilot; redesign

**The bug.** `ls_pack.py` used `Law2_ScGeom_FrictPhys_CundallStrack()` with Yade's default `sphericalBodies=True`, which
computes contact torques from radii: a NORMAL force never torques the grain. Right for spheres, wrong for every other
shape — Yade's own `examples/levelSet/levelSetBody.py` says `sphericalBodies=False` "is required in the general case".
Found because the verified literature disagreed with the pilot: Jiao, Stillinger & Torquato, PRE 81, 041304 (2010),
doi:10.1103/PhysRevE.81.041304 (full text read by the literature agent this session) report that frictionless MRJ
superball φ RISES monotonically as p = 1/eps moves away from 1 (cusp at the sphere; one-sided slope a₊ ≈ 0.229 per
unit p, their ref. 43; φ(p) values only in Fig. 6, **[UNVERIFIED]** as numbers) with Z plateauing near 8.15, well below
the isostatic 12; Yuan et al., Soft Matter 15 (2019), doi:10.1039/C9SM01932D (soft particles, quasi-static compression,
N 400) agree with Jiao at p = 1.25–2. The pilot had φ₀ FALLING by 0.026 at eps 0.5 (p 2) with Z₀ ≈ 6.3 — exactly the
translational isostatic count, as expected if rotations were frozen.
**Proof** (`scripts/ls_facepair.py`, two frictionless eps-0.5 grains face to face, B tilted 5°, pushed onto A):
tilt after 20000 steps = 5.00° (single-point law, sphericalBodies True), 5.00° (multi-point, True), **0.02°**
(single, False), **0.03°** (multi, False). So with the default, frictionless grains could not rotate at all; the
single-point contact itself holds a flush face contact once torques are live.
**Consequence:** every level-set record before this amendment (`ls_pilot`, `ls_ladder`, `ls_nodes`, `ls_repl`,
`ls_n200`, `ls_extra`, and the killed `method` check) is a pack of ORIENTATION-FROZEN grains, not of rounded cubes.
The pilot's "convex facets lower the frictionless edge φ₀ by 0.026 (t 17.5)" is withdrawn as a method artefact.
Clump packs (sphere members) are unaffected: a normal force on a sphere member acts through its centre.

**Redesign (the registered tests F1–F5 keep their bars; the cells change):** `SPHERICAL_BODIES = False`, recorded
per record. Launcher `scripts/ls_facet2.sh`, 20 workers, queue in this order:
- **Block 1** N 100, seeds 1–3, eps {1.0, 0.9, 0.8, 0.7, 0.5} × μ_p {0, 0.3}, `--post-iters 20000` (30 cells, seed-major).
- **Block 2** N 200, seeds 1–5, eps {1.0, 0.5} × μ_p {0, 0.3}, no post iterations (20 cells).
F1/F2/F3 are graded on Block 2 (as registered: N 200, 5 seeds) and repeated on Block 1 as a second size; F4 on
Block 1; F5 on both. W0 is dropped (no valid prior record exists to reproduce).
**New anchor A-lit (method validity, registered):** frictionless φ₀ must RISE from eps 1.0 to 0.9 and to 0.8
(seed-paired, Δ > 3 se_eff), as Jiao and Yuan report for superballs. If it does not, the level-set method is not
validated for this question and F1–F3 are reported as method-limited.
**Primary value:** φ, Z are read at the registered stop rule (`post[0]` where post iterations ran). The drift over
the 20000 post iterations is reported; if |Δφ| > PHI_SCALE on eps < 1 cells, the stop rule (force balance only,
no torque criterion) is flagged as inadequate for rotating grains.

## Amendment 2 (2026-09-27 ~11:35 CDT; one cell landed, disclosed; nothing graded) — health at the stop rule
The record-level `converged` in `ls_pack.py` is evaluated after the 20000 post iterations, where the stress servo keeps
the unbalanced force jittering at 0.03–0.06; the first landed cell (N 100, eps 0.5, μ_p 0.3, seed 2) met the registered
stop rule at iteration 8081 (unbalanced 0.0194) and reads `converged False` at the end. As Amendment 1 already defines
the primary value at the stop rule, the health rule is read there too (`post[0]`: unbalanced < 0.02 before the cap);
cells without post iterations use `converged` as before. `grade_facet.healthy()`.

## Amendment 3 (2026-09-27 ~12:40 CDT; 3 superball cells landed, none graded) — the sphere reference keeps the torque-free law
Probe (repo `ls_pack.py` + `LS_PROGRESS`, N 100, eps 1.0, μ_p 0.3, seed 1): with `sphericalBodies=False` a level-set
SPHERE's unbalanced force plateaus at 0.12–0.13 through 4500 iterations (the frozen version reached the stop rule at
3416). Cause: the contact point is the deepest surface NODE, laterally off the line of centres by up to ~half a node
spacing, so a normal force acquires a spurious lever arm (0.03–0.09 r on the census smoke cell) and spins spheres,
which have no restoring orientation. For a sphere `sphericalBodies=True` is the exact physics, so:
- **eps 1.0 cells use `--spherical-bodies 1`.** The eps-1.0 facet2 cells (rotations live) were killed at 12:18 before
  landing (rc 143) and a watcher kills any queued one; none is graded.
- **The pilot's eps-1.0, 1000-node records are VALID sphere references** (same code otherwise, same seeds, P0-passed):
  N 100 seeds 1–5 (`ls_ladder`, `ls_nodes`, `ls_repl`, `ls_extra`) and N 200 seed 1 (`ls_n200`). New: N 200 seeds 2–5 ×
  μ_p {0, 0.3}, tag `facet2ref`, census on (8 cells, `scripts/ls_facet2ref.sh`).
- The same node noise acts on superball contacts; unlike spheres they have a restoring orientation, and two μ_p 0.3
  superball cells met the stop rule. Whether FRICTIONLESS superballs reach it is being watched (four seed-1 cells
  restarted identically with progress logging, `scripts/ls_facet2_vis.sh`). An unconverged cell stays excluded.

## Amendment 4 (2026-09-27 ~12:55 CDT) — frictionless superballs do not reach the stop rule; a low-friction edge replaces them
Landed and disclosed at registration (seed-paired, N 100, stop-rule φ): μ_p 0.3: eps 0.9 s1 0.6003; eps 0.5 s1, s2
0.6214 ± 0.0029; sphere references (Am. 3) μ_p 0 0.6387 ± 0.0025, μ_p 0.3 0.5827 ± 0.0090 (seeds 1–3). No μ_p 0
superball cell has landed.
**Observation (progress logs, `facet2vis_*`, seed 1):** frictionless superballs keep the unbalanced force at
0.26 (eps 0.5), 0.32 (0.7), 0.45 (0.8), 0.60 (0.9) at 10 000 iterations, flat since ~6000, while φ creeps upward
(eps 0.5: 0.731 → 0.733). Same cause as Am. 3: node-discretisation torque noise (lever arm ~0.03–0.09 r) with only
Cundall damping to remove it; the nearer the sphere, the weaker the restoring orientation and the worse the plateau.
At the registered bar these cells would all end unconverged at the 300 000 cap (~10 h), i.e. excluded.
**Action:** the 12 queued/running frictionless superball cells without logging are stopped (N 100 seeds 2–3 and all
N 200); the four logged seed-1 cells run on and are reported DESCRIPTIVELY (φ and unbalanced vs iteration), never graded.
A-lit, F2, F4 and every F1 built on φ₀ are therefore **INCOMPLETE by method**, not failed.
**Registered replacement (before any of these cells runs):** a low-friction edge **μ_p 0.1**, the lowest friction at
which the Coulomb torque capacity μ_p r F_n exceeds the noise lever arm (~0.06 r F_n on the census smoke cell), N 100,
seeds 1–3, eps {0.9, 0.8, 0.7, 0.5} rotations live + eps 1.0 with the torque-free law (Am. 3), no post iterations,
census on (15 cells, tag `facet2mu01`).
- **F6 (frictional branch, graded on the μ_p 0.3 Block-1 cells):** Δ_m = φ_m(eps) − φ_m(1.0), seed-paired over seeds 1–3,
  se_eff as registered: RISE / FALL / UNRESOLVED at 3 se_eff with the sign on 3/3 seeds, per eps.
- **F7 (low-friction edge):** Δ_0.1 = φ(eps, μ 0.1) − φ(1.0, μ 0.1), same rule, per eps. Literature expectation
  (sign only): superballs pack denser than spheres at low friction (Jiao 2010, frictionless) and lose density faster
  with friction (Salerno et al., PRE 98, 050901(R) (2018), doi:10.1103/PhysRevE.98.050901).
- **F8 (the window proxy):** W_0.1 = 1 − φ(μ 0.3)/φ(μ 0.1), D = W_0.1(eps) − W_0.1(1.0) seed-paired: WIDENS / NARROWS /
  UNRESOLVED by F1's rule (3 se_eff, sign on all 3 seeds). W_0.1 is NOT the registered window (its low edge is not
  frictionless); it is reported as "the window between μ_p 0.1 and 0.3", and the frictionless edge stays open.
