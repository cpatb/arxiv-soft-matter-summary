# PREREG — the shape axis, 2026-09-18 (design registered; no cell run)

The author's call on the design fork (2026-09-18): **clumps now, LevelSet as the
control** — build the cheap axis, then measure the interlocking confound instead of
caveating it. This file registers the parameterisation, what it holds fixed, the two
defects a capability probe found before anything was built on them, and the gates.

**Nothing has been run.** `shape.py` and `test_shape.py` (19 anchors) are the geometry;
`scripts/probe_clump.sh` is a capability probe that writes no shard.

---

## 1. The parameterisation, written out — because a family is a PATH

The 2026-09-13 lesson: a one-parameter family is a path through a multi-parameter plane,
and that is exactly how this project manufactured its "bidisperse is non-monotone in CV"
negative — the family was 50:50 **by number**, so its fines carried 0.500 → 0.036 of the
**volume** as CV rose and it never visited the composition where the effect peaks. So the
shape axis is written out before it is run.

A grain is a rigid clump of two identical spheres of diameter `ds`, centres `s` apart,
with `u = s/ds` and **aspect ratio α = 1 + u**, α ∈ [1, 2]. The union volume is
`V = ds³·g(u)` with

    g(u) = π/3 − (π/12)(2+u)(1−u)²   (u ≤ 1),      g(u) = π/3   (u ≥ 1)

so the sub-sphere diameter that makes the grain's **volume-equivalent diameter** exactly
`d_eq` is `ds = d_eq·((π/6)/g(u))^(1/3)`.

**HELD FIXED along the α axis, by construction and not by hope:** `d_eq` is drawn from the
*same* `psd.make_psd`, so CV, `d32`, `d43`, skewness, span and every other PSD descriptor
are **identical at every α**. `mu_p`, `N`, the pressure ladder and the protocol are
unchanged. **MOVED:** α, and nothing else about the size distribution.

That makes **(CV, α) a grid, not a path**. The registered grid is
`α ∈ {1.0, 1.25, 1.5, 1.75, 2.0}` × the existing CV ladder, and **α = 1 is the first rung
deliberately**: at α = 1 the construction returns `ds = d_eq` exactly — the grain *is* the
sphere the whole project has run — so the shape campaign must reproduce the existing
sphere results at its own α = 1 before anything at α > 1 is believed.

**MOVED AS A SIDE EFFECT, and this is the confound:** a clump of overlapping spheres has a
bumpy surface. The waist depth relative to the equivalent radius is
`eps(u) = (ds/d_eq)·[1 − √(1−u²)]`, which runs **0 at α = 1 to 2^(−1/3) = 0.7937 at α = 2**,
where the neck radius is exactly zero. This is not a caveat, it is a second knob: at the
top of the ladder the grain's surface is as corrugated as a two-sphere body can be.
`eps` is recorded on every cell.

**THE CONTROL THAT TURNS IT INTO A MEASUREMENT.** A LevelSet body — smooth, matched α,
matched volume, `eps = 0` by construction — run at the same (CV, α) points. The
clump-minus-LevelSet difference *is* the interlocking term. Until that control exists,
**no result from the clump axis may be attributed to elongation rather than to
interlocking**, and the artifacts must say so in those words.

---

## 2. What the capability probe found, before anything was built on it

`scripts/probe_clump.sh`, α = 1.5, 60 clumps in a periodic cell, 50 steps. Four questions,
two of them answered badly, which is why the probe existed.

1. **Clumps work in the protocol.** 60/60 built in a periodic cell; the collider accepts
   them with `allowBiggerThanPeriod = False`; `PeriTriaxController`-style integration ran.
2. 🔴 **`utils.porosity()` DOUBLE-COUNTS the member overlap.** Measured: φ from the
   analytic union volume 0.06135923, φ from `utils.porosity()` **0.07272205**, ratio
   **1.185185**, which is exactly `2(ds/d_eq)³`. **This ratio is α-dependent and runs
   2.000 at α = 1 down to 1.000 at α = 2** — so an uncorrected clump campaign would have
   reported φ **falling by a factor of two along the aspect ladder from the bookkeeping
   alone**, in a project whose central observable is φ. `packing.measure` must compute φ
   from the analytic union volume for clump runs; `utils.porosity()` must not be called.
3. **No internal contacts.** 61 real interactions, 0 of them between members of the same
   clump, so `n_contacts`, `z` and the fabric need no filter.
4. 🔴 **Mass carries the same double count.** Clump mass 0.62056151 against a union mass
   of 0.52359878 — the same 1.185185. The material density must be scaled by
   `(π/6)/(2(ds/d_eq)³·(π/6))` per α, which is a constant independent of grain size, so
   one density per α fixes it.

Both defects are of the same kind, and neither would have raised an error.

---

## 3. Gates, registered

**S0 — α = 1 REPRODUCTION, runs first, exits non-zero.** At α = 1 the clump path must
return `phi`, `z_all`, `z_nr` and `rattler_frac` **bit-identical** to the existing sphere
records at the same (family, CV, seed, `p_hat`, N). `phi` is the gate (continuous
channel). This is the strongest available check on the whole construction: it tests the
volume bookkeeping, the density scaling, the insertion RNG order and the collider path in
one number, and an off-by-anything in §2 breaks it.

**S1 — the PSD must not move.** The realised `cv_realised`, `d32_realised`, `d43_realised`
recorded at α > 1 must equal those at α = 1 to 1e-12, at every (family, CV, seed). If the
construction is volume-preserving this is exact; it is checked rather than assumed.

**S2 — the cell must still be legal.** The longest extent is `span_factor(α)·d_eq`, which
is `α·((π/6)/g(α−1))^(1/3)` and reaches **2^(2/3) = 1.5874** at α = 2. N = 1500 clears the
periodic collider at α = 1 and **does not automatically clear it at α = 2**; the
clearance check is re-run per α and the required N is reported, never worked around by
tightening the truncation (which would move `d32`, the one moment this project turns on).

**S3 — the primary claim.** Per family, the paired CV 0 → 0.5 shift in φ_J at each α. The
registered question: **does the shape-of-the-distribution effect (the existing 4.47×
across families) survive, grow or vanish as grains stop being spheres?** No direction is
predicted. The grid is 2-D by construction, so the answer is read as a surface and not as
a path.

**S4 — the confound, quantified not caveated.** Every α > 1 result is reported with its
`eps`, and every S3 number is labelled **"elongation + interlocking"** until the LevelSet
control exists. A single-number claim about "aspect ratio" from clumps alone is
pre-registered as **not quotable**.

---

## 4. What this design cannot decide

Nothing about **absolute size** — the static DEM is proved scale-free, and a size axis
means a different model. Nothing about **concave or faceted** grains: a two-sphere clump
is a fixed shape family with one parameter, and α > 2 is a different family (the spheres
separate) needing more sub-spheres. And nothing about the **direction** of any fabric
anisotropy: orientations are drawn uniformly on the sphere (`shape.random_axes`, checked
against a θ–φ sampler that scores `n·aniso²` above 50 where uniform scores under 6), so
the pack has no built-in preferred axis to find.


---

## 5. AMENDMENT, 2026-09-18 -- written BEFORE any production cell, after a second probe

**Provenance, so this cannot be read as a post-hoc loosening.** At the moment this block
was written the only shards in `out/shape/` were the three `shpgate` cells of 5.2 below;
`out/shape/` contained no record with a production tag, and `grade_shape.py` refuses to
admit `shpgate` or `shptime` records into any S3 number. The file's own rule -- do not
edit after the first production cell -- is therefore satisfied, and the reason for the
amendment is a MEASUREMENT, not a result anybody liked or disliked.

`scripts/probe_clump_s0.sh` asked the mirror of the question 2 asked. 2 asked whether
clump members contact EACH OTHER (they do not). It did not ask whether an EXTERNAL
neighbour contacts BOTH members of one clump. It does, and that single fact moves three
things in this file.

### 5.1 THE CONTACT-MULTIPLICITY CONFOUND -- a third knob on this axis, not a caveat

Two clumps in contact can raise up to **four** member-member interactions. Measured
(`probe_clump_s0.sh` Q2, 200 clumps, seed 7, 200 steps):

| alpha | real interactions | distinct grain pairs | duplicates | `dup_frac` |
|---|---|---|---|---|
| 1.00 | 100 | 25 | 75 | **0.7500** |
| 1.25 | 28 | 22 | 6 | 0.2143 |
| 1.50 | 25 | 23 | 2 | 0.0800 |
| 1.75 | 30 | 26 | 4 | 0.1333 |
| 2.00 | 21 | 20 | 1 | 0.0476 |

**These are a LOOSE scene (z = 0.2-1.0) and are not the production numbers.** They
establish that the effect exists and is alpha-dependent; its size in a jammed pack is a
per-cell measurement, recorded as `dup_contact_frac` on every record, and is read from
the production shards and nowhere else.

Two consequences, and they are different in kind.

- **The COUNTING is repaired.** `packing.contact_graph` now builds the graph over
  DISTINCT GRAIN PAIRS, so `n_contacts`, `z_all`, `z_nr`, the rattler set and the
  backbone are all grain-level quantities and mean what their names say. The raw
  interaction count survives as `n_contacts_raw`. For a sphere pack the two are equal
  body for body, which gate S0a checks rather than assumes (measured: 2762 / 2308 / 2947
  raw = distinct on the three gate cells).
- **The MECHANICS cannot be repaired.** Two overlapping members touching the same
  neighbour are two springs in parallel, so the contact stiffness between that grain pair
  is doubled no matter how the graph is counted. **This axis therefore moves elongation,
  surface corrugation AND effective contact multiplicity together**, where 1 named only
  the first two. S4's label is extended accordingly: every alpha > 1 result is reported as
  **"elongation + interlocking + contact multiplicity"** until the LevelSet control
  exists, and `dup_contact_frac` is printed beside every S3 number.

### 5.2 GATE S0 WAS UNACHIEVABLE AS WRITTEN -- restated, with the original bar kept

3's S0 required a clump built at alpha = 1 to reproduce the sphere records
**bit-identically**. No correct implementation could have met that, for two independent
measured reasons:

1. **The body is not the sphere.** At alpha = 1 the two members COINCIDE, so every grain
   pair raises four interactions (`dup_frac` = 0.7500 above). A pack of coincident-member
   clumps has four times the contacts and four springs per contact; it is not the pack
   this project has run.
2. **The timestep moves.** The mass correction scales the material density by
   `1/(2(ds/d_eq)^3)`, and `PWaveTimeStep` reads the MEMBER density. Measured (Q3), as a
   ratio to the equal-volume sphere's: **0.7071 / 0.7449 / 0.7715 / 0.7879 / 0.7937** at
   alpha = 1 / 1.25 / 1.5 / 1.75 / 2. A different `dt` is a different trajectory. (It
   moves in the safe direction at every rung, and `dt` is recorded on every cell.)

A third measurement decides where the bit-identity bar can live at all: **the analytic
union volume does not reproduce `1 - utils.porosity()` bitwise even for spheres** --
Yade sums body volumes in its own order, and Q4 measured a 1-ulp gap (5.55e-17 at
phi = 0.171). So a bar read on the analytic `phi` would have been unachievable for a
third reason having nothing to do with clumps.

**THE CONSTRUCTION IS CORRECTED, NOT THE BAR.** The union of two coincident spheres *is*
one sphere, so **alpha = 1 builds a SINGLE SPHERE**. The ladder then has a **SEAM** at its
first rung rather than a limit -- `dup_frac` jumps 0 -> 0.75 the moment two members exist
-- and **S3's surface must be read with that seam marked**: the alpha = 1 -> 1.25 step is
the one step that mixes shape with the onset of contact multiplicity, and no "trend in
alpha" may be quoted across it without the alpha >= 1.25 sub-trend printed beside it.

S0 is replaced by three gates, all of which can fail:

**S0a -- the solver must not move.** The three cells of `PREREG_FABRIC_2026-09-18.md` F1
(bidisperse CV 0.3, uniform CV 0.3, lognormal CV 0.4; seed 1, `mu_p` = 0.3,
`p_hat` = 1e-3, N = 1500), re-run at alpha = 1 **through the shape runner**, return
`phi_porosity` **bit-identical** to the committed record. `phi_porosity` and not `phi`,
because `phi` is now the analytic union and the ulp above is real. This is a genuine
bit-identity gate on the trajectory: it can fail on the insertion RNG order (the
orientation axes are drawn from a separate seeded stream precisely so that a mistake
there is visible), on the material, on the engines, on the timestep and on the
grain-level rewrite of the contact graph.
**OUTCOME: S0a PASSES 3/3.** `0.6120403164804106` / `0.6101301588958694` /
`0.6339122825028567`, all sixteen digits, and `z_all`, `z_nr`, `rattler_frac` and
`n_contacts` matched too and are printed rather than read as the gate.

**S0b -- the two volume routes agree where they can be compared.** At alpha = 1,
`|phi_union - phi_porosity| <= 1e-13`. Measured: **1.876e-14 / 3.331e-16 / 6.661e-16**
(bidisperse / uniform / lognormal). The bidisperse value is ~85 ulp rather than 1 --
two discrete radii make the two summation orders diverge further than the probe's
continuous PSD did -- which is why the bar is set from the measurement and the three
numbers are printed. For scale, the seed scatter in phi is 2.2e-4, eleven orders larger.

**S0c -- the clump body is the body it claims to be, AT EVERY alpha.** This is the gate
that tests what S0 was supposed to test and, with alpha = 1 a sphere, no longer can. At
each alpha > 1 on one cell per family: `phi_porosity / phi_union` equals
**`2(ds/d_eq)^3`** -- the alpha-dependent double count the 2 probe measured -- to 1e-9,
and the clump mass read back from Yade equals the union mass to 1e-12. Both are
PREDICTIONS with values (2.0000, 1.4629, 1.1852, 1.0449, 1.0000 at the five rungs), not
observations, so a wrong density, a wrong sub-sphere diameter or a wrong `g(u)` breaks
them.

### 5.3 S1, S2, S4 as amended

S1 and S2 stand. S2's clearance check has been **run on the Windows side for the whole
grid before launching** (5.4). S4 stands with its label extended per 5.1.

### 5.4 S2 -- THE CLEARANCE RESULT, and the grid it forces

The longest extent is `span_factor(alpha)*d_max`, so the minimum N scales as
`span_factor^3`: **1.000 / 1.4285 / 2.0000 / 2.8005 / 4.0000**. Computed from the realised
PSDs at `truncate = (0.2, 5.0)`, worst over seeds 1-12, the required N is:

| family | CV | a=1.00 | a=1.25 | a=1.50 | a=1.75 | a=2.00 |
|---|---|---|---|---|---|---|
| lognormal | 0.3 | 410 | 585 | 819 | 1147 | 1638 |
| lognormal | 0.4 | 871 | 1244 | 1741 | 2437 | 3482 |
| lognormal | 0.5 | 1345 | 1921 | 2690 | 3766 | **5379** |
| uniform | 0.5 | 85 | 121 | 169 | 237 | 338 |
| bidisperse | 0.5 | 42 | 60 | 83 | 117 | 166 |

**The whole constraint is the lognormal tail**; uniform and bidisperse clear at N <= 338
everywhere. **N = 1500 does not cover the grid**: it fails lognormal CV 0.5 at every
alpha >= 1.25 and CV 0.4 at alpha >= 1.5 -- i.e. exactly the corner that carries this
project's effect. A ragged grid there is not acceptable, and **N may not vary across the
grid**, because phi_J drifts with N and an N that moved with alpha would wear alpha's
label.

### 5.5 THE GRID, RE-REGISTERED -- written before the first production cell

`N = 5400` for the whole grid (the clearance demands 5379 at the binding cell; 5400 is
that rounded up). N is a SINGLE value across every cell: phi_J drifts with N, so an N
that moved with alpha would wear alpha's label.

**What one cell costs, measured rather than estimated** (`scripts/run_shape.sh time`,
lognormal, seed 1, `mu_p` = 0.3, wall-clock seconds on one core):

| N | alpha | CV | `p_hat` | iterations | wall | `dup_contact_frac` | `z_all` |
|---|---|---|---|---|---|---|---|
| 1700 | 1.00 | 0.0 | 1e-3 | 10531 | 82.5 s | 0.000 | 4.969 |
| 1700 | 1.00 | 0.0 | 3e-5 | 51976 | 239.7 s | 0.000 | 4.248 |
| 1700 | 2.00 | 0.0 | 1e-3 | 15671 | 232.7 s | **0.130** | **6.341** |
| 5400 | 1.00 | 0.0 | 1e-3 | 19396 | 328.4 s | 0.000 | 4.921 |

Three things follow and each one moves the design.

1. **The bottom of the pressure ladder is the expensive rung, not the top**: 3e-5 costs
   2.9x 1e-3 at the same cell, because a softer target needs five times the iterations.
2. **A clump cell costs ~2.8x its sphere twin** at the same N and pressure -- twice the
   bodies and a longer relaxation.
3. **N = 5400 costs ~4x N = 1700** at the same corner (328.4 / 82.5 = 3.98 against an N
   ratio of 3.18).

Together those put the *whole* registered grid (5 alpha x 6 CV x 3 families x 3 seeds x
4 pressure rungs = 1080 cells) at roughly 17 core-hours per worker-equivalent, i.e. far
beyond one session on 20 workers. **The grid is therefore SAMPLED, and what was cut is
recorded here rather than discovered later from a ragged artifact.**

**REGISTERED PRODUCTION GRID, in two waves.**

- **Wave 1** (tag `shp`): `alpha in {1.0, 2.0}` x `CV in {0.0, 0.5}` x 3 families x
  seeds {1, 2, 3} x `p_hat in {1e-3, 3e-4, 1e-4, 3e-5}` = **144 cells**.
- **Wave 2** (same tag, resumed): `alpha = 1.5` on the same CV / family / seed /
  pressure set = **72 cells**, giving a three-rung alpha ladder.

**WHAT THE SAMPLING GIVES UP, stated before the result is known.**

- **CV is sampled at its two ENDPOINTS only.** S3's registered quantity is the paired
  CV 0 -> 0.5 shift, which this measures directly -- but the descriptive sweep's own
  correction (`REPORT_SWEEP_2026-09-18.md`: bidisperse's phi_J peaks at CV = 0.3 and
  falls back, taking a 4.47x endpoint ratio to 2.89x on the peak) means an endpoint
  shift can misstate a non-monotone family. **No claim of the form "the CV effect at
  alpha X is Y" may be quoted from this grid without that limitation attached**, and a
  non-monotonicity in CV is invisible to it BY CONSTRUCTION.
- **Three alpha rungs, one of which is the SEAM.** alpha = 1 is a sphere and alpha > 1
  is a clump, so the 1 -> 1.5 step carries the onset of contact multiplicity. The
  clump-only sub-trend (1.5 -> 2.0) is the part free of the seam and is reported
  separately. With wave 2 that sub-trend is two points, which can show a direction and
  not a curvature.
- **Three seeds.** A slope inside 2 se is UNDETERMINED at this replication, not absent;
  the grader prints the achieved resolving power beside every null.

**If wave 2 does not run, the campaign reports a two-rung alpha axis and says so.** A
partial grid is a real result provided the holes are named; a partial grid presented as
a complete one is not.

### 5.6 THE PRESSURE LADDER, RE-REGISTERED -- the bottom rung is not a measurement at N = 5400

**This is a defect found by running, not a budget cut, and the evidence is the reason.**
The first production cells at `N = 5400` show that the softest rung of the registered
pressure ladder does not reach a jammed state at the top of the CV ladder within this
project's iteration cap. Measured (lognormal, `mu_p` = 0.3, `alpha` = 1, N = 5400):

| CV | `p_hat` | iterations | `phi` | `z_all` | converged | wall |
|---|---|---|---|---|---|---|
| 0.0 | 1e-3 | 19396 | 0.5977 | 4.921 | True | 482 s |
| 0.0 | 3e-4 | 24601 | 0.5896 | 4.583 | True | 875 s |
| 0.0 | 1e-4 | 33786 | 0.5820 | 4.314 | True | 1444 s |
| 0.0 | 3e-5 | 74871 | 0.5770 | 4.264 | True | 1970 s |
| 0.5 | 1e-3 | 119346 | 0.6532 | 2.726 | True | 2017 s |
| **0.5** | **3e-5** | **200000 (CAP)** | **0.4191** | **0.023** | **False** | 2086 s |
| **0.5** | **3e-5** | **200000 (CAP)** | **0.4302** | **0.023** | **False** | 1213 s |

`z_all` = 0.023 is not a slow packing, it is **no packing at all** -- the cell is still a
gas when the cap stops it. The same rung converges at `N = 1500`, so this is a
consequence of the N the collider clearance forces, not of the pressure.

**Why it cannot simply be left to the converged filter.** The failure is not spread over
the grid: it is at **CV = 0.5 and not at CV = 0.0**, i.e. it censors the ladder
PREFERENTIALLY AT ONE END OF THE AXIS UNDER TEST. A phi_J extrapolated from four rungs at
CV 0 and from three at CV 0.5 is two different estimators wearing one name, and the
paired difference D would carry that difference rather than the physics. This project has
already recorded that failure mode twice (the frictionless drop rate running 0.92 -> 0.33
across CV; the fabric campaign's seven drops all at CV = 0.5).

**REGISTERED: the shape campaign's pressure ladder is `p_hat in {1e-3, 3e-4, 1e-4}`**,
MATCHED across every cell of the grid. `grade_shape.py` filters to that declared ladder
rather than to whatever happens to be on disk, so a stray 3e-5 record cannot lengthen one
cell's ladder and not another's.

**WHAT THIS COSTS, stated rather than discovered later.**

- **The free-exponent anchor is GONE.** `fit_phi_J` needs four points to fit
  `phi = phi_J + c*p_hat^b`, and `b` is the contact law's own prediction (`b` near 1 for
  Yade's linear `CundallStrack` law) -- the check that the ladder is in the critical
  regime at all. With three rungs it returns NaN. **Every phi_J in this campaign is
  therefore a linear extrapolation with NO independent exponent check**, and that is a
  weaker measurement than the sphere campaign's. It is comparable ACROSS this grid, which
  is what S3 needs, and it is not comparable in kind to a four-rung phi_J elsewhere in
  the project.
- The ladder spans one decade (1e-3 to 1e-4) instead of 1.5, so the extrapolation to
  `p_hat` -> 0 is over a shorter lever arm. The linear fit's own standard error is
  reported per cell and propagated into D.
- Three points leave one residual degree of freedom, so the fit can report a standard
  error but cannot diagnose its own functional form.

**The alternative that was NOT taken:** raising `PJ_MAX_ITER` would make the capped cells
converge eventually, at more than 200 000 iterations each and an unknown multiple of an
already 2000-second cell. That is a change to the protocol's convergence criterion made
to rescue a rung, which is the shape of a bar being widened; a matched, shorter ladder
that every cell can actually reach is the stricter choice.

### 5.7 WHAT A REFEREE FOUND IN THIS DESIGN, before the campaign finished

An independent referee saw sections 1-4, this amendment, `shape.py`, `packing.py`,
`scripts/probe_clump_s0.sh` and `grade_shape.py`, and was asked to refute six specific
claims. Four findings survived checking. All four are recorded here with what was done
about each; **none of them is repaired by widening a bar.**

**(i) GATE S0c IS NOT THE ANCHOR ITS OWN TEXT CLAIMED.** S0c compares
`phi_porosity/phi_union` against `2*(ds/d_eq)^3` -- and BOTH sides call
`shape.sub_sphere_diameter`, while `density_scale` is DEFINED as the reciprocal of the
same quantity. The identity therefore holds for any value that function returns, so a
wrong coefficient in `g(u)` would move the construction and the prediction together and
**S0c would still pass at 1e-9**. This is the two-anchors-sharing-one-identity failure
(lessons.md 2026-09-05), and the gate's text said the opposite.
**Done:** the gate KEEPS its bar and its numbers -- it does decide something real, namely
that Yade built members of the radius the code computed, that the clump assembly
preserved them, that the density was applied and that `utils.porosity()` sums member
volumes and nothing else -- and its text now states what it cannot decide. The closed
form is anchored OUTSIDE the pipeline by quadrature over the union's cross-sections
(`test_shape_campaign.py::test_clump_volume_against_independent_quadrature`, seven alpha
values, agreement better than 1e-10), which imports none of `shape.py`'s constants and is
SHOWN ABLE TO FAIL on a 1 % perturbation of `g(u)`'s coefficient.

**(ii) A FOURTH CONFOUND: THE CONTACT STIFFNESS SCALE.** Yade's contact law never sees
`d_eq`; it sees the MEMBER radius `ds/2`, and `ds/d_eq` runs **1.0000 -> 0.7937** along
the ladder. The normal stiffness of a contact therefore falls with alpha at fixed `p_hat`,
so `p_hat = p/E` does not hold the distance from the rigid-grain limit fixed across the
axis.
**Why phi_J is nonetheless first-order safe, and how that is CHECKED rather than argued:**
`ds/d_eq` is a CONSTANT within a rung (it does not vary from grain to grain), so the
effect is exactly a constant rescaling of the pressure axis. For the linear extrapolation
`phi = phi_J + c*p_hat` a rescaling `p_hat -> k*p_hat` changes `c -> c/k` and **leaves the
intercept unchanged**. The prediction that follows is sharp and is measured rather than
asserted: **the fitted slope `c` should rise by a factor `(d_eq/ds) = 1.2599` from
alpha = 1 to alpha = 2** if the rescaling is the whole story. `grade_shape.py` reports the
measured slope ratio beside that prediction. A slope ratio far from 1.2599 means the
stiffness change is doing something the rescaling argument does not cover, and the S3
numbers would then have to be read with that open.
**What is NOT recovered:** the check that the ladder sits in the near-jamming regime
EQUALLY WELL at both alphas is the free exponent `b`, and amendment 5.6 removed the
fourth rung that `b` needs. So this confound is bounded by an argument plus a slope test,
not by the exponent anchor. Recorded as a limit, not resolved.

**(iii) THE CLUMP'S MOMENT OF INERTIA IS WRONG WHERE THE MEMBERS OVERLAP.**
`density_scale` matches the clump's MASS to the union mass, but Yade assembles a clump's
inertia as the sum of the members' own inertias plus parallel-axis terms -- which weights
the overlap lens twice. After mass matching, the referee's quadrature puts the transverse
inertia **up to 15.6 % low, worst near alpha ~ 1.3, and 11.9 % low at alpha = 1.5**.
**The wave-1 grid is UNAFFECTED**: at alpha = 1 the grain is a single sphere and at
alpha = 2 the members are TANGENT, so the overlap lens is empty and the assembly is exact
at both rungs (ratio 1.0000). **Wave 2's alpha = 1.5 is exactly where the distortion
peaks**, which is the one rung that would supply curvature. Registered here BEFORE wave 2
runs: **wave 2 may not be quoted without either an inertia override or a paired
corrected-vs-uncorrected check at alpha = 1.5.**

**(iv) SECTION 5.6's EVIDENCE DOES NOT COVER THE RUNGS IT RETAINED.** The table there
shows CV = 0.5 / N = 5400 converging at `p_hat` = 1e-3 and failing at 3e-5. It shows
NOTHING at 3e-4 or 1e-4 -- the two intermediate rungs the retained ladder depends on --
and the 1e-3 cell already used 119 346 of the 200 000-iteration cap. **If either
intermediate rung also fails at CV = 0.5, the shortened ladder reproduces the exact
pathology it was introduced to remove**, and `grade_shape.py` would drop those cells
(it requires EXACTLY the registered rungs) rather than flag them as the same defect.
**Done:** the grader now reports the converged fraction BROKEN DOWN BY (CV, alpha, rung)
rather than pooled, and refuses to issue an S3 verdict while any (CV, alpha) cell of the
grid is missing its full ladder for that reason. The rungs are being measured by the
production run itself; the claim in 5.6 is narrowed to what was actually shown.

**Three claims the referee tried and could not break**, recorded because a design's
survivals are as informative as its defects: the alpha = 1 single-sphere construction and
the seam's characterisation (A); the PSD's invariance along alpha, confirmed by reading
the code path -- `radii` are drawn before alpha enters anything (B); and the paired
per-seed estimator, confirmed to share bit-identical `d_eq` draws AND insertion positions
across alpha, with CV 0 vs CV 0.5 correctly treated as unpaired, which is the
conservative direction (D). The referee also independently verified `g_of_u` by Monte
Carlo and quadrature at six values of u.

### 5.8 THE CV = 0 RUNG IS ONE CELL, NOT THREE -- and that is a correlation before it is a saving

Found while the campaign was running, from the records already on disk. At CV = 0 every
PSD family is the SAME delta distribution, so `psd.make_psd` returns bit-identical radii
for all three (`diag_shape_cv0.py`: sum 2700.0, min = max = 0.5, arrays bit-equal), and
the **nine CV = 0 cells that had been run in two families returned `phi` BIT-IDENTICALLY,
max |dphi| = 0.000e+00**. The discriminating half holds too: at CV = 0.5 the three
families genuinely differ (lognormal sum 2682.14 / min 0.1001 / max 2.1943, uniform
2695.78 / 0.0671 / 0.9328, bidisperse 2700.0 / 0.25 / 0.75), so the identity at CV = 0 is
a fact about the MODEL and not the family flag failing to do anything.

**Registered as gate S5**, which requires both halves -- identity at CV = 0 AND
difference at CV > 0 -- and which fails if either is missing. The grader aliases the
CV = 0 rung into the families that did not run it, and prints how many values it aliased.

**THE SAVING.** A third of the registered wave-1 grid was the same cell run three times.
The second and third families are therefore run at **CV = 0.5 only**: 18 cells each
instead of 36, and the nine already-run duplicates are kept as the gate's evidence rather
than discarded.

🔴 **THE CONSEQUENCE THAT IS NOT A SAVING, and it applies to results already in this
project.** Because all three families share their CV = 0 reference point EXACTLY,
`D = phi_J(CV_hi) - phi_J(CV 0)` for the three families is **not three independent
measurements**: they share a common subtrahend. Not merely "correlated" -- it is the same
floating-point number, so the shared component is perfectly correlated.

⚠ **THE FIRST VERSION OF THIS PARAGRAPH HAD THE SIGN BACKWARDS** and is corrected here
rather than quietly replaced. It said a cross-family comparison "understates its own
uncertainty if it treats the three as independent". **It overstates it.** One line of
algebra, which is what should have been written instead of an intuition about positive
correlation. Write `D_A = X_A - C` and `D_B = X_B - C` with `C` the shared reference:

    D_A - D_B = X_A - X_B          -- C cancels EXACTLY
    Var(D_A - D_B)          = Var(X_A) + Var(X_B)
    Var(D_A) + Var(D_B)     = Var(X_A) + Var(X_B) + 2*Var(C)      (the naive form)

so assuming independence ADDS `2*Var(C)` that is not there. Numerically, on the referee's
worked case, 1.338 against a correct 0.667 -- **a factor of two too wide**. The same sign
holds for the heterogeneity F statistic (a constant common to all three families does not
move the between-family variance but does inflate every denominator) and, to first order,
for a RATIO of two `D`s, since the covariance is positive there too.

Two readings follow.

- **S3's per-family paired `Delta` is unaffected**, and for a stronger reason than
  correlation bookkeeping: it differences the same reference point out on BOTH alpha
  rungs, so the shared term cancels identically rather than statistically. ⚠ A
  CROSS-FAMILY comparison of `Delta` -- the natural next question once "does the 4.47x
  survive" is asked of more than one family -- carries the same shared term as `D` does,
  and inherits the same caveat. Nothing in this campaign computes one; if a later one
  does, it must.
- **The descriptive sweep's 4.47x across-family spread in phi_J rests on three numbers
  whose CV = 0 arm is ONE measurement wearing three labels.** The direction of the effect
  is CONSERVATIVE, so no verdict there becomes less safe; but the reported uncertainties
  are wider than the data supports, and the exact size depends on the formula behind that
  particular statistic, which this prereg has not read. Recorded as a question for the
  sweep's own artifact, not settled here.
