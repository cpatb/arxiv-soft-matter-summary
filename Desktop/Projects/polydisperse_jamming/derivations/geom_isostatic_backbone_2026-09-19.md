# Item D — Constraint counting and the backbone

**VERDICT: SUPPORTED as the reason Z_nr is nearly flat, and REFUTED as an exact
statement.** The isostatic count contains no length, so polydispersity *cannot* move
z_iso — the flatness is forced, not a coincidence. The residual CV drift is resolved
(t = −7.25) in the one arm where the data is strongest, and it maps onto a fall in the
contact mobilisation fraction, which is measurable and unmeasured.

**What this anchor can and cannot decide.** The count decides what z_iso must be *given*
marginal rigidity and a mobilised fraction, and decides that grain size cannot enter it.
It CANNOT decide whether the measured drift is a real change in x_m or a failure of
marginal rigidity — **both predict the same z_iso.** The output that separates them, the
mobilised contact fraction, is not recorded.

Numbers printed by `geom_isostatic_backbone_2026-09-19.py` →
`geom_isostatic_backbone_2026-09-19.txt`.

---

## 1. z_iso for a frictional sphere packing, derived

Count **static determinacy**: contact forces are determined by balance alone exactly when
independent unknowns equal independent equations.

**Equations.** Each rigid grain gives `d` force-balance plus `d(d−1)/2` torque-balance
equations, so `d(d+1)/2` per grain → `N·d(d+1)/2`. At d = 3: **6N**.

**Unknowns.** Each contact carries `d` force components (1 normal, `d−1` tangential). At a
**fully mobilised** contact the tangential magnitude is pinned to `μ f_n`, removing one
unknown. With `x_m` the mobilised fraction: `N_c(d − x_m)`.

Setting them equal with `z = 2N_c/N` (SymPy, §1):

> **z_iso = d(d+1)/(d − x_m)**, and in d = 3, **z_iso = 12/(3 − x_m).**

**Limits, asserted:**

| x_m | z_iso (d=3) | meaning |
|---|---|---|
| 0 | **4 = d + 1** | nothing sliding, friction fully effective |
| 1 | **6 = 2d** | every contact at the cone, i.e. μ → 0 |

2D check: `d(d+1)/(d−x_m)` at d = 2 gives 3 (= d+1) at x_m = 0 and 6 at x_m = 1 —
correctly reproducing the 2D frictional and frictionless bounds.

This reproduces van Hecke's stated interval (*J. Phys.: Condens. Matter* **22**, 033101
(2010), doi:10.1088/0953-8984/22/3/033101, §IV.1.1 verbatim: *"frictional spheres can
attain a range of contact numbers: d+1 ≤ z_c ≤ 2d"*) and supplies the interpolant.

⚠ **One caveat on the monotonicity, from the source.** The review says z_c reaches d+1
*"only in the limit of μ→∞ and very slow equilibration"*, and its Fig. 17(a) explicitly
marks an arched region of (z_c, μ) combinations that *"while they are not reached in these
numerics, are perfectly possible."* So "z_c decreases with friction" is a property of a
preparation protocol, not a law. What the count above says without qualification is the
weaker and safer thing: **z_c is a one-to-one function of x_m and of nothing else**.
Whether x_m falls with μ is protocol-dependent; that it is the only handle is not.

⚠ **A convention that must travel with this.** van Hecke states the clean scalings hold
only if rattlers are **excluded from the contact count and included in the packing
fraction** — which is exactly what `packing.py:58` does (`phi` INCLUDES rattlers, `Z`
excludes them). The project's convention matches the literature's.

### The point for this project

Every symbol above is a count. Grain diameters appear **nowhere**. So:

> **PREDICTION P1: z_iso is independent of the size distribution at fixed μ.**

Fact 3 ("Z_nr is flat while the rattler fraction triples") is not a coincidence — it is
*forced*, provided (a) the pack is marginally rigid and (b) x_m does not itself depend on
CV. Proviso (b) is the whole item.

---

## 2. z_nr at a pressure is not z_iso

The contact law here is Cundall–Strack **linear** springs
(`packing.py: Law2_ScGeom_FrictPhys_CundallStrack`), so `p ~ (φ − φ_J)` and the standard
marginal scaling `z − z_c ~ (φ − φ_J)^{1/2}` gives `z_nr(p) = z_iso + a√p̂`. That is the
rule `measure_z_iso.py` uses; here it is refitted **per seed** so the drift carries an
error bar at the design's own 2 dof.

| family | μ_p | CV 0 | CV 0.1 | 0.2 | 0.3 | 0.4 | 0.5 |
|---|---|---|---|---|---|---|---|
| lognormal | 0.0 | 5.8531 | 5.8241 | 5.8686 | 5.8712 | 5.8419 | 5.8656 |
| lognormal | 0.3 | 4.5367 | 4.5356 | 4.5096 | 4.4375 | 4.4564 | 4.3878 |
| uniform | 0.0 | 5.8551 | 5.8656 | 5.8568 | 5.8773 | 5.9120 | 5.8352 |
| uniform | 0.3 | 4.5395 | 4.5445 | 4.5222 | 4.4673 | 4.4890 | 4.4634 |
| bidisperse | 0.3 | 4.5395 | 4.5389 | 4.5188 | 4.4785 | 4.4798 | 4.5651 |

(seed-mean; worst per-seed r² is 0.839, all fits on 3–5 rungs; see §2 of the .txt for se's)

---

## 3. Is the drift resolved? — and the answer is structured

Seed-paired CV 0 → 0.5, t_crit = 4.303 at 2 dof:

| family | μ_p | Δz_iso | se | t | verdict |
|---|---|---|---|---|---|
| lognormal | **0.0** | +0.0125 | 0.0304 | +0.41 | NOT resolved, P1 stands |
| uniform | **0.0** | −0.0199 | 0.0340 | −0.58 | NOT resolved, P1 stands |
| lognormal | 0.3 | −0.1489 | 0.0367 | −4.06 | NOT resolved (just) |
| **uniform** | **0.3** | **−0.0761** | **0.0105** | **−7.25** | **RESOLVED, P1 REFUTED** |
| bidisperse | 0.3 | +0.0256 | 0.0254 | +1.01 | NOT resolved, P1 stands |

**This is the structurally informative result, and it is exactly what §1 predicts.** At
μ_p = 0 the drift is zero to within 0.02 in both families — because at x_m ≈ 0.95 the
packing is against the frictionless ceiling z = 2d and there is nothing for CV to move.
The drift appears **only at finite friction**, where x_m is a free number. The count says
the *only* quantity through which polydispersity can reach z_iso is x_m, and x_m only
exists when μ > 0. The data agrees.

### The frictionless half is externally confirmed; the frictional half is untested in the literature

van Hecke §III.3, verbatim: *"Δz ~ (Δφ)^{0.50±0.03}, where z_c, the critical contact
number, is within error bars equal to the isostatic value 2d. **Note that this result is
independent of dimension, interaction potential or polydispersity**"* — sourced to O'Hern,
Silbert, Liu & Nagel, *Phys. Rev. E* **68**, 011306 (2003),
doi:10.1103/PhysRevE.68.011306, whose rattler-excluded Table 1 gives Z_c = 2d for 3D mono-
and bidisperse and 2D bidisperse; and confirmed for *continuously* polydisperse spheres by
Ozawa, Berthier & Coslovich, *SciPost Phys.* **3**, 027 (2017),
doi:10.21468/SciPostPhys.3.4.027 (*"the isostatic condition, Z = 6, indeed holds over the
entire J-line"*). **My μ_p = 0 rows are a third, independent confirmation on a CV ladder
rather than on two distributions** — O'Hern's evidence is two size distributions, not a
width sweep.

**The frictional half is a clean gap.** A literature pass on 2026-09-19 found no paper
that varies the PSD, measures z_c(μ) at finite friction, and reports whether d+1 shifts.
What exists is the counting argument itself (van Hecke App. A; Torquato & Stillinger,
*Rev. Mod. Phys.* **82**, 2633, doi:10.1103/RevModPhys.82.2633), which is manifestly
radius-free but is written for congruent spheres and never stated as extended; and
Danisch, Jin & Makse, *Phys. Rev. E* **81**, 051303 (2010),
doi:10.1103/PhysRevE.81.051303, which *applies* z = 6 and z = 4 inside an arbitrary-P(R),
arbitrary-friction model — an assumption, not a test. Silbert (2010) has the z_c(μ) table
but is monodisperse in 3D.

So: **the derivation in §1 is standard but its extension to a polydisperse frictional pack
appears to be assumed everywhere and tested nowhere, and §3 is a test of it.** The honest
counterweight is Maher & Torquato, *Phys. Rev. E* **110**, 064605 (2024),
doi:10.1103/PhysRevE.110.064605: *"rapid increase in the rattler concentration and rapid
drop-off in the fraction of isostatic configurations for ensembles with β > 2.0"* — the
count is radius-free, but whether a protocol *attains* it degrades with size ratio. That
is reading (b) in §4, and it is the reason §4 cannot be settled from z_iso alone.

---

## 4. Inverting the drift for the mobilised fraction

`x_m = 3 − 12/z_iso`. An x_m outside [0, 1] would refute marginal rigidity rather than
measure anything: **0 of 30 cells is inadmissible.**

| family | μ_p | x_m(CV=0) | x_m(CV=0.5) | change |
|---|---|---|---|---|
| lognormal | 0.0 | +0.9498 | +0.9542 | +0.0044 |
| uniform | 0.0 | +0.9505 | +0.9435 | −0.0070 |
| lognormal | 0.3 | +0.3549 | +0.2652 | **−0.0898** |
| uniform | 0.3 | +0.3565 | +0.3114 | **−0.0451** |
| bidisperse | 0.3 | +0.3565 | +0.3714 | +0.0148 |

So μ_p = 0.3 corresponds to roughly a third of contacts sitting at the Coulomb cone, and
μ_p = 0 to 95 % (as it must — a frictionless contact is trivially "mobilised" since it
carries no tangential force to constrain).

**The sign is the puzzle, and it is a real one.** A naive reading runs the other way: a
broader size distribution puts many small grains into the coarse skeleton's pore space;
those contacts carry little normal force; a contact with small `f_n` is *easier* to drive
to the cone. So x_m should RISE with CV and z_iso with it. **The measured drift is the
opposite sign.** Two readings survive, and one output separates them:

* **(a)** x_m really does fall — the lightly loaded small-grain contacts are *not*
  preferentially mobilised, i.e. the extra small grains are being held, not slipping.
* **(b)** the pack is not marginally rigid at CV = 0.5, so z_iso is not a count there and
  the inversion is meaningless.

Yade exposes the mobilised state of a contact directly (`|f_t| ≥ μ f_n`), so **(a) is
settled by one extra scalar in the measure step.** Reading (a) would also be the natural
companion to item B's required size selectivity S ≈ 0.8: if small grains are being held
rather than slipping, they are not rattling, which is exactly the "rattlers are strongly
but not perfectly size-selective" picture.

---

## 5. Three candidate artefacts

**(i) The extrapolation?** No. The raw per-rung drift is present at every rung and does
not shrink towards the lowest:

| family | p̂=1e−3 | 3e−4 | 1e−4 | 3e−5 | extrapolated |
|---|---|---|---|---|---|
| lognormal | −0.189 | −0.200 | −0.172 | −0.150 | −0.149 |
| uniform | −0.081 | −0.090 | −0.082 | −0.071 | −0.076 |
| bidisperse | −0.117 | −0.024 | +0.004 | −0.025 | +0.026 |

**(ii) The rattler threshold?** No. 217 records pass every filter; the drift barely moves
between k = 2 and k = 3 (lognormal at 1e−3: −0.1891 → −0.1404; uniform: −0.0812 →
−0.0500). ⚠ **k = 4 is UNDEFINED at μ_p = 0.3**: `z_nr_k4` is non-finite in 216 of 217
records, because iterated removal of grains with fewer than 4 contacts empties a backbone
whose mean coordination is only ~4.5. That is a property of the frictional arm, not a bug,
and it means the k-sweep is a two-point sweep here.

**(iii) Finite size?** `d_max/L` rises with CV by construction and the shards flag it.
**NAMED, NOT CHECKED** — it needs the N = 3000 block, which exists on disk under a
different tag and was not loaded.

---

## 6. ⚠ A discrepancy with `out/z_iso_measured.txt`

The committed artifact gives **z_iso = 4.609** at CV = 0 for all three families at
μ_p = 0.3 and 4.377 for lognormal at CV = 0.5 — a drift of **−0.232**, hardcoded as a
per-cell lookup table at `analyse_phase8.py:165-169` and quoted in five NEXT_STEPS
documents. The matched-ladder per-seed refit gives:

| family | CV 0 | CV 0.5 | drift |
|---|---|---|---|
| lognormal | 4.5367 ± 0.0089 | 4.3878 ± 0.0422 | **−0.1489** |
| uniform | 4.5395 ± 0.0155 | 4.4634 ± 0.0076 | −0.0761 |
| bidisperse | 4.5395 ± 0.0155 | 4.5651 ± 0.0155 | +0.0256 |

i.e. the lognormal drift is **36 % smaller** than the committed one, and the CV = 0
intercept is 0.07 lower.

The obvious explanation — an unmatched ladder — is **refuted by its own sign.** Exactly
one cell (lognormal, CV = 0, μ_p = 0.3) carries a fifth rung at p̂ = 3e−3 on one seed, and
it *is* the reference from which the drift is measured; but including it moves the
intercept 4.5750 → 4.5537, i.e. DOWN, which would make the drift smaller, and the
committed value is higher still.

What is established: on a ladder identical at every CV, the numbers above; and
`out/z_iso_measured.txt` reports "5 rungs" and 433 cells where this loader finds 4 rungs
and 395 cells, so its selection differs. **Cause not traced.** The candidate worth
checking first is the 2026-09-17 `_it{cap}` duplication
(`sweep_descriptive.py:243-262`): if `measure_z_iso` still counts one physical cell twice,
it duplicates rungs *inside* a single fit. Flagged, not asserted.

`z_iso = 4.535` does **not** appear in `paper/main.tex`; the manuscript quotes `Z_nr` at a
pressure (5.111 at CV = 0, falling by 0.189 / 0.081 / 0.117), so this affects the
strain-clock / shear line of work rather than the submitted manuscript.
