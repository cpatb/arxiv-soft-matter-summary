# Shape, size and the jammed ball: a geometric answer

2026-09-25. For the author's two questions: (1) why do polydisperse silica spheres fail to jam solidly where
monodisperse ones do? (2) would filtered, size-matched cornstarch jam like a ball, or does the sphericity of
Stöber silica matter? No new bench data. DEM numbers are from this project's Yade packs; every number is
printed by a script named beside it. Derivation: `derivations/geometric_window_2026-09-25.md`.
Registered campaign: `PREREG_ROLLCHI_2026-09-25.md` → `out/grade_rollchi.txt`.

## 1. The geometric picture in one relation

A random pack jams when its contacts supply enough constraints to make it rigid. Two geometric facts set the
packing fraction at which that happens:

- **Counting** sets how many contacts per grain are needed, Z_iso = 2 n_f / n_c: 6 for frictionless spheres,
  10–12 for frictionless non-spheres (their rotations are real degrees of freedom), and 4 for frictional grains
  of ANY shape (friction constrains rotation for everyone).
- **Filling** sets how much space each contact costs. Song, Wang & Makse (Nature 453, 629 (2008),
  doi:10.1038/nature06981) derive φ = Z/(Z + 2√3) for monodisperse spheres from Voronoi volumes. This project's
  packs follow it with a population-dependent filling factor w:

  **φ_J = Z / (Z + 2√3 w)**

  Monodisperse spheres need no fitted number (w = 1.00 frictionless, 0.99 at μ_p 0.3). A size distribution lowers
  w (lognormal CV 0.5: 0.76); an elongated grain raises it (doublet: 1.37); friction barely changes it
  (median spread 0.026 across friction). w drifts upward with the DEM pressure rung, so it is compared only at a fixed rung (p̂ 3e-5). As a prediction it is approximate: φ_m is predicted within ±0.014,
  biased +0.006 (`derivations/filling_factor_posthoc_2026-09-25.txt`, POST HOC).

The shear-jamming window is χ ≤ x < 1 with χ = φ_m/φ₀ and x = φ/φ₀. From the relation,
χ = Z_m(Z₀ + 2√3 w)/[Z₀(Z_m + 2√3 w)], so the window is opened by anything that lowers Z_m/Z₀ and narrowed
by better filling (∂χ/∂w < 0).

## 2. Question 1 — why polydisperse spheres stop jamming

**Size is a filling knob, not a counting knob.** Small grains sit in the voids of large ones and add solid
volume without joining the rigid network; at μ_p 0.6 the load-bearing backbone's filling factor is flat across
CV (1.400 / 1.395 / 1.406 at CV 0 / 0.2 / 0.5, POST HOC), i.e. the whole CV effect there is rattler volume
(consistent with rattling ∝ a², `derivations/rattler_closure_2026-09-24.md`).

At a fixed weighed φ that has two consequences, both against jamming:
1. **x falls.** w 1.00 → 0.76 (CV 0 → 0.5) raises φ₀ 0.632 → 0.694, so x drops by 9 %.
2. **The window narrows.** ∂χ/∂w < 0: measured χ 0.9109 → 0.9170 (N 1500, lognormal, p̂ 3e-5).

A monodisperse suspension at x = 0.95 (inside [0.911, 1)) becomes x ≈ 0.864 at CV 0.5, below its χ of 0.917. It can
then only thicken, however hard it is hit. The DNS drop shows the consequence (`cvmap`, EPS_MIN 0.03): at fixed
weighed φ, CV 0 jams and arrests by t/τ ≈ 0.1 while CV 0.5 cannot jam and spreads 35–67 % further before
retracting; the window carries 99.5–101.8 % of the change. ⚠ At μ_m 0.3 the claim "no weighed φ lets CV 0.5
jam" rests on a 0.007 margin, 1.2× its uncertainty (ANSWER_2026-09-23 §3b).

## 3. Question 2 — sphericity and cornstarch

### 3a. Sphericity matters — but the sign depends on WHICH non-sphericity
- **The geometric rule (counting).** The shear-jamming window is the gap between how many contacts a grain needs to
  be locked without friction and with it. A frictionless sphere needs 6: its rotations are free, since a normal force
  through its centre exerts no torque. With friction every shape needs 4. Shape can move the frictionless count
  either way:
  - **it can add rotational freedom** that frictionless point contacts must pin: elongation, a concave neck. The
    frictionless count rises, the gap grows, and **the window widens**;
  - **it can let a single contact pin rotation by itself:** a flat face resists tipping even without friction, like
    a box on a table. The frictionless pack then jams "frictional-like", the gap shrinks, and **the window narrows**.
- **Measured, elongation/concavity (two-sphere clumps and doublets, μ_p 0.3).** At the jamming rung (p̂ 3e-5, 3 seeds)
  1 − χ = 0.092 (sphere) → 0.131 (α 1.1) → 0.131 (α 1.75) → 0.138 (α 2); pure doublets ×1.49. Most of the widening
  is already there at 10 % elongation. ⚠ At α 1.1, frictionless neighbours nestle into the clump's shallow neck (32 %
  double contacts), so this family mixes elongation with concavity.
- **Measured, convex facets (level-set rounded cubes; PILOT: N 100, 2–3 seeds, p̂ 1e-3).** Window 0.07–0.10 (sphere)
  → 0.08–0.09 (ε 0.7, within scatter) → 0.04–0.05 (ε 0.5). The frictionless φ₀ falls by 0.027 while φ_m holds.
  Faceting narrows the window. This reverses the expectation I registered, and the earlier stand-in (spheres with a
  contact moment limit, which lower φ_m) does not describe convex facets.
- **The onset.** Mild convex elongation makes grains pack *denser* on both branches (φ_m 0.575 → 0.599 at α 1.5;
  frictionless φ₀ up to ~0.70, as for Donev's ellipsoids and Baule's dimers). A closing neck lowers φ_m continuously
  (0.583 → 0.541 from α 1.75 to 2, registered N1 GRADUAL), because filling worsens.
- **So Stöber spheres are not at an extreme of the window: they sit between the elongated/concave family (wider) and
  the faceted family (narrower).** Cornstarch granules are faceted AND irregular, so on geometry alone its window
  could go either way. The published same-paper fits (cornstarch χ 0.76, 1.5 µm silica spheres χ 0.83; Han et al. 2018,
  Royer et al. 2016) say wider, but they differ in lab, size, friction and φ convention and cannot isolate shape.

### 3b. What filtering cornstarch would change
- **Narrowing CV** raises w toward the monodisperse value: φ₀ falls, x rises at fixed φ, and the window
  widens slightly. By the relation, a filtered cornstarch shear-jams *more* readily than the native powder.
- **Mean size** barely moves under fractionation (the native range is ~5–20 µm; the specific granule-size
  and shape numbers for cornstarch are [UNVERIFIED], the citation check found no single primary source).

### 3c. Why a jammed cornstarch drop is still not a silica ball: drainage and contact stiffness, not geometry
`derivations/ball_hybrid_scaling_2026-09-25.md` (a scaling estimate) puts the two lab materials in different
corners of drainage Deborah number De and dilatant-strength ratio S: 800 nm silica De ≈ 7×10⁵, S ≈ 0.004;
10 µm cornstarch De ≈ 7, S ≈ 4. `derivations/cornstarch_fraction_2026-09-25.txt` then asks what a bench can do:
- **size fraction at fixed chemistry** (5–20 µm): De 29 → 1.8 (×20 either way across the particle-modulus
  bracket), **S unchanged at 4.1** (Hertz K_p is size-free at fixed confinement);
- **solvent viscosity** (water → 1 Pa s at 10 µm): De 7 → 7×10³, S unchanged;
- **sphericity** enters drainage only as ψ² ≥ 0.64 (Kozeny–Carman), against ×100–220 from size.
Neither handle reaches the silica corner, because S ≪ 1 needs a soft (double-layer-like) contact network,
i.e. surface chemistry. **So: both shape mechanisms this project can simulate open the window, and neither has been shown to transfer
to a real convex faceted grain; on the drainage/stiffness estimate, cornstarch will not form the same undrained,
weak-network ball unless its contacts are softened.**
⚠ This half is a scaling estimate on the project's poroelastic model; it is not a DEM or DNS result.

## 4. What the bench could do to separate the two halves
1. Fractionate cornstarch to CV < 0.15 and measure φ₀, φ_m on the rheometer: the relation predicts χ below the
   sphere value (window wider than Stöber's) and φ₀ lower than the native powder.
2. The same fraction in glycerol–water (µ_w ×10–100): De rises to ~70–700 — firmly undrained, though still far
   below silica's 7×10⁵ — while S does not move. If the drop then forms a ball, S is not the discriminant and the
   (De, S) reading is wrong.

## 5. Bottom line for the two questions

1. **Polydisperse spheres fail to jam solidly** because a size distribution is a *filling* change, not a
   *constraint* change: small grains fill voids and rattle, raising φ₀ (w 1.00 → 0.76 at CV 0.5) and slightly
   narrowing the window, so a fixed weighed φ lands below χ. They can still thicken; they cannot shear-jam.
2. **Filtered cornstarch:** narrowing its CV raises x, which favours jamming (the size result, measured on
   spheres). On shape the project's two DEM families point OPPOSITE ways — elongation/concavity widens the
   window, convex faceting narrows it (pilot) — so whether a size-matched cornstarch shear-jams more or less
   readily than silica is not decided by geometry here; it needs the rheometer's χ on the fractionated
   powder. What would stop it looking like the
   silica *ball* is not shape but drainage and contact stiffness: at ~10 µm it is near the drained/undrained
   crossover with a stiff, strong contact network, and fractionation barely moves either. The sharpest bench test
   is fractionated cornstarch in glycerol–water, which moves drainage toward silica's range (De ×10³ at 1 Pa s, still ~2 decades short of silica's 7×10⁵) and
   leaves the network strength where it is.

## 6. The registered DEM test (`PREREG_ROLLCHI_2026-09-25.md`, `out/grade_rollchi.txt`)
Partly complete (shared machine): R3 PASS, R7 (p̂ 1e-3) NOT CLUMP-LIKE, N1 GRADUAL, R2 FAIL; the rest INCOMPLETE.
Details, caveats and every number in `derivations/geometric_window_2026-09-25.md` §6. Figure R1
(`make_rollchi_figures.py`, `FIGURES_ROLLCHI.md`).
