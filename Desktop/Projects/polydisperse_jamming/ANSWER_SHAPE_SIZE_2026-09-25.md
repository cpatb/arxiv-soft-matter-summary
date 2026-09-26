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

### 3a. Sphericity does not make jamming easier; it makes the window as narrow as it can be
[filled from the registered campaign, §5]

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
i.e. surface chemistry. **So: size-matched cornstarch should jam at least as readily as silica (geometry),
but on this estimate it will not form the same undrained, weak-network ball unless its contacts are softened.**
⚠ This half is a scaling estimate on the project's poroelastic model; it is not a DEM or DNS result.

## 4. What the bench could do to separate the two halves
1. Fractionate cornstarch to CV < 0.15 and measure φ₀, φ_m on the rheometer: the relation predicts χ below the
   sphere value (window wider than Stöber's) and φ₀ lower than the native powder.
2. The same fraction in glycerol–water (µ_w ×10–100): De moves into the silica range, S does not. If the drop
   then forms a ball, S is not the discriminant and the (De, S) reading is wrong.

## 5. The registered test of the shape half — [filled from out/grade_rollchi.txt]
