# Beach sand dilatancy: field sampling to lab sphere-impact protocol

**Status:** proposal, 2026-08-06
**Origin:** field observation at Lake Michigan. Wet sand hardened underfoot on flat,
fine-grained stretches and did not harden at all on steeper, coarser stretches.
**Design:** collect sand from stations along a beach, characterise each sample, then run
steel-sphere drop impacts on reconstituted saturated beds under high-speed video.

**Reconciled 2026-08-09** against `DESIGN_FINDINGS_2026-08-07.md`, every number re-derived
from `dilatancy_design.py` (22 tests passing). Changed: (1) porosity now written `n`, `φ`
reserved for packing fraction; (2) impact speed range corrected from 0.5–4 m/s to
**0.2–2.5 m/s**, because the old range never reached the footstep regime it was meant to
model; (3) Q4 promoted from an afternoon's curiosity to a named aim, reframed as *measuring*
ψ; (4) §4 now requires the full grading curve and a Sauter diameter, not three percentiles;
(5) ψ, D and the Kozeny-Carman constant are swept rather than quoted. The §2 table also
picks up ρ_s = 2650 kg/m³ instead of an implied 2000, and the footstep-scale sentence in §2
was wrong and is now correct.

### Notation

`φ` is the **packing (solid) fraction** and only that: φ₀ for the as-placed value, φ_c for
the critical value, φ ≈ 0.60 for a natural sand pack. `n` is the **porosity** (void
fraction), n = 1 − φ, so n ≈ 0.40 for the same pack. The two are never interchangeable:
Kozeny-Carman carries n³/(1−n)², so substituting one for the other changes k by
(0.6/0.4)⁵ = 7.6× and moves every derived speed with it. This cost exactly that error once
already, while writing `dilatancy_design.py`.

---

## 1. What is already known

The mechanism is textbook and I am not proposing to rediscover it.

Dense saturated sand must dilate to shear (Reynolds 1885). Dilation opens pore volume that
water has to fill. If the water cannot arrive fast enough, pore pressure goes negative and
Terzaghi effective stress σ′ = σ − p rises, so the pack stiffens. Whether that happens is a
race between the loading time and the drainage time, and drainage time scales as d⁻²
through permeability.

Two lab papers have already isolated the two obvious control variables, each one at a time:

| Paper | Varied | Held fixed | Result |
|---|---|---|---|
| Jerome, Vandenberghe & Forterre, *PRL* **117**, 098003 (2016) | initial packing fraction φ₀ (0.560 vs 0.604) | grain size | sharp liquefaction / transient solidification transition |
| Strader, Causley, Dijksman & Clark, *PRF* **7**, 054304 (2022) | grain size (5 bins, 53 to 850 μm), fluid viscosity (1 to 1412 mPa·s) | φ₀ ≈ 0.59, not varied | a_max ∝ d^−0.8 η^0.4 V^1.6; theory fails above η ≈ 10 cP |

Supporting work: Zhao et al., *PRE* **91**, 022202 (2015) on penetration depth against
saturation; Iverson, *JGR* **110**, F02015 (2005) for the same dilatancy and pore-pressure
feedback in landslides; Rondon, Pouliquen & Aussillous, *Phys. Fluids* **23**, 073301 (2011)
for the loose-versus-dense split in immersed granular collapse.

On the field side, Bascom (1951) established that beach-face slope grows with grain size, so
slope and grain size are confounded on any real beach. The two things I thought I was
observing separately are largely one variable. Capillary fringe height goes as h_c ∝ 1/d,
roughly 50 cm for fine sand and 7 cm for coarse, so the width of the near-saturated firm
band is W ≈ h_c/tanβ ∝ 1/(d·tanβ). Both factors shrink it as grains coarsen, which is why
the transition across a beach feels abrupt rather than gradual.

### The gap

Both lab studies used monodisperse spherical glass beads at a single, lab-prepared packing
state. Real beach sand is polydisperse, angular, and sits at a packing fraction set by wave
loading rather than by a pluviation protocol. Wu et al., *Granular Matter* **23**, 88 (2021)
showed by X-ray CT that naturally-deposited beach sand has a measurably different
microstructure from air-pluviated lab sand at comparable porosity: narrower
coordination-number distribution, inclined rather than horizontal grain orientations. The
27-author perspective *Soft matter physics of the ground beneath our feet* (*Soft Matter*
**20**, 5859, 2024) names this lab-to-field translation as an outstanding problem and states
that monodisperse bead experiments cannot capture heterogeneous natural sediment.

So the open questions this project can actually reach:

**Q1, the best one.** Which grain diameter enters the Darcy-Reynolds law for polydisperse
sand? Permeability in graded sand is controlled by the fines (Hazen: k ∝ d₁₀²), not by the
median. So a_max ∝ d^−0.8 should collapse on d₁₀ or on a Kozeny-Carman effective diameter
rather than on d₅₀, and the discrepancy should grow with the uniformity coefficient
C_u = d₆₀/d₁₀. This is untestable with monodisperse beads by construction, and beaches hand
you a natural spread of sorting for free.

The question is comfortably well posed. For a lognormal grading, the predicted a_max using
d₁₀ against d₅₀ differs by 13% at C_u = 1.2, 59% at C_u = 2.0, and 108% at C_u = 3.0. At 15%
measurement scatter the discrimination threshold is C_u ≈ 1.23, which essentially every real
sand exceeds. But the Kozeny-Carman effective diameter is the **Sauter** mean d₃₂, and d₃₂
sits far closer to the median than d₁₀ does (0.903 d₅₀ against 0.561 d₅₀ at C_u = 2). So the
informative contrast is **Hazen (d₁₀) against Kozeny-Carman (d₃₂)** — a factor 1.46 in a_max
at C_u = 2 — not either one against d₅₀. This is why §4 demands the whole grading curve:
d₃₂ cannot be read off three percentiles.

**Q2.** Does the liquefaction/solidification boundary in the (φ₀, d) plane collapse onto one
dimensionless group, or are the two axes independent? Forterre varied φ₀, Clark varied d,
and nobody has crossed them.

**Q3.** Does grain angularity shift the critical packing fraction φ_c, and therefore the
sign of dilatancy at the packing fraction the beach actually delivers?

**Q4, a named aim, not an afternoon.** Suction is capped near 1 atm by cavitation and air
entry, so the V^1.6 law has to saturate at high enough impact speed. The ceiling sits at
V_cav = k(p_atm − p_vap)/(ηψD), which scales as d². With ψ = 0.1 and D = 20 mm:

| d₃₂ | V_cav | drop height |
|---|---|---|
| 0.10 mm | 0.49 m/s | 12 mm |
| 0.15 mm | 1.10 m/s | 62 mm |
| 0.20 mm | 1.96 m/s | 0.19 m |
| 0.30 mm | 4.40 m/s | 0.99 m |
| 0.50 mm | 12.2 m/s | 7.6 m |
| 1.00 mm | 48.9 m/s | 122 m |

Fine sand should therefore show the departure inside an ordinary bench drop, coarse sand
cannot, and the boundary should move as d². The contrast is the test, and it needs a shelf
rather than a tower — the §5 sweep to 2.5 m/s already covers d ≲ 0.23 mm.

What raises Q4 from curiosity to aim is that the simple model demonstrably *fails* here.
Strader et al. report an unbroken a_max ∝ d^−0.8 η^0.4 V^1.6 up to 2.5 m/s across 53 to
850 μm. Two requirements follow and they contradict each other: the fine end must not have
cavitated (ψD ≤ 1.10 × 10⁻⁴ m) while the coarse end must still have been pore-pressure
controlled (ψD ≥ 4.73 × 10⁻³ m). **The window is empty by a factor of 43.** Since k ∝ d²
the window closes as (d_fine/d_coarse)², so any wide grain-size sweep does this. Taken alone
the cavitation constraint bounds our own loosest number: ψ ≤ 0.0055 for D = 20 mm, at least
18× below the assumed 0.1.

Read that as a statement about the estimate, not about their measurements. A single constant
ψ, with the intruder diameter as the only length scale, cannot span a 16× range in grain
size. Candidate resolutions, all testable: ψ depends on the packing state rather than being
a constant; the relevant length is the dilating shear-band thickness rather than D; or
partial cavitation still transmits stress, so the ceiling is soft rather than hard. So Q4 is
reframed: **measure ψ(φ₀, d)** by locating the departure speed, rather than confirm a
saturation. The beach hands you the φ₀ spread for free. See §7.

Honest expectation: Q1 most likely confirms Hazen, which is a useful measurement and a
modest paper rather than a discovery. Q4 is now the one that can do better than that — a
measured ψ(φ₀, d) that closes the 43× gap would be a real contribution, and unlike Q1 it
cannot come out merely confirmatory, because the incumbent model is already known to be
inconsistent across the range. Venue realism is *Granular Matter*, or *Am. J. Phys.* if
framed pedagogically, or *PRF* if Q1 or Q4 comes out clean. Not a PRL.

---

## 2. Regime estimate, worth doing before buying anything

This section is computed, not asserted: every number below comes from `dilatancy_design.py`
and is covered by a test. Dilation at rate ~ψV over intruder scale D demands a Darcy flux,
requiring suction

Δp ≈ η ψ V D / k,  with k = n³d₃₂²/(180(1−n)²) ≈ 9.9 × 10⁻⁴ d₃₂² (Kozeny-Carman at n = 0.40)

n is the porosity, not the packing fraction — see Notation. ψ is the dilation-angle factor,
nominally 0.1, and it is by a wide margin the loosest number here. Compare Δp against the
inertial grain stress ρ_s V², using the pycnometer grain density (quartz, 2650 kg/m³, until
§4 measures it).

| d₃₂ | k (m²) | Δp at V = 2 m/s, D = 20 mm | ρ_s V² | Δp/ρ_s V² | Behaviour |
|---|---|---|---|---|---|
| 0.2 mm | 4.0 × 10⁻¹¹ | ~101 kPa | 10.6 kPa | 9.6 | suction dominates, hardens |
| 1.0 mm | 9.9 × 10⁻¹⁰ | ~4 kPa | 10.6 kPa | 0.4 | drained, foot goes through |

The Δp column is the suction the dilation *demands*, uncapped. For the fine case that demand,
101 kPa, has already overshot the cavitation ceiling of p_atm − p_vap = 99.0 kPa, so the
suction the pack can actually deliver is 99 kPa and the true ratio is 9.3 rather than 9.6.
The fine-sand row is therefore sitting on the ceiling, not comfortably below it. That is what
Q4 is about, and it is why the row is quoted as an upper bound.

**Only the group D/V matters.** Rearranging, Δp/(ρ_s V²) = ηψD/(k ρ_s V), so at fixed sand
the intruder enters *only* through D/V. Matching the field phenomenon means matching D/V,
not matching V. A footstep has D ≈ 0.10 m and V ≈ 1 m/s, hence D/V = 0.1 s. §5 is specified
to reach it.

**The footstep-scale estimate, done properly.** At D = 0.1 m and V = 1 m/s the uncapped
ratio is 95 for 0.2 mm sand and 3.8 for 1.0 mm sand, and the crossover to draining sits at
d ≈ 2.0 mm. So with ψ = 0.1 the model says *every* sand finer than granule gravel should
harden underfoot — which contradicts the field observation this project started from. The
contradiction is informative rather than fatal, because it is a statement about ψ:

| ψ | grain size at which a footstep stops hardening |
|---|---|
| 0.10 | 1.95 mm |
| 0.03 | 1.07 mm |
| 0.010 | 0.62 mm |
| 0.0055 (Strader bound, §1 Q4) | 0.46 mm |
| 0.003 | 0.34 mm |

At ψ ≈ 0.005 the crossover lands in fine-to-medium sand, where the beach actually puts it.
The field observation and Strader's unbroken power law therefore point the same way, at a ψ
roughly 20× below the textbook 0.1, from two completely independent directions. That is the
strongest argument for Q4 as a named aim.

**Sweep, do not quote.** ψ, D and the Kozeny-Carman constant (180) are the three numbers
this whole estimate rests on, and none is measured here. Every function in
`dilatancy_design.py` takes them as arguments for that reason. Report results as functions
of ψ and D — never as a single number with these buried inside it.

---

## 3. Field campaign

One cross-shore transect per beach, 6 to 10 stations from dune toe to swash. Three or four
beaches spanning grain size. Lake Michigan is fresh water, which removes salinity and clay
flocculation as confounds, a real advantage over an ocean beach.

Work a falling tide or a stable lake level. Record the time and the lake level from the
nearest NOAA gauge, since the water-table position lags.

Per station, record:

1. GPS position, distance from the waterline, time.
2. Beach-face slope β, from a phone inclinometer on a 1 m board, three readings.
3. In-situ packing fraction φ₀, using a thin-walled ring sampler (a cut length of 50 mm
   tube, volume measured by water fill). Press in vertically, excavate, trim flush both
   ends, seal. This is the most important field number and the one that cannot be recovered
   later.
4. Gravimetric water content, giving saturation S. Weigh the sealed ring wet, oven-dry at
   105 °C, reweigh.
5. Depth to water table. Dig a small pit, let it equilibrate 10 min, measure.
6. Field firmness proxy: a fixed-mass rod dropped from a fixed height, penetration read with
   a ruler, 10 repeats at fresh spots. Crude, but it is the direct analogue of the lab
   measurement and it anchors the whole thing. **Record the rod diameter and drop height**,
   not just the penetration: §2 shows the regime is set by D/V, and without those two numbers
   the rod's D/V = D/√(2gH) cannot be computed and the field points cannot be put on the same
   axis as the lab ones. A 10 mm rod from 0.1 m gives D/V = 0.0071 s, deep in the drained
   regime and 14× short of a footstep, so do not expect the rod to agree with a boot; that is
   a reason to record D/V, not a reason to change the rod.
7. Qualitative note on whether the sand whitens and hardens underfoot, plus the paced width
   W of the firm band.
8. Bulk sample, about 2 kg, sealed bag, labelled.

Do not try to preserve the field packing in transit, because it is not recoverable. The
design deliberately measures φ₀ in the field and reproduces it in the lab, sweeping φ₀
around that value. That turns the destroyed-packing problem into part of the experiment
rather than a flaw in it.

Gear: ring sampler tubes, trowel, sealable bags, kitchen scale reading to 0.1 g, phone with
inclinometer, 1 m board, tape measure, drop rod and ruler, sharpie, notebook, cooler. Under
$100 excluding the oven.

---

## 4. Sample characterisation

Per sample, before any impact work:

1. Full grain-size distribution by sieve stack. **Retain and archive the whole grading
   curve — the per-sieve retained masses — not just the percentiles derived from it.**
   Report d₁₀, d₅₀, d₆₀, C_u = d₆₀/d₁₀, Folk-Ward graphic sorting σ, *and* the Sauter
   (surface-area) mean

   d₃₂ = 1 / Σᵢ (wᵢ / dᵢ),

   with wᵢ the mass fraction retained in bin i and dᵢ the geometric mean of that bin's
   sieve openings. d₃₂ is the diameter Kozeny-Carman actually wants, because permeability
   is set by surface area per unit volume rather than by a length percentile, and it is the
   quantity Q1 tests against Hazen's d₁₀. It cannot be reconstructed from d₁₀/d₅₀/d₆₀ after
   the fact, so a discarded grading curve costs the sample. Use enough sieves to resolve the
   curve: successive openings in a ratio of √2 or finer. (Note that the sedimentological
   grain-size scale is also conventionally called "phi"; it is not used in this document,
   precisely to keep φ unambiguous.)

   Check C_u before committing a sample to Tier-2 work. At 15% scatter in a_max the
   hypotheses only separate above C_u ≈ 1.23 for d₁₀ against d₅₀, and above C_u ≈ 1.25 for
   the contrast that actually matters, d₁₀ against d₃₂ (§1 Q1). Both thresholds are
   independent of d₅₀ — sorting is what buys the discrimination, not grain size. A sample
   below C_u ≈ 1.25 is not discriminating, whatever else it is.
2. Grain shape, from a flatbed scan or stereoscope, with `skimage.measure.regionprops` for
   circularity and aspect ratio on at least 500 grains. Standard particle-analysis pipeline.
   Watch the usual traps: segment against a calibrated scale bar, exclude border-touching
   grains, verify against a known-size reticle.
3. Grain density by pycnometer. Do not assume 2650 kg/m³. φ is computed from it, so an error
   propagates straight into the packing state.
4. φ_min and φ_max by the standard loosest and densest protocols, so the field φ₀ can be
   reported as a relative density D_r. D_r is the quantity that compares across sands. Raw φ
   is not.

---

## 5. Lab impact rig

Container: cylindrical, inner diameter at least 10 times the sphere diameter and depth at
least 5 times the largest expected penetration. Confinement measurably changes penetration
(Seguin et al., arXiv:0805.4189) and an undersized container silently biases everything.

Saturation is the number one failure mode. Fill from below with de-aired water, slowly,
letting the front rise through the bed. Trapped air destroys the incompressibility the whole
mechanism depends on, and a bed that looks saturated from above usually is not. Use water
that has been boiled and cooled, or vacuum-degassed. Verify by weighing: measured bulk
density must match the value computed from φ and measured grain density.

Packing preparation: set φ₀ by controlled tapping for dense beds, or by fluidising and
letting the bed settle for loose ones. Measure φ₀ every time by weighing fluid and grains,
exactly as Strader et al. did. Never assume the protocol delivered the target. Sweep φ₀
across a range straddling the field value for that sample.

Free surface: water level exactly at the bed surface. If the sphere has to enter water
first, you have added a water-entry problem on top of the granular one.

Impact: steel spheres, two or three diameters, released from an electromagnet at heights
giving V = √(2gH). **The speed range is 0.2 to 2.5 m/s for the 20 mm sphere, i.e. drop
heights of 2 mm to 0.32 m** — not the 0.5 to 4 m/s an earlier draft of this protocol
specified. The derivation, since the choice is not free:

- *Slow end, set by the field.* §2 shows the regime is controlled by D/V alone, and a
  footstep is D/V = 0.10 m / 1.0 m/s = 0.1 s. The matched lab speed is therefore
  V = D/(0.1 s), which for D = 20 mm is **0.20 m/s**, a 2.0 mm drop. The old range topped
  out at D/V = 0.020/0.5 = 0.04 s, so at its *slowest* it was still 2.5× short of a
  footstep, and at 2 m/s it was 10× short. The rig as previously specified sat an order of
  magnitude further into the drained regime than the observation that motivated the project.
- *Fast end, set by the literature.* §6 requires reproducing Strader et al., whose reported
  range runs to 2.5 m/s. **2.5 m/s**, a 0.32 m drop, is therefore already required for the
  positive control and costs nothing extra here.

Together that sweeps D/V from 0.100 s down to 0.008 s, a factor of 12.5, with the footstep
at the top end. Sweep it deliberately and plot against D/V rather than against V: the whole
mechanism predicts that hardening is controlled by that one group, so the crossover should
collapse on it. The range is well placed for both sands of interest — for 1.0 mm sand the
predicted drained/hardening crossover falls at 0.76 m/s, inside the range, and for 0.2 mm
sand the Q4 cavitation ceiling falls at 1.96 m/s, also inside it.

For other sphere diameters the matched speed is V = D/(0.1 s) and the matched drop height is
H = (10D)²/2g ≈ 5.1 D², so it grows as D²: 0.5 mm for a 10 mm sphere, 2.0 mm for 20 mm,
8.2 mm for 40 mm. The largest sphere is the easiest to release cleanly at the matched speed,
which is a reason to make the largest sphere the workhorse rather than the exception.

One honest limit on the match. Fixing D/V reproduces the *drainage* group, but the
cavitation demand goes as V·D, which is 0.1 m²/s for a footstep and 0.004 m²/s at the
matched lab point — 25× lower. Matching both groups at once forces D = 0.1 m and V = 1 m/s,
i.e. full scale, so a bench rig cannot do it. A real footstep on 0.2 mm sand is
cavitation-limited; the matched lab drop is not. Say so in any write-up, and treat the
cavitation group as the separate axis it is, which is Q4.

Never drop twice into the same spot, because the bed is compacted after one impact.
Re-prepare between drops.

Imaging comes in two tiers.

*Tier 1, final penetration depth δ against drop height H.* A ruler and a still camera.
Robust, cheap, and directly comparable to Zhao et al. (*PRE* 2015). Do this first. It will
work.

*Tier 2, z(t) and a_max from high-speed video.* This is the measurement that touches the
Darcy-Reynolds law. At the top of the corrected range, 2.5 m/s, a 240 fps camera advances
the sphere 10 mm per frame, which cannot resolve a 1 to 3 cm stopping event; 960 fps is
marginal; 5000 fps makes it easy, so use the lab camera if it is available.

The slow end of the sweep inverts the problem rather than worsening it. At the
footstep-matched 0.20 m/s the sphere advances 0.83 mm per frame at 240 fps and only 0.04 mm
per frame at 5000 fps, which for any realistic px/mm is sub-pixel and will bury a_max in
differentiation noise. So do not fix one frame rate for the whole sweep. Choose fps per
speed so the sphere advances a few pixels per frame,

fps ≈ V × (px/mm × 10³) / (2 to 5 px),

and record the fps actually used with every recording. Calibrate px/mm from a target in the
impact plane every session and record the value with the data.

Record every session: sphere mass and diameter, drop height, water temperature (η varies
with T and enters the law as η^0.4), φ₀, camera fps, px/mm, ambient pressure.

---

## 6. Positive control, which is not optional

Before touching beach sand, buy two or three monodisperse glass bead bins spanning roughly
100 to 600 μm and reproduce a_max ∝ d^−0.8 V^1.6 in water. Run the control over Strader's
own 1 to 2.5 m/s, which is the top half of the §5 range; the slow end down to 0.2 m/s is
outside their published window and is ours to establish, so a disagreement there is not a
rig failure.

Without this, a null or odd result on natural sand is uninterpretable. It could equally be
trapped air, a wall effect, a tracking error, or real physics, and there is no way to tell
them apart after the fact. The reference curve from Strader et al. is the acceptance
criterion. If the rig cannot reproduce it on beads, the rig is not ready.

---

## 7. Analysis

1. **Q1.** Fit a_max ∝ d^−0.8 with d taken successively as d₅₀, d₁₀, and the Kozeny-Carman
   effective diameter d₃₂ from §4.1. Compare collapse quality and plot residuals against
   C_u. A residual that trends with sorting is the result. Report the d₁₀-against-d₃₂
   contrast as the primary one; d₅₀ is the straw man, and d₃₂ tracks it too closely
   (1.09× in a_max at C_u = 2) to be the interesting comparison.
2. **Q2.** For each sand, locate the φ₀ at which behaviour flips from sink-through to abrupt
   arrest. Plot that boundary in (φ₀, d) and test whether one dimensionless group collapses
   it.
3. **Q3.** Plot φ_c against measured circularity across samples.
4. **Q4, a named aim.** Push drop height until a_max departs from V^1.6. Saturation is
   predicted when ηψVD/k reaches p_atm − p_vap ≈ 99 kPa. Analyse it as a **measurement of
   ψ**, not as a yes/no confirmation:

   - A departure at V_dep gives ψ = k(p_atm − p_vap)/(η V_dep D) directly, per sample and
     per φ₀. That is a measured ψ(φ₀, d), which is the quantity §2 currently has to assume
     and the quantity the 43× inconsistency in §1 Q4 demands.
   - No departure up to 2.5 m/s is also a result: it bounds ψ ≤ k(p_atm − p_vap)/(η·2.5·D),
     which at D = 20 mm is ψ ≤ 0.078 for 0.2 mm sand and ψ ≤ 0.044 for 0.15 mm. Either way
     the sweep returns a number, which is exactly why this is an aim and not a gamble.
   - Fine sand is where it is reachable. With ψ = 0.1 the ceiling falls inside the §5 range
     for d ≲ 0.23 mm and is out of reach above 0.5 mm; if the true ψ is nearer the 0.0055
     bound the ceiling moves up by ~18× and only the finest fraction will show it. **Run the
     finest samples first**, and treat coarse sand as the control that must *not* depart.
   - The confound is real and must be excluded, not waved at: a departure could equally be a
     change in penetration-depth regime. The available discriminator is the d² scaling — the
     cavitation ceiling must move as d² across samples, a depth-regime change need not — so
     Q4 needs at least three grain sizes to mean anything, and a single departure curve
     proves nothing. (A partial vacuum over the bed would move the ceiling and nothing else,
     which is the clean test, but that is hardware this protocol does not specify. Day-to-day
     barometric variation is only a few percent and is far too small to use.)
   - Speeds above 2.5 m/s are an extension, not part of the standard sweep. The §5 range
     covers d ≲ 0.23 mm at ψ = 0.1; reaching 0.3 mm needs 4.4 m/s and a 0.99 m drop. Add
     height only for samples that the standard sweep has already failed to break.
5. **Field closure.** For each station, is the field φ₀ on the solid or the fluid side of
   that sample's own boundary? If the walkability classification falls out of the lab phase
   diagram, the loop is closed. Separately, plot field W against h_c/tanβ on log axes.

State the confound explicitly in any write-up: slope and grain size co-vary on beaches
(Bascom 1951). Any apparent slope effect has to be tested at fixed d₅₀ before it can be
claimed as independent of grain size.

---

## 8. Sequencing

1. Build and validate the rig on glass beads (§6), over Strader's 1 to 2.5 m/s. Nothing else
   matters until this passes.
2. Extend the bead rig down to 0.2 m/s and sweep D/V across the full 0.100 to 0.008 s of §5.
   This is cheap, uses beads already bought, and is what puts the rig in the footstep regime;
   do it before any beach sand, because it is also where the slow-end imaging settings get
   worked out.
3. One beach, one transect, full field protocol (§3). Characterise (§4), grading curves
   retained.
4. Tier-1 penetration depths on all samples. Check they order sensibly with d₁₀.
5. Tier-2 high-speed on a subset spanning the sorting range. This answers Q1.
6. Q4 on the finest samples, over the same §5 range, as a ψ measurement. Fine samples first,
   at least three grain sizes, coarse as the must-not-depart control.
7. φ₀ sweeps for the phase boundary. This answers Q2 and Q3, and feeds ψ(φ₀) back into Q4.
8. Additional beaches only after the first transect has produced a clean Q1 answer.

---

## References

- Reynolds, O. (1885). *Phil. Mag.* **20**, 469. Original dilatancy paper, and it uses the
  beach observation.
- Bascom, W. N. (1951). *Trans. AGU* **32**, 866. Sand size against beach-face slope.
- Jerome, J. J. S., Vandenberghe, N. & Forterre, Y. (2016). *PRL* **117**, 098003.
- Strader, J., Causley, N., Dijksman, J. A. & Clark, A. H. (2022). *PRF* **7**, 054304.
  arXiv:2107.12584.
- Zhao, S.-C. et al. (2015). *PRE* **91**, 022202. arXiv:1403.8124.
- Iverson, R. M. (2005). *JGR Earth Surf.* **110**, F02015.
- Rondon, L., Pouliquen, O. & Aussillous, P. (2011). *Phys. Fluids* **23**, 073301.
- Wu, Y. et al. (2021). *Granular Matter* **23**, 88. Natural against pluviated beach sand,
  X-ray CT.
- *Soft matter physics of the ground beneath our feet* (2024). *Soft Matter* **20**, 5859.
  arXiv:2308.00162.
- Seguin, A. et al. (2008). arXiv:0805.4189. Confinement effects on granular penetration.
