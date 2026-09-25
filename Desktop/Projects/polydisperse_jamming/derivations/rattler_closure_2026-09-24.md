# Who rattles: the rattler fraction from the contact-pair geometry

2026-09-24. Checks: `derivations/rattler_closure_check.py` → `out/rattler_closure.txt` (every number below
is printed there). Companion to `sigma_star_geometry_2026-09-24.md`, whose eq. 13 leaves the rattler
fraction f_r as the only non-geometric input to σ\*(CV).

> **Result (revised after referee pass 1 and a per-grain test, same night).** Give each grain a Poisson
> number of contacts with a size-dependent mean normalised to the measured Z. **The rattler fraction is
> predicted best, with no fitted parameter, by a mean proportional to the grain's SURFACE AREA,
> λ(a) = Z a²/⟨a²⟩**: out of sample on 702 isotropic packs (three PSD families, μ_p 0–3) rms 0.016,
> against 0.041 for the pair-closure form below and 0.128 for a size-blind null. A power sweep
> λ ∝ aᵖ has a sharp optimum at p = 2 (rms 0.114 / 0.079 / 0.016 / 0.092 / 0.158 at p = 0.5 / 1 / 2 / 3 / 4; `derivations/rattler_shape_compare.py`).
> A direct per-grain dump (§3b) agrees on the rattler probability P(rattler | a) (a²: rms 0.017 and 0.035
> at CV 0.2 and 0.5); for the MEAN ⟨k | a⟩ the rms ranking flips by cell (closure C better at CV 0.2, a² at
> CV 0.5), so the mean degree law is not a clean power.
> **Small grains rattle because a grain's contact count scales with its surface area, not because they fit
> in holes.** In a pair-resolved dump (one cell: CV 0.5, μ_p 0.3, isotropic), who touches whom carries no
> further size structure: ⟨aᵢaⱼ⟩_c = 0.532 against 0.536 for a configuration model with the MEASURED degrees,
> (Σkᵢaᵢ/Σkᵢ)², size assortativity r = −0.04 (§3b); the sheared dry-DEM and lubricated LF_DEM contacts agree
> to ≤ 4 % (σ\* derivation §5b; §3b below). ⚠ An earlier line here inferred a DISASSORTATIVE network by comparing
> sheared-pack pair statistics with isotropic-pack degrees; the same-cell test refutes it.
>
> First version (kept for the record): the closure-C form λ(a) = Z·a(a⟨a⟩+⟨a²⟩)/(2⟨a⟩⟨a²⟩) — the degree law
> of the pair weight that fits ⟨aᵢaⱼ⟩_c — gives rms 0.022 on the 52 sheared packs (null 0.100) and 0.041 out
> of sample. It was claimed as "closure C specifically"; the referee's sweep shows the OOS test does not
> pin C (on the sheared set C narrowly beats a², 0.022 vs 0.027; out of sample a² wins 2.5×).

## 1. Frame
- Data: `phase9_actshape_*` healthy records (sheared, constant pressure; `rattler_frac_ss`, `z_ss`,
  `z_nr_ss`), and `phase2*` records (isotropic compression; `rattler_frac`, `z_all`) for §4.
- Rattler (the project's primary definition, `measure.RATTLER_MIN_CONTACTS = 2`): a grain with fewer than
  2 load-bearing contacts, removed ITERATIVELY. The model below is NOT iterated, so iteration pushes it to
  under-predict f_r; the Poisson assumption (wider than the real, sub-Poissonian count distribution) pushes
  the other way. The CV-0 rows show the second bias alone (+0.004 in all 6 groups, referee m4).
- Assumption P: given its radius, a grain's contact count is Poisson. This is not derived; it is tested
  (it cannot be exact — contact counts in a jammed packing are narrower than Poisson at fixed Z — and the
  CV-0 rows measure how wrong it is when size plays no role).

## 2. Derivation
1. Under closure C the number of contacts between size classes a and b is ∝ p(a)p(b)·ab(a+b).
2. The expected contact count of one grain of radius a is then ∝ ∫p(b)·ab(a+b)db = a(a m₁ + m₂) (m_n =
   ⟨aⁿ⟩ over the NUMBER PSD).
3. Normalising so the PSD average is Z: ⟨a(a m₁ + m₂)⟩ = m₂m₁ + m₁m₂ = 2m₁m₂, so
   λ(a) = Z a(a m₁ + m₂)/(2 m₁ m₂) (eq. 3). Monodisperse: λ = Z.
4. P(rattler | a) = P(k ≤ 1) = e^{−λ}(1+λ); f_r = ⟨e^{−λ(a)}(1+λ(a))⟩ (eq. 4).
5. Backbone closure (ii): Z_nr = E[k; k ≥ 2]/P(k ≥ 2) with E[k; k≥2] = λ − λe^{−λ}; solve for the Z that
   reproduces the measured Z_nr, then eq. 4.

Dimensions: λ is a count; a(a m₁+m₂)/(m₁m₂) is L²/L² — dimensionless. Limits: CV → 0 gives the
size-blind Poisson value e^{−Z}(1+Z); a → 0 gives λ → 0 and P(rattler) → 1 (dust always rattles);
a → ∞ gives λ ∝ a², so large grains never rattle.

## 3. Results on the sheared packs (the set closure C was found on)
52 healthy `actshape` records, f_r (iterated, k < 2) from 0.009 to 0.682.

| model | input | mean error | rms | max abs |
|---|---|---|---|---|
| closure C + Poisson | measured Z | −0.009 | **0.022** | 0.070 |
| size-blind Poisson | measured Z | −0.062 | 0.100 | 0.231 |
| closure C + Poisson | measured Z_nr only | −0.038 | 0.067 | 0.157 |
| size-blind Poisson | measured Z_nr only | −0.125 | 0.203 | 0.440 |

From Z_nr alone, by CV: mean error +0.007 (CV 0), −0.012 (CV 0.2), −0.117 (CV 0.5); the null gives +0.007,
−0.050, −0.359. So closure C accounts for ~¾ of the CV-0.5 rattler excess that a size-blind count misses.
The predicted Z from Z_nr is 9 % high on average (rms 15 %).
⚠ This set is the one closure C was FOUND on (for ⟨aᵢaⱼ⟩_c, a different observable), so §4 is the real test.

## 3b. Direct per-grain test (`zbysize_dump.py` → `out/zbysize/*.npz`, `derivations/zbysize_check.py` → `out/zbysize_check.txt`)
Two isotropic cells rebuilt with the production `build_and_compress` (lognormal, μ_p 0.3, seed 1, p̂ 3e-5,
cap 1e6), dumping per grain (diameter, contact count, iterated rattler flag). **Wiring:** the CV-0.5 cell
reproduces the existing record's φ = 0.6380369863320716, z_all 2.216 and f_r 0.5173 exactly. Eight
quantile bins in d; each λ shape normalised to the measured Z:

| cell | ⟨k|d⟩ rms: C / a¹ / a² | P(rattler|d) rms: C / a¹ / a² | f_r pred C / a² (meas) |
|---|---|---|---|
| CV 0.2 (Z 3.89) | **0.18** / 0.48 / 0.38 | 0.058 / 0.096 / **0.017** | 0.137 / 0.163 (0.169) |
| CV 0.5 (Z 2.22) | 0.52 / 1.03 / **0.40** | 0.088 / 0.167 / **0.035** | 0.472 / 0.534 (0.517) |

**Pairs (CV-0.5 cell, rebuilt with pairs; φ again bit-identical to the record).** 1662 contacts:
⟨aᵢaⱼ⟩_c = 0.5319; configuration model with the measured degrees (Σkᵢaᵢ/Σkᵢ)² = 0.5360 (−0.8 %);
closure C m₃/m₁ = 0.4612 (+15 % below measured); uncorrelated a²b², (m₃/m₂)² = 0.5714; Newman size
assortativity r = −0.041. So in this isotropic pack the pair statistics follow from the degree law alone.
(For the SHEARED packs, the σ\* derivation §5b runs the id-free version of this test on load-bearing
contact ends: 0.973–0.982 at CV 0.5 against a finite-size null of 0.9999 — a small disassortative tilt. The
degree law k(a) under shear still needs a dump with grain ids.)

**The same test on LUBRICATED sheared suspensions (post hoc; `lfdem/decomp_lfsig.py` → `out/decomp_lfsig.txt`,
12 LF_DEM cells of PREREG_LFSIG, 1–3 int_ snapshots each in strain [1.5, 4], contacts = touching surfaces).**
⟨aᵢaⱼ⟩_c / configuration model = 0.999–1.003 at CV 0.2 (6 cells) and 0.960–0.999 at CV 0.5 (6 cells), with
size assortativity r between −0.15 and +0.09: the configuration-model statement holds to ≤ 4 %, with a
slight disassortative tilt at CV 0.5 and high stress. Closure C is within 2 % at CV 0.2 and 3–7 % HIGH at
CV 0.5, σ ≥ 30 σ₀. Few snapshots per cell; not a registered test.

The smallest CV-0.5 bin (d/⟨d⟩ 0.42) has ⟨k⟩ = 0.17 — below even a² (0.32): small grains are close to
contact-free, a threshold rather than a smooth power. The mean-degree law is not a clean power; the
rattler probability is, to within these errors, the a² one.

## 4. Out of sample: isotropic packs, three PSD families
702 converged `phase2` records (isotropic compression, a different protocol), lognormal / uniform /
bidisperse, μ_p 0–3, all p̂ rungs, group-mean f_r 0.012–0.563. Closure (i) from each record's own z_all:

- **closure C: rms 0.041 (mean −0.028); size-blind null: rms 0.128 (mean −0.094)** — three times better,
  with no parameter adjusted.
- Every group's mean residual is ≤ +0.005, and the misses grow with CV: worst −0.096 (bidisperse μ 0,
  CV 0.5), −0.075 (lognormal μ 0, CV 0.5), −0.063 (uniform μ 0, CV 0.5). That is the sign a non-iterated
  model must have. The misses SHRINK with friction (lognormal CV 0.5: −0.075, −0.047, −0.022, −0.008,
  −0.003 at μ 0, 0.3, 0.6, 1, 3).
- Bidisperse (two sizes only, where a Poisson count per size is least plausible) is the worst family and
  still ~3× better than the null.
- **Shape sweep (referee pass 1, reproduced in `derivations/rattler_shape_compare.py` →
  `out/rattler_shape_compare.txt`):** λ ∝ a² gives rms 0.016 (mean +0.009), best of every form tried.
- **Attrition (referee M2; `derivations/rattler_attrition_check.py` → `out/rattler_attrition.txt`).**
  Convergence falls with CV (97.6 % at CV 0 → 75.0 % at CV 0.5, deduplicated). The bias is FAMILY-SPECIFIC:
  uniform CV 0.5 keeps 35 of 60 and its non-converged cells rattle more (f_r 0.599 vs 0.411), so the uniform
  CV-0.5 f_r here is survivor-biased LOW; lognormal (82/108; 0.499 vs 0.491) and bidisperse (45/48; 0.452 vs
  0.460) show no such difference. (The first version said "41 of 84", carrying an arithmetic slip from the
  referee's pass-1 report; corrected by pass 2 and now printed by the script.) The model
  comparison is not: on the 97 excluded records alone, closure C still beats the null by > 3× (rms 0.047 vs
  0.157, referee's computation).

## 4b. The rattler VOLUME share, and backbone invariance (closes an UNDETERMINED item, model-conditionally)
`geom_rattler_decomposition_2026-09-19.md` could not decide whether polydispersity only adds rattlers to an
unchanged backbone (φ = φ_B/(1−ψ), φ_B fixed) because no record carries the rattler VOLUME share ψ. The
count-validated model supplies it: ψ = ⟨a³ P(k<2|λ(a))⟩/⟨a³⟩. At the p̂ 3e-5 rung (seed means):

| family, μ_p | Δφ (CV 0 → 0.5) | Δφ_B | f_r → | ψ_model → |
|---|---|---|---|---|
| lognormal, 0.3 | +0.062 | +0.007 | 0.08 → 0.52 | 0.07 → 0.15 |
| lognormal, 1 | +0.053 | −0.019 | 0.14 → 0.60 | 0.13 → 0.24 |
| lognormal, 3 | +0.052 | −0.022 | 0.14 → 0.60 | 0.14 → 0.25 |
| uniform, 0.3 | +0.026 | −0.009 | 0.08 → 0.47 | 0.07 → 0.13 |
| bidisperse, 0.3 | +0.012 | −0.023 | 0.08 → 0.50 | 0.07 → 0.13 |
| lognormal, 0 | +0.063 | +0.041 | 0.02 → 0.35 | 0.02 → 0.05 |

Reading, with the model's known bias stated first: the model UNDER-predicts f_r (by up to ~0.1 at CV 0.5,
μ 0), so ψ_model is low and φ_B is over-estimated at high CV; correcting it would push every Δφ_B DOWN.
- **Frictional (μ_p ≥ 0.3): the backbone does not densify.** |Δφ_B| ≤ 0.023 against Δφ up to 0.062; the
  CV-driven rise in φ is carried by rattler volume sitting in the backbone's pores. (With the bias, the
  frictional backbone if anything LOSES solid.) Rattlers are ~half the grains but only 13–25 % of the solid.
- **Frictionless: the backbone densifies** (lognormal Δφ_B +0.041 of +0.063 with closure-C ψ; +0.030, i.e.
  48 % of the rise, with the better a² ψ — the fraction is closure-dependent; the CV-0.5 cell is n = 1
  converged seed). Removing it would need ψ ≈ 0.108 at CV 0.5, twice the C model's 0.049. The uniform and
  bidisperse frictionless rows densify at CV 0.3–0.4 (Δφ_B +0.013, +0.012) but the ladder is not
  monotone: **bidisperse μ 0 CV 0.5 (n = 3) gives Δφ_B = −0.002**, i.e. no densification (referee m3).
- With the a² ψ the frictional conclusion holds and tightens: |Δφ_B| ≤ 0.015 on every frictional row
  (`derivations/backbone_psi_compare.py` → `out/backbone_psi_compare.txt`, which also prints the C column;
  it reproduces the referee's recomputation). With a² the frictionless uniform/bidisperse densification
  shrinks to +0.002 … +0.006.
- This is conditional on the model's ψ; a per-particle dump (id, radius, contacts) would measure ψ directly.

## 5. What it means for σ\* (written for closure C; the Z_nr → Z step has not been redone with a²)
With eq. 13 of the σ\* derivation, σ\*/κ ≈ Zφ/(8π⟨a⟩)·C/q, and Z is now a function of (Z_nr, PSD): the
whole CV-dependence of σ\* at fixed ⟨a⟩ reduces to the PSD, φ, the isostatic Z_nr (length-free) and two
O(5 %) tail corrections — up to the ¼ of the rattler excess closure (ii) does not capture.

## 6. Open
- Assumption P is the weak link; the per-particle contact-count distribution at fixed size is not in
  any dump. A per-particle dump (id, radius, k) would test eq. 3 directly rather than through f_r.
- Iterated vs non-iterated rattlers: the model has no cascade.
