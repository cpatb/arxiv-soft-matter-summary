# Who rattles: the rattler fraction from the contact-pair geometry

2026-09-24. Checks: `derivations/rattler_closure_check.py` → `out/rattler_closure.txt` (every number below
is printed there). Companion to `sigma_star_geometry_2026-09-24.md`, whose eq. 13 leaves the rattler
fraction f_r as the only non-geometric input to σ\*(CV).

> **Result.** Give each grain a Poisson number of contacts whose mean follows the pair weight that
> reproduces the contact statistics (closure C, w ∝ aᵢaⱼ(aᵢ+aⱼ)):
>
> **λ(a) = Z · a(a⟨a⟩ + ⟨a²⟩)/(2⟨a⟩⟨a²⟩)**, f_r = ⟨P_Poisson(k < 2 | λ(a))⟩_PSD.
>
> With the measured mean contact number Z and the PSD as the only inputs (no fitted parameter), this
> predicts the rattler fraction of 52 sheared packs (f_r 0.009–0.68) at rms 0.022; a size-blind null
> (λ = Z for every grain) misses at rms 0.100 and under-predicts CV 0.5 by 0.14–0.23.
> **Small grains rattle because a contact's probability scales with the partner sizes, not because they
> fit in holes.** Driven by the backbone count Z_nr alone (the length-free isostatic quantity), the same
> model recovers ~¾ of the CV-0.5 rattler excess (misses by −0.12 on average, null −0.36).
> OUT OF SAMPLE on 702 isotropic packs (three PSD families, μ_p 0–3): rms 0.041 vs null 0.128 (§4).

## 1. Frame
- Data: `phase9_actshape_*` healthy records (sheared, constant pressure; `rattler_frac_ss`, `z_ss`,
  `z_nr_ss`), and `phase2*` records (isotropic compression; `rattler_frac`, `z_all`) for §4.
- Rattler (the project's primary definition, `measure.RATTLER_MIN_CONTACTS = 2`): a grain with fewer than
  2 load-bearing contacts, removed ITERATIVELY. The model below is NOT iterated, so it should
  under-predict f_r (iteration can only add rattlers) — the sign of every residual is a check.
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
- **Frictionless: the backbone densifies** (lognormal Δφ_B +0.041 of +0.063, 2/3 of the rise; the CV-0.5
  cell is n = 1 converged seed). Removing that would need ψ ≈ 0.108 at CV 0.5, twice the model's 0.049,
  while the model's f_r miss there is −0.075 of 0.35 — so the bias shrinks the frictionless Δφ_B but does
  not plausibly remove it. The uniform and bidisperse frictionless rows (Δφ_B +0.013, +0.012 at CV 0.3–0.4)
  say the same with more seeds.
- This is conditional on the model's ψ; a per-particle dump (id, radius, contacts) would measure ψ directly.

## 5. What it means for σ\*
With eq. 13 of the σ\* derivation, σ\*/κ ≈ Zφ/(8π⟨a⟩)·C/q, and Z is now a function of (Z_nr, PSD): the
whole CV-dependence of σ\* at fixed ⟨a⟩ reduces to the PSD, φ, the isostatic Z_nr (length-free) and two
O(5 %) tail corrections — up to the ¼ of the rattler excess closure (ii) does not capture.

## 6. Open
- Assumption P is the weak link; the per-particle contact-count distribution at fixed size is not in
  any dump. A per-particle dump (id, radius, k) would test eq. 3 directly rather than through f_r.
- Iterated vs non-iterated rattlers: the model has no cascade.
