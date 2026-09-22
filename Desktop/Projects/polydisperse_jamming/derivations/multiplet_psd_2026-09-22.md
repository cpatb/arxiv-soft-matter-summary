# Aggregation as a packing variable: multiplets, trapped fluid, and the sign of the jamming shift

**Date** 2026-09-22 · **Status** unrefereed at time of writing (referee pass logged in §9)
**Checks** `derivations/multiplet_psd_2026-09-22_check.py` → `..._check.txt`
**Numbers** `multiplet_levers.py` → `out/multiplet_levers.txt`
**Code** `psd.MultipletPSD`, family string `multiplet:f<D_F>:n<N1>/<N2>/...`
**Tests** `test_multiplet.py` (17 anchors)

---

## RESULT BOX

> For a suspension of one primary particle that has aggregated into multiplets with
> number fractions `n_k` and mass-fractal dimension `d_f` (write `gamma = 3/d_f`):
>
> | quantity | closed form | fitted? |
> |---|---|---|
> | number-mean diameter | `d10 = d1 · M_{γ/3}` | no |
> | width | `CV² = M_{2γ/3}/M_{γ/3}² − 1` | no |
> | Sauter ratio | `d32/d10 = M_γ/(M_{2γ/3} M_{γ/3})` | no |
> | **trapped-fluid factor** | **`Λ = φ_eff/φ_solid = M_γ/M_1 ≥ 1`** | **no** |
>
> with `M_s = ⟨k^s⟩` under the number fractions. `Λ = 1` **exactly** at `d_f = 3`.
>
> **Aggregation moves a suspension toward jamming iff `Λ > G`**, where
> `G = φ_J(multiplet PSD)/φ_J(mono)` is the polydispersity gain. Since `Λ → ∞` as
> `d_f → 0` while `G` is bounded, **a critical `d_f*` must exist**, and it is found at
>
> **`d_f* = 2.89 – 2.95` across `⟨k⟩ = 1.1 – 8`** — essentially independent of *how much*
> the suspension has aggregated.
>
> **The closure-free form, which is the one to quote:** this project's *widest measured*
> pack reaches only `G = 1.0911`. `Λ` alone passes that ceiling at `d_f = 2.27 – 2.87`,
> and at the diffusion-limited `d_f = 1.8` it is **1.22 – 4.11**. So for any realistic
> colloidal aggregate (`d_f ≈ 1.8 – 2.1`) **the trapped-fluid lever wins by a factor of
> several, and no width effect this project has ever measured can offset it.**
>
> **Bench consequence.** The sign of the aggregation effect is set by the aggregate's
> *openness*, not by the amount of aggregation — so **two suspensions at the same measured
> CV can sit on opposite sides of the jamming shift**, and a CV measurement alone cannot
> predict the rheology.

---

## 1. Frame

**Sought.** How the jamming state of a suspension moves when its primary particles
aggregate, expressed in variables an experiment can measure.

**Why this family and not another.** Every other PSD family in `psd.py` chooses a *shape*
for the size distribution and asks what packing does. A real suspension has no such
freedom: one synthesis makes one primary particle, and the width of its distribution is a
**readout of the aggregation state**. Salt, pH and shear history move aggregation, not CV.
So "raise CV" and "aggregate more" are different experiments, and this derivation is about
the second one. (`test_multiplet.py::test_multiplet_cv_is_not_free_to_reach_arbitrary_values`
makes that concrete: singlets-plus-doublets at `d_f = 3` cannot exceed `CV = 0.14`, so the
`CV = 0.5` the continuous families run at is *unreachable* by compact doublet formation.)

**Assumptions**, each one load-bearing and each one named:

- **A1 — effective sphere.** A `k`-mer is replaced by a single sphere of diameter `d_k`.
  Real multiplets are non-spherical; §7 states what this costs and how to measure it.
- **A2 — fractal scaling.** `d_k/d_1 = k^(1/d_f)`, one exponent for all `k`. This is the
  standard mass-fractal relation. It is exact in the two limits it interpolates between
  (`d_f = 3`, compact; `d_f = 1`, straight chain) and is an assumption between them.
- **A3 — no interpenetration.** Aggregates are treated as impenetrable at their effective
  diameter. This makes `Λ` an **upper bound**; see §7.
- **A4 — rigid multiplets.** Aggregates do not break up or restructure under load. The
  derivation is about a state, not about its history.
- **A5 — the gauge.** All `φ_J` are `mixture_model`'s channel: `p̂ = 1e-3`, `μ_p = 0.3`,
  `N = 1500`, monodisperse `β = 0.5947` (measured, not fitted). `G` is a ratio *within*
  that channel, which is why `G` and not `φ_J` is the quantity carried forward.

**Dimensionless groups.** `k` is a pure count; `d_f` and `γ = 3/d_f` are pure exponents;
`Λ`, `G`, `CV` and both moment ratios are dimensionless. The only dimensional input is
`d_1`, and it enters every diameter exactly once (checked with `pint`, §5).

---

## 2. Derivation, every step

Let the monomer have diameter `d_1` and solid volume `v_1 = (π/6) d_1³`. Let `n_k` be the
**number** fraction of `k`-mers, `Σ n_k = 1`, and write `M_s ≡ ⟨k^s⟩ = Σ_k n_k k^s`.

**(1)** A `k`-mer contains exactly `k` monomers, so its **solid** volume is

&nbsp;&nbsp;&nbsp;&nbsp;`v_k^solid = k v_1`&nbsp;&nbsp;— exact, by definition of a `k`-mer.

**(2)** By **A2**, its **effective** diameter is `d_k = k^(1/d_f) d_1 = k^(γ/3) d_1`, so
its **effective** (excluded) volume is

&nbsp;&nbsp;&nbsp;&nbsp;`v_k^eff = (π/6) d_k³ = k^γ v_1`.

The whole derivation turns on (1) and (2) differing: solid goes as `k`, excluded volume as
`k^γ`. They coincide iff `γ = 1`.

**(3) Number-mean diameter.** `d10 = Σ n_k d_k = d_1 Σ n_k k^(γ/3) = d_1 M_{γ/3}`.
Inverting gives the monomer size from a measured mean: `d_1 = d10 / M_{γ/3}`.

**(4) Width.** `⟨d²⟩ = d_1² M_{2γ/3}`, so

&nbsp;&nbsp;&nbsp;&nbsp;`CV² = ⟨d²⟩/⟨d⟩² − 1 = M_{2γ/3} / M_{γ/3}² − 1`.

**(5) Sauter ratio.** `⟨d³⟩ = d_1³ M_γ`, hence
`d32/d10 = (⟨d³⟩/⟨d²⟩)/⟨d⟩ = M_γ / (M_{2γ/3} M_{γ/3})`.

**(6) De Brouckere ratio.** `d43/d10 = M_{4γ/3} / (M_γ M_{γ/3})`.

**(7) The trapped-fluid factor.** Summing (1) and (2) over the population, the ratio of
total excluded volume to total solid volume is

&nbsp;&nbsp;&nbsp;&nbsp;**`Λ ≡ φ_eff/φ_solid = Σ n_k k^γ v_1 / Σ n_k k v_1 = M_γ / M_1`.**

**(8) The bound.** For `γ ≥ 1` and `k ≥ 1`, `k^γ ≥ k` **pointwise**, so `M_γ ≥ M_1` and
`Λ ≥ 1`, with equality iff `γ = 1` or all weight sits at `k = 1`. No convexity argument is
needed; the pointwise inequality is enough.

**(9) The two levers.** A suspension mixed to solid fraction `φ_s` presents the packing
with `φ_eff = Λ φ_s`. Its multiplets are polydisperse, and polydispersity raises the
jamming fraction: `φ_J = G · φ_J^mono` with `G ≥ 1`. Distance to jamming, relative to the
*same suspension unaggregated* (where `Λ = 1` and `G = 1`):

&nbsp;&nbsp;&nbsp;&nbsp;**`J/J₀ = (Λ φ_s / G φ_J^mono) / (φ_s / φ_J^mono) = Λ / G`.**

Aggregation **promotes** jamming iff `Λ > G`.

**(10) Existence of `d_f*`.** At `γ = 1`, `Λ = 1` exactly while `G > 1` for any
non-degenerate `n_k` (compact multiplets still differ in size, `d_k = k^{1/3} d_1`), so
`J/J₀ < 1`: **compact aggregation can only suppress jamming.** As `γ → ∞`, if any weight
`w` sits at some `K > 1` then `Λ ≥ w K^γ / M_1 → ∞`, while `G = φ_J/φ_J^mono < 1/φ_J^mono`
because `φ_J < 1`. `Λ` is continuous and increasing in `γ`. **So a crossing `γ*` exists.**

**(11) Why `d_f*` hardly moves — first order about `γ = 1`.** With
`L₁ ≡ ⟨k ln k⟩/⟨k⟩ = dΛ/dγ |_{γ=1}`:

&nbsp;&nbsp;&nbsp;&nbsp;`Λ ≈ 1 + (γ−1) L₁`, &nbsp;&nbsp; `G ≈ G₀ + (γ−1) G₀'`,

&nbsp;&nbsp;&nbsp;&nbsp;`γ* − 1 = (G₀ − 1)/(L₁ − G₀') ≈ (G₀ − 1)/L₁`.

Both `G₀ − 1` and `L₁` grow with aggregation, so their **ratio** is nearly invariant —
which is the explanation for the flat `d_f*` column. Checked numerically in §6: the
first-order `d_f*` tracks the exact one to 0.005 and lies **below** it throughout, the
sign expected from dropping `G₀' > 0`.

---

## 3. SymPy check — `C1`

Symbolic `K = 3` population, number fractions carried as free symbols, normalisation
imposed symbolically. Each closed form is compared against a direct symbolic sum.

```
   [PASS] C1a  d10 = d1 * M_{gamma/3}
   [PASS] C1b  CV^2 = M_{2g/3}/M_{g/3}^2 - 1
   [PASS] C1c  d32/d10 = M_g/(M_{2g/3} M_{g/3})
   [PASS] C1d  d43/d10 = M_{4g/3}/(M_g M_{g/3})
   [PASS] C1e  Lambda = M_gamma / M_1
   [PASS] C1f  Lambda == 1 identically at gamma = 1 (any n)
```

`C1f` is the one worth naming: `Λ = 1` at `γ = 1` is asserted **identically in the
symbols**, not for a particular composition, which is what step (7) claims.

---

## 4. Limiting cases — `C2`

| limit | expected | got |
|---|---|---|
| `γ = 1` (`d_f = 3`, compact) | `Λ = 1` exactly | `1.0000000000000000` |
| monomers only, any `d_f` | `CV = 0`, `Λ = 1` | `0.00e+00`, `1.0000000000000000` |
| all particles are `K`-mers | `CV = 0`, `Λ = K^{γ−1}` | `3.9045287771` vs `3.9045287771` |
| `d_f = 1` (straight chain) | `d_k = k d_1`, `d10/d_1 = ⟨k⟩` | `1.5` at `⟨k⟩ = 1.5` |
| monodisperse | `d32/d10 = d43/d10 = 1`, matching `psd.py` at `CV = 0` | PASS |

The `K`-mer case is the discriminating one: the inverse ordering would give `K^{1−γ}`, and
it is the check that fixes which of `M_γ`, `M_1` is the numerator.

---

## 5. Dimensions — `C3`

`pint`: `d_k` is a length (`1.4909 µm` from `d_1 = 0.8 µm`, `k = 3`, `d_f = 1.7647`),
`v_k` a volume, `Λ` dimensionless (`2.1577`). By hand: `k` is a pure count and `γ` a pure
exponent, so `k^{γ/3}` is dimensionless and every `d_k` carries `d_1`'s metres exactly
once; `M_s` is a sum of (dimensionless × dimensionless).

---

## 6. Numeric spot-check against an independent route — `C4`, `C6`

The **anchor that does not pass through the algebra**: 2 × 10⁶ particles are drawn one at a
time, each converted to a diameter and a solid volume individually, and reduced with plain
means. No moment is ever formed. Four populations × five quantities:

| case | `γ` | worst relative error |
|---|---|---|
| DLCA-like, mostly singlets | 1.667 | 9.2 × 10⁻⁴ |
| RLCA-like, heavier tail | 1.429 | 2.1 × 10⁻⁴ |
| compact doublets only | 1.000 | 2.8 × 10⁻⁶ (`Λ` exactly 0) |
| strong aggregation | 2.000 | 2.4 × 10⁻⁴ |

Worst over all twenty comparisons: **9.16 × 10⁻⁴**, consistent with Monte-Carlo error at
this sample size. `Λ` in the compact case agrees to **0.00e+00** — it is an identity, and
the sampler reproduces it exactly rather than approximately.

`C5a`: `Λ ≥ 1` on 20 000 random `(n, K, γ ∈ [1,3])` draws, **0 violations**.

---

## 7. What this cannot decide, and what would decide it

1. **A1/A3 — the effective-sphere and no-interpenetration assumptions make `Λ` an UPPER
   BOUND.** Two open fractal aggregates can interpenetrate, which reduces the excluded
   volume below `(π/6) d_k³`. ⚠ This matters *most* where `Λ` is largest (low `d_f`).
   **But it matters least where the crossing is**: `d_f* ≈ 2.9` is a nearly-compact
   aggregate, where interpenetration is small. So the *location of the crossing* sits in
   the most defensible part of the model, while the *magnitude* at `d_f = 1.8` is an upper
   bound rather than an estimate. Stated as a bound throughout, never as a value.
2. **A1 — shape.** `G` here is the **size** consequence of aggregation only. A doublet is a
   dumbbell, and this project's shape axis measures separately that non-sphericity changes
   packing. `PREREG_MULTIPLET_2026-09-22.md` registers the DEM cells that separate them:
   the *same* `n_k` run as effective spheres and as rigid clumps. Until those run, the
   shape term is **unmeasured, not small**.
3. **A2 — one exponent for all `k`.** For small `k` the fractal relation is an
   interpolation, not a law: a physical doublet has a *specific* geometry, not a `d_f`.
   ⚠ The reported `d_f ≈ 1.8` (DLCA) and `≈ 2.1` (RLCA) are **[UNVERIFIED]** — not
   resolved to a DOI in this repository. This is why `d_f` is **swept** and never assumed,
   and why the result is stated as a *threshold in `d_f`* rather than as a value at one.
4. **`G` is the only fitted half.** §7 of `multiplet_levers.py` therefore re-runs the
   whole argument with `G` replaced by its **measured ceiling** (`1.0911`, the widest pack
   this project has ever made), and the conclusion survives without the closure.
5. **Not addressed here: the onset stress `σ*`.** Whether aggregation moves `σ*` as well as
   `φ_J` is a separate question with its own gauge problem, and `SIZE_THRESHOLD_2026-09-08.md`
   and `onset_size_degeneracy_2026-09-19.md` already own it. **No claim about `σ*` is made
   in this file.** It is taken up in `MULTIPLET_WC_2026-09-22.md`, which reads those two
   first.

---

## 8. Open questions

- Does the DEM reproduce `G` for a *multiplet* PSD? The closure was calibrated on binaries
  and transfers to three continuous families; multiplet PSDs are a fourth shape and no
  point of that kind is in the fit.
- Is `Λ` measurable in a DEM at all? It is a bookkeeping ratio between two definitions of
  volume, so a sphere pack cannot test it — only a **clump** pack, where the union volume
  and the enclosing-sphere volume genuinely differ, can. ⚠ `utils.porosity()` in Yade
  double-counts clump overlap, which this project has already recorded.
- At what `⟨k⟩` does the singlet/doublet/triplet description stop being adequate? §4 of
  `multiplet_levers.py` shows the truncation is harmless only once `K_max` is large: at
  `⟨k⟩ = 5, d_f = 1.8` the first cutoff tried (`K_max = 24`) was wrong by 3 %.

---

## 9. Referee

Pass 1: see §10 below (appended after the referee returns).
