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
> **`d_f* ≈ 2.89 – 2.95` over `⟨k⟩ = 1.1 – 8`**, stable to ±0.01 across the two refit
> calibrations and ±0.1 to the textbook one (§7 item 6).
>
> 🔴 **THE FLATNESS IS NOT ESTABLISHED AT THAT PRECISION, and an earlier version of this
> box claimed it was.** The plateau is only **0.053** wide, and one `PHI_SCALE` of
> protocol noise moves the *root* by
>
> | `⟨k⟩` | 1.1 | 1.2 | 1.5 | 2.0 | 3.0 | 5.0 | 8.0 |
> |---|---|---|---|---|---|---|---|
> | `Δ(d_f*)` | 0.081 | 0.042 | 0.019 | 0.011 | 0.007 | 0.005 | 0.004 |
> | as % of the plateau | **152 %** | **80 %** | **36 %** | 22 % | 14 % | 10 % | 8 % |
>
> So at low `⟨k⟩` the individual roots are less well pinned than the whole spread they are
> being compared across, and their agreement with the high-`⟨k⟩` rows **is not evidence of
> flatness**. ⚠ The first version of this box filtered rows on
> `|φ_J(d_f*) − β| ≥ PHI_SCALE`, which asks whether `G ≠ 1` at the crossing — **the wrong
> question**; the right one is the sensitivity of the root, above.
>
> **What survives without a precision argument is ORDINAL and is the claim to carry:**
> every row's `d_f*` lies near **2.9** — far below the compact limit `d_f = 3` and far
> above the reported colloidal values 1.8–2.1 — on all three calibrations. That is enough
> for the conclusion, because the conclusion is about which *side* of the crossing a real
> aggregate sits on.
>
> **The closure-free form, which is the one to quote:** this project's *widest measured*
> pack reaches only `G = 1.0911`, so the polydispersity lever is capped near **+9 %,
> measured**. Asking where the closed-form `Λ` alone exceeds that ceiling needs no fit:
>
> | `d_f` | smallest `⟨k⟩` at which `Λ` reaches `G_max` |
> |---|---|
> | 1.8 | 1.08 |
> | 2.1 | 1.14 |
> | 2.5 | 1.37 |
>
> **Above those thresholds, no width effect this project has ever measured can offset the
> trapped-fluid lever** — and `Λ` is then **at most** 1.22 – 4.11 over `⟨k⟩ = 1.2 – 5` at
> `d_f = 1.8` (an upper bound; see §7 item 1). **Below them the closure-free test is
> uninformative**, and the fitted closure still gives the same *sign* but by a margin of
> a per cent, not a factor of several.
>
> ⚠ An earlier version of this box quoted `Λ = 1.22 – 4.11` as a *bound over any
> realistic colloidal aggregate*. It is the range over the five tabulated `⟨k⟩` values,
> and a lightly aggregated suspension is equally realistic. Corrected 2026-09-22 after a
> referee pass.
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
which is the explanation for the flat `d_f*` column.

**The sign.** `G₀' > 0`, so dropping it makes the denominator *larger*, hence
`γ*_approx < γ*_exact`; and `d_f = 3/γ`, so **`d_f*_approx > d_f*_exact`** — the
first-order estimate sits **closer to `d_f = 3`, i.e. ABOVE**. Measured (printed by
`multiplet_levers.py` §6):

| `⟨k⟩` | 1.1 | 1.2 | 1.5 | 2.0 | 3.0 | 5.0 | 8.0 |
|---|---|---|---|---|---|---|---|
| first-order − exact | +0.0034 | +0.0038 | +0.0049 | **+0.0058** | +0.0011 | −0.0003 | −0.0003 |

**Above on 5 of 7 rows, worst deviation 0.0058.**

> 🔴 **CORRECTED 2026-09-22 after a referee pass.** This paragraph previously read *"the
> first-order `d_f*` tracks the exact one to 0.005 and lies **below** it throughout"*.
> Both halves are wrong and **the table this very file prints refutes them**: the
> estimate is above on 5 of 7 rows and the worst deviation is 0.0058. The *mechanism*
> sentence that followed ("it puts the crossing too close to `d_f = 3`") was right —
> "too close to 3" means *larger*, i.e. above — so the prose contradicted itself as well
> as its own artifact. The two rows that land marginally below (−0.0003) are second
> order and are reported rather than absorbed into the argument.

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
   in this file.** It is taken up in **`multiplet_wc.py` §4** (a module, not a document —
   an earlier version of this line pointed at a `.md` that does not exist), which states
   the two candidate conventions and refuses to choose between them.

6. **`G` is the only fitted half, and the fit's FORM is on record as too rigid.** Three
   calibrations are now carried side by side in `multiplet_levers.py` §2: textbook
   (1.020, 1.500), the ratio-2-and-3 refit (0.911, 2.122), and the **wider** ratio-1.25-to-3
   refit (0.961, 2.214), which `out/mixture_closure.txt` E2 records as the **better**
   held-out transfer (family `χ² = 34.9` against 63.0) and which the first version of this
   work did not use. `d_f*` moves by **at most ≈ 0.01** between the two refits, so no
   conclusion turns on the choice — but the same source file's E1 verdict is that the
   two-exponent **form** is too rigid (joint `χ²` 3.55 × the sum of per-ratio fits), so
   **neither pair is "the" answer** and `G` should be read as an O(10 %) quantity with a
   form error, not as a calibrated constant.

7. **A4 (rigid, non-restructuring multiplets) is an assumption with NO evidence behind it
   here, and it is the one most likely to fail in the regime of interest.** Shear and
   compression are known to restructure open flocs toward compactness, i.e. to raise
   `d_f` *during* the very loading this model is about — and this project's own
   `dilatancy_sign_mechanism` work established that the packing response is
   history-dependent rather than a function of `φ`. Restructuring would move a suspension
   **up** through `d_f*` as it is loaded, which is a *time-dependent* crossing this
   static model cannot represent at all. ⚠ This item was **absent** from the first version
   of this list, which named A4 in §1 and then never returned to it.

8. **The effective-diameter convention is a free O(1) factor that is NOT swept.** A2 fixes
   `d_k/d_1 = k^{1/d_f}`, normalised so that `k = 1` returns the monomer exactly. But
   *which* diameter of an aggregate sets its excluded volume is a convention — enclosing
   sphere, radius-of-gyration-equivalent, or hydrodynamic radius — and the aggregation
   literature writes the mass-fractal relation with a structural prefactor `k_f`
   (`N = k_f (R_g/a)^{d_f}`) that is commonly away from 1. **[UNVERIFIED — not resolved
   to a DOI in this repository.]** Any constant prefactor `c` on `d_k` multiplies every
   effective volume by `c³` and therefore multiplies `Λ` uniformly, which shifts `d_f*`
   by an O(1)-factor's worth. **This is a real degree of freedom the `d_f` sweep does not
   explore**, it was missing from the first version of this list, and it is the single
   largest un-bounded uncertainty in the *location* of `d_f*` — though not in the
   *existence* of a crossing, which needs only that `Λ` diverge.

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

### Pass 1 — 2026-09-22. **One FATAL, six Material, three Minor. All acted on.**

The referee recomputed every number independently (its own SymPy re-derivation, direct
sweeps, `brentq` cross-checks) rather than reading the check output.

| # | finding | action |
|---|---|---|
| **F1 FATAL** | Step 11's sign claim (*"lies below throughout"*, *"tracks to 0.005"*) is **refuted by the table this file prints**: above on 5/7 rows, worst 0.0058. The prose also contradicted its own following sentence. | **Accepted in full.** §2 step 11 rewritten with the measured deviations and the corrected direction; `multiplet_levers.py` §6 prints the same table and the retraction. |
| **F2** | The flat-`d_f*` headline included `⟨k⟩ = 1.1`, where the predicted `φ_J` shift is 0.89 × `PHI_SCALE` — **below the project's own resolution bar**, which §3 states but never applied to §2's own table. | **Accepted.** §2 now prints the shift in `PHI_SCALE` units per row, flags unresolvable rows, and the quoted range is computed from resolvable rows only (`⟨k⟩ = 1.2 – 8`). |
| **M1** | *"wins by a factor of several for any realistic aggregate"* is the range over the five tabulated `⟨k⟩`, not a bound; the closure-free test is **uninformative** below `⟨k⟩ ≈ 1.15`. | **Accepted.** The `⟨k⟩` thresholds are now computed and printed per `d_f`, and the claim is scoped to above them. |
| **M2** | `Λ` is stated as a value ("it **is** 1.22–4.11", "by 87.0 %") while §7 claims it is used as a bound throughout. | **Accepted.** Result box now says "at most" and points at §7 item 1. |
| **M3** | `PHI_J_MONO = 0.5748` in the check script is a **different rung** (`p̂ = 3e-5`) from the `β = 0.5947` gauge every other number uses. | **Accepted.** Set to `0.5947`; a new check `C5b'` asserts it agrees with `mixture_model`, so the two gauges cannot drift apart again. (The referee's "matches none of the gauges" is because it checked one rung; 0.5748 *is* a real project number, at the other rung. The inconsistency stands either way.) |
| **M4** | The fractal structural prefactor `k_f` is an un-enumerated convention; it scales `Λ` as `k_f³`. | **Accepted.** New §7 item 8. |
| **M5** | A4 (rigid multiplets) is named in §1 and has no limitations item. | **Accepted.** New §7 item 7, with the restructuring mechanism and why it matters here. |
| **M6** | A better-transferring calibration (0.961, 2.214; `χ² = 34.9` vs 63.0) exists in the cited source and was not tried. | **Accepted.** Carried as a third column; `d_f*` moves ≤ 0.01. New §7 item 6 also records the source's own "form is too rigid" verdict. |
| **m1** | No test exercises `multiplet_levers.py`'s own logic. | **Accepted** — `test_multiplet_levers.py` added. |
| **m2** | The §4 convergence gate covered `⟨k⟩ ∈ {1.2…5}` but §2 and §6 also quote 1.1 and 8.0. | **Accepted.** The gate now covers the union (42 points, all converged). |
| **m3** | Dangling reference to a non-existent `MULTIPLET_WC_2026-09-22.md`. | **Accepted** — it is `multiplet_wc.py` §4. |

**What the referee confirmed independently** (criteria A–D, F): all five closed forms
exact in its own SymPy; the `Λ ≥ 1` pointwise proof correct; **no double counting** — `Λ`
and `G` are orthogonal computations on the same object in one volume convention; the
existence argument for `d_f*` logically sound regardless of F2/M3; and the closure-free
section genuinely closure-free (no `PA`/`PB` enters). It also verified the committed
artifact is not stale.

### Pass 2 — 2026-09-22. **No FATAL. Five Material, three Minor. All acted on.**

Pass 2 was asked to check pass 1's *corrections* as well as look for new defects, and it
found that **two of the four corrections were half-fixes with the wrong criterion.**

| # | finding | action |
|---|---|---|
| **P2-1** | **The resolvability filter added by F2 asks the wrong question.** It tests `\|φ_J(d_f*) − β\| ≥ PHI_SCALE` ("is `G ≠ 1` at the crossing"), not whether the *root's location* is resolvable — which needs `Δ(d_f*) = (PHI_SCALE/β)/\|Λ′ − G′\|` by the implicit function theorem. At `⟨k⟩ = 1.2`, kept as resolvable by the old filter, noise alone moves `d_f*` by **80 % of the entire plateau**. | **Accepted in full.** `df_star_sensitivity()` added and printed per row; **the flatness headline is demoted to an ORDINAL claim** in the result box. Both columns stay printed. A test anchors the estimator against a direct re-solve of the perturbed root. |
| **P2-2** | **The admissibility ceiling `Λ ≤ 1/φ` is a trivial space-filling bound.** Rigid non-interpenetrating spheres cannot exceed *jamming*. **3 of 5 rows breach the right ceiling, not 1.** | **Accepted.** `multiplet_wc.py` now uses the multiplet's **own** `φ_0^eff = G₀·PHI0 = 0.6861`. ⚠ Note the referee proposed a fixed RCP ≈ 0.64, which would **over-flag**: the effective pack is polydisperse and its own jamming point is higher. The bound has to be the multiplet's own. **This DEMOTES the bench prediction** — the `d_f = 1.8` rows are outside the model's domain. |
| **P2-3** | **The `hypot` justification for the friction comparison is wrong on this repo's own terms**: `packing.py` seeds the radius draw independent of `μ_p`, so the arms are paired (measured `r = +0.56 … +0.84`). | **Accepted.** Replaced by a per-seed paired difference with the floor propagated through all four records. The verdict **survives** (0 of 4 resolved) but the margins move: CV 0.3 `1.73 → 2.73`, CV 0.5 `1.43 → 2.56`. Both estimators printed. |
| **P2-4** | **Block C's pair selection correlates with the sign.** Of 80 candidates **29 (36 %) are positive**; of the 13 resolvable, 9 neg / 4 pos — and **every positive one sits at the smallest `d_f` gap**, which the rule sorts against. | **Accepted.** `PREREG_MULTIPLET` **Amendment 1.2** withdraws the "uniformly negative" framing; the four predictions and tolerances are unchanged and still blind. The full population is now printed by `predict_multiplet.py`. Block C's *discriminating* content is undamaged: a CV-only closure predicts zero for every pair, either sign. |
| **P2-5** | **Block B does NOT discriminate a CV-only closure**, though `PREREG` M3 says it does: ordering its four points by `CV` and by `φ` gives the **identical** order, so any increasing function of `CV` reproduces the turnover. | **Accepted.** **Amendment 1.1.** The per-point predictions stand; the interpretive claim is withdrawn. This is the repo's own "a discriminator that passes BOTH candidates has discriminated nothing" lesson — caught *before* the cells ran. |
| P2-6 | The step-11 mechanism is verified only as a *linear-vs-semi-linear* statement and is invoked as a general one; it predicts the wrong direction for 3 of 7 rows. | **Accepted** — §2 step 11 already reports the two negative rows as unexplained; the general wording is the part to read narrowly. |
| P2-7 | The covariance claim is **confirmed against the solver source** (`eps` is scale-invariant, and the clip, `EPS_SOFT` and `MP_EXPO` all act on `eps`) — but is unconditionally worded. `DIL_COEFF > 0` adds a term built from a **fixed** `PHI_REF = 0.49`, which is not covariant. | **Accepted.** `multiplet_wc.py` §1 now scopes the claim to `DIL_COEFF = 0` and says why. No printed number changes (`DIL_COEFF` defaults to 0). |
| P2-8 | The "same sign below the threshold" check is near-trivial, since all tabulated `d_f` sit well below `d_f*`. | **Noted**, no change; recorded as a weak check rather than presented as a strong one. |

**Confirmed clean by pass 2**: the step-11 sign correction reproduced exactly from a fresh
script; the `β = 0.5947` gauge fix consistent throughout; and no FATAL anywhere.

⚠ **Pass 3 has not run**, and pass 2 changed four files. The loop rule is "repeat until a
pass returns nothing"; neither pass has yet returned nothing.
