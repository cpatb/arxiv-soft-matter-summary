# PREREG — `actshape`: does the size distribution set the SHARPNESS of the Wyart–Cates activation?

Written 2026-09-22, third block, **before launch**. Analysis: `actshape.py` (pure functions,
written before any cell ran). Grader: `grade_actshape.py` (written before the cells land).

## 0. Why this campaign, from the question down

The question is which particle property decides whether an impacting drop **shear-jams**
or stays **inertial-dilatant**. In the reduced framework (`three_knobs.py`) the drop can
only jam inside the window `χ ≤ x < 1`, and it jams there once the frictional fraction
reaches

&nbsp;&nbsp;&nbsp;&nbsp;`f_c = (1 − x)/(1 − χ)`.

Whether the impact stress gets `f̂` to `f_c` depends on the **activation law**. The solver
writes it `f̂_ss = exp(−(σ*/σ)^δ)` — a Weibull survival function of shape `δ = DELTA`.
Two of five particle axes have been placed (friction and shape move `χ`; CV moves `x`).
**Size and potential enter only through the activation law, and nobody has measured it.**

The previous handoff proposed the isotropic `tail_slope` as a lead and named its gap:
isotropic packs, not sheared. Checked at session start, **before this prereg**: the
sheared records already on disk (`phase8_commonlad`, `phase9_rm`) carry `force_moments`,
and the flattening **survives under shear with the same sign but about 3× weaker**
(constant-pressure, rate-matched, `p̂` = 3e-4: `tail_slope` −1.31 / −1.24 / −1.00 at CV
0 / 0.2 / 0.5, against the isotropic −2.02 / −1.78 / −0.94). ⚠ The handoff's premise that
`shear.py` "needs only `force_moments` recorded" was false — it already records them.

**Why that summary is not enough, and what this campaign adds.** `f/⟨f⟩` is the wrong
variable. For every surface force in play (screened electrostatics, steric, DLVO) the
Derjaguin approximation makes the threshold scale with the pair's reduced radius,
`F_rep = κ R*`. So a contact is frictional iff `f/R* > κ`, and in a polydisperse pack
`f` and `R*` are correlated: big grains carry big forces and have big thresholds. The
summary on disk cannot form `f/R*`. This campaign dumps the per-contact arrays
`(f_n, f_t, R₁, R₂)` so it can.

## 1. Grid

Constant-pressure shear, the rate-matched protocol of `PREREG_CONSTP_RATEMATCH`
(`cp-shear-v1`, gain 0.15, `γ_max` = 3, `N` = 1000, `φ_start` = 0.55, lognormal), tag
`actshape`, dumps of 10 snapshots evenly spaced over the steady half of the strain:

&nbsp;&nbsp;&nbsp;&nbsp;CV ∈ {0, 0.2, 0.5} × `μ_p` ∈ {0, 0.3} × `p̂` ∈ {1e-4, 3e-4, 1e-3} × seeds {1,2,3} = **54 cells**

| CV | steps at `p̂` 3e-4, 1e-3 | steps at `p̂` 1e-4 |
|---|---|---|
| 0.0 | 102 000 | 204 000 |
| 0.2 | 200 000 | 400 000 |
| 0.5 | 480 000 | 960 000 |

`γ̇` is matched across CV exactly as in the rate-matched ladder, and **halved at
`p̂` = 1e-4** because `I ∝ γ̇/√p` would otherwise cross `I_MAX` = 1e-2 there.

**Why `μ_p` = 0 is in the grid.** In Wyart–Cates the contacts BELOW threshold are
lubricated, i.e. frictionless. So the force network that decides the **onset** (small `f̂`)
is the frictionless one, and the network near **full activation** is the frictional one.
The truth lies between the two arms; this brackets it instead of assuming one.

## 2. Estimators, fixed now

Per cell: pool the contacts with `f_n > 0` over the 10 snapshots. With `σ = p` (the cell's
own sampled `p̂` at each snapshot):

- `g = f_n/(R* σ)` (`norm="derjaguin"`, **primary**); `"constant"` (`f_n/σ`, size-blind) and
  `"area"` (`f_n/(R*² σ)`, a diagnostic bracket only) are reported beside it.
- **`δ₁₀/₉₀`** = 3.0843 / ln(Q₉₀/Q₁₀) (**primary shape**) and `δ₂₅/₇₅` = 1.5727 / ln(Q₇₅/Q₂₅).
- **`σ*/κ = 1/Q_{1−1/e}(g)`**: the onset, in the solver's own σ* convention.
- **`σ_jam/κ = 1/Q_{1−f_c}(g)`** at `f_c` ∈ {0.1, 0.25, 0.5, 0.75, 0.9}. **The primary output**,
  because it needs no σ* convention and no Weibull assumption.

Seeds are averaged; CV and `μ_p` contrasts are **seed-paired**, because `packing.py` seeds
the radius draw from the seed regardless of CV and `μ_p` (lessons.md, 2026-09-22). With
n = 3, `t_crit(df = 2) = 4.303`.

⚠ **Gauge.** `σ = p`. The `τ = μ_M p` gauge is printed beside every onset number and is not
the registered one.

⚠ **Length.** The DEM unit is `d₁₀ = 1`. "At similar particle size" has no unique meaning
across CV, so `σ*` comparisons are printed at matched `d₁₀`, `d₅₀` and `d₃₂`. Only `d₁₀` is
graded, because it is the one the data sets without a choice.

## 3. Gates

**W1 — IDENTITY. Primary; it gates everything.** The 9 cells at (`μ_p` 0.3, `p̂` 3e-4) are
the rate-matched ladder's own cells, re-run with the dump on. The dump is read-only, so
**`p_hat_series`, `z_series` and the `force` block must come back BIT-IDENTICAL** to the
`phase9_rm_*` records. `phi_series` is compared too; if it differs at ≤ 1e-12 while the
other three are identical, that is the known cross-era `porosity()` characterisation and is
**reported, not failed**.

**W2 — THE DUMP IS THE SAME NETWORK.** On every cell, `measure.force_moments` recomputed from
the LAST snapshot of the dump must equal the record's `force` block exactly.

**W3 — THE ANCHOR THAT CAN FAIL.** At CV = 0 every radius is equal, so `R*` is a constant and
**`δ` must be identical under all three norms** (to 1e-12), while at CV = 0.5 it must
**differ** between `"constant"` and `"derjaguin"` by more than 1e-3 (otherwise the radii are
not reaching the estimator).

**W4 — HEALTH.** `steady`, `servo_held`, `quasistatic` per cell. Failing cells are PRINTED
with the reason and excluded from graded contrasts; the counts are printed.

## 4. Blind predictions (committed before any cell ran)

- **H1 — sharpness moves with CV, and Derjaguin does not remove it.** At `μ_p` 0.3, `p̂` 3e-4:
  `δ₁₀/₉₀(derjaguin)` at CV 0.5 is LOWER than at CV 0, seed-paired `|t| ≥ 4.303`, with ratio
  **in [0.60, 0.90]**.
- **H1b — the ordering of the norms.** At CV 0.5: `δ(constant) < δ(derjaguin) < δ(area)`.
  Reason: across a big–small pair the force is limited by the small grain's area, so
  `f ∝ p r_s²` while `R* ≈ r_s`; `f/R* ∝ p r_s` still carries one power of size, `f/R*²`
  carries none.
- **H2 — the onset shifts at fixed `d₁₀`.** `σ*/κ` at CV 0.5 is LOWER than at CV 0
  (larger grains at fixed `d₁₀` carry larger `f/R*`). Direction at fixed `d₅₀`/`d₃₂`: not
  predicted.
- **H3 — friction.** `δ(μ_p = 0)` vs `δ(μ_p = 0.3)`: no direction predicted, reported.
- **H4 — is sharpness a MATERIAL property?** At each (CV, `μ_p`), the span of `δ₁₀/₉₀` over
  the three `p̂` rungs is **less than half** of `|δ(CV 0.5) − δ(CV 0)|` at the same `μ_p`.
  If this fails, the activation shape is pressure-dependent and "CV is a sharpness knob" is
  not quotable as a material statement.
- **H5 — is the solver's Weibull form adequate?** `δ₂₅/₇₅ > δ₁₀/₉₀` on the CV-0 cells (a
  force distribution peaked near its mean with an exponential tail is narrower in its middle
  than a Weibull that matches its outer quantiles). Reported with the ratio.
- **H6 — the composite.** At fixed `d₁₀` and `κ`, the curves `σ_jam/κ` vs `f_c` for CV 0 and
  CV 0.5 **cross** inside `f_c` ∈ [0.1, 0.9]: polydispersity lowers the stress needed to reach
  small `f_c` (deep in the window) and raises it near the window edge.

**Anchors outside the code under test.** Exponential `u` gives `δ = 1` from both quantile
pairs exactly; a Weibull sample of known `δ` is recovered (both in `test_actshape.py`).
The solver's calibrated presets carry `DELTA` = 1.385 (`si_fit`, withdrawn) and 1.317
(`cs_fit`); the DEM `δ` is printed beside them and **not graded against them** — a
different material and a different fitting route.

## 5. What is NOT claimed by this campaign

`κ` itself (the potential), the solver's `τ_up`/`τ_relax`, and anything about relaxation.
The DEM has no repulsive force law: `f̂` here is the fraction that WOULD be frictional
under (1) given the force network of a pack that is uniformly frictional (`μ_p` 0.3) or
uniformly frictionless (`μ_p` 0). A pack whose friction switches contact by contact is a
capability build and is not attempted.

---

## Amendment 1 — the finite-size control, written 2026-09-22 after grading, BEFORE its cells ran

**Why.** The graded campaign passed W1 (9/9 bit-identical), W2′/W2″ (54/54), W3, H1 (ratio 0.811,
t = −30.2), H2, H4 and H6. But **every CV = 0.5 cell carries `FINITE SIZE RISK: d_max/L_dense`
= 0.34–0.38 > 0.25**, including the `μ_p` = 0.3 cells the headline contrasts rest on, and all
six W4 exclusions are frictionless CV = 0.5 cells. The CV 0 → 0.2 contrasts are free of the
flag. So the CV 0.5 numbers are **not quotable** until this control returns.

**Grid.** The same protocol at **N = 4000** (d_max/L ≈ 0.36·(1/4)^{1/3} ≈ 0.23), CV = 0.5,
`μ_p` ∈ {0, 0.3}, `p̂` = 3e-4, seeds {1,2,3}, 480 000 steps, dumps on. 6 cells, tag
`actshape_fs`.

**Registered comparison.** Seeds at different N are different radius draws, so this is
**unpaired** (Welch). For `μ_p` = 0.3:
- **F1** `|δ₁₀/₉₀(N=4000) − δ₁₀/₉₀(N=1000)|` < 0.05 (the within-cell sd is 0.01–0.04): the
  N = 1000 CV 0.5 value is not a finite-size artefact.
- **F2** H1 recomputed with CV 0.5 at N = 4000 against CV 0 at N = 1000 (Welch t, df by
  Welch–Satterthwaite): `δ` falls, `|t| ≥ t_crit`, ratio in [0.60, 0.90].
- **F3** `σ*/κ` (derjaguin, `d₁₀`) ratio CV 0.5/CV 0 recomputed the same way: < 1.
For `μ_p` = 0 the question is whether the servo holds at all at N = 4000 (W4); if it does,
the frictionless CV 0.5 numbers are printed, not graded.
