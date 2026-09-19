# Item E — Why does elongation reduce the CV effect to ~0.6×?

**VERDICT: UNDETERMINED on the magnitude, SUPPORTED on the sign, and two of the three
natural geometric accounts are REFUTED.**

* **E1** — Donev's "singular at the sphere value" is a **frictionless** statement and
  does not act at μ_p = 0.3: **SUPPORTED** (a derivation, plus a confrontation the
  measurement could have failed).
* **E2** — free volume sets the size of the CV effect: **REFUTED**, predicts 1.015 where
  0.615–0.648 is measured, and does not even get the sign.
* **E3** — depopulation carries it: **REFUTED**, the rattler-count CV response is
  unchanged by elongation (0.996) while the density response falls to 0.62.
* **E4** — an effective size ratio set by the member diameter: **SIGN ONLY**, and it is a
  heuristic, not a derivation.

**What this anchor can and cannot decide.** It can reject accounts that make a
parameter-free prediction for the ratio, and it rejected two. It CANNOT identify the
mechanism: a two-sphere clump moves elongation, corrugation and contact multiplicity
together — `CLAIMS_SHAPE_2026-09-18.md` says so in its own opening paragraph, and the
registered LevelSet control is what would separate them.

Numbers printed by `geom_shape_2026-09-19.py` → `geom_shape_2026-09-19.txt`. Shape
shards: 126 records → 91 after (converged, not failed, N = 5400, tag ∈ {shp, shpcap}) →
91 after collapsing to one record per physical cell. The measured ratios
`D(α=2)/D(α=1)` are read from `out/grade_shape.txt` S3: bidisperse 0.648 ± 0.132,
lognormal 0.622 ± 0.022, uniform 0.615 ± 0.054.

---

## E1. The Donev singularity is frictionless, and does not act here

Generalise §1 of `geom_isostatic_backbone` to a rigid body with `n_f` **effective**
degrees of freedom and `n_c` constraints per contact:

> **z_iso = 2 n_f / n_c**

**FRICTIONLESS (`n_c = 1`, only the normal gap constrains):**

| body | n_f | why | z_iso |
|---|---|---|---|
| sphere | **3** | rotations are exact null modes — a smooth sphere spins with no contact penalty | **6 = 2d** |
| spheroid / two-sphere clump | **5** | spin about the symmetry axis maps the body to itself, the other two rotations do not | **10** |
| general ellipsoid | **6** | no symmetry axis | **12** |

The jump 6 → 10 at infinitesimal asphericity **is** Donev's singularity (Donev *et al.*,
*Science* **303**, 990 (2004), doi:10.1126/science.1093010): the DOF count changes
discontinuously while the shape changes continuously.

**FRICTIONAL (`n_c = 3 − x_m` on average):**

* sphere: rotations are **not** null modes — tangential friction resists rolling — so
  `n_f = 6` already;
* two-sphere clump: rotation about the long axis moves every off-axis contact point
  tangentially, so it is resisted too; `n_f = 6`.

> **z_iso = 12/(3 − x_m) for both. The frictional count is shape-blind.**

**PREDICTION E1: at μ_p = 0.3 there is no DOF-driven jump in z from α = 1 to α = 2. Any
measured rise must come from contact multiplicity or corrugation.**

Confrontation (z_nr per rung, mean over seeds, μ_p = 0.3, N = 5400):

| family | CV | p̂ | z_nr(α=1) | z_nr(α=2) | rise | dup_frac(α=2) |
|---|---|---|---|---|---|---|
| lognormal | 0.0 | 1e−3 | 5.128 | 6.383 | +1.254 | 0.1363 |
| lognormal | 0.0 | 1e−4 | 4.701 | 5.741 | +1.040 | 0.1147 |
| lognormal | 0.5 | 1e−3 | 4.968 | 5.470 | +0.502 | 0.1901 |
| lognormal | 0.5 | 1e−4 | 4.568 | 4.961 | +0.393 | 0.1754 |
| uniform | 0.5 | 1e−3 | 5.038 | 5.869 | +0.832 | 0.1624 |
| bidisperse | 0.5 | 1e−3 | 5.087 | 6.046 | +0.959 | 0.1517 |

**A DOF jump would be a fixed +4 with no CV and no pressure dependence.** The measured
rise is a quarter of that (+0.39 to +1.25) and moves with both. E1's prediction survives,
and the printed `dup_contact_frac` (0.11–0.19 of clump pairs carrying two contacts) names
the channel that does carry the rise.

*This is the one item E result I would defend in print:* the literature's standard
argument for "shape changes jamming through the contact number" is a frictionless DOF
argument, and it does not transfer to a frictional production grid.

**Three things from the source (citation pass 2026-09-19) that sharpen this.**

1. Donev's own quotes, verified verbatim: *"the contact number appears singular at the
   sphere value and rises sharply for small deviations"*; Z goes 6 → ~10 (9.82 measured
   on M&Ms, 9.80 simulated) and then **saturates** — consistent with `2n_f/n_c` = 10 for a
   body of revolution, not 12.
2. **Donev's φ is NON-monotone**: a cusp-like *minimum* at α = 1, maxima φ ≈ 0.71 at
   α ≈ 1.5 (prolate), then falling. ⚠ This project measures φ(CV = 0) *falling* from
   0.5978 (α = 1) to 0.5702 (α = 2) — i.e. **below** the sphere value, which Donev's
   frictionless curve does not do anywhere. The systems differ in three ways at once
   (friction, a two-sphere clump rather than a smooth spheroid, and a compression protocol
   rather than a Lubachevsky–Stillinger growth), so this is not a contradiction of Donev;
   it is a reason not to import his φ(α) curve into this project's interpretation. Note
   also that the project's own α = 2 is *past* the frictionless maximum at α ≈ 1.5.
3. Donev states the frictional counting himself — *"If friction is strong, then less
   contacts are needed, Z = f + 1"* — and the authors hedge the singularity throughout
   ("putative", "appears", "not discontinuous"). The companion analysis is Donev,
   Connelly, Stillinger & Torquato, *Phys. Rev. E* **75**, 051304 (2007),
   doi:10.1103/PhysRevE.75.051304, which explains why the isocounting conjecture fails for
   nonspherical particles.

---

## E2. Free volume — REFUTED, and it does not even get the sign

In the Furnas ideal limit (`derivations/furnas_ideal_limit.md`) a binary mixture reaches
`φ_max = β(2−β)`, so the filling gain over the monodisperse skeleton is (SymPy asserts)

    φ_max − β = β(1 − β).

If the CV effect is void filling and nothing else changes, its size scales with that.

> **PREDICTION E2: D(α=2)/D(α=1) = [β₂(1−β₂)] / [β₁(1−β₁)].**

| p̂ | β(α=1) | β(α=2) | β(1−β) α=1 | β(1−β) α=2 | E2 ratio |
|---|---|---|---|---|---|
| 1e−3 | 0.59780 | 0.56960 | 0.24044 | 0.24516 | 1.0196 |
| 3e−4 | 0.59034 | 0.56627 | 0.24184 | 0.24561 | 1.0156 |
| 1e−4 | 0.58220 | 0.56380 | 0.24324 | 0.24593 | 1.0110 |

**E2 prediction: 1.0154.** Measured 0.615–0.648, refuted at t = −17.9 (lognormal) and
t = −7.4 (uniform); the bidisperse error bar (±0.132) admits it at t = −2.78, which is a
statement about that family's precision, not about E2.

The reason is instructive: `β(1−β)` is flat near β ≈ 0.58 because it sits near its
maximum at β = ½, so **a 0.028 fall in the monodisperse packing fraction buys essentially
no change in available void.** Free volume predicts *no reduction at all*, so it fails on
the sign as well as the size.

---

## E3. Depopulation — REFUTED, and it fails the same way friction did

If the CV effect is made by rattler production, the CV response of the rattler *number*
fraction should shrink by the same 0.62.

| family | p̂ | Δq(α=1) | Δq(α=2) | ratio |
|---|---|---|---|---|
| lognormal | 1e−3 | +0.4072 | +0.4119 | 1.0117 |
| lognormal | 3e−4 | +0.4404 | +0.4335 | 0.9844 |
| lognormal | 1e−4 | +0.4443 | +0.4551 | 1.0244 |
| uniform | 1e−3 | +0.3736 | +0.3601 | 0.9637 |

**Mean 0.996 over four (family, rung) cells, range 0.964–1.024.** Against a density
response of 0.615–0.648: refuted at t = −17.0 and −7.1.

> **The depopulation COUNT is essentially unchanged by elongation while the DENSITY
> effect falls to 0.62.**

This is the *same shape* as item A's friction result (§6 there: the rattler-count CV
response is 1.3–2.6× larger at μ_p = 0.3 than at μ_p = 0 while φ_J's is 1.07× larger).
Twice now, on two independent perturbations, the rattler-count channel and the φ channel
have failed to scale together. **Depopulation counted in grains is not the carrier of
either effect.** Under item B's decomposition the carrier has to be the conditional
rattler VOLUME — the same quantity item B identifies as unmeasured.

---

## E4. An effective size ratio — the only account with the right sign

**Stated as the heuristic it is.** Whether a fine grain can enter a coarse pore is set by
the pore *throat*, and for a two-sphere clump the throat is plausibly set by the MEMBER
diameter `d_s`, not the equivalent-sphere diameter `d_eq`. `CLAIMS_SHAPE_2026-09-18.md`
reports `d_s/d_eq = 0.7937` at α = 2, so a fine grain of size `u·d_eq` sees an effective
size ratio

    u_eff = u · (d_eq/d_s) = u / 0.7937 = 1.2599 u,

i.e. **less** size separation, hence less filling. Right sign by construction; the
question is the size.

Evaluated on the bidisperse family (the only one where the linear packing model applies
directly): CV = 0.5 gives atoms 0.5 and 1.5, so `u = 1/3` and the fines volume fraction
is `X_s = 0.5³/(0.5³+1.5³) = 0.035714`. Exponents PA = 1.02, PB = 1.50
(`paper/numbers.tex` `\fillPa`, `\fillPb`), evaluated through the shipped
`mixture_model.binary_density`:

| p̂ | u = 0.3333, β(α=1) | gain | u_eff = 0.4200, β(α=2) | gain | ratio |
|---|---|---|---|---|---|
| 1e−3 | 0.5978 | +0.00906 | 0.5696 | +0.00715 | **0.7892** |
| 3e−4 | 0.5903 | +0.00895 | 0.5663 | +0.00711 | **0.7945** |
| 1e−4 | 0.5822 | +0.00882 | 0.5638 | +0.00708 | **0.8021** |

Measured bidisperse ratio 0.648 ± 0.132; the prediction 0.79 sits 1.1 se away.

**What this is and is not.** The map `u → u·d_eq/d_s` is **asserted, not derived** —
nothing here shows the throat of a clump pore scales with the member diameter. It is the
only one of the three accounts that predicts a reduction at all, and it lands in the right
region for the one family where the linear packing model applies. **That is a sign result
with a plausible magnitude, not a prediction of 0.62.**

The honest summary the brief invited: *only the sign.*

---

## What would settle it

The registered LevelSet control in `PREREG_SHAPE_2026-09-18.md` separates elongation from
corrugation and multiplicity, which is the blocker. Short of that, E4 makes one cheap,
falsifiable prediction that the current campaign does not test: **the reduction factor
should depend on the size ratio `u`, not just on α.** E4 predicts the effect is strongest
where `a(u)` is steepest — so running the shape axis at a *second* binary size ratio (say
r = 5, where u = 0.2 → u_eff = 0.252) would give the map a second point, and E2 and E3
both predict no ratio dependence at all. That is a two-cell discriminator.
