# Item C, part 1 — Which class of functionals does the matched binary pair exclude?

**VERDICT: SUPPORTED.** The level set `{mean, CV, support}` on two-point PSDs is exactly
the pair of roots of one quadratic; every standardized moment of a two-point law is a
function of the composition alone; therefore a two-point PSD cannot hold CV and *any*
standardized moment fixed at two compositions. The manuscript's assertion at
`paper/main.tex:202` is here **proved**, and the class it excludes is stated precisely.

**What this anchor can and cannot decide.** It is pure moment algebra, so it decides
exactly which closures the pair excludes with no reference to any packing measurement. It
CANNOT say which functional is the right one — indeed it shows the pair is structurally
incapable of identifying one. That is a statement about the design, not about packing.

Numbers printed by `geom_psd_functional_class_2026-09-19.py` →
`geom_psd_functional_class_2026-09-19.txt`.

---

## 1. The level set is a quadratic, and its two roots are the manuscript's A and B

Atoms `d_s = s`, `d_l = r s` with `r > 1`; `q` = number fraction of the **large** atom.

    ⟨d⟩   = s(qr − q + 1)                        → fix ⟨d⟩ = 1: s = 1/(qr − q + 1)
    ⟨d²⟩  = s²(qr² − q + 1)
    CV²   = ⟨d²⟩/⟨d⟩² − 1 = −q(q−1)(r−1)²/(qr−q+1)²

Setting CV² = c² and clearing denominators gives a **quadratic in q** (SymPy asserts the
difference from the canonical form is exactly 0):

> **(1+c²)(r−1)² q² + (r−1)[2(1+c²) − r − 1] q + c² = 0**

Vieta, exact:

    q₁ q₂   = c² / [(1+c²)(r−1)²]
    q₁ + q₂ = (r − 1 − 2c²) / [(1+c²)(r−1)]
    discriminant = (r−1)²(r² − 2r + 1 − 4c²r)

So **two compositions realise any admissible (mean, CV, ratio)**, and no more. The level
set is a *pair*, not a continuum — the manuscript's pair exhausts what can be matched this
way. Equal roots is the tangency where the experiment becomes impossible.

**At r = 3, CV = 1/2 the quadratic is `5q² − 3q + 1/4 = 0`, with roots q = 1/10 and
q = 1/2** — reproducing the manuscript's pair as an *output*:

| root | n_small | d_s | d_l | realised mean | realised CV | skewness | X_s (fines volume) |
|---|---|---|---|---|---|---|---|
| q = 1/10 | 0.900 | 0.833333 | 2.500000 | 1.0000000000 | 0.5000000000 | **+2.666667** | **0.250000** |
| q = 1/2 | 0.500 | 0.500000 | 1.500000 | 1.0000000000 | 0.5000000000 | **+0.000000** | **0.035714 = 1/28** |

Against `paper/numbers.tex`: `\pairNsmallA` 0.500 / `\pairNsmallB` 0.900, `\pairXsA`
0.0357 / `\pairXsB` 0.25, `\pairSkewA` 0.00 / `\pairSkewB` +2.67. Reproduced to the
printed precision.

---

## 2. Skewness of a two-point PSD depends only on the composition

For atoms `a < b` with `P(b) = q`, SymPy gives

    var = q(1−q)(b−a)²
    μ₃  = q(1−q)(b−a)³(1−2q)
    γ₁  = μ₃/var^{3/2}

and the claim `γ₁ = (1−2q)/√(q(1−q))`, **free of a and b**, checks to 4.4e−15 over a
72-point grid in (a, b, q). Since

    dγ₁/dq = −√(1−q) / [2 q^{3/2} (q−1)²]  < 0 on (0,1),

`γ₁` is strictly decreasing and therefore **injective**. Two distinct compositions give
two distinct skewnesses. **A two-point PSD cannot hold CV and skewness fixed at two
compositions. QED.**

**Corollary, and it is the stronger statement.** Any two-point law is an affine image of
Bernoulli(q), and standardized moments are affine-invariant, so *every* standardized
moment is a function of q alone. Verified for the fourth: kurtosis
`= (1 − 3q(1−q))/(q(1−q))`, worst deviation 1.6e−14 over 24 points.

---

## 3. The excluded class, stated precisely

Let `L = { P : ⟨d⟩ = 1, CV = ½, supp(P) = {d_s, d_l} }`. The manuscript's pair is
`{P_A, P_B} ⊂ L`, and φ_A ≠ φ_B by 0.0635 ± 0.0032.

**EXCLUDED — every functional constant on L:**

* any function of the first two NUMBER moments alone (⟨d⟩, CV, variance, σ/μ, span);
* any function of the SUPPORT alone (size ratio r, d_max/d_min, d_max − d_min, species
  count K);
* any function of the two jointly;
* hence every closure of the form `φ_J = f(CV)`, `f(CV, r)`, `f(span)`.

**NOT EXCLUDED — and this is where "no descriptor of the size distribution" overstates
it.** Anything reading the third or higher moment survives (skewness differs by 8/3
across the pair), and so does any VOLUME-weighted moment, because volume weighting is
moments 3 and 4:

| | d₃₂/d₁₀ | d₄₃/d₁₀ | X_s | volume-weighted CV |
|---|---|---|---|---|
| q_large = 0.50 (A) | 1.400000 | 1.464286 | 0.035714 | 0.126735 |
| q_large = 0.10 (B) | 1.666667 | 2.083333 | 0.250000 | 0.346410 |

All four differ. So the exact statement is:

> **The pair excludes every closure that is a function of the number-weighted mean, the
> number-weighted CV and the support. "No descriptor of its WIDTH works" is exact; "no
> descriptor of the size distribution works" is too strong.**

**And the limit of the pair as evidence.** By §2, within two-point PSDs, matching
(mean, CV, support) and varying composition necessarily varies *every* higher standardized
moment simultaneously. The pair therefore cannot attribute the 0.0635 to the third moment
rather than the fourth, or to any single moment at all. It is a proof of the insufficiency
of width, not an identification of the sufficient descriptor.

---

## 4. The design that would identify the moment

A K-atom PSD on a fixed support has K−1 free composition parameters; fixing ⟨d⟩ and CV
leaves K−3.

| K | free parameters after (mean, CV) | what it buys |
|---|---|---|
| 2 | 0 | the two Vieta roots — two skewnesses and nothing between them |
| 3 | 1 | a **curve** of PSDs at matched mean, CV and support, along which skewness varies continuously |
| 4 | 2 | a surface — skewness and kurtosis can be varied **independently**, the first genuine discriminator between them |

Constructed and measured (§5), support `{1, u, u²}` rescaled to ⟨d⟩ = 1 at CV = 0.500 ±
0.002:

| support | solutions found | skewness range |
|---|---|---|
| {1, 2, 4} | 4067 | +2.900 to +4.565 (range 1.665) |
| {1, 3, 9} | 4005 | +0.010 to +15.334 (range 15.325) |

Against two isolated skewnesses (0.000 and +2.667) from the binary pair. The project
already has the machinery: `psd.make_psd` accepts a `"discrete:d1/2/4:y.../..."` spec and
derives CV from the composition, so this is a composition grid at fixed number-CV and no
new code.

*(Caveat for a follow-up: the {1,3,9} range is dominated by compositions that are almost
degenerate — nearly all the weight on one atom — and a usable ladder needs a minimum
occupancy per atom. That is a design decision, not a result of this file.)*
