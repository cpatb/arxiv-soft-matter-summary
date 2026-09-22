# The stress an impact needs to JAM a suspension, and which particle property enters where

Written 2026-09-22 (third block). Checks: `jam_stress_2026-09-22_check.py` →
`jam_stress_2026-09-22_check.txt` (**11/11 PASS**). Code: `actshape.py`.
Status: **derived and checked; referee pass 1 returned 1 FATAL + 2 MATERIAL, all in the
interpretive layer (§0, §5), all acted on below; the algebra was independently re-derived and
stands. Pass 2 not yet run.**

> **RESULT.** A suspension inside the shear-jamming window jams wherever the local stress exceeds
>
> &nbsp;&nbsp;&nbsp;&nbsp;**`σ_jam = κ_m · d^(m−2) / Q_{1−f_c}(ĝ_m)`**, &nbsp;&nbsp;`f_c = (1 − x)/(1 − χ)`, &nbsp;&nbsp;(8)
>
> with `x = φ/φ₀`, `χ = φ_m/φ₀`, `d` the particle size, `F_rep = κ_m R*^m` the repulsive
> threshold, and `ĝ_m` the dimensionless contact-force variable of §2 whose distribution the
> DEM measures. **The factors are:** packing fraction through `x`; friction and shape
> through `χ` **and, unmeasured, through `ĝ`**; surface potential through `κ_m`; size through
> `d^(m−2)` **only if `ĝ`'s distribution is scale-free (A4b, untested)**; the size
> DISTRIBUTION through `x` and through the distribution of `ĝ` (its onset scale and shape). For the solver's
> Weibull activation (8) reduces to `σ_jam = σ*[ln(1/f_c)]^(−1/δ)`, and the shape `δ` has a
> **sign-reversing** effect: sharper activation needs MORE stress to jam when `f_c < 1/e`
> (deep in the window) and LESS when `f_c > 1/e` (near its edge).

## 0. Frame

**Sought:** the stress above which a material point of a Wyart–Cates suspension is jammed,
as a function of every particle property, so that "does the drop jam or stay
inertial-dilatant" becomes a comparison of the impact stress with one number.

**Assumptions.**
- (A1) Wyart–Cates with the solver's linear interpolation `φ_J(f̂) = φ₀ − (φ₀ − φ_m) f̂`.
  Project convention (`three_knobs.py`): `x = φ/φ₀`, `χ = φ_m/φ₀`.
- (A2) **Quasi-steady activation**: `f̂ = f̂_ss(σ)`. The solver's measured clock relaxes `f̂`
  over a strain `γ₀ = 0.2`; an impact accumulates O(1) strain, so this is plausible but
  **not derived** — the `jamwin` DNS is what tests it, with the clock on.
- (A3) A contact is frictional iff its normal force exceeds `F_rep = κ_m R*^m`,
  `R* = R₁R₂/(R₁+R₂)`. `m = 1` is Derjaguin (any flat-plate law); `m = 0` is a
  size-independent threshold, which is what the one measured size series implies (Guy,
  Hermes & Poon, PRL 115, 088304 (2015), doi:10.1103/PhysRevLett.115.088304: `σ* ∝ d^−1.9`,
  sourced and discussed — as an OPEN problem with four live caveats — in
  `onset_size_degeneracy_2026-09-19.md` §6.2). **`m` is not settled; both are carried.**
- (A4) At fixed packing and friction, contact forces scale with the applied stress times a
  size-squared: `f = σ d² ĝ`, `ĝ` a dimensionless random variable whose distribution depends on
  the packing (CV, shape, `μ_p`, `x`) but **not** on `σ` (the "rescaling" assumption; TO BE
  tested by `PREREG_ACTSHAPE` H4 over a decade of `p̂` — which, in a pressure-servoed cell,
  moves the pack's compression together with `σ`, so a pass cannot separate the two).
- (A4b) **`ĝ_m`'s distribution does not depend on the absolute size scale**, only on the
  distribution's shape. Needed for "size enters only as `d^(m−2)`". **Untested for SHEARED
  packs**: `test_scale_invariance.py` covers static isotropic packs, and lessons.md records
  that a static scale-invariance does not transfer to the dynamics for free. In a real
  suspension Brownian and hydrodynamic forces break it outright (they carry their own size
  dependence), so (A4b) is a statement about the frictional-contact network only.
- (A6) **`f̂` is the UNWEIGHTED fraction of contacts, by count**, pooled over all pair types.
  This is a modelling choice, not a consequence of (A3)–(A4): it matches Wyart–Cates'
  constraint-counting motivation, but `φ_J` is a volume fraction and a size- or
  coordination-weighted fraction is a defensible alternative. Untested against one.
- (A5) Gauge: `σ` is the pack's mean pressure in the DEM. The solver's `σ` is a shear-stress
  magnitude; `τ = μ_M p` relates them in steady flow. Named, not assumed away.

## 1. The window (C1)

Jamming is `φ ≥ φ_J(f̂)`:

`φ ≥ φ₀ − (φ₀ − φ_m) f̂`   (1)
`⇔ x ≥ 1 − (1 − χ) f̂`     (2)  [divide by φ₀ > 0; substitute φ_m = χφ₀]
`⇔ f̂ ≥ (1 − x)/(1 − χ) ≡ f_c`   (3)  [1 − χ > 0 for χ < 1]

So `x ≥ 1` gives `f_c ≤ 0` (jammed at rest); `x < χ` gives `f_c > 1` (unreachable, since
`f̂ ≤ 1`); the window is `χ ≤ x < 1`. **C1: SymPy solves (1) for `f̂` and returns (3).**

## 2. Activation from the contact-force distribution (general form)

By (A3) and (A4), a contact is frictional iff `σ d² ĝ > κ_m (d r̂)^m`, with `r̂ = R*/d`:

`σ > κ_m d^(m−2) r̂^m / ĝ`.   (4)

Define `ĝ_m = ĝ / r̂^m` (per contact; in DEM units `d₁₀ = 1` this is `actshape.g_values` with
`norm = "constant"` for `m = 0` and `"derjaguin"` for `m = 1`). Then the fraction of
frictional contacts is

`f̂(σ) = P( ĝ_m > κ_m d^(m−2)/σ )`,   (5)

a survival function evaluated at `κ_m d^(m−2)/σ`. Setting (5) equal to `f_c` and inverting the
survival function gives the quantile:

`κ_m d^(m−2)/σ_jam = Q_{1−f_c}(ĝ_m)`  ⇒  (8).

(8) needs **no Weibull assumption and no σ\* convention** (**C6**: brute-force crossing of the
empirical survival of a deliberately non-Weibull sample against (8), worst 2.5e-5).
**Size exponent (C7b):** SymPy on (4) gives `d ln σ_thr/d ln d = m − 2`, i.e. `σ* ∝ 1/d` for
Derjaguin and `∝ 1/d²` for a fixed force. **Dimensions (C7a, pint):** `g = f/(R*σ)` is a
length and `κ/Q(g)` a pressure.

## 3. The solver's Weibull form

The solver writes `f̂_ss = exp(−(σ*/σ)^δ)`, i.e. `ĝ_m` Weibull with shape `δ`, and
`σ* ≡` the stress at which `f̂ = 1/e`. Setting `f̂_ss = f_c`:

`(σ*/σ_jam)^δ = ln(1/f_c) ≡ L`   (6)  [take −ln of both sides; 0 < f_c < 1 so L > 0]
`σ_jam = σ* L^(−1/δ)`.   (7)

**C2** (SymPy, with `f_c = e^−L`, `L > 0`): substituting (7) back returns `f_c` exactly.
**C2b:** `d f̂_ss/dσ = (δ/σ)(σ*/σ)^δ f̂_ss > 0`, so `f̂_ss ≥ f_c` exactly for `σ ≥ σ_jam`.
**C5** (brute-force threshold on a 2e6-point grid, 4 parameter points): worst 6.7e-6.

**The sign of δ's effect (C3).** From (7), `ln(σ_jam/σ*) = −(1/δ) ln L`, so

`∂ ln(σ_jam/σ*)/∂δ = ln L / δ²`,   (9)

positive iff `L > 1` iff **`f_c < 1/e`**. SymPy returns (9) exactly; evaluated at
`δ = 1.3`: +0.49 (`f_c` = 0.1), +0.19 (0.25), 0 (1/e), −0.74 (0.75).

⚠ **The crossing at exactly `1/e` is a property of the σ\* convention, not of physics.**
σ\* is *defined* as the stress at `f̂ = 1/e`, so two activation laws compared at matched σ\*
cross there by construction. The physical content is the comparison at matched `κ` and `d`,
(8), where the σ\*-shift and the shape both enter and no crossing is guaranteed.
`PREREG_ACTSHAPE` H6 asks the DEM whether the matched-`κ` curves cross.

## 4. Limits (C4, C8)

| limit | (7) gives | meaning |
|---|---|---|
| `δ → ∞` | `σ_jam → σ*` for every `f_c` | step activation: all contacts switch together |
| `δ → 0`, `f_c < 1/e` | `σ_jam → 0` | infinitely broad: some contacts are always active |
| `δ → 0`, `f_c > 1/e` | `σ_jam → ∞` | …and some never are |
| `f_c → 0` (`x → 1`) | `σ_jam → 0` | jammed at rest |
| `f_c → 1` (`x → χ`) | `σ_jam → ∞`, as `(1 − f_c)^(−1/δ)` (C8) | the lower window edge; a broader activation (smaller δ) diverges FASTER |

All six recovered by SymPy.

## 5. What this says about the question

A drop jams where its local stress exceeds (8). Reading the particle properties off (8):

| property | enters through | kind of effect |
|---|---|---|
| packing fraction `φ` | `x` → `f_c` | position in the window |
| friction `μ_p`, shape `α` | `χ` → `f_c`, and `ĝ` (unmeasured; `actshape` H3 measures friction's) | whether a window exists and how wide (measured: `three_knobs.py`) |
| surface potential | `κ_m` | pure stress scale; cannot change WHETHER, only at what stress |
| size `d` | `d^(m−2)`, if (A4b) | pure stress scale, exponent `m − 2` ∈ {−2, −1} |
| **size distribution CV** | `x` (φ₀ rises) **and** the distribution of `ĝ_m` | position; and, through `ĝ_m`, onset scale and activation shape |

The size distribution is the only property **measured** to enter (8) in more than one place.
Friction and shape may too (through `ĝ`); that is unmeasured, and `actshape` H3 measures the
friction half. ~~"the only property that enters (8) in more than one place"~~ was the v1 wording
and was withdrawn on referee pass 1 (FATAL: it contradicted the result box's own "weakly through `ĝ`",
which was itself unsupported). **Whether its shape
effect is large enough to matter is what `PREREG_ACTSHAPE` measures; whether shape matters
for a drop at all is what `PREREG_JAMWIN` measures.** Neither result is in this document.

## 6. Open

- (A2) quasi-steadiness: tested only by `jamwin`.
- (A4) rescaling in `σ`: `actshape` H4 (confounded with compression, see A4).
- (A4b) size-scale invariance of the SHEARED `ĝ`: no test planned yet.
- (A6) count-weighted `f̂`: no alternative computed yet.
- **The solver is rate-driven, not stress-controlled.** `wc_relax` activates on `η(f̂)γ̇`, so
  in a DNS the threshold is set by a bottleneck `f*` below `f_c`, and (7)'s reversal at
  `f_c = 1/e` moves to `f* = 1/e`: `jam_rate_driven_2026-09-22.py`, `PREREG_JAMWIN`
  Amendment 1. (8) describes a stress-controlled material point; reading a DNS with it needs
  that restatement.
- `m`: unsettled; Derjaguin says 1, the one size series says ≈ 0. (8) carries both.
- The gauge (A5) turns a DEM `p` into a solver `τ` through `μ_M`, which moves with CV and `μ_p`.
