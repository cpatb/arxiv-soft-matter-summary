# PREREG — `shapechi`: is SHAPE a window-knob or a position-knob?

Written 2026-09-22 **before launch**. Framework: `three_knobs.py` / `out/three_knobs.txt`.

## 1. The question, in the reduced variables

The Wyart–Cates closure the drop-impact DNS solves is **exactly covariant** under a common
rescaling of `(φ, φ₀, φ_m)`, so it depends on only

&nbsp;&nbsp;&nbsp;&nbsp;`x = φ/φ₀`  (position) `χ = φ_m/φ₀` (contrast) and `σ*` (onset).

**`χ` decides whether a suspension can shear-jam at all**: the window is `χ ≤ x < 1`, of
width `1 − χ`. So for every particle property the question "does it change the
*thickening*, or only *where you sit*?" is the single question **"does it move `χ`?"**

- **Friction moves `χ` by construction** — `φ₀` is the frictionless point, so `μ_p` acts on
  `φ_m` alone. Measured: `χ = 0.8872 ± 0.0042` at `μ_p` 0 → 3, window `1 − χ = 0.113`.
- **The size distribution does NOT move `χ`**: `Δχ = +0.0019 ± 0.0049`, `|t| = 0.40`
  against `t_crit` 4.303 across CV 0 → 0.2. It is a pure **x**-knob.

> **SHAPE is the remaining candidate `χ`-knob, and this repository cannot currently
> measure it: the entire shape campaign (126 records) ran at `μ_p = 0.3` ONLY**, so it
> measures neither `φ₀` nor `φ_m`. This campaign fixes exactly that, and nothing else.

**Why shape should be a `χ`-knob** (the hypothesis, not a result): elongation and
interlocking constrain particle ROTATION, which is mechanically what friction does, and
friction acts only on `φ_m`. If that is right, `α > 1` lowers `φ_m` more than `φ₀` and the
window **widens**.

## 2. Grid — deliberately the smallest that can answer it

`family = lognormal`, **`CV = 0`** (monodisperse, so the distribution axis is switched off
and shape is alone), `N = 1700`, seeds {1,2,3}, the standard 4-rung `p̂` ladder,
`PJ_MAX_ITER = 1000000`, tag **`shapechi`**:

&nbsp;&nbsp;&nbsp;&nbsp;`α ∈ {1.0, 2.0}` × `μ_p ∈ {0.0, 0.3, 3.0}` × 3 seeds × 4 rungs = **72 records**

⚠ **`N = 1700`, not the shape campaign's 5400, and that is a deliberate simplification.**
Timing on disk: `α = 2` costs a median **5497 s at N = 5400** and **233 s at N = 1700** —
24× for 3.2× the grains. `χ` is a **ratio at fixed `N`**, and the comparison registered
below is `χ(α=2)` against `χ(α=1)` **at the same `N`**, so the leading `N`-dependence
cancels. The price is that this campaign cannot be compared cell-for-cell against the
`N = 5400` shape block, and no gate below does so.

## 3. Gates

### W1 — WIRING. **Primary; gates everything.**

Three `N = 1700, CV = 0, μ_p = 0.3, seed = 1` records already exist and are re-run here.

> **REGISTERED: they must come back BIT-IDENTICAL.**
> `α = 1.0, p̂ = 1e-3 → φ = 0.595694`; `α = 1.0, p̂ = 3e-5 → φ = 0.575051`;
> `α = 2.0, p̂ = 1e-3 → φ = 0.566148`.
>
> ⚠ Compared on `phi` to full precision **and** on `iters`/`n_contacts`. This session
> established that cross-code-era records differ at ~1e-14 while the trajectory is
> identical; if that appears here it is **reported**, and the gate is read on whether
> `iters` and `n_contacts` match exactly.

### W2 — HEALTH.

> **REGISTERED: 72 of 72 converge at `1e6`.** Any cell that caps is reported, not
> re-raised. ⚠ The `μ_p = 0` arm is the one at risk — frictionless packs are what still
> cap at `1e6` elsewhere in this project, though only at high CV, and this campaign is at
> `CV = 0` where they converged 3/3.

### S1 — **PRIMARY. Does shape move the window?**

`χ(α) = φ(μ_p=3.0, α)/φ(μ_p=0, α)`, seed-paired, converged only, at the matched rung
`p̂ = 3e-5`. Bar: `max(seed-paired s.e., PHI_SCALE·√(1+χ²)/(φ₀·√n))`, `t_crit = 4.303`
at 2 dof.

> **REGISTERED, two-sided, and both outcomes are results:**
>
> **(W) WINDOW-KNOB** — `χ(α=2) < χ(α=1)` and resolved at `|t| > 4.303`. Shape widens the
> shear-jamming window; elongation acts like friction. **This is the hypothesis.**
>
> **(P) POSITION-KNOB** — `|Δχ|` not resolved. Shape joins the size distribution: it moves
> `φ₀` and `φ_m` together and changes only where you sit. **This would be the more
> surprising and the more useful answer**, because it would leave *friction as the only
> measured window-knob* and make the whole programme a friction-and-onset story.
>
> A resolved `χ(α=2) > χ(α=1)` — shape NARROWING the window — is registered as the third
> outcome and would refute the rotational-constraint argument outright.

### S2 — REPORTED, NOT GATED: the `α` effect on each branch separately.

`φ₀(α)` and `φ_m(α)` at the quoted rung, so that a `Δχ` can be read as *which branch
moved*. ⚠ The existing `μ_p = 0.3` records already show `α` is **NON-MONOTONE**
(`φ = 0.5957 / 0.6289 / 0.5661` at `α = 1.0 / 1.5 / 2.0`, `p̂ = 1e-3`), so a two-rung `α`
ladder can show a DIRECTION and not a SHAPE, and **`α = 1.5` may sit on the other side of
a maximum from both points measured here.** No curvature claim may be made.

### S3 — REPORTED: `χ` across all four rungs.

`χ` is **not** rung-independent — measured drift 0.016–0.025 across the ladder, *larger*
than `φ₀`'s own 0.005–0.010, because the two friction arms have different
compressibilities. So S1's comparison is **at a fixed rung by construction**, and S3
prints whether the `α` verdict is the same at every rung. **A verdict that flips between
rungs is reported as rung-dependent, not resolved.**

## 4. What would make this campaign worthless

- **W1 fails** → the runner is not reproducing its own records and nothing may be read.
- **The `μ_p = 0` arm caps** → `χ` is unmeasurable at that `α` and S1 is undecidable there.
- **S1 resolves but flips sign between rungs** (S3) → the ladder systematic is larger than
  the effect and the verdict is rung-dependent.
- ⚠ **A clump moves three things at once** — elongation, surface corrugation and contact
  multiplicity — and `CLAIMS_SHAPE_2026-09-18.md` already records that no clump-only run
  separates them. **So a (W) verdict says "non-sphericity widens the window", NOT
  "elongation does".** The registered LevelSet control is what would separate them, and it
  is not run here.

## 5. What is NOT being run

- No `CV > 0`: the distribution axis is deliberately switched off.
- No third `α`: two rungs is what a direction needs, and the non-monotonicity in S2 means
  a third would need to be chosen to bracket the maximum, which is a different campaign.
- No `N` ladder: `χ` is a fixed-`N` ratio (§2).
- No `σ*` claim: shape does not enter the onset stress in this framework.
