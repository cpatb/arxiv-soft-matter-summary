# The binary filling maximum: closed forms, limits, and what the finite ratio does

Checked by `furnas_ideal_limit_check.py` → `furnas_ideal_limit_numbers.txt` (this
directory). Every equation below is either printed by that script from SymPy or is an
identity it asserts to zero. Written 2026-09-14.

## Setup

Two size classes, coarse (1) and fine (2); fines volume fraction `x`; size ratio
`u = d₂/d₁ ≤ 1`; monodisperse packing density `β`. The linear packing model
(`mixture_model.packing_density`) reduces for two classes to

    γ₁ = β / (1 − (1 − a) x)                       coarse-dominant branch
    γ₂ = β / (1 − [1 − β + b β (1 − 1/β)] (1 − x))  fine-dominant branch
    φ  = min(γ₁, γ₂)

with the loosening `a(u)` and wall `b(u)` interaction functions in [0, 1],
`a(0) = b(0) = 0`, `a(1) = b(1) = 1`.

**Step 1 — the fine-branch bracket.** Expand `b β (1 − 1/β) = b β − b`, so
`1 − β + bβ − b = (1 − β) − b(1 − β) = (1 − β)(1 − b)`. SymPy asserts the difference is
zero. Hence

    γ₂ = β / (1 − (1 − β)(1 − b)(1 − x)).

Reading: on the fine-dominant branch the coarse grains are inclusions whose volume
`(1 − x)` displaces a void fraction `(1 − β)` of the fine matrix, reduced by the wall
effect `b`.

**Step 2 — monotonicity.** `γ₁`'s denominator `1 − (1 − a)x` decreases in `x` (for
`a < 1`), so `γ₁` increases in `x`. `γ₂`'s denominator `1 − (1 − β)(1 − b)(1 − x)`
increases in `x`, so `γ₂` decreases in `x`. The minimum of an increasing and a
decreasing function is maximised at their crossing.

**Step 3 — the crossing.** Setting `γ₁ = γ₂` gives one root (SymPy `solve`, one
solution):

    x* = (1 − b)(1 − β) / (2 − β − a − b(1 − β))

and the value there

    φ_max = β (2 − β − a − b(1 − β)) / (1 − aβ − ab(1 − β)).

(These are the script's factored forms `-(b-1)(β-1)/(a - bβ + b + β - 2)` and
`-β(a - bβ + b + β - 2)/(abβ - ab - aβ + 1)`, with the signs collected.)

## The ideal limit `a = b = 0` (perfect size separation, `u → 0`)

    x*     = (1 − β) / (2 − β)
    φ_max  = β (2 − β)

SymPy asserts both against the general expressions. Interpretation: the fines fill the
coarse void `(1 − β)` at density `β`, so the added solid is `β(1 − β)` per unit of the
coarse pack, and `φ_max = β + β(1 − β) = β(2 − β)`; the fines' share of the solid is
`β(1 − β)/[β + β(1 − β)] = (1 − β)/(2 − β)`.

Limits (script §4):
- `β → 1`: `x* → 0`, `φ_max → 1` — nothing to fill.
- `β → 0`: `x* → 1/2`, `φ_max/β → 2` — everything is void, each class fills the other.
- `a = b = 1` (`u = 1`, one size): `γ₁ = γ₂ = β` at every `x`, so `φ ≡ β`.

Numbers at this project's `β` (script §6, against the shipped `furnas_ideal`, agreement
to 5e−6 on a 100 001-point grid):

| β | x* | φ_max |
|---|---|---|
| 0.5947 (`μ_p` = 0.3) | 0.2884 | 0.8357 |
| 0.6413 (`μ_p` = 0)   | 0.2640 | 0.8713 |
| 0.55                 | 0.3103 | 0.7975 |

## Finite ratio

Both partials of `φ_max` are negative on the open unit square (script: max over a
60 × 60 grid is −9.8e−8 in `a` and −1.5e−6 in `b` at `β` = 0.5947, likewise at the other
two `β`), and `φ_max ≥ β` everywhere. So any interaction — loosening or wall — lowers the
achievable maximum, and the maximum is monotone in the size ratio whenever `a(u)`, `b(u)`
are monotone in `u`, which the shipped forms are.

Special cases (script §4 iv):

    φ_max(a, b = 0) = β (2 − β − a) / (1 − aβ)          loosening only
    φ_max(a = 0, b) = β (2 − β − b(1 − β))              wall effect only

The wall effect enters linearly; the loosening enters through a denominator and so bites
harder at large `a`.

With the shipped exponents `(PA, PB) = (0.961, 2.214)` the general crossing reproduces
the argmax of `binary_density` to `|Δφ| ≤ 6e−6`, `|Δx| ≤ 2e−5` at `r` = 2, 3, 4, 5
(script §6):

| r | a | b | x* | φ_max |
|---|---|---|---|---|
| 2 | 0.698 | 0.785 | 0.224 | 0.638 |
| 3 | 0.568 | 0.593 | 0.277 | 0.675 |
| 4 | 0.492 | 0.471 | 0.297 | 0.700 |
| 5 | 0.439 | 0.390 | 0.306 | 0.718 |
| ∞ | 0 | 0 | 0.288 | 0.836 |

The crossing composition is not monotone in `r`: it rises from 0.224 towards ~0.31 and
must come back to 0.288 as `a, b → 0`, because the two interaction functions do not
vanish at the same rate (`b` falls faster than `a` under these exponents, and `x*`
decreases in `b` and increases in `a`: `∂x*/∂a > 0`, `∂x*/∂b < 0` on the open unit
square, script §4 v, 60 × 60 grid at three `β`).

## What this does and does not establish

- It is the model's own algebra, so it certifies the code's branch structure, not the
  physics. The physics test is the extension block (`PREREG_MIXTURE_EXT_2026-09-14.md`).
- The measured `r = 3` maximum, 0.6763 at `X_s` = 0.30 (seed mean, N = 1500), sits
  0.0009 above the model's 0.6754 — the model's crossing at 0.277 and the grid's 0.30
  are within one grid step, so the location is not resolved by the r = 3 data.
