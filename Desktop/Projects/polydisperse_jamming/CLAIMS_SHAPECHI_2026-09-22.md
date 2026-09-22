# CLAIMS — `shapechi`: **shape is a window-knob.** All gates pass.

Graded 2026-09-22 against `PREREG_SHAPECHI_2026-09-22.md`, written before launch.
Artifact: `out/grade_shapechi.txt`. **72 of 72 records delivered, 72 of 72 converged.**

---

## The result

> **Non-sphericity nearly DOUBLES the shear-jamming window.**
>
> `χ = φ_m/φ₀` falls **0.8856 ± 0.0026 → 0.7884 ± 0.0043** from `α = 1` to `α = 2`.
> `Δχ = −0.0972 ± 0.0050`, **|t| = 19.3** against `t_crit = 4.303`.
> The window `1 − χ` goes **0.1144 → 0.2116, a factor 1.85.**

This is registered outcome **(W)**, the hypothesis: elongation constrains rotation, which
is mechanically what friction does, and friction acts only on `φ_m`.

## The mechanism is visible in the two branches separately

| | `α = 1` | `α = 2` | change |
|---|---|---|---|
| `φ₀` (frictionless) | 0.63269 | 0.62722 | **−0.9 %** |
| `φ_m` (frictional) | 0.56033 | 0.49453 | **−11.7 %** |

**The frictionless branch barely moves; the frictional one drops 13× more.** That is the
signature the hypothesis predicted, and it is why `χ` moves at all — a shape change that
lowered both branches equally would have been a pure **x**-knob, like the size
distribution.

## Robustness

| gate | result |
|---|---|
| **W1a** wiring | **3/3 BIT-IDENTICAL** to the pre-existing records (`\|Δφ\| = 0`, same `iters`) |
| **W1b** the prereg's transcribed literals | 3/3 within their own 6-dp precision |
| **W2** health | **72/72 converged** at `1e6` |
| **S3** rung consistency | **resolved at ALL FOUR rungs, same sign**: `Δχ` = −0.094 / −0.075 / −0.085 / −0.097, `\|t\|` = 25.8 / 20.4 / 24.0 / 19.3 |

S3 matters because `χ` is **not** rung-independent (it drifts 0.016–0.025 across the
ladder). The effect here is ~0.09, roughly **four times** that systematic, and it is the
same sign at every rung — so unlike the CV comparison, this verdict is not a fixed-rung
artefact.

> 🔴 **W1 FAILED ON ITS FIRST READING, AND THE GATE WAS THE PROBLEM.** The prereg quotes
> the three reference values to six decimals, and the gate demanded `|Δ| == 0.0` against
> those transcriptions — **a bar no correct run could pass**, since the true values are
> not 6-dp numbers. Measured against the full-precision records on disk the re-run is
> **bit-identical**; the ~3 × 10⁻⁷ discrepancy was exactly the rounding. Replaced by a
> **stricter** test (bit-identity on `phi`, `iters` and `n_contacts`), with the literal
> kept as a separate weaker check. The bar was not loosened.

## The bench consequence

At fixed concentration, a suspension sits at `x = φ/φ₀`. Between the two `χ` values it
**cannot shear-jam as spheres and can as clumps**:

&nbsp;&nbsp;&nbsp;&nbsp;**`0.7884 ≤ x < 0.8856`**, i.e. at `φ₀ ≈ 0.633`, **`φ` from 0.499 to 0.560**.

> **A 0.06-wide band in absolute volume fraction in which making the particles
> non-spherical — at the same loading, the same chemistry, the same everything else —
> flips the material from "thickens continuously" to "shear-jams".** At the midpoint
> `x = 0.837`, spheres give a finite `Θ` (8.8 at `n = 2`, 231 at `n = 5`) and clumps
> diverge.

That is directly testable, and it is the sharpest prediction this framework has produced.

## What this does NOT establish

- ⚠ **A clump moves three things at once** — elongation, surface corrugation and contact
  multiplicity. `CLAIMS_SHAPE_2026-09-18` already records that no clump-only run separates
  them. **So this is "non-sphericity widens the window", NOT "elongation does."** The
  registered LevelSet control is what would separate them and was not run.
- ⚠ **Two `α` rungs show a DIRECTION, not a shape.** The records already on disk have
  `α` **non-monotone** (`φ = 0.5957 / 0.6289 / 0.5661` at `α = 1.0 / 1.5 / 2.0`,
  `μ_p = 0.3`, `p̂ = 1e-3`), so `α = 1.5` may sit on the other side of a maximum from both
  points measured here. **No curvature claim.**
- ⚠ `N = 1700`, not the shape campaign's 5400. `χ` is a fixed-`N` ratio and the comparison
  is at one `N` throughout, but **no cell here is comparable to the `N = 5400` block**.
- ⚠ `μ_p = 3.0` stands in for `φ_m`. That stand-in was cross-checked elsewhere against an
  independent `φ_∞` extrapolation and agreed to under `PHI_SCALE`, but it remains a
  stand-in.

## Where this leaves the axis map

**Two window-knobs are now measured — friction and shape — and one position-knob,
the size distribution.** The remaining axes (size, surface potential) both act through
`σ*`, and **neither has a measurement**. That is now unambiguously the next campaign.
