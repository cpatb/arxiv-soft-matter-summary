# The programme: five particle axes, three numbers, one collapse test

**The question.** How do a particle's **size**, **shape**, **size distribution**, **surface
potential** and **material properties** change the way a suspension shear-thickens or jams
when a drop of it hits a surface?

**The claim of this document.** That question has a two-line answer, it is testable, and
most of what this repository currently contains is not needed to test it.

---

## 1. The reduction — why the answer is three numbers

The Wyart–Cates closure the drop-impact DNS solves,

&nbsp;&nbsp;&nbsp;&nbsp;`η = η_s (1 − φ/φ_J)⁻²`, `φ_J(f̂) = φ₀ − (φ₀ − φ_m) f̂`, `f̂ = f(σ/σ*)`,

is **exactly covariant** under `(φ, φ₀, φ_m) → c·(φ, φ₀, φ_m)` (SymPy-verified).

> ⚠ **SCOPE: at `DIL_COEFF = 0`.** The solver's Reynolds-dilatancy / Darcy term uses a
> `perm_factor(φ)` built on a **fixed `PHI_REF = 0.49`**, which does not rescale, so with
> dilatancy on the covariance is broken and this reduction is out of scope. Every preset
> on record uses the default `DIL_COEFF = 0`. This caveat was already written into
> `multiplet_wc.py` the same day and was not carried here until a referee found it twice.

So those are not three independent inputs. Writing `x = φ/φ₀` and `χ = φ_m/φ₀`:

&nbsp;&nbsp;&nbsp;&nbsp;**`η/η_s = (1 − x/[1 − (1−χ) f̂])⁻²`** — three numbers: **`x`, `χ`, `σ*`.**

The phase diagram is then exact and two-dimensional:

| condition | state |
|---|---|
| `x < χ` | never jams — continuous thickening only |
| **`χ ≤ x < 1`** | **shear-jamming window**, width `1 − χ` |
| `x ≥ 1` | jammed at rest — a solid, not a suspension |

> **`χ` decides *whether* a suspension can shear-jam and over how wide a concentration
> range. `x` decides *where in that range* it sits. `σ*` decides *at what stress*, and
> never whether.**

A drop rebounds as a solid iff **both** `σ_imp/σ* ≫ 1` (contacts activated) **and**
`χ ≤ x < 1` (the activated state is jammed). Two independent conditions.

---

## 2. The axis map — what each particle property actually does

| axis | knob | status |
|---|---|---|
| **material: friction `μ_p`** | **χ** | ✅ **measured**: `χ = 0.8872 ± 0.0042`, window `1 − χ = 0.113` (`μ_p` 0 → 3) |
| **distribution: CV** | **x** | ✅ **measured**: `Δχ = +0.0019 ± 0.0034` PAIRED (r = +0.86), `\|t\| = 0.56` — does **not** move the window |
| **shape: `α`** | **χ** | ✅ **measured**: `χ` 0.8856 → 0.7884, window **×1.85**, `\|t\| = 19.3`, same sign at all 4 rungs |
| **size `a`** | **σ\*** | 🔴 **degenerate**: a size ladder cannot separate the two grain-scale stresses (`onset_size_degeneracy_2026-09-19`) |
| **potential `F_rep`** | **σ\*** | 🔴 **not measured here** — it is a DEM *input*, not an output |
| **material: stiffness `E`** | — | ✅ the static packing is scale-free (`test_scale_invariance.py`) |
| *(aggregation)* | **x**, via `Λ` | ✅ derived: `Λ = ⟨k^γ⟩/⟨k⟩`; the model's own domain fails before the interesting corner |

**Two axes change the window (friction, shape) — both now MEASURED, and both work the
same way: they leave `φ₀` almost alone and pull `φ_m` down. Two change the onset stress
(size, potential) — and NEITHER is measured. One changes position (the distribution).**

> 🔑 **The bench consequence of the shape result:** between `χ(sphere) = 0.886` and
> `χ(clump) = 0.788` there is a band — at `φ₀ ≈ 0.633`, **`φ` from 0.499 to 0.560** — in
> which the same material at the same loading **shear-jams as clumps and cannot as
> spheres**. That is the sharpest prediction this framework has produced.

So the next campaign is `σ*`, not another packing study.

⚠ **`χ` is not rung-independent** — it drifts 0.016–0.025 across the pressure ladder,
*more* than `φ₀` does, because the two friction arms have different compressibilities. So
`χ` must always be quoted at a named rung. The tempting "it's a ratio, two cells measure
it" shortcut is **dead**, and it was tested rather than assumed.

---

## 3. The experimental half — and why it is easy

All three numbers are **rheometer-measurable with no new instrument** — but it is a
**concentration series, not one curve**, and that distinction matters:

- `σ*` — the onset stress, read off the thickening knee of a **single** curve.
- `φ₀`, `φ_m` — a Wyart–Cates fit to `η(σ)` **at several `φ`**. ⚠ It cannot be done from
  one curve: `η_low/η_s = (1−x)^{−n}` needs both `η_s` and the exponent `n`, and **neither
  is a datasheet number** — the calibrated `η_s` values (0.05 … 0.0039 Pa·s) are 4–50×
  any plausible carrier viscosity, so `η_s` is itself a fit output, and `n` is the
  `MP_EXPO` the solver's own notes call **degenerate with `φ`**.
- `x`, `χ` — ratios of those.

### 🔑 The falsifiable experiment: a COLLAPSE

> **Two suspensions with the same `(x, χ, σ_imp/σ*)` must impact IDENTICALLY, whatever
> their particles are made of, how big they are, or how wide their size distribution is.**

That is the strongest form of the claim and the cheapest to refute. It needs:

1. **two suspensions matched in `(x, χ, σ*)` and deliberately mismatched in particle
   size or PSD width** — the distribution result says this is *constructible*, because CV
   moves `x` without moving `χ`, so CV can be traded against concentration at fixed `χ`;
2. one drop-impact run each, at matched `We` and `Re`;
3. the spreading/rebound observable this lab already extracts.

**A collapse supports the reduction. A split refutes it**, and the direction of the split
says which axis carries information the three numbers do not — which is exactly the
result worth having either way.

⚠ **Which stress gauge `σ_imp` is remains unchosen**, and `lessons.md` records that
`ρ k U₀²` and the deviatoric wall shear give *opposite* verdicts on whether the activation
half is saturated. **The collapse test must name its gauge before it runs**, and the
matched pair should be designed so the verdict does not depend on that choice.

> 🔴 **THE COLLAPSE TEST MUST BE RUN ON THE BENCH, NOT IN THE SOLVER.** At `DIL_COEFF = 0`
> the DNS has **no particle-scale length at all** — size enters only through `σ*`, by the
> axis map's own logic — so two matched triples would collapse **trivially and
> uninformatively** there. The DNS pairs (§5 item 3) test something different: whether the
> phase diagram predicts the drop outcome. Do not conflate them.

> ⚠ **Matching `We` and `Re` is not enough when SIZE is the axis being varied.** A bench
> collapse must also control the **particle-scale** groups — Stokes number, density ratio
> `ρ_p/ρ_f`, and particle-size-to-film-thickness — or a failure to collapse cannot
> distinguish *"particle inertia / a boundary-scale effect"* from *"the three numbers are
> insufficient"*, which is the whole question. **Design the pair so those are matched too,
> or state which one is deliberately not and why.**

---

## 4. What to PRUNE

The reduction says which of this repository's ~45 pre-registrations bear on the question.

**Load-bearing — keep and finish:**
- the friction ladder (`mustarcap`) — it *is* the `χ` measurement;
- `shapechi` — the only untested window-knob;
- the WC mapping (`multiplet_wc.py`) and the covariance argument;
- the protocol scale `PHI_SCALE` and the convergence discipline — every number depends on
  them, and this session showed twice what happens when a cap is wrong.

**Not load-bearing for *this* question — finish or park, do not extend:**
- the `φ_J`-vs-PSD-descriptor programme (ternary, K-ladder, `u`-gap, skewness closure,
  additivity). It is a good result about **packing** — width alone does not determine
  `φ_J` — but in the reduction it is all **x**-knob detail, and `x` is also set by
  concentration, which an experimenter controls directly and far more precisely.
- the multiplet/aggregation work **beyond** the `Λ` derivation: its own admissibility
  bound shows the interesting corner is outside the model.
- the fabric/`Z` invariants, the strain clock, the first-passage work — none of them enter
  `(x, χ, σ*)`.

> **The blunt version: this project has been measuring `x` with great care, and `x` is the
> one number an experimentalist can already set by weighing powder.** The scientific value
> is in `χ` and `σ*`, and `χ` has exactly one measurement (friction) with a second running.

---

## 5. The four things that would close the question

In cost order, each with what it decides:

1. ~~`shapechi`~~ — **DONE**: shape IS a second window-knob (`CLAIMS_SHAPECHI`).
2. **`σ*` from the potential** — the one axis with no measurement at all. DEM with an
   explicit repulsive force, sweeping `F_rep`, checking `σ* ∝ F_rep/a²` and **which
   area convention** the simulation actually realises. This settles the size axis too,
   because size and potential enter `σ*` together and the existing derivation says a size
   ladder alone cannot separate them.
3. **The DNS pairs** — two `(x, χ, σ*)` triples run in `wc_drop_impact.c`, one inside the
   window and one outside, to check the phase diagram predicts the drop outcome.
4. **The bench collapse** (§3).

⚠ **Nothing in 1–4 needs the `N = 6000` campaign, the LevelSet control, or a third
`α` rung.** Those answer packing questions, not this one.
