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

is **exactly covariant** under `(φ, φ₀, φ_m) → c·(φ, φ₀, φ_m)` (SymPy-verified). So those
are not three independent inputs. Writing `x = φ/φ₀` and `χ = φ_m/φ₀`:

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
| **distribution: CV** | **x** | ✅ **measured**: `Δχ = +0.0019 ± 0.0049`, `\|t\| = 0.40` — does **not** move the window |
| **shape: `α`** | **χ?** | ⏳ **running** (`PREREG_SHAPECHI`) — the only untested window-knob |
| **size `a`** | **σ\*** | 🔴 **degenerate**: a size ladder cannot separate the two grain-scale stresses (`onset_size_degeneracy_2026-09-19`) |
| **potential `F_rep`** | **σ\*** | 🔴 **not measured here** — it is a DEM *input*, not an output |
| **material: stiffness `E`** | — | ✅ the static packing is scale-free (`test_scale_invariance.py`) |
| *(aggregation)* | **x**, via `Λ` | ✅ derived: `Λ = ⟨k^γ⟩/⟨k⟩`; the model's own domain fails before the interesting corner |

**Two axes can change the window (shape, friction). Two change the onset stress (size,
potential). The rest change position.** That is the whole answer at framework level, and
it is why the next campaign is the shape one and not another distribution one.

⚠ **`χ` is not rung-independent** — it drifts 0.016–0.025 across the pressure ladder,
*more* than `φ₀` does, because the two friction arms have different compressibilities. So
`χ` must always be quoted at a named rung. The tempting "it's a ratio, two cells measure
it" shortcut is **dead**, and it was tested rather than assumed.

---

## 3. The experimental half — and why it is easy

All three numbers are **rheometer-measurable from one steady-state flow curve**:

- `σ*` — the onset stress, read directly off the thickening knee.
- `φ₀`, `φ_m` — a Wyart–Cates fit to `η(σ)` at several `φ` returns both branches.
- `x`, `χ` — ratios of those.

So the bench half needs **no new instrument**, and the DEM half predicts how `x` and `χ`
move when you change the particles.

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

1. **`shapechi`** (running, 72 cells) — is shape a second window-knob?
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
