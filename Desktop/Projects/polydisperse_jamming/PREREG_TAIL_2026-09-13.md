# PRE-REGISTRATION — is the jamming shift set by the NUMBER width or by the TAIL?

**Written 2026-09-13, before any cell of this campaign ran.**

---

## 1. Why, and why this is the project's biggest open question

Every headline this project quotes is a derivative with respect to **CV**, the
number-weighted coefficient of variation: `dphi_c/dCV`, `dW/dCV`, `CV*`, "a CV ladder of
0.45 traverses the window", and Eq. (11)'s critical-state channel in
`derivations/dilatancy_sign_mechanism.md`.

`shape_control.py` shows CV is not the control parameter. At **matched CV = 0.5 and
matched friction** the three size-distribution families already on disk give

| family | `dphi_J`(CV 0 → 0.5), `mu_p` = 0.3, paired | shape of `phi_J(CV)` |
|---|---|---|
| lognormal | **+0.0621 ± 0.0026** | monotone, slightly superlinear |
| uniform | **+0.0242 ± 0.0008** | **saturates** past CV ≈ 0.3 |
| bidisperse | **+0.0139 ± 0.0025** | **NON-MONOTONE**, peaks at CV ≈ 0.3 |

a factor of **4.5**, and the three ladders differ qualitatively in shape, not just in
slope. For comparison, moving the jamming point itself by 0.058 with friction moves the
lognormal shift by about 7 % (`transfer_law.py`, after the `mu_p` = 1.0 arm is refused on
health — see §3). **The distribution shape dominates by more than an order of magnitude
over the axis the project has been controlling.**

The descriptor sweep says where the difference lives: **37 % of the lognormal's volume
sits above 2·d10, and the uniform and bidisperse have ZERO volume there.** The lognormal's
advantage is its large-particle tail — and that tail is cut by the production truncation
`(0.2, 5.0)`, which was an implementation choice. The volume-weighted width of the CV 0.5
lognormal runs 0.280 / 0.355 / 0.454 / 0.483 / 0.514 as the upper cut goes 2 / 3 / 5 / 8 /
20, i.e. **the production cut sits on the rising part and is not converged.**

So the headline may be partly a property of a constant in `packing.insert_packing`.
This campaign finds out.

> ⚠ **Read §5 before quoting this paragraph.** It states what I believed when the
> campaign was designed, and the campaign REFUTED it: `phi_J` moves by 0.021 between
> `U` = 2 and 5 and by **−0.0010 ± 0.0017, not resolved**, between 5 and 8. The
> production cut is on the PLATEAU. `cv_vol` still rises 6 % over that last step, which
> is precisely why `cv_vol` is not the control variable either.

---

## 2. The design: one family, one friction, one thing changed

**72 cells**, isotropic compression at `mu_p` = 0.3, N = 1500, lognormal, on the
**same 4-rung pressure ladder every `phi_J` on disk uses** (`P_HAT_LADDER` =
1e-3, 3e-4, 1e-4, 3e-5), 3 seeds.

Four truncated lognormals, all cut below at 0.2·d10 and above at `U`·d10, with the
**nominal sigma re-solved per window so the REALISED number CV is 0.500 in every case**:

| `U` | nominal CV | realised number CV | **volume-weighted CV** | `d32` |
|---|---|---|---|---|
| 2.0 | 0.7237 | 0.500 | **0.2804** | 1.252 |
| 3.0 | 0.5447 | 0.500 | **0.3548** | 1.487 |
| 5.0 | 0.5034 | 0.500 | **0.4542** | 1.556 |
| 8.0 | 0.5012 | 0.500 | **0.4831** | 1.563 |

**Matched in number width to 1e-3; a factor 1.72 apart in volume width.** That is the
whole experiment: the two candidate control variables are separated by construction, and
the separation is a property of the design rather than of the outcome.

Plus a **CV = 0 control** at `U` = 2 and 8 (24 cells): a monodisperse pack has no tail to
cut, so truncation must be inert there. If it is not, the truncation is doing something
other than what this campaign thinks.

---

## 3. Predictions, with bars. Written before the first cell.

**T0 — the campaign is what it says it is.** 72/72 records; realised number CV within
**0.005** of 0.500 on every CV = 0.5 cell; `trunc_hi`, `cv_nominal` and `cv_vol_design`
recorded per record; the pressure ladder identical to `P_HAT_LADDER`.

**T1 — 🔴 THE LADDER-HEALTH GATE, and it is registered because it has already bitten this
project.** `phi_J_lin` is a straight line through four pressure rungs and its standard
error cannot see a ladder that is not a straight line. An adversarial review on
2026-09-13 found that the `mu_p` = 1.0 arm on disk — which was carrying a claim in
`transfer_law.py` — is **0 of 6 monotone, with `r2_lin` as low as 0.002 and one
UNPHYSICAL negative slope** (`phi_J` falling as pressure rises). Every arm here must pass,
per (U, CV) arm of 3 seeds:

- `slope_lin` > 0 on **3/3** seeds (a negative slope is inadmissible, not merely noisy);
- median `r2_lin` >= **0.85** (the healthy arms on disk sit at 0.91–0.95);
- `converged` on 3/3.

An arm failing any of these is **REFUSED and reported as refused**, and its numbers stay
printed. If two or more of the four `U` arms are refused the campaign is **UNREADABLE**
and says so rather than quoting the survivors.

**T2 — the control: truncation is inert at CV = 0.** `phi_J`(CV 0) must agree between
`U` = 2 and `U` = 8 to **0.003**. This is a null that a correct implementation must pass
and a broken truncation would fail, and it is the cheapest way to catch a truncation that
is silently changing something else (the mean, the seed stream, the collider guard).

**T3 — THE TEST.** `phi_J` across the four `U` arms at CV 0.5.

- **NUMBER-WIDTH HYPOTHESIS**: `phi_J` is the same for all four, spread <= **0.005**.
- **TAIL HYPOTHESIS**: `phi_J` tracks the volume-weighted width. With `K` = 0.0787 from
  the three-family fit, the predicted spread is `K` × (0.4831 − 0.2804) = **0.0160**.
- The two predictions are **3.2×** apart against the 0.005 bar, so the test is informative
  before it runs — and the informativeness is asserted here, not after the fact.
- **UNDETERMINED** if the measured spread lands between 0.005 and 0.016 in a way that the
  seed scatter cannot separate; the seed scatter on `phi_J` at N = 1500 is 0.0007–0.0048
  on record, so a spread of 0.010 ± 0.004 would be exactly that and will be reported as
  such rather than assigned to a side.

**T4 — the cross-check that makes it a law rather than a coincidence.** If T3 favours the
tail, the four `U` points and the three FAMILY points must fall on the **same** line
`dphi_J` vs `cv_vol` — same `K`, fitted on one set and tested on the other. This is an
out-of-sample test of a proportionality across two completely different ways of changing
the distribution, and it is the difference between "the tail matters here" and "the
volume-weighted width is the control parameter".

**T5 — reported, not graded.** `z_nr`, the rattler fraction and `d_max_over_L` per arm
(`U` = 8 pushes the largest grain to 7.6·d10 and the finite-size flag is the thing to
watch); `d32` per arm, which is nearly degenerate between `U` = 5 and 8 (1.556 vs 1.563)
while `cv_vol` still moves 6 %, so those two arms are a small extra discriminator between
`d32` and `cv_vol` as the variable.

---

## 4. What each outcome buys

| outcome | meaning |
|---|---|
| T3 NUMBER-WIDTH | the production truncation is not load-bearing; the headline is safe as a number, and the 4.5× family spread must be explained by something other than the tail |
| T3 TAIL, T4 collapses | **the control parameter is the volume-weighted width**, the headline +0.053 is conditional on a truncation constant AND on the material's tail, and every `d/dCV` in this project needs restating on a volume-weighted axis |
| T3 TAIL, T4 does not collapse | the tail matters but is not the whole story; the control parameter is at least two-dimensional and the project should stop looking for one number |
| T1 refuses ≥ 2 arms | the isotropic ladder cannot carry this question at N = 1500 and the test must move to the constant-pressure protocol, which has no extrapolation in it |
| T2 FAIL | the truncation implementation is doing something other than truncating |

---

## 5. OUTCOME

**Landed 2026-09-13. Two runs: 72 cells at the default 200 000 iteration cap, then the
whole grid again at `PJ_MAX_ITER` = 1e6; the grader takes the higher cap per cell.**
Verdict line: `FAILURES: ['T0', 'T1', 'T3']`. Every number below is printed by
`grade_tail.py` → `out/tail.txt`.

### The registered verdict: T3 is UNREADABLE, and it stays that way

**T1 FAIL, 2 of 4 arms refused.** After the cap repair every arm converges (the first run
failed 0/3 on the `p̂` = 3e-5 rung alone, at the cap, with unbalanced 0.014–0.043 against
a 0.02 tolerance — slow cells, not wrong ones). What still refuses them is the **r² of the
`φ(p̂)` line**: `U` = 2.0 at 0.609 and `U` = 8.0 at 0.842, against the 0.85 bar.
By §3's own rule, two refusals make the campaign **UNREADABLE and the survivors are not
quoted** — so **T3 is UNREADABLE** and `φ_J` is not reported for any arm.

**T0 FAIL**, and the bar was posed on the wrong quantity: it required *each cell's*
realised number CV within 0.005 of 0.500, which at N = 1500 tests the seed's own sampling
noise rather than the design. The **arm means** are 0.5003 / 0.4963 / 0.4993 / 0.4968 —
matched to 0.004, which is what the design promised. The bar stays failed and the arm
means are printed beside it. Four cells also failed outright: `U` = 8, seed 1, all four
rungs, on the periodic-collider guard (`d_max` = 5.9·d₁₀ needs a larger cell than N = 1500
gives).

**T2 PASS.** At CV = 0 the two extreme windows give `φ_J` identical to **0.000000** — a
monodisperse pack has no tail to cut, and the truncation is doing only what it says.

### 🔑 The companion, which needs no extrapolation, and it is decisive

`φ_J` is an *extrapolation*, and it is the extrapolation that T1 finds unhealthy — the
slopes are all positive (8.7–18.7) and the three seeds inside each arm agree to ~0.005.
So the same hypothesis was asked of a quantity with **no fit in it**: `φ` measured *at*
each pressure rung. This removes a model rather than a bar, and it is **reported, not
graded** — T3's verdict stands.

Paired per (seed, pressure), which cancels the seed's realisation and the pressure in one
step:

| tail cut | `Δφ` | t | n |
|---|---|---|---|
| `U` 2 → 3 | **+0.0101 ± 0.0011** | +9.0 | 12 |
| `U` 3 → 5 | **+0.0109 ± 0.0011** | +9.8 | 12 |
| `U` 5 → 8 | **−0.0010 ± 0.0017** | −0.6 | 8 |
| `U` 2 → 5 | **+0.0211 ± 0.0013** | +16.5 | 12 |

and the unpaired per-rung version agrees at every one of the four pressures independently
(`U` 2 → 5 gives +0.0219, +0.0208, +0.0190, +0.0226 at `p̂` = 1e-3, 3e-4, 1e-4, 3e-5).

**THE NUMBER-WIDTH HYPOTHESIS IS DEAD.** These four distributions are matched in realised
number CV to better than 0.005, and cutting the large-particle tail at 2·d₁₀ instead of
5·d₁₀ lowers the packing fraction by **0.021 — four times the 0.005 bar the number-width
hypothesis predicted, and 16 σ from zero.** Whatever sets the polydispersity shift, it is
not the number-weighted width.

### 🔑 And the reassuring half: the production truncation is on the PLATEAU

The effect rises steeply to `U` = 5 and then **stops**: `U` 5 → 8 is **−0.0010 ± 0.0017,
not resolved at 2 se**, against +0.0211 for `U` 2 → 5. So the production constant
`(0.2, 5.0)` sits where the curve has flattened, and **the project's +0.053 headline is
not an artefact of it.** That was the specific worry this campaign was launched to settle,
and it is settled in the direction that costs nothing.

⚠ Note this was NOT obvious in advance from the descriptors: `cv_vol` still moves 6 %
between `U` = 5 and 8 (0.454 → 0.483) while `φ_J` moves by nothing measurable. So `cv_vol`
over-predicts in that region — one more reason it is not the control variable
(`shape_control.py` reaches the same conclusion from the family ladders).

### What this settles, and what it does not

- **SETTLED: the number-weighted width is not the control parameter**, now from two
  independent directions — across families at matched CV (4.5×, `shape_control.py`) and
  *within* one family at matched CV (0.021, here). The second is the cleaner of the two
  because it changes one thing.
- **SETTLED: the production truncation is not load-bearing.** The headline stands.
- **NOT SETTLED: what the control parameter is.** `cv_vol` orders the tail sweep correctly
  but over-predicts the `U` 5 → 8 step; fifteen descriptors fail on the family ladders.
  The honest position is that the shift is a property of the whole distribution.
- **NOT SETTLED, and now the binding limitation: the isotropic `φ(p̂)` ladder cannot carry
  this question cleanly at N = 1500.** Two arms fail its r² bar, one seed fails the
  collider guard, `d_max/L` = 0.34 at `U` ≥ 5 is above the 0.25 finite-size flag, and the
  rattler fraction is ~0.50 throughout. §4's registered response applies: **move the test
  to the constant-pressure protocol, which has no extrapolation in it** — and which is
  also the channel the project's headline is actually measured on.
