# Polydispersity translates the Wyart-Cates φ-axis; it does not change the frictional window

> ## 🔴 CORRECTED THE SAME DAY. Read `wc_friction_endpoint_2026-09-19.md` with this page.
>
> Two independent errors below were **cancelling**, and correcting either alone would have
> changed the conclusion. Both are arithmetic on data already committed, verified before
> being written down.
>
> 1. **`φ_J(μ_p = 0.3)` is NOT Wyart-Cates's `φ_m`.** The `μ_p = 1.0` and `μ_p = 3.0` blocks
>    have been on disk since before this session and no analysis had read them.
>    `μ_p = 0.3` captures only **77.2 %** of the `μ_p = 0 → 3.0` drop; the remaining
>    **22.8 %** is resolved at **\|t\| = 14.29**. `μ_p = 1.0` is converged (\|t\| = 2.78).
>    Everything below that says "`φ_m`" means "`φ_J` at `μ_p = 0.3`".
> 2. **The flatness test below used the wrong yardstick**, as referee pass 3 suspected. A
>    RANGE over six rungs against twice a TWO-POINT propagated s.e. does not cancel the
>    shared seed draw. On the per-seed **trend**, the `μ_p = 0.3` window is
>    `dW/dCV = −0.01217 ± 0.00155`, **\|t\| = 7.87 — RESOLVED, it is not flat** (10.5 % of
>    `W(0)` over the ladder). **The "range 0.00646 = 1.75 s.e." line in the table below is
>    therefore withdrawn as a test**, though the numbers in it are correct.
> 3. **With the converged endpoint the conclusion returns**: `dW/dCV = +0.00398 ± 0.00178`,
>    \|t\| = 2.24, not resolved — flat, at 2.7 % over the ladder.
> 4. **And the mechanism of the apparent deformation is identified**: the residual friction
>    gap above `μ_p = 0.3` **grows with CV** (+0.01550 → +0.02136 from CV 0 to 0.5, growth
>    +0.00586 ± 0.00077, \|t\| = 7.65). A wider size distribution needs more friction to
>    reach its frictional limit.
>
> **What survives unchanged, and it is the load-bearing claim: the CV shift is
> friction-independent** across `μ_p` = 0 / 0.3 / 1.0 (worst pairwise \|t\| = 1.40 against
> 4.303), while `φ_J` itself moves 13× the spread in the shift over the same range.
>
> ⚠ Everything below is left as written, with these corrections at the top rather than
> edited in, so the two versions can be compared.

**2026-09-19.** Every number is printed by `wc_endpoints_move_together_2026-09-19.py` into
`wc_endpoints_move_together_2026-09-19.txt` (rc 0, ALL CHECKS PASSED). No new cells: this
runs entirely on `out/phase2_fits.csv`, which has held both friction endpoints at matched CV
since before this session.

> ## RESULT
>
> In Wyart-Cates, `φ_J(σ)` interpolates between `φ_0` (frictionless) and `φ_m` (frictional).
> **Both endpoints exist on disk at matched CV for two of the three families**, on the full
> six-rung CV ladder with three seeds each. Their CV dependences are not resolved apart:
>
> | family | `Δφ_0` (CV 0→0.5) | `Δφ_m` | difference | s.e. | \|t\| | vs `t_crit` = 4.303 |
> |---|---|---|---|---|---|---|
> | lognormal | +0.05803 ± 0.00268 | +0.06207 ± 0.00255 | −0.00404 | 0.00370 | 1.09 | not resolved apart |
> | uniform | +0.02385 ± 0.00162 | +0.02420 ± 0.00080 | −0.00035 | 0.00181 | 0.20 | not resolved apart |
>
> and the frictional window `W = φ_0 − φ_m` is flat across the whole ladder:
>
> | family | CV 0.0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | range |
> |---|---|---|---|---|---|---|---|
> | lognormal | 0.05784 | 0.05938 | 0.05734 | 0.05433 | 0.05292 | 0.05380 | 0.00646 (1.75 s.e.) |
> | uniform | 0.05723 | 0.05703 | 0.05756 | 0.05649 | 0.05840 | 0.05688 | 0.00191 (1.06 s.e.) |
>
> **Reading: polydispersity TRANSLATES the Wyart-Cates φ-axis rather than deforming it.**
> The gap the model interpolates across is the same width at CV = 0.5 as at CV = 0, while
> both its ends move up together by `Δφ_m`.
>
> ⚠ **Verdict: SUPPORTED AS A BOUND, not as an equality, and NOT pre-registered.** At three
> seeds the design could only have resolved a difference larger than **0.0159 (lognormal,
> 25.7 % of the shift)** or **0.0078 (uniform, 32.2 %)**. The supportable sentence is *"the
> two endpoints move together to within a quarter to a third of the shift"*, never *"they
> move identically"*.

---

## 1. Why this matters for the drop

Two consequences follow, and neither needs the size-exponent machinery that
`onset_size_degeneracy_2026-09-19.md` found undetermined.

**(a) The shape of the DST diagram in `φ` is CV-independent.** Wyart-Cates writes
`φ_J(σ) = φ_0 f(σ/σ*) + φ_m [1 − f(σ/σ*)]` with `f` the frictionless fraction. If `φ_0` and
`φ_m` both shift by the same `δ`, then `φ_J(σ) → φ_J(σ) + δ` for every `σ`. The curve
translates. Nothing about its steepness, its stress dependence, or the width of the
thickening transition changes.

**(b) `σ*` is untouched, so the whole CV effect lands on the `φ` axis.** `σ*` is set by the
repulsive force scale between two grains, and Wyart-Cates makes it independent of `φ`.
Nothing in this measurement bears on it. So the prediction is clean: **widening the size
distribution at fixed `φ` moves the suspension down in `φ/φ_m` without moving the stress at
which thickening begins.** The packing fraction at which DST appears rises by `Δφ_m`, which
is 0.124 per unit CV for a lognormal family and 0.048 for a uniform one.

That is the same sign as the author's bench observation of 2026-08-31 (raising CV at similar
particle size moves the drop from dilatant toward spreading), and it reaches it without any
of the contested conventions of the companion page.

## 2. Checks that could have failed

- **`φ_0 > φ_m` at every (family, CV)** — friction lowers the jamming fraction, which is the
  sign the literature reports. 0 violations of 12.
- **The window is a real gap, not seed noise** — `W = 0.05784` at CV = 0, which is **21.6×**
  the seed s.e. of the shift.
- **Both shifts are positive at both frictions** — 4 of 4.
- **The grid is printed, not assumed.** `μ_p ∈ {0, 0.3, 1.0, 3.0}` are all present, but
  `μ_p = 1.0` has three CV levels and one family and `μ_p = 3.0` has one cell, so neither
  enters. The check asserts that bidisperse is *absent* from the frictionless block, so it
  fails if a future campaign adds one and the "two of three families" wording goes stale.

## 3. Caveats that travel with every number

1. 🔴 **NOT PRE-REGISTERED.** The session brief asked for the two candidate outcomes and the
   resolving power to be written down *before* computing. They were not; the endpoints were
   found on disk and the contrast computed in the same step. This is a post-hoc comparison
   and the resolving power in the result box is stated after the fact.
2. ⚠ **`paper/main.tex` already caveats `φ_0`**: *"The frictionless ladders are not monotone
   in `p̂` at this seed count, so this anchor is a check of the protocol's scale and not a
   result of the paper."* That caveat attaches to every `φ_0` number here.
3. ⚠ **Two families, not three — and the missing one is the awkward one.** Bidisperse, the
   non-monotone family, has no frictionless block. The family least like the others is the
   one not tested.
4. ⚠ **At CV = 0 the generator ignores the family label** (found 2026-09-19), so the two
   families' `φ_0` values share a cell. A test treating the families as independent
   over-states the variance of their difference, i.e. errs conservatively.
5. ⚠ **`μ_p = 0.3` is one friction, not the frictional limit.** `φ_m` in Wyart-Cates is the
   `μ_p → ∞` endpoint; 0.3 is a finite value the campaign chose. The `μ_p = 1.0` block has
   three CV levels and could bound how much that matters — **it was not used here**, and
   using it is the obvious next step.

## 4. What would change the verdict

- **More seeds.** The bound is 26–32 % of the effect at n = 3; six seeds would roughly halve
  it, and the cells are cheap at `μ_p = 0`.
- **A frictionless bidisperse block**, which would test the family the comparison currently
  skips and costs six cells.
- **Using the `μ_p = 1.0` ladder** to check whether `φ_m` at 0.3 stands in for the frictional
  limit, or whether the window keeps widening with friction.
- A *resolved* difference between `Δφ_0` and `Δφ_m` would kill the translation reading and
  say that polydispersity changes the shape of the DST diagram, not just its position.

## 5. Referee

⚠ **Unrefereed.** This page was written in the last hour of the session, after the companion
page had two headlines withdrawn by two referee passes. Its verdict is a bound rather than a
claim, which is the safer failure mode, but the "translation" reading in §1(a) is an
interpretation of the numbers and has not been attacked.
