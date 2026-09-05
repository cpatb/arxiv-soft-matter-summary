# The grid study, and a referee pass that took most of it back

**2026-09-03.** Written to close the item the 09-02 venue addendum opened — ASME JFE's
numerical-accuracy policy needs a discretisation-error estimate and the paper had two
grids and a percentage difference. Four grids were run. **A physics-referee pass then
falsified the study's central explanation and undermined its headline, and the corrected
version is below.** The first draft of this document claimed `C_P` was converged; that
claim does not survive and is retracted here rather than softened.

Artifacts: `cfd_openfoam/gci.py`, `convergence/gci_Re10000_H6_*.md`, `yplus_profile.py`,
`closure_sensitivity.py`, `closure_sensitivity.json`.

---

## 1. What survives, and it is worth having

**The iterative-convergence finding stands, and it is the most solid thing here.** Of the
grids in the original study only the coarsest two met `residualControl`; the paper's fine
grid stopped at the 3000-iteration cap. Restarted (`extend_case.sh`) it converges at
iteration **5105**, and `C_P` moves **1.5178 → 1.3928, −8.2 %**.

> **The manuscript's *"`C_P` changes 12.28 % between grids"* is really 3.0 %. Two thirds of
> the headline grid sensitivity was an unconverged solve.** That has to come out of the text
> regardless of everything below.

And the direction of that discovery is worth carrying: **correcting the iterative error made
the convergence study look worse.** At the capped state the same three grids gave a tidy
`p` = 1.70 and `GCI` = 6.06 % on `C_P`. **A clean-looking GCI is not evidence that the grids
converged.**

## 2. What the four grids give

Family `fine_conv`: refine **4 / 2 / 1**, `r` = 2, cells 491 520 / 122 880 / 30 720, all
meeting `residualControl` at iterations **11 552 / 5 105 / 2 261**.

| quantity | refine 4 | refine 2 | refine 1 | verdict |
|---|---|---|---|---|
| `n_t` | 1.9338 | 1.8922 | 1.9023 | **non-monotone** |
| **`C_P`** | **1.5809** | **1.3928** | **1.3518** | 🔴 **DIVERGING**, \|ε21/ε32\| = 4.59 |
| `τ_w`(15 mm) | 6.1374 | 6.2415 | 5.8658 | **non-monotone** |
| peak `τ_w` | 486.30 | 441.17 | 389.52 | `p` = 0.195, outside [0.5, 2.5] |

**Read as raw numbers this looks reassuring**: three of the four move by only 1.7–2.2 %
across a 16× cell-count range. Read as a convergence study it is not, for five reasons —
and the first one was a bug in my own tooling.

## 3. 🔴 Five findings, each verified here independently

### 3z. `C_P` is DIVERGING under refinement, and my GCI script reported |p| for it

**This is the biggest one and it came from a code review, not the referee.**
`observed_order` implemented Celik et al. Eq. (2) including its `abs()`, and guarded only
the SIGN of `ε32/ε21` — never the magnitude. Celik classify a triplet by
`R = ε21/ε32`: `0 < R < 1` is monotonic convergence, `−1 < R < 0` oscillatory, and
**`|R| ≥ 1` is DIVERGENCE, which the procedure does not cover.**

For `C_P` on this family, `ε21` = −0.1881 and `ε32` = −0.0410, so **`|R|` = 4.59**: the gap
between the two *finest* grids is four and a half times the gap between the two coarsest.
Refining takes `C_P` 1.3518 → 1.3928 → 1.5809 **with the steps growing**. The sequence is
running away, and the `p = 2.198` that was printed is `|−2.198|`.

**Both rows in the whole study that carried an "OK" verdict were the divergent ones.** So
`GCI_fine` = 4.14 % and the extrapolation to 1.6333 came from a diverging sequence and are
**void**. `closure_sensitivity.py`'s extrapolation row is deleted.

It also *explains* §3c rather than leaving it a puzzle: `p` = 1.703 / 0.315 / 2.198 across
the three families is |−1.703| / **+0.315** / |−2.198| — the two extremes are sign-flipped
divergences and the middle one is the only genuine observed order in the study.

`gci.py` now refuses both cases and names which, and asserts the downstream invariant that
would have caught it independently: with a true negative `p`, `r21^p − 1 < 0`, so the GCI
would come out **negative**. A bound that cannot be negative only because an `abs()`
upstream removed the sign is not a bound.

The full picture after the fix, twelve rows across three families and four quantities.
⚠ **Corrected 2026-09-04 after an independent verification pass — this paragraph said
"exactly one row is a genuine converging order" and there are TWO**, and the second was
printed in this document's own §2 table two pages above. `gci.py` returns a computed order
on three rows:

| family | quantity | `p` | in band? | iteratively sound? |
|---|---|---|---|---|
| `coarse_conv` | `C_P` | 0.315 | no | **yes** — all three grids converged |
| `fine_conv` | peak `τ_w` | 0.195 | no | **yes** — all three grids converged |
| `coarse` | `n_t` | 0.657 | yes | **no** — its finest grid is capped at 3000 it. |

So **two rows are genuine converging orders and both are outside the asymptotic band**,
and the one row that lands inside the band is disqualified by iterative non-convergence.
The miscount did not flatter the study — `fine_conv`'s peak `τ_w` carries `GCI` = 80.2%
and `e_ext21` = 39.1%, so adding it makes the grid-convergence picture *worse*, not
better. (The disqualified row's grid is the **finest** of its own three-grid triplet, not
its "middle" one; it is only middle relative to the four-grid superset.)

### 3a. The y⁺ mechanism is refuted — it has the wrong sign in the right place

The first draft explained the coarse family's failure by `y⁺` straddling 1 (medians 1.52 /
0.78 / 0.41), since `omegaWallFunction` blends the viscous and log layers there and the
quantity being converged is the wall shear stress. **The median is the wrong statistic and
the mechanism does not survive resolving it radially** (`yplus_profile.py`):

| case | median | r/d 2–5 | 5–10 | 10–14 | **14–20** | 20–25 |
|---|---|---|---|---|---|---|
| refine 4 | 0.211 | 0.206 | 0.100 | 0.058 | **0.038** | 0.025 |
| refine 2 | 0.406 | 0.398 | 0.190 | 0.113 | **0.079** | 0.058 |
| refine 1 | 0.781 | 0.776 | 0.372 | 0.219 | **0.152** | 0.112 |
| refine 0.5 | 1.517 | 1.517 | 0.733 | 0.436 | **0.299** | 0.219 |

**The `y⁺` ≈ 1.5 "straddle" is entirely a near-stagnation effect (r/d 2–5). At r/d = 15 —
where `τ_w`(15 mm) lives, where the non-monotonicity is, and where the design criterion
evaluates — `y⁺` is 0.04 to 0.30, well under 1 on *every* grid including the coarsest.**
Wall-function blending therefore cannot be the mechanism for the non-monotonicity it was
invoked to explain.

**Consequence for §2.** The change from `coarse_conv` (`p` = 0.315) to `fine_conv`
(`p` = 2.198) is then **a change of triplet, not the repair of a mechanism.** The first
draft presented it as the latter. It is not.

### 3b. `p` = 2.198 exceeds the formal order of the discretisation

`system/fvSchemes` specifies `div(phi,U) bounded Gauss linearUpwind grad(U)` (≈2nd order)
but **`div(phi,k)` and `div(phi,omega)` are `bounded Gauss upwind` — first order.** `nut`,
and hence the wall shear through the wall function, depends on `k` and `omega`. An observed
order **above the formal order of a first-order link in its own causal chain** is a signature
of accidental cancellation on one triplet, not of second-order accuracy.

### 3c. `p` is unstable across overlapping triplets from the same mesh family

| family | grids | as printed before | **corrected** |
|---|---|---|---|
| `coarse` | 2 / 1 / 0.5, refine-2 **capped** | 1.703 "OK" | **diverging**, \|R\| = 3.26 |
| `coarse_conv` | 2 / 1 / 0.5, refine-2 converged | 0.315 outside band | **`p` = 0.315**, the only genuine one |
| `fine_conv` | 4 / 2 / 1 | 2.198 "OK" | **diverging**, \|R\| = 4.59 |

A factor of seven from the same handful of grids — and §3z says why: two of the three were
sign-flipped divergences. **Cross-triplet comparison was the check that exposed it**, which
is the reason to run more than one triplet even when the first looks clean. A three-point
Richardson fit has zero residual degrees of freedom and cannot self-diagnose.

### 3d. 🔴 The designed independent cross-check cannot run, and its grid is genuinely stuck

`gci.py --family midfine` (refine 2 / 1.414 / 1, `r` = √2) exists precisely so the observed
order can be measured from a second, disjoint triplet. The refine = 1.414 grid **ran 20 000
iterations and did not converge**: every field reaches ~10⁻⁹ **except `Uz`, which sits at
1.8 × 10⁻²** — seven orders above the others and three above the target. `Uz` is the
out-of-plane component on a one-cell-thick wedge and should be ~0; a limit cycle there is a
defect of that case, not slowness. (`Uz` is 1.6 × 10⁻⁶ to 4.9 × 10⁻⁶ on the three grids that
did converge, so this is specific to refine = 1.414.)

**So the second, disjoint triplet is still outstanding, and the grid it needs is broken
rather than merely slow.** §3z has since answered the narrower question — `p` = 2.198 was an
artefact — but an independent order estimate is still the thing that would tell us whether
*any* order exists on this mesh family.

## 4. 🔴 And the design propagation extrapolates outside the fit band

This is the finding with the largest consequence and it is independent of the grid study.

`bidet_model.optimum()` solves for the nozzle diameter `d*`, then evaluates the closure at
the fixed cleaning radius `r_max` = 15 mm. Because `d*` shrinks as `C_P` grows:

| `C_P` from | `C_P` | `d*` [mm] | **`r_max/d*`** | in fit band [2, 20]? | `Q*` [L/min] |
|---|---|---|---|---|---|
| base grid — what the paper quotes | 1.3518 | 0.6535 | **22.95** | **NO** | 0.1432 |
| refine 2, converged | 1.3928 | 0.6397 | **23.45** | **NO** | 0.1372 |
| refine 4, converged | 1.5809 | 0.5844 | **25.67** | **NO** | 0.1145 |
| Richardson extrapolation | 1.6333 | 0.5710 | **26.27** | **NO** | 0.1093 |

**Every design point sits outside the band the closure was fitted over**, and the CFD domain
itself ends at `Rdom` = 25 mm — `r/d` = 25 at the simulated `d` = 1 mm — so the two finest
cases evaluate the power law **at or beyond the edge of the computation that produced it**.
`gen_case.py`'s own comment already says the `τ_c` = 1 Pa contour lies outside the domain.

Worse, the extrapolation is **not uniform across the comparison**: a larger `C_P` gives a
smaller `d*` and therefore a *worse* extrapolation, so the "corrected" rows are the ones
standing on the thinnest ice.

> **Therefore the 23.6 % fall in `Q*` is not a result this study can quote.** The arithmetic
> is right — `Q* ∝ C_P^{−2/(n_t−½)} = C_P^{−1.426}`, verified independently, an amplifying
> exponent — but the inputs are a closure evaluated 15–30 % past its fitted range at a nozzle
> diameter 40 % smaller than the one ever simulated, and that extrapolation error has never
> been estimated. Quote the direction, not the number.

⚠ Note that `run_bidet.py`'s own `[Result 1]` avoids this: it evaluates at `d` = 1 mm, i.e.
`r/d` = 15, in-band. The discipline exists in the codebase; it is the `optimum()` call that
escapes it.

### 4a. And `τ_w`(15 mm) is exactly the quantity that does not converge

`r_max` = 15 mm is literally the `tau_at_15mm_Pa` row of §2 — the one flagged non-monotone.
So the design criterion's own evaluation point is the location with no demonstrable
convergence, while `C_P`, the quantity that does have a `p`, is a **regression statistic**:
an intercept at `r/d` = 1, extrapolated past the *low* edge of the fit band, whose
statistical weight sits in the inner, best-resolved part of the profile. Convergence of that
intercept is not evidence that the field is converged where the design uses it.

## 5. What the paper should now say

1. 🔴 **Replace "`C_P` within 12 %"** — between converged grids it is 3.0 %, and see 2.
2. 🔴 **State that `C_P` DIVERGES under refinement.** Not "unconverged" — the steps grow:
   1.3518 → 1.3928 → 1.5809, `|ε21/ε32|` = 4.59. There is no extrapolated value and no
   GCI-backed error bar, and the raw spread (11.9 % between the two finest) is a **lower
   bound** on the discretisation error rather than an estimate of it, because the sequence
   is not bracketing anything.
3. 🔴 **Quote `n_t` and `τ_w`(15 mm) as bands** and say they are non-monotone in the grid.
   **Remove any peak `τ_w`** (GCI 80 %).
4. 🔴 **State that the design optimum evaluates the closure outside its fit band**
   (`r/d` ≈ 23–26 against a fitted [2, 20], and past the domain edge at 25). This is a
   pre-existing modelling gap the grid study happened to expose; it is not caused by it.
5. ✅ **Keep the sensitivity exponent** `Q* ∝ C_P^{−2/(n_t−½)}` — it is exact, it amplifies,
   and it states how much closure accuracy a design model of this shape needs. Quote it as
   a sensitivity, not as a corrected `Q*`.
6. ⚠ **The JFE numerical-accuracy policy is NOT satisfied, and the reason is now specific
   rather than vague.** The study has bounded the *iterative* error; it has shown the
   discretisation spread is a few per cent on three of four quantities; and of twelve
   (family × quantity) rows, **exactly two are genuine converging orders and both sit
   outside the asymptotic band** — `coarse_conv`'s `C_P` at `p` = 0.315 and `fine_conv`'s
   peak `τ_w` at `p` = 0.195, the latter carrying `GCI` = 80.2%. `C_P` on every family
   containing the refine-4 grid **diverges**. (Corrected 2026-09-04; this said "exactly
   one" and undercounted, in the direction that made the study look better than it is.)

## 6. What would settle it, in order of cost

1. **Fix or replace the refine = 1.414 grid** so `midfine` can run. Its `Uz` limit cycle is
   the cheapest open question in the study and it gates the only independent order estimate.
   Suspect the odd cell counts (`nx` = 169) against the collapsed wedge axis.
2. **Extend the domain to `Rdom` ≈ 50 mm and refit.** That removes the extrapolation in §4
   entirely rather than bounding it, and the 08 AUG card already identified it as the thing
   that would settle the cleaning-radius claim. It is the single highest-value CFD run left
   in this project.
3. **Second-order schemes for `k` and `omega`**, then repeat. Until then no observed order
   above 1 is interpretable.
4. **A second `Re_d`.** The `Re_d^{−1/2}` in the closure has never been tested against
   anything; the design runs at `Re_d` ≈ 7 400–8 400 against a single simulated 10 000.

## 7. Independent verification, and a process note

**Every number in this document has now been independently re-run and reproduced.** The
numerics verifier executed each claim itself rather than reading mine:

- The Celik arithmetic was **hand-recomputed without calling `gci.py`** — `p` = 2.19780,
  `GCI_fine` = 4.1454 %, `φ_ext` = 1.633328, agreeing to four significant figures.
- The three convergence iteration counts (**11 552 / 5 105 / 2 261**) were grepped from
  `log.simpleFoam` directly, as was `r1_ext`'s failure: 0 hits on "SIMPLE solution
  converged", `Uz` initial residual **1.85e-2** against ~1e-9 on every other field and
  3000–11 000× the `Uz` of the three grids that did converge.
- The `y⁺` bins reproduced exactly, and the binning was judged a fair reading of the same
  raw sample `fit_closure.py` uses (median per bin, true radial coordinate).
- `Q* ∝ C_P^{−1.426}` was recovered independently by calling `optimum()` at two `C_P`
  values and taking the log-ratio: **−1.4262**.
- The Birge inflation was recomputed and the intervals do overlap.

⚠ **One methodological catch worth keeping.** `latexmk`'s *intermediate* pass reports a
transient 3.9095 pt overfull box that disappears once `bibtex` resolves the references.
**Only the final `.log` is quotable** — which is the same shape as this document's own
headline finding, that an intermediate state can look like a result.

Two reviewers (code review, and the frisbee physics referee) were still outstanding when
this was written; three of the original five had died on a session rate limit and were
re-dispatched.
