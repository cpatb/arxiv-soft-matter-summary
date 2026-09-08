# Pre-registration — the extended-domain refit, and the §4 extrapolation

Written **2026-09-08, before the refit was run.** Not rewritten; errors go in an addendum.

## Why

`GRID_STUDY_2026-09-03.md` §4 records a defect that is not about the mesh at all:

> `optimum()` solves for `d*` then evaluates at fixed `r_max` = 15 mm, so
> **`r_max/d*` = 22.95 / 23.45 / 25.67 / 26.27 against a fit band of [2, 20]** — and the
> domain ends at `r/d` = 25, so the two finest evaluate *past the edge of the CFD that
> produced the fit*. The extrapolation is **worse for larger `C_P`**, i.e. worst for the
> "corrected" rows. **So the 23.6 % on `Q*` is a direction, not a number.**

The Status card names the repair as "the highest-value CFD run left": *extend the domain to
`Rdom` ≈ 50 mm and refit, which removes the §4 extrapolation rather than bounding it*.

**That run already exists and is converged.** `cfd_openfoam/case_Re10000_H6_R50` has
`Rdom = 0.05`, `refine = 1.0`, and its log ends `SIMPLE solution converged in 2689
iterations`. It has never been fitted, and `cfd_fit_summary.json` has never seen it. This
item does the fit and propagates it.

## What is NOT being claimed

- **This does not establish an observed order of convergence.** The JFE numerical-accuracy
  policy remains unsatisfied; `C_P` still diverges under mesh refinement (`|R|` = 4.59) and
  that finding is untouched here. Domain extent and mesh refinement are different axes.
- The wide case is at `refine = 1` only. Nothing here is a grid statement.

## The gates

### GB1 — the positive control, and it comes first

On the **matched** band `r/d ∈ [2, 20]` the wide case must reproduce the narrow one. If it
does not, the wide case is not a superset of the narrow one, the outflow boundary was
reaching into the interior, and nothing below means anything.

**Bar: `|ΔC_P|/C_P ≤ 5 %` and `|Δn_t|/n_t ≤ 5 %`**, matching the 5 % bar
`compare_domain_extent.py` already uses on the `τ_w` overlap (whose interior median is
0.691 % and max 1.791 %).

### GB2 — does one power law still describe the wall jet out to `r/d` = 40?

Refit over `[2, 40]`, which **covers every point `optimum()` evaluates**
(`r_max/d*` = 22.95–26.27) with margin.

**Bar: `r² ≥ 0.99`.** Chosen because the narrow-band fit returns 0.9968 and the project's
own reliability note already treats `r² ~ 0.94` as "good" and `0.77` as unusable; 0.99 is
the level at which the existing closure is quoted. **If it fails, the answer is that a
single power law does not reach `r/d` = 40 — which is itself the finding, and the fit band
must then be stated as a limitation rather than widened.**

### GB3 — registered predictions, so they can be scored

1. **`n_t` rises** on the wider band. Reason, computed in advance: the narrow fit
   (`C_P` = 1.3518, `n_t` = 1.9023) extrapolates to `τ_w` = **0.628 Pa** at `r/d` = 50,
   while the wide case *computes* **0.455 Pa** there — 27 % lower, so the far field decays
   faster than the `[2, 20]` power law.
2. **The paper's full-band extrapolated cleaning radius (39.14 mm) is too large**, and the
   outer-half one (34.75 mm) is close. Already corroborated by
   `compare_domain_extent.py` (computed 34.79 mm), and restated here so GB2's fit is
   consistent with it.
3. **`r_max/d*` under the refitted closure still lands above 20**, so the *narrow* band
   would still have been extrapolating — i.e. the defect was real and not an artefact of
   which `C_P` row was used.
4. **`Q*` moves by less than the 23.6 % the grid study attributes to `C_P`.** Reason: `Q*`
   depends on the closure through `d*`, and `d*` is set at `r_max` = 15 mm — inside BOTH
   bands — so a refit that mainly changes the far-field slope should move it little. **If
   this is refuted, the extrapolation was doing more work than §4 assumed and the §4
   sentence understates the problem.**

### GB4 — the honest denominator

Report the design numbers under **both** closures side by side, and state which band each
one's evaluation point falls in. No number from this item may be quoted without its band.
