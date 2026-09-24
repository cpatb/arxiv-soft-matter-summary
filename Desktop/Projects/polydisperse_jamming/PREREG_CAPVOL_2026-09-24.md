# PREREG — `capvol`: does the CAPPED VOLUME, not the jam class, set the drop's spread?

Written 2026-09-24, before any cell (commit time is the record). Follows `PREREG_STIFFNESS_2026-09-24.md`
(H1 PASS: peak viscosity orders the spread; H2 FAIL: at matched peak viscosity, jammed drops spread 3–17 % less).

## Idea
In the solver, jamming (f̂ ≥ f_c) and reaching the viscosity cap (ε ≤ EPS_MIN, i.e. f̂ ≥ f̂_cap =
(1 − x/(1−EPS_MIN))/(1−χ)) are the same kind of event, the floor on ε binding. A jammed cell reaches the cap
at f̂_cap < f_c. An unjammable cell with χ(1−EPS_MIN) ≤ x < χ reaches it at f̂_cap < 1 without jamming.
**Hypothesis HV:** the spread is set by how much of the drop sits at the cap (V), whatever the class. **Rival
HJ:** jamming matters beyond V.
V is recorded with no change to the dynamics: `jam_frac` counts the volume with f̂ > JAM_FHAT and is
diagnostic only (`wc_drop_impact.c:1480`). With JAM_FHAT = f̂_cap it is the capped-volume fraction.

## Cells (`out/capvol_cells.txt`; launcher `scripts/capvol_2026-09-24.sh`, G0-gated binary db413e87dd21)
All at EPS_MIN 0.03, L9, same settings as their source campaigns.
- 24 RERUNS: every `cvmapE03` grid cell and `kmermapE03` cell with f̂_cap < 1. The other 12 have f̂_cap ≥ 1,
  so V ≡ 0 by construction; they enter from their original runs with V = 0.
- 6 NEW, out of sample, x < χ so unjammable, yet f̂_cap < 1: singlets x 0.89 / 0.90 (χ 0.9083) and doublets
  x 0.86 (χ 0.8631), each at U₀ 2 and 4.
V = max over t/τ ≤ 0.5 of jam_frac(JAM_FHAT = f̂_cap). Class as before: J+ iff inside the window and J > 0.01
(J from the ORIGINAL run at JAM_FHAT = f_c). The new cells are J0 by construction.

## Predictions
- **G0 (wiring):** each rerun reproduces its original spread.dat EXACTLY on every column except jam_frac (7)
  and u_rms_jam (20). A mismatch in any other column means JAM_FHAT is not diagnostic-only, and V is void.
- **G1 (construction check):** V ≥ J(original) on every J+ rerun, because f̂_cap < f_c.
- **H3:** in each U₀ group with ≥ 6 cells, Spearman ρ(V, ln Rg_max) ≤ −0.7.
- **H4 (class-blind, discriminating):** in each U₀ group, fit an isotonic DECREASING ln Rg_max(V) to the J+
  cells. Every J0 cell with V > 0 whose V lies inside the J+ cells' V range gives a residual r = ln Rg_max − fit.
  Pooled over groups: **PASS if |median r| ≤ 0.05.** HJ predicts median r > +0.05 (unjammed at-cap drops
  spread more than jammed ones with the same capped volume). Informative only if ≥ 3 interior J0 cells pooled.
- Health: rc 0, t/τ ≥ 5.9, admissible (after t/τ 0.5: umax ≤ U₀, zc ≤ 0.5, h ≤ 1). Deaths are reported and
  never graded.

## Scope
One closure, one cap (0.03), L9. If HV passes, then in this DNS "jammed" has NO dynamical content beyond the
volume of fluid at the viscosity cap. The closure has no yield stress and no elasticity, so a bench signature of
jamming as a SOLID (rather than a very viscous liquid) cannot come from this DNS.
