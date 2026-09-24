# PREREG — `stiffness`: is "jammed" a distinct dynamical state in the drop DNS, or just the viscosity cap?

Written 2026-09-24, before any stiffness value was computed (commit time is the record). Post-processing
only, of DNS cells already on disk. Script `stiffness_collapse.py` → `out/stiffness_collapse.txt`, written
and committed with this file.

## Why
`cvmap` found that no drop observable (footprint, peak spread, apex recoil) separates a jammed drop from a
thickened one near the window edge. The solver's own law explains why that could be structural rather than a
missing observable (`wc_drop_impact.c:487-500`):

  η = η_s · max(ε, EPS_MIN)^(−2),  ε = 1 − φ/φ_J(f̂),  φ_J(f̂) = φ₀ − (φ₀ − φ_m) f̂.

"Jammed" (f̂ ≥ f_c, ε ≤ 0) means only that the ε floor binds. The viscosity is continuous through the
jamming point, and a thickened region can reach it without jamming. So the **hypothesis H** is that the drop
responds to the peak viscosity it reaches, not to whether the closure calls it jammed. The **rival H0** is
that jammed drops respond differently (spread less) at matched peak viscosity. Consequence if H holds: this
closure cannot represent a jammed SOLID (no yield stress, no elasticity), only a cap on viscosity. So a
jam-vs-thicken observable cannot be found in this DNS, and must come from the bench or a different constitutive
model.

## Variables (defined before computing)
- For each time row, the peak viscosity in the drop is η at the row's `fhat_max` (η is monotone in f̂), with
  the cell's own x = φ/φ₀, χ = φ_m/φ₀ and EPS_MIN from its banner, MP_EXPO 2.
  **S = ln(η_peak/η_s)**, the max over rows with t/τ ≤ 0.5 (the impact phase).
- Secondary: **S̄**, the same but at `fhat_mean` (a bulk measure), max over t/τ ≤ 0.5.
- Response: **ln Rg_max** (primary; peak spread, the one-signed channel in `cvmap`); S2 (recoil,
  `grade_jamsig.sig`) reported as secondary.
- Class: J+ if J > 0.01 AND the banner's JAM_FHAT equals the cell's f_c = (1−x)/(1−χ) to 1e-3; J0 if the cell
  is outside its window (x < χ). A cell inside its window whose JAM_FHAT ≠ f_c is excluded (counted).

## Cells
Every cell under `out/<campaign>/` that is: L9 (`_L9_` tag), MP_EXPO 2 (no `_mpe` tag), not EPS_SOFT
(jamsig3 excluded), cs_fit preset, rc 0 with t/τ ≥ 5.9, and ADMISSIBLE (after t/τ 0.5: umax ≤ U₀, zc ≤ 0.5 D₀,
h ≤ D₀; the gate added after `cvmap`'s 27th cell). N is printed before and after each filter.

## Predictions (per U₀ group ∈ {1, 2, 4} with ≥ 6 cells)
- **H1 (collapse):** Spearman ρ(S, ln Rg_max) ≤ −0.7.
- **H2 (class-blindness, the discriminating test):** fit an isotonic DECREASING curve ln Rg_max(S) to the J0
  cells alone. For each J+ cell whose S lies inside the J0 cells' S range (interpolation only), take the
  residual r = ln Rg_max − fit. **PASS if |median r| ≤ 0.05** (the `cvmap` Δ floor, set from kmerdecomp's
  ≤ 0.03 activation channel). **H0 predicts median r < −0.05.** H2 is informative only if at least 3 interior
  J+ cells exist in some U₀ group. Its verdict is the pooled one over all informative groups (all must pass).
- The same H1/H2 on S̄ are reported beside, not graded.
- **Scope:** this is one solver and one closure. H passing says the DNS cannot tell the two states apart
  BY CONSTRUCTION of its constitutive law. It says nothing about real suspensions.
