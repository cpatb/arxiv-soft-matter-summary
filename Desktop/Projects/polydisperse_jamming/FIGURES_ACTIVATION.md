# Activation-shape figures (2026-09-22, block 3)

Built by `python make_activation_figures.py` (`--check` rebuilds and compares md5).

## J1 — the jamming threshold across the window, two ways

`figures/J1_threshold_across_window.{pdf,png,csv}`. `chi = 0.80`, `EPS_MIN = 0.01`, `MP_EXPO = 2`,
`delta` in {0.7, 1.3169, 2.6}.

- **(a) stress-controlled:** `sigma_jam/sigma* = [ln(1/f_c)]^(-1/delta)`, `f_c = (1-x)/(1-chi)`.
  Every curve passes through 1 at `f_c = 1/e` (dotted, x = 0.926), and the order in `delta`
  reverses there. ⚠ The crossing POINT is the sigma* convention (sigma* is defined at
  fhat = 1/e); the physical statement is only that the order reverses somewhere.
- **(b) rate-driven, the solver's own drive** (`sigma = eta(fhat) gd`):
  `gd_c eta_s/sigma* = max_{f<f_c} [ln(1/f)]^(-1/delta) / (eta(f)/eta_s)`. The viscosity
  feedback moves the controlling fraction from `f_c` down to a bottleneck `f*`, and the
  order in `delta` then reverses only where `f_c > ~0.95` — within ~0.01 of the lower window
  edge, ~5 % of the window instead of ~63 % — and is non-monotone in delta there. The kink
  near x = 0.985 is where the unactivated gap `1 - x` falls to the `EPS_MIN` floor.

So the sign of the activation shape's effect on WHETHER a point jams depends on whether the
point is stress- or rate-controlled, and the solver is rate-controlled. `PREREG_JAMWIN`
Amendment 1 tests which reading describes the solver's drop. Derivation:
`derivations/jam_stress_2026-09-22.md`, `derivations/jam_rate_driven_2026-09-22.py`.

## J2 — the measured activation, by CV (PREREG_ACTSHAPE)

`figures/J2_activation_by_cv.{pdf,png,csv}`. Constant-pressure shear, `p̂ = 3e-4`, N = 1000,
W4-healthy cells, Derjaguin threshold (`F_rep = κR*`), lengths in `d₁₀`.

- **(a)** `δ₁₀/₉₀` per (μ_p, CV), seed mean ± sd. Both friction and polydispersity SOFTEN the
  activation (lower δ). CV = 0.5 is drawn HOLLOW: every CV 0.5 cell is finite-size flagged
  (`d_max/L` 0.34–0.38 > 0.25) and the N = 4000 control is pending. Grey band: the solver's
  rheometry-fitted DELTA (cs_fit 1.3169, si_fit 1.3852) — shown for scale, NOT a test
  (different materials and a different fitting route), though it does sit inside the DEM's
  frictional–frictionless bracket at CV = 0.
- **(b)** `f̂(σ/κ)` pooled over seeds at μ_p = 0.3. The curves differ modestly, which is why
  the jam_map decomposition finds the activation change second-order to the window shift.
