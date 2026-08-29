# Start here

State as of 2026-08-28, late evening. Full context is the vault card,
`Documents\Obsidian Vault\Status\dst-fit — Suspension Constitutive Fitting.md`.
Read its sections 1 to 3 and come back.

Two sessions worked on this package the same night — this one on G11, the Drop
Impact session on G12 and on a package-wide audit — so several things below are
corrections to claims that were in this file yesterday.

## Verify nothing rotted

```
cd C:\Users\compl\Desktop\Projects\dst_fit
python -m pytest tests -q     # 205 passed
python scripts/gate_n1.py     # 0
python scripts/gate_n2.py     # 0
python scripts/gate_n3.py     # 0
python scripts/gate_n4.py     # 1  <- correct, G3 fails on the 2023 corn starch
python scripts/gate_n5.py     # 1  <- correct, the joint fit loses by 55 AIC
python scripts/gate_n6.py     # 1  <- correct, G7a fails on Richards' corn starch
python scripts/gate_n7.py     # 0  <- ~3 min
python scripts/gate_g4.py     # 1  <- correct, G4a and G4b pass, G4c fails
python scripts/gate_g5.py     # 0
python scripts/gate_g10.py    # 0  <- sigma* vs a measured contact fraction
python scripts/gate_g11.py    # 1  <- NEW. Only G11c fails; see below
python scripts/gate_g12.py    # 0  <- NEW. Lin's solvent viscosities, and the 179x
python scripts/gate_n8.py     # 1  <- NEW. Steady-subset restriction; N8c fails
python scripts/audit_verdicts.py   # 0, prints a table. Not a gate
python scripts/probe_gauge_dependence.py   # prints two tables. Not a gate
python scripts/probe_two_clock.py          # prints a table. Not a gate
```

The gates needing the DataShare deposits are n6, n7, g4, g10, g11, g12 and the
audit. `DATA_INVENTORY.md` says where to get them.

## What is new since yesterday

### G11: the strain clock, the one constant a flow curve cannot contain

Full write-up in **`STRAIN_CLOCK_2026-08-28.md`**. In one paragraph: at steady
state the relaxation law's left-hand side is zero, so Richards Eq. 4's clocks
`gamma_0` and `t_r` are not merely ill-determined by a flow curve, they are absent
from it. Lin's reversal transients contain them. Inverting those records through
the same Wyart-Cates fit G10 validated gives `gamma_0`.

- **`gamma_0` is identifiable where `phi_m` is not, and provably so.** Writing
  `u = phi_J/phi`, the inversion is `f = A + B u` — affine in `u` for *any*
  packing triple — and a relaxation fit with both amplitudes free absorbs any
  affine map. So the clock does not depend on `phi`, `phi_m` or `phi_0` at all.
  Measured: bit-identical across four triples spanning `phi/phi_0` = 0.32 to 0.87.
  It depends only on `n`, and by 1.8% over `n` = 1.8 to 2.5.
- **The answer.** `gamma_0` = 0.54 to 1.59 on the 2 um silica and 0.71 to 1.64 on
  the 1.4 um PMMA, bracketed over three estimators. **Every one of 21 estimates
  exceeds the 0.2 the Basilisk solver borrows from Richards' corn starch**, by
  2.7x to 8.2x. Quote the ordering, not a value.
- **Wyart-Cates cannot reach the reversal dip.** The post-reversal minimum sits at
  0.56 to 0.74 of the model's own frictionless floor on every record, and Lin's
  PMMA at `gammadot = 0.05` is below it everywhere. No `f` reproduces that; the
  model has no microstructure variable. A WC closure will under-predict the
  softening after any flow reversal, which in a drop impact happens at the rim.
- **G11c fails on the silica and is untestable on the PMMA.** Silica clocks fall
  monotonically with rate (2.00x over 4 rates), which an additive time term cannot
  produce. PMMA clocks are *not* monotone (r² = -3.54), so the sign there is
  scatter. Do not quote this as "refuted on two materials".
- **And the single exponential is a blend.** Fitting two components jointly beats
  one by 19x in SSR for two parameters, and a shared-slow/per-rate-fast model beats
  that by 32 AICc. The rate trend survives the correction, shrinking 1.90x to
  1.70x. So the right sentence is "the clock falls with rate and correcting for a
  second component does not remove it", not "Richards Eq. 4 is refuted".

### G12: Lin's solvent viscosities exist, and the 179x is gone

`SOLVER_HANDOFF.md` §2 is rewritten. `Fig4a` is absolute viscosity and `Fig4b`
relative for the same sample, so `eta_s` is their ratio: **36.3 mPa.s** for the
silica and **1.770 Pa.s** for the PMMA, the latter agreeing to 0.29% with a direct
measurement in `FigS3`. Reduced stress cancels `eta_s` only *within* a dataset,
which is what §2 got wrong. The corrected gap is 24.6x, and the deposit's two
silica anchors disagree with *each other* by 180x, so this data cannot impeach the
solver's `sigma*` at all.

### The specification verdict now discriminates

`assess()` fired MISSPECIFIED on fourteen of fifteen real curves, which made it a
statement about data quality. `systematic_misfit()` adds the *size* of the
structure — the rms of a cubic through the ordered residuals, in log viscosity —
and the model is rejected only above `SYSTEMATIC_TOL = 0.05`, set from rheometer
reproducibility before the audit was run. That takes rejection from 15 of 23 real
curves to **5**. The two it rejects on Guy's 3770 nm series are 56% and 58%, which
are exactly the two the paper itself calls unsteady with massive wall slip.

`Assessment.verdict()` also now takes `profile=("phi_m",)` and reports "bounded
below by X, unbounded above" where that is the truth.

## The next job, in order

**A worked session plan with effort, ordering and decision points is in
`PLAN_2026-08-29.md`.** Read that instead if you are about to start work; the list
below is the same items with more background and less sequencing.

**Plan items 1, 2 and 3 are DONE as of 2026-08-29.** Item 3 fired its kill
criterion: see `STEADY_SUBSET_2026-08-29.md` and `scripts/gate_n8.py`. The
headline survives but the dramatic version of it does not, so **the spec must be
written around the weaker claim** in that file's section 3, not the old one.

### 1. Fix `joint.py`. It is rank-deficient by exactly one and prints nonsense.

**Verified twice tonight, not a suspicion.** `fit_series` frees `phi_0`, `phi_m`
and a per-curve `phi` simultaneously, and Wyart-Cates is homogeneous of degree one
in all three, so there is an exactly flat direction. Scaling the triple by `c`
changes every predicted viscosity by 2e-15; the residual Jacobian's smallest
singular value is 1.5e-5 against a largest of 119, with null vector proportional
to `(phi_0, phi_m, 0, phi)`.

`np.linalg.inv` does not raise on that. **That is why N5 prints
`phi_0 = 0.6912 +/- 924.4` and N6 `0.6778 +/- 1110`** — arithmetic noise, not wide
error bars. `n_params` at `joint.py:131` also counts the flat parameter, so every
joint AIC is 2 too high.

The conclusions survive (the joint fit loses by 55 and 568), but those numbers
must not be quoted. Fix: pin `phi_0` *or* `phi_m`, never both, and replace `inv`
with `pinvh` plus a rank check that raises.

### 2. Purge the same identity from the README

The covariance identity is quoted as a measurement in at least two more places:
"`sigma*` is orthogonal to `phi_0`, identical to four figures" and N6's "`phi_0`
moves the level of `phi_m` and leaves the reduced chi-squared bit-for-bit
unchanged, so candidate (1) is dead". Both are algebra and cannot fail. Check:
0.4634 x 0.64/0.546 = 0.5432, exactly the `phi_m` the gate reports.

The true statement is stronger: **`phi_0` is a pure gauge choice for a flow-curve
fit** — not hard to determine, not determined at all — and only `phi_m/phi_0` is
ever measured. Sweep `beta` instead, which moves `sigma*` by 15% over 0.8 to 1.2;
`beta = 0.8` fits Lin's silica 1.4x better than the assumed 1.

### 3. Three more from the audit, in descending severity

- **N4's "the fitter recovers the label, slope 0.970" is gauge-dependent** and so
  is not the validation it is presented as. Run on Guy's data, where `phi` is a
  genuinely known volume fraction, the same test gives slope **1.20 / 1.23 with a
  non-zero intercept**, which no `phi_0` removes. Promote the gauge-invariant
  version; it is a stronger check and it fails.
- **A bias the block bootstrap cannot see.** Refitting each curve on the top and
  bottom three-quarters of its own stress branch moves `phi_m` by +0.7 to +12.4
  Wald sigma, **always positive**. No resampling of residuals estimates a bias, so
  N7's margin against the needed 9.3-14.1x is a factor of ~1.5, not ~7. The
  degeneracy result is untouched.
- **N7 compares two different rulers.** `Profile.wald_width` uses
  `2 sqrt(f_crit) se` and `BootstrapInterval.wald_width` uses `2 x 1.95996 se`,
  differing by 4.9% at dof 26 and 13% at dof 5 — and that is the comparison the
  argument turns on. Also `block_bootstrap` refits with `n_starts=1` from a generic
  start rather than warm-starting, which is why its iid control comes back
  *narrower* than Wald.

### 4. Restrict every series to the curves its own authors call steady

Independently worth doing, and now supported from two directions. Richards'
0.46/0.48/0.50 curves are colour-coded in their own deposited `Figure1.py` as the
shear-jamming/oscillatory regime — time averages of an oscillating rate, not steady
flow curves. Restricting to the steady subset takes Richards from reduced chi² 221
to **11.0** and 0/7 to 3/4 within 2 sigma, and Guy from 344 to **36.9**, and on
Guy's steady subset the `phi_m` drift stops being monotone.

`systematic_misfit` reaches the same exclusion from the instrument side. Two
independent routes to the same subset is worth saying explicitly, because the
"monotone drift on three datasets" argument currently leans on the non-steady
curves.

### 5. Still open from G11

- **`Fig3` is the control on the method and is unused.** The simulation records
  evaluate contact and hydrodynamic stress *directly*, so `f` can be had without
  inverting anything. `read_lin_transients` already reads `Fig_3a`; pair it with a
  steady fit to `Fig_3c_total.csv` (26 rates) and the inversion is tested rather
  than assumed.
- ~~Copy the two-clock prototype into `scripts/`.~~ **Done** —
  `scripts/probe_two_clock.py`, numbers verified identical after the move. Its own
  output also refutes the follow-up guess that the fast component's share tracks
  `f_ss`: the share falls 0.465 to 0.357 while `f_ss` rises 0.116 to 0.885.
- **What form does fit the rate dependence.** G11c rules one out and proposes
  nothing.

## Housekeeping

- **Nothing is committed.** 44 files would be staged; `.gitignore` correctly
  excludes the 158 MB of deposits, the `DS_*.zip` and `logs/`. Verified with
  `git add -An`.
- **158 MB of literature data still lives inside this folder.** Consider moving
  the extracted trees to a sibling `dst_fit_data/` before packaging;
  `io_literature.py` takes paths, so nothing breaks.
- `figures/fig_phim_drift.md` carries a redraw brief, and the figure now
  understates its own result.
- A caution that cost time twice: a scripted documentation edit can silently match
  nothing. **Assert the replacement count before writing the file**, the same way
  you would check a numerical result.
