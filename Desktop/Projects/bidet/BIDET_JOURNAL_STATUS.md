# BIDET — Journal-readiness status (CFD-validated drive)

**Date:** 2026-07-25 · **Driver:** `FINISH_BIDET_JOURNAL.md`, executed end to end.
Environment: WSL2 + OpenFOAM v2606, Python (numpy/scipy/matplotlib), MiKTeX.

## One-paragraph summary
The one hard blocker for a journal — that the turbulent wall-shear closure (C_P, n_t)
was *calibrated to a single commercial datum, not computed* — is removed. An
axisymmetric OpenFOAM (`simpleFoam` + k-ω SST) simulation of a round jet impinging on a
plate now **computes** the closure. The computation confirms the qualitative reversal
away from the laminar r⁻⁵ law, **independently corroborates** the commercial calibration
(so the "circularity" caveat is gone), and — importantly — **overturns one headline
lever**: standoff is a *weak* lever for a coherent jet, not the 4.6× the analytic
derivation claimed. That is a finding, not a failure, and it is reported as such.

## What the CFD found
- **Computed turbulent closure** (Phares form, developed wall-jet band r/d∈[2,20], at the
  nominal seated point Re_d=10⁴, H/d=6): **n_t = 1.90, C_P = 1.35, r² = 0.997.**
  Far from the laminar exponent 5 (reversal holds), but higher than the Phares
  developed-wall-jet 1–1.3 the paper had assumed.
- **Wall shear at the cleaning radius:** τ_w(r_max=15 mm) = **5.9 Pa** (direct CFD) for the
  commercial 1 mm nozzle at Re_d=10⁴ ≈ the assumed τ_c = 5 Pa → **commercial nozzles are
  near-optimal**, and the commercial-datum calibration is corroborated to ~15%.
- **Standoff is a weak lever (CORRECTION):** τ_w(r_max) is nearly standoff-independent
  (5.57, 5.87, 6.76, 6.16 Pa at H/d = 1, 6, 20, 60), so the flow needed to clean r_max
  varies only ~1.14× over H/d 1–20. The paper's analytic "standoff = 4.6× primary lever"
  is retracted. Physical reason: far-field wall shear is set by the jet's conserved
  momentum flux, not centerline velocity decay. Caveat: single-phase CFD; jet breakup at
  large standoff (not modelled) could restore a penalty.
- **Weak Re dependence** of (C_P, n_t) at fixed H/d — validates the Phares Re^{1/2} collapse.

## Validation performed
- **Free-jet decay:** B = 6.3 at H/d=20 (target 5–6, Hussein-Capp-George) — passes.
  H/d=60 over-predicts (domain confinement); large-standoff fits flagged unreliable.
- **Grid independence** (×2 refine at nominal): n_t within 0.5%, τ_w(r_max) within 6.5%,
  C_P within 12%.
- **Convergence:** the three H/d=6 cases meet the strict 10⁻⁵-on-all-fields criterion; the
  rest reach residuals ≤ 4×10⁻⁴ at the 3000-iteration cap (τ_w settled).

## Which numbers changed vs the estimates
| Quantity | Old (estimate/calibrated) | New (CFD-computed) |
|---|---|---|
| decay exponent n_t | 1.05 (assumed Phares) | **1.90** (computed) |
| prefactor C_P | calibrated to 1 commercial datum | **1.35** (computed, corroborated) |
| standoff lever | 4.6× (analytic, "primary lever") | **weak, ~1.14×** (CFD) — retracted |
| pressure exponent | ≈ ΔP⁻² (n=1.05 est.) / flat (laminar) | **Q\* ∝ ΔP⁻⁰·⁵⁷** |
| oscillation saving | muddled 55% / 73% / 80% | **~52%** (κ=2/π) |
| τ_w(r_max), commercial | asserted below laminar bound | **5.9 Pa ≈ τ_c** (near-optimal) |

## τ_c (remaining honest gap — improved)
A *genuine* wall-detachment shear measurement now supports τ_c: Nejadnik et al. (2008),
hydrodynamic wall shear at bacterial attachment/detachment balance = **1–3 Pa**, bracketed
above by Gião & Keevil (2013), ~24–144 Pa for consolidated biofilm. τ_c is still labelled
an explicit assumption (not fecal-specific); levers are ratios with Q* ∝ τ_c^{1/3}, weakly
sensitive. The Woolley bulk-rheology cite is no longer presented as the wall-detachment shear.

## Citations verified (Crossref)
Bhagat & Wilson 2016 (Chem. Eng. Sci. 152, 606) and Wilson et al. 2014 (Chem. Eng. Sci.
109, 183) confirmed; the `% FLAG FOR AUTHOR: verify` tags removed.

## Recommended venue (revised after the review pipeline)
The multi-tier review (professor lens) judged **Physics of Fluids unrealistic**: the fluid
mechanics is a single-phase RANS of a canonical round-jet-on-plate flow recovering a
Phares-type scaling, applied to a product-design question — the contribution is the
application and the CFD-anchored design reading, not new jet physics. Honest recommendation:
**post to arXiv `physics.flu-dyn` now; submit the journal version to European Journal of
Physics (or American Journal of Physics)** as a worked pedagogical fluids case study. A
primary fluids journal would additionally require the multiphase/free-surface computation
and a full grid-convergence (GCI) study, scoped as future work. The manuscript now states
this venue framing in its conclusions.

## Remaining author decisions
- Whether to commission the multiphase/VOF + breakup CFD (Tier B extension) to settle the
  standoff question for a *real* (possibly breaking-up) jet — the single-phase result only
  settles the coherent-jet case.
- Final venue choice (PoF vs EJP).
- τ_c: whether to run the microfluidic wall-detachment measurement on an actual soft
  biological soil, or accept the Nejadnik bacterial-adhesion proxy.

## Provenance / reproducibility
- CFD: `cfd_openfoam/` — `gen_case.py` (case generator), `run_sweep.sh`, `fit_closure.py`,
  `validate_freejet.py`, `cfd_fit_summary.json`, `sweep_results.jsonl`, per-case `tau_w_r.csv`.
- Model: `model/` — `cfd_closure()` added; `run_bidet.py` recomputes all levers; 13/13 tests.
- Figures: `fig_cfd_tau_w.png`, `fig_cfd_closure_vs_Hd.png`, `model/fig_tau_closures.png`,
  `model/fig_Qstar_vs_n.png`.
- Manuscript: `bidet_optimization_1.JOURNAL.tex` (CORRECTED/REVIEWED preserved).
- Git: baseline → Stage 1 (CFD) → Stage 2 (model) → Stage 3 (manuscript) commits, scoped
  to `bidet/` only.

## Note on the driver document
`FINISH_BIDET_JOURNAL.md` referenced `bidet_optimization_1.CORRECTED.tex` and a
`CORRECTION_DOSSIER`; neither exists. The newest, dossier-confirmed working manuscript is
`bidet_optimization_1.REVIEWED.tex` (with blockers B1–B6 flagged inline), used as the
JOURNAL base; the blockers were resolved here using the CFD + model results.
## Stage 5 — review pipeline outcome
Ran the three substantive review tiers (validation, postdoc rigor, professor significance) as
independent subagents; consolidated into `bidet_optimization_1_REVIEW_DOSSIER_2026-07-25.md`.
The review caught **three verifiable errors**, all now fixed:
1. **τ_c sensitivity used the laminar exponent.** The robustness claim quoted Q* ∝ τ_c^{1/3}
   (2.15× per decade, the *laminar* branch) while all numbers use the CFD closure, under which
   **Q* ∝ τ_c^{1.43}** (~27× per decade). τ_c is a *strong* lever — corrected throughout. This
   also makes the direct-τ_c measurement more strongly motivated.
2. **τ_w(r_max)=5.9 Pa was mis-attributed to ΔP=0.07 MPa.** It is the Re_d=10⁴ (≈40 kPa) value;
   at 0.07 MPa it is ≈8.6 Pa. "Near-optimal" is now presented as **τ_c-contingent**: the
   commercial nozzle cleans to r≈25/21/16/12 mm at τ_c=1/3/5/10 Pa (within ~2× of r_max across
   the band; over-provisioned at the lower, better-supported τ_c).
3. **Oblique a*(θ) exponent** had a spurious factor of 2 (fixed; Q*(θ) exponent was already right).
Also: standoff weak-lever claim scoped to the validated H/d≤20 (H/d≈60 flagged domain-confined);
CFD methods augmented (domain, BCs, schemes, residuals, cell count, SST stagnation caveat).

**Final compile:** MiKTeX pdflatex ×2 → **12 pages, 0 undefined references, 0 errors**, no
placeholders / `a^{?}` / FLAG comments. Model tests: **13/13**.

## Definition-of-done checklist
- [x] `cfd_openfoam/` runs, produces fitted (C_P, n_t) with validation; `fig_cfd_tau_w.png` exists.
- [x] `model/` uses the CFD closure by default; 13/13 tests; levers recomputed.
- [x] `bidet_optimization_1.JOURNAL.tex` compiles, presents CFD as the evidence, survived the
      review pipeline; citations verified; τ_c honestly handled.
- [x] `BIDET_JOURNAL_STATUS.md` (this file): CFD findings, changed numbers, venue, author decisions.
- [x] Commits separable per stage (baseline → CFD → model → manuscript → review/finalize),
      scoped strictly to `bidet/`.
