# NEXT_STEPS.md — how to continue `biorheo` in a fresh session (any model)

Read this, then `README.md`, then the plan ledger at
`../biofluids_venture/docs/superpowers/plans/2026-09-05-phase0-foundation.md`. Do not
reconstruct state from chat history. State lives in the Obsidian card
`Documents\Obsidian Vault\Status\Biofluids Venture.md` and in this file.

## Before touching anything

```
cd C:\Users\compl\Desktop\Projects\biorheo
python -m pytest -q -m "not slow"      # must print: 137 passed, 1 deselected, 1 xfailed
git -C C:\Users\compl log --oneline -3  # top should be 5bc322d or later
```
If the count differs, find out why before adding work. Commit from `C:\Users\compl`
with `git add Desktop/Projects/biorheo && git -c core.safecrlf=false commit -m "..."`;
the CRLF warnings are noise.

## Rules that are specific to this repo

- Records under `src/biorheo/registry/records/` are GENERATED. Change the script
  (`scripts/records_from_carbopol.py`, `scripts/records_from_horner.py`), rerun, commit
  both. The loader rejects a record value outside its sidecar interval.
- Only `biorheo.identify.write_identifiability` writes `*.identifiability.yaml`. Never
  hand-write one; the loader rejects it unless `written_by: biorheo.identify`.
- Every number in a `.md` under `data/` or `case_studies/` is printed by a script. If you
  change prose that quotes a number, rerun the script and paste what it prints.
- `dst_fit` (`../dst_fit`) is consumed unchanged. `tests/test_gate_0_1_dstfit_unchanged.py`
  runs its suite; do not edit dst_fit to make biorheo convenient.
- No parameter value from memory. A record needs a DOI that resolves (offline cache in
  `src/biorheo/registry/doi_cache.json`; `resolve_doi(..., allow_network=True)` adds one).
- Writing a patch script through a bash heredoc breaks on `'''`; use the Write tool.
- Pre-registered bars are never loosened. A failed bar is recorded as failed with the
  number, in the test docstring, and a sharper test may be added beside it.

## Open work, in order

### 1. DONE 2026-09-05: the Horner protocol is settled from the thesis
Open-access PDF on udspace (URL + md5 in `data/horner2020/SOURCE.md`, verbatim quotes
there). At or below 2 1/s the steady value is the MAXIMUM of a 40-50 s transient; above
2 1/s the mean of the last 15 samples; 300 1/s preshear for 30 s. The decay after the
peak in the step sheet is attributed by the authors to syneresis (wall depletion), so
the earlier "structural family refuted on shape" claim was withdrawn. The text extract
`data/horner2020/thesis_text.txt` (gitignored, regenerate with pypdf from the PDF) is
greppable; sections 2.2, 3.3.1, 6.2.1 and 6.4 are the ones that matter.

### 2. Fit the structural family to the step RISE plus the steady curve (medium)
The step sheet's bulk signal, on the authors' reading, is the rise from the first
sample to the maximum (1-6 s at 1 1/s; donor D 16 s). `constitutive.structural`
(first-order kinetics, monotone build-up) is the right shape for that part and must
not be fitted to the decay. Pre-registered gate: per donor, fit steady curve (16
points, log-eta residuals) and rise samples (t <= t_max) jointly; the fitted
build-up time at 1 1/s must put the model's 95 % approach to steady state at 3-8 s on
at least 16/20 donors, with alpha and t_r each `measurable` or `bounded_only` (the
identify gate in `tests/test_structural.py` says they are degenerate from the steady
curve alone; the transient is what should break the degeneracy -- if it does not, that
is the result, record it). Before any of this, derive the constant-rate solution's
95 % time in `derivations/structural_kinetics.py` (it is a one-liner from fact 3) and
spot-check it numerically.
An independent test of the syneresis attribution is worth more than the fit: the decay
amplitude (1 - last/max, tabulated in `step_response.md`) against hematocrit and
fibrinogen from `extra` across the 20 donors. Pre-register the sign each mechanism
predicts before computing rho; n = 20 so quote the CI, not just the sign.

### 3. L3 solver bridge (plan unit 0.3-ish, not started)
Export a registry record to an OpenFOAM `viscosityModel` block (HB and Casson first)
and read it back. Gate: the round trip reproduces `eta_at` on 20 rates inside the
envelope to 1e-10, and refuses to export a record whose parameters are `degenerate`
or `not_identifiable` without an explicit override flag. Look at
`Desktop/Projects/bidet/` for an existing OpenFOAM case layout to copy from; do not
write a new one.

### 4. Second injectability surrogate on blood (small)
`biorheo.apps.injectability` is HB-only. Adding Casson needs the Casson needle-flow
closed form (Buckingham–Reiner for Casson), derived in `derivations/` first with the
same n=1 and yield→0 limits as `buckingham_reiner.py`. Then run the 20 blood records
through it and report how many are `extrapolated` at a 27G needle at 1 mL/min.

### 5. Referee items deliberately left undone
The finding-1 referee pass (2026-09-05) was absorbed: the authors'-residual runs test
and the run-count histogram print from `records_from_horner.py`; its synthetic-Casson
control (Cross/CY fit true-Casson + 1 % noise 2.5x WORSE than Casson, 20 trials) is
NOT yet a test in the repo and would be a good `tests/test_identify.py` addition.
Reviewer items deliberately left undone from the earlier pass: loader unit check
compares `repr(unit)` strings; `assert not base.fixed` missing in `regularise`; envelope
bounds not rounded outward; no note on dst_fit multi-start guess collision;
Wald-corr-at-bound caveat undocumented. None affects a shipped number.

## What NOT to do

- Do not fit any model to the post-peak decay of the step sheet (syneresis, per the thesis).
- Do not fit the LAOS sheets before item 2 has a model that reproduces the step rise.
- Do not add a record from a paper's table without its raw curve; the Horner case shows
  a published yield stress can be 71 % away from the same data under another objective.
- Do not put `si_fit` or any `impact_stress/materials.py` constant into the registry.
- Do not change `MAX_REPRO_DEV`, `RUNS_Z`, gate bars, or `LAMINAR_RE` to make a check pass.
