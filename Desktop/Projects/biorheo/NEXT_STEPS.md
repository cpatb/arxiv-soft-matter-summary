# NEXT_STEPS.md — how to continue `biorheo` in a fresh session (any model)

Read this, then `README.md`, then the plan ledger at
`../biofluids_venture/docs/superpowers/plans/2026-09-05-phase0-foundation.md`. Do not
reconstruct state from chat history. State lives in the Obsidian card
`Documents\Obsidian Vault\Status\Biofluids Venture.md` and in this file.

## Before touching anything

```
cd C:\Users\compl\Desktop\Projects\biorheo
python -m pytest -q -m "not slow"      # must print: 137 passed, 1 deselected, 1 xfailed
git -C C:\Users\compl log --oneline -3  # top should be 39b0988 or later
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

### 1. Settle the Horner steady-shear protocol (small, high value)
`data/horner2020/step_response.md` shows the steady sheet's 1 1/s stress is the 4–6 s
peak of the 300→1 1/s step on 15/20 donors. What protocol produced the steady curve is
UNVERIFIED. Horner, Armstrong, Wagner, Beris, J. Rheol. 63, 799 (2019),
doi 10.1122/1.5108737, is titled "steady and transient shear"; its methods section
should say. Quote it verbatim into `SOURCE.md` and replace the UNVERIFIED tags. If it
is not accessible, say so in `SOURCE.md` with the date tried. Do not guess.

### 2. A thixotropic model that can overshoot (medium)
The structural family (`src/biorheo/constitutive/structural.py`) is monotone at constant
rate and every donor overshoots after the step-down. Candidates, each to be derived in
`derivations/` with SymPy checks before any code, in the pattern of
`derivations/structural_kinetics.py`:
- two-structure kinetics (fast breakdown, slow rouleau growth with its own rate),
- a viscoelastic stress with a structure-dependent modulus (e.g. a Maxwell element
  whose G and lambda depend on lam).
Acceptance gate, pre-registered here: fitted to the steady curve AND the step transient
jointly, it must put the stress maximum at 3–8 s and the 50 s decay at 6–31 % of the
maximum on at least 18/20 donors, with every parameter `measurable` or `bounded_only`
under `fit_rate_curve`-style profiling. If it cannot, that is the result: record it.
Data: `data/horner2020/step_300_to_1_human.csv` (1 s samples, 40–50 s, 20 donors) and
`steady_shear_human.csv`.

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

- Do not fit the LAOS sheets before item 2 has a model that survives the step test.
- Do not add a record from a paper's table without its raw curve; the Horner case shows
  a published yield stress can be 71 % away from the same data under another objective.
- Do not put `si_fit` or any `impact_stress/materials.py` constant into the registry.
- Do not change `MAX_REPRO_DEV`, `RUNS_Z`, gate bars, or `LAMINAR_RE` to make a check pass.
