# NEXT_STEPS.md — how to continue `biorheo` in a fresh session (any model)

Read this, then `README.md`, then the plan ledger at
`../biofluids_venture/docs/superpowers/plans/2026-09-05-phase0-foundation.md`. Do not
reconstruct state from chat history. State lives in the Obsidian card
`Documents\Obsidian Vault\Status\Biofluids Venture.md` and in this file.

## Before touching anything

```
cd C:\Users\compl\Desktop\Projects\biorheo
python -m pytest -q -m "not slow"      # must print: 179 passed, 3 deselected, 1 xfailed
git -C C:\Users\compl log --oneline -3  # top should be 8258662 or later
```
If the count differs, find out why before adding work. Commit from `C:\Users\compl`
with `git add Desktop/Projects/biorheo && git -c core.safecrlf=false commit -m "..."`;
the CRLF warnings are noise.

The three slow tests are the `dst_fit` suite gate, an OpenFOAM re-run of the plane-Couette
anchor (~8 min, needs WSL), and the real-curve half of the model-comparison control.
Run them with `-m slow` before believing anything about the L3 bridge.

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
- Writing a patch script through a bash heredoc breaks on `'''` and on triple-quoted
  blocks; use the Write or Edit tool. This bit twice on 2026-09-05.
- Pre-registered bars are never loosened. A failed bar is recorded as failed with the
  number, in the test docstring, and a sharper test may be added beside it. Three bars
  failed on 2026-09-05 and all three were replaced by STRICTER tests with the originals
  and their numbers kept (`openfoam_couette_anchor.py`'s A1 and A5, and the case study's
  ISO comparison, which was retracted outright).
- OpenFOAM case directories live in `data/openfoam/runs/` and are gitignored; the JSON
  and `.md` beside them are the committed artifacts. `BIORHEO_OF_OUT` redirects the
  anchor's output so a test run cannot overwrite them.

## Open work, in order

### 1. DONE 2026-09-05: the Horner protocol is settled from the thesis
Open-access PDF on udspace (URL + md5 in `data/horner2020/SOURCE.md`, verbatim quotes
there). At or below 2 1/s the steady value is the MAXIMUM of a 40-50 s transient; above
2 1/s the mean of the last 15 samples; 300 1/s preshear for 30 s. The decay after the
peak in the step sheet is attributed by the authors to syneresis (wall depletion), so
the earlier "structural family refuted on shape" claim was withdrawn. The text extract
`data/horner2020/thesis_text.txt` (gitignored, regenerate with pypdf from the PDF) is
greppable; sections 2.2, 3.3.1, 6.2.1 and 6.4 are the ones that matter.

### 2. DONE 2026-09-05 (result: shape gate FAILED): structural family on the step rise
`scripts/structural_step_horner.py` -> `data/horner2020/structural_step.md`, pinned in
`tests/test_structural_step.py`. G1 (t95 in 3-8 s, tau determined on >= 16/20) PASSED
20/20; G2 (rise rmse <= 2x tail noise on >= 16/20) FAILED 6/20, and two post-hoc variants
do not rescue it (7/20, 5/20). The residual has one sign pattern on 15/20: the rise is
less concave than a single exponential. Timescale right, shape wrong, and the shape
failure is in the kinetics, not inherited from the steady misfit. Open: what shape the
rise IS. A 4-6 sample rise cannot discriminate a lag from a second timescale from
strain-dependent kinetics; the LAOS sheets might, and are still unused. Do not loosen G2.

### 2b. CLOSED 2026-09-05 as not runnable: syneresis correlation test
Both candidate mechanisms predict the SAME sign of the fibrinogen dependence, by the
thesis's own account (p. 168 and its Eq. 8.1), so a correlation over 20 donors would pass
either. Recorded in `SOURCE.md` with the one available order-of-magnitude check. A
discriminating measurement needs a gap sweep or roughened walls; not in the workbook.

### 3. DONE 2026-09-05 (all gates PASS): L3 OpenFOAM bridge
`src/biorheo/solver/openfoam.py` exports a record to an OpenFOAM v2606
`transportProperties` viscosity model (HerschelBulkley, Casson) and parses it back.
Three anchors outside the module: `foamDictionary -precision 17` recovers every
coefficient bit-exactly; `derivations/couette_uniform_shear.py` derives that plane Couette
imposes `gd = U0/H` exactly (on monotone TAU, not eta -- every law here is
shear-thinning); `scripts/openfoam_couette_anchor.py` runs `simpleFoam` on 17 cases and
compares OpenFOAM's own `nu` field and `wallShearStress` FO with `record.eta_at`,
touching none of this repo's transcription of `calcNu`. Gates: A1' 2.7e-13 (bar 1e-9),
A2 9.6e-13, A3 5.5e-14, A4 3.0e-14 (bars 1e-6), A5' 8.7e-12 (bar 1e-9). Pinned in
`tests/test_openfoam_export.py` and `tests/test_openfoam_couette.py`.

The findings, in `data/openfoam/couette_anchor.md`: the export is exact but the SOLVE is
not. `simpleFoam`'s outer iteration converges at a rate set by the local Bingham number
`B = tau_y/(k gd^n)`; at B = 29 it returns a viscosity 172 % wrong with `rc = 0`, no
warning and a written result. That B and not the viscosity magnitude is the controller
was MEASURED (a 2x2 at fixed shear rate: 100x in eta at constant B changes the error in
the tenth digit, against 5.6e4x across B). Under-relaxing harder is worse, but there are
two regimes -- 0.99 dies with SIGFPE, so the shipped 0.95 sits between a slow regime and
a divergent one. **The record's refusal to extrapolate does not survive the export**: a
dictionary has no envelope and a run started from rest has every cell below the lower
bound. The wall traction imbalance is a detector needing no reference answer (Spearman
0.98 with the true error over 14 decades).

### 4. DONE 2026-09-05: Casson needle flow and 20 blood donors through a 27G needle
`derivations/casson_tube_flow.py` derives the Casson tube-flow equation and anchors it
three ways, including a velocity-profile route that carries the plug EXPLICITLY (the
quadrature could not have caught a dropped plug term -- it shares the same
integration-by-parts identity). `injectability.flow_rate_casson` ships the FACTORED
bracket; the expanded one goes negative and is 6e2x wrong within 1e-6 of yield.
`x = 1` is a TRIPLE root, so `Q ~ (1-x)^3/3`.
`case_studies/blood_injectability/results.md`: all 20 donors are `extrapolated`, all by
the same 1.42 decades (the wall rate is `4Q/(pi R^3)` for an effectively Newtonian fluid,
so the extrapolation distance is a property of the DEVICE); the yield stress is invisible
in this geometry; and because `F ~ R^-4`, a 10.5 % bore error reproduces the whole
donor-to-donor spread. Pinned in `tests/test_casson_needle.py`.

**One claim was retracted here the day it was written** and the retraction is on the page:
"the ISO 9626 tolerance for regular-wall 27G moves the force more than the donor spread".
The band mixed a regular-wall figure with a thin-wall one and the claim REVERSES under
the narrower reading. ISO 9626:2016 cannot be checked offline, so the section was rebuilt
on the `R^-4` scaling, which needs no standard. **Do not reinstate the ISO comparison
without the actual standard in hand.**

### 5. PARTLY DONE 2026-09-05: referee items
DONE: the synthetic-Casson control for the "Cross/CY fit 5x better" claim is now
`tests/test_model_comparison_control.py` -- on true-Casson data the richer families do
NOT win (>= 18/20, median ratio >= 1), and on the real curves the ordering reverses by
more than 2x. The contrast is the argument; neither half alone is.

STILL UNDONE, none of which affects a shipped number: the loader's unit check compares
`repr(unit)` strings; `assert not base.fixed` is missing in `regularise`; envelope bounds
are not rounded outward; there is no note on dst_fit's multi-start guess collision; the
Wald-correlation-at-a-bound caveat is undocumented.

### 6. NEXT ITEM, not started: L3 bridge for the models that are NOT yield-stress
Written for zero context. `src/biorheo/solver/openfoam.py` currently exports two of the
registry's seven constitutive families (`herschel_bulkley` -> `HerschelBulkley`,
`casson` -> `Casson`) and raises `NotImplementedError` for the rest. OpenFOAM v2606 ships
`CrossPowerLaw` and `BirdCarreau` in the same directory
(`/usr/lib/openfoam/openfoam2606/src/transportModels/incompressible/viscosityModels/`,
readable in WSL). Neither is biorheo's `cross` or `carreau_yasuda`: read the two `calcNu`
bodies FIRST and write the parameter map out before coding, exactly as the module
docstring does for the two existing families, because the names collide and the forms do
not. In particular OpenFOAM's `BirdCarreau` has no `a` exponent (it is fixed at 2) while
biorheo's Carreau-Yasuda does, so the map is partial and the export must REFUSE a record
whose `a` is not 2 rather than silently dropping it -- that refusal is the point of the
task, and it is the same shape as the envelope finding in item 3.

Gate to pre-register before writing code, by analogy with item 3 (do not weaken it):
the round trip reproduces `eta_at` on 20 rates inside the envelope to 1e-10; the export
refuses a `degenerate` / `not_identifiable` parameter without the override; the export
refuses a Carreau-Yasuda record with `a != 2` with a message naming `a`; and the
plane-Couette anchor in `scripts/openfoam_couette_anchor.py` is extended to one record of
each new family and must pass A1', A2, A3 at their existing bars. `RECORD_IDS` in that
script is the list to add to. Expect the Cross records to converge easily -- the Bingham
number finding in item 3 is about yield stress, and these families have none, which is
itself worth reporting as the control that shows the B story is not just "hard fluid".

### 7. Also not started, smaller: is the blood step rise a lag or a second timescale?
Item 2 left this open and said a 4-6 sample rise cannot discriminate. Before designing
anything, check whether the LAOS sheets in the workbook have enough samples per cycle to
constrain a second timescale; if they do not, close this the way 2b was closed rather
than fitting something the data cannot separate.

## What NOT to do

- Do not fit any model to the post-peak decay of the step sheet (syneresis, per the thesis).
- Do not fit the LAOS sheets to the structural family as it stands: item 2 showed it
  does not reproduce the step rise's shape. A LAOS fit needs a model that does.
- Do not add a record from a paper's table without its raw curve; the Horner case shows
  a published yield stress can be 71 % away from the same data under another objective.
- Do not put `si_fit` or any `impact_stress/materials.py` constant into the registry.
- Do not change `MAX_REPRO_DEV`, `RUNS_Z`, gate bars, `LAMINAR_RE`, or any `BAR_*` in
  `scripts/openfoam_couette_anchor.py` to make a check pass.
- Do not reinstate the ISO 9626 needle-tolerance comparison without the standard in hand.
- Do not quote a flow rate from `flow_rate_casson` within about a part in 1e9 of yield:
  the formula is fine there, but `1 - tau_y/tau_w` is limited by how well the inputs are
  known, which is two or three figures.
