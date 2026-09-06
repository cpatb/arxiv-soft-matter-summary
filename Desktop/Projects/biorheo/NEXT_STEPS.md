# NEXT_STEPS.md — how to continue `biorheo` in a fresh session (any model)

Read this, then `README.md`, then the plan ledger at
`../biofluids_venture/docs/superpowers/plans/2026-09-05-phase0-foundation.md`. Do not
reconstruct state from chat history. State lives in the Obsidian card
`Documents\Obsidian Vault\Status\Biofluids Venture.md` and in this file.

## Before touching anything

```
cd C:\Users\compl\Desktop\Projects\biorheo
python -m pytest -q -m "not slow"      # must print: 206 passed, 4 deselected, 1 xfailed
git -C C:\Users\compl log --oneline -3  # top should be 5aba072 or later
```
If the count differs, find out why before adding work. Commit from `C:\Users\compl`
with `git add Desktop/Projects/biorheo && git -c core.safecrlf=false commit -m "..."`;
the CRLF warnings are noise.

The four slow tests are the `dst_fit` suite gate, an OpenFOAM re-run of the plane-Couette
anchor (~25 min, needs WSL), the real-curve half of the model-comparison control, and a
re-run of the oscillatory-capability census (needs the gitignored workbook). Run them
with `-m slow` before believing anything about the L3 bridge.

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
- **Do not write a patch script through a bash heredoc.** It breaks on `'''`, on
  triple-quoted blocks, and on regex backslashes. Use the Write or Edit tool. This has
  now bitten four times, twice on 2026-09-05 alone, once silently: a `.replace()` whose
  target string did not match applied nothing, the script printed no error, and the
  anchor's JSON went out missing a whole key while its `.md` had the section.
- Pre-registered bars are never loosened. A failed or MIS-SPECIFIED bar is recorded with
  its number and a *stricter* test is added beside it. Five gates have gone this way:
  `openfoam_couette_anchor.py`'s A1 -> A1' and A5 -> A5'; the case study's ISO comparison,
  retracted outright; the step-rise shape gate G2, recorded FAILED; and item 6's G7 ->
  G7' (see below), where the defect was in the pre-registration, not the code.
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
less concave than a single exponential. Timescale right, shape wrong. Do not loosen G2.
What shape the rise IS was item 7, now closed — see there.

### 2b. CLOSED 2026-09-05 as not runnable: syneresis correlation test
Both candidate mechanisms predict the SAME sign of the fibrinogen dependence, by the
thesis's own account (p. 168 and its Eq. 8.1), so a correlation over 20 donors would pass
either. Recorded in `SOURCE.md` with the one available order-of-magnitude check.

### 3. DONE 2026-09-05 (all gates PASS): L3 OpenFOAM bridge, yield-stress families
`src/biorheo/solver/openfoam.py` exports a record to an OpenFOAM v2606
`transportProperties` model and parses it back. Three anchors outside the module:
`foamDictionary -precision 17` recovers every coefficient bit-exactly;
`derivations/couette_uniform_shear.py` derives that plane Couette imposes `gd = U0/H`
exactly (on monotone TAU, not eta); `scripts/openfoam_couette_anchor.py` runs `simpleFoam`
and compares OpenFOAM's own `nu` field and `wallShearStress` FO with `record.eta_at`,
touching none of this repo's transcription of `calcNu`. Pinned in
`tests/test_openfoam_export.py` and `tests/test_openfoam_couette.py`.

The findings, in `data/openfoam/couette_anchor.md`: the export is exact but the SOLVE is
not. `simpleFoam`'s outer iteration converges at a rate set by the local Bingham number
`B = tau_y/(k gd^n)`; at B = 29 it returns a viscosity 172 % wrong with `rc = 0`, no
warning and a written result. **That B and not the viscosity magnitude is the controller
rests on the constant-`B` 2x2 in `couette_diagnosis.json`** — at fixed shear rate, 100x
in `eta` at constant B changes the error in the tenth digit against 5.6e4x across B — and
NOT on any comparison between records, which do not control for magnitude (see item 6).
**The record's refusal to extrapolate does not survive the export**: a dictionary has no
envelope and a run started from rest has every cell below the lower bound. The wall
traction imbalance is a detector needing no reference answer (Spearman 0.98 with the true
error over 14 decades).

### 4. DONE 2026-09-05: Casson needle flow and 20 blood donors through a 27G needle
`derivations/casson_tube_flow.py` derives the Casson tube-flow equation and anchors it
three ways, including a velocity-profile route that carries the plug EXPLICITLY.
`injectability.flow_rate_casson` ships the FACTORED bracket; the expanded one goes
negative and is 6e2x wrong within 1e-6 of yield. `x = 1` is a TRIPLE root, so
`Q ~ (1-x)^3/3`. `case_studies/blood_injectability/results.md`: all 20 donors are
`extrapolated`, all by the same 1.42 decades; the yield stress is invisible in this
geometry; and because `F ~ R^-4`, a 10.5 % bore error reproduces the whole donor spread.
Pinned in `tests/test_casson_needle.py`.

**One claim was retracted here the day it was written** and the retraction is on the page:
the ISO 9626 needle-tolerance comparison. **Do not reinstate it without the standard.**

### 5. DONE 2026-09-05: referee items
- The synthetic-Casson control for the "Cross/CY fit 5x better" claim is
  `tests/test_model_comparison_control.py`.
- **The loader's unit check was a STRING comparison and that was a silent ACCEPT.**
  `str()` of a pint dimensionality rounds each exponent to 6 significant figures, and
  Herschel-Bulkley's `k` carries the record's own fitted `n`, so
  `Pa*s**0.3333333` and `Pa*s**0.33333334` stringify identically. `units.same_dimension`
  compares the `UnitsContainer` objects; the loader uses it. `tests/test_units.py`.
- **`regularise` now refuses a base model with `fixed` constants** for the yield-term
  schemes, which rebuild the law from parameters and would have dropped them silently.
  The clamp wraps `base.func` and is unaffected. `tests/test_regularise.py`.

STILL UNDONE, none of which affects a shipped number: envelope bounds are written with
`:.6g` and are NOT rounded outward, so a record's own extreme measured rate can fall
just outside its envelope (seen on the Carbopol records: raw max 199.9612229703406
written as 199.961); there is no note on dst_fit's multi-start guess collision; the
Wald-correlation-at-a-bound caveat is undocumented.

### 6. DONE 2026-09-05: L3 bridge for the families that are NOT yield-stress
`cross -> CrossPowerLaw`, `carreau_yasuda -> BirdCarreau`, both `calcNu` bodies quoted
verbatim in the module docstring. 40 new blood records (both families, all 20 Horner
donors) written by `scripts/records_from_horner.py`, which already fitted them for
`model_comparison.md`. Gate pre-registered in `data/openfoam/nonyield_bridge_prereg.md`,
committed at `2561346` BEFORE any code, with an addendum recording what the run found.

**Item 6's stated premise was FALSE and is recorded, not worked around.** It said
v2606's `BirdCarreau` has no `a` exponent and the export must refuse `a != 2`. It has one
(`getOrDefault("a", dimensionedScalar("a", dimless, 2))`, BirdCarreau.C:81-88, used in
`calcNu`), so no such refusal is shipped, and G4 decided that by running `simpleFoam`
rather than by reading C++. The hazard points the other way: `a` is OPTIONAL with a silent
default, so the export always WRITES it, the parser fills a MISSING one with 2, and a
MALFORMED one raises — dropping the semicolon after `a` used to read as absent and give a
law 11 % wrong with no error anywhere.

**The Cross names collide with the meanings SWAPPED**: biorheo's `lam` (a time) is
OpenFOAM's `m` (`dimTime`), biorheo's `m` (an exponent) is OpenFOAM's `n` (`dimless`), and
a name-matched map raises nothing because both slots take a plain scalar. G3 carries that
wrong map as a discriminator.

**`cross` and `carreau_yasuda` REFUSE an `eta_max`** rather than ignoring one: their
`eta_0` is a fitted parameter, not a regularisation, so `nu0` IS the record.

**G7 passed and was UNINFORMATIVE, and the defect was in the pre-registration.** G7
compared the new families at `relax_U = 0.95, endTime = 20000` against the 172 % figure
taken at `0.7, 2000`; at the anchor's own settings the Herschel-Bulkley record passes A2
too, so both arms passed and the comparison decided nothing. G7's number is kept and
**G7' runs all four records at matched `(0.7, 2000)`**, reporting itself UNINFORMATIVE
unless a yield-stress record actually fails there. G7' still does not separate yield from
viscosity magnitude — the plateaux differ by decades — and the report says so; the
magnitude control is item 3's constant-`B` 2x2.

Tests: `tests/test_openfoam_export_nonyield.py` (G1/G2/G3/G5/G8) and the item-6 block of
`tests/test_openfoam_couette.py` (G4/G6/G7/G7').

### 7. CLOSED 2026-09-05 as not runnable: is the blood step rise a lag or a second timescale?
`scripts/oscillatory_capability_horner.py` -> `data/horner2020/oscillatory_capability.md`,
pinned in `tests/test_oscillatory_capability.py`. Five criteria fixed before any was
graded; four FAILED.

**The one that PASSED is the question item 7 itself asked.** The LAOS sheets carry 20 to
1571 samples per cycle over 3 to 64 cycles, against a 4-6 sample rise — sampling was never
the constraint, so the item could not be closed on the ground it proposed. What closes it:
the workbook stores a **storage modulus and nothing conjugate to it** (no loss modulus, no
phase angle, no tan(delta)), and a lag and a second timescale differ chiefly in phase;
**linearity is never demonstrated** (the amplitude sweep starts at 9.95 % strain, at or
above the frequency sweep's own operating strain, with G' still rising there in both ramps
on 9/15 donors — and the Read Me's "a strain amplitude of 10" has two defensible units,
both carried, neither guessed); and the **population differs** — 12/20 donors on any LAOS
sheet, 10/20 on every sheet, 13/20 on the frequency sweep, against a step-rise result on
all 20, with A, B, C, D, F, Q, R and T having no LAOS data at all.

Two data-structure findings worth carrying: **both sweep sheets are two concatenated
ramps** (26+26 and 20+20), so an argmax over a column mixes two runs; and the axis and
donor columns must be paired row-by-row, not filtered independently. The two ramps differ
by a median 83 % at the same strain, which is too large to call them repeats — the
workbook does not say what the second one is.

Do not reopen this by re-asking about sampling. A discriminating measurement needs a phase
channel on the same donors, an amplitude sweep reaching into the linear region, or a step
sampled an order of magnitude faster on the rise.

### 8. NEXT ITEM, not started: the envelope is written rounded INWARD
Written for zero context, and it is small and self-contained.

`scripts/records_from_carbopol.py` and `scripts/records_from_horner.py` write the measured
envelope with `{gmin:.6g}` / `{gmax:.6g}`. Six-significant-figure rounding goes whichever
way the digits fall, so it can round a MAXIMUM DOWN and a MINIMUM UP, which puts the
record's own extreme measured rate OUTSIDE its own envelope. Measured 2026-09-05 on the
Carbopol records: a raw maximum of `199.9612229703406` is written as `199.961`, and
`FluidRecord.eta_at` raises `EnvelopeError` for anything above that. The blood records are
unaffected (0.1 and 700 are exact at 6 figures), so no shipped number moves — this is a
latent defect, not a live one.

The fix is in the GENERATORS, not the loader or the records (records are generated). Round
OUTWARD at 6 significant figures: floor the minimum, ceil the maximum, e.g.

```python
def round_out(x: float, up: bool, sig: int = 6) -> float:
    """Round to `sig` significant figures AWAY from the data, so a written envelope can
    never exclude a rate that was measured."""
    if x == 0.0:
        return 0.0
    q = 10.0 ** (math.floor(math.log10(abs(x))) - (sig - 1))
    return (math.ceil(x / q) if up else math.floor(x / q)) * q
```

Gate to pre-register before writing code (do not weaken it): after regenerating, for
EVERY record, the raw minimum and maximum rate of the source curve must satisfy
`record.envelope.contains(...)` — checked in the generator itself, which exits non-zero
otherwise, and pinned in `tests/test_registry.py` against the two source CSVs. Report how
many of the 75 records move and by how much; expect the 15 Carbopol ones to move in the
7th figure and the 60 blood ones not at all. If a record does NOT move, say so rather than
implying the fix was needed everywhere.

### 9. Also not started, smaller: the two remaining item-5 referee notes
No shipped number depends on either. (a) Document dst_fit's multi-start guess collision —
two starts that collapse to the same guess make the multi-start report a spurious
agreement. (b) Document the caveat that a Wald correlation computed at a parameter bound
is not interpretable, which is the situation for the `bounded_only` verdicts.

## What NOT to do

- Do not fit any model to the post-peak decay of the step sheet (syneresis, per the thesis).
- Do not reopen item 7 by asking again whether the LAOS sheets are densely enough sampled.
  They are; that was measured and it is not the constraint.
- Do not fit the LAOS sheets to the structural family as it stands: item 2 showed it does
  not reproduce the step rise's shape.
- Do not add a record from a paper's table without its raw curve.
- Do not put `si_fit` or any `impact_stress/materials.py` constant into the registry.
- Do not change `MAX_REPRO_DEV`, `RUNS_Z`, gate bars, `LAMINAR_RE`, `BAR_*`, `G7P_RELAX_U`
  or `G7P_END_TIME` in `scripts/openfoam_couette_anchor.py` to make a check pass.
  `G7P_*` are the configuration being MATCHED, not bars.
- Do not reinstate the ISO 9626 needle-tolerance comparison without the standard in hand.
- Do not quote a flow rate from `flow_rate_casson` within about a part in 1e9 of yield.
- Do not add an `a != 2` refusal to the BirdCarreau export. v2606 honours `a`; that was
  decided by running the solver, and the refusal item 6 asked for would refuse something
  OpenFOAM can do.
