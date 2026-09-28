# NEXT_STEPS.md — how to continue `biorheo` in a fresh session (any model)

Read this, then `README.md`, then the plan ledger at
`../biofluids_venture/docs/superpowers/plans/2026-09-05-phase0-foundation.md`. Do not
reconstruct state from chat history. State lives in the Obsidian card
`Documents\Obsidian Vault\Status\Biofluids Venture.md` and in this file.

## Before touching anything

```
cd C:\Users\compl\Desktop\Projects\biorheo
python -m pytest -q -m "not slow"      # must print: 576 passed, 11 deselected, 1 xfailed (2026-09-28; 575 measured at 308 s before the mAb pin test was added -- the 576 is arithmetic, re-measure and correct this line if it differs)
git -C C:\Users\compl log --oneline -3  # top should be the item-29 commit ("G4c retired") or later
```
If the count differs, find out why before adding work. Commit from `C:\Users\compl`
with `git add Desktop/Projects/biorheo && git -c core.safecrlf=false commit -m "..."`;
the CRLF warnings are noise.

The seven slow tests are the `dst_fit` suite gate, an OpenFOAM re-run of the plane-Couette
anchor (~45 min, needs WSL), the real-curve half of the model-comparison control, a re-run
of the oscillatory-capability census, re-runs of the two seed-scatter studies (~20 min
and ~3 min), and a re-run of the item-12 injectability cross-check. All of them RESTORE the artifact they overwrite before asserting, so a
failure cannot destroy the reference it was checking against.

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
  its number and a *stricter* test is added beside it. Six gates have gone this way:
  `openfoam_couette_anchor.py`'s A1 -> A1' and A5 -> A5'; the case study's ISO comparison,
  retracted outright; the step-rise shape gate G2, recorded FAILED; item 6's G7 -> G7'
  (see below), where the defect was in the pre-registration, not the code; and item 11's
  G11.4b, whose `rc != 0` bar cannot tell "OpenFOAM rejected the dictionary" from "the
  case never got that far", so a stricter check is reported BESIDE it and the original
  verdict is kept as `pass_as_preregistered`.
- **A pre-registration is not rewritten.** Errors in one are corrected in an addendum at
  the foot of the file plus a one-line marker where the error is, never edited out --
  `newtonian_powerlaw_prereg.md` §0.2 is the worked example.
- **Adding a record to `RECORD_IDS` must not re-decide a closed gate.** It did, twice,
  silently: `g7_prime` iterated that list and `item6["A2_worst_all_records"]` was a max
  over every row on the page. Both now read a frozen list / one record's own rows. If you
  append to `RECORD_IDS`, grep for every other use of it first.
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

### 8. DONE 2026-09-06: the envelope was written rounded INWARD
`%.6g` rounds whichever way the digits fall, so a MAXIMUM could round down and a MINIMUM
up, putting a record's own extreme measured rate outside its own envelope --
`FluidRecord.eta_at` then refuses to evaluate where the data is. **Measured before the
fix: 12 of the 15 Carbopol records; 0 of 60 blood records (0.1 and 700 are exact at six
figures). After: 0 of 15.**

`registry.schema.round_out` rounds outward, working on the ROUNDED DECIMAL rather than on
`ceil(x/q)*q` -- the arithmetic form widens an exact 0.1 to 0.100001. One edge case the
property test caught: when `%.6g` crosses a power of ten (999999.9 -> 1e+06) the step
quantum must come from `x`'s exponent, not the rounded value's. Both generators now carry
the gate (write, read the record BACK, require its envelope to contain both extremes of
its own source curve, exit non-zero otherwise). `tests/test_registry.py` has a NEGATIVE
CONTROL that rebuilds the old envelopes and requires exactly 12 of 15 to fail -- if that
ever reads 0 the gate is inert.

**Two corrections to what the first commit of this item claimed, both found by review:**

1. It said BOTH generators had the gate. `records_from_horner.py` did not -- a patch had
   silently not applied and I trusted its "patched" print instead of grepping the file.
   Fixed and verified on disk. The blood records were never at risk (0.1 and 700 are
   exact) but the generator would not have caught it if they had been.
2. It said the anchor re-run left "every gate value BIT-IDENTICAL, so no shipped number
   moves". The six pass/fail scalars ARE bit-identical -- but they are MAXIMA over 33
   rows, and a max cannot see the other rows move. **7 of the 33 rows changed** (all
   seven Carbopol cases), and `verdict.A5prime_cauchy.worst_contraction_ratio` moved
   0.004368 -> 0.001856, a factor of 2.35, from a 5e-6 change in one input rate. This is
   the byte-identity-on-an-extremal trap from `lessons.md`, made by the person who wrote
   that lesson. The honest statement: the six verdicts and their scalars are unchanged;
   seven rows and one A5' diagnostic are not. Worth noting separately that
   `worst_contraction_ratio` -- the diagnostic added to catch a stalled A5' -- moves 2.35x
   under a 5e-6 rate change, so it is not carrying convergence information.

### 9. DONE 2026-09-06 (all three criteria PASS, after a RETRACTION): multi-start seed scatter
The referee item was "document dst_fit's multi-start guess collision". It was measured
instead: all 20 blood curves refitted five times under four models, changing nothing but
the optimiser's seed (`scripts/multistart_seed_scatter.py` ->
`data/horner2020/multistart_scatter.md`).

**Read the retraction in that file before using any of this.** The first version divided
the seed-to-seed span by `|mean|` and reported G9.2 and G9.3 as FAILED on one record --
Herschel-Bulkley, donor G, "4.149, i.e. 415 %, in `tau_y`". A docstring claim, a pinned
test and a new `strict_identifiability` API flag were built on that before a code review
caught it. The five seeds return `tau_y` = 1.61e-19, 9.81e-19, 5.00e-22, 3.27e-20,
5.94e-21 Pa against a smallest measured stress of 1.583e-03 Pa, with an identical profile
interval on all five (width 5.488e-04; quoted here as 0.01953 until item 13, and that
number was itself a bisection-grid point): **the optimiser converged the parameter to zero five
times out of five and the ratio was dust over dust.** Had `least_squares` returned exactly
0.0 the criterion would have passed. This repo's own lesson -- "a relative sensitivity
whose denominator is itself shrinking measures the shrinkage" -- applied to its own
instrument one session after being written down.

The metric now measures a span in units of the parameter's OWN 95 % profile interval; the
retracted ratio is still computed and printed so the number stays visible. Bars unchanged.
Results: **G9.1, G9.2 and G9.3 all PASS** -- 54 all-`measurable` records with a worst
spread of 2.972e-05, 63 export-accepted with the same worst. The two real movers are
Carreau-Yasuda donors B and L, both carrying `degenerate` / `not_identifiable`, which is
the positive content the criteria were written for. G9.3 is G9.2's predicate at a lower
bar over the SAME set, so they are one nested criterion, not two results.

The objective moves by at most 6.5e-10 in `rmse_log` (1.3e-9 in the least-squares cost),
so this is flatness and not a failed search.

**A real finding the review predicted:** the verdict SET is seed-dependent on 3 of the 80
pairs (Carreau-Yasuda A, B, L). Both `all_measurable` and `has_bad_verdict` are now taken
over all five seeds; the sidecar still records seed 0's only.

Also from this item: `check_exportable` gained `strict_identifiability=True` (off by
default). **Its original justification is withdrawn** -- what survives is that a genuinely
degenerate parameter whose covariance row is NaN is labelled `bounded_only` and accepted
by default, which is the gap the flag closes. Two other review fixes landed here: the
Wald-correlation-at-a-bound caveat is documented, and `verdict_for` now filters to FINITE
correlations before taking a max (`max([(nan,'a'),(0.995,'b')])` returns the NaN, so a
real degeneracy was found or missed depending on column order).

### 10. DONE 2026-09-06 (all three PASS): the same question, asked of Carbopol
`scripts/multistart_seed_scatter_carbopol.py` -> `data/carbopol_multistart_scatter.md`,
on an independent set: different material, instrument, 54 points per curve, and a
DIFFERENT FITTER for the shipped values. G10.1 and G10.2 pass at a worst spread of
5.409e-06; **G10.3 is the one specific to this set** -- because the shipped values come
from `carbopol/fit_hb.py`, seed scatter moves the SIDECAR rather than the record, and
`write_identifiability`'s containment rule had only ever been checked at seed 0. All 225
(seed, parameter, record) intervals contain the shipped value, and the verdict sets are
seed-stable on 15/15.

### 11. DONE 2026-09-06 (every gate PASS, and G11.6 REFUTED this item's own premise)
The L3 bridge now covers SIX of the registry's seven families: `newtonian -> Newtonian`
and `power_law -> powerLaw` are exported and parsed back, and `structural_steady` is
refused by a named `ThixotropyNotRepresentable` rather than a bare `NotImplementedError`.

Gate pre-registered in `data/openfoam/newtonian_powerlaw_prereg.md`, committed at
`e9c1ace` BEFORE any code (that commit contains that file and nothing else); the code
landed in `bcf84bf`. **Results are in the addendum at the foot of that file**, in
`data/openfoam/newtonian_powerlaw_gates.md` (written by `scripts/item11_gate_report.py`)
and in the item-11 section of `data/openfoam/couette_anchor.md`. Tests:
`tests/test_openfoam_export_null.py` and the item-11 block of
`tests/test_openfoam_couette.py`.

**This item's stated reason for refusing `structural_steady` was FALSE, and G11.6 was
written to find that out before the reason was put in an error message.** It said the
family "is thixotropic and has no counterpart there". It has an exact counterpart: the
steady flow curve of first-order structure kinetics IS `CrossPowerLaw` with n = 1, to
1.386e-16 relative -- `STRUCTURAL_STEADY.limits` already said so. Refusing on "no
counterpart" would have refused something OpenFOAM can do, which is item 6's mistake
repeated. The real reason is narrower and worse: `viscosityModel` is a function of the
instantaneous `strainRate()` alone, so the STRUCTURE EQUATION cannot be carried and the
export would be right in steady state and silently wrong in every transient. The message
says that and gives the deliberate substitution -- **per family**, because
`structural_steady_separate` carries `alpha` and `t_r` apart and has no `alpha_tr`.

**Three things the OpenFOAM source said that this item did not.** All were read before
coding, as the item required.
- `Newtonian` is the ONLY viscosityModel with no `<Model>Coeffs` sub-dictionary: it reads
  `nu` from the top level (`Newtonian.C:71`). The export writes no block and the parser
  REFUSES one rather than ignoring it.
- `powerLaw`'s `nuMin` and `nuMax` are REQUIRED (`readEntry`), not optional like
  `BirdCarreau`'s `a`, and `eta = k gd^(n-1)` has no plateau at either end -- so both are
  the caller's regularisation. New keyword-only `eta_min`, required for `power_law` and
  refused for the other five, each with its own reason.
- `powerLaw` floors `SMALL = 1e-15` on the DIMENSIONLESS `1[s]*strainRate()`, not
  `VSMALL`. (The prereg's §0.2 then says "the other four models use VSMALL", which is
  itself false -- only HB and Casson do. Corrected in the addendum, not edited out.)

**Records: 75 -> 115.** `scripts/records_from_horner.py` writes `newtonian` and
`power_law` per donor. They are NULL MODELS and say so in `notes`; the item forbids
inventing a record, and these are the same generator, source curve, DOI and envelope as
the four existing families. `model_comparison.md` gains the null columns and reports what
they say: **Casson beats a one-parameter constant on 20/20 donors but a two-parameter
power law on only 14/20**, at the same parameter count, and the `residual_structure` flag
does NOT separate a fit from a null -- 20/20 for both nulls AND for Herschel-Bulkley; the
only forms it clears are Cross and Carreau-Yasuda.

**What the gates do NOT establish.** Both G11.5 wrong-map separations are CONSTANTS:
G11.5a is `rho - 1` = 1059 on all 20 records and G11.5b is `1/gammadot_min - 1` = 9 on
all 20. They show the slot is USED and no more; neither is sensitive to the fitted values.
G11.5b's attribution survives its own control (the clamp moved 10x -> 100x leaves it at 9).

**Two code reviews, nine + four findings, all acted on** (`5f6ac8c`, `04332dd`). The one
that was a silently WRONG number: the Newtonian parser stripped sub-dictionaries with a
regex, and `[^}]*` cannot match balanced braces, so everything after a NESTED block leaked
out as top level -- it returned a buried `nu = 999`, and returned it even with the real
entry deleted. Four more of the same family as defects this repo has hit before:
- `G7'` iterated `RECORD_IDS`, and `item6["A2_worst_all_records"]` was a max over every
  row on the page, so appending two records here would have re-decided a CLOSED item-6
  gate -- and G7's yield/no-yield split keys on `NO_CLAMP_MODELS`, which lists neither
  null family, so both NULLS would have counted as yield-stress arms. Both now read a
  frozen list / the HB record's own rows.
- The anchor's exit code did not carry G11.7 or the refusal probes (rc = 0 on a failure).
- The per-family aggregation used builtin `max()`, which DROPS a NaN.
- Two docstring sentences said G11.3b/G11.4b "were decided by running the solver" while
  the anchor was still running.
- G11.4b could pass on a blockMesh failure with a BLANK message. The pre-registered bar
  (`rc != 0`) is unchanged and still recorded as `pass_as_preregistered`; a stricter
  check is reported beside it.
- The `foamDictionary -entry nu` test was cited as a second check of the Newtonian
  layout. It is not -- both candidates pass it. G11.4b is the only discriminator.

### 12. DONE 2026-09-06 (the prediction was REFUTED, and three of the write-up's own claims were too)
`case_studies/blood_injectability/crosscheck.py` -> `crosscheck.md` / `crosscheck.json`,
pinned in `tests/test_injectability_crosscheck.py`. Gate pre-registered in
`crosscheck_prereg.md`, committed at `1cd46de` BEFORE any of it was computed.

**Item 12 predicted the two forces would agree. They do not.** The worst
Cross-vs-Casson plunger force difference is **30.18 %**, median **-11.19 %**;
Carreau-Yasuda 33.02 %. That is ten times the yield term the case study calls
negligible and thirty times the calculation's own ~1 % resolution. The item's premise
(*the yield term contributes under 3 %*) is TRUE and was reproduced first; the conclusion
does not follow from it.

**The case study's ordinal claim SURVIVES** -- 30 % is well inside the 128 % donor spread,
so which patient the blood came from still matters more than which family was fitted, and
the yield stress is still invisible in that geometry. **Its precision does not**, and
`results.md` now carries that caveat, printed by `run.py` from `crosscheck.json`.

**The mechanism is MISFIT, not extrapolation, and the first write-up said the opposite.**
Casson's `rmse_log` is worse than Cross's on 20/20 donors (1.5x to 12.3x) and
`residual_structure` flags 19/20 against 2/20; the force difference tracks that misfit at
Spearman +0.559 (p = 0.010). Two checks killed the extrapolation story: the same fits
differ by up to 42 % at rates INSIDE the envelope (worse at the bottom than at the
needle), and re-running the whole comparison at a needle sized so the wall rate lands in
the measured range (ID 0.6236 mm) leaves **99.5 % of the median difference intact**.

**Do not describe this as "family choice moves the force".** It is not symmetric: the
two-parameter family the case study ships is the worse-fitting one on every donor.

Three claims were shipped and then corrected on the page (`439410f`), all found by review:
- *"both families predict a consistently LOWER force"* -- FALSE, 4/20 Cross donors are
  positive, up to +14.19 %. Use the paired sign test (16/20, p = 0.0118).
- *"the plateau residual is the size of the yield term"* -- FALSE and it is a SIGNED
  claim, so the unsigned bound that tested it could not see that 6/20 Cross and 7/20 CY
  donors have the wrong sign. Withdrawn.
- G12.3's Newtonian null does not license what it was used for: for a Newtonian fluid
  `dp = 8 mu L Q/(pi R^4)` exactly, so its 271 % is `(mu/eta_c)/(1+y)`, an identity
  restating "blood shear-thins 2-4x over the measured range". `power_law`, with no plateau
  at all, lands in the SAME band as `cross`, so the four-band verdict is coarser than the
  claim. **The bar was not changed**; both are recorded as findings.

**The null that does discriminate** (supplied by the review): Cross with `eta_inf`
replaced by that donor's Casson `eta_c`, keeping `eta_0`, `lam`, `m` fitted -- not Casson
in disguise. It collapses 30.18 % to **2.72 %**, inside the yield band. That, and not the
plateau-ratio arithmetic, is what the mechanism rests on.

A fourth control is on the page labelled as NOT independent evidence: a Newtonian fluid at
`eta_c` is Casson with `tau_y = 0`, so it recomputes the sibling study's
`yield_share_max`. It agrees to 1e-5, which makes it a cross-implementation check of that
number and nothing more.

### 13. DONE 2026-09-06 (all seven gates PASS; G13.9, added afterwards, FAILED): the profile intervals ended on a bisection grid

`dst_fit.profiling.profile_likelihood` sizes BOTH its outward step floor (`span/5000`)
and its ten bisections from the parameter's declared BOX BOUND. For a viscosity bounded
at `1e6` the first step is 200 on a parameter of value 0.0038, and the returned endpoint
is `best + m*span/5.12e6` for an integer `m` -- a bisection-grid point, not a likelihood
crossing.

Gate pre-registered in `data/profile_interval_prereg.md`, committed at `05abaee` BEFORE
any fix code (that commit contains that file and nothing else). Addendum A, at `eac96db`,
records what the pre-registration itself got wrong and pre-registers G13.9. **A
pre-registration is never rewritten**; sections 0-6 stand as committed.

**The item's own numbers were UNDERSTATED and the reproduction is what found that.**
`scripts/profile_interval_diagnosis.py` reads each model's declared bound rather than
three guessed candidates: **201 of 296 closed upper ends** on the grid (this file used to
say 117), worst **1639x the best-fit value** above it (used to say 82x).

**G13.2 is the gate that could have made this cosmetic.** Evaluating the profile F
statistic AT the shipped endpoint over the 70-end pre-registered sample: **median 710x
`f_crit`, worst 8130x**. A confidence limit has `F = f_crit` there, by definition.

The fix is `src/biorheo/identify/profiling.py`, called by `fit_rate_curve`. `dst_fit` is
untouched (G13.8 checks a pathspec-scoped `git status` plus the pinned figure md5). It
imports `dst_fit`'s own private `_nuisance_minimiser`, so **only the WALK differs** and a
difference between the two routines cannot be an objective difference; that coupling is
pinned by `tests/test_profiling_wrapper.py`.

Gates (`scripts/profile_interval_gates.py` -> `data/profile_interval_gates.md`):

| gate | result |
|---|---|
| G13.2 shipped endpoint is not the crossing | PASS, 70/70, median 710x `f_crit` |
| G13.3 corrected endpoint IS the crossing | PASS, worst 2.8e-4 against a 1e-2 bar |
| G13.4 matches a CLOSED FORM using none of this code | PASS, worst 2.4e-9 against 1e-6 |
| G13.5 the correction can only NARROW | PASS, 98.3 %, worst narrowing 1227x |
| G13.6 items 9 and 10 re-graded, bars unchanged | PASS, all six |
| G13.7 loader containment survives | PASS, 115 records, 325 params, 0 failures |
| G13.8 `dst_fit` untouched | PASS |

**G13.4 is the load-bearing one.** G13.3 grades this module with this module's own
machinery. The `newtonian` closed form does not: for a one-parameter fit in log space the
profile is exactly `SSR(u) - SSR_min = n (u - ybar)^2`, so the crossing is
`exp(ybar +/- s sqrt(f_crit/n))`, and it passes through neither the walk nor the root
finder. Derived in the pre-registration, not reverse-engineered from the answer.

**Three things the run found that the pre-registration did not anticipate.**

1. **On two blood Carreau-Yasuda records the profile is NOT A FUNCTION of the
   parameter.** The nuisance re-minimisation lands on different local optima at
   neighbouring values, so `F` moves 0.239 -> 10.06 between points 0.1 % apart. A sign
   change exists and a root finder returns a number from it; that number is a
   DISCONTINUITY. Such an end is now rejected, reported OPEN, and carries the value and
   `F` that were rejected so the rejection is auditable. 5 ends across donors B and L;
   11 further apparent crossings were refuted mid-walk by a second nuisance start.
2. **The first version of that check rejected a GOOD end** (donor C, `F = f_crit` to
   3e-4) because the flanking test compared without a noise band, and just outside a
   genuine crossing `F - f_crit` is below the inner solver's own noise. Over-rejection
   discards good data as surely as under-rejection keeps bad data.
3. **The prediction that intervals narrow "dramatically in most cases" is wrong.** The
   MEDIAN narrowing is 1.0002x. The effect is a TAIL -- 201 of 296 upper ends -- not a
   shift.

**Registry effect.** Verdicts move ONLY in the `carreau_yasuda` column, on 4 donors, and
only toward MORE flagging (`eta_0` measurable 20 -> degenerate 2 / measurable 18). No
fit, rmse or verdict changed for `casson`, `herschel_bulkley`, `cross`, `newtonian` or
`power_law`. Item 9's worst spread moved 2.972e-05 -> 2.973e-05 and its `casson` column
moved 8.9x (3.51e-07 -> 3.14e-06) -- the predicted direction, from a number 300x under
its bar. Item 10's `.md` is bit-identical. Every number in the blood injectability case
study is unchanged; only `git_head` moved in its `results.json`.

**G13.9 (addendum A.2, `eac96db`) FAILED both arms, against a prediction of 20/20 on
both.** It asked step 4's question: is item 12's Cross-vs-Casson force disagreement
resolved at 95 %?

- **G13.9a FAIL, 9/20 excluded against a bar of 11/20.**
- **G13.9b FAIL, 12/20 contained by the old hull against a bar of 19/20.**

**The cause is a sign error in the reasoning, not in any code.** The bug inflated the
UPPER ends only, but the force difference is SIGNED -- Cross below Casson on 16/20 -- so
the comparison turns on the LOWER end of the hull, which was never corrupted. Measured:
the binding end is the lower one on 8 donors and the upper one on 1; on the old intervals
lower on 8, upper on 0. **Exactly one donor of twenty flips.** G13.9b was additionally
close to undecidable by construction, since narrowing an interval can only shrink a
push-through range, so `old_contains >= new_contains` always.

**What IS established, and it is the useful part** (`data/profile_interval_g139.md`):
the Casson push-through hull is a median **29.1 %** of the force at the 27G needle, the
covariance-sampled joint interval **25.2 %**, against a median family difference of
**12.5 %**. **A record's own parameter uncertainty there is LARGER than the family-choice
effect.** So item 12's precision caveat is STRENGTHENED, not lifted, and `crosscheck.md`
now says so: do not quote a needle force from that study to better than about 30 %.
Item 12's mechanism is untouched -- Casson misfits on 20/20 donors in `rmse_log`, which
is not a statement about parameter uncertainty.

Tests: `tests/test_profiling_wrapper.py` (16), including a NEGATIVE CONTROL that asserts
`dst_fit`'s endpoint is still on the grid and still past `f_crit` -- if that ever stops
being true the wrapper is unnecessary, which is worth failing a test over.

### 14. DONE 2026-09-06 (four gates PASS, G14.5 FAILED as predicted): the case study now carries uncertainty, and it cannot rank donors

`case_studies/blood_injectability/uq.py` -> `uq.md` / `uq.json`, with a section added to
`results.md` by `run.py`. Gate pre-registered in `uq_prereg.md`, committed at `3d267db`
BEFORE any of it was computed. Pinned in `tests/test_injectability_uq.py`.

Phase 0's unit 0.6 (`../biofluids_venture/05_ROADMAP.md`) requires this study to carry
"full UQ from the L2 posterior" and to state which parameters the force is insensitive
to. Both are now on the page.

**The intervals are READ from `data/profile_interval_g139.json`, not recomputed.** That
artifact already pushed each Casson parameter through its corrected interval at this
exact needle. G14.2 and G14.3 exist to check that borrowing is legitimate.

| gate | result |
|---|---|
| G14.2 the borrowed intervals belong to these forces | PASS, worst 1.2e-06 relative over 20 donors |
| G14.3 both modules describe the SAME device | PASS, all four constants exactly equal |
| G14.4 an open interval is never printed as a number | PASS, 0 donors |
| **G14.5 donors can still be ranked at 95 %** | **FAIL, 36/190 pairs = 18.9 % against a 50 % bar** |
| G14.6 insensitivity by an independent route | PASS, `tau_y` insensitive 20/20 |

**G14.5 is the result and it FAILED.** Only **36 of 190** donor pairs have disjoint 95 %
force intervals, and ranked by force **two donors must be 17 places apart out of 20** to
be reliably separated. The pre-registration predicted a failure at 35-50 %; the measured
18.9 % is worse, so the prediction was right in direction and wrong in size. Both are
recorded, and the bar was not moved.

**Why that gate and not the obvious one.** Re-checking "the 128 % donor spread beats the
12.5 % family effect" cannot fail — adding a 29 % uncertainty does not change that
arithmetic. The comparison that CAN fail is the pairwise one, because **the hull is PER
DONOR and the spread is ACROSS donors**. A customer does not use the spread; they ask
whether this patient needs more force than that one.

**A sentence on the same page had to be corrected**, and it is the reason this item was
worth doing. The crosscheck section said *"which patient the blood came from still
matters more than which family was fitted"* — a range-versus-effect-size comparison read
as a per-patient claim — twenty lines below the new section saying donors cannot be
ranked pairwise. The RANGE comparison is kept because it is true; the pairwise reading is
withdrawn on the page, with the old wording quoted so the correction is visible.

**Sizes, and they are different objects.** The one-at-a-time push-through hull is a
median **29.1 %** of the force; the covariance-sampled joint interval is **25.2 %**; the
donor spread is **128 %**. The hull is WIDER than the joint region it is drawn from (the
joint region is an ellipse inscribed in that box); the sample IS joint but assumes the
Wald normal, which is the approximation a profile likelihood exists because it distrusts.
**Neither is simply the better number and the page names which one it is quoting, every
time.**

**G14.6 is a third independent route to a claim the page already made twice.** `results.md`
argues the yield stress is invisible from `x = tau_y/tau_w <= 1.7e-4` and from zeroing
`tau_y` (under 2.9 %). `sensitivities()` is a finite-difference `d ln F / d ln p` on the
force solver and shares neither construction: `|d ln F / d ln tau_y| <= 1.48e-02` on
20/20, while `|d ln F / d ln eta_c| = 0.995` — the force is very nearly proportional to
the plateau viscosity and to nothing else. All three agree. The test asserts the SECOND
half too, because a sensitivity route that called everything insensitive would pass the
first half and mean nothing.

**Not added into these intervals:** the family term (median 12.5 %, worst 30.2 %), the
needle bore (`F ~ R^-4`, still [UNVERIFIED]), and any error in the Casson FORM, which
misfits these curves on 20/20 donors. A correct interval on a wrong model is still a
wrong uncertainty.

### 15. DONE 2026-09-06 (three gates PASS; the ITEM'S OWN PREMISE was false): how far item 13 reached the Carbopol study

`case_studies/injectability_surrogates/uq_check.py` -> `uq_check.md` / `uq_check.json`,
with a section added to `results.md` by `run.py`. Gate pre-registered in that study's
`uq_prereg.md`, committed at `850a647` before anything was recomputed. Pinned in
`tests/test_surrogates_uq_check.py`.

**This item's stated premise was FALSE and one grep settled it.** It said the Carbopol
study "ships point forces too" and asked for a UQ column. `run.py` has called
`force_with_uq` since it was written, and `results.md` has carried a profile-envelope
column, a covariance column, a per-parameter sensitivity line and an `Insensitive to:`
line all along. **Fourth time in this repo** (items 6, 11, 13 and now 15).

**What is actually true is worse and more useful:** those intervals were drawn by
`dst_fit.profiling` before item 13, so the page's shipped error bars were bisection-grid
artefacts. The item became a CORRECTION to a shipped page rather than an extension.

| gate | result |
|---|---|
| G15.1 the correction reaches this page | PASS, worst endpoint move 2.94 % against a 1 % bar |
| G15.2 the corrected intervals NARROW | PASS, 2/2 |
| G15.3 no end is silently open | PASS, 0 open ends |
| G15.4 is the parameter interval still dominant | nothing flipped; ratios 7.5x and 2.3x |

**THE PARAMETER INTERVALS MOVED AND THE FORCE DID NOT.** Endpoints move up to 2.94 %;
the push-through hull moves at most **0.0081 %**. The mechanism was measured rather than
argued: item 13's step floor is `box_bound/5000`, so **the defect's size scales with the
declared box**. Herschel-Bulkley's boxes are not alike — `tau_y` at `1e5`, `k` at `1e6`,
the flow index `n` at O(1) — so `k` moved 2.94 %, `tau_y` 1.87 % and **`n` 0.0011 %,
three orders less**. And the force is set by `n` on BOTH gels
(`d ln F / d ln n` = 4.87 and 7.44).

**That is the same shape as item 13's own G13.9 failure, one item earlier.** A correction
reaches a downstream number only through the parameters that number is sensitive to;
there it was the wrong END of an interval, here the wrong PARAMETER. Neither was
predicted. The pre-registration's SPLIT prediction (flip on the 0.5 wt% gel, not the
0.1 wt%) has its flip half **refuted**, and the reasoning behind it was backwards — it
assumed a grid inflates narrow intervals most, when the grid's size comes from the box.

**Two of my own claims were corrected mid-item, both written before reading the rows.**
The report's narrative was drafted assuming G15.1 would FAIL and it passed. And both the
report and `profiling.py`'s docstring called `step_rule = "value"` a fallback for an
unusable Wald error; it is not. The step is `max(wald/4, |best| * 1e-3)` and the floor
wins whenever it is larger, which it does on **2 of these 6 parameters, all of which have
finite positive Wald errors**. Behaviour right, description wrong. The description is
fixed and the test now asserts `step_rule != "box_fallback"` rather than the one I had
assumed.

**CORRECTED 2026-09-06 by item 16.** This item read that assertion as *"the box never
entered, which is the direct statement that the defect could not act here"*. It is not:
`box_fallback` belongs to the CORRECTED routine and needs both `wald/4` and `|best|·1e-3`
to be zero, while the DEFECT acted whenever `span/5000 > wald/4` — true on **five of
these six parameters**. The conclusion survives and the reason does not. `n`'s box is
`2.0`, so its `(span/5000)/(wald/4)` is **1.09** — the floor was applied to the flow index
too and happened to be the right size — against **11 770** for `k`. The ratio explains
the three orders of magnitude this item could only assert.

**G15.4 puts two error terms in the same units for the first time here.** `rmse_log` is a
relative viscosity error and `d ln F / d ln k` is ~1, so they are comparable as
percentages: the hull is 7.5x the misfit on the 0.5 wt% gel and only **2.3x** on the
0.1 wt%. On that gel the error bar and the form error are the same order, and the page now
says the quoted interval is a **parameter** interval, not a total error.

**Two gates were named UNFALSIFIABLE in the pre-registration and are reported without
being counted:** the ordering reversal survives (the hulls were already disjoint by a
factor of three against a 3.9x effect, and narrowing can only help) and the nominal forces
are unchanged.

**Still untested:** the "an open interval is never printed as a number" branch has now
failed to fire in BOTH case studies. It is unexercised code, recorded as such.

### 16. DONE 2026-09-06 (all four gates PASS, after G16.1 failed TWICE on my own instruments): which numbers were exposed, and by how much

`scripts/exposure_scan.py` -> `data/exposure_scan.md` / `.json`, plus
`scripts/profile_interval_regrade.py` -> `data/profile_interval_prefix.md` / `.json`.
Gate pre-registered in `data/exposure_prereg.md`, committed at `1e33b56` before the scan
was run. Pinned in `tests/test_profiling_wrapper.py`.

**The question.** The same defect was **large** in the registry, **irrelevant** to the
Carbopol needle forces (0.008 %), and **irrelevant** to the blood family comparison for a
third and different reason. Three studies, one bug, three answers, and what decides which
you get was not written anywhere a caller would look.

**The answer is a RATIO, not a boolean:**

```
exposure = (span / 5000) / (wald / 4),   span = |bound - best|
```

`exposure = 1.09` means the floor was applied and happened to be the right size;
`exposure = 11 770` means the first step overshot by four orders of magnitude. **Both
ends, always** — a viscosity at `best = 21` has `span_hi/5000 = 200` against
`span_lo/5000 = 0.0043`.

| gate | result |
|---|---|
| G16.1 the grid signature implies the box floor won | PASS, **0** counterexamples over 609 ends (201 on the grid) |
| G16.2 the converse is NOT claimed | 8 inflated but not detectably on the grid — a limit of the DETECTOR |
| G16.3 the ratio tracks the overreach | Spearman **0.784** (p = 5e-43, n = 201) |
| G16.4 exposure surfaced where a reader will find it | every sidecar carries `wald_stderr`, `exposure_hi`, `exposure_lo` |

**G16.1 FAILED TWICE and both failures were in my own measuring code.**

1. **15 counterexamples.** The scan computed one ratio per PARAMETER from the upper bound
   and compared it against grid points detected on either end. 14 of the 15 were
   lower-end points graded against the upper ratio. Fixed by working per END.
2. **18 counterexamples, and this one was upstream.** The grid detector asked only whether
   `m = delta * 5.12e6 / span` is an INTEGER. With `_REFINEMENTS = 10` a bisection can
   only return `m` in `[1, 1024]`; `m > 1024` is not expressible. It therefore accepted
   any `delta/span` that happened to be round — **including every end with `delta == span`,
   i.e. `m = 5 120 000`, the endpoint sitting AT the box bound and not at a crossing at
   all.** Five of the eighteen were exactly that.

**It corrects item 13's own headline: 201 of 296, not 203.** The worst case is unchanged
at 1639x (`m = 251`) and **G13.2 is untouched** — still `n = 70`, still a median 710x
`f_crit` — because every rejected end was on a `cross` or `carreau_yasuda` record.
Addendum B of `data/profile_interval_prereg.md` records it; the pre-registration itself is
not edited.

**This is a change to a measuring instrument that makes a failing gate pass**, so the
justification has to stand alone, and it does: the range comes from `_REFINEMENTS`, not
from the gate, and would be the right bound if G16.1 had never been written. The
corroboration is the rejected set itself.

**A reproducibility defect surfaced with it.** `profile_interval_diagnosis.py` reads the
committed sidecars, and item 13 regenerated those sidecars **in the same commit that first
wrote the artifact** — so from `d38bc77` onward nobody could reproduce it, and running it
today finds zero. The rows are preserved verbatim as
`data/profile_interval_prefix_rows.json` and re-graded by
`scripts/profile_interval_regrade.py`; the live scan is relabelled as what it now is, a
**regression check that must report zero**. Both consumers were re-pointed, because
reading the live scan would have handed them an empty sample — a gate that passes by
vacuum.

**It also corrects an item-15 sentence.** Item 15 said `box_fallback` firing on none of
the Carbopol parameters is "the direct statement that the defect could not act here".
`box_fallback` is a property of the CORRECTED routine (it needs both `wald/4` and
`|best|·1e-3` to be zero); the DEFECT acted whenever `span/5000 > wald/4`, which was true
on **five of the six**. Item 15's conclusion survives and its reason does not: `n`'s box
is `2.0`, so its `exposure` is **1.09** while `k`'s is **11 770**. The ratio explains the
three orders of magnitude that the boolean could only assert.

### 18. DONE 2026-09-07 (all seven gates PASS; G18.6' added beside a bar that was close to unfalsifiable): the needle-bore exponent

`scripts/needle_bore.py` -> `data/needle_bore.md` / `.json`. Gate pre-registered in
`data/needle_bore_prereg.md`, committed at `8f4fe9c` before any code. Library code in
`src/biorheo/apps/design_space.py`. Pinned in `tests/test_design_space.py`.

**The quantity is `b = d ln F / d ln R` at fixed flow rate** — the exponent a device group
implicitly assumes is `-4` when it sizes a needle from a Newtonian calculation. Two closed
forms anchor it without touching the force solver: `newtonian` gives `-4` exactly and
`power_law` gives `-(3n+1)` exactly, and the registry has both null families on the same
20 donor curves as the fitted ones (item 11's records earning their keep).

| gate | result |
|---|---|
| G18.1 `b = -4` on the 20 `newtonian` records | PASS 20/20, worst 4.66e-10 |
| G18.2 `b = -(3n+1)` on the 20 `power_law` records | PASS 20/20, worst 2.66e-10 |
| G18.3 the discriminating half — not `-4` for everything | PASS, 36 records past 0.1, worst 1.867 |
| G18.4 no record steeper than Newtonian | PASS 115/115 |
| G18.5 SPLIT prediction on the Newtonian rule's error | PASS, both halves as pre-registered |
| G18.6 derating factor below 1 | PASS 91/91 — **and nearly unfalsifiable, see below** |
| G18.6' zero-width intervals give exactly 1, real ones less | PASS |
| G18.7 three error terms ranked | reported |

**G18.5 was pre-registered as a SPLIT and came out that way.** A Newtonian `R^-4` sizing
rule over-predicts the force at a 30 % finer bore by **under 1.24 % on 20/20 blood donors**
and by **47.9–93.9 % on 15/15 Carbopol records**. Whether `R^-4` is safe is a property of
the fluid AND the operating point: blood at that wall rate sits on its infinite-shear
plateau, so it is Newtonian there however shear-thinning it is at rheometer rates. The
blood half's FAILURE was predicted from item 14's `|d ln F / d ln eta_c| = 0.995`.

**G18.6 as pre-registered can essentially only be failed by a fit that locally maximises
the force**, because the upper envelope is a max over profile ends and the nominal force is
not one of them. A routine returning a constant 0.5 would have passed it. The bar and its
number stand; **G18.6' is the stricter test beside it** and pins both ends the bar cannot.

**The identifiability split is a TAIL, not a shift, and the first draft of that section
said the opposite.** Medians 0.796 (all-measurable) vs 0.866 (has a non-measurable
verdict); worst 0.123 vs 7.169e-05. An identifiability verdict does not predict how wide a
typical error bar is — it predicts whether the error bar can be catastrophic. The first
version restricted to closed intervals, which empties the very group being compared, and
had prose written around the empty result.

**A venture constraint surfaced here and it is not an item detail: ALL 115 records are
`commercial_use: false`** — 15 Carbopol under CC BY-NC-SA, 100 blood from a workbook that
states no licence at all. Not one curve in the registry can back a paid study. This does
not block a paid study on a CLIENT's curve; it blocks shipping a public demo on a named
fluid. `tests/test_report.py` pins it so it fails the day it stops being true.

### 19. DONE 2026-09-07: the deliverable — the page a client receives

`src/biorheo/apps/report.py` (`DeviceSpec`, `build`, `render`) plus four worked pages in
`case_studies/deliverable_example/`. Tests: `tests/test_report.py` (20 fast, 1 slow),
`tests/test_deliverable_example.py` (9 fast, 1 slow), `tests/test_offer_claims.py` (13).

It computes **nothing new**. Every number comes from a function anchored elsewhere. Its job
is ordering and **refusing**, and each refusal changes the HEADLINE rather than being
appended — a caveat in a footer is a caveat nobody reads:

| flag | how the headline changes |
|---|---|
| `no_sidecar` | "a POINT ESTIMATE with **no uncertainty**" |
| `open_interval` | "a **BOUND, not a limit**", naming the parameters |
| `hull_unbounded_above` | the quoted interval is a **LOWER BOUND** on the uncertainty |
| `non_commercial` | stamped **above** the answer |
| `extrapolated` | the distance in decades, above the fold |

**Two defects in already-shipped code, both found by the new tests.**

1. **`sensitivities()` stepped OUTSIDE the model box and raised** for a parameter sitting
   on its bound — `carreau_yasuda`'s `a` is fitted at 10 on several donors and a 1e-4
   relative step put it at 10.001. It now clips both offsets into the box and divides by
   **the span actually spent**, degrading to a one-sided difference at a bound. The 40
   shipped sensitivity values move by at most **1.061e-12** (`p_up - p_dn` is not exactly
   `2h` in floating point) and no verdict moves. That was checked by recomputing against
   `uq.json` with the comparison COUNT asserted — the first attempt compared an empty set
   because the key is `sensitivities`, not `sensitivity`, and printed a clean 0.0.
2. **An interval end can be a fluid nothing can deliver.** Donor B's `carreau_yasuda`
   `eta_inf` sits on its 1e6 Pa s bound, where 1 mL/min needs more than 1e12 Pa. Dropping
   that end would NARROW the interval — the direction that flatters the answer — so the
   hull end is reported open.

**The clamp branch has never fired on this registry** and there is a test that says so
rather than a test that pretends to exercise it.

`test_the_four_pages_demonstrate_four_different_things` first asserted four distinct FLAG
sets and failed at three. The criterion was wrong, not the pages: `01_clean` and
`04_thinning` differ on the Newtonian-rule verdict, not on flags. Corrected to
(flags, bore verdict); not weakened to three.

`../biofluids_venture/OFFER.md` is the sales document built on this, and
`tests/test_offer_claims.py` recomputes **every number in it** from the artifacts.

### 17. THE CASE IS WRITTEN 2026-09-08; THE DECISION IS THE FOUNDER'S (D5) AND IS STILL OPEN

`scripts/dstfit_upstream_proposal.py` -> `data/dstfit_upstream_proposal.md`, pinned by
`tests/test_dstfit_proposal.py` (8 tests). **`dst_fit` is untouched** — a pathspec-scoped
`git status` on it is empty and `test_gate_0_1_dstfit_unchanged.py` still passes.

**The number four items had deferred: the blast radius on `dst_fit`'s OWN results is
ZERO, and that is measured, not argued.** Item 16's exposure ratio says the box floor wins
exactly when `wald < span/1250`, which depends only on a parameter's declared box and the
precision the fit reports — no run needed.

- **No shipped `dst_fit` script passes `profile=` to `assess()`** (0 of them). The only
  shipped consumer of the profile walk is `scripts/gate_n7.py`, and it profiles **`phi_m`
  only**. `sigma_star` is profiled in a *test*, and a test is not a published number.
- **`phi_m`'s box is `[0.3, 0.639]`** and `gate_n7.py` quotes `phi_m = 0.5876 ± 0.0013`.
  Exposure is **0.177** at the lower end and **0.0316** at the upper — below 1 by factors
  of 6 and 32, so the box floor never sets the first step.
- **`sigma_star`'s box is `[1e-4, 1e5]`**, so `span/5000` is 20 Pa on a parameter of order
  44 Pa — a first step of ~45 % of the parameter itself. Nothing ships that today. **That
  is the risk the fix buys against: not a wrong number now, but a silently wrong one the
  first time anyone profiles a wide-boxed parameter.**

Recommendation on the page: **fix it, but not because anything is currently wrong** — the
argument is entirely about the next caller and the failure mode is a wrong answer with no
error. If it proceeds: pre-register "unchanged" per parameter via the exposure ratio (the
table IS that pre-registration and predicts no movement),
`test_gate_0_1_dstfit_unchanged.py` moves in the same commit, and **biorheo's wrapper is
not deleted** — it carries the endpoint verification, the second nuisance start and the
`step_rule`/`exposure` reporting, none of which is in the one-line fix.

**Rejecting it is a defensible position and should be recorded on that page** so item 17
stops being reopened. What is not defensible is leaving it undecided for a sixth item.

### 17b. The original framing, kept for the record: `dst_fit` still ships the defect

Written for zero context.

**The state of play.** `biorheo` routes around the defect and `dst_fit` still has it.
Every item from 13 to 16 has deferred the upstream question with the same sentence — *"if
the fix goes upstream it needs its own justification, not this one"* — and four items
later that justification now exists in measured form and has still never been written
down as a decision.

**What is known, and it is enough to decide on.**

- The endpoints are not confidence limits: a median **710x `f_crit`** over a 70-end
  pre-registered sample (item 13, G13.2).
- The corrected routine matches a **closed form** that touches none of its code, to
  **2.4e-9** on 20 records (G13.4).
- The blast radius is computable per parameter and is now a field in every sidecar
  (item 16, G16.4).
- `dst_fit`'s own suite is green and its one figure is byte-identical; nothing in it
  depends on the profile step floor being sized from the box (`tests/test_gate_0_1_dstfit_unchanged.py`).

**What to do, and the first step is NOT to write code.**

1. **Decide whether this is a decision an agent may make.** `dst_fit` is the founder's own
   package and the repo rule says it is consumed unchanged. Item 17 is a *proposal* to be
   read, not a patch to be landed. **Write the case; do not merge it.**
2. If it proceeds: the fix is one line in `dst_fit.profiling.profile_likelihood` (size the
   floor from the parameter rather than its box) plus the endpoint verification item 13
   added, and `tests/test_gate_0_1_dstfit_unchanged.py` **must be updated in the same
   commit** — it pins a figure hash that will move, and a gate whose whole purpose is to
   detect change cannot be left to fail silently.
3. **Pre-register what "unchanged results" means before touching it.** `dst_fit`'s own
   published numbers (`phim_drift`, the steady-subset study) were computed with the
   defective walk. The honest question is not whether they change but by how much, and
   the exposure ratio computes that per parameter without running anything.
4. **Do not delete biorheo's wrapper if the upstream fix lands.** It carries the endpoint
   verification, the second nuisance start and the `step_rule` / `exposure` reporting, none
   of which is part of the one-line fix, and all of which caught something real.

**The counter-argument, which is not weak.** `dst_fit` is consumed unchanged precisely so
that biorheo cannot quietly reshape it, and the wrapper already gives biorheo correct
intervals. Patching upstream buys correctness for OTHER callers of `dst_fit` — of which
there is currently one, `dst_fit` itself. Write that down too.

### 21. DONE 2026-09-07 (nine gates PASS, G21.9 FAILED on a mis-specified bar): a temperature axis, without inventing one

`src/biorheo/constitutive/shift.py` + `scripts/temperature_shift.py` ->
`data/temperature_shift.md` / `.json`. Gate pre-registered in
`data/temperature_shift_prereg.md`, committed at `6059c80` BEFORE the module existed.
Pinned in `tests/test_shift.py` (70 tests).

`../biofluids_venture/OFFER.md` named "no temperature dependence" as the single largest
missing capability for the first product. This closes it **narrowly**, and the narrowness
is the point: the layer carries a shift the CALLER supplies and refuses to guess one.

```
eta(gd; T)  =  a * eta_ref(a * gd)        a > 1 is the colder, more viscous fluid
```

**The pre-registration derives the parameter map for all six families by hand, before any
code**, and records the consequence nobody chose: **the yield stress comes out INVARIANT in
both yield-stress families.** Casson's `sqrt(a)` outside cancels the `sqrt(1/a)` inside;
Herschel-Bulkley's `a` cancels the `1/a` in `tau_y/(a gd)`. That is the physically right
answer -- a yield stress is a stress scale set by structure, not a viscosity -- but it is a
*consequence* here rather than an imposition, and the opposite would have been a reason to
distrust the construction. G21.6 checks it with its discriminating half beside it (the
viscous parameter must move, or invariance is being reported by doing nothing).

| gate | result |
|---|---|
| G21.1 map reproduces `a*eta(a*gd)`, six families | PASS, 24 cases, worst 6.66e-16 |
| G21.2 `a = 1` is bitwise identity | PASS 6/6 |
| G21.3 the shift composes | PASS, worst 2.22e-16 |
| G21.4 discriminating half -- the map is not the identity | PASS, smallest move 14.3 % |
| G21.5 the envelope moves and the sidecar is dropped | PASS |
| G21.6 `tau_y` invariant, viscous parameter not | PASS |
| G21.7 Newtonian force scales exactly as `a` | PASS, worst 1.84e-12 |
| G21.8 `structural_steady` refused for the RIGHT reason | PASS |
| **G21.9 collapse residual < 1e-10 on a synthetic pair** | **FAILED at 1.442e-05** |
| **G21.9' that residual CONVERGES at the interpolation order** | PASS, order 1.832 |

**G21.9's bar was MIS-SPECIFIED, not failed, and it has not been moved.** It demanded exact
recovery from a routine that compares against a LINEAR INTERPOLATION of the reference curve
in `(log gd, log eta)`, so the residual is bounded below by the interpolation error of
whatever grid it is handed -- **no correct implementation can beat it.** The number stays on
the page. G21.9' is stricter in kind rather than looser in size: refine the grid 8x and the
residual must fall at the interpolation order (1.4e-05 -> 3.2e-07, observed order **1.832**
against 2). A wrong parameter map, a broken objective or a genuine shape mismatch all leave
a residual that does not fall with `h` at all, so it cannot be passed by being merely small.

**Two things a caller must not miss, both enforced rather than documented.**
- **The envelope moves to `[gd_min/a, gd_max/a]`.** The shifted curve needs `a*gd` inside
  the reference envelope. A shift that left the envelope alone would silently widen the
  record's claim -- the same failure item 3 records for the OpenFOAM export.
- **The identifiability sidecar is DROPPED, not mapped.** A profile interval belongs to a
  fit against measured data and there is no data at the shifted temperature. A shifted
  record therefore reports as a POINT ESTIMATE in `apps/report.py`, which is checked.

`collapse_shift` returns the best `a` **and its residual and `n_used`**, because
time-temperature superposition is an assumption: a best-fitting `a` always exists, and two
curves differing only in SHAPE give one that means nothing (rms log residual 0.943 against
1.4e-05 for a real shift). No default activation energy is shipped for any fluid and a test
asserts `arrhenius_shift`'s `E_a` has no default.

### 22. DONE 2026-09-08, **THEN CORRECTED THE SAME DAY BY ITEM 23**: the registry's first licensable records

> **Read item 23 first.** This item's headline -- "G22.1' PASS, the strongest EXTERNAL
> anchor here, worst 0.0299" -- was graded against the paper's FIGURE CAPTION indices
> (its COMSOL inputs), not its Table 2 fit results. Against Table 2 the same bar FAILS
> at 0.0599. The records themselves did not move. Everything below is left as written.

`scripts/records_from_dmem.py` -> eight records + sidecars, `data/dmem_sints2024/records.md`
/ `.json`. Gate pre-registered in `data/dmem_record_prereg.md`, committed at `5825e58`
BEFORE the generator. Pinned in `tests/test_dmem_records.py` (16 tests).

**Eight `power_law` records, CC BY 4.0, `commercial_use: true`** -- the first in this
registry that may appear in a public demo, a webinar figure or a paid deliverable. The
other 115 remain `commercial_use: false`. Registry is now **123 records**.

**G22.1' is the strongest external anchor this repo has.** The source paper's authors fit
the same Ostwald-de Waele power law and published their flow indices in 2024. Fitting the
raw deposited files through `biorheo.identify` reproduces them on 4/4 conditions:

| condition | pooled `n` | paper `n` | dev |
|---|---|---|---|
| 10 % FBS, 20 C | 0.3195 | 0.32 | **0.0005** |
| 5 % FBS, 20 C | 0.3310 | 0.31 | 0.0210 |
| 1 % FBS, 20 C | 0.5699 | 0.54 | 0.0299 |
| 10 % FBS, 37 C | 0.5120 | 0.50 | 0.0120 |

Every other anchor in this repo compares biorheo against biorheo, a closed form, or a
solver. **This one compares it against a number somebody else published before this code
existed.**

**G22.1 FAILED at 0.065 against a 0.05 bar and the bar was NOT moved.** It compared each
INTERVAL's index against a value the paper publishes once per CONDITION; the per-interval
spread WITHIN a condition reaches 0.076, larger than the bar itself, so no fit could have
met it. G22.1' removes the interval-choice excuse instead of widening the tolerance. A
defect in the grader is recorded rather than quietly fixed: the first pooled run keyed on
the published index as well as the condition, splitting the 37 C condition and pooling only
its 5-point half.

**G22.6 was pre-registered as a PREDICTED FAILURE of a collapse, and it failed as
predicted.** A one-parameter shift maps `k -> a^n k` with `n` FIXED (item 21), and the
paper's indices are 0.32 at 20 C against 0.5 at 37 C for the same medium -- the two curves
differ in a parameter the shift cannot touch. A SMALL residual would have meant item 21's
shift or these fits were wrong. A positive control sits beside it.

**Item 11's failure was pre-empted, not repeated.** The pre-registration enumerated by grep,
BEFORE any record was written, every consumer depending on the registry's contents. Five
broke; four were on the list. **The fifth was not**: `test_registry` selected blood records
as "every record that is not `herschel_bulkley`", which was the same set only while blood
was the only non-HB family. `tests/conftest.py` now supplies a frozen 115-record registry
for gates CLOSED against it, and no existing gate changed its value.

**These records cannot carry an injectability claim** -- rates stop at 1 1/s against
1e4-1e5 in a needle -- and they say so in their own `notes`. `OFFER.md` says it too, so
they cannot be oversold as an injectability demo. What is still missing is a licensable
curve that reaches rheometer-to-needle rates.

### 23. DONE 2026-09-08 (G23.1 and G23.9 FAILED as predicted; no bar moved): the item-22 anchor was graded against the wrong published numbers

`scripts/dmem_anchor_recheck.py` -> `data/dmem_sints2024/anchor_recheck.md` / `.json`.
Pre-registered in `data/dmem_anchor_correction_prereg.md` at `d539ceb`, before the code, with
**every gate labelled `RESTATEMENT` or `BLIND`** because a review had already measured most
of them -- a restatement is not evidence the way a blind gate is, and the file says which is
which. Pinned in `tests/test_dmem_anchor_recheck.py` (15 tests).

**Found by four independent reviewer passes** (code review, physics referee, units audit,
numerics re-run). Two found the mis-anchoring without seeing each other; the numerics pass
then re-derived the corrected number from a from-scratch numpy log-log fit that never touches
`biorheo.identify`, matching to six decimals.

**The defect.** `SPEC` graded against `0.32, 0.31, 0.54, 0.50`, which are the paper's
**figure-caption** indices -- the values it fed to COMSOL, rounded and shared across
conditions. The paper also carries **Table 2, "Fit results for `K` and `n` values in (4) for
the experimental rheometry data shown in Fig. 1"**, one row per sample for exactly the four
deposited conditions, and it publishes `K` too, which item 22 never compared:

| condition | `SPEC` used | Table 2 says |
|---|---|---|
| 10 % FBS 20 C | 0.32 | 0.31 |
| 5 % FBS 20 C | 0.31 | **0.33** |
| 1 % FBS 20 C | 0.54 | **0.51** |
| 10 % FBS 37 C | 0.50 | **0.54** |

**All four disagree.** This is `lessons.md`'s 2026-09-01 entry -- *read the sample table
before comparing to it* -- for the second time.

**G23.1 FAILED at the unchanged 0.05**, worst 0.0599 on 1 % FBS. **G23.2**: neither single
correction rescues the old pass -- dropping the out-of-range 37 C interval alone gives 0.0651,
regrading against Table 2 alone gives 0.0599. **G23.3**: the flagship 0.0005 was a pooling
coincidence, that pooled index (0.3195) lying OUTSIDE the span of its own two interval fits
[0.2945, 0.3134], so the two intervals are mutually inconsistent under one power law.

**What survives is better anchored than what it replaces.** `K` and `n` trade off along the
data, so a deviation in `n` alone is uninterpretable. **G23.4 PASS**: the two parameter sets
compared **as viscosity curves** over each condition's measured range agree to a median
1.0-7.1 %, worst **19.03 %** -- inside the **20 %** repeatability the authors state for their
own measurements in S2 Fig. That is a scale the SOURCE declares, not one this repo invented.
`K` agrees to 6.7 % or better on all four and 3 of 4 round to the published value at its
printed precision.

**G23.9, added after the first run and labelled so, is stricter in kind than the bar beside
it**: does the published `n` fall inside biorheo's OWN profile interval on `n`? A fixed 0.05
bar treats a 5-point fit and a 30-point fit alike. Inside on 3 of 4; the exception is 1 % FBS,
[0.5380, 0.6017] against 0.51 -- **the only RESOLVED disagreement in the dataset**.

**Two of my own predictions were REFUTED and are in addendum A.** I predicted G23.1 and G23.4
while pooling every interval, then specified a pooling restricted to 0.01-1 1/s -- the same
knob mismatch the item exists to correct, committed in the file correcting it. Resolved by
reading what Table 2 says it fitted, *"the entire ensemble of data points"*, so the
restriction was item 22's invention and not the paper's. **G23.1b** reports its size: at 37 C
the three data choices give `n` = 0.4349 (5 points), **0.5393** (the 20-point interval alone,
0.0007 from the published 0.54) and 0.5120 (both). The restriction, not the fit, is what
disagreed there.

**And a test written to be able to fail caught a third error of mine**: "three of the four
values disagree" was in six files and it is FOUR of four. Two of them are the right number
against the WRONG ROW, which is how a total mismatch reads as a partial one. Addendum A.4.

**Five more defects fixed**, each with its original number kept: G22.1's `len(graded) == 6`
against 7 (could not pass whatever the fits returned; verdict unmoved, reason now true); the
spread argument's 0.1044 taken over a condition with one graded interval (0.0753 on fully
graded ones, both over the bar); verdicts read from the sidecar BEFORE it was written (blank
on a clean run, stale on a re-run); three `max()` calls a NaN can poison; and
`scripts/item11_gate_report.py` calling `load_all()`, so item 22's records would have entered
a CLOSED item-11 gate -- **G23.7: frozen, it reproduces the committed artifact BIT-FOR-BIT**,
so the artifact was right and the next run would have drifted silently.

`needle_bore.md`, `OFFER.md` and `arrhenius_shift` also corrected; see commit `4711dba`.

**G23.5**: no record moved, 0 of 8. **G23.6**: no file under `case_studies/` or `apps/`
references a DMEM record, so the correction reaches no shipped number.

### 20. DONE 2026-09-08 (superseded by item 22): find a licensable curve

**Why this matters more than it looks.** All 115 records are `commercial_use: false`
(item 18), so all four worked deliverable pages are stamped NOT LICENSABLE. A paid study on
a client's own curve is unaffected; a public demo, a webinar figure, a short course and a
conference slide are all blocked.

**Found: `data/dmem_sints2024/` — DMEM + FBS rheology, Zenodo 13828992,
`10.5281/zenodo.13828992`, CC BY 4.0.** Anton Paar MCR 502, cone-plate CP50-2, four
samples (1 %, 5 %, 10 % FBS at room temperature; 10 % at 37 °C), eight intervals, shear
rate / shear stress / viscosity / torque columns. Source paper arXiv:2409.13650. The 5 KB
zip is COMMITTED with per-file md5s in `SOURCE.md`.

**Read `data/dmem_sints2024/SOURCE.md` before touching it.** Three traps are recorded
there: each `.txt` holds MORE THAN ONE interval (an importer reading the first table
discards half the data — the same shape as item 7's concatenated ramps); one file opens at
2.24e-7 1/s with an apparent viscosity of 8710 Pa s, which is the instrument's zero-rate
point and not a measurement; and **"room temperature" is not a number**, so the envelope's
`temperature_C` cannot be filled without inventing one.

**A hypothesis of mine was REFUTED by the scan and is recorded rather than deleted.**
`scripts/dmem_candidate_scan.py` -> `candidate_scan.md` asked whether the affine-looking
stress column was a Bingham intercept masquerading as shear thinning. It is not: the
two-parameter power law fits better on **7 of 8 intervals**, same parameter count, same
error measure, `n` between 0.30 and 0.70. The paper's shear-thinning reading is the one the
data supports. What remains odd is the MAGNITUDE — the room-temperature intervals sit
14–32x water for a fluid that is ~99 % water — and the paper's own abstract says its
observed stresses "significantly exceed literature reports".

**What is left, in order.** (1) Read arXiv:2409.13650. (2) Resolve the DOI into
`doi_cache.json`. (3) Fit through `biorheo.identify` with intervals and verdicts, not the
scan's three hand fits. (4) Pre-register the gate BEFORE the generator, as every other
record has been. (5) Then decide whether this becomes the registry's first
`commercial_use: true` record.

**Even if it does, it cannot carry an injectability claim**: the rates stop at 1 1/s
(10 1/s on one interval) against 1e4–1e5 in a fine needle, i.e. four to five decades of
extrapolation against 1.42 for blood. Its value is that it is licensable, that it is a
biofluid, and that it is **the only data in the project with the same fluid at two
temperatures** — the one thing `constitutive.shift.collapse_shift` (item 21) has never been
run on outside a synthetic pair.

Other candidates found in the same search, all CC BY 4.0 and unexamined:
Zenodo 11389467 (2 wt% Carbopol 974p in propylene glycol — flow curve, stress relaxation;
a licensable counterpart to the NC Carbopol records), 8334762 (pipe rheology of aqueous
foams, with CSVs — pipe flow, so directly comparable to needle flow), 21868871 (620
water-based drilling-fluid formulations at six temperatures), 18961871 (detergent foams,
solution flow curves).

**How to search Zenodo:** its API returns HTTP 504 to a default `python-urllib`
User-Agent and works with an explicit one. Use `size=20` with backoff; `size=50` errored.

Do NOT digitise a figure to get around any of this, and do not add a record from a paper's
table without its raw curve. Both are standing rules.

### 24. DONE 2026-09-08 (G24.1 PASS -- the shift layer's first external anchor; G24.3 and G24.4 REFUTED my own predictions, no bar moved): the licensable yield-stress curve, at three temperatures

`scripts/records_from_carbopol_pg.py` -> two records + sidecars,
`data/carbopol_pg_vinutha2024/records.md` / `.json`. Pre-registered in
`data/carbopol_pg_prereg.md` at `b36c0fa`, before the generator, every gate BLIND. Pinned in
`tests/test_carbopol_pg_records.py` (12 tests). Read `SOURCE.md` for the four traps.

**THE RESULT: G24.1 PASSES, and it is the shift layer's first external anchor.**
`collapse_shift` fits the shift factor from the two flow curves alone and never sees a
viscosity; the deposit publishes the medium viscosities. `a_fit` = **0.13984** against
`a_pred` = 0.012/0.102 = **0.11765**, ratio **1.1887** against a bar of 1.30 set before the
number was known. Item 21 had only ever been exercised on a synthetic pair and on item 22's
collapse, which was pre-registered to FAIL.

**G24.5 and G24.6 confirm the rest of item 21's derived parameter map on real data**: `n`
invariant (0.5389 [0.5339, 0.5438] and 0.5325 [0.5277, 0.5373]) and `k` moving as `a^n`
(ratio 0.3461 against a predicted 0.3237). **G24.2** shows the test is INFORMATIVE rather
than assuming it: the high-rate residual rises **16.6x** when `a` is doubled, bar 3.

**G24.3 and G24.4 were pre-registered to PASS and FAILED. Neither bar was moved.**

**G24.3** predicted the low-rate half of the curve would be blind to the shift factor; it
moved **69 %** against a 10 % bar. The §1.1 algebra is right and my SPLIT was wrong -- it
divided at the geometric mean of the overlap, where this material is far above yield.
**G24.3ʺ** measures the real thing: the plateau (`sigma < 1.2 tau_y`) ends at **0.00398 1/s**
while the shiftable window starts at **0.0143 1/s**, so **the plateau lies entirely outside
the window a shift can act on**. The trap I pre-empted is real algebra and **inoperative on
this dataset**, which is also why G24.2's separation is so large. *A nuisance you pre-empt
may not be able to act, and whether it can is a property of the data's range.* A second
defect found while diagnosing: the three probes ran on DIFFERENT point sets (n = 38, 45, 30),
so the ratio was partly a change of sample; both gates now also report the intersection and
both verdicts are unchanged.

**G24.4** predicted the yield stress would be invariant, as item 21 derives and the authors
state. It is not: **8.855 [8.682, 9.027] against 9.176 [9.049, 9.303] Pa**, intervals missing
by 0.022 Pa, point estimates differing by 3.56 %. **The scale is reported beside it**: a
profile interval is statistical and carries no protocol-to-protocol term, and G24.8 measures
that term at **2.93 %** -- comparable to the difference being called significant. Consistent
as a material claim, inconsistent as a statistical one, and the page says both.

**G24.8 is the project's FIRST INDEPENDENT REPLICATE for any record.** The cessation files
carry a second, unannounced steady flow curve (`preshear-stress` at 16 rates per
temperature); it agrees with `Fig 1(a)` to a median 2.5 % at 10 C and 2.9 % at 50 C. G15.4
had to approximate this comparison with a model misfit because none existed.

**G24.7 refuses the 20 C block on four RECOMPUTED signals, and caught a defect of mine
first**: numpy 2's `repr(np.float64(1330.0))` is `"np.float64(1330.0)"`, so the
significant-figure counter counted the wrapper and returned 0/67 instead of 62/67. The gate
failed for exactly the reason it exists.

**G24.9** refuses to grade a Herschel-Bulkley `tau_y` against the published 8.9 Pa, which is
a fit to the three-component model this package does not implement.

⚠ **The pre-registration's own consumer enumeration was incomplete, for the third time in
this repo.** §4 tabulated five and missed a sixth, because the grep looked for record counts
and for `load_all` while that one selects on `commercial_ok`. It broke loudly. Fixed in kind:
`test_registry.py` now selects the Carbopol dataset **by id** rather than as "every
`herschel_bulkley` record" (item 22's defect from the other side), `test_openfoam_export.py`'s
closed item-3 gate is frozen, and the licensable set is asserted EXACTLY over both families.

Registry **123 -> 125 records, 10 licensable**.

### 24 (background, kept): why this dataset was taken first

**Zenodo 11389467, CC BY 4.0**, 2 wt% Carbopol 974p in propylene glycol, Anton Paar MCR 300,
1 deg cone, **both surfaces roughened**, Peltier hood, flow curves at **10 / 20 / 50 deg C**.
Paper arXiv:2312.00251. All five experimental files committed with md5s.

**Why it is worth more than the DMEM records.** It is a yield-stress microgel -- the same
family as the 15 CC BY-NC-SA Carbopol records the injectability case study ships, so it is
the licensable counterpart to the family this project actually demonstrates on -- and it
reaches **1e3 1/s**, one to two decades from a needle wall rate against the DMEM records'
four to five.

**And it is the first real test bed for `constitutive/shift.py`.** The authors chose
propylene glycol because its viscosity moves about tenfold over 10-50 deg C, so the material
is the same at all three temperatures and only the suspending medium changes. Two numbers
follow, both published in the deposit's own named columns and corroborated by the paper's
Methods:

| | 10 deg C | 20 deg C | 50 deg C |
|---|---|---|---|
| medium viscosity, Pa s | 0.102 | 0.055 | 0.012 |
| yield stress, Pa | 8.9 | 8.9 | 8.9 |

So **the shift factor is PREDICTED, not fitted**: if temperature acts only through the medium
viscosity then `a = eta_s(T)/eta_s(T_ref)`, which is `shift.shift_from_two_viscosities`, while
`shift.collapse_shift` fits `a` from the flow curves alone and never sees `eta_s`. Comparing
the two is a check of a FITTED quantity against a number from outside this repo. Item 21's
derivation and the authors' own sentence -- *"The yield stress is independent of the solvent
viscosity"* -- are the same statement reached independently, and G22.6 tested only the
NEGATIVE case. This is the positive one.

**The four things that had to be settled before the generator, and were:**

1. **Do not build a record from the 20 deg C block.** `SOURCE.md` trap 2: it is the only one
   whose rates are off the programmed logarithmic grid (3 of 67 against 37 of 48 and 50 of
   60), the only one whose **rate goes backwards** inside a decreasing ramp (32.4 -> 33.6 and
   24.4 -> 26.4 1/s), the only one rounded to 3 significant figures (62 of 67 points), and the
   only one whose independent replicate is offset one-sidedly. Four signals, no one of which
   would be enough alone. The shift test runs 10 -> 50 deg C, which is the larger viscosity
   contrast anyway.
2. **Use the second flow curve nobody has noticed.** Trap 3: both cessation files carry a
   `preshear-stress` column at 16 rates at each temperature -- a steady flow curve of the
   same material on the same instrument under a different protocol. It agrees with
   `Fig 1(a)` to a median 2.5 % at 10 deg C and 2.9 % at 50 deg C. **That is the only
   independent replicate this project has for any record**, and it is the scale a fitted
   parameter's confidence interval should be compared against -- the thing G15.4 had to
   approximate with a model misfit.
3. **Do not grade a Herschel-Bulkley `tau_y` against the published 8.9 Pa.** That is a fit to
   the *three-component model* (`sigma = sigma_y + k1 gd^0.5 + eta_bg gd`), which this repo
   does not implement. A yield stress is a fit-objective convention before it is a material
   constant -- Horner's Casson `tau_y`, and now item 23's flow index, both say so.
4. **Pre-register the collapse test's DISCRIMINATING subset.** Below yield both curves tend
   to `sigma_y`, so `eta -> sigma_y/gd` and `a*eta_ref(a*gd) = sigma_y/gd` for **any** `a`:
   the low-rate half of the curve contributes zero residual whatever the shift factor is. A
   residual over the whole curve is therefore dominated by points that cannot discriminate.
   Report the residual on the high-rate half as well, and report how much it moves when `a`
   is wrong by a factor of two, or the test will pass for a reason that is not the physics.

Rate ranges, in case a later session needs them without re-reading: 10 deg C, 48 points,
2.0e-3 to 1.0e2 1/s; 50 deg C, 60 points, 1.26e-3 to 1.0e3 1/s.

### 25. DONE 2026-09-08 (six gates PASS, G25.1 FAILED as predicted): THE WEDGE FLUID

`case_studies/mab_injectability/run.py` -> `results.md` / `.json`. Pre-registered in
`data/mab_injectability_prereg.md`; provenance and the source's own words in
`data/mab_prasnikar2025/SOURCE.md`. Pinned in `tests/test_mab_injectability.py` (12 tests).

**The first result in this project computed on the fluid the venture is actually about.**
Everything else here is human blood, a Carbopol microgel, or cell-culture medium.

Source: **Zenodo 14409284, CC BY 4.0**, supporting Table S1 of Prašnikar et al.,
*Int. J. Pharm.* **674** (2025) 125451 (open access). **58 formulations of one model
monoclonal antibody at 159.3–212.2 mg/mL**, viscosity 13.57–47.07 mPa·s each with an SD,
measured on a RheoSense m-VROC at 25 °C.

**THE HEADLINE.** At **matched** mAb concentration the best and worst excipient differ by
**7.33x** in viscosity, which through the concentration law is worth **62.6 mg/mL of extra
mAb at the same viscosity**. The raw spread across the table is 3.47x and 39.1 mg/mL and it
is the **wrong number** — these formulations span 159–212 mg/mL and viscosity is exponential
in concentration, so the raw spread mixes the excipient lever with the concentration lever.
The confound was pre-registered and controlled, and controlling it moved the answer 1.6x in
the direction that matters: **the naive spreadsheet number understates what the excipient
buys.** Viscosity doubles every **21.8 mg/mL**.

**G25.4: the needle extrapolates only 0.74 DECADES** past what the source measured — the
least of any fluid here, against blood's 1.42 and the DMEM records' four to five. The source
ran a rate sweep at **960–3329 1/s and categorized the formulations as Newtonian**, so the
one-parameter model is *its* verdict rather than our convenience, and an m-VROC measures near
needle rates by construction. That is why this dataset was worth chasing.

**G25.3: the median formulation needs 7.1x the plunger force of the median blood donor**
(1.98 N against 0.28 N) on the **same** needle, syringe and flow rate — the blood study's
device, read from its module rather than retyped, so the wedge fluid and the demonstration
fluid finally sit on one axis. **No injectable-force limit is quoted from any standard.**

**G25.6: 1525 of 1653 pairs (92.3 %) separate at 2 SD**, against item 14's 36 of 190 (18.9 %)
for the blood donors. A formulation decision made on these numbers is supported by them in a
way a patient-ranking decision on the blood data was not.

**G25.1 FAILED as predicted, and that is the point.** The source's own
`viscosity reduction [%]` column does not reproduce: back-solving its Mooney baseline and
fitting one exponential gives **R² = 0.4575** on the 49 rows of the main series, where a
single calibration curve would collapse to numerical precision. So the study uses the
**measured viscosities** and never that column. What survives is the **slope** (0.03131 and
0.03233 mL/mg on the two series, agreeing to 3.2 %), and it is what makes the mg/mL
conversion possible. There are **two series** and only the labels say so — nine rows carry a
`200 mAb` suffix; pooling them is what takes the back-solve from R² = 0.945 to 0.458.

**G25.7 demonstrates a SCHEMA BOUNDARY rather than asserting one.** These cannot become
registry records: a one-parameter Newtonian model on a single viscosity has zero residual
degrees of freedom and `fit_rate_curve` raises `ValueError: 1 points cannot determine 1
parameters`. The registry is untouched at 125 records and a test asserts it. **State that
boundary plainly rather than hide it — a table of formulation viscosities is what most
arriving client data looks like**, and the study reads the deposited table directly on
purpose, which is the path a paid study would actually take.

**A defect a type checker caught before it ran**: `import run` from one case study imports
**itself**, because both case studies are called `run.py` and the first on `sys.path` wins.
Every device constant would have resolved to this module's own unset names. Loaded by path.

**`OFFER.md` also still carried item 22's corrected anchor claim in two places** and now
states the corrected result. Its gap 3 is rewritten as what actually remains — the word
**curve** — rather than deleted.

### 26. DONE 2026-09-08: the deliverable page for SCREEN data, and the first LICENSABLE worked example

`biorheo.apps.report.build_from_screen()`; worked example
`case_studies/deliverable_example/05_screen.md`. Pinned in `tests/test_report_screen.py`
(11 tests) and `tests/test_deliverable_example.py`.

**Why.** Item 25 found the registry cannot hold the shape most client data arrives in — one
measured viscosity per formulation, no flow curve — so the report app could not produce a
page for the majority of realistic inbound work, and all four shipped worked examples were
on non-commercial data. Both are now fixed.

The page differs from the fitted one in the three places where it would otherwise lie:

1. the interval is a **MEASUREMENT** interval from the replicate SD, named as such and
   explicitly **not** a parameter interval. A Newtonian force is exactly proportional to
   viscosity, so the relative error passes through unchanged and no propagation
   approximation is involved. **With no SD there is no interval at all** and the headline
   says so rather than showing a zero-width one.
2. the flag **replaces** `no_sidecar` instead of joining it. `no_sidecar` means somebody
   forgot to fit; screen data has nothing to fit, and the two call for different responses.
3. the prose stops saying **"fitted"**. The page had said "the fitted law", "the fitted
   envelope" and "Fluid record X" — all false here, and it would have shipped in the
   **client-facing** document.

The licence comes from the **caller**: whether a client's own data may appear in a public
artifact is the client's answer, never this package's. A test asserts both directions.

**Two of the deliverable-set tests broke loudly and both were true when written.**
`test_every_page_is_stamped_non_commercial` carried the note *"this fails the day a
licensable record is added, which is the day the demo set should be rebuilt on it"* — today
was that day; it is now a licence **split** test. And `test_every_page_is_an_extrapolation`
required more than a decade on every page, while the screen fluid extrapolates **0.74** —
the smallest in the set and the whole reason that dataset was chosen — so the bar is restated
as "every page is an extrapolation" plus an assertion of the screen page's **rank**, rather
than kept or dropped.

### 27. DONE 2026-09-08 (G27.2 and G27.4 PASS; G27.1 and G27.6 FAILED, no bar moved; FIVE of seven registered predictions REFUTED): the repo's first EXTERNAL tube-flow anchor

`scripts/foam_pipe_anchor.py` -> `data/foam_koponen2023/foam_pipe_anchor.md` + `.json`.
Pre-registration `data/foam_koponen2023/PREREG_item27.md` (written before the generator,
with addendum A). Source and traps `data/foam_koponen2023/SOURCE.md`. 24 tests in
`tests/test_item27_foam_pipe_anchor.py`, one of them slow.

**Why this and not another record.** `apps/injectability.py` names its own uncovered gap in
its module docstring -- *"ALSO NOT included, and not inertial: wall slip"* -- and that
sentence had never carried a number, because **every tube anchor in this repo was
analytic** (Buckingham-Reiner, the Casson closed form, Poiseuille, a Couette OpenFOAM run
against an analytic profile). Nothing had ever been compared against a measured
pressure-driven flow. The claim `OFFER.md` sells is that a model identified from one
measurement predicts a different geometry, and that claim had no experiment behind it.

**The data.** Zenodo `10.5281/zenodo.8334762`, **CC BY 4.0**, Koponen et al.,
*Chem. Eng. Sci.* (2023) `10.1016/j.ces.2023.119282`. `Set_2` is a full factorial: 3 pipe
radii (4.075 / 5.95 / 9.50 mm) x 3 PVOH concentrations x 5 foam densities x ~15 flow rates
= 660 rows, and every one of the 15 cells exists at all three radii. Each row carries the
pressure gradient AND **a slip velocity measured independently by high-speed imaging**,
which is what makes the item possible. It is a FOAM: PREREG section 4 lists what may not be
claimed from it, and a test greps the page for unqualified mentions of blood or a mAb.

**A correction the deposit forces, found before any gate was graded.** The flow-rate column
is measured at the **absolute pressure of column 8** (115-408 kPa, at the generator outlet)
while the pressure gradient and the sample density are at **~1 atm** in the measurement
span. Foam is compressible, so `Q_local = Q / [(1-phi) + phi P_atm/P_abs]` -- a factor of
**1.091-3.120, median 1.574**. It was not settled by argument but by a check that can fail:
**a flow with a plug cannot slip faster than its own mean velocity**, which the raw column
violates on **11 of 660 rows** (worst 1.403) and the corrected one on **0** (worst 0.916).

**Results.**
- **G27.1 FAILED, 2 of 5.** Against the paper's **Table 2** (its fit table, not a caption --
  item 23's lesson applied in advance), the pooled flow index agrees at phi = 0.70 and 0.85
  and misses by 3.1-5.7 of the paper's OWN sigma at the other three. Bar not moved.
  **G27.1'** compares the same thing as VISCOSITY CURVES using both published parameters
  (`B` and `m` trade off exactly as `k` and `n` do): median **8.1 %** over 660 points.
- **G27.1 grades the REDUCTION, not `identify`.** After Weissenberg-Rabinowitsch,
  `eta = tau_w/gd_w` with `gd_w = gd_a (3N+1)/(4N)` and `N` the log-log slope of the same
  points, so a power law fitted to `(gd_w, eta)` returns `n = N` to printed precision. Said
  on the page, because an anchor that shares an identity with its target cannot check it.
- **G27.2 PASS at 18.74 %** (bar 20 %), 210 rows: identify on the 19 mm pipe, predict the
  8.15 mm pipe through `pressure_for_flow_rate`. **G27.4** (same pipe) is 2.40 %, so the
  fit describes its own data and G27.2's error is extrapolation rather than misfit.
- **G27.3 REFUTED its own registered prediction and it is the consequential result.**
  Deleting the slip term made the prediction BETTER (10.58 % vs 18.74 %), separation
  **0.56x** against a 1.6 bar, so **G27.2's pass is recorded NOT INFORMATIVE about slip**.
- **G27.3' measures why, and it is arithmetic.** `4 v_s/R = c gd_a` with `c = v_s/v_avg`, so
  the slip correction MULTIPLIES the apparent rate by `1-c`, and a power law rescaled by a
  constant in `gd` has the same exponent. Measured at matched `tau_w`, `c` is
  **0.361 / 0.267 / 0.331** across the three pipes -- **not even monotone** -- and the part
  that survives, `(1-c)_big/(1-c)_small`, is **1.068**: 6.8 % against an 18.7 % prediction
  error. **Slip is sub-dominant for cross-diameter transfer HERE, by about 3x**, which is a
  property of this data and not a general licence.
- **G27.6 FAILED at 140 %** (bar 50 %) and **G27.6'** says how: the Mooney three-diameter
  slip construction runs **2.20x HIGH** and returns a **negative slip velocity on 25 of 144
  points** -- inadmissible, not merely wrong -- because the three points are not monotone in
  `1/R`. The paper reached the same conclusion by a different route and imaged the slip
  instead. **The number this adds is how far apart the two standard methods are on one
  dataset**, which matters because nobody images slip inside a needle.
- **G27.7**: `Set_1` and `Set_2` differ by a median **6.4 %** at 19 shared conditions, so
  they look like repeats. Registered as un-graded, and no gate was re-graded on it.
- **15 of 225 rows could not be predicted and 1 cell was lost**, counted on the page. Three
  cells reduce to `n` at or near its 0.05 bound, and `sigma = k gd^0.06` reaches the needed
  stress only past `gd = 1e12`. They are the 0.5 % PVOH foams at the lowest densities --
  the corner the paper itself flags for partial phase separation.

**Route note.** `NEXT_STEPS` item 20 recorded the Zenodo REST API as the search route. It
now returns **HTTP 504 to an explicit User-Agent as well**, from two independent clients.
**DataCite answers** and indexes the same records with the same licence metadata, plus
figshare, Dryad, 4TU and Mendeley Data. `scripts/datacite_candidate_search.py` is the
replacement; it found **138 unique dataset DOIs, 107 of them under a licence this repo may
build a public demo on**, written to `data/candidate_search_2026-09-08.md`. Zenodo's
landing pages and file downloads still work; only `/api/` is down.

### 28. DONE 2026-09-08 (P1, P2, P5 PASS; P3 and P4 REFUTED, no bar moved): the term that defaults to zero

`injection_force` takes `friction: float = 0.0` and **no case study in this repo supplies
it**, so every force ever published here -- including the 1.87 N headline on the
CLIENT-FACING `05_screen.md` -- is the hydrodynamic term alone. That is labelled correctly
on every page and it is not the finding. The finding is which shipped CONCLUSIONS are
properties of the assumption that the friction is zero.

`case_studies/device_friction/` (`prereg.md`, `run.py`, `results.md`, `results.json`,
`tests/test_device_friction.py`, 8 tests). It reads only the shipped artifacts and refits
nothing, so it cannot disturb what it audits.

**The algebra decides which questions are live, and it is stated in the prereg before any
gate.** With `F = h + f` and `f` constant: ordering, force DIFFERENCES and
deterministic-interval separability are invariant BY CONSTRUCTION. They are checked as
arithmetic and are deliberately NOT gates. Two things are not invariant.

- **Q1 -- a RATIO is.** G25.3's *"the median formulation needs 7.1x the plunger force of
  the median blood donor"* is **7.05x at f = 0, 2.33x at f = 1 N, 1.52x at f = 3 N**, and
  first falls below 2.0 at **f = 1.44 N**. The statement that survives any friction is the
  DIFFERENCE, **1.70 N more than blood**, not the multiple. P1 and P2 PASS.
- **Q2 -- separability is destroyed by friction VARIABILITY, which a constant cannot
  touch, and no page mentions it.** G25.6's 92.3 % of formulation pairs separable at 2 SD
  falls below 50 % at a syringe-to-syringe friction SD of **s\* = 0.1397 N** and below
  10 % at 0.3254 N. The blood cohort's 36/190 reaches **0 %** at s = 0.1 N (P5 PASS).
  **The deliverable is the THRESHOLD, not a verdict at an assumed friction** -- a threshold
  needs no external number to be correct, and a datasheet only says where on the curve a
  device sits. No friction magnitude is quoted anywhere in the study, and a test greps the
  page for that.

**P3 REFUTED**: predicted `s* < 0.05 N`, measured 0.1397 N -- right in direction, 2.8x the
bar. **P4 REFUTED, and the defect was in my own pre-registration**: it required
`mab_growth > 3` (holds, 4.56x) and `blood_growth < 1.3` (fails, 3.33x), and it specified
the two arms under two DIFFERENT operations -- quadrature on the mAb side, a linear `2s`
on each hull side for blood. Same shape as item 6's G7 and item 23's pooling. The graded
verdict stands; a matched-operation recomputation is reported beside it, labelled POST-HOC:
blood 2.54x against mAb 4.56x, so the mAb cohort is the more exposed by **1.80x**, about
half the separation the bar implied.

**Why the mAb screen is the exposed one**: its interval is a pure measurement repeatability,
a median **0.42 %** of the force (7.8 mN), while the blood interval already carries a
median **33.4 %** parameter hull. *The tighter the fluid measurement, the more of the
answer the device owns.*

**C1 is the reason any of it is believable**: at f = 0, s = 0 the study reproduces all
three shipped numbers exactly -- 7.052465, 1525/1653 and 36/190. **C4, the null that would
have made the item pointless (`s* > 1 N`, i.e. friction is irrelevant here), did not fire**
and is recorded either way.

**Where a real device sits, from a DOI-resolved source and nowhere else.** A literature
survey run the same day (`data/injection_force_terms_lit_2026-09-08.md`) resolved
**breakloose ~3 N silicone-free to ~9 N siliconised** for a 1 mL glass prefilled syringe
with a staked needle (Funke et al., *J. Pharm. Sci.* 2016, doi:10.1002/jps.23945) and a
**1.3-1.5 N glide** bracket from a trade article with **no DOI**, carried as UNVERIFIED.
Glide is the comparable term, because the case studies compute a steady flow and
breakloose is the start transient the module docstring already excludes. **At 1.3-1.5 N
the shipped 7.05x ratio is 2.07-1.95x.** That syringe is not the one either case study
specifies -- their 2.3 mm barrel is a CHOICE -- so it places a device on the curve and
does not correct any number.

**And the Q2 side has no comparator at all, which is a result.** The survey resolved
friction MEANS and **no friction SD or syringe-to-syringe repeatability for any device**.
The quantity that decides whether a formulation screen is decidable is unpublished, so
`s* = 0.1397 N` stands as a threshold with nothing to compare it against. It is also the
cheapest bench measurement in this whole programme: N empty-barrel injections on one
syringe lot.

**Open, and a measurement rather than a calculation**: friction is treated as additive and
FLUID-INDEPENDENT. That assumption is what makes the ordering claims survive. If plunger
friction depends on the fluid through the lubricant film, even the ordering is exposed and
nothing here can tell.

### 29–40. OPEN 2026-09-26: the in-silico work-order sheet — take them in wave order

Founder decision 2026-09-26: **no bench, in silico only** (`../biofluids_venture/09_INSILICO_PLAN_2026-09-26.md`).
Briefs, gates and `owns` sets are in **`WORK_ORDERS_2026-09-26.md`**; the still-open 09-08
packages B1–B6 are folded in there. Take ONE per session, in this order, and mark it here
DONE/FAILED/CLOSED with the numbers exactly as items 1–28 are:

| item | package | one line | status |
|---|---|---|---|
| 29 | D1 | retire dst_fit gate G4c with the reason (founder D5 taken); update `test_gate_0_1_dstfit_unchanged.py` expectation | **DONE 2026-09-28** — see below |
| 30 | B3 | audit every external anchor against its source's FIT TABLE (read-only) | **DONE 2026-09-28** — see below |
| 31 | N1 | `biorheo` CLI (`fit` / `report` / `screen` / `export`), wheel, CI; byte-identical page from a client-shaped CSV | OPEN |
| 32 | J1 | freeze φ_J = Z/(Z+2√3 w) as `polydisperse_jamming/phij_model.py` with a trained-domain guard; held-out seeds chosen before fitting | OPEN |
| 33 | B5 | referee items 24, 25, 26 (read-only) | OPEN |
| 34 | B1 | per-cohort friction-repeatability budget `s*` on the deliverable page | OPEN |
| 35 | B4 | omitted terms sized for the mAb cohort; reconcile the 3× Hagenbach vs end-correction disagreement | OPEN |
| 36 | N2 | L5 credibility manifest (V&V 40) as code; wrap mAb study and the bidet grid study; every factor defaults MISSING | OPEN |
| 37 | B2 | one licensable protein FLOW CURVE (or the figshare mucus deposit) as records — exclusive on `RECORD_IDS` | OPEN |
| 38 | N4 | device terms as explicit inputs with sensitivities; Arrhenius input for the shift; existing pages byte-identical at zero | OPEN (after 34, 35) |
| 39 | N3 | Fluent + COMSOL exporters with parse-back round trips and the regularisation-error curve | OPEN |
| 40 | J2 | `biorheo.apps.jamming` formulation-risk page (after 31, 32) | OPEN |
| 41 | L2 | capillary-rheometer analysis pipeline (Rabinowitsch, Bagley, Mooney with refusal, uncertainty budget) on the Koponen control | OPEN |
| 42 | L1 | high-shear rig BOM from vendor list prices (research only) | OPEN |
| then | B6, N5 (repeated), D2, J3 | see the sheet's wave 4 | OPEN |

Rules that hold for all of them are in the sheet's header and in "Rules that are specific to
this repo" above. Nothing on this list asks for a measurement; anything that would need one
is an input with a sensitivity, never a default.

#### 29 (D1). DONE 2026-09-28: G4c retired, its number kept, exit code pinned

`dst_fit/scripts/gate_g4.py` was run BEFORE the edit (rc=1, 6 s, last lines
`G4a PASS   G4b PASS   G4c FAIL` / `G4c: 1 of 6 within 2 sigma, needed 4`) and AFTER
(rc=0, `G4a PASS   G4b PASS   G4c RETIRED (would have read FAIL)`, with
`bar was 4 of 6; reads 1 of 6 within 2 sigma`). G4a (+0.0003) and G4b (-0.0151) are
unchanged. The check is still computed and printed; only the exit code no longer depends
on it. The bar constants (`G4_MIN_CONSISTENT = 4`, `G4_N_SIGMA = 2.0`, `G4_TOLERANCE =
0.02`) are untouched and a static test asserts they stay so. Reason, recorded in the
script docstring, `dst_fit/README.md` "Gate G4, as run" and `dst_fit/NEXT.md`: the rule
G4c encodes (cross-curve consistency of `phi_m`) is not satisfiable in any material,
because a flow curve determines only `phi/phi_J` (N4, N6, G4, exact gauge covariance);
it was added on top of the pre-registered G4a/G4b ahead of its run and tested the
degeneracy rather than the fitter. Founder decision D5, 2026-09-26.

`tests/test_gate_0_1_dstfit_unchanged.py` gained two fast tests: a static one (the
RETIRED marker is in the source AND the three bar constants are verbatim) and a dynamic
one (runs the gate, ~6 s: rc 0, the retained "1 of 6 ... would have been FAIL" line
present, "GATE G4 FAIL" absent). Both halves are needed: a later edit that made G4c pass
by touching the data or the bar would fail the second; one that dropped the number would
fail it too. The `README.md` "Status" sentence "eight pass, five fail" and the table
row "FAIL, 1 of 6" are kept verbatim with a RETIRED note beside them.

Not done: `dst_fit/NEXT.md` line 23 (the verify block's `# 1  <- ... G4c fails`
comment) still reads the OLD exit code — the edit was refused by the session's
permission classifier, not by any rule; the new "What is new" section in the same file
states the new exit code. One-line fix for the next session. No reviewer was run: no
physics or numerics changed, and the diff is a gate's exit-code policy plus prose.

#### 30 (B3). DONE 2026-09-28: every external comparison graded; three have NO fit table

`data/external_anchor_audit_2026-09-08.md`, ten rows, read-only. Sources read this
session: Prašnikar 2025 accepted PDF (Ljubljana repository, 13 pp), Vinutha 2024 arXiv v1
full text, Biermann 2025 full text (Europe PMC `PMC12251366`); MDPI and Elsevier both
refuse fetches. Verdicts: DMEM (item 23), Koponen and Horner are GRADED AGAINST A FIT
TABLE; the Carbopol-PG shift factor compares a MEASURED solvent viscosity; the
Carbopol-PG 8.9 Pa is CAPTION-ONLY (three-component `k1`, `eta_bg` unpublished) and was
correctly never graded; the mAb calibration `η = a·e^(b·c)` is fitted in R with a CI band
and **`a`, `b`, R² are printed nowhere**, so the 21.8 mg/mL doubling is unanchorable (prose
points imply 0.03785 vs the back-solved 0.03131, 1.21×, on "~150" rounding — a note, not a
grade); the 15 Carbopol ETD2020 records rest on the sibling `carbopol/fit_hb.py` range
check ("τ_y max 41.1 at 0.5 wt%, k 3.4–15.9, n 0.40–0.53 — reproduced") and the paper has
**no per-condition table**, only heatmaps A1/A2, so a range check is the ceiling; ISO 9626
withdrawn; OpenFOAM anchors analytic. **One defect found outside the anchors**:
`WORK_ORDERS_2026-09-26.md` N4.2 cites "the Funke 2016 glide bracket (doi:10.1002/jps.23945,
1.3–1.5 N)" — the 1.3–1.5 N is the no-DOI trade article's glide (UNVERIFIED on the page),
Funke publishes 3–9 N breakloose. **Before N4 runs, its owner must re-cite N4.2** (keep the
number, fix the source, or re-register on Funke's breakloose which `results.json` also
carries). Two one-sentence fixes ranked below it: the mAb page should say the authors'
`b` is unpublished; the Carbopol records' `notes` should carry the range check and the
sibling's "reproduced" should read "inside the published ranges". No number moved.

**Acted on the same day (commit `3f7798d`):** N4.2 got its marker plus Addendum A1 at the
foot of `WORK_ORDERS_2026-09-26.md` (the pre-registration itself is not rewritten; N4's
owner picks option (a) or (b) before running). The mAb page now carries the "no external
check" paragraph; `results.json` differed only in its `git_head` stamp, and
`test_the_page_says_the_calibration_exponent_is_unpublished` pins the sentence (52 tests
in the three affected files passed, 472 s). **Still open from the audit:** the Carbopol
records' `notes` range check (touches 15 generated YAMLs via `records_from_carbopol.py`
— regenerate and confirm byte-identity elsewhere) and the sibling `carbopol/fit_hb.py`
wording; neither is a biorheo package, both are one-line.


## What NOT to do

- **Do not quote a shipped injection force as the force a user feels.** Every one of them
  is the hydrodynamic term; `friction` defaults to 0 and no case study supplies it. Item 28
  measures what that costs: the G25.3 ratio is 7.05x at zero friction and 2.33x at 1 N.
- **Do not quote the 7.1x (or any force RATIO) without the friction caveat.** Quote the
  DIFFERENCE instead -- 1.70 N -- which is invariant under any constant friction.
- **Do not read item 28 as a claim about any syringe.** It quotes a threshold (s* = 0.1397 N)
  and no friction magnitude at all; a test greps the page to keep it that way.


- **Do not quote item 27's G27.2 pass as evidence that biorheo handles wall slip.** Its own
  pre-registered discriminator refused it: deleting the slip term makes the SAME prediction
  better (10.58 % vs 18.74 %, separation 0.56x against a 1.6 bar). G27.3' measures why --
  the slip fraction at matched wall stress is 0.361/0.267/0.331 across the three pipes, not
  even monotone, and the surviving part is 6.8 % against an 18.7 % error.
- **Do not read that as "slip does not matter".** It is `(1-c)` cancelling between
  identification and prediction because `c` happens not to depend on diameter here. The
  moment it does, the same arithmetic makes slip the leading term.
- **Do not use the Mooney multi-diameter slip construction without checking monotonicity in
  `1/R` first.** G27.6': 2.20x high on the median and 25 of 144 points NEGATIVE. The
  deposit's authors imaged the slip instead and said why in their section 2.
- **Do not reduce the Koponen deposit without the compressibility factor.** The flow rate is
  measured at 115-408 kPa and the pressure gradient at ~1 atm; the factor is a median 1.574
  and omitting it makes the imaged slip velocity exceed the mean velocity on 11 rows.
- **Do not add a foam record to the registry from item 27.** The reduced curve is DERIVED
  (slip subtraction + Weissenberg-Rabinowitsch), and the registry has never held a derived
  curve. A test asserts no `*foam*` record exists.
- **Do not describe item 27 as a rheometer-to-needle anchor.** The deposit's rheometer file
  is the PVOH SOLUTION, not the foam, so there is no rheometer curve of the flowing fluid
  anywhere in it. It is pipe-to-pipe: the same extrapolation in kind, a smaller one in
  degree.
- **Do not search Zenodo through its REST API.** It returns 504 to an explicit User-Agent
  as of 2026-09-08. Use `scripts/datacite_candidate_search.py`; Zenodo's landing pages and
  file downloads are fine.

- **Do not compare a screen page's interval with a fitted page's.** One is the measurement's
  repeatability, the other is a parameter hull. `build_from_screen` names both and they are
  not interchangeable.
- **Do not let the report call anything "fitted" on the screen path.** Nothing was. A test
  greps the rendered page for `fitted law`, `fitted envelope`, `identifiability verdict` and
  `profile push-through`.
- **Do not default a licence in `build_from_screen`.** The caller supplies it, because
  whether a client's data may appear in public is the client's answer.

- **Do not quote item 25's raw excipient spread (3.47x, 39 mg/mL).** It is confounded: the
  58 formulations span 159-212 mg/mL and viscosity is exponential in concentration. The
  answer is the concentration-normalised 7.33x / 62.6 mg/mL, and the raw number is printed
  beside it only to show the size of the confound.
- **Do not use the source's `viscosity reduction [%]` column.** G25.1: it does not reproduce
  from a single Mooney curve (R2 = 0.4575), so nothing may be built on it. The SLOPE does
  survive and is what the mg/mL conversion uses.
- **Do not try to make the mAb formulations into registry records.** One viscosity cannot
  determine a one-parameter model; G25.7 shows the exception. This is a schema boundary and
  also a description of most client data.
- **Do not quote an injectable-force limit** in the mAb study or anywhere else. The scale is
  this repo's own blood study on the same device. The retracted ISO 9626 comparison is the
  precedent.
- **Do not `import run` from inside a case study.** Both `case_studies/*/run.py` share a
  name, so a bare import silently imports the importing file. Load by path.

- **Do not use the Carbopol-in-propylene-glycol records for an injectability number.** They
  stop at 1e3 1/s against 1e4-1e5 at a needle wall, the medium is propylene glycol, and the
  sample is a rheology model system. They exist because they are LICENSABLE and because they
  carry a temperature axis. Their own `notes` say so.
- **Do not read G24.4's failure as "the shift layer is wrong".** `tau_y` moves 3.56 % between
  10 and 50 C while the independent replicate's own scatter is 2.93 %. The profile interval
  is a STATISTICAL interval on one curve and excludes the protocol term; the failure is real
  at that level and not at the level of the material. Quote both.
- **Do not reinstate G24.3 as posed.** Its "low-rate half" was the power-law region, not the
  plateau, and the plateau for this material lies entirely OUTSIDE the range a shift can act
  on (0.00398 vs 0.0143 1/s). The trap is real algebra and inoperative here; G24.3ʺ is the
  measurement that settles it.
- **Do not count significant figures off a numpy scalar's `repr`.** Under numpy 2 it is
  `"np.float64(1330.0)"` and the count includes the wrapper. Convert with `float()` first.
- **Do not trust a consumer enumeration that greps only for record counts and `load_all`.**
  Item 24's missed a sixth consumer that selects on `commercial_ok`. Grep for the PROPERTY a
  new record changes, not only for the shapes of the last defect.

- **Do not quote item 22's anchor numbers.** "Reproduces the published flow indices on 4/4
  conditions, worst 0.0299, 0.0005 on the flagship" is CORRECTED (item 23). Against the
  paper's Table 2 the same bar fails at 0.0599, and the 0.0005 was a pooling coincidence --
  that pooled index lies outside the span of the two interval fits it pools. What may be
  quoted is `data/dmem_sints2024/anchor_recheck.md`: agreement as a VISCOSITY CURVE to a
  median 1.0-7.1 %, worst 19.03 %, inside the 20 % repeatability the authors declare.
- **Do not grade against a number from a figure caption when the paper has a fit table.**
  A caption states what the authors FED a simulation; a fit table states what they measured.
  They differed on all four conditions here, and twice two of them were the right number
  against the wrong row, which is how a total mismatch reads as a partial one.
- **Do not compare a fitted `n` alone against a published `n`.** `K` and `n` trade off along
  the data, so a deviation in one is uninterpretable without the other. Compare the parameter
  sets as viscosity curves over the measured range, and scale the answer against a
  repeatability the SOURCE declares rather than a tolerance invented here.
- **Do not let a generator that writes a CLOSED gate's artifact call `load_all()`.** Item 22
  froze the TESTS and missed `scripts/item11_gate_report.py`; eight new records would have
  entered a gate whose page says 115. Use `CLOSED_GATE_FAMILIES` from `tests/conftest.py`.
- **Do not call `arrhenius_shift` with celsius.** It now refuses below 100 K. The only
  per-record temperature in the schema is `envelope.temperature_C`, so the wrong call is the
  natural one.

- Do not use the DMEM records for an injectability number. Their rates stop at 1 1/s
  against 1e4-1e5 in a needle: four to five decades. They exist because they are
  LICENSABLE, not because they suit the wedge.
- Do not grade a per-interval fit against a per-condition published value. That is G22.1's
  recorded failure; the spread WITHIN a condition is larger than the bar.
- Do not point a CLOSED gate at `load_all()`. Use `closed_gate_records` from
  `tests/conftest.py`, which is frozen at the 115 records of 2026-09-07. Item 11's failure
  was a shared list growing under a decided gate; item 22 pre-empted it and still found a
  fifth consumer nobody listed, because it selected blood as "not herschel_bulkley".
- Do not quote `b = d ln F / d ln R` as a material constant. It is a property of the
  fluid AND the operating point: blood returns -3.98 in a fine needle because it sits on
  its infinite-shear plateau THERE, not because it is nearly Newtonian.
- Do not read item 18's derating factor being below 1 as a result. That bar is close to
  unfalsifiable (the upper envelope is a max over profile ends); the SIZE is the result
  and G18.6' is what discriminates.
- Do not delete the clamp in `apps/report.py` because its test says the branch never
  fires. The test records untested code; it does not license removing it.
- Do not edit a number in `../biofluids_venture/OFFER.md` without rerunning
  `tests/test_offer_claims.py`. Every figure in that document is pinned to an artifact.
- Do not describe item 21 as "biorheo does temperature". It carries a shift factor the
  CALLER supplies. There is no temperature model, no activation energy for any fluid, and
  no record measured at more than one temperature.
- Do not carry an identifiability sidecar through a temperature shift. There is no fit at
  the shifted temperature, so there is no interval; a shifted record is a point estimate
  and `apps/report.py` says so.
- Do not "fix" G21.9 by widening its bar. It is recorded FAILED with its number because the
  bar demanded something no correct implementation can deliver; G21.9' is the replacement
  and it is stricter, not looser.

- Do not fit any model to the post-peak decay of the step sheet (syneresis, per the thesis).
- Do not reopen item 7 by asking again whether the LAOS sheets are densely enough sampled.
  They are; that was measured and it is not the constraint.
- Do not fit the LAOS sheets to the structural family as it stands: item 2 showed it does
  not reproduce the step rise's shape.
- Do not add a record from a paper's table without its raw curve.
- Do not put `si_fit` or any `impact_stress/materials.py` constant into the registry.
- Do not change `MAX_REPRO_DEV`, `RUNS_Z`, gate bars, `LAMINAR_RE`, `BAR_*`, `G7P_RELAX_U`
  or `G7P_END_TIME` in `scripts/openfoam_couette_anchor.py` to make a check pass.
  `G7P_*` are the configuration being MATCHED, not bars. Same for `G91_BAR`, `G92_BAR`,
  `G93_BAR`, `G101_BAR`, `G102_BAR` in the two seed-scatter scripts.
- Do not form a RELATIVE spread without checking what the denominator is. Item 9's
  retraction is the whole argument: a parameter converged to 1e-19 gave a 415 % move.
  Measure against a declared scale -- the parameter's own confidence interval, or a
  physical quantity -- before believing any ratio.
- Do not read a byte-identity or bit-identity claim off a MAXIMUM. Item 8's second
  correction: six gate scalars were bit-identical while 7 of 33 underlying rows moved.
- Do not reinstate the ISO 9626 needle-tolerance comparison without the standard in hand.
- Do not quote a flow rate from `flow_rate_casson` within about a part in 1e9 of yield.
- Do not add an `a != 2` refusal to the BirdCarreau export. v2606 honours `a`; that was
  decided by running the solver, and the refusal item 6 asked for would refuse something
  OpenFOAM can do.
- Do not cite `foamDictionary -entry nu` as evidence about OpenFOAM's Newtonian LAYOUT.
  It reads a top-level key out of a file biorheo wrote with a top-level key, using a
  generic dictionary reader; it would pass unchanged if the layout claim were false. It
  anchors bit-exactness, which is a different claim. G11.4b is the only discriminator.
- Do not add an `a != 2` refusal to the BirdCarreau export, and do not add a
  "structural_steady has no counterpart" message to the structural refusal. Both are
  premises this repo tested and refuted (items 6 and 11).
- Do not describe item 12's finding as "which constitutive family you pick moves the
  force". It is not symmetric: Casson fits these curves WORSE on 20/20 donors, and the
  effect tracks the misfit. The shipped family is the worse-fitting one.
- Do not read a sidecar `hi` as a limit without checking `hi_at_bound`. Item 13 made
  MORE ends open, not fewer: an end the walk could not verify as a crossing is now
  reported at the box bound, so it is a `bounded_only` statement and not a number.
- Do not quote a `force_with_uq` push-through hull as a 95 % confidence interval. It
  is a ONE-AT-A-TIME box hull and is wider than the joint region it is drawn from;
  the covariance-sampled interval beside it is a joint region but assumes the Wald
  normal. Name which one you are quoting. (The blanket ban on `force_with_uq` was
  lifted when item 13 closed; the reason it existed -- 201 of 296 closed upper ends
  were bisection-grid artefacts -- is fixed, not waived.)
- Do not expect a corrected interval to have changed item 12. It did not, and G13.9
  measured why: the bug inflated only the UPPER ends and the force difference is
  signed downward, so exactly one donor of twenty flips.
- Do not read "the donor spread is bigger than the family effect" as "which patient it is
  matters more". The first is a range against an effect size; the second is a claim about
  telling two patients apart, and item 14 measured that one directly: 36 of 190 pairs
  separate at 95 %, and two donors must be 17 places apart out of 20. A shipped page said
  the second for a day; the correction is on it, with the old wording quoted.
- Do not quote a `force_with_uq` interval without naming the object. The push-through hull
  is one-at-a-time and WIDER than the joint region; the covariance sample IS joint but
  assumes the Wald normal. Neither is simply better.
- Do not restate item 13's size without its condition. "The bug moved endpoints by 1639x"
  and "the bug moved these forces by 0.008 %" are both true, of the registry and of the
  Carbopol needle study. What decides which you get is the parameter's declared BOX width
  (the floor is `box/5000`), the number of points on the curve (a large Wald error means
  `wald/4` wins the `max`), and which parameter the downstream quantity is sensitive to.
- Do not call `step_rule = "value"` a fallback for an unusable Wald error. The step is
  `max(wald/4, |best| * 1e-3)` and that floor wins whenever it is larger -- on 2 of the 6
  Carbopol parameters, all with finite positive Wald errors.
- Do not read `box_fallback` never firing as "the defect could not act here". That branch
  belongs to the CORRECTED routine and needs BOTH terms to be zero; the defect acted
  whenever `span/5000 > wald/4`, which was true on five of those same six parameters. The
  quantity that answers it is `exposure_hi` / `exposure_lo` in the sidecar, and it is a
  ratio: 1.09 for a flow index bounded at 2.0, 11 770 for a consistency bounded at 1e6.
