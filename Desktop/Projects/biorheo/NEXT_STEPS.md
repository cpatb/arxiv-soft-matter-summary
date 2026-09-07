# NEXT_STEPS.md — how to continue `biorheo` in a fresh session (any model)

Read this, then `README.md`, then the plan ledger at
`../biofluids_venture/docs/superpowers/plans/2026-09-05-phase0-foundation.md`. Do not
reconstruct state from chat history. State lives in the Obsidian card
`Documents\Obsidian Vault\Status\Biofluids Venture.md` and in this file.

## Before touching anything

```
cd C:\Users\compl\Desktop\Projects\biorheo
python -m pytest -q -m "not slow"      # must print: 311 passed, 7 deselected, 1 xfailed
git -C C:\Users\compl log --oneline -3  # top should be c771970 or later
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
three guessed candidates: **203 of 296 closed upper ends** on the grid (this file used to
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
   MEDIAN narrowing is 1.0002x. The effect is a TAIL -- 203 of 296 upper ends -- not a
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
fixed and the test now asserts the invariant that matters — **the BOX never entered
(`step_rule != "box_fallback"`)** — rather than the one I had assumed.

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

### 16. NEXT ITEM, not started: nothing in this repo tells a caller whether ITS numbers were affected by item 13

Written for zero context.

**The problem, stated from what items 13, 14 and 15 measured.** The same bug was:

- **large** in the registry — 203 of 296 closed upper ends on the grid, worst 1639x the
  best-fit value, shipped endpoints at a median 710x `f_crit`;
- **irrelevant** to the Carbopol needle forces — 0.0081 %;
- **irrelevant** to the blood family comparison, for a third and different reason (it
  moved the upper ends and that comparison binds on the lower ones).

Three studies, one bug, three answers, and **the thing that decides which you get is not
written down anywhere a caller would look**. It is: how wide the parameter's declared box
is (the floor is `box/5000`), how many points the curve has (a large Wald error means
`wald/4` wins the `max`), and which parameter the downstream quantity is sensitive to.

**What to do.**

1. **Do not build a general "impact estimator".** That is the over-engineering trap this
   repo's rules warn about. Build the smallest thing that answers the question a caller
   actually asks, which is *"is the number I am about to quote affected?"*
2. The cheapest honest answer already exists per-profile and is not surfaced:
   `ProfileInterval.step_rule`. **`box_fallback` never fires and `value`/`wald` both mean
   the box was not consulted**, so a profile whose rule is not `box_fallback` is one the
   defect could not have acted on. Consider exposing that in the sidecar, where a reader
   can see it without re-running anything.
3. **Pre-register, and pick a gate that can fail.** Items 14 and 15 both had to have their
   obvious gate discarded as unfalsifiable. A candidate that can fail: *does `step_rule`
   actually predict which registry parameters were on the grid?* The diagnosis artifact
   (`data/profile_interval_diagnosis.json`) lists all 219 grid points by record and
   parameter, and the profiler can report its rule for each. **If they do not agree, the
   explanation on three pages of this repo is wrong** — and that is worth knowing before
   it is quoted again.
4. Note the asymmetry before predicting: `wald/4 > box/5000` means the box did not set the
   step, but a parameter can still land on the grid through a LATER step, because the step
   grows geometrically from wherever it starts. So the prediction is one-directional and
   should be stated that way.

**Do not** patch `dst_fit` as part of this, and do not restate any of the three numbers
above without its condition — "the bug moved things by 1639x" and "the bug moved things by
0.008 %" are both true and both useless without the clause that follows.

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
  lifted when item 13 closed; the reason it existed -- 203 of 296 closed upper ends
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
  Carbopol parameters, all with finite positive Wald errors. `box_fallback` is the branch
  that means the box was consulted, and it has never fired.
