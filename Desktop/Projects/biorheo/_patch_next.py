"""One-shot: rewrite NEXT_STEPS items 8-12 and the resume block. Deleted after it runs."""
from pathlib import Path

p = Path("NEXT_STEPS.md")
t = p.read_text(encoding="utf-8")

t = t.replace("# must print: 207 passed, 4 deselected, 1 xfailed",
              "# must print: 223 passed, 5 deselected, 1 xfailed")
t = t.replace("# top should be 50ff1c3 or later", "# top should be 3e82547 or later")
t = t.replace(
    "The four slow tests are the `dst_fit` suite gate, an OpenFOAM re-run of the plane-Couette\n"
    "anchor (~25 min, needs WSL), the real-curve half of the model-comparison control, and a\n"
    "re-run of the oscillatory-capability census (needs the gitignored workbook). Run them\n"
    "with `-m slow` before believing anything about the L3 bridge.",
    "The five slow tests are the `dst_fit` suite gate, an OpenFOAM re-run of the plane-Couette\n"
    "anchor (~25 min, needs WSL), the real-curve half of the model-comparison control, a re-run\n"
    "of the oscillatory-capability census (needs the gitignored workbook), and a re-run of the\n"
    "blood seed-scatter study (~20 min). Run them with `-m slow` before believing anything\n"
    "about the L3 bridge.")

BODY = '''### 8. DONE 2026-09-06: the envelope was written rounded INWARD
`%.6g` rounds whichever way the digits fall, so a MAXIMUM could round down and a MINIMUM
up, putting a record's own extreme measured rate outside its own envelope --
`FluidRecord.eta_at` then refuses to evaluate where the data is. **Measured before the
fix: 12 of the 15 Carbopol records, not the handful the item implied; 0 of 60 blood
records (0.1 and 700 are exact at six figures). After: 0 of 15.**

`registry.schema.round_out` rounds outward, working on the ROUNDED DECIMAL rather than on
`ceil(x/q)*q` -- the arithmetic form widens an exact 0.1 to 0.100001, because the float
0.1 sits a hair above the decimal. One edge case the property test caught: when `%.6g`
crosses a power of ten (999999.9 -> 1e+06) the step quantum must come from `x`'s exponent,
not the rounded value's. Both generators carry the gate (read the record BACK, require its
envelope to contain both extremes of its own source curve, exit non-zero otherwise), and
`tests/test_registry.py` has a NEGATIVE CONTROL that rebuilds the old envelopes and
requires exactly 12 of 15 to fail -- if that reads 0 the gate is inert.

The anchor was re-run because its rates are `geomspace` over the envelope: the top Carbopol
rate moved 199.972 -> 199.973 and every gate value came back BIT-IDENTICAL, so "no shipped
number moves" is verified rather than asserted.

### 9. DONE 2026-09-06 (G9.1 PASS, G9.2 and G9.3 FAILED): multi-start seed scatter
The referee item was "document dst_fit's multi-start guess collision". Prose about it says
nothing about whether it matters, so it was measured: all 20 blood curves refitted five
times under four models, changing nothing but the optimiser's seed
(`scripts/multistart_seed_scatter.py` -> `data/horner2020/multistart_scatter.md`).

- **G9.1 PASSED**: 54 all-`measurable` records, worst seed spread 5.997e-05 (bar 1e-3).
- **G9.2 FAILED**, on the held-out half: Herschel-Bulkley on donor G moves 4.149 (415 %)
  in `tau_y` while that parameter is graded `bounded_only` -- neither of the verdicts it
  asked about, and NOT one `check_exportable` refuses.
- **G9.3**, the stricter replacement added beside it with G9.2's number kept: nothing the
  export ACCEPTS may move more than 1e-3. Fails over all 80 pairs on the same record;
  passes over the 43 export-accepted records in the shipped families.

The objective moves by at most 6.5e-10 across all 80 pairs, so this is FLATNESS and not a
failed search. `check_exportable` gained `strict_identifiability=True` (off by default) to
refuse `bounded_only`; it refuses 7 of the 58 records the default accepts, all
Carreau-Yasuda blood. The Wald-correlation-at-a-bound caveat is documented in
`identify/fit.py`: the `|corr| > 0.99` signature is read off a local quadratic
approximation and is consulted ONLY when the profile hit a bound, i.e. exactly where that
approximation is weakest.

**Do not read `bounded_only` as "determined, with a one-sided interval".**

### 10. DONE 2026-09-06 (all three PASS): the same question, asked of Carbopol
Item 9 covered blood only. `scripts/multistart_seed_scatter_carbopol.py` ->
`data/carbopol_multistart_scatter.md` closes that on an independent set -- different
material, instrument, 54 points per curve, and a DIFFERENT FITTER for the shipped values.
G10.1 and G10.2 pass at a worst spread of 5.409e-06; **G10.3 is the one specific to this
set**: because the shipped values come from `carbopol/fit_hb.py`, seed scatter moves the
SIDECAR rather than the record, and `write_identifiability`'s containment rule had only
ever been checked at seed 0. All 225 (seed, parameter, record) intervals contain the
shipped value, and the verdict sets are seed-stable on 15/15.

So the blood failure is one record in one family that does not ship, not a property of the
identify layer.

### 11. NEXT ITEM, not started: the L3 bridge for `newtonian` and `power_law`
Written for zero context. `solver/openfoam.py` exports four of the registry's seven
constitutive families and raises `NotImplementedError` for `newtonian`, `power_law` and
`structural_steady`. OpenFOAM v2606 ships `Newtonian` and `powerLaw` in the same directory
as the four already mapped
(`/usr/lib/openfoam/openfoam2606/src/transportModels/incompressible/viscosityModels/`,
readable in WSL). `structural_steady` is thixotropic and has no counterpart there; that
refusal should be made EXPLICIT, with a message saying why, rather than left as a generic
`NotImplementedError`.

**Read both `calcNu` bodies FIRST and write the parameter map out before coding**, exactly
as the module docstring does for the four existing families. Item 6 is the reason it is
worth the ten minutes: its stated premise about `BirdCarreau` turned out to be false, and
`CrossPowerLaw`'s parameter names collide with biorheo's with the meanings SWAPPED.
Assume nothing about `powerLaw`'s names, or about whether it carries `nuMin`/`nuMax`
limiters -- check.

Gate to pre-register before writing code, by analogy with item 6 (do not weaken it):
the round trip reproduces `eta_at` on 20 rates inside the envelope to 1e-10; the export
refuses a `degenerate` / `not_identifiable` parameter without the override; if `powerLaw`
has viscosity limiters, DECIDE AND DOCUMENT whether they are a clamp the caller chooses
(as for Herschel-Bulkley) or fitted parameters (as for Cross), and make `coefficients`
require or refuse `eta_max` accordingly -- getting that wrong silently is exactly what G8
exists to prevent; and a discriminating negative case in the style of G3, some plausible
WRONG map that must fail the round trip by more than 1e-3, reported beside the verdict.

There is currently no `newtonian` or `power_law` RECORD in the registry, so this needs one
of each before the plane-Couette anchor can be extended. Do not invent one: either derive
it from a source already on disk with a DOI that resolves, or state plainly that the export
is tested by round trip only and put that limitation in the gate. Note the anchor already
has a hand-made Newtonian control (`newtonian_control()`), which is NOT a record and does
not exercise the export path.

### 12. Also not started: injectability for the families with no yield stress
`injectability` has closed forms for Herschel-Bulkley and Casson only. The registry now
holds 20 Cross and 20 Carreau-Yasuda blood records fitted to the SAME curves, and the
needle case study's central finding -- that the yield stress is invisible in this geometry
-- makes a prediction that can be pre-registered and tested: the plunger force computed
from a donor's Cross record should agree with the force from that donor's Casson record to
within the donor-to-donor spread, because the two describe the same curve over the same
envelope and the yield term contributes under 3 %.

That is a real cross-check of the whole case study rather than an extension of it, and it
can fail. Cross and Carreau-Yasuda have no closed-form tube solution, so it needs the
numerical quadrature route rather than a new derivation. Beware the trap item 4 already
hit: `flow_rate_casson`'s value near yield is limited by the CALLER's conditioning in
`1 - tau_y/tau_w`, so compare forces at the operating point, not near yield.

'''

start = t.index("### 8. NEXT ITEM, not started: the envelope is written rounded INWARD")
end = t.index("## What NOT to do")
t = t[:start] + BODY + t[end:]

old = ("- Do not change `MAX_REPRO_DEV`, `RUNS_Z`, gate bars, `LAMINAR_RE`, `BAR_*`, `G7P_RELAX_U`\n"
       "  or `G7P_END_TIME` in `scripts/openfoam_couette_anchor.py` to make a check pass.\n"
       "  `G7P_*` are the configuration being MATCHED, not bars.")
new = ("- Do not change `MAX_REPRO_DEV`, `RUNS_Z`, gate bars, `LAMINAR_RE`, `BAR_*`, `G7P_RELAX_U`\n"
       "  or `G7P_END_TIME` in `scripts/openfoam_couette_anchor.py` to make a check pass.\n"
       "  `G7P_*` are the configuration being MATCHED, not bars. Same for `G91_BAR`, `G92_BAR`,\n"
       "  `G93_BAR`, `G101_BAR`, `G102_BAR` in the two seed-scatter scripts: G9.2 and G9.3\n"
       "  FAILED and their bars stay where they were written.\n"
       "- Do not read `bounded_only` as \"determined\". It is not one of the verdicts the export\n"
       "  refuses by default, and such a parameter has been measured moving 415 % under nothing\n"
       "  but a change of the fitter's seed.")
assert t.count(old) == 1
t = t.replace(old, new)

p.write_text(t, encoding="utf-8", newline="\n")
print("ok")
