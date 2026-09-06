# Pre-registration — NEXT_STEPS item 11: the L3 bridge for `newtonian` and `power_law`

Written **2026-09-06, before any line of code**, and committed on its own so a later
session cannot claim a bar was chosen after the answer. Same discipline as
`nonyield_bridge_prereg.md` (item 6), which is worth re-reading first: item 6's stated
premise about `BirdCarreau` was FALSE, and `CrossPowerLaw`'s parameter names collide with
biorheo's *with the meanings swapped*. So nothing below is assumed from the family names.

## 0. What was read first, and what it says

> **Correction, see the addendum at the foot of this file:** §0.2's last sentence is FALSE.
> Only `HerschelBulkley` and `Casson` use `VSMALL`; `CrossPowerLaw`, `BirdCarreau` and
> `Newtonian` have no floor at all. The text is left as written because a pre-registration
> is not rewritten; the corrected statement is in `openfoam.py` beside `SMALL`.

`/usr/lib/openfoam/openfoam2606/src/transportModels/incompressible/viscosityModels/`,
read in WSL on 2026-09-06.

### `Newtonian/Newtonian.C`

```
calcNu() (lines 48-57):   returns a uniform field of nu0_
constructor (line 71):    nu0_("nu", dimViscosity, viscosityProperties_)
read()      (line ~103):  viscosityProperties_.readEntry("nu", nu0_)
```

**`Newtonian` is the only one of the six that does NOT use a `<Model>Coeffs`
sub-dictionary.** Every other model in that directory constructs from
`viscosityProperties.optionalSubDict(typeName + "Coeffs")`; `Newtonian` reads the key
`nu` straight out of the top-level `transportProperties`. A `NewtonianCoeffs { nu ...; }`
block is therefore not a Newtonian dictionary — OpenFOAM would look for a top-level `nu`,
not find one, and stop. There are no other coefficients and no limiters.

### `powerLaw/powerLaw.C`

```
calcNu() (lines 53-73):
    max(nuMin_, min(nuMax_, k_*pow(max(1[s]*strainRate(), SMALL), n_.value() - 1)))
constructor (lines 88-91):
    k_("k", dimViscosity, powerLawCoeffs_)      n_("n", dimless, powerLawCoeffs_)
    nuMin_("nuMin", dimViscosity, powerLawCoeffs_)
    nuMax_("nuMax", dimViscosity, powerLawCoeffs_)
read() (lines 118-121): readEntry on all four
```

Three things that are NOT what the item assumed and are settled by reading the source:

1. **`nuMin` and `nuMax` are REQUIRED, not optional.** They are constructed by the
   `dimensionedScalar(name, dims, dict)` form and re-read with `readEntry`, both of which
   raise on a missing key — unlike `BirdCarreau`'s `a`, which uses `getOrDefault`. So
   there is no silent-default hazard here, and the export has no choice: it must write
   both.
2. **The floor `SMALL` is `1.0e-15`, not `VSMALL = 1.0e-300`** (doubleScalar.H:62 vs 64),
   and it is applied to the *dimensionless product* `1[s]*strainRate()`, not to the rate
   itself. The other four models this bridge already covers use `VSMALL`. Transcribing
   `VSMALL` here would be wrong.
3. **The exponent slot carries `n`, and OpenFOAM subtracts the 1 itself**
   (`n_.value() - scalar(1)`). biorheo's `power_law` is `eta = k gd^(n-1)`, written with
   the `-1` explicit. A reader who matches the two expressions term by term without
   noticing that OpenFOAM's `pow` already has `n-1` in it writes `n_of = n_bio - 1`. That
   is the discriminating wrong map for this family (G11.5b).

`k` goes into a slot OpenFOAM declares `dimViscosity`, and for `n != 1` that is
dimensionally false in exactly the way it already is for `HerschelBulkley` — `k_si` is
Pa·s^n so `k/rho` is m²·s^(n-2). It is numerically right because OpenFOAM carries the
`1[s]` factor and everything here is SI. Do not "fix" it by writing a dimension bracket.

### Parameter map, written out before coding

```
newtonian    eta = mu                        ->  Newtonian    (TOP LEVEL, no Coeffs block)
    nu = mu/rho
power_law    eta = k gd^(n-1)                ->  powerLaw
    k = k/rho    n = n    nuMin = eta_min/rho    nuMax = eta_max/rho
```

## 1. The clamp decision, stated before it is coded

Item 11 requires this to be DECIDED AND DOCUMENTED, because getting it wrong silently is
what G8 exists to prevent. The decision:

**`nuMin` and `nuMax` for `powerLaw` are a two-sided CLAMP the caller chooses, not fitted
parameters, so `coefficients` REQUIRES both `eta_max` and `eta_min` for `power_law`.**

The argument is the one already used for `HerschelBulkley` and `Casson` and it is
stronger here. `power_law` has no plateau at either end: for `n < 1` it diverges as
`gd -> 0` and goes to zero as `gd -> inf`, and for `n > 1` (the bounds allow up to 2.0)
the divergence is at the top instead. A record's `k` and `n` are the fit; neither `nuMin`
nor `nuMax` is anywhere in the record, and there is no floor that "never changes the law"
the way Casson's `eta_c` is. So both are regularisation and both are the caller's.

Consequences, all of which must be enforced rather than documented:

- `newtonian` **refuses** both `eta_max` and `eta_min`: a constant viscosity has nothing
  to regularise and OpenFOAM's `Newtonian` has no limiter slot to put one in.
- `herschel_bulkley` **refuses** `eta_min` (OpenFOAM's `HerschelBulkley` has only the
  one-sided `min(nu0, ...)`; there is no floor slot) and still requires `eta_max`.
- `casson` **refuses** `eta_min`: its `nuMin` is set to `eta_c`, the record's own
  infinite-rate limit, which is the largest floor that cannot change the law
  (`eta > eta_c` at every finite rate). That is determined by the record, not chosen.
- `cross` and `carreau_yasuda` continue to refuse both.

A new API argument is added: `eta_min`, keyword-only, default `None`.

## 2. Gates. These are the bars. They are not to be moved.

Every one is graded by a script and pinned in a test. A FAILED gate is recorded with its
number and a *stricter* test is added beside it; it is never loosened.

**G11.1 — round trip, both new families.** For every `newtonian` and `power_law` record
in the registry: export -> parse -> `nu_openfoam` reproduces `record.eta_at(gd)/rho` at
20 log-spaced rates inside the envelope, wherever neither clamp binds.
**Bar: max relative error <= 1e-10.**

**G11.2 — the verdict refusal reaches the new families.** A record carrying a
`degenerate` or `not_identifiable` verdict is refused by `export_transport_properties`
without `allow_unidentified=True`, and accepted with it. **Bar: pass/fail, and the
refusal must name the offending parameter.** (For `newtonian`, whose single parameter is
a median, a refusal may never fire on the shipped set; the gate is then graded on a
constructed record and says so.)

**G11.3 — the clamp decision is enforced, and "required" is decided by RUNNING the
solver.**
  a. `coefficients` raises for: `power_law` without `eta_min`; `power_law` without
     `eta_max`; `eta_min` passed to any of the other five families; `eta_max` passed to
     `newtonian`, `cross`, `carreau_yasuda`. **Bar: all eight raise, each with a message
     naming the family and the reason.**
  b. A `powerLaw` dictionary with the `nuMin` entry DELETED is handed to `simpleFoam`.
     **Predicted: it fails (non-zero rc) with a message naming `nuMin`.** If it instead
     runs, `nuMin` is optional after all, the premise of §1 is wrong for that key, and
     that is recorded as the result. **Bar: rc != 0 and "nuMin" appears in the log.**

**G11.4 — the Newtonian dictionary structure, also decided by running the solver.**
  a. The export writes a top-level `nu` and NO `NewtonianCoeffs` block.
  b. The same record written with `nu` inside a `NewtonianCoeffs` block is handed to
     `simpleFoam`. **Predicted: it fails (non-zero rc).** **Bar: rc != 0.** If it runs,
     the sub-dictionary form is legal too and the finding in §0 is wrong for `Newtonian`;
     that is recorded.

**G11.5 — discriminating wrong maps, in the style of G3.** Each must FAIL the round trip
by more than 1e-3 relative, and the separation is reported beside the verdict so that
G11.1 passing at 1e-10 is evidence about the MAP and not only about arithmetic.
  a. `newtonian`: `nu = mu` (the dynamic viscosity, `/rho` forgotten).
  b. `power_law`: `n_of = n_bio - 1` (the `-1` applied twice — see §0.3).
  **Bar: both separations > 1e-3.** If either comes in below 1e-3 on the shipped records
  the gate reports itself UNINFORMATIVE for that family rather than passing.

**G11.6 — the `structural_steady` refusal, and its stated reason is CHECKED, not
asserted.** Item 11 says `structural_steady` "is thixotropic and has no counterpart
there". Before writing that sentence into an error message it is tested, because
`STRUCTURAL_STEADY.limits` already claims the opposite:
`"any gd": "equals Cross(eta_0 = eta_inf (1 + beta), eta_inf, lam = alpha_tr, m = 1)"`.
  **Pre-registered numeric claim:** for 20 log-spaced rates over four decades and three
  parameter sets, `structural_steady(gd; eta_inf, beta, alpha_tr)` equals
  `cross(gd; eta_0 = eta_inf(1+beta), eta_inf = eta_inf, lam = alpha_tr, m = 1)` to
  **<= 1e-12 relative**.
  - If it HOLDS, the steady flow curve *does* have an exact OpenFOAM counterpart, the
    item's stated reason is wrong, and the refusal message must give the true reason
    instead: OpenFOAM's `viscosityModel` interface is a function of the instantaneous
    strain rate only, so the STRUCTURE EQUATION cannot be carried; an export would be
    right in steady state and silently wrong in every transient. The message must also
    give the exact substitution a caller can make deliberately.
  - If it FAILS, the item's wording stands and the number is recorded.
  **Bar: <= 1e-12, and whichever way it goes the refusal is an explicit, named error
  (not a bare `NotImplementedError`) carrying the reason that the test supports.**

**G11.7 — the plane-Couette anchor, extended (slow, needs WSL).** One `newtonian` and one
`power_law` record run through `simpleFoam` at `N_RATES` wall speeds over the whole
measured envelope. OpenFOAM's own `nu` field and its `wallShearStress` function object
are compared with `record.eta_at`, touching none of this repo's transcription of
`calcNu`. **Bars: the existing ones, unchanged — A1' `<= 1e-9`, A2 `<= 1e-6`,
A3 `<= 1e-6`.**

## 3. The records, and the limitation if there are none

There is no `newtonian` and no `power_law` record in the registry, and G11.1/G11.7 need
one of each. **No record will be invented.** The registry's rule is that a record is a
fit to a curve from a source on disk whose DOI resolves.

The plan is to fit both families to the **Horner 2020 human blood curves already on
disk** (`DOI_PAPER` in `scripts/records_from_horner.py`; the same 16 points per donor that
the four existing families are fitted to) by adding them to `records_from_horner.py`'s
`COMPARE` and `EXTRA_FAMILIES`. That is not an invention: it is the same generator, the
same source, the same envelope.

It is also not a claim that blood is Newtonian or a power-law fluid. Both are NULL
MODELS — `newtonian` is described in `constitutive/gnf.py` as "the null model every other
one must beat" — and their records will say so in `notes`, will carry their own
`rmse_log` beside Casson's, and will be flagged `structured` by
`residual_structure` if they misfit, which is the mechanism the L3 layer already has for
exactly this. Adding them also gives `model_comparison.md` the null column it currently
lacks, which is a result in its own right: **it is worth knowing whether the published
two-parameter Casson fit beats a one-parameter constant on every donor.** That comparison
is reported whichever way it falls.

If for any reason the records cannot be written, the fallback stated in item 11 is taken
instead: the export is tested by round trip only, G11.7 is recorded NOT RUN, and the
limitation goes in the gate report rather than being worked around.

## 4. What would make this item's work wrong, stated in advance

- Writing a `NewtonianCoeffs` block because every other family has one (G11.4).
- Applying the `-1` to the power-law exponent twice (G11.5b).
- Transcribing `VSMALL` instead of `SMALL`, or applying the floor to the rate rather
  than to `1[s] * rate` (§0.2).
- Accepting an `eta_min` for a family with no floor slot and dropping it silently — the
  failure mode G8 was written for.
- Writing "no counterpart exists" into the `structural_steady` refusal without checking
  it (G11.6), which is item 6's mistake repeated.

---

## Addendum, 2026-09-06, after the run and the code review

**Nothing above this line was edited.** The commit that carries it, `e9c1ace`, contains
this file and nothing else, and is the proof of order; the code it grades landed in
`bcf84bf`. No bar was moved. What follows is what the run and the review found.

### The gates

| gate | verdict | number | bar | where |
|---|---|---|---|---|
| G11.1 `newtonian` | PASS | 0.000e+00 (exact) | 1e-10 | `newtonian_powerlaw_gates.md` |
| G11.1 `power_law` | PASS | 2.392e-16 | 1e-10 | same |
| G11.2 | PASS | pass/fail | — | `tests/test_openfoam_export_null.py` |
| G11.3a | PASS | 8 refusals, each with its own reason | — | same |
| G11.3b | **PASS, prediction confirmed** | control rc 0, broken rc 1, `nuMin` named | rc != 0 | `couette_anchor.md` |
| G11.4a | PASS | pass/fail | — | `tests/test_openfoam_export_null.py` |
| G11.4b | **PASS, prediction confirmed** | control rc 0, broken rc 1 | rc != 0 | `couette_anchor.md` |
| G11.5a | PASS | 1059 minimum separation | > 1e-3 | `newtonian_powerlaw_gates.md` |
| G11.5b | PASS | 9 minimum separation | > 1e-3 | same |
| G11.6 | **PASS, and it REFUTES item 11's premise** | 1.386e-16 | 1e-12 | same |
| G11.7 | PASS | `newtonian` err_nu 0.000e+00, `power_law` 6.370e-15 | A2 1e-6 | `couette_anchor.md` |

OpenFOAM's own words on the two dictionaries it was predicted to reject:
`Entry 'nuMin' not found in dictionary ".../transportProperties/powerLawCoeffs"` and
`Entry 'nu' not found in dictionary ".../transportProperties"`. The second is the more
informative of the two: OpenFOAM looked for `nu` at the TOP LEVEL and did not find it,
i.e. the `NewtonianCoeffs` block was not read at all — which is the claim §0 made from
`Newtonian.C:71` and the reason the parser refuses such a block rather than ignoring it.

### G11.6: the item's stated reason was wrong, and the refusal is kept with a new one

§2 said this in advance and it happened. `structural_steady`'s steady flow curve equals
`CrossPowerLaw` with n = 1 to 1.386e-16, so "no counterpart exists" would have refused
something OpenFOAM can do — item 6's mistake repeated. `ThixotropyNotRepresentable` now
gives the true reason: `viscosityModel` is a function of the instantaneous `strainRate()`
alone, so the structure equation cannot be carried and the export would be right in
steady state and silently wrong in every transient.

### What §0 got wrong, recorded rather than edited out

**§0.2's last sentence is FALSE.** It says "The other four models this bridge already
covers use `VSMALL`". Only two do — `HerschelBulkley.C:65` and `Casson.C:69`, which
divide by the rate. `CrossPowerLaw`, `BirdCarreau` and `Newtonian` have **no floor at
all**, because nothing in them is singular at zero rate. No number moves: `nu_openfoam`
transcribes each model as written and the two that have no floor never had one. The
corrected statement is in the module docstring beside `SMALL`. Found by code review, not
by the run.

### What the gates do NOT establish, measured after the fact

**Both G11.5 separations are constants and neither is a property of the fit.** G11.5a is
`rho - 1` = 1059 on all 20 records (spread 2.3e-13) because forgetting `/rho` is a pure
scaling; G11.5b is `1/gammadot_min - 1` = 9 on all 20 (spread 7.1e-15) because the
doubled `-1` multiplies the law by `gd^-1` and the largest such factor sits at the bottom
of the envelope. So they establish that the slot is USED, and no more. A discriminator
sensitive to the fitted values would need a different candidate. That second
attribution is not free -- the clamp sits 10x outside the envelope and the wrong map's
value at `gammadot_min` lands on `nuMax`, so 9 could equally have been the clamp factor
minus 1 -- and it survives its own control: `item11_gate_report.py` re-measures it with
the clamp 100x out and prints the number beside the verdict.

### Five rigor defects the review found in the code that implements these gates

All fixed in `5f6ac8c`; listed because each could have produced a passing verdict from a
broken measurement.

1. **The Newtonian parser was FAIL-OPEN on nested sub-dictionaries** and returned a
   silently wrong `nu` — a regex cannot match balanced braces, so everything after a
   nested block leaked out as top level. It returned the buried value even with the real
   entry deleted, i.e. it succeeded on a dictionary OpenFOAM rejects.
2. **The anchor's exit code did not carry G11.7 or the refusal probes.** They were in the
   JSON only, so a FAILED gate left `rc = 0`.
3. **The per-family aggregation used builtin `max()`, which DROPS a NaN**
   (`max([5.0, nan]) == 5.0`), so a family could report PASS on a case that produced no
   number.
4. **`G7'` iterated `RECORD_IDS`**, so appending two records here would have re-decided
   item 6's CLOSED gate — and its yield/no-yield split keys on `NO_CLAMP_MODELS`, which
   lists neither null family, so both NULLS would have been counted as yield-stress arms.
   `G7'` now has its own frozen list, and its `worst_yield` is unchanged at 1.7202393803.
5. **Two docstring sentences claimed G11.3b and G11.4b had been decided by running the
   solver** while the anchor was still running. They now say the claim is read from the
   C++ and point here for the verdict — which, as of this addendum, agrees.

### Effect on the numbers already on this page's sibling artifacts

Of `couette_anchor.json`'s 33 pre-existing rows, **every field of all 33 is bit-identical**
except `contraction`, which is NaN in both (NaN compares unequal to itself). Every item-6
gate scalar is unchanged. **A3's headline moved 5.53e-14 -> 6.28e-14 and no old row
moved**: A3 is a maximum over every case on the page, and a NEW `newtonian` row is now the
largest. That distinction is the point of reporting the item-11 families in their own
table.
