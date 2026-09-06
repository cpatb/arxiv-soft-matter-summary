# Pre-registration: the L3 bridge for the families that have no yield stress

Written 2026-09-05, **before any code was written**, and committed on its own so the
history proves the order. NEXT_STEPS item 6. Nothing below may be loosened; a bar that
fails is recorded as FAILED with its number and a *stricter* test may be added beside it.

## 0. A premise of item 6 is wrong, and this is what replaces it

Item 6 as written says:

> OpenFOAM's `BirdCarreau` has no `a` exponent (it is fixed at 2) while biorheo's
> Carreau-Yasuda does, so the map is partial and the export must REFUSE a record whose
> `a` is not 2 rather than silently dropping it -- that refusal is the point of the task.

**That is false for OpenFOAM v2606.** `BirdCarreau.C` lines 81-88 and 120-124 read

```
    a_
    (
        BirdCarreauCoeffs_.getOrDefault
        (
            "a",
            dimensionedScalar("a", dimless, 2)
        )
    ),
```

and `calcNu` (lines 55-59) uses it:

```
    nuInf_ + (nu0_ - nuInf_)*pow(scalar(1) + pow(k_*strainRate(), a_), (n_ - 1.0)/a_)
```

so `a` is a real, honoured coefficient whose only special property is that it is
**optional with a silent default of 2**. The map is therefore TOTAL, not partial, and an
`a != 2` refusal would be refusing something the solver can do.

The hazard of the same shape is still there and it points the other way: because `a` is
optional, a dictionary that OMITS it is a valid `BirdCarreau` that runs as a Carreau
fluid. So the silent-drop risk is in the WRITE (an export that forgets `a` becomes a
different fluid with no error) and in the READ (a parser that treats a missing `a` as an
error disagrees with OpenFOAM). Gates G4 and G5 below test exactly that, and G4 is the
one that decides the premise empirically rather than by reading C++.

Also recorded before coding, because it is the trap item 6 did warn about correctly: the
Cross parameter names **collide across the two codebases with the meanings swapped**.

| biorheo `cross` | meaning | OpenFOAM `CrossPowerLaw` |
|---|---|---|
| `eta_0` | zero-shear plateau, Pa s | `nu0` (kinematic) |
| `eta_inf` | infinite-shear plateau, Pa s | `nuInf` (kinematic) |
| `lam` | time constant, s | **`m`** (`dimTime`) |
| `m` | exponent, dimensionless | **`n`** (`dimless`) |

biorheo `eta = eta_inf + (eta_0 - eta_inf)/(1 + (lam gd)^m)`;
OpenFOAM `nu = (nu0 - nuInf)/(1 + (m sr)^n) + nuInf`. Same form, and a map written by
matching NAMES would put an exponent into a time slot. G3 exists to prove the round trip
can detect that.

Carreau-Yasuda -> `BirdCarreau` is a straight rename with no collision:
`eta_0 -> nu0`, `eta_inf -> nuInf`, `lam -> k` (`dimTime`), `a -> a`, `n -> n`.

## 1. The gates

Bars are fixed here. `rho` is the caller's density, as for the two existing families;
records carry none.

**G1 round trip.** For every `cross` and `carreau_yasuda` record in the registry, the
chain record -> exported text -> `parse_transport_properties` -> `nu_openfoam` reproduces
`rho * record.eta_at` on 20 log-spaced rates spanning the record's envelope to a relative
error <= **1e-10**. Same bar as the existing families.

**G2 refusal on verdict.** `check_exportable` refuses a record any of whose parameters is
graded `degenerate` or `not_identifiable` unless `allow_unidentified=True`. Reported as a
count over the 20 donors per family; the count itself is the result, not a bar. The bar is
that the refusal fires on **every** such record and on no other (no false positives, no
false negatives against the sidecar).

**G3 discriminator: the wrong Cross map must fail G1.** A deliberately name-matched map
(biorheo `m` written into OpenFOAM's `m` slot, biorheo `lam` into `n`) evaluated through
the same `nu_openfoam` must differ from `rho * record.eta_at` by **more than 1e-3
relative on at least one of the 20 rates**, for the shipped Cross record. If the wrong map
also passes at 1e-10, G1 has discriminated nothing and G1's verdict is void. Report the
candidate-separation factor beside the verdict.

**G4 the `a` premise, decided by a running solver, not by reading C++.** Pre-registered
prediction, written now: **OpenFOAM v2606 honours `a`,** so a plane-Couette case built
from a Carreau-Yasuda record whose `a` is far from 2 will reproduce `record.eta_at`
through OpenFOAM's own `nu` field at the existing `BAR_A2 = 1e-6`.

The test is only informative if the two candidates are far apart, so the separation is
pre-registered too: for the record chosen, `max |nu(a fitted) / nu(a = 2) - 1|` over the
20 anchor rates must be **>= 1e-2**, i.e. >= 1e4 times `BAR_A2`. If the record's fitted
`a` is too close to 2 to give that separation, the record is unsuitable and a synthetic
`a` far from 2 is used for this gate alone, stated as such. Outcomes:

- passes A2 at 1e-6 -> `a` is honoured, item 6's premise is refuted, no refusal is written;
- fails A2 by about the separation -> `a` is ignored, item 6's premise stands, and the
  export must refuse `a != 2`.

**G5 the omission hazard.** A `BirdCarreau` dictionary with the `a` entry deleted must
parse to `a = 2` (mirroring `getOrDefault`) rather than raising, and `nu_openfoam` on that
parsed dictionary must equal the `a = 2` law and therefore differ from the shipped one by
the G4 separation. And the exporter must ALWAYS write `a`, checked by asserting the token
is present in the exported text.

**G6 the anchor, at the existing bars.** `scripts/openfoam_couette_anchor.py` is extended
with one `cross` and one `carreau_yasuda` record added to `RECORD_IDS`, and must pass
**A1' <= 1e-9, A2 <= 1e-6, A3 <= 1e-6** -- the bars already in that script, unchanged.

**G7 the control that item 6 is really about.** These families have no yield stress, so
the Bingham number `B = tau_y/(k gd^n)` that controls `simpleFoam`'s convergence for the
Herschel-Bulkley record is identically zero. Pre-registered claim: the new-family cases
converge without the pathology of item 3, quantified as **worst `err_nu` over all new
cases <= `BAR_A2` = 1e-6**, against the 172 % that the Herschel-Bulkley record returns at
B = 29. Iteration counts are reported beside it. If any new case exceeds 1e-6 the control
is recorded FAILED and the "B story is about yield, not about hard fluid" reading of
item 3 is weakened in writing.

**G8 no clamp knob.** For `herschel_bulkley` and `casson`, `nu0`/`nuMax` is a
regularisation the caller chooses and `eta_max` is a required argument. For `cross` and
`carreau_yasuda` the zero-shear plateau `eta_0` is a FITTED PARAMETER and the law is
bounded without any clamp, so there is no `eta_max` to choose. The export must therefore
**refuse an `eta_max` for these families** rather than accept and ignore it; silently
ignoring a caller's regularisation request is the same silent-drop failure this item
exists to prevent. Checked by a test that asserts the raise and names `eta_max`.

## 2. What would make me wrong

- G3 passing the wrong map: then G1 proves nothing and both are reported void.
- G4's separation coming out below 1e-2: the gate is uninformative and is re-run with a
  synthetic `a`, said so in the results.
- G7 failing on any new case: the control fails and item 3's reading is weakened, not the
  bar.

## 3. Records this needs

The registry has no `cross` or `carreau_yasuda` record. `scripts/records_from_horner.py`
already FITS both families to all 20 donor curves for `model_comparison.md` but writes
only the Casson record. It is extended to write all three families, keeping the existing
Casson ids unsuffixed so nothing that pins them breaks, with the new ids suffixed
`_cross` and `_carreau_yasuda`. That asymmetry is ugly and is documented in the script.

All 20 donors are written, not one, because **G2's refusal count is only a measurement
over the whole set** -- one record cannot demonstrate a refusal that depends on
identifiability.

---

## Addendum, 2026-09-05, after the run and the code review

Nothing above was edited except a wrong date in the first line (it said 2026-09-06; the
file was written on 2026-09-05, and the commit that carries it, `2561346`, is the proof
of order). No bar was moved. What follows is what the run and the review found.

**G1, G2, G3, G5, G6, G8: PASS.** Numbers in `tests/test_openfoam_export_nonyield.py`
and in `couette_anchor.md`.

**G4: PASS, and the premise recorded in §0 is REFUTED by a running solver**, which is
what §0 said it would take. The Carreau-Yasuda record's fitted `a` is far enough from 2
that ignoring the entry would move `nu` by 1.1e-1; A2 on those cases lands far below its
1e-6 bar instead. No `a != 2` refusal is shipped.

**G7: PASS as written, and UNINFORMATIVE — the gate was mis-specified, by me, here.**
Its number stands and is not withdrawn. But G7 compares the new families, run at the
anchor's `relax_U = 0.95, endTime = 20000`, against the 172 % figure from
`couette_diagnosis.json`, which was taken at `relax_U = 0.7, endTime = 2000`. **At the
anchor's own settings the Herschel-Bulkley record passes A2 as well.** Both arms pass, so
the comparison decides nothing about yield. This is the same error the lessons file calls
"a discriminator that passes BOTH candidates has discriminated nothing", and it was in
the pre-registration rather than in the code.

The replacement is STRICTER and is added rather than substituted, following this repo's
A1/A1' and A5/A5' pattern:

**G7' matched settings.** Every record in `RECORD_IDS` — both yield-stress families and
both no-yield ones — is run once at its own lowest envelope rate at the SAME
`(relax_U, endTime) = (0.7, 2000)`, the configuration in which the Herschel-Bulkley
record is known to fail. Pre-registered before it ran: **the comparison is informative
only if at least one yield-stress record fails `BAR_A2` there**, because with no failing
arm there is nothing to contrast; if none fails, G7' is reported UNINFORMATIVE and no
claim about yield may be quoted from the file.

**And G7' still does not separate yield from viscosity magnitude.** The no-yield record's
zero-shear plateau is three to four decades below the Herschel-Bulkley case's viscosity
at the failing rate, so a "hard fluid" reading is not excluded by that table and the
report says so. The magnitude control is the constant-`B` 2x2 already in
`couette_diagnosis.json` — at fixed shear rate, 100x in `eta` at constant `B` changes the
error in the tenth digit against 5.6e4x across `B` — and that, not the new records, is
what the "about yield, not magnitude" reading rests on. The earlier draft of the report
asserted the opposite from these records alone; that sentence is withdrawn.
