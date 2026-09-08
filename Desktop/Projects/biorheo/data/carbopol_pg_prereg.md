# Item 24 pre-registration — a licensable yield-stress curve, and the first real test of the temperature shift

Written 2026-09-08, **before any generator code**, as every record in this repo has been.
Data downloaded and documented first: `data/carbopol_pg_vinutha2024/SOURCE.md`.

Unlike item 23's, every gate here is **blind**. Nothing below has been computed. The only
numbers already in hand are the ones printed in `SOURCE.md`, which are structural (point
counts, rate ranges, the four signals that disqualify the 20 °C block) and the two the
deposit itself publishes.

---

## 0. What this item is for, and the one thing it is not for

**Zenodo 11389467, CC BY 4.0**, 2 wt% Carbopol 974p in propylene glycol, Anton Paar MCR 300,
1 deg cone, **both surfaces roughened to minimise wall slip**, Peltier plate and hood,
T to 0.1 deg C. Flow curves at 10, 20 and 50 deg C. Paper arXiv:2312.00251.

Three reasons it is worth more than the DMEM records of item 22:

1. It is a **yield-stress microgel** — the same family as the 15 CC BY-NC-SA Carbopol
   records the injectability case study actually ships, so it is the licensable counterpart
   to the family this project demonstrates on.
2. It reaches **1e3 1/s**, one to two decades from a needle wall rate, against the DMEM
   records' four to five.
3. It is the **first real test bed for `constitutive/shift.py`**, and the test is an anchor
   rather than a demonstration. See §1.

**It is not an injectability record.** 1e3 1/s is not 1e4–1e5, the suspending medium is
propylene glycol rather than water, and the material is a model system for a rheology study,
not a formulation anyone injects. Any force computed from it must say all three.

---

## 1. Why the shift test here is an anchor and not a demonstration

The authors chose propylene glycol because its viscosity moves about tenfold over
10–50 deg C, *"which enables us to assess the importance of the viscous stresses
(`eta_s gd`) independently of the shear rate"*. **The material is the same at all three
temperatures and only the suspending medium changes.** The deposit publishes both relevant
numbers in named columns of `Fig 1(c) and Fig S3 residual stress_exp.txt`:

| | 10 deg C | 20 deg C | 50 deg C |
|---|---|---|---|
| `medium viscosity [Pa s]` | 0.102 | 0.055 | 0.012 |
| `yield stress [Pa]` | 8.9 | 8.9 | 8.9 |

So the shift factor is **predicted, not fitted**. If temperature acts only through the
medium viscosity then

    a = eta_s(T) / eta_s(T_ref)

which is exactly `shift.shift_from_two_viscosities`. Meanwhile `shift.collapse_shift` fits
`a` from the two flow curves alone and **never sees `eta_s`**. Comparing them checks a
fitted quantity against a number from outside this repo.

The algebra, to be re-derived in the script rather than trusted here. Under
`eta(gd;T) = a·eta_ref(a·gd)` the stress transforms as `sigma(gd;T) = sigma_ref(a·gd)`, so
for Herschel-Bulkley `sigma = tau_y + k gd^n`:

    tau_y -> tau_y     (invariant)
    k     -> a^n · k
    n     -> n         (invariant)

The invariance of `tau_y` is item 21's derived consequence, written down before this dataset
was found. The authors state the same thing from their own data: *"The yield stress is
independent of the solvent viscosity, which denotes that sigma_y is exclusively set by the
intrinsically solid characteristics of the dense ensemble of soft spheres."* Two independent
routes to one statement, and this item tests it.

**Item 21 has only ever been exercised on a synthetic pair and on item 22's collapse, which
was pre-registered to FAIL and did. This is the positive case, and it has never been run.**

## 1.1 The trap that would make this test pass for the wrong reason

**Below yield the collapse residual is identically zero for every shift factor.** As
`gd -> 0` both curves tend to `tau_y`, so `eta -> tau_y/gd`, and

    a · eta_ref(a·gd) = a · tau_y/(a·gd) = tau_y/gd = eta_other(gd)

for **any** `a`. The low-rate half of the curve therefore contributes nothing to the
objective whatever the answer is, and `collapse_shift`'s rms over the whole curve is
dominated by points that cannot discriminate. A small residual would then mean "most of my
points are on the plateau", not "the curves collapse".

This is the "a discriminator that passes both candidates has discriminated nothing" failure,
and it is pre-empted rather than discovered: **G24.2 and G24.3 below measure it directly**,
and the item reports itself uninformative if the separation is too small.

---

## 2. Gates

All blind. `a_pred = 0.012/0.102 = 0.117647...` throughout, with 10 deg C as the reference
and 50 deg C as the other curve — the pair with the largest viscosity contrast (8.5x) and
the two blocks whose provenance is certified (§3).

### G24.1 — the fitted shift factor reproduces the published viscosity ratio

`collapse_shift` on the 10 -> 50 deg C flow curves, at their common rates, versus
`a_pred`. **Bar: `|ln(a_fit / a_pred)| < ln(1.30)`**, i.e. agreement to 30 %.

Predicted: PASS. If it fails, either the temperature does not act only through the medium
viscosity or the shift layer is wrong, and which one is decided by G24.4–G24.6.

30 % is chosen because it is far tighter than the 8.5x range `a` is being picked out of and
far looser than the 3 % scale on which the medium viscosities are quoted. It is set here,
before the number is known, and it will not be moved.

### G24.2 — and the test must be shown INFORMATIVE, or G24.1 decides nothing

Report the collapse residual on the **high-rate half** of the overlapping range (rates above
the geometric mean of the overlap) at `a_fit`, at `2·a_fit` and at `a_fit/2`.

**Bar: the residual at `2·a_fit` must exceed the residual at `a_fit` by a factor of at least
3.** Below that, this item reports G24.1 as **uninformative** and does not quote the
agreement, exactly as the frisbee discriminator now refuses to call a test informative below
a separation of 1.6x.

Predicted: PASS, and the separation is expected to be large because the high-rate half is
where the viscous stress lives.

### G24.3 — the plateau's insensitivity is measured, not argued

Same three shift factors, residual on the **low-rate half** instead.

**Bar: that residual changes by less than 10 % between `a_fit` and `2·a_fit`.**

Predicted: PASS, i.e. the low-rate half really is blind to `a`, which is what makes G24.2
necessary. This is the positive control for §1.1: if the low-rate residual *does* move, my
algebra is wrong and G24.2's framing must be rewritten before G24.1 is quoted.

### G24.4 — the yield stress comes out invariant, as item 21 derived and the authors state

Fit Herschel-Bulkley independently at 10 and 50 deg C. **Bar: the two `tau_y` profile
intervals overlap.**

Predicted: PASS. Report both intervals and both point estimates whatever happens; an
overlap of two wide intervals is a weak pass and must be labelled one, so also report the
ratio of the difference to the mean interval width.

### G24.5 — the flow index comes out invariant, because no shift can move it

**Bar: the two `n` profile intervals overlap.** Predicted: PASS.

This is the discriminating companion to G24.1: item 22's collapse failed precisely because
its two curves differed in `n`. If `n` differs here too, no shift factor is meaningful and
G24.1's number must not be quoted regardless of what it says.

### G24.6 — the consistency index moves by the predicted power of the shift factor

**Bar: `|ln( (k_50/k_10) / a_pred^n ) | < ln(1.40)`**, with `n` the mean of the two fitted
indices. Predicted: PASS.

This is a second, independent route to G24.1 that uses the fitted parameters rather than the
curve, so it can disagree with G24.1 and that disagreement would be informative.

### G24.7 — the 20 deg C block is refused, and the refusal is recomputed not asserted

The script must recompute all four disqualifying signals from `SOURCE.md` trap 2 and print
them: the fraction of rates on the programmed logarithmic grid, the count of rate
reversals, the significant-figure histogram of the stress column, and the sign of the offset
against the independent replicate.

**Bar: no record is written from the 20 deg C block, and all four signals reproduce.**

Predicted: PASS. If any signal does not reproduce, the exclusion is reconsidered in writing
rather than kept out of habit.

### G24.8 — put the parameter uncertainty and the measurement scatter in the same units

The cessation files carry a second, unannounced steady flow curve (`preshear-stress` at 16
rates per temperature), i.e. **the only independent replicate this project has for any
record**. Compare it with `Fig 1(a)` at matched rates, then compare that scatter with the
width of the fitted parameters' profile intervals.

**No bar; this is a measurement.** Report which is larger. G15.4 had to approximate this
comparison with a model misfit because no replicate existed; here one does.

Predicted, and recorded so it can be wrong: the replicate scatter will be the **larger**
term, as it was for the Carbopol case study's 0.1 wt% gel.

### G24.9 — the published 8.9 Pa is NOT graded against a Herschel-Bulkley `tau_y`

The 8.9 Pa is the authors' fit to the **three-component model**
(`sigma = sigma_y + k1 gd^0.5 + eta_bg gd`), which this repo does not implement. A yield
stress is a fit-objective convention before it is a material constant — Horner's Casson
`tau_y` and item 23's flow index both say so.

**Bar: no gate compares an HB `tau_y` against 8.9 Pa.** The difference is reported as a
model-form comparison with both forms named, and the fact that biorheo cannot reproduce the
authors' number *because it does not implement their model* is stated rather than presented
as a disagreement.

---

## 3. Which blocks may become records

**10 deg C and 50 deg C only.** Both are certified instrument output: rates on the
programmed 1/10-decade grid (37 of 48 and 50 of 60 exactly, the rest to 1e-4), log step
-0.1000 with a spread under 0.001, strictly monotone, stress carried to 7-8 and 4-5
significant figures.

**The 20 deg C block is excluded** on the four independent signals in `SOURCE.md` trap 2. It
cannot be certified as a raw curve at the grade of the other two, and this repo does not add
a record it cannot certify. The cost is the reference temperature the medium viscosity is
quoted at; it is not a cost to the shift test, which runs 10 -> 50 deg C.

Two records, `carbopol_pg_vinutha2024_10C` and `..._50C`, family `herschel_bulkley`,
`commercial_use: true`, envelope rounded OUTWARD per item 8, sidecars written only by
`biorheo.identify.write_identifiability`.

## 4. Consumers that depend on the registry's contents

Enumerated by grep BEFORE any record is written, as item 22's pre-registration required and
as item 23 found item 22 had still missed one — a *generator*, not a test. **This list
therefore covers `scripts/` as well as `tests/`.**

Adding `herschel_bulkley` records is the case item 22 did not exercise. Item 22 added a
family that existed only in the null models; this adds to a family with fifteen incumbents,
so every consumer that treats "the `herschel_bulkley` records" as a synonym for "the
Carbopol ETD2020 records" is now wrong.

| consumer | how it selects | effect of +2 `herschel_bulkley` | action |
|---|---|---|---|
| `tests/test_registry.py:397` | **`model.name == "herschel_bulkley"`, then asserts 15 and regex-matches each id against `carbopol_etd2020_*`** | ⚠ **the same defect item 22 found from the other side** — it names one family and means one dataset | select by `id.startswith("carbopol_etd2020_")`, as line 206 in the same file already does |
| `tests/test_registry.py:206` | `id.startswith("carbopol_etd2020_")` | none | nothing; this is the correct idiom |
| `tests/test_openfoam_export.py:38` | `model.name in ("herschel_bulkley", "casson")`, asserts 35 | a CLOSED item-3 gate whose sample would grow to 37 | freeze to `closed_gate_records` |
| `tests/test_dmem_records.py:102` | `len(load_all()) == 123` | 125 | restate as the count of DMEM records plus the frozen 115, or update with the reason |
| `tests/test_offer_claims.py:101` | `(n_nc, len(recs)) == (115, 123)` and the string `"8 of the 123 records"` in `OFFER.md` | 125, and **`OFFER.md`'s own sentence goes stale** | update both together; `OFFER.md` is pinned by this test on purpose |
| `tests/test_openfoam_export_nonyield.py:321`, `test_openfoam_export_null.py:189`, `test_dmem_records.py:221` | already frozen at 115 | none | nothing |
| `tests/test_design_space.py` | reads the artifact's own frozen rows | none | nothing; and `needle_bore.md`'s licence sentence was already scoped by item 23 |
| `tests/test_profiling_wrapper.py:369` | `>= 115` | none | nothing; written as an inequality |
| `tests/test_casson_needle.py`, `test_model_comparison_control.py`, `case_studies/blood_injectability/*` | `model.name == "casson"` | none | nothing |
| `scripts/item11_gate_report.py` | frozen by item 23 | none | nothing |
| `scripts/openfoam_couette_anchor.py:1110` | `load_all()` then indexes an explicit `RECORD_IDS` | none | nothing; the explicit list is the sample |

**Every breakage above is LOUD** — a count assertion fires — which is luck rather than
design, since item 22's fifth consumer was silent. The selectors are therefore being fixed
in kind, not just made to pass: a test that means "the Carbopol dataset" must say so by id,
so that the next family addition cannot depend on an assertion happening to fire.

Rule, unchanged: **no existing gate may change its VALUE.** Any gate whose sample would grow
gets the frozen list from `tests/conftest.py`.

## 5. What this item does NOT do

- No injectability number and no force. §0.
- No temperature MODEL and no activation energy. The shift factor is supplied, and here it
  is supplied from the deposit's own published viscosities.
- No record from the 20 deg C block, and no digitisation of anything.
- No three-component model. Adding one is a separate item with its own pre-registration.
- No change to `dst_fit`, and no change to any existing record.
