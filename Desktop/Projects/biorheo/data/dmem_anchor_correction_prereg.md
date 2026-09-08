# Item 23 pre-registration — the item-22 anchor was graded against the wrong published numbers

Written 2026-09-08, before any change to `scripts/records_from_dmem.py`, after three
independent reviewer passes on item 22 (a code review, a physics referee, and a numerics
re-run). Two of the three found the same thing without seeing each other's reports.

**This pre-registration is unusual and the difference is stated here rather than buried.**
Every other pre-registration in this repo was written before its numbers were known. This
one is written after a review measured most of them. Section 0 says which gates are
**restatements** of a number already in hand and which are genuinely **blind**. A
restatement is not evidence in the way a blind gate is, and calling one the other would be
the exact failure this file exists to correct.

`data/dmem_record_prereg.md` is **not edited**. An addendum at its foot points here.

---

## 0. What is being corrected

Item 22 shipped this claim, and it is in the Status card, `NEXT_STEPS.md`, `OFFER.md` and
the project's memory index:

> The source paper's authors fit the same Ostwald-de Waele power law and published their
> flow indices in 2024. Fitting the raw deposited files through `biorheo.identify`
> reproduces them on 4/4 conditions -- 0.0005 on the 10 % FBS 20 C condition, worst 0.0299.
> Every other anchor here compares biorheo against biorheo, a closed form, or a solver;
> this one compares it against a number somebody else published before the code existed.

Four defects, in descending order of how much they move the conclusion.

### 0.1 The graded values are not the paper's fit-results table (FATAL to the claim as worded)

`SPEC` in `scripts/records_from_dmem.py` grades against `0.32, 0.31, 0.54, 0.50`. Those
come from **figure captions** -- Fig. 4, Fig. 5 and S1 Fig -- which state the flow index
used as the **input to the COMSOL simulations**, rounded and **shared across conditions**
(S1 Fig's caption gives one `n = 0.31` for "10 % FBS *and* 5 % FBS at room temperature" and
one `n = 0.54` for "10 % FBS at 37 C *and* 1 % FBS at room temperature").

The paper also has **Table 2, "Fit results for `K` and `n` values in (4) for the
experimental rheometry data shown in Fig. 1"** -- one row per sample, resolving exactly the
four conditions the deposit contains, for exactly the equation biorheo fits:

| sample | `K` | `n` |
|---|---|---|
| DMEM + 10 % FBS | 0.028 | 0.31 |
| DMEM + 5 % FBS | 0.037 | 0.33 |
| DMEM + 1 % FBS | 0.016 | 0.51 |
| DMEM + 10 % FBS (37 C) | 0.005 | 0.54 |

Three of the four values `SPEC` uses disagree with this table. This is the same failure
already in `lessons.md` from 2026-09-01 -- *"Read the sample table before comparing to
it"*, written after a contradiction with Singh et al. was built entirely out of comparing
against an abstract instead of the paper's Table 2. It has now happened twice.

**Table 2 also publishes `K`, and item 22 never compared it.**

### 0.2 The replacement gate pooled an interval the same file excludes

`SPEC` leaves `10pct_37C_b` unbarred with a reason stated in a comment: the paper's indices
are for the 0.01-1 1/s range and this is the only interval running to 10 1/s. G22.1' then
concatenates that same interval into the 37 C pooled fit, where it supplies 20 of the 25
points. The exclusion rule and the gate that replaces it contradict each other.

### 0.3 G22.1's pass expression could not be true

`"pass": bool(len(graded) == 6 and worst < G221_TOL)` with `len(graded) == 7`. The gate
reads FAILED in the artifact for a reason that has nothing to do with any fit. The recorded
narrative attributes the failure to the tolerance.

### 0.4 The number offered as proof that G22.1 was mis-specified is from the wrong condition

`worst_spread` is a max over all four conditions including 10 % FBS 37 C, only one of whose
intervals is graded -- so it cannot illustrate "no fit could have met the bar on **both**
intervals of that condition". The generator's own comment says 0.076 while the code prints
0.1044.

---

## 1. Gates

`RESTATEMENT` = the number was measured during the review and is recorded here so the
artifact and the prose cannot drift. `BLIND` = not yet computed at the time of writing.

### G23.1 `RESTATEMENT` — the corrected anchor FAILS the original bar, and the bar is not moved

Grade the four per-condition pooled fits against **Table 2**, pooling only intervals that
lie in the paper's stated rate range, at the **original `G221_TOL = 0.05`**.

Predicted: FAIL on exactly one condition, 1 % FBS 20 C, at `|dn| = 0.0599`; worst over the
four = 0.0599. The other three pass at 0.0095, 0.0010, 0.0280.

**The bar stays at 0.05 and the failure is reported.** Widening it to 0.06 to collect a
pass is the move this repo does not make.

### G23.2 `RESTATEMENT` — the old pass does not survive either correction, separately

Two independent single-change reruns of the ORIGINAL G22.1' (still against `SPEC`, not
Table 2):

- drop the out-of-range 37 C interval from the pooling: worst goes 0.0299 -> **0.0651**, FAIL.
- keep the pooling, grade against Table 2: worst goes 0.0299 -> **0.0599**, FAIL.

Both must be reported. A claim that survives neither of two independent corrections was not
resting on the data.

### G23.3 `RESTATEMENT` — the flagship 0.0005 is a pooling coincidence, not a precision

The 10 % FBS 20 C pooled index is 0.3195 while its own two interval fits are 0.2945 and
0.3134. **A pooled estimator outside the span of its inputs** means the two intervals are
mutually inconsistent under one power law, so four-figure agreement with any external
number is a coincidence of combining them. Predicted: `pooled_inside_span` is False for
this condition and True for the other two two-interval conditions.

This number must stop being quoted as a precision anywhere in the repo.

### G23.4 `RESTATEMENT` — what actually survives is a comparison item 22 never made

Compare the two **parameter sets as viscosity curves** over each condition's own measured
rate range. This is immune to the `K`-`n` tradeoff that makes a deviation in `n` alone
uninterpretable, and it is the quantity a client uses.

Predicted: median 2.04 / 0.99 / 7.06 / 4.86 %, worst 4.16 / 1.23 / 19.03 / 11.80 %.

**Scale it against a declared number, not against a tolerance this repo invented.** The
paper's S2 Fig states its own repeatability: *"A 20 % deviation from the fit line is
included as a measure of repeatability of results."* Predicted: all four conditions agree
to inside that 20 %, the worst being 19.03 %.

`K` alone, predicted: 0.1 / 0.8 / 6.7 / 5.9 % relative, and 3 of 4 fitted values round to
the published value at the precision the paper prints it.

### G23.5 `BLIND` — does the correction move any record?

Predicted **no**: the fits are untouched and only the grading changes, so all eight
`dmem_*.yaml` and their sidecars must be byte-identical after regeneration. If any record
moves, this item has a second defect and the prediction is wrong.

Read on a continuous channel, not on a verdict: compare the YAML bytes, not the verdicts.

### G23.6 `BLIND` — does the correction reach any shipped number?

Predicted **no**: no DMEM record enters a force, needle or design-space calculation.
Settled by grep over `case_studies/` and `apps/`, not by argument. If it does reach one,
say which and by how much.

This is the item-13/15 shape -- a correction reaches a downstream number only through the
parameters that number uses -- and the prediction is that the reach here is zero.

### G23.7 `BLIND` — the closed item-11 gate must be regenerable at its own sample

`scripts/item11_gate_report.py` calls `load_all()`, so it now recomputes a CLOSED gate over
123 records where its committed artifact says 115. Predicted: pointing it at the frozen
115-record list reproduces `data/openfoam/newtonian_powerlaw_gates.md` **bit-for-bit**. If
it does not, the artifact was already stale before item 22 and that is a separate finding.

### G23.8 `BLIND` — is `K` identifiable here at all?

A comparison against a published `K` means nothing if `K` is unidentified. Predicted: every
DMEM record's `k` carries a closed profile interval and is verdict `identified`. If any is
not, the `K` comparison is reported as `bounded_only` for that condition.

---

## 2. What is NOT being changed

- **`G221_TOL` stays 0.05.** Every failure above is reported at that bar.
- **`data/dmem_record_prereg.md` is not edited**, per the repo's standing rule; it gets an
  addendum.
- **G22.1 and G22.1' keep their original numbers in the artifact**, including the 0.0299,
  labelled as superseded and with the reason. A retraction that deletes the number it
  retracts cannot be audited -- that is the 2026-09-06 lesson about item 9's retraction
  quoting a figure that was itself the artefact.
- **The eight records stay.** The defect is in the grading, not in the fits, which is
  exactly what G23.5 tests.

## 3. What the honest headline becomes

Not *"reproduces the paper's published flow indices on 4/4 conditions, worst 0.0299"*.

Instead, and only if the gates come out as predicted:

> Fitting the raw deposited curves through `biorheo.identify` reproduces the source paper's
> own Table 2 `(K, n)` as a **viscosity curve** to a median of 1.0-7.1 % over each
> condition's measured range, worst 19.0 %, inside the 20 % repeatability the authors
> state for their own measurements. On the flow index alone it agrees to 0.03 on three of
> four conditions and **misses the fourth by 0.0599 against a pre-registered 0.05**. `K`
> agrees to 6.7 % or better on all four, and 3 of 4 round to the published value.

Weaker than the claim it replaces, better anchored, and it names its own failure.

## 4. A limit on what any version of this anchor can show

Recorded now so it is not discovered later and presented as new.

1. **The fit method is shared.** For a power law, minimising rms `log(eta)` residuals is
   ordinary least squares of `log eta` on `log gd`. That is the near-universal way anyone
   extracts `(K, n)`, and the paper gives no sign of doing anything else. Agreement
   therefore certifies this repo's **parsing** -- the interval split, the `gd > 1e-3` drop
   rule, the unit handling -- and does **not** exercise the multi-start search or the
   profile likelihood, which this problem is too well conditioned to strain.
2. **The paper's own scatter is the same size as the agreement.** Its figure captions and
   its Table 2 differ by 0.01-0.04 on the same conditions, and its Table 2 fit is over an
   ensemble of three samples, refills and different days, while the deposit is one file per
   condition. So the reference number is not a point.
3. **The dataset's own noise floor exceeds the agreement on half the conditions.** Repeat
   intervals of the same condition move `n` by 0.075 (5 % FBS) and 0.061 (1 % FBS), larger
   than the deviations being called agreement there. Only the 10 % FBS 20 C condition's
   agreement is well inside its own interval-to-interval scatter.

None of this makes the anchor worthless -- it is still the only comparison in this repo
against a number computed by someone else, on data this repo did not process. It makes it
an anchor on the **parse and the units**, at the resolution of the source's own
repeatability, which is what it will be called from now on.
