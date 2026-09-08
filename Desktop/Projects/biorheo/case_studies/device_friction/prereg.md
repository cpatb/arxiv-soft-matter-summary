# Item 28 pre-registration — the term that defaults to zero

Written 2026-09-08, **before any number in this study was computed**. Not to be edited;
corrections go in an addendum at the bottom, as in `newtonian_powerlaw_prereg.md` §0.2 and
`profile_interval_prereg.md` addendum A.

## 0. The gap

`biorheo.apps.injectability.injection_force` computes

    F = dp * pi R_barrel^2 + friction,        friction: float = 0.0

`dp` is the fully developed pressure drop in the needle cannula. `friction` is a
**caller-supplied** term with a **default of zero**, and it is described everywhere in the
repo as "a property of the device, not the fluid" — correctly. Every shipped page says so:
the five `deliverable_example/*.md` pages each carry the sentence *"Syringe friction is a
property of the device, not the fluid; it is added only if the caller supplied it."*

**No case study supplies it.** `grep -rn friction case_studies/ --include=*.py` returns
one hit, and it is a comment in `injectability_surrogates/run.py` telling the reader to
replace the geometry with a datasheet. So every force this repo has ever published —
0.287 N for blood donor A, 1.868 N for the `Arg_25` mAb formulation, the 1.87 N headline
on `05_screen.md`, which is the **client-facing** page — is the hydrodynamic term alone.

That is honestly labelled and it is not the finding. The finding this item is looking for
is the one nobody has checked: **which of the shipped CONCLUSIONS survive a non-zero
friction term, and which are properties of the assumption that it is zero.**

The companion omission, the entrance (Hagenbach) loss, has already been sized —
`blood_injectability/results.md` puts it at 0.43–0.99 % of `dp`, the same size as the
yield-stress contribution that study calls negligible. Friction has never been sized at
all, and unlike the entrance loss it is not small: it enters as an **additive constant in
newtons**, not as a fraction of `dp`, so its importance is set by the size of `dp * A`
rather than by any dimensionless group.

## 1. What is being tested

Two shipped claims and one shipped uncertainty:

- **G25.3** — *"the median formulation needs 7.1x the plunger force of the median blood
  donor (1.98 N vs 0.28 N)"*, `mab_injectability/results.md`. A RATIO of two forces.
- **G25.6** — *"92.3 % of formulation pairs separate at 2 SD, against item 14's 18.9 % for
  blood"*, i.e. a formulation decision is supported by the data where a patient-ranking
  decision was not. A PAIRWISE separability.
- **G14.5** — *"only 36 of 190 blood donor pairs have disjoint 95 % force intervals"*.

## 2. The algebra, stated first so no gate is written that cannot fail

`F_i = h_i + f` with `h_i` the hydrodynamic term and `f` a constant. Then:

- **Any ordering claim is invariant.** `F_i > F_j` iff `h_i > h_j`. Spearman rank
  correlation before/after is exactly 1. Testing this would decide nothing, so it is
  **not a gate**; it is recorded here as the thing that survives by construction.
- **Any RATIO claim is not.** `(h_1 + f)/(h_2 + f) -> 1` monotonically as `f` grows.
- **Any DIFFERENCE claim is invariant**: `F_i - F_j = h_i - h_j`.
- **Pairwise separability under a DETERMINISTIC `f` is invariant**, because both interval
  ends shift by the same `f`. It is **not** invariant under a friction that VARIES from
  syringe to syringe, which is what a real device does. That variability, and not the mean
  friction, is the quantity that can destroy G25.6 — and it is the one no page mentions.

So the study has exactly two live questions, and they are separated on purpose:

- **Q1 (mean friction):** how far does a constant `f` move the RATIO claims?
- **Q2 (friction repeatability):** at what syringe-to-syringe friction SD `s` does the
  formulation screen stop being decidable?

## 3. Pre-registered predictions

Written before running. Each can fail.

- **P1.** The G25.3 ratio `F_mAb_median / F_blood_median = 7.1` falls **below 3.0 at
  f = 1.0 N** and **below 2.0 at f = 3.0 N**.
- **P2.** The value of `f` at which that ratio first falls below 2.0 is **under 4 N**.
- **P3 (the deliverable).** Adding a syringe-to-syringe friction SD `s` in quadrature to
  each formulation's own force SD, the fraction of the 1653 mAb formulation pairs that
  separate at 2 SD falls **below 50 % at some s* < 0.05 N**.
- **P4 (the directional one, with a mechanism).** The SAME absolute `s` damages the mAb
  screen far more than the blood ranking, because the mAb interval is a pure measurement
  repeatability (median force SD 0.0078 N, 0.42 %) while the blood interval carries a
  parameter uncertainty item 14 measured at a median 28.8 % (0.08 N absolute). Concretely:
  at `s = 0.05 N` the mAb median interval half-width grows by **more than 3x** and the
  blood one by **less than 1.3x**.
- **P5.** The blood pairwise separation `36/190` falls below **10 %** at `s = 0.1 N`.

## 4. Controls, and what each can and cannot decide

- **C1 (reproduction before extension).** At `f = 0` and `s = 0` this study must reproduce
  the shipped numbers EXACTLY: G25.3's ratio, G25.6's `disjoint/1653`, G14.5's `36/190`.
  If any of the three disagrees, the study stops and that disagreement is the finding.
  This is the only check that certifies the re-implementation; every number below is
  meaningless without it.
- **C2 (the two cohorts use different yardsticks, and that is not fixed here).** Blood
  intervals are parameter push-through hulls; mAb intervals are `+/- 2 SD` of a measured
  repeatability. They are NOT the same quantity, and this study does not make them one. It
  applies the same `f` and the same `s` to both and reports each against its OWN shipped
  baseline. A cross-cohort comparison of the two separability numbers is therefore
  reported as it already is on the pages, and inherits their existing caveat.
- **C3 (what the sensitivity curve can decide without any literature).** The deliverable is
  **s\*, the threshold**, not a verdict at an assumed `s`. A threshold needs no external
  number to be correct. A literature or datasheet value for real syringe friction only
  says WHERE on the curve a given device sits, and if no value resolves to a DOI, the
  threshold still stands and the verdict is withheld. **No friction magnitude will be
  quoted in this study from memory, and none from an unresolved source.**
- **C4 (the null that would make this item pointless).** If the answer to Q2 is
  `s* > 1 N` — i.e. the screen survives any plausible device — then the finding is that
  friction is irrelevant to these conclusions, and that is a real verdict which will be
  reported as such rather than being made to look like a defect.

## 5. What may NOT be claimed from this study

- That friction is or is not large for any particular syringe. Nothing here measures a
  syringe.
- That the shipped forces are wrong. They are the hydrodynamic term, correctly computed and
  correctly labelled. What is at issue is the CONCLUSIONS drawn from them.
- That friction is independent of the fluid. It is a plunger/barrel/lubricant interaction;
  treating it as fluid-independent is an ASSUMPTION stated here, and it is the assumption
  that makes `f` additive and therefore makes the ordering claims survive. If friction
  depended on the fluid, even the ordering would be at risk, and this study cannot tell.
- That a "total injection force" computed here is what a patient feels. Needle-hub losses,
  the start transient and tissue back-pressure are all still absent.

---

## Addendum A (2026-09-08, after the run; the body above is NOT edited)

**Two of the five predictions were REFUTED and no bar was moved.**

- **P3 REFUTED.** Predicted `s* < 0.05 N` for the mAb screen to fall below 50 % separable
  pairs; measured **`s* = 0.1397 N`**, 2.8x the bar. The direction is right and the size is
  not. The finding stands on the measured threshold, not on the prediction.
- **P4 REFUTED, and the defect is in this pre-registration rather than in the result.**
  P4 required `mab_growth > 3` (holds, 4.56x) AND `blood_growth < 1.3` (fails, 3.33x).
  Worse than a failed bar: **the two arms were specified under two different operations** —
  the mAb arm adds the device term in quadrature and the blood arm adds `2s` linearly to
  each hull side — so the comparison was not like-for-like, which is the same defect this
  repo recorded on 2026-09-05 (item 6's G7) and again on 2026-09-08 (item 23's pooling).
  The graded verdict and its numbers stand. A post-hoc recomputation under the SAME
  operation on both arms is reported beside it: blood **2.54x**, mAb **4.56x**, so the mAb
  cohort is the more exposed by **1.80x** — the direction P4 asserted, at roughly half the
  separation the bar implied. It is labelled POST-HOC in `results.json` and on the page.

**P1, P2 and P5 passed**, and **C1 reproduced all three shipped numbers exactly**
(7.052465, 1525/1653, 36/190), so the arrays audited are the shipped ones.

**C4 (the null) did not fire**: `s*` is 0.1397 N, well under the 1 N at which the honest
verdict would have been *friction is irrelevant to these conclusions*.

**What this item did NOT establish, and a later item should.** The friction term is
treated here as additive and fluid-independent (prereg §5). If real plunger friction
depends on the fluid — through the lubricant film the plunger rides on — then even the
ordering claims are exposed, and nothing here can tell. That is a measurement, not a
calculation.
