# Pre-registration — NEXT_STEPS item 12: does a donor's Cross record give the same needle force as its Casson record?

Written **2026-09-06, before any of the comparison was computed**, and committed on its
own so a later session cannot claim a bar was chosen after the answer. Same discipline as
`data/openfoam/nonyield_bridge_prereg.md` (item 6, `2561346`) and
`data/openfoam/newtonian_powerlaw_prereg.md` (item 11, `e9c1ace`). Both of those items had
a stated premise that turned out to be FALSE, so §0 checks this one's premise first.

## 0. The premise, reproduced before anything is built on it

Item 12's stated premise is *"the yield term contributes under 3 %"*. That is on the case
study's own page. Re-run of `case_studies/blood_injectability/run.py` on 2026-09-06
reproduces it: setting `tau_y` to zero and refitting the pressure moves the plunger force
by **0.92 % to 2.94 %** across the 20 donors (`yield_share_min`, `yield_share_max` in
`results.json`), and `x = tau_y/tau_w` runs 1.63e-05 to 1.68e-04. **The premise is TRUE.**

But the same page carries a second number that bounds what this comparison can mean, and
item 12 does not mention it. **The entrance (Hagenbach) loss, which the model omits
entirely, is 132.5 Pa — 0.43 % to 0.99 % of the pressure drop.** Fahraeus-Lindqvist is a
second omitted correction of unquantified size pointing the other way. So:

> **The calculation's own unresolved systematic is of order 1 % in force.** Any agreement
> between two models that lands below about 1 % is agreement below the resolution of the
> model, and must be reported as such rather than as a corroboration.

That is the single most important thing this pre-registration adds, and it is fixed here
before the numbers exist.

## 1. What is compared, and why it can fail

For each of the 20 Horner donors, the plunger force is computed at the SAME operating
point — `Q` = 1 mL/min, the nominal 27G needle (ID 0.210 mm, L 12.7 mm), the same syringe
— from that donor's records in five families, all fitted to the SAME 16 points:

| family | route | high-rate behaviour |
|---|---|---|
| `casson` | closed form (`flow_rate_casson`) | -> `eta_c`, a fitted plateau |
| `cross` | quadrature (`flow_rate_gnf`) | -> `eta_inf`, a DIFFERENT fitted plateau |
| `carreau_yasuda` | quadrature | -> `eta_inf`, a third fitted plateau |
| `newtonian` | quadrature | constant `mu`, fitted over the WHOLE curve |
| `power_law` | quadrature | `k gd^(n-1)`, no plateau: keeps thinning |

**This is not a formality and it can fail.** The needle's wall rate is about 1.8e4 1/s,
**1.42 decades above the 700 1/s top of every record's envelope**, so every family is
being read where none of them was measured. Casson, Cross and Carreau-Yasuda all approach
a constant there, but they approach **three different constants**, each of which is partly
an extrapolation of a fit that stops at 700 1/s. Whether those three constants agree is a
real question with a real answer, and "they fit the same curve" does not settle it.

`newtonian` and `power_law` are the NULLS and they are the reason this comparison is
worth running at all (see §3).

## 2. Gates. These are the bars. They are not to be moved.

Every number is printed by `case_studies/blood_injectability/crosscheck.py` and pinned in
`tests/test_injectability_crosscheck.py`. A failed gate is recorded with its number and a
*stricter* test added beside it; it is never loosened.

The metric throughout is the **relative force difference against the donor's own Casson
force**, `F_family / F_casson - 1`. The denominator is a declared, non-collapsing scale:
Casson forces run 0.223 to 0.509 N over the 20 donors (`results.json`), so this is not
item 9's divide-by-dust. The raw forces are stored in the artifact so a reader sees
numbers and not only a ratio.

**Three scales are named in advance**, all from the case study's own artifact:

- **1 %** — the model's unresolved systematic (entrance loss, §0). Below this, two models
  are not distinguishable BY THIS CALCULATION.
- **2.94 %** — the largest yield-stress contribution, the effect the case study calls
  negligible. An agreement at or under this is "no larger than the term we already
  neglect".
- **128 %** (a factor 2.2792) — the donor-to-donor force spread, the quantity item 12
  proposes as the tolerance. **This is a nearly vacuous bar** and is recorded as such:
  almost any pair of fits to the same curve would pass it.

**G12.1 (the item's own prediction, stated at its own bar).** For every donor,
`|F_cross/F_casson - 1|` is below the donor-to-donor spread. **Bar: < 1.2792** (the
spread, 2.2792, minus 1). **Pre-registered as WEAK**: this is the bar item 12 proposes,
it is reported for that reason, and passing it establishes very little.

**G12.2 (the informative bar).** The worst `|F_cross/F_casson - 1|` over the 20 donors is
compared with the three scales above and the verdict is one of four, decided by which band
it falls in — the bands, not the verdict, are what is fixed here:

| band | verdict |
|---|---|
| < 1 % | AGREE BELOW RESOLUTION — the two models are not distinguishable by this calculation, and that is the honest statement, not "they agree" |
| 1 % to 2.94 % | AGREE TO WITHIN THE NEGLECTED YIELD TERM |
| 2.94 % to 128 % | DISAGREE BY MORE THAN THE EFFECT THE CASE STUDY CALLS NEGLIGIBLE, but less than the donor spread — the case study's ranking survives, its precision does not |
| > 128 % | FAILS — the choice of constitutive family matters more than which patient the blood came from, and the case study's central claim needs restating |

Same table applied independently to `carreau_yasuda`.

**G12.3 (the informativeness control — the gate that decides whether G12.2 means
anything).** A comparison in which every candidate passes has discriminated nothing;
item 6's G7 is this repo's worked example. So, in advance:

> **G12.2 is INFORMATIVE only if at least one of the two NULL families
> (`newtonian`, `power_law`) lands in a DIFFERENT band from `cross`.**

If `newtonian` — a single constant fitted to the whole 0.1-700 1/s curve, which is not the
high-rate plateau — also reproduces the Casson force to under 1 %, then the needle has
washed out every difference between every model, the agreement of Cross and Casson says
nothing about their plateaux, and G12.2 is reported UNINFORMATIVE regardless of its number.
**Bar: the null must separate.** The separation ratio (worst null deviation / worst Cross
deviation) is printed beside the verdict, and the comparison is called uninformative below
**1.6x**, the same floor this programme has used since 2026-08-30.

**G12.4 (identifiability split, decided before it is looked at).** 7 of the 20
Carreau-Yasuda records carry a `bounded_only` `a`, and some Cross records carry a
`degenerate` `lam`. All 20 are run, and the comparison is reported TWICE: over all 20, and
over the subset whose every parameter is `measurable`. **Bar: the verdict of G12.2 is the
same on both subsets.** If it is not, the headline is the all-`measurable` one and the
difference is the finding. Counts before and after the split are printed.

**G12.5 (nothing is silently dropped).** `flow_rate_gnf` refuses a non-monotone flow curve
(`_assert_monotone_stress`), which a Cross record with `m > 1` can produce, and
`pressure_for_flow_rate` can run out of bracket. Every such failure is CAUGHT, counted,
named by donor and family, and printed. **Bar: N in and N out are both reported for every
family, and any family losing a donor is reported as a partial result rather than
compared as if complete.** A dropped donor is data about the fit, not a nuisance.

## 3. What this cannot establish, written in advance

- **It cannot validate the extrapolation.** All five families are read 1.42 decades past
  the data on every donor. If they agree, they agree about an extrapolation; the
  measurement that would settle it is tube viscometry at this diameter, which the Horner
  dataset does not contain. No amount of agreement here is evidence that any of them is
  RIGHT at 1.8e4 1/s.
- **It cannot resolve anything at the 1 % level** (§0), and the entrance and
  Fahraeus-Lindqvist corrections point in OPPOSITE directions, so they do not even
  cancel predictably.
- **A pass of G12.1 is close to vacuous** and is labelled that way on the page.

## 4. What would make this item's work wrong, stated in advance

- Reporting "the models agree" for a number below 1 % without saying it is below the
  calculation's own resolution.
- Grading only against the donor spread (G12.1) because that is the bar item 12 wrote.
- Running the nulls and not letting them decide informativeness (G12.3).
- Dropping a donor whose flow curve the quadrature refuses, instead of counting it.
- Forming `F_family/F_casson - 1` without checking the denominator's scale — it is fine
  here (0.223-0.509 N) and that is checked, not assumed.
