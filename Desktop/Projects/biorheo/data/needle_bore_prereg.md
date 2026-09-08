# Item 18 pre-registration — the needle-bore exponent, and whether a force budget can be met at 95 %

**Committed before any of this was computed.** Bars are not moved. A bar that turns out to
be the wrong bar is recorded with its number and a *stricter* test is added beside it. A
pre-registration is never rewritten; corrections go in an addendum at the foot of this file.

## 0. Why this item exists, and it is a commercial reason

`../../biofluids_venture/docs/market_research_2026-09-06/00_SYNTHESIS.md` ranks the
venture's first product as an **injectability calibration study**: a fitted law with a
confidence interval, an identifiability verdict, and a needle force with an error bar, sold
as a fixed-price deliverable. Report 01 of that set records two dated, URL-sourced failures
that the deliverable is supposed to address:

- *"Deviations from theoretical models were identified in drug product concentrations above
  230 mg/mL and in thin needles (27G, 30G) potentially due to ... insufficient
  characterization of high shear viscosity behavior"* — JHU thesis,
  <https://jscholarship.library.jhu.edu/items/32605470-29f4-42c5-a309-a530919d5d64>
- *"using a Newtonian model will overestimate injection time should shear thinning occur"*;
  a narrower needle *"increases dynamic pressure for Newtonian fluids but decreases for
  shear-thinning"* — Roche, *Med Devices (Auckl)* 2018,
  <https://www.tandfonline.com/doi/full/10.2147/MDER.S151727>

Both statements are about the same derivative. This item measures it.

**The quantity.** At fixed volumetric flow rate `Q`, needle length `L` and syringe barrel,

```
b  =  d ln F / d ln R
```

by central difference in `ln R`. With `friction = 0` this is `d ln dp / d ln R`; the
friction term is a property of the syringe, not the fluid, and is excluded so that `b` is a
property of the fluid and the operating point alone.

**Two closed forms bound it and neither passes through the force solver.**

- Newtonian: `dp = 8 mu L Q / (pi R^4)`, so **`b = -4` exactly**, for every `mu`, `Q`, `L`.
- Power law (`sigma = k gd^n`): `Q = (pi n / (3n+1)) R^3 (tau_w/k)^(1/n)` with
  `tau_w = dp R / (2L)`, so `Q ∝ R^(3 + 1/n) dp^(1/n)`, hence `dp ∝ Q^n R^-(3n+1)` and
  **`b = -(3n+1)` exactly**. At `n = 1` this returns `-4`; for a shear-thinning fluid
  (`n < 1`) it is **greater than** `-4`, i.e. less steep.

So the physical claim under test is one-sided: *a shear-thinning fluid pays a smaller force
penalty for a finer needle than the Newtonian `R^-4` rule says it does.* The registry
contains `newtonian` and `power_law` records for all 20 Horner donors (item 11), written as
null models, which makes both closed forms available as anchors on the same curves the
fitted families use.

## 1. Geometry is declared, not standardised

**The bore axis is DIMENSIONLESS.** This repo retracted an ISO 9626 comparison once
(item 4) because the standard was quoted from memory and mixed a regular-wall inner
diameter with a thin-wall one; the retraction is on `case_studies/blood_injectability/results.md`
and the standing rule in `NEXT_STEPS.md` is not to reinstate it without the standard in hand.
That has not changed and nothing here needs it.

The bore axis is therefore the ratio

```
rho = R / R_ref,     R_ref = 1.05e-4 m
```

`R_ref` is the *same declared value* the blood case study already runs at (0.210 mm inner
diameter, carried there with its own `[UNVERIFIED]` flag). Every result of this item is a
statement about `rho`; the absolute radius cancels out of `b`, and the one place an absolute
force appears it is labelled as being at the declared `R_ref` and inherits that flag.
`rho = 0.7` is used as the "much finer needle" corner. **It is not called 30G anywhere.**

Needle length `L = 12.7e-3` m and syringe barrel radius are inherited from the blood case
study unchanged, so that any force printed here is comparable with one already shipped.

## 2. Flow-rate axis

Two declared rates:

- `Q1 = 1 mL/min` — the rate the shipped blood study runs at.
- `Q2 = 1 mL / 10 s` — sourced: *"a 10-second injection of 1 mL through a 27 G needle
  presents shear rates up to 10^5 s-1"*, <https://pmc.ncbi.nlm.nih.gov/articles/PMC13119207/>
  (report 01). This is 6x `Q1`.

## 3. Force budget

The one budget figure available with a source is a **rule of thumb, not a standard**, and is
labelled as such everywhere it is printed: *"27-gauge thin walled needle, the liquid viscosity
should be maintained at or below 20 cP in order to not exceed the glide force of 25 Newton"*
— JHU thesis, URL in §0. **`F_BUDGET = 25.0` N.** No claim is made that any regulator or
standard uses it; it is a scale against which a derating factor can be quoted, and if it is
wrong every derating number scales with it in a way the report states explicitly.

## 4. Gates

All are graded over the registry records that carry a fitted flow curve and an
identifiability sidecar. The **frozen record list** for this item is the 20 Horner donors in
each of six families plus the 15 Carbopol Herschel-Bulkley records, read from a list private
to this item's script — **not** from `RECORD_IDS`, because appending to a shared list has
silently re-decided a closed gate in this repo before (item 11).

| gate | statement | bar | prediction |
|---|---|---|---|
| **G18.1** | positive control: `b` on the `newtonian` records equals `-4` | `|b + 4| < 1e-3` on **20/20** | PASS |
| **G18.2** | closed-form anchor: `b` on the `power_law` records equals `-(3n+1)` from that record's own fitted `n`, computed without the force solver | relative error `< 1e-3` on **20/20** | PASS |
| **G18.3** | the discriminating half: the routine does not return `-4` for everything | at least **one** record with `|b + 4| > 0.1` | PASS (the Carbopol gels have `n ≈ 0.4`) |
| **G18.4** | the one-sided physical claim: no record is steeper than Newtonian | `b >= -4 - 1e-6` on **every** record graded | PASS |
| **G18.5** | **SPLIT PREDICTION.** Using the Newtonian `R^-4` rule instead of the fitted law over-predicts the force at `rho = 0.7` by more than **10 %** | **FAIL on blood** (< 10 % on ≥ 15/20 Casson donors) and **PASS on Carbopol** (> 10 % on 2/2 gels) | as stated |
| **G18.6** | derating: the flow rate at which the 95 % upper force bound hits `F_BUDGET` is lower than the one at which the nominal force does | ratio `Q_max(95%) / Q_max(nominal) < 1` on **every** record with a closed interval | PASS |
| **G18.7** | the three error terms are put in one table and ranked, per record: bore tolerance (`|b| * delta_R/R` at a declared 5 % bore error), parameter uncertainty (the item-14 push-through hull), and constitutive family choice (item 12's Cross-vs-Casson difference) | reported, not barred — but the report must name a **winner per record** and state how often each wins | no prediction |

### 4.1 Why G18.5 is split, and what would refute each half

The blood half is predicted to FAIL *because of a result this repo already has*: at the 27G
needle the Casson force is proportional to the plateau viscosity and to nothing else
(`|d ln F / d ln eta_c| = 0.995`, item 14 G14.6). A fluid sitting on its infinite-shear
plateau is Newtonian there, so shrinking the bore cannot thin it further and `b` must be
close to `-4`. If the blood half instead comes back above 10 %, either that sensitivity
result or this reasoning is wrong, and finding out which is the item.

The Carbopol half is predicted to PASS because those records are Herschel–Bulkley with a
fitted `n` near 0.4 and no plateau at all, so `b ≈ -(3n+1) ≈ -2.2` and the Newtonian rule
over-predicts by roughly `0.7^(-4) / 0.7^(-2.2) ≈ 2.1`, i.e. 110 %. If it comes back under
10 % the closed form in §0 is wrong or the gels are not where the fits say they are.

**Both halves passing would make G18.5 uninformative** — it would then be a statement that
every fluid behaves the same way, which is the outcome this repo has repeatedly had to
label "both candidates pass". If that happens it is recorded as UNINFORMATIVE and a
sharper test is added beside it, not a moved bar.

### 4.2 What G18.6 is not

`Q_max` is where a *model* crosses a *rule-of-thumb* budget with the omissions listed in
`injectability.py`'s docstring still omitted — no barrel pressure drop, no entrance loss, no
wall slip, no needle friction. It is a derating **factor** that is meaningful as a ratio and
an absolute number that is not. The report says so beside every printed value.

## 5. What this item does NOT claim

- Nothing about any real needle gauge, any standard, or any real device. §1.
- Nothing about a mAb or any protein formulation. The registry has none; the wedge fluid is
  not in the data and the report says the demonstration runs on blood and a microgel.
- Nothing about temperature. There is no temperature axis in any record.
- No new record is added. Item 18 writes no record and fits nothing new.

## 6. Artifacts

- `src/biorheo/apps/design_space.py` — `bore_exponent`, `newtonian_rule_force`,
  `max_flow_rate`, `DeratingResult`. Library code, so a client study can call it.
- `scripts/needle_bore.py` -> `data/needle_bore.md` / `.json`.
- `tests/test_design_space.py` — pins the gate verdicts and asserts G18.3, the
  discriminating half, separately from the anchors.
