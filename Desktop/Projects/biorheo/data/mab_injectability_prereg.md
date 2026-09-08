# Item 25 pre-registration — the wedge fluid, at last: a real high-concentration mAb

Written 2026-09-08, before any analysis code. Data downloaded and its provenance recorded in
`data/mab_prasnikar2025/SOURCE.md`.

**This is the first time anything in this project has been demonstrated on the fluid the
venture is actually about.** `OFFER.md` §6 names it as gap 3: *"No protein or
high-concentration formulation curve, so nothing here has been demonstrated on the wedge's
own fluid. The examples are human blood and a microgel."*

Source: **Zenodo 14409284, CC BY 4.0**, supporting Table S1 of Prašnikar et al.,
*"Additive effects of the new viscosity-reducing and stabilizing excipients for monoclonal
antibody formulation"*, Int. J. Pharm. **674** (2025) 125451, `10.1016/j.ijpharm.2025.125451`
(open access; the accepted PDF is in the University of Ljubljana repository).

**58 formulations of one model monoclonal antibody at 159.3–212.2 mg/mL**, each with a
measured viscosity (13.57–47.07 mPa·s) and its standard deviation, plus the excipient and
its concentration, plus SEC stability data. The paper's own keywords are *"subcutaneous
administration, highly concentrated mAb formulations, viscosity-reducing agent"*.

---

## 0. Two facts from the source that decide what may be done, and one that decides what may not

**0.1 The fluid is Newtonian, and the authors established that rather than assuming it.**
Methods §2.2.4: *"we measured the viscosity of highly concentrated protein solutions at
shear rates of 960–3329 s⁻¹ … categorized them as Newtonian. Consequently, samples were
measured at 25 °C and shear rates within the described range to ensure a constant flow
rate."* So a one-parameter Newtonian model is the source's own verdict over a **stated rate
range**, not a convenience.

**0.2 That rate range is needle-adjacent.** 960–3329 1/s against 1e4–1e5 at the wall of a
fine needle. Every other fluid in this project is decades away: the DMEM records stop at
1 1/s (four to five decades), the Carbopol-in-propylene-glycol records at 1e3 1/s, and the
blood donors extrapolate 1.42 decades. **G25.4 measures where this one lands** and it is
expected to be the least-extrapolated fluid in the registry's history.

**0.3 These cannot become registry records, and G25.7 demonstrates that rather than
asserting it.** A record needs a fit and an identifiability sidecar written by
`biorheo.identify`. A Newtonian model has one parameter and each formulation here has **one
viscosity**: zero residual degrees of freedom, no profile interval, nothing for
`write_identifiability` to write. This is a real boundary of the schema and it is worth
stating plainly, because **it is also what a client's data usually looks like** — a table of
formulation viscosities, not a flow curve. The case study therefore reads the deposited
table directly, which is exactly the arriving-data path a paid study would take.

---

## 1. Gates

`RESTATEMENT` = measured while assessing the dataset, recorded here so prose and artifact
cannot drift. `BLIND` = not computed at the time of writing.

### G25.1 `RESTATEMENT` — the source's `reduction %` column does NOT reproduce, so nothing is built on it

The paper defines it: *"the viscosity reduction was calculated as a relative decrease in
viscosity (%) of the sample containing the test compound compared to the starting
formulation without the test compound at the same mAb concentration, which was calculated
using the exponential calibration curve"* — i.e. `red = 100·(1 − η/η_base(c))` with
`η_base = a·exp(b·c)` (Mooney).

Back-solving `η_base = η/(1 − red/100)` and fitting one exponential gives
**R² = 0.4575 on the 49 rows of the main series** (rms 0.215 in ln η). A single calibration
curve would collapse this to numerical precision. Either more than one curve was used, or
the column carries processing the Methods do not describe.

**Consequence, and it is the whole point of the gate: this study uses the MEASURED
viscosities and never the `reduction %` column.** A number that cannot be reproduced cannot
be extended. Predicted: R² < 0.7, i.e. FAIL to reproduce, recorded rather than worked around.

### G25.2 `RESTATEMENT` — but the Mooney SLOPE is robust, and it is the lever G25.5 needs

Fitted independently on the two label series, the exponent is **0.03131 mL/mg** (49 rows,
c = 159–184) and **0.03233 mL/mg** (9 rows, c = 186–212). **Bar: the two agree to within
20 %.** The scatter that sinks G25.1 is in the intercept, not the slope.

Predicted PASS, and it licenses one statement: **viscosity doubles every `ln2/b` mg/mL**,
about 22 mg/mL.

### G25.3 `BLIND` — the force, on the same device as the blood study

Compute the plunger force for all 58 formulations through **exactly the needle, syringe and
flow rate the blood case study uses**, changing nothing about the device. Reusing them
introduces no new unverified constant and makes the two studies directly comparable — and
the comparison is the point, because **this repo has never put its wedge fluid and its
demonstration fluid on the same axis.**

**No external injectable-force limit is used.** This project retracted an ISO 9626
comparison once for quoting a standard from memory, and a human-factors force limit would be
the same mistake. The scale is the blood donors' own forces through the same needle.

Predicted: every mAb formulation needs more force than every blood donor, because 13.6–47.1
mPa·s against blood's few mPa·s at needle rates. **Bar: report the ratio; if any mAb
formulation needs LESS force than the median blood donor, this prediction is refuted and
said so.**

### G25.4 `BLIND` — how far the device extrapolates past what was measured

`wall_shear_rate` at the delivery flow rate, against the source's own 960–3329 1/s.

**Bar: report the extrapolation in decades and compare with the three fluids already in the
project** (blood 1.42, Carbopol-in-PG, DMEM). Predicted: **under 1.5 decades, and the
smallest of any fluid here.** If it is larger than blood's 1.42 the headline of §0.2 is
withdrawn.

This is item 5's finding — *how far a device extrapolates its own rheology is a property of
the DEVICE* — applied to the fluid the venture is for.

### G25.5 `BLIND` — the headline: which lever moves the force more, the excipient or the concentration?

The question a formulation scientist actually asks. Both levers are measurable here:

* the **excipient** lever, from the spread of measured viscosity across formulations at
  comparable mAb concentration;
* the **concentration** lever, from G25.2's Mooney slope.

Express the excipient's effect in the concentration's units: **how many mg/mL of extra mAb
concentration does the best excipient buy at equal viscosity?** That is
`ln(η_worst/η_best)/b`.

**No bar; this is a measurement, and it is the number the study exists to produce.**
Predicted: tens of mg/mL, i.e. the two levers are comparable rather than one dominating. If
the excipient lever turns out negligible against a few mg/mL of concentration drift, that is
a more useful finding than the reverse and it must be reported that way.

⚠ **Confound to control, or the answer is arithmetic rather than physics.** The formulations
do NOT sit at one mAb concentration — they span 159–212 mg/mL, and viscosity is exponential
in concentration, so the raw viscosity spread mixes both levers. **The excipient lever must
be measured at matched concentration**, by normalising each formulation to a reference
concentration through the Mooney slope before taking the spread. Report both the raw and the
normalised spread; if they differ materially, the raw one is the confounded quantity and the
normalised one is the answer. This is the pooled-vs-within-cell lesson.

### G25.6 `BLIND` — is the measurement or the formulation the larger term?

Each viscosity carries a measured SD. Propagate it to the force (for a Newtonian fluid the
force is exactly proportional to viscosity, so the relative error passes straight through)
and compare with the formulation-to-formulation spread.

Predicted: the formulation spread is far larger, i.e. the measurement is not the limiting
term — the opposite of item 14's blood finding, where the parameter uncertainty was 29 % of
the force and the donors could not be ranked pairwise. **Bar: report the fraction of the
`58·57/2` formulation pairs whose force intervals are disjoint at 2 SD.** Item 14's blood
answer was 36 of 190 pairs; predicted here: a much higher fraction.

### G25.7 `BLIND` — demonstrate the schema boundary rather than asserting it

Attempt to build a registry record from one formulation and record what fails and where.
Predicted: `fit_rate_curve` cannot fit a one-parameter model to a single point with zero
residual degrees of freedom, so no identifiability sidecar exists and the loader would
reject the record.

**Bar: the failure must be shown, with the exception, not argued.** If it unexpectedly
succeeds, the registry gains 58 licensable wedge-fluid records and §0.3 is wrong.

---

## 2. What this item does NOT do

- **No claim that any formulation is or is not injectable by a human.** No force limit is
  quoted, from a standard or otherwise. The scale is this repo's own blood study.
- **No use of the `reduction %` column** (G25.1), and no use of the SEC stability columns —
  stability is not rheology and this package has nothing to say about it.
- **No new needle, syringe or flow rate.** The blood study's device is reused verbatim.
- **No registry records** (G25.7), and therefore no change to any existing gate's sample.
- **No shear-thinning model.** The source verified Newtonian over its own measured range and
  that is the model used; whether the fluid is still Newtonian at a needle wall is exactly
  what G25.4's extrapolation number is for, and it is not assumed either way.
