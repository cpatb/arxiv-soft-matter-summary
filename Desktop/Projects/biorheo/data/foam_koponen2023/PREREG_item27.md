# Pre-registration — item 27: an EXTERNAL tube-flow anchor, and the first measurement of the wall-slip gap

Written **2026-09-08, before the generator**, as every record in this repo has been.
A pre-registration is never rewritten; errors in it are corrected in an addendum at the
foot of this file plus a one-line marker where the error is.

---

## 0. Why this item exists

`src/biorheo/apps/injectability.py` names its own uncovered gap in its module docstring:

> ALSO NOT included, and not inertial: **wall slip**. The closed form assumes no slip at
> the needle wall; microgel suspensions (Carbopol) and many biological fluids slip at wall
> stresses of the kPa order this geometry produces, which lowers the true force below the
> number here and is not covered by the extrapolation flag either (it is a boundary
> condition, not a rate).

That sentence has never been attached to a number, because nothing in this repo has ever
been compared against a **measured pressure-driven tube flow**. Every anchor the tube
chain has is analytic: the Buckingham–Reiner limit, the Casson closed form, Poiseuille,
and a plane-Couette OpenFOAM run against an analytic profile. The claim the OFFER sells —
*a model identified from one measurement predicts the pressure drop in a different
geometry* — has therefore never been tested against an experiment at all.

This dataset can test it, and can measure the slip term, because it contains **both**.

## 1. The data

**Zenodo `10.5281/zenodo.8334762`, CC BY 4.0** (SPDX `cc-by-4.0`, confirmed from the
DataCite record, not from the landing page). Koponen, Viitala, Tanaka, Prakash, Laukkanen,
Jäsberg, *Pipe rheology of wet aqueous application foams*, **Chem. Eng. Sci.** (2023),
`10.1016/j.ces.2023.119282`; the DataCite record carries that DOI as `IsReferencedBy`.
File md5s are in `SOURCE.md`.

`Set_2_Data_PVOH_foam_rheology.csv` is a **full factorial**: three pipe radii
(hydrodynamic `R` = 4.075, 6.0, 9.5 mm from the `D hydrodynamic` column) [**ERROR — see addendum A.1**], three PVOH
concentrations (0.5, 1.0, 3.0 %), five foam density setpoints (100–300 g/l), ~15 flow
rates each — 660 rows. **Every one of the 15 (PVOH, density) cells exists at all three
radii.** Each row carries the volumetric flow rate, the pressure gradient, the sample
density, the Sauter mean bubble diameter, and — this is what makes the item possible —
the **slip velocity, measured independently by high-speed video imaging of bubbles at the
wall**, which the paper states to better than 3 %.

The paper's own declared measurement errors, which are the only scales here this repo did
not invent: **slip velocity < 3 %**, **flow rate ≤ 5 %** (§2, the splitter experiment).

`Set_1` is a different, unbalanced campaign (24 cells, three radii, four PVOH
concentrations) and is NOT used for the gates. It is loaded only for G27.7.

## 2. The reduction, stated before it is run

For each row, with `R` the hydrodynamic radius, `Q` the volumetric flow rate and
`dP/dL` the measured pressure gradient:

[**INCOMPLETE — see addendum A.2: `Q` needs a compressibility factor before this line**]

```
v_avg  = Q / (pi R^2)
v_bulk = v_avg - v_s                  (v_s = measured slip velocity; Eq. 8 of the paper)
gd_a   = 4 v_bulk / R                 (apparent shear rate,          Eq. 9)
tau_w  = R (dP/dL) / 2                (exact for stationary pipe flow, Eq. 7)
```

The Weissenberg–Rabinowitsch correction is applied **per cell** by fitting
`log tau_w` against `log gd_a` and taking the local slope `N = d log tau_w / d log gd_a`:

```
gd_w = gd_a (3 N + 1) / (4 N)
eta  = tau_w / gd_w
```

This is the paper's own reduction. Nothing here is a choice this repo made.

## 3. Gates

Every bar below is fixed now. A failed bar is recorded with its number and a **stricter**
test is added beside it; no bar is loosened.

### G27.1 — reproduce the source's own fit table before extending it

The paper publishes **Table 2**: a regression of Eq. (19),
`mu/mu_l = B(phi) (gd D mu_l / sigma_l)^m(phi)`, performed **at each air content `phi`
separately, pooling over PVOH concentration, pipe radius and bubble size**, with
`n = 1 + m` stated in the text (line 472). The published values, with the paper's own
standard errors:

| `phi` | density setpoint (g/l) | `m(phi)` | `n = 1 + m` |
|---|---|---|---|
| 0.70 | 300 | −0.639 ± 0.014 | 0.361 ± 0.014 |
| 0.75 | 250 | −0.637 ± 0.006 | 0.363 ± 0.006 |
| 0.80 | 200 | −0.691 ± 0.004 | 0.309 ± 0.004 |
| 0.85 | 150 | −0.735 ± 0.005 | 0.265 ± 0.005 |
| 0.90 | 100 | −0.762 ± 0.006 | 0.238 ± 0.006 |

(`phi = (rho_l - rho_f)/rho_l` with `rho_l = 1000 kg/m^3`, which maps the five setpoints
onto the five table rows exactly.)

**Gate.** Fit `power_law` through `biorheo.identify.fit_rate_curve` on the reduced
`(gd_w, eta)` of every cell, pool `n` across PVOH and radius **at fixed density setpoint**
— the same pooling the paper states — and compare the pooled median against `1 + m`.

**Bar: `|n_biorheo - n_paper| <= 3 x sigma_paper(m)` on at least 4 of the 5 air contents.**
Chosen as three of the SOURCE's own published standard errors, not a tolerance invented
here (the item-23 lesson). The bar is deliberately not tighter, for a reason stated in
advance: Eq. (19) collapses the bubble-size and liquid-viscosity dependence onto the same
exponent as the shear rate, an approximation the paper itself quantifies at "the average
difference was below 4 %" (line 456), while a per-cell fit varies only the shear rate.
**A per-cell `n` and a pooled `m` are therefore not the same estimator**, and the gate is
a consistency check between two estimators, not an identity.

If this fails, nothing downstream is trustworthy and the item stops here.

### G27.2 — THE HEADLINE. Identify big, predict small: the direction the business sells

A rheometer measures at low wall rates; a needle runs at high ones. The analogue here is
identify on the **largest** pipe (`R` = 9.5 mm) and predict the **smallest** (`R` = 4.075
mm), a factor of **2.33 in radius**.

For each of the 15 cells: fit `power_law` on the `R` = 9.5 mm rows only; then for every
`R` = 4.075 mm row of the same cell, take that row's **measured** `Q` and **measured**
`v_s`, form `Q_bulk = Q - v_s pi R^2`, and solve
`biorheo.apps.injectability.pressure_for_flow_rate` for the pressure drop that carries
`Q_bulk` through a tube of that radius. Report `(dP/dL)_pred` against `(dP/dL)_meas`.

**Bar: median over all predicted rows of `|rel err|` in `dP/dL` <= 20 %.**
Justification, stated in advance: the source declares 5 % on the flow rate, and
`dP/dL ~ Q^n` with `n ~ 0.3` propagates that to ~1.5 %; 20 % is deliberately **four times
looser than the propagated measurement error**, because this is a cross-geometry
extrapolation and not a repeat measurement. The full distribution (median, p90, worst,
per-cell) is reported whatever the verdict.

### G27.3 — the discriminator. Without it, G27.2 decides nothing

Repeat G27.2 **with the slip term deleted everywhere**: `v_bulk := v_avg` in the
identification and `Q_bulk := Q` in the prediction. This is the calculation
`injectability.py` performs today.

**Prediction, registered now: this FAILS**, because the paper states the slip carries most
of the volumetric flow at low flow rates (abstract, line 22).

**Bar: the separation ratio** `median|rel err|_noslip / median|rel err|_slip` **must
exceed 1.6** before G27.2's verdict may be called informative. That threshold is this
lab's standing rule for a discriminator (`lessons.md`, 2026-08-30): a test whose two
candidates predict the same thing has discriminated nothing, however it comes out.

### G27.4 — the positive control

Identify on `R` = 9.5 mm and predict **the same pipe's own rows**. This is not a
prediction; it is the fit's residual expressed in the same units as G27.2, and it bounds
how much of G27.2's error could ever be removed by a better geometry model.

**Bar: median `|rel err|` <= 5 %.** If the fit cannot reproduce its own data to that,
`power_law` does not describe these curves and G27.2's number is a misfit reading.

### G27.5 — the reverse direction

Identify on `R` = 4.075 mm, predict `R` = 9.5 mm. No bar: reported for the asymmetry.
Registered prediction: **the reverse direction is WORSE**, because the small pipe's rows
sit at higher wall rates and the extrapolation is then downward into the region where a
power law with `n ~ 0.3` diverges as `gd -> 0`.

### G27.6 — a slip channel that never touches the constitutive fit

Mooney's construction: at fixed wall shear stress, `gd_a,total = 4 Q/(pi R^3)` plotted
against `1/R` has slope `4 v_s` and intercept the true apparent rate. With three radii
this is a two-parameter fit at each `tau_w`, and it uses **no constitutive model at all**.

**Gate: the Mooney slip velocity against the video-imaged slip velocity**, interpolated to
the same `tau_w`. **Bar: median `|rel err|` <= 50 %.** Loose on purpose and said so: the
Mooney slope is a three-point regression over a 2.33x lever in `1/R`, so it is a weak
estimator, and the point of the gate is that it is **independent of everything else in
this item** — it is the anchor that does not pass through the code under test.

### G27.7 — is `Set_1` a repeat of `Set_2`?

The deposit ships two sets and says nothing about their relationship. 12 (radius, PVOH,
density) combinations occur in both. Interpolate each set's `dP/dL` onto the other's flow
rates within the overlapping range and report the median absolute relative difference.

**No bar.** This is a measurement of the deposit, registered because if the two sets ARE
repeats the number is a repeatability that could scale every gate above, and if they are
not, quoting it as one would be an invention. **Whichever it is, it is reported and no
gate above is re-graded on it**, because the gates' bars are fixed before this is known.

## 4. What this item may NOT claim

- **Not a biofluid.** PVOH foam is a foam. What is being anchored is the **chain** —
  reduce a tube measurement, identify a model with an interval, predict a different tube —
  which is fluid-agnostic. No number here may be quoted about blood, a mAb formulation, or
  a microgel.
- **Not an injectability number.** The wall rates here are pipe-flow rates, orders below a
  27G needle. The standing rule about extrapolation flags applies unchanged.
- **Not a claim that biorheo models slip.** It does not. This item measures the size of the
  error the omission causes in one dataset, at one set of geometries.
- **No record enters the registry from this item** unless a separate decision says so.
  The reduced curve is derived through a correction (WR + slip subtraction), not measured,
  and the registry has never held a derived curve.

## 5. Registered predictions, so they can be scored

1. G27.1 passes on 4 or 5 of 5. — *scored in §6 of the report*
2. G27.2 passes.
3. G27.3 fails, with a separation ratio above 3.
4. G27.4 passes.
5. G27.5 is worse than G27.2.
6. G27.6 passes.
7. `Set_1` and `Set_2` are NOT repeats (different campaigns), so G27.7's number is large.

---

# Addendum A (2026-09-08, after the run)

A pre-registration is not rewritten. These are the errors found in the text above and the
gates added after a pre-registered one failed to decide.

## A.1 — the middle pipe's radius is wrong in §1

§1 says the three hydrodynamic radii are 4.075, **6.0** and 9.5 mm. The `D hydrodynamic`
column gives 8.15, **11.90**, 19.00 mm, so the middle radius is **5.95 mm**, not 6.0. The
paper's own methods section says "6.0 mm", which is where the slip came from. It affects
only G27.6, which is the one gate that uses all three pipes; G27.2-G27.5 use the outer two,
whose values were right. The script reads the column and never used the wrong number.

## A.2 — §2's reduction is INCOMPLETE, and the missing term is large

§2 writes `v_avg = Q / (pi R^2)` directly from the deposit's flow-rate column. That is
wrong, and the deposit's own column descriptions say so: the flow rate is measured **at the
absolute pressure of column 8** (115-408 kPa, at the foam-generator outlet), while the
pressure gradient and the sample density are at **~1 atm**, in the measurement span. Foam is
compressible, so with `phi` the air content at 1 atm,

    Q_local = Q * 1 / [ (1 - phi) + phi * P_atm / P_abs ] .

**Measured: the factor runs 1.091-3.120, median 1.574.** Omitting it understates the flow in
the measurement pipe by about 57 %.

**It was not found by argument.** The check that settles it can fail: a flow with a plug
cannot have a wall slip velocity above its own mean velocity. On the raw column that bound
is violated on **11 of 660 rows** (worst 1.403); with the factor applied it is violated on
**0**, and the largest slip fraction in the deposit is 0.916. The correction was applied
before any gate was graded, and the bars in §3 are unchanged.

Not corrected and named rather than modelled: the measurement span is not isobaric (0.3 m
at up to ~60 kPa/m is ~18 % of an atmosphere in the smallest pipe).

## A.3 — three gates added after the fact, each STRICTER than the one it sits beside

None replaces a bar. Every original verdict and number stays on the report.

- **G27.1'** (after G27.1 FAILED, 2 of 5). G27.1 compares one fitted exponent against one
  published exponent, and `B` and `m` trade off along the data exactly as `k` and `n` do.
  G27.1' evaluates Eq. (19) with BOTH published parameters at every reduced point and
  compares as viscosity curves over the measured range. No bar is invented for it; it is
  reported as a size (median 8.1 %).
- **G27.3'** (after G27.3 REFUTED its own prediction: deleting slip made the prediction
  BETTER, separation 0.56 against a 1.6 bar). It does not compare two predictions; it
  measures whether they could ever have differed. `4 v_s/R = c * gd_a` with `c = v_s/v_avg`,
  so the slip correction is a multiplication by `1 - c`, and a power law rescaled by a
  constant in `gd` has the SAME exponent. If `c` does not depend on diameter the rescaling
  cancels between identification and prediction and slip is degenerate with the fit.
  Measured: `c` is 0.361 / 0.267 / 0.331 across the three pipes at matched `tau_w`, **not
  even monotone**, and the surviving part `(1-c)_big/(1-c)_small` is 1.068 -- 6.8 % against
  an 18.7 % prediction error.
- **G27.6'** (after G27.6 FAILED at 140 % against a 50 % bar). An unsigned error hides both
  the sign and the admissibility. Mooney runs HIGH by a median factor of 2.20, and **25 of
  144 points return a NEGATIVE slip velocity** -- not a large error but an inadmissible
  answer, from the same non-monotonicity G27.3' measures.

## A.4 — scoring §5

Five of the seven registered predictions are REFUTED: 1, 3, 5, 6 and 7. Predictions 2 and 4
are confirmed. Prediction 3's refutation is the consequential one, because it removes the
ground for reading G27.2's pass as evidence about wall slip.
