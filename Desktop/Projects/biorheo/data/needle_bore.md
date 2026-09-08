# Item 18 — the needle-bore exponent, and what a force budget costs once the fit has an error bar

Written by `scripts/needle_bore.py` from `data/needle_bore_prereg.md`. biorheo 0.0.1. Every number on this page is printed by that script.

## What is measured

```
b  =  d ln F / d ln R        at fixed flow rate, needle length and syringe
```

A Newtonian sizing calculation assumes `b = -4`. For a shear-thinning fluid it is not: a narrower bore raises the wall shear rate, which thins the fluid, so the force rises less steeply. **`b` is a property of the fluid AND the operating point**, not a material constant — a fluid already on its infinite-shear plateau at the wall returns `-4` however shear-thinning it is at rheometer rates, and that is exactly what the blood records do.

Geometry is DECLARED, not standardised. The bore axis is the ratio `rho = R/R_ref` with `R_ref = 0.000105 m`, the same declared value the shipped blood study runs at and carrying that study's `[UNVERIFIED]` flag. **No needle gauge is named anywhere on this page**; this repo retracted an ISO 9626 comparison once and the standard is still not in hand.

## Gates

| gate | statement | result |
|---|---|---|
| G18.1 | positive control, `b = -4` on the 20 `newtonian` records | **PASS**, 20/20, worst |b+4| = 4.66e-10 |
| G18.2 | closed form `b = -(3n+1)` on the 20 `power_law` records, computed without the force solver | **PASS**, 20/20, worst relative error 2.66e-10 |
| G18.3 | the discriminating half — the routine does not return `-4` for everything | **PASS**, 36 records with |b+4| > 0.1, worst 1.867 (`carbopol_etd2020_0p3wt_50uL`) |
| G18.4 | no record is steeper than Newtonian | **PASS**, 115/115, steepest -4 |
| G18.5 | SPLIT prediction on the Newtonian rule's over-prediction at `rho = 0.7` | **PASS** (see below) |
| G18.6 | the derating factor is below 1 on every closed-interval record | **PASS**, 91/91, strongest derating 0.123 |
| **G18.6′** | the stricter test beside it: zero-width intervals must give a factor of exactly 1, real ones strictly less | **PASS**, 1 and 0.8476 |
| G18.7 | three error terms ranked per record | reported, not barred |

**G18.6 as pre-registered is close to unfalsifiable, and that is recorded rather than fixed by moving it.** The upper force envelope is a maximum over profile ends and the nominal force is not one of them, so the factor is below 1 unless the fitted parameters locally maximise the force in every direction at once. A routine that returned a constant 0.5 would have passed it. The bar and its number stand; **G18.6′** is the discriminating test added beside it, and it pins both ends the bar cannot.

## G18.5, the split prediction, and it came out as pre-registered

The pre-registration predicted the blood half would **FAIL** its 10 % threshold and the Carbopol half would **PASS** it, for a stated reason in each case.

- **Blood: 20 of 20** Casson donors are BELOW 10 %; the worst is 1.24 %. The Newtonian rule is very nearly right for blood in this needle, because at the wall the fluid is on its infinite-shear plateau — the same fact item 14 measured a different way (`|d ln F / d ln eta_c| = 0.995`).
- **Carbopol: 2 of 2** named gels are ABOVE 10 %, and 15 of 15 of all Carbopol records are. The over-prediction runs 47.9 % to 93.9 %.

**The commercial reading, and it is the point of the item.** Whether a Newtonian sizing rule is safe is a property of *your* fluid at *your* operating point, and the two fluids in this registry sit on opposite sides of the question by an order of magnitude. A shop that quotes `R^-4` without measuring `b` cannot tell which case it is in — and the measurement costs one extra evaluation of a curve that has already been fitted.

## The three error terms, in the same units

All three are relative force errors at `Q1`, so they can be compared directly:

- **bore tolerance** — `|b| x 5 %`, a DECLARED bore error, not a standard's tolerance. It scales linearly with whatever error the caller declares.
- **parameter uncertainty** — the one-at-a-time profile push-through hull, which is *wider* than the joint 95 % region it is drawn from.
- **family choice** — item 12's Cross-vs-Casson difference, READ from that study's artifact and not recomputed. Defined only for the Casson blood records.

| winner | records |
|---|---|
| parameter_uncertainty | 87 of 115 |
| bore_tolerance | 28 of 115 |

On the 20 records where all three terms are defined (the Casson blood donors), the ranking is **parameter_uncertainty** 17, **bore_tolerance** 3.

## Derating against a force budget

`Q_max` is the flow rate at which the plunger force reaches **25 N**. That figure is a **rule of thumb** quoted in the JHU thesis, not a standard, and the absolute `Q_max` inherits every omission in `injectability.py` (no barrel pressure drop, no entrance loss, no wall slip, no hub friction, no stiction). **The ratio is the usable number**, because those omissions largely cancel between the two crossings.

Over 91 records with fully closed profile intervals, the derating factor `Q_max(95 % upper) / Q_max(fit)` runs from 0.12304 to 0.9718. **The parameter uncertainty alone costs up to a factor of 8.13 in deliverable flow rate**, before any of the physics the force model leaves out.

**The split by identifiability is a TAIL, not a shift, and the medians do not show it.**
Records whose every parameter is `measurable` derate by a median 0.7959 (91 records); records carrying at least one non-`measurable` verdict derate by a median 0.8664 (24 records) — **the same, to within the spread**. The difference is entirely at the bad end: the worst derating among fully-measurable records is 0.123, a factor of 8.13, against 7.169e-05 — a factor of 1.39e+04 — among those with a non-measurable verdict. So an identifiability verdict does not predict how wide a typical error bar will be; what it predicts is whether the error bar can be catastrophic. The point forecast looks identical either way.

## Per-record table

| record | model | `b` at Q1 | `b` at Q2 | closed form | Newtonian-rule over-prediction at rho = 0.7 | hull | derating | all measurable? | dominant term |
|---|---|---|---|---|---|---|---|---|---|
| `blood_horner2020_donor_A` | casson | -3.981 | -3.9922 | - | 0.526 % | 33.4 % | 0.8476 | yes | parameter_uncertainty |
| `blood_horner2020_donor_A_cross` | cross | -3.9201 | -3.9639 | - | 2.3 % | 193 % | 0.82087 | NO | parameter_uncertainty |
| `blood_horner2020_donor_A_carreau_yasuda` | carreau_yasuda | -3.8939 | -3.9471 | - | 3.15 % | 1.5e+04 % | 0.0055245 | NO | parameter_uncertainty |
| `blood_horner2020_donor_A_newtonian` | newtonian | -4 | -4 | -4 | 2.03e-10 % | 80.6 % | 0.67506 | yes | parameter_uncertainty |
| `blood_horner2020_donor_A_power_law` | power_law | -3.2077 | -3.2077 | -3.2077 | 32.7 % | 45.1 % | 0.5875 | yes | parameter_uncertainty |
| `blood_horner2020_donor_B` | casson | -3.9709 | -3.988 | - | 0.809 % | 15.1 % | 0.92694 | yes | bore_tolerance |
| `blood_horner2020_donor_B_cross` | cross | -3.978 | -3.9927 | - | 0.578 % | 358 % | 0.80629 | NO | parameter_uncertainty |
| `blood_horner2020_donor_B_carreau_yasuda` | carreau_yasuda | -3.978 | -3.9926 | - | 0.578 % | 6.32e+05 % | 0.00015635 | NO | parameter_uncertainty |
| `blood_horner2020_donor_B_newtonian` | newtonian | -4 | -4 | -4 | -3.35e-12 % | 104 % | 0.60752 | yes | parameter_uncertainty |
| `blood_horner2020_donor_B_power_law` | power_law | -3.0207 | -3.0207 | -3.0207 | 41.8 % | 105 % | 0.26415 | yes | parameter_uncertainty |
| `blood_horner2020_donor_C` | casson | -3.9587 | -3.983 | - | 1.15 % | 14.2 % | 0.93105 | yes | bore_tolerance |
| `blood_horner2020_donor_C_cross` | cross | -3.9831 | -3.9951 | - | 0.426 % | 156 % | 0.92653 | NO | parameter_uncertainty |
| `blood_horner2020_donor_C_carreau_yasuda` | carreau_yasuda | -3.9846 | -3.9956 | - | 0.387 % | 2.08e+03 % | 0.39921 | NO | parameter_uncertainty |
| `blood_horner2020_donor_C_newtonian` | newtonian | -4 | -4 | -4 | -3.19e-12 % | 126 % | 0.55179 | yes | parameter_uncertainty |
| `blood_horner2020_donor_C_power_law` | power_law | -2.8301 | -2.8301 | -2.8301 | 51.8 % | 125 % | 0.1607 | yes | parameter_uncertainty |
| `blood_horner2020_donor_D` | casson | -3.9556 | -3.9817 | - | 1.24 % | 18.6 % | 0.9107 | yes | bore_tolerance |
| `blood_horner2020_donor_D_cross` | cross | -3.9935 | -3.9985 | - | 0.154 % | 483 % | 0.88546 | NO | parameter_uncertainty |
| `blood_horner2020_donor_D_carreau_yasuda` | carreau_yasuda | -3.9912 | -3.9978 | - | 0.215 % | 14.7 % | 0.93063 | NO | bore_tolerance |
| `blood_horner2020_donor_D_newtonian` | newtonian | -4 | -4 | -4 | -5.33e-13 % | 132 % | 0.53898 | yes | parameter_uncertainty |
| `blood_horner2020_donor_D_power_law` | power_law | -2.7915 | -2.7915 | -2.7915 | 53.9 % | 142 % | 0.12304 | yes | parameter_uncertainty |
| `blood_horner2020_donor_F` | casson | -3.9772 | -3.9906 | - | 0.633 % | 25.3 % | 0.8816 | yes | parameter_uncertainty |
| `blood_horner2020_donor_F_cross` | cross | -3.9667 | -3.9877 | - | 0.902 % | 221 % | 0.87063 | NO | parameter_uncertainty |
| `blood_horner2020_donor_F_carreau_yasuda` | carreau_yasuda | -3.9619 | -3.9853 | - | 1.04 % | 2.06e+04 % | 0.0028439 | NO | parameter_uncertainty |
| `blood_horner2020_donor_F_newtonian` | newtonian | -4 | -4 | -4 | -3.21e-12 % | 89.7 % | 0.64735 | yes | parameter_uncertainty |
| `blood_horner2020_donor_F_power_law` | power_law | -3.1347 | -3.1347 | -3.1347 | 36.2 % | 73.3 % | 0.40598 | yes | parameter_uncertainty |
| `blood_horner2020_donor_G` | casson | -3.9862 | -3.9943 | - | 0.383 % | 46 % | 0.79772 | yes | parameter_uncertainty |
| `blood_horner2020_donor_G_cross` | cross | -3.9874 | -3.9967 | - | 0.309 % | 32.3 % | 0.87885 | yes | parameter_uncertainty |
| `blood_horner2020_donor_G_carreau_yasuda` | carreau_yasuda | -3.9241 | -3.9643 | - | 2.21 % | 28.2 % | 0.87624 | NO | parameter_uncertainty |
| `blood_horner2020_donor_G_newtonian` | newtonian | -4 | -4 | -4 | -2.81e-12 % | 68.4 % | 0.71474 | yes | parameter_uncertainty |
| `blood_horner2020_donor_G_power_law` | power_law | -3.3261 | -3.3261 | -3.3261 | 27.2 % | 45.7 % | 0.60999 | yes | parameter_uncertainty |
| `blood_horner2020_donor_H` | casson | -3.9715 | -3.9883 | - | 0.791 % | 43.1 % | 0.80829 | yes | parameter_uncertainty |
| `blood_horner2020_donor_H_cross` | cross | -3.9902 | -3.9977 | - | 0.235 % | 16.5 % | 0.92632 | yes | bore_tolerance |
| `blood_horner2020_donor_H_carreau_yasuda` | carreau_yasuda | -3.9767 | -3.9927 | - | 0.602 % | 16.9 % | 0.9218 | yes | bore_tolerance |
| `blood_horner2020_donor_H_newtonian` | newtonian | -4 | -4 | -4 | -3.02e-12 % | 102 % | 0.61232 | yes | parameter_uncertainty |
| `blood_horner2020_donor_H_power_law` | power_law | -3.0196 | -3.0196 | -3.0196 | 41.9 % | 75.2 % | 0.37301 | yes | parameter_uncertainty |
| `blood_horner2020_donor_I` | casson | -3.9722 | -3.9886 | - | 0.773 % | 39 % | 0.82453 | yes | parameter_uncertainty |
| `blood_horner2020_donor_I_cross` | cross | -3.9898 | -3.9974 | - | 0.249 % | 20 % | 0.9134 | yes | parameter_uncertainty |
| `blood_horner2020_donor_I_carreau_yasuda` | carreau_yasuda | -3.9738 | -3.9913 | - | 0.687 % | 11.9 % | 0.94364 | NO | bore_tolerance |
| `blood_horner2020_donor_I_newtonian` | newtonian | -4 | -4 | -4 | -3.56e-12 % | 101 % | 0.6159 | yes | parameter_uncertainty |
| `blood_horner2020_donor_I_power_law` | power_law | -3.0335 | -3.0335 | -3.0335 | 41.2 % | 78.4 % | 0.3618 | yes | parameter_uncertainty |
| `blood_horner2020_donor_J` | casson | -3.9718 | -3.9884 | - | 0.783 % | 28.8 % | 0.8666 | yes | parameter_uncertainty |
| `blood_horner2020_donor_J_cross` | cross | -3.9943 | -3.9987 | - | 0.136 % | 18.7 % | 0.9183 | yes | bore_tolerance |
| `blood_horner2020_donor_J_carreau_yasuda` | carreau_yasuda | -3.9894 | -3.997 | - | 0.265 % | 20.8 % | 0.90674 | NO | parameter_uncertainty |
| `blood_horner2020_donor_J_newtonian` | newtonian | -4 | -4 | -4 | 3.71e-12 % | 102 % | 0.61336 | yes | parameter_uncertainty |
| `blood_horner2020_donor_J_power_law` | power_law | -3.0362 | -3.0362 | -3.0362 | 41 % | 97.4 % | 0.29738 | yes | parameter_uncertainty |
| `blood_horner2020_donor_K` | casson | -3.9744 | -3.9895 | - | 0.712 % | 31.6 % | 0.85492 | yes | parameter_uncertainty |
| `blood_horner2020_donor_K_cross` | cross | -3.9311 | -3.9712 | - | 1.94 % | 165 % | 0.8671 | NO | parameter_uncertainty |
| `blood_horner2020_donor_K_carreau_yasuda` | carreau_yasuda | -3.9439 | -3.9791 | - | 1.53 % | 1.43e+03 % | 0.26589 | NO | parameter_uncertainty |
| `blood_horner2020_donor_K_newtonian` | newtonian | -4 | -4 | -4 | -3.36e-12 % | 96.1 % | 0.62908 | yes | parameter_uncertainty |
| `blood_horner2020_donor_K_power_law` | power_law | -3.0696 | -3.0696 | -3.0696 | 39.4 % | 62.9 % | 0.4445 | yes | parameter_uncertainty |
| `blood_horner2020_donor_L` | casson | -3.9587 | -3.983 | - | 1.15 % | 25.8 % | 0.87893 | yes | parameter_uncertainty |
| `blood_horner2020_donor_L_cross` | cross | -3.9579 | -3.9857 | - | 1.11 % | 727 % | 0.62283 | NO | parameter_uncertainty |
| `blood_horner2020_donor_L_carreau_yasuda` | carreau_yasuda | -3.9579 | -3.9857 | - | 1.12 % | 6.66e+04 % | 7.1691e-05 | NO | parameter_uncertainty |
| `blood_horner2020_donor_L_newtonian` | newtonian | -4 | -4 | -4 | -6.5e-11 % | 125 % | 0.55442 | yes | parameter_uncertainty |
| `blood_horner2020_donor_L_power_law` | power_law | -2.8234 | -2.8234 | -2.8234 | 52.1 % | 95.9 % | 0.23626 | yes | parameter_uncertainty |
| `blood_horner2020_donor_M` | casson | -3.9574 | -3.9824 | - | 1.19 % | 43.5 % | 0.80648 | yes | parameter_uncertainty |
| `blood_horner2020_donor_M_cross` | cross | -3.9905 | -3.9978 | - | 0.226 % | 33.8 % | 0.86122 | yes | parameter_uncertainty |
| `blood_horner2020_donor_M_carreau_yasuda` | carreau_yasuda | -3.9758 | -3.9929 | - | 0.616 % | 19.1 % | 0.91261 | NO | bore_tolerance |
| `blood_horner2020_donor_M_newtonian` | newtonian | -4 | -4 | -4 | 5.54e-11 % | 127 % | 0.54958 | yes | parameter_uncertainty |
| `blood_horner2020_donor_M_power_law` | power_law | -2.8092 | -2.8092 | -2.8092 | 52.9 % | 104 % | 0.20761 | yes | parameter_uncertainty |
| `blood_horner2020_donor_N` | casson | -3.9643 | -3.9853 | - | 0.993 % | 38.7 % | 0.82557 | yes | parameter_uncertainty |
| `blood_horner2020_donor_N_cross` | cross | -3.9378 | -3.9774 | - | 1.68 % | 11.3 % | 0.9461 | yes | bore_tolerance |
| `blood_horner2020_donor_N_carreau_yasuda` | carreau_yasuda | -3.9266 | -3.9718 | - | 2.02 % | 15.5 % | 0.94165 | yes | bore_tolerance |
| `blood_horner2020_donor_N_newtonian` | newtonian | -4 | -4 | -4 | -4.45e-12 % | 115 % | 0.57829 | yes | parameter_uncertainty |
| `blood_horner2020_donor_N_power_law` | power_law | -2.8963 | -2.8963 | -2.8963 | 48.2 % | 64.9 % | 0.40125 | yes | parameter_uncertainty |
| `blood_horner2020_donor_O` | casson | -3.9732 | -3.989 | - | 0.744 % | 28.8 % | 0.86653 | yes | parameter_uncertainty |
| `blood_horner2020_donor_O_cross` | cross | -3.9754 | -3.992 | - | 0.642 % | 14.5 % | 0.93503 | yes | bore_tolerance |
| `blood_horner2020_donor_O_carreau_yasuda` | carreau_yasuda | -3.9724 | -3.9907 | - | 0.729 % | 257 % | 0.86572 | NO | parameter_uncertainty |
| `blood_horner2020_donor_O_newtonian` | newtonian | -4 | -4 | -4 | -3.38e-12 % | 98.4 % | 0.62237 | yes | parameter_uncertainty |
| `blood_horner2020_donor_O_power_law` | power_law | -3.0547 | -3.0547 | -3.0547 | 40.1 % | 77.4 % | 0.36934 | yes | parameter_uncertainty |
| `blood_horner2020_donor_P` | casson | -3.9715 | -3.9883 | - | 0.791 % | 25 % | 0.8828 | yes | parameter_uncertainty |
| `blood_horner2020_donor_P_cross` | cross | -3.9606 | -3.9855 | - | 1.06 % | 164 % | 0.88947 | NO | parameter_uncertainty |
| `blood_horner2020_donor_P_carreau_yasuda` | carreau_yasuda | -3.9594 | -3.9849 | - | 1.1 % | 1.48e+04 % | 0.004004 | NO | parameter_uncertainty |
| `blood_horner2020_donor_P_newtonian` | newtonian | -4 | -4 | -4 | -3.38e-12 % | 102 % | 0.61253 | yes | parameter_uncertainty |
| `blood_horner2020_donor_P_power_law` | power_law | -3.0242 | -3.0242 | -3.0242 | 41.6 % | 82.4 % | 0.34365 | yes | parameter_uncertainty |
| `blood_horner2020_donor_Q` | casson | -3.9691 | -3.9873 | - | 0.86 % | 34.6 % | 0.84234 | yes | parameter_uncertainty |
| `blood_horner2020_donor_Q_cross` | cross | -3.9778 | -3.9934 | - | 0.567 % | 12.5 % | 0.94243 | yes | bore_tolerance |
| `blood_horner2020_donor_Q_carreau_yasuda` | carreau_yasuda | -3.9637 | -3.9874 | - | 0.966 % | 14.6 % | 0.92549 | yes | bore_tolerance |
| `blood_horner2020_donor_Q_newtonian` | newtonian | -4 | -4 | -4 | -3.38e-12 % | 107 % | 0.60002 | yes | parameter_uncertainty |
| `blood_horner2020_donor_Q_power_law` | power_law | -2.9789 | -2.9789 | -2.9789 | 43.9 % | 78.1 % | 0.34544 | yes | parameter_uncertainty |
| `blood_horner2020_donor_R` | casson | -3.966 | -3.986 | - | 0.947 % | 31.5 % | 0.85483 | yes | parameter_uncertainty |
| `blood_horner2020_donor_R_cross` | cross | -3.9858 | -3.9962 | - | 0.351 % | 13 % | 0.94043 | yes | bore_tolerance |
| `blood_horner2020_donor_R_carreau_yasuda` | carreau_yasuda | -3.9785 | -3.9934 | - | 0.552 % | 13.8 % | 0.93397 | yes | bore_tolerance |
| `blood_horner2020_donor_R_newtonian` | newtonian | -4 | -4 | -4 | -2.3e-12 % | 112 % | 0.58487 | yes | parameter_uncertainty |
| `blood_horner2020_donor_R_power_law` | power_law | -2.9332 | -2.9332 | -2.9332 | 46.3 % | 92.3 % | 0.27829 | yes | parameter_uncertainty |
| `blood_horner2020_donor_S` | casson | -3.9657 | -3.9859 | - | 0.955 % | 29.3 % | 0.86415 | yes | parameter_uncertainty |
| `blood_horner2020_donor_S_cross` | cross | -3.9838 | -3.9955 | - | 0.405 % | 10.9 % | 0.9494 | yes | bore_tolerance |
| `blood_horner2020_donor_S_carreau_yasuda` | carreau_yasuda | -3.9778 | -3.9932 | - | 0.57 % | 12.4 % | 0.93901 | yes | bore_tolerance |
| `blood_horner2020_donor_S_newtonian` | newtonian | -4 | -4 | -4 | 2.86e-12 % | 113 % | 0.58341 | yes | parameter_uncertainty |
| `blood_horner2020_donor_S_power_law` | power_law | -2.9288 | -2.9288 | -2.9288 | 46.5 % | 93.9 % | 0.27357 | yes | parameter_uncertainty |
| `blood_horner2020_donor_T` | casson | -3.9686 | -3.9871 | - | 0.874 % | 26.7 % | 0.87548 | yes | parameter_uncertainty |
| `blood_horner2020_donor_T_cross` | cross | -3.9752 | -3.9921 | - | 0.644 % | 14.6 % | 0.9341 | yes | bore_tolerance |
| `blood_horner2020_donor_T_carreau_yasuda` | carreau_yasuda | -3.9671 | -3.9886 | - | 0.874 % | 20.5 % | 0.93079 | NO | parameter_uncertainty |
| `blood_horner2020_donor_T_newtonian` | newtonian | -4 | -4 | -4 | -3.39e-12 % | 108 % | 0.59712 | yes | parameter_uncertainty |
| `blood_horner2020_donor_T_power_law` | power_law | -2.9743 | -2.9743 | -2.9743 | 44.2 % | 88.3 % | 0.30351 | yes | parameter_uncertainty |
| `blood_horner2020_donor_U` | casson | -3.9692 | -3.9874 | - | 0.856 % | 26.4 % | 0.8766 | yes | parameter_uncertainty |
| `blood_horner2020_donor_U_cross` | cross | -3.9514 | -3.9817 | - | 1.32 % | 73.1 % | 0.9471 | NO | parameter_uncertainty |
| `blood_horner2020_donor_U_carreau_yasuda` | carreau_yasuda | -3.9506 | -3.9813 | - | 1.35 % | 3.63e+04 % | 0.00055105 | NO | parameter_uncertainty |
| `blood_horner2020_donor_U_newtonian` | newtonian | -4 | -4 | -4 | -8.99e-13 % | 107 % | 0.60042 | yes | parameter_uncertainty |
| `blood_horner2020_donor_U_power_law` | power_law | -2.9814 | -2.9814 | -2.9814 | 43.8 % | 80.2 % | 0.34045 | yes | parameter_uncertainty |
| `carbopol_etd2020_0p1wt_30uL` | herschel_bulkley | -2.7334 | -2.7357 | - | 57.1 % | 75.7 % | 0.42533 | yes | parameter_uncertainty |
| `carbopol_etd2020_0p1wt_40uL` | herschel_bulkley | -2.9033 | -2.9038 | - | 47.9 % | 68.1 % | 0.50258 | yes | parameter_uncertainty |
| `carbopol_etd2020_0p1wt_50uL` | herschel_bulkley | -2.7757 | -2.7763 | - | 54.7 % | 42.7 % | 0.61936 | yes | parameter_uncertainty |
| `carbopol_etd2020_0p3wt_30uL` | herschel_bulkley | -2.223 | -2.2483 | - | 87.9 % | 13 % | 0.7758 | yes | parameter_uncertainty |
| `carbopol_etd2020_0p3wt_40uL` | herschel_bulkley | -2.1751 | -2.2044 | - | 91 % | 9.72 % | 0.82696 | yes | bore_tolerance |
| `carbopol_etd2020_0p3wt_50uL` | herschel_bulkley | -2.1327 | -2.1625 | - | 93.9 % | 11.3 % | 0.79586 | yes | parameter_uncertainty |
| `carbopol_etd2020_0p5wt_30uL` | herschel_bulkley | -2.4912 | -2.5055 | - | 71 % | 1.87 % | 0.9718 | yes | bore_tolerance |
| `carbopol_etd2020_0p5wt_40uL` | herschel_bulkley | -2.3127 | -2.3373 | - | 82 % | 2.68 % | 0.95777 | yes | bore_tolerance |
| `carbopol_etd2020_0p5wt_50uL` | herschel_bulkley | -2.2252 | -2.2563 | - | 87.6 % | 4.3 % | 0.92613 | yes | bore_tolerance |
| `carbopol_etd2020_0p7wt_30uL` | herschel_bulkley | -2.5294 | -2.5422 | - | 68.7 % | 5.48 % | 0.92203 | yes | bore_tolerance |
| `carbopol_etd2020_0p7wt_40uL` | herschel_bulkley | -2.4203 | -2.44 | - | 75.2 % | 4.52 % | 0.934 | yes | bore_tolerance |
| `carbopol_etd2020_0p7wt_50uL` | herschel_bulkley | -2.3617 | -2.3865 | - | 78.8 % | 4.57 % | 0.93381 | yes | bore_tolerance |
| `carbopol_etd2020_0p9wt_30uL` | herschel_bulkley | -2.5706 | -2.5859 | - | 66.1 % | 16.7 % | 0.79375 | yes | parameter_uncertainty |
| `carbopol_etd2020_0p9wt_40uL` | herschel_bulkley | -2.4219 | -2.4425 | - | 75.1 % | 8 % | 0.88266 | yes | bore_tolerance |
| `carbopol_etd2020_0p9wt_50uL` | herschel_bulkley | -2.3848 | -2.4087 | - | 77.3 % | 3.01 % | 0.95484 | yes | bore_tolerance |

## What this page does NOT say

- Nothing about any real needle gauge, any standard, or any real device.
- Nothing about a monoclonal antibody or any protein formulation. **The registry contains none.** The wedge fluid for the venture's first product is not in this data, and this page demonstrates the method on human blood and a microgel.
- Nothing about temperature; no record carries a temperature axis.
- **EVERY record on this page is `commercial_use: false` — all 115 of 115 rows cannot appear in a commercial deliverable.** The 15 Carbopol records are CC BY-NC-SA; the 100 blood records come from an author-posted workbook that states no licence at all, which is not permission. **This is a standing constraint on the venture, not a detail of this item**: the method demonstrated here is licensable, and not one curve currently in the registry is. A paid study needs either the client's own data or a curve whose licence permits commercial use.
- No record was added, no curve was fitted, and `dst_fit` was not touched.
