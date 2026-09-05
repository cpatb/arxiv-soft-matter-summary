# Horner 2020 human blood rheology (via Farrington, Jariwala, Armstrong 2023)

Downloaded 2026-09-05. Not a licensed dataset: the workbook is author-posted with a
citation request and no licence statement, so every record built from it carries
`commercial_use: false` and is for validation and demonstration only.

| file | origin | md5 |
|---|---|---|
| `Horner-Thesis-Data.xlsx` (16 MB, gitignored) | https://sites.udel.edu/wagnergroup/files/2016/06/Horner-Thesis-Data.xlsx | `9142d5ac31e70da5e63aeb97d961e407` |
| `HornerData_CassonFitTable.xlsx` (gitignored) | https://raw.githubusercontent.com/FarringtonSM/Hemostatistics_BloodRheology_and_MachineLearning/main/DATA/HornerData_CassonFitTable.xlsx | `92003f2d9fed31a142300aeed2280093` |
| `DonorOutliers.xlsx` (gitignored) | same repository, `DATA/DonorOutliers.xlsx` | `f2a3b0506aa9ee41014aa120e50e5265` |
| `steady_shear_human.csv` | the `Steady Shear` sheet of the workbook, human donors only, stress in Pa | `0d53329ee629efa583655f26f88b89a1` |
| `Horner_udel_0060D_14251.pdf` (8.5 MB, gitignored; `thesis_text.txt` is its pypdf extract, also gitignored) | https://udspace.udel.edu/server/api/core/bitstreams/e1aa1af1-569c-4909-a422-022eda350f56/content (DSpace item b79dda1f-cad2-43ff-a995-9265b9f86ab8, open access, downloaded 2026-09-05) | `27efa964c26779821aa68b99d5cd357f` |
| `physiology_and_casson_farrington2023.csv` | the Casson fit table: per-donor Casson yield stress and viscosity (mPa, mPa s, with the authors' standard errors), hematocrit, fibrinogen, sex | `f1ef42aa0ce1525641e372d6ae378d7d` |

Paper: Farrington, Jariwala, Armstrong, "Hemostatistics: blood rheology and machine
learning", Rheologica Acta (2023), doi 10.1007/s00397-023-01402-2 (resolved; in
`doi_cache.json`). Thesis: Horner, J. S., "An Experimental and Theoretical
Investigation of Blood Rheology", University of Delaware dissertation (2020); the
workbook's ReadMe asks that this be cited by any publication using the data.

## What the workbook's ReadMe says (verbatim points)

- ARES-G2 rheometer, double wall concentric cylinder geometry.
- Human measurements at 37 C; animal blood at the species' body temperature.
- All samples anticoagulated with EDTA and measured within ~4 hours of withdrawal.
- Sheets: Steady Shear; Step from 300 s-1 to 1 s-1; Amplitude Sweep (2 Hz); Frequency
  Sweep (strain amplitude 10); UD-LAOS and LAOS at (strain amplitude, rad/s) pairs.

The steady-shear sheet builds the records. The step-shear sheet (300 -> 1 1/s, 1 s
sampling, 40-50 s) is extracted to `step_300_to_1_human.csv`
(md5 `956171e5c1a772fdcbed512b693e5ccc`) and read by `scripts/step_shear_horner.py`;
the LAOS sheets are untouched.

## The measurement protocol, from the thesis itself (verified 2026-09-05)

The workbook's ReadMe does not state the protocol; the thesis does. Verbatim, from
Horner (2020) section 3.3.1 ("Experimental Tests and Protocol"):

> "The first test conducted on the sample upon loading into the rheometer and
> temperature equilibration was a steady shear curve at shear rates ranging from
> 0.1 s-1 to 700 s-1. For each shear rate, the sample was measured for 40 s or a minimum
> strain of 50. For shear rates above 2 s-1, the steady shear stress value was taken as
> the average of the final 15 transient stress measurements obtained on the sample. As
> for shear rates at or below 2 s-1, the steady shear stress response was taken as the
> maximum of the transient data. This method follows previously established guidelines
> for performing steady shear measurements on human blood and is implemented to account
> for the effects of syneresis, which can result in significantly lower measured
> stresses particularly at low shear rates."

Section 2.2 says the same for the Chapter 2 donors ("For measurements below 2 s-1 on
these steady state flow curves, the points were obtained by taking the maximum of the
time dependent data sampled over an interval of 1 s at the constant shear rate") and
adds the preshear: "Between tests, the sample was sheared at 300 s-1 for 30 s to remove
any effects of the previous test." Section 6.2.1: "Prior to each measurement, a preshear
of 300 s-1 was applied for 30 s."

On the step sheet itself, section 6.4 and Fig. 6.6 ("Transient stress measurements at
1 s-1 after an instantaneous step down from 300 s-1"): "An initial increase in the
stress is observed for species with rouleaux formation followed by a longer decay. This
is due to the initial build-up of structure followed by the formation of a depletion
layer. The observation of syneresis in only those species whose RBCs form rouleaux
lends credence to the hypothesis that syneresis is a result of" [the sentence continues
on the next page with the excluded-volume argument].

Consequences for what is built here:

- Six of the 16 steady points (0.1, 0.2, 0.35, 0.7, 1, 2 1/s) are BY DESIGN the maximum
  of a 40-50 s transient after a 300 1/s preshear; the other ten are 15-sample
  averages at the end of the hold. That is what `step_response.md` found empirically
  for the 1 1/s point (bit-identical to the transient's maximum on 14/20 donors). Donor
  C's steady value is the 14 s sample rather than the 6 s maximum, and donor D's is not
  in the step sheet at all, so the stated rule was not followed exactly on 2 of 20.
- The authors attribute the decay AFTER the peak to syneresis (a wall depletion layer),
  not to bulk thixotropy. The decay is therefore not a bulk-rheology signal for a
  constitutive model to fit, and the "structural family is refuted on shape" reading
  in the first version of `step_response.md` (2026-09-05, commit 39b0988) was wrong:
  a monotone rise to a maximum is exactly what first-order structure build-up gives,
  and the non-monotone part is, on the authors' interpretation, a measurement artefact.
  The thesis says "lends credence to the hypothesis": it is their attribution, not a
  demonstrated mechanism, and this repository has no independent test of it.
- The Chapter 2 Casson table (thesis Table 2.2, donor A: 5.09 mPa, 3.02 mPa s) equals the
  Farrington 2023 table to the printed digits and is captioned "Casson plot and linear
  fit", which is the sqrt-space line reproduced below.

## The steady 1 1/s point is a transient peak (step_response.md)

On 14 of 20 donors the steady sheet's 1 1/s stress is bit-identical to the maximum of
the step-down transient at 4-6 s (donor C: identical to the 14 s sample; 4 more within
1e-4 of a sample; donor D 3 % off). That is the stated protocol above, seen in the data. The value is
7-46 % above the stress 40-50 s later. Refitting Casson with that one point replaced by
the last sample moves tau_y by only -1 to -7 % (a lower bound: the five other peak-hold
points cannot be replaced, no step exists for them). Every donor's stress after the step
rises to 4-6 s and then decays; the rise is structure build-up, the decay is attributed
by the authors to syneresis.

## Steady shear, human: 20 donors, 16 rates in [0.1, 700] 1/s

Donor letters A-U without E (20 donors). `DonorOutliers.xlsx` lists which donors the
paper's own analysis treated as outliers; it is not applied here, every donor gets a
record and the flag is carried in `extra`.

## The authors' Casson values are a convention, reproduced exactly

The Casson table's (tau_y, eta_c) per donor are reproduced to max |dev| = 1.2e-8 by an
UNWEIGHTED least-squares line through (sqrt(gd), sqrt(tau)) over all 16 points:
`sqrt(tau) = sqrt(tau_y) + sqrt(eta_c) sqrt(gd)`. That objective weights the high-rate
end (sqrt(700) = 26 against sqrt(0.1) = 0.3), so the "yield stress" it returns is the
intercept of a line dominated by the 100-700 1/s points. `biorheo.identify` fits
log-viscosity residuals (equal weight per decade) and returns tau_y 0.4-71 % lower
(median 33 %) and eta_c 0.3-35 % higher on the same data. Neither is wrong; the Casson
form does not describe the curve: on 18 of 20 donors the log-residuals ordered by rate
fall into 3 runs against 9 expected (z = -3.1; the minimum possible at n = 16 is -3.62
for 2 runs; donor D is the exception at 8 runs, z = -0.5), so its parameters depend on
the objective by that much. The same test on the authors' own sqrt-space residuals
flags 18 of 20, so the misfit is the form's, not the objective's (referee check,
2026-09-05).
Cross and Carreau-Yasuda fit the same curves with rmse_log 0.03 against Casson's 0.15,
and are flagged structured on 2 and 1 donors; the price is that lam and a are degenerate
or bounded on a third of donors. Full table: `model_comparison.md`.
`scripts/records_from_horner.py` writes both: the biorheo fit as the record, the
authors' values and the reproduction deviation in `extra`.
