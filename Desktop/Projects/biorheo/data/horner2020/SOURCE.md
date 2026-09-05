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
the LAOS sheets are untouched. The measurement protocol (time per steady point,
peak-hold or average, preconditioning) is NOT in the workbook's ReadMe. The papers the
data belongs to, resolved via Crossref on 2026-09-05 and not yet read (paywalled):
Horner, Armstrong, Wagner, Beris, J. Rheol. 62 (2018) doi 10.1122/1.5017623 and
J. Rheol. 63, 799 (2019) doi 10.1122/1.5108737 ("steady and transient shear"). Reading
the 2019 methods section is the next step for the protocol question below.

## The steady 1 1/s point is a transient peak (step_response.md)

On 15 of 20 donors the steady sheet's 1 1/s stress is bit-identical to the maximum of
the step-down transient at 4-6 s (donor C: identical to the 14 s sample; 4 more within
1e-4; donor D 3 % off). So the lowest-rate end of the "steady" curve is a peak-hold
value, 7-46 % above the stress 40-50 s later. Refitting Casson with that one point
replaced moves tau_y by only -1 to -7 %, so the records survive it; what does not
survive is the structural family: every donor's stress after the step OVERSHOOTS
(rise to 4-6 s, then decay), and a scalar structure with first-order kinetics is
monotone at constant rate.

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
form does not describe the curve: on 19 of 20 donors the log-residuals ordered by rate
fall into 2 runs against 9 expected (z = -3.1, which is the floor for n = 16; donor D is
the exception at z = -0.5), so its parameters depend on the objective by that much.
Cross and Carreau-Yasuda fit the same curves with rmse_log 0.03 against Casson's 0.15,
and are flagged structured on 2 and 1 donors; the price is that lam and a are degenerate
or bounded on a third of donors. Full table: `model_comparison.md`.
`scripts/records_from_horner.py` writes both: the biorheo fit as the record, the
authors' values and the reproduction deviation in `extra`.
