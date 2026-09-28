# External anchor audit — is each comparison graded against the source's FIT TABLE?

Work order B3 (`WORK_ORDERS_2026-09-08.md`), NEXT_STEPS item 30. **FIRST PASS, written
2026-09-28 in a 20-minute slice; NOT complete.** Read-only: nothing here fixes anything.
Every row below was built from the deposit's `SOURCE.md` and the anchor's own page, both
already on disk; no source paper was re-opened in this pass, so any cell that would need
the PDF is marked UNVERIFIED rather than guessed.

The question, from item 23: a figure caption states what the authors FED a model (rounded,
shared across conditions); a fit table states what they MEASURED. Grading against the first
is not an external anchor. Two of item 23's four caption values were the right number
against the WRONG ROW, which is how a total mismatch read as a partial one.

Verdict vocabulary, as the brief fixes it: GRADED AGAINST A FIT TABLE / GRADED AGAINST A
CAPTION / NO FIT TABLE EXISTS / UNVERIFIED.

## The five anchors the brief names

| # | anchor (repo page) | source | the compared number came from | source publishes a fit table? used? | verdict | rank to fix |
|---|---|---|---|---|---|---|
| 1 | DMEM/FBS power-law records, G22.1 -> G23.x (`dmem_sints2024/anchor_recheck.md`) | Šints et al. 2024, arXiv:2409.13650, Zenodo 13828992 | **Item 22: figure captions** (n = 0.31/0.32/0.5/0.54, the COMSOL inputs). **Item 23: Table 2** (`K`, `n` per deposited sample, "fit results ... for the experimental rheometry data shown in Fig. 1"), compared AS CURVES over the measured range | YES (Table 2). Not used by item 22; used by item 23 | was GRADED AGAINST A CAPTION; **now GRADED AGAINST A FIT TABLE** | done (item 23). Nothing further |
| 2 | Carbopol 974p / propylene glycol shift anchor, G24.1 (`carbopol_pg_prereg.md`, `carbopol_pg_vinutha2024/`) | Vinutha et al. 2024, arXiv:2312.00251, Zenodo 11389467 | The deposit's own `medium viscosity [Pa s]` column of `Fig 1(c) and Fig S3 residual stress_exp.txt` (0.102 / 0.055 / 0.012 Pa s), corroborated by the paper's Methods. `a_pred = 0.11765` is a ratio of two MEASURED solvent viscosities, not of fitted parameters | The paper fits a three-component model (`sigma = sigma_y + k1 gd^0.5 + eta_bg gd`, ref. [18]); the deposit column carries only `sigma_y = 8.9 Pa`. **Whether the paper tabulates `k1` and `eta_bg` per temperature is UNVERIFIED** (not read in this pass). G24.9 deliberately does NOT grade a Herschel-Bulkley `tau_y` against the 8.9 Pa (different model form) | G24.1: **NO FIT TABLE EXISTS** for the quantity compared (a measured viscosity ratio) — the right kind of anchor. Yield stress: UNVERIFIED whether a fuller fit table exists | **medium**: if the paper publishes the full three-component parameters, the item-23 style check (compare the two model forms AS CURVES over the measured range, scaled by the 2.5-2.9 % replicate scatter the deposit gives) is available and has never been run |
| 3 | Koponen foam pipe anchor, G27.1-G27.7 (`foam_koponen2023/foam_pipe_anchor.md`) | Koponen et al. 2023, CES 10.1016/j.ces.2023.119282, Zenodo 8334762 | **Table 2 of the paper**, read off the PDF: `B(phi)`, `m(phi)` with the authors' own sigma and R² per air content — a regression table. G27.1 compared one fitted exponent against one published exponent (FAILED 2 of 5); G27.x then compared the published parameters at every reduced point over the measured range (the curve comparison) | YES (Table 2), and it IS what was used. Caveat recorded in `SOURCE.md` trap 4: the paper's Table 1 solution viscosity (1.32 mPa s at 0.5 % PVOH) disagrees with the deposit's own rheometer file (~1.56 mPa s); Table 1 was used because Eq. (19) was regressed with it, and the discrepancy is quantified, not resolved | **GRADED AGAINST A FIT TABLE** | low. Open only: the Table 1 vs deposit-file disagreement is a source-internal inconsistency, already recorded |
| 4 | Horner blood Casson parameters (`horner2020/SOURCE.md`, `records_from_horner.py`, `physiology_and_casson_farrington2023.csv`) | Horner 2020 thesis (Table 2.2) and Farrington, Jariwala, Armstrong 2023, Rheol. Acta 10.1007/s00397-023-01402-2 (`HornerData_CassonFitTable.xlsx`) | The per-donor Casson fit table (tau_y, eta_c with the authors' standard errors). Reproduced to max dev 1.2e-8 by an unweighted line in (sqrt gd, sqrt tau); biorheo's own log-viscosity fit differs by up to 71 % in tau_y, and that difference is the recorded finding (a yield stress is a fit-objective convention) | YES (a fit table, per donor, with errors), and it IS what was used; the thesis Table 2.2 agrees to the printed digits | **GRADED AGAINST A FIT TABLE** | none |
| 5 | mAb concentration law, the "63 mg/mL" and "doubles every 21.8 mg/mL" headline (`case_studies/mab_injectability/results.md`, `mab_prasnikar2025/SOURCE.md`) | Prašnikar et al. 2025, Int. J. Pharm. 674, 125451, Zenodo 14409284 | The slope 0.03131 mL/mg (main series) / 0.032334 (higher-concentration series) is **back-solved by this repo** from the deposit's `viscosity reduction [%]` column via the Methods' definition (`100 (1 - eta/eta_base(c))`, `eta_base` an exponential calibration curve). The intercept does NOT reproduce (R² 0.4575 on 49 rows, 0.945 on series B alone) | **UNVERIFIED whether the paper publishes its exponential calibration curve's parameters** (a fit the authors made and the repo re-derived). If it does, the headline's slope was graded against nothing external, when a published number existed | **UNVERIFIED** | **highest**: the case-study headline goes through this slope, and the one comparison that could anchor it (the authors' own calibration parameters) has not been looked for. Check the open-access PDF (Ljubljana repository record 167794) §2.2.4 and any figure of viscosity vs concentration for the fitted exponential; if present, compare as curves over 159-212 mg/mL and scale by the deposit's per-formulation SD |

## What this pass did NOT do (the brief says "every external comparison in the repo")

Not enumerated, not graded — the next session must sweep these before the item is closed:

- the OpenFOAM plane-Couette / Poiseuille anchors (`data/openfoam/`): analytic, so the
  "fit table" question is moot, but the record they reproduce should be named;
- the Funke 2016 glide-force bracket (doi:10.1002/jps.23945, 1.3-1.5 N) quoted in item 28
  and in work order N4.2 — a literature NUMBER, where in the paper it comes from is not
  recorded on disk as far as this pass saw;
- the Roche 2018 "<12 %" model and the Hagenbach/entrance-loss literature in
  `data/injection_force_terms_lit_2026-09-08.md` (B4's territory; overlap to be resolved
  with B4, not duplicated here);
- the retracted ISO 9626 comparison (already withdrawn; should appear in the table as the
  worked example of a number quoted from memory);
- every `doi` in `registry/doi_cache.json` that is cited by a record `extra` field — a grep
  of `extra:` blocks for `published` / `source_value` keys is the mechanical way to find
  comparisons this pass missed.

Suggested next command for the sweep (read-only):

```
cd C:\Users\compl\Desktop\Projects\biorheo
rg -n "published|Table [0-9]|caption" src/biorheo/registry/records data case_studies --glob "*.md" --glob "*.yaml" -l
```

## Ranking (what needs fixing, highest first)

1. **mAb calibration slope (row 5)** — find out whether the source publishes the
   exponential calibration parameters; if yes, the headline has an external anchor it never
   used. If no, record NO FIT TABLE EXISTS and the row closes.
2. **Carbopol-PG three-component fit (row 2)** — same question for `k1`, `eta_bg`; enables a
   curve-vs-curve comparison across model forms that G24.9 currently forgoes.
3. Rows 1, 3, 4 need nothing.

No bar was touched, no number was changed, and no record was regenerated in this pass.
