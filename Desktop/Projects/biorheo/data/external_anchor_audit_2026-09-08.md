# External anchor audit — is each comparison graded against the source's FIT TABLE?

Work order B3 (`WORK_ORDERS_2026-09-08.md`), NEXT_STEPS item 30. First pass 2026-09-28
08:30 (five named anchors, from files on disk only); **completed 2026-09-28 09:20** by
reading the three sources that were UNVERIFIED (Prašnikar 2025 accepted PDF from the
Ljubljana repository; Vinutha 2024 arXiv v1 full text; Biermann 2025 full text via Europe
PMC `PMC12251366`) and by sweeping the rest of the repo for external comparisons (every
DOI outside `raw/` and `candidate_search_*.json`, every `published|Table|caption` hit under
`data/` and `case_studies/`, every `extra.*` key in `registry/records/`). Read-only:
nothing here fixes anything. Every number quoted from a source is quoted with where it sits
in that source.

The question, from item 23: a figure caption states what the authors FED a model (rounded,
shared across conditions); a fit table states what they MEASURED. Grading against the first
is not an external anchor. Two of item 23's four caption values were the right number
against the WRONG ROW, which is how a total mismatch read as a partial one.

Verdict vocabulary, as the brief fixes it: GRADED AGAINST A FIT TABLE / GRADED AGAINST A
CAPTION / NO FIT TABLE EXISTS / UNVERIFIED. One more was needed and is used sparingly:
**NOT AN ANCHOR** — the repo never compared a number against this source, so there is
nothing to grade.

## Every external comparison in the repo

| # | anchor (repo page) | source | the compared number came from | does the source publish a fit table? was it used? | verdict | rank to fix |
|---|---|---|---|---|---|---|
| 1 | DMEM/FBS power-law records, G22.1 → G23.x (`dmem_sints2024/anchor_recheck.md`) | Šints et al. 2024, arXiv:2409.13650, Zenodo 13828992 | Item 22: figure captions (n = 0.31/0.32/0.5/0.54, the COMSOL inputs). Item 23: **Table 2** (`K`, `n` per deposited sample, "fit results ... for the experimental rheometry data shown in Fig. 1"), compared AS CURVES over the measured range | YES (Table 2). Not used by item 22; used by item 23 | was GRADED AGAINST A CAPTION; **now GRADED AGAINST A FIT TABLE** | none (item 23) |
| 2a | Carbopol 974p / propylene glycol **shift factor**, G24.1 (`carbopol_pg_prereg.md`) | Vinutha et al. 2024, arXiv:2312.00251, Zenodo 11389467 | The deposit's own `medium viscosity [Pa s]` column (0.102 / 0.055 / 0.012), a MEASURED solvent viscosity, corroborated by the paper's text ("η_s = 0.055 Pa s", "η_s = 0.012 Pa s" in the Fig. 1(b) caption). `a_pred = 0.11765` is a ratio of two measurements, not of fitted parameters | No fit table is involved and none is needed | **NO FIT TABLE EXISTS** (the right kind of anchor: a measured quantity) | none |
| 2b | Carbopol 974p / PG **yield stress 8.9 Pa**, deliberately NOT graded (G24.9) | same | **A CAPTION value.** arXiv v1, Fig. 1(a) caption, verbatim: *"Lines through the data are fits to the three component model [18]; the resulting yield stress is σ_y = 8.9 Pa."* Searching the full text for "three component", "Herschel" and "Caggioni" (ref. 18) returns only this caption and the reference list: `k_1` and `η_bg` of the three-component fit are **not published** in v1 | NO. One number, in a caption, one value for all three temperatures | **NO FIT TABLE EXISTS**; the only published number is caption-grade, and G24.9's refusal to grade a Herschel–Bulkley `tau_y` against it was right for a second reason (it is not even a per-temperature fit result) | none. The item-23-style curve-vs-curve comparison across model forms suggested in the first pass is **not available**: the three-component parameters are not published |
| 3 | Koponen foam pipe anchor, G27.1–G27.7 (`foam_koponen2023/foam_pipe_anchor.md`) | Koponen et al. 2023, CES 10.1016/j.ces.2023.119282, Zenodo 8334762 | **Table 2 of the paper**, read off the PDF: `B(φ)`, `m(φ)` with the authors' σ and R² per air content — a regression table. G27.1 compared one exponent against one published exponent (FAILED 2 of 5); the follow-up compared the published parameters at every reduced point over the measured range | YES, and it IS what was used. Recorded caveat (`SOURCE.md` trap 4): Table 1's solution viscosity (1.32 mPa s at 0.5 % PVOH) disagrees with the deposit's own rheometer file (~1.56); Table 1 used because Eq. (19) was regressed with it | **GRADED AGAINST A FIT TABLE** | low; the Table 1 vs deposit-file disagreement is source-internal and already recorded |
| 4 | Horner blood Casson parameters (`horner2020/SOURCE.md`, `records_from_horner.py`, `extra.authors_casson_*` on all 20 records) | Horner 2020 thesis (Table 2.2) and Farrington, Jariwala, Armstrong 2023, Rheol. Acta 10.1007/s00397-023-01402-2 (`HornerData_CassonFitTable.xlsx`) | The per-donor Casson fit table (τ_y, η_c with standard errors). Reproduced to 1.7e-9 … 2.1e-9 per record (`authors_casson_reproduced_to`); biorheo's log-viscosity fit differs by up to 71 % in τ_y and that difference is the finding (a yield stress is a fit-objective convention) | YES (per donor, with errors), and it IS what was used; thesis Table 2.2 agrees to the printed digits | **GRADED AGAINST A FIT TABLE** | none |
| 5 | mAb concentration law: "doubles every 21.8 mg/mL", "63 mg/mL" (`case_studies/mab_injectability/results.md`, `mab_prasnikar2025/SOURCE.md`) | Prašnikar et al. 2025, Int. J. Pharm. 674, 125451, Zenodo 14409284; accepted PDF read 2026-09-28 (13 pages + refs) | The slope 0.03131 mL/mg (main series) / 0.032334 (200-mAb series) is **back-solved by this repo** from the deposit's `viscosity reduction [%]` column via the Methods definition. The paper: §2.2.4 *"the starting formulation samples were prepared with mAb concentrations of 110–200 mg mL⁻¹, and their viscosity was measured. An exponential function η = a·e^(b·c) was fitted to the data"*; §2.2.7 *"analyzed using non-linear regression in R software"*; **Fig. 1** plots the fitted curve with a 95 % CI band for the starting formulation and the aqueous mAb. **`a`, `b` and R² are printed nowhere** — not in the caption, the text, the tables, or the deposited Table S1 | **NO.** The fit exists (Fig. 1, with a CI band) but its parameters are not published, so the repo's back-solve is the ONLY route to `b` and it cannot be graded against the source | **NO FIT TABLE EXISTS** — and this row is *unanchorable*, not merely unanchored | **medium**: the case-study page should say in one sentence that the authors' fitted `b` is not published, so the 21.8 mg/mL doubling is a repo-derived number with no external check. The prose points the paper does give are rounded ("~150 mg mL⁻¹ … 33 mPa s", "190 mg mL⁻¹ … 150 mPa s", aqueous 20 / 98 mPa s): they imply `b` = ln(150/33)/40 = 0.03785 (script-printed below), 1.21× the repo's 0.03131, but "~150" cannot carry a 21 % comparison, so this is NOT a grade, only a note that the two are the same order |
| 6 | **Carbopol ETD2020, all 15 registry records** (`records_from_carbopol.py` ← `../carbopol/hb_fits_etd2020.csv` ← `../carbopol/fit_hb.py`), and everything built on them: `injectability_surrogates`, item 18 (`needle_bore.md`), the OpenFOAM anchors, the multistart study | deposit 10.25835/hjn2qcvy (Biermann & Mesecke 2024, CC BY-NC-SA); paper Biermann et al., *Materials* 18, 3164 (2025), 10.3390/ma18133164, full text via Europe PMC PMC12251366 | biorheo itself makes NO comparison — the records carry the paper DOI as `role: paper` and nothing else. The sibling `carbopol/fit_hb.py` (line 88) prints *"paper reports tau_y max 41.1 Pa at 0.5 wt%, k 3.4-15.9, n 0.40-0.53 … reproduced"*, and `CARBOPOL_RHEOLOGY_LITERATURE.md` line 38 repeats it: a comparison of one maximum and two RANGES | **NO.** Europe PMC full text: *"yield stress τ_y, determined by Herschel–Bulkley modeling, ranged from below 1 Pa at 0.1 wt% to above 41 Pa at 0.5 wt%"*, *"consistency index k between 3.4 and 15.9, and a flow index n in the range of 0.40 to 0.53"*, *"Parameter fitting was performed using nonlinear regression across the full shear rate range"*; per-formulation values appear only as **heatmaps (Figs. A1, A2)** and representative curves (Fig. 4). No table of the 15 (τ_y, k, n) exists | **NO FIT TABLE EXISTS.** The sibling's range check is the strongest comparison the source allows; it is a RANGE comparison (a range is not a pairwise comparison — one max and two intervals cannot fail on any single condition except the 0.5 wt% maximum, which matched: 41.11 vs "above 41") | **medium**: (i) biorheo's records should carry the range check in `notes` (they currently cite the paper and say nothing about agreement); (ii) the sibling's word "reproduced" overstates a range check and should say "the 15 fits fall inside the paper's published ranges and reproduce its stated maximum". Digitising Figs. A1/A2 is forbidden by the repo rule and would not be a fit table anyway |
| 7 | Device friction placement, item 28 / `case_studies/device_friction/results.md` and `.json` | Funke et al. 2016, J. Pharm. Sci., doi:10.1002/jps.23945 (RESOLVED) and a *BioProcess International* trade article with **no DOI** | `results.json`: `"glide_note": "1.3-1.5 N glide, from a trade article with NO DOI … treated as UNVERIFIED"`; `breakloose 3.0-9.0 N` from Funke (silicone-free / siliconised). The survey (`injection_force_terms_lit_2026-09-08.md` §1) lists the glide as 1.28–1.45 N from the trade article and Funke's numbers as ~3 N / ~9 N breakloose | Not a fit table question (force measurements, not constitutive fits); the correct question is *which source*, and `results.json` gets it right | **UNVERIFIED**, correctly labelled on the page. **But `WORK_ORDERS_2026-09-26.md` N4.2 mis-attributes it**: *"with the Funke 2016 glide bracket (doi:10.1002/jps.23945, 1.3–1.5 N)"* — the 1.3–1.5 N is the trade article's, not Funke's, and Funke publishes breakloose, not glide | **high, and cheap**: a pre-registered bar (N4.2) names the wrong source for its number. Fix belongs to N4's owner (`WORK_ORDERS_2026-09-26.md`), not to B3: change the citation to "trade-article glide bracket, no DOI, UNVERIFIED" and keep the 1.3–1.5 N; or re-register N4.2 on Funke's breakloose 3–9 N, which `results.json` also carries (ratio 1.18–1.52×) |
| 8 | Injection-force omitted terms (`data/injection_force_terms_lit_2026-09-08.md`, consumed by item 28 and by B4) | Boger, Gupta & Tanner 1978 (Bagley `e(n)`), Cilurzo 2011, Wu 2024, Allmendinger 2014/2015, Rathore 2012, Krayukhina 2020, Verwulgen 2018, MacDonald 2017, Astin 1999 | Literature NUMBERS, each self-graded RESOLVED / PARTIAL / UNVERIFIED in the file, several read from search-engine summaries (Boger's `e = 0.834 … 2.917`, Allmendinger's back-pressure ranges, Verwulgen's 80/40 N) | Not fit tables. None of these numbers is graded against a biorheo computation as a gate; they size omitted terms | **UNVERIFIED** where the file says so (PARTIAL rows); nothing hidden | low for B3; B4 owns the Hagenbach 0.43–0.99 % vs Bagley 1.4–4.6 % reconciliation, which is a *model* disagreement, not a source-reading one |
| 9 | ISO 9626 27G inner-diameter band (`blood_injectability/results.md` "What was retracted", `needle_bore.md`) | ISO 9626:2016, never in hand | 0.184–0.241 mm mixed a regular-wall and a thin-wall figure; **retracted 2026-09-05**, both readings tagged `[UNVERIFIED]`, section rebuilt on `F ~ R^-4` (measured −3.9723), which needs no standard | not applicable | **withdrawn** (the worked example of a number from memory) | none |
| 10 | OpenFOAM plane-Couette and Poiseuille anchors (`data/openfoam/*.md`, `derivations/`) | analytic (Couette line, Buckingham–Reiner, Casson tube flow) + a hand-computed Newtonian control | OpenFOAM's `nu` and `wallShearStress` against `record.eta_at`; A4 against a hand value | no external publication involved | **NOT AN ANCHOR** in B3's sense (internal/analytic); listed so the sweep is complete | none |

Script-printed arithmetic for row 5 (`python -c`, 2026-09-28):

```
starting formulation: b = ln(150/33)/40 = 0.03785 per mg/mL, doubling every 18.3 mg/mL
aqueous mAb:          b = ln(98/20)/40  = 0.03973 per mg/mL, doubling every 17.4 mg/mL
repo back-solved slopes 0.03131 / 0.032334 -> doubling every 22.1 / 21.4 mg/mL; ratio prose/repo = 1.209
```

(The repo page prints "doubles every 21.8 mg/mL" from the series-A slope with its own
rounding; 22.1 here is ln 2 / 0.03131 to the CSV's precision. Not a discrepancy to chase:
the paper's prose points are rounded to two figures and "~150" is not a measurement.)

## Out of B3's scope, listed so nobody re-discovers them

- **`dst_fit`'s own anchors** (Guy, Hermes & Poon 2015 `phi_m = 0.558(5)` in gate G4a; Lin
  et al. 2015 contact fraction in G10; Richards et al. 2019 in N6). `dst_fit` is consumed
  unchanged and audited in its own README. G4a is fit-table grade: it refits the authors'
  deposited `eta_2(phi)` branch and reads their `phi_m - phi` column, and G4c was retired
  2026-09-26 with its number kept (item 29).
- **Roche 2018 "<12 % error" injection-time model** — quoted from a URL in
  `biofluids_venture/docs/market_research_2026-09-06/01_*.md` as market evidence; never
  compared against a biorheo number. Not an anchor.
- **Overcashier 2006 "20 cP / 25 N"** — quoted in `needle_bore_prereg.md` as a design rule
  the field uses; marked UNVERIFIED (trade journal, no DOI) in the survey. Not an anchor.

## What the sweep found, ranked

1. **N4.2's pre-registered bar cites the wrong source** (row 7). The number is fine and
   already labelled UNVERIFIED on the page it came from; the work order welds it to a DOI
   that publishes a different quantity. Fix in the work order before N4 runs — a bar
   graded later against "Funke 2016" would be graded against 3–9 N breakloose and fail for
   a reason that has nothing to do with the physics.
2. **The mAb headline is unanchorable from its source** (row 5): the authors fitted
   η = a·e^(b·c) in R with a 95 % CI band and published neither `a`, `b` nor R². The
   deliverable page should say so in one sentence. The prose points give the same order
   (1.21×) and no more.
3. **The 15 Carbopol records rest on a range check made in a sibling project** (row 6);
   the paper has no per-condition table, so a range check is the ceiling. Carry it into the
   records' `notes`, and soften "reproduced" to what it is.
4. Rows 1, 3, 4 are fit-table grade and need nothing. Rows 2a/2b need nothing: the shift
   anchor compares a measurement, and the 8.9 Pa was correctly not graded and is now known
   to be caption-only.

Three of the four "candidates" the brief expected to be caption-graded turned out to have
**no fit table to grade against at all** (rows 2b, 5, 6). That is a different failure
mode from item 23's: not "graded against the wrong published number" but "there is no
published number, and the page does not say so". The fix in every case is a sentence, not
a computation.

No bar was touched, no number was changed, and no record was regenerated in this audit.
