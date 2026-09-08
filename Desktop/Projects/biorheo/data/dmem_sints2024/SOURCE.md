# DMEM + FBS rheology (Sints et al. 2024) — the project's first COMMERCIALLY LICENSABLE dataset

Downloaded 2026-09-07/08. **CC BY 4.0.** This matters more than the data does: every one of
the 115 records in the registry today is `commercial_use: false`, so nothing in it can
appear in a public demo, a webinar figure, a short course or a paid deliverable. This is
the first candidate that can (NEXT_STEPS item 20).

| field | value |
|---|---|
| Zenodo record | 13828992 |
| DOI | `10.5281/zenodo.13828992` |
| licence | **CC BY 4.0** — commercial use permitted with attribution |
| source paper | Šints, Cīmurs, Birjukovs, Driķis, Goluba, Jēkabsons, Parfejevs, Riekstiņa, Mozoļevskis, Rimša, Kitenbergs, *"Physical model of serum supplemented medium flow in organ-on-a-chip systems"*, arXiv:2409.13650 |
| instrument | Anton Paar MCR 502, SN81094851, cone-plate **CP50-2** SN26597, `d = 0.208 mm` |
| measured | 2023-03 to 2023-10 |

The download is committed (5 KB) rather than gitignored, because it is small enough that
byte-level reproducibility is cheaper than a re-download.

| file | md5 |
|---|---|
| `DMEM_FBS_rheology.zip` | `5685d4deea3279e1428ef4da1c310319` |
| `DMEM_FBS_rheology/DMEM_and_10pct_FBS.txt` | `65d100041ff039ef171871d1da819476` |
| `DMEM_FBS_rheology/DMEM_and_10pct_FBS_37C.txt` | `aea4de63345b60df67d2fa785c1829d1` |
| `DMEM_FBS_rheology/DMEM_and_5pct_FBS.txt` | `5bf525af95f60c02282233dd2d3b6e77` |
| `DMEM_FBS_rheology/DMEM_and_1pct_FBS.txt` | `4333c4c1681d87bc04758ea2dc6ac445` |
| `DMEM_FBS_rheology/readme.txt` | `8554f9e503b7fd0506f12793f9dd6ebb` |

## What the deposited readme says, verbatim

> All data sets contain results of rheological measurements for Dulbecco's Modified Eagle
> Medium (DMEM) supplemented with Fetal Bovine Serum (FBS). All measurements performed with
> Anton Paar MCR 502 rheometer in a cone-plate mode.
>
> DMEM_and_10pct_FBS: DMEM supplemented with 10% FBS at room temperature
> DMEM_and_5pct_FBS: DMEM supplemented with 5% FBS at room temperature
> DMEM_and_1pct_FBS: DMEM supplemented with 1% FBS at room temperature
> DMEM_and_10pct_FBS_37: DMEM supplemented with 10% FBS at T=37 C

## Structure, and two traps in it

Each `.txt` is an Anton Paar RHEOPLUS export holding **more than one interval** — a header
block, a `Meas. Pts.` table, then another header and another table. Eight intervals across
the four files. **An importer that reads the first table only silently discards half the
data**, which is the same shape as the two-concatenated-ramps trap already recorded for the
Horner workbook (NEXT_STEPS item 7).

Columns: `Meas. Pts. | Shear Rate [1/s] | Shear Stress [Pa] | Viscosity [Pa·s] |
Speed [1/min] | Torque [µNm] | Status`.

**`DMEM_and_10pct_FBS.txt` interval 1 opens at a shear rate of 2.24e-7 1/s with an apparent
viscosity of 8710 Pa s** and a torque of 0.0641 µN·m. That is the instrument's own
zero-rate point, not a measurement of anything. `scripts/dmem_candidate_scan.py` excludes
it by a stated rule (rate below 1e-3 1/s) and prints the count.

## Two things that are NOT in this dataset, and both bind

1. **"Room temperature" is not a number.** Three of the four samples are labelled only
   "room temperature". The envelope schema wants `temperature_C`, and guessing 20 or 25
   would be inventing a measurement. The RT/37 °C pair can still give an *empirical* shift
   factor (`collapse_shift` needs no temperature at all), but it cannot give an activation
   energy.
2. **The rates stop at 1 1/s**, 10 1/s on one interval. A fine needle imposes 1e4–1e5 1/s,
   so anything injectability-shaped built on this would be four to five decades outside the
   measured range, against 1.42 decades for the blood records.

## Due diligence already done — and a hypothesis of mine that it refuted

`scripts/dmem_candidate_scan.py` -> `candidate_scan.md`. **It creates no record and ships
no verdict.**

The paper claims DMEM + FBS *"requires a more complex treatment (i.e., is a non-Newtonian
fluid)"*. Reading the raw file, the **shear stress** column looked close to affine in shear
rate with a non-zero intercept, and a viscosity formed as `tau/gd` from such a series falls
like `tau_0/gd + eta_p` — which looks like shear thinning whatever `tau_0` is. So the scan
asked whether a two-parameter Bingham form described the curves as well as a two-parameter
power law.

**It does not. The power law wins on 7 of the 8 intervals**, same parameter count, same
error measure, with `n` between 0.30 and 0.70. The paper's shear-thinning reading is the
one the data supports and my offset hypothesis was wrong. It is recorded rather than
deleted, in the scan's docstring and on its page.

What remains open is the **magnitude, not the shape**: the room-temperature intervals sit
14–32x water and the 37 °C ones 2–5x, for a fluid that is ~99 % water. The paper's own
abstract says its *"observed shear stress values significantly exceed literature reports"*,
so the authors are aware of it. **Read the paper before building a record.**

## What must happen before a record is written

Nothing here is a record and none of the numbers in `candidate_scan.md` may be quoted as
one. In order:

1. Read arXiv:2409.13650 — what do the authors attribute the magnitude to, what temperature
   is "room temperature", and do they fit a constitutive form?
2. Resolve the DOI into `doi_cache.json`.
3. Fit the registry's own families through `biorheo.identify`, with profile intervals and
   an identifiability verdict, rather than the three hand fits in the scan.
4. Pre-register the gate BEFORE writing the generator, as every other record has been.
5. Only then decide whether this becomes the registry's first `commercial_use: true`
   record — and note that even then it cannot carry an injectability claim, for the
   four-to-five-decade reason above. Its value is that it is **licensable**, that it is a
   **biofluid**, and that it is the only data in the project with the **same fluid at two
   temperatures**, which is the one thing `constitutive.shift.collapse_shift` (item 21) has
   never been run on outside a synthetic pair.
