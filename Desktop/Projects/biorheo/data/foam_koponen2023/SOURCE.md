# Source — PVOH application-foam pipe rheology (Koponen et al. 2023)

**Deposit.** Zenodo `10.5281/zenodo.8334762`, *Pipe rheology of wet aqueous application
foams*. Landing page `https://zenodo.org/records/8334762`.
**Licence: CC BY 4.0** (`cc-by-4.0`, SPDX), read from the **DataCite** record
`https://api.datacite.org/dois/10.5281/zenodo.8334762`, whose `rightsList` gives
`Creative Commons Attribution 4.0 International` + `info:eu-repo/semantics/openAccess`.
Creators: Antti Koponen, Janika Viitala, Atsushi Tanaka, Baranivignesh Prakash,
Olli-Ville Laukkanen, Ari Jäsberg (VTT Technical Research Centre of Finland).

**Paper.** Chemical Engineering Science (2023), `10.1016/j.ces.2023.119282`, carried on the
DataCite record as `IsReferencedBy`. The accepted manuscript is IN the deposit and is the
copy read here.

**Route note, worth carrying.** `NEXT_STEPS.md` item 20 recorded the Zenodo REST API as the
way to search, with the caveat that it "returns HTTP 504 to a default `python-urllib`
User-Agent and works with an explicit one". On 2026-09-08 it returns **504 to an explicit
User-Agent too**, from two independent clients. **DataCite answers**, indexes the same
Zenodo records with the same licence metadata, and additionally covers figshare, Dryad,
4TU and Mendeley Data. Zenodo's *landing pages and file downloads* work normally; it is
only `/api/` that is down. `scripts/datacite_candidate_search.py` is the replacement route.

## Files (raw/, downloaded 2026-09-08)

| file | bytes | md5 |
|---|---|---|
| `Set_1_Data_PVOH_foam_rheology.csv` | 13639 | `4c9b9f8db41ff6fba8cf224711dc8d60` |
| `Set_2_Data_PVOH_foam_rheology.csv` | 29718 | `f3f63abba394b1c4543098bbc3828dfc` |
| `Data_PVOH_foam_rheology.xlsx` | 103543 | `f7684d3b2bfe82a8646c84913b86b167` |
| `Data_column_descriptions.csv` | 905 | `235809792efc3cfcd15a072d35a8d3f5` |
| `PVOH shear viscosity data.xlsx` | 14179 | `7f58b8c041e36c92cb41f7ce08c232d0` |
| `Koponen_pipe_rheology_foams.pdf` | 1363628 | `448ce1130442439aa75c1cbc12d4ccc0` |

`paper_text.txt` is `pdfplumber`'s extraction of the PDF, kept so a quotation can be
checked without re-running the extractor. Equation glyphs in it are mangled; **every
equation used in this project was read off the PDF page, not off that text file.**

## Columns (from `Data_column_descriptions.csv`, verbatim)

`D nominal` · `D hydrodynamic` (the real pipe diameter determined by hydrodynamic testing)
· `PVOH concentration` [% by weight of resin in the foaming solution] · `Foam density
setpoint` [g/l] · `Sample density` [g/l] (measured at 1 atm at the discharge) · `Flow rate`
[l/min] (at the absolute pressure of column 8) · `Pressure gradient` [kPa/m] (measured at
pressure close to 1 atm) · `Absolute pressure` [kPa] · `Slip velocity` [m/s] (**measured
with video imaging**) · `Sauter mean diameter` [µm].

## Traps found before any analysis

1. **`D nominal` is not the pipe.** The three pipes are 8.15, **11.90** and 19.00 mm
   hydrodynamic, i.e. R = 4.075, **5.95**, 9.50 mm. Using the nominal 12 mm for the middle
   pipe is a 0.8 % error in `R` and, through `Q ~ R^4` for a Newtonian fluid, ~3.4 % in a
   flow rate. Use `D hydrodynamic`.
2. **`Sample density` is not `Foam density setpoint`.** At setpoint 200 g/l the samples
   run 182–194 g/l in the first four rows alone. The paper's Table 2 is indexed by air
   content derived from the NOMINAL density, so a gate that compares against Table 2 must
   bin by setpoint and say so; anything computed per row should use the sample density.
3. **`Set_1` and `Set_2` are different campaigns and the deposit never says how they
   relate.** `Set_1` is unbalanced (24 cells, 4 PVOH concentrations, mostly the 12 mm
   pipe); `Set_2` is a clean 3 x 3 x 5 factorial. Item 27 uses `Set_2` for every gate and
   measures the `Set_1`/`Set_2` difference as an explicitly un-graded observation (G27.7).
4. **The deposit's own rheometer file disagrees with the paper's Table 1.**
   `PVOH shear viscosity data.xlsx` is the *solution* (not foam) flow sweep and gives
   ~1.56 mPa·s for 0.5 % PVOH at 1 1/s where Table 1 prints 1.32 mPa·s. Item 27 uses
   Table 1, because Table 1 is what Eq. (19) was regressed with; the discrepancy is
   quantified in the report rather than resolved.
5. **The rheometer file is the SOLUTION, not the foam.** There is no rheometer flow curve
   of the foam anywhere in the deposit, so this dataset cannot test a rheometer-to-pipe
   round trip. What it can test is pipe-to-pipe, which is the same extrapolation in kind.

## What the paper declares about its own errors

- slip velocity: **better than 3 %** (§2, from the scatter within each trial point);
- flow rate: **maximum 5 %** (§2, the two-pipe splitter experiment);
- surface tension: 0.3 mN/m; dilatational modulus: 0.02 mN/m.

These are the only accuracy scales in item 27 that this repo did not invent.

## Table 1 (PVOH rows), read off the PDF, 22 °C

| PVOH % (w/w) | viscosity [mPa·s] | surface tension [mN/m] | dilatational modulus [mN/m] |
|---|---|---|---|
| 0.5 | 1.32 | 52.7 | 6.16 |
| 1.0 | 1.59 | 50.5 | 4.07 |
| 3.0 | 3.45 | 47.8 | 2.34 |
| 6.0 | 9.70 | 46.2 | 2.35 |

## Table 2, read off the PDF

Eq. (19) is `mu/mu_l = B(phi) (gd D mu_l / sigma_l)^m(phi)`, regressed **at each air
content separately**, pooling PVOH concentration, pipe radius and bubble size. The text
states `n = 1 + m` where `n` is the flow index.

| `phi` | `B(phi)` | `m(phi)` | R² |
|---|---|---|---|
| 0.70 | 0.466 ± 0.039 | −0.639 ± 0.014 | 0.94 |
| 0.75 | 0.636 ± 0.023 | −0.637 ± 0.006 | 0.98 |
| 0.80 | 0.659 ± 0.016 | −0.691 ± 0.004 | 0.99 |
| 0.85 | 0.718 ± 0.021 | −0.735 ± 0.005 | 0.99 |
| 0.90 | 0.883 ± 0.022 | −0.762 ± 0.006 | 0.99 |

`phi = (rho_l - rho_f)/rho_l` with `rho_l = 1000 kg/m^3` maps the five density setpoints
300/250/200/150/100 g/l onto these five rows exactly.
