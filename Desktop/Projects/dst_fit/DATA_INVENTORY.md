# What is in the three downloaded datasets

Surveyed 2026-08-27. All three are Edinburgh DataShare deposits under the DataShare
depositor agreement, so acknowledge the depositors and cite the papers.

The zips are still here alongside the extracted folders: `DS_10283_834.zip` (78 KB),
`DS_10283_886.zip` (2.6 MB), `DS_10283_3357.zip` (37.5 MB). Safe to delete once the
readers are written.

---

## 1. Guy, Hermes & Poon, PRL 115, 088304 (2015)

`Towards a unified description of the rheology of hard-particle suspensions/`
arXiv:1506.09041. **This is the G4 reference.**

### Correction: it is PMMA, not silica

`DATA_SOURCES.md` called this a silica paper. It is not. The system is PHSA-stabilised
poly-methylmethacrylate spheres, d = 268 to 3770 nm, in a density-matching CHB/decalin
mixture. Silica appears only as literature comparison points in panel (b). PMMA is still
hard, monodisperse and non-swelling, so it serves the purpose, and it is better
characterised than most silica work. But the label was wrong.

### Correction, 2026-08-27 evening: two solvent viscosities, not one

An earlier draft of this file recorded a single **eta_f = 2.4 mPa.s**. The paper says
"the solvent viscosities were eta_f = 2.83 mPa.s for the large and 2.4 mPa.s for the small
particles at 19 degC". Using 2.4 throughout makes every d = 3770 nm relative viscosity 18%
too high, and that error goes into phi_m.

Both values are now confirmed from the deposit itself rather than from the sentence:
dividing each flow curve's lowest-stress viscosity by the matching `eta_1` entry in
`Fig1c.xlsx` returns the solvent viscosity, and gives **2.35 to 2.41 mPa.s for `Fig1a`**
and **2.66 to 2.81 mPa.s for `Fig1b`**. `verify_guy_size_assignment()` performs that check.

### Correction: sigma* = 200 Pa belongs to d = 404 nm, not d = 3770 nm

The paper: sigma* ~= 200 Pa for the 404 nm particles, ~= 2 Pa for the 3770 nm particles,
with a least-squares size scaling **sigma* = B d^alpha, alpha = -1.9(1)** rather than an
exact d^-2. Fitting both series here recovers a ratio of **68x** against the 70x that
alpha = -1.9 predicts, which is an independent confirmation of the size assignment below.

### Files

| File | Sheets | Columns |
|---|---|---|
| `Fig 1a/Fig1a.xls` | 40, 47, 53, 54.5, 56, 58, 59, 61 % | shear rate (1/s), shear stress (Pa), viscosity (Pa.s) |
| `Fig 1b/Fig1b.xlsx` | 19, 41, 47, 51, 54, 56, 58 % | same |
| `Fig 1c/Fig1c.xlsx` | literature compilations, `d=3770 nm eta_1`, `d=3770 nm eta_2`, `d=404 nm eta_1` | phi, relative viscosity, and for eta_2 a `phi_m - phi` column |
| `Fig 1d/Fig1d.xlsx` | 39, 45 % | same |

Row 0 of each sheet is the header, data starts at row 1. Volume fractions are true volume
fractions, not weight fractions, which is what makes this usable where the corn starch was
not.

**The file names do not match the panel letters, and this is now settled.** `Fig1b.xlsx`
holds flow curves; `Fig1c.xlsx` holds the phi-axis compilation. The published caption puts
the eta(phi) compilation and its eta_2 fit in panel **(b)** and the d = 3770 nm flow curves
in panel **(c)**. So `Fig1b.xlsx` is panel (c) and `Fig1c.xlsx` is panel (b). Confirmed
twice: by the caption, and by matching each workbook's viscosity plateaus to the `eta_1`
columns (`Fig1a` -> `d=404 nm eta_1`, `Fig1b` -> `d=3770 nm eta_1`). `Fig1d.xlsx` is still
unchecked and panel (d) is a schematic, so do not use it.

**Two curves per panel are not steady measurements.** "The flow in both (a) and (c) was
unsteady for phi >= 0.56, and points represent temporal averages", plus wall slip at all
stresses below phi_m and massive wall slip above it. Those are the 56% and 58% curves.

From the paper's Figure 1 caption: panel (a) is **d = 404 nm**, stress controlled, downward
rate sweeps below 10 Pa and upward stress sweeps above it. Panel (c) is **d = 3770 nm**,
upward stress sweeps except phi = 58 % which is a downward sweep. Panel (d) is a
*schematic* for d = 350 um calculated from theory, not measurement, so `Fig1d.xlsx` needs
checking before use. The file-to-panel mapping is not verified beyond (a) and (c); confirm
against the figure before trusting a label.

**Flow was unsteady for phi >= 0.56 in both (a) and (c), and those points are temporal
averages.** That matters for a fit that treats each point as an independent steady
measurement.

### The numbers G4 compares against

From the caption, the authors' own least-squares fit to the upper viscosity branch:

> eta_2(phi) = A (1 - phi/phi_m)^(-n),  with  A = 0.20(9),  **phi_m = 0.558(5)**,  n = 2.2(2)

Independently confirmed from the data: the `d=3770 nm eta_2` sheet carries a `phi_m - phi`
column, and phi + (phi_m - phi) = **0.560 exactly** on all six rows, consistent with 0.558(5)
rounded.

Also stated: the onset stress is phi-independent within each size, at **sigma\* ~= 200 Pa
for d = 404 nm** and **~= 2 Pa for d = 3770 nm**, scaling as **d^-1.9(1)** across sizes.

The fitted Wyart-Cates sigma* comes out about 3 to 5x higher than both published figures.
That is a definition difference, not a disagreement: the paper reads sigma* off where
thickening visibly begins, while the fitted sigma* is where f = 1/e. The near-constant
factor across both sizes is what says so.

### Three things dst-fit has pinned that are wrong for this data

This is the useful part, and it has to be fixed before G4 can run honestly.

1. **The exponent.** `models.py` hardcodes `(1 - phi/phi_J)^-2`. The reference fit used
   **n = 2.2(2)**. Comparing a phi_m obtained at n = 2 against one obtained at n = 2.2 is
   not a like-for-like test, and the difference will land in phi_m.
2. **The solvent viscosity.** `ETA_S_WATER = 8.9e-4` is the default. This system is
   **2.4e-3 Pa.s**. Nearly a factor of three.
3. **phi_0.** Pinned at 0.64. The frictionless branch eta_1 here diverges near 0.61, which
   the `d=3770 nm eta_1` sheet supports: relative viscosity 2297 at phi = 0.5984 implies
   phi_0 = 0.611 under an exponent of 2.

**Done, 2026-08-27.** All four (`eta_s`, `phi_0`, `n`, `beta`) are now arguments bound per
dataset with `Model.with_fixed()`; the module defaults are unchanged and pinned by test.
Fitting the frictionless branches at the same exponent as phi_m gives
**phi_0 = 0.6076(22) at d = 404 nm and 0.6123(10) at d = 3770 nm**, agreeing to 0.005 across
a 9x size range. At n = 2 the same branches give 0.6103 and 0.6007, which is why the
exponent has to be settled before the poles are compared.

Sensitivity, measured rather than assumed: moving phi_0 from 0.6076 to 0.64 moves the
weighted-mean phi_m from 0.539 to 0.568 and leaves the reduced chi-squared of the
consistency test **identical at 344**. So phi_0 sets the level of phi_m and has no effect at
all on whether phi_m holds still.

---

## 2. Lin, Guy, Hermes, Ness, Sun, Poon & Cohen, PRL 115, 228304 (2015)

`Hydrodynamic and Contact Contributions to Continuous Shear Thickening.../`
arXiv:1509.02750. **This is the actual silica.**

Shear-reversal experiments separating contact from hydrodynamic stress.

| Group | System | Contents |
|---|---|---|
| `Fig2a`, `Fig2b` | **2 um silica** | `rate_sweep`, `steady_state`, `contact`, `hydro` as shear rate vs **relative viscosity**; plus instantaneous viscosity vs strain after reversal at four rates |
| `Fig2c`, `Fig2d` | 1.4 um PMMA | same structure |
| `Fig3` | simulation | viscosity vs strain after reversal, contact/hydro decomposed by direct evaluation |

CSVs are two-column and headerless: shear rate, relative viscosity.

**Limitation for G4: this is a single volume fraction per system, not a phi-series.** So it
cannot test phi_m consistency, which is what G3 and G4 are about.

**Used 2026-08-27 for gate G10, and it passed.** Its value was exactly as anticipated and
larger than expected. Everything is *relative* viscosity against shear rate, so the reduced
stress `Sigma = eta_r * gamma_dot = sigma / eta_s` cancels the solvent viscosity, which
neither system states. Fitting Wyart-Cates to `Fig_2a_rate_sweep` and comparing the
recovered `sigma*` against where `Fig_2a_contact` overtakes `Fig_2a_hydro` gives:

| system | fitted sigma* | measured crossover | ratio |
|---|---|---|---|
| `Fig_2a`, 2 um silica | 8.70 +/- 0.48 | 11.85 | **0.73x** |
| `Fig_2c`, 1.4 um PMMA | 3.89 +/- 0.20 | 2.12 | **1.83x** |

The decomposition closes exactly (`contact + hydro = total` to 1e-10), and the contact
fraction saturates near 0.7 rather than 1, because hydrodynamic stress never switches off.
`sigma*` is **identical to four figures** as the unmeasured phi_0 sweeps 0.58 to 0.68, so
the phi/phi_m degeneracy does not touch it. The 2 um silica fit is the only well-specified
fit on real data anywhere in this project (runs test p = 0.092).

**Still unused, and both are worth a session.** `Fig2b` and `Fig2d` hold instantaneous
viscosity against accumulated strain after reversal at four rates each: that is the
*relaxation* of the contact network, i.e. the `gamma_0` of Richards et al.'s Eq. 4, which
is fitted nowhere in this package. `Fig3` is the same decomposition from simulation, where
contact and hydrodynamic stress are evaluated directly rather than inferred from reversal,
which is the control on G10's own method.

---

## 3. Richards, Royer, Liebchen, Guy & Poon, PRL 123, 038004 (2019)

`Competing Timescales Lead to Oscillations in Shear-Thickening Suspensions/`
arXiv:1902.07655. Already in memory as the `WC_CLOCK` source.

**This is corn starch in 50 wt% glycerol-water, and it is the most directly useful thing in
the download for the N4 question.**

`Figure1/Plotting/` holds a seven-point series: `38wt%.csv` through `50wt%.csv` in steps of
2, each 31 to 36 rows with a real header:

```
Shear-rate,Error,Stress
```

Shear rate in 1/s, its error in 1/s, external stress in Pa. Time-averaged flow curves,
stress controlled. `Figure1/Data/` and `Figure3/Data/cs50glyc/` hold the underlying
time series (time, strain, strain rate, torque, strain acceleration) at 47 and 48 wt% for
every applied stress, roughly 20 stresses each.

### Why this matters for N4

N4 found phi_m sliding from 0.492 to 0.580 across your own 2023 corn-starch series and left
open whether that is a property of corn starch or of that data. **This is an independent
corn-starch series, different lab, different solvent, different rheometer.** Running it
through the same pipeline answers the question directly:

- if phi_m drifts here too, the drift is corn starch and the model is being applied outside
  its assumptions, which is a publishable observation rather than a data problem;
- if phi_m holds still here, the 2023 series has something wrong with it and the fault is
  local.

Either way it is a cleaner test than any of the four candidates currently in the spec, and
it costs one afternoon. **Do this before the literature cross-check.**

One conversion needed: these are **weight** fractions, so a volume fraction requires the
corn-starch and solvent densities. The wt%-to-vol% conversion is exactly the ambiguity
candidate (3) in the spec named, so carry the conversion explicitly rather than folding it
into a constant.

---

## Suggested order of work

1. Make `eta_s`, `phi_0` and the divergence exponent per-dataset parameters instead of
   module constants. Nothing below is honest until this is done.
2. Run the Richards corn starch. It answers the open N4 question with a like-for-like
   comparison.
3. Run Guy `Fig1a` (d = 404 nm) and `Fig1c` (d = 3770 nm) against phi_m = 0.558(5) for G4,
   using n = 2.2 and eta_f = 2.4 mPa.s.
4. ~~Use Lin `Fig2a` as the single-curve silica check.~~ **Done 2026-08-27, gate G10,
   passed.** Next in this deposit: `Fig2b`/`Fig2d` for the contact-network relaxation
   strain, and `Fig3` as the simulation control on the reversal decomposition itself.
