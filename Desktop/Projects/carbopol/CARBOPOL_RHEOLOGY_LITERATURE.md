# Carbopol rheology: where the usable flow data is

Literature scan, 17 Aug 2026. Written against the bath experiment sketched in
`../carbopol_advisor_meeting_2026-08-03.md`, which needs a yield-stress ladder and a
defensible Herschel-Bulkley parameterisation.

Sorted by how much work it takes to get numbers out. Tier 1 is a file you can load
today; Tier 2 is a table you retype; Tier 3 is a figure you digitise.

---

## Tier 1: raw flow curves you can download

### Biermann & Mesecke (2024), Carbopol ETD2020

**DOI [10.25835/hjn2qcvy](https://doi.org/10.25835/hjn2qcvy)**, CC BY-NC-SA 4.0.
Companion paper: Biermann et al., *Materials* **18**, 3164 (2025),
[10.3390/ma18133164](https://doi.org/10.3390/ma18133164), open access.

Downloaded and staged under `data/etd2020_biermann2024/`. It is a full factorial:
5 carbomer concentrations (0.1, 0.3, 0.5, 0.7, 0.9 wt%) by 3 triethanolamine doses
(30, 40, 50 µL per 100 g) by 3 replicate stress sweeps, 54 points per sweep, 20 °C.
That gives 15 flow curves over 10⁻³ to 200 s⁻¹ with replicates, so you can see both the
scatter and the trend. That is unusual for deposited rheology.

Two traps in the file format, both handled in `fit_hb.py`:

- The `Shear rate` column is the rounded set point: three significant figures, and only
  two decimals, so every point below 0.005 s⁻¹ is written `0,00`. Taken at face value you
  lose the bottom decade, the one that sets τ_y. The viscosity column is not rounded, so
  the real rate comes back exactly as `γ̇ = τ/η`. The script validates that recovery
  against the set point before using it.
- German locale: semicolon-delimited, decimal comma, latin-1 encoding.

`fit_hb.py` fits `τ = τ_y + K γ̇ⁿ` to all 15 and writes `hb_fits_etd2020.csv`;
`fig_etd2020_flowcurves.png` overlays the curves. Log-log axes, because the rate spans
five decades and a linear axis would collapse the yield plateau onto zero. The 15 fits
fall inside the paper's published ranges (τ_y max 41.1 Pa at 0.5 wt%, k 3.4 to 15.9,
n 0.40 to 0.53) and reproduce its stated maximum. That is a range check, not a
per-condition one: the paper publishes no table of the 15 fitted (τ_y, k, n), only the
ranges in prose and heatmaps (Figs. A1, A2), so nothing stronger is possible without
digitising a figure (checked against the Europe PMC full text, PMC12251366, 2026-09-28;
`biorheo/data/external_anchor_audit_2026-09-08.md` row 6).

The ladder, in Pa:

| wt% \ TEA | 30 µL | 40 µL | 50 µL |
|---|---|---|---|
| 0.1 | 2.6 | 0.8 | 0.8 |
| 0.3 | 19.3 | 32.2 | 39.1 |
| 0.5 | 7.9 | 23.8 | 41.1 |
| 0.7 | 6.7 | 16.5 | 28.4 |
| 0.9 | 8.8 | 16.8 | 26.5 |

Read the rows before the columns. τ_y is not monotonic in carbomer concentration. It
peaks at 0.3 to 0.5 wt% and falls again, because at a fixed neutraliser dose more polymer
means less complete neutralisation, so the gel comes out weaker rather than stronger. The
neutraliser is the control knob, not the concentration. A ladder built by varying wt% at
fixed NaOH will not be monotonic in τ_y, and the crater-depth trend will look like noise.

One caveat limits how far this data can be pushed. The geometry was smooth 50 mm parallel
plates at 0.5 mm gap, not roughened or serrated, so apparent slip is uncontrolled and
every τ_y above is a lower bound. Fine for picking a ladder, not for a quoted
constitutive parameter.

### rheodata (github.com/ManonMarchand/rheodata)

Flow curves digitised from published figures, CSV plus YAML metadata, one file per
figure. The README names Carbopol as a target material, but there is no Carbopol in it
yet: the eight flow curves present are emulsions, PNIPAM microgels and PMMA suspensions
(Mason 1996, Dinkgreve 2015, Seth 2011, Ghosh 2019, Paredes 2013, Petekidis 2004,
Pamvouxoglou 2021). Useful as soft-glass comparison, not as Carbopol.

---

## Tier 2: papers with Herschel-Bulkley parameters tabulated

### Jalaal, Stoeber & Balmforth, *Spreading of viscoplastic droplets*, JFM (2021)

[arXiv:2010.02894](https://arxiv.org/abs/2010.02894). Carbopol Ultrez 21, neutralised
with triethanolamine. Anton Paar MCR-302, sand-blasted PP25-S parallel plates (4 µm
roughness), so slip is actually controlled here. Table 1 gives all five samples plus the
elastic modulus, which almost nobody reports:

| sample | n | K (Pa sⁿ) | τ_y (Pa) | G (Pa) |
|---|---|---|---|---|
| 1 | 0.541 | 0.834 | 1.01 | 9 ± 2 |
| 2 | 0.437 | 2.455 | 5.03 | 33 ± 2 |
| 3 | 0.448 | 4.579 | 10.82 | 57 ± 1 |
| 4 | 0.410 | 8.354 | 16.60 | 64 ± 1 |
| 5 | 0.427 | 10.080 | 25.46 | 66 ± 1 |

Density taken as water's. They did not measure surface tension and used 63 mN/m as the
average of literature values spanning 51 to 69 mN/m. Worth knowing, since crater and
splash arguments are sensitive to it and there is no good measurement to lean on.

### Divoux, Barentin & Manneville, *From transient fluidization to Herschel-Bulkley behaviour*, Soft Matter (2011)

[arXiv:1012.0693](https://arxiv.org/abs/1012.0693). Carbopol ETD 2050 at C = 0.5, 0.7, 1
and 3 wt%, rough plate-plate (glued sandpaper, roughness about 46 to 60 µm), 1 mm gap,
with a deliberate smooth-versus-rough comparison. Yield stresses, from the Fig. 4
caption: 0.5 wt% gives 17.3 Pa, 0.7 wt% gives 30.3 Pa, 1 wt% gives 30 Pa, 3 wt% gives
56.5 Pa. The preparation protocol in section II.A.1 is the most careful one in this list
and is worth copying: dispersion at 50 °C, neutralisation, centrifugation to clear
bubbles, and an explicit statement that microgel size depends on grade, concentration,
final pH, and the stirring speed during neutralisation.

The wider point of the paper matters for a bath too. Under smooth boundary conditions
Carbopol creeps and shear-bands transiently before reaching steady state. Sheared and
then left to sit, "the same fluid" is not the same fluid immediately afterwards.

### Rinkens, Verhoosel, Alicke, Anderson & Jaensson, *Bayesian model selection for complex flows of yield stress fluids* (2026)

[arXiv:2601.10115](https://arxiv.org/abs/2601.10115). Carbopol 980 at 0.5 wt%,
NaOH-neutralised, two independent batches, 50 mm rough cross-hatched plates *and* smooth
plates at three gaps, MCR501, 0.1 to 100 s⁻¹. The full preparation recipe is in section
4.1 and is the closest published match to a 940/980 plus NaOH bath. No data availability
statement in the preprint, so the curves would have to be digitised or requested.

Read its headline before fitting any bath data: HB parameters inferred from rheometry did
not transfer to their squeeze-flow experiment on the same material. A τ_y measured on a
rheometer is not automatically the τ_y that governs a different flow geometry. A crater
is a different flow geometry.

---

## Tier 3: reference points to digitise or cite

- **Roberts & Barnes**, *New measurements of the flow-curves for Carbopol dispersions
  without slip artefacts*, Rheol. Acta **40**, 499 (2001),
  [10.1007/s003970100178](https://doi.org/10.1007/s003970100178). Complete steady-state
  flow curves for Carbopol 980 across 0.045 to 1.0 wt%, using a vane-and-basket geometry
  specifically to kill slip. This is the canonical concentration series for the grade the
  bath experiment plans to use. Paywalled; the figures are the deliverable.
- **Piau**, *Carbopol gels: elastoviscoplastic and slippery glasses made of individual
  swollen sponges*, J. Non-Newt. Fluid Mech. **144**, 1 (2007). The review. Compiles
  rheology across grades and concentrations and explains why slip is the default failure
  mode.
- **Dinkgreve, Fazilati, Denn & Bonn**, *Carbopol: from a simple to a thixotropic yield
  stress fluid*, J. Rheol. **62**, 773 (2018),
  [10.1122/1.5016034](https://doi.org/10.1122/1.5016034). Same grade, different stirring
  history, different physics: prolonged strong stirring fragments the polymer, creates a
  depletion attraction, and turns a simple yield-stress fluid thixotropic. Directly
  relevant, because a bath gets mixed, poured, and then impacted repeatedly.
- **Luu & Forterre**, *Drop impact of yield-stress fluids*, JFM **632**, 301 (2009).
  Already in `Desktop/Papers/Soft Matter & Fluids/Drop Impact/`, with its own HB
  characterisation of the drops.

---

## Two things this scan changes about the plan

**1. The 125 Pa in the proposal has no support in anything found here.** Thrust 3 quotes
0.1 to 1 % w/w giving τ_y of roughly 3 to 125 Pa. At or below 1 wt%, the four independent
studies above report 0.8 to 41 Pa (ETD2020), 1.0 to 25.5 Pa (Ultrez 21), and 17 to 30 Pa
(ETD 2050 at 0.5 to 1 wt%); ETD 2050 needs 3 wt% to reach 56.5 Pa. Nothing here reaches
125 Pa below 1 wt%. Carbopol 940 and 980 are more strongly crosslinked than the ETD and
Ultrez grades and do sit higher at equal weight, so the number is not impossible, but it
needs a source before it is defended in front of a panel. Roberts & Barnes (980, 0.045 to
1.0 wt%) is where to find one.

**2. Build the ladder on neutraliser, not on carbomer.** The ETD2020 factorial says so
directly, and it is the cheaper axis anyway: one stock dispersion, split and titrated to
different pH, rather than five separate batches whose hydration histories differ.

## What is still missing

No open dataset for Carbopol 940 or 980 with NaOH neutralisation, which is what the bath
will actually be. The ETD2020 set is the right *shape* of data on the wrong grade with
the wrong neutraliser. It can size the ladder and check the rheometer, but the bath still
needs its own flow curves, measured on roughened or vane geometry, since every
smooth-plate τ_y in this list is a lower bound.
