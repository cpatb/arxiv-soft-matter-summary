# Pre-registration: a second protocol on silica — constant-stress start-up

**Registered 2026-09-27 ~23:55 CDT, before any clock was fitted to these records.** Seen
before registering: the raw rate at the first/last samples of 4 files (inertial rate peak
0.22-0.99 1/s at gamma ~0.05; final rates 0.065-0.107 1/s, matching the steady curve).

## Data

Richards et al. 2019 deposit (DataShare 10283/3357), Figure3/Data/silica87glyc: 4 um
silica in 87 wt% glycerol-water, phi = 0.574, 17 constant-stress steps 1-11 Pa, 20 s each,
same direction, separated by ~14 s at zero stress (during which the sample recoils ~0.06
in strain). Steady curve `gw_phi_0p574_flow.txt` (relative viscosity given, so
eta_s = eta/eta_r from the same file). Stress-controlled, where Lin is rate-controlled.

## Method (fixed now)

eta_r(gamma) = eta/eta_s. Start of fit = the maximum of the rate within gamma < 0.3 (end of
rotor spin-up; analogue of the reversal minimum). From there, the SAME package pipeline as
the reversal clock: `log_bin`, u = 1/(1 - eta_r^(-1/n)) (n = 2), `fit_relaxation`
(affine-invariant, so no packing triple is needed). Report per stress: gamma_eff, profile
interval, chi2_r, and the viscosity rise eta_end/eta_start. Onset of thickening from the
steady curve: the lowest stress where steady eta_r exceeds 2x its low-stress plateau.

## Predictions

- **P1.** At stresses above onset, the start-up clock is below 0.30 on every step, i.e. at
  least 1.8x below the smallest reversal clock on Lin's silica (0.54). Fails if any
  above-onset step gives >= 0.30.
- **P2 (the control that can embarrass P1).** Below onset, where Wyart-Cates has almost no
  frictional fraction to rebuild, a transient still exists (seen in the raw rate at 1 Pa).
  If its clock is within 2x of the above-onset median, the start-up transient is NOT
  primarily the f-relaxation of Eq. 1 and P1 cannot be read as a measurement of gamma_0;
  it is then a statement about the microstructure. Reported either way.
- **P3 (sanity).** The final viscosity of each step agrees with the steady curve at that
  stress to within 30 % (the steady file's own replicate spread reaches ~20 %).

## Addendum 1 (~23:40; git commit da49f2e time is authoritative; after P1-P3 were graded, before any recoil strain was computed)

Graded result, `logs/probe_startup_clock.txt`: **P1 FAIL** (above-onset clocks 0.027-0.885,
chi2_r 73-11343, 2 fits at the grid edge; neighbouring stresses differ up to 30x -- a
single strain exponential does not describe these records, so no clock is measured).
P2 reported: medians 0.045 below vs 0.159 above onset. **P3 FAIL** at 3.0 and 9.5 Pa.

Plotting the records (not a fit) shows why: above onset the rate COLLAPSES abruptly at a
nearly stress-independent strain gamma ~0.07-0.09, and below onset eta rises around the same
strain. The rest before each step recoils the sample by ~0.06. Hypothesis M (strain memory):
the network re-forms once the restart has undone the recoil.

- **P4 (fit-free).** Define gamma_j per step = the strain at the steepest fall of the rate
  (max -d(rate)/d(gamma) on the raw record, gamma < 0.3), and R = |net strain of the
  preceding Relax file| (final minus initial). M predicts gamma_j - R is small and
  stress-independent: |gamma_j - R| < 0.03 on at least 80 % of above-onset steps, and
  gamma_j tracks R across steps (Spearman rho > 0.5). Competing reading: gamma_j is a
  constant (~0.08) unrelated to R, i.e. rho ~ 0 -- then the jam strain is a material
  strain scale, not memory. Both are reported; the first step (0_Stress_1, no preceding
  Relax) is excluded by construction.

## What this cannot do

Different particle size (4 vs 2 um), phi and solvent from Lin's silica, so it is the same
material CLASS, not the same sample. Rest periods of ~14 s with recoil mean each step starts
from a partly reversed structure, not from rest. Both stated wherever the result is quoted.
