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

## What this cannot do

Different particle size (4 vs 2 um), phi and solvent from Lin's silica, so it is the same
material CLASS, not the same sample. Rest periods of ~14 s with recoil mean each step starts
from a partly reversed structure, not from rest. Both stated wherever the result is quoted.
