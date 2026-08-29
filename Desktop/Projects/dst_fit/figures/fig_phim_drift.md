# fig_phim_drift

Corn starch, six volume fractions, TA flow sweeps from August 2023. Regenerate with
`python scripts/make_fig_phim_drift.py`.

## What each panel shows

**(a)** Measured viscosity against stress, open circles, one colour per nominal volume
fraction. Solid lines are Wyart-Cates fits. Each line covers only the stress range it
was fitted over, which is the thickening branch, from the viscosity minimum upward.

The gap between the circles and the lines at low stress is not a bad fit. Wyart-Cates
has no thinning branch at all: its viscosity increases monotonically with stress. Corn
starch thins by a decade or more before it thickens, so that part of every curve is
outside what the model can describe and is excluded rather than fitted. Where the cut
falls changes the recovered `sigma*`, so the fit records it.

**(b)** The `phi_m` recovered from each of those six fits, against the nominal volume
fraction. Grey horizontal line is the inverse-variance weighted mean, 0.5697. Dashed
line is `phi_m = phi_nom`, for reference only.

`phi_m` is the frictional jamming fraction. It is a property of the particles, so six
fits of the same corn starch have to return the same number. It climbs from 0.492 to
0.580.

**(c)** The same six values in units of their own standard errors, deviation from the
weighted mean. Grey band is plus or minus two sigma.

## Why (b) and (c) both exist

Panel (b) makes the drift look minor, and in absolute terms it is: 0.09 on a quantity
bounded near 0.64. The error bars in (b) are drawn but invisible, because they run from
0.0017 to 0.0099 on an axis spanning 0.68.

That invisibility is the point. Panel (c) rescales by those error bars and the same
six points span minus 7.8 to plus 6.2 sigma, with none inside the two-sigma band.
Reduced chi-squared for a single common `phi_m` is 32.4 on 5 degrees of freedom.

Small and insignificant are different things, and a figure showing only (b) would
suggest agreement while a figure showing only (c) would suggest a large effect. Neither
alone is honest.

## Axis choices

Panel (a) is log-log. The stress axis spans more than five decades and the viscosity
axis more than three, so linear axes would collapse everything below 10 Pa into the
origin.

Panels (b) and (c) are linear and start at zero rather than cropping to the data range.
Cropping (b) to 0.48 to 0.59 would fill the panel with the drift and make a 0.09 change
look like the whole story. The zero-based axis is what makes (c) necessary, which is the
correct order: show the real magnitude, then show the significance separately.

## What the figure does not settle, and what has since been settled elsewhere

**This figure now understates its own result, and should be redrawn before it is
published.** It shows one corn-starch series, and reads as a statement about corn starch.
It is not.

Two published series were run through the same pipeline on 2026-08-27:

- **Richards et al. 2019 corn starch**, a different lab, solvent, rheometer and operator,
  drifts harder: reduced chi-squared **221 on 6 dof** against the 32.4 shown here.
- **Guy et al. 2015 PMMA** — monodisperse, non-porous, non-swelling, density matched, the
  system Wyart-Cates was built for — fails the same test: **1 of 6** within two sigma,
  reduced chi-squared **344**.

In both cases the weighted-mean `phi_m` lands within 0.016 of the authors' own published
value, so panel (b)'s *level* is right and only its *spread* fails.

The mechanism is identifiability, not swelling. `phi` and `phi_m` enter Wyart-Cates only
through `phi/phi_J`, so a curve approaching jamming pins the ratio and says nothing about
either number: `r(phi, phi_m)` runs to +0.99 and `phi_m - phi` collapses to 0.0013. Pinning
`phi` at an independently known volume fraction takes pairwise agreement from 1 of 15 to
11 of 15 — and makes the residuals three times worse, which is the real finding.

`phi_0` is dead as an explanation for a second time: moving it changes the level of `phi_m`
and leaves the reduced chi-squared bit-for-bit identical, on both materials. Granule
swelling cannot be the explanation either, since PMMA does not swell.

**Redraw brief:** the honest version of this figure is panel (c) repeated for three
datasets side by side, with the published `phi_m` marked on each. See the repository
README, gates N4, N5, N6 and G4.
