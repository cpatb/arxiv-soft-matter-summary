# F1 — the jamming packing fraction against polydispersity

**What is plotted.** `φ_J`, the packing fraction at which the assembly first supports
stress, against the *realised* coefficient of variation of the particle diameter
distribution. One series per interparticle friction coefficient `μ_p`. Error bars are the
standard deviation over three independent seeds. The dotted horizontal line is the
published frictionless random-close-packing value, 0.64.

**How `φ_J` was obtained.** Never from a single run. For each `(family, CV, μ_p, seed)` a
ladder of reduced pressures `p̂ = p/E` was compressed and `φ(p̂)` extrapolated to `p̂ → 0`.
The linear extrapolation is what a harmonic contact law predicts; the free-exponent fit is
carried alongside as a check and is reported in `out/phase2_verdict.txt`.

**Why the y-axis starts at zero.** House rule — a linear axis begins at zero so the eye
reads the *size* of the effect, not just its ordering. The CV-induced shift is genuinely
small on that scale (a few percent of `φ_J`), and the figure should say so rather than
magnify it by cropping.

**What it shows.** `φ_J` rises monotonically with CV at both friction values, and the two
`μ_p` curves are separated by far more than the CV effect. The friction axis is the larger
one.

**Anchors this figure passes.**
- The `μ_p = 0`, CV = 0 point must sit at 0.64 ± 0.02 (published RCP).
- The ordering `φ_J(μ_p = 0) > φ_J(μ_p = 0.3)` reproduces Silbert, *Soft Matter* **6**,
  2918 (2010), who find 0.64 → 0.55 over decades of friction.

**Caveats.** 19% of cells did not reach the unbalanced-force criterion within the
iteration cap, and the failures correlate with CV (5% at CV = 0, 46% at CV = 0.5).
Dropping them shifts `φ` by at most 0.004 and by a median of 0.0006, so the selection is a
compute-budget artefact rather than a bias — but the check is printed by
`analyse_phase2.py` on every run and should be read, not assumed.
