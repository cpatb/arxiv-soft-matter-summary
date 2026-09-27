# Figure FC1 — convex facets, rotations live (`make_facet_figures.py`, `figures/FC1_facet_window.{pdf,png,csv}`)

Regenerate: `python make_facet_figures.py --check` (PNG regenerates byte-identically). Inputs: `out/lsdem/ls_facet2*.jsonl`
and the progress logs `out/lsdem/facet2mu01_*.log`, `facet2vis_*.log` (synced by `scripts/sync_facet.sh`), plus the
pilot's 1000-node sphere records (torque-free law, valid for a sphere; `PREREG_FACET` Am. 3). Level-set superellipsoids
|x|^(2/ε)+|y|^(2/ε)+|z|^(2/ε) ≤ 1, volume-equivalent to d = 1, N 100, p̂ 1e-3, 1000 surface nodes. The x axis is the
superball exponent p = 1/ε (1 sphere, 2 a rounded cube). The CSV carries every plotted point with its seed and a
`converged` flag. Registration and verdicts: `PREREG_FACET_2026-09-27.md`, `out/grade_facet.txt`,
`FACET_RESULT_2026-09-27.md`.

**(a) Packing fraction at the target pressure.** Orange: μ_p 0.3. Blue: μ_p 0.1. One point per seed. Filled = met the
stop rule (unbalanced force < 0.02); open = never met it, read as the mean over the last quarter of its progress log
(`PREREG_FACET` Am. 6, descriptive; these creep upward while the run continues). Dark diamonds: frictionless — the
sphere's stop-rule value (filled) and the seed-1 superball trajectories at their last log (open; they never converge,
Am. 4). Grey crosses: the WITHDRAWN frozen-rotation pilot (frictionless), drawn to show how large the artefact was.
Linear axis from zero.

**(b) W₀.₁ = 1 − φ(μ_p 0.3)/φ(μ_p 0.1) per seed**, where both edges exist; open if either edge is unconverged. This is
the window between μ_p 0.1 and 0.3, not the registered frictionless window, which this method cannot measure.

**(c) Face–face contact enrichment at μ_p 0.3**: the fraction of contacts that are face-on (m = max|nᵢ| > 0.95 in the
body frame) for BOTH grains, divided by its uniform-surface-area null (`facet_census.py`). 1 = no preference. The sphere
point is one N-200 cell with ~6 face–face contacts, i.e. scatter.

⚠ Node-discretisation torque noise densifies every rotating-grain pack (a live-torque sphere: +0.018 at 1000 nodes,
+0.0099 at 4000 on the registered readout, which did not pass its halving bar); the superball points carry it and the sphere reference does not. At ε 0.5 the frictional gap to
the sphere (+0.04–0.05) exceeds it; at ε 0.9 it does not.
