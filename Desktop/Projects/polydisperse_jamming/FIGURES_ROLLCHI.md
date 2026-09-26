# Figure R1 — the geometric window (`make_rollchi_figures.py`, `figures/R1_geometric_window.{pdf,png,csv}`)

Regenerate: `python make_rollchi_figures.py --check` (PNG regenerates byte-identically). Inputs: `out/rollchi/*.jsonl`
(PREREG_ROLLCHI, synced with `scripts/sync_rollchi.sh`) and `out/shape/phase2_{kmer,shapechi}_*.jsonl`. N 1700, CV 0
unless open symbols; seed means, bars = seed standard deviation; clump cells pass the K4 overlap gate; kmer μ_p 3
multiplet cells are excluded (K4-failed in `grade_kmer.txt`). The figure is regenerated from whatever records have
landed, so its n per point is in the CSV.

**(a) φ against the non-rattler coordination Z at p̂ 3e-5.** Black circles: spheres across friction μ_p 0–3. Blue
squares: spheres with a contact moment limit η > 0 (filled CV 0, open CV 0.5). Green triangles up: two-sphere clumps
(α 1.1–2.0). Pink triangles down: singlet/doublet mixtures (kmer). Grey curves: φ = Z/(Z + 2√3 w) — solid w = 1
(Song, Wang & Makse 2008, no free parameter), dotted w of pure doublets at μ_p 0 (printed on the panel). Friction and
rolling resistance slide a population ALONG its curve (they change Z); size distribution and a closing neck move it
to a different curve (they change w). ⚠ SWM state the relation for Z ∈ [4, 6]; points below Z = 4 (rolling) are an
extrapolation that the data happen to follow.

**(b) The frictional jamming fraction φ_m(α) at μ_p 0.3.** Orange circles p̂ 1e-3, green triangles p̂ 3e-5. Dashed
lines: spheres with a saturated moment limit (η 1) at the same rung — the level pure doublets also reach. The curve
rises for mild asphericity (more neighbours at unchanged filling) and falls continuously as the neck closes toward
α = 2 (PREREG_ROLLCHI N1: GRADUAL). Some interior points have n = 1–2 (CSV).

**(c) Window width 1 − χ = 1 − φ(μ_p 0.3)/φ(μ_p 0) against α at p̂ 1e-3** — the rung at which both branches exist for
every α (the frictionless clump cells at 3e-5 did not finish in the session). Spheres 0.072; every clump 0.10–0.11.
The window opens by α 1.1 and then stays flat. ⚠ Linear axis from zero; α 1.2 and 1.5 rest on single seeds.
