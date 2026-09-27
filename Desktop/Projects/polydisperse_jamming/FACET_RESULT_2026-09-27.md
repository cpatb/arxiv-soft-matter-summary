# Convex facets and the shear-jamming window — the pilot was an artefact; what the corrected runs say

2026-09-27. Registration: `PREREG_FACET_2026-09-27.md` (+ Amendments 1–6, all committed before the cells they govern
landed, with every disclosure in place). Grader: `python grade_facet.py` → `out/grade_facet.txt` (every number below
is printed there unless marked). Figure FC1: `python make_facet_figures.py --check` (`FIGURES_FACET.md`).

## 1. What went wrong with the level-set pilot (method, settled)

- **Frozen rotations.** `ls_pack.py` used Yade's `Law2_ScGeom_FrictPhys_CundallStrack()` default `sphericalBodies=True`,
  which computes contact torques from radii: a normal force never torques the grain. Right for a sphere, wrong for any
  other shape (Yade's own `examples/levelSet/levelSetBody.py`: `sphericalBodies=False` "is required in the general
  case"). Two-body test (`scripts/ls_facepair.py`, frictionless ε 0.5 grains face to face, 5° initial tilt): tilt after
  20 000 steps 5.00° (single-point law, default), 5.00° (multi-point law, default), 0.02° (single, False), 0.03° (multi,
  False). So every pilot "rounded cube" was an orientation-frozen grain. The pilot's frictionless Z₀ ≈ 6.3 is the
  translational isostatic count — the signature of frozen rotations, which the pilot had read as "a face pins tilt".
- **How it was caught:** the pilot's frictionless φ₀ FELL by 0.026 at ε 0.5; the verified literature says it RISES
  away from the sphere (Jiao, Stillinger & Torquato, PRE 81, 041304 (2010), doi:10.1103/PhysRevE.81.041304: MRJ
  superballs, φ rises monotonically with |p − 1|, cusp at p = 1, Z near 8 against the isostatic 12; Yuan et al., Soft
  Matter 15 (2019), doi:10.1039/C9SM01932D agree for soft superballs under quasi-static compression). Figure values of
  both papers were not read; only the sign and the text statements are used here.
- **Withdrawn:** "convex facets lower the frictionless edge φ₀ by 0.026 (t 17.5)" in `ANSWER_SHAPE_SIZE_2026-09-25.md`
  §3a, `REGIME_MAP_2026-09-26.md`, `derivations/geometric_window_2026-09-25.md` §7 and `PREREG_LSPILOT` (flags in place).
  Clump packs are unaffected (sphere members; a normal force on a sphere acts through its centre).

## 2. What live rotations exposed: node-discretisation torque noise (method, settled)

The level-set contact point is the deepest surface NODE, laterally off the ideal point by up to about half a node
spacing, so with rotations live a normal force acquires a spurious lever arm (median 0.047 r on sphere contacts in the
census). Consequences, each measured:
- **A live-torque sphere never reaches the stop rule** (unbalanced 0.12–0.13 through 4500 iterations, against 0.004 for
  the torque-free sphere) and is **densified** by the agitation: φ 0.598 vs 0.580 at μ_p 0.3, seed 1, 1000 nodes.
  At 4000 nodes the excess is +0.0056 at a matched ~4000–4500 iterations (vs +0.018), so node spacing drives it; but the
  unconverged live sphere keeps compacting (+0.0099 by 9000 iterations), so the Am.-5 check graded on the log's last
  quarter is time-dependent and sits at its half-way bar (0.009) — read it as "smaller, not gone". The sphere reference therefore keeps the torque-free law, which is exact for a sphere (Am. 3).
- **Frictionless superballs never reach the stop rule** (unbalanced 0.26–0.60, flat from ~6000 to 12 000 iterations;
  worse the closer to a sphere). The frictionless edge φ₀ is NOT measured by this method (Am. 4). Their trajectories
  are reported descriptively only.
- **Frictional superballs (μ_p 0.3) do reach it**, but always by a fluctuation dipping just under the bar (unbalanced
  0.018–0.019 at the stop, after 8k–82k iterations); φ drifts ≤ 0.0007 over 20 000 further iterations, so the reading
  is stable. **μ_p 0.1** sits at unbalanced 0.04–0.16 at 11–12k iterations.
- **Node refinement of the superball values (Am. 5 N1): FAIL at both ε as registered** (|Δφ| bar 0.0042): ε 0.5 gives
  0.6193 / 0.6056 (unconverged) / 0.6299 at 1000 / 2000 / 4000 nodes, ε 0.9 falls 0.6003 → 0.5911 (converged at 4000 nodes). The ε-0.5
  sequence is non-monotone — a different node set is a different trajectory, as the pilot also found — so single-seed
  φ at N 100 carries ~±0.01 from the node set alone, comparable to the sphere's seed scatter (sd 0.009 at μ_p 0.3).

## 3. Results (grade of 2026-09-27 14:40 CDT; cells still running are listed in `out/grade_facet.txt`)

**Registered verdicts.** A-lit, F1–F4 (frictionless edge): **INCOMPLETE by method** — no frictionless superball
reaches the stop rule (§2). F5 (census, registered on N 200 frictionless): INCOMPLETE. F6 (μ_p 0.3 branch), F7 (μ_p 0.1
edge), F8 (window between them): **INCOMPLETE** at every ε — at least one of three seeds had not met the stop rule when
the session ended. N1 (node robustness): **FAIL** at ε 0.5 and 0.9 (§2). Nothing below is a graded verdict.

**What the landed data do support (descriptive; single N 100, 1–3 seeds; sphere = torque-free reference, seeds 1–3):**

| | sphere (ε 1) | ε 0.9 (p 1.11) | ε 0.8 (p 1.25) | ε 0.7 (p 1.43) | ε 0.5 (p 2) |
|---|---|---|---|---|---|
| φ at μ_p 0.3, stop rule, 1000 nodes | 0.5827 ± 0.0090 (3) | 0.6003 (1) | 0.6130 ± 0.0107 (2) | — | 0.6214 ± 0.0029 (2) |
| same, 4000 nodes, seed 1 | (live torques: 0.587 u, noise) | 0.5911 (stop rule) | — | — | 0.6299 (stop rule) |
| φ at μ_p 0.1 (u = unconverged, log readout) | 0.6167 / 0.6173 / 0.6244 | 0.638–0.644 u | 0.640–0.657 u | 0.658–0.661 u | 0.680–0.687 u |
| W₀.₁ = 1 − φ(0.3)/φ(0.1) | 0.059 (0.051–0.068) | 0.060 (1 seed, u) | 0.055, 0.067 (u) | — | 0.092, 0.083 (u) |
| frictionless, seed 1, last log (u) | 0.6387 ± 0.0025 (stop rule) | 0.665 | 0.683 | 0.700 | 0.736 |
| face–face enrichment E_ff at μ_p 0.3 | 0.56 (N 200, 1 seed) | 2.01 | 1.74 | — | 1.14 (N 100), 1.00 (N 200) |

1. **Convex faceting RAISES the jamming density on every branch that could be read — the opposite of the withdrawn
   pilot.** At ε 0.5 the frictional φ_m is above the sphere on every seed and every node count tried: seed-paired
   +0.040 / +0.048 (1000 nodes), +0.050 at 4000 nodes (converged), and even the 2000-node unconverged value (0.6058) is
   +0.026 above its seed's sphere — larger than the whole node-noise densification of a sphere (+0.018 at 1000 nodes).
   The low-friction edge rises more (+0.06 at μ_p 0.1, upper-biased because unconverged cells keep creeping upward under
   the noise), and every frictionless trajectory sits far above the sphere, in the direction Jiao et al. report.
2. **Near the sphere the effect is not resolved.** At ε 0.9 the frictional +0.02 shrinks to +0.011 at 4000 nodes (0.5911 vs the seed's 0.5797 sphere),
   the size of the residual node noise, and W₀.₁ at ε 0.9 and 0.8 equals the sphere's within scatter.
3. **The window between μ_p 0.1 and 0.3 is ~1.5× the sphere's at ε 0.5 (0.083–0.092 vs 0.051–0.068), on two seeds
   with the μ_p 0.1 edge unconverged.** That is suggestive, not established: both superball edges carry node noise and
   the direction of its bias on a RATIO of the two is not measured. The frictionless window (the registered question)
   is not measured at all.
4. **The frictional rise is a counting effect, not a filling one (post hoc).** In φ = Z/(Z + 2√3 w) with Z = z_nr:
   sphere Z 4.89, w 1.01; ε 0.8 Z 5.53, w 1.01; ε 0.5 Z 5.55, w 0.98 — faceted grains jam frictionally with ~0.6 more
   load-bearing contacts at unchanged filling. (z_nr is read at the end of the post-relaxation, not at the stop.)
5. **Faces do not pair up at ε 0.5.** The fraction of contacts that are face-on for both grains matches the
   uniform-area null (E_ff 1.00–1.14; P(m > 0.95) 0.55–0.58 vs 0.567). Weakly faceted grains (ε 0.8–0.9) DO prefer
   their flattest patches (E_ff 1.7–2.0). The sphere's E_ff 0.56 at μ_p 0.3 is one N-200 cell (6 face–face contacts
   expected ~11); treat as scatter.

## 4. What this does to the author's question 2 (sphericity and cornstarch)

- The only mechanism this project had found that pointed the other way — "convex facets lower the frictionless edge
  and may narrow the window" — was a frozen-rotation artefact. **No shape family simulated in this project now
  narrows the window**: elongation and concavity (clumps, doublets) widen it ×1.4–1.7 (registered, replicated), and
  convex faceting raises both edges, the low-friction one more (descriptive, ε 0.5 only).
- On geometry alone the sphere is therefore the NARROW-window extreme among the shapes tested, and Stöber silica's
  sphericity is load-bearing in that sense. A faceted and irregular grain like cornstarch is pushed wider by both
  features, consistent in sign with the published χ 0.76 (cornstarch) vs 0.83 (1.5 µm silica) — which still mixes lab,
  size, friction and φ convention, so it is not a test.
- Unchanged from `ANSWER_SHAPE_SIZE` §3c: whether a jammed cornstarch drop looks like a silica BALL is set by drainage
  and contact stiffness (De, S), not by geometry.

## 5. What would close it
1. A frictionless measurement needs a contact that is torque-exact for convex shapes in a periodic cell: either many
   more surface nodes (the sphere's noise excess fell from 0.018 to 0.0056 at matched iterations going 1000 → 4000 nodes, and keeps creeping;
   cost ∝ nodes),
   or a polyhedral/superquadric code with analytic contact (not available in this Yade build). Until then φ₀(ε) rests
   on Jiao 2010 / Yuan 2019 (sign only).
2. The running μ_p 0.1 / 0.3 cells can be graded when they land: `wsl -e bash ~/sync_facet.sh` then
   `python grade_facet.py`; F6–F8 grade automatically once all three seeds of a cell have met the stop rule.
3. For the bench: a lower χ (wider window) for fractionated cornstarch than for Stöber silica at matched CV is now the
   direction BOTH shape mechanisms predict; the facet-specific signature the withdrawn pilot proposed (a lower
   frictionless φ₀) is withdrawn with it.
