# Figure M1 — the ball-vs-hybrid plane (`make_regime_figures.py`, `figures/M1_ball_hybrid_plane.{pdf,png,csv}`)

Regenerate: `python make_regime_figures.py --check` (PNG byte-identical). A SCALING ESTIMATE on the project's poroelastic
model (`derivations/ball_hybrid_scaling_2026-09-25.md`, `impact_stress/pore_relaxation.py`), not a measurement.

Axes: dilatant-strength ratio S = K_pΔφ/(φρU₀²) (how much dilatant suction the contact network can carry, relative to
the impact stress) against drainage Deborah number De = τ_pore/τ_impact (≫ 1: undrained). Grey dashed: S = 1, De = 1.
- Blue circles: 800 nm Stöber silica, φ 0.44, lab impact (D₀ 2.5 mm, U₀ 2.3 m/s), across the lab salt range
  λ_D = 30 → 10 → 3 → 1 nm (K_p = σ* a/6λ_D, σ* held fixed; with the bench's onset shift folded in, the path is shorter
  — `derivations/potential_ladder_2026-09-26.txt`).
- Pink squares: silica at fixed chemistry, particle size 0.4 → 20 µm (Derjaguin σ* ∝ a⁻¹, so K_p and S are size-free and
  only De moves, ∝ a⁻²).
- Orange triangles up: 10 µm-class cornstarch, φ 0.42, [M26] impact (D₀ 4 mm, U₀ 5.1 m/s), size fractions 5 → 20 µm
  (Hertz K_p recomputed; size-free at fixed confinement).
- Green triangles down: 10 µm cornstarch, solvent viscosity 10⁻³ → 1 Pa s (De ∝ μ_f; S unchanged).

Reading: the materials are separated HORIZONTALLY — by the contact network's modulus (double layer vs Hertzian contact,
×1870 in K_p) — and every size, salt, fractionation or viscosity knob moves them VERTICALLY. None crosses S = 1.
⚠ The silica K_p rule ties P_conf to σ* (silica_poro's declared choice, the most model-laden line); σ* for silica is a
literature bracket, the repo value is a retracted fit kept for comparability; cornstarch's particle-modulus bracket
(E_p 10⁷–10⁹ Pa) spans S = 0.89–19 (`cornstarch_fraction_2026-09-25.txt`), so at its soft end cornstarch sits just
below S = 1 — the plotted default (4.1) is not a sharp position, and S ~ 1 is a crossover, not a boundary. Silica's
S stays ≤ 0.04 across every knob plotted: the separation between the materials is 1.4 decades even at the closest ends.
