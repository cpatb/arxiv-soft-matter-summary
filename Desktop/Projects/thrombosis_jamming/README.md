# thrombosis-jamming — SIPA as a stress-activated jamming transition

A Tier-0 closure toy for **shear-induced platelet aggregation (SIPA) / arterial thrombosis**,
the third fluid in the "one state-variable jamming closure, many fluids" program
(silica drop impact → blood rheology → thrombosis). Built 2026-07-16 (overnight, with Claude)
from a research agent's spec.

## The idea
In high-shear stenotic arteries (~10⁴–10⁵ s⁻¹), clots form by shear-induced, VWF-mediated
platelet aggregation — not the biochemical cascade. Structurally this is a stress-activated
jamming transition: shear unfolds von Willebrand factor → switches on adhesive platelet
contacts → a bonded network percolates → flow arrests. This reuses the **Wyart–Cates skeleton**
with the microphysics relabelled.

## What's novel (honest scope)
Continuum thrombosis models exist (Fogelson; Ku/Mehrabadi two-phase; Antaki VWF-field). They
produce arrest by *deposition kinetics* (smooth, monotonic speed-up with shear). **None treat
occlusion as a bulk rheological jamming divergence.** This toy's distinctive, falsifiable
prediction is a **critical occlusion threshold**: no occlusion below a critical wall shear rate,
and occlusion time *diverging* as `t_occ ∝ (γ̇ − γ̇_crit)^(−α)` above it — a genuine critical
point the incumbent models don't predict. (A 2026 Wyart–Cates closure for epithelial *tissue*,
arXiv 2602.20886, is genre precedent — verify.)

## The closure (see `thrombosis_jamming.py`)
- Activation by VWF unfolding, set by a **critical shear rate** `γ̇_c ≈ 5000 s⁻¹` (Schneider et al. 2007 — verify):
  `f(γ̇) = 1 − exp(−(γ̇/γ̇_c)^β)`, f = fraction of VWF-bridged adhesive contacts.
- Migrating jamming point: `φ_J(f) = φ_m·f + φ0·(1−f)` (φ_m < φ0; adhesion rigidifies looser).
- Effective viscosity diverges: `η_eff = η_plasma·(1 − φ/φ_J)^(−2)`.
- 1-D throat densification → **occlusion time** vs shear rate; diverges at the threshold.

## Result (validated, signature PASS)
`run_thrombosis.py` → `thrombosis_occlusion.png`: occlusion times ~0.4–100 s (physiologically
plausible), no occlusion below γ̇_crit, monotone decreasing + critical divergence above. The
validation checks the *signature* that distinguishes jamming from deposition kinetics.

## Run
```
python run_thrombosis.py     # signature validation + figure
```
Requires numpy, matplotlib.

## First-paper path (wet-lab-free) & risks
- **Falsifiable test:** digitize published occlusion-time-vs-shear (Ku group: *Sci. Rep.* 2020,
  *Biomicrofluidics* 2019) and replot vs (γ̇ − γ̇_crit) — jamming collapses to a power law,
  deposition doesn't. Fix γ̇_c and the modulus scale from *independent* literature so occlusion
  is a genuine prediction.
- **Risks (from the research brief):** prior art closer than it looks (must show fewer params +
  a real critical point, not relabeling); occlusion-vs-shear data may be too sparse to test an
  exponent — *scout the exact figures before writing equations*; real occlusion also involves
  shear-*gradient* effects, platelet contraction, embolization — this is a deliberately minimal model.
- All parameters here are ILLUSTRATIVE. See the vault: *Frontiers/Thrombosis as Jamming — From Analogy to First Paper*.
