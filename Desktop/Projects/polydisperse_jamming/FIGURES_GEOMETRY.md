# Figures — contact geometry of σ\* and of rattling (2026-09-24)

`python make_geometry_figures.py` → `figures/G1_sigma_star_geometry.{pdf,png,csv}`.

## G1
**(a)** The CV → σ\* effect at fixed number-mean diameter, split into the exact factors of the Love–Weber
identity (`derivations/sigma_star_geometry_2026-09-24.md` eq. 7): contact density Zφ, the contact-pair size
product ⟨aᵢaⱼ⟩_c/⟨a³⟩, the force–size correlation C, and the inverse tail-shape factor 1/q; black is the
measured σ\* ratio (their product, to ≤ 1.1 %). Sheared constant-pressure Yade packs, μ_p 0.3, seed-paired
means (n = 3), at p̂ 3e-4 and 1e-3 (the 1e-4 cells unjam intermittently and are left out). Read: only Zφ
moves far from 1, and C and 1/q push the other way.

**(b)** Per-grain probability of being a rattler (iterated, < 2 contacts) against size, in eight quantile
bins, for two isotropic packs (μ_p 0.3, p̂ 3e-5, seed 1; CV 0.2 blue, 0.5 orange). Lines: Poisson prediction
from the pack's own mean contact number with a per-grain mean ∝ a² (solid) or ∝ the closure-C form
a(a⟨a⟩+⟨a²⟩) (dotted); neither has a fitted parameter. The a² form is the better one in both packs.
`derivations/rattler_closure_2026-09-24.md` §3b.

The CSV header carries a generation timestamp (house `save()`), so it is not byte-identical across runs;
the PDF/PNG and the data columns are.

## G2 (`python make_geometry_figures_g2.py`)
Out-of-sample tests on binary (135 records, size ratio 1.25–5) and ternary 1 : 2 : 4 (36 records) isotropic
packs, μ_p 0.3 (`derivations/binary_contact_check.py`; rattler_closure §4d).
**(a)** Per-class rattler fraction predicted from the pack's own mean contact number with a per-grain Poisson
mean ∝ a² (coloured; no fitted parameter) against the measured one; grey is the size-blind Poisson(Z)
prediction for the same classes. Identity line. The ternary points fall below the line (rattling
under-predicted by ~0.08).
**(b)** Number of small–large contacts divided by the configuration-model value D_S D_L/(D_S+D_L) from the
measured class degrees, one point per binary record, against size ratio. 1 = random pairing given degrees.
