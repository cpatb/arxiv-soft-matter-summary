# glymphatic-transport — advection-vs-diffusion regime map

A transparent **dimensional-analysis** artifact for the central controversy in brain-fluid
clearance: does solute transport go by bulk flow (advection) or diffusion? Built 2026-07-16
(overnight, with Claude) to operationalize the "Fluid Dynamics & the Brain" frontier note.

**This is not a simulation** — it is a Péclet-number regime map (`Pe = U·L/D`) from
literature-cited numbers, showing where each regime dominates.

## Result (validated against the field synthesis)
`glymphatic_peclet.py` → `glymphatic_peclet.png`. Crossover lengths (Pe=1):
- Perivascular, small molecule: L* ≈ 33 µm; amyloid-β: ≈ 7 µm → **advective** at PVS-segment scale (~250 µm, Pe≈7.5).
- Deep parenchyma (low Darcy velocity, tortuosity-corrected D): **diffusive** at ECS scale (~20 µm, Pe≈0.05).

This reproduces the current synthesis (Kelley & Thomas ARFM 2023; Holter PNAS 2017): **advection in
the fast perivascular spaces, diffusion deep in the low-permeability parenchyma.** The map makes the
crossover quantitative and solute-size-dependent.

## Inputs (literature-cited — verify before formal use)
- Perivascular CSF speed ~10–20 µm/s (Mestre et al., Nat. Commun. 2018)
- Parenchymal Darcy speed ~0.1–1 µm/s (Holter et al., PNAS 2017)
- D: small molecule ~5e-10, amyloid-β monomer ~1e-10 m²/s; tissue D_eff = D/λ², λ≈1.6 (Nicholson)
- Oscillatory (Taylor–Aris) dispersion can raise effective D tens-of-% to several-fold (Asgari 2016) — an uncertainty band, not modeled precisely here.

## Run
```
python glymphatic_peclet.py   # regime-map validation + figure
```
Requires numpy, matplotlib.

## Where it goes (from the frontier note)
The natural next builds — and the ones the field is actually bottlenecked on — are (1) a 1-D
pulsatile-annulus **Taylor–Aris dispersion** calculator for realistic (eccentric) PVS cross-sections,
and (2) a **PTV/velocimetry + uncertainty pipeline** for perivascular imaging (your differentiating
skill). See vault: *Frontiers/Fluid Dynamics & the Brain* (best entry: Kelley/Rochester, Tithof/Minnesota).
