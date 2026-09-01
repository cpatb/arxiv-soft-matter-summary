# polydisperse_jamming

A hydrodynamics-free DEM toy model asking whether the **width of a particle size
distribution** can move an impacting dense suspension between a jammed, rebounding solid
and an inertial-dilatant hybrid, and, if so, through which mechanism.

**Read `SPEC_2026-08-31.md` first.** It carries the physics, the pre-registered tests, the
anchors, and a §11 log of every design change made during execution. This file is only the
map.

---

## The question in one paragraph

The author observes on the bench that at *similar particle size*, raising the coefficient
of variation `CV` moves an impacting drop from purely dilatant to spreading and
inertial-dominated. "Similar size" means matched **number mean**. But Carman–Kozeny
permeability reads the **Sauter mean** `d₃₂ = ⟨d³⟩/⟨d²⟩`, and for a lognormal at fixed
number mean `d₃₂/d₁₀ = (1+CV²)²` exactly, so raising CV raises permeability, and the
lab's own inertial–dilatant number `ID ∝ d₃₂²` with it, *in the observed direction with no
packing physics at all*. This project exists to separate that null (M1) from two
structural candidates (M2, M3).

## The three mechanisms

| | mechanism | where it lives | status |
|---|---|---|---|
| **M1** | permeability via `d₃₂` | the interstitial fluid | **excluded** over most of the range: `ID` spans 7.2e-5 to 1.8e-3 against a boundary at 1 (SPEC §1.2) |
| **M2** | distance to jamming, via `φ_J(CV)` | the packing | Phase 2 |
| **M3** | force-network structure at matched mechanical state | the packing | Phase 3 |

**A DEM packing is scale-invariant, so it cannot see `d₃₂` at all**, which means every
result here is free of M1 by construction, not by matching (SPEC §11.7).

## Layout

| file | what it does |
|---|---|
| `psd.py` | size distributions and their moments. `LogNormalPSD`, and `UniformPSD` (skewness ≡ 0) to separate CV from skewness |
| `packing.py` | Yade packing generation + isotropic stress-controlled compression |
| `measure.py` | pure functions of a contact graph: rattlers, backbone, fabric, force stats |
| `phase2_phiJ.py` | the `φ_J(CV, μ_p)` sweep, resumable, shardable |
| `phase4_impact.py` | a carved aggregate in a free / capillary / caged boundary, driven into a wall |
| `analyse_phase2.py` | φ_J extrapolation and the M2 verdict |
| `analyse_phase3.py` | structure at matched mechanical state, the M3 verdict |
| `analyse_phase4.py` | boundary-condition comparison |
| `make_figures.py` | every figure, with a byte-identical regeneration check |
| `test_psd.py`, `test_measure.py` | 30 engine-independent anchor tests |
| `smoke_yade.py` | proves Yade *simulates*, not just imports |
| `diag_ladder.py` | prints the raw φ(p̂) ladders, for checking the extrapolation model |

## Running it

Yade is installed **userspace** (no root on this machine). `scripts/env.sh` documents the
three things that each silently broke the import, read it before debugging anything.

```bash
# engine-independent tests, from Windows          -> 30 passed
python -m pytest test_psd.py test_measure.py -q

# Yade actually simulates                          -> SMOKE_PASS
wsl bash -lc '~/rl.sh smoke_yade.py'

# the phi_J sweep, 12 single-threaded workers, resumable
wsl bash -lc '~/lp2.sh 12'

# everything remaining, chained and unsupervised
wsl bash -lc '~/rr.sh'

# analysis + figures + sync results back to Windows
wsl bash -lc '~/aa.sh'
```

Code flows Windows → WSL (`~/pj`); results flow back to `out/` and `figures/`.

## Things that will bite you

- **The rattler threshold is a definition that changes the answer.** Iterated removal of
  grains with `< k` contacts is k-core decomposition, and at frictional mean `Z` ≈ 4.9 the
  **k = 4 core is empty**, the whole pack reads as rattlers, while one non-iterated pass
  removes 15%. All of `k` = 2/3/4 are recorded; `k` = 2 is primary.
- **`N` is bounded below by the periodic collider**, and the bound binds at the *final*
  density, not the initial one. At `CV` = 0.5 that means `N ≥ 1299`.
- **The CV ladder stops at 0.5** and everything is reported by *realised* CV: the
  truncation needed to keep grains inside the cell puts a nominal 0.8 at a realised 0.72
  with `d₃₂` 21% low.
- **Yade body forces must come from a `PyRunner`.** `ForceResetter` zeroes `O.forces` every
  timestep, so a force added from Python between `O.run()` calls acts once in N steps. This
  is what made the first Phase-4 version return COR = 2.6e37.
- **Many single-threaded workers beat few multi-threaded ones.** Idle OpenMP threads spin;
  8 workers × 2 threads gave load 42 on 16 physical cores.
- **`si_fit` is withdrawn**, do not calibrate anything here against the lab's silica
  closure parameters.

## Provenance and limits

- The CV dataset motivating this is **not on disk**; every CV is a bench-note figure.
- Phase 4 runs at a granular inertial number `I ~ d/D_agg ~ 0.08` against `4e-4` for a real
  drop, **~200× too inertial per grain**. It compares boundary conditions and predicts no
  experimental regime boundary.
- Literature anchors: Silbert, *Soft Matter* **6**, 2918 (2010), arXiv:1108.0012 (φ_J falls
  0.64 → 0.55 with friction, `Z` 6 → 4); Desmond & Weeks, *PRE* **90**, 022204 (2014),
  arXiv:1303.4627 (skewness can outweigh polydispersity for φ_RCP).

State card: `Documents\Obsidian Vault\Status\Drop Impact — Polydispersity Jamming Model.md`
