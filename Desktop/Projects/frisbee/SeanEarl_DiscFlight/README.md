# Sean Earl Disc Flight

Project 1 codebase: coupling analytical disc aerodynamics with site-specific CFD
over Loyola University's Sean Earl Field.

## Repository layout

```
SeanEarl_DiscFlight/
├── P1_AnalyticalModel/      6-DOF disc trajectory integrator (Python)
├── P2_Geometry/             Sean Earl + campus block STL pipeline
├── P3_CFD/                  OpenFOAM case template + 13-direction sweep
├── P4_Coupling/             CFD-to-trajectory coupling, Monte Carlo, gust gen
└── P5_Synthesis/            Throwability metric, figures
```

## Phase status

| Phase | What it produces | Status |
|---|---|---|
| P0 | OpenFOAM install, environment | **You.** |
| P1 | `disc_6dof.py` validated against Hummel (2003) | **Drafted + sandbox-tested.** |
| P2 | `sean_earl_buildings.stl`, ground DEM, tree porous-zone | Scaffold + workflow doc. |
| P3 | 13 RANS cases, sampled flow on 0.5 m grid 1–6 m AGL | Case template + generator. |
| P4 | Monte-Carlo landing distributions per (x, θ_wind, throw) | Coupling module + MC driver. |
| P5 | `Q(x, θ_wind)` maps, paper figures | Metric calculator + plot stub. |

Each subdirectory has its own `README.md` with running instructions and a
checklist of remaining decisions.

## Dependencies

- Python 3.10+ with numpy, scipy, matplotlib (P1, P4, P5)
- pyvista (P2 STL handling — optional, can use trimesh)
- OpenFOAM v2406 or 11 (P3)
- snappyHexMesh and ParaView (P3 post-processing)
- requests + rasterio (P2 data fetching — optional)

## How to test P1 immediately

```bash
cd P1_AnalyticalModel
python -m pip install numpy scipy matplotlib
python validation.py                  # reproduces Hummel-style trim glide
python -m pytest test_disc_6dof.py    # 4 sanity tests
```

## How to work the rest

P2 → P3 → P4 → P5 in order; each phase's outputs feed the next. See each
phase README for prerequisites and outputs.
