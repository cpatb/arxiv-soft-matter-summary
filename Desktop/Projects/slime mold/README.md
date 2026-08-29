# Growing the Chicago 'L' with Slime-Mold Physics

A constraint-aware, demand-weighted model that uses the *Physarum polycephalum*
(slime mold) network-formation physics of Tero et al. (Science, 2010) to **grow**
a candidate expansion of the Chicago Transit Authority 'L' rail system, benchmarked
against classical optimization baselines.

## What it does
- Builds a realistic Chicago candidate graph (678 nodes, 1,937 edges) on a 1.6 km
  triangular lattice over the buildable land area, with **Lake Michigan as a hard
  no-build constraint** and a **construction-cost field** (right-of-way reuse discount,
  downtown tunnelling premium) folded into the flow physics as an effective tube length.
- Loads demand as a convex blend of **population density (0.4)**, **activity centers (0.3)**,
  and **transit-desert severity (0.3)** across 22 target neighborhoods + 16 existing anchors.
- Evolves the Tero current-reinforcement dynamics to a self-organized network.
- Compares against **minimum spanning tree**, **Steiner tree** (Mehlhorn 2-approx),
  and **hub-and-spoke** baselines using the Tero cost-benefit metrics: total length (TL),
  mean inter-terminal distance (MD), and fault tolerance (FT).

## Headline result (gamma = 1.3, balanced regime)
| Network            | TL (km) | MD (km) | FT   |
|--------------------|--------:|--------:|-----:|
| Minimum span. tree |   147.0 |   26.95 | 0.00 |
| Steiner tree       |   145.2 |   29.29 | 0.00 |
| Hub-and-spoke      |   132.9 |   19.23 | 0.00 |
| **Physarum 1.3**   | 247.5+/-11 | 18.69+/-0.34 | 0.88+/-0.02 |

The slime-mold network spends ~1.7x the MST length to cut mean trips ~30% AND raise
fault tolerance from 0 to ~0.89. A single exponent gamma traces the entire
cost-vs-robustness Pareto front (mesh at gamma<1, tree at gamma>1).

## Layout
```
slime mold/
  code/   data_chicago.py     # real coordinates, demand weights, lake polygon, cost field
          physarum.py         # Tero solver, candidate graph, baselines, metrics
          run_experiments.py  # drives runs + figures (stage 1 and stage 2)
  figures/ fig1..fig6 .png    # constraints, evolution, grown net, baselines, Pareto, gamma sweep
  data/    results.json/.csv  # all metrics; multiseed.json, gamma_sweep.json
  doc/     chicago_slime_mold.tex / .pdf   # full write-up (9 pp)
  REVIEW.md                   # 5x verification log
```

## Reproduce
```
pip install numpy scipy matplotlib networkx --break-system-packages
cd code
python run_experiments.py 1   # headline + baselines + figs 1-4 + results.csv
python run_experiments.py 2   # multi-seed + gamma sweep + figs 5-6 + results.json
```
Compile the paper: `cd doc && pdflatex chicago_slime_mold.tex` (run twice).

## Key references
- Tero et al., *Rules for biologically inspired adaptive network design*, Science 327:439-442 (2010).
- Nakagaki et al., *Maze-solving by an amoeboid organism*, Nature 407:470 (2000).
- Tero, Kobayashi, Nakagaki, J. Theor. Biol. 244:553-564 (2007).
- Bonifaci, Mehlhorn, Varma, *Physarum can compute shortest paths*, J. Theor. Biol. 309:121-133 (2012).
- Watanabe et al., *Biosystems* 105:225-232 (2011).

---

## Extension: protected bike-lane gap analysis (`bike/`)

A second study reusing the same Physarum engine on **real City of Chicago open data**
to find where to build protected bike lanes next.

- **Real inputs:** facility mileage by type and named protected corridors (Bike Routes
  dataset `hvv9-38ut`); 17,677 cyclist crashes / 40 deaths / 1,575 serious injuries
  (Traffic Crashes `85ca-t3if`). Pulled via the City's Socrata API.
- **Finding:** only 68.7 mi (15%) of mapped bikeways are protected. The grown ideal
  low-stress network is 229 km, of which **37% (84.5 km) has no protected lane**, gaps
  concentrated on the West and South sides (Austin–Belmont Cragin, North Lawndale,
  Englewood, Roseland).
- **Deliverable:** `bike/doc/chicago_bike_gap_brief.pdf` — a 6-page policy brief with a
  ranked build/upgrade list. Code in `bike/code/`, figures in `bike/figures/`,
  verification in `bike/REVIEW_bike.md`.

Reproduce: `cd bike/code && python bike_run.py` (needs `physarum.py` from `../code` on the
path; the script imports it directly).
