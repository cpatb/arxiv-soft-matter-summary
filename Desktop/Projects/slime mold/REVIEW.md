# 5x Review Log

Five independent verification passes were run over the approach, numerics, and write-up.
All issues found were fixed; the items below record what was checked and the outcome.

## Pass 1 - Model fidelity to the literature (independent review agent + web search)
- [OK] Flux law Q=(D/L)(p_i-p_j), Kirchhoff/Laplacian solve, and adaptation
  dD/dt = f(|Q|) - mu*D with f(Q)=Q^g/(1+Q^g) match Tero et al. (2010). The mu
  coefficient (=1 in the original) is a correct generalization.
- [OK] gamma>1 -> tree, gamma<1 -> mesh, gamma~1 -> balanced/Tokyo regime: correct direction.
- [OK] Cost-benefit metrics TL / MD / FT correctly identified and defined.
- [FIXED, was MAJOR] Watanabe citation year was inconsistently shown as 2014; corrected
  everywhere to Biosystems 105:225-232 (2011).
- [FIXED, was MINOR] Steiner-tree claim was over-attributed to Bonifaci et al. (2012),
  who proved only single source-sink shortest-path convergence. Softened to
  "empirically observed to approximate the minimum Steiner tree."
- [FIXED, was MINOR] Added a caveat that the gamma~1 balance is approximate and that the
  original Tokyo study used a tuned feedback form.

## Pass 2 - Numerical correctness (automated assertions)
- [OK] Kirchhoff residual: max |L p - b| at non-ground nodes = 0.000 (exact solve).
- [OK] MST baseline is a genuine tree (114 edges, 115 nodes) and FT = 0.000, as it must be.
- [OK] Physarum network connects ALL demand terminals after the feasibility repair.
- [OK] Multi-seed spread is tight (FT 0.88+/-0.02, MD 18.69+/-0.34 km over 8 seeds):
  results are not a random-seed artifact.

## Pass 3 - Geography / constraints
- [OK] Lake Michigan no-build: a test point east of the shoreline is correctly OUTSIDE
  the land polygon; 100% of candidate-edge midpoints lie on land (no lake crossings).
- [OK] Station coordinates spot-checked: O'Hare (-87.904, 41.980), Loop (-87.630, 41.879),
  95th/Dan Ryan (-87.624, 41.722), Midway (-87.738, 41.787) - all correct.
- [OK] Demand weights: transit deserts (Englewood, Austin, Roseland) carry the highest
  underserved scores; dense areas (Lakeview, Belmont Cragin) carry the highest population.

## Pass 4 - Figures vs. reported numbers
- [OK] Table 1 and the abstract match results.json exactly (MST TL 147.0; Physarum
  247.5+/-11.1; MD 18.69; FT 0.88).
- [OK] Pareto figure: baselines sit on FT=0; Physarum family forms the dominating front.
- [OK] Gamma sweep is monotone (TL 1734->117, FT 1.0->0.0 as gamma 0.7->2.2).

## Pass 5 - Reproducibility & internal consistency
- [OK] All three code modules compile and run end-to-end from the delivered copies
  (verified by re-import and a fresh candidate-graph build: N=678, E=1937).
- [OK] Seeds are fixed; stages documented in README; LaTeX compiles with zero undefined
  citations or references (9 pages).
- [OK] No dead cross-references; bibliography complete.

## Residual limitations (documented in the paper, Sec. 7)
Schematic existing-line geometry and lattice abstraction; normalized (not ridership-
calibrated) demand; structural (graph) fault tolerance rather than an operations-level
reliability simulation. None affect the central, seed-stable finding.
