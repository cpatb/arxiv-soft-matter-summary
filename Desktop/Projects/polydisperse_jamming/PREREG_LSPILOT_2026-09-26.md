# PREREG_LSPILOT_2026-09-26 — PILOT: does a convex faceted grain (rounded cube) widen the shear-jamming window?

Written 2026-09-25 ~23:25 CDT before these cells ran. **A pilot, not a measurement**: one seed, N = 100, one rung.
Script `ls_pack.py` (level-set superellipsoids, eps_e = eps_n = eps, volume-equivalent to d = 1; volume formula
anchored at eps = 1, 1e-3 and by Monte Carlo at 0.5); launcher `scripts/ls_pilot.sh`.
Disclosed shakeout (`out/lsdem/ls_smoke.jsonl`, seed 99, N 60, NOT converged): eps 0.5, mu_p 0.3 → phi 0.635 at
p_hat 1.15e-3.

Cells: eps {1.0 (sphere, same level-set method), 0.5 (rounded cube)} x mu_p {0, 0.3}, seed 1, N 100, p_hat 1e-3,
max 300000 iterations, grid spacing 0.05 of the extent diameter, 102 surface nodes.

**P1 (the question):** window width W(eps) = 1 − phi(mu .3)/phi(mu 0). Registered expectation from the counting
argument (a non-sphere's rotations must be pinned frictionlessly, not frictionally): W(0.5) > W(1.0). Rival: facets
act only through the frictional branch like a contact moment limit, lowering phi_m, which ALSO widens W — so the
sign alone does not discriminate them; the pilot reports W, phi_0 and phi_m separately.
**P0 (method control):** the level-set sphere reproduces the FrictMat sphere at the same N and rung within 0.01 on
each branch (a loose bar: N 100, one seed). If P0 fails, P1 is reported as method-limited.
Health: converged (unbalanced < 0.02 before the cap); an unconverged cell is reported, never graded.

## Amendment 1 (2026-09-25 ~23:55, after P0 FAILED on the three landed cells; registered before these runs)
P0 failed: level-set sphere phi 0.6576 / 0.5985 vs FrictMat 0.6332 / 0.5800 (mu 0 / 0.3, N 100, p 1e-3), +0.024 / +0.019.
Two-body probe (`scripts/ls_stiffness_probe.py`): kn uses the extent radius (0.499 sphere, 0.435 eps 0.5 — a 13 %
stiffness difference, negligible for phi at fixed p_hat), but for a true overlap of 0.010 the level-set SPHERE reports a
penetration of 0.0035 (102 surface nodes, node spacing ~0.17 d, the deepest point falls between nodes) while a face-on
eps-0.5 contact reports 0.0116. The artefact is SHAPE-DEPENDENT, so P1's comparison is method-limited as registered.
**N-test (node convergence):** eps 1.0, mu 0.3, same cell, with 400 and 1000 surface nodes. Prediction if the artefact is
node resolution: phi falls monotonically toward the FrictMat 0.5800 (|phi(1000) - 0.5800| < |phi(102) - 0.5800| / 2).
If it does, eps 0.5 at mu 0 / 0.3 is re-run at the node count that passes P0 (within 0.01) and P1 is graded on those.

## Amendment 2 (2026-09-26 ~00:10, after the N-test; before these cells run)
N-test RESULT: eps 1.0, mu 0.3 phi = 0.5985 / 0.5929 / 0.5797 at 102 / 400 / 1000 surface nodes vs FrictMat 0.5800:
monotone, and |phi(1000) - 0.5800| = 0.0003 < 0.0185/2 -> **PASS; the P0 failure was node resolution.**
All further level-set cells use **1000 nodes**. Ladder (seed 1, N 100, p 1e-3, max 300000 iterations): eps {1.0, 0.7,
0.5, 0.35} x mu_p {0, 0.3} (eps 1.0 mu 0.3 at 1000 nodes already exists). P0 is re-graded at 1000 nodes on both branches
(bar 0.01 vs FrictMat 0.6332 / 0.5800). P1 as registered, on these cells: W(eps) = 1 - phi(mu .3)/phi(mu 0); expectation
W(eps < 1) > W(1). Also reported: phi_0(eps), phi_m(eps) and the (G) decomposition (Z, w) per cell, and whether
the window at eps 0.7 (a mildly faceted near-sphere) already carries most of the eps-0.35 opening. One seed: every
difference is reported against the FrictMat seed scatter at N 100 (not measured here) — i.e. as a pilot.
