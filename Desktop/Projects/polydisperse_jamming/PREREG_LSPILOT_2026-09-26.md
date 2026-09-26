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
