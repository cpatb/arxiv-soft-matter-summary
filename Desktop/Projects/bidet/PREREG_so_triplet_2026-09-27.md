# Pre-registration: second-order turbulence advection, r = 1.5 triplet

**Registered 2026-09-27 ~21:38 CDT, before any case below was generated or run** (lanes
launched 21:41; `run_so_lane.sh`, staged in WSL `~/bidet_so`).
Written against `GRID_STUDY_2026-09-03.md` §3b/§3d/§6 items 1 and 3.

## Question

`C_P` diverges on the `fine_conv` family (refine 4/2/1, |R| = 4.59) and no quantity has an
observed order inside [0.5, 2.5]. `div(phi,k)` and `div(phi,omega)` are first-order
`bounded Gauss upwind`, so the effective diffusion in the free shear layer of the jet
falls with h at first order and that can change the impingement profile, and with it
`C_P` and `n_t`. **Does replacing them with a second-order bounded TVD scheme give
converging orders?**

## Design

Re_d = 1e4, H/d = 6, Rdom = 25 mm. Everything else unchanged from `gen_case.py`.

- **SO** family: `div(phi,k)` and `div(phi,omega)` = `bounded Gauss limitedLinear 1`
  (the standard OpenFOAM second-order TVD choice for turbulence scalars), refine
  **1 / 1.5 / 2.25**, so r = 1.5 exactly and every block count is an integer with no
  truncation (nx 120/180/270, nrj 16/24/36, nro 240/360/540).
- **FO** control family: identical, first-order upwind as shipped, refine 1.5 and 2.25.
  Its refine-1 member is the existing converged `case_Re10000_H6`.
- endTime 30 000, `residualControl` unchanged (1e-5). A case that does not meet it is
  **excluded, not graded** (health separate from the bar).

## Gates (graded by a script reading `tau_w_r.csv` through `fit_closure.py`)

Per family and per quantity (`C_P`, `n_t`, `tau_w(15 mm)`, peak `tau_w`), Celik R = eps21/eps32:

- **S1 (primary).** SO `C_P`: 0 < R < 1 and p in [0.5, 2.5].
- **S2.** SO `n_t` and `tau_w(15 mm)`: 0 < R < 1 (monotone converging). p reported.
- **S3 (discriminator — must be able to fail).** The attribution "the first-order scheme
  caused the divergence" is licensed only if FO on the SAME 1/1.5/2.25 grids is NOT
  monotone-converging for `C_P` (R <= 0 or |R| >= 1). If FO converges too, the scheme
  contrast is **NOT INFORMATIVE**, whatever S1 says, and the r = 2 divergence is then a
  property of the refine-4 grid, not of the scheme.
- **S4 (Uz stall).** FO refine 1.5 meets residualControl: the refine-1.414 limit cycle
  was specific to that truncated grid. Fails if the Uz residual stalls above 1e-4.
- **S5 (anchor that does not pass through the triplet).** SO refine 1 vs FO refine 1:
  predicted |Delta C_P| / C_P > 1 % (the scheme matters at the base grid). If it is below
  1 %, the scheme is not a lever and S1's outcome cannot be attributed to it.

## Predictions (mine, before data)

- S1 PASS, S3 FO diverges or oscillates, S4 PASS, S5 PASS (scheme moves `C_P` by several %).
- Competing reading: the top-hat inlet's lip discontinuity or the stagnation region is the
  cause, in which case SO diverges too and S1 FAILS.

## Amendment 1 (21:55, before any refine-2.25 case reached convergence; no gate value seen)

Measured cost on the shared VM (load ~20/22): ~5 s per iteration for refine 2.25 on 2
ranks, so a cold start (~6000 iterations, by analogy with refine 2's 5105) would not
finish before ~06:00. Both refine-2.25 cases are therefore RESTARTED from a mapFields
initialisation of the converged refine-1 case of the same scheme (FO from
`case_Re10000_H6`, SO from `case_Re10000_H6_so_r1` once it converges). Only the initial
condition changes; every gate reads the converged steady state, which a unique steady
solution does not remember. Check that goes with it (S6, reported): the mapped-start FO
refine-2.25 case is graded exactly as the others — if its converged `C_P` were an
artefact of the start, the FO r=1.5 triplet would show it as a non-smooth jump, which S3
reports either way. The cold-start runs are killed at ~iteration 100 and their partial
output is discarded.

## What will NOT be done

No bar is loosened after data. p outside [0.5, 2.5] is reported, not re-banded. If the
triplet cannot finish tonight, the verdict is INCOMPLETE.
