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

## Amendment 2 (~23:15, before any gate quantity of a refine-1.5 or 2.25 case was computed)

(Operational note, 23:20: FO refine 1.5 is also started by mapFieldsPar from
`case_Re10000_H6`, on 2 ranks, instead of cold on the serial lane, whose driver was
detached; SO refine 1.5 continues cold and serial and is finalised by hand.)

**Finding that forces it.** Both refine-2.25 cases stalled with the worst initial residual
~1e-3, always `Uz`, while every other field fell. `uz_magnitude.py` on the written fields:
FO refine 2.25 at iteration 1000 has max|Uz| = 3.4e-13 m/s against max in-plane 9.0 m/s
(ratio 3.7e-14); the old refine-1.414 case (`r1_ext`, the "broken" midfine grid of
GRID_STUDY §3d) has max|Uz| = 3.9e-14. The converged cases have 2.1e-10 to 2.6e-10. `Uz` is
the out-of-plane component on a one-cell wedge and is round-off everywhere; OpenFOAM
normalises each component's residual by that component's own scale, so a round-off field
returns an O(1e-3 to 1e-2) "residual" that no amount of iteration reduces. Whether a
case's `Uz` happened to dip under 1e-5 (1.6e-6 to 4.9e-6 on the converged ones) is luck.
`residualControl` on the vector U therefore cannot fire on some grids for a reason
unrelated to the solution.

**Health rule, replacing "meets residualControl":** a case is converged when (H1) every
initial residual except `Uz` is below 1e-5 at the same iteration (the registered target,
unchanged), and (H2) max|Uz| / max|U_inplane| < 1e-10 in the final field (the component
is round-off, so excluding it hides nothing). Both are printed per case. Operationally a
watcher (`so_watch.sh`) sets `stopAt writeNow` when H1 holds, since simpleFoam will not.
The original criterion's verdict is reported beside it for every case.

**S6 (Cauchy check on the gate quantities, registered 23:30 before any refine-1.5/2.25
gate value was printed).** The mapped-start cases reached H1 within a few hundred
iterations and then sat on a residual plateau (~5e-6), and a plateaued normalised residual
does not prove the field has stopped moving. So each mapped-start case (FO/SO 2.25, FO 1.5)
is CONTINUED 400 further iterations from its final state (`so_continue.sh`), and every
gate quantity (`C_P`, `n_t`, `tau_w(15 mm)`, peak `tau_w`) must move by less than 0.5 %
relative. A case that fails is NOT CONVERGED on the quantities and its family is graded
INCOMPLETE; both values are printed. 0.5 % is set against the smallest grid step the gates
must resolve (the refine-1 FO/SO difference in `n_t` was 1.2 %).

**Consequence registered now, graded later:** under H1/H2 the refine-1.414 `midfine`
triplet (§3d, "cannot run") may be gradeable after all. It is graded as a SEPARATE,
first-order family and does not enter S1-S5.

## Amendment 3 (23:35, before any SO refine-1.5/2.25 gate value was printed): S7

The first-order study's trouble appears only at its finest step (refine 2 → 4; the midfine
family shows C_P merely wandering over refine 1–2), so a triplet ending at 2.25 could
pass S1 without testing the regime that failed. **S7 (reported, not part of S1-S5):** SO
refine 3.375 (nx 405, nrj 54, nro 810; r = 1.5 again, all integer), started by
mapFieldsPar from SO refine 2.25, health H1/H2 and the S6 Cauchy continuation as above;
the second SO triplet 3.375 / 2.25 / 1.5 is graded with the same classification. If
S1 passes and S7's triplet is not monotone-converging, S1's pass is reported as not
robust to further refinement.

## Amendment 4 (~23:45 27 Sep; git commit time is authoritative): convergence judged on the gate quantities themselves

**Why.** H1 is not a convergence test for a mapped start: SO refine 3.375, mapped from
2.25, met H1 at iteration 26, which cannot be a converged 350k-cell solution. A good
initial guess makes the normalised residual small from the outset. (The original
`residualControl` has the same weakness; it was only meaningful for cold starts, which
needed thousands of iterations to reach it.) Meanwhile the cold SO refine-1.5 run has
plateaued at 2.2e-5 and may never meet H1. S6's 400-iteration window can miss a slow drift.

**Rule C1, replacing H1 as the health criterion for every refine > 1 case of both
families and for S7** (refine 1 cases are cold-started and met the original criterion):
each case is continued 2400 iterations from its current final state, writing every 400
(`so_trajectory.sh`); the wall-shear sample is regenerated at each written time by
`simpleFoam -postProcess` (verified to reproduce an in-run record to 4e-9 relative) and
fitted by `fit_closure.py` with `FIT_TIME`. A case is converged when, for every gate
quantity, the relative change between the last two writes is < 0.1 % and between the
+1600 and +2400 writes is < 0.2 %. The graded value is the +2400 one. H2 (Uz round-off)
still applies. A case that fails C1 makes its family INCOMPLETE. S6 (400-iteration
Cauchy) is superseded, since C1 starts from the same state over a longer window; the one
S6 value already on screen (FO refine 1.5 continued: `C_P` 1.3821, printed by its
launcher) is disclosed here. The SO refine-2.25 S6 run is stopped unfinished.

## Amendment 5 (00:40 28 Sep; commit time authoritative): C1 caught a drifting case

**What happened.** FO refine 1.5 — stopped by H1 at iteration 189 after its mapped start —
was continued 2400 iterations (`convergence/`, grader) and is **still drifting**: `C_P`
1.3808 → 1.4060 (+1.8 %), peak `tau_w` 394.8 → 417.3 (+5.7 %), `tau_w`(15 mm) 5.848 → 5.658
(−3.3 %), still ~0.8 %/400 iterations on `C_P` at the end. It FAILS C1. (S6 would also have
failed it, on `tau_w`(15 mm) 0.64 % and peak 1.12 %, while `C_P` alone moved 0.27 %.)
**So a mapped start stopped on a residual is not converged, and every gate value would have
been read off a transient.**

**Operational changes, no bar moved.** (a) A case failing C1 is continued (`TOTAL=4800`
env in `so_trajectory.sh`, chained as `<case>_traj_traj`), and C1 is applied to the end of
the concatenated trajectory with the same 0.1 %/0.2 % bars. (b) S7 (refine 3.375, mapped
from iteration 29 of refine 2.25) is stopped at iteration 563 to free its threads; it is
reported NOT RUN unless restarted after 03:00. (c) Registered now, before any legacy
trajectory exists: the cold-started refine-1 FO case (`case_Re10000_H6`, certified only by
the original residual criterion) is continued 2400 iterations the same way (**L1**). If its
gate quantities move by more than the C1 bars, every `C_P` in GRID_STUDY and in the
manuscript carries an iterative error that the residual criterion did not see.

## L1 result (01:45 28 Sep) and L2 registered

**L1:** the cold-started refine-1 FO case (`case_Re10000_H6`, the manuscript's grid, which
met `residualControl` at iteration 2261) continued to 4661
(`convergence/L1_legacy_r1_trajectory.txt`): `C_P` 1.3518 → 1.3663 (**+1.08 %**), `n_t`
1.9023 → 1.9102 (+0.42 %), `tau_w`(15 mm) 5.866 → 5.792 (**−1.26 %**), peak 389.5 → 390.1;
flat from ~3600 (last-800 window < 0.02 %). **The paper's own numbers carry ~1 % iterative
error that the residual criterion did not see.**

**L2 (registered now, run when threads allow):** the same continuation on the legacy
refine-2 case (`case_Re10000_H6_r2_ext`, met the criterion at 5105). Reported with the same
C1 bars. If refine 2 and refine 4 drift by different amounts, the `fine_conv` "divergence"
(|R| = 4.59) is partly an iterative artefact and must be re-graded on converged values;
refine 4 (491k cells) is out of reach tonight and is stated as such.

## What will NOT be done

No bar is loosened after data. p outside [0.5, 2.5] is reported, not re-banded. If the
triplet cannot finish tonight, the verdict is INCOMPLETE.
