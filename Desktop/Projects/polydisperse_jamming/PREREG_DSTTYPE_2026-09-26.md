# PREREG_DSTTYPE_2026-09-26 — does the steady-shear TYPE (CST vs DST) decide what an impacting drop does?

Written 2026-09-25 ~23:10 CDT (the commit timestamp is authoritative), before any cell ran. Solver: `wc_drop_impact.c`, the G0-gated binary md5 db413e87dd21
(same as cvmap), L9, EPS_MIN 0.03, TEND 6, CFL 0.10, THETA0 5, WC_CLOCK 1, PHI0 0.62, MP_EXPO 2 (n = 2). Launcher
`scripts/dsttype_2026-09-26.sh` (via `scripts/run_dsttype.sh`); cells `out/dsttype_cells.txt`; grader `grade_dsttype.py`.

## Why
`derivations/thickening_type_2026-09-25.md` gives the steady-shear type in closed form: CST for x < x_DST(χ, nδ), DST
for x_DST ≤ x < χ, SJ above. POST HOC on the 24 cvmap drops (`derivations/cvmap_by_type_2026-09-26.txt`), the two
"thickened but unjammed" drops the project could not tell from jammed ones are exactly the two DST-band drops
(Rg_max 0.376/0.402 vs SJ 0.361–0.398 vs CST 0.488–0.670). Two readings predict that: **TYPE** (an S-shaped flow
curve makes the drop respond like a solid) and **Θ** (stiffness grows smoothly with the thickening ratio
Θ = [(1−x)/(1−x/χ)]^n, which is large near χ anyway). This test separates them with pairs matched in Θ but not in type.

## Cells (20; JAM_FHAT 1.0 everywhere: every cell is below χ, so none can jam)
Input set C5 (the cvmap CV 0.5 inputs: PHIM 0.56855, DELTA 1.0733, SIGMA* 220.5; χ 0.91702, nδ 2.147, x_DST 0.89854):
x = 0.8735 (Θ 7.1, CST), 0.8835 (10.2, CST), 0.8925 (16.2, CST), 0.9035 (43.1, DST), 0.9105 (160.3, DST).
Input set C0 (the cvmap CV 0 inputs: PHIM 0.56459, DELTA 1.3169, SIGMA* 300.0; χ 0.91063, nδ 2.634, x_DST 0.87919):
x = 0.8542 (5.5, CST), 0.8642 (7.1, CST), 0.8732 (9.5, CST), 0.8842 (15.9, DST), 0.8912 (26.0, DST).
φ = 0.62 x (5 dp); each at U₀ = 2 and 4.

## Registered tests (Rg = the cell's Rg_max from spread.dat column 11; health: rc 0 and t/τ ≥ 5.9; a dead cell is
EXCLUDED and reported, never a FAIL)
- **D1 (primary, matched Θ, different type).** Pair P16 = (C0 x 0.8842, DST, Θ 15.9) vs (C5 x 0.8925, CST, Θ 16.2);
  control pair P10 = (C0 x 0.8732, CST, Θ 9.5) vs (C5 x 0.8835, CST, Θ 10.2). Δ = ln Rg(C5 member) − ln Rg(C0 member).
  **TYPE** if Δ16 > 0.10 AND Δ16 > 3|Δ10| at both U₀; **Θ** if |Δ16| < 0.05 at both U₀; otherwise UNDETERMINED.
- **D2 (a jump at x_DST).** Per input set and U₀, on ln Rg against x: the step between the two cells bracketing x_DST
  exceeds 2× the largest step between adjacent same-type cells. Graded on the 4 (set, U₀) ladders; PASS if ≥ 3/4.
- **D3 (recoil).** DST cells S2 ≤ 0.2 and CST cells S2 ≥ 0.3 (the gap seen in cvmap; S2 as `grade_jamsig.sig`). Report
  the count.
- **D4 (ordering, report only).** Is Rg monotone decreasing in Θ across ALL 10 cells at each U₀? (Θ reading predicts yes.)

## Amendment 1 (2026-09-25 ~23:45 CDT, before ANY cell landed; referee 2026-09-26 MATERIAL #2)
D1's control pair spans two input sets that also differ in delta (23 %) and sigma* (36 %), and its Theta match is 6.7 %
(P16: 1.5 %), so a D1 "TYPE" could come from unequal sensitivity to the residual mismatch. Added, graded beside D1:
**D1' (within-set, no cross-set confound).** Per input set and U0, fit ln Rg = a + b ln Theta on that set's three CST
cells, extrapolate to the Theta of its FIRST DST cell (C0: x 0.8842, Theta 15.9, extrapolation x1.7 in Theta from
9.5; C5: x 0.9035, Theta 43.1, x2.7 from 16.2), and take the residual r = ln Rg(DST) - prediction.
**TYPE** if r < -0.10; **THETA** if |r| < 0.05; otherwise UNDETERMINED. Primary for D1' is set C0 (shorter
extrapolation); C5 is reported. The registered D1 is unchanged and still graded.
