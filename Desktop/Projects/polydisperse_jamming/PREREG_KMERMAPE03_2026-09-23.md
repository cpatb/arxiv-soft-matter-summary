# PREREG — `kmermapE03`: the kmermap regime map at a softer viscosity cap (EPS_MIN 0.03)

Written 2026-09-23 before any cell; commit time is the record. Launcher `scripts/kmermapE03_2026-09-23.sh`
(the kmermap grid, identical cells and binary `db413e87dd21`, with `EPS_MIN=0.03` added to every cell;
η_max 1.1e3 η_s instead of 1e4). Grader `grade_kmermapE03.py`, written before the data.

## Why
At the default cap, 6/18 kmermap cells died at L9 (dt-collapse under high viscosity), which left R2
and R3 INCOMPLETE. At EPS_MIN 0.03 all 4 `epssig` cells survived, including the x 0.88 U₀ 4 doublet
cell that died at the default. `epssig` also showed that J's MAGNITUDE depends on the cap while the
CLASS held on the two pairs tested. So this is a cap-robustness arm for the ORDERING and CLASS claims,
not a replacement for the default-cap numbers. Its J magnitudes are not comparable with kmermap's.

## Predictions (the kmermap sentences, re-registered for this cap)
- **R1′** x 0.88: doublets jam at ≥ 1 U₀ (J > 0.01); singlets J ≡ 0 (by construction).
- **R2′** x 0.92: J(w₂ 1) > J(w₂ 0) at EVERY U₀.
- **R3′** inside a window, J is non-decreasing in U₀.
- **R4′** β_max(w₂ 1) ≤ β_max(w₂ 0) + 1/64 at every (x, U₀).
- **R5′** health: rc 0, t/τ ≥ 5.9, banner eps_min = 0.03; else EXCLUDED (a death is never a jam).
- **X1 (cross-cap class):** every (w, x, U₀) healthy at BOTH caps has the same class (J+ if J > 0.01,
  J0 if J ≡ 0) at 0.03 as at 0.01. Cells in the window at J ≤ 0.01 are reported as ambiguous.
- Not tested here: the recoil signature (its contrast collapses at this cap by `epssig` G2).

## RESULT (graded 2026-09-23; commit time is the record; `grade_kmermapE03.py` → `out/grade_kmermapE03.txt`)
18/18 cells healthy, provenance 18/18 (eps_min 0.03 verified in every banner). **Every registered
gate passes:**
| gate | verdict | numbers (EPS_MIN 0.03) |
|---|---|---|
| R1′ doublets jam at x 0.88 | **PASS** | J 0.041 / 0.357 / 0.673 at U₀ 1 / 2 / 4; singlets 0 |
| R2′ doublet > singlet at every U₀, x 0.92 | **PASS 3/3** | 0.435 vs 0.058, 0.739 vs 0.338, 0.928 vs 0.646 |
| R3′ J non-decreasing in U₀ | **PASS 3/3** | all three in-window rows |
| R4′ β(doublets) ≤ β(singlets) + 1/64 | **PASS 9/9** | |
| X1 same class at both caps | **PASS 12/12** | every cell healthy at both caps |

- **Control, exact:** the five singlet cells at x 0.85/0.88 are BYTE-IDENTICAL to the default-cap run
  (spread.dat md5). Their ε = 0.064 / 0.031 exceeds the 0.03 floor, so the cap cannot bind there.
  Cells where it binds differ (w₁ x 0.85: ε 0.015; w₀ x 0.92). The override acted only where the algebra
  says it should.
- **The footprint reading also holds at this cap**: ln Rg cut (doublets vs singlets, U₀ 4) −0.329 at
  x 0.85 (nothing jams), −0.070 at x 0.88, −0.005 at x 0.92.
- J magnitudes are 1.5–6× the default-cap values, as `epssig` warned; only classes and orderings
  transfer between caps.

**Reading.** The orderings that the L9 deaths left INCOMPLETE at the default cap (R2, R3) are complete
here and pass. Every class seen at both caps is the same. Within the tested cap range
(1.1e3–1e4 η_s), the regime map is cap-robust in class and ordering. The softer cap is also the one
at which the recoil gap falls to a quarter (`epssig` G2), so recoil is NOT read from this run.

## Amendment 1 — a grid check of the headline (registered before its cells; commit time is the record)
R1 (doublets jam at x 0.88 where singlets cannot) has no grid check. At the default cap, L10 jammed
cells die at impact (`jamsig` L10_w0.50). At EPS_MIN 0.03 every L9 cell survived and the class is
cap-robust (X1 12/12). Two LEVEL-10 cells, otherwise identical to `kmermapE03` w1_x0.88_U2.0/U4.0
(campaign `kmermapE03L10_2026-09-23`). The singlet arm at x 0.88 is J ≡ 0 by construction (JAM_FHAT
1.0) and is not re-run.
- **R1L10:** J > 0.01 in ≥ 1 healthy L10 doublet cell (L9 values 0.357 / 0.673).
- Reported, no bar: J(L10)/J(L9) per U₀ (J is an f̂ quantity, which lessons.md records as
  grid-sensitive in this solver family, so only the class is registered).
- Health: rc 0, t/τ ≥ 5.9, banner tag `_L10_` and eps_min 0.03.

### Amendment 1 — RESULT (commit time is the record)
Both L10 cells healthy (tag `_L10_`, eps_min 0.03 verified). **R1L10 PASS**: J = 0.386 (U₀ 2) and
0.672 (U₀ 4) against L9 0.357 / 0.673, ratios 1.083 / 0.998. The headline class holds at L10, and at
this cap J is grid-converged to within 8 % (U₀ 2) and 0.2 % (U₀ 4). That is better than the f̂-family
grid sensitivity recorded for this solver elsewhere, but it is measured on these two cells only.

## Amendment 2 — does the switch survive the DEM pressure-rung convention? (registered before its cells)
Every DNS input so far takes χ from the p̂ 3e-5 rung (singlets 0.9083, doublets 0.8631), which puts the
switch window at x ∈ [0.863, 0.908). At the p̂ 1e-3 rung (`grade_kmer.txt`: 0.9284 / 0.8978, both
clean) the window is [0.898, 0.928). "Doublets jam at x 0.88" is therefore RUNG-SPECIFIC by algebra:
at 1e-3, x 0.88 lies outside both windows. The rung-free claim is that a switch window EXISTS. Test it
at the 1e-3 rung. Doublet cells at x 0.91: φ 0.5642, φ_m 0.55664, JAM_FHAT = f_c = 0.8806, σ\* 235.8 Pa,
δ 1.1163. Singlets at x 0.91 are outside their 1e-3 window (0.91 < 0.9284), J ≡ 0 by construction.
Campaign `kmermapRG_2026-09-23`: `rg_d_x0.91_U2.0_e03`, `rg_d_x0.91_U4.0_e03` (EPS_MIN 0.03) and
`rg_d_x0.91_U2.0_e01`, `rg_d_x0.91_U4.0_e01` (default, reported).
- **RG1:** at EPS_MIN 0.03, J > 0.01 in ≥ 1 healthy doublet cell. FAIL if both are healthy with J ≤ 0.01.
- Default-cap cells reported only (deaths expected).

### Amendment 2 — RESULT (commit time is the record)
| cell (doublets, x 0.91, χ from p̂ 1e-3) | health | J |
|---|---|---|
| U₀ 2, EPS_MIN 0.03 | healthy | 0.318 |
| U₀ 4, EPS_MIN 0.03 | healthy | 0.624 |
| U₀ 2, default | healthy | 0.061 |
| U₀ 4, default | died t/τ 2.6 | ≥ 0.122 |
**RG1 PASS.** With χ from the other DEM pressure rung, doublets jam inside that rung's switch
window, where singlets cannot (0.91 < 0.9284). A healthy default-cap cell agrees (J 0.061). **The
existence of a doublet-opened switch window is rung-robust; only its position in x is a convention.**
