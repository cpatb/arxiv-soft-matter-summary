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
