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
