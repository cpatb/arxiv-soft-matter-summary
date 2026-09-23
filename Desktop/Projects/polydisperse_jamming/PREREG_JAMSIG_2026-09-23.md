# PREREG — `jamsig`: is there a camera-visible DNS signature of JAMMING, as opposed to thickening?

Written 2026-09-23, committed 00:07 CDT (git commit time), BEFORE any of the observables below was computed on any cell.
Post-processing only, on existing healthy L9 cells. Grader `grade_jamsig.py`, written with this file.

## Why
`kmermap` showed the footprint (Rg, β) does not discriminate: the largest cut is where nothing
jams. The bench test for the multiplet hypothesis needs an observable that does. This asks whether
the DNS itself carries one among the columns a camera could see.

## Cells (healthy, provenance-passed, L9, same binary `db413e87dd21`)
kmermap (12 healthy), kmerdecomp (5 healthy L9), kmerdrop (w0.00…w1.00 + all_w0.50). A cell is
**J+** if its J_max > 0.01, **J0** if J_max ≡ 0 (outside its window). No cell is re-labelled.

## Candidate observables (all from spread.dat; camera-visible)
- **S1 late motion**: mean of `u_rms` over t/τ ∈ [5, 6], divided by U₀ (the solver's velocity unit
  is not assumed; the ratio is compared only WITHIN one U₀).
- **S2 recoil**: max of h/D0 after its first minimum minus that minimum (apex rebound).
- **S3 late creep**: Rg_r(t/τ = 6) − Rg_r(t/τ = 3), relative to Rg_r(t/τ = 3).

## Test (per observable, WITHIN each U₀ that has both J+ and J0 cells)
- **Informative** only if the J+ and J0 ranges do NOT overlap at that U₀, AND the group medians
  differ by ≥ 1.6× (the project's discriminator bar). Otherwise NOT INFORMATIVE at that U₀.
- **J1 (prediction)**: S1 separates at every U₀ tested, with LOWER late motion for J+ (a jammed
  region carries the divergent viscosity). **J2**: S2 and S3 do NOT separate (both are dominated by
  the thickened bulk, like the footprint). No direction is registered for S2/S3.
- **Confound, stated before:** at fixed U₀ the J+ cells sit at higher x than the J0 cells (x is how
  they got into the window), so any separation may be an x effect, not a jamming effect. The
  discriminating subset is the one where x is shared: kmerdrop (all x 0.85, U₀ 2, J0 at w₂ 0/0.25,
  J+ at w₂ ≥ 0.5) and kmerdecomp/kmermap x 0.85 U₀ 4 (all J0) vs x 0.88 U₀ 4. **J3**: an observable
  counts as a jamming signature only if it also separates inside the kmerdrop ladder (shared x).
- Report every value; a separation that holds only across x is labelled "x-confounded".

## RESULT (committed 2026-09-23 00:08 CDT (git commit time), `grade_jamsig.py` → `out/grade_jamsig.txt`; 23 healthy cells)

| | verdict | numbers |
|---|---|---|
| **J1** S1 late motion separates, J+ lower | **FAIL** | overlaps at U₀ 2 (ratio 1.05) and U₀ 4 (1.32) |
| **J2** S2, S3 do NOT separate | **FAIL** | S2 separates at U₀ 2 and 4; S3 separates at U₀ 4 only |
| **J3** signature also at shared x (kmerdrop ladder) | **S2 only** | S2 J+ 0.016–0.078 vs J0 0.172–0.469 (no overlap, ×6.8) |

**S2 (apex recoil) separates everywhere it was tested**: after the apex height reaches its minimum,
J0 drops (thickened, cannot jam) rise again by h/D0 0.125–0.586, J+ drops by 0–0.078. (S2 = 0 means
no rise of ≥ 1 grid quantum was resolved, not a measured zero. Referee: in S2 = 0 cells h is pinned
for 235 samples while Rg and J keep evolving, so the solver is live and the apex genuinely does not
move; J0 traces show a real multi-quantum rise-then-relax, e.g. 30 quanta. So it is an arrested
vs oscillating apex, not a flicker.) Reading
(HYPOTHESIS): a thickened liquid still retracts under surface tension; a drop with a jammed core
does not. This is the first candidate camera-visible discriminator between jamming and thickening
in this DNS.

⚠ **Why it is a candidate and not a result.**
- (a) h/D0 is GRID-QUANTIZED (1/64 at L9). The separation is J+ ≤ 5 quanta against J0 ≥ 8 quanta,
  a 3-quantum gap, and several J+ values are exactly 0 or 1 quantum.
- (b) J2 was registered the other way, so S2 was not predicted to separate. The observable and its
  test share one set of 23 cells.
- (c) J+ and J0 also differ in footprint: J+ drops spread less. The closest shared-x pair is
  kmerdrop w0.25 (J0, Rg 0.3875, S2 0.172) vs w0.50 (J+, Rg 0.3631, S2 0.078).

Registered follow-up: Amendment 1 below.

## Amendment 1 — held-out and grid test of S2 (committed 2026-09-23 00:08 CDT (git commit time), before its cells)
Campaign `jamsig_2026-09-23` (launcher `scripts/jamsig_2026-09-23.sh`, same binary and settings as the
kmerdrop ladder; x 0.85, σ\* 300 Pa, δ 1.3169, χ(w₂) as `kmerdrop`). The threshold is FIXED here:
**S2\* = 0.1** (6.4 quanta, between the observed ≤ 5 and ≥ 8).
- **K1 held out (L9, U₀ 4, cells never seen):** kd_w0.25_U4.0, kd_w0.50_U4.0, kd_w0.75_U4.0. Every
  healthy J+ cell (J > 0.01) has S2 < 0.1 and every healthy J0 cell has S2 ≥ 0.1. J status comes from
  the run, not from this file. If all healthy cells fall in one class, K1 is UNINFORMATIVE on the
  other side and is reported as such.
- **K2 grid (L10, U₀ 2):** L10_w0.25_U2.0 and L10_w0.50_U2.0. The pair keeps its L9 classes
  (J0 / J+) AND its S2 sides of 0.1. A class change at L10 is reported and makes K2 undecidable. At
  L10 one quantum is 1/128, so the threshold is 12.8 quanta.
- Health as D5 of `PREREG_KMERDECOMP`; a dead cell decides nothing.

### Amendment 1 — RESULT (committed 2026-09-23 00:46 CDT (git commit time))
| cell | health | J | class | S2 | side of 0.1 |
|---|---|---|---|---|---|
| kd_w0.25_U4.0 (held out) | **died** (SIGFPE, t/τ 3.6) | 0 | J0 | ≥ 0.078 (partial) | undecided |
| kd_w0.50_U4.0 (held out) | healthy | 0.132 | J+ | 0.063 | predicted ✓ |
| kd_w0.75_U4.0 (held out) | healthy | 0.544 | J+ | 0.000 | predicted ✓ |
| L10_w0.25_U2.0 | healthy | 0 | J0 | 0.156 | predicted ✓ |
| L10_w0.50_U2.0 | **died** (rc 2, t/τ 0.075, at impact) | ≥ 0.087 | J+ | — | undecided |

**K1: UNINFORMATIVE** — 2/2 healthy held-out cells on the predicted side, but both are J+; the only
held-out J0 cell died, so the threshold's other side was never tested out of sample. **K2:
INCOMPLETE** — the J0 half holds at L10 (0.156 ≥ 0.1, i.e. 20 quanta at 1/128), the J+ half died at
impact. Every healthy new cell fell on its predicted side (3/3); no cell contradicts S2.
**Verdict: S2 (apex recoil) stays a CANDIDATE signature**, consistent with every healthy cell seen
(23 + 3), but its out-of-sample test is undecided on the J0 side and its grid test on the J+ side.
The deaths are the L9/L10 high-viscosity failure again; they are what blocks this, as they block
kmermap R2/R3 (Status §3 item 3). ⚠ Survivorship: the J0 death came late (t/τ 3.6, inside the
recoil window) and the J+ deaths at impact. If large recoils crash the solver preferentially,
surviving cells are biased toward gentle cases, so K1/K2 must not later be read as 'resolved by
attrition' without addressing that.

## Amendment 2 — out-of-sample test on SINGLETS across x (registered before its cells; commit time is the record)
Amendment 1's held-out test failed to inform because its J0 cell died. This one uses regimes that
SURVIVED at the default cap (the singlet kmermap cells at x 0.85/0.88/0.92, U₀ 2/4 were all
healthy). Eight never-seen singlet cells (χ 0.9083: φ_m 0.56315, σ\* 300 Pa, δ 1.3169, default
EPS_MIN 0.01, L9), campaign `jamsig2_2026-09-23`:
| label | x | φ | JAM_FHAT (f_c) | expected class |
|---|---|---|---|---|
| s_x0.86_U2.0 / _U4.0 | 0.86 | 0.5332 | 1.0 (outside window) | J0 by construction |
| s_x0.87_U2.0 / _U4.0 | 0.87 | 0.5394 | 1.0 (outside window) | J0 by construction |
| s_x0.93_U2.0 / _U4.0 | 0.93 | 0.5766 | 0.7634 | J+ expected, decided by the run |
| s_x0.94_U2.0 / _U4.0 | 0.94 | 0.5828 | 0.6543 | J+ expected, decided by the run |
- **K3 (the out-of-sample test, S2\* = 0.1 unchanged):** every healthy J+ cell has S2 < 0.1 and every
  healthy J0 cell has S2 ≥ 0.1. It is informative only if ≥ 1 healthy cell of EACH class exists.
- One particle type throughout, so composition, σ\* and δ cannot carry the separation; only x and
  the jam state differ.
- Health as before; deaths are reported, never graded.

### Amendment 2 — RESULT (commit time is the record; `out/grade_jamsig.txt`)
| cell | health | J | class | S2 |
|---|---|---|---|---|
| s_x0.86_U2.0 / U4.0 | healthy / healthy | 0 / 0 | J0 | 0.461 / 0.500 |
| s_x0.87_U2.0 / U4.0 | healthy / healthy | 0 / 0 | J0 | 0.391 / 0.391 |
| s_x0.93_U2.0 | died (t/τ 1.25, dt-collapse) | ≥ 0.126 | J+ | — |
| s_x0.93_U4.0 | healthy | 0.329 | J+ | **0.000** |
| s_x0.94_U2.0 / U4.0 | died (rc 134 at t/τ 1.23; t/τ 3.45) | ≥ 0.453 / ≥ 0.452 | J+ | — |

**K3 PASS: 5/5 healthy cells on the predicted side of the pre-fixed S2\* = 0.1, both classes present.**
This is the first registered OUT-OF-SAMPLE pass for apex recoil, on one particle type, so composition,
σ\* and δ cannot carry it. ⚠ **Its J+ side is n = 1**: 3 of 4 jammed cells died at the L9 floor,
exactly the survivorship pattern flagged under Amendment 1. S2 is therefore upgraded from
"candidate" to **"passed one out-of-sample test, J0 side n = 4, J+ side n = 1"**, not to validated.

## Amendment 3 — can EPS_SOFT rescue the L9 deaths? (registered before its cells; commit time is the record)
The deaths are the solver's HARD dt floor (dt < 1e-7, not scaled by DTFLOOR_FACT) tripping after
Poisson non-convergence: stop dt 6.8e-8 to 5.8e-11, lowdt mostly 1. So DTFLOOR_FACT cannot help. The
solver documents EPS_SOFT for exactly this: the hard clip at EPS_MIN puts a C0 corner "at exactly the
state the jammed plug occupies". EPS_SOFT = 0.002 makes the jammed-state viscosity
η_s/(EPS_MIN + EPS_SOFT/2)² = 0.83 η_max (a 17 % softer cap; classes survived 10× in `epssig`). For
unjammed cells far from the floor (ε − EPS_MIN ≫ EPS_SOFT) the change is ≈ EPS_SOFT²/(4d), 0.03 % at
x 0.86. Campaign `jamsig3_2026-09-23`, singlet cells as Amendment 2, default EPS_MIN, EPS_SOFT 0.002:
- rescue: `s_x0.93_U2.0`, `s_x0.94_U2.0`, `s_x0.94_U4.0` (died at default);
- controls: `s_x0.93_U4.0` (J+, S2 0 at default), `s_x0.86_U4.0` (J0, S2 0.500 at default).
Predictions:
- **T1 (rescue):** ≥ 2 of the 3 rescue cells are healthy.
- **T2 (controls keep class and S2 side of 0.1).**
- **K3′:** every healthy J+ cell in this campaign has S2 < 0.1 and every healthy J0 cell S2 ≥ 0.1. It
  counts as out-of-sample for the RESCUED cells only (the controls were seen).
- Health: rc 0, t/τ ≥ 5.9, and the tag must carry the EPS_SOFT marker.
