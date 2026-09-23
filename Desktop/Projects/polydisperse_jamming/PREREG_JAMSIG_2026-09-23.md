# PREREG — `jamsig`: is there a camera-visible DNS signature of JAMMING, as opposed to thickening?

Written 2026-09-23 ~01:40 CDT, BEFORE any of the observables below was computed on any cell.
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

## RESULT (graded 2026-09-23 ~01:55 CDT, `grade_jamsig.py` → `out/grade_jamsig.txt`; 23 healthy cells)

| | verdict | numbers |
|---|---|---|
| **J1** S1 late motion separates, J+ lower | **FAIL** | overlaps at U₀ 2 (ratio 1.05) and U₀ 4 (1.32) |
| **J2** S2, S3 do NOT separate | **FAIL** | S2 separates at U₀ 2 and 4; S3 separates at U₀ 4 only |
| **J3** signature also at shared x (kmerdrop ladder) | **S2 only** | S2 J+ 0.016–0.078 vs J0 0.172–0.469 (no overlap, ×6.8) |

**S2 (apex recoil) separates everywhere it was tested**: after the apex height reaches its minimum,
J0 drops (thickened, cannot jam) rise again by h/D0 0.125–0.586, J+ drops by 0–0.078. Reading
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

## Amendment 1 — held-out and grid test of S2 (written 2026-09-23 ~02:05 CDT, before its cells)
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
