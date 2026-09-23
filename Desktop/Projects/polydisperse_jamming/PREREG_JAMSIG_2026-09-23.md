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
