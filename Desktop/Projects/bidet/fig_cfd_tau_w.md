# fig_cfd_tau_w.png

CFD wall shear on the plate for a round jet impinging normally, axisymmetric
k-ω SST, Re_d = 10⁴, H/d = 6, V₀ = 8.9 m/s, radial domain R = 50 mm.
Case `cfd_openfoam/case_Re10000_H6_R50`. Regenerate with
`python cfd_openfoam/fig_cfd_tau_w.py case_Re10000_H6_R50`.

## What is plotted

- **CFD** — τ_w(r) at 501 sample radii, from `tau_w_r.csv`. Peak 389.4 Pa at
  r = 0.92 mm; 0.455 Pa at the outer boundary.
- **Turbulent fit** — the Phares closure τ_w/(ρV₀²)·Re_d^{1/2} = C_P (r/d)^{−n_t}
  over its own fit band r/d ∈ [2, 20] (shaded), n_t = 1.898, C_P = 1.353,
  r² = 0.99688. Drawn only over the band it was fitted on.
- **Laminar** — the Watson Region-II r^{−5} law, for contrast. It is far steeper
  than the data, which is the point: this flow is turbulent well before r_max.
- **τ_c = 1 Pa** — the cleaning threshold best supported by the literature
  (Nejadnik et al. 2008 measure 1.0 ± 0.2 Pa for *S. aureus* and 1.1 ± 0.2 Pa for
  *P. aeruginosa* in genuine wall detachment). Crossed at **r = 34.79 mm**, marked
  with a star.
- **R = 50 mm** (solid) — the outer radial boundary of this computation.
- **earlier R** (dotted, 25 mm) — where the case used to end, with τ_w still at
  1.96 Pa. The τ_c = 1 Pa contour then lay *outside* the domain and could only be
  bounded (r > 25 mm) or extrapolated.

## Why the domain was widened, and what it settled

The paper's most literature-supported cleaning threshold was the only one its
computation could not resolve. Re-running at R = 50 mm puts the contour at
34.79 mm, and two extrapolations can now be scored against it:

| method | r at τ_c = 1 Pa | error |
|---|---|---|
| **computed, R = 50 mm** | **34.79 mm** | — |
| outer-half refit, r/d ∈ [10, 24.8], n_t = 2.10 | 34.75 mm | −0.1% |
| full-band closure extrapolated past r/d = 20 | 39.14 mm | +12.5% |

The far field is steeper than the whole-band fit, so extrapolating a
single-power-law fit through a profile whose local slope is still steepening
overestimates the reach. A local refit does not. That is worth more than the
number itself.

**Validation.** The widened case must contain the narrow one, or the comparison is
worthless. Over r ≤ 22.5 mm the two agree to a median 0.7% and a maximum 1.8% in
τ_w. They diverge only in the outer 10% of the narrow domain, reaching 5.7% at its
last cell — the expected signature of an outflow boundary reaching back into the
profile it was meant to terminate, and exactly the artifact the re-run removed.
Run `python cfd_openfoam/compare_domain_extent.py` to reproduce.

**Mesh.** 60 120 cells against 30 720, with the first radial cell held at 7.0512 µm
so the near-impingement peak is resolved identically; max aspect ratio 106.6 against
105.9 and identical non-orthogonality. Converged in 2689 SIMPLE iterations,
y⁺ = 0.047–1.39 (median 0.45).

## Axes

Both axes are logarithmic, against the standing preference for zero-anchored linear
axes, and the reason is that τ_w spans 0.455 to 389 Pa (2.9 decades) over a radius
range of 0.02 to 50 mm (3.4 decades). A linear axis from zero would compress
everything past r ≈ 2 mm — which is the entire region the paper's design claims live
in — into the bottom-left corner. The power-law behaviour that the figure exists to
show is only visible as a straight line on log-log.
