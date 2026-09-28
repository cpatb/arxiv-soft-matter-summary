# The grid study, continued: the "broken" grid was never broken, and a scheme test

**2026-09-27 overnight.** Pre-registration: `PREREG_so_triplet_2026-09-27.md` (registered
before launch, two amendments and S6 registered before any fine-grid gate value).
Grader: `cfd_openfoam/grade_so_triplet.py` → `convergence/so_triplet_grade.md`.
Status of the gates: see §3 (filled in when the last case lands).

---

## 1. GRID_STUDY §3d is answered: the Uz "limit cycle" is a residual normalisation artefact

GRID_STUDY §3d recorded the refine-1.414 grid as *genuinely stuck* — 20 000 iterations,
every field ~1e-9 except `Uz` at 1.8e-2 — and §6 made fixing it the cheapest open item.
Tonight both refine-2.25 cases showed the same symptom (`Uz` ~1e-3, flat), on grids with
no truncated cell counts, so the suspected cause (odd `nx` = 169) was wrong.

**The field says what the residual cannot.** `cfd_openfoam/uz_magnitude.py` on the
written `U`:

| case | max\|Uz\| (m/s) | max\|Uz\| / max\|U_in-plane\| | `Uz` residual | met residualControl? |
|---|---|---|---|---|
| refine 1 (converged) | 2.6e-10 | 2.8e-11 | 1.6e-6 | yes |
| refine 2 `_ext` (converged) | 2.1e-10 | 2.3e-11 | 3.3e-6 | yes |
| refine 1.414 `_ext` ("stuck") | 3.9e-14 | 4.3e-15 | **1.8e-2** | **no** |
| FO refine 2.25, it. 1000 ("stalled") | 3.4e-13 | 3.7e-14 | 7.8e-4 | no |

`Uz` is the out-of-plane component on a one-cell wedge and is **round-off everywhere,
including on the "stuck" grid** — where it is *smaller* than on the converged ones.
OpenFOAM normalises each component's residual by that component's own scale, so a
round-off field returns an O(1e-3–1e-2) residual that iteration cannot reduce, and
whether a case happened to see `Uz` dip below 1e-5 (the converged ones sit at 1.6–4.9e-6)
was luck. **`residualControl` on the vector `U` can fail for a reason unrelated to the
solution.** The refine-1.414 case's other fields reached 9.7e-9.

Health rule used from here on (Amendment 2): every initial residual except `Uz` below
1e-5 (H1) and `Uz` at round-off in the final field (H2), with the original verdict
printed beside it.

## 2. With the grid admitted, the independent order check runs — and C_P oscillates

The `midfine` family (refine 2 / 1.414 / 1, first-order `k`/`omega`, the shipped
schemes; r from cell counts 1.419 / 1.409):

| quantity | refine 2 | refine 1.414 | refine 1 | R = ε21/ε32 | class |
|---|---|---|---|---|---|
| `C_P` | 1.3928 | 1.4098 | 1.3518 | −0.293 | **oscillatory** |
| `n_t` | 1.8922 | 1.9237 | 1.9023 | −1.469 | non-monotone |
| `τ_w`(15 mm) | 6.2415 | 5.7726 | 5.8658 | −5.034 | non-monotone |
| peak `τ_w` | 441.17 | 417.2 | 389.52 | +0.866 | converging, p = 0.473 (outside band) |

Read together with `fine_conv` (refine 1 / 2 / 4), the five first-order grids give
`C_P` = 1.3518, 1.4098, 1.3928, 1.5809 at refine 1, 1.414, 2, 4 (and refine 0.5 from
`coarse_conv`). **Over refine 1–2 `C_P` wanders by ±2–4 % without a direction; the
"divergence" of `fine_conv` (|R| = 4.59) is carried by the single step to refine 4
(+13.5 %).** That changes what the refine-4 grid is: not the last point of a runaway
sequence but the one outlier of a noisy one. It does not make the study converge — no
quantity has an in-band order on any first-order family — but it relocates the problem
from "the sequence runs away" to "something happens between refine 2 and 4", which is
where the next check belongs (§4).

Caveat: the refine-1.414 case's saved `U` is at iteration 3000 (the 20000 directory holds
only the post-processed fields), so H2 there is read on the 3000 field.

## 3. The registered scheme test (S1–S7)

*Pending — filled from `convergence/so_triplet_grade.md`.* Health is now judged by rule C1
(Amendment 4, QoI trajectory over +2400 iterations), because H1 turned out to be no test
at all for a mapped start (refine 3.375 met it at iteration 26).

## 4. Exploratory, post hoc: where refine 4 departs, and a model-form fact

`convergence/explore_band_fo_ladder.txt`. **Not a registered gate; it re-reads the
first-order grids already on disk and none of the gated SO fine grids.**

**Radially resolved, refine 4 does not shift the profile, it tilts it.** Median Phares
LHS per r/d bin relative to refine 2 (`_ext`):

| r/d | 1–2 | 2–3 | 3–5 | 5–8 | 8–12 | 12–16 | 16–20 | 20–25 |
|---|---|---|---|---|---|---|---|---|
| refine 4 | +5.6 % | +5.3 % | +7.5 % | +9.6 % | +7.8 % | +1.2 % | −6.7 % | **−23.8 %** |
| refine 1.414 | −2.6 % | −1.9 % | −3.2 % | −3.4 % | −4.5 % | −6.4 % | −8.4 % | −8.5 % |
| refine 1 | −6.3 % | −4.0 % | −4.5 % | −3.6 % | −4.5 % | −5.6 % | −6.2 % | −5.8 % |

The refine-4 far field falls a quarter below its neighbour right against the outlet
(r/d = 25), and the [2, 20] fit's intercept is dragged by that tilt.

**The closure is not a single power law, on any grid.** Refitting each grid on sub-bands:

| band (r/d) | `n_t` on refine 0.5 / 1 / 1.414 / 2 / 4 |
|---|---|
| [2, 8] | 1.716 / 1.732 / 1.749 / 1.732 / 1.688 |
| [5, 12] | 2.003 / 2.017 / 2.037 / 2.000 / 2.050 |
| [12, 20] | 2.102 / 2.104 / 2.133 / 2.086 / **2.391** |

The local exponent rises outward by ~0.4 on every grid, so `C_P` and `n_t` are properties
of the fit band as much as of the flow. The design evaluates the closure at r/d ≈ 21
(`DOMAIN_REFIT_2026-09-08.md`), in the region where the exponent is largest. And refine 4's
anomaly is confined to the outer band (2.39 against 2.09–2.13 on the other four grids),
while the inner band [2, 8] gives the only monotone `C_P` triplet in the study
(refine 4/2/1: R = +0.34, p = 1.54, in band) — with `n_t` still diverging there, so this
is not a rescue.

**Two candidates for refine 4's outer band, and the literature does not pick.**
(i) An outflow-boundary effect that reaches inward only at high resolution (the R = 50 mm
domain test, GB1, cleared the boundary at refine 1 only). (ii) The finer grid is *closer*
to the self-similar turbulent wall jet. `LIT_WALLJET_DECAY_2026-09-28.md` (DOIs resolved):
the far-field law Phares et al. 2000 themselves adopt (their Eq. 2.37, from Poreh, Tsuei &
Cermak 1967) is tau_w ∝ (r/H)^-2.3, i.e. **tau/(rho U0²) = 0.0367 (r/d)^-2.3** at
Re_d = 1e4, H/d = 6, ~15 % on magnitude. So a local exponent rising outward toward ~2.3 is
*expected*, and refine 4's 2.39 is nearer that asymptote than the other grids' 2.09–2.13.

**An external magnitude anchor the paper does not yet use** (`convergence/explore_poreh_farfield.txt`),
CFD / Poreh-law, median per r/d bin:

| grid | 5–8 | 8–12 | 12–16 | 16–20 | 20–25 | 25–35 | 35–48 |
|---|---|---|---|---|---|---|---|
| refine 0.5 | 0.826 | 0.938 | 1.012 | 1.064 | 1.095 | – | – |
| refine 1 | 0.827 | 0.934 | 1.007 | 1.058 | 1.097 | – | – |
| refine 1.414 | 0.830 | 0.929 | 0.993 | 1.036 | 1.069 | – | – |
| refine 2 | 0.863 | 0.982 | 1.070 | 1.130 | 1.165 | – | – |
| refine 4 | 0.949 | 1.055 | 1.073 | 1.050 | 0.902 | – | – |
| SO refine 1 | 0.828 | 0.931 | 1.002 | 1.050 | 1.088 | – | – |
| R = 50 mm, refine 1 | 0.834 | 0.946 | 1.023 | 1.076 | 1.121 | 1.178 | 1.238 |

On every grid the wall-jet shear is within 0.83–1.17 of the published law over r/d 5–25,
mostly inside its ~15 %: **the first external check on the magnitude of the closure's
wall-jet region.** Two cautions. The CFD decays more slowly than r^-2.3 on every grid
except refine 4 (the ratio climbs outward, to 1.24 at r/d 35–48 on the wide domain), so the
agreement is in magnitude over a band, not in exponent. And Poreh's own r/d range and
Reynolds range are unverified (paper paywalled), so Re_d = 1e4 may lie outside it.
Deciding (i) against (ii) needs refine ≥ 3 on the 50 mm domain: if the refine-4 bend at
r/d 20–25 moves outward with the boundary, it is (i).

## 5. What would settle the rest

*Pending the SO grades.*
