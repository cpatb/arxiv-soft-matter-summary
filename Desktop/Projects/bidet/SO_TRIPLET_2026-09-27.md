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

### 3a. What C1 found first: a residual-certified case still drifting by percents

FO refine 1.5, mapped from refine 1 and stopped when every non-`Uz` residual fell below
1e-5 (iteration 189), continued 2400 iterations with the gate quantities read every 400:

| iteration | `n_t` | `C_P` | `tau_w`(15 mm) | peak `tau_w` |
|---|---|---|---|---|
| 191 (the H1 stop) | 1.9089 | 1.3784 | 5.866 | 391.6 |
| 400 | 1.9109 | 1.3808 | 5.848 | 394.8 |
| 1200 | 1.9158 | 1.3771 | 5.764 | 402.9 |
| 2000 | 1.9259 | 1.3944 | 5.681 | 412.3 |
| 2591 | 1.9304 | 1.4060 | 5.658 | 417.3 |

**Still moving at the end** (`C_P` ~0.8 % per 400 iterations): NOT CONVERGED under C1, and
continued (Amendment 5). The residuals sat below 1e-5 the whole time. Three consequences:

1. **Every value read at a residual stop after a mapped start is a transient.** The
   registered S6 (+400, 0.5 %) would have failed this case on two of four quantities, so the
   pre-registration was not blind to it, but `C_P` alone would have passed.
2. **The drift is the size of the effects being measured.** The refine-1 SO/FO difference
   in `C_P` is 5.2 % and grid steps are 1–4 %; a 1.8 %-and-counting iterative drift is not
   small against either.
**Continued to 7200 (`convergence/fo_r1p5_trajectory.txt`), it is not monotone:** `C_P`
rises to 1.4267 at 4800 and then falls (1.4177 at 7200), `tau_w`(15 mm) bottoms at 5.652
(3200) and climbs (5.759 at 7200), still ~0.1 %/400 at the end — a very slow mode, possibly
a limit cycle, of amplitude ~0.5 % in `C_P`. By contrast the cold-started refine-1 FO case
(L1) and SO refine 1.5 both go flat by ~3600 iterations. Not converged under C1.

3. **The legacy grids were certified by the same residual criterion** (cold starts, which
   take longer to get there, so the error is likely smaller but is unmeasured). L1 in
   Amendment 5 continues the refine-1 case to find out. If it drifts, GRID_STUDY's
   iterative-convergence finding (§1 there) was itself incomplete.

### 3b. 🔴 The legacy grids were not iteratively converged either (L1, L2)

Registered in Amendment 5 before either ran: continue the legacy first-order grids that
GRID_STUDY certified by `residualControl`, and read the gate quantities every 400 iterations.

| case (GRID_STUDY status) | residual stop → last read | `C_P` | `n_t` | `tau_w`(15 mm) | state at the end |
|---|---|---|---|---|---|
| refine 1, `case_Re10000_H6` ("converged at 2261", the manuscript's grid) | 2261 → 4661 | 1.3518 → **1.3663** (+1.08 %) | 1.9023 → 1.9102 | 5.866 → **5.792** (−1.26 %) | **flat** from ~3600 |
| refine 2, `r2_ext` ("converged at 5105") | 5105 → 7505 | 1.3928 → **1.4296** (+2.64 %) | 1.8922 → 1.9243 | 6.241 → **5.865** (−6.03 %) | **still moving** (~0.3 %/400 in `C_P`) |

`convergence/L1_legacy_r1_trajectory.txt`, `L2_legacy_r2_trajectory.txt`.

**This retracts GRID_STUDY §1's "the iterative-convergence finding stands, and it is the
most solid thing here."** That finding moved the refine-2 grid from a capped state to a
`residualControl` stop and called it converged; the residual stop is itself ~6 % short on
`tau_w`(15 mm). Consequences, each now to be re-read:

- **"`tau_w`(15 mm) is non-monotone in the grid"** (6.1374 / 6.2415 / 5.8658 on refine
  4 / 2 / 1) rested on refine 2's 6.2415, which is already 5.865 and still falling;
  converged refine 1 is 5.792. The non-monotonicity is not established.
- **The `C_P` "divergence"** (1.3518 → 1.3928 → 1.5809, |R| = 4.59) now has refine 1 at
  1.3663 (converged) and refine 2 at ≥ 1.4296 and rising; refine 4 (491k cells, certified
  the same way) is unmeasured tonight. With the two coarser members moved, the ratio that
  made it a divergence is no longer known.
- **The manuscript's `C_P` = 1.3518 carries a +1.1 % iterative error**; the design `Q*`
  moves with `C_P^-1.426`, i.e. about −1.5 %.

**Independently verified** (numerics-verifier, 28 Sep ~03:40, own reader of the raw
`wallShearStress_*.raw` samples with no shared code): every number above reproduces
exactly; the continuations are clean restarts (identical mesh by md5, byte-identical
`fvSchemes`, relaxation factors and boundary conditions, start fields bit-identical to the
base cases' final fields; only `residualControl` and I/O cadence differ); both bases
genuinely logged "SIMPLE solution converged"; and every residual, `Uz` included, stayed
below 1e-5 throughout the continuations. The drift is real iterative non-convergence under
a criterion that was genuinely met, not an artefact of how it was continued.

The mechanism is the same as §3a: a normalised residual reaches 1e-5 while a slow mode in
the wall-jet solution is still evolving. It shows up in `tau_w`(15 mm) and `C_P` (the
wall-jet quantities) and hardly in peak `tau_w` (stagnation: 441.2 → 440.2).

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
