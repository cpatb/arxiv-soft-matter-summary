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

## 3. The registered scheme test (S1–S6)

*Pending — filled from `convergence/so_triplet_grade.md`.*

## 4. What would settle the rest

*Pending.*
