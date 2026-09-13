# Adversarial review of the transfer-law claim — 2026-09-13

Run through `physics-referee` in a fresh context, seeing only `transfer_law.py`,
`out/transfer_law.txt`, `dst_window.py`'s docstring and
`derivations/dilatancy_sign_mechanism.md` §7–10, plus the acceptance criteria it was
asked to attack. It recomputed the χ² and exponent fit independently (agreeing to
rounding) and pulled the raw `phase2_fits.csv` rows behind each arm.

**The claim it was given:** over the friction range `phi_J`(CV 0) = 0.5590–0.6324, the
CV 0 → 0.5 shift is constant within error; free-volume is disfavoured at p = 0.014 and
relative at p = 0.029 against p = 0.32 for absolute; therefore the shift is a property of
the size distribution rather than of the packing, and +0.053 can be carried to a material
that jams lower.

---

## Findings, and what was done about each

**1 — MATERIAL. The `mu_p` = 1.0 arm's ladders are not trustworthy, and it was carrying
the claim.** `phi_J_lin` is a straight line through four pressure rungs; its standard
error is computed as though that model were right and cannot see a ladder that is not a
line. Checking the raw rows:

| `mu_p` | non-monotone | `r2_lin` range |
|---|---|---|
| 0.0 | 2/6 | 0.66–0.94 |
| 0.3 | 3/6 | 0.77–0.98 |
| **1.0** | **6/6** | **0.002–0.56** |

and one `mu_p` = 1.0 seed has `slope_lin` = **−7.06**, i.e. `phi_J` *falling* as pressure
rises — inadmissible for an equation of state, not merely noisy. This bites specifically
because the feature carrying ABSOLUTE's χ² win was that the shift is **non-monotonic** in
`phi_J` (0.058 → 0.062 → 0.056), a shape no single-sign power law can produce — and that
non-monotonicity is set entirely by the direction of the 0.3 → 1.0 step.

✅ **ACTED ON.** `transfer_law.health()` now refuses any arm with a non-positive slope or
median `r2` < 0.85, applied *before* any number is read; `mu_p` = 1.0 is REFUSED with its
numbers still printed. **The conclusion changed:** on the two healthy arms, ABSOLUTE
(χ² = 1.19, p = 0.28) and FREE-VOLUME (χ² = 1.61, p = 0.21) are **not separable** —
free-volume's rejection came entirely from the bad arm. Only RELATIVE is still rejected
(χ² = 7.02, p = 0.008). The lever is 2 points and 0.058 wide, not 3 and 0.073, the
non-monotonicity is withdrawn, and the exponent is no longer fitted at all (2 points,
2 parameters, zero residual dof); the admissible range is the span of the surviving named
forms, n ∈ [0, 1]. The same health gate was added to `shape_control.py`, where **all
family arms pass** (median `r2` 0.916–0.953, no negative slopes).

**2 — MATERIAL. Friction is collinear with the contact number.** `Z_iso` falls with
`mu_p` over the same range as `phi_J`, so this test cannot separate "the shift tracks
packing density" from "the shift tracks contact topology". Both are equally "a property
of the packing it sits in". The suggested clean lever: move `phi_J` by **protocol** at
fixed `mu_p`, where isostaticity pins `Z`.

✅ **ACTED ON** as a recorded limitation: the claim is narrowed to "rather than of
`phi_J` specifically". The protocol lever is written into the next-steps.

**3 — MATERIAL. The law is fitted on the wrong channel.** The exponent is fitted entirely
on ISOTROPIC ladders and applied multiplicatively to the SHEARED critical state, which
exists at exactly one friction. At the one friction where both exist they differ by
**+0.0088 ± 0.0032, 2.75σ** (isotropic +0.0621 ± 0.0026 vs sheared +0.0533 ± 0.0018 —
independently reproduced). There is no lever at all on the channel the transfer is needed
for.

✅ **ACTED ON** as a recorded limitation. The cheapest fix, already built, is to run
`PREREG_CONSTP_RATEMATCH`'s protocol at `mu_p` = 0 and 1.0 as well as 0.3.

**4 — Minor-material. The χ² table and the exponent bound answer different questions.**
The table compares three *named, motivated* forms at equal parameter count and that
comparison is quotable; the bound is a bound and "n = 0.10" must never be quoted alone.
They are not in tension, because FREE-VOLUME is more than "n = 1".

✅ **ACTED ON**: both are printed with an explicit note, and with two arms the exponent is
not fitted at all.

**5 — Minor. The literature framing was a category error.** Farr & Groot
(doi:10.1063/1.3276799) and Desmond & Weeks (doi:10.1103/PhysRevE.90.022204) give
`phi_RCP` as a convex, accelerating function of polydispersity — but **random close
packing is frictionless by construction**, so nothing in that literature can predict `n`
for a test run *across friction*. The referee also flagged its own read of Desmond &
Weeks' `c2` as **[UNVERIFIED]** (0.857 vs 0.0857 ambiguous; the 0.857 reading extrapolates
to `phi_rcp` > 1 by CV = 0.5, outside the paper's own fitted range) and said not to use it.

✅ **ACTED ON**: the citation is withdrawn as support for the exponent hypothesis and kept
only as qualitative corroboration of the `mu_p` = 0 ladder's convexity — this project's own
secant slope steepens 0.051 → 0.178 per unit CV over CV [0, 0.1] → [0.4, 0.5], a 3.5×
steepening. **No coefficient from either paper is quoted anywhere.**

**6 — Minor. Rolling resistance is not the clean second lever I assumed.** It shares
finding 2's `Z`–`phi_J` collinearity; standard rolling-resistance laws use
`M_max = eta * R_eq * F_n`, introducing an **explicit particle-size coupling** that sliding
friction does not have — a new confound on the one axis friction was chosen to avoid; and
`phase6_cohesion.py` exercises `CohFrictMat` + `momentRotationLaw` only for *bonded*
cohesion with `etaRoll` never touched, so "already proves this codebase can drive it"
overstated what was shown.

✅ **ACTED ON, and the scoping agrees.** The path was implemented behind three gates
(`check_roll_path.py`: default bit-identical; `CohFrictMat` with `etaRoll` = 0 reproduces
`FrictMat` to 0.0042; `etaRoll` = 0.3 lowers `phi` by 0.021, so not inert) — and then a
ladder probe showed **the knob saturates by `etaRoll` ≈ 0.1 and buys only 0.022 in `phi`**,
far short of the `phi_J` ≈ 0.50 assumed. Combined with finding 6's size coupling, rolling
resistance is **not pursued as the transfer-law lever**. The code and its gates are kept
because they are verified and cheap to reuse.

---

## What survives, in the referee's own terms

- **Quotable as stated:** the χ² comparison among the three named forms — but now on two
  healthy arms, where it rejects only RELATIVE.
- **Weaker than written:** "the escape clause is supported over the friction range" — it
  is supported by two arms spanning 0.058, not three spanning 0.073.
- **Not established:** that `phi_J` (rather than `Z`) is the right axis to bound on.
- **Not touched at all:** whether the transfer holds on the sheared channel.

## What this redirected the work to

The review's cumulative effect is that the **friction lever is weak** — two arms, a
collinearity it cannot resolve, and the wrong channel. Meanwhile `shape_control.py`, which
the referee did not see, finds that at **matched CV and matched friction** the three size
distributions give shifts spanning **4.47×** while friction moves the same shift by
**1.07×**. The size-distribution shape is the larger effect by a wide margin, it is
measured on arms that all pass the health gate, and it is free of every confound above
because it varies the independent variable of the project itself.

`PREREG_TAIL_2026-09-13.md` is the follow-up: four lognormals **matched in realised number
CV to 1e-3** and differing only in where the large-particle tail is cut, separating the
number-weighted from the volume-weighted width by 1.72× *by construction*. Its §3 carries
this review's health gate as a registered bar.
