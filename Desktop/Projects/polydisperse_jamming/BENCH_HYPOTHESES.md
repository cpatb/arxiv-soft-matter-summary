# BENCH_HYPOTHESES — what the jamming model now says to expect at the bench

Written 2026-09-14, **from the graded verdicts**, not alongside them. Every number below
is printed by a script named beside it; nothing here is quoted from a plan or a status
card. Where a statement is an inference rather than a measurement, it says so.

Scale convention throughout: `PHI_SCALE` = **2.1e−3**, the protocol's own reproducibility
(`mixture_model.py`), NOT a seed standard error. "20 scale" means Δφ = 0.042.

---

## 1. The four things that are measured

| # | Statement | Evidence | Artifact |
|---|---|---|---|
| A | **The control parameter is a pairwise-filling functional of the PSD, not any descriptor of it.** Two packs matched in number CV (to 5.6e−17), size ratio, `d_max/d_min`, `d95/d05`, span and grain count — differing only in which atom carries the volume — pack at 0.6042 and 0.6677. **Δφ = +0.0635 ± 0.0032, 30 scale, 19.7 se**, where every such descriptor predicts 0. | 123 cells, `FAILURES: []` | `out/mixture.txt` |
| B | **The composition effect is friction-independent.** It survives at `mu_p` = 0 and at two pressures, so the interaction functions are geometric and `beta` carries the friction. | F1–F3 | `out/mixture_mu0.txt` |
| C | **The closure is calibrated at `u` ≥ 1/3 and is not a predictor below it.** Blind at ratios 4 and 5 it under-predicts 9 of 9 compositions, rms 7.50 against a bar of 4.0. | X1 FAIL, 39 cells | `out/grade_mixture_ext.txt` |
| D | 🔑 **The functional is NOT ADDITIVE beyond two classes.** Six ternaries at `d` = 1, 2, 4 — whose three pairs sit at `u` = 1/2, 1/2, 1/4, *inside* the calibration range — are under-predicted **6 of 6**, mean **+7.53**, rms **7.83** against a pre-registered bar of 2.8. **Confirmed at doubled N**: 6/6 positive again, mean **+8.67**. | T1 FAIL + NA1 NON-ADDITIVE, 42 cells, 0 errors | `out/grade_ternary.txt`, `out/grade_ternary_N.txt` |

**D is the new one and it is the reason this document exists.**

### D survived the confound that would have killed it

Grading T1 surfaced a defect the original design could not break: at fixed N, composition
and coarse-phase **sample size** are collinear by construction (number fraction goes as
`y_i/d_i³`), and the six residuals were **perfectly** rank-correlated with both of the
sample-size variables — `spearman(residual, n_L) = −1.0000` and
`spearman(residual, d_max/L) = +1.0000`, on 6 of 6. That is exactly what a finite-size
artifact looks like, and it made T1's verdict ambiguous.

`PREREG_TERNARY_N_2026-09-14.md` breaks it by holding every composition **fixed** and
doubling N to 3000, which doubles `n_L` and multiplies `d_max/L` by 2^(−1/3) while moving
no composition variable at all. NA1 was written so it could **withdraw** D.

| | N = 1500 | N = 3000 |
|---|---|---|
| mean residual | +7.86 | **+8.67** |
| sign structure | 6/6 positive | **6/6 positive** |
| `spearman(res, n_L)` *within* block | −1.0000 | −1.0000 |
| `spearman(res, n_L)` **pooled over both** | — | **−0.3077** (p 0.33) |
| `spearman(res, d_max/L)` **pooled** | — | **+0.3077** (p 0.33) |

* **The residual does not shrink when the coarse phase doubles — it grows slightly.** A
  finite-size artifact must shrink.
* **The pooled correlations collapse from ±1.0000 to ±0.31.** The N = 3000 points sit at
  the *same* residuals with *double* the `n_L`, which is what destroys them. `n_L` and
  `d_max/L` are excluded as explanations at the registered bar of 0.6.
* **`y_S` survives** at +0.878 (p 0.021) within the N = 3000 block alone. The residual
  tracks the **fines volume fraction**, which is composition.
* **NA2: 5 of 6** compositions agree between the two N within 2 se of the difference. The
  exception, (0.55, 0.35, 0.10), moved **+1.81 against an se of 0.48** — reported, not
  explained; it is the one composition where a residual N-dependence is not excluded.

*(Reconciling the two artifacts: `grade_ternary_N.py` grades each N against that block's
own monodisperse control, `beta` = 0.593985 at N = 1500, while `grade_ternary.py` uses the
isotropic lognormal CV 0 `beta` = 0.594687. The 0.33-scale difference between the two
betas is the whole difference between +7.53 and +7.86. Both choices are stated where they
are made; the N-block re-reads `beta` per block deliberately, because using one `beta`
across two N would fold a finite-size shift in the reference into the quantity under
test — that shift is +0.64 scale and is printed.)*

## 2. Why D is about additivity and not about calibration

The pre-registration named two explanations for the closure's trouble on continuous
PSDs — that `a(u)` is mis-calibrated, or that the pair *sum* is the wrong combination
rule — and the ternary was built to separate them. It does:

* **`a(1/2)` would have to be 0.589–0.661**, against a binary-measured **0.7854 with
  profile interval [0.7500, 0.7999]**. All six lie outside, by 4–7 interval widths. This
  interaction enters both the L–M and M–S terms, sits in the binding branch, and moves φ
  across its whole range (0.701 → 0.644 as `a` goes 0.3 → 0.9), so its root is a real
  constraint.

⚠ **A second inversion was withdrawn on referee review (2026-09-14).** An earlier draft
called "no admissible `a(1/4)` reproduces the measurement, 6/6" *the stronger half of the
result*. **That was backwards.** Asked what φ actually does as `a(1/4)` moves, the answer
is **nothing**: φ is **bit-identical** (0.65617225…) for every value from −3 to +0.6,
because the L–S pair term never enters the binding branch at these compositions. "No
root" therefore says `a(1/4)` is **not identifiable here** — a statement about branch
structure, not about filling. A degenerate parameter cannot testify either way, and the
`a(1/2)` half alone carries the evidence.

So no recalibration of the one interaction that *is* identifiable can absorb the residual.
(Honest limit: this inverts one interaction at a time, so a joint move of both is not
excluded. Evidence, not proof.)

Three further checks that could have killed D and did not:

* **Wiring.** Emptying one class of the K = 3 spec reproduces the corresponding binary to
  **1e−12 on 14 of 14** cases, across both ratios and all three emptyings. The residual is
  not plumbing. (My first version of this check was itself wrong — it asked a
  middle-emptied spec to match a ratio-2 binary when emptying the middle leaves ratio 4.)
* **Ceiling.** All six sit **between** the additive prediction and the ideal-separation
  limit, only **4.6–8.3 %** of the way up that gap. The filling picture is intact; it is
  the combination rule that is at issue.
* **Protocol.** The monodisperse control through the new PSD code path lands **−0.33
  scale** from `beta` (T0 PASS, bar 2.0).

**Size, against the effect it is a residual on:** the additive model captures **73–87 %**
of each ternary's composition effect; the miss is **13–27 %** of it. D is a real
correction, not a refutation of the whole picture.

**Direction:** the model always under-predicts. Three classes acting together fill
**better** than the sum of their separate pairwise effects — the same sign as the
continuous-family failure that prompted the block (−6.9 scale at CV 0.5, held-out
χ² 108.1 against J's 34.9).

**The truth is BRACKETED, and that is a cleaner statement than "sequential wins".**
A sequential (nested-binary) construction — coarse pair packs first, fines fill the
composite — using the *same* measured interactions and **no new free parameter** gives
rms **3.80** against additive **7.83**, ratio **2.06**, above the pre-registered
informativeness bar of 1.6.

But rms only says which candidate is *nearer*, not that the winner is right, and the sign
structure says more: **the additive rule is below the measurement on 6 of 6 and the fully
sequential `d43` rule is above it on 6 of 6.** The measurement is bracketed between them
on every composition. So the defensible reading is **not** "the combination rule is
sequential" — the rms winner still carries the *opposite* bias, of comparable size, to the
rms loser. The truth is an **intermediate, partial-synergy** rule: additive under-fills,
fully sequential over-fills. (Raised by the physics-referee pass, 2026-09-14; an earlier
draft read the rms win as identifying the mechanism.)

## 3. What this predicts at the bench

These follow from A–D. Each is stated so it could come out the other way.

**H1 — Two powders sieved to the same CV can jam 0.06 apart in φ.** This is measurement
A, not an inference. At the drop-impact programme's working φ ≈ 0.44, a 0.06 shift in
φ_J is larger than any φ tolerance the experiment holds. *Bench consequence:* quoting a
suspension's polydispersity as a single width does not specify its jamming point, and two
batches matched on a Mastersizer CV are not matched.

**H2 — Blending to raise φ_J must target the fines' VOLUME fraction, not their number
fraction.** Number fraction goes as `y_i/d_i³`, so at a 4:1 ratio a 60 % *volume* coarse
fraction is 5 % of the *grains*. The composition maximum sits near `X_s` ≈ 0.25–0.30 by
volume, which is ~97 % of the grains. *Falsifiable:* a blend targeted by number count
will land far off the optimum in a predictable direction.

**H3 — A real (continuous) powder will jam HIGHER than a binary-calibrated estimate, by
of order 10–25 % of the composition effect.** This is D read forward. It is the sharpest
bench-testable consequence of non-additivity and it has a sign. *Falsifiable:* if a
measured continuous-PSD φ_J came in *at or below* the pairwise estimate, D would not
generalise past ternaries.

**H4 — The load is carried by the coarse skeleton and the fines are largely spectators.**
In the ternaries, **87–96 % of the small grains are rattlers** and `z_S` = 0.37–0.97,
while the coarse class has `z_L` = 8.0–12.0 and a rattler fraction of **0.000–0.008**;
the middle class is genuinely intermediate in both, on 6 of 6 (`z_S < z_M < z_L`,
`rat_S > rat_M > rat_L`). *Bench consequence, and this one is an INFERENCE not a
measurement:* a rheological property set by the force network (yield stress, the onset of
thickening) should track the coarse sub-pack, while φ — a volume count — reads the fines
too. **The two need not move together, and a bench test that assumed they did would be
measuring the wrong thing.** This project has not measured any rheological quantity on
these packs and cannot settle it.

**H5 — The closure under-predicts Δφ_c on BOTH protocols by 13–34 %, and it does not
know which protocol it is on.** Re-basing it on the constant-pressure shear protocol's own
`beta` (0.591127) captures **87 % and 75 %** of the measured Δφ_c at CV 0.2 and 0.5.

⚠ **The isotropic comparator was at the WRONG PRESSURE in an earlier draft** (found by
the physics-referee pass, 2026-09-14). It used the p_hat = 1e−3 ladder — the pressure the
closure is *calibrated* at — while the shear protocol runs at p\* ≈ 3e−4, and a
pressure-matched isotropic ladder already exists in this campaign
(`grade_constp_prod.phi_iso_at_p(3e-4)`, used by `grade_constp_rm.py`'s own R4). On the
matched ladder the captured fraction is **77 % and 68 %**, not the 66 % and 73 % quoted
before, which changes the CV 0.2 comparison from a 21-point protocol gap to a 10-point
one.

The matched ladder also runs five rungs, and the fraction captured is **strikingly
constant**: 73 / 77 / 66 / 71 / 68 % at CV 0.1 / 0.2 / 0.3 / 0.4 / 0.5. **The closure
under-predicts Δφ_c by a roughly constant ~30 % across the whole isotropic ladder** —
which is a more useful statement than any single rung, and the same direction as the
ternary's non-additivity residual.

⚠ **That is a weaker result than it first reads, and the weakening is measured.** The
predicted Δφ is almost independent of `beta`: `d(Δ_pred)/d(beta)` = **+0.003 and +0.006**,
so moving `beta` by the whole 0.0036 between protocols changes the prediction by
**1e−5**, which is **0.02 of the measurement's own standard error**. The comparison
therefore **cannot fail through the re-basing**. What it actually is: **one fixed
prediction against two different measurements** (isotropic Δ = +0.0144 vs shear +0.0109
at CV 0.2), with the prediction sitting between them and below both. The differing
"fraction captured" is a property of the two *denominators*, not of the model.

*Bench consequence:* a closure calibrated by isotropic compression can be used to RANK
shear-protocol blends, not to predict their absolute φ_c — and a real test of whether the
closure's SHAPE is protocol-dependent has to move something the prediction is sensitive
to, which `beta` is not. (`out/transfer_shear.txt`. A transfer **check**, not a
pre-registered gate — both datasets were on disk before it was written.)

## 4. What is NOT established, and would be quoted wrongly if it were

* **How the non-additivity grows past K = 3.** One ladder, one friction, one pressure.
  The size of the miss at K = 3 is measured; nothing here says it accumulates.
* ⚠ **The design is 2-LEVEL in the fines fraction.** `y_S` takes only **two** values
  (0.10 and 0.15) across the six admitted compositions, inside a corner with
  `y_L` ∈ [0.55, 0.70]. The residual tracks `y_S` and the whole headline therefore rests
  on a two-level contrast in the one variable that survives NA3. More fines giving
  proportionally more extra filling is what a genuine three-body synergy predicts, so
  this does not contradict the finding — but "the functional is not additive" is stated
  more generally than a two-level design in a narrow corner supports, and the next block
  should widen `y_S` before the claim is generalised. (Raised by the physics-referee
  pass, 2026-09-14.)
* **That the sequential rule is correct.** It wins under one of its two conventions and
  loses informativeness under the other. It is the better of two named candidates on six
  points, which is not a calibrated closure.
* **Anything at `u` < 1/4.** A `d` = 1, 2, 8 ternary would reach there and would
  re-confound additivity with extrapolation, which is why it was deliberately not run.
* **Ratio 5 vs ratio 4**, still unresolved at n = 3 (X4): the predicted separation sits
  inside the three-seed standard error at two of three points.
* **Any rheological consequence.** See H4. The project measures packing, not flow.
* **A residual N-dependence on one composition.** NA1 excludes finite size as *the*
  explanation, but NA2 finds (0.55, 0.35, 0.10) moved +1.81 between the two N against an
  se of 0.48. Five of six agree; that one does not, and it is not explained.
* **The `d_max/L` margin**, though it is no longer the open question it looked like. The
  solver's own record gives 0.1997–0.2256 at N = 1500 and 0.1586–0.1791 at N = 3000, all
  below the 0.25 flag, with no cell carrying a note in either block — and the residual
  does not track it once both N are pooled (+0.31, p 0.33).
* **N > 3000.** Both blocks reached N = 3000 inside the 600 s agent cap, which the
  project's own notes had assumed impossible above 1500; a ternary cell is cheap. Nothing
  here bounds behaviour above 3000.

## 5. The cheapest next experiments, in the order that buys the most

1. **A four-class pack on the same ladder** (`d` = 1, 2, 4 plus one more inside
   `u` ≥ 1/4). If the miss grows with the number of classes, the combination rule is
   wrong in a way that matters for real powders; if it saturates at the K = 3 value, the
   correction is a constant and can be absorbed.
2. **Calibrate the sequential rule properly**, with its diameter convention swept rather
   than chosen, and re-run the continuous-family transfer. That is the direct test of
   whether D explains the χ² 108 vs 35 that started this.
3. **A continuous PSD at a composition the filling model puts far from the families'
   locus** — the existing lognormal/uniform/bidisperse ladders all run through the
   fines-poor corner, which is what hid the mechanism for a month.
