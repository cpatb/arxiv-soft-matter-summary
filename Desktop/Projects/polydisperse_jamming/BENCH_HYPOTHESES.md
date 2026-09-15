# BENCH_HYPOTHESES — what the jamming model now says to expect at the bench

Written 2026-09-14 and extended 2026-09-15 with the K-ladder block (row E, and the
four items it closed in sections 4 and 5), **from the graded verdicts**, not
alongside them. Every number below
is printed by a script named beside it; nothing here is quoted from a plan or a status
card. Where a statement is an inference rather than a measurement, it says so.

Scale convention throughout: `PHI_SCALE` = **2.1e−3**, the protocol's own reproducibility
(`mixture_model.py`), NOT a seed standard error. "20 scale" means Δφ = 0.042.

---

## 1. The five things that are measured

| # | Statement | Evidence | Artifact |
|---|---|---|---|
| A | **The control parameter is a pairwise-filling functional of the PSD, not any descriptor of it.** Two packs matched in number CV (to 5.6e−17), size ratio, `d_max/d_min`, `d95/d05`, span and grain count — differing only in which atom carries the volume — pack at 0.6042 and 0.6677. **Δφ = +0.0635 ± 0.0032, 30 scale, 19.7 se**, where every such descriptor predicts 0. | 123 cells, `FAILURES: []` | `out/mixture.txt` |
| B | **The composition effect is friction-independent.** It survives at `mu_p` = 0 and at two pressures, so the interaction functions are geometric and `beta` carries the friction. | F1–F3 | `out/mixture_mu0.txt` |
| C | **The closure is calibrated at `u` ≥ 1/3 and is not a predictor below it.** Blind at ratios 4 and 5 it under-predicts 9 of 9 compositions, rms 7.50 against a bar of 4.0. | X1 FAIL, 39 cells | `out/grade_mixture_ext.txt` |
| D | 🔑 **The functional is NOT ADDITIVE beyond two classes.** Six ternaries at `d` = 1, 2, 4 — whose three pairs sit at `u` = 1/2, 1/2, 1/4, *inside* the calibration range — are under-predicted **6 of 6**, mean **+7.53**, rms **7.83** against a pre-registered bar of 2.8. **Confirmed at doubled N**: 6/6 positive again, mean **+8.67**. | T1 FAIL + NA1 NON-ADDITIVE, 42 cells, 0 errors | `out/grade_ternary.txt`, `out/grade_ternary_N.txt` |
| E | 🔑 **The miss does NOT accumulate with the number of classes, and the closure's SIGN in K is wrong.** On a geometric ladder of fixed span 4 filled in at K = 2, 3, 4, 5 at matched fines fraction, measured φ **rises** with K where the additive closure says it **falls**, on 3 of 3 levels (2 of 3 resolved against seed error, 1 undetermined). The residual steps once between three and four classes and does not step again — though whether that is SATURATION is undetermined, see qualification 2 below: `R4−R3` = +3.4/+7.7/+8.0, `R5−R4` = +0.2/−1.6/−1.8, i.e. increment ratios −0.22 to +0.14 against the 1.00 / 1.33 / 2.00 that per-class, per-pair and per-triple growth predict. **All three growth laws refuted, and no-K-dependence refuted too.** | G1, G7, G1c; 135 runs at two N, 135 admitted, 0 errors; refereed 2026-09-15 | `out/grade_kladder.txt` |

| F | 🔑 **The filling interaction `a(u)` IS measurable in the gap the K-ladder reads — at `u` = 0.70711, and NOT at `u` = 0.8.** `a(0.70711) = 0.9481 [0.8799, 0.9699]`, `b = 0.9376 [0.8849, 0.9599]`, against an identifiability bar of 0.25 quoted verbatim from the closed block. What made it work is **where in composition the block samples**, not how many seeds: the profile bar is Δχ² = 1 *per point*, so the same point set at n = 8, 16, 32 returns w(a) = 0.0400 three times identically. At `u` = 0.8 the filling signal (max(φ−β) = **2.21 scale**) falls below what this design can resolve (**2.57**), so **no experiment at this design could determine `a(0.8)`** — an inference from a null, putting `a(0.8) ≳ 0.975`. | U0 12/12 bit-identical, U1/U1b/U4/U5 PASS, U1c/U3 FAIL; 108 runs at two N, 108 admitted, 0 errors | `out/grade_ugap.txt` |

**D, E and F are the new ones. E closed §5's first two items; F closes the third and is
what makes E's qualification 2 a measurement rather than a flag.**

⚠ **E must always be quoted with three qualifications**, all printed by
`out/grade_kladder.txt` (sections G1 and G1c) rather than argued:

1. **Most of the step is the MODEL.** `R = (φ − ADD)/scale`, and at `y_S` = 0.05 and 0.08
   **79–93 % of the K3→K4 step is the additive closure's own prediction dropping by
   ~6.4 scale**, not the packing rising (φ moves only +0.5 to +1.6). Only at `y_S` = 0.03
   is it mostly a real rise (+2.5, +3.6 against se 1.2, 0.4).
2. **"And then saturates" is UNDETERMINED — and as of 2026-09-16 that is MEASURED, not
   merely flagged.** The gap above `u` = 0.667 now has a knot (row F), and re-running the
   comparison with it gives `R5−R4` = **+2.70 / +0.69 / −0.26** at N = 1500 and
   **+3.13 / +1.47 / +0.30** at N = 3000 — ≤ 0 on **1 of 3** and **0 of 3**, against the
   pre-registered bar of ≥ 2 of 3 at both, so **the registered gate U3 FAILS**. But swept
   over the new knot's **own profile interval** the saturation verdict survives on
   **33 of 87 jointly admissible corners (37.9 %)**. So the honest verdict is
   **undetermined, leaning against saturation**: it is no longer quotable as established,
   and it is not refuted either. (A referee corrected a first version that swept the
   MARGINAL box and reported 110/169 = 65 %; 82 of those corners are (a, b) pairs the
   joint likelihood excludes. Both numbers stay printed.) K = 2, 3, 4 pairs all sit at `u` ≤ 0.630, so **the K3→K4 rise is
   unaffected** except through interpolant leakage (`out/grade_ugap.txt`, prereg C-iii).
3. **A step in ADD across K is not by itself a filling statement.** `packing_density` is a
   minimum over branches, and the argmin reassigns somewhere on the K axis at **3 of 3**
   `y_S` levels.

**The decomposition-free half of E is the sign test**, which needs no residual
bookkeeping and no bar: measured φ rises with K where the closure says it falls, at every
level and both N — **2 of 3 levels resolved against their own seed error, 1 undetermined**
(`y_S` = 0.08 is +0.80 against se 0.85). That is the actionable defect.

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

* **How the non-additivity grows past K = 3.** **HALF closed by E**: it does not grow in
  the way any of the three registered laws predicts — per-class, per-pair and per-triple
  growth are all refuted, and so is no-K-dependence. ⚠ **Whether it then SATURATES is
  UNDETERMINED and stays open** (E qualification 2): with the U-GAP knot in hand the
  registered gate U3 FAILS, while the same knot's own profile interval leaves saturation
  standing on 65 % of its corners. Still one friction and one pressure, still only
  measured to K = 5 — and **K = 8 or 20 is now known to be out of reach on a span-4
  ladder**, because those rungs read the closure at `u` ≥ 0.82 where §5 item 1 shows it
  cannot be calibrated.
* ~~⚠ **The design is 2-LEVEL in the fines fraction.**~~ **CLOSED 2026-09-15 by G2**:
  widened to **five** levels spanning **6.0×** (`y_S` 0.03 → 0.18, reaching above the old
  maximum of 0.15 and 3.3× below its minimum). The trend survives at K = 3 —
  spearman +0.700 (p 0.188) at N = 1500 and **+0.900 (p 0.037)** at N = 3000, span +8.5
  and +11.0 scale — and it is **absent and of the opposite sign at K = 2**
  (−0.600, −0.500), where the closure is in-sample. The N = 1500 correlation lands
  *exactly* on its own 0.70 bar, so the claim rests on the N = 3000 replicate and must be
  quoted that way.
* **That the sequential rule is correct** — and 2026-09-15 put a number on how far from
  correct it is. With its diameter convention **swept** over all four moment ratios
  rather than chosen, the best nested rule beats the additive one by **1.93×** on the
  K-ladder, against a **1.96× spread across the conventions themselves**: the undeclared
  choice is worth as much as the result. On the held-out continuous families the spread
  is **11.9×** (χ² 42.8 to 509.3), the best convention closes +0.89 of the additive-to-J
  gap and the worst opens it by 5.48×. **The comparison does not settle.** And the
  registered BRACKET sub-clause FAILED at 74 %/74 % against ≥ 80 %.
* **Anything at `u` < 1/4.** A `d` = 1, 2, 8 ternary would reach there and would
  re-confound additivity with extrapolation, which is why it was deliberately not run.
* **Ratio 5 vs ratio 4**, still unresolved at n = 3 (X4): the predicted separation sits
  inside the three-seed standard error at two of three points.
* **Any rheological consequence.** See H4. The project measures packing, not flow.
* ~~**A residual N-dependence on one composition.**~~ **CLOSED 2026-09-15 by G5, and it
  closes the opposite way round from how it was posed.** Re-run at **six** seeds at both
  N, (0.55, 0.35, 0.10) shifts +2.19 scale with a standard error *of the difference* of
  0.29 — **7.6 se**, Welch t = 7.64, p = 0.0001 — so the registered prediction that six
  seeds would show it as seed noise **FAILED; the shift is real.** But it is **not
  anomalous**: its own shift is **+1.89 against a block-wide mean of +1.87 over 18 cells**
  (rms 2.04, range +0.14 to +3.84). **Every cell in this campaign moves by about +1.9
  scale when N doubles**, so the quantity that needed the second N was the block-wide
  offset, not that composition. Seeds 1–3 reproduced the closed ternary block
  **bit-identically** (6/6 including the monodisperse control, `diff = +0.000e+00`), so
  it is the same packing and not a near one.
* **The residual's absolute level at any single N**, precisely because of the previous
  point: R shifts by ~+1.9 scale on 18 of 18 cells between N = 1500 and 3000. Increments
  and contrasts at fixed N are the robust quantities; a single R is not. (An artefact
  would *shrink* with N, which is why the direction is reassuring rather than alarming.)
* 🔴 **The grain-count "null" is a CONFOUNDED DIAGNOSTIC and its 0.6 bar is WITHDRAWN**
  (physics-referee, 2026-09-15). `n_min` is collinear with **both** design axes by
  construction — spearman(`y_S`, `n_min`) = **−1.0000 on all 8 (K, N) combinations**, and
  `n_min` also falls with K at fixed `y_S` — so a pooled spearman(R, `n_min`) is
  guaranteed large and negative whether or not any leak exists. The measured **−0.4694
  (p 0.0067)** is a restatement of the grid and testifies to nothing; an earlier draft
  read it as a partial residual leak. `d_max/L` is not collinear that way and sits at
  +0.1400 (p 0.445). **The valid finite-size control is the gated N-doubling comparison**,
  which passes 3 of 3, plus the direction of the N-shift.
* **The `d_max/L` margin**, though it is no longer the open question it looked like. The
  solver's own record gives 0.1997–0.2256 at N = 1500 and 0.1586–0.1791 at N = 3000, all
  below the 0.25 flag, with no cell carrying a note in either block — and the residual
  does not track it once both N are pooled (+0.31, p 0.33).
* **N > 3000.** Both blocks reached N = 3000 inside the 600 s agent cap, which the
  project's own notes had assumed impossible above 1500; a ternary cell is cheap. Nothing
  here bounds behaviour above 3000.

## 5. The cheapest next experiments, in the order that buys the most

~~1.~~ **DONE 2026-09-15** — the K-ladder block (K = 2–5 at fixed span 4, `y_S` widened
   6×, 135 runs at two N). Answer: the miss does **not accumulate** with K — it steps once
   between three and four classes, and the closure's **sign** in K is wrong. ⚠ **Whether
   it then SATURATES is undetermined** (E qualification 2; this line said "saturates"
   until 2026-09-16 and contradicted that qualification in its own document).
   See E and `CLAIMS_KLADDER_2026-09-15.md`.
~~2.~~ **DONE 2026-09-15** — the sequential rule's convention swept, and the transfer
   re-run: **it does not settle**, because the convention spread (11.9×) exceeds the
   effect. See `out/transfer_kladder.txt`.

~~1.~~ **DONE 2026-09-16, and it came back HALF YES** — the U-GAP block, 108 runs at two
   N (`CLAIMS_UGAP_2026-09-16.md`). `a(0.70711)` is measured (row F); `a(0.8)` is **not
   measurable on this protocol at any seed count**, because the filling signal there is
   below the design's resolving power. The K-axis blocker is therefore **lifted at
   `u` = 0.70711 only**, which is exactly what K = 5 needs — and **not** what K = 8 or 12
   would need.

1. 🔴 **THE K AXIS IS UNBLOCKED AT K = 5 AND STILL BLOCKED ABOVE IT.** A span-4 ladder's
   adjacent-pair ratio goes 4^(−1/7) = **0.82** at K = 8 and 4^(−1/11) = **0.88** at
   K = 12, and the U-GAP block shows the closure cannot be calibrated at `u` ≳ 0.8 by this
   method — the signal is gone. So **extending K on a span-4 ladder walks into a region
   that is now known to be unmeasurable**, and item 2's old plan is not available as
   written. The two ways out, in cost order: **(a) widen the SPAN** so that K = 8 keeps
   its adjacent pair near 0.7 (span 4^(7/4) ≈ 11 puts it at 0.707, at the price of
   breaking the `u` ≥ 1/4 constraint row C says the closure needs — so span and K trade
   against each other and that trade is the design question); **(b) a different observable
   at `u` → 1** with more signal than φ − β, which is what fails there.
2. **Settle `a(u)` between 0.70711 and 1 by a means that is NOT a binary composition
   curve**, or accept that the closure is an extrapolation there and quote it as one.
   `out/diag_ugap.txt` shows that region still moves the K = 5 prediction by up to
   **9.59 scale**, so it is not a cosmetic gap — it is just not reachable the way the
   knots below it were.
3. **A continuous PSD at a composition the filling model puts far from the families'
   locus** — the existing lognormal/uniform/bidisperse ladders all run through the
   fines-poor corner, which is what hid the mechanism for a month. (Unchanged, and now
   the highest-value item that is not a ladder extension.)
4. **A block designed to separate the descriptors the post-hoc scan could not.**
   `out/diag_kladder_descriptor.txt` ranks CV (−0.88), `fill_weight` (+0.88),
   `pair_weight` (+0.87) and `triple_weight` (+0.87) above the sample-size null
   (0.70–0.73), while K, `n_pairs` and `n_triples` score *below* it. But
   Σ_{i<j} y_i y_j = (1 − Σ y_i²)/2 exactly, so the top four are **one quantity — how
   evenly the volume is spread across classes — seen four ways**, and this grid cannot
   tell them apart. That is a hypothesis selected on the cells that motivated it, not a
   result; the test places cells where the candidates disagree, which is exactly the move
   that turned the bidisperse negative into the filling functional.
5. **The closure's wrong sign in K is the actionable defect** (E). Any repair should be
   graded on the sign test first, because it needs no residual bookkeeping and no bar.
