# CLAIMS — the U-GAP block, 2026-09-16

Pre-registration `PREREG_UGAP_2026-09-16.md`, committed `8267379` **before any cell ran**,
with blind predictions in `out/predict_ugap.txt`. Graded by `grade_ugap.py` →
`out/grade_ugap.txt`, **exit code 1**,
`FAILURES: ['U1c: arm B came back UNDETERMINED', 'U3: 1/3 and 0/3 levels <= 0']`.

**108 runs, 108 admitted, 0 refused, 0 errors** — 81 at N = 1500 across three ratios and
27 at N = 3000 on arm A's identical cell list. Tag `ug` / `ug3k`, matching none of
`mixr*`, `extr*`, `tern*`, `kl*`; `grade_ugap.assert_scope()` checks that structurally.

---

## C1 — **a(u) IS measurable at u = 0.70711, and it is now measured**

**a(0.70711) = 0.9481, profile interval [0.8799, 0.9699], width 0.090**
**b(0.70711) = 0.9376, profile interval [0.8849, 0.9599], width 0.075**

against the registered bar of 0.25, quoted verbatim from
`grade_mixture_ext.tabulated_closure` and asserted unmoved by the grader. **U1 PASS.**
The gap in the closure above u = 0.667 — the one `CLAIMS_KLADDER_2026-09-15.md` C1a(ii)
flagged as the place where K = 5 reads a curve with no data — now has a measured knot at
exactly the u that K = 5's four adjacent pairs sit at.

**This is the answer to the question the block was commissioned for: (a) measurable with a
better-designed binary block — but only at u ≈ 0.707, not at u ≈ 0.8.**

**What made it work was not more seeds.** `out/diag_uident.txt` §3 shows the registered
rule cannot see them: the profile bar is Δχ² = 1 **per point**, so the same r = 1.5 point
set at n = 8, 16 and 32 returns w(a) = **0.0400** three times, identically. Width is set by
curvature, i.e. by **where in composition the block samples**. The design moved the grid to
straddle the model's branch crossover (x\* ≈ 0.356, four points below and five above)
because `a` enters only the coarse branch and `b` only the fine one. The old block put six
points below and two above.

⚠ **A PASS here is not a claim about precision** (caveat C-iv): the bar is in units of
`PHI_SCALE`, not of the seed error.

## C2 — **the block is on the old protocol, and that is checked by bit-identity**

**U0: 12 of 12 runs reproduce the `mixr` block with `diff = +0.000e+00`.** Four (r, x)
cells — (1.25, 0.05), (1.25, 0.20), (1.5, 0.05), (1.5, 0.20) — were put in the grid for
this purpose; `phase12_mixture.cell_key` carries no tag, so an identical key must give an
identical number. Not "close": exactly zero, three seeds each.

## C3 — **the attribution control holds: the grid did not move the answer**

**U1b PASS.** Arm C re-measures the *existing* u = 0.6667 knot on the *new* grid:
a **0.9249 → 0.9336** (Δ +0.0087), b **0.8985 → 0.8972** (Δ −0.0013), both values inside
the union of the two profile intervals. So arm A's number can be attributed to the
**ratio**, not to the composition grid. Without arm C, C1 would not be attributable.

## C4 — 🔴 **u = 0.8 is NOT measurable on this protocol, and the reason is not the one the prereg gave**

**U1c FAIL** — arm B came back UNDETERMINED again (a [0.0001, 0.9949], b [0.0001, 0.9949]),
against a registered prediction of IDENTIFIABLE.

⚠ **Two corrections to the pre-registration, both recorded rather than resolved in the
direction that flatters it.**

1. The prereg predicted the mechanism would be the *escape branch* (a = 1 or b = 1, where
   the closure collapses to the constant φ = β). **It is not.** The new grid **closed**
   that escape at arm B: χ² = **14.12 against a bar of 10.41**, where the old grid had it
   open at 7.09 against 8. The registered sentence was right about the verdict and wrong
   about the reason — and the reason is what the bound rests on.
2. The prereg's "evidence that a(0.8) > 0.98" took the 0.98 from **arm A's** robustness
   table, which is the wrong arm.

**What is supportable, argued from the measured signal and not from a point estimate**
(arm B's fitted b = 0.5003 is 1/(1+e⁰), the Nelder–Mead *start* — on a flat surface the
optimiser never moved, and the 201-grid argmin lands elsewhere; a meaningless point
estimate cannot carry a bound):

* arm B's measured max(φ − β) = **2.21 scale**;
* the noise-free ceiling for this design resolves down to a signal of **2.57 scale** and
  fails at **2.55**;
* so the filling signal at u = 0.8 sits **below what this design can resolve**, and **no
  experiment at this design could have determined a(0.8), whatever it is.** Mapping the
  signal through the ceiling table puts **a(0.8) ≳ 0.975**.

**Status of that number: an inference FROM A NULL**, via the design's own resolving power.
It is not a measurement, it has no interval, and it must never be quoted with one. It is
consistent with the measured a(0.70711) = 0.9481 under monotonicity in u.

## C5 — **the shipped extrapolation over-predicts filling at u\*, marginally**

**U2 PASS at exactly its bar** — the shipped closure's residual at arm A is negative on
**7 of 9** compositions (bar ≥ 7), mean **−0.70 scale**, rms 1.04, max 2.30. Fitted
a(u\*) = 0.9481 against the shipped 0.9398 (**+0.0083**); b 0.9376 vs 0.9177 (**+0.0199**).

⚠ **Three things weaken this and all three are printed beside it.**

1. **7/9 is exactly the bar.** A marginal pass, labelled as one.
2. **The registered sentence was internally ambiguous.** It read "≥ 7 of 9 (sign-test
   p ≤ 0.09)"; 7/9 gives one-sided p = **0.0898** and two-sided p = **0.1797**, so it
   clears the parenthetical only on a reading the sentence did not specify. The prereg
   explicitly named the **count** as the gate, so the count is what is graded — but the
   ambiguity is recorded as a defect, not resolved in the passing direction.
3. **The interval contains the shipped value.** [0.8799, 0.9699] ∋ 0.9398, so U2's second
   half — "the fitted a(u\*) comes out above 0.9398" — is true of the point estimate and
   **not resolved**. This is a trend, not a measured disagreement.

The blind direction prediction came from held-out data: r = 1.25 failed identifiability, so
`tabulated_closure` never put a knot there, and the shipped closure over-predicts it on
**8 of 8** compositions (mean −1.03 scale) against a fitted-knot control at r = 1.5 of
+0.11, 3/8. That prediction was made in `out/predict_ugap.txt` before any cell ran.

## C6 — 🔑🔴 **"AND THEN SATURATES" IS UNDETERMINED, NOT CONFIRMED — and the registered gate FAILS**

This is what the block existed to settle, and the two sentences below are different
sentences. Both are recorded.

**The registered gate: U3 FAIL.** With the measured knot inserted, R(5) − R(4) ≤ 0 on
**1 of 3** `y_S` levels at N = 1500 and **0 of 3** at N = 3000, against a registered bar of
≥ 2 of 3 at both. The K-ladder's C1a sub-clause does not survive at the point estimate.

| N | y_S | R5−R4 old | R5−R4 new | ΔR5 | ΔR4 | Δφ(5−4) | ΔADD(5−4) |
|---|---|---|---|---|---|---|---|
| 1500 | 0.03 | +0.20 | **+2.70** | +1.66 | −0.84 | +0.15 | −2.55 |
| 1500 | 0.05 | −1.58 | **+0.69** | +1.58 | −0.69 | +0.59 | −0.10 |
| 1500 | 0.08 | −1.77 | **−0.26** | +0.99 | −0.52 | −0.72 | −0.46 |
| 3000 | 0.03 | +0.63 | **+3.13** | +1.66 | −0.84 | +0.58 | −2.54 |
| 3000 | 0.05 | −0.80 | **+1.47** | +1.58 | −0.69 | +1.38 | −0.09 |
| 3000 | 0.08 | −1.21 | **+0.30** | +0.99 | −0.52 | −0.17 | −0.47 |

φ is fixed data here, so **every change from old to new is the model**, and it is split
between the two arms because caveat C-iii says both move: ΔR5 = **+0.99 to +1.66** and
ΔR4 = **−0.52 to −0.84**, i.e. the new knot pushes the two arms in **opposite** directions
and roughly two thirds of the swing is K = 5's own.

**The scientific statement: UNDETERMINED.** Sweeping a and b over arm A's **own profile
intervals** — the verdict read with its error bar rather than at a point — saturation
**survives on 110 of 169 corners (65.1 %)**, declining monotonically with a:

```
a0.880:13/13  a0.887:12/13  a0.895:12/13  a0.902:11/13  a0.910:10/13  a0.917:9/13  a0.925:9/13
a0.932:8/13   a0.940:7/13   a0.947:6/13   a0.955:5/13   a0.962:4/13   a0.970:4/13
```

So the measured knot **is not precise enough to decide saturation either way**. The
supportable claim is that C1a's "and then SATURATES" is **undetermined**, not that it is
refuted — and, equally, it can no longer be quoted as established.

⚠ **Not an artefact of arm C.** Recomputed with **arm A's knot only**, keeping the existing
u = 0.6667 knot rather than replacing it with arm C's re-measurement, the verdict is
**0/3 and 0/3** — the two readings agree, so it does not turn on having silently
re-measured an existing knot.

⚠ **K = 2 and K = 3 are untouched**, and so are the K3→K4 step and the G7 sign result: those
rungs' pairs all sit at u ≤ 0.630. They are affected only through the interpolant leakage
of caveat C-iii.

## C7 — **U4 PASS: the interaction functions stay monotone**

| u | a | b | source |
|---|---|---|---|
| 0.00000 | 0.0000 | 0.0000 | proved limit |
| 0.25000 | 0.3010 | 0.4796 | existing |
| 0.33333 | 0.5210 | 0.5505 | existing |
| 0.50000 | 0.7854 | 0.7471 | existing |
| 0.66667 | 0.9336 | 0.8972 | **re-measured (arm C)** |
| 0.70711 | **0.9481** | **0.9376** | **NEW (arm A)** |
| 1.00000 | 1.0000 | 1.0000 | proved limit |

## C8 — **U5 PASS: the identifiability is not a finite-size accident**

At N = 3000 arm A is **IDENTIFIABLE**: a 0.9259 [0.8599, 0.9549], b 0.9114 [0.8449,
0.9349]. Δa(N) = −0.0223, Δb(N) = −0.0262, **not graded** — the prereg declared in advance
that a knot fitted at N = 3000 is not comparable to a closure calibrated at N = 1500,
because the K-ladder's G5 measured a block-wide offset when N doubles. This block's own
figure on the identical cell list is **+1.51 scale**, against the K-ladder's +1.9.

## C9 — **U6, POST-HOC: a knot above u = 0.667 does NOT rescue the transfer, and makes it worse**

Declared post-hoc in the prereg and searched before the block ran
(`out/diag_transfer_ugap.txt`). The continuous families do put most of their pair weight in
the gap — lognormal 0.45–1.00 and uniform 0.80–1.00 of it above u = 0.6667 — so the knot
*could* have mattered. It does not help:

* held-out family χ²: J **34.9**, ADDtab **108.1 → 147.1** with the measured knot — **worse**;
* nested rule with the measured knot: d10 316.8, d21 108.3, d32 **14.9**, d43 430.1 —
  spread **28.95×** against the standing 11.90×.
* swept over the whole admissible knot range beforehand, the convention spread runs
  **2.00× to 29.03×** and never collapses: its narrowest point (2.00×) has a *best* member
  of **427.9**, i.e. the spread narrows because every member got worse.

**The spread narrowing while every member degrades is not a collapse**, and both the ratio
and the best member are printed so it cannot read as one. **G8 stays unsettled**, and the
reason is now ruled out explicitly rather than by argument.

---

## Caveats, carried verbatim from `out/predict_ugap.txt` (asserted by the grader)

* **C-i** One friction (`mu_p` = 0.3), one pressure rung (`p_hat` = 1e-3), one N (1500).
  Nothing here says a(u) is pressure- or friction-independent.
* **C-ii** Arm A's knot is at u = 0.70711 **exactly**. It does **not** fill the gap: the
  interval (0.70711, 1) still has no measured point.
* **C-iii** Inserting a knot moves a(u), b(u) at **neighbouring** u through the
  interpolant's derivatives — up to **3.49 scale** on the K = 4 cells. Any K = 4 number
  recomputed after this block is **not** independent of it.
* **C-iv** The identifiability bar is in units of `PHI_SCALE`, not of the seed error, so a
  PASS is **not** a claim about precision.

## What is NOT established

* That the two-parameter (a, b) **form** is right at high u. Only its parameters were fitted.
* That a(u) between 0.70711 and 1 is anything in particular. It is still PCHIP between a
  measured knot and a proved limit, and `out/diag_ugap.txt`'s sweep shows that region moves
  the K = 5 prediction by up to 9.59 scale.
* Any statement about pressure or friction dependence (C-i).
* That saturation is **refuted**. C6 says undetermined; the registered gate says fail. Do
  not let the gate's wording become the claim.
