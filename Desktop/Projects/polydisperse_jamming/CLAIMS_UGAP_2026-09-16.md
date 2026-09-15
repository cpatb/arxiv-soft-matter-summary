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
a **0.9249 → 0.9336** (Δ +0.0087), b **0.8985 → 0.8972** (Δ −0.0013).

🔴 ⚠ **The registered wording for this gate was VACUOUS and has been replaced** (referee,
2026-09-16). It asked whether both values lie "within the **union** of the two profile
intervals" — but `identifiable()` already guarantees each best fit lies inside its *own*
interval, and each interval is a subset of the union, so **any two individually
identifiable fits pass, even with disjoint intervals an order of magnitude further apart
than these.** It was a gate that could not fire. The vacuous result stays printed beside a
replacement that is stricter in kind: the two intervals must actually **intersect**, and
|Δ| must not exceed the narrower interval's width. Measured: a overlap **+0.0300**
(narrower interval 0.0400), b overlap **+0.0650** (0.0650); |Δa| 0.0087 vs 0.0400, |Δb|
0.0013 vs 0.0650.

**And a third test, stronger in kind again, added on the referee's second pass**: interval
overlap is shape-blind — two intervals can intersect while each fit's optimum sits at a
far edge — so the gate now also runs a **nested-model test** asking directly whether *one*
(a, b) explains *both* grids. Pooling the 8 old and 9 new r = 1.5 points costs
**Δχ² = 0.77 for 2 extra degrees of freedom**, against the project's own 1-per-dof
convention (bar 2.0). **PASS.**

So the conclusion is unchanged — arm A's number is attributable to the **ratio**, not to
the composition grid — but it now rests on a test that could have failed, and on the
strongest of the three.

## C4 — 🔴 **u = 0.8 is NOT measurable on this protocol, and the reason is not the one the prereg gave**

**U1c FAIL** — arm B came back UNDETERMINED again (a [0.0001, 0.9949], b [0.0001, 0.9949]),
against a registered prediction of IDENTIFIABLE.

⚠ **Two corrections to the pre-registration, both recorded rather than resolved in the
direction that flatters it.**

1. The prereg predicted the mechanism would be the *escape branch* (a = 1 or b = 1, where
   the closure collapses to the constant φ = β). **It is not.** The new grid **closed**
   that escape at arm B: χ² = **14.12 against a bar of 10.41** (= χ²_min 1.41 + n 9, the
   same bar the profile itself uses), where the old grid had it open at 7.09 against 8.
   The registered sentence was right about the verdict and wrong about the reason — and
   the reason is what the bound rests on.
   ⚠ That 10.41 was **numerically right but traceable to nothing**: the grader printed
   only "bar = n = 9" and the first claims draft quoted a bar no artifact contained
   (referee, 2026-09-16). `grade_ugap.py` now prints χ²_min, n and the bar in one row for
   every arm, so the figure traces.
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
4. **The 9 compositions are not 9 independent trials** (referee, 2026-09-16). They share
   the same three seed values; measured, the per-seed residual columns correlate at
   r = **+0.726** (seeds 1, 2), −0.225 (1, 3) and −0.003 (2, 3) across compositions. One
   pair is strongly correlated, so the binomial p-values above are **optimistic**. No
   verdict changes — U2 was already marginal on three other grounds — but the count must
   not be quoted as if it were nine independent draws.

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

**The scientific statement: UNDETERMINED.** The verdict read with its error bar rather
than at a point — sweeping a and b over arm A's own intervals on a 13 × 13 grid.

🔴 ⚠ **The first version of this swept the MARGINAL box and the number was wrong**
(referee, 2026-09-16). Two wide marginals enclose corners the **joint** likelihood
excludes — the very point `out/diag_uident.txt` §7 makes about arm B, not applied here
until the referee asked. Of the 169 marginal corners only **87** satisfy
χ² ≤ χ²_min + n = 12.41; the other 82 are combinations the data actively excludes.

| sweep | saturation survives | status |
|---|---|---|
| marginal box (first version) | 110 of 169 — **65.1 %** | **SUPERSEDED**, printed so the correction is visible |
| **joint likelihood region** | **33 of 87 — 37.9 %** | **the number the verdict rests on** |

Survival by a(u\*), jointly admissible corners only:

```
a0.880:2/2  a0.887:2/2  a0.895:3/3  a0.902:2/4  a0.910:2/5  a0.917:2/6  a0.925:2/6
a0.932:2/7  a0.940:3/9  a0.947:3/10  a0.955:4/12  a0.962:3/11  a0.970:3/10
```

So the measured knot **is not precise enough to decide saturation either way**. The
supportable claim is that C1a's "and then SATURATES" is **undetermined**, not that it is
refuted — and, equally, it can no longer be quoted as established. ⚠ At 37.9 % rather than
65 % the balance leans **toward** the graded FAIL, and the honest reading is "undetermined,
leaning against saturation", not "a coin flip".

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
(`out/diag_transfer_ugap.txt`). It does not help:

* held-out family χ²: J **34.9**, ADDtab **108.1 → 147.1** with the measured knot — **worse**;
* nested rule with the measured knot: d10 316.8, d21 108.3, d32 **14.9**, d43 430.1 —
  spread **28.95×** against the standing 11.90×.
* swept over the whole admissible knot range beforehand, the convention spread runs
  **2.00× to 29.03×** and never collapses: its narrowest point (2.00×) has a *best* member
  of **427.9**, i.e. the spread narrows because every member got worse.

**The spread narrowing while every member degrades is not a collapse**, and both the ratio
and the best member are printed so it cannot read as one. **G8 stays unsettled**, and the
reason is now ruled out explicitly rather than by argument.

🔑 **And the mechanism is legible, which makes this the block's most consequential
by-product.** A knot at u\* = 0.70711 sits at the **bottom** of the gap; the continuous
families put their weight at the **top** of it. `mm.psd_classes` discretises each family
into **160 equal-count classes**, and adjacent classes of a continuous PSD are nearly the
same size, so:

| family | CV | w(u > 0.6667) | w(u > 0.70711) | **w(u > 0.90)** | w(u > 0.95) |
|---|---|---|---|---|---|
| lognormal | 0.1 | 0.9965 | 0.9865 | **0.5428** | 0.2797 |
| lognormal | 0.5 | 0.4456 | 0.3809 | **0.1075** | 0.0445 |
| uniform | 0.1 | 1.0000 | 0.9999 | **0.5471** | 0.2973 |
| uniform | 0.5 | 0.7990 | 0.7466 | **0.3346** | 0.1757 |

**Over half the lognormal and uniform pair weight at CV = 0.1 sits above u = 0.90.** The
new knot therefore does not *cover* that weight — it only changes how PCHIP **guesses** on
[u\*, 1], which is why the transfer moved, and moved the wrong way.

🔴 ⚠ **AND THE CONCLUSION I FIRST DREW FROM THAT TABLE IS WITHDRAWN.** The first version of
this section read "**the continuous-family transfer needs a(u) at u > 0.9**". A referee
refuted it on the second pass, and the check is one column — put the weight beside the
residual it is supposed to explain:

| family | CV | w(u > 0.90) | residual (φ_meas − ADDtab)/scale |
|---|---|---|---|
| lognormal | 0.1 | **0.5428** | **−0.1** |
| lognormal | 0.3 | 0.1915 | +3.1 |
| lognormal | 0.5 | **0.1075** | **+6.9** |
| uniform | 0.1 | **0.5471** | **−0.2** |
| uniform | 0.5 | 0.3346 | +2.4 |

**spearman(|residual|, w(u > 0.90)) = −0.939 (p 0.000)** over the 10 continuous cells.
The cells that lean hardest on u > 0.9 are **already fit to within noise**; the badly
misfit ones (CV 0.4–0.5, up to +6.9 scale) have the **least** weight there. Sensitivity and
misfit are **anti-correlated**, so resolving u > 0.9 would improve the predictions that do
not need improving and do essentially nothing for the ones that drive "the transfer does
not work."

**A weight table says where φ is SENSITIVE to a(u). It does not say where the model
FAILS**, and here the two sit at opposite ends of the CV axis. What survives is the
narrower statement: a knot at u\* cannot cover the families' weight, which is *why* U6's
verdict is what it is — not a prescription for where to measure next. (The same caution
the section's last paragraph already applies to the nested rule's convention had to be
extended to the plain-ADDtab argument.)

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
* **C-v** ⚠ **ADDED after the referee pass — it is not in the blind file.** Every number
  here rests on `beta`, and `mm.measured()` prints a standing **non-determinism** warning:
  15 of 82 duplicated keys disagree, worst |Δφ| = **5.03e−3 = 2.4 × PHI_SCALE**, resolved
  by **filename sort order**. It is tracked in `PREREG_DST_WINDOW_2026-09-09.md` and is not
  introduced by this block — but `beta`'s own precision carries that exposure, and neither
  the prereg nor the first claims draft said so.

## What is NOT established

* That the two-parameter (a, b) **form** is right at high u. Only its parameters were fitted.
* That a(u) between 0.70711 and 1 is anything in particular. It is still PCHIP between a
  measured knot and a proved limit, and `out/diag_ugap.txt`'s sweep shows that region moves
  the K = 5 prediction by up to 9.59 scale.
* Any statement about pressure or friction dependence (C-i).
* That saturation is **refuted**. C6 says undetermined; the registered gate says fail. Do
  not let the gate's wording become the claim.
