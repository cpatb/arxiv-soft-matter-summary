# CLAIMS_KLADDER_2026-09-15 — what the K-ladder block establishes, and what it does not

Written FROM `out/grade_kladder.txt`, `out/transfer_kladder.txt`,
`out/diag_kladder_descriptor.txt` and `out/predict_kladder.txt`. Pre-registration:
`PREREG_KLADDER_2026-09-15.md` (committed 7b8fcf0, before any cell ran).

Campaign: **135 runs, 135 admitted, 0 errors, 0 non-converged** — 63 at N = 1500 and 72
at N = 3000 on the identical cell list plus three N = 3000-only corners.

---

## C1 (HEADLINE) — the non-additive miss does NOT accumulate with K: it STEPS ONCE
between three and four classes. Whether it then SATURATES is **undetermined** (C1a ii).

At matched fines fraction `y_S`, with `R(K)` the residual (measured φ − additive
tabulated closure) in PHI_SCALE = 2.1e−3 units:

| `y_S` | R(2) | R(3) | R(4) | R(5) | **R4−R3** | **R5−R4** | ratio |
|---|---|---|---|---|---|---|---|
| 0.03 | −0.85 | +1.68 | +5.09 | +5.29 | **+3.42** (se 1.26) | **+0.20** (se 0.77) | 0.06 |
| 0.05 | +0.18 | −0.39 | +7.34 | +5.76 | **+7.73** (se 0.76) | **−1.58** (se 1.06) | −0.20 |
| 0.08 | −0.40 | +1.24 | +9.27 | +7.50 | **+8.03** (se 1.12) | **−1.77** (se 0.80) | −0.22 |

and at N = 3000 on the same cells: +4.54/+0.63, +7.44/−0.80, +7.01/−1.21 (ratios 0.14,
−0.11, −0.17).

**The second increment is zero or negative on 6 of 6 (K, N) combinations.** The
registered hypotheses predict increment ratios of **1.00** (`H_K`, per class), **1.33**
(`H_P`, per pair) and **2.00** (`H_T`, per triple); the measured ratios run **−0.22 to
+0.14**. All three growth laws are refuted, and `H_sat` (no K dependence at all) is
refuted by the first increment, which is +3.4 to +8.0 against a bar of 2.02.

**None of the four registered hypotheses describes the data** -- with the important qualification of C1a below, which was found while writing this file and changes how the step must be read. The registered
sub-clause — R(5) − R(3) ≥ 2.02 on ≥ 2 of 3 levels — passed **3 of 3**, so "the residual
grows with K" is correct as far as it goes; the finer structure is that all of the growth
is in the single step from three classes to four.

### C1a — ⚠ MOST OF THE STEP IS THE CLOSURE FALLING, NOT THE PACKING RISING, AND C1 MUST BE READ THROUGH THAT

`R = (φ − ADD)/scale`, so any step in R has two parts. Decomposed
(`ΔR = Δφ − ΔADD`, scale units):

| N | `y_S` | step | ΔR | **Δφ (measured)** | **ΔADD (model)** | model's share |
|---|---|---|---|---|---|---|
| 1500 | 0.03 | K3→K4 | +3.42 | **+2.51** | −0.91 | 27 % |
| 1500 | 0.05 | K3→K4 | +7.73 | +1.63 | **−6.09** | **79 %** |
| 1500 | 0.08 | K3→K4 | +8.03 | +1.52 | **−6.51** | **81 %** |
| 3000 | 0.03 | K3→K4 | +4.54 | **+3.61** | −0.93 | 20 % |
| 3000 | 0.05 | K3→K4 | +7.44 | +1.32 | **−6.12** | **82 %** |
| 3000 | 0.08 | K3→K4 | +7.01 | +0.52 | **−6.49** | **93 %** |

**At `y_S` = 0.05 and 0.08, 79–93 % of the "step" is the additive closure's own
prediction DROPPING by 6.1–6.5 scale when the fourth class is added** — the measurement
moves only +0.5 to +1.6. Only at `y_S` = 0.03 is the step mostly a genuine rise in φ
(+2.51 and +3.61, against seed se of 1.2 and 0.4). The same applies to the saturation —
but **not uniformly, and the three `y_S` levels are not interchangeable**: at K4→K5 the
model comes back UP by +2.17 (`y_S` 0.05) and +1.05 (0.08), which is most of why ΔR goes
negative there, while at **`y_S` = 0.03 it does not come back up at all** (ΔADD = −0.05 at
both N). At 0.03 the flat ΔR (+0.20, ratio 0.06) is therefore mostly a genuine
deceleration in φ's own rise (Δφ = +0.15 of ΔR = +0.20); at 0.05 and 0.08 it is almost
entirely the closure bouncing. An earlier draft wrote "the model comes back UP by +1.0 to
+2.2" without that exception.

So C1 is a correct statement **about the residual** and it does refute the three growth
laws, which are hypotheses about the residual. It is **not** a statement that the packing
does something dramatic between three and four classes: measured φ rises smoothly and
modestly across the whole K axis. **The robust, decomposition-free version of the finding
is C4** — the closure's sign in K is wrong — and C1 should always be quoted next to it.

Two further qualifications, both from the physics-referee pass of 2026-09-15, both now
printed by the grader as **G1c** rather than argued here:

**(i) the yardstick moves by BRANCH REASSIGNMENT, not only by interaction values.**
`packing_density` is a MINIMUM over which class is dominant, so its value can change
because a different branch wins. The argmin index (0 = coarsest) runs:

| `y_S` | K = 2 | K = 3 | K = 4 | K = 5 |
|---|---|---|---|---|
| 0.03 | 0 | 0 | **1** | 1 |
| 0.05 | 0 | 0 | **1** | 1 |
| 0.08 | 0 | **1** | 1 | **2** |

**The branch changes somewhere on the K axis on 3 of 3 rows.** The referee's stronger
reading — that the K3→K4 drop *coincides* with the reassignment — holds at `y_S` = 0.05
but **not** at 0.08, where the −6.51 drop happens entirely *within* branch i = 1; and at
`y_S` = 0.03 the branch does switch while the drop is only −0.91. So reassignment is
present everywhere and is **not** in one-to-one correspondence with the size of the drop.
The supportable statement is the weaker one: **a step in ADDtab across K is not by itself
a filling statement, because part of it is the min-over-branches structure reacting to
how volume is redistributed when a class is added.**

**(ii) K = 5 is evaluated where the closure has NO DATA, so the SATURATION IS NOT
ATTRIBUTABLE.** The tabulated interactions are a PCHIP through the identifiable binary
ratios at `u` = 0.250, 0.333, 0.500, 0.667 plus the *proved but unmeasured* limits at
`u` = 0 and 1. Pair coverage:

| K | pairs | `u` values (count) | above the top measured knot 0.667 |
|---|---|---|---|
| 2 | 1 | 0.2500 ×1 | none |
| 3 | 3 | 0.2500 ×1, 0.5000 ×2 | none |
| 4 | 6 | 0.2500 ×1, 0.3968 ×2, 0.6300 ×3 | none |
| 5 | 10 | 0.2500 ×1, 0.3536 ×2, 0.5000 ×3, **0.7071 ×4** | **4 of 10** |

**K = 5's most numerous pair type sits above every measured ratio**, in the widest gap in
the whole construction. The rung that decides "saturation" is therefore evaluated on
interpolation, and at `y_S` = 0.05 the K4→K5 "recovery" in ADDtab (+2.17, no branch
change) is a statement about how PCHIP guesses `a(0.707)`, not a measurement. **C1's
"and then saturates" must be read as UNDETERMINED by this design.**

What survives
untouched is the K3→K4 rise (every pair on both those rungs is inside the measured
bracket) and the sign test C4, which compares data with model and does not care why the
model does what it does. This is the symmetric twin of the caveat already in "What this
block does NOT establish" about `u` < 1/4 — the block says nothing trustworthy about
`u` ∈ (0.667, 1) either, and it was missing until the referee asked.

> **ADDENDUM 2026-09-16 — the gap has now been MEASURED, and this caveat was right.**
> `CLAIMS_UGAP_2026-09-16.md` puts a knot at `a(0.70711) = 0.9481 [0.8799, 0.9699]`.
> Re-running the comparison with it, the registered saturation gate **FAILS** (R5−R4 ≤ 0
> on 1 of 3 levels at N = 1500 and 0 of 3 at N = 3000) — but swept over that knot's own
> profile interval, saturation survives on **33 of 87 jointly admissible corners
> (37.9 %)**. So "and then saturates" remains **UNDETERMINED**, now for a measured reason
> rather than a structural one, and it is **not** refuted either — though at 37.9 % the
> balance leans *against* it. The K3→K4 rise and the sign test C4 are unaffected,
> **exactly as this paragraph predicted.**

### C1b — the registered *ranking* half of G1 is NOT INFORMATIVE, and the grader says so
Every hypothesis in the registered family is a MULTIPLE of R(3), and R(3) landed **inside
the 2.02 bar of zero on 5 of the 6** (K, N) rows (|R3| = 0.39 to 1.81). A multiple of a
quantity consistent with zero predicts nothing. The ranking is printed with
`<-- R3 INSIDE THE BAR: ranking DEGENERATE` beside it and is listed under NOT INFORMATIVE.
The absolute-increment table above is the statistic that decides, it is reported BESIDE
the registered one, and **no bar was moved.**

---

## C2 — it is NOT class counting: a degenerate fourth class does almost nothing, where a genuine one does everything (G3)

Splitting the middle class of a K = 3 pack into two classes at d_mid × 0.95 and × 1.05,
each taking half its volume — K goes 3 → 4, pairs 3 → 6, the distribution barely moves:

| | `y_S` 0.05 | `y_S` 0.08 |
|---|---|---|
| R(split) − R(K = 3), N = 1500 | **+0.83** (se 1.37) | **+1.45** (se 0.84) |
| R(split) − R(K = 3), N = 3000 | **−0.31** (se 0.98) | **+2.22** (se 0.80) |
| the genuine K 3 → 4 increment at the same `y_S`, N = 3000 | +7.44 | +7.01 |

🔴 **The registered sentence — "|Δ| < 2.02 on 2 of 2" — does not name an N, and the code
scoped it to N = 1500, which is a NARROWER question than the one registered.** Widened
after the referee pass; `grade_kladder.py` now grades **both** N, and **G3 therefore FAILS
at N = 3000, `y_S` = 0.08 (+2.22 against the 2.02 bar)**. The bar is unchanged and the
failure is in `FAILURES`. (lessons.md 2026-09-13: a pre-registered rule is a SENTENCE, and
the code can compute a different quantity than the sentence names.)

**The decisive comparison is the RATIO, and it must be quoted with all four values, not
three.** The near-duplicate class buys, as a fraction of what a genuine one buys:

| | `y_S` 0.05 | `y_S` 0.08 |
|---|---|---|
| N = 1500 | +11 % | +18 % |
| N = 3000 | **−4 %** | +32 % |

An earlier draft quoted "11–32 %", which silently dropped the one **negative** value
(−0.31 / +7.44 at N = 3000, `y_S` = 0.05 — the split moved the residual the *wrong* way,
though at 0.3 se it is indistinguishable from zero). **Honest restatement: on 3 of 4
measurements a near-duplicate fourth class buys 11–32 % of what a genuinely intermediate
one buys, and on the fourth it buys nothing.** The conclusion — the residual is a
functional of the size distribution, not a count of classes — is supported by all four,
since none of them approaches the +7.0 to +8.0 a genuine class delivers; but the range as
printed was not the whole set.

---

## C3 — the `y_S` dependence survives a six-fold widening, and it is a K ≥ 3 phenomenon (G2)

Five levels of `y_S` spanning 6.0× (0.03 → 0.18), against the closed ternary block's two
levels spanning 1.5×:

| K | spearman(R, `y_S`), N = 1500 | N = 3000 | R(0.18) − R(0.03) |
|---|---|---|---|
| 3 | **+0.700** (p 0.188) | **+0.900** (p 0.037) | +8.49 / +10.95 |
| 2 | −0.600 (p 0.285) | −0.500 (p 0.391) | −1.21 / −0.63 |

Registered bars: spearman ≥ +0.70 and span ≥ 2.02. **PASS at K = 3 at both N** — but the
N = 1500 correlation lands **exactly on** the bar with p = 0.188, so the claim rests on
the N = 3000 replicate, and that has to be said.

⚠ **The K = 2 row is NOT a clean control, and "control" overstates it** (referee pass,
2026-09-15). Three things weaken it: its spearman is **not significant at either N**
(p = 0.285, 0.391); `K2_yS0.13` gives R = **−2.97**, which *exceeds* TAB's own in-sample
rms of 2.02 and so breaks the R(2) ≈ 0 anchor at that composition (it is still inside the
in-sample **max** of 5.74); and `K2_yS0.13` and `K2_yS0.18` carry the **two largest
N-doubling swings in the whole 18-cell campaign** (+3.84 and +2.12 against a campaign rms
of 2.04). **Honest restatement: the K = 2 trend is qualitatively opposite in direction to
K = 3's and is not resolved.** That is still useful — a `y_S` dependence created by the
binary calibration leaking in would have to appear at K = 2 with the *same* sign, and it
does not — but it is a direction, not a measured null.

---

## C4 — the additive closure has the WRONG SIGN in K (G7) — **2 of 3 levels resolved, 1 undetermined**

At matched `y_S`, the additive tabulated closure predicts φ **falls** as classes are added
inside a fixed span. Measured, it **rises** — on 3 of 3 levels by point estimate, of which
**2 of 3 clear their own standard error**:

| `y_S` | measured φ, K = 3, 4, 5 | K5−K3 (scale) | ADDtab says | ADDtab K5−K3 |
|---|---|---|---|---|
| 0.03 | 0.6325, 0.6377, 0.6381 | **+2.65** (se 1.22) | 0.6290, 0.6271, 0.6269 | −0.97 |
| 0.05 | 0.6432, 0.6467, 0.6479 | **+2.22** (se 0.99) | 0.6441, 0.6313, 0.6358 | −3.93 |
| 0.08 | 0.6543, 0.6575, 0.6560 | **+0.80** (se 0.85) | 0.6517, 0.6380, 0.6402 | −5.45 |

This needs no residual bookkeeping and no closure-specific bar: it is a sign test, and the
model is on the wrong side of it everywhere. 🔴 **But "3 of 3" was the wrong headline**
(referee pass): at `y_S` = 0.08 the measured rise is **+0.80 against se 0.85**, a point
estimate smaller than its own error and consistent with zero. **Report it as "2 of 3
resolved, 1 undetermined."** The two resolved levels (+2.65 / se 1.22 and +2.22 / se 0.99)
are what carry C4, and the model's own prediction is negative and *large* on all three
(−0.97, −3.93, −5.45), so the model-vs-data disagreement is resolved at 3 of 3 even where
the measured rise on its own is not. `test_kladder.py` asserts `resolved == 2`, so the
honest count is pinned and cannot drift back to 3.

---

## C5 — the finite-size nulls are excluded, and the residual GROWS with N rather than shrinking (G4)

Every cell was run at N = 1500 and N = 3000 with composition, span, CV and every pair
ratio identical, moving every class count ×2 and `d_max/L` ×2^(−1/3).

* the graded quantity, R(5) − R(3), moves by 1.55, 0.49 and 0.46 against 2×se of 2.74,
  3.10 and 2.08 — **PASS 3 of 3**;
* 🔴 **the `n_min` correlation is a CONFOUNDED DIAGNOSTIC and its 0.6 bar is WITHDRAWN as
  meaningless** (referee pass, 2026-09-15). `n_min` is collinear with **both** design axes
  by construction: spearman(`y_S`, `n_min`) = **−1.0000 on all 8 (K, N) combinations**,
  and `n_min` also falls monotonically with K at fixed `y_S` (because n ∝ y/d³ and adding
  classes dilutes every count). A pooled spearman(R, `n_min`) is therefore **guaranteed**
  to come out large and negative whether or not any finite-size leak exists — it is a
  restatement of the grid, not a test of it. The measured **−0.4694 (p 0.0067)** is
  reported for the record and **must not be read as "a partial residual finite-size
  leak."** An earlier draft did read it that way. `d_max/L` is not collinear in the same
  way and sits at +0.1400 (p 0.445), which is a real if weak reassurance. **The valid
  finite-size control is the N-doubling comparison above**, which is the gated one and
  passes 3 of 3;
* the grader previously printed "bar 0.6" beside a number that **never fed `failures`** —
  a bar that cannot fail. It now prints the collinearity that makes the number
  uninformative instead. Nothing that *was* gated changed;
* the residual **rises** with N on 18 of 18 cells (mean **+1.87**, rms 2.04 scale). A
  finite-size artefact must shrink as the box grows.

`d_max/L` never exceeded 0.237 (N = 1500) or 0.196 (N = 3000), below the 0.25 flag.

---

## C6 — the NA2 "exception" was never an exception (G5)

The closed ternary block left (0.55, 0.35, 0.10) as its one unexplained cell: +1.81 scale
between N = 1500 and 3000 on three seeds against an se of 0.48, where the other five
agreed. Re-run at **six** seeds at both N:

* the shift is **REAL**: **Δφ/scale = +2.19**, se of the difference 0.29, **7.64 se**,
  Welch t = 7.64, p = 0.0001, df 6.8. **The registered prediction (that six seeds would
  show it as seed noise) FAILS.**
* and it is **NOT AN OUTLIER OF THIS GRID**: its **ΔR = +1.89** against a mean **ΔR** of
  **+1.87** over the 18 cells (rms 2.04, range +0.14 to +3.84).
  ⚠ **The two numbers are different quantities and must be labelled** (referee pass):
  +2.19 is Δφ/scale, the quantity actually significance-tested; +1.89 is ΔR, which is net
  of the model's own `beta` recalibration between the two N and is the right thing to
  compare against a pooled ΔR mean. Both are defensible; quoting them adjacently without
  labels was not.
  ⚠ And the comparator is **the entire enumerated design grid, not a sample from a
  population**, spanning +0.14 to +3.84 — a 27× range. "Near the mean" therefore
  establishes only that NA2 is **not an extreme outlier of this grid**, which is the
  conclusion drawn; it is not a distributional typicality claim and "mean ± rms" language
  should not be read as one.

So the open item closes the other way round from how it was posed: there is nothing
special about that composition, and the quantity that needed the second N was the
block-wide offset, not one cell.

Seeds 1–3 reproduced the closed block **bit-identically** (6 of 6 including the
monodisperse control, `diff = +0.000e+00`), so this is the same packing, not a near one.

---

## C7 — the nested rule beats the additive one by 1.93×, and its own undeclared convention moves it by 1.96× (G6, FAILED sub-clause)

rms over the 19 admitted multi-class cells at N = 1500:

| model | rms | max\|res\| | mean signed | above | below |
|---|---|---|---|---|---|
| SEQ`d21` | **2.70** | 6.14 | +1.00 | 13 | 6 |
| SEQ`d32` | 3.05 | 8.14 | −1.19 | 6 | 13 |
| ADD`j` | 4.43 | 10.17 | +3.26 | 16 | 3 |
| SEQ`d10` | 4.46 | 8.82 | +3.08 | 15 | 4 |
| ADD`tab` | 5.22 | 10.17 | +3.39 | 14 | 5 |
| SEQ`d43` | 5.28 | 13.82 | −3.30 | 5 | 14 |
| IDEAL | 101.79 | 162.81 | −86.49 | 0 | 19 |

ADD`tab` / SEQ`d21` = **1.93** (bar 1.6) — but the **spread across the four composite
diameter conventions is 1.96**, i.e. **the undeclared convention is worth as much as the
result.** The convention was swept, not chosen, which is why this is visible.

**The registered BRACKET sub-clause FAILED:** ADD`tab` sits below the measurement on 74 %
of cells and SEQ`d43` above on 74 %, against a registered ≥ 80 %. The bar stays and the
failure is reported. Both additive closures (TAB and the shipped J) under-fill on the
majority of cells, which is the ternary block's finding reproduced on a disjoint grid.

---

## C8 (POST-HOC, DECLARED) — swapping the sum for a nested rule does NOT settle the continuous-family transfer

`out/transfer_kladder.txt`. The standing numbers reproduce exactly: held-out χ² over the
15 family points is **34.9** for the shipped J and **108.1** for the additive tabulated
closure.

| rule | χ² |
|---|---|
| ADD`j` | 34.9 |
| SEQ`d21` | **42.8** |
| SEQ`d32` | 48.0 |
| ADD`tab` | 108.1 |
| SEQ`d10` | 186.3 |
| SEQ`d43` | 509.3 |

The best convention closes **+0.89** of the ADD`tab` − J gap; the worst **opens it by
5.48×**. The two straddle zero, the convention spread is **11.9×**, and the registered
reading applies: **this comparison does not settle whether non-additivity explains the
transfer gap.** What does survive is the direction — ADD`tab` under-predicts the lognormal
monotonically with CV (+0.13 → −6.85 scale) and SEQ`d43` over-predicts it monotonically
(+1.49 → +17.91), so the measurement is bracketed there too.

Two limits are printed rather than argued: the families' pair weight below `u` = 0.25,
where the interactions are interpolated rather than measured, is ≤ 2.8 % (lognormal
CV 0.5) and zero on 11 of 15 points; and the bidisperse family arrives at K = 2, where all
six columns are identical — the K = 2 wiring gate appearing inside the transfer itself.

---

## C9 (POST-HOC, DECLARED, AND IT IS A SEARCH NOT A TEST) — what orders the residual

`out/diag_kladder_descriptor.txt`, 18 cells at each N. Ranked by |spearman| against R,
with the sample-size nulls in the same table:

| descriptor | N = 1500 | N = 3000 |
|---|---|---|
| CV | **−0.882** | **−0.887** |
| `fill_weight` = Σ y_i y_j (1 − u_ij) | **+0.878** | **+0.889** |
| `pair_weight` = Σ y_i y_j | **+0.874** | **+0.878** |
| `triple_weight` = Σ y_i y_j y_k | **+0.873** | **+0.869** |
| `n_min` (NULL, sample size) | −0.702 | −0.726 |
| K = n_pairs = n_triples (rank-identical) | +0.680 | +0.659 |
| `d_max/L` (NULL) | +0.371 | +0.404 |
| `y_S` alone | +0.173 | +0.207 |

**Caveats that matter more than the ranking.** The top four are not four candidates: for a
volume vector `y`, Σ_{i<j} y_i y_j = (1 − Σ y_i²)/2, so `pair_weight` is a measure of how
EVENLY the volume is spread across classes, `triple_weight` and `fill_weight` are close
relatives of it, and CV falls as that evenness rises — **they are one quantity seen four
ways, and the block cannot separate them.** K, `n_pairs` and `n_triples` are
rank-identical by construction and score **below the sample-size null**, which is the
same verdict C1 reaches by a different route. And `y_S` scores near zero **because the
design deliberately uncrosses it from K**, so a pooled correlation cannot see it — C3 is
where the `y_S` effect is measured, at fixed K.

A descriptor that wins on the cells that motivated the search has been **selected, not
tested**. The next block is the one that places cells where these candidates disagree.

---

## What this block does NOT establish

* **"Per class" vs "per pair" (`H_K` vs `H_P`) was never decidable here** and was recorded
  as such before the data existed: P = K(K−1)/2 is a deterministic function of K, and at
  K ≤ 5 the two laws differ by 1.25–1.50×, below the 1.6 bar. In the event both are
  refuted, so the limit did not bind — but it was a real limit and it is still registered.
* Nothing about `u` < 1/4. Span 4 puts every pair on every rung at `u` ≥ 1/4 by
  construction, which is what keeps additivity unconfounded with extrapolation — and it
  means the block says nothing about the regime the shipped closure was already known to
  fail in (X1, `u` ≤ 1/3).
* 🔴 **And nothing about `u` ∈ (0.667, 1) either — the symmetric caveat, missing until
  the referee asked.** The tabulated interactions have no measured knot above `u` = 0.667,
  and **4 of K = 5's 10 pairs sit at `u` = 0.7071**. So the K = 4 → 5 rung, the one that
  would decide saturation, is evaluated on interpolation toward a proved-but-unmeasured
  endpoint. The K = 2 → 3 → 4 rungs are entirely inside the measured bracket and are not
  affected.
* **That any step in the ADDITIVE prediction across K is a filling statement.** It is a
  minimum over branches, and the argmin reassigns somewhere on the K axis at 3 of 3 `y_S`
  levels (C1a i). The residual is still the registered quantity and the growth laws are
  still refuted by it, but the *mechanism* behind the yardstick's own motion is not
  established here.
* **A resolved sign test at `y_S` = 0.08** (C4): +0.80 against se 0.85.
* Nothing about friction: every cell is `mu_p` = 0.30.
* Nothing about a continuous PSD directly. C8 is the transfer, and C8 does not settle.
* The saturation in C1 is measured at K = 4 → 5 only, **and is undetermined even there**
  (C1a ii). 🔴 **Extending the same ladder to K = 8 or 20 — which is where a continuous
  PSD lives — would make this WORSE, not better**, and an earlier draft of this bullet
  said the opposite. On a span-4 geometric ladder the adjacent-pair ratio is
  `4^(−1/(K−1))`: 0.707 at K = 5, **0.82 at K = 8, 0.88 at K = 12**, so every class added
  piles more pair weight into the stretch where `a(u)` and `b(u)` have no measured knot.
  **The blocker is therefore not K — it is binary calibration above `u` = 0.667**, which
  costs what a ternary block costs and unlocks the whole axis. Until it exists, no ladder
  can answer "does the miss saturate?"
