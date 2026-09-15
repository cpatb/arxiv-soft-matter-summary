# PREREG_KLADDER_2026-09-15 — does the non-additive residual ACCUMULATE with K, or is it a `y_S` effect already saturated at three classes?

**Written before any `kl*` cell existed.** Blind predictions are printed by
`predict_kladder.py` → `out/predict_kladder.txt`, which reads only the closed binary
(`mixr*`, `extr*`) and ternary (`tern*`) blocks. `grade_kladder.py` reproduces those
numbers off disk before grading anything against them.

Grid, models and admissibility all live in **`kladder.py`**, which is the single source
of truth: the prediction script, the runner and the grader import it and none of them
rebuilds the grid.

---

## 1. What is being tested, and why the last block could not test it

`grade_ternary.py` (closed, 2026-09-14) established that the pairwise-filling functional
is **not additive** beyond two classes. Six ternaries at `d` = 1, 2, 4 are under-predicted
**6/6**, mean **+7.53** scale, rms **7.83** against a bar of 2.8; confirmed at N = 3000
(mean **+8.67**, 6/6) with `n_L` and `d_max/L` excluded as confounds (pooled spearman
±0.3077, p 0.33).

Two things weaken it, and they are this block:

1. **`y_S` took only TWO values** — 0.10 and 0.15 — across the six admitted cells, with
   `y_L` ∈ [0.55, 0.70]. The residual tracks `y_S` at spearman **+0.878**, so the
   headline rests on a two-level contrast in a narrow corner.
2. **Nothing says whether the miss grows with the number of classes.** If it saturates at
   the K = 3 value the correction is a constant and can be absorbed. If it grows, every
   continuous-PSD prediction this project makes needs the new rule — and a continuous PSD
   is K → ∞. The held-out transfer that started this arm (tabulated closure χ² **108.1**
   against J's **34.9**, under-predicting the lognormal monotonically with CV to −6.9
   scale) is exactly what an accumulating miss would look like.

---

## 2. The design

A **geometric diameter ladder of FIXED SPAN 4**, filled in with more classes:

| K | diameters | pairs P | triples T | pair `u` values |
|---|---|---|---|---|
| 2 | 1, 4 | 1 | 0 | 0.2500 |
| 3 | 1, 2, 4 | 3 | 1 | 0.2500, 0.5000 |
| 4 | 1, 4^⅓, 4^⅔, 4 | 6 | 4 | 0.2500, 0.3969, 0.6300 |
| 5 | 1, 4^¼, …, 4 | 10 | 10 | 0.2500, 0.3536, 0.5000, 0.7071 |

**Every pair on every rung has `u` ≥ 1/4 by construction**, because the span is 4 and no
two classes are further apart than the extremes. That is the constraint that refused
`d` = 1, 2, 8 for the ternary block, and span 4 satisfies it identically at every K — so
additivity stays unconfounded with extrapolation on the whole ladder. `d_min`, `d_max`,
the size ratio and the span are the SAME on every rung: **the K axis adds classes without
widening the distribution.**

Composition is a one-parameter family `y_i ∝ d_i^q`, but the grid is addressed by the
**fines volume fraction** `y_S` = y(d_min): for each (K, `y_S`) the exponent `q` is
*solved* so the finest class carries exactly `y_S`. Levels
**`y_S` ∈ {0.03, 0.05, 0.08, 0.13, 0.18}** — five levels spanning **6.0×**, against the
closed block's two levels spanning 1.5×, and reaching **above** its maximum (0.18 > 0.15)
as well as 3.3× below its minimum.

### 2.1 The two axes are UNCROSSED, and that is the whole point

`predict_kladder.py` §1c, printed before the run:

* balanced sub-block (12 cells, K = 2,3,4,5 × `y_S` = 0.03, 0.05, 0.08):
  **corr(K, ln `y_S`) = +0.000000** — zero *by construction*, not by luck;
* whole core grid (16 cells; K = 2 and 3 additionally carry `y_S` = 0.13 and 0.18):
  corr = −0.3214.

This is the direct repair of the defect the ternary block had to break with a second N.
There, `n_L` and `y_S` were collinear *by construction* because n ∝ y/d³. Here K and the
composition axis are orthogonal by construction.

`n_L` is **not** held fixed and is not claimed to be — it cannot be, for the same n ∝ y/d³
reason. What breaks that confound is the **N = 3000 block, which runs the IDENTICAL cell
list** (§2.3) and moves every count by exactly 2× and `d_max/L` by 2^(−1/3).

### 2.2 The cells

Admissibility is applied BEFORE the run: `MIN_PHASE_COUNT` = 50 grains on **every** class
at **N = 1500**, inherited verbatim from `phase12_mixture.py`. Refused, not run and
dropped. Full table with `n_i`, `d_max/L`, CV, `d32`, `d43` and the family string per cell
is `out/predict_kladder.txt` §1b.

* **arm `T0`** — the monodisperse protocol control, written as the 3-class spec with two
  empty classes, i.e. the string `phase12_ternary.y_string` emitted for its own T0. At
  seeds 1–3 it is therefore the **bit-identical packing**.
* **arm `K`** — 16 cells: K = 2 and 3 at all five `y_S`; K = 4 and 5 at `y_S` = 0.03,
  0.05, 0.08. **REFUSED at N = 1500:** (K4, 0.13) at 47 grains, (K4, 0.18) at 31,
  (K5, 0.13) at 30, (K5, 0.18) at 19.
* **arm `SPLIT`** — 2 cells. The K = 3 cell at `y_S` = 0.05 and 0.08 with its **middle
  class split into two nearly-degenerate classes** at d_mid × (1 ∓ 0.05), each taking half
  its volume. K goes 3 → 4 and P goes 3 → 6 while the distribution barely moves: the
  additive prediction shifts **0.07 and 1.41 scale**, both inside the bar.
* **arm `NA2`** — the unexplained exception carried over from
  `PREREG_TERNARY_N_2026-09-14.md`: (`y_L`,`y_M`,`y_S`) = (0.55, 0.35, 0.10) at **six**
  seeds. Its family string is byte-for-byte the ternary block's, so seeds 1–3 are a
  reproduction gate and seeds 4–6 are the new information.
* **arm `KX`** — N = 3000 ONLY: the three (K, `y_S`) corners N = 1500 cannot afford
  ((K4, 0.13), (K4, 0.18), (K5, 0.13)). Labelled separately so **no N-comparison can
  silently include a cell that exists at one N only**. (K5, 0.18) is refused at both.

### 2.3 Both N, on the SAME cell list

N = 1500 (tag `kl`, 63 runs) and N = 3000 (tag `kl3k`, 72 runs = the same 63 plus the 9
`KX` runs). The core grid's admissibility is judged at N = 1500 *whatever N is asked
for*, so the two blocks contain the same cells and the N-doubling is a comparison of the
same points — lessons.md 2026-09-08, evaluating two probes of a sensitivity on different
samples makes the difference partly a change of sample.

`beta` is **re-measured at N = 3000** from that block's own monodisperse control, as
`PREREG_TERNARY_N_2026-09-14.md` did, so a finite-size shift in the reference is not
folded into the quantity under test.

### 2.4 Wiring gates — nothing reduced "by inspection"

`PREREG_TERNARY_2026-09-14` §3 asserted that its sequential rule "reduces exactly to the
binary closure at K = 2". It raised `IndexError` there. Every reduction this block relies
on is therefore **run**, in `predict_kladder.py` §0, before any prediction is quoted:

| gate | claim | result |
|---|---|---|
| W1 | `kladder.ideal_ceiling` at K = 2 == `mixture_model.furnas_ideal` | **12/12 to 1e−12** |
| W2 | `kladder.sequential_k` at K = 2 == the additive binary closure, both closures | **24/24 to 1e−12** |
| W3 | `kladder.sequential_k` at K = 3 == `grade_ternary.sequential`, both conventions, both closures | **24/24, worst 1.1e−16** |
| W4 | every family string round-trips through `parse_discrete_family` to the SAME floats | **23/23 exact** |
| W5 | the `mono` and `NA2_x6` strings are byte-identical to `phase12_ternary.y_string`'s | **True, True** |

W1 earned its keep immediately: the first `ideal_ceiling` took the a = b = 0 limit
*numerically* (pa = pb = 1e−13, giving a ≈ 3e−7) and missed `furnas_ideal` by up to 3e−8
on 6 of 12 binaries. Rewritten with the zeros substituted exactly. **A limit evaluated
numerically is not the limit** — lessons.md 2026-09-14, match the arithmetic, not the
algebra.

---

## 3. Bars, and where they come from

Re-derived and asserted in `predict_kladder.py` §2; **nothing in this block refits
anything.**

* TAB (the tabulated closure, primary) on all 41 binaries: rms **2.02**, max **5.74**.
* J (the shipped two-exponent closure, control) on its own 32 points: rms 2.83, max 6.69.

**BARS: rms 2.02, max 5.74 — TAB's own in-sample scatter.** A ternary or quaternary must
be predicted as well as the binaries the closure was calibrated on; anything else is a bar
invented for this block.

**Separation rule.** A comparison is declared **NOT INFORMATIVE** if its candidates'
predictions are separated by less than **1.6 ×** the bar over the cells it is made on.
Printed in §4 of the prediction file, before the data exists:

| comparison | mean \|sep\| (scale) | verdict |
|---|---|---|
| ADDtab vs SEQ`d32` | 5.28 | informative |
| ADDtab vs SEQ`d43` | 7.09 | informative |
| ADDtab vs SEQ`d10` | 2.05 | **NOT INFORMATIVE** (below 1.6 × 2.02 = 3.23) |
| ADDtab vs ADDj | 3.91 | informative |
| ADDtab vs IDEAL | 89.03 | informative |

---

## 4. The three hypotheses, and what this design can and CANNOT separate

Let `R(K)` be the mean residual (measured φ − additive TAB prediction) in PHI_SCALE units
at matched `y_S`. `R(2) = 0` is an **anchor, not a blind prediction**: K = 2 *is* TAB's
calibration set, so its residual is the in-sample noise floor by construction. Stated here
so no reader takes it for a success.

| hypothesis | meaning | R(3) | R(4) | R(5) |
|---|---|---|---|---|
| `H_sat` | composition only; K irrelevant | R₃ | R₃ | R₃ |
| `H_K` | one unit per extra **class** | R₃ | 2 R₃ | 3 R₃ |
| `H_P` | one unit per extra **pair**, P = K(K−1)/2 | R₃ | 2.5 R₃ | 4.5 R₃ |
| `H_T` | one unit per **triple**, T = C(K,3) | R₃ | 4 R₃ | 10 R₃ |

Pairwise separation ratios (printed in §5 of the prediction file):

| pair | at K = 4 | at K = 5 | verdict |
|---|---|---|---|
| `H_sat` vs `H_K` | 2.00 | 3.00 | separable |
| `H_sat` vs `H_P` | 2.50 | 4.50 | separable |
| `H_sat` vs `H_T` | 4.00 | 10.00 | separable |
| **`H_K` vs `H_P`** | **1.25** | **1.50** | **NOT SEPARABLE by this design** |
| `H_K` vs `H_T` | 2.00 | 3.33 | separable |
| `H_P` vs `H_T` | 1.60 | 2.22 | separable |

**REGISTERED IN ADVANCE: `H_K` and `H_P` cannot be told apart here.** P = K(K−1)/2 is a
deterministic function of K, so no design that moves only K can decouple them, and at
K ≤ 5 the two laws differ by 1.25–1.50×, below the 1.6 bar. **No verdict of this block may
claim to have distinguished "per class" from "per pair."** That is a limit of the design,
recorded before the data, not a discovery afterwards.

What the design *can* decide: **saturation vs growth**, **linear-in-K vs three-body**,
and — through the `y_S` axis at fixed K and the `SPLIT` control — **whether the residual
is a property of the size distribution at all.**

R₃ is recalibrated at grade time from **this block's own K = 3 cells at matched `y_S`**,
not from the closed block's 7.53 (which was measured at `y_S` = 0.10–0.15, off this
grid). The 7.53 version is printed for reference only.

---

## 5. Registered gates

Applied in code by `grade_kladder.py`. **A failure is the model's range, not a licence to
refit.** No bar below may be moved; if one turns out to be the wrong test, the replacement
must be *stricter in kind* and the original number stays printed.

### Admission (registered, applied in code)
A cell is admitted iff: no `error`; `converged` true; `n_requested` == the block's N;
`mu_p` == 0.30; and every design class count ≥ 50. `N` before and after the filter is
printed. Errors and non-converged cells are **counted**, never silently dropped.

---

**G0 — protocol control + cross-block reproduction.**
`mono` lands within **2.0 scale** of `beta`; AND at seeds 1–3 φ is **bit-identical** to
the ternary block's T0, same for `NA2_x6`.
*Prediction: PASS, with the bit-identity exact.* A non-zero difference is reported as a
finding (it would mean the WSL tree or the solver has moved), and the 2.0-scale bar is
then the binding one.

**G1 (PRIMARY) — the K profile.**
At each `y_S` ∈ {0.03, 0.05, 0.08} independently, R(K) for K = 2, 3, 4, 5 is compared with
the four hypothesis curves (R₃ taken from that `y_S`'s own K = 3 cell). Reported: the rms
distance of each hypothesis, the winner, and whether the winner separates from the runner
up by ≥ 1.6×.
*Sub-clause that can fail on its own:* **R(5) − R(3) ≥ 2.02 (the bar) on at least 2 of the
3 `y_S` levels.**
*Registered prediction: **the residual GROWS with K** — the sub-clause passes 3/3, and
`H_K` (or its unseparated twin `H_P`) beats both `H_sat` and `H_T`.*
*What failure means:* if the sub-clause fails on ≥ 2 of 3, the miss is **already saturated
at three classes**, the correction is a constant, and the continuous-family transfer
problem is NOT explained by accumulation.

**G2 — the `y_S` axis, at fixed K.**
Over the five `y_S` levels at K = 3 (and, reported beside it, at K = 2): spearman(R, `y_S`)
**≥ +0.70**, AND R(`y_S` = 0.18) − R(`y_S` = 0.03) **≥ 2.02**.
*Registered prediction: PASS at K = 3.* The closed block's +0.878 came from two levels; a
six-fold sweep either reproduces it or reveals it as a fluke of that corner.
*What failure means:* the composition reading of T1 fails and the residual is driven by
something else on the grid.

**G3 — the `SPLIT` control (the sharpest single test here).**
R(`SPLIT_yS`) − R(`K3_yS`) at `y_S` = 0.05 and 0.08. The split moves K 3 → 4 and P 3 → 6
while leaving the distribution essentially unchanged (the additive prediction moves 0.07
and 1.41 scale).
*Registered prediction: **the split does NOT move the residual** — |Δ| < 2.02 on 2 of 2.*
*What failure means:* if |Δ| is instead comparable to the genuine K 3 → 4 increment from
G1, the residual is **counting classes**, which would make it an artefact of how a PSD is
discretised rather than physics.

> **G1 and G3 are NOT in tension and can both pass.** G1 manipulates K with genuinely
> distinct sizes; G3 manipulates K with a near-duplicate. Under a literal class-counting
> law G1 passes and G3 fails; under saturation G1 fails and G3 passes; under "genuine size
> diversity" both pass. All four combinations are reachable. This paragraph exists because
> the datacenter block once registered two gates that could not both pass
> (lessons.md 2026-09-08) and the tension, not either verdict, was the error.

**G4 — N-doubling, on the identical cell list.**
Every count ×2, `d_max/L` ×2^(−1/3), no composition variable moved. The quantity graded is
the **K-profile increment** R(5) − R(3) at each `y_S`: |Δ_3000 − Δ_1500| must be within
**2 × the pooled seed standard error** of the increment. Also reported: pooled
spearman(R, `n_min`) and spearman(R, `d_max/L`) over both N, against the closed block's
registered 0.6 bar.
*Registered prediction: PASS — the increment survives the doubling.*
*What failure means:* the K growth is finite-size, exactly as the closed block's residual
might have been.

**G5 — the NA2 exception.**
(0.55, 0.35, 0.10) at six seeds and both N. Registered quantity: φ(N = 3000) − φ(N = 1500)
in scale units, with its **standard error of the DIFFERENCE** (not the per-seed scatter —
lessons.md 2026-09-02), and a Welch t against n = 6 per arm.
*Registered prediction: with six seeds the +1.81 shift falls below 2 standard errors of
the difference, i.e. it was seed noise.* A shift that survives at n = 6 is a real finding
about that composition and is reported as such.

**G6 — the SEQUENTIAL rule, its convention SWEPT not chosen.**
rms over the whole admitted grid for ADDtab and for SEQ at all four composite-diameter
conventions (`d10`, `d21`, `d32`, `d43`). Reported with the sign structure (how many cells
each sits above/below the measurement), because the closed block's T4 showed the rms
winner carries the OPPOSITE bias of comparable size to the loser — "sequential beats
additive" must not be read as "the rule is sequential".
*ADDtab vs SEQ`d10` is pre-declared NOT INFORMATIVE (mean separation 2.05 < 3.23).*
*Registered prediction: the measurement is BRACKETED — ADD below on ≥ 80 % of cells and
SEQ`d43` above on ≥ 80 % — as it was at K = 3.*

**G7 — a directional prediction the additive closure can fail outright.**
At matched `y_S`, **the additive closure predicts φ FALLS with K** (e.g. at `y_S` = 0.08:
0.6517 → 0.6380 → 0.6402 for K = 3, 4, 5). Physically, more intermediate sizes inside a
fixed span should fill better.
*Registered prediction: the MEASURED φ RISES with K at matched `y_S` on ≥ 2 of the 3
levels, contradicting the closure's own sign.* This is a sign test that needs no residual
bookkeeping and no bar beyond the seed error.

**G8 — transfer to the continuous families (declared POST-HOC, not blind).**
Both datasets are already on disk, so this is **not** a pre-registered gate and will not be
reported as one. With the sequential rule's convention swept rather than chosen, the
held-out χ² of the lognormal/uniform/bidisperse families is recomputed for ADD and for all
four SEQ conventions against the standing ADD 108.1 / J 34.9. Reported as a description of
whether non-additivity can account for that gap, with the convention's own spread printed
beside it.

---

## 6. What would make this block report NOTHING

Stated in advance so the answer is not chosen afterwards:

* if G0's bit-identity fails, every cross-block statement is withdrawn and the block is
  reported as measuring a different protocol;
* if the K = 2 cells do not sit at TAB's in-sample noise floor, the K axis has no zero and
  G1 is reported as uncalibrated;
* if the winner of G1 does not separate from the runner-up by 1.6×, G1 is reported **NOT
  INFORMATIVE** regardless of which curve is nearest;
* `H_K` vs `H_P` is not reportable at all (§4).
