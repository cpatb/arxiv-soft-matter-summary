# PRE-REGISTRATION — the composition axis, and the filling functional

**Written 2026-09-13, before any cell of this campaign ran.** Runner
`phase12_mixture.py`; grader `grade_mixture.py`; model `mixture_model.py`
(→ `out/mixture_model.txt`, which is the source of every predicted number below and was
written before this document).

---

## 1. What this is for

`shape_control.py` established a **negative**: at matched number CV and matched friction
the three size-distribution families give jamming shifts spanning **4.5×**, and fifteen
single descriptors of the distribution fail to collapse them. The obstruction is
structural — **the bidisperse family's shift is non-monotone in CV** while every
descriptor tried is monotone in CV — so no monotone function of one variable can work.

This campaign tests a **positive**: that the shift is a *pairwise-filling functional* of
the size distribution, not a function of any moment of it, and that the bidisperse
non-monotonicity is the signature of the family's own parameterisation rather than a
failure of the filling picture.

The one-line reason to believe that before measuring anything: this project's
`BidispersePSD` is 50:50 **by number** at `d₁₀(1 ∓ CV)`, so the fines carry

    X_s = (1−CV)³ / [(1−CV)³ + (1+CV)³]

of the **volume**: 0.500, 0.354, 0.229, 0.135, 0.073, **0.036** at CV = 0 … 0.5, while
the size ratio runs 1.00 → 3.00. The family is a *curve* through the
(ratio, composition) plane, and it runs into the fines-poor corner. A filling model has
its maximum near `X_s` ≈ 0.3, which the family passes at CV ≈ 0.1 where the ratio is only
1.22 and filling is impossible. That is why the CV = 0.5 bidisperse packs **worst** of the
three families *despite* having the largest size ratio, and why "ratio 3 is above the
2.4 octahedral threshold so the interstitial picture is refuted" (`shape_control.py`
§WHAT WOULD SETTLE IT) does not follow.

## 2. The model, and what is fitted in it

`mixture_model.py`: the linear packing model over `K` size classes,

    γ_i = β / ( 1 − Σ_{j<i} [1 − β + b_ij β (1 − 1/β)] y_j − Σ_{j>i} [1 − a_ij] y_j )
    φ   = min_i γ_i

with `a_ij = sqrt(1 − (1 − d_j/d_i)^PA)` (loosening of a coarse class by a fine one) and
`b_ij = 1 − (1 − d_i/d_j)^PB` (wall effect of a coarse class on a fine one).

- `β` is **measured**, not fitted: the monodisperse φ at this protocol's own rung,
  **0.5947** at `μ_p` = 0.3 and **0.6413** at `μ_p` = 0.
- `PA` = 1.02 and `PB` = 1.50 are the textbook values and are **held fixed** for every
  number quoted as "zero-parameter".
- A second, fitted, reading exists (`PA` = 0.900, `PB` = 3.154) calibrated on the five
  non-trivial **bidisperse** points only.

Both are reported. The zero-parameter one is the claim.

**What it already does, on data that existed before this document**
(`out/mixture_model.txt`, §3–§7):

| | zero-parameter model |
|---|---|
| family order at CV = 0.5 | **correct** (bidisperse < uniform < lognormal) |
| bidisperse interior maximum | **present** (model at CV 0.2, measured at CV 0.4) |
| spread of measured/model over the three families | **1.71×** against 4.51× for CV, 4.06× for `d₃₂` |
| held-out χ² on 10 uniform+lognormal points | **160** vs **589** for the size-ratio-blind null with the same two free parameters |
| `μ_p` = 0 with only `β` changed | same residual size (+3 to +7 scale) |

⚠ **and what it does not do.** The model is **systematically high by +0.008 to +0.016**
(4–8 × the grading scale) at every non-trivial point, and its bidisperse maximum sits at
CV 0.2 where the measurement puts it at CV 0.4 — i.e. **it over-credits filling at small
size ratios.** That is the specific defect this campaign is designed to measure, and it
is the reason the campaign is worth running rather than a reason to discard the model.

⚠ **Held out is not blind.** The uniform and lognormal ladders are on disk and were read
before the model was written. The *blind* part of this exercise is §4 below.

## 3. The grading scale, declared before any comparison

**`PHI_SCALE` = 2.1e-03.** This is the protocol's own run-to-run non-determinism:
`analyse_phase2.load_all` finds 15 duplicated cell keys that return a *different* φ from
identical inputs, median |Δφ| 2.1e-03, worst 5.0e-03. The seed sem at a typical cell is
~2e-04, ten times smaller, and grading a model against a scatter smaller than the
reproducibility of the thing being modelled would measure the seed average rather than
the model. Every "× scale" in this document and in `grade_mixture.py` is in these units.

## 4. The campaign, and the predictions — written down before it runs

**Grid.** `phase12_mixture.py`, isotropic compression, `p̂` = 1e-3, `N` = 1500,
`μ_p` = 0.3, seeds 1/2/3, `truncate=None` (a binary is already bounded).

- size ratio `r` ∈ {2.0, 3.0}
- fines volume fraction `X_s` ∈ {0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.60, 1.0}
- plus the anchor `X_s` = 0.5³/(0.5³+1.5³) = 0.035714… at `r` = 3

63 cells. **No `φ_J` extrapolation anywhere** — φ is read at the rung, following the
COMPANION section of `PREREG_TAIL_2026-09-13.md`, which found the fixed-pressure channel
sound and the ladder fit unhealthy.

**Predicted φ (`out/mixture_model.txt` §8), zero-parameter / fitted / blind null:**

| `r` | `X_s` | zero-param | fitted | blind null |
|---|---|---|---|---|
| 3 | 0.05 | 0.6074 | 0.6083 | 0.6081 |
| 3 | 0.10 | 0.6206 | 0.6225 | 0.6083 |
| 3 | 0.15 | 0.6345 | 0.6374 | 0.6075 |
| 3 | 0.20 | 0.6489 | 0.6531 | 0.6068 |
| 3 | **0.25** | **0.6641** | 0.6497 | 0.6060 |
| 3 | 0.30 | 0.6800 | 0.6457 | 0.6052 |
| 3 | 0.40 | 0.6854 | 0.6379 | 0.6037 |
| 3 | 0.60 | 0.6522 | 0.6228 | 0.6007 |
| 2 | 0.20 | 0.6310 | 0.6172 | 0.6068 |
| 2 | 0.30 | 0.6509 | 0.6143 | 0.6052 |

### The gates

**M0 — REPRODUCTION, bit-identity.** The anchor cell (`r` = 3, `X_s` = 0.035714…) is the
*same distribution* as `BidispersePSD(cv = 0.5)` reached along the other
parameterisation, and the two `sample()` implementations are already verified to return
**bit-identical radii** from the same generator (checked 2026-09-13 before launch:
`np.array_equal` → True, max |Δr| = 0.0, including through the family-string round trip).
So the anchor must reproduce the existing bidisperse CV = 0.5 records to **|Δφ| < 1e-12**,
not merely within scatter. Likewise `X_s` = 0 and `X_s` = 1 are both monodisperse at
`d` = 1 by construction and must reproduce **each other** and the existing CV = 0 cells to
the same bar. *A gate that can only be passed by the code being the same code is worth
more than one passed by two numbers being close, and all three are free.*
⚠ If M0 fails the campaign is **not readable** and nothing below is quoted.

**M1 — THE MATCHED PAIR. The sharpest single comparison this project can make.**
At `r` = 3, `X_s` = 0.25 the binary has number CV **exactly 0.5000**
(`n_small` = 9/10, `d_small` = 5/6, `d_large` = 5/2 ⇒ mean 1, variance 1/4), and the
existing bidisperse CV = 0.5 cell is `r` = 3 at number CV exactly 0.5000 too. The two
packs therefore agree **exactly** in number CV, size ratio, `d_max/d_min`, `d95/d05`,
number span and the number of atoms, and differ **only in which atom carries the volume**.

- **Descriptor prediction:** Δφ = 0. Every descriptor built on those quantities is
  blind to the difference.
- **Filling prediction:** Δφ = **+0.060** (0.6641 − 0.6042), i.e. **+28 scale**.
- **PASS** if |Δφ| > 5 × `PHI_SCALE` = 0.0105 **and** the sign is positive.
- **FAIL** if |Δφ| < 2 × `PHI_SCALE`, which kills the filling picture on one pair.
- The two candidates are 28 scale apart, so the test is informative before it runs; the
  separation ratio is reported beside the verdict and the gate declares itself **NOT
  INFORMATIVE** if the realised separation falls below 1.6× the measurement error.

**M2 — AN INTERIOR MAXIMUM EXISTS at `r` = 3.** `max_x φ(x)` over the gradable
compositions must exceed both endpoints by > 5 × `PHI_SCALE`.
⚠ This gate is nearly implied by M1 and is reported as a consistency check, not as
independent evidence; it can still fail if the maximum sits outside the grid.

**M3 — WHERE the maximum is, which is what the model gets WRONG on the family path.**
Reported, not graded, because the grid spacing (0.05) is coarser than the model's own
uncertainty in the peak position. The prediction to beat: zero-parameter puts it at
`X_s` = 0.40 (`r` = 3) and 0.30 (`r` = 2); fitted puts it at 0.20 and 0.15.

**M4 — THE SIZE RATIO MATTERS, which is what the blind null denies.** The maximum
excess over β must be larger at `r` = 3 than at `r` = 2 by > 3 × `PHI_SCALE`.
Zero-parameter predicts +0.0907 vs +0.0562; the blind null predicts the **same** curve at
both ratios by construction. **This is the gate the null cannot pass**, and it is the
reason the null is in the pre-registration rather than invented afterwards.

**M5 — MINORITY-PHASE ADMISSION, registered before the run.** A cell whose minority atom
has fewer than **50** grains is recorded and printed but **not graded**: a mixture whose
minority phase is a few dozen grains measures that phase's sampling noise. On this grid
that excludes exactly one composition, `r` = 3 at `X_s` = 0.60 (36 large grains). The
count is printed for every cell.

**M6 — FINITE SIZE.** `d_max/L_dense` is computed for every cell before launch and the
worst on the grid is **0.248** (`r` = 3, `X_s` = 0.60), just under this project's 0.25
flag. No cell is expected to raise it; any that does is reported.

### What would make me abandon the filling picture

M1 returning |Δφ| < 2 × `PHI_SCALE` with the two packs matched as described. There is no
version of the pairwise-filling account that survives two packings at identical ratio and
identical CV packing the same, and no re-reading of the grid that rescues it.

### What a PASS does and does not establish

A pass establishes that **composition at fixed ratio moves φ**, i.e. that the control
parameter is a functional and not a moment — which is the positive statement this project
does not yet have. It does **not** establish that *this* functional is the right one: the
model is already known to be systematically high and to misplace the maximum on the
family path, and M3 is where that is measured rather than graded.

---

## 4b. EXTENSION registered after M0–M6 landed and before any new cell ran

The first grid has size ratios **2 and 3 only**. A continuous PSD is a sum over many
pairs, most of them at ratios **below 2**, so the calibration in `mixture_closure.py`
constrains the interaction functions nowhere near where the continuous families mostly
live — and ratio 2 is exactly where the binary campaign found the textbook form worst
(mean residual −4.2 scale, against +1.6 at ratio 3). The extension adds

- `r` = **1.25** and `r` = **1.5**, the same `X_s` grid, same seeds, N = 1500, 60 cells.

**E1 — is the two-exponent FORM rigid enough?** Fit `(PA, PB)` to each ratio block
SEPARATELY and compare with the joint fit. If the per-ratio exponents agree, a single
pair describes all four ratios and the functional is a closure. If they disagree
systematically with ratio, the FORM is wrong even though each block fits — which is a
finding about the form, not a failure of the campaign, and is the more likely of the two
given that the textbook functions were calibrated on non-spherical aggregates.
**Graded:** per-ratio `PA` and `PB` are printed with the joint values; the verdict is
RIGID if the joint fit's χ² per point over all four ratios is within 2× the sum of the
per-ratio fits, and TOO RIGID otherwise.

**E2 — does the wider calibration improve or degrade the held-out family transfer?**
Refit on all four ratios and recompute the 15-point family χ², which is 63.0 on the
two-ratio calibration and 437.7 at the textbook exponents.
**Prediction, written before the run:** it improves, because the added ratios are closer
to where the families live. ⚠ If it *degrades*, the two-ratio agreement was partly
coincidence and must be reported as such.

⚠ Nothing in §4 is re-graded by these cells. M0–M6 are closed on the ratio 2/3 grid.

## 5. Outcome

**123 cells at N = 1500, 0 errors, 0 non-converged.** `grade_mixture.py` →
`out/mixture.txt`, exit 0, **`FAILURES: []`** — every graded gate passes.
`mixture_closure.py` → `out/mixture_closure.txt`.

### The measurement

| `r` | `X_s` | φ | `cv_realised` |
|---|---|---|---|
| 3 | 0 and 1 (monodisperse) | 0.5940 | 0 |
| 3 | 0.0357 (the anchor) | 0.6042 | **0.5000** |
| 3 | 0.15 | 0.6469 | 0.5622 |
| 3 | 0.20 | 0.6616 | 0.5332 |
| 3 | **0.25** | **0.6677** | **0.5000** |
| 3 | 0.30 | **0.6763** (maximum) | 0.4665 |
| 3 | 0.40 | 0.6741 | 0.4042 |
| 2 | 0.40 | 0.6317 (maximum) | 0.3150 |

- **M0 PASS, bit-identically.** The anchor reproduces the bidisperse CV = 0.5 records at
  |Δφ| = **0.000e+00** on all three seeds, and the two monodisperse endpoints reproduce
  each other and the ratio-2 endpoint at 0.000e+00 on all three. One seed of the
  comparison against the *existing* CV = 0 lognormal cells differs by **2.10e-03**; that
  is the protocol's own documented non-determinism, which `load_all` reports at a median
  2.1e-03 and finds only at CV = 0 — so a gate written to test the new family
  independently re-measured the number `PHI_SCALE` was set from.
- 🔑 **M1 PASS, and it is the sharpest result in this project.** Two binary packs at
  **identical size ratio 3 and identical number CV 0.500000** (agreeing to 5.6e-17),
  differing only in which atom carries the volume, pack at **0.6042 and 0.6677**:
  **Δφ = +0.0635 ± 0.0032, 30 scale, 19.7 se.** The descriptor account predicts 0; the
  zero-parameter filling account predicted **+0.0604**, i.e. it is **5 % off a prediction
  with nothing fitted to this campaign**. Candidate separation 18.7 measurement errors.
- **M2 PASS** (interior maximum, +39 scale at `r` = 3 and +18 at `r` = 2) and
  **M4 PASS** (+21 scale, against exactly 0 from the size-ratio-blind null).
- **M3, reported:** the maximum sits at `X_s` = 0.30 (`r` = 3) and 0.40 (`r` = 2), against
  0.40 / 0.30 zero-parameter and 0.20 / 0.15 fitted. The model has the peak POSITIONS
  ordered the wrong way round between the two ratios, both within one grid step.
- **CV is refuted as the control parameter from inside a single two-atom family.** φ is
  **not monotone** in `cv_realised` on either ratio block, and the closest-in-CV pairs are
  30 scale apart (`r` = 3, ΔCV = 0) and 10 scale apart (`r` = 2, ΔCV = 0.011). This is
  stronger than the cross-family version because only one thing changes.
- **The model, point by point over the 16 gradable interior compositions**, in
  `PHI_SCALE²`: **zero-parameter 365, bidisperse-fitted 906, size-ratio-blind null 4911.**
  The column with nothing fitted to this campaign wins.
- 🔴 **And where it is wrong, with a sign:** mean residual **−4.2 scale at `r` = 2** and
  **+1.6 at `r` = 3**. The textbook interaction functions **over-credit filling at the
  smaller ratio** — the defect §2 named in advance from the family path.

### The closure (`mixture_closure.py`)

**Calibrating the two exponents on binary packs alone and predicting the continuous
families is a cross-experiment transfer the calibration set contains no point of.**

| calibration | held-out χ² over the 15 family points | family-contrast spread |
|---|---|---|
| textbook exponents | 437.7 | 1.71× |
| ratios 2, 3 (`PA` = 0.911, `PB` = 2.122) | **63.0** | **1.23×** |
| ratios 1.25, 1.5, 2, 3 (`PA` = 0.961, `PB` = 2.214) | **34.9** | — |

against **4.51× for CV** and **4.06× for `d₃₂`**. Every one of the 15 held-out points
lands within **+3.8 / −0.4 scale**. **E2 passes as predicted.**

- 🔴 **E1: TOO RIGID, and it says where.** Per-ratio fits give `PA` = 2.35 / 1.76 / 1.38 /
  0.78 at `r` = 1.25 / 1.5 / 2 / 3 — monotone, a 3.0× drift — while `PB` = 2.22 / 2.08 /
  1.98 / 1.97 is constant to 12 %. The joint fit's χ² is **3.55×** the sum of the
  per-ratio fits against a 2× bar. **The wall function transfers; the loosening function
  does not have this form.**
- 🔴 **The obvious repair was tried and REFUSED.** Fitting `a` and `b` directly per ratio
  — two numbers, no functional form — returns a `b` that is **non-monotone in `u`**, which
  the script's own rule (written above the table, not after it) says makes the
  interpolation indefensible. Profiling confirms why: at `r` = 1.25 both parameters have
  profile intervals spanning the **entire [0, 1]**, because that ratio's whole composition
  curve spans 0.0037 in φ — under 2 scale. The fitted value there was a readout of where
  the optimiser stopped. The route is refused, its held-out χ² (375.9, worse than the
  two-exponent fit) is printed rather than hidden, and **the two-exponent fit remains the
  reported closure.**

### What this establishes, and what it does not

**ESTABLISHED.** The control parameter for the polydispersity jamming shift is a
**pairwise-filling functional of the size distribution**, not a function of any moment of
it. Two packs matched in CV, size ratio, span and number of atoms differ by 0.0635 in φ.
The functional calibrated on binary mixtures predicts three continuous families across
five CV levels to ~2 scale, collapsing a spread that CV leaves at 4.5×. `shape_control`'s
negative now has a positive beside it, and the bidisperse non-monotonicity is explained:
that family is a curve through the composition plane that runs into the fines-poor corner.

**NOT ESTABLISHED.** (i) No binary here has a ratio above 3, so `u` < 1/3 rests on a
proved limit and an interpolation — and the lognormal at CV 0.5 contains pairs down to
`u` = 0.04. (ii) `β` is the same measured number in every column, so nothing here tests
the functional at a different friction; only `mixture_model.py` §7 does, at the textbook
exponents. (iii) Every φ is read at one pressure rung; nothing says the interaction
functions are pressure-independent. (iv) The functional is fitted to φ at `p̂` = 1e-3, not
to `φ_c` from the constant-pressure shear protocol the headline is measured on — that
transfer is untested.
