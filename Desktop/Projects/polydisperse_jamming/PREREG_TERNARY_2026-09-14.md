# PREREG_TERNARY_2026-09-14 — is the pairwise-filling functional ADDITIVE?

Written **before any ternary cell existed**. Every number in §3 is printed by
`predict_ternary.py`, which loads only the binary blocks (`phase12_mixr*`,
`phase12_extr*`) and the monodisperse `beta`; the ternary tag is `tern`, which
neither of those globs matches. The artifact is `out/predict_ternary.txt`.

Runner `phase12_ternary.py`, launcher `scripts/run_ternary.sh` → `scripts/launch_ternary.sh`,
grader `grade_ternary.py` → `out/grade_ternary.txt`.

---

## 1. The question, and why a binary cannot answer it

`grade_mixture.py` (123 cells) established that the control parameter for the jamming
shift is a **pairwise-filling functional** of the size distribution, not a descriptor of
it. `grade_mixture_ext.py` (39 cells) extended the binary calibration to size ratios 4
and 5 and measured the two interaction functions directly at
`u = d_fine/d_coarse` = 0.667, 0.500, 0.333, 0.250.

The closure built from those **measured** interactions then transfers to the continuous
families **worse** than the two-exponent form it was meant to improve on: held-out
χ² **108.1** against J's **34.9**, under-predicting the lognormal monotonically with CV
to **−6.9 scale at CV 0.5** (measured φ *above* predicted φ).

Two explanations are exactly confounded in every measurement this project has:

* **(i) calibration.** `a(u)` is wrong at small `u`, where a continuous PSD's widest
  pairs live and where the binaries constrain least.
* **(ii) additivity.** The model sums an interaction over every ordered pair. If three
  or more classes acting together fill *better* than the sum of their separate effects,
  an additive evaluation under-predicts — which is the sign observed.

**A binary pack has exactly one pair.** No binary campaign, at any size ratio, with any
number of seeds, can distinguish (i) from (ii), because the sum has one term. A ternary
is the cheapest object whose prediction depends on how the terms combine.

## 2. Design

**Diameters** `d` = 1, 2, 4 (relative; the family string carries them exactly). The three
ordered pairs are at `u` = 1/2 (1–2), 1/2 (2–4) and 1/4 (1–4).

Both values sit **inside** the binary calibration range (`u` = 0.20 … 0.80), so the
interaction *values* are interpolated, not extrapolated. The only new thing asked of the
model is whether three simultaneous pairs combine the way it says. That is the whole
design.

The ladder is bracketed from both sides and is the only one that meets both constraints:
classes far enough apart for filling to be on (at `u` = 1/2 the measured `a` = 0.785, so
a ratio-2 pair still interacts strongly), and close enough that the coarse class is a
phase rather than a sampling fluctuation, since number fraction goes as `y_i/d_i³`.

**Compositions**, as (`y_L`, `y_M`, `y_S`) — VOLUME fractions of the `d` = 4, 2, 1
classes. Admissibility (`MIN_PHASE_COUNT` = 50 grains, inherited verbatim from
`phase12_mixture.py`) is registered here and applied in the runner, not after the fact:

| `y_L` | `y_M` | `y_S` | `n_L` | `n_M` | `n_S` | CV | `d32` | `d_max/L` |
|---|---|---|---|---|---|---|---|---|
| 0.60 | 0.25 | 0.15 | 74 | 246 | 1180 | 0.5443 | 1.7941 | 0.210 |
| 0.55 | 0.30 | 0.15 | 66 | 287 | 1147 | 0.5246 | 1.7281 | 0.212 |
| 0.55 | 0.35 | 0.10 | 84 | 431 | 985 | 0.5259 | 1.6645 | 0.195 |
| 0.65 | 0.20 | 0.15 | 82 | 203 | 1215 | 0.5649 | 1.8654 | 0.208 |
| 0.70 | 0.20 | 0.10 | 121 | 276 | 1103 | 0.5983 | 1.8710 | 0.187 |
| 0.60 | 0.30 | 0.10 | 96 | 383 | 1021 | 0.5487 | 1.7279 | 0.192 |
| ~~0.50~~ | ~~0.30~~ | ~~0.20~~ | **48** | 229 | 1223 | 0.4923 | 1.6864 | 0.228 | **REFUSED** |

`d_max/L` at the dense state is 0.187–0.212 on all six, **below** `packing.py`'s own
0.25 finite-size flag, computed at `beta` (a lower bound on the ternary φ, so this is the
worst case). The refused composition is refused *before* the run, not dropped after it.

> **CORRECTION appended 2026-09-14, after the run; the original wording is kept above.**
> "a lower bound on the ternary φ, so this is the worst case" has the direction
> **backwards**. `L³ ∝ V_solid/φ`, so a *lower* φ gives a *larger* cell and therefore a
> *smaller* `d_max/L` — evaluating at `beta` is the BEST case, not the worst. The solver's
> own recorded `d_max/L` is **0.1997–0.2256**, above my estimate, exactly as the corrected
> reasoning requires. **All six are still below the 0.25 flag and no cell carries a
> finite-size note (21 of 21), so the admission decision is unchanged**; only my
> justification for it was wrong.
>
> Related, and it matters more: `d_max/L` turns out to be rank-correlated with the T1
> residual at **spearman = +1.0000**, as `n_L` is at −1.0000. Both are collinear with
> `y_S` by construction at fixed N. `PREREG_TERNARY_N_2026-09-14.md` breaks both at once,
> since doubling N doubles `n_L` and multiplies `d_max/L` by 2^(−1/3) = 0.794.

**Protocol**, byte-for-byte the production one: `build_and_compress`, `p_hat` = 1e−3,
`N` = 1500, `mu_p` = 0.3, seeds 1/2/3, `truncate=None`. Plus a **monodisperse control**
through the same `DiscretePSD` code path. 7 compositions × 3 seeds = **21 cells**,
12 shards.

**Tag `tern`**, chosen NOT to match `phase12_mixr*` or `phase12_extr*`. M0–M6, F1–F3 and
X0–X5 are closed and this block cannot re-decide any of them.

**`DiscretePSD` gate, passed before launch** (`check_discrete_psd.py`, 34 checks,
0 failures): a two-class `DiscretePSD` returns the **same radii array**, element for
element (`np.array_equal`), as the `BinaryPSD` it is the same distribution as, on 8
(ratio, composition) cases, and all six derived scalars compare **exactly equal**. The
first implementation failed this at 1e−16 on 3 of 8 cases — the same trap
`phase12_mixture.py` records for a rounded family string — and was rewritten to use
`BinaryPSD`'s own arithmetic (weight `y·(d_max/d_i)³`, coarsest class taking
`1 − Σ rest`). A three-class case is checked against hand-computed closed forms, because
bit-identity at K = 2 is degenerate in every way the K = 3 path could be wrong.

## 3. Blind predictions

`beta` = **0.594687** (measured monodisperse φ, `mu_p` 0.3, `p_hat` 1e−3, N 1500;
n = 3, sem 9.7e−4). `PHI_SCALE` = **2.1e−3**, the protocol's own reproducibility,
declared in `mixture_model.py` and NOT the seed sem.

**J** (the shipped two-exponent closure, fitted on the ratio ≤ 3 binaries):
PA = 0.9614, PB = 2.2144.
**TAB** (the tabulated/nonparametric closure, `grade_mixture_ext.tabulated_closure`):
identifiable (a, b) at `u` = 0.2500 → (0.3010, 0.4796), 0.3333 → (0.5210, 0.5505),
0.5000 → (0.7854, 0.7471), 0.6667 → (0.9249, 0.8985), plus the proved limits (0,0) at
`u` = 0 and (1,1) at `u` = 1, joined by PCHIP.

| `y_L` | `y_M` | `y_S` | **J** | **TAB** | ideal | J−β (scale) | TAB−β (scale) |
|---|---|---|---|---|---|---|---|
| 0.60 | 0.25 | 0.15 | 0.6591 | **0.6562** | 0.9072 | +30.7 | +29.3 |
| 0.55 | 0.30 | 0.15 | 0.6560 | **0.6525** | 0.9072 | +29.2 | +27.5 |
| 0.55 | 0.35 | 0.10 | 0.6452 | **0.6449** | 0.8783 | +24.1 | +23.9 |
| 0.65 | 0.20 | 0.15 | 0.6623 | **0.6599** | 0.9072 | +32.2 | +31.1 |
| 0.70 | 0.20 | 0.10 | 0.6545 | **0.6558** | 0.8496 | +28.5 | +29.1 |
| 0.60 | 0.30 | 0.10 | 0.6483 | **0.6485** | 0.9054 | +25.5 | +25.6 |
| 1.00 | 0.00 | 0.00 | 0.5947 | 0.5947 | 0.5947 | +0.0 | +0.0 |

`ideal` is the zero-parameter perfect-separation ceiling (all interactions off). It is
not a candidate; it is printed so the two candidates can be read against the most any
filling picture could give.

### The non-additive alternative, so the test has two arms

**SEQUENTIAL (nested-binary) filling.** The two coarsest classes pack first; the finest
then fills the voids of that **composite** as if it were one medium. It uses the SAME
measured binary closure, adds **no free parameter**, and reduces **exactly** to the
binary closure at K = 2 — which is precisely why the binary calibration cannot tell it
from the additive model, and why a ternary can.

> **CORRECTION appended 2026-09-14, after the run; the original wording is kept above.**
> "reduces **exactly** to the binary closure at K = 2" was asserted **by inspection and
> was FALSE of the code as shipped**. `sequential()` indexed the composite as the
> hardcoded pair `(d[1], d[2])`, but `DiscretePSD` **drops** a zero-volume class, so a
> two-class input arrived with two diameters and raised `IndexError` — the construction
> was not merely wrong at K = 2, it was **undefined** there. A physics-referee pass noted
> that this function had no wiring gate (unlike J and TAB); adding that gate is what
> surfaced the crash, and `grade_ternary.py` had been **exiting rc = 1** while leaving a
> stale artifact on disk that still read PASS. The composite is now written as "every
> class coarser than the finest", with a single-coarse-class branch that packs at `beta`
> and presents its own diameter, so the reduction holds and is **tested** rather than
> asserted.
>
> **No T4 number moves.** `sequential()` is evaluated in T4 only on the six genuine
> three-class compositions, where the old and new forms are algebraically identical; the
> defect was reachable only from the K = 2 gate that did not previously exist.

It needs one convention: what diameter the composite presents to the fines. That
convention is a knob, so **both** defensible choices (`d32`, `d43`) are computed and the
test is declared not informative if they bracket the additive prediction. They do not —
both lie **above** ADD on all six compositions:

| `y_L` | `y_M` | `y_S` | ADD (TAB) | SEQ `d32` | SEQ `d43` | SEQ−ADD (scale) |
|---|---|---|---|---|---|---|
| 0.60 | 0.25 | 0.15 | 0.6562 | 0.6863 | 0.6961 | +14.4 … +19.0 |
| 0.55 | 0.30 | 0.15 | 0.6525 | 0.6851 | 0.6949 | +15.6 … +20.2 |
| 0.55 | 0.35 | 0.10 | 0.6449 | 0.6644 | 0.6704 | +9.3 … +12.1 |
| 0.65 | 0.20 | 0.15 | 0.6599 | 0.6816 | 0.6902 | +10.3 … +14.4 |
| 0.70 | 0.20 | 0.10 | 0.6558 | 0.6608 | 0.6659 | +2.4 … +4.8 |
| 0.60 | 0.30 | 0.10 | 0.6485 | 0.6704 | 0.6767 | +10.4 … +13.4 |

### Candidate separation, stated before the run

* **ADD vs SEQ** (the arms that matter), taking each composition's LARGER separation:
  min **4.8**, median **13.9**, max **20.2** scale. The smallest separation under the
  `d32` convention is **2.4** scale, at (0.70, 0.20, 0.10). Against a three-seed standard error of ~0.6–1.2 scale this is resolvable on
  five of six compositions. **(0.70, 0.20, 0.10) separates by only 2.4 scale and is
  registered here as unable to decide on its own.**
* **J vs TAB**: rms **1.04** scale, max 1.66. **This comparison is DECLARED NOT
  INFORMATIVE BEFORE THE RUN.** J and TAB are two spellings of the same additive
  assumption and they agree to within the seed error on this ladder; grading one against
  the other would be a test that both candidates pass. It is printed, and it decides
  nothing. TAB is the **primary** candidate (it is the only one calibrated at `u` = 0.25,
  which this ladder contains); J is the **control**.

## 4. Gates, with bars declared

Bars come from the **in-sample residuals on the binary calibration set** — a ternary must
be predicted as well as the binaries the closure was calibrated on. Anything else is a
bar invented for this block.

* J on its own calibration set (32 points, ratios ≤ 3): rms **2.83**, max **6.69**
  (this reproduces `PREREG_MIXTURE_EXT`'s recorded 2.83 / 6.7 — a reproduction check on
  the whole prediction path, run before anything new was computed).
* J on all 41 binaries including ratios 4 and 5: rms 4.31, max 14.19.
* TAB on all 41 binaries: rms **2.02**, max **5.74**.

| Gate | Statement | Bar |
|---|---|---|
| **T0** | Monodisperse control through the `DiscretePSD` path lands on `beta` | \|dev\| ≤ **2.0** scale (X0's bar, verbatim) |
| **T1** | TAB (additive, primary) predicts the six ternaries | rms ≤ **2.8** scale (1.4 × its in-sample 2.02, the same 1.4 X1 used) AND max \|res\| ≤ **5.74** (its in-sample max, verbatim) |
| **T2** | J (additive, control) predicts the six ternaries | rms ≤ **4.0** scale (X1's bar, verbatim) AND max \|res\| ≤ **6.69** |
| **T3** | **Sign.** If additive, the six signed TAB residuals are mixed | ADDITIVE ⇔ NOT 6/6 same sign, AND \|mean signed\| ≤ **2.8** scale. NON-ADDITIVE ⇔ 6/6 positive AND mean signed ≥ **+2.8** |
| **T4** | **ADD vs SEQ discrimination.** rms residual ratio loser/winner | declared **NOT INFORMATIVE below 1.6**; the winner is named either way |
| **T5** | J vs TAB | **pre-declared NOT INFORMATIVE** (separation 1.04 scale). Printed, decides nothing |
| **T6** | *Exploratory, not in FAILURES.* Three-class structure: per-design-class `z` and rattler fraction, which no binary can show | recorded |

**Registered prediction, written before the run.** The continuous-family under-prediction
is real and sizeable (−6.9 scale at CV 0.5, χ² 108 vs 35), so I expect **T1 to FAIL with
6/6 positive residuals** (additivity is the defect) with probability I would put a little
above even. I expect **T2 to FAIL** regardless, because this ladder contains `u` = 0.25
and X1 already established J under-predicts there by +7.5 rms — T2's failure carries no
information about additivity and is registered as such now, so it cannot be quoted later
as if it did.

**A fail is the model's range of validity, not a refit.** No exponent, no interaction
value and no bar is touched after the data lands. If T1 fails, the statement is that the
pairwise-filling functional is **not additive beyond two classes**, with the measured
size of the failure. If T1 passes, the statement is that additivity holds on this ladder
and the continuous-family residual is a **calibration** problem at small `u` — which
names a different next experiment.

## 5. What is NOT tested here

* One diameter ladder, one friction, one pressure. A ternary at `d` = 1, 2, 8 would put
  a pair at `u` = 1/8, outside the calibration range, and would confound additivity with
  extrapolation again. That is deliberately not this block.
* Four or more classes. If additivity fails at K = 3 the size of the failure at K = 3 is
  what is measured; nothing here says how it grows.
* The SEQUENTIAL alternative is a named, parameter-free competitor, **not** a proposed
  replacement closure. If it wins, the finding is that the combination rule is
  sequential-like on this ladder, and it needs its own calibration campaign before it is
  a closure.
