# PRE-REGISTRATION — the first-passage clock, and whether the approach to `phi_c` is CV-dependent

**Written 2026-09-13, before any cell of this extension ran.** It is the Status card §3
item 1: close the gap `PREREG_STRAINCLOCK_2026-09-13.md` §7 left open — *whether `gamma_c`
diverges at the critical line is undecided, not answered* — by **fixing the observable,
not the ladder**.

---

## 1. What is already known, and which part of it is EXPLORATORY

The strain-clock campaign registered `gamma_c`, the strain to reach **half of the cell's
own `p_ss`**, and that observable is **not well posed** at `x` <= +0.10: 0 of 16 jammed
cells readable, first and last half-crossings 5–15x apart. Its G4 verdict (COLLAPSE)
therefore rests on `x` >= 0.35 and says nothing about the approach to `phi_mu`.

After that failure I looked for an observable that does not normalise by `p_ss` — because
`p_ss` is itself collapsing as `phi -> phi_c`, so `gamma_c` chases a moving target — and
found one already computable from the 45 committed records:

> **`t_fp(p_thr)`** — the strain at which `p_hat` **first** reaches a **FIXED** threshold
> `p_thr`, by linear interpolation, with **`frac_above`** = the fraction of post-passage
> samples still at or above `p_thr`.

**Everything in this section is EXPLORATORY.** The estimator was chosen after seeing the
registered one fail, and the numbers below are in-sample. They are written down here so
the out-of-sample test in §3 can be stated against them, not to be quoted as results.

Exploratory finding, at `p_thr` = 1e-04 (`out/firstpassage_explore.txt`):

- It **is** well posed where `gamma_c` is not: at `x` = +0.10, `frac_above` = 0.93–1.00
  on 9 of 9 cells, i.e. the pack reaches the threshold and stays there.
- It rises steeply toward the critical line where `gamma_c` was flat. At
  `p_thr` = 3e-04 the geometric means run 0.0184 -> 0.0354 -> 0.1558 (CV 0) over
  `x` = 0.60 -> 0.35 -> 0.10, an 8.5x rise.
  ⚠ **Corrected the same day, before any new cell landed: most of that rise is NOT
  kinetic.** A fixed threshold read against a `p_ss` that is itself collapsing measures
  the collapse — this project's own standing lesson, in mirror image. Measured:
  `p_ss ~ Delta phi^b` with `b` = 1.033 / 1.114 / 1.218 (CV 0 / 0.2 / 0.5), against
  `a` = 1.219 / 1.262 / 1.446, so **`a/b` = 1.180 / 1.133 / 1.188** and
  `d ln t_fp / d ln p_ss` = -1.180 / -1.131 / -1.194. A universal build curve
  `p(gamma) = p_ss f(gamma/gamma*)` with `gamma*` independent of `phi` predicts exactly
  -1; the measured -1.13 to -1.19 is the small departure from it, and it is the SAME at
  every CV. The first draft of this section said the estimator "shows the divergence
  `gamma_c` was blind to", which overstates it: the two observables are two readings of
  one picture, and F6 below is what separates the kinetic part from the pressure scale.
- Over a decade in `Delta phi` = `phi - phi_mu` it is a **power law**,
  `t_fp ~ Delta phi^(-a)`, with max |residual| in ln of 0.06–0.32:

| CV | `a` at `p_thr` = 1e-04 | `a` at 3e-04 | worst seed sd(ln) |
|---|---|---|---|
| 0.0 | 1.219 | 1.138 | 0.18 |
| 0.2 | 1.262 | 1.010 | 0.30 |
| 0.5 | 1.446 | 1.814 | 0.65 |

- 🔑 **And the CV 0.5 arm does not extrapolate.** Fitting on the three far rungs
  (`x` = 0.35, 0.60, 0.85) and predicting `x` = 0.10 gives, **on the primary threshold
  1e-04**, obs/pred = **1.14 (CV 0), 0.94 (CV 0.2), 0.42 (CV 0.5)**; on 3e-04 the same
  three are 1.28 / 0.90 / **0.30**. Either the broad pack's approach to critical really is
  different — which would put **M4 back, in exactly the region G4 had no data** — or CV
  0.5's three-seed scatter (sd(ln) up to 0.65) makes its far-rung fit meaningless. Three
  seeds cannot tell those apart. That is what this campaign is for.

---

## 2. The grid, fixed here

**54 new cells**, same protocol, same driver, same rate, same `n_samples` = 1500, N = 1000,
`mu_prep` = 0: `phase10_strainclock.py --x-rungs ... --seeds ... --tag fp`. The cell key
carries `phi` and `seed`, so these cannot collide with the 45 strain-clock records, and
the 45 are **reused unchanged** — same protocol string, same everything.

| what | rungs `x` | seeds | cells |
|---|---|---|---|
| a NEW rung nearer critical | +0.05 | 1–6 | 18 |
| more seeds at the existing rungs | +0.10, +0.35, +0.60, +0.85 | 4–6 | 36 |

giving **6 seeds at every rung from `x` = +0.05 to +0.85** when pooled with the committed
records. `phi(x = 0.05)` = 0.5932 / 0.6042 / 0.6467 at CV 0 / 0.2 / 0.5.

⚠ **`x` = +0.05 is `Delta phi` = 0.0021–0.0023, which at CV 0.5 is INSIDE that level's own
`phi_c` seed scatter (sd 0.0043).** So the CV 0.5 point at `x` = 0.05 may not be where the
label says it is, and no verdict below rests on it alone; it is a check on the law's shape,
and the CV 0 / 0.2 points (whose `phi_c` sd is 0.0010 / 0.0013, i.e. half the offset) carry
the out-of-sample test. Stated now rather than discovered later.

⚠ `x` = +0.85 had 2 of 9 cells **prepared jammed** in the strain-clock campaign
(`p_hat_prep` up to 1.96e-04, the `phi_0(CV 0.5)` uncertainty). Those cells are excluded
from any fit by the admission rule below, and the count is printed.

---

## 3. Predictions, with bars. Written before the first new cell.

Primary threshold **`p_thr` = 1e-04** (chosen because it is the one where all 9 cells at
`x` = +0.10 already satisfy `frac_above` > 0.9, i.e. the estimator is admissible on the
whole ladder; `p_thr` = 3e-04 = `p*` is reported beside it throughout and its CV 0.5
`x` = 0.10 refusal stands as registered).

**F0 — the campaign is what it says it is.** 54/54 records, `protocol` =
`cv-shear-v3-fprep`, `friction_switch_ok` on all, `mu_prep` = 0, `n_samples` = 1500,
N = 1000, `gdot` spread across CV below 1.2x. Pooled with the committed records: 6 seeds
at each of the 5 rungs at each of 3 CVs, 90 cells.

**F1 — the admission rule, fixed here.** A cell enters a fit only if
`frac_above` >= **0.90** and `p_hat_prep` < 1e-06. Cells failing either are printed with
their numbers and excluded. PASS requires >= 85 % of cells at `x` >= +0.05 admitted; below
that the estimator is not usable near critical and the campaign says so.

**F2 — THE OUT-OF-SAMPLE TEST.** Fit `ln t_fp = ln A - a ln(Delta phi)` on
`x` = 0.35, 0.60, 0.85 **only**, per CV, on the 6-seed data. Predict `x` = **+0.05**, a
rung no fit has seen.

- **PASS** if `obs/pred` at `x` = 0.05 lies in **[0.5, 2.0]** at CV 0 and CV 0.2 (the two
  levels whose `phi_c` scatter is well inside the offset). Two of two required.
- CV 0.5 is **reported, not graded** — §2's `Delta phi` caveat.
- A factor 2 is the bar because the extrapolation is 0.85 -> 0.05 in `x`, a factor 7 in
  `Delta phi`, and the in-sample seed sd(ln) is already 0.18–0.30, i.e. a factor 1.2–1.35
  per cell mean.

**F3 — does the EXPONENT depend on CV?** This is G4's question asked where G4 had no data,
and it is the M2a/M4 discriminator in its strongest form: if polydispersity acts only
through `phi_c`, then plotted against `Delta phi` the approach law must be the SAME curve.

- **COLLAPSE (M2a)** if `a`(CV 0.2) and `a`(CV 0.5) each lie within **2 se** of `a`(CV 0),
  with `se` from the per-rung seed scatter, AND the prefactor ratio `A(CV)/A(0)` lies in
  [0.5, 2.0].
- **DOES NOT COLLAPSE (M4 lives near critical)** if either exponent differs by more than
  2 se and the difference is larger than the threshold-to-threshold spread of `a` within
  one CV (which is 0.08–0.37 on the exploratory data, and is a systematic of the
  estimator, not of the material — so it is the floor below which no exponent difference
  means anything).
- **UNDETERMINED** otherwise, and this is a real verdict.

**F4 — does the CV 0.5 extrapolation failure survive six seeds?** Pre-registered,
two-way, on the number that motivated the campaign: the exploratory obs/pred at
`x` = 0.10, CV 0.5 was **0.42** at the primary threshold and **0.30** at 3e-04, both
outside [0.5, 2.0].

> ⚠ **Corrected 2026-09-13 before any new cell landed.** The first draft of §1 and of this
> paragraph quoted "1.14 / 0.94 / 0.30", which SPLICES two thresholds — the first two are
> 1e-04 and the third is 3e-04. Caught by running `grade_firstpassage.py --explore`, which
> computes both. The splice did not change the verdict (0.42 is outside the band too), and
> it is recorded because a number welded from two settings is exactly the defect this
> project has been bitten by before.

- **It was scatter** if the 6-seed obs/pred at `x` = 0.10, CV 0.5 lands in [0.5, 2.0].
- **It is real** if it stays below 0.5.
- Either way the CV 0 and CV 0.2 values (1.14, 0.94) must stay in [0.5, 2.0] or the
  estimator itself is in question and F2–F4 are all withdrawn.

**F6 — THE PROPERLY NORMALISED COLLAPSE, added 2026-09-13 before any new cell landed
and after the exploratory check in §1.** `a` inherits the CV dependence of `b`, so
comparing `a` across CV (F3) partly re-measures `b`. The scale-free statistic is

    `a/b`  equivalently  `d ln t_fp / d ln p_ss`,

fitted per CV on the same `x` = 0.35, 0.60, 0.85 rungs, with `b` from `p_ss` on those
same cells so the two exponents share their sample exactly.

- **COLLAPSE** if `a/b` at CV 0.2 and CV 0.5 each lie within **2 se** of CV 0's, se by the
  same seed bootstrap, at BOTH thresholds.
- **DOES NOT COLLAPSE** if either differs by more than 2 se at both thresholds.
- **UNDETERMINED** if the thresholds disagree.
- Reported beside it: `b` itself per CV, because `b`(CV) is a statement about the
  critical-state pressure scale and is worth having whichever way F6 falls.

⚠ F3 is **kept and graded as registered** — it is the question as it was asked, and its
number stays printed. F6 is the sharper test beside it, not a replacement, and if the two
disagree that disagreement is the result.

**F5 — reported, not graded.** `t_fp` at both thresholds side by side; `frac_above` per
rung; the ratio `t_fp(3e-04)/t_fp(1e-04)` per cell, which is a shape statistic and must be
CV-independent if the pressure builds by the same trajectory everywhere; and `t_fp`
against `x` as well as against `Delta phi`, since the two coordinates differ by `W(CV)`
and the collapse could hold in one and not the other.

---

## 4. What each outcome buys

| outcome | meaning |
|---|---|
| F2 PASS, F3 COLLAPSE, F4 "scatter" | the approach to critical is one law, CV enters it only through `phi_c`, and the strain-clock G4 verdict extends all the way to the critical line. M4 is dead on both observables. |
| F2 PASS, F3 NOT COLLAPSE | the *steady* state collapses (G4) and the *approach* does not — M4 lives, but only near critical, which is a sharper claim than either campaign could have made alone and is directly testable on the bench as a `phi`-dependent onset time |
| F4 "real", F3 UNDETERMINED | the CV 0.5 arm is anomalous and not yet explained; the next step is seeds, not rungs |
| F1 FAIL | the first-passage estimator is not usable near critical either, and the honest end state is that this DEM cannot resolve the approach to `phi_c` at this N and rate |

---

## 5. OUTCOME

**Landed 2026-09-13. 54 new cells, 0 failures; pooled with the 45 committed strain-clock
records that is 99 cells and 6 seeds at every rung from `x` = +0.05 to +0.85.** Every
number below is printed by `grade_firstpassage.py` -> `out/firstpassage.txt`, whose
verdict line is `FAILURES: ['F2']`.

### 🔑 F6 — the properly normalised clock COLLAPSES, at every CV and both thresholds

| | `a` | `b` (`p_ss ~ dphi^b`) | `a/b` |
|---|---|---|---|
| CV 0.0 | 1.229 ± 0.120 | 1.021 | **1.203 ± 0.131** |
| CV 0.2 | 1.357 ± 0.131 | 1.048 | **1.295 ± 0.116** |
| CV 0.5 | 1.464 ± 0.377 | 1.063 | **1.378 ± 0.191** |

(primary threshold; at 3e-04 `a/b` = 1.055 / 1.059 / 1.184.) Against CV 0 the differences
are **0.5, 0.8, 0.0 and 0.7 se** across the four (CV, threshold) combinations — COLLAPSE
on all four, and the verdict is the same at both thresholds for both CVs.

**So the approach to the critical line is one law, and CV enters it only through the
pressure scale.** This is G4's question asked where G4 had no data, and it gives G4's
answer. **M4 is now excluded on both observables and on both sides of the ladder.**

Two things make `a/b` the statistic rather than `a`, and both are measured:

- **`a` alone is convention-dependent by the size of the effect it is being asked to
  resolve.** Subtracting this protocol's own `phi_c` (§below) instead of the
  constant-pressure `phi_mu` moves `a` by **0.16–0.30** — comparable to the CV
  differences F3 tests — and moves `a/b` by **0.000–0.004**. Both exponents are
  `d ln(.)/d ln(dphi)` on the same abscissa, so moving the origin rescales them by the
  same local factor and it cancels.
- **`b` is 1.02 / 1.05 / 1.06**, i.e. `p_ss` is linear in `dphi` at every CV. The
  CV-dependence `a` appears to have is mostly inherited from `b` and from the origin.

### F3 — as registered: COLLAPSE at CV 0.2, UNDETERMINED at CV 0.5

`da` = +0.128 ± 0.178 and +0.032 ± 0.114 at CV 0.2 (both thresholds COLLAPSE); at CV 0.5
the two thresholds disagree (UNDETERMINED at 1e-04 on the prefactor band, COLLAPSE at
3e-04), so the registered floor rule returns **UNDETERMINED**. F3 is kept and printed; F6
is the sharper test beside it and they do not conflict — F6 resolves what F3 could not,
which is what a normalisation is for.

### 🔴 F2 FAIL — and the failure LOCATES a 0.004 disagreement between two `phi_c` estimates

Fitting on `x` = 0.35/0.60/0.85 and predicting `x` = +0.05 — a factor 7 in `dphi`, a rung
no fit saw — gives obs/pred = **0.51 (CV 0), 0.49 (CV 0.2), 0.27 (CV 0.5)** at the primary
threshold. The bar was [0.5, 2.0] at CV 0 and 0.2, and CV 0.2 misses it by 0.01. **FAIL as
registered.**

The cause is measured, not argued. This campaign has two estimates of `phi_c`: the
constant-pressure `phi_mu` the ladder was placed on, and the crossing this protocol itself
returns (`PREREG_STRAINCLOCK` §7 G5), which is **0.0035–0.0046 lower** — at `x` = +0.05
a factor **2.5–3.0** in `dphi`. Re-running the same extrapolation with that origin:

| | obs/pred at `x` = +0.05, `phi_mu` | on the protocol's own `phi_c` |
|---|---|---|
| CV 0.0 | 0.51 | **1.06** |
| CV 0.2 | 0.49 | **1.04** |
| CV 0.5 | 0.27 | 0.68 |

**The power law extrapolates a factor 7 in `dphi` to within 6 % at CV 0 and 0.2 once the
origin is the one this protocol measures.** So F2 did not find a saturation of the approach
law; it found that the two `phi_c` estimates differ, and it is sensitive to that difference
because `x` = +0.05 sits 2–3 `phi_c`-discrepancies from critical.

⚠ **How independent is that?** The G5 crossing is read from `p_ss` on these same cells, so
it is not independent of `b`. The part that is: `t_fp` is a *different observable* from
`p_ss`, and its extrapolation prefers the *same* origin.

### F4 — the number that motivated the campaign was mostly three-seed scatter

The exploratory obs/pred at `x` = +0.10, CV 0.5 was 0.42 (primary) / 0.30 (3e-04) on three
seeds. On six: **0.86 / 0.79 / 0.66** at the primary threshold (all "scatter" by the
registered rule) and 1.08 / 0.95 / **0.48** at 3e-04, where CV 0.5 is marginally "REAL".
So the anomaly that prompted this campaign largely dissolved with three more seeds — which
is the outcome the two-way registration existed to allow, and the reason it was registered
two-way.

### F0, F1 — the campaign is what it says it is, and the estimator reaches the new rung

99 cells, all `cv-shear-v3-fprep`, all with `friction_switch_ok`, `gdot` spread 1.016
against 1.2, 6 seeds at every rung at every CV. **F0 PASS.**
Admission at the primary threshold: **84/90 (93 %)** of the cells at `x` >= +0.05, against
an 85 % bar — **F1 PASS**. The estimator works within `dphi` = 0.002 of the critical line,
where the half-rise clock was unreadable on 16/16. The six refusals are 5 `x` = 0.85 cells
at CV 0.5 prepared JAMMED (the `phi_0` uncertainty §2 flagged) and one `frac_above` = 0.90
boundary case.
At 3e-04 admission is 81 %, below the bar — because `p_ss` at `x` = +0.05 is ~5e-04 and a
3e-04 threshold is then 60 % of steady, so `frac_above` bites. That is why 1e-04 was named
primary in §3 *before* the run, and the 3e-04 column is reported throughout rather than
graded.

### What this changes

1. **The Status card §3 item 1 is discharged.** Whether the clock diverges at the critical
   line was undecided; it does diverge, as `dphi^(-1.2)`, and **`a/b` = 1.06–1.38 shows
   that divergence is the pressure scale plus a small, CV-independent kinetic excess over
   the universal-curve value of 1.**
2. **M4 is excluded near critical as well as far from it.** The strain clock collapsed at
   `x` >= 0.35; the normalised first-passage clock collapses down to `dphi` = 0.002.
3. **A new number for the next campaign:** the two `phi_c` estimates differ by
   0.0035–0.0046 and F2 is the first test sensitive enough to care. The four-pressure
   constant-pressure ladder (Status card §3 item 2) is what settles which is right.
