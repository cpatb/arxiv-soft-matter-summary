# CLAIMS_ADDITIVITY_2026-09-17 — the additivity repair, and the K = 6 interpolation

Graded against `PREREG_ADDITIVITY_2026-09-17.md`, committed at `1c90457` **before
`grade_additivity.py` existed**. **No new cells were run.** Every number below is printed
by `out/grade_additivity.txt` or `out/diag_k6_interp.txt`; `test_additivity.py` (20 tests)
asserts that none of them can drift from its artifact.

`FAILURES: ['P0 calibration arms disagree (0.5214 vs 0.4969)', 'P2b mix 56.5 does not
beat SEQd21 42.8', 'P3 held-out 56.5 >= J 34.9']`

Sign convention: **residual = MEASURED − MODEL**; positive = the model under-predicts φ.
`transfer_kladder.py` uses the opposite sign and the grader **asserts** the two agree
after a flip (−1.6946 both ways).

---

## A. The result: one λ describes both DISCRETE ladders and does not carry to continuous PSDs

This is the whole item in one table, at the registered primary convention `d43`:

| set | λ | ±1 | cells | with leverage |
|---|---|---|---|---|
| ternary (calibration) | **0.5214** | 0.0274 | 6 | 6 |
| K-ladder (calibration) | **0.4969** | 0.0228 | 19 | 14 |
| both ladders (the primary fit) | **0.5069** | 0.0175 | 25 | 20 |
| continuous families (**HELD OUT**, oracle) | **0.3033** | 0.0313 | 15 | 10 |

* ternary vs K-ladder: |Δ| = 0.0245, **t = 0.69**
* ladders vs families: |Δ| = 0.2036, **t = 5.67** — **8× the between-ladder gap**

Two datasets built on different principles — a three-class composition simplex and a
geometric size ladder at four class counts — want the **same** amount of nesting. The
continuous PSDs want a **different** one, by eight times the disagreement between the two
that agree. That is the supportable statement, and it does not depend on any χ² bar.

⚠ **Read the table only with the leverage column.** At K = 2 the nested rule reduces
**exactly** to the additive binary closure, so `S − A ≡ 0` and λ has literally no leverage
there: **5 of the 15 held-out points (the entire bidisperse family) are λ-independent**,
and **5.4 of ADDtab's 108.1 is a floor no mixing model can go below.**

## B. The bracket premise holds; the transfer does not

* **P1 PASS.** λ̂ = **0.5069 [0.4894, 0.5245]** at `d43` — strictly interior to (0, 1)
  with its interval inside it, fitted **unbounded** so a truth outside the bracket could
  have shown itself. The ADD/SEQd43 bracket contains the answer on the ladders.
* **P2 PASS**, held-out χ² **56.5 < 108.1**. Registered as **weak** and it is: `SEQd21`
  already manages 42.8 with no fitting at all.
* **P2b FAIL**, the discriminating bar. The convention that wins on the **calibration**
  set — selected without touching one held-out point — is `SEQd21` (calibration rms 2.64
  against d43's 5.80), and its held-out χ² is **42.8**. The fitted mix's **56.5 does not
  beat it.** One free parameter bought less than a convention choice that costs none.
* **P3 FAIL**, the headline: **56.5 ≥ J's 34.9**.
* **P4 PASS**, 3 of 4 conventions clear P2 (`d10` does not, at 178.9).

### The failure is a TRANSFER failure, not a FORM failure — and the oracle is what says so

The **oracle** λ — the value that minimises χ² **on the held-out points**, which is *not a
result* and must never be quoted as one — reaches **14.2** at `d43`, comfortably below J's
34.9. So no property of the functional form stops this rule from beating J. What stops it
is **where the calibration puts λ**: the ladders want 0.5069 ± 0.0175 and the families
want 0.3033, a gap of **11.6 half-widths**. Refitting is not the answer, and neither is a
better ladder dataset.

## C. Mixing made the composite-diameter convention WORSE, not better

| | spread of held-out χ² | factor |
|---|---|---|
| unmixed SEQ (recorded, `transfer_kladder.txt`) | 42.8 → 509.3 | **11.90×** |
| **the fitted mix** | **11.6 → 178.9** | **15.38×** |

λ̂ by convention: d10 +0.926, d21 +0.898, d32 +0.671, d43 +0.507.

`d32` mixed reaches **11.6**, three times better than J — and it is **not quotable**,
because it is the best of four undeclared conventions and the spread is larger than the
effect. The undeclared composite-diameter choice has now been worth 11.9×, 28.95× (U6)
and 15.38× on three separate occasions; it is the dominant term in every nested result
this project has produced.

## D. Where the improvement actually is, per family

| family | n | ADDtab | J | MIX (d43) | SEQd21 | MIX mean res |
|---|---|---|---|---|---|---|
| lognormal | 5 | 84.4 | **6.7** | 47.8 | 32.0 | −2.35 |
| uniform | 5 | 18.3 | 9.8 | **3.3** | 5.4 | +0.31 |
| bidisperse | 5 | 5.4 | 18.4 | 5.4 (floor) | 5.4 | +0.15 |

The mix essentially **solves the uniform family** and leaves the lognormal's monotone
drift with CV intact (−0.82, −0.09, −1.77, −3.40, −5.70 across CV 0.1→0.5). A
one-parameter mix **cannot** cure a monotone drift, because both endpoints of the bracket
drift the same way. **J and MIX fail on complementary families** — J is three times worse
than ADDtab on the bidisperse pack and seven times better on the lognormal — so "J wins"
is a statement about the pooled χ², not about a family.

## E. P0 failed on the sentence that was registered, and it stays failed

P0 required each arm's point estimate to lie inside the other's Δχ² = 1 interval **and**
the intervals to intersect. The intervals intersect by +0.0258; the ternary's 0.5214 sits
**0.0017 outside** the K-ladder's upper end. **FAIL, recorded.**

Two sharper statistics are **printed beside it and neither replaces it**, and both point
the loosening way, which is the direction to distrust:

* |Δλ| = 0.0245 against the se of the **difference** 0.0357 → **t = 0.69**. The registered
  sentence compares a point to the *other* arm's interval, which is a stricter and less
  standard test than the two-sample one.
* χ²_cal = **50.3 on 24 dof**, i.e. χ²_red = **2.10**: the fit is dominated by **model
  misfit, not seed noise**, so a Δχ² = 1 interval in `PHI_SCALE` units is "one scale unit
  of misfit", **not a confidence interval**. Inflated by √χ²_red = 1.45× the two arms are
  [+0.4817, +0.5611] and [+0.4639, +0.5299] and overlap comfortably.

The honest reading is that the two ladder arms **agree**, and that **the gate's wording,
not the data, produced the failure** — which is why section A states the result as a
t-statistic and not as a P0 pass.

## F. The reproduction gate, and what it cleared

14 of 14 rows OK before anything else ran, exiting non-zero on a mismatch: the held-out
χ² of ADDtab (108.0717 vs 108.1), J (34.8649 vs 34.9) and all four SEQ conventions;
T1's TAB mean +7.5292 and 6/6 positive; G6's ADDtab +3.3931, SEQd43 −3.2955 and the
**74 %/74 % bracket — which is G6's own FAILED ≥ 80 % bar** and is re-printed rather than
read in the passing direction.

N before/after every filter: ternary 12 shard files, 21 records → 18 admitted over 6
cells; K-ladder 20 files, 63 records at N = 1500 → 63 admitted, 19 surviving the G6 core
rule. **Calibration 25 cells, held out 15 points.**

The design-time back-of-envelope recorded in the prereg said ~60; computed **56.5**, a
6 % deviation. That is a **flag on the script**, not evidence for the model — the envelope
was made from the same tables.

## G. The weight correlation reproduces, and weakens on the mix residual

`spearman(|residual|, w(u > 0.90))` on ADDtab: pooled (lognormal + uniform, n = 10)
**−0.939**, lognormal **−1.000**, uniform **−0.900** — all three recovered exactly. On the
**MIX** residual it weakens to **−0.636**. The pooled figure overstates independence
(within each family both axes are monotone in CV alone, effective N nearer 2 than 10), and
the within-group split is printed beside it every time. **A weight table says where φ is
sensitive, not where the model fails.**

---

# H. The K = 6 forward content is convention-dominated: UNRESOLVED, not established

`out/diag_k6_interp.txt`, post-hoc, **not a gate**. Written after a physics-referee pass on
`derivations/k_axis_bound.md` reported that the document's 2.99 comes from a **two-point**
log-log interpolation that **discards the third measured point**, and that an equally
defensible form through all three lands below the bar. **Reproduced from the `ug` shards,
not from the referee's report**, with arm B as the control (recorded 2.21, computed
**2.2056**) and the document's own number recovered by its own form (**2.9937** vs a
printed 2.99). Exits non-zero if either control fails.

Measured signals: u = 0.66667 → **7.0470**, u = 0.70711 → **4.4273**, u = 0.80000 →
**2.2056**, all as `max(φ − β)/PHI_SCALE`, **the same units as the 2.57 bar** — not the U1
table's "φ span" column, whose arm B is 2.6, and the referee confirmed that substituting
the wrong column flips K = 7 from no to yes.

| form | pars | anchored at (1, 0) | s(K = 6) | margin | s(K = 7) | margin |
|---|---|---|---|---|---|---|
| log-log, 2 bracketing pts **[the document's]** | 2 | no | 2.994 | 1.165 | 2.306 | 0.897 |
| linear in (u, s), 2 pts | 2 | no | 3.214 | 1.250 | 2.356 | 0.917 |
| quadratic in (u, s), 3 pts | 3 | no | **2.558** | **0.995** | 2.189 | 0.852 |
| quadratic in (ln u, ln s), 3 pts | 3 | no | 2.858 | 1.112 | 2.280 | 0.887 |
| power law in (1−u), 2 pts | 2 | **yes** | 3.128 | 1.217 | 2.334 | 0.908 |
| power law in (1−u), 3 pts LSQ | 3 | **yes** | 3.253 | 1.266 | 2.297 | 0.894 |
| quadratic in (ln(1−u), ln s), 3 pts | 3 | **yes** | 2.758 | 1.073 | 2.248 | 0.875 |

"Anchored" = forced through the **exact** physical point s(u = 1) = 0: at u = 1 the two
classes are identical, the mixture **is** monodisperse, φ = β and the filling signal
vanishes. That constraint is free and is not a fitted parameter.

* **K = 6 (u = 0.75786): UNRESOLVED.** The signal spans **2.558 to 3.253** (a factor
  **1.27**) and **the 2.57 bar is inside that span**; the verdict is "attemptable" on
  **6 of 7** forms. The document's 2.99 is one draw from that range. **"The K axis closes
  at 6" is not established by this estimate** — though the balance leans toward
  attemptable, and the single dissenting form misses by 0.5 %.
* **K = 7 (u = 0.79370): CLOSED, robustly.** Span **2.189 to 2.356**, factor 1.08, bar
  **outside** it, **0 of 7** forms reach it. The recorded 2.31 is safe.
* **No two-parameter form fits the three points.** Fit each to the two outer points and
  predict the middle one: log-log misses by **9.4 %**, linear by **26.0 %**, power-law in
  (1−u) by **18.6 %**. The interpolation is genuinely **under-determined by the data**, so
  no form among them is privileged. That is a stronger statement than "the margin is
  thin", and it is what turns the K = 6 estimate from marginal into unresolved.

**Nothing here moves a bar and nothing here can decide K = 6.** What can is unchanged and
is what the document already asks for: recompute the noise-free identifiability ceiling at
**r = 4^(1/5) = 1.31951**, at several assumed truths, before any pre-registration.

---

## What this block can and cannot decide

**CAN.** Whether one mixing parameter fitted on discrete ladders carries to continuous
PSDs (**no**, by 5.67σ); whether the two discrete ladders want the same amount of nesting
(**yes**, t = 0.69); whether the P3 failure is a form or a transfer failure (**transfer**,
by the oracle); and whether the K = 6 forward content survives its own interpolation
choice (**no**).

**CANNOT.** Which nesting is right — the 15.38× convention spread is the honest width.
Anything about u < 0.25, where `a(u)`, `b(u)` are interpolated toward a proved limit rather
than measured. And it **cannot be quoted as a parameter-count win over J**: J's two
exponents saw **binaries only**, while λ saw the ternary and the K-ladder as well, so the
mix has seen strictly more data.
