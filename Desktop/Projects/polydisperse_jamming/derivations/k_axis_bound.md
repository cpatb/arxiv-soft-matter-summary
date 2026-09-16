# The K-axis bound for calibrated size ladders

**Written 2026-09-17.** Checks: `derivations/k_axis_bound_check.py`, transcript
`derivations/k_axis_bound_check.txt`. Every number below is printed by that script
except where it is quoted from a named artifact and said to be.

> **RESULT.** A ladder of `K` size classes on which *every* pair ratio lies in a
> calibration window `[u_lo, u_hi] ⊂ (0,1)` exists **if and only if**
> `u_lo ≤ u_hi^(K−1)`, i.e.
>
> **K − 1 ≤ κ ≡ ln(u_lo) / ln(u_hi)**,   so   **K_max = 1 + ⌊κ⌋**,
>
> for **any** ladder, geometric or not. The span is then confined to
> `u_hi^−(K−1) ≤ d_max/d_min ≤ 1/u_lo`, and **when κ is an integer and K = K_max the
> ladder is unique**: the geometric one of ratio `u_hi` and span exactly `1/u_lo`.

---

## 1. Frame

**Object.** A ladder is a finite set of diameters `d_1 < d_2 < … < d_K`, `K ≥ 2`.
The *pair ratio* of classes `i < j` is `u_ij = d_i/d_j ∈ (0,1)` (small over large —
this is the project's convention, `mixture_closure.tabulated_closure` keys its knots on
`1/ratio`; §5 checks the bound is the same under the other convention). Write the
*adjacent ratios* `v_k = d_k/d_{k+1}`, `k = 1 … K−1`.

**Calibration window.** `[u_lo, u_hi] ⊂ (0,1)` is the interval of `u` over which the
pairwise-filling closure carries *measured* knots. A ladder is **calibrated** when every
one of its `K(K−1)/2` pair ratios lies in that window — i.e. the closure is asked only
to interpolate, never to extrapolate, on any pair of that ladder.

**Assumptions**, all structural and none physical:
1. The closure is *pairwise* — it reads the PSD through its pair ratios. This is the
   form established in `mixture_model.packing_density` / `mixture_closure`; the bound
   says nothing about a closure of another form.
2. Classes are distinct, so every `v_k ∈ (0,1)` strictly.
3. Interpolation between measured knots is admissible and extrapolation beyond them is
   not. This is the whole content of "calibrated" and is a *convention*, not a theorem.
   §6 reports what the bound gives under two different, defensible readings of it.

**Dimensionless groups.** There is only one: `u`, a ratio of lengths. `K` is a count.
Nothing here has a dimension, and §5 checks that explicitly.

---

## 2. Derivation

**(1) Telescoping.** For `i < j`, by cancelling the intermediate diameters,

    u_ij = d_i/d_j = (d_i/d_{i+1})(d_{i+1}/d_{i+2}) ⋯ (d_{j-1}/d_j) = Π_{k=i}^{j−1} v_k .

In particular the extreme ratio is the product of **all** adjacent ratios:

    u_1K = Π_{k=1}^{K−1} v_k .                                                      (1)

**(2) The extremes are adjacent and extreme.** Every `v_k ∈ (0,1)`, so multiplying a
product by one more such factor strictly decreases it. `u_ij` is a product over the
index set `{i, …, j−1}`; the smallest such sets are the singletons (`j = i+1`) and the
largest is the whole set (`i=1, j=K`). Hence

    max_{i<j} u_ij = max_k v_k ,      min_{i<j} u_ij = u_1K .                       (2)

**(3) The constraint.** "Every pair ratio in `[u_lo, u_hi]`" is therefore *exactly* the
pair of conditions

    max_k v_k ≤ u_hi          (no pair too close in size)                           (3a)
    u_1K      ≥ u_lo          (no pair too far apart)                               (3b)

— nothing is lost by dropping the other `K(K−1)/2 − K` pairs, because (2) says they are
bracketed by these two.

**(4) Necessity.** Combining (1), (3a) and (3b):

    u_lo  ≤  u_1K  =  Π_{k=1}^{K−1} v_k  ≤  u_hi^(K−1) ,                            (4)

the middle step being (1) and the last being (3a) applied `K−1` times. Taking logs and
dividing by `ln(u_hi) < 0` (which **reverses** the inequality):

    K − 1  ≤  ln(u_lo) / ln(u_hi)  ≡  κ .                                           (5)

**(5) Sufficiency — the bound is attained, not merely implied.** Suppose
`u_lo ≤ u_hi^(K−1)`. Take the geometric ladder `v_k = u_lo^{1/(K−1)}` for every `k`.
Then `u_1K = u_lo`, satisfying (3b) with equality, and
`max_k v_k = u_lo^{1/(K−1)} ≤ u_hi` exactly when `u_lo ≤ u_hi^(K−1)`, which is the
hypothesis. Every intermediate pair has `u_ij = u_lo^{(j−i)/(K−1)} ∈ [u_lo, u_lo^{1/(K−1)}]
⊆ [u_lo, u_hi]`. So a calibrated `K`-class ladder exists. Hence

    **K_max = 1 + ⌊κ⌋ ,   κ = ln(u_lo)/ln(u_hi) .**                                 (6)

**(6) The span is confined, and at the top rung it is pinned.** With
`S = d_K/d_1 = 1/u_1K`, (3b) gives `S ≤ 1/u_lo` and (1)+(3a) give `S ≥ u_hi^−(K−1)`:

    u_hi^−(K−1)  ≤  S  ≤  1/u_lo ,                                                  (7)

an interval that is non-empty exactly under (5). **If κ is an integer and `K = 1 + κ`,
the two ends of (7) coincide** at `S = 1/u_lo = u_hi^−(K−1)`, so the span is *forced*.
More: equality in (4) forces `Π v_k = u_hi^(K−1)` with every `v_k ≤ u_hi`, which can
only hold if `v_k = u_hi` for every `k`. **The top-rung ladder is unique and geometric.**

---

## 3. Checks

All four ran; the command is `python3 derivations/k_axis_bound_check.py`, transcript in
`derivations/k_axis_bound_check.txt`, `FAILURES: []`, exit 0.

| # | Check | Result |
|---|---|---|
| 3 | **SymPy.** Step (1) verified symbolically on all 15 pairs at `K = 6` from a ladder built from the *definition* `v_k = d_k/d_{k+1}`; the telescoping to `u_1K = Π v_k` verified separately; the inversion `u_lo = u_hi^n ⟹ n = ln u_lo/ln u_hi` verified by `sympy.solve`. | PASS |
| 3b | **Step (2), independently.** `max(pairs) == max(adjacents)` and `min(pairs) == d_1/d_K`, required as **exact float equality**, on 2000 random ladders of `K = 2…8`. 0 violations. | PASS |
| 4 | **Limits.** `u_hi → u_lo⁺` (window a single point) → `κ → 1`, `K_max = 2`: only a binary at that exact ratio fits. `u_hi → 1⁻` → `κ → ∞`: unbounded, as it must be, since a ladder can be subdivided arbitrarily finely with no upper constraint. `u_lo → 0⁺` → `κ → ∞`: unbounded, since with no lower constraint an arbitrarily long ladder fits inside one adjacent step. `K = 2` on `[0.25, 0.70711]`: `K_max = 5 ≥ 2`, any `u` in the window works. | 4/4 PASS |
| 5 | **Dimensions / invariance.** Rescaling every diameter by `10^±6` leaves every pair ratio bit-unchanged — the bound is a statement about ratios only, and has no dimension to get wrong. Under the *other* convention (`U = large/small`, window `[1/u_hi, 1/u_lo]`) `κ` is unchanged to 10 decimals: `4.0000000000` both ways. | PASS |
| 6 | **Numeric, brute force.** Feasibility decided **without the formula** — 4000 random ladders (Dirichlet log-gaps × random span) plus a 4000-point deterministic sweep of geometric ladders — over a 4 × 4 × 6 grid of (`u_lo`, `u_hi`, `K`). **96 of 96 cells agree** with (6). | PASS |

**Anchor outside the derivation.** The shipped `kladder.ladder(K)`, which the derivation
never saw, at span 4:

| K | span | min pair | max pair | in `[0.25, 0.70711]`? |
|---|---|---|---|---|
| 2 | 4.0000 | 0.2500000000 | 0.2500000000 | YES |
| 3 | 4.0000 | 0.2500000000 | 0.5000000000 | YES |
| 4 | 4.0000 | 0.2500000000 | 0.6299605249 | YES |
| 5 | 4.0000 | 0.2500000000 | **0.7071067812** | YES |
| 6 | 4.0000 | 0.2500000000 | 0.7578582833 | **NO** |

The bound says `K_max = 5` on that window, and the shipped ladder crosses out of it at
exactly `K = 6`. The `K = 6` max pair is `0.7578582833 = 0.25^{1/5}` to ten figures,
which is (7)'s lower end — i.e. the anchor reproduces the *mechanism* of the bound, not
only its verdict.

**Prefactors.** There are none. Every constant in this document is either a window end
quoted from `grade_mixture_ext.tabulated_closure`'s knot set, or `⌊·⌋`.

---

## 4. Is the zero slack a coincidence, an enforced constraint, or a theorem?

The shipped `K = 5` span-4 ladder sits on **both** window edges with zero slack
(`min pair 0.2500000000`, `max pair 0.7071067812`), and `κ = 4.0000` exactly. Three
separate things are going on and they were worth separating before the tightness was
remarked on.

**(a) The span being pinned is a THEOREM, not luck.** §2 step (6): whenever `κ` is an
integer and `K = K_max`, (7) collapses to a point and the ladder is unique. So *given*
`κ = 4`, a `K = 5` calibrated ladder **has** to be the geometric span-4 one, with zero
slack at both ends. Nothing could have been otherwise. `kladder.py` does not enforce
this — it chose span 4 for a stated and different reason ("every pair has `u ≥ 1/4` by
construction") — and arrived at the only ladder available.

**(b) `κ = 4.0000` exactly is ENFORCED by a chain of two design choices, and is not a
coincidence.** The chain, each link documented in the file that made it:

1. `u_lo = 1/4` is *measured*: `r = 4` is the lowest binary ratio that passes
   `tabulated_closure`'s identifiability rule (`r = 5` fails, `out/diag_uident.txt`
   §1), so `1/4` is the closure's lowest knot.
2. `kladder.py` then set `SPAN = 4 = 1/u_lo` deliberately, so that the extreme pair of
   every rung sits exactly at the lowest knot.
3. `ugap.py` then set `U_STAR = SPAN^{−1/4} = 4^{−1/4}`, the adjacent ratio of the
   `K = 5` rung of *that* ladder.

Composing 2 and 3, `u_hi = u_lo^{1/4}` **identically**, for any span — so
`κ = ln(u_lo)/ln(u_lo^{1/4}) = 4` is an identity of the design, not an arithmetical
accident.

**(c) Therefore the bound is NOT independent evidence that `K = 5` is special.** It
says the calibration window was extended to exactly the point that licenses the shipped
ladder and not one step further. That is a weaker claim than "the K axis is closed at 5
by a theorem", and it is the one this document supports. What the bound *does* supply
that is new is forward content: §6.

**A consistency check that the bound passes and that matters.** *Before* the U-GAP block
the window was `[0.25, 0.66667]`, giving `κ = 3.41902` and `K_max = 4`. So the `K = 5`
rung of `PREREG_KLADDER_2026-09-15.md` was **already outside the calibrated window when
it ran** — which is exactly the defect gate G1c raised independently, from the knot
table rather than from any bound. The bound reproduces a known, separately-found result
it was not built to reproduce.

---

## 5. The window's lower end: two readings, two bounds, both reported

The brief asked whether the honest lower end is `1/4` or `1/3`, noting that
`BENCH_HYPOTHESES.md` row C says the closure "is calibrated at `u ≥ 1/3` and is not a
predictor below it". **The two are about different objects and both are right**, so both
bounds are reported and each is labelled with the claim it governs.

| Window | Referent | κ | K_max | span at K_max |
|---|---|---|---|---|
| `[0.25, 0.66667]` | ADDtab's knots **before** the U-GAP block | 3.41902 | **4** | [3.375, 4.000] |
| `[0.25, 0.70711]` | **ADDtab's knots now** — the closure the K-ladder actually reads | **4.00000** | **5** | **[4.000, 4.000] — FORCED** |
| `[0.33333, 0.70711]` | the *parametric* closure's validated range, row C | 3.16993 | **4** | [2.828, 3.000] |

**Which applies where.**

* The K-ladder's residual `R = (φ − ADDtab)/scale` reads the **tabulated** closure, whose
  `u = 0.25` entry is a *fitted, identifiable knot* (`r = 4`, 6 points, `w(a) = 0.045`,
  `w(b) = 0.025`, `out/diag_uident.txt` §1), not an extrapolation. For that object
  `u_lo = 0.25` and **K ≤ 5**.
* Row C's failure — "blind at ratios 4 and 5 it under-predicts 9 of 9, rms 7.50 against
  a bar of 4.0", X1 FAIL — is a failure of the **parametric** closures `J` / `EXT`, which
  *extrapolate* to `u = 0.25` from calibration at ratios 2 and 3. For any claim resting
  on that form, `u_lo = 1/3` and **K ≤ 4** — so the shipped `K = 5` rung is outside its
  honest window.

**The consequence, stated plainly rather than resolved in the convenient direction.**
`kladder.py`'s docstring justifies the design with "*every* pair on *every* rung has
`u ≥ 1/4` **by construction** … so additivity stays unconfounded with extrapolation on
the whole ladder". That sentence is true of ADDtab and **false of the parametric
closure**, and the file does not say which it means. Any K-ladder statement graded
against `J` or `EXT` at `K = 5` is outside the calibrated window; every statement graded
against ADDtab is inside it. The graded gates G1, G1c, G7 and U3 all use ADDtab, so they
stand — but the general sentence in the docstring is too strong and is flagged in
§7 as the one thing this derivation asks to be changed.

---

## 6. What the bound is actually good for: the forward content

The bound turns "can we push K past 5?" from a budget question into an arithmetic one
with a single number attached. To license `K` on a span-4 ladder the window must reach

    u_hi ≥ u_lo^{1/(K−1)} = 0.25^{1/(K−1)} .

The U-GAP block measured the *filling signal* `max(φ − β)/scale` at three ratios, and
measured the resolving power of that design (**2.57 scale**, quoted from
`out/grade_ugap.txt`'s arm-B ceiling table, which fails at 2.55). Log-interpolating the
measured signal gives, per rung:

| K | `u_hi` needed | signal est. | margin | attemptable? |
|---|---|---|---|---|
| 5 | 0.70711 | 4.43 | 1.72 | already measured |
| **6** | **0.75786** | **2.99** | **1.16** | **yes** |
| 7 | 0.79370 | 2.31 | 0.90 | no |
| 8, 10, 12 | 0.820, 0.857, 0.882 | < 2.21 | < 0.86 | no — below arm B's measured 2.21, which FAILED |

The signals are computed in the check script from the `ug` shards, and its **control is
exact**: arm B comes out at 2.21, the number `out/grade_ugap.txt` records. Rows above
`u = 0.80` are *not* extrapolated — the signal is monotone decreasing in `u` and
vanishes at `u = 1`, so they are bounded above by arm B's measured 2.21, which already
failed. That is an argument from monotonicity, not a fit past its range.

🔴 **This qualifies the Status card's §3 item 2** ("do not push K past 5 on a span-4
ladder"). Its two exhibits — adjacent ratios 0.82 at `K = 8` and 0.88 at `K = 12` — are
correctly excluded, but *by a measurement* (arm B's 2.21 failing at `u = 0.80`) rather
than by the argument given, and **`K = 6` was never evaluated**. It needs a knot at
`u = 0.75786`, whose estimated signal clears the bar by 1.16×.

**Status of that number: an ESTIMATE from a three-point log interpolation, not a
measurement, and the margin is thin.** It says `K = 6` is worth attempting where
`a(0.8)` was not. It does not say it will succeed.

---

## 7. Open questions and what this asks to be changed

1. **`kladder.py`'s docstring overstates its own constraint.** "`u ≥ 1/4` by
   construction … unconfounded with extrapolation" is true against ADDtab and false
   against the parametric closure. It should name which. (§5)
2. **The `K = 6` estimate is not a measurement.** It rests on a three-point log
   interpolation of the filling signal and a resolving-power figure taken from arm B's
   ceiling table at arm B's ratio, not at `u = 0.7579`. A proper design would recompute
   the ceiling at `r = 4^{1/5} = 1.31951` before any cell. (§6)
3. **The bound assumes a strictly pairwise closure.** The ternary block (BENCH row D)
   found the functional is **not** additive beyond two classes, so a repaired closure may
   read triples — for which "every *pair* inside the window" is the wrong feasibility
   condition and this whole bound would need restating. That is the one assumption that
   is under active attack elsewhere in the project.
4. The bound is a statement about *what can be asked of the closure*, never about
   physics. Nothing here says a `K = 6` pack behaves differently from a `K = 5` one.
