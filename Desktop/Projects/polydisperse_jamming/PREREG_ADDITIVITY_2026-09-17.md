# PREREG_ADDITIVITY_2026-09-17 — the additivity repair

**Committed before `grade_additivity.py` was written or run. No new cells.** Every cell
this reads is already on disk and was read months (families, binaries) or days (ternary,
K-ladder) ago.

---

## 0. HONESTY DECLARATION, first, because it changes how everything below reads

**This is NOT a blind prediction.** Both datasets are on disk and their residual tables
are printed in `out/transfer_kladder.txt`, `out/grade_ternary.txt` and
`out/grade_kladder.txt`, and I read all three while designing this. What is genuinely
registered here is the **rule** — which cells calibrate, which are held out, which
functional form, which convention is primary, and what counts as a win — fixed before the
number exists. That is a weaker claim than blindness and is labelled as such, exactly as
`mixture_closure.py`'s own docstring labels its transfer.

**And I did a back-of-envelope while designing.** Mixing the published mean residuals at
λ ≈ 0.52 gives a held-out χ² of **roughly 60** — i.e. it suggested P2 passes and P3
fails. That estimate is written down *here* so that (i) the result cannot later be read
as confirming a blind prediction, and (ii) a large disagreement between 60 and the
computed value is itself a flag that the script is computing something other than what
this document describes.

---

## 1. What is being repaired, and why this is the useful half of a withdrawn claim

`transfer_kladder.py` (G8, declared post-hoc) asked whether replacing the additive pair
sum with a fully nested sequential rule closes the held-out transfer gap. It does not:
the answer moved by **11.90×** across the four undeclared composite-diameter conventions,
which is more than the effect, so G8 stayed unsettled. The stronger claim drawn from that
table — "the continuous-family transfer needs `a(u)` at `u` > 0.9" — was **WITHDRAWN**
by a referee, because `spearman(|family residual|, w(u > 0.90)) = −0.939`: the cells
leaning hardest on that region are already fit to within noise.

What survives is a **bracket**, and it is the same bracket in all three datasets:

| set | n | ADDtab (additive) | SEQd43 (fully nested) |
|---|---|---|---|
| ternary, T1/T4 | 6 | under-fills **6/6**, mean **+7.53** | over-fills **6/6**, mean **−6.47** |
| K-ladder, G6 (N = 1500 multi-class) | 19 | under on **74 %**, mean **+3.39** | over on **74 %**, mean **−3.30** |
| continuous families, held out | 15 | under on 11/15, mean **+1.69** | over on 12/15, mean **−2.89** |

Sign convention throughout **this** document and `grade_additivity.py`:
**residual = MEASURED − MODEL**, the `grade_ternary` / `grade_kladder` convention.
Positive = the model under-predicts φ. `transfer_kladder.py` uses the opposite sign
(`model − φ`); the script **asserts** the two agree after a flip, because this project has
already shipped one artifact that contradicted its own table over exactly this.

⚠ **The bracket is 6/6 on the ternary and only 74 %/74 % on the K-ladder — G6's own
registered bar was ≥ 80 % and G6 FAILED it.** That failure is a known, recorded weakness
of the premise this model is built on, and it is not being read in the passing direction.
It is re-printed by the grader beside the result.

---

## 2. The model — ONE parameter

For every cell, with `A` = the additive tabulated prediction `packing_density_tabulated`
and `S_c` = `kladder.sequential_k` at composite-diameter convention `c`:

```
    phi_mix(lambda, c) = (1 - lambda) * A + lambda * S_c
```

`lambda` is the **only** new free parameter. The interaction functions `a(u)`, `b(u)`
inside both `A` and `S_c` are the shipped tabulated closure, calibrated on **binaries
only** and **not refitted here** — asserted unmoved by the grader.

Fitted by **unweighted least squares in units of `PHI_SCALE` over the calibration set**,
which for a linear model has the closed form

```
    lambda_hat = sum (m - A)(S - A) / sum (S - A)^2
```

The grader computes it **both** in closed form and by `scipy.optimize.minimize`, and
asserts they agree to 1e-9. A Δχ² = 1 profile interval on `lambda` is reported.

**Primary convention: `c = d43`**, because that is the arm T4 registered as the bracket's
upper edge ("SEQd43 over-fills 6/6"). **All four conventions (`d10`, `d21`, `d32`, `d43`)
are fitted and reported**, and the spread of the held-out χ² across them is printed beside
the primary answer — the undeclared composite-diameter choice has twice been worth more
than the effect it was being used to measure (11.90× here, 28.95× in U6).

### What λ CANNOT do, stated in advance

At **K = 2** the nested rule reduces **exactly** to the additive binary closure (proved in
`sequential_k`'s docstring, wiring-gated in `predict_kladder.py`). So `S - A ≡ 0` there
and **λ has literally zero leverage on every K = 2 cell**. The grader **counts and prints**
the leverage split in both sets. Expected: 5 of 19 K-ladder cells and **5 of 15 held-out
points (the entire bidisperse family)** are λ-independent, so the held-out test has
**10 effective points, not 15**, and 5/15 of the held-out χ² is a constant floor equal to
ADDtab's contribution there. Any statement about the held-out χ² that ignores this floor
is wrong.

---

## 3. The sets — fixed here, not at read time

**CALIBRATION (λ is fitted to these and to nothing else):**

* **Ternary**, all 6 compositions of `phase12_ternary.COMPOSITIONS`, N = 1500, μ_p = 0.3,
  through `grade_ternary.load()` and `grade_ternary.admissible()` verbatim.
* **K-ladder**, `grade_kladder`'s own G6 core at N = 1500: `arm ∈ {K, SPLIT, NA2}` and
  `K > 1`, i.e. 19 cells. The monodisperse T0 control (K = 1) is excluded because it is
  β by construction and carries no interaction at all.
* **N = 3000 is excluded from the primary fit** and reported as a replicate only.
  Reason: the families are measured at N = 1500 and the K-ladder measured its own
  block-wide N offset at +1.51 scale, so pooling two N would put a known offset into λ.

Expected total: **25 calibration cells, 20 of them carrying leverage.**

**HELD OUT (no point enters the fit):** the **15** continuous-family points —
`{lognormal, uniform, bidisperse} × CV {0.1, 0.2, 0.3, 0.4, 0.5}` — read through
`mixture_model.measured(mu_p=0.3)` and `mixture_model.psd_classes`, at
`mixture_model.K_CLASSES = 160` equal-count classes, i.e. the identical objects
`transfer_kladder.py` scores. The grader **reproduces ADDtab = 108.1 and J = 34.9 on
them before doing anything else and exits non-zero on a mismatch.**

---

## 4. The competitors, and the bars

Standing held-out χ² over the 15 points, quoted **verbatim** from
`out/transfer_kladder.txt` and **re-derived by the grader**:

| model | held-out χ² | free parameters beyond the binary calibration |
|---|---|---|
| **ADDtab** | **108.1** | 0 |
| **J** (two-exponent parametric) | **34.9** | 0 (its 2 exponents are binary-fitted) |
| SEQd10 / d21 / d32 / d43 | 186.3 / 42.8 / 48.0 / 509.3 | 0, but the convention is a choice |
| **MIX** (this document) | to be computed | **1** (λ) |

⚠ **The parameter count is not the whole story and the grader says so.** J's two
exponents saw binaries only; MIX's λ saw the ternary and the K-ladder as well. MIX has
therefore seen strictly more data than J, and beating J is correspondingly less impressive
than the "1 vs 0" column makes it look. Recorded here so it cannot be quoted the other way.

### Registered predictions — each one able to FAIL

* **P0 (the physics content, and the one that matters most).** λ̂ fitted on the **ternary
  alone** and λ̂ fitted on the **K-ladder alone** agree: each point estimate lies inside
  the other's Δχ² = 1 interval, **and** the two intervals intersect. *Fails if the two
  discrete-ladder families want different amounts of nesting, in which case λ is a
  property of the dataset and not of the combination rule, and the transfer below means
  nothing whatever it scores.*
* **P1 (the bracket premise).** On the primary convention d43, λ̂ ∈ (0, 1) strictly, and
  its Δχ² = 1 interval is contained in [0, 1]. *Fails if the truth is not between the two
  rules — λ is fitted UNBOUNDED precisely so this can fail.*
* **P2 (the weak win).** Held-out χ² at λ̂ < **108.1**. *Nearly free — SEQd21 already
  manages 42.8 with no fitting at all — and it is registered as weak for that reason.*
* **P2b (the DISCRIMINATING win).** Held-out χ² at λ̂ < the held-out χ² of the **single
  SEQ convention that wins on the CALIBRATION set**. That competitor is selected without
  touching one held-out point, costs no free parameter, and is the honest thing a reader
  would do instead. *This is the bar the exercise can most plausibly fail, and on my
  back-of-envelope it probably does.*
* **P3 (the headline win).** Held-out χ² at λ̂ < **34.9**, i.e. the one-parameter mix of
  two measured-interaction rules beats the smoothed parametric closure. *My envelope says
  no (~60).*
* **P4 (convention robustness).** P2 holds on **≥ 3 of the 4** conventions. *Fails if the
  win is a property of the composite-diameter choice rather than of the mixing.*

**None of these bars moves.** If one fails it is reported as a failure and the model's
range is what the failure says it is. No refit, no reweighting, no dropped cell.

### Registered NULLS — outcomes that are real results, not non-results

* **λ̂ ≈ 0** (interval containing 0): the discrete ladders want **no** nesting admixture,
  the ternary/K-ladder under-fill is *not* explained by the combination rule, and the
  additivity story is dead. Publishable.
* **λ̂ interior but P2b fails**: the bracket is real, one λ describes both discrete
  families, and it still does not buy anything a plain convention choice does not — i.e.
  **the non-additivity correction and the composite-diameter convention are degenerate on
  this data**, which is a statement about the FORM. Publishable, and the most likely
  outcome.
* **P0 fails**: λ is not transferable and nothing downstream of it is quotable.

---

## 5. What the grader must print (rule, not decoration)

1. **Reproduction first, exit non-zero on mismatch**: ADDtab 108.1, J 34.9, SEQ×4 on the
   held-out 15; T1's TAB mean +7.53 and 6/6; G6's ADDtab mean +3.39 and the 74 %/74 %
   bracket. *A diagnostic that cannot reproduce the number it audits has not audited it.*
2. **N before and after every filter**, on both sets.
3. The **leverage split** (§2) in both sets, and the λ-independent χ² floor.
4. The **within-group split beside every pooled figure** — the pooled
   `spearman(|residual|, w(u>0.90)) = −0.939` overstates independence because within each
   family both axes are monotone in CV alone, so the effective N is nearer 2 than 10.
   Report **lognormal −1.000 / uniform −0.900** beside it, and recompute all three.
5. **The nested-convention spread beside the primary answer**, always.
6. Per-family held-out χ², so a win concentrated in one family cannot read as a general one.
7. The **closed-form vs optimiser** agreement, and the seed-se-weighted λ̂ as a robustness
   figure that is **not** the primary.
8. `FAILURES: [...]` as the last line, empty only if every P above passed.

No number in `CLAIMS_*` or the Status card may be typed by hand; all are printed by this
script.

---

## 6. Scope guard

`grade_additivity.py` writes only `out/grade_additivity.txt`. It re-decides **no** closed
gate: it imports `grade_ternary` and `grade_kladder` **as loaders only** and calls neither
`main()`. `test_additivity.py` asserts that structurally, on the import lines.
