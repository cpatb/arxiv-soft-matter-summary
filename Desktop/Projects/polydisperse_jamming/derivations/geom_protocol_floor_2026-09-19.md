# The geometry batch's significances, re-graded against the protocol's own scale

**2026-09-19.** Numbers printed by `geom_protocol_floor_2026-09-19.py` into
`geom_protocol_floor_2026-09-19.txt` (rc 0, ALL CHECKS PASSED). Corrects
`geom_separability_2026-09-19.md` (item A) and `geom_skewness_closure_2026-09-19.md`
(item C2), both of which carry a banner pointing here.

> ## RESULT
>
> A physics-referee pass returned one FATAL: every significance in items A and C2 was
> computed from the **seed-paired standard error alone**, and this repository documents in
> its own printed artifact that that is the wrong bar.
>
> `mixture_model.py`: *"grading scale PHI_SCALE = 2.1e-3: the protocol's own run-to-run
> non-determinism, NOT the seed sem (~2e-04), which is 10× smaller than the reproducibility
> of the thing being modelled."* `test_mixture.py` enforces it for a structurally identical
> claim. **And `\phiScale` = 2.1e-3 is already a macro in the manuscript** — so this is
> grading against a scale the project declares, and then not using it.
>
> | statistic | t, seed s.e. | t, protocol floor | verdict |
> |---|---|---|---|
> | item A additive, lognormal | 27.70 | **1.67** | 🔴 REFUTED → **UNDETERMINED** |
> | item A multiplicative, lognormal | 56.97 | **4.03** | 🔴 REFUTED → **UNDETERMINED** |
> | item A additive, uniform | 0.28 | 0.15 | unchanged (not resolved either way) |
> | item C2 zero-skewness pair | resolved **3 of 5** CV | resolved **1 of 5** | conclusion holds, robustness does not |

---

## 1. The propagation is derived, because it changes verdicts

`PHI_SCALE` is the run-to-run scatter of **one** cell's `φ_J`. A statistic built from `L`
independent cells and averaged over `n` seeds therefore has a floor

```
    se_proto = √L · PHI_SCALE / √n
```

and **`L` counts the cells that do not cancel.** That clause is load-bearing:

- **Item A** combines `(μ,CV)` = (0, 0), (0, 0.5), (0.3, 0), (0.3, 0.5) — four distinct
  cells, none shared. `L = 4`, `se_proto = 0.00242`.
- **Item C2** compares two families' CV-shifts, and uniform and bidisperse share their
  CV = 0 cell **bit-identically**, so that term cancels exactly and the statistic reduces to
  `φ_u(CV) − φ_b(CV)`. `L = 2`, `se_proto = 0.00171`.

**Using `L = 4` for item C2 would give `t = 4.25` at CV = 0.5 and flip it to "not
resolved"**; the correct `L = 2` gives `t = 6.01`. The referee's own figure (1 of 5, with
CV = 0.5 "marginal" at 4.91) used a third propagation. The count agrees; the margin does
not, and the derivation above is why this page's is the one quoted.

The bar used is `max(seed s.e., se_proto)` — the protocol value is a **floor**, not a
replacement. Where the seed scatter is larger it is the honest number, and at item C2's
CV = 0.5 the two coincide to three figures.

## 2. Item A: both forms become UNDETERMINED

This is a **stronger** correction than the referee proposed. It re-graded the additive form
(27.70 → 1.92 on its propagation) and left the multiplicative one resolved at 4.45, because
it computed no separate floor for a **ratio**. A ratio of packing fractions carries a
*relative* floor `PHI_SCALE/φ`, which at `φ ≈ 0.63` is 0.0033 per leg; propagated,
`se_proto = 0.00404` and `t = 4.03 < 4.303`.

**So neither separable form is excluded at three seeds.** What survives is an ordering: the
additive residual (−0.00404) is smaller than the multiplicative one (−0.01628) on every
family. That is a descriptive comparison of two one-parameter forms and needs no
significance test — and it is all "additive to ~6.5 % of the shift" was ever entitled to say.

## 3. What this does NOT touch

- **Items C1, B, and the `z_iso` and Furnas derivations are ALGEBRA.** They carry no
  t-test, cannot fail on data, and the referee re-derived all four independently — including
  SymPy confirmation that `φ = φ_B/(1−ψ)` is exact, that the two-point level set is the root
  pair of `(1+c²)(r−1)²q² + (r−1)[2(1+c²)−r−1]q + c²` with roots 1/10 and 1/2 at
  `r = 3, CV = ½`, and that `z_iso = d(d+1)/(d−x_m)` reproduces both limits.
- **Item B's verdict already leans "not resolved"**, and under-counting noise cannot
  manufacture that, so UNDETERMINED is unaffected.
- 🔴 **Item D's isostatic drift (t = −7.25) uses the same machinery and CANNOT be re-graded
  here**, because there is no `z_iso` analogue of `PHI_SCALE` on record. The honest
  statement is that its significance is **unaudited**, not that it is wrong. Establishing a
  run-to-run scale for `Z_nr` is a cheap, unregistered gap.

## 4. Where else this bar belongs

`paper/extract_numbers.py` and `paper/main.tex` were corrected to it the same day: the
manuscript now states the bar, the propagation and its reason, and prints **both** counts
(`\zeroSkewNResolved` = 1 and `\zeroSkewNResolvedSeedOnly` = 3), saying plainly that the
seed-only reading is not supported by the reproducibility of the quantity compared.

⚠ **The general form of the defect is worth carrying**: a seed-paired standard error is the
scatter of the *sampling*, not of the *protocol*. Wherever a project has measured its own
re-run reproducibility — and this one has, twice, in `PHI_SCALE` and in the conflicting
duplicate records — that number is a floor under every error bar built from seeds.
