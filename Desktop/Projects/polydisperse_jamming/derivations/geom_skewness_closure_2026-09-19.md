# Item C, part 2 — Is the family effect a SKEWNESS effect? An external test

**VERDICT: PARTIALLY SUPPORTED, AND REFUTED AS A COMPLETE DESCRIPTION.** Skewness at
fixed width is the right leading variable — an external two-parameter closure, fitted in
2013 to 81 other packings, predicts this manuscript's headline gap to within 2 se. But
(δ, S) is **not sufficient**: two S = 0 families at identical δ are measured apart at
t = +6.02.

**What this anchor can and cannot decide.** It is external and could have failed: the
coefficients are someone else's and the matched pair is not in their data set. It CANNOT
decide the friction question (Desmond & Weeks are frictionless; the pair was measured at
μ_p = 0.3), it CANNOT decide the *level* (only increments), and at S = 2.667 it is a
6.7 % extrapolation past their own posted range.

Numbers printed by `geom_skewness_closure_2026-09-19.py` →
`geom_skewness_closure_2026-09-19.txt`.

---

## 0. The closure, and a typo that matters by a factor of ten

K. W. Desmond and E. R. Weeks, *Influence of particle size distribution on random close
packing of spheres*, **Phys. Rev. E 90, 022204 (2014)**,
doi:10.1103/PhysRevE.90.022204 (preprint **arXiv:1303.4627**). Equation (4):

> **φ_RCP = φ*_RCP + c₁ δ + c₂ S δ²**,  φ*_RCP = 0.634, c₁ = 0.0658, **c₂ = 0.0857**

with, their Eqs. (1)–(2), on the **number-weighted** radius distribution `P(R)`:

    δ = √⟨ΔR²⟩ / ⟨R⟩            ΔR = R − ⟨R⟩
    S = ⟨ΔR³⟩ / ⟨ΔR²⟩^{3/2}

Both δ and S are scale-invariant ratios, so radius versus diameter is immaterial; the
weighting is not — they are number moments, which this project's families are.

### ⚠ c₂ = 0.0857 or 0.857? Settled here, not by a lookup

**Two independent literature passes on 2026-09-19 returned different values of c₂ from
the same paper** — 0.0857 (reported as the published PRE value) and 0.857 (which is what
arXiv:1303.4627v2 prints, in both the Fig. 2 caption and the body). A factor of ten in
the skewness term decides everything below, so it is settled by arithmetic on quantities
the paper states *in words*, in §0 of the script:

| check | c₂ = 0.0857 | c₂ = 0.857 |
|---|---|---|
| (a) D&W's own worked example: *"changing S to ±1 shifts φ_RCP by ≈±0.02"* at δ = 0.40 → `c₂·1·0.16` | **0.01371** ✓ | 0.13712 — 7× their stated shift |
| (b) admissibility at their own largest posted (δ, S) = (0.413, 2.50); FCC limit is π/(3√2) = 0.740480 | **0.69772**, admissible | **1.02662 — above FCC, impossible** |
| (c) predicted matched-binary gap (§5) against a measured +0.0635 ± 0.0032 | **+0.05713** | +0.57133 |

All three agree on **c₂ = 0.0857**, and (b) is decisive on its own: the preprint value
makes the closure return a packing fraction above the densest possible sphere packing
*inside its own fitted range*. The δ term is consistent with either
(`c₁ × 0.40 = 0.02632` against the same "≈0.02").

Supporting evidence from the lookups, not relied on: 0.0857 reproduces D&W's stated
residual on their posted 81-packing data set (mean |res| 0.00105 against their stated
0.001; 0.857 gives 0.0332, 33× worse), and it is the value Weeks uses when replotting the
same data in arXiv:2404.02316 Eq. (24). **Cite and code from the published version.**

No uncertainties are published for any of the three coefficients. What is published is
residual accuracy: ±0.002 for "nearly all" points (88 % on their posted data), mean
absolute 0.001.

**Fitted box**, taken from their posted data because the paper prints none:
δ ∈ [0.049, 0.413]; |S| ≤ 2.50, but |S| above 1.17 comes **only from their binary
packings**. Kurtosis: "K seems not to be a useful fit parameter for φ_RCP".

**Frictionless?** The word *friction* does not appear in their paper. It is an inference
from their cited algorithm (Xu, Blawzdziewicz & O'Hern, PRE **71**, 061306 (2005),
*"…frictionless disks"*). Treat as unverified-in-source, near-certain from provenance. The
like-for-like arm here is therefore **μ_p = 0**, which this project has for lognormal and
uniform over the full CV ladder.

---

## 1. The three families in D&W coordinates

Closed forms, each verified against a 10⁶-sample draw (worst |diff| 0.0094):

| family | δ | S |
|---|---|---|
| lognormal | CV | **3·CV + CV³** |
| uniform | CV | **0** (symmetric on an interval) |
| bidisperse, 50:50 by number | CV | **0** (symmetric two-point) |

At CV = 0.5 the lognormal reaches S = 1.625 — inside D&W's binary range, above their
continuous one. CV = 0.5 is δ = 0.5, which is 21 % past their largest δ; flagged, not
hidden.

---

## 2. The frictionless arm — increments, not levels

The absolute levels are not comparable (different protocol, and 0.634 is D&W's own fitted
intercept). What is comparable is the **CV increment**, because Eq. (4) is additive in
that intercept — which is precisely the structure item A measured. So compare
`D_meas(CV) = φ_J(CV) − φ_J(0)` against `D_DW(CV) = c₁δ + c₂Sδ²`.

**lognormal, μ_p = 0, n = 3 seeds:**

| CV | S | D_DW | D_meas | se | diff | t | ratio |
|---|---|---|---|---|---|---|---|
| 0.1 | +0.301 | +0.00684 | +0.00512 | 0.00025 | −0.00171 | −6.91 | 0.749 |
| 0.2 | +0.608 | +0.01524 | +0.01302 | 0.00044 | −0.00222 | −5.05 | 0.854 |
| 0.3 | +0.927 | +0.02689 | +0.02595 | 0.00122 | −0.00094 | −0.77 | 0.965 |
| 0.4 | +1.264 | +0.04365 | +0.04025 | 0.00083 | −0.00340 | −4.10 | 0.922 |
| 0.5 | +1.625 | +0.06772 | +0.05803 | 0.00268 | −0.00969 | −3.62 | 0.857 |

**uniform, μ_p = 0, n = 3 seeds:**

| CV | S | D_DW | D_meas | se | diff | t | ratio |
|---|---|---|---|---|---|---|---|
| 0.1 | 0 | +0.00658 | +0.00361 | 0.00015 | −0.00297 | −19.99 | 0.549 |
| 0.2 | 0 | +0.01316 | +0.01184 | 0.00035 | −0.00132 | −3.81 | 0.900 |
| 0.3 | 0 | +0.01974 | +0.02002 | 0.00054 | +0.00028 | +0.53 | 1.014 |
| 0.4 | 0 | +0.02632 | +0.02347 | 0.00035 | −0.00285 | −8.13 | 0.892 |
| 0.5 | 0 | +0.03290 | +0.02385 | 0.00162 | −0.00905 | −5.58 | 0.725 |

D&W **over-predicts the absolute increment** in 9 of 10 cells, by 10–45 %, resolved in 5
of 10. That is a protocol-level offset and it is expected: different algorithm, different
finite-size treatment, δ = 0.5 outside their box.

---

## 3. The discriminating statistic — and it works

A closure that returns the same number for both families would have agreed with neither
family's shape, so the informative statistic is the **lognormal/uniform ratio**, which in
D&W is set *entirely* by the `c₂Sδ²` term since uniform has S = 0 exactly:

| CV | D&W ratio | measured ratio | se | t | verdict |
|---|---|---|---|---|---|
| 0.3 | 1.362 | 1.296 | 0.070 | −0.94 | agree |
| 0.4 | 1.659 | 1.715 | 0.044 | +1.30 | agree |
| 0.5 | 2.058 | 2.433 | 0.200 | +1.88 | agree |

**A single published coefficient, c₂ = 0.0857, fitted to other people's packings,
reproduces this project's lognormal-versus-uniform contrast at three CVs.** That is the
strongest external corroboration in this session: the family effect *is*, to leading
order, a skewness effect.

---

## 4. The refuting half — uniform against bidisperse

Uniform and bidisperse both have **S = 0 exactly** and the same δ, so Eq. (4) predicts
them **identical at every CV**. Only μ_p = 0.3 has both. They share their CV = 0 cell
(`paper/main.tex:149-151`), so the reference cancels exactly in the difference.
Seed-paired, seeds {1,2,3}, t_crit = 4.303:

| CV | D_unif | D_bidi | difference | se | t | verdict |
|---|---|---|---|---|---|---|
| 0.1 | +0.00382 | +0.00140 | +0.00242 | 0.00049 | **+4.90** | RESOLVED → D&W refuted |
| 0.2 | +0.01151 | +0.01204 | −0.00053 | 0.00239 | −0.22 | not resolved |
| 0.3 | +0.02077 | +0.02148 | −0.00070 | 0.00080 | −0.87 | not resolved |
| 0.4 | +0.02230 | +0.01896 | +0.00335 | 0.00042 | **+7.94** | RESOLVED → D&W refuted |
| 0.5 | +0.02420 | +0.01389 | +0.01031 | 0.00171 | **+6.02** | RESOLVED → D&W refuted |

**Resolved apart at 3 of 5 CVs.** So (δ, S) is not a complete description of a packing's
φ_J: two distributions matched in both are measured 0.0103 apart at CV = 0.5 — 43 % of the
uniform family's whole CV effect.

⚠ **This is a test `paper/main.tex:222` says has not been done.** The manuscript writes:
*"the only test that could exclude it is the two zero-skewness families against each other
(\dPhiJUnif against \dPhiJBidi), and that pair has not been tested at the calibrated
threshold."* Run seed-paired at the design's own 2 dof, it separates them at t = +6.02.
The manuscript's caveat about a *calibrated* threshold applies to a max-versus-min pair
SELECTED from the data; this pair is named in advance by the hypothesis, so a plain t-test
is the right instrument. **Recommend the manuscript be updated.**

---

## 5. The matched-binary pair — D&W's hardest case, and the headline number

Both mixtures have δ = 0.500. From `geom_psd_functional_class`, A has S = 0 and B has
S = +8/3.

| mixture | δ | S | φ_DW |
|---|---|---|---|
| A (n_small = 0.500, X_s = 1/28) | 0.500000 | +0.000000 | 0.666900 |
| B (n_small = 0.900, X_s = 0.25) | 0.500000 | +2.666667 | 0.724033 |

> **D&W predicted gap = c₂ · S_B · δ² = 0.0857 × 2.666667 × 0.25 = +0.057133**
>
> **Manuscript measured gap = +0.0635 ± 0.0032** (`\pairDPhi`, `\pairDPhiSe`; p̂ = 1e−3,
> μ_p = 0.3, N = 1500)
>
> **difference = +0.00637 = +1.99 se → AGREE at 2 dof**

For scale, the project's own two-exponent filling closure predicts +0.0604
(`\pairFillingPred`), which is +0.0033 from D&W's number. **A closure that knows nothing
about size ratio, composition or the linear packing model lands within 10 % of a
purpose-built one, on the manuscript's headline result.**

**Four caveats, all load-bearing.**

1. It is a genuine *external* prediction — D&W's constants were fitted before this pair
   existed, to distributions this project did not simulate.
2. S = +2.667 is **outside** their posted range (max |S| = 2.50) by 6.7 %, and that end of
   their range is carried by binary packings only — which is at least the same kind of
   distribution as this pair.
3. It compares **increments**, not levels. D&W's absolutes (0.6669, 0.7240) are
   frictionless; the measurement is at μ_p = 0.3 (0.6042, 0.6677). The comparison is
   licensed only by item A's additive separability, and inherits its ~6 % residual.
4. **One agreement is not a validation.** The refuting test is §4, and it fires.

---

## 6. What this leaves

The picture that survives both halves:

* **Skewness at fixed width is the leading missing variable.** The two-parameter closure
  gets the family ordering, the lognormal/uniform ratio at three CVs, and the headline
  binary gap.
* **It is not the whole story.** Something separates uniform from bidisperse at matched
  (δ, S). Those two differ in their *support*: at CV = 0.5 the uniform runs on
  [0.134, 1.866] (size ratio 13.9) while the bidisperse has exactly two atoms at ratio 3.0.
  So the residual variable is plausibly the support/size-ratio axis that the project's own
  filling functional reads and D&W does not — which is exactly the pair of descriptors
  `geom_psd_functional_class` §3 shows the binary pair could never separate.
