# PRE-REGISTRATION — the extension block: size ratios 4 and 5, and the structure of the branches

**Written 2026-09-14 00:05 CDT, before any cell of this block ran.** Runner
`phase12_mixture.py` (tags `extr4`, `extr5` — chosen NOT to match the `phase12_mixr*`
glob that `grade_mixture.py` and `test_mixture.py` load, so the closed gates M0–M6 and
F1–F3 cannot be re-decided by this data); grader `grade_mixture_ext.py`. Every predicted
number below was printed by `mixture_model.py` / `mixture_closure.py` from the data on
disk at the time of writing and is reproduced by the grader's `--predict` mode.

Grading scale throughout: **`PHI_SCALE` = 2.1e−03**, the protocol's own run-to-run
non-determinism (PREREG_DST_WINDOW §3c), not the seed sem.

---

## 1. Why this block exists

The pairwise-filling closure (`PREREG_MIXTURE_2026-09-13.md`, M0–M6 and F1–F3, all
passed) rests on binaries at size ratio ≤ 3. Its own audit found the two-exponent FORM
**TOO RIGID** (E1): the per-ratio loosening exponent drifts

    r = 1.25  1.50  2.00  3.00
    PA = 2.348 1.760 1.384 0.781        PB = 2.222 2.083 1.983 1.972

while the in-sample residuals of the joint fit (PA 0.961, PB 2.214) are **signed by
ratio** — every r = 1.5 and r = 2 point at X_s ≤ 0.3 is over-predicted, every r = 3
point under-predicted, up to 6.7 scale (rms 2.83, mean |res| 2.29, n = 32). So
everything the closure says about `u = d_fine/d_coarse < 1/3` is an extrapolation of a
form already known to bend. This block measures there.

Two things a bench experiment needs that no cell on disk provides:

1. **Whether the maximum keeps rising and moving with ratio** the way the model says
   (φ_max 0.676 measured at r = 3 → predicted 0.700 at r = 4, 0.718 at r = 5, both at
   X_s ≈ 0.30), or saturates — that is the difference between "use ratio 3" and "use the
   largest ratio you can suspend".
2. **A size-resolved structural signature.** Every rattler number so far is pooled and
   number-weighted, so it reads the small class only. The mechanism makes a claim about
   *which class* carries load on each branch; that has never been recorded.

## 2. Design

Protocol identical to the production block: isotropic compression to `p̂` = 1e−3,
N = 1500, `μ_p` = 0.3, seeds 1/2/3, `BinaryPSD(d10 = 1, ratio, x_small)`, no anchor.

    ratio 4:  X_s ∈ {0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 1.0}   8 × 3 = 24 cells
    ratio 5:  X_s ∈ {0, 0.05, 0.10, 0.15, 1.0}                     5 × 3 = 15 cells

The grids stop where `MIN_PHASE_COUNT` = 50 large grains stops admitting compositions
(r = 4: n_large = 88/67/**53** at 0.20/0.25/0.30; r = 5: 100/**65** at 0.10/0.15, 46 at
0.20 refused). X_s = 0 and 1 are monodisperse controls at d = d_large and d = d_small:
both must reproduce β = 0.5947 ± 2 scale or the block is not comparable to the
production one (gate X0).

**New observable, recorded on every cell:** `measure.rattler_by_class` — per-class
(small/large, split at the geometric mean of the two design diameters, exact for a
binary) rattler fraction at k = 2/3/4 by number, rattler fraction by volume, per-class
coordination, and the SS/SL/LL contact partition. Hand-checked on a 4-grain graph
before launch (this session); it is a pure function of the contact graph and does not
touch `packing.py`.

Admissibility (same as F-gates): no `error`, `converged`, `n_requested` = 1500,
`n_small_design` ≥ 50 and `n_large` ≥ 50 realised. Non-admitted cells are listed, never
silently dropped.

## 3. Predictions (blind — no cell above ratio 3 exists)

Three candidates, all with β = 0.5947 measured and no free parameter at grading time:

- **J** — joint calibration on all 32 ratio ≤ 3 points, (PA, PB) = (0.961, 2.214). *The
  closure as shipped.*
- **R3** — nearest-ratio calibration, (0.781, 1.972) from the 8 ratio-3 points.
- **EXT** — the drift extrapolated: PA(u) = −0.287 + 3.226 u, PB(u) = 1.759 + 0.531 u,
  linear least squares on the four per-ratio fits, giving (0.520, 1.892) at r = 4 and
  (0.358, 1.866) at r = 5. *What "the form bends this way" predicts.*
- Ideal Furnas (no loosening, no wall effect) as the upper bound.

|   r |  X_s |      J |     R3 |    EXT |  ideal | J–EXT sep. (scale) |
|----:|-----:|-------:|-------:|-------:|-------:|-------------------:|
| 4.0 | 0.05 | 0.6102 | 0.6115 | 0.6139 | 0.6260 |  1.8 |
| 4.0 | 0.10 | 0.6265 | 0.6294 | 0.6345 | 0.6608 |  3.8 |
| 4.0 | 0.15 | 0.6438 | 0.6483 | 0.6565 | 0.6996 |  6.0 |
| 4.0 | 0.20 | 0.6620 | 0.6684 | 0.6800 | 0.7434 |  8.6 |
| 4.0 | 0.25 | 0.6813 | 0.6898 | 0.7053 | 0.7929 | 11.4 |
| 4.0 | 0.30 | 0.6997 | 0.7087 | 0.7119 | 0.8302 |  5.8 |
| 5.0 | 0.05 | 0.6118 | 0.6131 | 0.6170 | 0.6260 |  2.5 |
| 5.0 | 0.10 | 0.6300 | 0.6326 | 0.6410 | 0.6608 |  5.2 |
| 5.0 | 0.15 | 0.6493 | 0.6535 | 0.6670 | 0.6996 |  8.4 |

Model maximum (J): r = 4 at X_s = 0.296, φ = 0.7003; r = 5 at X_s = 0.306, φ = 0.7178.
Ideal-limit closed forms (r → ∞): X_s* = (1−β)/(2−β) = 0.2883, φ_max = β(2−β) = 0.8357
(derivation `derivations/furnas_ideal_limit.md`, this session).

## 4. Gates

**X0 controls.** φ(X_s = 0) and φ(X_s = 1) at both ratios within **2 scale** of
β = 0.5947 (seed mean). Fails → the block is not on the production protocol; stop.

**X1 blind transfer of the shipped closure (J).** Over the 9 non-trivial compositions,
seed-mean residual `(φ_meas − J)/scale`: **rms ≤ 4.0** (1.4× the in-sample 2.83) and
**every |res| ≤ 6.7** (the in-sample maximum). A pass means the closure calibrated at
u ≥ 1/3 predicts u = 1/4 and 1/5 no worse than it fits its own data. A fail is reported
as the closure's range of validity, not repaired by refitting on this block.

**X2 which way does the form bend?** On the three compositions where J and EXT are
≥ 8 scale apart — (4, 0.20), (4, 0.25), (5, 0.15) — the candidate with the smaller rms
residual wins. Informative only if the winner's rms is below the loser's by ≥ 1.6× (the
lab's discriminator bar); otherwise "not informative", both numbers printed. The
prediction on record: **EXT wins** (the r = 3 under-prediction is the first sign of the
same drift). If J wins, the drift was not monotone in u and the per-ratio PA sequence
was over-read.

**X3 the maximum keeps rising with ratio.** Measured, model-free: the largest seed-mean
φ over the r = 4 grid exceeds the largest over the r = 3 grid (0.6763 at X_s = 0.30)
by **≥ 3 scale**, and every 0.05 step in X_s from 0 to 0.30 at r = 4 is **positive by
≥ 3 scale** (J predicts every step ≥ +7 scale; so does EXT). Fail → the maximum
saturates below the model at ratio 4 and the bench recommendation is "ratio 3 is enough".

**X4 ratio ordering at fixed composition** (model-free, the mechanism's own sign: both
loosening and wall effect weaken as u falls). At each shared X_s ∈ {0.05, 0.10, 0.15}:
φ(r = 5) ≥ φ(r = 4) − 1 scale and φ(r = 4) ≥ φ(r = 3) − 1 scale, seed means. 6 of 6
inequalities must hold.

**X5 size-resolved structure — exploratory, recorded not load-bearing.** Predictions
from the mechanism, graded ordinally at the primary k = 2:
 (a) `rat_S_k2 > rat_L_k2` at every X_s ≤ 0.30 at both ratios (the small class is the
     interstitial one on the coarse branch);
 (b) `rat_L_k2 ≤ 0.10` at every X_s ≤ 0.30 at both ratios (the large skeleton carries
     load on the whole coarse branch);
 (c) `z_LL` (large–large contacts per large grain) is non-increasing in X_s at r = 4 and
     falls from X_s = 0.05 to 0.30 by more than the seed sd (the skeleton is pried apart
     as the crossing is approached);
 (d) the pooled `rattler_frac_k2` peaks at an X_s BELOW the density maximum at r = 4, as
     it did at r = 3 (0.15 against 0.30) — because it is number-weighted and reads the
     small class.
Any (a)–(d) failing is reported as a limit of the structural picture; none re-decides
X1–X4.

## 5. Analysis registered for after unblinding (not gates)

- Direct per-ratio (a, b) fit (`mixture_closure.fit_ab_direct`) now at six u values
  {0.8, 0.667, 0.5, 0.333, 0.25, 0.2}: the tabulated a(u), b(u) becomes the closure's
  nonparametric form; its identifiability (profile interval width) is printed per u.
- Re-transfer of the tabulated closure to the three continuous families (held-out χ²
  against the 63.0 / 437.7 already on record).
- The bench-hypothesis document `BENCH_HYPOTHESES.md` is written FROM the graded
  verdicts, not alongside them.

## 6. What this block cannot decide

Nothing above ratio 5, nothing at `μ_p` = 0 (the frictionless transfer at r ≤ 3 is on
record; extending it is a separate block), and nothing about the shear protocol — the
transfer of the closure to the constant-pressure `φ_c(CV)` is a separate analysis
(`out/transfer_shear.txt`, this session) on data that already exists.
