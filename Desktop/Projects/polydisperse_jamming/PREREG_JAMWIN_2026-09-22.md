# PREREG — `jamwin`: does the ACTIVATION SHAPE decide whether an impacting drop jams, and does it do so with the sign the reduced framework predicts?

Written 2026-09-22, third block, **before any cell runs**. Solver:
`Github/Drop Impact/basilisk-drop-impact/wc_drop_impact.c`. Companion DEM campaign:
`PREREG_ACTSHAPE_2026-09-22.md` (measures the shape `δ` this campaign varies).

## 0. What is being tested

The reduced framework (`three_knobs.py`) says a suspension at `x = φ/φ₀` can jam only in
the window `χ ≤ x < 1`, and does so where the frictional fraction reaches
`f_c = (1 − x)/(1 − χ)`. With the solver's activation `f_ss = exp(−(σ*/σ)^δ)`, a material
point in steady state jams iff its stress exceeds

&nbsp;&nbsp;&nbsp;&nbsp;`σ_jam = σ* · [ln(1/f_c)]^(−1/δ)`.                                    (1)

(1) is a closed form (checked in `test_actshape.py` against a Weibull sample). It makes a
**sign-reversing** prediction, which is what makes it testable rather than illustrative:

- **Deep in the window, `f_c < 1/e`** (`ln 1/f_c > 1`): `σ_jam/σ*` RISES with `δ`, so a
  SHARPER activation jams LESS.
- **Near the window edge, `f_c > 1/e`**: `σ_jam/σ*` FALLS with `δ`, so a sharper activation
  jams MORE.

The companion DEM campaign asks whether polydispersity moves `δ`. This one asks whether
`δ` moves the drop — in the full DNS, with inertia, capillarity, a finite state clock
(`WC_CLOCK = 1`, `γ₀ = 0.2`) and a spatially varying stress, none of which (1) contains.

## 0b. A defect in the observable, found while designing this

`wc_drop_impact.c` reports `jam_frac` as the volume fraction with **`f̂ > 0.5`**, a fixed
threshold that knows neither `x` nor `χ`. The closure's own jamming condition is
`φ ≥ φ_J(f̂)`, i.e. `f̂ ≥ f_c`, and for `x < χ` it is **unreachable at any `f̂`**. The
shipped `cs_fit` preset has `x = 0.671`, `χ = 0.781`, so it can never jam, and the
programme's "arrest" readings of `jam_frac > 0.5` on it are readings of ACTIVATION, not of
jamming. A new env override **`JAM_FHAT`** (default 0.5, output byte-identical when unset)
lets `jam_frac` be read at the closure's `f_c`. Arm N below exists to show the conventional
reading reporting "jammed" on a suspension that provably cannot jam.

## 1. Grid

Base preset `cs_fit` (`η_s` = 0.003939 Pa·s, `EPS_MIN` = 0.01, `MP_EXPO` = 2,
`DIL_COEFF` = 0, `WC_CLOCK` = 1), `U₀` = 2.0 m/s, `θ₀` = 5°, `TEND` = 6, `CFLNUM` = 0.10,
`SHAPE_MODE` = 1, **L9**, overridden by env:

`PHI0` = 0.62, `PHIM` = 0.496, so **`χ` = 0.80** (window width 0.20), and

| arm | `x` | `PHI` | `f_c` | `JAM_FHAT` | role |
|---|---|---|---|---|---|
| **A** | 0.95 | 0.589 | 0.25 (< 1/e) | 0.25 | deep in window |
| **B** | 0.85 | 0.527 | 0.75 (> 1/e) | 0.75 | near the edge |
| **N** | 0.75 | 0.465 | unreachable | 0.5 | outside: cannot jam |

× `δ` ∈ {0.7, 1.3169, 2.6} (the `deltastiff` ladder) × `σ*` ∈ {30, 300, 3000} Pa = **27 cells**.

The `σ*` ladder spans the impact's stress scales: `ρU₀²` = 4800 Pa; the unactivated
viscous stress `η(f̂=0) U₀/D₀` ≈ 1050 Pa (A), 115 Pa (B), 42 Pa (N). `δ` can only matter
where `σ ~ σ*`, so a single `σ*` could sit where every `δ` saturates or none activates.

Plus the identity gate (§2, G0): **2 cells** reproducing `deltastiff_2026-09-03/DS_d1.3169_L9`.

## 2. Gates

- **G0 — IDENTITY (primary).** The new binary, run on `DS_d1.3169_L9`'s exact configuration,
  must reproduce its `spread_*.dat` **byte-for-byte (md5)** with `JAM_FHAT` unset (G0a) and
  with `JAM_FHAT = 0.5` set explicitly (G0b). The build source differs from the one that
  made `DS_d1.3169_L9` only by the `JAM_FHAT` lines (verified: the old build source is
  md5-identical to repo HEAD).
- **G1 — PROVENANCE.** Every cell's stderr banner must print the `phi`, `phi0`, `phiM`,
  `sigma*` and `delta` it was asked for, and `JAM_FHAT=` when set. The preset block runs
  before the env overrides; this repo corrupted an anchor on 2026-08-31 by trusting that.
- **G2 — HEALTH.** `rc = 0` and the run reaches `t/τ ≥ 5.9`. Failing cells are printed and
  excluded; the count is printed. ⚠ **A cell that dies is not a jam.**

## 3. Observable and predictions

**Primary observable:** `J = max over t of jam_frac` (volume fraction with `f̂ > f_c`),
per cell. Secondary: `f̂_mean` at its maximum, `β_max`, and `u_rms_jam/u_rms` late.

A `(arm, σ*)` pair is **INFORMATIVE** only if `J` is neither ≤ 1e-3 nor ≥ 0.99 at all three
`δ` (otherwise every `δ` is saturated or silent and the pair decides nothing). Only
informative pairs are graded; the count of informative pairs is printed.

- **P1 — arm A (f_c = 0.25):** on every informative `σ*`, `J` is **non-increasing** in `δ`,
  and `J(0.7) > J(2.6)`.
- **P2 — arm B (f_c = 0.75):** on every informative `σ*`, `J` is **non-decreasing** in `δ`,
  and `J(2.6) > J(0.7)`.
- **P3 — the reversal is the claim.** P1 and P2 must BOTH hold on at least one informative
  `σ*` each. A monotone response of the same sign in both arms (e.g. `δ` just speeding the
  whole activation) **refutes** the reduced framework's reading of `δ`, and would be
  reported as such.
- **P4 — arm N:** `J` (at the conventional 0.5) exceeds 0.1 at `σ*` = 30 Pa for all `δ`,
  while the closure guarantees the jamming gap there never falls below
  `1 − x/χ = 0.0625` > `EPS_MIN`. That is the demonstration of §0b. **Expected, and graded
  only as a demonstration; it is not a test of the physics.**

**What would make this uninformative:** at L9 the `f̂` family is grid-sensitive
(`deltastiff`: `f̂_mean` 21.9× between L9 and L10 at `cs_fit`). The predictions here are
**SIGNS within one grid**, which a grid error can still flip if it is `δ`-dependent (and
`deltastiff`'s D2 was designed precisely to ask whether it is). If P1–P3 pass, an L10 rung
on the informative pairs is the registered follow-up before anything is quoted.

---

## Amendment 1 — written 2026-09-22 ~13:50 CDT, BEFORE any grid cell was read

Only `G0a_unset`/`G0b_half` had been pulled (`out/jamwin/` held those two and the reference).

**§0's prediction (1) is a STRESS-controlled threshold. The solver is not stress-controlled.**
`wc_relax` activates on `σ = η(f̂)·γ̇` — the lagged viscosity times the local strain rate —
and `η(f̂)` rises steeply towards `f_c`. So activation feeds back on itself, and a material
point relaxing from `f̂ = 0` reaches `f_c` iff `f̂ < F(f̂) = exp(−(σ*/(η(f̂)γ̇))^δ)` on all of
`[0, f_c)`. The threshold is then set by the BOTTLENECK `f*`, not by `f_c`:

&nbsp;&nbsp;&nbsp;&nbsp;`γ̇_c · η_s/σ* = max_{0<f<f_c} (ln 1/f)^(−1/δ) / (η(f)/η_s)`,

and `∂γ̇_c/∂δ` has the sign of `ln ln(1/f*)` — the reversal sits at `f* = 1/e`, not `f_c = 1/e`.
Computed for this grid (`derivations/jam_rate_driven_2026-09-22.py`, brute-force relaxation
check 6/6 PASS): `f*` = 0.011–0.050 (arm A) and **0.053–0.295 (arm B), below 1/e in both
arms**, and `γ̇_c` RISES with `δ` in both: A 3.41e-4 / 7.59e-4 / 1.29e-3, B 7.03e-3 / 9.43e-3 /
1.31e-2 (in `σ*/η_s`, at δ 0.7 / 1.32 / 2.6).

**So two predictions are now registered, and the DNS discriminates between them:**

- **(S) stress-thresholded** (§3 as written): A falls with `δ`, **B rises** — reversal (P1–P3).
- **(R) rate-driven bottleneck** (this amendment): **A and B BOTH fall with `δ`** — no reversal.
  P1 is shared; **P2 flips**; P3 is predicted to FAIL.

The gates, bars and informative-pair rule are unchanged. The grader reports which of (S)
and (R) the arm-B sign matches. If B rises, the rate-driven reading of the solver is wrong;
if B falls, the stress-controlled derivation does not describe the solver's drop, and eq.
(1)/(8) must be restated in rate-driven form before it is used to read a DNS.

---

## Amendment 2 — written 2026-09-22 ~13:55 CDT, BEFORE any arm-B or arm-N cell had finished

**Health finding (arm A, 3/3 finished cells):** every arm-A cell dies within `t/τ ≤ 0.1`
(`dt`-collapse to ~2e-8, or `umax` blow-up), and in each the rows before death already show
`jam_frac(f̂ > f_c)` of 0.70–0.86. The solver cannot integrate a drop whose jammed region
sits at `η_max = η_s/EPS_MIN²` = 10⁴ `η_s`. Under G2 these cells are EXCLUDED and their
pre-death values are **not** the registered observable. (They were seen while diagnosing
the deaths; arm A carries only P1, which (S) and (R) share, so nothing that discriminates
was seen.)

**Contingency, fixed now.** If arm B's cells also fail G2, P1–P3 are recorded
**UNDETERMINED (observable unmeasurable)**, and the follow-up campaign `jamwin2` uses:

- **Outcome, robust to death-after-jam:** JAMMED if `jam_frac(f̂ > f_c) ≥ 0.5` at any written
  row; NOT JAMMED if the run reaches `t/τ ≥ 5.9` without that; UNDETERMINED if it dies with
  `max jam_frac < 0.5`. ⚠ A death is still never a jam; only a written row above 0.5 is.
- **Forcing ladder:** arm B only (it alone separates (S) from (R)), `σ*` on a factor-√10 ladder
  {10, 31.6, 100, 316, 1000, 3162, 10000} Pa × `δ` {0.7, 1.3169, 2.6} = 21 cells, L9,
  everything else as §1.
- **Observable:** `σ*_c(δ)` = the largest `σ*` rung that JAMS (the jam/no-jam boundary), read
  only if the outcomes are monotone in `σ*` at that `δ` (else UNDETERMINED for that `δ`).
- **Predictions:** (S) `σ*_c` RISES with `δ` in arm B (sharper activation jams at lower
  forcing near the window edge); (R) `σ*_c` FALLS with `δ`. A tie across all three `δ` at
  this ladder spacing is UNDETERMINED, not support for either.
