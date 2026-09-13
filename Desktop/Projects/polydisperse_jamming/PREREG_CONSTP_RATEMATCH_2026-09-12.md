# PRE-REGISTRATION — the RATE-MATCHED constant-pressure ladder

**Written 2026-09-12, after `PREREG_CONSTP_PROD_2026-09-12.md` §6 and before any cell of
this run launched.** Grader: `grade_constp_rm.py` → `out/constp_ratematch.txt`. Every bar
below is fixed here.

---

## 1. What this fixes

The production ladder left `Δφ_c`(CV 0 → 0.5) **UNRESOLVED** for one reason: `γ̇` is an
OUTPUT of the packing (`γ̇ = γ_max/(steps·dt)` and `dt` is set by the *smallest* grain), so
it ran **4.70×** across the CV ladder and every CV 0.5 cell landed at `I` = 1.4–1.5e-02
against the 1e-02 quasi-static bar.

This run matches the rate instead of caveating it. **`shear_steps` is set per CV** from
the measured `dt` of the production cells, so all three levels shear at the same `γ̇`:

| CV | measured `dt` | `shear_steps` here | `γ̇` | predicted `I` | est. wall |
|---|---|---|---|---|---|
| 0.0 | 0.250000 | **102 000** | 1.1765e-04 | 6.8e-03 | ~410 s |
| 0.2 | 0.127959 | **200 000** | 1.1747e-04 | 6.7e-03 | ~765 s |
| 0.5 | 0.053180 | **480 000** | 1.1751e-04 | 6.1e-03 | ~1410 s |

Target `γ̇` = **1.1747e-04**, which is CV 0.2's *existing* rate — chosen because it makes
CV 0.2's step count unchanged and therefore a **control** (§3 R3), and because all three
`I` values then clear the bar with ≥ 30 % margin. It is not a value picked to make a
verdict come out: the selection criterion is margin on `I`, stated before the run.

Everything else is `PREREG_CONSTP_PROD` §2 unchanged: lognormal, `μ_p` = 0.3, N = 1000,
`p*` = 3e-04, `K` = 0.15, `γ_max` = 3, `φ_start ∈ {0.55, 0.64}`, seeds {1, 2, 3}, and the
same two absolute preparations at every CV so **no `φ_J` enters the grid**.

## 2. Admission

`PREREG_CONSTP_PROD` §3 unchanged, with one repair it names in its own §6:

- the `z` floor becomes **relative** — `z_win_min / z_ss ≥ 0.5` — because the absolute 2.0
  refused CV 0.5 for having a rattler fraction of 0.42–0.52, i.e. for being polydisperse
  (`z_series` is `z_ALL`). Read relatively the production run gave 0.77–1.00 on all 18
  cells. The genuine unjamming signature is `p_win_min/p*` ~ 3e-05, so **`p_win_min ≥ 0.1 p*`
  is kept unchanged** and is what actually catches unjamming.
- ⚠ This repair is **not** expected to change any verdict: the production refusals were
  decided by `I`, which every cell here is predicted to pass.

## 3. Predictions and bars

**R1 — the rate is matched.** The between-CV spread of realised `γ̇`
(`max/min` over the three CV means) must be **≤ 1.20** — the project's own
`GDOT_MATCH_TOL` — against **4.70** in the production run. And every cell must satisfy
`I ≤ 1e-02`. *Predicted: both hold.* If the spread misses, `dt` moved between runs and the
step counts need recomputing from this run's own `dt`, not the last one's.

**R2 — 🔑 `Δφ_c`(CV 0 → 0.5), and the direction of the move is predicted.**
Bar: within **±0.015** of the common ladder's midpoint **+0.0525**, i.e. in
[+0.0375, +0.0675].
*Predicted: PASSES, and `Δφ_c` moves UP from the unmatched +0.0412.* The direction is not
a guess — the constant-volume campaign's own 30-cell rate-matched control moved `Δφ_c` by
**at most +0.0107 and in the conservative direction**, i.e. matching the rate *increased*
the contrast there. If it moves DOWN instead, the two protocols disagree about what the
rate confound does, which is itself a finding and must be reported as one rather than
absorbed.
⚠ A miss on the high side (> +0.0675) is as informative as a miss on the low side. Do not
treat agreement with the existing number as the success condition.

**R3 — the bit-identity control, and it can fail.** CV 0.2 runs at **exactly the 200 000
steps the production cells used**, same seeds, same everything, under a new tag. Its six
cells must come back **bit-identical** in `φ_ss` to the production run's.
*Predicted: identical.* This is the first determinism check on the constant-pressure
protocol, and it matters because the *isotropic* protocol was found NON-deterministic at
CV 0 (15 duplicated keys returning a different `φ`, worst 5.0e-03). If these six differ,
the servo loop carries the same defect and every constant-pressure number needs a spread
attached.
⚠ Read on `φ_ss`, `p_ss_hat` and `α` — continuous channels. `z_ss` is quantized and
`φ_ss` is an average, so identity on those alone would not settle it.

**R4 — the sign flip, now on admissible cells.** `φ_ss − φ_iso`(3e-04) ran
+0.0099 → +0.0064 → **−0.0072** across CV 0 / 0.2 / 0.5, the last on refused cells.
*Predicted: it stays negative at CV 0.5.* If it does, **the shear-minus-compression offset
changes sign with polydispersity** — a new statement, and the first candidate explanation
for H-3 (the E4/E4m split that `PREREG_ZISO_2026-09-12.md` §8 left unexplained). ⚠ The
`φ_iso` values are N = 1500 and these are N = 1000, so the magnitude mixes a finite-size
difference in; the SIGN at CV 0.5 is −0.0072 against an isotropic seed scatter of
~0.0015, i.e. ~5σ, which finite size would have to be implausibly large to invent.

**R5 — is `α` a coefficient or a rate readout?** CV 0's `α` was measured at
`γ̇` = 6.0e-05 (dilating 6.21 ± 1.14, contracting 3.49 ± 0.68). Here CV 0 runs at
1.1765e-04, **1.96× faster**, with everything else identical. A Reynolds dilatancy
coefficient must not depend on `γ̇`. Bar: **each branch's `α` within 2 seed scatters of its
production value**. *Predicted: within — but this is the prediction most likely to fail,
and a failure means `α` is not a material coefficient at this rate and the derivation
cannot use it.*
And the branch asymmetry itself must survive: between/within ≥ **2.0** (production gave
2.98 and 3.62).

## 4. What would make this uninformative

If R1 misses, nothing downstream is a rate-matched measurement and R2 is quoted as
UNRESOLVED again rather than passed. If R3 fails, every `φ_ss` in the project acquires a
run-to-run spread that has to be measured before any contrast is quoted — that becomes the
result of this run and R2 waits.

## 5. OUTCOME

**Graded 2026-09-12** by `grade_constp_rm.py` → `out/constp_ratematch.txt`. 18 of 18 cells
landed and were admitted (no refusals; `z_win_min/z_ss` 0.83–1.00, servo within ±4.3 %).

| bar | predicted | measured | verdict |
|---|---|---|---|
| R1 rate spread ≤ 1.20, all `I` ≤ 1e-02 | holds | spread **1.001** (production 4.70); max `I` **6.98e-03** | **PASS** |
| R2 `Δφ_c`(0→0.5) in [+0.0375, +0.0675], UP from +0.0412 | passes, moves up | **+0.0533**, up by +0.0121 | **PASS** |
| R3 CV 0.2 bit-identical to production | identical | **6 of 6** identical on `φ_ss`, `p_ss_hat`, `α` | **PASS** |
| R4 `φ_c − φ_iso(p*)` < 0 at CV 0.5 | negative | **+0.0028** (CV 0: +0.0078, CV 0.2: +0.0063) | **FAIL** |
| R5 `α`(CV 0) within 2 sd of production, asymmetry ≥ 2.0 | within | within on both starts; asymmetry **1.12** | **FAIL** as registered |

`φ_c` = 0.5911 ± 0.0010 (CV 0), 0.6020 ± 0.0013 (CV 0.2), 0.6444 ± 0.0043 (CV 0.5); sd
over the six pooled cells per level.

**R2 — the first rate-matched `Δφ_c`.** +0.0533, in band, and moved up from the unmatched
number in the direction the constant-volume control predicted. The two protocols agree on
what the rate confound does. ⚠ **The move is not linear in CV**: +0.0109 over 0→0.2 then
+0.0425 over 0.2→0.5, i.e. +0.054 then +0.142 per unit CV. That is reported, not a bar;
a 4-level ladder would be needed to call it curvature rather than two slopes.

**R3 — the constant-pressure protocol is deterministic**, unlike the isotropic one. Every
constant-pressure `φ_ss` in the project is a single-valued function of its key.

**R4 — the sign flip was the rate.** The −0.0072 at CV 0.5 came from cells at `I` = 1.5e-02;
at matched `I` = 6e-03 the offset is +0.0028 against a combined seed error of ~0.0023, so
it is **indistinguishable from zero and not negative**. The offset `φ_c − φ_iso(p*)` falls
monotonically with CV (+0.0078 → +0.0063 → +0.0028) but does not change sign. **The "H-3
candidate" in R4's prose is withdrawn**. The offset *appears* to shrink with
polydispersity (0.0050 over the ladder), but that is roughly its combined seed error
(`φ_iso` scatter ~0.0015 at each end, `φ_c` 0.0004–0.0018), so the trend is not resolved
either.

**R5 — graded as registered and then set aside.** After this pre-registration was written,
`α` was retracted (`derivations/alpha_from_dem.md`): the fitted rate reads the servo
clamp on the `φ₀` = 0.64 side and the servo compressing an unjammed pack on the
`φ₀` = 0.55 side, and this file's branch labels are swapped (0.55 is *below* `φ_c` and
contracts). Both starts land within 2 sd of production — as a servo readout would, since
the clamp is the same — and the asymmetry collapses to 1.12 at the doubled rate. The bar
stays FAIL; the quantity it grades is not a material coefficient, and the failure is not
evidence about one.

**Reported, not a bar — the rate sensitivity `a = −dφ_ss/dI` from the two runs:**
CV 0: `a` = 0.63 ± 0.15 (I 6.9e-03 → 3.5e-03, `Δφ` = +0.0021 ± 0.0005); CV 0.5:
`a` = 1.16 ± 0.27 (I 6.1e-03 → 1.5e-02, `Δφ` = −0.0100 ± 0.0024). The polydisperse pack
is ~1.9× more rate-sensitive in `I`, at ~1.7σ from equal. Two-point differences on three
seeds; a finding to test, not one to quote.

**Decision:** `Δφ_c`(CV 0 → 0.5) = **+0.053** at matched `I` ≈ 6–7e-03 is the project's
constant-pressure number. `PREREG_CONSTP_PROD` P1 is RESOLVED by this run.
