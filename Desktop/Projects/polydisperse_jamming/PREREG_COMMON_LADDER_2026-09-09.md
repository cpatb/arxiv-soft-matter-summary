# Pre-registration — the common absolute `phi` ladder

**Written 2026-09-09, BEFORE any cell was launched.** The launcher
(`scripts/launch_commonlad.sh`) and the grading script (`grade_common_ladder.py`) are
named below; neither existed when this file was written. This is
`NEXT_STEPS_2026-09-10.md` §2 item 1(a), the run that file calls "the single
highest-value run left".

---

## 1. The defect this run exists to settle

Every `phi_c` ladder in this project is centred on the project's own published
**isotropic** `phi_J` for that cell (`phase8_criticalstate.PHI_J_ISO`), with fixed
offsets. So the fitting window MOVES with CV by exactly `delta phi_J_iso`:

| cv | ladder centre `phi_J_iso` | prod2 rungs (absolute `phi`) | downlad/downlad2 rungs |
|---|---|---|---|
| 0.0 | 0.574533 | 0.5745 … 0.6945 | 0.5145 … 0.5945 |
| 0.2 | 0.588050 | 0.5881 … 0.7081 | 0.5280 … 0.6080 |
| 0.5 | 0.636854 | 0.6369 … 0.7569 | 0.5769 … 0.6569 |

`PREREG_ROOT_PINNING_2026-09-09.md` showed E2's root follows the bottom of its window
(response 1.15 down, 0.06 up), so for E2 the contrast is inherited **by construction**.
E4 (the `z_nr` isostatic crossing) interpolates *inside* the ladder and is structurally
immune to that, and it gives `delta phi_c` = +0.0505 (E4) / +0.0410 (E4m) against the
isotropic +0.0623. **But E4's ladders are still PLACED by `phi_J`**, so the immunity is
an argument, not a measurement.

**A ladder whose rungs are identical absolute numbers at every CV cannot transmit
`delta phi_J` into `delta phi_c` at all.** That is the whole design.

## 2. The design

**Grid (identical at every CV, 6 rungs, spacing 0.025):**

    phi = 0.545, 0.570, 0.595, 0.620, 0.645, 0.670

The ladder TARGETS are realised to 1e-6 at all three CV levels via per-CV `--offsets`
(checked before launch by running `build_grid`). The `phi` the constant-volume
controller actually achieves differs from its target by up to **8.2e-05**; §6 makes
that a bar, not an assumption. 3 CV (0.0, 0.2, 0.5) x 3 seeds (1,2,3) x 6 rungs = **54 cells**, tag
`commonlad`, glob `phase8_commonlad*.jsonl`, which matches no existing analysis.

**Why this window.** It must bracket the root at BOTH ends of the CV axis with rungs on
both sides, using the E4/E4m numbers already on disk (`out/full_ladder_verdict.txt`):

| cv | E4 `phi_c` | E4m `phi_c` | rungs below | rungs above |
|---|---|---|---|---|
| 0.0 | 0.5763 | 0.5948 | 2 (0.545, 0.570) | 3–4 |
| 0.2 | 0.5872 | 0.5971 | 2 | 3–4 |
| 0.5 | 0.6268 | 0.6358 | 4 | 1–2 |

Every root is interior at every CV. The narrowest margin is cv = 0.5 under E4m
(0.6358, one rung of 0.025 below the top rung 0.670), which is declared here as the
run's weakest point.

**Every other knob is held at the centred ladders' setting, deliberately.**
`N` = 1500, `mu_p` = 0.3, family lognormal, `gamma_max` = 1.5,
`shear_steps` = 200 000, same protocol. In particular the run is **NOT rate-matched**:
`gdot` still varies 4.78x across the CV axis exactly as it does in `prod2`/`downlad`,
because the comparison being made is *common window vs centred window* and the two arms
must differ in one knob only. The rate confound is separately measured
(`out/rate_matched_verdict.txt`: +0.0002 to +0.0107, always **positive**, i.e.
conservative) and is applied as a known additive correction, never re-derived here.

## 3. The two hypotheses, and what each predicts

Reference value, from the centred ladders, the estimator that does not inherit the
window: **`delta phi_c`(cv 0 -> 0.5) = +0.0410 (E4m) to +0.0505 (E4)**, midpoint
**+0.0458**.

**H_inherit — the centred-ladder contrast is the ladder placement, not the material.**
Removing the CV-dependence of the window removes the contrast:

    delta phi_c(common)  in  [-0.015, +0.015]     on E4 AND E4m

**H_real — `phi_c` genuinely rises with CV.** The window has nothing to do with it:

    delta phi_c(common)  within 0.015 of the centred value, i.e. in [+0.026, +0.066]
    on E4 AND E4m

**Separation.** The two predicted intervals are disjoint and their midpoints are
**0.0458 apart**. The seed scatter of E4's contrast on the centred ladder is 0.0046
(cv 0.5) and 0.0013 (cv 0), giving a contrast standard error of about 0.0048/sqrt(3)
combined ~ 0.0028. Separation / scatter = **16**, well above the repo's 1.6x
informativeness floor. The test is informative.

**Neither (the third branch, declared in advance).** If `delta phi_c(common)` lands in
(+0.015, +0.026) — or if E4 and E4m land on opposite sides of +0.015 — then the window
carries part but not all of the contrast, and the honest statement is a **reduced**
effect quoted with the common-ladder number, not the centred one. That is a real
outcome and will not be reported as either hypothesis passing.

**Refusal is also an outcome.** If E4 returns NaN or an endpoint at any CV level, the
estimator has failed on this ladder and the run reports that rather than a number.

## 4. The gates and the estimators, fixed now

Reported for **all four gate sets side by side** (ungated / quasistatic / steady /
both), never one alone, exactly as `analyse_phase8.py` already does. Primary estimator
**E4**, secondary **E4m**; E2 is reported but **not quoted** (it is window-dependent by
the previous pre-registration, and this ladder does not repair that). E3 and dense-E1
are printed for completeness.

## 5. The three secondary results this run also produces, declared now

**(a) The window test, on the common ladder itself.** Drop the bottom rung (0.545) and
re-run every estimator. Pre-registered bar: the E4 root moves **< 0.010** at every CV.
If it moves more, E4 is window-dependent too and the whole `phi_c` programme needs a
different observable, not a different ladder.

**(b) P3 becomes EXACT.** On the centred ladders the matched-`phi` pressure contrast
pairs cells that differ by 0.0023 in `phi`; here the three CV levels are given the same six
`phi` TARGETS, so the residual mismatch is the constant-volume controller's own
error rather than a ladder offset. ⚠ Measured on the first 30 cells that landed, the
worst |realised - target| is **8.2e-05**, not zero -- 28x smaller than the centred
ladders' 0.0023, which is the point, but it is not exact and is not quoted as exact. Pre-registered:
the ratio is > 1 (higher CV, lower pressure) at every rung at or above 0.595, and the
range over the six rungs is quoted, never a single rung. No prediction is made about
0.545 and 0.570, where both arms may be below their own root and the ratio is a ratio of
two near-zeros — those rungs are reported and explicitly excluded from the headline.

**(c) A fit-free bracket at every CV.** `p_ss` and `z_nr` are read at six shared `phi`,
so the critical state is bracketed by the two rungs it falls between at each CV with no
estimator at all. This is the number to quote if every estimator refuses.

## 6. What would make this run uninterpretable, said in advance

- More than 3 of the 54 cells failing, or any CV level losing a seed.
- The steadiness gate failing on the two lowest rungs at cv = 0.5 (they are 0.09 below
  that cell's isotropic `phi_J` and may never reach a steady state). This is *expected*
  at 0.545 and is why four rungs sit below the cv = 0.5 root: the estimators must
  survive losing the bottom one.
- A realised `phi` differing from the target by more than 1e-4 at any rung, which would
  mean the ladders are not in fact common.

## 7. The commands

```
# launch (WSL, against the copy in ~/pj; nothing else may be running)
wsl bash -lc 'bash ~/pj/scripts/launch_commonlad.sh'

# grade
python analyse_phase8.py --pattern "phase8_commonlad*.jsonl" --out common_ladder_verdict.txt
python grade_common_ladder.py            # the comparison against the centred ladders,
                                         # the window test, and exact-P3
```

---

## 8. OUTCOME — graded 2026-09-09, `out/common_ladder_grade.txt`

**54 cells, 0 failures, 3 seeds at every CV, worst |realised − target φ| = 8.2e-05
(bar 1e-04).** Every §6 integrity bar passed.

### §3 — **H_real, unanimously**

| estimator | gate | Δ`φ_c` (CV 0 → 0.5) | call |
|---|---|---|---|
| **E4** | ungated | **+0.0554** | H_real |
| E4 | steady | +0.0521 | H_real |
| **E4m** | ungated | **+0.0497** | H_real |
| E4m | quasistatic | +0.0506 | H_real |
| E4m | steady | +0.0494 | H_real |

No estimator and no resolved gate set gives a different call. Per-level means
(E4, ungated): `φ_c` = 0.5745 ± 0.0015 / 0.5850 ± 0.0010 / 0.6299 ± 0.0048, n = 3 each.

🔑 **On a fitting window that is IDENTICAL at every CV, Δ`φ_c` = +0.049 to +0.055 —
squarely on top of the centred ladders' +0.041 to +0.051.** The contrast was not
inherited from the ladder placement. This repairs the strongest doubt raised on
2026-09-08/09 and makes Δ`φ_c` a material measurement again.

⚠ The **quasi-static gate removes almost everything** (0 resolved groups for E4, 1 seed
per level for E4m), because the two lowest rungs are far below jamming and run at
`I` ~ 0.1. The evidence is therefore ungated + steady, and that is how it is quoted.
⚠ **E2 is unresolved on this ladder at every gate set** — which is the previous
pre-registration's finding reproduced, not a new failure: a power law cannot represent
rungs below its own root.

### §5a — the window test **PASSED AND IS VACUOUS**; §5a′ is the one that could fail

`E4` came back **bit-identical (move = 0.0000) on all 9 groups** when the bottom rung was
dropped. That exact zero is the tell: E4 interpolates between the two rungs bracketing
`z_nr` = 4, which are rungs 2–3 at CV 0/0.2 and 4–5 at CV 0.5, so deleting rung 1 cannot
touch it and **any correct implementation must pass** — as would a routine returning a
constant. The number stays printed and the bar stays; a stricter test was added beside
it rather than the bar being loosened.

**§5a′, the placement test.** Same cells, ladders at genuinely different absolute `φ`
(common 0.545–0.670 against centred 0.515–0.695 at CV 0 and 0.577–0.757 at CV 0.5):

```
   worst |E4(common) - E4(centred)| = 0.0051 over 9 cells   (bar 0.010)  PASS
   range 0.0003 to 0.0051 -- NOT zero, so the test had room to fail
```

Moving the ladder by up to 0.06 in absolute `φ` moves E4 by at most 0.0051. **E4's
window-independence was an argument until now; it is a measurement now.**

### §5b — P3 at matched `φ`, and it is sharper than before

| `φ` | `p_ss`(CV 0) | `p_ss`(CV 0.5) | ratio | `Z_nr` | rattlers |
|---|---|---|---|---|---|
| 0.545 | 4.70e-07 | 1.39e-06 | 0.3× | 3.00 → 2.39 | 0.386 → 0.944 |
| 0.570 | 2.29e-06 | 2.87e-06 | 0.8× | 3.82 → 2.55 | 0.159 → 0.863 |
| 0.595 | 4.60e-04 | 8.02e-06 | **57.5×** | 4.83 → 3.11 | 0.048 → 0.715 |
| 0.620 | 7.02e-03 | 3.96e-05 | **177.5×** | 5.83 → 3.62 | 0.007 → 0.562 |
| 0.645 | 1.35e-02 | 8.14e-04 | 16.6× | 6.30 → 4.59 | 0.002 → 0.447 |
| 0.670 | 2.03e-02 | 6.43e-03 | 3.2× | 6.70 → 5.45 | 0.001 → 0.288 |

**Quote 3.2×–177.5×** over the rungs at or above 0.595, all > 1 as pre-registered. The
two lowest rungs are ratios of two near-zeros with both arms below their own root and
were excluded in advance. ⚠ A defect in the grader's own first version keyed the match
on `round(phi, 4)` and **silently dropped the 0.595 row** (realised 0.59505 → 0.5951);
fixed to tolerance matching, and the headline range is unchanged because 57.5× falls
inside it.

### §5c — the fit-free bracket, now on a common grid

`z_nr` crosses 4 between **0.570 and 0.595** at CV 0 *and* CV 0.2, and between **0.620
and 0.645** at CV 0.5. That is the whole result with no estimator at all: **at least one
full grid interval of separation, read off the same six numbers at every CV.**

### H-3 (`PREREG_DST_WINDOW` §5) — **SPLIT, and the split is diagnostic**

| estimator | offset `φ_c(shear) − φ_J(iso, μ_p = 0.3)` at CV 0 / 0.2 / 0.5 | span | verdict |
|---|---|---|---|
| E4 | +0.0000 / −0.0030 / −0.0069 | 0.0069 | constant to <0.010 — **passes** |
| E4m | +0.0157 / +0.0113 / +0.0031 | 0.0126 | **fails** |

The two differ *only* in the isostatic target, and `Z_ISO_MEASURED` is the one that
drifts with CV (4.609 → 4.377). **So E4m's failure is that drift restated, not new
information about protocols**, and H-3 is not settled either way until the isostatic
count is. It is one more item that waits on `Z_ISO`, which `NEXT_STEPS_2026-09-10.md` §0
already names as structural.
