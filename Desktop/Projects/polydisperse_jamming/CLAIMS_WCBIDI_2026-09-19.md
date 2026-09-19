# The frictionless bidisperse block: the branch DOES inherit the turnover, and the campaign fails its own worthlessness criterion

**2026-09-19.** Graded by `grade_wcbidi.py` into `out/grade_wcbidi.txt` (rc 1, 3 gates
failed). Registered in `PREREG_WCBIDI_2026-09-19.md` before launch, with two amendments both
written while the campaign was in flight and before any of its records were read.
72 records delivered, `wcbidi` tag, `μ_p = 0`, bidisperse, CV 0.0–0.5, 3 seeds, 4 rungs.

> ## RESULT
>
> **G3, the primary gate: (T) TRANSLATION — the frictionless branch DOES inherit the
> bidisperse family's non-monotonicity.** argmax of `φ_0(CV)` is at **CV = 0.4** on *both*
> channels (extrapolated and fit-free), 1.9 and 2.6 s.e. clear of the runner-up. The
> registered alternative (D), a monotone frictionless branch with argmax at CV = 0.5, is not
> what the data show.
>
> 🔴 **But the verdict is PROVISIONAL, by the prereg's own §5**: *"More than 2 of 72 cells
> cap → G3's argmax is being read off a censored curve."* **14 of 72 capped**, and the
> argmax rung itself has **0 of 3** cells with a fully converged ladder. The criterion was
> written down in advance precisely so this could not be argued away afterwards.

---

## The gates, as registered

| gate | registered | delivered | |
|---|---|---|---|
| **W1** wiring | new CV = 0 cells bit-identical to existing | 39 field mismatches | 🔴 **FAIL — but mis-specified, see below** |
| **W1′** physics | *(added on grading)* | fresh run reproduces an existing record **bit-for-bit at 12 of 12** (seed, rung) pairs | ✅ |
| **W1″** protocol | *(added on grading)* | the tie-break picks the reproduced value on **20 of 24**; 4 disagreements, **all at the bottom rung** | ⚠ |
| **G1** health | 72 of 72 converge | **58 of 72** | 🔴 **FAIL** |
| **W2** (Amd 2) readability | report span/sd | **7.3 – 23.0**, all six groups readable | ✅ |
| **G2** | registered NOT INFORMATIVE in advance | \|t\| = 0.55; `t_crit·se_diff` = 82 % of the effect | — printed, not read |
| **G3** primary | argmax ∈ {0.2,0.3,0.4} = T; = 0.5 = D | **T on both channels**, argmax CV = 0.4 | ✅ **provisional** |
| **G4** blind | `Δφ_0 ∈ [+0.0125, +0.0146]` | **+0.01406** fit-free (inside), **+0.01534** extrapolated (outside) | ⚠ split |
| **G5a** | `W(0) = 0.05723` | 0.05719 | ✅ |
| **G5b** | window range < 2 × paired s.e. | **0.00497 vs 0.00168** | 🔴 **FAIL** |

## 1. W1 was mis-specified, and the campaign is what showed it

The gate asked for bit-identity against "the existing record". **There is no such thing at
CV = 0**: several of those cells carry 2–4 *conflicting duplicate* records, and the gate
silently compared against whichever the loader happened to keep last.

Split into the two claims it was conflating, both pass in the useful sense:

- **The physics.** The fresh bidisperse run reproduces *some* existing record **exactly**,
  at **12 of 12** (seed, rung) pairs. So the generator really does ignore the family label at
  CV = 0, at `μ_p = 0` as well as at 0.3.
- 🔑 **The protocol, and this is new.** Because the fresh run is an independent repetition,
  it *identifies which member of each conflicting pair is reproducible*. The project's
  tie-break ("highest `max_iter`, then first in file order") picks that member on **20 of
  24** multi-record cells. The 4 it misses are **all at the bottom rung `p̂ = 3e-5`**, which
  carries the most leverage in the extrapolation:

  | family | seed | rung | tie-break | reproduced | Δ |
  |---|---|---|---|---|---|
  | lognormal | 2 | 3e-5 | 0.630897017 | 0.630649425 | +2.48e-04 |
  | lognormal | 3 | 3e-5 | 0.633235891 | 0.633183843 | +5.20e-05 |
  | uniform | 2 | 3e-5 | 0.630897017 | 0.630649425 | +2.48e-04 |
  | uniform | 3 | 3e-5 | 0.633235891 | 0.633183843 | +5.20e-05 |

  **Effect on the published anchor `\phiJFrictionless`: −0.00004**, i.e. below its printed
  precision. The rule is not reproducible at the bottom rung and it does not matter here —
  both halves are worth recording, because "a reproduction criterion" is a better tie-break
  than "first in file order" and this campaign shows one is available.

## 2. G1 failed, and the two questions that raises are both answered

**58 of 72 converged.** Per the prereg, reported as an amendment and **not** re-run at a
raised cap.

- **Does the non-convergence correlate with the axis under test?** No. Per CV it is
  10/11/10/10/8/9 of 12, slope **−0.333** per unit CV over a 0.67–0.92 range. It cannot
  manufacture the argmax.
- **Does G3 survive dropping the non-converged records?** **It cannot be tested.** Only 5 of
  18 cells have all four rungs converged, and **none at CV = 0.4 or 0.5** — including the
  argmax rung. This is why the verdict is provisional rather than clean.

## 3. What W2 adds that was not asked for

The frictionless *bidisperse* ladders are **healthy**: span/sd 7.3–23.0 on every CV group.
That is a direct contrast with the frictionless *lognormal* CV = 0.5 group, which
`ladder_dynamic_range_2026-09-19.md` measures at **1.3**. So the ladder degeneracy found
earlier today is family-dependent, not a blanket property of `μ_p = 0`.

## 4. G4: the blind prediction split between channels

Registered **[+0.0125, +0.0146]** from the ratio the other two families give. Delivered
**+0.01406** on the fit-free matched rung (inside) and **+0.01534** on the extrapolated
channel (outside, by 0.0007). The band is not widened. Amendment 2 registered the fit-free
channel as the one to quote where the two disagree, and it does.

## 5. G5b failed

The window's range across the ladder is 0.00497 against a bar of 0.00168 — and the CV = 0.4
rung alone moves it by +0.00433 at \|t\| = 30.3. ⚠ Two reasons not to over-read this: the
CV = 0.5 row is **absent** from the table (the frictional bidisperse cells do not supply a
complete four-rung ladder there after filtering), so the range is over five rungs, not six;
and the same censoring that makes G3 provisional applies here.

## 6. What would settle it

1. **Re-run the 14 capped cells at a raised `MAX_ITER`, under a new tag**, and grade G3 on
   converged-only data. That is the single thing standing between "provisional" and a
   result, and it is ~14 cells.
2. The `μ_p = 1.0` bidisperse block, which does not exist, would say whether the turnover
   survives to the frictional limit.

## 7. Referee

⚠ **Unrefereed.** The gates are as registered and the two mis-specifications found during
grading (W1's reference, and the `mu_p`-less loader key that made W1 report "0 comparisons"
on its first run) are recorded in the grader's own comments rather than edited out.
