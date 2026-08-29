# The direction ordering is only half resolved

**2026-08-07.** Result of `q_significance.py`, which is new. Regenerate with
`python q_significance.py` from this directory; it writes `significance_generated.tex`
(35 macros) and `q_significance.json`.

## The claim under test

`FRISBEE_STATUS.md` §6 and `results_generated.tex` report

> Best of the three computed = 210° (SSW). Worst = 270° (westerly).

with `\QbestHeading{210}` and `\QworstHeading{270}` as generated macros. The three means
are Q̄ = 0.0443 / 0.0483 / 0.0356, sitting inside heavily overlapping Q_min–Q_max ranges,
and nothing in the pipeline had asked whether that ordering survives its own noise.

## Two questions that were being conflated

**Is the reported Q range an uncertainty?** No. Q_min–Q_max is the spread across the 12
launch points, which is spatial structure. The Monte Carlo contribution is separate and
much smaller:

| heading | Q̄ | spatial sd | Monte Carlo se of the mean | ratio |
|---|---|---|---|---|
| 90° | 0.0443 | 0.0087 | 0.00099 | 9× |
| 210° | 0.0483 | 0.0127 | 0.00093 | 14× |
| 270° | 0.0356 | 0.0094 | 0.00078 | 12× |

So the sweep is adequately sampled and the Q map is real structure. **But the paper must
not present Q_min–Q_max as an error bar**, because it is not one.

**Is the ordering resolved?** Partly. The 12 launch points are common to every
direction, so the right comparison is paired, which is far more powerful than comparing
means against the spatial spread. Paired t plus Wilcoxon, Holm-corrected across the three
comparisons, with a paired bootstrap CI on each difference:

| comparison | mean ΔQ | 95% CI | p (t, Holm) | p (W, Holm) | points favouring | verdict |
|---|---|---|---|---|---|---|
| 90 vs 210 | −0.0040 | [−0.0095, +0.0018] | 0.206 | 0.301 | 4/12 | **not resolved** |
| 90 vs 270 | +0.0086 | [+0.0035, +0.0148] | 0.031 | 0.004 | 11/12 | resolved |
| 210 vs 270 | +0.0127 | [+0.0056, +0.0205] | 0.026 | 0.019 | 10/12 | resolved |

A pair counts as resolved only if both the parametric and the rank test survive
correction, which is the conservative reading and the right one at n = 12.

## What this changes

**Supported:** 270° (westerly) is the worst of the three computed, resolved against both
others, in 10–11 of 12 launch points, with a CI that excludes zero.

**Not supported:** that 210° is the best. Its 95% CI against 90° spans zero and only 4 of
12 launch points favour 90°, so the two are indistinguishable on this evidence.

Concretely, in the manuscript:

1. `\QbestHeading` should not be asserted. Either withdraw the "best" claim or state it
   as "90° and 210° are not separable; both are better than 270°."
2. The H3 hypothesis verdict, currently **SUPPORTED** ("best conditions at a light SSW
   ~210°"), should become **PARTIALLY SUPPORTED**: 210° is among the better directions
   but is not resolved from 90°.
3. Wherever Q_min–Q_max appears, say it is the spatial range across launch points, and
   quote the Monte Carlo standard error separately if an uncertainty is wanted.

None of this weakens the paper's central result. The headline is the sheltering
finding — disc-altitude wind at 24–28% of the 10 m reference in every direction — and
that is unaffected. What changes is a secondary ordering claim that the data cannot
carry, and stating it correctly is cheaper than having a referee find it.

## A provenance trap this analysis walked into, and the guard added for it

`throwability.py` supports `--exclude-throws`, and its own source records that the
shipped `Q.npz` was built with `--exclude-throws hammer` while `runDownstream.sh` passed
no such flag, so the committed script does not reproduce the committed artifact. The
shipped `Q.npz` also predates the provenance keys, so it does not say which throw subset
it used.

Averaging the binomial error over all three throw types when the artifact used two
inflates the Monte Carlo error by 1.5× and would have made the sweep look under-sampled
when it is not. `resolve_throw_subset()` therefore uses the recorded provenance if
present and otherwise recovers the subset by matching the stored success rate, reporting
which route it took. On the shipped artifact it recovers
`right_backhand, right_forehand` exactly, confirming the hammer exclusion independently.

**Action:** regenerate `Q.npz` with the current `throwability.py` so the provenance keys
are recorded, and make `runDownstream.sh` pass `--exclude-throws hammer` so the script
reproduces the artifact. Until then `\QsigThrowProvenance` reads "inferred by matching",
which is honest but should not be permanent.

## Tests

`test_q_significance.py`, 13 checks: Holm against a hand computation, bootstrap
reproducibility and coverage, the throw-subset recovery, and the three headline
conclusions pinned so a future edit cannot quietly change them.
