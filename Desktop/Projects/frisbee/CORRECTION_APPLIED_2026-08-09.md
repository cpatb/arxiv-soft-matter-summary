# Frisbee: the best-direction claim is withdrawn

**2026-08-09.** The correction had been written and tested and never applied. It is applied
now, in the text, in the figure and in the verdict JSON.

## What the data actually supports

`q_significance.py` compares the three headings pairwise at the **same 12 launch points**,
which is what makes the test possible at all: the across-field standard deviation of Q is
**9 to 14 times** the Monte Carlo standard error at every heading, so an unpaired
comparison is swamped by spatial variation. Bootstrapping over launch points, 20 000
resamples, α = 0.05:

| comparison | paired difference | 95% CI | p | favoured at | verdict |
|---|---|---|---|---|---|
| 90° vs 210° | −0.0040 | [−0.0095, +0.0018] | 0.206 | 4/12 | **not resolved** |
| 90° vs 270° | +0.0086 | [+0.0035, +0.0148] | 0.031 | 11/12 | resolved |
| 210° vs 270° | +0.0127 | [+0.0056, +0.0205] | 0.026 | 10/12 | resolved |

`verdict: {"best_supported": null, "worst_supported": 270}`.

**270° is the worst of the three, resolved against both others. No best heading is
resolvable.** The SSW does have the highest mean Q — but the gap to the easterly is
smaller than the uncertainty on it, and the *sign* is not even consistent across the
field: the SSW wins at only 4 of 12 launch points.

## What changed

1. **`significance_generated.tex` is now `\input`.** It existed, carrying 35 macros
   including `\QbestSupported{none}`, and the manuscript never pulled it in — so every one
   of those macros was inert. Added with a `\MISSINGNUM` fallback block matching the
   file's existing convention, so an absent file still compiles and renders visibly.
2. **§Throwability** now reports the paired test rather than the raw mean ordering, states
   why pairing is necessary (the 9–14× spatial-to-MC ratio), and closes on the worst
   heading rather than the best.
3. **§The hypotheses, tested** — "Best conditions at a light SSW --- SUPPORTED" is now
   **PARTIALLY SUPPORTED**, matching `hypothesis_check.json` verbatim rather than
   inventing a third wording. The old caveat said three samples can rank candidates but
   not locate an optimum; the corrected text notes that was too generous — at 12 launch
   points they cannot reliably rank the top two either.
4. **§Conclusions** now leads with the westerly being the least throwable and states
   explicitly that no best heading is claimed.
5. **The count line** was wrong independently of all this: it read "three falsifiable
   predictions. Two are supported, two of the sub-claims are contradicted" — two and two
   out of three. Now: five sub-claims, one supported, two contradicted, one partially
   supported, one not testable, matching the five entries in `hypothesis_check.json`.

## The figure

`fig4_throwability` previously printed "best" into the image. **A figure is pixels, so no
macro guard in the manuscript could see it** — the same gap that let three withdrawn
claims survive in `make_stf_figure.py` on the simulation side. `_significance_caption()`
now builds the subtitle from `q_significance.json`'s verdict, so the figure cannot assert
a ranking the test declines to support. Regenerated.

## Provenance

`Q.npz` was regenerated with the current `throwability.py` — it previously predated the
analysis code by ~11 h and was missing provenance keys. No headline number moved as a
result.

## Verification

`latexmk` → **26 pp, exit 0, 0 undefined control sequences, 0 undefined references, 0
`\MISSINGNUM` renders.** No occurrence of the old "--- SUPPORTED" string for the
best-direction claim remains.

**Still needs the author:** the author line.
