# The fringe numbers were 1.566x too small

**2026-08-09.**

## The defect

`fringe_analysis.py` was the one script in the pipeline that did not use the *measured* fry
body length. When `body_length.median_px` was absent from a clip's metadata it fell back to
`working_body_px` — the scale the segmentation filters were built at, 21.70 px against a
measured 13.86 px — and it did so **silently**. Every length reported in body-length units
was therefore 21.70/13.86 = **1.566× too small**, and nothing in the pipeline flagged it,
because a silent fallback produces a plausible number rather than an error.

## The fix

1. `body_length` added to the clip metadata, so the measured value is available.
2. `fringe_analysis.py` now **raises** when the measured length is missing, naming the file
   to regenerate, instead of substituting the working scale. A fallback that produces a
   wrong answer is worse than one that stops.
3. The output records `body_px_measured` and `working_body_px` **as separate fields**, so
   the two can never be confused again by a reader of the JSON.
4. Re-run. `fig5_fringe_individuals.png` regenerated from the corrected values.
5. `digest.py` read `xi_bl` where the pipeline produces `xi_decay_bl`; fixed, and the digest
   now runs clean end to end.

## What moved

| Clip | tracks (BL/s) | flow (BL/s) | ratio | NND (BL) |
|---|---|---|---|---|
| | was → now | was → now | was → now | was → now |
| 0544 | 2.76 → **3.96** | 1.54 → **2.41** | 1.79 → 1.64 | 0.85–0.95 → **1.33** |
| 0545 | 4.41 → **5.53** | 3.24 → **4.54** | 1.36 → 1.22 | → **1.33** |
| 0546 | 3.06 → **5.09** | 1.57 → **2.38** | 1.95 → 2.14 | → **1.38** |
| 0547 | 3.33 → **5.78** | 2.37 → **3.79** | 1.40 → 1.52 | → **1.39** |

Track counts also moved (975–1,294 → 1,025–1,300), because the re-run routes through the
same camera-motion correction as the rest of the pipeline rather than bypassing it.

## One conclusion reversed

RESULTS.md §7 previously read *"Fringe NND is 0.85 to 0.95 BL, so under one body length"* and
built on it. The corrected value is **1.33 to 1.39 BL — over one body length.** That specific
claim is **withdrawn**; it was an artifact of dividing by a body length 57% too long. The
surrounding argument survives unchanged: NND is still an upper bound, since detection
recovers roughly half the individuals and a missed neighbour can only make the nearest
surviving one look further away.

**What did NOT move is the section's actual conclusion.** Individuals swim ~1.6× faster than
their local group velocity, so the fluctuation component is comparable to the mean. That is a
*ratio of two speeds sharing the same length unit*, so the calibration error cancels out of
it exactly. The absolute speeds were wrong; the comparison was not.

## Two audit items that did not hold up

`LIMBO_AUDIT_2026-08-08.md` flagged two things here that are **not** defects:

- **fig6.** The audit said it shows a 55 s fission while the docs describe a 5.4 s fusion.
  It does not. `FIGURES.md` calls fig6 "a merger", and RESULTS.md says "the longest single
  episode is a fusion rather than a fission". Those agree. The 5.4 in RESULTS.md is from an
  unrelated table.
- **RESULTS.md §8's claim** that every kinematic number is post-correction was false *only*
  because of the §7 defect above. Now that §7 is re-run it is true, and the sentence carries
  a note recording that it was not true when written — the general point being that a
  pipeline-wide correctness claim is only as good as its least-checked script, and here the
  script that skipped the check was also the one that did not announce it.

## Not done

The folder is ~5,900 files and ~2 GB, unversioned. It is not a git repository and I did not
make it one. If it becomes one, the frames and .MOV files need a `.gitignore` first.
