# Grid independence: the old test was invalid, the mesh error is small, and the directional claim still does not survive

**2026-08-30.** Grid independence at disc altitude was the project's one failed validation
and, per `FRISBEE_STATUS.md` §8, "the only thing gating per-location claims". It was
reported as **FAIL, 57.1% of local speed / 10.3% of u_ref**.

Two conclusions, and they point in opposite directions.

**The old test was not a grid-independence test**, for three independent reasons, and its
headline number was an artefact. A proper four-grid convergence study puts the production
mesh's discretisation error at **~4.5%** in the quantity the paper reports — not 57%, with
an observed order of 2.25 and a 1.19% GCI one level finer. §1–6.

**But the paper's one resolved directional claim does not survive refinement.** Re-running
the entire pipeline — all three headings at twice the resolution, 97,200 fresh trajectories
— **"270° is the worst of the three computed" is not reproduced**: at gridL4 no heading is
separable from any other, and hypothesis H2a reverses sign. The sheltering magnitude and
the Campus Towers finding do survive. §8.

So the old FAIL was pointing at something real, but not at what it measured, and not where
it was applied: the restriction it motivated was placed on per-*location* claims while the
per-*heading* claim, which is exactly as grid-sensitive, was published.

---

## 1. Why the old test was not a grid-independence test

`validate_cfd.test_grid_independence` compared the production mesh against one coarser
mesh at five hand-picked stations and reported the worst as a fraction of the **local**
speed.

**(a) Two grids give a difference, not an error.** `|f_fine − f_coarse|` says the meshes
disagree. It says nothing about how far the fine mesh is from the continuum. An observed
order of convergence needs at least three geometrically similar grids.

**(b) Normalising by the local speed.** Disc-altitude speeds over this pitch run
0.02–1.8 m/s. Dividing by a local speed that approaches zero in a recirculation makes the
metric diverge exactly where the flow is least resolved. The 57.1% is **one station**
sitting in a wake; the same station is **10.3% of u_ref**. The verdict was set by whichever
of five points happened to land in the worst place.

**(c) Five stations under-sampled the field.** Over the whole pitch at 2 m, the
fine-vs-coarse difference is **median 4.8%, p95 14.3%, max 21.2% of u_ref** — worse
coverage than the five stations implied. At the 12 launch points where Q is actually
evaluated it runs 0.34%–13.2%.

The verdict now grades on the u_ref-normalised measure (`validate_cfd.py`), both numbers
are still printed, and the test states in its own output that it cannot establish grid
independence and points at `grid_convergence.py`.

## 2. The root cause: the coarse mesh does not resolve the buildings

Nothing had checked how many cells span the buildings. Campus Towers — the 45 m block nine
metres off the eastern sideline that the paper's headline claim is about — is **24 m across
its narrow horizontal dimension**:

| grid | fieldBox cell | cells across Campus Towers | COST 732 / AIJ (≥10) |
|---|---|---|---|
| gridL1 | 8.0 × 2.0 m | **3.0** | ✗ far below |
| coarse | 4.0 × 1.0 m | **6.0** | ✗ below |
| **production** | 2.0 × 0.5 m | **12.0** | ✓ marginal |
| gridL4 | 1.0 × 0.25 m | **24.0** | ✓ comfortable |

Five of the six largest buildings are below the guideline at h = 4 m.

**So the original test compared a resolved mesh against an unresolved one.** A gap between
those two measures the unresolved mesh, not the resolved one.

## 3. A proper four-member family

The existing ladder's `fine` entry is unusable for this: it moves `dxy`, `near` **and**
`field`, so it is a different mesh family and would corrupt a Richardson fit. Two entries
were added that move **only** the fieldBox level, keeping the meshes geometrically similar
at exactly r = 2 (`generateCases.py`, `gridL1` and `gridL4`).

Every member was run through the identical pipeline: solve to 1400 iterations, restart from
the 1250 write, cycle-average over 800 iterations, sample onto the same 201×141×6 grid.

> **An undocumented protocol detail.** The production cases' real `endTime` is **1400**,
> not the 2500 default; they terminate cleanly there. `FRISBEE_STATUS.md`'s reproduce block
> omits that flag, so following it as written produces a different run. The new members
> were matched to 1400 by hand.

| grid | cells | max skewness |
|---|---|---|
| gridL1 | 260,399 | 6.96 |
| coarse | 276,466 | — |
| production | 446,408 | 7.93 |
| gridL4 | 1,687,455 | 7.97 |

Solver health is comparable across members (gridL4's final residuals are slightly lower
than production's; both bound k on ~99% of iterations, a pre-existing characteristic).

## 4. Result: it converges, at the resolved end

Pitch-mean horizontal speed at 2 m, heading 90°:

| grid | h | mean \|U_h\| |
|---|---|---|
| gridL1 | 8.0 m | 0.9456 |
| coarse | 4.0 m | 0.9930 |
| production | 2.0 m | 1.2057 |
| gridL4 | 1.0 m | 1.2504 |

**Coarse triplet (production/coarse/gridL1): p = −2.17, refused.** Refining 8→4 m moved the
answer 0.047; refining 4→2 m moved it 0.213, four times *more*. Differences growing under
refinement is divergence, and `grid_convergence.py` declines to emit a GCI rather than
quoting a meaningless one. That is exactly what a family straddling *building unresolved →
building resolved* should look like: Richardson extrapolation across a regime change is
undefined.

**Fine triplet (gridL4/production/coarse): p = 2.25** — close to the scheme's formal second
order. Richardson extrapolation `f_ext = 1.2623 m/s`; **GCI on gridL4 = 1.19%**.

Against that extrapolated continuum value:

- **production is 4.5% low** — 24.11% of u_ref against an extrapolated **25.25%**
- gridL4 is 0.9% low

## 5. The quantity that matters is not the one that was tested

The sheltering headline rests on `U_mean`. The **270-is-worst verdict rests on the
turbulence factor (1 − I)**, so a study that convergence-tests only the speed says nothing
about the claim the paper actually resolves. All three:

| quantity | ratio | p | f_ext | err(production) | err(gridL4) |
|---|---|---|---|---|---|
| `U_mean` | 4.759 | **2.25** | 1.2623 | 4.5% | 0.9% |
| `k_mean` | 2.098 | 1.07 | 0.1754 | 18.3% | 8.7% |
| `I_median` | 1.286 | **0.36** | 0.1809 | *extrapolation refused* | — |

`I_median`'s fitted order is far below the scheme's formal order, so the family is **not in
the asymptotic range** for it: Richardson then projects a large change out of small
grid-to-grid differences, and the resulting "25.7%" is an artefact of the tiny exponent
rather than a measurement. The script flags it and reports what *is* measured instead — the
**production→gridL4 change of 4.5%**.

> **A correction worth recording.** On the three coarsest grids the turbulence intensity
> looked converged: it moved 26% from level 1→2 and only −0.6% from 2→3. That was a
> coincidence. gridL4 shows it continuing to drift. Three grids would have produced a
> confident and wrong "converged" verdict here; only the fourth grid exposed it.
>
> A second correction: the first aggregation used was the *mean* of the per-cell intensity,
> which has exactly the local-normalisation fragility this document criticises in the old
> test — it swings 0.33 → 0.56 across the family while the median moves 0.233 → 0.217,
> because a handful of near-stagnant cells dominate a ratio whose denominator approaches
> zero. The median is what is graded.

Cellwise over the pitch (three finest grids): 58.9% of cells converge monotonically, 25.6%
sit in a usable asymptotic range (up from 9.7% on the coarse triplet), median observed
order 1.71, median GCI 4.28%, p95 56.4%. The field is not uniformly resolved and that is
reported as a distribution rather than collapsed to a max.

## 6. What this changes for the paper

- **"Grid independence FAILED at 57.1%" should be withdrawn.** It is not a measurement of
  the production mesh's error. The defensible statement is that the production mesh carries
  a **~4.5% discretisation error** in the reported disc-altitude speed, with an observed
  order of 2.25 and a 1.19% GCI on the mesh one level finer.
- **The sheltering headline moves slightly.** Heading 90 goes from 24.11% to an extrapolated
  25.25% of u_ref. The published "24–28%" band would need re-stating once all three headings
  are extrapolated; the physical claim (heavily sheltered, ~a quarter of the reference) is
  untouched.
- **The per-location gating should be re-examined.** `FRISBEE_STATUS.md` withholds
  per-launch-point claims because of this failed test. The basis for that restriction was
  the invalid 57.1%; whether it should still stand is now a judgement about a 4.5% error,
  not a 57% one.

## 7. Still open

The heading **comparison** — the paper's one resolved directional claim — is not settled by
a single-heading study. 210° and 270° are running at gridL4 now (1,746,021 cells for 210°),
after which the ranking can be tested at a resolution where every major building clears the
guideline. That result will be appended here.

A formal GCI restricted to *resolved* meshes only is not reachable on this hardware: it
would need a fifth member at h = 0.5 m (~13M cells), against 15 GB of RAM.

## Reproduce

```
cd SeanEarl_DiscFlight/P3_CFD
python3 generateCases.py --headings 90 --resolution gridL1 --tag _gridL1 --nproc 11
python3 generateCases.py --headings 90 --resolution gridL4 --tag _gridL4 --nproc 11
# set endTime 1400 / writeInterval 1250 in each system/controlDict to match production
bash runGridStudy.sh cases/dir_090_gridL1 11 800
bash runGridStudy.sh cases/dir_090_gridL4 11 800
python3 grid_convergence.py --heading 90 --out grid_convergence.json
```

---

# 8. The heading comparison at gridL4: the directional claim does not survive

210° and 270° were run at gridL4 (1,746,021 and 1,709,091 cells), then the **entire
downstream pipeline** was re-run on the three-heading gridL4 field: 97,200 trajectories
(202.6 min), the same launch grid, throws, bearings, seed and `sweeping` gust model, with
**only the CFD field changed**. Every driving parameter in `meta_json` — `u_ref`, rotation,
`z0_inlet`, `z0_ground`, `u_star`, `k_inlet_ref` — is identical between the two fields; only
`effective_field_cell_m` differs.

## The result

| comparison | production (2 m) | gridL4 (1 m) |
|---|---|---|
| 90 vs 210 | p = 0.229 — not resolved | p = 0.602 — not resolved |
| 90 vs 270 | **p = 0.023 — RESOLVED** | p = 0.602 — **not resolved** |
| 210 vs 270 | **p = 0.0052 — RESOLVED** | p = 0.166 — **not resolved** |

```
verdict, production : worst_supported = 270,  beats = {90: [270], 210: [270], 270: []}
verdict, gridL4     : worst_supported = None, beats = {90: [],    210: [],    270: []}
```

**The paper's one resolved directional claim — "270° is the worst of the three computed" —
is not reproduced at the finer mesh.** At gridL4 no heading is separable from any other.

## Why: the between-heading turbulence contrast collapses

The verdict rests on the turbulence factor, and refinement compresses the very quantity it
depends on. Mean turbulence intensity at the 12 launch points:

| | 90° | 210° | 270° | spread |
|---|---|---|---|---|
| production | 0.2578 | 0.3701 | 0.4417 | **0.184** |
| gridL4 | 0.4006 | 0.4633 | 0.4657 | **0.065** |

A **2.8× compression**. The headings become nearly indistinguishable in turbulence
intensity once the buildings are resolved to 24 cells rather than 12. Consistent with this,
the pitch-mean `I_median` gap between 270 and 210 shrinks 7.5×, from 0.0527 to 0.0070, and
the 210/270 mean-speed ordering flips outright.

## It is the field, not the metric

At gridL4 the turbulence factor `(1 − I)` saturates to zero at **2 of 12 launch points**
(one each for 90° and 270°), which `throwability.py`'s own docstring warns is the wrong
functional form for sheltered sites. That is a real caveat, but it is **not** what removes
the significance — dropping the saturated points does not restore it:

| 210 vs 270, gridL4 | n | diff | p_t | p_w |
|---|---|---|---|---|
| all launch points | 12 | +0.0130 | 0.055 | 0.064 |
| saturated points excluded | 10 | +0.0113 | 0.120 | 0.131 |

(uncorrected; Holm over three comparisons makes both weaker still)

## A confound nobody had checked

The production **270° mesh is the lowest-quality of the three**: max skewness **16.68**
against 7.93 (90°) and 7.34 (210°). Refinement halves it to 8.79. So the heading the paper
declares worst was computed on the one production mesh carrying a quality outlier, and
refining that mesh both fixes the skewness and removes the separation. This does not prove
the verdict was an artefact, but it is a second mechanism pointing the same way.

## What survives and what does not

| claim | production | gridL4 | verdict |
|---|---|---|---|
| Sheltering ~24–28% of u_ref | 24.11 / 28.07 / 27.76% | 25.01 / 25.06 / 27.60% | **SURVIVES** (band becomes 25–28%) |
| Campus Towers wake, not the lake IBL, controls the flow | — | — | **SURVIVES** (a magnitude claim, and refinement only sharpens the wake) |
| **270° is the worst of the three** | RESOLVED | **NOT RESOLVED** | **DOES NOT SURVIVE** |
| 90 vs 210 not separable | not resolved | not resolved | SURVIVES |
| H1 easterly Q gradient | SUPPORTED | SUPPORTED | SURVIVES |
| H2a westerly Q patchier | SUPPORTED (CV .255/.245/.302) | **CONTRADICTED** (CV .494/.362/.376) | **FLIPS** |
| H2b forehand benefits most under westerly | CONTRADICTED | CONTRADICTED | SURVIVES |
| H3 best = 210 | PARTIALLY SUPPORTED | PARTIALLY SUPPORTED | SURVIVES (weaker) |

`FRISBEE_STATUS.md` §9 argued that what survives the grid-independence failure is "the
*magnitude* of the sheltering, because it is a statement both meshes agree on", and that
what does not is "the spatial detail of Q within the field — so no ranking of individual
launch positions is claimed". That reasoning is vindicated, and it should have been applied
one step further: **the ranking of headings is exactly as grid-sensitive as the ranking of
locations, and it was claimed anyway.**

## What this means

1. **Withdraw "270° is the worst of the three computed"**, or state explicitly that it is
   not reproduced at twice the resolution. It is currently the paper's only resolved
   directional result.
2. **H2a must be re-stated**; it reverses sign under refinement.
3. **The sheltering magnitude and the Campus Towers finding stand.** They are magnitude
   statements and they move by ~1 point of u_ref.
4. The old restriction on per-launch-point claims should stay — for the right reason now.

## Cost, for the record

| stage | wall clock |
|---|---|
| gridL1 mesh + solve + average | ~31 min |
| gridL4 90° mesh + solve + average | ~11 min + 91 min + ~55 min |
| gridL4 210° end to end | 5,433 s (91 min) |
| gridL4 270° end to end | 5,264 s (88 min) |
| Monte Carlo on gridL4, 97,200 trajectories | 202.6 min |

About 9 hours on 11 cores, all of it unattended.
