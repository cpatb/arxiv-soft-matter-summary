# Results

Four clips, 5,780 frames, one afternoon at Montrose Point. All lengths are in body lengths
(BL); see `docs/METHODS.md` §2 for what that unit is. The figures are described in
`figures/FIGURES.md`.

---

## 1. What was measured, and how much of it worked

The school was found in 99.6% of frames (5,734 of 5,780). The misses are real: in IMG_0545
the school swims out of shot for about a second, and the detector reports it absent rather
than fitting an outline to gravel.

Fry body length came out at 13.9, 14.6, 14.5 and 14.3 px for the four clips, with a median
silhouette aspect ratio of 2.6 to 2.8. Four schools filmed independently over four minutes
giving the same fry size to within 5% is a useful internal check. It is what you would
expect from one cohort of one species, and it would be hard to produce by accident from a
segmentation artefact.

## 2. The schools

| | 0544 | 0545 | 0546 | 0547 |
|---|---|---|---|---|
| Diameter *D* (BL) | 32.2 | 38.3 | 33.5 | 33.4 |
| Aspect ratio | 1.76 | 1.53 | 2.13 | 2.34 |
| Circularity | 0.53 | 0.57 | 0.58 | 0.65 |
| Solidity | 0.87 | 0.90 | 0.89 | 0.94 |
| Boundary roughness | 1.38 | 1.33 | 1.32 | 1.24 |
| Projected packing φ | 0.78 | 0.84 | 0.82 | 0.86 |
| Speed (BL s⁻¹) | 2.4 | 4.5 | 2.5 | 3.6 |

The schools are 30 to 40 body lengths across, which is a few hundred to a few thousand fish,
though counting them is exactly what this imagery cannot do. They are mildly elongated
rather than round, and their boundaries run 25 to 40% longer than a circle of the same area
would. Ragged edges, not smooth ones.

## 3. All three collective states appear, within four minutes

This is the result I did not expect.

| Clip | Polarisation Φ | Rotation *M* | State |
|---|---|---|---|
| 0544 | 0.51 | 0.23 | swarm |
| 0545 | 0.69 | 0.08 | polarised |
| 0546 | 0.50 | 0.48 | milling |
| 0547 | 0.77 | 0.09 | polarised |

Schooling fish are known to occupy at least three dynamically stable collective states: a
disordered swarm, a rotating mill, and an aligned polarised group (Tunstrøm et al., PLoS
Comput Biol 2013). All three are here, in four clips filmed at one spot over about four
minutes.

Polarisation alone would have missed it. Clips 0544 and 0546 have almost identical Φ (0.51
and 0.50) and are doing completely different things: one is an internally disordered swarm,
the other is rotating coherently. Only the rotational order parameter separates them.

The state thresholds are conventions rather than measurements, so the raw Φ and *M* values
are in `results/collective_state.json` for anyone who wants to draw the lines elsewhere.

## 4. The interior is more ordered than the edge

In all four clips, without exception, the core of the school is more polarised than the rind.

| Clip | Φ rind | Φ core |
|---|---|---|
| 0544 | 0.44 | 0.56 |
| 0545 | 0.63 | 0.71 |
| 0546 | 0.45 | 0.53 |
| 0547 | 0.76 | 0.78 |

Packing shows the same pattern far more strongly, rising from 0.55-0.66 at the boundary to
0.94-0.97 in the core.

The packing result needs its caveat stated plainly. φ is projected occupancy, not a
three-dimensional volume fraction. The camera flattens a water column, and that column is
deeper in the middle of a lens-shaped school, so part of the rise is geometry rather than
crowding. The polarisation difference carries no such caveat, because it is an angle rather
than an area, and it points the same way in every clip.

## 5. Coordination is real but short-ranged, and not scale-free

Velocity fluctuations are genuinely correlated in space. The test is a permutation null that
shuffles velocities among the same positions, preserving the velocity distribution, the
point geometry, the school outline and the mean-subtraction constraint, and destroying only
the spatial arrangement. The measured correlation length beats that null in every clip at
p < 0.001.

| Clip | ξ (BL) | shuffled null | ratio | flow-resolution floor |
|---|---|---|---|---|
| 0544 | 5.4 | 1.6 | 3.4x | 1.8 |
| 0545 | 2.2 | 1.4 | 1.5x | 1.7 |
| 0546 | 8.8 | 1.5 | 5.8x | 1.7 |
| 0547 | 3.1 | 1.6 | 1.9x | 1.8 |

The shuffled null lands almost exactly on the optical-flow resolution floor, which is a
reassuring sign that both are measuring the same thing, namely nothing.

Two qualifications.

For the two polarised clips the result is marginal. A ξ of 2.2 or 3.1 BL against a floor of
1.7 to 1.8 BL is only 1.3 to 1.7 times the resolution. Those should be read as "short, at or
near the limit of what this footage can resolve", not as measurements.

And ξ does not grow with school size. Binning all frames by school diameter and averaging
C(r) within each bin gives ξ = 4.3, 5.5, 4.8, 4.1, 3.6 BL for *D* = 26 up to 47 BL. That is
a log-log slope of −0.43, not the +1 that scale-free correlation predicts. Over this range,
coordination has a roughly fixed reach of a few body lengths rather than one that scales
with the group.

The negative result is real but weakly powered. The accessible range of *D* spans a factor
of 1.8, and fitting an exponent across less than a decade is not strong evidence for any
particular value. What can be said is that there is no sign of the ξ ∝ *L* growth that
Cavagna and colleagues found in starling flocks across far wider size ranges.

One caveat cuts the other way. The correlation is computed on fluctuations about the
school's mean velocity. For a milling school the collective motion is rotational, and
subtracting a linear mean does not remove it, so some of what is measured as "correlated
fluctuation" in clip 0546 is the mill itself. It is probably not a coincidence that the
milling clip has the longest ξ at 8.8 BL while the two polarised clips have the shortest at
2.2 and 3.1 BL. The method was designed for polarised flocks and it behaves differently on a
mill.

## 6. Splitting and rejoining

Eleven separation episodes were detected across the four clips, after discarding runs shorter
than five frames as boundary flicker.

Every one of them ended in reunion. Eleven out of eleven. Across roughly three minutes of
usable footage, no piece of a school was ever permanently lost.

The separations are substantial rather than a few strays wandering off. At the deepest point
of each episode the largest piece holds only 58 to 74% of the total, so these are genuine
divisions into comparable parts. Individual episodes last from 0.17 s to 55 s (per-clip
medians 0.33, 0.88, 5.44 and 28.2 s), and clip 0544 spends 72% of its 79 seconds as two
distinct schools before they merge.

The longest single episode is a fusion rather than a fission: two separate schools converging
and merging over 5.4 s, with a visible bridge of individuals streaming between them before
the outlines join (fig6). Direction is determined from the data rather than assumed. Calling
every multi-component episode a "split" would have described this one backwards.

A school does not measurably neck before it splits. The obvious hypothesis was that solidity,
occupied area over convex hull, should dip as a waist forms and before the pieces separate.
It does not: median drop 0.014, present in only a third of events, Wilcoxon p = 0.73 across
nine events with a usable onset. Either these splits are not preceded by gradual necking, or
nine events is too few to see it. I would not claim the first from this evidence.

## 7. Individual tracking: what worked and what did not

> **Corrected 2026-08-09.** Every number in this section moved. `fringe_analysis.py` was
> the one script in the pipeline that did not use the *measured* body length: it fell back
> silently to `working_body_px`, the scale the segmentation filters were built at, which is
> 21.70 px against a measured 13.86 px. Every length in body-length units was therefore
> **1.566× too small**, and because the fallback was silent nothing flagged it. The script
> now raises rather than substituting, and the numbers below are the re-run values. One
> conclusion reversed as a result — see the nearest-neighbour paragraph.

Between 1,025 and 1,300 short tracks per clip were recovered at the school fringe. They
serve one useful purpose and one limited one.

The useful purpose is cross-validation. All the kinematics elsewhere come from dense optical
flow, a single method. Individual tracking measures the same quantity a different way. The
two agree in ordering across clips, with individual speeds running 1.2 to 2.1 times higher
than the flow field.

| Clip | tracks (BL s⁻¹) | flow (BL s⁻¹) | ratio |
|---|---|---|---|
| 0544 | 3.96 | 2.41 | 1.64 |
| 0545 | 5.53 | 4.54 | 1.22 |
| 0546 | 5.09 | 2.38 | 2.14 |
| 0547 | 5.78 | 3.79 | 1.52 |

The gap is expected and informative rather than a contradiction. Optical flow averages over a
window of about 1.8 BL, so it measures the local group velocity, while a track follows one
fish including its motion relative to its neighbours. That individuals move roughly 1.6 times
faster than their local group velocity says the fluctuation component is comparable to the
mean. These fish are actively swimming relative to each other, not being carried along.
That conclusion is unchanged by the correction, because the ratio is between two speeds
that share the same length unit — the calibration error cancels out of it exactly. It is the
absolute speeds that moved, not the comparison.

The limited purpose is nearest-neighbour spacing, **and this is where the correction changes
the answer rather than the number.** Fringe NND is 1.33 to 1.39 BL — *over* one body length,
where the uncorrected calibration put it at 0.85 to 0.95 and this document previously
concluded "under one body length". The surrounding logic survives: it is still an upper
bound, because detection recovers roughly half the individuals and a missed neighbour can
only make the nearest surviving one look further away, so true spacing is tighter than
1.33 BL. But the specific claim that fringe fry sit closer together than their own body
length is **withdrawn** — it was an artifact of dividing by a body length that was 57% too
long. Complete individual tracking through the school is not possible from this footage, and
no amount of further processing changes it. The fry occupy a volume and the camera records
one projection of it.

## 8. Camera motion was the largest single threat to these results

The footage is handheld, and a pan moves every point in the frame coherently, which is
exactly what an order parameter is built to detect.

Across 240 sampled frames, removing camera motion changes median polarisation only modestly,
from 0.61 to 0.57. But 32.5% of frames were inflated by more than 0.2, and the worst case is
severe. IMG_0545 frame 1001, with the camera panning at 40 px per frame, reads Φ = 0.998
uncorrected, a school in essentially perfect collective alignment. With the photographer's
motion removed it reads Φ = 0.171.

Every kinematic number in this document is post-correction. Any analysis of handheld footage
that skips this step will find spectacular collective order that is not there.

⚠ **That sentence was not true when it was written, and the exception is worth recording
rather than quietly repairing.** §7's fringe-tracking numbers were produced by the one
script in the pipeline that bypassed both this correction *and* the measured body-length
calibration. It was corrected on 2026-08-09 and §7 now carries re-run values, so the
sentence is true as of that date. The general lesson is the same one this section is about:
a pipeline-wide claim of correctness is only as good as its least-checked script, and the
script that skipped the check was also the one that did not announce it.

---

## What I would say with confidence, and what I would not

Confident: body length is consistent across all four schools; the schools are 30 to 40 BL
across, mildly elongated, with ragged boundaries; all three canonical collective states are
present within four minutes; the core is more polarised than the rind in all four clips;
velocity fluctuations are spatially correlated well beyond chance; every separation episode
reunited; camera correction is essential and was applied.

Suggestive but not established: the absolute value of ξ for the two polarised clips, which
sits near the resolution floor; the absence of scale-free scaling, which rests on a 1.8x
range in school size; the rind-to-core packing rise, which is entangled with projection
geometry.

Not established at all: any count of fish; any metric length or speed; complete individual
trajectories; the presence or behaviour of a guarding parent. No adult fish appears in any
frame I examined, so the parental-care literature is context here rather than an observation.

## Obvious next steps, in order of value per effort

1. Photograph a fry against a ruler. Five minutes on a future visit, and it retro-calibrates
   every number in this document into centimetres.
2. Film from a fixed mount. A phone clamped to the pier rail removes the single largest error
   source and would make the marginal ξ measurements decisive.
3. Get a wider range of school sizes. The scale-free test is limited by a 1.8x span in *D*,
   not by anything about the analysis. More schools of more different sizes would settle it
   either way.
4. Record for longer. Eleven separation episodes is enough to say they all reunited. It is
   not enough to test what precedes a split.
