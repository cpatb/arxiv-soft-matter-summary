# Methods, and the decisions behind them

Written so that a reader can tell which choices are principled, which are pragmatic, and
which were forced by something going wrong. The failures are included because they are
the part most likely to matter to anyone reusing this.

---

## 1. Converting the video

The four `.MOV` files are H.264, 1920×1080, 30000/1001 fps, from an iPhone 11. Three
carry a −90° display rotation in the container; one does not.

Frames were extracted with

```bash
ffmpeg -v error -i IMG_0544.MOV -vsync 0 -q:v 2 data/frames/IMG_0544/%05d.jpg
```

Two deliberate choices.

**Rotation is baked into the pixels.** ffmpeg applies the container's display rotation by
default. Working from extracted frames rather than decoding the `.MOV` in OpenCV sidesteps
a genuine trap: how `cv2.VideoCapture` treats QuickTime rotation metadata has varied
across versions, and a silently transposed frame would corrupt every direction-dependent
result. Orientation, polarisation and velocity would all be wrong while every scalar still looked fine.

**JPEG at `-q:v 2`, not PNG.** The source is already lossy at 15 Mbit/s, so near-maximum-
quality JPEG adds negligible degradation relative to what the camera threw away, and it
turns 15 GB of PNG into 1.9 GB. Extracted frame counts match the container's `nb_frames`
exactly for all four clips (2369, 1674, 804, 933), which is the check that the extraction
dropped nothing.

## 2. The length unit

There is no ruler and no calibration target, the camera height is unknown, and it drifts.
Metric scale is unrecoverable and is not reported.

The unit is the fish body length, BL, measured per frame from silhouettes verified to be
isolated. Beyond removing units this cancels the drift: move the camera closer and school
and fish both grow in pixels by the same factor, so every ratio is preserved.

**Only isolated silhouettes are used.** A merged pair of overlapping fish has a longer
major axis than either fish. Because merging is worse in dense frames, admitting merged
clumps would make apparent body size correlate with school density, manufacturing a
relationship out of a segmentation artefact.

Isolation is enforced by dilating each candidate silhouette by 0.35 BL and requiring that
it touches no other. This is strict: typically only one to five fish per frame qualify, so
body length is pooled over a ±4 s rolling window rather than computed per frame. That is legitimate, because body length is a property of the
population rather than of a frame, and the window still tracks slow camera drift.

**A wrong assumption, corrected.** The first version also demanded solidity ≥ 0.88, on the
reasoning that a merged pair is lumpy while a single fish is convex. That rejected every
silhouette in every frame. A catfish fry is comma-shaped, with a broad head and a curved
tapering tail, so a *correctly* segmented single fish has solidity around 0.65. The
criterion was testing for the wrong shape. Solidity is now only a loose sanity floor, and
merging is detected directly by the clearance test, which is what actually distinguishes
one fish from two.

**A circular definition that ran away.** The first version set the silhouette acceptance
window as a fraction of the *assumed* body length squared. Because body length is what
was being measured, this was circular, and not benignly so. A low starting guess narrows
the window, which excludes genuinely large single fish, which lowers the estimate, which
narrows the window further. On IMG_0546 the pooled estimate fell from 21.9 px to 15.9 px
between the pre-scan and the batch pass, purely from this feedback.

Every criterion is now either an absolute noise floor (area ≥ 25 px, far below a fish) or
a scale-free shape ratio (aspect 2-6, eccentricity ≥ 0.80), and the isolation clearance is
scaled to each silhouette's *own* length rather than to a global assumption. The estimate
is now stable at ~14.5 px on IMG_0546 whether the iteration starts at 22, 32 or 60 px.

**Precision is good; absolute accuracy is not established.** Two different questions, and
conflating them would be the easy mistake.

*Precision.* Sweeping the background threshold from the 98th to the 99.9th percentile over
120 sampled frames per clip moves the estimate by only 2-3% (14.3-15.0 px on IMG_0546,
13.9-14.3 px on IMG_0547; see `results/bodylength_summary.json`). The measurement is
reproducible and largely insensitive to where the threshold is set, because the accepted
population is selected by the isolation test rather than by the threshold. An earlier
five-frame probe suggested ±15%; that was small-sample noise, and the larger sample
supersedes it.

*Accuracy.* Whether ~14.5 px equals the fish's true total length is a separate question
this imagery cannot settle. The quantity measured is the major axis of the thresholded
dark silhouette, and a faint tail below threshold would not be counted. An attempt to
quantify the shortfall by re-growing accepted silhouettes at a much lower threshold
failed: at a threshold loose enough to capture faint tails, silhouettes merge with their
neighbours, so the comparison measures clumping rather than tail length, and only two
usable pairs survived. That number is therefore not reported.

Projection is a second reason the value runs short. The fry occupy a water column and are
not all swimming level; a fish tilted relative to the camera projects shorter than it is.
That is visible in the spread: the accepted silhouettes run from about 12 to 18 px
between the quartiles against a median near 14.5. It can only bias the median downward,
never upward. Taking a high percentile instead of the median would be an equally
defensible convention and would give a larger unit; the median is used throughout, and
consistency, not the particular choice, is what the results depend on.

The honest description of the unit is: **the median major axis of an isolated, thresholded
fry silhouette**, a consistent, reproducible biological length scale, which is probably
somewhat shorter than a specimen measured in the hand. Settling the difference needs a
photograph of a fry against a ruler, which is a five-minute job on a future visit and
would retro-calibrate every number here.

Three consequences, and they are not equal:

- **Dimensionless quantities are unaffected.** Circularity, solidity, aspect ratio,
  boundary roughness, packing fraction, polarisation and rotation involve no length unit
  at all. These are the most trustworthy numbers in the project.
- **Ratios of two lengths are unaffected.** In particular ξ/D, the quantity the scale-free
  test actually turns on, has BL in both numerator and denominator and cancels exactly.
- **Absolute lengths in BL inherit any scale error.** School diameter and ξ in body
  lengths, and areas in BL², would all shift together if the unit is off.

Because the same BL divides every frame of a clip, an error in it is a *scale* error, not
noise: it multiplies values together and cannot create, destroy, or reverse a trend. Every comparison made in the results survives it intact: between clips, over time, between
rind and core, and ξ against school size.

**An indicative conversion, clearly flagged as such.** Published accounts put bullhead fry
at the guarded-school stage at roughly 20-30 mm total length. If that holds here, a school
measuring 28 BL across corresponds to something like 60-80 cm. That is an illustration for
a reader who wants a physical picture, derived from the literature and not from these
recordings, and no result in this project depends on it.

## 3. Segmenting the school

Three cues, each rejecting a confounder the others cannot.

- **Darkness.** The fish are dark on a pale bottom. Referenced to the 60th percentile of
  L\*, which sits above the school and inside the substrate, so it tracks illumination
  without being dragged down by the fish.
- **Not green.** Filamentous algae are also dark. In CIELAB they sit at strongly negative
  a\* while the fish sit near neutral, so subtracting a normalised greenness term
  separates them where lightness alone does not.
- **Texture.** The school casts a diffuse shadow on the lake bed which is dark *and*
  colour-neutral, so the first two cues cannot reject it. But the school is finely
  textured at fish scale, alternating silhouette and gap, while a shadow is smooth.
  Local standard deviation of L\* separates them cleanly.

The score is the product of a clipped darkness-and-hue term and a normalised texture term,
so a region must satisfy all three to score highly. Threshold by Otsu with hysteresis:
the Otsu level seeds the school, and the mask grows into connected material above half
that level, which recovers the sparse fringe without leaking into the shadow.

**The failure that produced this design.** The first version thresholded at three robust
standard deviations above the background median of a darkness field. It returned an empty
mask on a frame that plainly contains a school (IMG_0547, frame 559). The cause: that
frame's background is genuinely heterogeneous, bright sand, dark shadow, dark water, so
the MAD was measuring scene structure rather than noise, and three MADs above the median
landed *above the field maximum*. The general lesson is worth stating: a gate calibrated
purely on spread becomes mis-calibrated exactly where the scene is most structured. The
same frame now segments correctly.

Two explicit guards remain: a frame is declared school-free if the 99.9th percentile of
the score falls below an absolute floor (frames where the school swims out of shot must
report *absent*, not a threshold fitted to gravel), and a mask covering more than 55% of
the frame is refused rather than emitted.

## 4. Individual detection, and what it cannot do

Individuals are found as band-pass maxima inside the envelope. The band-pass scales matter:
an early attempt used a black top-hat with a structuring element of 2.5 BL, which by
construction can only find dark objects *smaller than the element and surrounded by
brighter background*. Inside a dense school there is no brighter background, so the
interior was invisible to it, the detector found only the rim. A difference of Gaussians
at 0.15 and 0.60 BL estimates the background at roughly one fish and does respond inside
the school.

**It is still incomplete, and this is not fixable.** Checked against hand inspection of
upscaled crops, recall is roughly half in moderately dense regions and worse in the core.
The reason is not the algorithm. The fry occupy a three-dimensional volume and the camera
records a single projection of it, so silhouettes genuinely superimpose, at 8×
magnification the "sparse" fringe still shows overlapping fish that no observer could
individuate with confidence.

The consequence is taken seriously rather than worked around:

- Density is reported as **packing fraction**, measured from occupied *area*, which needs
  no individuation and is therefore trustworthy.
- Nearest-neighbour distances are reported for the fringe only, and flagged as **upper
  bounds**, a missed neighbour can only make the nearest one look further away.
- Kinematics come from **dense optical flow**, which measures the motion of the pattern
  and does not require knowing which fish is which. This is the standard approach for
  dense active-matter systems whose constituents are unresolvable, not an improvisation.

## 5. Removing the camera

The phone was handheld, so apparent motion of the school is fish plus photographer, and
the photographer's contribution is coherent across the whole frame, precisely what an
order parameter is designed to detect.

Camera motion is estimated from Lucas-Kanade tracking of corner features on the lake bed
*outside* a generously dilated school mask, fitted with a RANSAC similarity transform
(translation, rotation, uniform scale). A full homography was rejected: the bottom is
roughly planar and roughly parallel to the sensor, so the extra projective freedom mostly
fits surface-ripple refraction noise.

**This step is load-bearing.** In IMG_0545 around frame 1000 the camera pans at about
34 px per frame. Uncorrected, that frame reports polarisation 0.996 and a median speed of
41 px per frame, a school in near-perfect collective alignment. With camera motion
removed, the same frame gives polarisation 0.279 and 6.4 px per frame. The apparent order
was entirely the photographer's arm. Any analysis of handheld footage that skips this
step will find spectacular collective order and it will be meaningless.

The estimator was validated against a synthetic translation and recovers an imposed
(+7.0, −4.0) px shift as (+7.01, −4.02).

## 6. Correlation analysis

Following the starling-flock method: subtract the school's mean velocity to obtain
fluctuations, then measure how far apart two points can be while their fluctuations still
align. The correlation length ξ is the first zero crossing of

    C(r) = ⟨δv_i · δv_j⟩ over pairs at separation r, divided by ⟨δv·δv⟩

Two things learned from validating the estimator on synthetic fields, both of which
changed the analysis:

**Per-frame ξ is too noisy to use.** On synthetic fields with a known correlation scale,
the single-realisation estimate has a spread comparable to its own value. C(r) is
therefore averaged across frames *first*, and ξ taken from the averaged curve. Averaging
per-frame ξ values would be dominated by noise in the tail of that distribution. For the
size-scaling test, frames are grouped into school-size bins and C(r) averaged within each
bin.

**A zero crossing exists even for pure noise.** Subtracting the sample mean constrains the
fluctuations to sum to zero over the sampled points, which forces C(r) negative somewhere
regardless of the underlying field. A correlation length is therefore not self-evidently
meaningful. Two safeguards: on synthetic white noise the spurious ξ stays small relative
to the domain (verified), and the real measurement is compared against a **permutation
null** in which velocity vectors are shuffled among the sampled positions. The shuffle
destroys spatial arrangement while preserving the velocity distribution, the point
geometry, the school outline, and the mean-subtraction constraint, so everything that
could produce a spurious ξ survives it, and only genuine spatial structure does not.

ξ is also checked against the largest separation actually sampled. A value approaching
that limit is a lower bound, not a measurement, and is flagged as saturated.

**There is also a floor, and it is not removed by the permutation null.** The Farnebäck
estimator uses a 25 px averaging window, which at the measured body length of roughly
21-22 px is about 1.2 BL. Neighbouring flow vectors therefore share input pixels and are
correlated by construction, so C(r) contains structure out to roughly one body length
that comes from the estimator rather than from the fish.

The permutation null does not control for this. Shuffling velocities among positions
destroys the flow-induced correlation along with everything else, so the null is a test
of spatial structure in general, not of structure beyond the smoothing scale.

The consequence is a hard interpretive rule, stated before the numbers were seen: a
correlation length of order 1 BL is indistinguishable from the measurement kernel and
means nothing. Only ξ substantially larger than ~1.2 BL is evidence of real coordination
between fish, and both numbers are reported together so the comparison is always visible.
A finer-grained flow estimate would lower this floor but would also become dominated by
the tail-beat motion of individual fish rather than the motion of the group.

## 7. Validation

`src/test_numerics.py` checks each routine against a case with a known answer: circularity
and aspect ratio of a digital disc and ellipse against their analytic values; the
body-length normalisation against its expected scaling; constructed splits and merges
against the fission/fusion detector; polarisation against perfectly aligned, isotropic,
and cone-distributed fields including the analytic sin(w)/w; the correlation length
against imposed smoothing scales; the zero crossing against a hand-computed
interpolation; camera motion against a synthetic shift; and the input guards against
out-of-range, NaN, and degenerate frames. All checks pass.

Two of them failed on first run and both failures were informative. One was a bug in the
test, not the code, the sign convention of `scipy.ndimage.shift` was backwards, and the
estimator had been right all along. The other was a genuinely underpowered test that
revealed the per-frame ξ noise problem described above, and changed the analysis.

## 8. Performance

The batch pass is ~2.4 s per frame per clip with four clips running concurrently, about
90 minutes for all 5780 frames on 22 cores.

Two optimisations were needed to get there from an initial 3 s per frame *single*-clip
rate. The CIELAB conversion and the large-σ Gaussians were being recomputed up to five
times per frame by separate functions; they are now computed once and shared, which alone
roughly halved the cost. Gaussian smoothing moved from scipy to OpenCV, 4.9× faster, after
confirming the two agree to 0.07 L\* units, 0.3% of the fish/background contrast.

One memory bug was caught before the long run: fission detection originally buffered every
frame's label image, which at 8 MB per frame would have needed roughly 19 GB on the
longest clip. It is now incremental against the previous frame only.
