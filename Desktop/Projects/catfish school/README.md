# Catfish fry-ball schooling, Montrose Point, Chicago, 8 August 2026

Analysis of four handheld iPhone videos of dense schools of bullhead catfish fry
(*Ameiurus* sp.) filmed from the J-hook of Montrose Pier, over the shallow water inside
the hook by the bird sanctuary.

## The recordings

| Clip | Duration | Frames | Orientation | What it shows |
|---|---|---|---|---|
| IMG_0544 | 79.0 s | 2369 | portrait | Two separate schools in frame together |
| IMG_0545 | 55.8 s | 1674 | landscape | One large crescent-shaped school; school leaves frame partway |
| IMG_0546 | 26.8 s | 804 | portrait | A school in the middle of splitting into two lobes |
| IMG_0547 | 31.1 s | 933 | portrait | Compact ellipse; the school's own shadow visible on the lake bed |

All four: 1920×1080, 30000/1001 fps (29.97), H.264, iPhone 11, recorded 12:55-12:59
local time. Three carry a −90° display rotation; one does not.

## There is no metric scale, by construction

No ruler, no calibration target, unknown and drifting camera height. Metric lengths and
speeds are **not recoverable** and are not reported anywhere in this project.

Instead the unit is the **body length (BL)**, measured from the fish themselves in the
same frame, from silhouettes verified to be isolated rather than merged. This is not
merely a substitute for centimetres; it is better behaved for this data. When the phone
drifts closer, the school and the fish both grow in pixels by the same factor, so every
ratio (school width in BL, speed in BL/s, correlation length in BL) is invariant to the
drift that cannot be measured.

Precisely, BL is the major axis of an isolated, thresholded fry silhouette. It is
*reproducible*: varying the segmentation threshold moves it by only 2-3%. Whether it
equals a fry's total length in the hand is not established here, since a faint tail below
threshold would not be counted. Any such error is a scale factor common to a whole clip,
so it multiplies values together and cannot create or reverse a trend. Dimensionless
quantities (circularity, packing fraction, polarisation) and ratios of two lengths
(notably ξ/D) are immune to it entirely.

## Layout

```
IMG_054*.MOV              originals, untouched
data/frames/<clip>/       extracted frames, rotation baked in (JPEG q=2, 1.9 GB)
data/meta/                data-inspection record
src/catfish/              analysis package
  frames.py               frame access, per-clip metadata
  segment.py              school envelope + individual detection
  calibrate.py            body length from isolated silhouettes
  topology.py             shape descriptors, fission/fusion
  flow.py                 camera-motion removal, velocity fields, C(r)
  plotstyle.py            minimal figure style
src/run_pipeline.py       batch pass over every frame -> results/*.csv
src/analyse_results.py    normalise to BL, make figures, summary tables
src/fringe_analysis.py    individual tracking at the fringe; speed cross-check
src/make_visual_figures.py image figures (gallery, fission, velocity field)
src/test_numerics.py      validation against analytic and synthetic cases
results/                  CSV/JSON output
figures/                  figures; FIGURES.md describes each
docs/LITERATURE_REVIEW.md what is known, and where the gap is
docs/RESULTS.md           findings
```

## Reproducing

```bash
# 1. extract frames from the originals (ffmpeg; bakes in display rotation)
#    see docs/METHODS.md for the exact command

# 2. validate the numerics first -- all checks must pass
python src/test_numerics.py

# 3. batch pass (~2.5 s/frame/clip; run the four clips in parallel)
python src/run_pipeline.py --clips IMG_0544

# 4. normalise, analyse, plot
python src/analyse_results.py
python src/fringe_analysis.py
python src/make_visual_figures.py
```

## What this analysis can and cannot do

**Measured reliably.** School shape and topology; fission and fusion; packing fraction
and its rind-to-core profile; velocity fields by optical flow; polarisation; the
velocity-fluctuation correlation function and its correlation length.

**Measured with stated caveats.** Individual positions and spacing at the sparse fringe
only. Detection recall there is roughly half, and missed neighbours can only inflate a
nearest-neighbour distance, so those values are upper bounds.

**Not possible.** Complete individual tracking. The fry occupy a three-dimensional
volume and the camera records one projection of it, so silhouettes genuinely overlap. In
the dense core they cannot be separated by eye at any magnification either. This is an
information limit of the recording, not a deficiency of the detector, and no amount of
additional processing removes it.

## Software

Python 3.12, numpy 2.4.2, scipy 1.17.1, scikit-image 0.26.0, opencv 4.13.0,
pandas 2.3.3, matplotlib 3.10.8, trackpy 0.7. Random seeds are fixed and recorded where
any sampling occurs.
