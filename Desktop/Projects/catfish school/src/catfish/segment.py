"""Segmentation of the school envelope and of individual fish.

Two products, deliberately separate because they need different methods and have
different reliability:

1. ``school_envelope`` -- the region occupied by the school, including the dense
   interior where silhouettes merge. Used for topology (area, shape, fission).

2. ``detect_individuals`` -- positions of individual fish. Reliable in sparse
   regions, incomplete in the dense core where fish overlap. Used for body
   length, nearest-neighbour statistics, and tracking. Its incompleteness is
   measured rather than assumed.

Discriminator
-------------
Three cues are needed, each rejecting a confounder the others cannot:

* **darkness** -- the fish are dark against a pale substrate.
* **not-green** -- filamentous algae are also dark, but strongly green. In
  CIELAB the algae sit at markedly negative a* and the fish near a* = 0.
* **texture** -- the school casts a diffuse shadow on the lake bed that is both
  dark and colour-neutral, so darkness and hue alone cannot reject it. The
  school is finely textured at fish scale (alternating silhouette and gap); a
  shadow is smooth. Local standard deviation of L* separates them.

The school score is the product of a clipped fishness term and a normalised
texture term, so a region must satisfy all three cues to score highly.

Threshold
---------
Otsu on the score, applied with hysteresis: the Otsu level seeds the school and
the mask is then grown into connected material above half that level, which
recovers the sparse fringe of individuals without leaking into the shadow.

An earlier version thresholded at a fixed number of robust standard deviations
above the background median. That failed on frames whose background is genuinely
heterogeneous (bright sand plus dark shadow plus dark water): there the MAD
measures scene structure rather than noise, and the threshold rose above the
field maximum, silently returning an empty mask on a frame that plainly contains
a school. The lesson generalises -- a gate calibrated purely on spread becomes
mis-calibrated exactly where the scene is most structured.
"""

from __future__ import annotations

from dataclasses import dataclass

import cv2
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import color, filters, morphology

# --- Parameters. Named constants, not magic numbers buried in expressions. -----

#: Nominal fish length in px at these camera distances, from visual inspection of
#: upscaled crops and confirmed by the measured distribution (~30 px). Only used
#: to set filter scales; results are normalised by the *measured* body length.
NOMINAL_BODY_PX = 32.0

#: Smoothing scale for both the fishness and texture fields, in body lengths.
FIELD_SMOOTH_BL = 0.5

#: Hysteresis low threshold as a fraction of the Otsu (high) threshold.
HYSTERESIS_RATIO = 0.5

#: A frame is declared to contain no school if the 99.9th percentile of the
#: score falls below this. Frames with a school reach ~0.4-0.5; frames without
#: one sit near 0, so the floor is far from both.
SCORE_PRESENCE_FLOOR = 0.05

#: Morphological closing radius for the envelope, in body lengths.
ENVELOPE_CLOSE_BL = 0.8

#: Discard envelope components smaller than this, in units of BL^2.
MIN_SCHOOL_AREA_BL2 = 12.0

#: A school occupying more than this fraction of the frame indicates the
#: threshold has collapsed onto background rather than fish.
MAX_SCHOOL_AREA_FRAC = 0.55

#: Band-pass scales for individual detection, in body lengths. The small scale
#: suppresses sensor/compression noise; the large scale estimates the local
#: background. Crucially the large scale is ~1 body length, not several: a wide
#: structuring element has no bright background to work against inside a dense
#: school and so cannot detect its interior at all.
DOG_SMALL_BL = 0.15
DOG_LARGE_BL = 0.60

#: Minimum separation between detected individuals, in body lengths.
PEAK_MIN_DIST_BL = 0.40

#: Accepted silhouette area range for an individual, in BL^2.
FISH_AREA_BL2 = (0.04, 0.60)

#: Minimum eccentricity for a detection to count as a clean, isolated single
#: fish when measuring body length. Merged pairs are rounder or lumpier.
FISH_MIN_ECCENTRICITY = 0.75


@dataclass(frozen=True)
class Envelope:
    """School occupancy mask and its connected components."""

    mask: np.ndarray  # bool (H, W)
    labels: np.ndarray  # int32 (H, W), 0 = background
    n_components: int
    score: np.ndarray  # the school score the mask was cut from
    present: bool  # False when no school is in frame


@dataclass(frozen=True)
class Fields:
    """Per-frame colour and score fields, computed once and reused.

    The CIELAB conversion and the large-sigma Gaussians together dominate the
    cost of this pipeline. Recomputing them inside each of fishness_field,
    texture_field, school_score and dog_response made a frame take ~3 s, which
    is 5 hours across this dataset; sharing them brings it under a second.
    """

    L: np.ndarray
    a: np.ndarray
    fishness: np.ndarray
    texture: np.ndarray
    score: np.ndarray
    body_px: float


def _blur(img: np.ndarray, sigma: float) -> np.ndarray:
    """Gaussian blur via OpenCV, which is markedly faster than scipy here."""
    if sigma <= 0:
        return img.astype(np.float64, copy=True)
    k = int(2 * round(3.0 * sigma) + 1)  # odd kernel covering +/- 3 sigma
    return cv2.GaussianBlur(
        img.astype(np.float32), (k, k), sigmaX=sigma, sigmaY=sigma,
        borderType=cv2.BORDER_REPLICATE,
    ).astype(np.float64)


def compute_fields(rgb: np.ndarray, body_px: float = NOMINAL_BODY_PX) -> Fields:
    """Compute every per-frame field in one pass."""
    _check_rgb(rgb)
    lab = color.rgb2lab(rgb)
    L, a = lab[..., 0], lab[..., 1]

    L_ref = float(np.percentile(L, 60))
    darkness = (L_ref - L) / max(L_ref, 1e-9)
    greenness = np.clip(-a, 0.0, None) / 20.0
    fish = darkness - greenness

    s = FIELD_SMOOTH_BL * body_px
    mu = _blur(L, s)
    mu2 = _blur(L * L, s)
    tex = np.sqrt(np.clip(mu2 - mu * mu, 0.0, None))

    fish_s = _blur(np.clip(fish, 0.0, None), s)
    tex_s = _blur(tex, s)
    tex_n = tex_s / max(float(np.percentile(tex_s, 99)), 1e-9)

    return Fields(L=L, a=a, fishness=fish, texture=tex, score=fish_s * tex_n,
                  body_px=body_px)


def _check_rgb(rgb: np.ndarray) -> None:
    if rgb.ndim != 3 or rgb.shape[2] != 3:
        raise ValueError(f"expected (H, W, 3) RGB, got {rgb.shape}")
    if not np.all(np.isfinite(rgb)):
        raise ValueError("non-finite values in input image")
    if rgb.min() < 0.0 or rgb.max() > 1.0:
        raise ValueError(f"expected RGB in [0, 1], got [{rgb.min()}, {rgb.max()}]")


# The three fields below are thin wrappers over compute_fields. They exist for
# readability at call sites and in tests; the formulae live in exactly one place
# so the batch path and the interactive path can never drift apart.


def fishness_field(rgb: np.ndarray) -> np.ndarray:
    """Dark-and-not-green score. Higher = more fish-like.

    Independent of body_px, so the default scale is harmless here.
    """
    return compute_fields(rgb, NOMINAL_BODY_PX).fishness


def texture_field(rgb: np.ndarray, body_px: float = NOMINAL_BODY_PX) -> np.ndarray:
    """Local standard deviation of L* at half a body length."""
    return compute_fields(rgb, body_px).texture


def school_score(rgb: np.ndarray, body_px: float = NOMINAL_BODY_PX) -> np.ndarray:
    """Combined dark + neutral + textured score. Zero or near-zero off-school."""
    return compute_fields(rgb, body_px).score


def school_envelope(
    rgb: np.ndarray | None = None,
    body_px: float = NOMINAL_BODY_PX,
    fields: Fields | None = None,
) -> Envelope:
    """Segment the region occupied by the school.

    Pass ``fields`` when they have already been computed for this frame; that
    avoids repeating the CIELAB conversion and the smoothing, which is most of
    the per-frame cost.
    """
    if fields is None:
        if rgb is None:
            raise ValueError("provide either rgb or fields")
        fields = compute_fields(rgb, body_px)
    score = fields.score
    empty = np.zeros(score.shape, dtype=bool)

    # Presence gate on an absolute level, not on a ratio or a spread: frames
    # where the school has swum out of shot must return "absent", not a
    # threshold fitted to background texture.
    if float(np.percentile(score, 99.9)) < SCORE_PRESENCE_FLOOR:
        return Envelope(empty, empty.astype(np.int32), 0, score, present=False)

    hi = float(filters.threshold_otsu(score))
    lo = HYSTERESIS_RATIO * hi
    mask = filters.apply_hysteresis_threshold(score, lo, hi)

    if mask.mean() > MAX_SCHOOL_AREA_FRAC:
        # Threshold has collapsed onto background; refuse rather than emit a
        # plausible-looking but wrong mask.
        return Envelope(empty, empty.astype(np.int32), 0, score, present=False)

    body_px = fields.body_px
    r = max(1, int(round(ENVELOPE_CLOSE_BL * body_px)))
    se = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r + 1, 2 * r + 1))
    m = cv2.morphologyEx(mask.astype(np.uint8), cv2.MORPH_CLOSE, se)
    m = ndi.binary_fill_holes(m.astype(bool))

    # scikit-image >= 0.26 deprecates min_size in favour of max_size, which
    # removes objects smaller than *or equal to* its value; hence the -1.
    min_area = int(MIN_SCHOOL_AREA_BL2 * body_px**2)
    m = morphology.remove_small_objects(m, max_size=max(min_area - 1, 0))

    if not m.any():
        return Envelope(empty, empty.astype(np.int32), 0, score, present=False)

    labels, n = ndi.label(m)
    return Envelope(m, labels.astype(np.int32), int(n), score, present=True)


def dog_response(
    rgb: np.ndarray | None = None,
    body_px: float = NOMINAL_BODY_PX,
    fields: Fields | None = None,
) -> np.ndarray:
    """Band-pass response to dark, fish-scale features.

    Difference of Gaussians on inverted lightness. Unlike a large-structuring-
    element black top-hat, this responds inside a dense school, because its
    background estimate is taken at ~1 body length rather than several.
    """
    if fields is None:
        if rgb is None:
            raise ValueError("provide either rgb or fields")
        fields = compute_fields(rgb, body_px)
    body_px = fields.body_px
    dark = -fields.L  # fish become bright maxima
    small = _blur(dark, DOG_SMALL_BL * body_px)
    large = _blur(dark, DOG_LARGE_BL * body_px)
    return small - large


def detect_individuals(
    rgb: np.ndarray | None = None,
    body_px: float = NOMINAL_BODY_PX,
    envelope: Envelope | None = None,
    fields: Fields | None = None,
) -> pd.DataFrame:
    """Detect individual fish as band-pass maxima within the school envelope.

    Returns
    -------
    DataFrame with columns y, x (centroid px), response, in_core (bool).
    Empty DataFrame with the same columns when no school is present.
    """
    cols = ["y", "x", "response", "in_core"]
    if fields is None:
        if rgb is None:
            raise ValueError("provide either rgb or fields")
        fields = compute_fields(rgb, body_px)
    body_px = fields.body_px
    if envelope is None:
        envelope = school_envelope(fields=fields)
    if not envelope.present:
        return pd.DataFrame(columns=cols)

    dog = dog_response(fields=fields)

    # Search slightly outside the envelope so fringe fish are not clipped.
    r = max(1, int(round(0.5 * body_px)))
    se = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r + 1, 2 * r + 1))
    search = cv2.dilate(envelope.mask.astype(np.uint8), se).astype(bool)

    # Null model from outside the school: there the band-pass sees only substrate
    # texture, so its upper tail sets a threshold that admits fish and rejects
    # gravel. Guard against frames where the school fills most of the view.
    outside = ~search
    ref = dog[outside] if outside.sum() > 10_000 else dog
    med = float(np.median(ref))
    mad = float(np.median(np.abs(ref - med)))
    robust_sd = 1.4826 * mad
    thresh = med + 3.0 * robust_sd if robust_sd > 0 else float(np.percentile(dog, 99))

    min_dist = max(2, int(round(PEAK_MIN_DIST_BL * body_px)))
    from skimage import feature  # local import keeps module import cheap

    coords = feature.peak_local_max(
        dog,
        min_distance=min_dist,
        threshold_abs=thresh,
        labels=search.astype(np.int32),
        exclude_border=False,
    )
    if len(coords) == 0:
        return pd.DataFrame(columns=cols)

    df = pd.DataFrame({"y": coords[:, 0].astype(float), "x": coords[:, 1].astype(float)})
    df["response"] = dog[coords[:, 0], coords[:, 1]]

    # Flag detections in the dense core (envelope eroded by 1.5 BL) so edge and
    # core statistics can be reported separately -- their reliability differs.
    r_core = max(1, int(round(1.5 * body_px)))
    se_core = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r_core + 1, 2 * r_core + 1))
    core = cv2.erode(envelope.mask.astype(np.uint8), se_core).astype(bool)
    df["in_core"] = core[coords[:, 0], coords[:, 1]]
    return df[cols]


def packing_fraction(
    rgb: np.ndarray | None, envelope: Envelope, fields: Fields | None = None
) -> float:
    """Fraction of the envelope's *projected area* covered by fish, not gap.

    This is a two-dimensional occupancy, not a three-dimensional volume
    fraction, and the distinction matters. The fry occupy a water column that
    the camera flattens into one plane, so a pixel counts as occupied if *any*
    fish lies above it at any depth. Adding more fish at greater depth therefore
    drives this quantity towards 1 without the school becoming any more crowded
    in three dimensions.

    Read it as optical density -- how completely the school blocks the view of
    the bottom beneath it -- and not as how tightly packed the fish are.
    Comparisons between frames and between clips remain meaningful; the absolute
    value is an upper bound on any 3D packing interpretation.
    """
    if not envelope.present or not envelope.mask.any():
        return float("nan")
    f = fields.fishness if fields is not None else fishness_field(rgb)
    outside = ~envelope.mask
    if outside.sum() < 1000:
        return float("nan")
    ref = float(np.percentile(f[outside], 99))
    return float((f[envelope.mask] > ref).mean())
