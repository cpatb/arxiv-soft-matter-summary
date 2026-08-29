"""The length unit: fish body length, measured from the fish themselves.

These recordings have no ruler, no calibration target, and an unknown camera
height that drifts as the phone is held over the water. No metric scale can be
recovered. What *can* be recovered is a biological scale -- the median fish body
length, BL -- measured in the same frame as everything else.

Normalising by a per-frame BL does more than remove units. If the camera drifts
closer, the school and the fish grow in pixels by the same factor, so every
ratio (school width in BL, nearest-neighbour distance in BL, speed in BL/s) is
invariant to that drift. The measurement is self-calibrating.

Only *isolated* silhouettes are used. A merged clump of two overlapping fish has
a longer major axis than either fish, so admitting clumps would bias BL upward
and, because clumping is worse in dense frames, would imprint a spurious
correlation between school density and apparent body size.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import measure, morphology

from . import segment

#: Absolute noise floor on silhouette area, in px. Fixed rather than scaled to
#: an assumed body length, and far below any fish, so it cannot influence the
#: measured value.
#:
#: An earlier version set the acceptance window as a fraction of the *assumed*
#: body length squared. That made the measurement circular in a way that did not
#: merely add noise but ran away: a low starting guess narrowed the window, which
#: excluded genuinely large single fish, which lowered the estimate further. On
#: IMG_0546 it drove the pooled estimate from 21.9 px down to 15.9 px. Every
#: criterion below is now either an absolute noise floor or a scale-free shape
#: ratio, so the result cannot depend on where the iteration started.
MIN_ABS_AREA_PX = 25.0

#: Upper bound on area, in px, purely to exclude the whole merged school (which
#: reaches ~3e5 px). Set far above a plausible single fish so it never binds on
#: real individuals; merged *pairs* are rejected by the clearance test, not here.
MAX_ABS_AREA_PX = 4000.0

#: Fish-like elongation, major/minor. A fry seen from above is a slender comma.
#: Scale-free, so it does not reintroduce the circularity.
ASPECT_RANGE = (2.0, 6.0)

#: Loose sanity floor on solidity only. It is tempting to demand a high value on
#: the reasoning that a merged pair is lumpy, but that is wrong here: a catfish
#: fry is comma-shaped, with a broad head and a curved tapering tail, so a
#: *correctly* segmented single fish has solidity ~0.65. Requiring 0.88 rejected
#: every silhouette in every frame. Merging is instead detected directly, by the
#: isolation-clearance test below, which is what actually distinguishes a lone
#: fish from a touching pair.
ISOLATED_MIN_SOLIDITY = 0.45

#: Fish are elongated; round blobs are debris, gaps, or head-on fish whose
#: apparent length is foreshortened and therefore not a body length.
ISOLATED_MIN_ECCENTRICITY = 0.80

#: An isolated fish must have no other silhouette within this distance, in
#: nominal body lengths, of its own boundary.
ISOLATION_CLEARANCE_BL = 0.35

#: Minimum number of accepted silhouettes for a frame's BL estimate to be used.
MIN_ISOLATED_PER_FRAME = 8


@dataclass(frozen=True)
class BodyLength:
    """Measured body length.

    ``body_px`` is the median major axis of the fitted ellipse, i.e. a chord.
    For a fish whose tail is curved this slightly underestimates the true
    along-body length; for the curvatures seen here the shortfall is a few per
    cent, and it is common to every frame, so ratios are unaffected.
    """

    body_px: float
    n_isolated: int
    iqr_px: tuple[float, float]
    per_fish: pd.DataFrame


def isolated_silhouettes(
    rgb: np.ndarray | None = None,
    body_px: float = segment.NOMINAL_BODY_PX,
    fields: segment.Fields | None = None,
    env: segment.Envelope | None = None,
) -> pd.DataFrame:
    """Find fish silhouettes that do not touch or nearly touch another.

    Pass ``fields`` and ``env`` when they are already available for this frame;
    recomputing them here was doubling the per-frame cost of the batch run.

    Returns a DataFrame with major_px, minor_px, area_px, solidity,
    eccentricity, y, x for each accepted silhouette.
    """
    if fields is None:
        if rgb is None:
            raise ValueError("provide either rgb or fields")
        fields = segment.compute_fields(rgb, body_px)
    body_px = fields.body_px
    if env is None:
        env = segment.school_envelope(fields=fields)

    cols = ["y", "x", "area_px", "major_px", "minor_px", "solidity", "eccentricity"]
    if not env.present:
        return pd.DataFrame(columns=cols)

    # Threshold the raw (unsmoothed) fishness against the off-school
    # distribution: this cuts silhouettes at their true edges rather than at a
    # smoothed contour, which matters when the quantity being measured is a
    # length.
    f = fields.fishness
    outside = ~env.mask
    if outside.sum() < 10_000:
        return pd.DataFrame(columns=cols)
    thr = float(np.percentile(f[outside], 99.5))

    # Search a collar around the school: the fringe, where background is
    # visible between fish. Inside the dense core silhouettes always touch.
    grown = ndi.binary_dilation(env.mask, iterations=int(round(1.0 * body_px)))  # search collar only
    fg = (f > thr) & grown

    fg = morphology.remove_small_objects(fg, max_size=int(MIN_ABS_AREA_PX) - 1)
    if not fg.any():
        return pd.DataFrame(columns=cols)

    labels, n = ndi.label(fg)
    if n == 0:
        return pd.DataFrame(columns=cols)

    # Isolation test: dilate each silhouette by a clearance proportional to its
    # own length and require that it does not meet any other label.
    rows = []
    for p in measure.regionprops(labels):
        if not (MIN_ABS_AREA_PX <= p.area <= MAX_ABS_AREA_PX):
            continue
        if p.axis_minor_length <= 0:
            continue
        aspect = p.axis_major_length / p.axis_minor_length
        if not (ASPECT_RANGE[0] <= aspect <= ASPECT_RANGE[1]):
            continue
        if p.solidity < ISOLATED_MIN_SOLIDITY:
            continue
        if p.eccentricity < ISOLATED_MIN_ECCENTRICITY:
            continue

        # Clearance scaled to this silhouette's own length, not to a global
        # assumed body length -- the whole point is that the criterion must not
        # know the answer in advance.
        clearance = max(2, int(round(ISOLATION_CLEARANCE_BL * p.axis_major_length)))

        # Local isolation: within a padded bounding box, dilate this label and
        # check whether the dilation touches any other label.
        r0, c0, r1, c1 = p.bbox
        pad = clearance + 2
        R0, C0 = max(0, r0 - pad), max(0, c0 - pad)
        R1, C1 = min(labels.shape[0], r1 + pad), min(labels.shape[1], c1 + pad)
        sub = labels[R0:R1, C0:C1]
        me = sub == p.label
        me_d = ndi.binary_dilation(me, iterations=clearance)
        others = (sub > 0) & (sub != p.label)
        if (me_d & others).any():
            continue

        y, x = p.centroid
        rows.append(
            {
                "y": float(y),
                "x": float(x),
                "area_px": float(p.area),
                "major_px": float(p.axis_major_length),
                "minor_px": float(p.axis_minor_length),
                "solidity": float(p.solidity),
                "eccentricity": float(p.eccentricity),
            }
        )

    return pd.DataFrame(rows, columns=cols)


def body_length(rgb: np.ndarray, body_px: float = segment.NOMINAL_BODY_PX) -> BodyLength:
    """Median major-axis length of isolated silhouettes, in pixels."""
    df = isolated_silhouettes(rgb, body_px=body_px)
    if len(df) < MIN_ISOLATED_PER_FRAME:
        return BodyLength(float("nan"), len(df), (float("nan"), float("nan")), df)
    v = df["major_px"].to_numpy()
    return BodyLength(
        body_px=float(np.median(v)),
        n_isolated=len(df),
        iqr_px=(float(np.percentile(v, 25)), float(np.percentile(v, 75))),
        per_fish=df,
    )


def converged_body_length(
    rgb: np.ndarray, start_px: float = segment.NOMINAL_BODY_PX, n_iter: int = 3
) -> BodyLength:
    """Iterate the measurement so it does not depend on the initial guess.

    The filter scales are set in body lengths, and body length is what is being
    measured, so the estimate is nominally circular. Re-running with the measured
    value breaks the circularity; in practice it converges in two iterations.
    """
    est = start_px
    out = body_length(rgb, body_px=est)
    for _ in range(n_iter - 1):
        if not np.isfinite(out.body_px):
            break
        if abs(out.body_px - est) / est < 0.02:  # converged to within 2%
            break
        est = 0.5 * (est + out.body_px)  # damped update, avoids oscillation
        out = body_length(rgb, body_px=est)
    return out
