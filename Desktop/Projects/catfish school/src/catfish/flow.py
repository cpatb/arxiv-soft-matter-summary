"""Velocity fields, camera-motion removal, and correlation analysis.

Individual fish cannot be resolved throughout these schools -- silhouettes
overlap in projection because the fish occupy a three-dimensional volume that
the camera sees as one plane. Rather than pretend otherwise, the kinematics are
measured as a *field*, by dense optical flow, exactly as collective motion is
measured in dense active-matter experiments where individual constituents are
likewise unresolvable.

Camera motion
-------------
The video is handheld. Apparent motion of the school therefore contains the
phone's motion, which must be removed before any velocity is believed. It is
estimated from feature correspondences on the lake bed *outside* the school and
applied as a similarity transform. Without this step a slow pan would appear as
coherent whole-school drift and would inflate the polarisation to near 1.

Correlation
-----------
Following Cavagna et al. (PNAS 2010) on starling flocks: subtract the school's
mean velocity to obtain fluctuations, then measure how far apart two points can
be while their fluctuations still align. The correlation length is where that
correlation first crosses zero. Reported in body lengths.
"""

from __future__ import annotations

from dataclasses import dataclass

import cv2
import numpy as np

#: Farneback parameters. The window must be large enough to average over several
#: fish (so the estimate is of the school's local motion, not one tail beat) yet
#: small compared with the school (so structure is not smeared away).
FARNEBACK = dict(
    pyr_scale=0.5,
    levels=4,
    winsize=25,
    iterations=3,
    poly_n=7,
    poly_sigma=1.5,
    flags=0,
)

#: Minimum number of RANSAC inliers for a camera-motion estimate to be trusted.
MIN_CAMERA_INLIERS = 25

#: Grid spacing for sampling the dense flow field, in body lengths. Sampling at
#: roughly one fish spacing keeps the correlation function from being dominated
#: by the optical-flow smoothing kernel.
SAMPLE_SPACING_BL = 1.0


@dataclass
class CameraMotion:
    """Similarity transform taking the previous frame onto the current one."""

    matrix: np.ndarray | None  # 2x3 affine, or None if estimation failed
    n_inliers: int
    ok: bool

    def velocity_at(self, ys: np.ndarray, xs: np.ndarray) -> np.ndarray:
        """Apparent (u, v) the camera alone would produce at these points."""
        if not self.ok or self.matrix is None:
            return np.zeros((len(ys), 2))
        pts = np.stack([xs, ys], axis=1).astype(np.float64)
        moved = pts @ self.matrix[:, :2].T + self.matrix[:, 2]
        return moved - pts  # (dx, dy)


def estimate_camera_motion(
    prev_gray: np.ndarray, cur_gray: np.ndarray, exclude: np.ndarray | None = None
) -> CameraMotion:
    """Estimate frame-to-frame camera motion from static background features.

    ``exclude`` masks out the school so that the fish -- which are moving
    independently -- cannot drive the estimate.
    """
    if prev_gray.shape != cur_gray.shape:
        raise ValueError("frame shape changed between frames")

    mask = None
    if exclude is not None:
        # Dilate generously: fish just outside the envelope would otherwise be
        # treated as static background.
        k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (81, 81))
        mask = (~cv2.dilate(exclude.astype(np.uint8), k).astype(bool)).astype(np.uint8) * 255

    pts = cv2.goodFeaturesToTrack(
        prev_gray, maxCorners=1200, qualityLevel=0.01, minDistance=12, mask=mask
    )
    if pts is None or len(pts) < MIN_CAMERA_INLIERS:
        return CameraMotion(None, 0, False)

    nxt, status, _ = cv2.calcOpticalFlowPyrLK(
        prev_gray,
        cur_gray,
        pts,
        None,
        winSize=(21, 21),
        maxLevel=3,
        criteria=(cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 30, 0.01),
    )
    if nxt is None:
        return CameraMotion(None, 0, False)

    good = status.ravel() == 1
    if good.sum() < MIN_CAMERA_INLIERS:
        return CameraMotion(None, int(good.sum()), False)

    src = pts[good].reshape(-1, 2)
    dst = nxt[good].reshape(-1, 2)

    # Similarity (translation + rotation + uniform scale) rather than a full
    # homography: the lake bed is roughly planar and roughly parallel to the
    # sensor, and the extra projective freedom mostly fits surface-ripple
    # refraction noise.
    M, inl = cv2.estimateAffinePartial2D(
        src, dst, method=cv2.RANSAC, ransacReprojThreshold=2.5, maxIters=3000
    )
    if M is None or inl is None or int(inl.sum()) < MIN_CAMERA_INLIERS:
        return CameraMotion(None, 0 if inl is None else int(inl.sum()), False)

    return CameraMotion(M.astype(np.float64), int(inl.sum()), True)


def dense_flow(prev_gray: np.ndarray, cur_gray: np.ndarray) -> np.ndarray:
    """Farneback dense optical flow. Returns (H, W, 2) array of (dx, dy) px/frame."""
    return cv2.calcOpticalFlowFarneback(prev_gray, cur_gray, None, **FARNEBACK)


def sample_school_velocity(
    flow: np.ndarray,
    mask: np.ndarray,
    body_px: float,
    camera: CameraMotion | None = None,
) -> dict[str, np.ndarray]:
    """Sample the flow field on a grid inside the school, camera motion removed.

    Returns arrays y, x (px), u, v (px/frame, fish motion only).
    """
    step = max(2, int(round(SAMPLE_SPACING_BL * body_px)))
    H, W = mask.shape
    ys, xs = np.mgrid[step // 2 : H : step, step // 2 : W : step]
    ys, xs = ys.ravel(), xs.ravel()
    keep = mask[ys, xs]
    ys, xs = ys[keep], xs[keep]
    if len(ys) == 0:
        return {k: np.array([]) for k in ("y", "x", "u", "v")}

    u = flow[ys, xs, 0].astype(np.float64)
    v = flow[ys, xs, 1].astype(np.float64)

    if camera is not None and camera.ok:
        cam = camera.velocity_at(ys, xs)
        u -= cam[:, 0]
        v -= cam[:, 1]

    return {"y": ys.astype(np.float64), "x": xs.astype(np.float64), "u": u, "v": v}


def order_parameters(u: np.ndarray, v: np.ndarray) -> dict:
    """Polarisation and speed statistics of a sampled velocity field.

    Polarisation is |<v_hat>|, the length of the mean unit velocity: 1 when every
    sample points the same way, ~0 when directions are uniformly spread. Points
    with speed below a small floor are dropped, because the direction of a
    near-zero vector is noise and would otherwise dilute the average.
    """
    if len(u) == 0:
        return {"polarisation": np.nan, "mean_speed_px": np.nan,
                "median_speed_px": np.nan, "n_samples": 0, "n_used": 0}

    speed = np.hypot(u, v)
    floor = max(0.05, 0.1 * float(np.median(speed)))
    ok = speed > floor
    if ok.sum() < 5:
        return {"polarisation": np.nan, "mean_speed_px": float(np.mean(speed)),
                "median_speed_px": float(np.median(speed)),
                "n_samples": len(u), "n_used": int(ok.sum())}

    ux, vy = u[ok] / speed[ok], v[ok] / speed[ok]
    pol = float(np.hypot(ux.mean(), vy.mean()))
    return {
        "polarisation": pol,
        "mean_speed_px": float(np.mean(speed)),
        "median_speed_px": float(np.median(speed)),
        "n_samples": int(len(u)),
        "n_used": int(ok.sum()),
    }


def correlation_function(
    y: np.ndarray,
    x: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    body_px: float,
    n_bins: int = 40,
    max_r_bl: float | None = None,
) -> dict:
    """Spatial correlation of velocity fluctuations, after Cavagna et al. (2010).

    The mean velocity of the school is subtracted first, so what remains is the
    fluctuation about collective motion. C(r) is the normalised dot product of
    fluctuations for pairs separated by r:

        C(r) = < dv_i . dv_j >_{|r_ij| = r}  /  < dv_i . dv_i >

    C(0) = 1 by construction. The correlation length xi is the first zero
    crossing: beyond it, fluctuations are on average anti-aligned, which is the
    signature of a finite-size collective mode rather than independent noise.

    Returns r in body lengths, C(r), pair counts per bin, and xi.
    """
    n = len(u)
    if n < 30:
        return {"r_bl": np.array([]), "C": np.array([]), "counts": np.array([]),
                "xi_bl": np.nan, "n_points": n}

    du = u - u.mean()
    dv = v - v.mean()
    norm = float(np.mean(du * du + dv * dv))
    if norm <= 0:
        return {"r_bl": np.array([]), "C": np.array([]), "counts": np.array([]),
                "xi_bl": np.nan, "n_points": n}

    # Pairwise separations and fluctuation dot products, upper triangle only.
    dy = y[:, None] - y[None, :]
    dx = x[:, None] - x[None, :]
    r = np.hypot(dy, dx) / body_px
    dot = du[:, None] * du[None, :] + dv[:, None] * dv[None, :]

    iu = np.triu_indices(n, k=1)
    r, dot = r[iu], dot[iu]

    if max_r_bl is None:
        # Beyond half the maximum separation the pair statistics are dominated
        # by the few widest-separated points and by the school's own shape.
        max_r_bl = float(np.percentile(r, 50))
    sel = r <= max_r_bl
    r, dot = r[sel], dot[sel]
    if len(r) < 30:
        return {"r_bl": np.array([]), "C": np.array([]), "counts": np.array([]),
                "xi_bl": np.nan, "n_points": n}

    edges = np.linspace(0.0, max_r_bl, n_bins + 1)
    idx = np.clip(np.digitize(r, edges) - 1, 0, n_bins - 1)
    counts = np.bincount(idx, minlength=n_bins)
    sums = np.bincount(idx, weights=dot, minlength=n_bins)

    with np.errstate(invalid="ignore", divide="ignore"):
        C = np.where(counts > 0, sums / np.maximum(counts, 1) / norm, np.nan)
    centres = 0.5 * (edges[:-1] + edges[1:])

    # Require a minimum number of pairs per bin before believing a value.
    good = counts >= 30
    centres, C, counts = centres[good], C[good], counts[good]

    return {
        "r_bl": centres,
        "C": C,
        "counts": counts,
        "xi_bl": _first_zero_crossing(centres, C),
        "n_points": n,
    }


def _first_crossing(r: np.ndarray, C: np.ndarray, level: float) -> float:
    """Linear interpolation of the first downward crossing of `level`.

    Returns NaN when C never reaches the level within the measured range, and
    the smallest measured separation when it is already below it at the start.
    Callers must distinguish those two cases -- they mean opposite things.
    """
    if len(r) < 2:
        return float("nan")
    finite = np.isfinite(C)
    r, C = r[finite], C[finite]
    if len(r) < 2:
        return float("nan")
    if C[0] <= level:
        # Already below at the smallest separation resolved: the correlation
        # length is shorter than the measurement can see, not undefined.
        return float(r[0])
    below = np.where(C < level)[0]
    if len(below) == 0:
        return float("nan")  # never reaches the level within range
    i = int(below[0])
    r0, r1, c0, c1 = r[i - 1], r[i], C[i - 1], C[i]
    if c0 == c1:
        return float(r1)
    return float(r0 + (r1 - r0) * (c0 - level) / (c0 - c1))


def _first_zero_crossing(r: np.ndarray, C: np.ndarray) -> float:
    """First downward zero crossing of C(r); NaN if it never crosses."""
    return _first_crossing(r, C, 0.0)


#: Level defining the correlation decay length, as a fraction of C(0) = 1.
DECAY_LEVEL = 1.0 / np.e


def correlation_decay_length(r: np.ndarray, C: np.ndarray) -> float:
    """Separation at which C(r) has fallen to 1/e of its value at zero.

    Preferred over the zero crossing for this dataset. The measured C(r) decays
    smoothly but stays positive out to the largest separation that can be
    sampled reliably, so a zero crossing frequently does not exist -- and where
    one does appear in a single noisy frame, it is set by noise rather than by
    the correlation. The 1/e level is reached well inside the measured range and
    is therefore a real measurement rather than an extrapolation.

    C(0) = 1 by construction, so the level is simply 1/e.
    """
    return _first_crossing(r, C, DECAY_LEVEL)
