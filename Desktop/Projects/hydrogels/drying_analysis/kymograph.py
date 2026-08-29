"""
Side-view kymograph and crack-vs-peel ordering (MOC section 5B).

The sharpest mechanism test in the program asks whether a vertical crack nucleates
at the free surface *first* and then deflects into a horizontal peel (Chung-
Mahadevan), or whether peel/blister comes first. The workhorse is a kymograph: a
single vertical pixel-line through a peel site, sampled across a high-speed frame
stack, producing a depth (z) vs time (t) image. Timestamps of the crack arrival
and the deflection are read from it.

This module:
  * builds the kymograph from a (T, H, W) frame stack given a column (or a small
    column band, averaged) -- `build_kymograph`;
  * extracts the crack-tip depth trajectory z(t) as the leading (shallowest)
    strongly-changed pixel per frame, using a frame-difference threshold --
    `crack_tip_trajectory`;
  * classifies ordering from the trajectory shape -- `classify_ordering`.

No physical model is imposed; these are geometric readouts of the profile. px/mm
and frame rate convert to physical depth and time downstream. Validated on a
synthetic descending-then-turning crack.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def build_kymograph(frames, col, band: int = 1) -> np.ndarray:
    """Depth-time image from a frame stack.

    Parameters
    ----------
    frames : (T, H, W) array (grayscale). If a color (T,H,W,C) stack is passed,
        it is averaged over channels.
    col : column index of the vertical sample line (depth axis = image rows).
    band : half-width; average over columns [col-band, col+band] inclusive
        (2*band+1 columns, clipped to the image). band=0 is a single column.

    Returns
    -------
    kymo : (H, T) array -- rows are depth z, columns are time.
    """
    f = np.asarray(frames)
    if f.ndim == 4:
        f = f.mean(axis=-1)
    if f.ndim != 3:
        raise ValueError(f"frames must be (T,H,W) or (T,H,W,C); got {f.shape}")
    if band < 0:
        raise ValueError(f"band must be >= 0; got {band}")
    T, H, W = f.shape
    if not (0 <= col < W):
        raise ValueError(f"col {col} out of range [0,{W})")
    lo = max(0, col - band)
    hi = min(W, col + band + 1)
    strip = f[:, :, lo:hi].mean(axis=2)   # (T, H)
    return strip.T                        # (H, T)


@dataclass
class TipTrajectory:
    t_index: np.ndarray       # frame indices where a tip was detected
    z_index: np.ndarray       # depth (row) of the leading changed pixel
    t_first: int              # first frame with any detection (onset)
    z_first: float            # depth at onset (surface if ~ surface_row)
    surface_row: int          # row treated as the free surface
    depth: np.ndarray = None  # |z_index - surface_row| = depth below surface


def crack_tip_trajectory(kymo, thresh=None, surface_row: int = 0,
                         going_down: bool = True, ref_col: int = 0) -> TipTrajectory:
    """Leading-edge depth per frame from frame-to-frame change.

    A pixel is 'changed' when |kymo(z,t) - kymo(z, ref_col)| exceeds `thresh`,
    referenced to column `ref_col` (default 0). The crack TIP is the changed
    pixel FARTHEST from the surface (the leading edge of the growing crack) each
    frame when going_down is True -- i.e. the deepest changed pixel below the
    surface. (`going_down=False` selects the shallowest, for a crack advancing
    toward the surface.)

    IMPORTANT: `ref_col` must be a PRE-EVENT frame (before any crack/peel is
    visible on the sample line) -- frame-differencing is blind to a feature that
    is already present in the reference frame. In practice, use the Phantom
    pre-trigger circular buffer so the stack starts before the event. The noise
    threshold is estimated from the reference column only, so a large event later
    in the stack does not inflate it.
    """
    k = np.asarray(kymo, dtype=float)
    if k.ndim != 2:
        raise ValueError(f"kymo must be 2-D (H,T); got {k.ndim}-D")
    H, T = k.shape
    if not (0 <= ref_col < T):
        raise ValueError(f"ref_col {ref_col} out of range [0,{T})")
    if not (0 <= surface_row < H):
        raise ValueError(f"surface_row {surface_row} out of range [0,{H})")
    ref = k[:, ref_col][:, None]
    change = np.abs(k - ref)
    if thresh is None:
        # noise floor from the (pre-event) reference column's frame-to-frame
        # scatter; robust MAD scaled to sigma. Constants named for traceability.
        NOISE_SIGMA = 5.0
        MAD_TO_SIGMA = 1.4826
        near_ref = change[:, max(0, ref_col - 1):ref_col + 2]
        mad = np.median(np.abs(near_ref - np.median(near_ref)))
        thresh = float(np.median(near_ref) + NOISE_SIGMA * MAD_TO_SIGMA * mad)
        if thresh <= 0:
            thresh = float(NOISE_SIGMA)  # fully flat reference: any real change
    changed = change > thresh

    t_idx, z_idx = [], []
    for t in range(T):
        rows = np.where(changed[:, t])[0]
        if rows.size == 0:
            continue
        # tip = leading edge = pixel FARTHEST from the surface (deepest) when the
        # crack grows downward; shallowest when it grows upward.
        dist = np.abs(rows - surface_row)
        tip = rows[np.argmax(dist)] if going_down else rows[np.argmin(dist)]
        t_idx.append(t)
        z_idx.append(int(tip))

    t_idx = np.asarray(t_idx, dtype=int)
    z_idx = np.asarray(z_idx, dtype=int)
    if t_idx.size == 0:
        return TipTrajectory(t_idx, z_idx, -1, np.nan, surface_row,
                             np.asarray([], dtype=float))
    depth = np.abs(z_idx - surface_row).astype(float)
    return TipTrajectory(t_idx, z_idx, int(t_idx[0]), float(z_idx[0]),
                         surface_row, depth)


def classify_ordering(traj: TipTrajectory, surface_tol_px: int = 3,
                      settle_frac: float = 0.8) -> str:
    """Classify a peel-site trajectory as one of:
       'crack->peel'  : tip starts at the surface and descends, then SETTLES
                        (reaches its max depth well before the end and holds) ->
                        a vertical crack that deflects into a horizontal peel.
       'peel-first'   : first detection is subsurface (not at the surface).
       'crack-only'   : descends and is still descending at the end (reaches max
                        depth only in the final frames) -> no deflection.
       'indeterminate': too few points.

    Classification is on depth-below-surface = |z - surface_row|, so it is
    independent of whether the surface is at the top or bottom of the frame.
    Kinematic only; it does not prove crack/peel continuity (that needs the full
    kymograph read).
    """
    if traj.t_index.size < 4:
        return "indeterminate"
    depth = traj.depth
    if depth is None or depth.size < 4:
        return "indeterminate"
    started_at_surface = abs(traj.z_first - traj.surface_row) <= surface_tol_px
    if not started_at_surface:
        return "peel-first"
    dmax = depth.max()
    if dmax <= surface_tol_px:
        # started at surface and never went anywhere -> not a descending crack
        return "indeterminate"
    # frame index (within the detected sequence) where depth first reaches
    # ~max: if that happens well before the end and depth then holds, the crack
    # stopped descending (deflected); if it only reaches max at the very end,
    # it was still propagating.
    i_settle = int(np.argmax(depth >= 0.95 * dmax))
    n = depth.size
    settled_early = i_settle <= settle_frac * (n - 1)
    tail_flat = (depth[i_settle:].max() - depth[i_settle:].min()) <= surface_tol_px
    if settled_early and tail_flat:
        return "crack->peel"
    return "crack-only"
