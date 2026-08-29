"""Motion-compensated temporal median for colorization chroma.

WHY THIS EXISTS
---------------
DDColor predicts colour per frame with no temporal state, so it flips between plausible
colour hypotheses from frame to frame (blue shirt / brown shirt). That flicker is
*multimodal*, which is why a temporal MEDIAN is the right statistic and a temporal MEAN
is actively wrong -- a mean of {blue, blue, brown, blue, brown} lands between the modes
and is wrong on every frame, while a median selects a mode and tolerates up to half the
window being wrong.

But a naive per-pixel temporal median asks "what colour was this PIXEL LOCATION over the
last N frames?" On moving content that is the wrong question. The right one is "what
colour was this SURFACE POINT?" The fix is to warp neighbouring frames' chroma into the
current frame's coordinate system with optical flow before taking the median.

This is a recombination of two standard primitives (flow-warped temporal consistency from
video stylization, and median filtering), not a novel method.

WHY IT IS CHEAP HERE
--------------------
DDColor computes chroma at a fixed 512x512 regardless of source resolution, and the
luminance channel never passes through the network -- full-resolution L is concatenated
back at the end. Two consequences:

  1. All temporal work happens on (2, 512, 512) arrays, not (3, 3840, 2160). The flow
     field only needs to be computed at 512x512, because the chroma carries no
     information above that scale. This turns flow-warping from expensive into nearly
     free, and it can run on the CPU alongside GPU inference.
  2. No amount of chroma smoothing can soften an edge, because sharpness lives entirely
     in L and L is untouched. The usual "temporal filtering costs you sharpness"
     objection does not apply.

Requires: numpy, opencv-python. No GPU, no torch.
"""

from __future__ import annotations

import numpy as np

try:
    import cv2
except ImportError as exc:  # pragma: no cover
    raise ImportError("temporal_chroma requires opencv-python (cv2)") from exc


# Forward-backward round-trip error above this many pixels marks a pixel as occluded
# or as a flow failure. 1.0 px is conservative; raise it if too much of the frame is
# being rejected on noisy footage, lower it if colour drags across occlusion edges.
DEFAULT_FB_THRESHOLD_PX = 1.0

# Odd window length for the temporal median. 5 is a reasonable start. Longer windows
# buy more stability and cost more lag wherever the flow fails. This is a taste knob.
DEFAULT_WINDOW = 5


def make_dis_flow(preset: str = "medium") -> "cv2.DISOpticalFlow":
    """Build a DIS optical flow estimator.

    DIS is the right default here: it is CPU-side, markedly better than Farneback, and
    cheap enough at 512x512 to run off the GPU critical path. RAFT is deliberately not
    used -- it needs the GPU VRAM that DDColor and Real-ESRGAN are already competing for.
    """
    presets = {
        "ultrafast": cv2.DISOPTICAL_FLOW_PRESET_ULTRAFAST,
        "fast": cv2.DISOPTICAL_FLOW_PRESET_FAST,
        "medium": cv2.DISOPTICAL_FLOW_PRESET_MEDIUM,
    }
    if preset not in presets:
        raise ValueError(f"preset must be one of {sorted(presets)}, got {preset!r}")
    return cv2.DISOpticalFlow.create(presets[preset])


def _as_flow_input(gray: np.ndarray) -> np.ndarray:
    """Coerce a luminance frame to the uint8 single-channel form DIS expects."""
    if gray.ndim != 2:
        raise ValueError(f"expected a 2-D luminance array, got shape {gray.shape}")
    if gray.dtype == np.uint8:
        return gray
    # Float input is assumed to be in [0, 1]; anything else is scaled by its own range
    # so that a 0-100 CIE L channel also works.
    finite = np.isfinite(gray)
    if not finite.all():
        gray = np.nan_to_num(gray, nan=0.0, posinf=0.0, neginf=0.0)
    lo = float(gray.min())
    hi = float(gray.max())
    if hi <= 1.0 and lo >= 0.0:
        scaled = gray * 255.0
    elif hi > lo:
        scaled = (gray - lo) * (255.0 / (hi - lo))
    else:
        scaled = np.zeros_like(gray)
    return np.clip(scaled, 0, 255).astype(np.uint8)


def warp_by_flow(image: np.ndarray, flow: np.ndarray) -> np.ndarray:
    """Sample ``image`` at the locations ``flow`` points to.

    ``flow`` must be the flow FROM the reference frame TO the frame that ``image``
    belongs to, so that ``flow[y, x]`` is the displacement of the surface point at
    (x, y) in the reference frame. The result is ``image`` resampled into the
    reference frame's coordinate system.
    """
    if flow.ndim != 3 or flow.shape[2] != 2:
        raise ValueError(f"flow must be (H, W, 2), got {flow.shape}")
    height, width = flow.shape[:2]
    if image.shape[:2] != (height, width):
        raise ValueError(f"image {image.shape[:2]} does not match flow {(height, width)}")

    grid_x, grid_y = np.meshgrid(
        np.arange(width, dtype=np.float32), np.arange(height, dtype=np.float32)
    )
    map_x = grid_x + flow[..., 0]
    map_y = grid_y + flow[..., 1]
    return cv2.remap(
        image, map_x, map_y,
        interpolation=cv2.INTER_LINEAR,
        borderMode=cv2.BORDER_REPLICATE,
    )


def fb_consistency_mask(
    flow_fw: np.ndarray, flow_bw: np.ndarray,
    threshold_px: float = DEFAULT_FB_THRESHOLD_PX,
) -> np.ndarray:
    """Boolean mask of pixels whose flow survives a forward-backward round trip.

    Skipping this check is the single most common way flow-warped smoothing goes wrong:
    without it, colour is dragged across occlusion boundaries and produces comet-tail
    artifacts that are worse than the flicker being removed.

    ``flow_fw`` is reference->other, ``flow_bw`` is other->reference.
    """
    if threshold_px <= 0:
        raise ValueError("threshold_px must be positive")
    warped_bw = warp_by_flow(flow_bw, flow_fw)
    residual = flow_fw + warped_bw
    error = np.sqrt(residual[..., 0] ** 2 + residual[..., 1] ** 2)
    return error <= threshold_px


def flow_warped_median(
    ab_stack: np.ndarray,
    gray_stack: np.ndarray,
    center: int,
    flow_engine=None,
    fb_threshold_px: float = DEFAULT_FB_THRESHOLD_PX,
    return_diagnostics: bool = False,
):
    """Motion-compensated temporal median of chroma for one frame.

    Parameters
    ----------
    ab_stack : (N, H, W, 2) float32
        Chroma for a window of N frames, in DDColor's native Lab a/b.
    gray_stack : (N, H, W)
        Matching luminance, used only to estimate flow.
    center : int
        Index within the window of the frame being filtered.
    flow_engine :
        A ``cv2.DISOpticalFlow``; one is created if omitted.
    fb_threshold_px : float
        Forward-backward round-trip tolerance in pixels.
    return_diagnostics : bool
        If True, also return a dict with the per-pixel valid-sample count.

    Returns
    -------
    (H, W, 2) float32 chroma for the centre frame.

    Pixels where every neighbour was rejected fall back to the centre frame's own
    chroma -- that is, we trust DDColor rather than inventing a value.
    """
    ab_stack = np.asarray(ab_stack, dtype=np.float32)
    gray_stack = np.asarray(gray_stack)

    if ab_stack.ndim != 4 or ab_stack.shape[3] != 2:
        raise ValueError(f"ab_stack must be (N, H, W, 2), got {ab_stack.shape}")
    if gray_stack.shape[:3] != ab_stack.shape[:3]:
        raise ValueError(
            f"gray_stack {gray_stack.shape[:3]} does not match ab_stack {ab_stack.shape[:3]}"
        )
    n_frames = ab_stack.shape[0]
    if not 0 <= center < n_frames:
        raise IndexError(f"center {center} outside window of {n_frames}")
    if not np.isfinite(ab_stack).all():
        raise ValueError("ab_stack contains NaN or inf")

    if flow_engine is None:
        flow_engine = make_dis_flow()

    center_gray = _as_flow_input(gray_stack[center])

    # Always include the centre frame itself, unwarped, as a valid sample.
    samples = [ab_stack[center]]
    valid = [np.ones(ab_stack.shape[1:3], dtype=bool)]

    for index in range(n_frames):
        if index == center:
            continue
        other_gray = _as_flow_input(gray_stack[index])

        # centre -> other, so warping `other` by this lands it in centre's frame.
        # Buffers are preallocated rather than passing None: cv2 accepts either, but an
        # explicit array keeps the call unambiguous and avoids a per-call allocation.
        flow_shape = (ab_stack.shape[1], ab_stack.shape[2], 2)
        flow_fw = flow_engine.calc(center_gray, other_gray, np.zeros(flow_shape, np.float32))
        flow_bw = flow_engine.calc(other_gray, center_gray, np.zeros(flow_shape, np.float32))

        mask = fb_consistency_mask(flow_fw, flow_bw, fb_threshold_px)
        samples.append(warp_by_flow(ab_stack[index], flow_fw))
        valid.append(mask)

    stacked = np.stack(samples, axis=0)          # (N, H, W, 2)
    mask_stack = np.stack(valid, axis=0)         # (N, H, W)

    # Median over valid samples only. Invalid entries become NaN so nanmedian skips
    # them; the centre frame is always valid, so no pixel can be all-NaN.
    masked = np.where(mask_stack[..., None], stacked, np.nan)
    result = np.nanmedian(masked, axis=0).astype(np.float32)

    if not np.isfinite(result).all():
        raise AssertionError("median produced non-finite output; this should be impossible")

    if return_diagnostics:
        return result, {"valid_counts": mask_stack.sum(axis=0)}
    return result


def naive_median(ab_stack: np.ndarray, center: int) -> np.ndarray:
    """Plain per-pixel temporal median, for A/B comparison against the above."""
    ab_stack = np.asarray(ab_stack, dtype=np.float32)
    if ab_stack.ndim != 4 or ab_stack.shape[3] != 2:
        raise ValueError(f"ab_stack must be (N, H, W, 2), got {ab_stack.shape}")
    if not 0 <= center < ab_stack.shape[0]:
        raise IndexError("center outside window")
    return np.median(ab_stack, axis=0).astype(np.float32)


def upsample_chroma(ab: np.ndarray, height: int, width: int) -> np.ndarray:
    """Resize chroma to full resolution with BILINEAR interpolation.

    DDColor's shipped pipeline calls ``F.interpolate`` with no ``mode``, and PyTorch's
    default is ``"nearest"``. Nearest-neighbour upsampling from 512x512 to 4K replicates
    each chroma sample into a ~7.5x4.2 px block with a hard edge; as content drifts
    sub-pixel between frames, output pixels SNAP between neighbouring chroma samples at
    different moments across the frame. That is chroma jitter generated by the resize
    itself, and it will otherwise be measured as, and blamed on, DDColor's flicker.

    Bilinear is also the perceptually correct choice: chroma is a smooth signal. (Nearest
    is mandatory only for label/mask images, where interpolating invents classes.)
    """
    if ab.ndim != 3 or ab.shape[2] != 2:
        raise ValueError(f"ab must be (H, W, 2), got {ab.shape}")
    if height <= 0 or width <= 0:
        raise ValueError("height and width must be positive")
    return cv2.resize(ab, (width, height), interpolation=cv2.INTER_LINEAR)
