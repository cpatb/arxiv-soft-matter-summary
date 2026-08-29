"""Shape and topology of the school envelope, and fission/fusion events.

Every length is reported in body lengths (BL) and every area in BL^2, so the
numbers are comparable across clips shot at different and drifting camera
heights. Dimensionless shape descriptors (circularity, solidity, aspect ratio)
need no normalisation at all and are the most trustworthy quantities here.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy import ndimage as ndi
from skimage import measure

from .segment import Envelope


@dataclass
class ShapeMetrics:
    """Descriptors of one connected school component."""

    label: int
    area_bl2: float
    perimeter_bl: float
    #: 4*pi*A/P^2. 1 for a disc, smaller for elongated or ragged outlines.
    circularity: float
    #: A / A_convex_hull. Measures how much the outline is indented -- a school
    #: that is starting to split has low solidity before it actually separates.
    solidity: float
    #: Ratio of the principal second moments of the occupied area. 1 = isotropic.
    aspect_ratio: float
    #: Orientation of the long axis, radians, image frame.
    orientation_rad: float
    #: Radius of gyration, the natural size scale for a diffuse object.
    gyration_bl: float
    #: Perimeter divided by that of a disc of equal area. 1 = smooth, larger =
    #: convoluted boundary. Reported instead of a fractal dimension, which is not
    #: identifiable over the ~1 decade of scale this imagery spans.
    boundary_roughness: float
    centroid_y: float
    centroid_x: float


def component_shapes(env: Envelope, body_px: float) -> list[ShapeMetrics]:
    """Shape descriptors for every connected component of the envelope."""
    if not env.present or env.n_components == 0:
        return []
    if not np.isfinite(body_px) or body_px <= 0:
        raise ValueError(f"body_px must be finite and positive, got {body_px}")

    out: list[ShapeMetrics] = []
    for p in measure.regionprops(env.labels):
        area_px = float(p.area)
        # Crofton perimeter is markedly less biased than the boundary-pixel count
        # for digitised outlines, which matters because circularity goes as P^2.
        perim_px = float(measure.perimeter_crofton(env.labels == p.label, directions=4))
        if perim_px <= 0:
            continue

        circ = 4.0 * np.pi * area_px / (perim_px**2)
        equiv_perim = 2.0 * np.sqrt(np.pi * area_px)  # perimeter of equal-area disc

        # Second moments of the occupied area give the principal axes.
        l1, l2 = p.inertia_tensor_eigvals  # l1 >= l2 >= 0
        aspect = float(np.sqrt(l1 / l2)) if l2 > 0 else float("inf")
        gyration_px = float(np.sqrt(l1 + l2))

        out.append(
            ShapeMetrics(
                label=int(p.label),
                area_bl2=area_px / body_px**2,
                perimeter_bl=perim_px / body_px,
                circularity=float(circ),
                solidity=float(p.solidity),
                aspect_ratio=aspect,
                orientation_rad=float(p.orientation),
                gyration_bl=gyration_px / body_px,
                boundary_roughness=float(perim_px / equiv_perim),
                centroid_y=float(p.centroid[0]),
                centroid_x=float(p.centroid[1]),
            )
        )
    return out


def frame_summary(env: Envelope, body_px: float) -> dict:
    """One row per frame: the whole-school view.

    When the school has fragmented, area is summed over components while shape
    descriptors describe the largest one, and ``n_components`` records the split.
    """
    shapes = component_shapes(env, body_px)
    if not shapes:
        return {"present": False, "n_components": 0}

    shapes.sort(key=lambda s: s.area_bl2, reverse=True)
    big = shapes[0]
    total = sum(s.area_bl2 for s in shapes)

    row = {"present": True, "n_components": len(shapes), "total_area_bl2": total}
    row.update({f"main_{k}": v for k, v in asdict(big).items()})
    # Fraction of the school in the largest component: 1.0 when intact, and it
    # drops towards 0.5 as a school splits into two comparable halves.
    row["largest_fraction"] = big.area_bl2 / total if total > 0 else np.nan
    return row


# --- fission / fusion --------------------------------------------------------


def match_components(
    prev_labels: np.ndarray, cur_labels: np.ndarray, min_overlap_px: int = 50
) -> dict[int, list[int]]:
    """Map each previous component to the current components it overlaps.

    Overlap-based matching is appropriate here because the school moves only a
    few pixels between consecutive frames at 30 fps, so a component's footprint
    at t and t+1 overlaps heavily. Camera motion is removed before this is
    called; without that, a fast pan would look like a fission event.
    """
    out: dict[int, list[int]] = {}
    prev_ids = [i for i in np.unique(prev_labels) if i != 0]
    for pid in prev_ids:
        m = prev_labels == pid
        hit = cur_labels[m]
        ids, counts = np.unique(hit[hit != 0], return_counts=True)
        out[int(pid)] = [int(i) for i, c in zip(ids, counts) if c >= min_overlap_px]
    return out


def detect_events(labels_series: list[np.ndarray], frame_indices: list[int]) -> list[dict]:
    """Detect fission (1 -> many) and fusion (many -> 1) between frames.

    Returns a list of event dicts with frame index, type, and the component
    counts either side. Consecutive duplicate reports of the same ongoing split
    are collapsed by the caller if desired.
    """
    events: list[dict] = []
    for k in range(1, len(labels_series)):
        prev, cur = labels_series[k - 1], labels_series[k]
        fwd = match_components(prev, cur)
        bwd = match_components(cur, prev)

        for pid, children in fwd.items():
            if len(children) > 1:
                events.append(
                    {
                        "frame": frame_indices[k],
                        "type": "fission",
                        "parent": pid,
                        "n_children": len(children),
                    }
                )
        for cid, parents in bwd.items():
            if len(parents) > 1:
                events.append(
                    {
                        "frame": frame_indices[k],
                        "type": "fusion",
                        "child": cid,
                        "n_parents": len(parents),
                    }
                )
    return events


def packing_profile(
    env: Envelope, fishness: np.ndarray, body_px: float, n_bins: int = 12
) -> dict:
    """Fish-area fraction as a function of depth into the school.

    Distance is measured inward from the envelope boundary via the Euclidean
    distance transform and expressed in body lengths. A flat profile means the
    school is uniformly packed; a rising one means a denser core; a falling one
    means a denser rind. This is the border-versus-interior question that
    Cavagna and co-workers examined in starling flocks, asked of a fish school.
    """
    if not env.present or not env.mask.any():
        return {"depth_bl": [], "packing": [], "n_px": []}

    dist_px = ndi.distance_transform_edt(env.mask)
    outside = ~env.mask
    if outside.sum() < 1000:
        return {"depth_bl": [], "packing": [], "n_px": []}
    thr = float(np.percentile(fishness[outside], 99))
    occupied = fishness > thr

    depth_bl = dist_px / body_px
    max_depth = float(depth_bl.max())
    if max_depth <= 0:
        return {"depth_bl": [], "packing": [], "n_px": []}

    edges = np.linspace(0.0, max_depth, n_bins + 1)
    centres, packing, counts = [], [], []
    for i in range(n_bins):
        sel = env.mask & (depth_bl >= edges[i]) & (depth_bl < edges[i + 1])
        n = int(sel.sum())
        if n < 200:  # too few pixels for a stable fraction
            continue
        centres.append(float(0.5 * (edges[i] + edges[i + 1])))
        packing.append(float(occupied[sel].mean()))
        counts.append(n)

    return {"depth_bl": centres, "packing": packing, "n_px": counts}
