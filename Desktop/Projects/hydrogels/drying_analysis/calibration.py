"""
Spatial calibration: pixels -> physical length, from a fiducial of known size.

The drying-fracture roadmap requires "fix px/mm per run" (no hardcoded scale).
Every measurement module in this package takes a `px_per_mm` argument; this
module is the one place that value is derived, from a fiducial (ruler tick,
calibration dot, or the mold edge) of known physical length.

Method: pure ratio of measured pixel distance to known physical distance.
No interpolation, no assumptions about image content.
"""
from __future__ import annotations

import numpy as np


def px_per_mm_from_fiducial(p0, p1, known_length_mm: float) -> float:
    """px/mm from two endpoint pixel coordinates of a fiducial of known length.

    Parameters
    ----------
    p0, p1 : (row, col) pixel coordinates of the fiducial endpoints.
    known_length_mm : the true physical distance between p0 and p1, in mm.

    Returns
    -------
    px_per_mm : float
    """
    if known_length_mm <= 0:
        raise ValueError(f"known_length_mm must be positive, got {known_length_mm}")
    p0 = np.asarray(p0, dtype=float)
    p1 = np.asarray(p1, dtype=float)
    dist_px = float(np.hypot(*(p1 - p0)))
    if dist_px < 1.0:
        raise ValueError(
            f"fiducial endpoints are <1 px apart ({dist_px:.3f} px) -- pick "
            "endpoints spanning a resolvable distance, or the scale is unreliable.")
    return dist_px / known_length_mm


def mm_per_px(px_per_mm: float) -> float:
    if px_per_mm <= 0:
        raise ValueError(f"px_per_mm must be positive, got {px_per_mm}")
    return 1.0 / px_per_mm
