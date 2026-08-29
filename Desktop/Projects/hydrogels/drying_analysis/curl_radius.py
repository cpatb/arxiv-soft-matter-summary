"""
Curl radius rho from a side-view arc (circle fit to the curled edge profile).

The side-view profile of a curling ("cinnamon-roll") layer is a circular arc.
Fitting a circle to the traced edge points gives the curl radius rho; the MOC
predicts rho ~ h / (3 eps_bar) (Style-Peppin mud-peel), so rho is compared to
h and the differential drying strain.

Method: algebraic (Kasa) circle fit -- solves a linear least-squares system for
the circle that best fits the points, then reports centre, radius, and the RMS
residual. Algebraic fitting is used (not iterative geometric) because it is
non-iterative, has no starting-guess sensitivity, and is accurate for arcs that
span a reasonable angular range. For very short arcs the fit is ill-conditioned;
the returned arc-span angle flags that case.

Points are (row, col) pixel coordinates; radius is returned in pixels and (if a
calibration is given) in mm.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class CircleFit:
    center_rc: tuple          # (row, col) in pixels
    radius_px: float
    radius_mm: float | None
    rms_residual_px: float
    arc_span_deg: float       # angular extent of the data about the centre


def fit_curl_radius(points_rc, px_per_mm: float | None = None) -> CircleFit:
    """Algebraic (Kasa) circle fit to side-view arc points.

    Parameters
    ----------
    points_rc : (N, 2) array of (row, col) pixel coordinates along the arc.
    px_per_mm : optional calibration; if given, radius_mm is populated.
    """
    pts = np.asarray(points_rc, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError(f"points must be (N,2); got {pts.shape}")
    if pts.shape[0] < 3:
        raise ValueError(f"need >=3 points to fit a circle; got {pts.shape[0]}")
    if px_per_mm is not None and px_per_mm <= 0:
        raise ValueError(f"px_per_mm must be positive if given; got {px_per_mm}")
    y = pts[:, 0]  # row
    x = pts[:, 1]  # col
    # Kasa: minimize sum (x^2+y^2 + Dx + Ey + F)^2  -> linear in (D,E,F)
    A = np.column_stack([x, y, np.ones_like(x)])
    b = -(x**2 + y**2)
    (D, E, F), *_ = np.linalg.lstsq(A, b, rcond=None)
    cx = -D / 2.0
    cy = -E / 2.0
    r = np.sqrt(cx**2 + cy**2 - F)
    if not np.isfinite(r) or r <= 0:
        raise ValueError("degenerate circle fit (collinear points?)")
    resid = np.hypot(x - cx, y - cy) - r
    rms = float(np.sqrt(np.mean(resid**2)))
    angles = np.arctan2(y - cy, x - cx)
    span = float(np.ptp(np.unwrap(np.sort(angles))))
    return CircleFit(
        center_rc=(float(cy), float(cx)),
        radius_px=float(r),
        radius_mm=(float(r / px_per_mm) if px_per_mm else None),
        rms_residual_px=rms,
        arc_span_deg=float(np.degrees(span)),
    )
