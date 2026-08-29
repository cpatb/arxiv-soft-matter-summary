"""
Top-view crack metrics: polygon (fragment) area A_p and crack spacing Lambda.

Given a binary crack mask (True on crack lines, from a top-view frame), the gel
fragments are the connected components of the *complement*. For each fragment we
report area A_p and an equivalent spacing. The MOC length ledger expects
A_p ~ h^{4/3} (Ma-Lowensohn-Burton) and Lambda ~ h^{2/3..1}, so these per-fragment
areas vs thickness are the scaling-law test.

Method: skimage.measure.label on ~crack_mask, then regionprops for area and
equivalent diameter. Fragments touching the image border are excluded by default
(they are truncated by the field of view, so their area is not physical).
Crack spacing Lambda is reported as the median equivalent diameter of interior
fragments (a robust central spacing), plus the mean sqrt(area).

All areas convert to mm^2 via px_per_mm^2; lengths via px_per_mm.

Reference: Ma, Lowensohn & Burton, Phys. Rev. E 99, 012802 (2019).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from skimage.measure import label, regionprops
from skimage.segmentation import clear_border


@dataclass
class CrackMetrics:
    n_fragments: int
    areas_mm2: np.ndarray
    equiv_diam_mm: np.ndarray
    median_area_mm2: float
    median_spacing_mm: float          # median equivalent diameter
    mean_spacing_mm: float            # mean sqrt(area)
    per_fragment: list = field(default_factory=list)


def crack_metrics(crack_mask, px_per_mm: float, min_area_px: int = 9,
                  exclude_border: bool = True) -> CrackMetrics:
    """Fragment areas and crack spacing from a binary top-view crack mask.

    Parameters
    ----------
    crack_mask : 2-D bool array, True on crack pixels.
    px_per_mm : spatial calibration.
    min_area_px : drop fragments smaller than this (noise / speckle).
    exclude_border : drop fragments touching the image edge (truncated by FOV).
    """
    mask = np.asarray(crack_mask)
    if mask.ndim != 2:
        raise ValueError(f"crack_mask must be 2-D; got {mask.ndim}-D")
    if px_per_mm <= 0:
        raise ValueError("px_per_mm must be positive")
    # Guard against a float probability/intensity map being passed by mistake:
    # a non-boolean input whose values are not a subset of {0,1} would make
    # ~mask.astype(bool) treat every nonzero pixel as crack, silently corrupting
    # the fragment segmentation.
    if mask.dtype != bool:
        uniq = np.unique(mask)
        if not np.all(np.isin(uniq, (0, 1))):
            raise ValueError(
                "crack_mask must be boolean or 0/1-valued; got values "
                f"{uniq[:6]}... -- threshold your intensity/probability map first.")
    fragments = ~mask.astype(bool)
    lbl = label(fragments, connectivity=1)
    if exclude_border:
        lbl = clear_border(lbl)

    px_area_mm2 = 1.0 / px_per_mm**2
    areas, diams, rows = [], [], []
    for p in regionprops(lbl):
        if p.area < min_area_px:
            continue
        a_mm2 = p.area * px_area_mm2
        d_mm = p.equivalent_diameter_area / px_per_mm
        areas.append(a_mm2)
        diams.append(d_mm)
        rows.append({"label": p.label, "area_mm2": a_mm2,
                     "equiv_diam_mm": d_mm, "centroid_rc": p.centroid})

    areas = np.asarray(areas, dtype=float)
    diams = np.asarray(diams, dtype=float)
    if areas.size == 0:
        return CrackMetrics(0, areas, diams, np.nan, np.nan, np.nan, rows)
    return CrackMetrics(
        n_fragments=int(areas.size),
        areas_mm2=areas,
        equiv_diam_mm=diams,
        median_area_mm2=float(np.median(areas)),
        median_spacing_mm=float(np.median(diams)),
        mean_spacing_mm=float(np.mean(np.sqrt(areas))),
        per_fragment=rows,
    )
