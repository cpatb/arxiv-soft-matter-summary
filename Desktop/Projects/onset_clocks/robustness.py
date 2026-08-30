"""Which conclusions survive the fact that the calibration behind them is withdrawn?

The problem
-----------
Every number in this project reaches the material through
``impact_stress.materials.lab_silica()``, which ships the ``si_fit`` preset from
``Github/Drop Impact/closure/params_fitted.yaml``. That preset was fitted on
2026-07-20 to ``E:\\rheology``, and:

- **those flow curves were retracted on 2026-07-29**
  (``Desktop/Papers/STATE_OF_THE_FIELD/OPEN_PROBLEMS_MASTER.md``: "the eleven
  flow curves were retracted 2026-07-29 ... No sigma* ground truth on any axis");
- **``si_fit`` itself was withdrawn on 2026-08-26**
  (``basilisk-drop-impact/WC_DNS_MASTER_PLAN.md``, entry A4 -> S1: "No cleared
  silica rheology exists (E:\\rheology RETRACTED; 26 FEB 26 in-house barred by
  the 2026-08-24 policy; nothing published on disk)");
- the named defect is that ``eta_s`` = 4.736 mPa.s is 4.736x water for a
  *water*-based suspension, and the thickening ceiling is exactly proportional
  to ``eta_s``, so the fit bought 4.736x more achievable thickening than water
  allows. Held at water, that entry reports ``phi_0`` falling 0.70 -> 0.5306.

So the parameters cannot be refitted (the data is gone) and cannot be trusted
(the fit is withdrawn). What can still be done is to ask **over what range of
those parameters each conclusion survives** -- which is a stronger statement
than a fit would have given, because it does not depend on one.

The ranges, and why each
------------------------
``eta_s``  1.0 to 4.736 mPa.s: water at one end, the withdrawn fit at the other.
``phi_0``  0.5306 to 0.70: the value the withdrawal entry derives when eta_s is
           held at water, up to the fitted one.
``phi_m``  0.44 to 0.58: from the lowest packing that arrests in the dataset up
           past random close packing for frictional spheres.
``sigma*`` 4.4 to 442 Pa: the fitted value times 0.1 and 10, which brackets both
           the definitional 3-5x offset between "where f = 1/e" and "where
           thickening visibly begins" and the 7.53x salt span.

A conclusion that holds on every corner of that box does not depend on the
withdrawn calibration at all.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from itertools import product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from clocks import wall_shear_gauge  # noqa: E402

ETA_S_RANGE = (1.0e-3, 4.736e-3)
PHI_0_RANGE = (0.5306, 0.70)
PHI_M_RANGE = (0.44, 0.58)
SIGMA_STAR_RANGE = (4.42, 442.3)
DELTA = 1.3852          # the activation sharpness; held, and swept separately


@dataclass(frozen=True)
class Params:
    """One corner of the parameter box."""

    eta_s: float
    phi_0: float
    phi_m: float
    sigma_star: float
    delta: float = DELTA

    @property
    def valid(self) -> bool:
        """phi_0 must exceed phi_m or the model thins instead of thickening."""
        return self.phi_0 > self.phi_m


def corners(n_per_axis: int = 3):
    """The parameter box, sampled on a grid. Invalid corners are skipped."""
    def grid(lo, hi):
        return np.linspace(lo, hi, n_per_axis)

    for e, p0, pm, ss in product(grid(*ETA_S_RANGE), grid(*PHI_0_RANGE),
                                 grid(*PHI_M_RANGE), np.geomspace(*SIGMA_STAR_RANGE,
                                                                 n_per_axis)):
        q = Params(float(e), float(p0), float(pm), float(ss))
        if q.valid:
            yield q


# --------------------------------------------------------------------------
# the conclusions, each as a predicate over (drops, params)
# --------------------------------------------------------------------------
def can_any_drop_jam(drops, q: Params) -> bool:
    """C1: does ANY packing in the dataset reach phi_J = phi at f = 1?

    The project's premise is that none can. Returns True if the premise fails.
    """
    return any((q.phi_0 - d.phi) / (q.phi_0 - q.phi_m) <= 1.0 for d in drops)


def f_range(drops, q: Params) -> float:
    """C2: dynamic range of the frictional fraction across the dataset."""
    f = np.array([np.exp(-((q.sigma_star / wall_shear_gauge(d)) ** q.delta))
                  for d in drops])
    f = f[np.isfinite(f) & (f > 0)]
    return float(f.max() / f.min()) if f.size else float("nan")


def all_above_ceiling(drops, q: Params, p: float = 5.0) -> tuple[bool, float]:
    """C3: does every arrested drop demand more viscosity than f = 1 allows?

    Returns (all above, the smallest gap). The ceiling is
    ``eta_s (1 - phi/phi_m)^-2``, infinite where phi >= phi_m.
    """
    gaps = []
    for d in drops:
        if d.phi >= q.phi_m:
            return False, 0.0          # the model jams outright; no ceiling
        ceiling = q.eta_s * (1.0 - d.phi / q.phi_m) ** -2.0
        rho = d.phi * 2000.0 + (1 - d.phi) * 1000.0
        implied = rho * d.u0 * d.d0 / d.beta_p**p
        gaps.append(implied / ceiling)
    gaps = np.array(gaps)
    return bool((gaps > 1).all()), float(gaps.min())


def marginal_jamming_offset(drops, q: Params, p: float = 5.0) -> float:
    """C4: the largest |phi_J(at arrest) - phi| over the drops. Should be small."""
    off = []
    for d in drops:
        rho = d.phi * 2000.0 + (1 - d.phi) * 1000.0
        eta = rho * d.u0 * d.d0 / d.beta_p**p
        r = eta / q.eta_s
        if r <= 1:
            continue
        gap = r**-0.5
        if gap >= 1:
            continue
        off.append(abs(d.phi / (1.0 - gap) - d.phi))
    return float(np.max(off)) if off else float("nan")
