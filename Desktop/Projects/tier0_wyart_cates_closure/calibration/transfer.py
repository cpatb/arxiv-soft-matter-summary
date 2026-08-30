"""What a steady flow curve does and does not hand a drop-impact solver.

The tier-0 closure has six parameters the README asks to be calibrated:
``phi``, ``phi_m``, ``sigma_star``, ``sigma_yield``, ``f_res`` and ``tau``. A
steady flow curve is the only measurement available for the job. This module
establishes, on three published datasets whose volume fractions are known
independently, which of the six survive the transfer.

Three protocols, and they are not equivalent
--------------------------------------------
A. **Free fit.** Fit ``(phi, phi_m, sigma_star)``. This is what everyone does.
B. **Pinned fit.** Supply the ``phi`` the sample was made at; fit
   ``(phi_m, sigma_star)``.
C. **Plateau inversion.** Read ``phi`` off the low-stress viscosity with no
   optimiser at all: ``phi = phi_0 (1 - eta_r^(-1/n))``. Only defined where the
   suspension is non-Brownian at the lowest sampled stress; see
   :mod:`calibration.brownian`.

A drop-impact solver advects the *real* volume fraction, the one the
experimentalist weighed out. So it can only use protocol A's ``phi_m`` if
protocol A's ``phi`` is that same number. It is not:

    phi_inferred = a * phi_true + b

Measured, with a curve-level bootstrap, on the two series where the plateau
route is admissible at all: ``a = 1.15 [1.06, 1.34]``, ``b = -0.101 [-0.206,
-0.060]`` on Guy's d = 3770 nm PMMA and ``a = 1.12 [1.04, 1.18]``,
``b = -0.048 [-0.074, -0.013]`` on Richards' corn starch. Neither interval
contains ``(1, 0)``, in two different labs on two different materials. The free
fit gives the same shape more strongly on Guy (``a = 1.24``) and marginally on
Richards.

**Guy's d = 404 nm series cannot test this**, and saying so is half the result:
its Peclet number at the lowest sampled stress runs 1e-4 to 3e-3, so its
low-stress viscosity is Brownian rather than the model's frictionless branch,
and the fit-free route is simply not defined there. Its free-fit interval is
correspondingly wide and does contain the identity.

The map is a property of the data, not of the fitter: on curves generated from
the model itself and fitted back with the same constants, the same code returns
slope 1.0000 and intercept 0.0000 at zero noise and stays inside 1.04 at 3%
log-normal noise.

The gauge, and what survives it
-------------------------------
Wyart-Cates is exactly invariant under ``(phi, phi_m, phi_0) -> c (phi, phi_m,
phi_0)``, because ``phi/phi_J`` is homogeneous of degree zero in the triple. So
``a`` and ``b`` both scale with the assumed ``phi_0`` and neither is meaningful
alone. Their ratio is not:

    phi_x = -b / a

is invariant under the gauge. It is the volume fraction at which the inferred
packing would agree with the real one, and the distance from it measures how far
the fit has drifted.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class MaterialSpec:
    """Everything known about a published system before any curve is fitted.

    ``radius`` and ``temperature`` are needed only for the Peclet admission and
    are allowed to be ``None`` when a source does not record them, in which case
    the plateau route is refused rather than guessed.
    """

    name: str
    eta_s: float                # Pa.s
    phi_0: float
    n: float
    beta: float
    phi_units: str              # "volume" or "weight"
    radius: float | None        # m
    temperature: float | None   # K
    phi_m_published: float | None
    source: str

    def as_fixed(self) -> dict[str, float]:
        return {"eta_s": self.eta_s, "phi_0": self.phi_0, "n": self.n, "beta": self.beta}


@dataclass(frozen=True)
class TransferMap:
    """The affine map from the sample's volume fraction to the inferred one.

    ``slope``/``intercept``/``phi_x`` are point estimates; the ``*_ci`` fields
    are percentile bootstrap intervals over curves, which is the right
    resampling unit because the curves are the independent measurements.
    """

    slope: float
    intercept: float
    r2: float
    n_curves: int
    phi_x: float
    slope_ci: tuple[float, float]
    intercept_ci: tuple[float, float]
    phi_x_ci: tuple[float, float]
    residual_max: float

    @property
    def is_identity(self) -> bool:
        """Does the map contain (1, 0) inside both bootstrap intervals?"""
        return (self.slope_ci[0] <= 1.0 <= self.slope_ci[1]
                and self.intercept_ci[0] <= 0.0 <= self.intercept_ci[1])

    def __str__(self) -> str:
        return (
            f"phi_inferred = {self.slope:.4f} phi {self.intercept:+.4f}   "
            f"(slope 95% [{self.slope_ci[0]:.3f}, {self.slope_ci[1]:.3f}], "
            f"intercept 95% [{self.intercept_ci[0]:+.4f}, {self.intercept_ci[1]:+.4f}])  "
            f"r2 {self.r2:.5f}, n={self.n_curves}, "
            f"phi_x = {self.phi_x:.4f} [{self.phi_x_ci[0]:.4f}, {self.phi_x_ci[1]:.4f}], "
            f"max residual {self.residual_max:.4f}"
        )


def _line(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    slope, intercept = np.polyfit(x, y, 1)
    return float(slope), float(intercept)


def transfer_map(phi_true, phi_inferred, *, n_boot: int = 4000, seed: int = 20260830,
                 alpha: float = 0.05) -> TransferMap:
    """Regress an inferred volume fraction on the known one, with a bootstrap.

    Resamples *curves*, not points: each curve contributes one (true, inferred)
    pair and the curves are what is independent. Requires at least three pairs
    and refuses rather than extrapolating below that.
    """
    x = np.asarray(phi_true, dtype=float)
    y = np.asarray(phi_inferred, dtype=float)
    if x.shape != y.shape:
        raise ValueError("phi_true and phi_inferred must have the same shape")
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 3:
        raise ValueError(f"need at least 3 curves for a transfer map, got {x.size}")
    if np.ptp(x) <= 0:
        raise ValueError("all curves share one volume fraction; no map is defined")

    slope, intercept = _line(x, y)
    r2 = float(np.corrcoef(x, y)[0, 1] ** 2)
    resid = y - (slope * x + intercept)

    rng = np.random.default_rng(seed)
    sl_b, ic_b, px_b = [], [], []
    for _ in range(n_boot):
        idx = rng.integers(0, x.size, x.size)
        xb, yb = x[idx], y[idx]
        if np.ptp(xb) <= 0:
            continue
        s, i = _line(xb, yb)
        if s == 0:
            continue
        sl_b.append(s)
        ic_b.append(i)
        px_b.append(-i / s)
    lo, hi = 100 * alpha / 2, 100 * (1 - alpha / 2)

    def ci(v):
        return (float(np.percentile(v, lo)), float(np.percentile(v, hi)))

    return TransferMap(
        slope=slope, intercept=intercept, r2=r2, n_curves=int(x.size),
        phi_x=float(-intercept / slope) if slope else float("nan"),
        slope_ci=ci(sl_b), intercept_ci=ci(ic_b), phi_x_ci=ci(px_b),
        residual_max=float(np.max(np.abs(resid))),
    )


def jamming_stress(phi, phi_m, phi_0, sigma_star, beta):
    """Stress at which the material shear-jams, or ``inf`` if it never does.

    Jamming is ``phi_J(f) = phi``, so the frictional fraction required is

        f_jam = (phi_0 - phi) / (phi_0 - phi_m),

    and inverting ``f = exp(-(sigma*/sigma)^beta)`` gives

        sigma_jam = sigma* (-ln f_jam)^(-1/beta).

    ``f_jam >= 1`` means the fully frictional state still flows at this packing,
    so the material never jams and the stress is ``+inf``. ``f_jam <= 0`` means
    it is jammed at rest and the stress is 0.

    ``f_jam`` is a ratio of differences of the packing triple, so it is
    invariant under the model's ``(phi, phi_m, phi_0) -> c(...)`` gauge, and so
    is ``sigma_jam``. This is the one calibrated number that a change of the
    assumed frictionless packing cannot move.
    """
    phi, phi_m, phi_0 = float(phi), float(phi_m), float(phi_0)
    if phi_0 <= phi_m:
        raise ValueError(f"phi_0 ({phi_0}) must exceed phi_m ({phi_m})")
    f_jam = (phi_0 - phi) / (phi_0 - phi_m)
    if f_jam >= 1.0:
        return float("inf")
    if f_jam <= 0.0:
        return 0.0
    return float(sigma_star * (-np.log(f_jam)) ** (-1.0 / beta))


def jams_under_stress(phi, phi_m):
    """Will a suspension at packing ``phi`` shear-jam at any stress?

    True exactly when ``phi > phi_m``: at ``f -> 1`` the jamming point is
    ``phi_m``, so anything denser has no flowing high-stress branch. This is the
    single yes/no a drop-impact solver asks of the rheology.
    """
    return bool(float(phi) > float(phi_m))


@dataclass(frozen=True)
class VerdictBand:
    """Range of packings over which one material's own curves disagree.

    ``lo`` and ``hi`` are the smallest and largest ``phi_m`` any admitted curve
    returned. A sample made between them is called jammed by one curve's
    calibration and flowing by another's, from the same material in the same
    paper.
    """

    lo: float
    hi: float
    labels: tuple[str, ...] = field(default_factory=tuple)

    @property
    def width(self) -> float:
        return self.hi - self.lo

    def contains(self, phi: float) -> bool:
        return self.lo < float(phi) < self.hi

    def __str__(self) -> str:
        return (f"phi_m in [{self.lo:.4f}, {self.hi:.4f}], width {self.width:.4f} "
                f"over {len(self.labels)} curves")


def verdict_band(phi_m_values, labels=None) -> VerdictBand:
    v = np.asarray(list(phi_m_values), dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 2:
        raise ValueError("need at least two phi_m values to form a band")
    return VerdictBand(float(v.min()), float(v.max()), tuple(labels) if labels else ())
