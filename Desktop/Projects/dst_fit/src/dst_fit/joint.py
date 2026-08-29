"""Fit a whole phi-series at once, sharing the parameters that ought to be shared.

Fitting each flow curve on its own asks less of the model than the model claims.
Wyart-Cates says ``phi_m`` is a particle property and ``sigma*`` is set by the
repulsive force scale, so neither should move when you put more starch in the
pot. Six independent fits give those two parameters six chances to absorb
whatever else is wrong, and the resulting spread then looks like a measurement.

The joint fit removes that freedom. One ``phi_m``, one ``sigma*``, and one
``phi`` per curve, fit to every point at once. The comparison against the
independent fits is then a real model-selection question rather than a matter of
taste, and it is settled with AIC so that the extra parameters have to be paid
for.

Why ``phi_0`` is fixed here, and why it used to be free
-------------------------------------------------------
Until 2026-08-29 this module freed ``phi_0`` alongside ``phi_m`` and a per-curve
``phi``, on the reasoning that 0.64 is random close packing for frictionless
monodisperse spheres, that corn starch granules are none of those things, and
that a wrong ``phi_0`` is absorbed by ``phi_m`` through
``phi_J = (1-f) phi_0 + f phi_m``. The reasoning is sound and the conclusion was
wrong, because freeing it cannot answer the question.

Wyart-Cates is **exactly covariant** under ``(phi, phi_m, phi_0) -> c (phi, phi_m,
phi_0)``: ``phi_J`` is homogeneous of degree one, so ``phi/phi_J`` -- the only
place any of the three enters -- is invariant pointwise. Freeing all three
therefore adds a **provably flat direction** to the cost surface. Measured: the
residual Jacobian's smallest singular value was 1.5e-5 against a largest of 119,
with a null vector proportional to ``(phi_0, phi_m, 0, phi)``, and scaling the
triple moved every predicted viscosity by 2e-15.

``np.linalg.inv`` does not raise on that. It returned finite numbers, and the
gates duly printed ``phi_0 = 0.6912 +/- 924.4`` and ``0.6778 +/- 1110`` as though
they were uncertainties. They were arithmetic noise off a singular normal matrix.
``n_params`` also counted the flat parameter, so every joint AIC was 2 too high.

So ``phi_0`` now joins ``eta_s``, ``n`` and ``beta`` as a *fixed constant*
supplied per dataset. It is not offered as a fitted parameter, because there is
no data anywhere that determines it: **for a flow-curve fit ``phi_0`` is a gauge
choice, and only ``phi_m/phi_0`` is ever measured.** The covariance is now formed
with ``pinvh`` behind an explicit rank check that raises rather than returning
numbers.

What this costs: nothing in the conclusions. The joint fit still loses to
independent fits by 55 AIC on the 2023 corn starch and 568 on Richards', and a
2-unit correction does not touch either. What it costs is a gate condition -- see
``scripts/gate_n5.py``, where G6c tested the fitted ``phi_0`` and is retired.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import pinvh
from scipy.optimize import least_squares

from dst_fit.models import (
    ETA_S_WATER,
    PHI_0_RCP,
    WC_BETA,
    WC_EXPONENT,
    divergence,
    frictional_fraction,
    jamming_fraction,
)

_PENALTY = 25.0

# Bounds on the fitted parameters. There is deliberately no PHI0_BOUNDS any more:
# phi_0 is a fixed constant, and a bound on it would imply it is being searched.
# The upper bound on phi_m and phi is derived from phi_0 inside fit_series, since
# both must sit below the frictionless packing fraction or the thickening inverts.
PHIM_BOUNDS = (0.30, 0.70)
SIGMA_STAR_BOUNDS = (1e-4, 1e5)
PHI_BOUNDS = (0.10, 0.70)


@dataclass
class JointFit:
    phi_0: float
    phi_m: float
    sigma_star: float | np.ndarray  # scalar when shared, one per curve otherwise
    phi: np.ndarray
    stderr: dict
    residuals: np.ndarray
    n_points: int
    n_params: int
    success: bool
    shared_sigma_star: bool
    metadata: dict = field(default_factory=dict)

    @property
    def rmse_log(self) -> float:
        return float(np.sqrt(np.mean(self.residuals**2)))

    @property
    def aic(self) -> float:
        """Gaussian AIC on the log-viscosity residuals."""
        n = self.n_points
        ssr = float(self.residuals @ self.residuals)
        return n * np.log(ssr / n) + 2 * self.n_params

    def summary(self) -> str:
        e = self.stderr
        lines = [
            f"joint fit: {'ok' if self.success else 'FAILED'}  "
            f"n={self.n_points} k={self.n_params} rmse={self.rmse_log:.4f} aic={self.aic:.1f}",
            f"    phi_0      = {self.phi_0:.4f}   (fixed, not fitted: a gauge choice)",
            f"    phi_m      = {self.phi_m:.4f} +/- {e['phi_m']:.4f}",
        ]
        if self.shared_sigma_star:
            lines.append(
                f"    sigma_star = {self.sigma_star:.4f} +/- {e['sigma_star']:.4f}  (shared)"
            )
        else:
            ss = np.atleast_1d(self.sigma_star)
            lines.append(
                "    sigma_star = " + ", ".join(f"{v:.3f}" for v in ss) + "  (per curve)"
            )
        lines.append("    phi        = " + ", ".join(f"{v:.4f}" for v in self.phi))
        return "\n".join(lines)


def _stderr_or_raise(jac: np.ndarray, s2: float, n_params: int) -> np.ndarray:
    """Standard errors from the Jacobian, refusing a rank-deficient design.

    ``np.linalg.inv`` happily inverts a numerically singular ``J^T J`` and
    returns enormous finite numbers. That is how this module came to print
    ``phi_0 = 0.6912 +/- 924.4`` for a parameter the data does not constrain at
    all: the caller cannot tell a wide interval from a meaningless one, and
    ``+/- 924`` on a packing fraction reads as the former.

    So the rank is checked explicitly, against ``numpy.linalg.matrix_rank``'s own
    tolerance, and a deficient design raises. ``pinvh`` is then used for the
    inverse because ``J^T J`` is symmetric positive semi-definite by construction
    and a pseudo-inverse is the numerically better-behaved choice even when the
    rank is full.
    """
    jtj = jac.T @ jac
    rank = int(np.linalg.matrix_rank(jtj))
    if rank < n_params:
        s = np.linalg.svd(jtj, compute_uv=False)
        raise np.linalg.LinAlgError(
            f"joint fit: the design is rank {rank} for {n_params} parameters, so at "
            f"least one direction is unconstrained and no standard error exists for "
            f"it. Singular values span {s[0]:.3g} to {s[-1]:.3g}. Refusing to invert; "
            f"a finite number here would be arithmetic noise, not an uncertainty."
        )
    return np.sqrt(np.abs(np.diag(s2 * pinvh(jtj))))


def _eta(sigma, phi, phi_m, sigma_star, phi_0, eta_s=ETA_S_WATER, n=WC_EXPONENT, beta=WC_BETA):
    f = frictional_fraction(sigma, sigma_star, beta=beta)
    phi_j = jamming_fraction(f, phi_0=phi_0, phi_m=phi_m)
    return divergence(np.broadcast_to(float(phi), np.shape(phi_j)), phi_j, eta_s=eta_s, n=n)


def fit_series(
    curves: list[tuple[np.ndarray, np.ndarray]],
    *,
    shared_sigma_star: bool = True,
    phi_0: float = PHI_0_RCP,
    eta_s: float = ETA_S_WATER,
    n: float = WC_EXPONENT,
    beta: float = WC_BETA,
    n_starts: int = 12,
    seed: int = 0,
) -> JointFit:
    """Fit ``[(stress, viscosity), ...]`` with a shared phi_m and optionally sigma*.

    Every curve gets its own ``phi``. Residuals are in log viscosity, as in the
    single-curve fitter, so the two are directly comparable.

    ``phi_0``, ``eta_s``, ``n`` and ``beta`` are fixed constants, the same ones
    the single-curve models take (see :mod:`dst_fit.models`), and must be set to
    the same values used for the independent fits or the AIC comparison is not
    like for like.

    ``phi_0`` is fixed rather than fitted, and that is not a default anyone should
    override by reaching into the parameter vector: the model is exactly covariant
    under a common rescaling of ``(phi, phi_m, phi_0)``, so freeing it adds a flat
    direction and returns noise. See the module docstring.

    Raises
    ------
    numpy.linalg.LinAlgError
        If the design is rank-deficient at the optimum, rather than returning
        standard errors computed from a singular matrix.
    """
    n_curves = len(curves)
    if n_curves < 2:
        raise ValueError("a joint fit needs at least two curves")

    sigmas = [np.asarray(s, dtype=float) for s, _ in curves]
    log_etas = [np.log(np.asarray(e, dtype=float)) for _, e in curves]
    n_points = int(sum(s.size for s in sigmas))

    if not np.isfinite(phi_0) or not (0.0 < phi_0 < 1.0):
        raise ValueError(f"phi_0 must be a packing fraction in (0, 1), got {phi_0!r}")

    n_ss = 1 if shared_sigma_star else n_curves
    # phi_m, the sigma* block, and one phi per curve. phi_0 is NOT counted: it is
    # a fixed constant, and counting it is what made every earlier joint AIC 2 too
    # high. See the module docstring.
    n_params = 1 + n_ss + n_curves

    phim_hi = min(PHIM_BOUNDS[1], phi_0 - 1e-3)
    if phim_hi <= PHIM_BOUNDS[0]:
        raise ValueError(
            f"phi_0 = {phi_0:g} leaves no room for phi_m above {PHIM_BOUNDS[0]:g}"
        )
    lo = np.array(
        [PHIM_BOUNDS[0]] + [SIGMA_STAR_BOUNDS[0]] * n_ss + [PHI_BOUNDS[0]] * n_curves
    )
    hi = np.array(
        [phim_hi] + [SIGMA_STAR_BOUNDS[1]] * n_ss + [min(PHI_BOUNDS[1], phim_hi)] * n_curves
    )

    def unpack(p):
        phi_m = p[0]
        ss = p[1 : 1 + n_ss]
        phis = p[1 + n_ss :]
        return phi_m, ss, phis

    def residual(p):
        phi_m, ss, phis = unpack(p)
        out = []
        for i in range(n_curves):
            s_star = ss[0] if shared_sigma_star else ss[i]
            with np.errstate(all="ignore"):
                pred = _eta(sigmas[i], phis[i], phi_m, s_star, phi_0,
                            eta_s=eta_s, n=n, beta=beta)
            r = np.full(pred.shape, _PENALTY)
            ok = np.isfinite(pred) & (pred > 0)
            r[ok] = np.log(pred[ok]) - log_etas[i][ok]
            out.append(r)
        return np.concatenate(out)

    rng = np.random.default_rng(seed)
    starts = [
        np.concatenate([[0.58], [5.0] * n_ss, np.full(n_curves, 0.50)]),
        np.concatenate([[0.52], [1.0] * n_ss, np.full(n_curves, 0.45)]),
    ]
    for _ in range(max(0, n_starts - len(starts))):
        u = rng.uniform(size=n_params)
        starts.append(lo + u * (hi - lo))

    best = None
    for x0 in starts:
        try:
            sol = least_squares(
                residual, np.clip(x0, lo + 1e-9, hi - 1e-9),
                bounds=(lo, hi), method="trf", max_nfev=20000,
            )
        except (ValueError, np.linalg.LinAlgError):
            continue
        if best is None or sol.cost < best.cost:
            best = sol
    if best is None:
        raise RuntimeError("joint fit: every start failed")

    dof = max(best.fun.size - n_params, 1)
    s2 = float(best.fun @ best.fun) / dof
    err = _stderr_or_raise(best.jac, s2, n_params)

    phi_m, ss, phis = unpack(best.x)
    e1, ess, ephi = unpack(err)
    return JointFit(
        phi_0=float(phi_0),
        phi_m=float(phi_m),
        sigma_star=float(ss[0]) if shared_sigma_star else np.asarray(ss),
        phi=np.asarray(phis),
        stderr={
            # phi_0 is fixed, so its uncertainty is not a small number, it is
            # not a number. NaN rather than 0.0: 0.0 would read as "measured
            # perfectly", which is the opposite of what is true.
            "phi_0": float("nan"),
            "phi_m": float(e1),
            "sigma_star": float(ess[0]) if shared_sigma_star else np.asarray(ess),
            "phi": np.asarray(ephi),
        },
        residuals=best.fun,
        n_points=n_points,
        n_params=n_params,
        success=bool(best.success),
        shared_sigma_star=shared_sigma_star,
        metadata={
            "n_starts": len(starts),
            "seed": seed,
            "fixed": {"eta_s": float(eta_s), "n": float(n), "beta": float(beta)},
        },
    )


def aic_from_residuals(residuals: np.ndarray, n_params: int) -> float:
    n = residuals.size
    ssr = float(residuals @ residuals)
    return n * np.log(ssr / n) + 2 * n_params
