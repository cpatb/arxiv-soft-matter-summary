"""Fit a whole phi-series at once, sharing the parameters that ought to be shared.

Fitting each flow curve on its own asks less of the model than the model claims.
Wyart-Cates says ``phi_m`` is a particle property and ``sigma*`` is set by the
repulsive force scale, so neither should move when you put more starch in the
pot. Six independent fits give those two parameters six chances to absorb
whatever else is wrong, and the resulting spread then looks like a measurement.

The joint fit removes that freedom. One ``phi_m``, one ``sigma*``, one
``phi_0``, and one ``phi`` per curve, fit to every point at once. The comparison
against the independent fits is then a real model-selection question rather than
a matter of taste, and it is settled with AIC so that the extra parameters have
to be paid for.

``phi_0`` is freed here because it was pinned at 0.64 in the single-curve fits.
That value is random close packing for frictionless monodisperse spheres. Corn
starch granules are polydisperse, angular and swollen, so 0.64 is a guess, and a
wrong ``phi_0`` is absorbed directly by ``phi_m`` through
``phi_J = (1-f) phi_0 + f phi_m``. Freeing it is the cheapest way to find out
whether the ``phi_m`` drift was ever about ``phi_m``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from dst_fit.models import (
    ETA_S_WATER,
    WC_BETA,
    WC_EXPONENT,
    divergence,
    frictional_fraction,
    jamming_fraction,
)

_PENALTY = 25.0

# Bounds on the shared parameters. phi_0 is allowed to roam between random loose
# packing and a little past FCC so that the fit can say something implausible if
# that is what the data wants; the gate, not the bounds, decides admissibility.
PHI0_BOUNDS = (0.50, 0.74)
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
            f"    phi_0      = {self.phi_0:.4f} +/- {e['phi_0']:.4f}",
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


def _eta(sigma, phi, phi_m, sigma_star, phi_0, eta_s=ETA_S_WATER, n=WC_EXPONENT, beta=WC_BETA):
    f = frictional_fraction(sigma, sigma_star, beta=beta)
    phi_j = jamming_fraction(f, phi_0=phi_0, phi_m=phi_m)
    return divergence(np.broadcast_to(float(phi), np.shape(phi_j)), phi_j, eta_s=eta_s, n=n)


def fit_series(
    curves: list[tuple[np.ndarray, np.ndarray]],
    *,
    shared_sigma_star: bool = True,
    eta_s: float = ETA_S_WATER,
    n: float = WC_EXPONENT,
    beta: float = WC_BETA,
    n_starts: int = 12,
    seed: int = 0,
) -> JointFit:
    """Fit ``[(stress, viscosity), ...]`` with shared phi_0, phi_m and optionally sigma*.

    Every curve gets its own ``phi``. Residuals are in log viscosity, as in the
    single-curve fitter, so the two are directly comparable.

    ``eta_s``, ``n`` and ``beta`` are the same fixed constants the single-curve
    models take (see :mod:`dst_fit.models`), and must be set to the same values
    used for the independent fits or the AIC comparison is not like for like.
    ``phi_0`` is *not* among them: freeing it is the whole point of this module.
    """
    n_curves = len(curves)
    if n_curves < 2:
        raise ValueError("a joint fit needs at least two curves")

    sigmas = [np.asarray(s, dtype=float) for s, _ in curves]
    log_etas = [np.log(np.asarray(e, dtype=float)) for _, e in curves]
    n_points = int(sum(s.size for s in sigmas))

    n_ss = 1 if shared_sigma_star else n_curves
    n_params = 2 + n_ss + n_curves

    lo = np.array(
        [PHI0_BOUNDS[0], PHIM_BOUNDS[0]]
        + [SIGMA_STAR_BOUNDS[0]] * n_ss
        + [PHI_BOUNDS[0]] * n_curves
    )
    hi = np.array(
        [PHI0_BOUNDS[1], PHIM_BOUNDS[1]]
        + [SIGMA_STAR_BOUNDS[1]] * n_ss
        + [PHI_BOUNDS[1]] * n_curves
    )

    def unpack(p):
        phi_0, phi_m = p[0], p[1]
        ss = p[2 : 2 + n_ss]
        phis = p[2 + n_ss :]
        return phi_0, phi_m, ss, phis

    def residual(p):
        phi_0, phi_m, ss, phis = unpack(p)
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
        np.concatenate([[0.64, 0.58], [5.0] * n_ss, np.full(n_curves, 0.50)]),
        np.concatenate([[0.58, 0.52], [1.0] * n_ss, np.full(n_curves, 0.45)]),
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
    try:
        cov = s2 * np.linalg.inv(best.jac.T @ best.jac)
        err = np.sqrt(np.abs(np.diag(cov)))
    except np.linalg.LinAlgError:
        err = np.full(n_params, np.nan)

    phi_0, phi_m, ss, phis = unpack(best.x)
    e0, e1, ess, ephi = unpack(err)
    return JointFit(
        phi_0=float(phi_0),
        phi_m=float(phi_m),
        sigma_star=float(ss[0]) if shared_sigma_star else np.asarray(ss),
        phi=np.asarray(phis),
        stderr={
            "phi_0": float(e0),
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
