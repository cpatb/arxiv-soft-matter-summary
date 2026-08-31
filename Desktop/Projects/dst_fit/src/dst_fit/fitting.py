"""Least-squares fitting of the models in :mod:`dst_fit.models`, with uncertainty.

Residuals are taken in log viscosity. A suspension flow curve spans decades, so
a linear residual is dominated entirely by the few highest-viscosity points and
the fit stops being a fit of the curve.

Parameter uncertainty comes from the Jacobian at the optimum,

    cov = s^2 * (J^T J)^-1,    s^2 = SSR / (n - p)

which is the standard linearised estimate. It is a local, Gaussian approximation
and it is optimistic wherever the model is strongly nonlinear in its parameters,
which for Wyart-Cates it is. The correlation matrix is returned alongside so a
caller can see when two parameters are trading off, rather than reading a small
standard error as though it meant the parameter was pinned. Turning that into an
actual verdict is N3; this module only supplies the numbers.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from dst_fit.models import Model, get_model

# Residual returned in place of a non-finite prediction. Large enough to push the
# optimiser away from jammed or unyielded parameter regions, small enough not to
# swamp the Jacobian with infinities.
_PENALTY = 25.0


@dataclass
class FitResult:
    model: Model
    params: np.ndarray
    stderr: np.ndarray
    cov: np.ndarray
    corr: np.ndarray
    residuals: np.ndarray
    n_points: int
    success: bool
    message: str
    branch: slice | None = None
    metadata: dict = field(default_factory=dict)

    @property
    def n_free(self) -> int:
        return self.model.n_params

    @property
    def dof(self) -> int:
        return self.n_points - self.n_free

    @property
    def rmse_log(self) -> float:
        """Root-mean-square residual in log viscosity, i.e. a relative error."""
        return float(np.sqrt(np.mean(self.residuals**2)))

    def as_dict(self) -> dict[str, tuple[float, float]]:
        return {
            name: (float(p), float(e))
            for name, p, e in zip(self.model.param_names, self.params, self.stderr)
        }

    def parameters_at_bound(self, tol: float = 1e-4) -> tuple[str, ...]:
        """Names of parameters sitting on one of their bounds.

        A bound-hitting fit is a *constrained* optimum. The cost surface was still
        descending when the optimiser ran out of room, so the returned value is a
        property of the bound, not of the data, and its standard error is
        meaningless: the quadratic approximation is taken at a point that is not a
        stationary point.

        This is not hypothetical. Guy et al.'s d = 3770 nm 19% curve fits to
        ``phi_m = 0.3000 +/- 2.3768``, which is exactly the lower bound with an
        error bar seven times the plausible range of the quantity. At ``phi``
        around 0.10 that suspension barely thickens, so the curve carries no
        information about ``phi_m`` and the fitter returned the bound instead of
        refusing. Counting it as a measurement dragged a cross-curve consistency
        test from 3-of-4 to 1-of-5.

        ``tol`` is a fraction of the parameter's own bound range, measured **in
        log space for a parameter whose bounds span more than two decades** and
        linearly otherwise. That distinction is not cosmetic. ``sigma_star`` is
        bounded over 1e-4 to 1e5, so a linear span-relative test calls anything
        below 10 "at the lower bound", which flags every ordinary fit. A packing
        fraction bounded over 0.3 to 0.64 has to be judged linearly, because zero
        is not a meaningful origin for it.
        """
        lo = np.asarray(self.model.bounds[0], dtype=float)
        hi = np.asarray(self.model.bounds[1], dtype=float)
        p = np.asarray(self.params, dtype=float)

        logarithmic = (lo > 0) & (hi > 0) & (hi / np.where(lo > 0, lo, 1.0) > 100.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            lo_c = np.where(logarithmic, np.log10(np.maximum(lo, 1e-300)), lo)
            hi_c = np.where(logarithmic, np.log10(np.maximum(hi, 1e-300)), hi)
            p_c = np.where(logarithmic, np.log10(np.maximum(p, 1e-300)), p)

        span = np.where(hi_c > lo_c, hi_c - lo_c, 1.0)
        near = (np.abs(p_c - lo_c) <= tol * span) | (np.abs(p_c - hi_c) <= tol * span)
        near &= np.isfinite(p_c)
        return tuple(n for n, flag in zip(self.model.param_names, near) if flag)

    def worst_correlation(self) -> tuple[str, str, float]:
        """The most strongly correlated parameter pair, by absolute value."""
        c = np.abs(self.corr - np.eye(self.n_free))
        i, j = np.unravel_index(int(np.argmax(c)), c.shape)
        names = self.model.param_names
        return names[i], names[j], float(self.corr[i, j])

    def summary(self) -> str:
        lines = [
            f"{self.model.name}: {'ok' if self.success else 'FAILED'}  "
            f"n={self.n_points} dof={self.dof} rmse(log eta)={self.rmse_log:.4f}"
        ]
        for name, (p, e) in self.as_dict().items():
            rel = f" ({100 * e / abs(p):.1f}%)" if p != 0 and np.isfinite(e) else ""
            lines.append(f"    {name:>12s} = {p:12.6g} +/- {e:<10.4g}{rel}")
        a, b, r = self.worst_correlation()
        lines.append(f"    worst correlation: {a} vs {b} = {r:+.4f}")
        return "\n".join(lines)


def _initial_guess(model: Model, sigma: np.ndarray, eta: np.ndarray) -> np.ndarray:
    """Heuristic start, clipped into bounds. Multi-start covers what this misses."""
    lo, hi = np.array(model.bounds[0]), np.array(model.bounds[1])
    name = model.name
    if name == "wyart_cates":
        # sigma* sits somewhere in the measured stress window; phi below phi_m.
        guess = np.array([0.50, 0.58, float(np.sqrt(sigma.min() * sigma.max()))])
    elif name == "wyart_cates_known_phi":
        # phi_m must exceed the supplied phi or the curve is jammed everywhere,
        # so start just above it rather than at a generic 0.58.
        known = float(model.fixed.get("phi", 0.50))
        guess = np.array([min(known * 1.05, hi[0]), float(np.sqrt(sigma.min() * sigma.max()))])
    elif name == "maron_pierce":
        guess = np.array([0.50, 0.58])
    elif name == "herschel_bulkley":
        guess = np.array([float(sigma.min()) * 0.5, float(np.median(eta)), 1.0])
    elif name == "ellis":
        guess = np.array(
            [float(eta.max()), float(eta.min()), float(np.median(sigma)), 1.0]
        )
    else:  # pragma: no cover - registry is closed
        guess = 0.5 * (lo + hi)
    return np.clip(guess, lo + 1e-12, hi - 1e-12)


def _residual_factory(model: Model, sigma: np.ndarray, log_eta: np.ndarray):
    def residual(params: np.ndarray) -> np.ndarray:
        with np.errstate(all="ignore"):
            pred = np.asarray(model(sigma, *params), dtype=float)
        out = np.full(pred.shape, _PENALTY)
        ok = np.isfinite(pred) & (pred > 0)
        out[ok] = np.log(pred[ok]) - log_eta[ok]
        return out

    return residual


def _covariance(jac: np.ndarray, residuals: np.ndarray, n_free: int):
    """Linearised covariance, or NaNs if the Jacobian is singular."""
    dof = residuals.size - n_free
    nan = np.full((n_free, n_free), np.nan)
    if dof <= 0:
        return nan, nan
    s2 = float(residuals @ residuals) / dof
    jtj = jac.T @ jac
    try:
        cov = s2 * np.linalg.inv(jtj)
    except np.linalg.LinAlgError:
        return nan, nan
    d = np.sqrt(np.diag(cov))
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = cov / np.outer(d, d)
    return cov, corr


def fit_model(
    model: Model | str,
    stress: np.ndarray,
    viscosity: np.ndarray,
    *,
    branch: slice | None = None,
    n_starts: int = 8,
    seed: int = 0,
    start: np.ndarray | None = None,
) -> FitResult:
    """Fit ``model`` to a measured (stress, viscosity) pair.

    Parameters
    ----------
    branch
        Optional slice restricting the fit, e.g. the thickening branch for
        Wyart-Cates. Recorded on the result so the choice travels with it.
    n_starts
        Multi-start count. The heuristic start is always tried first; the rest
        are drawn log-uniformly inside the bounds with a fixed seed, so the
        result is reproducible. Single-start fits of these models are known to
        be initial-guess sensitive.
    start
        An extra starting point tried *before* the heuristic one and counted on
        top of ``n_starts``. Its use is bootstrap resampling: a resample is a
        perturbation of the data the full fit was made on, so the full fit's
        parameters are the best available guess, and starting a resample from
        the generic heuristic instead throws that away. Doing so made the iid
        control interval come back *narrower* than Wald -- an optimiser
        artefact, since a resample that converges to a worse optimum reports a
        parameter closer to wherever the generic start led, not a real spread.
    """
    mdl: Model = get_model(model) if isinstance(model, str) else model
    sigma = np.asarray(stress, dtype=float)
    eta = np.asarray(viscosity, dtype=float)
    if sigma.shape != eta.shape:
        raise ValueError(f"shape mismatch: stress {sigma.shape} vs viscosity {eta.shape}")

    if branch is not None:
        sigma, eta = sigma[branch], eta[branch]

    good = np.isfinite(sigma) & np.isfinite(eta) & (sigma > 0) & (eta > 0)
    sigma, eta = sigma[good], eta[good]
    if sigma.size <= mdl.n_params:
        raise ValueError(
            f"{mdl.name}: {sigma.size} usable points for {mdl.n_params} parameters; "
            "refusing to fit an underdetermined system"
        )

    log_eta = np.log(eta)
    residual = _residual_factory(mdl, sigma, log_eta)
    lo, hi = np.array(mdl.bounds[0]), np.array(mdl.bounds[1])

    rng = np.random.default_rng(seed)
    starts = [] if start is None else [np.asarray(start, dtype=float)]
    starts.append(_initial_guess(mdl, sigma, eta))
    for _ in range(max(0, n_starts - 1)):
        u = rng.uniform(size=mdl.n_params)
        starts.append(np.exp(np.log(np.maximum(lo, 1e-12)) + u * (np.log(hi) - np.log(np.maximum(lo, 1e-12)))))

    best = None
    for x0 in starts:
        try:
            sol = least_squares(
                residual, np.clip(x0, lo + 1e-12, hi - 1e-12),
                bounds=(lo, hi), method="trf", max_nfev=5000,
            )
        except (ValueError, np.linalg.LinAlgError):
            continue
        if best is None or sol.cost < best.cost:
            best = sol

    if best is None:
        raise RuntimeError(f"{mdl.name}: every start failed to converge")

    cov, corr = _covariance(best.jac, best.fun, mdl.n_params)
    return FitResult(
        model=mdl,
        params=best.x,
        stderr=np.sqrt(np.diag(cov)),
        cov=cov,
        corr=corr,
        residuals=best.fun,
        n_points=int(sigma.size),
        success=bool(best.success),
        message=str(best.message),
        branch=branch,
        metadata={"n_starts": len(starts), "seed": seed},
    )


@dataclass
class DivergenceFit:
    """One fit of eta_r = A (1 - phi/phi_c)^-n to a branch measured in phi."""

    A: float
    phi_c: float
    n: float
    stderr: dict[str, float]
    residuals: np.ndarray
    n_points: int
    n_free: int
    free_exponent: bool

    @property
    def rmse_log(self) -> float:
        return float(np.sqrt(np.mean(self.residuals**2)))

    def summary(self) -> str:
        e = self.stderr
        exp = f"{self.n:.3f} +/- {e['n']:.3f}" if self.free_exponent else f"{self.n:.3f} (fixed)"
        return (
            f"eta_r = A (1 - phi/phi_c)^-n   n={self.n_points} k={self.n_free} "
            f"rmse(log eta)={self.rmse_log:.4f}\n"
            f"    A     = {self.A:.4f} +/- {e['A']:.4f}\n"
            f"    phi_c = {self.phi_c:.4f} +/- {e['phi_c']:.4f}\n"
            f"    n     = {exp}"
        )


def fit_divergence(
    phi: np.ndarray,
    eta_rel: np.ndarray,
    *,
    free_exponent: bool = True,
    n: float = 2.0,
    n_starts: int = 12,
    seed: int = 0,
) -> DivergenceFit:
    """Fit a jamming divergence measured as relative viscosity against phi.

    This is the form published papers quote their ``phi_m`` from: a whole branch
    of a flow-curve series collapsed to one viscosity per volume fraction, then
    fit with ``eta_r = A (1 - phi/phi_c)^-n``. It is *not* the stress-space
    Wyart-Cates fit the rest of the package does, and the two are different
    measurements of the same quantity, which is exactly why comparing them is
    worth doing.

    Residuals are in log relative viscosity, as everywhere else here, so an
    rmse is directly comparable with :class:`FitResult.rmse_log`.

    ``free_exponent=False`` pins ``n``; a ``phi_c`` obtained at one exponent
    cannot be compared with one obtained at another, so the exponent used is
    recorded on the result rather than assumed by the reader.
    """
    x = np.asarray(phi, dtype=float)
    y = np.asarray(eta_rel, dtype=float)
    if x.shape != y.shape:
        raise ValueError(f"shape mismatch: phi {x.shape} vs eta_rel {y.shape}")
    good = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x, y = x[good], y[good]
    order = np.argsort(x)
    x, y = x[order], y[order]

    n_free = 3 if free_exponent else 2
    if x.size <= n_free:
        raise ValueError(
            f"fit_divergence: {x.size} usable points for {n_free} parameters; "
            "refusing to fit an underdetermined system"
        )

    log_y = np.log(y)
    # phi_c must exceed every measured phi, else the branch contains its own pole.
    lo = np.array([1e-4, float(x.max()) + 1e-4] + ([0.5] if free_exponent else []))
    hi = np.array([1e3, 0.95] + ([6.0] if free_exponent else []))
    if lo[1] >= hi[1]:
        raise ValueError(
            f"fit_divergence: largest phi = {x.max():.4f} leaves no room for phi_c below "
            f"{hi[1]:g}"
        )

    def residual(p):
        A, phi_c = p[0], p[1]
        exponent = p[2] if free_exponent else n
        with np.errstate(all="ignore"):
            pred = A * (1.0 - x / phi_c) ** -exponent
        out = np.full(pred.shape, _PENALTY)
        ok = np.isfinite(pred) & (pred > 0)
        out[ok] = np.log(pred[ok]) - log_y[ok]
        return out

    rng = np.random.default_rng(seed)
    starts = [np.array([1.0, min(float(x.max()) * 1.05, 0.94)] + ([n] if free_exponent else []))]
    for _ in range(max(0, n_starts - 1)):
        starts.append(lo + rng.uniform(size=n_free) * (hi - lo))

    best = None
    for x0 in starts:
        try:
            sol = least_squares(residual, np.clip(x0, lo + 1e-12, hi - 1e-12),
                                bounds=(lo, hi), method="trf", max_nfev=5000)
        except (ValueError, np.linalg.LinAlgError):
            continue
        if best is None or sol.cost < best.cost:
            best = sol
    if best is None:
        raise RuntimeError("fit_divergence: every start failed to converge")

    cov, _ = _covariance(best.jac, best.fun, n_free)
    err = np.sqrt(np.abs(np.diag(cov)))
    names = ["A", "phi_c"] + (["n"] if free_exponent else [])
    stderr = {k: float(v) for k, v in zip(names, err)}
    stderr.setdefault("n", 0.0)
    return DivergenceFit(
        A=float(best.x[0]),
        phi_c=float(best.x[1]),
        n=float(best.x[2]) if free_exponent else float(n),
        stderr=stderr,
        residuals=best.fun,
        n_points=int(x.size),
        n_free=n_free,
        free_exponent=free_exponent,
    )


def synthetic_curve(
    model: Model | str,
    params,
    sigma: np.ndarray,
    *,
    noise: float = 0.0,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate (stress, viscosity) from a model with multiplicative lognormal noise.

    ``noise`` is the standard deviation in log viscosity, so 0.02 is roughly 2%
    relative scatter. Points where the model returns a non-finite viscosity
    (jammed or unyielded) are dropped rather than clipped.
    """
    mdl: Model = get_model(model) if isinstance(model, str) else model
    mdl.check_params(params)
    sigma = np.asarray(sigma, dtype=float)
    eta = np.asarray(mdl(sigma, *params), dtype=float)
    keep = np.isfinite(eta) & (eta > 0)
    sigma, eta = sigma[keep], eta[keep]
    if noise > 0:
        rng = np.random.default_rng(seed)
        eta = eta * np.exp(rng.normal(0.0, noise, size=eta.shape))
    return sigma, eta
