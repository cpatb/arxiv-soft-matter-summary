"""
Swelling / drying kinetics -> collective diffusivity D, surface transfer k_m, Bi.

Implements the model-free extractions from the hydrogel MOC section 5A for a
gravimetric curve F(t) = M_t / M_inf of a plane sheet (slab) of full thickness L,
losing/absorbing solvent through both faces (symmetric, 1-D).

Three model-free routes (do all; agreement is the Fickian self-consistency check):

  1. Early-time sqrt(t)  [F <~ 0.6]:   F = (4/L) sqrt(D t / pi)
        -> D = pi (slope * L / 4)^2                     (Crank 1975, eq. 4.20)
  2. Long-time single-exponential:      F = 1 - (8/pi^2) exp(-pi^2 D t / L^2)
        -> ln(1 - F) = const - (pi^2 D / L^2) t         (Tanaka-Fillmore 1979)
  3. Power law  F ~ t^n over the early rise:
        n ~ 0.5 Fickian, n ~ 1 Case II (moving front), intermediate = anomalous.

For drying in air the surface flux is finite (Robin / evaporative BC); a naive
sqrt(t) fit then returns too-small D. `fit_D_km_robin` fits (D, k_m) together by
forward-simulating the 1-D diffusion equation with a Robin boundary condition and
least-squares matching F(t); it returns Bi = k_m (L/2) / D as a byproduct (the
Biot number on the diffusion length = half-thickness). The Robin simulator uses
the SAME full-thickness-L, two-face convention as Routes 1-2, so the D it returns
is directly comparable to theirs (this consistency is asserted in validate.py).

Assumptions: 1-D through-thickness transport, sheet thin vs lateral size (or rim
sealed), constant D over the fitted range (a single number is an average over the
phi-range traversed -- report early and long-time D as a bracket if they differ).

References: Crank, *The Mathematics of Diffusion* (1975); Tanaka & Fillmore,
J. Chem. Phys. 70, 1214 (1979).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares


# --------------------------------------------------------------------------
# Analytic Crank series (used both as a forward model and to validate the fits)
# --------------------------------------------------------------------------
def crank_fixed_surface(t, D, L, n_terms: int = 200):
    """F(t)=M_t/M_inf for a plane sheet of thickness L, fixed surface conc,
    diffusion through both faces (Crank 1975, eq. 4.18).

    F = 1 - sum_{n=0}^inf 8/((2n+1)^2 pi^2) exp(-(2n+1)^2 pi^2 D t / L^2)
    """
    t = np.atleast_1d(np.asarray(t, dtype=float))
    if D <= 0 or L <= 0:
        raise ValueError("D and L must be positive")
    n = np.arange(n_terms)
    coef = 8.0 / ((2 * n + 1) ** 2 * np.pi**2)
    expo = -((2 * n + 1) ** 2) * np.pi**2 * D * t[:, None] / L**2
    F = 1.0 - (coef[None, :] * np.exp(expo)).sum(axis=1)
    return np.clip(F, 0.0, 1.0)


# --------------------------------------------------------------------------
# Route 1: early-time sqrt(t)
# --------------------------------------------------------------------------
@dataclass
class SqrtFit:
    D: float
    slope: float          # slope of F vs sqrt(t)
    r2: float
    n_points: int


def fit_D_early_sqrt(t, F, L, f_max: float = 0.6) -> SqrtFit:
    """Fit D from the early-time linear F-vs-sqrt(t) region (F <= f_max)."""
    t = np.asarray(t, dtype=float)
    F = np.asarray(F, dtype=float)
    if L <= 0:
        raise ValueError("L must be positive")
    mask = (F <= f_max) & (t > 0) & np.isfinite(F)
    if mask.sum() < 2:
        raise ValueError(f"need >=2 points with F<={f_max}; got {mask.sum()}")
    x = np.sqrt(t[mask])
    y = F[mask]
    # force through origin? No -- allow intercept, but report slope only.
    A = np.vstack([x, np.ones_like(x)]).T
    (slope, intercept), *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ np.array([slope, intercept])
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    D = np.pi * (slope * L / 4.0) ** 2
    return SqrtFit(D=float(D), slope=float(slope), r2=r2, n_points=int(mask.sum()))


# --------------------------------------------------------------------------
# Route 2: long-time exponential
# --------------------------------------------------------------------------
@dataclass
class LongTimeFit:
    D: float
    slope: float          # slope of ln(1-F) vs t  (= -pi^2 D / L^2)
    r2: float
    n_points: int


def fit_D_longtime(t, F, L, f_min: float = 0.6, f_max: float = 0.999) -> LongTimeFit:
    """Fit D from the long-time slope of ln(1 - F) vs t."""
    t = np.asarray(t, dtype=float)
    F = np.asarray(F, dtype=float)
    if L <= 0:
        raise ValueError("L must be positive")
    mask = (F >= f_min) & (F < f_max) & np.isfinite(F)
    if mask.sum() < 2:
        raise ValueError(f"need >=2 points with {f_min}<=F<{f_max}; got {mask.sum()}")
    x = t[mask]
    y = np.log(1.0 - F[mask])
    A = np.vstack([x, np.ones_like(x)]).T
    (slope, intercept), *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ np.array([slope, intercept])
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    D = -slope * L**2 / np.pi**2
    return LongTimeFit(D=float(D), slope=float(slope), r2=r2, n_points=int(mask.sum()))


# --------------------------------------------------------------------------
# Route 3: power-law exponent (Fickian vs Case II diagnostic)
# --------------------------------------------------------------------------
def fit_exponent_n(t, F, f_max: float = 0.6) -> float:
    """Fit n in F ~ k t^n over the early rise (log-log slope). ~0.5 Fickian,
    ~1.0 Case II."""
    t = np.asarray(t, dtype=float)
    F = np.asarray(F, dtype=float)
    mask = (F > 0) & (F <= f_max) & (t > 0) & np.isfinite(F)
    if mask.sum() < 2:
        raise ValueError(f"need >=2 points with 0<F<={f_max}; got {mask.sum()}")
    x = np.log(t[mask])
    y = np.log(F[mask])
    n = np.polyfit(x, y, 1)[0]
    return float(n)


# --------------------------------------------------------------------------
# Route 4: (D, k_m) from a Robin-BC forward model -> Bi
# --------------------------------------------------------------------------
def _simulate_robin(t_eval, D, k_m, L, n_nodes: int = 41):
    """F(t) for a plane sheet of FULL thickness L drying with a finite surface
    flux, matching the two-face full-thickness-L convention of `crank_fixed_surface`
    and Routes 1-2.

    By symmetry the two-face sheet is modelled on the half-domain x in [0, L/2]
    with a no-flux (symmetry) plane at the mid-plane x=0 and a Robin
    (evaporative) BC at the surface x=L/2:
        -D dC/dx|_{L/2} = k_m (C_surf - C_env),  C_env = 0.
    (This half-domain problem is also *exactly* the one-sided drying of a
    thickness-L/2 gel with a sealed base -- so for a physically one-sided run of
    thickness H, pass L = 2H.)

    In the limit k_m -> inf the surface concentration is pinned to 0 and this
    reduces to `crank_fixed_surface(t, D, L)` (asserted in validate.py).

    Returns F(t) = fractional mass loss. Because the initial field is uniform
    C0 = 1, <C>(0) = 1 exactly, so F = 1 - <C>(t) -- referenced to t=0, not to
    the first evaluated time (which may be > 0 for real gravimetric series).
    Method of lines with a second-order finite-difference Laplacian.
    """
    half = 0.5 * L
    x = np.linspace(0.0, half, n_nodes)
    dx = x[1] - x[0]
    C0 = np.ones(n_nodes)

    def rhs(_t, C):
        d2 = np.empty_like(C)
        d2[1:-1] = (C[2:] - 2 * C[1:-1] + C[:-2]) / dx**2
        # symmetry (no-flux) at x=0: ghost node C[-1]=C[1] -> d2[0]=2(C1-C0)/dx^2
        d2[0] = 2.0 * (C[1] - C[0]) / dx**2
        # Robin at surface x=L/2: ghost C_{N} = C_{N-2} - 2 dx k_m/D C_{N-1}
        cghost = C[-2] - 2.0 * dx * k_m / D * C[-1]
        d2[-1] = (cghost - 2 * C[-1] + C[-2]) / dx**2
        return D * d2

    sol = solve_ivp(rhs, (0.0, float(t_eval[-1])), C0, t_eval=t_eval,
                    method="BDF", rtol=1e-8, atol=1e-10)
    if not sol.success:
        raise RuntimeError(f"solve_ivp failed: {sol.message}")
    mean_C = np.trapezoid(sol.y, x, axis=0) / half   # <C> over the half-domain
    F = 1.0 - mean_C                                 # <C>(0) = 1 exactly
    return F


@dataclass
class RobinFit:
    D: float
    k_m: float
    Bi: float
    cost: float


def fit_D_km_robin(t, F, L, D0=None, km0=None, n_nodes: int = 41) -> RobinFit:
    """Jointly fit (D, k_m) to a drying-in-air curve; returns Bi = k_m L / D.

    Uses an early-time sqrt fit for D0 if not supplied (that D0 is biased low by
    the finite BC, but is a fine starting guess). Fits in log-space to keep both
    parameters positive.
    """
    t = np.asarray(t, dtype=float)
    F = np.asarray(F, dtype=float)
    good = np.isfinite(F) & (t >= 0)
    t, F = t[good], F[good]
    if D0 is None:
        try:
            D0 = fit_D_early_sqrt(t, F, L).D
        except ValueError:
            D0 = L**2 / (t[-1] if t[-1] > 0 else 1.0)
    if km0 is None:
        km0 = D0 / L  # Bi ~ 1 starting guess

    def resid(logp):
        D, k_m = np.exp(logp)
        Fmodel = _simulate_robin(t, D, k_m, L, n_nodes=n_nodes)
        return Fmodel - F

    # Bounded trust-region fit in log-space (LM is unbounded and can wander to
    # unphysical D on the shallow-gradient drying curve). Bound each parameter to
    # a wide but physical window around the starting guess.
    p0 = np.log([D0, km0])
    lo = p0 - np.log(1e3)
    hi = p0 + np.log(1e3)
    res = least_squares(resid, p0, method="trf", bounds=(lo, hi),
                        x_scale="jac", max_nfev=400)
    D, k_m = np.exp(res.x)
    # Bi uses the diffusion length = half-thickness L/2 (mid-plane to surface),
    # the physical transport length of the symmetric two-face sheet.
    Bi = k_m * (0.5 * L) / D
    return RobinFit(D=float(D), k_m=float(k_m), Bi=float(Bi), cost=float(res.cost))
