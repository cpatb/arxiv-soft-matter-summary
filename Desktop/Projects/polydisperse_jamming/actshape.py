"""actshape.py -- the ACTIVATION SHAPE of Wyart-Cates, read off sheared contact forces.

Pre-registration: PREREG_ACTSHAPE_2026-09-22.md. Grader: grade_actshape.py.
Pure functions only; runs on Windows without Yade.

THE CHAIN, in one place
-----------------------
Wyart-Cates calls a contact frictional when its normal force exceeds a repulsive
threshold F_rep, and fhat is the FRACTION of such contacts. For the surface forces a
dense colloid carries (screened electrostatics, steric brushes, DLVO), the Derjaguin
approximation gives a threshold proportional to the pair's REDUCED RADIUS,

    F_rep = kappa R*,        R* = R1 R2 / (R1 + R2),                         (1)

with kappa a force per length set by the potential alone. So contact c is frictional
iff f_c / R*_c > kappa. Writing the DEM force in units of the applied stress,

    g_c = f_c / (R*_c sigma)          (dimensionless; length unit d10 = 1),     (2)

    fhat(sigma) = P( g > kappa / sigma ).                                      (3)

Everything the solver needs is in the distribution of g:

  * the ONSET scale: sigma* := the stress at which fhat = 1/e, i.e. the solver's own
    definition (fss = exp(-(sigma*/sigma)^delta) equals 1/e at sigma = sigma*):
        sigma* = kappa / Q_{1-1/e}(g),                                         (4)
    Q_p the p-quantile. At fixed kappa and d10, Q_{1-1/e}(g) IS the size-distribution
    dependence of sigma*.
  * the SHAPE: for a Weibull survival exp(-(u/lambda)^delta) the quantiles satisfy
        delta = ln[ ln(1-p_hi) / ln(1-p_lo) ] / ln(Q_hi/Q_lo),                 (5)
    so delta is read from two quantiles with no fit. Two pairs (10/90, 25/75) agreeing
    is the test that the Weibull form -- the solver's form -- describes the data at all.
  * the JAMMING STRESS: in the reduced framework (three_knobs.py) a suspension at
    x = phi/phi0 inside the window chi <= x < 1 jams once fhat >= f_c, with
        f_c = (1 - x) / (1 - chi).                                             (6)
    From (3),  sigma_jam = kappa / Q_{1-f_c}(g).                               (7)
    (7) needs NO convention for sigma* and NO Weibull assumption: at fixed potential
    and size it compares two size distributions directly. It is the primary output.

⚠ A CONVENTION TRAP, stated so nobody quotes it as physics: if every distribution is
first normalised by its own Q_{1-1/e} (i.e. compared at MATCHED sigma*), all activation
curves cross at fhat = 1/e BY CONSTRUCTION, and "a flatter tail helps jamming above
fhat = 1/e and hurts below" is then a tautology of the normalisation. The comparison at
matched KAPPA, (7), is the one that is not.

GAUGE. sigma in (2) is the mean PRESSURE p of the sheared pack. The solver's sigma is a
shear-stress magnitude; in the steady state tau = mu_M p with mu_M the macroscopic
friction, which differs between the mu_p = 0 and 0.3 arms and moves with CV. Both
gauges are reported (`gauge="p"` / `"tau"`); neither is silently assumed.
"""

from __future__ import annotations

import math

import numpy as np

# Quantile pairs for (5), and the Weibull constant each implies.
Q_PAIRS = ((0.10, 0.90), (0.25, 0.75))
INV_E = math.exp(-1.0)


def weibull_const(p_lo: float, p_hi: float) -> float:
    """ln[ln(1-p_hi)/ln(1-p_lo)]: the numerator of (5). 3.0843 for (0.1, 0.9)."""
    return math.log(math.log(1.0 - p_hi) / math.log(1.0 - p_lo))


def delta_from_quantiles(u: np.ndarray, p_lo: float, p_hi: float) -> float:
    """Weibull shape from two quantiles, eq. (5). Fit-free; scale-invariant in u."""
    u = np.asarray(u, dtype=float)
    u = u[np.isfinite(u) & (u > 0)]
    if u.size < 50:
        return float("nan")
    q_lo, q_hi = np.quantile(u, [p_lo, p_hi])
    return weibull_const(p_lo, p_hi) / math.log(q_hi / q_lo)


def reduced_radius(r1: np.ndarray, r2: np.ndarray) -> np.ndarray:
    r1 = np.asarray(r1, dtype=float)
    r2 = np.asarray(r2, dtype=float)
    return r1 * r2 / (r1 + r2)


def g_values(fn, r1, r2, sigma: float, norm: str = "derjaguin") -> np.ndarray:
    """Per-contact g, eq. (2), for a threshold law named by `norm`.

    "derjaguin": F_rep proportional to R*          -- the physical case, eq. (1)
    "constant" : F_rep the same for every contact  -- size-blind; what f/<f> assumes
    "area"     : F_rep proportional to R*^2        -- a DIAGNOSTIC bracket only; no
                 surface force this project maps onto scales this way
    Contacts with fn <= 0 are dropped (they carry no load and cannot be frictional).
    """
    fn = np.asarray(fn, dtype=float)
    rs = reduced_radius(r1, r2)
    ok = np.isfinite(fn) & (fn > 0)
    if norm == "derjaguin":
        g = fn[ok] / rs[ok]
    elif norm == "constant":
        g = fn[ok]
    elif norm == "area":
        g = fn[ok] / rs[ok] ** 2
    else:
        raise ValueError(f"unknown norm {norm!r}")
    return g / float(sigma)


def f_crit(x: float, chi: float) -> float:
    """Eq. (6): the frictional fraction at which phi_J(fhat) reaches phi.

    Defined only inside the window chi <= x < 1; outside it the suspension either can
    never jam (x < chi, returns inf -- no fraction <= 1 suffices) or is jammed at rest
    (x >= 1, returns 0).
    """
    if not (0.0 < chi < 1.0):
        raise ValueError(f"chi must lie in (0, 1), got {chi}")
    if x >= 1.0:
        return 0.0
    if x < chi:
        return math.inf
    return (1.0 - x) / (1.0 - chi)


def sigma_jam_over_kappa(g: np.ndarray, fc: float) -> float:
    """Eq. (7): sigma_jam / kappa = 1 / Q_{1-fc}(g). inf when fc >= 1 is unreachable."""
    if fc <= 0.0:
        return 0.0
    if fc >= 1.0:
        return math.inf
    g = np.asarray(g, dtype=float)
    return 1.0 / float(np.quantile(g, 1.0 - fc))


def sigma_star_over_kappa(g: np.ndarray) -> float:
    """Eq. (4): sigma*/kappa = 1/Q_{1-1/e}(g), the solver's own sigma* convention."""
    return sigma_jam_over_kappa(g, INV_E)


def weibull_sigma_jam_ratio(fc: float, delta: float) -> float:
    """sigma_jam/sigma* for the solver's Weibull activation: ln(1/fc)^(-1/delta)."""
    if fc <= 0.0:
        return 0.0
    if fc >= 1.0:
        return math.inf
    return math.log(1.0 / fc) ** (-1.0 / delta)


def eta_rel(f, x: float, chi: float, eps_min: float = 0.01, mp_expo: float = 2.0):
    """eta/eta_s exactly as wc_drop_impact.c eta_wc(): eps^-MP_EXPO, eps = 1 - x/(1-(1-chi)f),
    floored at eps_min (EPS_SOFT, off by default in the solver, is not modelled)."""
    f = np.asarray(f, dtype=float)
    eps = np.maximum(1.0 - x / (1.0 - (1.0 - chi) * f), eps_min)
    return eps ** -mp_expo


def rate_threshold(stress_at, x: float, chi: float, eps_min: float = 0.01,
                   mp_expo: float = 2.0, n: int = 200_001) -> tuple[float, float]:
    """RATE-driven jamming threshold, the solver's own drive (sigma = eta(fhat) gd).

    `stress_at(f)` is the stress at which the activation reaches fraction f (vectorised):
    sigma*.(ln 1/f)^(-1/delta) for the Weibull law, or kappa/Q_{1-f}(g) for a measured g.
    A point relaxing from fhat = 0 reaches f_c iff eta(f) gd > stress_at(f) on all of
    [0, f_c), so
        gd_c * eta_s = max_{0<f<f_c} stress_at(f) / eta_rel(f).
    Returns (gd_c * eta_s, f*), f* the bottleneck. derivations/jam_rate_driven_2026-09-22.py.
    """
    fc = f_crit(x, chi)
    if not (0.0 < fc < 1.0):
        raise ValueError(f"rate threshold defined only inside the window; f_c = {fc}")
    f = np.linspace(1e-6, fc * (1 - 1e-9), n)
    r = np.asarray(stress_at(f)) / eta_rel(f, x, chi, eps_min, mp_expo)
    i = int(np.argmax(r))
    return float(r[i]), float(f[i])


def activation_curve(g: np.ndarray, s: np.ndarray) -> np.ndarray:
    """Empirical fhat at sigma/kappa = s: P(g > 1/s), eq. (3)."""
    g = np.sort(np.asarray(g, dtype=float))
    thr = 1.0 / np.asarray(s, dtype=float)
    return 1.0 - np.searchsorted(g, thr, side="right") / g.size
