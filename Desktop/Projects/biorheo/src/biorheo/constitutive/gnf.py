"""Generalised-Newtonian viscosity models, eta(gammadot), as ``RateModel`` instances.

Every model here is single-valued in shear rate and is written exactly as the solver
bridge will export it. Units: eta in Pa s, gammadot in 1/s, stresses in Pa. Each
model's ``limits`` are derived in ``derivations/gnf_limits.md`` and exercised
numerically in ``tests/test_gnf.py``; the numbers in the tests were computed by hand
from the formula, not by running this file.

Parameter bounds are wide, physical boxes: they exist so the fitter cannot wander
into negative viscosities, not to encode a prior. A parameter sitting on one of them
after a fit is a ``bounded_only`` verdict in ``biorheo.identify``, never a value.
"""

from __future__ import annotations

import numpy as np

from .base import RateModel

__all__ = [
    "NEWTONIAN",
    "POWER_LAW",
    "HERSCHEL_BULKLEY",
    "CASSON",
    "CARREAU_YASUDA",
    "CROSS",
    "MODELS",
    "get_model",
]

# n must stay below 1 for "thinning", but the box allows mild thickening so that a fit
# can *report* n > 1 rather than being pinned at 1 by construction.
_N_BOUNDS = (0.05, 2.0)


# --- Newtonian ---------------------------------------------------------------------

def _newtonian(gd: np.ndarray, mu: float) -> np.ndarray:
    return np.full_like(gd, mu, dtype=float)


def _guess_newtonian(gd, eta):
    return np.array([np.median(eta)])


NEWTONIAN = RateModel(
    name="newtonian",
    param_names=("mu",),
    func=_newtonian,
    bounds=((1e-6,), (1e6,)),
    describes_thickening=False,
    notes="eta = mu; the null model every other one must beat",
    guess_fn=_guess_newtonian,
    limits={"all": "eta = mu at every rate"},
    units={"mu": "Pa*s"},
)


# --- power law ---------------------------------------------------------------------

def _power_law(gd: np.ndarray, k: float, n: float) -> np.ndarray:
    return k * gd ** (n - 1.0)


def _guess_power_law(gd, eta):
    # slope of log eta against log gd is n - 1; np.polyfit is the least-squares line
    lg, le = np.log(gd), np.log(eta)
    slope, intercept = np.polyfit(lg, le, 1) if lg.size >= 2 else (0.0, le.mean())
    return np.array([np.exp(intercept), slope + 1.0])


POWER_LAW = RateModel(
    name="power_law",
    param_names=("k", "n"),
    func=_power_law,
    bounds=((1e-6, _N_BOUNDS[0]), (1e6, _N_BOUNDS[1])),
    describes_thickening=False,
    notes="eta = k gd^(n-1); diverges as gd -> 0 for n < 1, so never a low-rate model",
    guess_fn=_guess_power_law,
    limits={
        "n = 1": "eta = k, Newtonian",
        "gd -> 0, n < 1": "eta -> inf as gd^(n-1)",
        "gd -> inf, n < 1": "eta -> 0 as gd^(n-1)",
    },
    units={"k": "Pa*s**n", "n": "dimensionless"},
)


# --- Herschel-Bulkley --------------------------------------------------------------

def _herschel_bulkley(gd: np.ndarray, tau_y: float, k: float, n: float) -> np.ndarray:
    # sigma = tau_y + k gd^n  =>  eta = tau_y/gd + k gd^(n-1); inf at gd = 0.
    with np.errstate(divide="ignore"):
        return tau_y / gd + k * gd ** (n - 1.0)


def _guess_herschel_bulkley(gd, eta):
    sigma = eta * gd
    tau_y = 0.5 * np.min(sigma)
    n = 0.5
    k = np.median((sigma - tau_y) / gd ** n)
    return np.array([tau_y, max(k, 1e-6), n])


HERSCHEL_BULKLEY = RateModel(
    name="herschel_bulkley",
    param_names=("tau_y", "k", "n"),
    func=_herschel_bulkley,
    bounds=((0.0, 1e-6, _N_BOUNDS[0]), (1e5, 1e6, _N_BOUNDS[1])),
    describes_thickening=False,
    notes="sigma = tau_y + k gd^n in rate space; tau_y = 0 is the power law, n = 1 is Bingham",
    guess_fn=_guess_herschel_bulkley,
    limits={
        "tau_y = 0": "power law, eta = k gd^(n-1)",
        "n = 1": "Bingham, eta = tau_y/gd + k",
        "gd -> 0": "sigma -> tau_y; eta ~ tau_y/gd -> inf",
        "gd -> inf, n < 1": "sigma ~ k gd^n; the yield term is a vanishing fraction tau_y/(k gd^n)",
    },
    units={"tau_y": "Pa", "k": "Pa*s**n", "n": "dimensionless"},
)


# --- Casson -------------------------------------------------------------------------

def _casson(gd: np.ndarray, tau_y: float, eta_c: float) -> np.ndarray:
    # sqrt(sigma) = sqrt(tau_y) + sqrt(eta_c gd)  =>  eta = (sqrt(tau_y/gd) + sqrt(eta_c))^2
    with np.errstate(divide="ignore"):
        return (np.sqrt(tau_y / gd) + np.sqrt(eta_c)) ** 2


def _guess_casson(gd, eta):
    # linear in sqrt space: sqrt(sigma) = sqrt(tau_y) + sqrt(eta_c) sqrt(gd)
    s, x = np.sqrt(eta * gd), np.sqrt(gd)
    slope, intercept = np.polyfit(x, s, 1) if x.size >= 2 else (np.sqrt(np.median(eta)), 0.0)
    return np.array([max(intercept, 0.0) ** 2, max(slope, 1e-3) ** 2])


CASSON = RateModel(
    name="casson",
    param_names=("tau_y", "eta_c"),
    func=_casson,
    bounds=((0.0, 1e-6), (1e5, 1e6)),
    describes_thickening=False,
    notes="sqrt(sigma) = sqrt(tau_y) + sqrt(eta_c gd); the classical whole-blood form",
    guess_fn=_guess_casson,
    limits={
        "tau_y = 0": "Newtonian, eta = eta_c",
        "gd -> 0": "sigma -> tau_y; eta ~ tau_y/gd -> inf",
        "gd -> inf": "eta -> eta_c from above, as eta_c + 2 sqrt(tau_y eta_c/gd) + tau_y/gd",
    },
    units={"tau_y": "Pa", "eta_c": "Pa*s"},
)


# --- Carreau-Yasuda ----------------------------------------------------------------

def _carreau_yasuda(gd: np.ndarray, eta_0: float, eta_inf: float, lam: float, a: float,
                    n: float) -> np.ndarray:
    return eta_inf + (eta_0 - eta_inf) * (1.0 + (lam * gd) ** a) ** ((n - 1.0) / a)


def _guess_carreau_yasuda(gd, eta):
    eta_0, eta_inf = np.max(eta), np.min(eta)
    # the transition sits where eta has fallen to the geometric mean of the plateaux
    mid = np.sqrt(eta_0 * max(eta_inf, 1e-12))
    i = int(np.argmin(np.abs(eta - mid)))
    lam = 1.0 / gd[i]
    return np.array([eta_0, eta_inf, lam, 2.0, 0.5])


CARREAU_YASUDA = RateModel(
    name="carreau_yasuda",
    param_names=("eta_0", "eta_inf", "lam", "a", "n"),
    func=_carreau_yasuda,
    bounds=((1e-6, 0.0, 1e-6, 0.1, _N_BOUNDS[0]), (1e6, 1e6, 1e6, 10.0, _N_BOUNDS[1])),
    describes_thickening=False,
    notes="eta = eta_inf + (eta_0 - eta_inf)[1 + (lam gd)^a]^((n-1)/a); a = 2 is Carreau",
    guess_fn=_guess_carreau_yasuda,
    limits={
        "gd -> 0": "eta -> eta_0",
        "gd -> inf": "eta -> eta_inf",
        "lam gd >> 1, eta_inf = 0": "power law with k = eta_0 lam^(n-1)",
        "n = 1": "eta = eta_0 at every rate (eta_inf drops out)",
    },
    units={"eta_0": "Pa*s", "eta_inf": "Pa*s", "lam": "s", "a": "dimensionless",
           "n": "dimensionless"},
)


# --- Cross --------------------------------------------------------------------------

def _cross(gd: np.ndarray, eta_0: float, eta_inf: float, lam: float, m: float) -> np.ndarray:
    return eta_inf + (eta_0 - eta_inf) / (1.0 + (lam * gd) ** m)


def _guess_cross(gd, eta):
    eta_0, eta_inf = np.max(eta), np.min(eta)
    mid = 0.5 * (eta_0 + eta_inf)
    i = int(np.argmin(np.abs(eta - mid)))
    return np.array([eta_0, eta_inf, 1.0 / gd[i], 0.7])


CROSS = RateModel(
    name="cross",
    param_names=("eta_0", "eta_inf", "lam", "m"),
    func=_cross,
    bounds=((1e-6, 0.0, 1e-6, 0.05), (1e6, 1e6, 1e6, 3.0)),
    describes_thickening=False,
    notes="eta = eta_inf + (eta_0 - eta_inf)/(1 + (lam gd)^m)",
    guess_fn=_guess_cross,
    limits={
        "gd -> 0": "eta -> eta_0",
        "gd -> inf": "eta -> eta_inf",
        "lam gd = 1": "eta = eta_inf + (eta_0 - eta_inf)/2, the half-way rate",
        "lam gd >> 1, eta_inf = 0": "power law with k = eta_0 lam^(-m), n = 1 - m",
    },
    units={"eta_0": "Pa*s", "eta_inf": "Pa*s", "lam": "s", "m": "dimensionless"},
)


from .structural import STRUCTURAL_STEADY  # noqa: E402  (steady state of the kinetics)

MODELS: dict[str, RateModel] = {
    m.name: m
    for m in (NEWTONIAN, POWER_LAW, HERSCHEL_BULKLEY, CASSON, CARREAU_YASUDA, CROSS,
              STRUCTURAL_STEADY)
}


def get_model(name: str) -> RateModel:
    """Look a model up by its registry name; raises ``KeyError`` naming the options."""
    try:
        return MODELS[name]
    except KeyError:
        raise KeyError(f"unknown model {name!r}; known: {sorted(MODELS)}") from None
