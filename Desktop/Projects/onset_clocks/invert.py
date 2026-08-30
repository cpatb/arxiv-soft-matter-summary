"""Ask the data what the closure would need, instead of asking whether it works.

The audit in ``run_audit.py`` tests candidate criteria and eliminates them. This
module does the complementary thing: it takes the closure the drop-impact repo
already fitted, treats ``phi_m`` as the one unknown, and asks what value each
drop would need for the model to reproduce that drop's own arrested footprint.

Why this is the sharper question
--------------------------------
``phi_m`` is supposed to be a **material constant**. It cannot depend on how
fast you threw the drop. So:

- if the required ``phi_m`` comes back constant, the closure is fine and the
  fitted 0.5176 was simply wrong -- the contradiction dissolves and the
  inversion has *measured* the right value;
- if the required ``phi_m`` drifts with ``U0`` at fixed packing and fixed salt,
  no value of it can work, and the closure is structurally wrong for impact.

That is a cleaner test than any correlation, and it needs no bench work.

What this is NOT
----------------
The implied viscosity below is a **reparameterisation of** ``beta_p``, not an
independent measurement of anything. Feeding it back through the model cannot
produce information the footprint did not already contain. The content of this
module is entirely in whether the *required* ``phi_m`` holds still -- a
statement about the model's internal consistency, not about the data.

The one modelling assumption, stated once
-----------------------------------------
Arrest is taken to be viscous: the drop stops when its kinetic energy has gone
into viscous dissipation, which gives the standard viscous maximum-spread
balance

    beta ~ Re^(1/p),   Re = rho U0 D0 / eta   =>   eta = rho U0 D0 / beta^p

with ``p = 5`` (Eggers et al. 2010, the viscous branch) as the primary and
``p = 4`` reported as a sensitivity. **The O(1) prefactor is unknown**, so
absolute viscosities carry it and only ratios and trends are quoted. Every
conclusion below is checked at both exponents.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from clocks import Drop, frictional_fraction_shear, material  # noqa: E402

# Exponent in beta ~ Re^(1/p). 5 is the viscous branch of the Eggers/Laan
# collapse; 4 is the cruder energy balance. Both are run.
SPREAD_EXPONENT_PRIMARY = 5.0
SPREAD_EXPONENTS = (4.0, 5.0)


@dataclass(frozen=True)
class Inversion:
    """What one drop's footprint demands of the closure."""

    label: str
    phi: float
    salt_mM: float
    u0: float
    beta_p: float
    eta_implied: float      # Pa.s, from the footprint
    eta_model_f1: float     # Pa.s, the closure's FULLY FRICTIONAL viscosity
    phi_j_required: float
    phi_m_required: float
    f_at_impact: float

    @property
    def viscosity_gap(self) -> float:
        """How many times more viscous the footprint demands than f = 1 allows."""
        return self.eta_implied / self.eta_model_f1


def implied_viscosity(d: Drop, p: float = SPREAD_EXPONENT_PRIMARY) -> float:
    """Effective viscosity the arrested footprint implies [Pa.s].

    ``eta = rho U0 D0 / beta^p``. Prefactor 1; see the module docstring.
    """
    m = material(d.phi)
    return float(m.rho * d.u0 * d.d0 / d.beta_p**p)


def required_phi_j(d: Drop, eta: float) -> float:
    """Jamming fraction the closure needs to produce ``eta`` at this drop's phi.

    Inverts ``eta = eta_s (1 - phi/phi_J)^-n`` with n = 2, the exponent this
    closure was fitted at:

        phi_J = phi / (1 - (eta/eta_s)^(-1/n))

    Returns NaN where ``eta <= eta_s``, which is outside the model.
    """
    m = material(d.phi)
    r = eta / m.eta_s
    if not np.isfinite(r) or r <= 1.0:
        return float("nan")
    gap = r ** (-0.5)           # n = 2
    if gap >= 1.0:
        return float("nan")
    return float(d.phi / (1.0 - gap))


def required_phi_m(d: Drop, phi_j: float, f: float) -> float:
    """phi_m the closure needs, given the required phi_J and the actual f.

    Inverts ``phi_J = phi_0 - (phi_0 - phi_m) f``:

        phi_m = phi_0 - (phi_0 - phi_J) / f

    With f near 1 this is close to phi_J itself; the correction matters where
    the activation is incomplete, which on the deviatoric gauge it is.
    """
    m = material(d.phi)
    if not np.isfinite(phi_j) or f <= 0:
        return float("nan")
    return float(m.phi0 - (m.phi0 - phi_j) / f)


def invert(d: Drop, p: float = SPREAD_EXPONENT_PRIMARY) -> Inversion:
    m = material(d.phi)
    eta = implied_viscosity(d, p)
    f = frictional_fraction_shear(d)
    pj = required_phi_j(d, eta)
    return Inversion(
        label=d.label, phi=d.phi, salt_mM=d.salt_mM, u0=d.u0, beta_p=d.beta_p,
        eta_implied=eta, eta_model_f1=float(m.eta_frictional),
        phi_j_required=pj, phi_m_required=required_phi_m(d, pj, f),
        f_at_impact=f,
    )


def drift(values, against, labels=("required phi_m", "u0")):
    """Spearman correlation and an OLS slope of one array against another.

    Returned rather than printed so the caller decides what to do with it.
    """
    from scipy.stats import spearmanr
    v = np.asarray(values, float)
    a = np.asarray(against, float)
    ok = np.isfinite(v) & np.isfinite(a)
    v, a = v[ok], a[ok]
    if v.size < 3:
        return None
    rho, pval = spearmanr(a, v)
    slope, intercept = np.polyfit(a, v, 1)
    return {"n": int(v.size), "rho": float(rho), "p": float(pval),
            "slope": float(slope), "intercept": float(intercept),
            "span": float(v.max() - v.min()),
            "x": labels[1], "y": labels[0]}


def prefactor_phi_m(d: Drop, p: float, C: float) -> float:
    """Required phi_m with an explicit prefactor on the spreading law.

    The O(1) coefficient in ``eta = C rho U0 D0 / beta^p`` is not known. This
    exists so the sensitivity to it is measured rather than assumed away: the
    divergence is steep, so a large implied viscosity forces phi_J close to phi
    whatever C is, and the answer should barely move.
    """
    m = material(d.phi)
    eta = C * m.rho * d.u0 * d.d0 / d.beta_p**p
    r = eta / m.eta_s
    if not np.isfinite(r) or r <= 1.0:
        return float("nan")
    gap = r ** -0.5
    if gap >= 1.0:
        return float("nan")
    phi_j = d.phi / (1.0 - gap)
    f = frictional_fraction_shear(d)
    return float(m.phi0 - (m.phi0 - phi_j) / f)
