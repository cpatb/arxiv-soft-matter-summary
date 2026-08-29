"""
conventions.py — "what IS the shear rate of an impact?", and why the answer
moves by a factor of ten.

THE PROBLEM
-----------
Every attempt to connect drop impact to rheology needs one number: the shear
rate the impact imposes. Four different answers are in print, three of them
from Cheng's or Driscoll's groups, and they span an order of magnitude:

  (A) gammadot = U0/D0
      Shah, Arora & Driscoll, Commun. Phys. 5, 222 (2022), Methods; also
      Boyer et al. PRF 1, 013901 (2016); also Shah & Driscoll, Soft Matter
      (2024). Argument: at t = 0 the drop's bottom is at rest and its apex
      still moves at U0, so U0/D0 is the maximum rate. Quoted 233-1333 1/s.
      THIS IS THE ABSCISSA OF THIS LAB'S OWN FITTED gammadot_c = 632 1/s.

  (B) gammadot = U0/(2 D0)
      Same paper, used to pick eta for Re_eff: "average of the two extremes".

  (C) gammadot = U0/l with l ~ 100 um, the minimum lamella thickness
      Same paper again, for Re_eff; gives ~1e4 1/s.

  (D) gammadot = U0/L_c with L_c = 0.986 Re^alpha(We) in units of D0
      Mobaseri, Kumar & Cheng, PNAS 122, e2500163122 (2025).
      alpha(We) = -0.4 + 0.06 exp(-0.03 We);  Re = rho U0 D0 / mu(1/L_c)
      Calibrated against the MEASURED viscous dissipation integral in >1200
      Basilisk simulations. It is the only one of the four that was validated
      rather than asserted, and the paper explicitly REFUTES (A): the measured
      L_c is always below 0.45 D0, so (A) understates the rate by >2x.

  (E) gammadot(r,t) from the Philippi field, tau/eta with tau from
      CSG22 Eq. [13] — a field, not a number. See stress.py.

WHY THIS MATTERS HERE
---------------------
This lab's headline number, the impact-solidification onset
gammadot_c = 632 / 479 / 195 1/s across 0/10/100 mM LiCl, is on convention
(A). Any comparison with a rheometer's critical shear rate is a comparison of
(A) against a real shear rate, and the conversion factor D0/L_c is ~5-10 for
this system. That factor is larger than the salt effect being measured.

Nothing here says the lab's number is wrong. It says the number is
convention-bound, and that the convention has to be stated every time it is
compared with anything.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

# Mobaseri, Kumar & Cheng, PNAS 122, e2500163122 (2025)
LC_C = 0.986          # +/- 0.002
LC_A0 = -0.4
LC_A1 = 0.06
LC_A2 = -0.03


def lc_exponent(We: float) -> float:
    """alpha(We) = -0.4 + 0.06 exp(-0.03 We).

    At We > ~150 the exponential is dead and alpha = -0.4 exactly, so
    L_c ~ Re^(-2/5) — the same exponent as the classical viscous
    maximum-spreading law, which is not a coincidence: both come from the
    boundary layer meeting the descending free surface.
    """
    return LC_A0 + LC_A1 * np.exp(LC_A2 * We)


def lc_over_D0(Re: float, We: float) -> float:
    """L_c / D0 = 0.986 Re^alpha(We)."""
    return LC_C * Re ** lc_exponent(We)


def solve_Lc(rho: float, U0: float, D0: float, gamma: float,
             eta_of_gammadot, eta_guess: float) -> dict:
    """Solve the PNAS-2025 definition self-consistently.

    L_c depends on Re, Re depends on mu evaluated AT the shear rate 1/L_c, and
    1/L_c depends on L_c. For a Newtonian drop it closes in one step; for a
    shear-thickening one it must be iterated, and it is iterated here by
    bisection on log(L_c) rather than by fixed-point (which oscillates once
    the viscosity is rising steeply with rate).

    Parameters
    ----------
    eta_of_gammadot : callable, shear rate [1/s] -> viscosity [Pa.s]
    eta_guess : starting viscosity, used only to bracket the root

    Returns dict with L_c [m], gammadot_c [1/s], Re, We, and the viscosity
    actually used.
    """
    We = rho * U0**2 * D0 / gamma

    def residual(log_Lc_over_D0):
        LcD = np.exp(log_Lc_over_D0)
        gd = U0 / (LcD * D0)
        mu = eta_of_gammadot(gd)
        Re = rho * U0 * D0 / mu
        return np.log(lc_over_D0(Re, We)) - log_Lc_over_D0

    lo, hi = np.log(1e-4), np.log(1.0)
    if residual(lo) * residual(hi) > 0:
        raise ValueError("L_c root not bracketed; check eta_of_gammadot")
    log_LcD = brentq(residual, lo, hi, xtol=1e-12)
    LcD = float(np.exp(log_LcD))
    gd = U0 / (LcD * D0)
    mu = eta_of_gammadot(gd)
    return {
        "L_c_over_D0": LcD,
        "L_c_m": LcD * D0,
        "gammadot_c_per_s": gd,
        "Re": rho * U0 * D0 / mu,
        "We": We,
        "eta_used_Pa_s": mu,
        "amplification_vs_U0_over_D0": 1.0 / LcD,
        "eta_guess_Pa_s": eta_guess,
    }


def convention_table(rho: float, U0: float, D0: float, gamma: float,
                     eta_of_gammadot, eta0: float) -> dict:
    """All four scalar conventions side by side, in 1/s."""
    sol = solve_Lc(rho, U0, D0, gamma, eta_of_gammadot, eta0)
    return {
        "(A) U0/D0  [Shah 2022; this lab's fitted axis]": U0 / D0,
        "(B) U0/2D0 [Shah 2022, for Re_eff]": U0 / (2 * D0),
        "(C) U0/l, l=100 um [Shah 2022, for Re_eff]": U0 / 100e-6,
        "(D) U0/L_c [Mobaseri 2025, validated]": sol["gammadot_c_per_s"],
        "_L_c_over_D0": sol["L_c_over_D0"],
        "_amplification_D_over_Lc": sol["amplification_vs_U0_over_D0"],
        "_Re": sol["Re"],
        "_We": sol["We"],
    }
