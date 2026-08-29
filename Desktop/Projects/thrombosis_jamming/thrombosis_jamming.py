"""
Thrombosis as a stress-activated jamming transition — a Tier-0 closure toy.

Third fluid in the "one state-variable jamming closure, many fluids" program
(silica drop impact -> blood rheology -> arterial thrombosis). Implements the
minimal Wyart-Cates-style closure for shear-induced platelet aggregation (SIPA)
proposed as a wet-lab-free first project, and produces its distinctive,
FALSIFIABLE prediction: a CRITICAL DIVERGENCE of occlusion time at a threshold
wall shear rate — a genuine critical point that the incumbent deposition-kinetics
thrombosis models (which give a smooth monotonic speed-up) do not predict.

Physics (all ILLUSTRATIVE parameters; the point is the qualitative signature)
---------------------------------------------------------------------------
- Activation is by VWF unfolding, set by a CRITICAL SHEAR RATE (not stress):
      f(gdot) = 1 - exp(-(gdot/gdot_c)^beta),   gdot_c ~ 5000 /s  (Schneider 2007)
  f = fraction of VWF-bridged (adhesive) platelet contacts.
- Migrating jamming point (Wyart-Cates skeleton):
      phi_J(f) = phi_m*f + phi0*(1-f),   phi_m < phi0
  bonded contacts let the network rigidify at looser packing (phi_m).
- The throat platelet packing densifies by deposition toward phi_max:
      phi(t) = phi0_base + (phi_max - phi0_base)*(1 - exp(-t/tau)),  tau ~ 1/f
- Occlusion when phi(t) reaches phi_J(f) (effective viscosity diverges):
      t_occ = -tau * ln(1 - (phi_J(f) - phi0_base)/(phi_max - phi0_base))
  If phi_J(f) > phi_max (low shear, little activation) -> NEVER occludes.
  As gdot -> gdot_crit+ (where phi_J = phi_max), t_occ -> infinity: CRITICAL DIVERGENCE.

Falsifiable test: replot published occlusion-time-vs-shear data against
(gdot - gdot_crit); a jamming mechanism collapses to a power-law divergence,
a deposition mechanism does not.
"""
from __future__ import annotations

import numpy as np

# --- illustrative constants ---
GDOT_C = 5000.0     # 1/s, VWF unfolding critical shear rate (Schneider et al. 2007)
BETA = 2.0
PHI0 = 0.64         # unbonded platelet-suspension jamming fraction
PHI_M = 0.40        # bonded/adhesive jamming fraction (looser, adhesion-stabilised)
PHI0_BASE = 0.30    # baseline throat cellular packing
PHI_MAX = 0.62      # max packing reachable by deposition
TAU0 = 1.0          # s, deposition timescale prefactor (tau = TAU0 / f)


def bonded_fraction(gdot, gdot_c=GDOT_C, beta=BETA):
    """f(gdot): fraction of VWF-bridged adhesive contacts, in [0,1)."""
    gdot = np.asarray(gdot, dtype=float)
    return 1.0 - np.exp(-((gdot / gdot_c) ** beta))


def jamming_fraction(f, phi0=PHI0, phi_m=PHI_M):
    return phi_m * f + phi0 * (1.0 - f)


def occlusion_time(gdot, *, phi0_base=PHI0_BASE, phi_max=PHI_MAX, tau0=TAU0):
    """Occlusion time vs wall shear rate. np.inf where no occlusion occurs.

    Returns (t_occ, gdot_crit). Vectorised over gdot.
    """
    gdot = np.asarray(gdot, dtype=float)
    f = bonded_fraction(gdot)
    phi_J = jamming_fraction(f)
    tau = tau0 / np.clip(f, 1e-6, None)             # faster deposition at higher shear
    # fraction of the deposition "budget" needed to reach the jamming point
    frac = (phi_J - phi0_base) / (phi_max - phi0_base)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(frac < 1.0, -tau * np.log(np.clip(1.0 - frac, 1e-12, None)), np.inf)
    t = np.where(frac < 0.0, 0.0, t)                # already jammed at baseline
    # critical shear rate where phi_J(f) == phi_max  (occlusion onset)
    #   phi_J = phi_max -> f_crit, then invert f(gdot)
    f_crit = (PHI0 - phi_max) / (PHI0 - PHI_M)
    f_crit = min(max(f_crit, 1e-6), 1 - 1e-6)
    gdot_crit = GDOT_C * (-np.log(1.0 - f_crit)) ** (1.0 / BETA)
    return t, gdot_crit
