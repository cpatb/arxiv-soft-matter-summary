"""
Blood as a structured suspension — a bridge from the drop-impact closure work
to biological fluids.

Thesis of the bridge
--------------------
The Tier-0 shear-thickening closure (silica, `tier0_wyart_cates_closure`) uses a
first-order-kinetics STATE VARIABLE for the fraction of frictional contacts f,
which RISES with stress and drives thickening / jamming. Blood is the mirror
image in the SAME mathematical framework: a state variable lambda for the
fraction of aggregated red cells (rouleaux) that FALLS with shear rate and drives
shear THINNING. One closure skeleton, opposite sign of the structure->viscosity
coupling. This file demonstrates that transferability concretely.

Method
------
Structural (thixotropic) rheology model of the Moore / Coussot lineage,
standard for blood:
  - Aggregation state lambda in [0, 1] (1 = fully aggregated rouleaux network).
  - Steady state:  lambda_ss(gamma_dot) = 1 / (1 + (gamma_dot / gamma_c))
    (build by Brownian/bridging, break by shear; monotone decreasing).
  - Structural viscosity:  eta(lambda) = eta_inf * (1 + (eta0/eta_inf - 1) * lambda)
    so eta -> eta0 (aggregated, low shear) and eta -> eta_inf (dispersed, high shear).
  - Casson-style yield from the rouleaux network at low shear (added as an
    illustrative offset; Casson is the textbook blood constitutive law).
Transient thixotropy:
  d(lambda)/dt = k_build*(1 - lambda) - k_break*gamma_dot*lambda
which gives the hysteresis loop blood shows under an up/down shear ramp -- the
biological analogue of the silica retention hysteresis in the drop-impact work.

Parameters are physiological-ballpark (Hct ~ 45%), illustrative not fitted:
  eta_inf ~ 3.5 mPa.s (high-shear plasma+cells), eta0 ~ 200 mPa.s (aggregated),
  yield ~ few mPa. To be calibrated to real viscometry next.

Assumptions/caveats: homogeneous shear, no Fahraeus-Lindqvist (no confinement),
no margination, single relaxation time. This is a 0-D constitutive toy, the
entry rung for a biological-fluids track.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

# --- physiological-ballpark constants (Hct ~ 45%) ---
ETA_INF = 3.5e-3     # Pa.s, high-shear (dispersed) blood viscosity
ETA0 = 200.0e-3      # Pa.s, low-shear (aggregated) blood viscosity
GAMMA_C = 5.0        # 1/s, shear rate where rouleaux break up
SIGMA_YIELD = 3.0e-3  # Pa, Casson-like yield from the rouleaux network


def lambda_ss(gamma_dot, gamma_c=GAMMA_C):
    """Steady aggregation fraction; monotone decreasing in shear rate, in [0,1]."""
    gamma_dot = np.asarray(gamma_dot, dtype=float)
    return 1.0 / (1.0 + gamma_dot / gamma_c)


def structural_viscosity(lam, eta_inf=ETA_INF, eta0=ETA0):
    """eta(lambda): eta0 when aggregated (lambda=1), eta_inf when dispersed (0)."""
    lam = np.asarray(lam, dtype=float)
    return eta_inf * (1.0 + (eta0 / eta_inf - 1.0) * lam)


def steady_flow_curve(gamma_dot, *, sigma_yield=SIGMA_YIELD, **kw):
    """Steady blood flow curve. Returns dict of arrays.

    eta_app is the apparent viscosity sigma/gamma_dot with a Casson-like yield.
    """
    gamma_dot = np.asarray(gamma_dot, dtype=float)
    lam = lambda_ss(gamma_dot, kw.get("gamma_c", GAMMA_C))
    eta = structural_viscosity(lam, kw.get("eta_inf", ETA_INF), kw.get("eta0", ETA0))
    sigma = sigma_yield + eta * gamma_dot          # structural + yield
    eta_app = sigma / gamma_dot
    return {"gamma_dot": gamma_dot, "lambda": lam, "eta_struct": eta,
            "sigma": sigma, "eta_app": eta_app}


def thixotropic_loop(gamma_max=100.0, t_ramp=5.0, n=600,
                     k_build=1.0, k_break=0.05, gamma_c=GAMMA_C):
    """Up-then-down linear shear-rate ramp; returns the (gamma_dot, sigma) loop.

    A finite loop area = thixotropic hysteresis (structure cannot rebuild
    instantly on the down-ramp). This is the blood analogue of the silica
    retention hysteresis. Returns dict with up/down branches.
    """
    t_up = np.linspace(0.0, t_ramp, n)
    t_dn = np.linspace(t_ramp, 2 * t_ramp, n + 1)[1:]  # drop duplicate apex point

    def gdot(t):
        if t <= t_ramp:
            return gamma_max * (t / t_ramp)
        return gamma_max * (1.0 - (t - t_ramp) / t_ramp)

    def rhs(t, y):
        g = gdot(t)
        return [k_build * (1.0 - y[0]) - k_break * g * y[0]]

    # start from the steady structure at gamma_dot = 0 (fully aggregated)
    sol = solve_ivp(rhs, (0.0, 2 * t_ramp), [1.0],
                    t_eval=np.concatenate([t_up, t_dn]),
                    method="LSODA", rtol=1e-8, atol=1e-10)
    if not sol.success:
        raise RuntimeError(sol.message)
    lam = sol.y[0]
    g = np.array([gdot(t) for t in sol.t])
    eta = structural_viscosity(lam, ETA_INF, ETA0)
    sigma = SIGMA_YIELD + eta * g
    half = len(t_up)
    return {
        "g_up": g[:half], "sigma_up": sigma[:half],
        "g_dn": g[half:], "sigma_dn": sigma[half:],
        "lambda": lam, "t": sol.t,
    }
