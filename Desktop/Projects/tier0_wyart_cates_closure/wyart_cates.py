"""
Tier-0 reduced-order closure for shear-thickening / shear-jamming suspensions.

Method
------
Wyart & Cates frictional-contact model of shear thickening
(M. Wyart & M. E. Cates, Phys. Rev. Lett. 112, 098302, 2014), in the
steady + first-order-relaxation ("state variable") form used across the
dense-suspension literature (Mari, Seto, Morris; Baumgarten & Kamrin).

Physical picture
----------------
- A stress-activated fraction of frictional contacts f(sigma) rises from 0
  (lubricated) to 1 (frictional) as the shear stress crosses an onset stress
  sigma_star.
- The jamming volume fraction interpolates between a frictionless value phi0
  (~ random close packing) at f=0 and a lower frictional value phi_m at f=1:
      phi_J(f) = (1 - f) * phi0 + f * phi_m
- The relative viscosity diverges as phi -> phi_J:
      eta_r(phi, f) = (1 - phi / phi_J(f)) ** (-2)
- Regimes at fixed packing phi:
      phi  <  phi_m        -> continuous thickening (CST); always flows
      phi_m < phi < phi0   -> discontinuous thickening / shear jamming (DST/SJ);
                              at high stress phi_J(f->1)=phi_m < phi so eta_r -> inf
                              and the flow curve back-bends.

Reversibility (Tier-0 state variable)
-------------------------------------
The frictional fraction is given memory via first-order kinetics toward a
possibly hysteretic target:
      df/dt = (f_target(sigma, f) - f) / tau
      f_target = max( f_ss(sigma), f_residual * S(f) )
where S(f) is a smooth on-switch (~0 below f_on, ~1 above). f_residual = 0 is a
fully reversible material: on stress removal f -> 0 (recovery, corn starch). A
positive f_residual encodes an irreversible frictional / DLVO-collapsed network
that does not fully release once formed -> the microstructure is retained even
after the stress is gone (a "stuck pancake", silica). The residual only engages
after the network has formed (S(f)~1), so it never spuriously activates at rest.
tau additionally sets how fast each approaches its target.

All parameters here are ILLUSTRATIVE, chosen to reproduce the qualitative
contrast as a cheap sanity gate. They are meant to be *calibrated* against
measured flow curves next (see README / Simulations MOC).

Assumptions / caveats
---------------------
- Homogeneous simple shear; no spatial structure, no shear-banding.
- Additive Bingham-style yield offset for adhesive systems (corn starch);
  the yield term is phenomenological, not derived.
- Retention is modelled here purely as timescale separation (large tau).
  True contact hysteresis (irreversible frictional network) is a Tier-0.5
  extension noted in the README.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp


# --- default microstructure constants (spheres; illustrative) ---
PHI0_DEFAULT = 0.65    # frictionless (lubricated) jamming fraction ~ RCP
PHI_M_DEFAULT = 0.56   # frictional jamming fraction (random loose ~ mu-dependent)


def frictional_fraction_ss(sigma, sigma_star, beta=1.0):
    """Steady-state fraction of frictional contacts, f_ss in [0, 1).

    f_ss(sigma) = exp(-(sigma_star / sigma) ** beta),  f_ss(0) := 0.
    Vectorised; safe at sigma = 0.
    """
    sigma = np.asarray(sigma, dtype=float)
    out = np.zeros_like(sigma)
    pos = sigma > 0
    out[pos] = np.exp(-((sigma_star / sigma[pos]) ** beta))
    return out


def jamming_fraction(f, phi0=PHI0_DEFAULT, phi_m=PHI_M_DEFAULT):
    """phi_J(f) interpolating frictionless -> frictional jamming."""
    f = np.asarray(f, dtype=float)
    return (1.0 - f) * phi0 + f * phi_m


def relative_viscosity(phi, f, phi0=PHI0_DEFAULT, phi_m=PHI_M_DEFAULT):
    """eta_r(phi, f) = (1 - phi/phi_J)^(-2); +inf once phi >= phi_J (jammed)."""
    phi_J = jamming_fraction(f, phi0, phi_m)
    phi = np.asarray(phi, dtype=float)
    gap = 1.0 - phi / phi_J
    with np.errstate(divide="ignore", invalid="ignore"):
        eta_r = np.where(gap > 0, gap ** (-2.0), np.inf)
    return eta_r


def steady_flow_curve(
    phi,
    sigma_star,
    *,
    phi0=PHI0_DEFAULT,
    phi_m=PHI_M_DEFAULT,
    eta_s=1.0e-3,
    sigma_yield=0.0,
    beta=1.0,
    n_points=400,
    sigma_hydro_range=(1.0e-2, 1.0e4),
):
    """Steady flow curve, parametrised by the hydrodynamic (WC) stress sigma_h.

    Returns dict with arrays: sigma_h, f, phi_J, eta_r, gamma_dot,
    sigma_total (= sigma_yield + sigma_h), eta_app (= sigma_total/gamma_dot).

    Parametrising by sigma_h avoids solving the implicit sigma(gamma_dot);
    the flow curve back-bends (gamma_dot non-monotonic) exactly in the DST/SJ
    regime, which is the signature we test for.
    """
    sigma_h = np.logspace(
        np.log10(sigma_hydro_range[0]), np.log10(sigma_hydro_range[1]), n_points
    )
    f = frictional_fraction_ss(sigma_h, sigma_star, beta)
    phi_J = jamming_fraction(f, phi0, phi_m)
    eta_r = relative_viscosity(phi, f, phi0, phi_m)
    gamma_dot = sigma_h / (eta_r * eta_s)          # -> 0 as eta_r -> inf (jam)
    sigma_total = sigma_yield + sigma_h
    with np.errstate(divide="ignore", invalid="ignore"):
        eta_app = np.where(gamma_dot > 0, sigma_total / gamma_dot, np.inf)
    return {
        "sigma_h": sigma_h,
        "f": f,
        "phi_J": phi_J,
        "eta_r": eta_r,
        "gamma_dot": gamma_dot,
        "sigma_total": sigma_total,
        "eta_app": eta_app,
    }


def is_discontinuous(flow_curve):
    """DST/SJ signature: the flow curve back-bends, i.e. gamma_dot is
    non-monotonic in sigma_h (there is a range where d(gamma_dot)/d(sigma_h) < 0).
    """
    gd = flow_curve["gamma_dot"]
    finite = np.isfinite(gd)
    gd = gd[finite]
    if gd.size < 3:
        return False
    return bool(np.any(np.diff(gd) < 0))


def relax_state(
    t_eval,
    sigma_of_t,
    tau,
    sigma_star,
    *,
    f0=0.0,
    beta=1.0,
    f_residual=0.0,
    f_on=0.5,
    width=0.05,
):
    """Integrate df/dt = (f_target(sigma(t), f) - f)/tau over t_eval.

    f_target = max(f_ss(sigma), f_residual * S(f)),  S = logistic on-switch.
    f_residual = 0 -> fully reversible (recovery); f_residual > 0 -> retention
    of the frictional network once formed (irreversibility / "stuck pancake").

    sigma_of_t : callable t -> shear stress (Pa).
    Returns f(t) sampled at t_eval (shape == t_eval).
    """
    # Two states: f (frictional fraction) and a in [0,1], an IRREVERSIBLE
    # activation latch that ratchets up once the network forms (S(f)~1) and
    # never decays. The retained floor is f_residual * a, so it persists after
    # the stress is removed. f_residual=0 -> no latch effect (reversible).
    k_on = 1.0 / tau  # latch engages on the same fast scale as f rises

    def rhs(t, y):
        f, a = y
        f_ss = float(frictional_fraction_ss(sigma_of_t(t), sigma_star, beta))
        switch = 1.0 / (1.0 + np.exp(-(f - f_on) / width))  # ~0 until network forms
        f_target = max(f_ss, f_residual * a)
        df = (f_target - f) / tau
        da = k_on * switch * (1.0 - a)   # monotonic; no deactivation term
        return [df, da]

    sol = solve_ivp(
        rhs,
        (t_eval[0], t_eval[-1]),
        [f0, 0.0],
        t_eval=t_eval,
        method="LSODA",
        rtol=1e-8,
        atol=1e-10,
        max_step=(t_eval[-1] - t_eval[0]) / 200.0,
    )
    if not sol.success:
        raise RuntimeError(f"state ODE failed: {sol.message}")
    return sol.y[0]
