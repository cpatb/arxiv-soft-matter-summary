"""poro_impact.py — 1D consolidation with a Reynolds dilatancy source, during impact.

WHY THIS EXISTS
---------------
`derivations/dilatancy_sign_mechanism.md` derives the pore-pressure ratio

    Pi = alpha mu_w |phi - phi_c| ell^2 / (rho kappa U0 D0)                     Eq. (6)

by balancing what a dilating skeleton demands against what Darcy can supply. That is a
QUASI-STEADY, SINGLE-LENGTH balance (assumptions A2 and A5), and it is exactly the
substitution `[M26]` flag with "~". Two things it cannot say:

  1. whether the balance is reached within the impact time at all, and
  2. the O(1) geometric prefactor, which the "~" leaves undetermined.

This module solves the transient problem instead, so both become measurements.

THE MODEL, DERIVED NOT ASSERTED
-------------------------------
Saturated pack, incompressible pore fluid, incompressible grains. Mixture
incompressibility plus Darcy gives, with `eps_v` the skeleton's volumetric strain rate
(POSITIVE IN DILATION),

    div(v_s) = eps_v = div[(kappa/mu_w) grad p]                                  (i)

The skeleton's volume change has an elastic part and a Reynolds part, with `K'` the
DRAINED bulk modulus and sigma' = sigma - p the effective mean stress (compression
positive):

    eps_v = -(1/K') d(sigma')/dt + alpha gammadot (phi - phi_c)                  (ii)

Eliminating sigma' = sigma - p and substituting into (i):

    dp/dt = c_v d2p/dx2 - K' alpha gammadot (phi - phi_c) + d(sigma)/dt          (iii)
    c_v   = kappa K' / mu_w          [m^2/s]  the consolidation coefficient

which is Terzaghi consolidation with a dilatancy source and an undrained-loading term.
SIGN: phi > phi_c makes the source term NEGATIVE, so `p` falls -- suction -- which is
what the derivation claims. phi < phi_c makes it positive: overpressure.

The packing fraction evolves as the skeleton dilates, which SHUTS THE SOURCE OFF:

    dphi/dt = -phi * alpha gammadot (phi - phi_c)                                (iv)

so `phi - phi_c` decays on a strain of order 1/alpha. This is the model's own account of
why the author's bench observation is *"inertial-dominated for some time"* rather than
permanently: the effect has a strain clock built into it.

NONDIMENSIONAL FORM (what the code actually integrates)
-------------------------------------------------------
    x~ = x/L,  t~ = t U0/D0,  p~ = p/(rho U0^2),  gammadot = (U0/D0) g(t~)

    dp~/dt~ = Lambda d2p~/dx~2 - S (phi - phi_c) g(t~) + dsigma~/dt~
    Lambda  = c_v D0 / (U0 L^2)         drainage number  (= c_v t_imp / L^2)
    S       = K' alpha / (rho U0^2)     dilatancy-stiffness number

and the QUASI-STEADY limit of this system reproduces Eq. (6) exactly, with a geometric
prefactor the closed form cannot supply:

    p~_peak = S (phi - phi_c) / (2 Lambda) = Pi/2       (impermeable base, drained top)

`test_poro.py::test_quasistatic_limit_recovers_Pi_over_two` pins that, and
`test_poro.py::test_terzaghi_series` pins the source-free solver against the published
analytic consolidation series -- an anchor that does not pass through this code.

BOUNDARY CONDITIONS. x~ = 0 is the substrate: impermeable, dp~/dx~ = 0. x~ = 1 is the
free surface: drained, p~ = 0 (gauge). `L` is the DRAINAGE LENGTH, and section 6 of the
derivation says nobody has measured it -- so it is an explicit input here and every
result is reported against it.

CAVITATION. The pore fluid cannot sustain more than about one atmosphere of suction.
With `p_cav_hat` set, `p~` is prevented from falling below -p_cav_hat and the binding is
REPORTED, because the derivation's own estimate exceeds that ceiling by ~459x.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

# --- named constants, with origin ---------------------------------------------------
P_ATM = 1.013e5          # Pa; the suction a pore fluid sustains before cavitating
ALPHA_DIL = 1.0          # tan(psi) = alpha (phi - phi_c); alpha ~ 1 per [M26]
N_NODES = 201            # method-of-lines grid; convergence is checked in test_poro.py


@dataclass
class PoroResult:
    """One transient solve. JSON-serialisable."""

    # inputs (all dimensionless)
    lam: float                    # Lambda, drainage number
    s_dil: float                  # S, dilatancy-stiffness number
    dphi0: float                  # phi - phi_c at t = 0  (SIGNED: the mechanism)
    phi0: float
    t_end: float
    alpha: float = ALPHA_DIL
    p_cav_hat: float = float("inf")
    n_nodes: int = N_NODES
    loading: bool = False         # include the undrained-loading term dsigma~/dt~
    # outputs
    p_peak_hat: float = float("nan")     # signed extremum of p~ over space and time
    t_at_peak: float = float("nan")
    p_base_final: float = float("nan")
    dphi_final: float = float("nan")
    strain_to_half_dphi: float = float("nan")
    cavitated: bool = False
    frac_time_cavitated: float = 0.0
    steady_ratio: float = float("nan")   # p_peak / the quasi-steady prediction Pi/2
    ok: bool = False
    p_final_field: list[float] = field(default_factory=list)
    t: list[float] = field(default_factory=list)
    p_base: list[float] = field(default_factory=list)
    dphi: list[float] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_json(self) -> str:
        return json.dumps(asdict(self), sort_keys=True)


def _shear_history(t: float) -> float:
    """g(t~), the shear-rate history, normalised so g(0) = 1.

    Constant by default: the impact delivers a roughly steady strain rate U0/D0 over the
    contact time. A shaped history is a knob this module deliberately does not turn --
    the mechanism's claims are about the SIGN, which no positive g(t) can change.
    """
    return 1.0


def solve(lam: float, s_dil: float, dphi0: float, phi0: float = 0.44,
          t_end: float = 2.0, alpha: float = ALPHA_DIL,
          p_cav_hat: float = float("inf"), n_nodes: int = N_NODES,
          loading: bool = False, dsigma_dt: float = 1.0,
          p_init: float = 0.0, freeze_phi: bool = False,
          rtol: float = 1e-8, atol: float = 1e-12) -> PoroResult:
    """Integrate (iii)+(iv) on 0 <= x~ <= 1, impermeable at 0, drained at 1.

    `p_init` seeds a uniform initial pore pressure -- used only by the Terzaghi anchor,
    where the source is switched off and the decay of that initial field is compared
    against the published series. `freeze_phi` holds (iv) still, which is what the
    quasi-steady comparison needs: with phi free the source shuts itself off before a
    steady state can be reached, so the two tests would be testing each other.
    """
    res = PoroResult(lam=lam, s_dil=s_dil, dphi0=dphi0, phi0=phi0, t_end=t_end,
                     alpha=alpha, p_cav_hat=p_cav_hat, n_nodes=n_nodes,
                     loading=loading)
    res.notes.append(f"p_init={p_init:g} freeze_phi={freeze_phi}")
    n = n_nodes
    dx = 1.0 / (n - 1)
    phi_c = phi0 - dphi0                        # held fixed; phi is what evolves

    def rhs(_t, y):
        p = y[:n]
        phi = y[n:]
        g = _shear_history(_t)
        # diffusion, second order, with dp/dx = 0 at x=0 (ghost node mirror) and
        # p = 0 enforced at x = 1 by holding that node's value fixed.
        d2 = np.empty_like(p)
        d2[0] = 2.0 * (p[1] - p[0]) / dx**2
        d2[1:-1] = (p[2:] - 2 * p[1:-1] + p[:-2]) / dx**2
        d2[-1] = 0.0
        src = s_dil * alpha * (phi - phi_c) * g
        dp = lam * d2 - src
        if loading:
            dp = dp + dsigma_dt
        dp[-1] = 0.0                            # Dirichlet p~(1) = 0
        # cavitation: the fluid cannot be pulled below -p_cav_hat
        if np.isfinite(p_cav_hat):
            at_floor = (p <= -p_cav_hat) & (dp < 0)
            dp = np.where(at_floor, 0.0, dp)
        dphi_dt = (np.zeros_like(phi) if freeze_phi
                   else -phi * alpha * (phi - phi_c) * g)
        return np.concatenate([dp, dphi_dt])

    y0 = np.concatenate([np.full(n, p_init), np.full(n, phi0)])
    y0[n - 1] = 0.0
    # Sparsity pattern of the Jacobian. dp_i touches p_{i-1}, p_i, p_{i+1} and phi_i;
    # dphi_i touches phi_i only. Handing this to BDF turns an O(n^3) dense factorisation
    # into a banded one -- with n = 401 the dense version was minutes per solve.
    from scipy.sparse import lil_matrix

    jac = lil_matrix((2 * n, 2 * n), dtype=int)
    for i in range(n):
        for j in (i - 1, i, i + 1):
            if 0 <= j < n:
                jac[i, j] = 1
        jac[i, n + i] = 1
        jac[n + i, n + i] = 1
    sol = solve_ivp(rhs, (0.0, t_end), y0, method="BDF", rtol=rtol, atol=atol,
                    dense_output=False, max_step=t_end / 50.0,
                    jac_sparsity=jac.tocsr())
    if not sol.success:
        res.notes.append(f"SOLVER FAILED: {sol.message}")
        return res

    P = sol.y[:n, :]
    PHI = sol.y[n:, :]
    res.ok = True
    res.t = [float(v) for v in sol.t]
    res.p_base = [float(v) for v in P[0, :]]
    res.dphi = [float(v) for v in (PHI[0, :] - phi_c)]

    k = int(np.argmax(np.abs(P)))
    res.p_peak_hat = float(P.ravel()[k])
    res.t_at_peak = float(sol.t[k % P.shape[1]])
    res.p_base_final = float(P[0, -1])
    res.dphi_final = float(PHI[0, -1] - phi_c)
    res.p_final_field = [float(v) for v in P[:, -1]]

    # strain at which |phi - phi_c| has halved. gammadot = 1 in these units, so t~ IS
    # the strain; NOT read off an assumed monotone decay -- (iv) is monotone in sign by
    # construction (dphi and its derivative have opposite signs), and the crossing is
    # taken as the FIRST one with the rule named here.
    d = np.abs(np.asarray(res.dphi))
    if d[0] > 0:
        below = np.flatnonzero(d <= 0.5 * d[0])
        res.strain_to_half_dphi = (float(sol.t[below[0]]) if below.size
                                   else float("nan"))

    if np.isfinite(p_cav_hat):
        hit = np.any(P <= -p_cav_hat * (1 - 1e-9), axis=0)
        res.cavitated = bool(hit.any())
        res.frac_time_cavitated = float(hit.mean())

    # the quasi-steady prediction this transient is being compared against
    pi_over_two = -s_dil * alpha * dphi0 / (2.0 * lam)
    res.steady_ratio = (float(res.p_peak_hat / pi_over_two)
                        if pi_over_two != 0 else float("nan"))
    return res


def quasistatic_peak(lam: float, s_dil: float, dphi0: float,
                     alpha: float = ALPHA_DIL) -> float:
    """The closed form this solver is being tested against: p~ at the impermeable base.

    Steady state of  Lambda p'' = S alpha dphi  with p'(0) = 0, p(1) = 0 gives
    p(x) = (S alpha dphi / (2 Lambda)) (x^2 - 1), hence p(0) = -S alpha dphi/(2 Lambda).
    The factor 1/2 is the geometric prefactor that Eq. (6)'s "~" does not supply.
    """
    return -s_dil * alpha * dphi0 / (2.0 * lam)


def terzaghi_series(x: np.ndarray, T: float, n_terms: int = 200) -> np.ndarray:
    """Analytic 1D consolidation, uniform initial p, SEALED at x=0, DRAINED at x=1.

        p/p0 = sum_m (2/M) sin(M (1 - x)) exp(-M^2 T),   M = (2m+1) pi / 2

    Standard result (Terzaghi 1943), written for the orientation this solver uses. The
    textbook form is usually printed with the coordinate measured from the DRAINED face,
    where the eigenfunctions are sin(M z); measured from the SEALED face they are
    cos(M x), and sin(M(1-x)) = (-1)^m cos(M x) is the same series with the sign
    absorbed.

    ⚠ The first version of this function used sin(M x) with x from the sealed face, i.e.
    the mirror image of the solver's geometry, and disagreed with it by 0.61 in a field
    of order 1. The anchor caught it on its first run, which is what an anchor outside
    the code under test is for. `test_terzaghi_orientation_is_not_mirrored` pins it.

    Used ONLY as an anchor; `solve` never calls it.
    """
    out = np.zeros_like(np.asarray(x, float))
    for m in range(n_terms):
        k = (2 * m + 1) * math.pi / 2.0
        out += (2.0 / k) * np.sin(k * (1.0 - np.asarray(x, float))) * math.exp(-(k**2) * T)
    return out


def physical_numbers(kappa: float, mu_w: float, k_drained: float, rho: float,
                     u0: float, d0: float, ell: float) -> dict:
    """Lambda and S from physical inputs, so the map has an address in SI."""
    c_v = kappa * k_drained / mu_w
    return {
        "c_v_m2_s": c_v,
        "t_imp_s": d0 / u0,
        "Lambda": c_v * d0 / (u0 * ell**2),
        "S": k_drained * ALPHA_DIL / (rho * u0**2),
        "p_cav_hat": P_ATM / (rho * u0**2),
    }
