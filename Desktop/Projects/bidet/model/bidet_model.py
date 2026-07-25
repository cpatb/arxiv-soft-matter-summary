"""
bidet_model.py — validated reduced-order wall-shear + nozzle-optimization model
for the "Optimal Bidet" paper. This is the Tier-A numerical solver the manuscript
promised, plus the turbulent-closure recompute the review demanded.

Core physics
------------
An axisymmetric jet (nozzle diameter d, exit speed V0) impinges on a surface and
spreads as a wall jet whose wall shear stress tau_w(r) must exceed the fecal
detachment shear tau_c out to a contamination radius r_max.

Two closures for tau_w(r), written in a common power-law form
    tau_w(r) = C_tau * rho * V0^2 * Re_d^p * (d/r)^n
with Re_d = V0 d / nu:
  * LAMINAR Watson Region II (the paper's tauII):  p = 1,   n = 5,  C_tau = c2/64
    (from tau_w = c2 rho V0^2 Re_a (a/r)^5, a=d/2, Re_a=Re_d/2 -> C_tau=c2/2/2^5=c2/64)
  * TURBULENT (Phares 2000 form): (tau_w/rho V0^2) Re_d^{1/2} = C_P (r/d)^{-nt}
    -> p = -1/2, n = nt (~1.05 developed turbulent), C_tau = C_P.
    C_P is CALIBRATED to one commercial operating point (the only datum available),
    which is the honest resolution of the laminar model over-predicting water use.

The laminar solver also integrates the film momentum ODE through the transition r0
so the absolute prefactor does not rely on the (order-100) asymptotic matching fudge.

All results are checked against analytic limits in test_bidet_model.py.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

# ---- constants ----
RHO = 1000.0          # water density [kg/m^3]
NU = 8.9e-7           # kinematic viscosity at 25 C [m^2/s]
SIGMA = 0.072         # surface tension [N/m]
CV = 0.97             # velocity coefficient (Torricelli)
CD = 0.62             # discharge coefficient (sharp-edged orifice) -- KEPT (fixes B5)
C2 = 27.0 / 50.0      # Watson Region-II prefactor (0.54)
C_DELTA = np.sqrt(280.0 / 39.0)  # cubic-profile bl constant (~2.68)


def V0_from_dP(dP):
    """Torricelli exit velocity from gauge supply pressure dP [Pa]."""
    return CV * np.sqrt(2.0 * dP / RHO)


def flow_rate(d, V0):
    """Volumetric flow through a nozzle of diameter d at exit speed V0, with C_d."""
    return CD * (np.pi / 4.0) * d**2 * V0        # [m^3/s]


# =====================================================================
# Closures: tau_w(r) = C_tau rho V0^2 Re_d^p (d/r)^n
# =====================================================================
class Closure:
    def __init__(self, C_tau, p, n, label):
        self.C_tau, self.p, self.n, self.label = C_tau, p, n, label

    def tau_w(self, r, V0, d):
        Re_d = V0 * d / NU
        return self.C_tau * RHO * V0**2 * Re_d**self.p * (d / r)**self.n


# Laminar Watson Region II, expressed in the common form
LAMINAR = Closure(C_tau=C2 / 64.0, p=1.0, n=5.0, label="laminar Watson II (r^-5)")


def laminar_tauII_direct(r, V0, a):
    """The paper's Eq.(tauII) directly, for cross-checking the common-form C_tau."""
    Re_a = V0 * a / NU
    return C2 * RHO * V0**2 * Re_a * (a / r)**5


# ---- transition radius (Region I -> II) ----
def r0_over_a(V0, a):
    Re_a = V0 * a / NU
    return (Re_a / (4.0 * C_DELTA**2))**(1.0 / 3.0)   # ~0.330 Re_a^{1/3}


# =====================================================================
# Tier-A numerical film solver (laminar): integrate the Region-II film
# momentum ODE outward from the transition and return tau_w(r) without the
# asymptotic matching constant. Validates the closed-form Region-II law.
# =====================================================================
def film_tau_numeric(r_eval, V0, a):
    """
    Depth-averaged laminar film: mass  r h ubar = Q/(2 pi) = a^2 V0 / 2.
    Region-II momentum (parabolic profile): d/dr( (6/5) r h ubar^2 ) = -3 nu r ubar / h.
    Integrate h(r) from r0 (initial h = a^2/(2 r0), the depth where the bl fills the film).
    tau_w = 2 mu U_s / h, U_s = (3/2) ubar, ubar = Q/(2 pi r h).
    """
    Q_over_2pi = a**2 * V0 / 2.0
    r0 = r0_over_a(V0, a) * a
    h0 = a**2 / (2.0 * r0)

    def rhs(r, h):
        h = max(h[0], 1e-9)
        ubar = Q_over_2pi / (r * h)
        # d/dr[(6/5) r h ubar^2] = -3 nu r ubar/h ; expand with r h ubar = const:
        # (6/5) d/dr( ubar * const ) = -3 nu r ubar/h  ->  (6/5) const ubar' = -3 nu r ubar/h
        # ubar = const/(r h) -> ubar' = -const (h + r h')/(r h)^2. Solve for h'.
        const = Q_over_2pi
        # momentum: (6/5) d/dr(r h ubar^2) = -3 nu r ubar / h
        # r h ubar^2 = const*ubar ; d/dr = const*ubar'
        # ubar' = -const (h + r h')/(r^2 h^2)
        # (6/5) const * [ -const (h + r h')/(r^2 h^2) ] = -3 nu r const/(r h * h) = -3 nu const/(h^2)...
        # solve for h':
        # (6/5) const (h + r h')/(r^2 h^2) = 3 nu / h^2   [signs: both sides positive]
        # (6/5) const (h + r h') / r^2 = 3 nu
        # h + r h' = 3 nu r^2 * 5 / (6 const) = (5 nu r^2)/(2 const)
        hp = ((5.0 * NU * r**2) / (2.0 * const) - h) / r
        return [hp]

    sol = solve_ivp(rhs, [r0, max(np.max(r_eval), r0 * 1.0001)], [h0],
                    t_eval=np.atleast_1d(r_eval), rtol=1e-8, atol=1e-12, dense_output=True)
    h = sol.y[0]
    ubar = Q_over_2pi / (np.atleast_1d(r_eval) * h)
    Us = 1.5 * ubar
    mu = RHO * NU
    return 2.0 * mu * Us / h, r0


# =====================================================================
# Cleaning-criterion inversion: solve tau_w(r_max) = tau_c for nozzle d*,
# then Q* = C_d (pi/4) d*^2 V0. General closure; geometry enters via r_max.
# With standoff H: use the impinging (spread-jet) velocity/width.
# =====================================================================
def optimum(closure, dP, r_max, tau_c, H=None, B_jet=6.0, theta=0.0):
    """
    Return dict with d* [m], V0 [m/s], Q* [m^3/s] (and L/min), Re_d, We.
    If H is None: idealized H/d -> small (no free-jet decay), impingement at exit speed.
    If H given (>~6 d): free jet spreads; centerline Vc = B_jet V0 d / H, effective
      impinging width d_eff = H/B_jet, momentum Vc d_eff = V0 d so Re_eff = Re_d.
      tau_w(r) = C_tau rho Vc^2 Re_d^p (d_eff/r)^n.
    Oblique angle theta: normal velocity component V0 cos(theta) (Region-II sec-theta law).
    """
    V0 = V0_from_dP(dP) * np.cos(theta)
    Ct, p, n = closure.C_tau, closure.p, closure.n

    if H is None:
        # tau_c = Ct rho V0^{2+p} nu^{-p} d^{p+n} r_max^{-n}
        # d* = [ tau_c r_max^n nu^p / (Ct rho V0^{2+p}) ]^{1/(p+n)}
        d_star = (tau_c * r_max**n * NU**p / (Ct * RHO * V0**(2 + p)))**(1.0 / (p + n))
    else:
        # spread jet: Vc = B V0 d / H, d_eff = H/B, Re_eff = Re_d = V0 d/nu
        # tau_w(r_max) = Ct rho Vc^2 Re_d^p (d_eff/r_max)^n
        #  = Ct rho (B V0 d/H)^2 (V0 d/nu)^p (H/(B r_max))^n
        #  = Ct rho V0^{2+p} nu^{-p} d^{2+p} B^{2-n} H^{n-2} r_max^{-n}
        # set = tau_c, solve d^{2+p}:
        d_star = (tau_c * r_max**n * NU**p * H**(2 - n) /
                  (Ct * RHO * V0**(2 + p) * B_jet**(2 - n)))**(1.0 / (2 + p))

    Q = flow_rate(d_star, V0_from_dP(dP))   # flow set by exit speed (not cos-reduced) & d*
    Re_d = V0_from_dP(dP) * d_star / NU
    We = RHO * V0_from_dP(dP)**2 * d_star / SIGMA
    return dict(closure=closure.label, dP=dP, d_star=d_star, V0=V0_from_dP(dP),
                Q=Q, Q_Lmin=Q * 6e4, Re_d=Re_d, We=We, H=H, theta=theta)


def calibrate_turbulent(nt, commercial, r_max, tau_c):
    """
    Calibrate the turbulent prefactor C_P so the model reproduces a commercial
    operating point: nozzle d_c at pressure dP_c cleans r_max at tau_c when placed
    at standoff H_c (H/d in the seated range). Returns a calibrated Closure.
    commercial = dict(d=..., dP=..., H=..., theta=...)
    """
    p = -0.5
    d_c, dP_c, H_c = commercial["d"], commercial["dP"], commercial.get("H", None)
    th = commercial.get("theta", 0.0)
    V0 = V0_from_dP(dP_c) * np.cos(th)
    Re_d = V0_from_dP(dP_c) * d_c / NU
    if H_c is None:
        geom = (d_c / r_max)**nt
        C_P = tau_c / (RHO * V0**2 * Re_d**p * geom)
    else:
        B = 6.0
        Vc = B * V0 * d_c / H_c
        d_eff = H_c / B
        C_P = tau_c / (RHO * Vc**2 * Re_d**p * (d_eff / r_max)**nt)
    return Closure(C_tau=C_P, p=p, n=nt, label=f"turbulent (r^-{nt:g}, calibrated)")


def cfd_closure(C_P, n_t, label=None):
    """
    Turbulent closure built DIRECTLY from CFD-fitted Phares coefficients -- no
    commercial-point calibration. The axisymmetric simpleFoam/k-omega-SST solve
    (cfd_openfoam/) fits, over the developed wall-jet band,
        (tau_w / rho V0^2) Re_d^{1/2} = C_P (r/d)^{-n_t},
    which is the model's common form with C_tau = C_P, p = -1/2, n = n_t. This is
    the honest resolution of the review blocker: C_P and n_t are COMPUTED, not
    calibrated to the single commercial datum.
    """
    return Closure(C_tau=C_P, p=-0.5, n=n_t,
                   label=label or f"turbulent (CFD-fitted, r^-{n_t:g})")


# ---- oscillation: principled uniformity factor from the sinusoidal dwell law ----
def kappa_sinusoid():
    """
    Sinusoidal sweep phi(t)=alpha sin(2 pi f t): the impingement point follows the
    arcsine distribution, density rho(s) = 1/(pi sqrt(A^2 - s^2)), A = H alpha.
    The least-served point is the centre s=0, density 1/(pi A), vs a uniform sweep's
    1/(2A). Ratio (worst-served / uniform) = 2/pi. So kappa = 2/pi ~ 0.637.
    """
    return 2.0 / np.pi
