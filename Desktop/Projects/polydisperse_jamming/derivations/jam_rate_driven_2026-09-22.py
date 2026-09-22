"""Rate-driven jamming threshold for the solver's own activation (amendment to PREREG_JAMWIN).

The solver activates on sigma = eta(fhat) * gd (lagged viscosity x local strain rate), not
on an imposed stress. Steady states solve  f = F(f) := exp(-(sigma*/(eta(f) gd))^delta),
and a material point relaxing from f = 0 reaches f_c iff F(f) > f on all of [0, f_c)
(no stable fixed point below f_c; F is increasing in f because eta is). Writing
g(f) = (ln 1/f)^(-1/delta), F(f) > f  <=>  eta(f) gd / sigma* > g(f), so

    gd_c / (sigma*/eta_s) = max_{0 < f < f_c}  g(f) / (eta(f)/eta_s).          (R1)

The maximiser f* is the BOTTLENECK. d gd_c / d delta has the sign of ln ln(1/f*)
(the same algebra as eq. 9 of jam_stress_2026-09-22.md, evaluated at f* instead of f_c),
so the reversal sits at f* = 1/e, not at f_c = 1/e.

Viscosity exactly as wc_drop_impact.c eta_wc(): eta_s * eps^-MP_EXPO, eps = 1 - PHI/phiJ(f),
floored at EPS_MIN, phiJ(f) = phi0 - (phi0 - phim) f. EPS_SOFT is off by default.
"""
import numpy as np

PHI0, PHIM, EPS_MIN, MP = 0.62, 0.496, 0.01, 2.0
ARMS = {"A": 0.589, "B": 0.527}
DELTAS = (0.7, 1.3169, 2.6)


def eta_rel(f, phi):
    phij = PHI0 - (PHI0 - PHIM) * f
    eps = np.maximum(1.0 - phi / phij, EPS_MIN)
    return eps ** -MP


def check_mumax_inert():
    """The solver's outer clamp min(e, mu_max_wc), mu_max_wc = min(eta_s/EPS_MIN^2, MUMAX_WC),
    is inert here iff eta_s_code/EPS_MIN^2 < MUMAX_WC. Anchor: the jamwin banners print
    'eta_max=5.47 code' (out/jamwin/*/banner.txt)."""
    eta_s, rho, u0, d0, mumax = 0.003939, 1200.0, 2.0, 3.0e-3, 1e9   # cs_fit, RHO_L, U0, D0
    eta_s_code = eta_s / (rho * u0 * d0)
    cap = eta_s_code / EPS_MIN ** 2
    print(f"mu_max_wc check: eta_s_code {eta_s_code:.4g}, eta_s/EPS_MIN^2 = {cap:.3g} code "
          f"(banner prints 5.47) < MUMAX_WC {mumax:g}: {'INERT' if cap < mumax else 'BINDS'}")
    assert abs(cap - 5.47) < 0.01 and cap < mumax


def main():
    print(__doc__.split("\n")[0])
    check_mumax_inert()
    print(f"phi0 {PHI0} phim {PHIM} chi {PHIM / PHI0:.3f}  EPS_MIN {EPS_MIN}  MP_EXPO {MP}")
    rows = {}
    for arm, phi in ARMS.items():
        x, chi = phi / PHI0, PHIM / PHI0
        fc = (1 - x) / (1 - chi)
        f = np.linspace(1e-6, fc * (1 - 1e-9), 2_000_001)
        print(f"\narm {arm}: x = {x:.4f}, f_c = {fc:.4f}")
        print("   delta   f*      gd_c*eta_s/sigma*   [stress-driven sigma_jam/sigma* = g(f_c)]")
        for d in DELTAS:
            g = np.log(1 / f) ** (-1 / d)
            r = g / eta_rel(f, phi)
            i = int(np.argmax(r))
            gfc = np.log(1 / fc) ** (-1 / d)
            rows[(arm, d)] = r[i]
            print(f"   {d:6.4f}  {f[i]:.4f}   {r[i]:.5e}          [{gfc:.4f}]")
        seq = [rows[(arm, d)] for d in DELTAS]
        trend = "RISES" if seq[0] < seq[1] < seq[2] else ("FALLS" if seq[0] > seq[1] > seq[2]
                                                          else "NON-MONOTONE")
        print(f"   rate-driven threshold gd_c {trend} with delta "
              f"=> J should {'FALL' if trend == 'RISES' else 'RISE' if trend == 'FALLS' else '?'}"
              " with delta at fixed forcing")
    # brute-force check of (R1): integrate df/dt = F(f) - f from 0 at gd just above/below
    print("\nbrute-force check: relax df/dt = F(f) - f from f = 0 (sigma* = eta_s = 1)")
    for (arm, d), gdc in rows.items():
        phi = ARMS[arm]
        fc = (1 - phi / PHI0) / (1 - PHIM / PHI0)
        res = []
        for fac in (0.98, 1.02):
            gd, f = gdc * fac, 0.0
            for _ in range(200_000):
                F = np.exp(-(1.0 / (eta_rel(f, phi) * gd)) ** d) if f < 1 else 1.0
                f += 0.01 * (F - f)
            res.append(f >= fc)
        print(f"   {arm} delta {d}: reaches f_c at 0.98 gd_c? {res[0]}   at 1.02 gd_c? {res[1]}"
              f"   {'PASS' if (not res[0] and res[1]) else 'FAIL'}")


if __name__ == "__main__":
    main()
