"""ball_hybrid_scaling.py -- scaling estimate for the ball -> inertia-dilatant-hybrid transition
(derivations/ball_hybrid_scaling_2026-09-25.md; PREREG_DILSIZE_2026-09-25.md).

    python derivations/ball_hybrid_scaling.py  ->  out/ball_hybrid_scaling.txt

Four candidate control groups, each computed with the EXISTING impact_stress modules (nothing retyped):
  G_act   = <tau_wall>(t*=1) / sigma*       shear activation (dilatancy.lambda_shear's quantity; the
                                             wall shear from stress.wall_shear_mean)
  De      = tau_p(L_contact) / (D0/U0)      drainage Deborah number (pore_relaxation.drainage_number)
  S       = p_dil_ceiling / (rho U0^2)      dilatant-strength ratio (capacity K_p dphi/phi vs inertia)
  ID      = [M26] Eq. 6                     inertial-dilatant number (inertial_dilatant.ID)
for (i) the two lab materials at their lab operating points, and (ii) a size ladder at FIXED chemistry,
under explicitly declared conventions (the 2026-09-19 lesson: the size exponents are convention-
dependent, so each convention is printed rather than chosen).
Viscosity for G_act is a declared BRACKET (the lab-silica eta_0 comes from a retracted fit, si_fit).
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
IMPACT = os.path.join(os.path.dirname(ROOT), "impact_stress")
sys.path.insert(0, IMPACT)
import inertial_dilatant as idm  # noqa: E402
import pore_relaxation as pr     # noqa: E402
import stress as st              # noqa: E402
from materials import Fluid      # noqa: E402
from dilatancy import SIGMA_STAR_LITERATURE  # noqa: E402

lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


# lab operating points (impact_stress/run_pore_relaxation.py; cornstarch = [M26])
LAB = {"silica": dict(D0=2.5e-3, U0=2.3), "cornstarch": dict(D0=4.0e-3, U0=5.1)}
ETA_BRACKET = (0.01, 0.1)          # Pa s, declared bracket for eta_0 of the unthickened suspension
# Pa: min / repo / max of dilatancy.SIGMA_STAR_LITERATURE, read not retyped. The middle 'repo' value is the
# retracted si_fit (failed fit gate); it is kept only because silica_poro uses it for P_conf.
_SS = SIGMA_STAR_LITERATURE
SIGMA_STAR_SILICA = (min(_SS.values()), next(v for k, v in _SS.items() if "repo" in k), max(_SS.values()))


def tau_wall(rho: float, eta: float, U0: float, D0: float) -> float:
    fl = Fluid(name="susp", rho=rho, c=1500.0, eta=eta, gamma=0.072, provenance="scaling bracket")
    imp = st.Impact(drop=fl, U=U0, R=D0 / 2, target=st.rigid_target())
    return float(st.wall_shear_mean(imp, 1.0))


def groups(p: pr.PoroSystem, D0: float, U0: float, K_p: float | None = None) -> dict:
    s = p.sys
    K = p.K_p if K_p is None else K_p
    L = pr.drainage_lengths(D0, s.a, s.phi)["contact"]
    tau_p = pr.pore_relaxation_time(L, s.kappa, K, s.phi)
    De = pr.drainage_number(tau_p, D0, U0)
    ceil = pr.dilatant_pressure_ceiling(s, K)
    rho_u2 = s.rho * U0**2
    return dict(a=s.a, phi=s.phi, dphi=s.dphi, kappa=s.kappa, K_p=K, De=De, ceil=ceil, S=ceil / rho_u2,
                rho_u2=rho_u2, ID=idm.ID(s, U0, D0),
                tw=[tau_wall(s.rho, e, U0, D0) for e in ETA_BRACKET])


def main() -> int:
    say("ball -> hybrid scaling: four candidate groups (impact_stress modules; see docstring)")
    say(f"eta_0 bracket for the shear-activation group: {ETA_BRACKET} Pa s (declared, not measured)")
    say("\n=== 1. the two lab materials at their lab operating points ===")
    mats = {"silica phi 0.44": (pr.silica_poro(phi=0.44), LAB["silica"]),
            "silica phi 0.50": (pr.silica_poro(phi=0.50), LAB["silica"]),
            "cornstarch phi 0.42": (pr.cornstarch_poro(phi=0.42), LAB["cornstarch"])}
    for name, (p, op) in mats.items():
        g = groups(p, **op)
        glo = groups(p, K_p=p.K_p_low, **op)
        ghi = groups(p, K_p=p.K_p_high, **op)
        say(f"  {name}: a {g['a']*1e6:.2f} um, dphi {g['dphi']:.3f}, kappa {g['kappa']:.2e} m^2, "
            f"K_p {g['K_p']:.3g} Pa [{p.K_p_low:.3g}, {p.K_p_high:.3g}]")
        say(f"     De (contact drainage) {g['De']:.3g} [{ghi['De']:.3g}, {glo['De']:.3g}]   "
            f"S = ceiling/rho U0^2 = {g['ceil']:.3g}/{g['rho_u2']:.3g} = {g['S']:.3g} "
            f"[{glo['S']:.3g}, {ghi['S']:.3g}]   ID {g['ID']:.3g}")
        say(f"     <tau_wall>(t*=1) = {g['tw'][0]:.3g} .. {g['tw'][1]:.3g} Pa over the eta bracket")
    say(f"  silica sigma* bracket {SIGMA_STAR_SILICA} Pa (dilatancy.SIGMA_STAR_LITERATURE); cornstarch sigma* is NOT "
        "on disk as a measured value, so G_act is printed only for silica:")
    g = groups(pr.silica_poro(phi=0.44), **LAB["silica"])
    for ss in SIGMA_STAR_SILICA:
        say(f"     sigma* {ss:6.1f} Pa: G_act = {g['tw'][0] / ss:.3g} .. {g['tw'][1] / ss:.3g}")
    say("  silica phi 0.44, De and S with P_conf (the pressure in K_p) set to each sigma* of the bracket:")
    p0 = pr.silica_poro(phi=0.44)
    for ss in SIGMA_STAR_SILICA:
        K = pr.network_modulus_edl(ss, p0.sys.a)
        pp = pr.PoroSystem(sys=p0.sys, K_p=K, K_p_low=K, K_p_high=K, K_p_origin="sigma* bracket")
        gg = groups(pp, **LAB["silica"])
        say(f"     P_conf {ss:6.1f} Pa: K_p {K:.3g} Pa, De {gg['De']:.3g}, S {gg['S']:.3g}")
    say("  NOTE: P_conf = sigma* equates a SHEAR onset stress with a particle PRESSURE; the two differ by an O(1)")
    say("  friction coefficient (dilatancy.py), which moves K_p, De and S by that O(1) factor, not by decades.")

    say("\n=== 2. size ladder at FIXED chemistry (silica anchor a0 = 0.8 um, phi 0.44, dphi 0.02; lab impact) ===")
    say("  conventions (each printed, none chosen):")
    say("    sigma*: [F] fixed repulsive force F* -> sigma* ~ a^-2 ;  [D] Derjaguin F* ~ a -> sigma* ~ a^-1")
    say("    K_p   : EDL K_p = P_conf a/(6 lambda_D) with P_conf = sigma* (silica_poro's own rule) -> K_p ~ a sigma*")
    say("    drainage length: contact patch 0.49 D0 (PORE_PHYSICS §3's choice for this question)")
    say("  a/um   kappa       |  [F]: sigma*  K_p     De        S         |  [D]: sigma*  K_p     De        S")
    base = pr.silica_poro(phi=0.44)
    a0, s0 = base.sys.a, SIGMA_STAR_SILICA[1]
    for a in (0.4e-6, 0.8e-6, 2e-6, 5e-6, 10e-6, 20e-6):
        sysd = idm.DilatantSystem(name="ladder", phi=base.sys.phi, phi_c=base.sys.phi_c, a=a, rho=base.sys.rho)
        row = []
        for exp in (-2, -1):
            ss = s0 * (a / a0) ** exp
            K = pr.network_modulus_edl(ss, a)
            p = pr.PoroSystem(sys=sysd, K_p=K, K_p_low=K, K_p_high=K, K_p_origin="ladder")
            gg = groups(p, **LAB["silica"])
            row.append((ss, K, gg["De"], gg["S"]))
        say(f"  {a*1e6:5.1f}  {sysd.kappa:.2e}  |  {row[0][0]:7.3g} {row[0][1]:7.3g} {row[0][2]:9.3g} {row[0][3]:9.3g} "
            f"|  {row[1][0]:7.3g} {row[1][1]:7.3g} {row[1][2]:9.3g} {row[1][3]:9.3g}")
    say("  exponents (d ln / d ln a): kappa +2; [F] sigma* -2, K_p -1, De = L^2 mu/(kappa K_p) tau_i^-1 ~ a^-1, "
        "S ~ K_p ~ a^-1;  [D] sigma* -1, K_p 0, De ~ a^-2, S ~ a^0")

    say("\n=== 3. CV at fixed number-mean size (d10), same chemistry ===")
    say("  kappa uses d32 (Carman-Kozeny): d32/d10 = (1+CV^2)^2 -> kappa x (1+CV^2)^4; sigma* x 0.735 at CV 0.5")
    say("  (sigma_star_geometry, fixed d10); K_p's CV dependence is not measured and is held fixed here.")
    for cv in (0.2, 0.5):
        k = (1 + cv**2) ** 4
        say(f"  CV {cv}: kappa x{k:.3f} -> De x{1/k:.3f}, ID x{k:.3f}; window: x = phi/phi_0 falls (cvmap)")
    with open(os.path.join(ROOT, "out", "ball_hybrid_scaling.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
