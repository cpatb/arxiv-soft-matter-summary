"""diag_kappa_sweep.py — is the pore pressure really independent of permeability?

THE CLAIM UNDER TEST
--------------------
`derivations/dilatancy_sign_mechanism.md` §7 says that where the network-modulus ceiling
binds, permeability leaves the problem entirely, so mechanism M1 is not merely
sub-dominant but ABSENT. That is a STATIC argument: Eq. (16), `p_dil^max = K_p dphi/phi`,
simply has no `kappa` in it.

A static argument is not enough, because the ceiling has to be REACHED. If the pressure
builds diffusively then early times are still `kappa`-dependent and M1 would come back in
through the transient. This script settles it dynamically: it integrates `poro_impact`
at the lab's own physical numbers and sweeps `kappa` over ten decades, holding everything
else fixed, then reports the pressure at a fixed impact strain.

It exists because the claim was first checked by a referee in a transcript, and a number
that lives only in a transcript is not reproducible. `python diag_kappa_sweep.py`.

WHAT TO READ
------------
  * `p(gamma=1)/ceiling` flat across `kappa` ⇒ M1 absent, dynamically as well.
  * the SIZE LADDER rows are the ones that matter for the experiment: 0.8 -> 8 um is
    100x in `kappa` and is the ladder E-4 proposes.
  * `Lambda` is printed beside every row, because the flat region is exactly the region
    where `Lambda << 1` and the honest statement is conditional on that.
"""

from __future__ import annotations

import os
import sys


HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "impact_stress")))

import poro_impact as PI                     # noqa: E402
import pore_relaxation as PR                 # noqa: E402

PHI, PHI_C = 0.44, 0.42
U0, D0 = 3.0, 2.5e-3
MU_W = 1.0e-3
GAMMA_READ = 1.0        # the strain at which the pressure is read: one impact strain


def row(kappa: float, k_p: float, label: str) -> dict:
    rho = 2000.0 * PHI + 1000.0 * (1 - PHI)
    n = PI.physical_numbers(kappa=kappa, mu_w=MU_W, k_drained=k_p, rho=rho,
                            u0=U0, d0=D0, ell=D0)
    r = PI.solve(lam=n["Lambda"], s_dil=n["S"], dphi0=PHI - PHI_C, phi0=PHI,
                 t_end=GAMMA_READ, n_nodes=201, rtol=1e-10, atol=1e-14)
    p_at = abs(r.p_base_final) * rho * U0**2          # back to Pa
    ceiling = k_p * (PHI - PHI_C) / PHI               # Eq. (16), Pa
    return {"label": label, "kappa": kappa, "Lambda": n["Lambda"], "S": n["S"],
            "p_Pa": p_at, "ceiling_Pa": ceiling, "ratio": p_at / ceiling,
            "dphi_final": r.dphi_final}


def main() -> int:
    ps = PR.silica_poro(phi=PHI, phi_c=PHI_C)
    k_p = ps.K_p
    kappa0 = ps.sys.kappa
    print(f"lab silica: phi={PHI}, dphi={PHI - PHI_C}, K_p={k_p:.4g} Pa, "
          f"kappa={kappa0:.4g} m^2, read at gamma={GAMMA_READ}")
    print(f"{'kappa/kappa0':>12s} {'kappa[m^2]':>12s} {'Lambda':>10s} "
          f"{'p(Pa)':>10s} {'ceiling':>9s} {'p/ceiling':>10s} {'dphi_end':>9s}")
    rows = []
    for e in range(0, 11):
        r = row(kappa0 * 10.0**e, k_p, f"x1e{e}")
        rows.append(r)
        print(f"{10.0**e:12.0e} {r['kappa']:12.3e} {r['Lambda']:10.2e} "
              f"{r['p_Pa']:10.4g} {r['ceiling_Pa']:9.4g} {r['ratio']:10.4f} "
              f"{r['dphi_final']:9.6f}")

    flat = [r for r in rows if r["Lambda"] < 0.1]
    if flat:
        lo = min(r["ratio"] for r in flat)
        hi = max(r["ratio"] for r in flat)
        print(f"\nover the {len(flat)} rungs with Lambda < 0.1 the ratio spans "
              f"{lo:.6f} to {hi:.6f}, i.e. {hi / lo:.6f}x "
              f"({(10.0 ** (len(flat) - 1)):.0e} in kappa)")

    print("\n--- the size ladder E-4 proposes, at fixed phi, dphi, U0, D0 ---")
    print(f"{'a[um]':>7s} {'kappa[m^2]':>12s} {'Lambda':>10s} {'p(Pa)':>10s} "
          f"{'p/p(0.8um)':>11s}")
    base = None
    for a_um in (0.8, 2.5, 8.0):
        a = a_um * 1e-6
        kap = (1 - PHI) ** 3 * a**2 / (180 * PHI**2)
        # K_p for the EDL contact goes as the particle size, so it moves too
        kp_a = PR.network_modulus_edl(PR.P_CONF_SILICA, a, PR.DEBYE_DEFAULT_M)
        r = row(kap, kp_a, f"a={a_um}um")
        base = r["p_Pa"] if base is None else base
        print(f"{a_um:7.1f} {kap:12.3e} {r['Lambda']:10.2e} {r['p_Pa']:10.4g} "
              f"{r['p_Pa'] / base:11.4f}")
    print("\nNOTE the size ladder moves K_p as well as kappa (K_p ~ a for a")
    print("double-layer contact), so its p column is NOT a pure permeability test --")
    print("that is what the kappa sweep above is for.")

    out = os.path.join(HERE, "out", "kappa_sweep.txt")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("kappa_ratio,kappa,Lambda,S,p_Pa,ceiling_Pa,ratio,dphi_final\n")
        for r, e in zip(rows, range(0, 11)):
            fh.write(f"{10.0 ** e:g},{r['kappa']:.6e},{r['Lambda']:.6e},"
                     f"{r['S']:.6e},{r['p_Pa']:.6e},{r['ceiling_Pa']:.6e},"
                     f"{r['ratio']:.6e},{r['dphi_final']:.6e}\n")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
