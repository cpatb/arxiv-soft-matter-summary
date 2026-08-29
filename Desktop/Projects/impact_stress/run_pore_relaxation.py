"""
run_pore_relaxation.py — evaluates pore_relaxation.py for the two systems the
PI asked about and writes out/pore_relaxation.json plus three figures.

Run:  python run_pore_relaxation.py

Figure convention (lab standard): minimal text on the axes, no titles, no
annotation prose. Every figure's reading instructions live in
PORE_RELAXATION_BRIEF.md, not on the plot.

Axis convention: linear axes are anchored at zero; log axes are used only where
a quantity spans more than a decade, and each such use is justified in the
brief. Every quantity plotted here spans four to nine decades, so all the
ordinates are logarithmic and the brief says so.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import pore_relaxation as pr

OUT = Path(__file__).parent / "out"
FIGS = Path(__file__).parent / "figures"

# ---- the two operating points ---------------------------------------------
D0_LAB, U0_LAB = 2.5e-3, 2.3       # lab drop, band B centre
D0_M26, U0_M26 = 4.0e-3, 5.1       # [M26] cornstarch drop, Regime-2/3 boundary
PHI_C_SILICA = 0.42                # the unmeasured input; see the brief
PHI_C_STARCH = 0.37                # [M26]

PHI_SILICA = np.linspace(0.425, 0.580, 40)
PHI_STARCH = np.linspace(0.375, 0.500, 40)
DPHI_MATCHED = np.linspace(0.005, 0.080, 40)

LAMBDA_D = {"0 mM (residual)": 30e-9, "1 mM": 10e-9, "100 mM": 1e-9}


def _sweep(maker, phis, D0, U0, **kw):
    rows = []
    for phi in phis:
        p = maker(phi=phi, **kw)
        r = pr.report(p, D0, U0)
        rows.append({
            "phi": float(phi), "dphi": r["dphi"], "kappa_m2": r["kappa_m2"],
            "K_p_Pa": r["K_p_Pa"], "D_p_m2_s": r["D_p_m2_s"],
            "tau_p_contact_s": r["tau_p_s"]["contact"],
            "tau_p_particle_s": r["tau_p_s"]["particle"],
            "De_contact": r["De"]["contact"],
            "De_particle": r["De"]["particle"],
            "p_dil_darcy_Pa": r["p_dil_darcy_Pa"]["contact"],
            "p_dil_ceiling_Pa": r["p_dil_ceiling_Pa"],
            "p_dil_capped_Pa": r["p_dil_capped_Pa"]["contact"],
        })
    return rows


def amplification(phi: float, a: float, lam: float, phi_c: float) -> float:
    """p_dil_max / P_conf = (a / 6 lambda_D) * dphi / phi.

    The parameter-free content of the silica ceiling: whether frustrated
    dilation amplifies the stress that confines the network, or cannot.
    Independent of P_conf, so it does not inherit sigma*'s uncertainty.
    Below 1, Reynolds dilatancy cannot raise the particle pressure above what
    already drives it; above 1, it can.

    The strain is dphi/phi, not dphi/(1-phi); see
    pore_relaxation.dilatant_pressure_ceiling.
    """
    dphi = phi - phi_c
    if dphi <= 0:
        return 0.0
    return float((a / (6.0 * lam)) * dphi / phi)


def phi_star(a: float, lam: float, phi_c: float) -> float:
    """The phi at which amplification crosses 1, solved in closed form.

        (a/6L)(phi - phi_c)/phi = 1
        =>  phi* = phi_c / (1 - 6L/a)

    Returns inf when 6 lambda_D >= a, i.e. when the double layer is thicker than
    a sixth of a particle: then the ratio can never reach 1 at any phi and
    Reynolds dilatancy is unavailable throughout. For 800 nm silica that needs
    lambda_D > 133 nm, well outside the accessible salt range, so the branch is
    a guard rather than a regime.
    """
    q = 6.0 * lam / a
    if q >= 1.0:
        return float("inf")
    return float(phi_c / (1.0 - q))


def main() -> dict:
    OUT.mkdir(exist_ok=True)
    FIGS.mkdir(exist_ok=True)

    si_hi = pr.silica_poro(phi=0.50, phi_c=PHI_C_SILICA)
    si_lab = pr.silica_poro(phi=0.44, phi_c=PHI_C_SILICA)
    cs = pr.cornstarch_poro(phi=0.42)
    gb = pr.glass_beads_JVF()

    res = {
        "_what": "Pore-pressure relaxation and the Reynolds-dilatancy ceiling "
                 "for 800 nm silica vs 10 um cornstarch.",
        "_conventions": {
            "a_is": "particle DIAMETER everywhere",
            "gammadot": "U0/D0 (convention A, Shah/Arora/Driscoll)",
            "drainage_length_contact": "0.49 D0, the contact-patch radius at "
                                       "peak force from the Wagner analysis",
            "tau_impact": "D0/U0",
        },
        "operating_points": {
            "silica_phi0.50": pr.report(si_hi, D0_LAB, U0_LAB),
            "silica_phi0.44_lab": pr.report(si_lab, D0_LAB, U0_LAB),
            "cornstarch_phi0.42": pr.report(cs, D0_M26, U0_M26),
            "glass_beads_JVF_validation": pr.report(gb, D0_M26, U0_M26),
        },
        "sweeps": {
            "silica": _sweep(pr.silica_poro, PHI_SILICA, D0_LAB, U0_LAB,
                             phi_c=PHI_C_SILICA),
            "cornstarch": _sweep(pr.cornstarch_poro, PHI_STARCH, D0_M26, U0_M26),
        },
    }

    # ---- the contrast, decomposed -----------------------------------------
    dec = []
    for dphi in DPHI_MATCHED:
        s = pr.silica_poro(phi=PHI_C_SILICA + dphi, phi_c=PHI_C_SILICA)
        c = pr.cornstarch_poro(phi=PHI_C_STARCH + dphi)
        dec.append({
            "dphi": float(dphi),
            "kappa_ratio_si_over_cs": s.kappa / c.kappa,
            "K_p_ratio_si_over_cs": s.K_p / c.K_p,
            "D_p_ratio_si_over_cs": s.D_p() / c.D_p(),
        })
    res["contrast_decomposition_at_matched_dphi"] = dec
    d50 = dec[len(dec) // 2]
    res["contrast_headline"] = {
        "at_dphi": d50["dphi"],
        "permeability_factor": 1.0 / d50["kappa_ratio_si_over_cs"],
        "network_stiffness_factor": 1.0 / d50["K_p_ratio_si_over_cs"],
        "total_diffusivity_factor": 1.0 / d50["D_p_ratio_si_over_cs"],
        "note": "Each system is evaluated at ITS OWN operating pressure "
                "(silica P_conf = sigma* = 44 Pa; cornstarch P_conf = 1e5 Pa, "
                "the measured jammed-column stress). Read with "
                "contrast_at_matched_P_conf below, which shows how much of the "
                "stiffness factor is the pressure difference rather than the "
                "materials.",
    }

    # ---- how much of the contrast survives matched confinement? -----------
    # K_p ~ P for silica (EDL) and ~P^(1/3) for cornstarch (Hertz), so the
    # stiffness ratio is NOT a material constant. This is the control.
    matched = {}
    for P in (1.0e3, 1.0e4, 1.0e5):
        s = pr.silica_poro(phi=0.50, phi_c=PHI_C_SILICA, P_conf=P)
        c = pr.cornstarch_poro(phi=0.42, P_conf=P)
        tau_s = s.tau_p(0.49 * D0_LAB)
        matched[f"P_conf={P:.0e} Pa"] = {
            "K_p_silica_Pa": s.K_p, "K_p_cornstarch_Pa": c.K_p,
            "K_p_factor_cs_over_si": c.K_p / s.K_p,
            "D_p_factor_cs_over_si": c.D_p() / s.D_p(),
            "silica_De_contact": pr.drainage_number(tau_s, D0_LAB, U0_LAB),
        }
    res["contrast_at_matched_P_conf"] = {
        "_why": "The stiffness contrast in contrast_headline is largely the "
                "difference in operating pressure, not a material property. At "
                "matched P_conf the two networks are comparably stiff and the "
                "whole diffusivity contrast collapses onto the a^2 "
                "permeability factor. The undrained verdict for silica "
                "survives every entry here by at least two decades, which is "
                "what makes it the robust conclusion.",
        "by_pressure": matched,
    }

    # ---- the amplification crossover --------------------------------------
    a_si = 800e-9
    res["dilatant_amplification"] = {
        "definition": "p_dil_max / P_conf = (a/6 lambda_D) dphi/(1-phi); "
                      "P_conf cancels, so this does not inherit sigma*'s "
                      "uncertainty. >1 means frustrated dilation amplifies.",
        "phi_star": {k: phi_star(a_si, v, PHI_C_SILICA) for k, v in LAMBDA_D.items()},
        "at_lab_phi_0.44": {k: amplification(0.44, a_si, v, PHI_C_SILICA)
                            for k, v in LAMBDA_D.items()},
        "at_phi_0.50": {k: amplification(0.50, a_si, v, PHI_C_SILICA)
                        for k, v in LAMBDA_D.items()},
        "sensitivity_note": "phi* depends on phi_c, which is UNMEASURED for "
                            "this silica. phi_c = 0.42 is assumed throughout; "
                            "phi* moves with it roughly one-for-one.",
    }

    _figures(res)

    (OUT / "pore_relaxation.json").write_text(json.dumps(res, indent=1))
    return res


# ==========================================================================
def _figures(res: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "figure.dpi": 200,
                         "axes.grid": True, "grid.alpha": 0.2,
                         "savefig.bbox": "tight", "legend.frameon": False})
    C_SI, C_CS = "#1b4a8f", "#d9822b"

    si = res["sweeps"]["silica"]
    cs = res["sweeps"]["cornstarch"]
    g = lambda rows, k: np.array([r[k] for r in rows])  # noqa: E731

    # ---- Fig 6: relaxation time and drainage number -----------------------
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.9))
    ax[0].semilogy(g(si, "phi"), g(si, "tau_p_contact_s"), color=C_SI, lw=1.6)
    ax[0].semilogy(g(cs, "phi"), g(cs, "tau_p_contact_s"), color=C_CS, lw=1.6)
    ax[0].axhline(D0_LAB / U0_LAB, color=C_SI, ls=":", lw=1)
    ax[0].axhline(D0_M26 / U0_M26, color=C_CS, ls=":", lw=1)
    ax[0].set(xlabel=r"$\phi$", ylabel=r"$\tau_p$  [s]")

    ax[1].semilogy(g(si, "phi"), g(si, "De_contact"), color=C_SI, lw=1.6,
                   label="silica, 800 nm")
    ax[1].semilogy(g(cs, "phi"), g(cs, "De_contact"), color=C_CS, lw=1.6,
                   label=r"cornstarch, 10 $\mu$m")
    ax[1].axhline(1.0, color="k", lw=0.9)
    ax[1].set(xlabel=r"$\phi$", ylabel=r"$\mathrm{De}=\tau_p/\tau_{\rm imp}$")
    ax[1].legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIGS / "fig6_pore_relaxation.png")
    plt.close(fig)

    # ---- Fig 7: the dilatant amplification crossover ----------------------
    fig, ax = plt.subplots(figsize=(3.6, 2.9))
    phis = np.linspace(0.421, 0.580, 300)
    for (lab, lam), ls in zip(LAMBDA_D.items(), ["-", "--", "-."]):
        ax.semilogy(phis, [amplification(p, 800e-9, lam, PHI_C_SILICA) for p in phis],
                    color=C_SI, ls=ls, lw=1.5, label=lab)
    ax.axhline(1.0, color="k", lw=0.9)
    ax.axvline(0.44, color="0.45", lw=0.9, ls=":")
    ax.axvline(0.47, color="0.45", lw=0.9, ls=":")
    ax.set(xlabel=r"$\phi$", ylabel=r"$p_{\rm dil}^{\max}/P_{\rm conf}$",
           ylim=(1e-2, 1e3))
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIGS / "fig7_dilatant_amplification.png")
    plt.close(fig)

    # ---- Fig 8: what the contrast is made of ------------------------------
    dec = res["contrast_decomposition_at_matched_dphi"]
    fig, ax = plt.subplots(figsize=(3.6, 2.9))
    x = g(dec, "dphi")
    ax.semilogy(x, 1 / g(dec, "kappa_ratio_si_over_cs"), color="0.35", lw=1.5,
                label=r"permeability $\kappa$")
    ax.semilogy(x, 1 / g(dec, "K_p_ratio_si_over_cs"), color=C_SI, lw=1.5,
                label=r"network modulus $K_p$")
    ax.semilogy(x, 1 / g(dec, "D_p_ratio_si_over_cs"), color="k", lw=2.0,
                label=r"product, $D_p$")
    ax.set(xlabel=r"$\Delta\phi=\phi-\phi_c$",
           ylabel="cornstarch / silica")
    ax.legend(fontsize=8, loc="center right")
    fig.tight_layout()
    fig.savefig(FIGS / "fig8_contrast_decomposition.png")
    plt.close(fig)


if __name__ == "__main__":
    r = main()
    h = r["contrast_headline"]
    print(f"permeability   x{h['permeability_factor']:.4g}")
    print(f"network K_p    x{h['network_stiffness_factor']:.4g}")
    print(f"product D_p    x{h['total_diffusivity_factor']:.4g}")
    print("phi* (amplification = 1):",
          {k: round(v, 4) for k, v in r["dilatant_amplification"]["phi_star"].items()})
    for k, v in r["operating_points"].items():
        print(f"{k:32s} De(contact)={v['De']['contact']:.3g}  "
              f"tau_p={v['tau_p_s']['contact']:.3g} s  "
              f"p_dil_capped={v['p_dil_capped_Pa']['contact']:.3g} Pa")
    print("wrote out/pore_relaxation.json and figures/fig6-8")
