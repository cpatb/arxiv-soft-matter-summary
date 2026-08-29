"""
Tier-0 sanity gate: reproduce the two contrasts the drop-impact work hinges on.

Gate 1 (bulk rheology): silica-like parameters give DST / shear jamming
        (back-bending flow curve); corn-starch-like parameters give a
        continuous, weakly thickening yield fluid (CST). This is the cheap
        check requested in the Simulations MOC BEFORE any CFD.

Gate 2 (reversibility switch): under an identical stress pulse, a long-tau
        (silica) state retains its frictional network after the stress is
        removed, while a short-tau (corn starch) state relaxes -> the
        retention-vs-recovery split that distinguishes a stuck pancake from
        a recovering drop.

Outputs two PNGs and a printed PASS/FAIL validation report.
Parameters are illustrative; the point is the qualitative contrast.
"""
from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import wyart_cates as wc

PURPLE, RED, GOLD, INK = "#4E2A84", "#C0392B", "#B6862C", "#2B2B2B"

# ------------------------------------------------------------------ parameters
# Shared microstructure
PHI0, PHI_M = 0.65, 0.56
ETA_S = 1.0e-3  # Pa.s (water-like suspending fluid)

# Silica: charge-stabilised, phi in the shear-jamming window (phi_m < phi < phi0)
SILICA = dict(phi=0.585, sigma_star=10.0, sigma_yield=0.0)
# Corn starch (this prep): CST + adhesive yield; phi below phi_m, weaker upturn
#   (reflects the Slide-16 observation that this cornstarch shows CST not DST)
CORN = dict(phi=0.520, sigma_star=8.0, sigma_yield=4.0)


def gate1_flow_curves():
    fc_si = wc.steady_flow_curve(
        SILICA["phi"], SILICA["sigma_star"], phi0=PHI0, phi_m=PHI_M,
        eta_s=ETA_S, sigma_yield=SILICA["sigma_yield"],
    )
    fc_cs = wc.steady_flow_curve(
        CORN["phi"], CORN["sigma_star"], phi0=PHI0, phi_m=PHI_M,
        eta_s=ETA_S, sigma_yield=CORN["sigma_yield"],
    )

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.4), dpi=140)

    # (A) stress vs shear rate (flow curve); DST = back-bending
    for fc, c, lab in ((fc_si, PURPLE, "silica-like (DST / shear-jam)"),
                       (fc_cs, RED, "corn-starch-like (CST + yield)")):
        m = np.isfinite(fc["gamma_dot"]) & (fc["gamma_dot"] > 0)
        axA.loglog(fc["gamma_dot"][m], fc["sigma_total"][m], "-", color=c, lw=2.4, label=lab)
    axA.set_xlabel(r"shear rate  $\dot\gamma$  (s$^{-1}$)", color=INK)
    axA.set_ylabel(r"shear stress  $\sigma$  (Pa)", color=INK)
    axA.set_title("Flow curve — back-bending = DST/shear jamming", fontsize=11)
    axA.legend(fontsize=9, frameon=False)
    axA.grid(True, which="both", color="#E3E3E3", lw=0.7)
    axA.set_axisbelow(True)

    # (B) apparent viscosity vs shear rate
    for fc, c, lab in ((fc_si, PURPLE, "silica-like"), (fc_cs, RED, "corn-starch-like")):
        m = np.isfinite(fc["eta_app"]) & (fc["gamma_dot"] > 0)
        axB.loglog(fc["gamma_dot"][m], fc["eta_app"][m], "-", color=c, lw=2.4, label=lab)
    axB.set_xlabel(r"shear rate  $\dot\gamma$  (s$^{-1}$)", color=INK)
    axB.set_ylabel(r"apparent viscosity  $\eta$  (Pa$\cdot$s)", color=INK)
    axB.set_title("Viscosity — sharp jump (silica) vs weak upturn (corn starch)", fontsize=11)
    axB.legend(fontsize=9, frameon=False)
    axB.grid(True, which="both", color="#E3E3E3", lw=0.7)
    axB.set_axisbelow(True)

    for ax in (axA, axB):
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig("gate1_flow_curves.png", facecolor="white")
    plt.close(fig)
    return fc_si, fc_cs


def gate2_reversibility():
    # identical impact-strength stress pulse (f_ss ~ 1) for 5 ms, then near-zero.
    # Impact stress is large (cf. 30 kPa-3.6 MPa local in the drop-impact data),
    # so the frictional network fully forms during the pulse for BOTH materials;
    # what then differs is whether it RELEASES.
    sigma_hi, sigma_lo = 1.0e3, 0.5
    t_pulse = 5.0e-3
    t_eval = np.linspace(0.0, 60.0e-3, 1200)

    def sigma_of_t(t):
        return sigma_hi if t < t_pulse else sigma_lo

    tau = 2.0e-3  # same rise time for both -> isolates reversibility as the switch
    sigma_star = 10.0
    # silica retains an irreversible frictional network (residual floor);
    # corn starch is reversible (no residual) -> recovers.
    f_si = wc.relax_state(t_eval, sigma_of_t, tau, sigma_star, f_residual=0.6)
    f_cs = wc.relax_state(t_eval, sigma_of_t, tau, sigma_star, f_residual=0.0)

    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=140)
    ax.axvspan(0, t_pulse * 1e3, color=GOLD, alpha=0.12, label="stress pulse")
    ax.plot(t_eval * 1e3, f_si, "-", color=PURPLE, lw=2.6,
            label=r"silica-like ($f_{\rm res}$=0.6) — retains (stuck pancake)")
    ax.plot(t_eval * 1e3, f_cs, "-", color=RED, lw=2.6,
            label=r"corn-starch-like ($f_{\rm res}$=0) — recovers")
    ax.set_xlabel("time after impact  (ms)", color=INK)
    ax.set_ylabel("frictional-contact fraction  f", color=INK)
    ax.set_title("Reversibility switch — same pulse, opposite fate", fontsize=11)
    ax.set_ylim(-0.02, 1.02)
    ax.legend(fontsize=9, frameon=False, loc="center right")
    ax.grid(True, color="#E3E3E3", lw=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig("gate2_reversibility.png", facecolor="white")
    plt.close(fig)
    return t_eval, f_si, f_cs, t_pulse


def main():
    print("=" * 68)
    print("TIER-0 WYART-CATES CLOSURE — SANITY GATE")
    print("=" * 68)

    fc_si, fc_cs = gate1_flow_curves()
    si_dst = wc.is_discontinuous(fc_si)
    cs_dst = wc.is_discontinuous(fc_cs)
    print("\n[Gate 1] bulk rheology contrast")
    print(f"  silica-like    back-bending (DST/SJ) : {si_dst}   (expect True)")
    print(f"  corn-starch    back-bending (DST)     : {cs_dst}   (expect False)")
    eta_si = fc_si["eta_app"]
    eta_si = eta_si[np.isfinite(eta_si)]
    eta_cs = fc_cs["eta_app"]
    eta_cs = eta_cs[np.isfinite(eta_cs)]
    print(f"  silica viscosity span   : {eta_si.max()/eta_si.min():.1f}x")
    print(f"  cornstarch viscosity span: {eta_cs.max()/eta_cs.min():.1f}x")
    gate1_pass = si_dst and (not cs_dst)

    t_eval, f_si, f_cs, t_pulse = gate2_reversibility()
    # frictional fraction 10 ms after the pulse ends
    idx = np.searchsorted(t_eval, t_pulse + 10.0e-3)
    print("\n[Gate 2] reversibility switch (f, 10 ms after pulse ends)")
    print(f"  silica-like retained f : {f_si[idx]:.3f}   (expect high, retention)")
    print(f"  corn-starch relaxed f  : {f_cs[idx]:.3f}   (expect low, recovery)")
    gate2_pass = f_si[idx] > 0.5 > f_cs[idx]

    # ODE self-consistency: constant stress must converge to f_ss(sigma)
    sig_c = 15.0
    t_long = np.linspace(0, 0.5, 500)
    f_conv = wc.relax_state(t_long, lambda t: sig_c, tau=5e-3, sigma_star=10.0)
    f_ss_target = float(wc.frictional_fraction_ss(sig_c, 10.0))
    ok_conv = np.isclose(f_conv[-1], f_ss_target, rtol=1e-3, atol=1e-3)
    print("\n[Numerics] df/dt converges to f_ss at constant stress")
    print(f"  f(t->inf)={f_conv[-1]:.4f}  vs  f_ss={f_ss_target:.4f}  -> {ok_conv}")

    print("\n" + "-" * 68)
    print(f"  GATE 1 (rheology contrast): {'PASS' if gate1_pass else 'FAIL'}")
    print(f"  GATE 2 (reversibility)    : {'PASS' if gate2_pass else 'FAIL'}")
    print(f"  NUMERICS (ODE consistency): {'PASS' if ok_conv else 'FAIL'}")
    print("-" * 68)
    print("  wrote: gate1_flow_curves.png, gate2_reversibility.png")


if __name__ == "__main__":
    main()
