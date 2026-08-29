"""
Blood bridge demo: one structural-closure framework, opposite responses.

Figure 1: apparent viscosity vs shear rate — blood SHEAR-THINS (RBC rouleaux
          break up) while a silica suspension SHEAR-THICKENS (frictional
          contacts form). Both are first-order structural state-variable
          closures; only the sign of structure->viscosity differs.
Figure 2: blood thixotropic hysteresis loop under an up/down shear ramp — the
          biological analogue of the silica retention hysteresis in the
          drop-impact work.

Prints a PASS/FAIL validation report.
"""
from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import blood_rheology as blood

PURPLE, RED, INK = "#4E2A84", "#C0392B", "#2B2B2B"


def silica_cst_curve():
    """Minimal inline Wyart-Cates CST (continuous thickening) curve for contrast.
    phi < phi_m so it thickens but keeps flowing (clean comparison to blood).
    """
    phi, phi_m, phi0, sig_star, eta_s = 0.54, 0.56, 0.65, 10.0, 1.0e-3
    sig = np.logspace(-1, 4, 400)
    f = np.where(sig > 0, np.exp(-sig_star / sig), 0.0)
    phi_J = (1 - f) * phi0 + f * phi_m
    eta_r = (1 - phi / phi_J) ** (-2.0)
    gdot = sig / (eta_r * eta_s)
    eta_app = sig / gdot                       # = eta_r * eta_s
    return gdot, eta_app


def main():
    print("=" * 66)
    print("BLOOD BRIDGE — one closure framework, opposite response")
    print("=" * 66)

    # --- blood steady curve ---
    gdot_b = np.logspace(-2, 3, 400)
    fc = blood.steady_flow_curve(gdot_b)
    eta_b = fc["eta_app"]
    lam = fc["lambda"]

    # --- silica contrast ---
    gdot_s, eta_s = silica_cst_curve()

    fig, ax = plt.subplots(figsize=(7.4, 4.6), dpi=140)
    ax.loglog(gdot_b, eta_b * 1e3, "-", color=RED, lw=2.6,
              label="blood — shear-THINS (rouleaux break up)")
    ax.loglog(gdot_s, eta_s * 1e3, "-", color=PURPLE, lw=2.6,
              label="silica — shear-THICKENS (frictional contacts)")
    ax.set_xlabel(r"shear rate  $\dot\gamma$  (s$^{-1}$)", color=INK)
    ax.set_ylabel(r"apparent viscosity  (mPa$\cdot$s)", color=INK)
    ax.set_title("Same state-variable closure, opposite sign of structure coupling",
                 fontsize=11)
    ax.legend(fontsize=9.5, frameon=False)
    ax.grid(True, which="both", color="#E7E7E7", lw=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig("fig1_blood_vs_silica.png", facecolor="white")
    plt.close(fig)

    # --- thixotropic hysteresis loop ---
    loop = blood.thixotropic_loop()
    fig2, ax2 = plt.subplots(figsize=(6.6, 4.6), dpi=140)
    ax2.plot(loop["g_up"], loop["sigma_up"] * 1e3, "-", color=RED, lw=2.4,
             label="increasing shear (structure breaking)")
    ax2.plot(loop["g_dn"], loop["sigma_dn"] * 1e3, "--", color=PURPLE, lw=2.4,
             label="decreasing shear (structure not yet rebuilt)")
    ax2.set_xlabel(r"shear rate  $\dot\gamma$  (s$^{-1}$)", color=INK)
    ax2.set_ylabel(r"shear stress  (mPa)", color=INK)
    ax2.set_title("Blood thixotropic hysteresis loop", fontsize=11)
    ax2.legend(fontsize=9.5, frameon=False)
    ax2.grid(True, color="#E7E7E7", lw=0.7)
    ax2.set_axisbelow(True)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    fig2.tight_layout()
    fig2.savefig("fig2_thixotropic_loop.png", facecolor="white")
    plt.close(fig2)

    # -------- validation --------
    thinning = np.all(np.diff(eta_b) < 1e-12)          # monotone decreasing
    ratio = eta_b.max() / eta_b.min()
    lam_ok = np.all((lam >= 0) & (lam <= 1)) and np.all(np.diff(lam) <= 1e-12)
    thickening = eta_s[-1] > eta_s[np.argmin(eta_s)]   # rises at high shear
    # thixotropy: down-branch stress below up-branch at matched rates -> loop area
    g_common = np.linspace(5, 90, 50)
    su = np.interp(g_common, loop["g_up"], loop["sigma_up"])
    sd = np.interp(g_common, loop["g_dn"][::-1], loop["sigma_dn"][::-1])
    d = su - sd  # up-branch minus down-branch stress at matched shear rates
    loop_area = float(np.sum((d[:-1] + d[1:]) / 2.0 * np.diff(g_common)))  # trapezoid
    hysteresis = loop_area > 0

    print(f"\n  blood shear-thinning (eta decreasing)   : {thinning}   (expect True)")
    print(f"  blood eta0/eta_inf ratio                : {ratio:.0f}x  (physiological ~ tens)")
    print(f"  aggregation lambda in [0,1], decreasing : {lam_ok}   (expect True)")
    print(f"  silica shear-thickening (contrast)      : {thickening}   (expect True)")
    print(f"  thixotropic loop area > 0 (hysteresis)  : {hysteresis}   (expect True)")
    print("\n  wrote: fig1_blood_vs_silica.png, fig2_thixotropic_loop.png")
    allpass = thinning and lam_ok and thickening and hysteresis and (ratio > 10)
    print("-" * 66)
    print(f"  VALIDATION: {'PASS' if allpass else 'FAIL'}")
    print("-" * 66)


if __name__ == "__main__":
    main()
