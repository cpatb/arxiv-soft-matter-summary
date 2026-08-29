"""
Thrombosis-as-jamming: the falsifiable signature.

Plots occlusion time vs wall shear rate and shows the CRITICAL DIVERGENCE at a
threshold shear rate (no occlusion below it; t_occ -> inf as gdot -> gdot_crit+).
This is the qualitative signature that distinguishes a jamming mechanism from the
incumbent deposition-kinetics thrombosis models (which give a smooth monotonic
speed-up with shear and NO critical point).

Prints a PASS/FAIL validation of that signature. Parameters are illustrative;
the falsifiable test on real data is: replot published occlusion-time-vs-shear
against (gdot - gdot_crit) and look for a power-law collapse.
"""
from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import thrombosis_jamming as tj

RED, PURPLE, GOLD, INK = "#C0392B", "#4E2A84", "#B6862C", "#2B2B2B"


def main():
    print("=" * 62)
    print("THROMBOSIS AS JAMMING — falsifiable critical-divergence test")
    print("=" * 62)

    gdot = np.logspace(3.0, 4.6, 700)          # ~1000 .. 40000 /s (spans the threshold)
    t_occ, gdot_crit = tj.occlusion_time(gdot)
    finite = np.isfinite(t_occ)

    print(f"\n  VWF activation shear rate gdot_c : {tj.GDOT_C:.0f} /s")
    print(f"  predicted occlusion threshold    : {gdot_crit:.0f} /s")
    occ = gdot[finite]
    print(f"  occlusion only above threshold   : {occ.min():.0f}-{occ.max():.0f} /s")
    if occ.size:
        print(f"  t_occ range (above threshold)    : {t_occ[finite].min():.2f}-{t_occ[finite].max():.1f} s")

    # -------- validation of the critical-divergence signature --------
    below = gdot < gdot_crit * 0.98
    no_occ_below = np.all(~np.isfinite(t_occ[below])) if below.any() else True
    # monotone decreasing above threshold
    ta = t_occ[finite]
    mono = np.all(np.diff(ta) <= 1e-9)
    # divergence: t_occ just above threshold is far larger than the typical value
    diverges = ta[0] > 5 * np.median(ta) if ta.size > 10 else False
    allpass = no_occ_below and mono and diverges

    print(f"\n  no occlusion below threshold     : {no_occ_below}   (expect True)")
    print(f"  t_occ monotone decreasing above  : {mono}   (expect True)")
    print(f"  t_occ diverges near threshold    : {diverges}   (expect True)")

    # figure
    fig, ax = plt.subplots(figsize=(7.6, 4.7), dpi=140)
    ax.axvspan(gdot.min(), gdot_crit, color="#BBBBBB", alpha=0.18)
    ax.text(gdot_crit * 0.75, 0.5, "no occlusion\n(flows)", ha="right", va="bottom",
            fontsize=9.5, color="#666666")
    ax.semilogy(gdot[finite], t_occ[finite], "-", color=RED, lw=2.8,
                label="occlusion time (jamming closure)")
    ax.axvline(gdot_crit, color=PURPLE, ls="--", lw=1.8,
               label=fr"$\dot\gamma_{{\rm crit}}\approx${gdot_crit:.0f} s$^{{-1}}$ (onset)")
    ax.axvline(tj.GDOT_C, color=GOLD, ls=":", lw=1.8,
               label=fr"$\dot\gamma_c$={tj.GDOT_C:.0f} s$^{{-1}}$ (VWF unfolding)")
    ax.set_xlabel(r"wall shear rate  $\dot\gamma$  (s$^{-1}$)", color=INK)
    ax.set_ylabel("occlusion time  (s)", color=INK)
    ax.set_title("SIPA as stress-activated jamming: a critical occlusion threshold",
                 fontsize=11)
    ax.set_xlim(gdot.min(), gdot.max())
    ax.legend(fontsize=9, frameon=False, loc="upper right")
    ax.grid(True, which="both", color="#EAEAEA", lw=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig("thrombosis_occlusion.png", facecolor="white")
    plt.close(fig)

    print("\n  wrote: thrombosis_occlusion.png")
    print("-" * 62)
    print(f"  SIGNATURE VALIDATION: {'PASS' if allpass else 'FAIL'}")
    print("  (Distinguishes jamming from smooth deposition-kinetics models.)")
    print("-" * 62)


if __name__ == "__main__":
    main()
