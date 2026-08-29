"""
Glymphatic clearance: an advection-vs-diffusion regime map (Peclet analysis).

This is NOT a simulation -- it is a transparent dimensional-analysis map of the
field's central controversy (does brain solute clearance go by bulk flow or
diffusion?), built from literature-cited numbers. The Peclet number
    Pe = U * L / D
compares advective to diffusive transport over a length L. Pe >> 1 -> advection
dominates; Pe << 1 -> diffusion dominates; Pe ~ 1 -> crossover.

Current synthesis in the field (Kelley & Thomas ARFM 2023; Holter PNAS 2017):
advection dominates in the fast perivascular spaces, diffusion dominates deep in
the low-permeability parenchyma. This map makes that quantitative and shows the
crossover length for solutes of different size (small molecule vs amyloid-beta).

Literature-cited inputs (verify against sources before formal use):
  - Perivascular CSF speed U ~ 10-20 um/s   (Mestre et al., Nat. Commun. 2018)
  - Parenchymal (Darcy) speed  U ~ 0.1-1 um/s or less (Holter PNAS 2017; low k)
  - D (small molecule, e.g. mannitol) ~ 5e-10 m^2/s
  - D (amyloid-beta monomer)          ~ 1e-10 m^2/s
  - tissue tortuosity lambda ~ 1.6 -> D_eff = D / lambda^2 ~ D / 2.6 (Nicholson)
Oscillatory (Taylor-Aris) dispersion can raise the effective D in pulsating
channels by tens-of-% to several-fold (Asgari 2016) -- shown as an uncertainty
band, not a precise value.
"""
from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- literature-cited constants (SI) ----
U_PVS = 15e-6          # m/s, perivascular CSF speed (Mestre 2018)
U_PARENCHYMA = 0.5e-6  # m/s, parenchymal Darcy speed (Holter 2017; illustrative low-k)
D_SMALL = 5e-10        # m^2/s, small-molecule diffusivity
D_ABETA = 1e-10        # m^2/s, amyloid-beta monomer
TORTUOSITY = 1.6       # -> D_eff = D / tortuosity^2 in tissue


def peclet(U, L, D):
    return U * L / D


def crossover_length(U, D):
    """L where Pe = 1."""
    return D / U


def main():
    print("=" * 62)
    print("GLYMPHATIC CLEARANCE — advection-vs-diffusion regime map")
    print("=" * 62)

    L = np.logspace(-6, -2, 400)  # 1 um .. 10 mm transport length

    cases = [
        ("perivascular, small molecule", U_PVS, D_SMALL, "#4E2A84", "-"),
        ("perivascular, amyloid-beta", U_PVS, D_ABETA, "#C0392B", "-"),
        ("parenchyma, small molecule", U_PARENCHYMA, D_SMALL / TORTUOSITY**2, "#4E2A84", "--"),
        ("parenchyma, amyloid-beta", U_PARENCHYMA, D_ABETA / TORTUOSITY**2, "#C0392B", "--"),
    ]

    fig, ax = plt.subplots(figsize=(8.0, 5.0), dpi=140)
    print("\n  crossover length L* (Pe=1):")
    for name, U, D, color, ls in cases:
        Pe = peclet(U, L, D)
        ax.loglog(L * 1e6, Pe, ls, color=color, lw=2.4, label=name)
        Lstar = crossover_length(U, D) * 1e6  # um
        print(f"    {name:<32}: L* = {Lstar:8.1f} um")

    ax.axhline(1.0, color="#555555", lw=1.6, ls=":")
    ax.text(1.3, 1.4, "Pe = 1  (advection = diffusion)", fontsize=9, color="#555555")
    ax.text(1.3, 8e2, "advection dominates", fontsize=9.5, color="#333333", style="italic")
    ax.text(1.3, 3e-3, "diffusion dominates", fontsize=9.5, color="#333333", style="italic")
    # annotate physiological length scales
    for Lann, lab in [(20, "ECS gap→cell"), (250, "PVS segment"), (2000, "cortical depth")]:
        ax.axvline(Lann, color="#CFCFCF", lw=1.0, zorder=0)
        ax.text(Lann * 1.05, 1e-4, lab, rotation=90, fontsize=8, color="#999999", va="bottom")

    ax.set_xlabel(r"transport length  $L$  ($\mu$m)")
    ax.set_ylabel("Péclet number  Pe = UL/D")
    ax.set_title("Where brain clearance is advective vs diffusive", fontsize=11)
    ax.legend(fontsize=8.5, frameon=False, loc="lower right")
    ax.grid(True, which="both", color="#EEEEEE", lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig("glymphatic_peclet.png", facecolor="white")
    plt.close(fig)

    # -------- validation (sanity, not simulation) --------
    # perivascular small-molecule should be advective at PVS-segment scale (~250 um)
    pe_pvs = peclet(U_PVS, 250e-6, D_SMALL)
    # deep parenchyma should be diffusive at ECS scale (~20 um)
    pe_par = peclet(U_PARENCHYMA, 20e-6, D_SMALL / TORTUOSITY**2)
    print(f"\n  Pe (perivascular, 250 um) = {pe_pvs:6.1f}   -> {'advective' if pe_pvs>1 else 'diffusive'}")
    print(f"  Pe (parenchyma,   20 um)  = {pe_par:6.3f}   -> {'advective' if pe_par>1 else 'diffusive'}")
    ok = pe_pvs > 1 and pe_par < 1
    print("\n  matches field synthesis (PVS advective, deep tissue diffusive):", ok)
    print("  wrote: glymphatic_peclet.png")
    print("-" * 62)
    print(f"  REGIME MAP VALIDATION: {'PASS' if ok else 'FAIL'}")
    print("-" * 62)


if __name__ == "__main__":
    main()
