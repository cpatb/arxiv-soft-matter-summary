"""What phi_m would each drop need? Run the inversion and report the drift.

Run:  python run_inversion.py
Writes out/inversion.json and out/inversion.png.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from clocks import material  # noqa: E402
from data import is_arrested, load_drops  # noqa: E402
from invert import (  # noqa: E402
    SPREAD_EXPONENTS,
    drift,
    invert,
    prefactor_phi_m,
)

OUT = HERE / "out"
# The packing tolerance a suspension is actually made to, and the yardstick a
# "constant" phi_m has to hold still against.
PHI_TOLERANCE = 0.002
# beta_p's own noise floor, measured two ways in this session: 6.7% within-cell
# replicate CV over 12 cells with n >= 3, and 20% between the independent side
# and under-view channels on the 21 drops that have both. The tighter one is
# used, so the "phi_m must hold still" bar is the demanding one.
BETA_NOISE_CV = 0.067


def main() -> int:
    drops, prov = load_drops()
    m44 = material(0.44)
    OUT.mkdir(exist_ok=True)

    print("=" * 100)
    print("WHAT phi_m WOULD THE FOOTPRINTS NEED?")
    print("=" * 100)
    print(f"\nThe closure as fitted: phi0 = {m44.phi0}, phiM = {m44.phiM}, "
          f"eta_s = {m44.eta_s} Pa.s, sigma* = {m44.sigma_star:.2f} Pa")
    print(f"Its fully frictional viscosity at phi = 0.44: "
          f"{m44.eta_frictional:.4f} Pa.s -- the most the model can ever give.")
    print(f"\nbeta_p noise floor used as the 'holds still' bar: {BETA_NOISE_CV:.1%} "
          f"(within-cell replicate CV)")

    payload = {"provenance": prov, "beta_noise_cv": BETA_NOISE_CV, "exponents": {}}

    # The inversion assumes viscous ARREST, so it is only meaningful on drops
    # that arrested. The 0 and 10 mM cells straddle the transition; applying an
    # arrest model to a drop that did not arrest is what produced the large
    # apparent velocity drift in the unrestricted version of this analysis.
    arrested = [d for d in drops if is_arrested(d)]
    print(f"Restricted to arrested drops (beta_p < 0.9 beta_eq): "
          f"{len(arrested)} of {len(drops)}")

    for p in SPREAD_EXPONENTS:
        invs = [invert(d, p) for d in arrested]
        ok = [i for i in invs if np.isfinite(i.phi_m_required)]
        tag = f"p={p:g}"
        print("\n" + "-" * 100)
        print(f"SPREAD EXPONENT {tag}   (beta ~ Re^(1/{p:g}))")
        print("-" * 100)

        eta = np.array([i.eta_implied for i in invs])
        gap = np.array([i.viscosity_gap for i in invs])
        print(f"  implied viscosity: {eta.min():.4g} to {eta.max():.4g} Pa.s")
        print(f"  the closure's ceiling (f = 1) at phi = 0.44: "
              f"{m44.eta_frictional:.4f} Pa.s")
        print(f"  gap, implied / ceiling: {gap.min():.2f}x to {gap.max():.0f}x; "
              f"{(gap > 1).sum()} of {len(gap)} drops are ABOVE the ceiling")

        pm = np.array([i.phi_m_required for i in ok])
        print(f"\n  required phi_m: {pm.min():.4f} to {pm.max():.4f}  "
              f"(span {pm.max() - pm.min():.4f} = "
              f"{(pm.max() - pm.min()) / PHI_TOLERANCE:.0f}x the phi tolerance)")
        print(f"  median {np.median(pm):.4f}, against the fitted {m44.phiM}")

        # The decisive test: does the required phi_m drift with U0 at fixed
        # packing AND fixed salt? phi_m is a material constant; it cannot.
        print("\n  DOES IT HOLD STILL? drift of required phi_m within each "
              "(phi, salt) cell:")
        cells = defaultdict(list)
        for i in ok:
            cells[(i.phi, i.salt_mM)].append(i)
        rows = []
        for (phi, salt), group in sorted(cells.items()):
            if len(group) < 4:
                continue
            u = [g.u0 for g in group]
            if max(u) / min(u) < 1.3:
                continue                      # no velocity lever in this cell
            v = [g.phi_m_required for g in group]
            dr = drift(v, u)
            if dr is None:
                continue
            dr |= {"phi": phi, "salt_mM": salt,
                   "u0_ratio": float(max(u) / min(u))}
            rows.append(dr)
            print(f"    phi {phi:.2f}, {salt:5.0f} mM, n={dr['n']:2d}, "
                  f"u0 x{dr['u0_ratio']:.2f}: rho {dr['rho']:+.3f} "
                  f"(p {dr['p']:.3f}), span {dr['span']:.4f} = "
                  f"{dr['span'] / PHI_TOLERANCE:5.1f}x tol, "
                  f"slope {dr['slope']:+.4f} per m/s")

        if rows:
            spans = np.array([r["span"] for r in rows])
            rhos = np.array([r["rho"] for r in rows])
            print(f"\n    {len(rows)} cells with a velocity lever. "
                  f"Required phi_m spans {spans.min():.4f} to {spans.max():.4f} "
                  f"within a cell,")
            print(f"    i.e. {spans.min() / PHI_TOLERANCE:.0f}x to "
                  f"{spans.max() / PHI_TOLERANCE:.0f}x the tolerance phi is set to, "
                  f"at FIXED phi and FIXED salt.")
            same_sign = int((rhos > 0).sum()), int((rhos < 0).sum())
            print(f"    sign of the drift: {same_sign[0]} positive, "
                  f"{same_sign[1]} negative out of {len(rhos)}")

        # And across salt at fixed phi, which the closure has no mechanism for
        # at all once f is saturated.
        print("\n  ACROSS SALT at phi = 0.44, matched on velocity where possible:")
        by_salt = defaultdict(list)
        for i in ok:
            if i.phi == 0.44:
                by_salt[i.salt_mM].append(i.phi_m_required)
        for salt in sorted(by_salt):
            v = np.array(by_salt[salt])
            print(f"    {salt:5.0f} mM: n={v.size:2d}  median phi_m required "
                  f"{np.median(v):.4f}  [{v.min():.4f}, {v.max():.4f}]")
        meds = [float(np.median(by_salt[s])) for s in sorted(by_salt)]
        if len(meds) >= 3:
            print(f"    span of the medians: {max(meds) - min(meds):.4f} = "
                  f"{(max(meds) - min(meds)) / PHI_TOLERANCE:.0f}x the tolerance")

        # How much of this survives the unknown O(1) prefactor in the
        # spreading law? Sweep it over four decades and report the median.
        print("\n  PREFACTOR SENSITIVITY  (eta = C rho U0 D0 / beta^p):")
        pref = {}
        for C in (0.01, 0.1, 1.0, 10.0, 100.0):
            v = [prefactor_phi_m(d, p, C) for d in arrested]
            v = np.array([x for x in v if np.isfinite(x)])
            pref[f"C={C:g}"] = float(np.median(v))
            print(f"    C = {C:6.2f}: median required phi_m = {np.median(v):.4f}")
        span = max(pref.values()) - min(pref.values())
        print(f"    over 4 decades of C the median moves {span:.4f}, against a "
              f"{m44.phiM - np.median(pm):.4f} gap to the fitted value.")

        payload["exponents"][tag] = {
            "prefactor_sensitivity": pref,
            "eta_implied_min": float(eta.min()), "eta_implied_max": float(eta.max()),
            "eta_model_ceiling": float(m44.eta_frictional),
            "n_above_ceiling": int((gap > 1).sum()), "n_total": len(gap),
            "phi_m_required_min": float(pm.min()),
            "phi_m_required_max": float(pm.max()),
            "phi_m_required_median": float(np.median(pm)),
            "within_cell_drift": rows,
            "salt_medians": {f"{s:g}": float(np.median(by_salt[s]))
                             for s in sorted(by_salt)},
            "curves": [asdict(i) for i in invs],
        }

    (OUT / "inversion.json").write_text(json.dumps(payload, indent=2),
                                        encoding="utf-8")
    print(f"\nwrote {OUT / 'inversion.json'}")
    figure(drops, payload)
    return 0


def figure(drops, payload) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p = SPREAD_EXPONENTS[-1]
    invs = [invert(d, p) for d in drops]
    m44 = material(0.44)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 4.2))
    colours = {0.0: "#4C72B0", 10.0: "#55A868", 50.0: "#DD8452", 100.0: "#C44E52"}

    # Panel 1: the implied viscosity against the closure's own ceiling.
    for salt, c in colours.items():
        sel = [i for i in invs if i.salt_mM == salt and i.phi == 0.44]
        if not sel:
            continue
        ax1.plot([i.u0 for i in sel], [i.eta_implied for i in sel], "o",
                 ms=4.5, color=c, mfc="none", label=f"{salt:g} mM")
    ax1.axhline(m44.eta_frictional, color="k", ls="--", lw=1.0)
    ax1.axhline(m44.eta_0, color="k", ls=":", lw=1.0)
    ax1.set_yscale("log")
    ax1.set_xlabel(r"impact speed $U_0$  (m/s)")
    ax1.set_ylabel(r"viscosity the footprint implies  (Pa$\cdot$s)")
    ax1.legend(frameon=False, fontsize=8, title=r"$\phi = 0.44$",
               title_fontsize=8, loc="upper right")
    ax1.tick_params(labelsize=8)

    # Panel 2: the required phi_m, which is meant to be one number.
    for salt, c in colours.items():
        sel = [i for i in invs
               if i.salt_mM == salt and i.phi == 0.44
               and np.isfinite(i.phi_m_required)]
        if not sel:
            continue
        ax2.plot([i.u0 for i in sel], [i.phi_m_required for i in sel], "o",
                 ms=4.5, color=c, mfc="none")
    ax2.axhline(m44.phiM, color="k", ls="--", lw=1.0)
    ax2.axhline(0.44, color="0.5", ls=":", lw=1.0)
    ax2.set_xlabel(r"impact speed $U_0$  (m/s)")
    ax2.set_ylabel(r"$\phi_m$ the closure would need")
    ax2.tick_params(labelsize=8)

    fig.tight_layout()
    dest = OUT / "inversion.png"
    fig.savefig(dest, dpi=200)
    plt.close(fig)
    print(f"wrote {dest}")


if __name__ == "__main__":
    raise SystemExit(main())
