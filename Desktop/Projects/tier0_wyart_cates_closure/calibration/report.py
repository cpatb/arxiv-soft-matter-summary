"""The tier-0 calibration verdict: what transfers into a solver, and what does not.

Reads ``out/calibration.json`` (written by ``calibrate.py``), runs the impact
propagation, prints the verdict tables and writes the figure.

Every table is reported twice: over every usable curve, and over the subset the
*source lab* calls steady (``calibration/steady.py``). The two differ enough
that quoting only one of them would be a choice about the answer.

Run:  python calibration/calibrate.py && python calibration/report.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from propagate_impact import (  # noqa: E402
    K_ADDED_MASS,
    PHI_TOLERANCE,
    RHO_SUSPENSION,
    U0_ONSET,
    band_consequence,
    flip_window,
    impact_stress,
    saturation,
)

# The packing the drop-impact programme actually works at, from its own Status
# card: phi = 0.44 silica, which is the salt axis's pivot.
PROGRAMME_PHI = 0.44

OUT = HERE / "out"
SUBSETS = (("", "all usable curves"), ("_steady", "the source lab's own steady subset"))


def load() -> dict:
    path = OUT / "calibration.json"
    if not path.is_file():
        raise SystemExit(f"{path} not found; run calibration/calibrate.py first")
    return json.loads(path.read_text(encoding="utf-8"))


def usable_curves(d: dict) -> list[dict]:
    return [c for c in d["curves"] if not c["at_bound_free"]]


def steady_curves(d: dict) -> list[dict]:
    labels = set(d.get("steady_labels", ()))
    return [c for c in usable_curves(d) if c["label"] in labels]


def saturation_boundary(phi_0=0.61, phi_m=0.56, band=0.03, beta=1.0,
                        lo=1e-3, hi=1e7, tol=1e-3, gauge="inertial") -> float:
    """Smallest sigma* at which the activation half stops being saturated.

    Bisection on a monotone predicate: the residual ``(1-f)(phi_0 - phi_m)``
    grows with sigma*, so there is exactly one crossing. Returned in Pa.
    """
    def not_sat(ss):
        return not saturation(ss, beta, phi_0, phi_m, band, gauge=gauge).saturated

    if not_sat(lo):
        return lo
    if not not_sat(hi):
        return float("inf")
    while hi / lo > 1 + tol:
        mid = float(np.sqrt(lo * hi))
        if not_sat(mid):
            hi = mid
        else:
            lo = mid
    return float(np.sqrt(lo * hi))


def _identity_verdict(tm: dict) -> str:
    inside = (tm["slope_ci"][0] <= 1 <= tm["slope_ci"][1]
              and tm["intercept_ci"][0] <= 0 <= tm["intercept_ci"][1])
    return "contains the identity" if inside else "EXCLUDES the identity"


def print_transfer(data: dict) -> None:
    print("1. Does one flow curve fix the volume fraction the solver advects?")
    print("   Bootstrap over curves. The identity is slope 1, intercept 0.\n")
    for suffix, header in SUBSETS:
        print(f"   --- {header} ---")
        for d in data.values():
            n = d.get("n_steady") if suffix else d.get("n_usable")
            print(f"   {d['material']}  ({n} curves)")
            tm = d.get(f"transfer_free{suffix}")
            if tm:
                print(f"      free fit : slope {tm['slope']:.3f} "
                      f"[{tm['slope_ci'][0]:.3f}, {tm['slope_ci'][1]:.3f}]   "
                      f"intercept {tm['intercept']:+.4f} "
                      f"[{tm['intercept_ci'][0]:+.4f}, {tm['intercept_ci'][1]:+.4f}]"
                      f"   -> {_identity_verdict(tm)}")
            else:
                print("      free fit : refused, fewer than 3 curves")
            tp = d.get(f"transfer_plateau{suffix}")
            if tp:
                print(f"      plateau  : slope {tp['slope']:.3f} "
                      f"[{tp['slope_ci'][0]:.3f}, {tp['slope_ci'][1]:.3f}]   "
                      f"intercept {tp['intercept']:+.4f} "
                      f"[{tp['intercept_ci'][0]:+.4f}, {tp['intercept_ci'][1]:+.4f}]"
                      f"   -> {_identity_verdict(tp)}")
            else:
                n_adm = d.get(f"n_plateau_admitted{suffix}", 0)
                print(f"      plateau  : refused, {n_adm} curve(s) non-Brownian at the "
                      f"lowest stress sampled")
        print()


def print_saturation(data: dict) -> None:
    print("2. Is the half of the model a flow curve DOES measure in play during impact?\n")
    print("   Reported on BOTH stress gauges, because they disagree and the choice is")
    print("   not a detail. rho k U0^2 is a TOTAL stress; under undrained loading the")
    print("   pore fluid carries the isotropic part, so the deviatoric (wall shear)")
    print("   gauge is the one that reaches the particle network. The lab's own")
    print("   impact_stress/STRESS_AND_DILATANCY_BRIEF.md establishes this for this")
    print("   geometry and measures the ratio at 5 to 10x.\n")
    for gauge in ("inertial", "deviatoric"):
        print(f"   --- {gauge} gauge ---")
        for d in data.values():
            rows = steady_curves(d) or usable_curves(d)
            band = d.get("verdict_band_free_steady") or d.get("verdict_band_free")
            if not rows or not band:
                continue
            ss = float(np.median([c["sigma_star_free"] for c in rows]))
            pm = float(np.median([c["phi_m_free"] for c in rows]))
            s = saturation(ss, d["constants"]["beta"], d["constants"]["phi_0"], pm,
                           band["width"], gauge=gauge)
            print(f"   {d['material'][:44]:<45s} {s}")
        # Computed from each material's OWN band. The earlier version called
        # this bare, with phi_0 - phi_m = 0.05 and a 0.03 band matching none of
        # the three, and then quoted the resulting particle size as a result.
        print("   boundary sigma* above which the activation half stops being saturated,")
        print("   from each material's own band:")
        for d in data.values():
            band = d.get("verdict_band_free_steady") or d.get("verdict_band_free")
            rows = steady_curves(d) or usable_curves(d)
            if not band or not rows:
                continue
            pm = float(np.median([c["phi_m_free"] for c in rows]))
            b = saturation_boundary(phi_0=d["constants"]["phi_0"], phi_m=pm,
                                    band=band["width"], beta=d["constants"]["beta"],
                                    gauge=gauge)
            print(f"      {d['material'][:40]:<41s} {b:12,.1f} Pa")
        print()


def print_bands(data: dict) -> None:
    print("3. What the phi_m band does to the impact prediction\n")
    print(f"   {'material':<32s} {'subset':>8s} {'protocol':>8s} "
          f"{'phi_m band':>18s} {'n':>3s} {'width':>8s} {'x tol':>6s}")
    for d in data.values():
        first = True
        for suffix, _ in SUBSETS:
            for proto in ("free", "pinned"):
                band = d.get(f"verdict_band_{proto}{suffix}")
                if not band:
                    continue
                _, _, width, xtol = flip_window(band["lo"], band["hi"])
                name = d["material"][:31] if first else ""
                first = False
                print(f"   {name:<32s} {('steady' if suffix else 'all'):>8s} "
                      f"{proto:>8s} [{band['lo']:.4f}, {band['hi']:.4f}] "
                      f"{band['n']:3d} {width:8.4f} {xtol:6.0f}")
    print(f"\n   'x tol' is the band width divided by the {PHI_TOLERANCE} volume-fraction "
          f"tolerance the retake campaign holds itself to.\n")


def print_flip(data: dict) -> None:
    print("4. Does the phi_m band straddle a packing anyone actually works at?\n")
    print("   Evaluating the model at the band's own midpoint makes the verdict flip by")
    print("   construction and tests nothing. So ask first whether the drop-impact")
    print("   programme's own packing lands inside any of these bands.\n")
    for d in data.values():
        for tag in ("", "_steady"):
            band = d.get(f"verdict_band_free{tag}")
            if not band:
                continue
            inside = band["lo"] < PROGRAMME_PHI < band["hi"]
            print(f"   {d['material'][:40]:<41s} {('steady' if tag else 'all'):>7s} "
                  f"[{band['lo']:.4f}, {band['hi']:.4f}]  phi = {PROGRAMME_PHI} is "
                  f"{'INSIDE' if inside else 'outside'}")
    print(f"\n   phi = {PROGRAMME_PHI} sits below every band. For this programme's own")
    print("   material the calibration does not straddle the verdict; the band bites for")
    print("   anyone working between 0.53 and 0.59 (PMMA) or 0.45 and 0.51 (corn starch).")
    print("\n   Band width against the null of a random subset of the same size:")
    for d in data.values():
        null = d.get("band_null_random4")
        band = d.get("verdict_band_free_steady")
        if not null or not band:
            continue
        print(f"      {d['material'][:40]:<41s} steady {band['width']:.4f}   "
              f"random-{null['k']}: mean {null['mean']:.4f}, min {null['min']:.4f}, "
              f"max {null['max']:.4f}")
    print("\n   The steady subset is the LOWEST-phi block, and phi_m rises with phi, so it")
    print("   is the tightest subset available whatever steadiness has to do with it --")
    print("   for two of the three it IS the minimum over all subsets of that size. The")
    print("   defensible statement is the null's floor: however four curves of one")
    print("   material are chosen, the band is at least ~0.012 and typically ~0.048,")
    print(f"   i.e. 6x to 24x the {PHI_TOLERANCE} tolerance phi is set to.")
    print("\n   The model at each band's midpoint, for completeness:\n")
    for d in data.values():
        band = d.get("verdict_band_free_steady") or d.get("verdict_band_free")
        rows = steady_curves(d) or usable_curves(d)
        if not band or not rows or band["width"] <= 0:
            continue
        ss = float(np.median([c["sigma_star_free"] for c in rows]))
        phi_mid = 0.5 * (band["lo"] + band["hi"])
        bc = band_consequence(phi_mid, band["lo"], band["hi"],
                              d["constants"]["phi_0"], ss, d["constants"]["beta"],
                              n=d["constants"]["n"])
        print(f"   {d['material']} at phi = {phi_mid:.4f}")
        print(f"      phi_m = {band['lo']:.4f} -> eta_r {bc.eta_r_at_lo:>10.4g}, "
              f"jams {bc.jams_at_lo}")
        print(f"      phi_m = {band['hi']:.4f} -> eta_r {bc.eta_r_at_hi:>10.4g}, "
              f"jams {bc.jams_at_hi}")
        print(f"      verdict flips: {bc.verdict_flips}   "
              f"over a window {band['width']:.4f} wide "
              f"({band['width'] / PHI_TOLERANCE:.0f}x the phi tolerance)\n")


def main() -> int:
    data = load()
    lo_s, hi_s = impact_stress()
    print("=" * 100)
    print("TIER-0 CALIBRATION VERDICT")
    print("=" * 100)
    print(f"\nImpact stress scale, Sigma = rho k U0^2 with rho = {RHO_SUSPENSION:.0f} kg/m^3, "
          f"U0 = {U0_ONSET} m/s, k = {K_ADDED_MASS[0]:g} to {K_ADDED_MASS[1]:g}:")
    print(f"    {lo_s:,.0f} to {hi_s:,.0f} Pa\n")
    print_transfer(data)
    print_saturation(data)
    print_bands(data)
    print_flip(data)
    figure(data)
    return 0


def figure(data: dict) -> None:
    """Two panels: the transfer residual, and the phi_m band against the flip."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 4.0))
    colours = {"guy_404": "#4C72B0", "guy_3770": "#C44E52",
               "richards_cornstarch": "#55A868"}
    short = {"guy_404": "PMMA 404 nm", "guy_3770": "PMMA 3770 nm",
             "richards_cornstarch": "corn starch"}

    # Panel 1: the residual, inferred minus made. Plotted as a residual rather
    # than as inferred-against-made because the house rule keeps zero on a
    # linear axis, and on a 0-to-0.65 axis a 0.04 offset is invisible. Zero here
    # is where a flow curve would sit if it reported the packing of the sample.
    for key, d in data.items():
        rows = usable_curves(d)
        if not rows:
            continue
        steady = {c["label"] for c in steady_curves(d)}
        x = np.array([c["phi_known"] for c in rows])
        y = np.array([c["phi_free"] for c in rows]) - x
        filled = np.array([c["label"] in steady for c in rows])
        ax1.plot(x[filled], y[filled], "o", color=colours[key], ms=5, label=short[key])
        if (~filled).any():
            ax1.plot(x[~filled], y[~filled], "o", color=colours[key], ms=5, mfc="none")
        p = np.array([c["phi_plateau"] for c in rows])
        keep = np.isfinite(p)
        if keep.any():
            ax1.plot(x[keep], p[keep] - x[keep], "^", color=colours[key], ms=5,
                     mfc="none", alpha=0.6)
        tm = d.get("transfer_free")
        if tm:
            xs = np.linspace(x.min(), x.max(), 20)
            ax1.plot(xs, tm["slope"] * xs + tm["intercept"] - xs, "-",
                     color=colours[key], lw=1.0, alpha=0.7)
    ax1.axhline(0.0, color="k", ls="--", lw=0.8)
    ax1.set_xlim(0.0, 0.65)
    ax1.set_xlabel(r"$\phi$ the sample was made at")
    ax1.set_ylabel(r"$\phi$ inferred $-$ $\phi$ made")
    ax1.legend(frameon=False, fontsize=8, loc="lower right")

    # Panel 2: the phi_m band per material, all curves against the steady subset
    ys, labels = [], []
    row = 0
    for key, d in data.items():
        for suffix, _ in SUBSETS:
            band = d.get(f"verdict_band_free{suffix}")
            if not band:
                continue
            ax2.plot([band["lo"], band["hi"]], [row, row], "-", lw=6,
                     color=colours[key], alpha=0.45 if suffix else 0.9,
                     solid_capstyle="butt")
            if d.get("phi_m_published") is not None:
                ax2.plot([d["phi_m_published"]], [row], "k|", ms=12, mew=1.4)
            ys.append(row)
            labels.append(f"{short[key]}, {'steady' if suffix else 'all'}")
            row += 1
    ax2.set_yticks(ys)
    ax2.set_yticklabels(labels, fontsize=8)
    ax2.set_xlim(0.40, 0.65)
    ax2.set_xlabel(r"$\phi_m$ returned by the material's own curves")
    ax2.invert_yaxis()

    for ax in (ax1, ax2):
        ax.tick_params(labelsize=8)
    fig.tight_layout()
    dest = OUT / "calibration_transfer.png"
    fig.savefig(dest, dpi=200)
    plt.close(fig)
    print(f"wrote {dest}")


if __name__ == "__main__":
    raise SystemExit(main())
