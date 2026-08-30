"""Run the elimination audit and write the table, the JSON and the figure.

Run:  python run_audit.py
"""
from __future__ import annotations

import json
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from audit import audit_clock  # noqa: E402
from clocks import (  # noqa: E402
    CLOCKS,
    K_ADDED_MASS_LO,
    frictional_fraction_at_impact,
    jamming_deficit,
    material,
)
from data import ARREST_RATIO, load_drops  # noqa: E402

OUT = HERE / "out"


def model_viscosity_span(drops) -> dict:
    """The closure's own predicted impact viscosity, split by what drives it.

    Evaluating the model the drop-impact repo actually calibrated at the stress
    each impact actually delivers gives a large spread -- but almost all of it
    comes from ``phi``, which the campaign varied deliberately, and almost none
    from the impact stress, which is the part the closure's stress activation is
    supposed to supply.

    So the span is reported three ways: over the whole dataset, over the
    velocity axis at fixed ``phi`` (the largest single-``phi`` group), and over
    the packing axis at fixed velocity. Quoting only the first would credit the
    stress clock with variation that belongs to the packing.
    """
    def eta(d):
        return float(material(d.phi).eta_of_f(
            frictional_fraction_at_impact(d, k=K_ADDED_MASS_LO)))

    etas = np.array([eta(d) for d in drops])
    fs = np.array([frictional_fraction_at_impact(d, k=K_ADDED_MASS_LO) for d in drops])

    phis = [d.phi for d in drops]
    main_phi = max(set(phis), key=phis.count)
    at_phi = [d for d in drops if d.phi == main_phi]
    e_phi = np.array([eta(d) for d in at_phi])
    u_phi = np.array([d.u0 for d in at_phi])
    b_phi = np.array([d.beta_p for d in at_phi])

    return {
        "eta_min": float(etas.min()), "eta_max": float(etas.max()),
        "eta_ratio": float(etas.max() / etas.min()),
        "f_min": float(fs.min()), "f_max": float(fs.max()),
        "fixed_phi": float(main_phi),
        "n_fixed_phi": len(at_phi),
        "u0_ratio_fixed_phi": float(u_phi.max() / u_phi.min()),
        "eta_ratio_fixed_phi": float(e_phi.max() / e_phi.min()),
        "beta_ratio_fixed_phi": float(b_phi.max() / b_phi.min()),
    }


def main() -> int:
    drops, prov = load_drops()
    beta = np.array([d.beta_p for d in drops])
    OUT.mkdir(exist_ok=True)

    print("=" * 104)
    print("WHAT ARRESTS AN IMPACTING SUSPENSION DROP?  An elimination audit")
    print("=" * 104)
    print(f"\nData: {prov['n_rows']} rows -> {prov['n_used']} used "
          f"({prov['n_bounded_beta_p_dropped']} dropped for a bounded beta_p: "
          f"{prov['bounded_kinds']}, {prov['n_non_finite_dropped']} non-finite)")
    print(f"  phi levels  {prov['phi_levels']}")
    print(f"  salt levels {prov['salt_levels']} mM")
    print(f"  u0  {prov['u0_range'][0]:.2f} to {prov['u0_range'][1]:.2f} m/s "
          f"({prov['u0_range'][1] / prov['u0_range'][0]:.2f}x)")
    print(f"  d0  {prov['d0_range_mm'][0]:.2f} to {prov['d0_range_mm'][1]:.2f} mm")
    print(f"  RESPONSE beta_p {prov['beta_p_range'][0]:.3f} to "
          f"{prov['beta_p_range'][1]:.3f}  "
          f"({prov['beta_p_range'][1] / prov['beta_p_range'][0]:.2f}x, "
          f"log-range {prov['beta_p_log_range']:.3f})")
    print(f"  arrested by beta_p < {ARREST_RATIO} beta_eq: "
          f"{prov['n_arrested']} of {prov['n_used']}")

    # ---------------------------------------------------------------- part 1
    print("\n" + "-" * 104)
    print("1. CAN THIS MATERIAL SHEAR-JAM AT ALL, BY THE REPO'S OWN CALIBRATION?")
    print("-" * 104)
    print("   Jamming needs phi_J(f) to reach phi, i.e. f_jam = (phi0-phi)/(phi0-phiM) < 1.")
    print(f"   {'phi':>6s} {'f_jam required':>16s} {'can jam?':>10s}   n")
    for phi in prov["phi_levels"]:
        d0 = next(d for d in drops if d.phi == phi)
        fj = jamming_deficit(d0)
        n = sum(1 for d in drops if d.phi == phi)
        print(f"   {phi:6.2f} {fj:16.4f} {'yes' if fj < 1 else 'NEVER':>10s}   {n}")
    m = material(0.44)
    print(f"\n   phi0 = {m.phi0}, phiM = {m.phiM} (Github/Drop Impact/closure/"
          f"params_fitted.yaml)")
    print("   Every packing in this dataset needs f > 1 to jam. The model says these")
    print("   drops cannot shear-jam at any stress. They demonstrably arrest.")

    # ---------------------------------------------------------------- part 2
    print("\n" + "-" * 104)
    print("2. THE MODEL'S OWN VISCOSITY, ACROSS THE WHOLE DATASET")
    print("-" * 104)
    span = model_viscosity_span(drops)
    print(f"   f at the impact stress:  {span['f_min']:.6f} to {span['f_max']:.6f}")
    print(f"   eta the closure predicts: {span['eta_min']:.5f} to "
          f"{span['eta_max']:.5f} Pa.s   ({span['eta_ratio']:.4f}x over the whole set)")
    print(f"   measured beta_p:          {prov['beta_p_range'][0]:.3f} to "
          f"{prov['beta_p_range'][1]:.3f}   "
          f"({prov['beta_p_range'][1] / prov['beta_p_range'][0]:.2f}x)")
    print(f"\n   But that {span['eta_ratio']:.1f}x is almost all phi, not stress. "
          f"At the single packing")
    print(f"   phi = {span['fixed_phi']:g} ({span['n_fixed_phi']} drops), over a "
          f"{span['u0_ratio_fixed_phi']:.2f}x range of impact speed:")
    print(f"      the closure's viscosity moves {span['eta_ratio_fixed_phi']:.4f}x")
    print(f"      the measured footprint moves  {span['beta_ratio_fixed_phi']:.2f}x")
    print(f"\n   Ratio: the response outruns the closure's stress response by "
          f"{np.log(span['beta_ratio_fixed_phi']) / np.log(span['eta_ratio_fixed_phi']):.0f}"
          f"x in log terms.")

    # ---------------------------------------------------------------- part 3
    print("\n" + "-" * 104)
    print("3. ELIMINATION AUDIT")
    print("-" * 104)
    print("   FALSIFIED: drops whose group agrees to 5% still span more than half the")
    print("              total log-range of beta_p.")
    print("   STRAINED : the group would need |d ln beta_p / d ln G| > 10.\n")
    verdicts = []
    for name, fn in CLOCKS.items():
        g = np.array([fn(d) for d in drops])
        v = audit_clock(name, g, beta)
        verdicts.append(v)
        print("   " + str(v))

    print("\n   survives:")
    survivors = [v for v in verdicts if v.verdict == "survives"]
    if survivors:
        for v in sorted(survivors, key=lambda x: -abs(x.spearman_rho)):
            print(f"      {v.name}   (rho {v.spearman_rho:+.3f}, "
                  f"p {v.spearman_p:.2g}, needs exponent {v.required_exponent:.2f})")
    else:
        print("      nothing")

    payload = {
        "provenance": prov,
        "arrest_ratio": ARREST_RATIO,
        "model_viscosity_span": span,
        "jamming_deficit": {
            f"{phi:g}": jamming_deficit(next(d for d in drops if d.phi == phi))
            for phi in prov["phi_levels"]
        },
        "verdicts": [asdict(v) | {"verdict": v.verdict, "g_range": v.g_range}
                     for v in verdicts],
    }
    (OUT / "audit.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"\nwrote {OUT / 'audit.json'}")
    figure(drops, beta, verdicts)
    return 0


def figure(drops, beta, verdicts) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 4.2))

    # Panel 1: the contrast that carries the result -- the strongest survivor
    # against the stress clock the audit falsifies. Both on their own x-axis,
    # because they live on incomparable scales; the point is the SHAPE, not a
    # common abscissa.
    best = max((v for v in verdicts if v.verdict == "survives"),
               key=lambda v: abs(v.spearman_rho))
    worst = min(verdicts, key=lambda v: v.g_range)
    for v, colour, marker, ax in ((best, "#55A868", "o", ax1),
                                  (worst, "#C44E52", "s", ax1.twiny())):
        g = np.array([CLOCKS[v.name](d) for d in drops])
        ax.plot(g, beta, marker, ms=4, color=colour, mfc="none")
        ax.set_xlabel(f"{v.name}   ({v.g_range:.3g}x, rho {v.spearman_rho:+.2f})",
                      fontsize=7.5, color=colour)
        ax.tick_params(axis="x", labelsize=7, colors=colour)
    ax1.set_ylabel(r"arrested footprint $\beta_p$")
    ax1.set_ylim(0.0, 2.1)
    ax1.tick_params(labelsize=8)

    # Panel 2: the audit itself -- dynamic range of each clock against the
    # dynamic range the response demands. Log axis because the clocks span
    # eight decades between them; the response's line is the thing to beat.
    names = [v.name for v in verdicts]
    ranges = np.array([v.g_range for v in verdicts])
    order = np.argsort(ranges)
    y = np.arange(len(order))
    resp = float(np.max(beta) / np.min(beta))
    colours = {"FALSIFIED": "#C44E52", "STRAINED": "#DD8452", "survives": "#55A868"}
    for yy, i in zip(y, order):
        v = verdicts[i]
        ax2.plot([1.0, max(ranges[i], 1.0)], [yy, yy], "-", lw=5,
                 color=colours[v.verdict], alpha=0.85, solid_capstyle="butt")
    ax2.axvline(resp, color="k", ls="--", lw=1.0)
    ax2.set_xscale("log")
    ax2.set_yticks(y)
    ax2.set_yticklabels([names[i] for i in order], fontsize=7)
    ax2.set_xlabel("dynamic range of the group over the 75 drops")
    ax2.tick_params(labelsize=8)
    ax2.invert_yaxis()

    fig.tight_layout()
    dest = OUT / "elimination_audit.png"
    fig.savefig(dest, dpi=200)
    plt.close(fig)
    print(f"wrote {dest}")


if __name__ == "__main__":
    raise SystemExit(main())
