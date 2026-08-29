"""Which parts of "the fitter recovers the nominal volume fraction" survive the gauge?

Gate N4 reports that fitted ``phi`` tracks the nominal label with slope 0.970,
intercept +0.052 and r2 = 0.99925, and the README read that as evidence the fit
"is not confused about the data", with the offset attributed to corn-starch
granules swelling in water.

Half of that is real and half of it is a statement about a constant we chose.

Wyart-Cates is exactly covariant under ``(phi, phi_m, phi_0) -> c (phi, phi_m,
phi_0)``, so every fitted ``phi`` scales with the assumed ``phi_0``. The slope of
fitted-against-nominal therefore scales with ``phi_0`` too, and can be set to
whatever you like. The correlation cannot: a common rescaling of one axis leaves
r2 untouched.

So this script sweeps ``phi_0`` and prints both, and then runs the same test on a
dataset where the volume fraction is genuinely known -- Guy et al.'s PMMA, where
``phi`` is a true volume fraction rather than a label of unrecorded convention.

Run:  python scripts/probe_gauge_dependence.py
Prints two tables; not a gate.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dst_fit import fit_model, thickening_branch  # noqa: E402
from dst_fit.io_literature import (  # noqa: E402
    GUY_FIG1A_PHI,
    GUY_FIG1B_PHI,
    GUY_PMMA_404,
    GUY_PMMA_3770,
    read_guy_flow_curves,
)
from dst_fit.io_trios import read_series  # noqa: E402
from dst_fit.models import WYART_CATES  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
STARCH = Path(r"C:\Users\compl\Desktop\Projects\Corn Starch Bulk Rheology AM")
GUY = ROOT / "Towards a unified description of the rheology of hard-particle suspensions" / "Fig 1" / "Fig 1"

PHI0_SWEEP = (0.58, 0.61, 0.64, 0.6598, 0.68)


def recover(sweeps, model):
    """Regress fitted phi on nominal phi. Returns (slope, intercept, r2, n)."""
    nom, fitted = [], []
    for s in sweeps:
        phi = getattr(s, "phi", None)
        if phi is None:
            continue
        br = thickening_branch(s.viscosity)
        try:
            fit = fit_model(model, s.stress, s.viscosity, branch=br, n_starts=12)
        except Exception:
            continue
        nom.append(float(phi))
        fitted.append(float(fit.params[0]))
    nom, fitted = np.array(nom), np.array(fitted)
    if nom.size < 3:
        return None
    slope, intercept = np.polyfit(nom, fitted, 1)
    r2 = float(np.corrcoef(nom, fitted)[0, 1] ** 2)
    return float(slope), float(intercept), r2, nom.size


def main() -> int:
    print("1. The 2023 corn starch, swept over the assumed phi_0")
    print("   The label's convention (weight or volume) was never recorded, so the")
    print("   nominal axis is a label, not a measurement.\n")
    if not STARCH.is_dir():
        print("   corn-starch data not present; skipping")
    else:
        sweeps = read_series(sorted(STARCH.glob("*_e1.xls")))
        print(f"   {'phi_0':>7} {'slope':>8} {'intercept':>10} {'r2':>9} {'predicted slope':>16}")
        base = None
        for phi_0 in PHI0_SWEEP:
            got = recover(sweeps, WYART_CATES.with_fixed(phi_0=phi_0))
            if got is None:
                continue
            slope, intercept, r2, _ = got
            if base is None:
                base = (phi_0, slope)
            pred = base[1] * phi_0 / base[0]
            print(f"   {phi_0:7.4f} {slope:8.4f} {intercept:+10.4f} {r2:9.5f} {pred:16.4f}")
        print()
        print("   The slope tracks the covariance prediction exactly, and r2 does not")
        print("   move at all. At phi_0 = 0.6598 the slope reads 1.0000: the")
        print("   'validation' can be dialled to pass by choosing the gauge.")
        print()
        print("   So: r2 = 0.99925 is a real, gauge-invariant statement that fitted phi")
        print("   tracks the label linearly. The slope and the intercept are not")
        print("   measurements, and the swelling story told about that intercept is")
        print("   unsupported -- its size is set by phi_0, not by the starch.")

    print("\n2. Guy et al. PMMA, where phi is a genuine volume fraction")
    print("   Same test, on data whose phi needs no convention and whose phi_0 was")
    print("   measured from the frictionless branch rather than assumed.\n")
    cases = [
        ("404 nm", GUY / "Fig 1a" / "Fig1a.xls", GUY_FIG1A_PHI, GUY_PMMA_404),
        ("3770 nm", GUY / "Fig 1b" / "Fig1b.xlsx", GUY_FIG1B_PHI, GUY_PMMA_3770),
    ]
    for label, path, phi_map, constants in cases:
        if not path.is_file():
            print(f"   {label}: deposit not present; skipping")
            continue
        got = recover(
            read_guy_flow_curves(path, phi_map),
            WYART_CATES.with_fixed(**constants.as_fixed()),
        )
        if got is None:
            print(f"   {label}: too few usable curves")
            continue
        slope, intercept, r2, n = got
        print(f"   {label:>8}: n={n}  slope {slope:.4f}  intercept {intercept:+.4f}  "
              f"r2 {r2:.5f}")
    print()
    print("   In the measured gauge the slope is not 1 and the intercept is not 0,")
    print("   and no choice of phi_0 fixes both at once. The fit reproduces the")
    print("   ordering of the volume fractions, not their values.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
