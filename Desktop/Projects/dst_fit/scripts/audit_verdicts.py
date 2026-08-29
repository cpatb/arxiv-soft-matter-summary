"""How often does the specification verdict fire, and how often should it?

``assess()`` fired MISSPECIFIED on fourteen of the fifteen real curves fitted
anywhere in this project -- every one except Lin's 2 um silica. A verdict that
almost always fires is barely a verdict, and the reason is that the Wald-Wolfowitz
runs test measures whether structure is *detectable*, which on a precise 25-point
sweep it nearly always is.

:func:`dst_fit.identifiability.systematic_misfit` separates that from *how big*
the structure is. This script runs both over every real curve on disk and
tabulates them, so the split between "the model cannot make this shape" and "the
model is slightly off and the data is precise" is a measured table rather than an
assertion.

The threshold it is judged against, ``SYSTEMATIC_TOL = 0.05``, was set from the
reproducibility of a rheometer flow curve before this script was first run, and
is not adjusted afterwards.

Run:  python scripts/audit_verdicts.py
Prints a table; always exits 0. This is an audit, not a gate.
"""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dst_fit import assess, fit_model, thickening_branch  # noqa: E402
from dst_fit.identifiability import SYSTEMATIC_TOL  # noqa: E402
from dst_fit.io_literature import (  # noqa: E402
    GUY_FIG1A_PHI,
    GUY_FIG1B_PHI,
    GUY_PMMA_404,
    GUY_PMMA_3770,
    RICHARDS_CORNSTARCH,
    read_guy_flow_curves,
    read_lin_reversal,
    read_richards_series,
)
from dst_fit.models import WYART_CATES  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RICHARDS_DIR = (
    ROOT
    / "Competing Timescales Lead to Oscillations in Shear-Thickening Suspensions"
    / "Figure1"
    / "Plotting"
)
GUY = ROOT / "Towards a unified description of the rheology of hard-particle suspensions" / "Fig 1" / "Fig 1"
LIN_DIR = ROOT / (
    "Hydrodynamic and Contact Contributions to Continuous Shear Thickening "
    "in Colloidal Suspensions"
) / "Fig2" / "Fig2"


def rows_from(sweeps, model, source):
    out = []
    for s in sweeps:
        br = thickening_branch(s.viscosity)
        try:
            fit = fit_model(model, s.stress, s.viscosity, branch=br, n_starts=16)
        except Exception as exc:  # a curve the fitter refuses is still a row
            out.append((source, s.sheet, None, str(exc)[:40]))
            continue
        out.append((source, s.sheet, assess(fit), ""))
    return out


def main() -> int:
    rows = []

    if RICHARDS_DIR.is_dir():
        model = WYART_CATES.with_fixed(**RICHARDS_CORNSTARCH.as_fixed())
        rows += rows_from(read_richards_series(RICHARDS_DIR), model, "Richards starch")

    fig1a = GUY / "Fig 1a" / "Fig1a.xls"
    fig1b = GUY / "Fig 1b" / "Fig1b.xlsx"
    if fig1a.is_file():
        rows += rows_from(
            read_guy_flow_curves(fig1a, GUY_FIG1A_PHI),
            WYART_CATES.with_fixed(**GUY_PMMA_404.as_fixed()),
            "Guy PMMA 404nm",
        )
    if fig1b.is_file():
        rows += rows_from(
            read_guy_flow_curves(fig1b, GUY_FIG1B_PHI),
            WYART_CATES.with_fixed(**GUY_PMMA_3770.as_fixed()),
            "Guy PMMA 3770nm",
        )

    if LIN_DIR.is_dir():
        model = WYART_CATES.with_fixed(eta_s=1.0, phi_0=0.64)
        for system, label in (("Fig_2a", "Lin 2um silica"), ("Fig_2c", "Lin 1.4um PMMA")):
            split = read_lin_reversal(LIN_DIR, system)
            br = thickening_branch(split.sweep_eta_r)
            fit = fit_model(
                model, split.sweep_sigma, split.sweep_eta_r, branch=br, n_starts=16
            )
            rows.append((label, system, assess(fit), ""))

    if not rows:
        print("no deposits present; nothing to audit")
        return 0

    print("Specification verdict on every real curve on disk.")
    print("runs_p asks whether structure is DETECTABLE; sys_amp asks how BIG it is,")
    print(f"in log viscosity. The model is blamed only above sys_amp = {SYSTEMATIC_TOL}.\n")
    print(f"{'dataset':<18} {'curve':<10} {'n':>4} {'rmse':>7} {'runs_p':>9} "
          f"{'sys_amp':>8} {'sys_frac':>8}  verdict")
    n_detect = n_blame = n_total = 0
    for source, label, a, err in rows:
        if a is None:
            print(f"{source:<18} {label:<10} {'':>4} {'':>7} {'':>9} {'':>8} {'':>8}  {err}")
            continue
        n_total += 1
        n_detect += not a.well_specified
        n_blame += a.structurally_wrong
        tag = (
            "MISSPECIFIED" if a.structurally_wrong
            else "slightly off" if not a.well_specified
            else "describes"
        )
        print(f"{source:<18} {str(label):<10} {a.fit.n_points:>4} {a.rmse_log:>7.4f} "
              f"{a.runs_p:>9.2g} {a.systematic_amplitude:>8.4f} "
              f"{a.systematic_fraction:>8.2f}  {tag}")

    print(f"\n  structure detected on {n_detect} of {n_total} curves "
          f"(this is what the old verdict reported)")
    print(f"  large enough to reject the model on {n_blame} of {n_total}")
    print(f"  the remaining {n_detect - n_blame} are the model being slightly wrong "
          f"and the data being good enough to notice")
    print()
    print("  An external check nobody arranged: on the d = 3770 nm series the two")
    print("  curves this rejects are 56% and 58%, which are exactly the two Guy et al.")
    print("  say were unsteady, with massive wall slip, at phi >= 0.56. The threshold")
    print("  was set from instrument reproducibility and never told that.")
    print()
    print("  The blind spot, stated because it is real: rejection needs BOTH a")
    print("  significant runs test and a large amplitude, so a short sweep cannot be")
    print("  rejected at all. Guy's d = 404 nm 59% and 61% curves carry amplitudes of")
    print("  0.077 and 0.059, above the threshold, but only 11 points, so the runs test")
    print("  reaches p = 0.058 and the conjunction never fires. The paper calls those")
    print("  unsteady too. A criterion for short sweeps is an open question, not")
    print("  something to fix by loosening RUNS_ALPHA after seeing this table.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
