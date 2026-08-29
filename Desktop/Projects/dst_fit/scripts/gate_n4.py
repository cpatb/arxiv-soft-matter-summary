"""Gate N4 / G3: fit the real corn-starch phi-series and ask whether phi_m holds still.

phi_m is the frictional jamming fraction. It is a property of the particles, not
of how many of them are in the pot, so six independent fits of the same corn
starch at six volume fractions have to return the same phi_m. If they do not,
either Wyart-Cates does not describe this material or the fit is not measuring
what it claims. Either answer is worth more than the package.

Pass condition, fixed in the spec on 2026-08-27 before any real data was fit and
NOT revised after seeing the result:

  G3  phi_m consistent across at least 4 of the 6 curves within 2 sigma

Made precise here, since the spec wording needed one choice: consistency is
measured against the inverse-variance weighted mean of all six phi_m estimates,
and a curve counts as consistent when |phi_m_i - weighted mean| <= 2 * sigma_i.
That choice is recorded rather than tuned; the alternative (all fifteen pairwise
comparisons) is stricter and is reported alongside for information.

Everything below the gate is observation, not a pass condition, because none of
it was pre-registered.

Run:  python scripts/gate_n4.py
Exit code is 0 only if G3 passes.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dst_fit import (  # noqa: E402
    assess,
    consistency_across_curves,
    fit_model,
    read_series,
    thickening_branch,
    weighted_mean,
)

DATA_DIR = Path(r"C:\Users\compl\Desktop\Projects\Corn Starch Bulk Rheology AM")

G3_MIN_CONSISTENT = 4
G3_N_SIGMA = 2.0


def main() -> int:
    sweeps = read_series(sorted(DATA_DIR.glob("*_e1.xls")))
    print(f"corn starch, {len(sweeps)} sweeps, Wyart-Cates on the thickening branch\n")

    rows = []
    for s in sweeps:
        br = thickening_branch(s.viscosity)
        try:
            fit = fit_model("wyart_cates", s.stress, s.viscosity, branch=br, n_starts=12)
        except (RuntimeError, ValueError) as exc:
            print(f"  phi_nom {s.phi:.2f}: refused to fit ({exc})")
            rows.append(None)
            continue
        a = assess(fit)
        rows.append({"sweep": s, "fit": fit, "assess": a, "branch": br})

        p = dict(zip(fit.model.param_names, fit.params))
        e = dict(zip(fit.model.param_names, fit.stderr))
        flags = []
        if not a.well_specified:
            flags.append("struct")
        if not a.identifiable:
            flags.append("degen")
        print(
            f"  phi_nom {s.phi:.2f}  n={fit.n_points:2d}  "
            f"phi_fit {p['phi']:.4f}+/-{e['phi']:.4f}   "
            f"phi_m {p['phi_m']:.4f}+/-{e['phi_m']:.4f}   "
            f"sigma* {p['sigma_star']:8.3f}+/-{e['sigma_star']:<7.3f}  "
            f"rmse {fit.rmse_log:.3f}  {','.join(flags) if flags else 'usable'}"
        )

    good = [r for r in rows if r is not None]
    if len(good) < G3_MIN_CONSISTENT:
        print(f"\nG3 FAIL: only {len(good)} curves fit at all")
        return 1

    phi_m = np.array([r["fit"].params[r["fit"].model.param_names.index("phi_m")] for r in good])
    sig = np.array([r["fit"].stderr[r["fit"].model.param_names.index("phi_m")] for r in good])

    check = consistency_across_curves(phi_m, sig, n_sigma=G3_N_SIGMA)
    consistent = check.n_consistent
    print("\nG3: phi_m " + check.report([f"phi_nom {r['sweep'].phi:.2f}" for r in good]))
    print(f"    (need {G3_MIN_CONSISTENT})")

    print("\nobservations (not gated, none of this was pre-registered):")
    print("    fitted phi vs nominal phi (corn-starch granules swell, so these are "
          "not expected to match exactly):")
    for r in good:
        f = r["fit"]
        pf = f.params[f.model.param_names.index("phi")]
        print(f"      nominal {r['sweep'].phi:.2f}  ->  fitted {pf:.4f}")
    sstar = np.array([r["fit"].params[r["fit"].model.param_names.index("sigma_star")] for r in good])
    print(f"    sigma* across the series: {np.min(sstar):.3f} to {np.max(sstar):.3f} Pa "
          f"(ratio {np.max(sstar) / np.min(sstar):.1f}x); "
          f"Wyart-Cates expects it phi-independent")

    # Does the fit recover the label? The nominal phi is the only external
    # reference available, so this is the closest thing to ground truth here.
    fitted = np.array([r["fit"].params[r["fit"].model.param_names.index("phi")] for r in good])
    nominal = np.array([r["sweep"].phi for r in good])
    slope, intercept = np.polyfit(nominal, fitted, 1)
    r2 = float(np.corrcoef(nominal, fitted)[0, 1] ** 2)
    print(f"    fitted phi vs nominal phi: slope {slope:.4f}, intercept {intercept:+.4f}, "
          f"r2 {r2:.5f}")

    # Is the failure a phi/phi_m degeneracy, or a real drift in phi_m?
    gaps, gerrs, corrs = [], [], []
    for r in good:
        f = r["fit"]
        i = f.model.param_names.index("phi")
        j = f.model.param_names.index("phi_m")
        gaps.append(f.params[j] - f.params[i])
        gerrs.append(np.sqrt(max(f.cov[j, j] + f.cov[i, i] - 2 * f.cov[i, j], 0.0)))
        corrs.append(f.corr[i, j])
    gaps, gerrs = np.array(gaps), np.array(gerrs)
    gmean, _ = weighted_mean(gaps, gerrs)
    gchi2 = float(np.sum(((gaps - gmean) / gerrs) ** 2))
    print(f"    distance to jamming phi_m - phi: "
          f"{gaps.min():.4f} to {gaps.max():.4f}, "
          f"chi2 for a common value {gchi2:.1f} on {len(good) - 1} dof")
    print("    r(phi, phi_m) per curve: " + ", ".join(f"{c:.3f}" for c in corrs))

    print()
    if consistent >= G3_MIN_CONSISTENT:
        print("GATE N4 PASS (G3)")
        return 0
    print(f"GATE N4 FAIL: G3 wanted {G3_MIN_CONSISTENT} of {len(good)} "
          f"within {G3_N_SIGMA:.0f} sigma, got {consistent}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
