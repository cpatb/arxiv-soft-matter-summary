"""Gate N6 / G7: the same test as N4, on somebody else's corn starch.

N4 found phi_m sliding from 0.492 to 0.580 across six fits of one 2023
corn-starch series, reduced chi-squared 32.4 on 5 dof, where Wyart-Cates says
phi_m is one number. N5 killed the cheapest explanation. What N4 could not say
is whether that drift is a property of corn starch or a property of that data.

This answers it. Richards, Royer, Liebchen, Guy & Poon, PRL 123, 038004 (2019),
deposited seven time-averaged flow curves of corn starch in 50 wt% glycerol-
water: different lab, different solvent, different rheometer, different
operator, seven weight fractions from 0.38 to 0.50. Same pipeline, same rule.

  - if phi_m drifts here too, the drift belongs to corn starch and Wyart-Cates
    is being used outside its assumptions, which is a result;
  - if phi_m holds still here, the 2023 series has a local fault.

The constants are theirs, not ours
----------------------------------
Their Figure 1 caption and the deposit's own plotting script both give
phi_m = 0.457, phi_rcp = 0.546, beta = 0.94, sigma* = 5.1 Pa, eta_s = 6 mPa.s,
exponent 2. Every one of those differs from this package's defaults, and running
at the defaults would compare a phi_m measured under one set of constants with a
phi_m published under another. Since 2026-08-27 they are arguments, so this gate
binds them.

**They use weight fractions and say why**: "Here we use weight fractions, phi,
due to the porosity of our main model system (cornstarch)". So there is no
weight-to-volume conversion to apply, and phi_0 = 0.546 rather than 0.64 is not
a fudge, it is the frictionless jamming point measured in those units. The
conversion is quantified below as an observation, because whether it is legitimate
is candidate (3) of the N4 open questions, but it is not applied.

Pass conditions, fixed before the data was fit and not revised afterwards:

  G7a  phi_m consistent across at least 4 of the 7 curves, each within 2 sigma
       of the inverse-variance weighted mean of all seven. This is G3's rule
       verbatim, via the same ``consistency_across_curves`` the N4 gate calls,
       so the two results are comparable and neither can be tuned separately.
  G7b  that weighted mean lies within 0.02 of the published phi_m = 0.457. This
       is G4's tolerance applied to a second published system.

G7a is the one that answers N4. G7b asks the different question of whether the
package reproduces a published number at all.

Run:  python scripts/gate_n6.py
Exit code is 0 only if both pass.
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
    thickening_branch,
    transition_coverage,
)
from dst_fit.io_literature import (  # noqa: E402
    RICHARDS_CORNSTARCH,
    mass_to_volume_fraction,
    read_richards_series,
)
from dst_fit.joint import aic_from_residuals, fit_series  # noqa: E402
from dst_fit.models import WYART_CATES  # noqa: E402

DATA_DIR = Path(__file__).resolve().parents[1] / (
    "Competing Timescales Lead to Oscillations in Shear-Thickening Suspensions"
) / "Figure1" / "Plotting"

G7_MIN_CONSISTENT = 4
G7_N_SIGMA = 2.0
G7B_TOLERANCE = 0.02

# Dry corn-starch granule density and 50 wt% glycerol-water density at 20 degC,
# g/cm3. Used once, for the observation that quantifies what a weight-to-volume
# conversion would do. Not used by either gate; see mass_to_volume_fraction.
RHO_CORNSTARCH_DRY = 1.55
RHO_GLYCEROL_WATER_50 = 1.126


def fit_series_at(constants, sweeps, *, phi_0=None, beta=None, n_starts=12):
    """Fit every sweep with Wyart-Cates under one set of fixed constants."""
    fixed = constants.as_fixed()
    if phi_0 is not None:
        fixed["phi_0"] = phi_0
    if beta is not None:
        fixed["beta"] = beta
    model = WYART_CATES.with_fixed(**fixed)

    rows = []
    for s in sweeps:
        br = thickening_branch(s.viscosity)
        try:
            fit = fit_model(model, s.stress, s.viscosity, branch=br, n_starts=n_starts)
        except (RuntimeError, ValueError) as exc:
            print(f"  phi {s.phi:.2f}: refused to fit ({exc})")
            continue
        rows.append({"sweep": s, "fit": fit, "branch": br})
    return model, rows


def phi_m_of(rows):
    i = rows[0]["fit"].model.param_names.index("phi_m")
    return (
        np.array([r["fit"].params[i] for r in rows]),
        np.array([r["fit"].stderr[i] for r in rows]),
    )


def main() -> int:
    failures: list[str] = []
    sweeps = read_richards_series(DATA_DIR)
    c = RICHARDS_CORNSTARCH

    print(f"{c.name}\n  {c.source}")
    print(f"  constants, as published: {c}")
    print(f"  {len(sweeps)} curves, weight fractions "
          f"{min(s.phi for s in sweeps):.2f} to {max(s.phi for s in sweeps):.2f}\n")

    model, rows = fit_series_at(c, sweeps)
    print(f"Wyart-Cates on the thickening branch, {model.describe_fixed()}\n")
    for r in rows:
        f, s = r["fit"], r["sweep"]
        p = dict(zip(f.model.param_names, f.params))
        e = dict(zip(f.model.param_names, f.stderr))
        a = assess(f)
        cov = transition_coverage(p["sigma_star"], s.stress[r["branch"]], beta=c.beta)
        flags = []
        if not a.well_specified:
            flags.append("struct")
        if not a.identifiable:
            flags.append("degen")
        if not cov.brackets:
            flags.append("no-transition")
        print(
            f"  phi {s.phi:.2f}  n={f.n_points:2d}  "
            f"phi_fit {p['phi']:.4f}+/-{e['phi']:.4f}   "
            f"phi_m {p['phi_m']:.4f}+/-{e['phi_m']:.4f}   "
            f"sigma* {p['sigma_star']:7.3f}+/-{e['sigma_star']:<6.3f}  "
            f"f {cov.f_min:.2f}-{cov.f_max:.2f}  "
            f"rmse {f.rmse_log:.3f}  {','.join(flags) if flags else 'usable'}"
        )

    if len(rows) < G7_MIN_CONSISTENT:
        print(f"\nG7a FAIL: only {len(rows)} curves fit at all")
        return 1

    phi_m, sig = phi_m_of(rows)
    check = consistency_across_curves(phi_m, sig, n_sigma=G7_N_SIGMA)
    print("\nG7a: phi_m " + check.report([f"phi {r['sweep'].phi:.2f}" for r in rows]))
    print(f"    need {G7_MIN_CONSISTENT} of {check.n_total}")
    if not check.passes(G7_MIN_CONSISTENT):
        failures.append(
            f"G7a: {check.n_consistent} of {check.n_total} within "
            f"{G7_N_SIGMA:.0f} sigma, needed {G7_MIN_CONSISTENT}"
        )

    delta = check.mean - c.phi_m_published
    print(f"\nG7b: weighted-mean phi_m {check.mean:.4f} +/- {check.mean_err:.4f} "
          f"against the published {c.phi_m_published:.3f}: {delta:+.4f} "
          f"(tolerance +/-{G7B_TOLERANCE})")
    if abs(delta) > G7B_TOLERANCE:
        failures.append(f"G7b: weighted-mean phi_m off the published value by {delta:+.4f}")

    # ----------------------------------------------------------------------
    # observations. None of this was pre-registered.
    # ----------------------------------------------------------------------
    print("\nobservations (not gated):")

    # The strongest available check, and it needs no fit at all: take the
    # authors' own five numbers, put the nominal weight fraction in for phi, and
    # ask how well their published curve describes their published data. If it
    # describes it well, then a per-curve phi_m that drifts is a statement about
    # per-curve freedom, not about the model.
    published = WYART_CATES.with_fixed(**c.as_fixed())
    print("    the authors' own parameters, forward-evaluated, no fitting:")
    print(f"      phi_m={c.phi_m_published}, sigma*={c.sigma_star_published} Pa, "
          f"phi = the nominal weight fraction")
    rmses = []
    for r in rows:
        s, br = r["sweep"], r["branch"]
        pred = published(s.stress[br], s.phi, c.phi_m_published, c.sigma_star_published)
        obs = s.viscosity[br]
        ok = np.isfinite(pred) & (pred > 0)
        rms = float(np.sqrt(np.mean((np.log(pred[ok]) - np.log(obs[ok])) ** 2)))
        rmses.append(rms)
        jam = int(np.count_nonzero(~ok))
        print(f"      phi {s.phi:.2f}   rmse(log eta) {rms:.3f}   "
              f"fitted-per-curve rmse {r['fit'].rmse_log:.3f}"
              + (f"   ({jam} of {ok.size} points predicted jammed)" if jam else ""))
    print(f"      published-parameter rmse across the series: "
          f"{np.mean(rmses):.3f} mean, {np.max(rmses):.3f} worst")

    sstar = np.array([r["fit"].params[2] for r in rows])
    print(f"    sigma* across the series: {sstar.min():.3f} to {sstar.max():.3f} Pa "
          f"(ratio {sstar.max() / sstar.min():.1f}x); published {c.sigma_star_published} Pa, "
          f"phi-independent")

    fitted = np.array([r["fit"].params[0] for r in rows])
    nominal = np.array([r["sweep"].phi for r in rows])
    slope, intercept = np.polyfit(nominal, fitted, 1)
    r2 = float(np.corrcoef(nominal, fitted)[0, 1] ** 2)
    print(f"    fitted phi vs nominal weight fraction: slope {slope:.4f}, "
          f"intercept {intercept:+.4f}, r2 {r2:.5f}")

    # What a weight-to-volume conversion would do, quantified but not applied.
    phi_v = mass_to_volume_fraction(nominal, RHO_CORNSTARCH_DRY, RHO_GLYCEROL_WATER_50)
    print(f"    if the weight fractions were converted to volume at "
          f"rho_p={RHO_CORNSTARCH_DRY} and rho_f={RHO_GLYCEROL_WATER_50} g/cm3, "
          f"{nominal.min():.2f}-{nominal.max():.2f} would become "
          f"{phi_v.min():.3f}-{phi_v.max():.3f}, a ratio running "
          f"{(phi_v / nominal).min():.3f} to {(phi_v / nominal).max():.3f}. "
          f"Not constant, so the conversion is not a change of units. Not applied.")

    # N5's control, applied here. Per-curve phi_m and phi are correlated, so the
    # drift could in principle be phi_m absorbing per-curve phi. Forcing one
    # phi_m and one sigma* on the whole series, with phi still free per curve,
    # is the best case for a shared phi_m; AIC then says whether it survives
    # paying for the parameters it saved.
    curves = [(r["sweep"].stress[r["branch"]], r["sweep"].viscosity[r["branch"]])
              for r in rows]
    indep_residuals = np.concatenate([r["fit"].residuals for r in rows])
    indep_aic = aic_from_residuals(indep_residuals, 3 * len(rows))
    print(f"\n    joint fit, one phi_m and one sigma* for all {len(rows)} curves "
          f"(phi_0 free, as in N5):")
    print(f"      independent: k={3 * len(rows)} n={indep_residuals.size} "
          f"rmse={float(np.sqrt(np.mean(indep_residuals ** 2))):.4f} AIC={indep_aic:.1f}")
    try:
        jf = fit_series(curves, shared_sigma_star=True, eta_s=c.eta_s, n=c.n,
                        beta=c.beta, n_starts=16)
    except RuntimeError as exc:
        print(f"      joint fit failed to converge: {exc}")
    else:
        for line in jf.summary().splitlines():
            print("      " + line)
        d_aic = jf.aic - indep_aic
        print(f"      AIC vs independent: {d_aic:+.1f} "
              f"({'joint wins' if d_aic < 0 else 'independent wins'})")
        j_slope, j_int = np.polyfit(nominal, jf.phi, 1)
        j_r2 = float(np.corrcoef(nominal, jf.phi)[0, 1] ** 2)
        print(f"      phi vs nominal: slope {j_slope:.4f}, intercept {j_int:+.4f}, "
              f"r2 {j_r2:.5f}")

    # How much of the answer is the constants rather than the data.
    print("\n    sensitivity of the weighted-mean phi_m to the fixed constants:")
    for label, kwargs in [
        ("as published (phi_0=0.546, beta=0.94)", {}),
        ("beta = 1.0, the original WC form", {"beta": 1.0}),
        ("phi_0 = 0.64, this package's old default", {"phi_0": 0.64}),
        ("phi_0 = 0.64 and beta = 1.0, i.e. the pre-2026-08-27 code",
         {"phi_0": 0.64, "beta": 1.0}),
    ]:
        _, alt_rows = fit_series_at(c, sweeps, **kwargs)
        if len(alt_rows) < 2:
            print(f"      {label:>46s}: too few curves fit")
            continue
        v, e = phi_m_of(alt_rows)
        alt = consistency_across_curves(v, e, n_sigma=G7_N_SIGMA)
        print(f"      {label:>46s}: phi_m {alt.mean:.4f} +/- {alt.mean_err:.4f}, "
              f"spread {v.min():.4f}-{v.max():.4f}, "
              f"{alt.n_consistent}/{alt.n_total} consistent, "
              f"reduced chi2 {alt.reduced_chi2:.1f}")

    print()
    if failures:
        print(f"GATE N6 FAIL ({len(failures)} problem(s)):")
        for f in failures:
            print("  " + f)
        return 1
    print("GATE N6 PASS (G7a, G7b)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
