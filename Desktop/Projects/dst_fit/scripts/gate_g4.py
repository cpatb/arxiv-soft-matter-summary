"""Gate G4: the literature cross-check, and the control the corn starch needs.

G4 was written into the spec on 2026-08-27 as

    "Literature cross-check, 2 or more published DST papers.
     Recovered phi_m within the paper's own value +/- 0.02"

**One substitution, recorded rather than hidden.** The spec said *silica*. No
silica phi-series exists in the three deposits on disk: the genuine silica (Lin
et al. 2015 Fig 2a, and the two curves in Richards et al.'s Figure 3) is a single
volume fraction per system, and a single curve cannot test whether phi_m holds
still. So G4 runs on the two published systems that do have phi-series: Guy,
Hermes & Poon 2015 (PMMA, here) and Richards et al. 2019 (corn starch, in
``gate_n6.py``). Two papers, two labs, two materials, as the spec asked; not the
material the spec named.

Why this gate is worth more than a cross-check
-----------------------------------------------
N4 and N6 both failed the phi_m-consistency rule on corn starch, in two
independent datasets. That is only interesting if the rule is satisfiable. Guy
et al.'s PMMA is the control: hard, monodisperse, non-porous, non-swelling
spheres in a density-matched solvent, the system Wyart-Cates was built for and
the one whose phi_m the paper actually quotes. If phi_m holds still here and
drifts on corn starch, the difference is the material. If it drifts here too,
the pipeline is the problem and the corn-starch result means nothing.

Pass conditions, fixed before the data was fit:

  G4a  fitting eta_r = A (1 - phi/phi_c)^-n to the authors' own deposited
       ``d=3770 nm eta_2`` branch recovers phi_c within 0.02 of their published
       phi_m = 0.558(5). This tests the fitter against a published fit of the
       same numbers, so a failure here is a failure of dst-fit, not of the data.
  G4b  Wyart-Cates fits to the d = 3770 nm flow curves give an inverse-variance
       weighted-mean phi_m within 0.02 of 0.558.
  G4c  at least 4 of the admitted curves lie within 2 sigma of that weighted
       mean. G3's rule verbatim, via the same function N4 and N6 call.

  Admission, fixed with the rest and applied in code rather than by eye: a curve
  enters G4b and G4c only if its measured stresses bracket the frictional
  transition and it keeps more than MIN_BRANCH_POINTS = 4 points after the
  thickening-branch cut, so that a three-parameter fit has at least two degrees
  of freedom. A curve that stops below sigma* returns a confident-looking sigma* that
  is not a measurement, and its phi_m inherits that. The criterion is a property
  of the stress window, not of phi_m, so it cannot select for the answer; but it
  is computed from the *fitted* sigma*, so it can only catch the cases where the
  fit itself puts the transition outside the window. Fewer than 4 admitted
  curves fails the gate rather than lowering the bar.

Run:  python scripts/gate_g4.py
Exit code is 0 only if all three pass.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dst_fit import (  # noqa: E402
    assess,
    consistency_across_curves,
    fit_divergence,
    fit_model,
    thickening_branch,
    transition_coverage,
)
from dst_fit.io_literature import (  # noqa: E402
    GUY_FIG1A_PHI,
    GUY_FIG1B_PHI,
    GUY_PMMA_404,
    GUY_PMMA_3770,
    read_guy_branch,
    read_guy_flow_curves,
    read_guy_phi_m_column,
    verify_guy_size_assignment,
)
from dst_fit.models import WYART_CATES, WYART_CATES_KNOWN_PHI  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
GUY = ROOT / "Towards a unified description of the rheology of hard-particle suspensions"
FIG1 = GUY / "Fig 1" / "Fig 1"
FIG1A = FIG1 / "Fig 1a" / "Fig1a.xls"
FIG1B = FIG1 / "Fig 1b" / "Fig1b.xlsx"
FIG1C = FIG1 / "Fig 1c" / "Fig1c.xlsx"

G4_TOLERANCE = 0.02
G4_MIN_CONSISTENT = 4
G4_N_SIGMA = 2.0
MIN_BRANCH_POINTS = 4

# Particle diameters, nm. Used only for the sigma* ~ d^-2 observation.
D_FIG1A, D_FIG1B = 404.0, 3770.0


def fit_curves(sweeps, model, beta):
    """Fit, assess and check transition coverage for every curve in a workbook."""
    rows = []
    for s in sweeps:
        br = thickening_branch(s.viscosity)
        n_branch = len(range(*br.indices(s.viscosity.size)))
        if n_branch <= MIN_BRANCH_POINTS:
            rows.append({"sweep": s, "fit": None, "reason": f"{n_branch} branch points"})
            continue
        try:
            fit = fit_model(model, s.stress, s.viscosity, branch=br, n_starts=16)
        except (RuntimeError, ValueError) as exc:
            rows.append({"sweep": s, "fit": None, "reason": str(exc)})
            continue
        p = dict(zip(fit.model.param_names, fit.params))
        cov = transition_coverage(p["sigma_star"], s.stress[br], beta=beta)
        rows.append({
            "sweep": s, "fit": fit, "branch": br, "cov": cov,
            "assess": assess(fit), "admitted": cov.brackets, "reason": "",
        })
    return rows


def print_rows(rows):
    for r in rows:
        s = r["sweep"]
        label = s.metadata.get("phi_label", s.sheet)
        phi = f"{s.phi:.5f}" if s.phi is not None else "unknown"
        if r["fit"] is None:
            print(f"  {label:>6s}  phi {phi:>7s}   not fit: {r['reason']}")
            continue
        f, a, cov = r["fit"], r["assess"], r["cov"]
        p = dict(zip(f.model.param_names, f.params))
        e = dict(zip(f.model.param_names, f.stderr))
        flags = []
        if not a.well_specified:
            flags.append("struct")
        if not a.identifiable:
            flags.append("degen")
        if not cov.brackets:
            flags.append("NO-TRANSITION")
        print(
            f"  {label:>6s}  phi {phi:>7s}  n={f.n_points:2d}  "
            f"phi_fit {p['phi']:.4f}+/-{e['phi']:.4f}  "
            f"phi_m {p['phi_m']:.4f}+/-{e['phi_m']:.4f}  "
            f"sigma* {p['sigma_star']:9.2f}+/-{e['sigma_star']:<8.2f} "
            f"f {cov.f_min:.2f}-{cov.f_max:.2f}  rmse {f.rmse_log:.3f}  "
            f"{','.join(flags) if flags else 'usable'}"
        )


def main() -> int:
    failures: list[str] = []
    c = GUY_PMMA_3770  # the size the published phi_m was measured on
    c404 = GUY_PMMA_404
    print(f"{c.source}")
    print(f"  gate size:  {c}")
    print(f"  other size: {c404}\n")

    # ----------------------------------------------------------------------
    # provenance: the file names do not match the panel letters, so the size
    # assignment is established from the numbers before anything is fit.
    # ----------------------------------------------------------------------
    print("workbook-to-size assignment, checked against the phi-axis branch sheets:")
    checks = [
        verify_guy_size_assignment(FIG1A, FIG1C, "d=404 nm eta_1", GUY_FIG1A_PHI,
                                   eta_s=c404.eta_s),
        verify_guy_size_assignment(FIG1B, FIG1C, "d=3770 nm eta_1", GUY_FIG1B_PHI,
                                   eta_s=c.eta_s),
    ]
    for chk in checks:
        print("  " + str(chk).replace("\n", "\n  "))
    if not all(chk.matched for chk in checks):
        print("\nrefusing to fit: the size assignment could not be confirmed")
        return 1
    print("  Fig1a is d = 404 nm and Fig1b is d = 3770 nm, from the data, not the")
    print("  file names. Note the published panel letters are the other way round:")
    print("  Fig1c.xlsx holds the phi-axis compilation, which is panel (b).\n")

    # ----------------------------------------------------------------------
    # G4a: reproduce the published branch fit from the published numbers
    # ----------------------------------------------------------------------
    phi, eta_r = read_guy_branch(FIG1C, "d=3770 nm eta_2")
    branch_fit = fit_divergence(phi, eta_r, free_exponent=True, n_starts=24)
    tabulated = read_guy_phi_m_column(FIG1C)
    print("G4a: the authors' own upper branch, refit here")
    print("    " + branch_fit.summary().replace("\n", "\n    "))
    print(f"    published:  A = 0.20(9)   phi_m = {c.phi_m_published}"
          f"({c.phi_m_published_err * 1000:.0f}e-3)   n = {c.n}(2)")
    print(f"    phi_m also readable without fitting, from the deposited "
          f"'phi_m - phi' column: {tabulated:.3f}")
    d_a = branch_fit.phi_c - c.phi_m_published
    print(f"    phi_c - published phi_m = {d_a:+.4f} (tolerance +/-{G4_TOLERANCE})")
    if abs(d_a) > G4_TOLERANCE:
        failures.append(f"G4a: branch-fit phi_c off the published value by {d_a:+.4f}")
    if abs(branch_fit.A - 0.20) > 3 * max(branch_fit.stderr["A"], 0.09):
        print(f"    NOTE: A = {branch_fit.A:.3f} does not reproduce the published 0.20(9). "
              f"phi_m and n do, and they are what G4 tests; the published A is most "
              f"likely from the multi-system collapse in panel (b) rather than from "
              f"these six rows alone.")

    # frictionless poles, measured at the same exponent, for the phi_0 used below
    print("\n    frictionless branches, fit at the same exponent so the two poles are"
          " comparable:")
    poles = {}
    for sheet in ("d=3770 nm eta_1", "d=404 nm eta_1"):
        p, e = read_guy_branch(FIG1C, sheet)
        m = p > 0.35  # the dilute end carries no information about the pole
        f0 = fit_divergence(p[m], e[m], free_exponent=False, n=c.n, n_starts=24)
        poles[sheet] = (f0.phi_c, f0.stderr["phi_c"])
        print(f"      {sheet:>16s} (phi > 0.35, {int(m.sum())} pts): "
              f"phi_0 = {f0.phi_c:.4f} +/- {f0.stderr['phi_c']:.4f}  rmse {f0.rmse_log:.4f}")
    print(f"      using phi_0 = {c.phi_0} for the fits below")

    # ----------------------------------------------------------------------
    # G4b, G4c: Wyart-Cates on the d = 3770 nm flow curves
    # ----------------------------------------------------------------------
    model = WYART_CATES.with_fixed(**c.as_fixed())
    print(f"\nWyart-Cates on the thickening branch, {model.describe_fixed()}")

    print(f"\nd = {D_FIG1B:.0f} nm, {FIG1B.name}  (the size the published phi_m came from)")
    rows_b = fit_curves(read_guy_flow_curves(FIG1B, GUY_FIG1B_PHI), model, c.beta)
    print_rows(rows_b)

    admitted = [r for r in rows_b if r["fit"] is not None and r["admitted"]]
    rejected = [r for r in rows_b if r["fit"] is not None and not r["admitted"]]
    print(f"\n  admitted {len(admitted)} of {len(rows_b)}; "
          f"{len(rejected)} rejected for not bracketing the transition"
          + (": " + ", ".join(r["sweep"].sheet for r in rejected) if rejected else ""))

    if len(admitted) < G4_MIN_CONSISTENT:
        print(f"\nG4b/G4c FAIL: only {len(admitted)} curves admitted, "
              f"need {G4_MIN_CONSISTENT}")
        failures.append(f"G4b/G4c: only {len(admitted)} admitted curves")
    else:
        i = model.param_names.index("phi_m")
        vals = np.array([r["fit"].params[i] for r in admitted])
        errs = np.array([r["fit"].stderr[i] for r in admitted])
        check = consistency_across_curves(vals, errs, n_sigma=G4_N_SIGMA)

        d_b = check.mean - c.phi_m_published
        print(f"\nG4b: weighted-mean phi_m {check.mean:.4f} +/- {check.mean_err:.4f} "
              f"against the published {c.phi_m_published:.3f}: {d_b:+.4f} "
              f"(tolerance +/-{G4_TOLERANCE})")
        if abs(d_b) > G4_TOLERANCE:
            failures.append(f"G4b: weighted-mean phi_m off the published value by {d_b:+.4f}")

        print("\nG4c: phi_m " + check.report([r["sweep"].sheet for r in admitted]))
        print(f"    need {G4_MIN_CONSISTENT} of {check.n_total}")
        if not check.passes(G4_MIN_CONSISTENT):
            failures.append(
                f"G4c: {check.n_consistent} of {check.n_total} within "
                f"{G4_N_SIGMA:.0f} sigma, needed {G4_MIN_CONSISTENT}"
            )

        # sensitivity to the one constant that was chosen rather than published
        print("\n    sensitivity of the weighted-mean phi_m to phi_0:")
        for label, p0 in [
            ("the d = 404 nm pole", poles["d=404 nm eta_1"][0]),
            ("the d = 3770 nm pole, used above", poles["d=3770 nm eta_1"][0]),
            ("0.64, frictionless RCP, this package's old default", 0.64),
        ]:
            alt_model = WYART_CATES.with_fixed(**{**c.as_fixed(), "phi_0": p0})
            alt = [r for r in fit_curves(read_guy_flow_curves(FIG1B, GUY_FIG1B_PHI),
                                         alt_model, c.beta)
                   if r["fit"] is not None and r["admitted"]]
            if len(alt) < 2:
                continue
            v = np.array([r["fit"].params[i] for r in alt])
            e = np.array([r["fit"].stderr[i] for r in alt])
            a = consistency_across_curves(v, e, n_sigma=G4_N_SIGMA)
            print(f"      phi_0 = {p0:.4f}, {label:>48s}: "
                  f"phi_m {a.mean:.4f} +/- {a.mean_err:.4f}, "
                  f"{a.n_consistent}/{a.n_total} consistent, "
                  f"reduced chi2 {a.reduced_chi2:.1f}")

    # ----------------------------------------------------------------------
    # observations. None of this was pre-registered.
    # ----------------------------------------------------------------------
    print("\nobservations (not gated):")

    # The control that says whether a failed G4c means "the model is wrong" or
    # "one flow curve cannot measure phi_m". phi and phi_m enter only through
    # phi/phi_J, so near jamming a curve pins the ratio and not the two numbers.
    # Guy et al. knew phi: they made the samples. Pin it and refit.
    print("\n  phi_m with phi pinned at the authors' own volume fraction, not fitted:")
    known_rows = []
    for r in rows_b:
        s = r["sweep"]
        if r["fit"] is None or s.phi is None or not r["admitted"]:
            continue
        km = WYART_CATES_KNOWN_PHI.with_fixed(phi=s.phi, **c.as_fixed())
        try:
            kf = fit_model(km, s.stress, s.viscosity, branch=r["branch"], n_starts=16)
        except (RuntimeError, ValueError) as exc:
            print(f"      {s.sheet:>6s}  refused: {exc}")
            continue
        known_rows.append((s, kf, r))
        p = dict(zip(kf.model.param_names, kf.params))
        e = dict(zip(kf.model.param_names, kf.stderr))
        free = dict(zip(r["fit"].model.param_names, r["fit"].params))
        print(f"      {s.sheet:>6s}  phi {s.phi:.5f} (pinned)   "
              f"phi_m {p['phi_m']:.4f}+/-{e['phi_m']:.4f}   "
              f"sigma* {p['sigma_star']:8.2f}+/-{e['sigma_star']:<7.2f} "
              f"rmse {kf.rmse_log:.3f}   "
              f"[phi free gave phi_m {free['phi_m']:.4f} at phi {free['phi']:.4f}]")
    if len(known_rows) >= 2:
        v = np.array([kf.params[0] for _, kf, _ in known_rows])
        e = np.array([kf.stderr[0] for _, kf, _ in known_rows])
        kchk = consistency_across_curves(v, e, n_sigma=G4_N_SIGMA)
        print("      " + kchk.report([s.sheet for s, _, _ in known_rows]).replace(
            "\n", "\n      "))
        print(f"      weighted mean {kchk.mean:.4f} vs published "
              f"{c.phi_m_published:.3f}: {kchk.mean - c.phi_m_published:+.4f}")

    # The same trade-off, stated directly: what the free fit actually determines
    # is the distance to jamming at the top of the sweep, not phi_m.
    print("\n  distance to jamming, phi_m - phi, from the free fits:")
    for r in admitted:
        f = r["fit"]
        i_phi = f.model.param_names.index("phi")
        i_pm = f.model.param_names.index("phi_m")
        gap = f.params[i_pm] - f.params[i_phi]
        gap_err = float(np.sqrt(max(
            f.cov[i_pm, i_pm] + f.cov[i_phi, i_phi] - 2 * f.cov[i_phi, i_pm], 0.0)))
        print(f"      {r['sweep'].sheet:>6s}  phi_m - phi = {gap:.4f} +/- {gap_err:.4f}   "
              f"r(phi, phi_m) = {f.corr[i_phi, i_pm]:+.4f}")
    print(f"\n  d = {D_FIG1A:.0f} nm, {FIG1A.name}  (same particles, 9.3x smaller)")
    rows_a = fit_curves(read_guy_flow_curves(FIG1A, GUY_FIG1A_PHI), model, c.beta)
    print_rows(rows_a)

    adm_a = [r for r in rows_a if r["fit"] is not None and r["admitted"]]
    if adm_a and admitted:
        i = model.param_names.index("sigma_star")
        ss_a = np.median([r["fit"].params[i] for r in adm_a])
        ss_b = np.median([r["fit"].params[i] for r in admitted])
        print(f"\n    median fitted sigma*: {ss_a:.1f} Pa at {D_FIG1A:.0f} nm, "
              f"{ss_b:.2f} Pa at {D_FIG1B:.0f} nm, ratio {ss_a / ss_b:.1f}x")
        print(f"      published onset stresses: {c404.sigma_star_published:.0f} Pa and "
              f"{c.sigma_star_published:.0f} Pa, ratio "
              f"{c404.sigma_star_published / c.sigma_star_published:.0f}x. The paper's "
              f"own least-squares size scaling is sigma* ~ d^-1.9(1), predicting "
              f"{(D_FIG1B / D_FIG1A) ** 1.9:.1f}x (a clean d^-2 would give "
              f"{(D_FIG1B / D_FIG1A) ** 2:.1f}x).")
        print(f"      so the *ratio* is reproduced to "
              f"{100 * abs(ss_a / ss_b / (D_FIG1B / D_FIG1A) ** 1.9 - 1):.0f}%, while each "
              f"absolute value comes out {ss_a / c404.sigma_star_published:.1f}x and "
              f"{ss_b / c.sigma_star_published:.1f}x the published one. A near-constant "
              f"factor on both is a definition difference, not a disagreement: the paper "
              f"reads sigma* off where thickening visibly begins, while the fitted "
              f"sigma* is where f = 1/e, i.e. where the transition is already a third "
              f"complete. Neither number entered the other, so the ratio is a real "
              f"check on the size assignment.")

    for name, rows in (("d = 404 nm", rows_a), ("d = 3770 nm", rows_b)):
        usable = [r for r in rows if r["fit"] is not None and r["sweep"].phi is not None]
        if len(usable) < 3:
            continue
        j = model.param_names.index("phi")
        fitted = np.array([r["fit"].params[j] for r in usable])
        nominal = np.array([r["sweep"].phi for r in usable])
        slope, intercept = np.polyfit(nominal, fitted, 1)
        r2 = float(np.corrcoef(nominal, fitted)[0, 1] ** 2)
        print(f"\n    {name}: fitted phi vs the authors' phi over {len(usable)} curves: "
              f"slope {slope:.4f}, intercept {intercept:+.4f}, r2 {r2:.5f}")

    print("\n  what the paper says about its own two highest-phi curves:")
    print("    'The flow in both (a) and (c) was unsteady for phi >= 0.56, and points")
    print("    represent temporal averages.' The paper also reports wall slip at all")
    print("    stresses below phi_m and massive wall slip above it. Those are the 56%")
    print("    and 58% curves, which are exactly the two carrying the largest phi_m")
    print("    above. They are kept here because the admission rule was fixed in")
    print("    advance and says nothing about steadiness; excluding them after seeing")
    print("    the result would be choosing the answer. Recorded so the next reader")
    print("    can pre-register a steadiness criterion rather than inherit this one.")

    print()
    verdicts = {"G4a": True, "G4b": True, "G4c": True}
    for f in failures:
        verdicts[f.split(":")[0].split("/")[0]] = False
    print("  " + "   ".join(
        f"{k} {'PASS' if v else 'FAIL'}" for k, v in sorted(verdicts.items())
    ))
    if failures:
        print(f"GATE G4 FAIL ({len(failures)} problem(s)):")
        for f in failures:
            print("  " + f)
        return 1
    print("GATE G4 PASS (G4a, G4b, G4c)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
