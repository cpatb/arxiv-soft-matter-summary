"""Gate G5: on data that cannot answer the question, does the tool refuse?

The spec's wording was "3-point curves, Newtonian data, monotone no-thickening
data; declines to fit, never extrapolates." Made concrete below. A fitting
package is dangerous exactly where the data is too thin to say anything, because
least squares always returns a number and always attaches an error bar to it.

The case worth the most here is the last one. Wyart-Cates puts all of its stress
dependence in ``f(sigma) = exp(-sigma*/sigma)``, which moves from 0 to 1 over
about two decades around ``sigma*``. Measure only below that and every larger
``sigma*`` fits equally well; measure only above and every smaller one does.
The fit still converges, still reports a standard error, and the number is a
bound rather than a measurement. That is a property of the experiment, so it is
reported separately and in a form the user can act on: measure further, and this
far.

Pass conditions, fixed before the run:

  G5a  fewer points than parameters: raises, does not fit
  G5b  Newtonian data (flat viscosity) fit with Wyart-Cates: not usable
  G5c  monotone thinning data fit with Wyart-Cates: not usable
  G5d  all stresses far below the transition: reported as not bracketing, and
       sigma* carries a relative error above 100%
  G5e  all stresses far above the transition: reported as not bracketing
  G5f  the control: data that does bracket the transition is reported as
       bracketing, with sigma* to better than 50% relative. Without this, a
       checker that always cried wolf would pass G5d and G5e

Run:  python scripts/gate_g5.py
Exit code is 0 only if every gate passes.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dst_fit import assess, fit_model, synthetic_curve  # noqa: E402
from dst_fit.identifiability import transition_coverage  # noqa: E402
from dst_fit.models import WYART_CATES  # noqa: E402

TRUTH = (0.52, 0.58, 5.0)
SIGMA_STAR = TRUTH[2]
NOISE = 0.02

REL_ERR_UNCONSTRAINED = 1.0  # G5d: >100% relative error
REL_ERR_CONSTRAINED = 0.5  # G5f: <50% relative error


def fit_window(lo_decade: float, hi_decade: float, seed: int = 0):
    sigma = np.logspace(lo_decade, hi_decade, 25)
    s, e = synthetic_curve(WYART_CATES, TRUTH, sigma, noise=NOISE, seed=seed)
    fit = fit_model(WYART_CATES, s, e, n_starts=8, seed=seed)
    i = WYART_CATES.param_names.index("sigma_star")
    cov = transition_coverage(fit.params[i], s)
    rel = float(fit.stderr[i] / fit.params[i]) if fit.params[i] > 0 else np.inf
    return fit, cov, rel


def main() -> int:
    failures: list[str] = []

    # ---- G5a ------------------------------------------------------------
    s, e = synthetic_curve(WYART_CATES, TRUTH, np.logspace(0, 1, 3))
    try:
        fit_model(WYART_CATES, s, e)
        failures.append("G5a: fitted a 3-point curve with a 3-parameter model")
        print("  G5a FAIL: fitted an underdetermined system")
    except ValueError as exc:
        print(f"  G5a PASS: refused, {str(exc)[:70]}")

    # ---- G5b ------------------------------------------------------------
    sigma = np.logspace(-2, 3, 30)
    rng = np.random.default_rng(0)
    flat = np.full(sigma.shape, 0.05) * np.exp(rng.normal(0, NOISE, sigma.shape))
    a = assess(fit_model(WYART_CATES, sigma, flat, n_starts=8))
    print(f"  G5b {'PASS' if not a.usable else 'FAIL'}: Newtonian data -> "
          f"{'not usable' if not a.usable else 'USABLE'} "
          f"(well_specified={a.well_specified}, identifiable={a.identifiable})")
    if a.usable:
        failures.append("G5b: Newtonian data came back usable under Wyart-Cates")

    # ---- G5c ------------------------------------------------------------
    thinning = 1.0 * sigma**-0.4 * np.exp(rng.normal(0, NOISE, sigma.shape))
    a = assess(fit_model(WYART_CATES, sigma, thinning, n_starts=8))
    print(f"  G5c {'PASS' if not a.usable else 'FAIL'}: thinning data -> "
          f"{'not usable' if not a.usable else 'USABLE'}")
    if a.usable:
        failures.append("G5c: monotone thinning data came back usable under Wyart-Cates")

    # ---- G5d / G5e / G5f ------------------------------------------------
    print(f"\n  transition coverage (true sigma* = {SIGMA_STAR:g} Pa)")
    cases = [
        ("G5d", "far below", -4.0, -2.0, False),
        ("G5e", "far above", 3.0, 5.0, False),
        ("G5f", "bracketing", -1.0, 3.0, True),
    ]
    for tag, label, lo, hi, should_bracket in cases:
        fit, cov, rel = fit_window(lo, hi)
        i = WYART_CATES.param_names.index("sigma_star")
        fitted_ss = float(fit.params[i])
        ok = cov.brackets == should_bracket
        print(f"    {tag} {label:>10s}: sigma {cov.sigma_min:.3g} to {cov.sigma_max:.3g} Pa, "
              f"sigma*_fit {fitted_ss:10.4g} (true {SIGMA_STAR:g}), "
              f"f {cov.f_min:.3f} to {cov.f_max:.3f}, "
              f"brackets={str(cov.brackets):5s} rel err {rel:6.1%}  "
              f"{'ok' if ok else 'MISREPORTED'}")
        if not ok:
            failures.append(
                f"{tag}: brackets reported {cov.brackets}, expected {should_bracket}"
            )
        if tag == "G5d" and rel < REL_ERR_UNCONSTRAINED:
            failures.append(
                f"G5d: sigma* relative error only {rel:.1%} on data that cannot constrain it"
            )
        if tag == "G5f" and rel > REL_ERR_CONSTRAINED:
            failures.append(
                f"G5f: sigma* relative error {rel:.1%} on data that should constrain it"
            )

    print("\n  message the user would see, data entirely below the transition:")
    _, cov, _ = fit_window(-4.0, -2.0)
    print(f"      {cov}")

    print()
    if failures:
        print(f"GATE G5 FAIL ({len(failures)} problem(s)):")
        for f in failures:
            print("  " + f)
        return 1
    print("GATE G5 PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
