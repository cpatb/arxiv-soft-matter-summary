"""Gate N7 / G9: is the phi_m spread real, or were the error bars wrong?

N4, N6 and G4 all rejected a common ``phi_m`` across a phi-series, on the
strength of standard errors taken from the Jacobian at the optimum. Those errors
come out absurdly small on the curves nearest jamming: one Guy et al. curve
returns ``phi_m = 0.5876 +/- 0.0013`` from a fit whose ``phi_m - phi`` is itself
0.0013. If those errors are too small, the whole disagreement is an artefact of
the error bar and there is nothing to explain.

Two things could be wrong with them, and they need separating.

**The shape of the cost surface.** The Wald error is a quadratic approximation.
Wyart-Cates diverges as ``phi -> phi_J``, so the surface should be strongly
asymmetric in ``phi_m`` and the approximation should fail. A profile likelihood
follows the real surface and would show it.

**The noise model.** Both the Wald and the profile interval assume the residuals
are independent draws. The runs test fires on every real curve in this project,
so on this data they are not. A moving-block bootstrap with blocks sized from
the observed run structure drops that assumption while keeping the rest.

The two are tested separately on purpose. If the first explains the spread, the
fix is a better interval; if the second does, the fix is a better model, and no
amount of interval arithmetic will help.

Pass conditions, fixed before the run:

  G9a  On synthetic data drawn from the model itself, the 95% profile interval
       covers the truth in at least 80 of 100 realisations. Same bar as G1b, on
       the same nominal 95%. Establishes that the machinery is not broken before
       any conclusion is drawn from it.
  G9b  On the same synthetic data, the 95% block-bootstrap interval also covers
       the truth in at least 80 of 100. A method that is wide enough to swallow
       real disagreement must first be shown not to be wide by construction.
  G9c  The verdict must be decidable: on the real Guy et al. series, the
       block-bootstrap intervals must either all overlap (the spread is an
       error-bar artefact) or fail to overlap (the spread is real). This gate
       passes as long as the intersection test returns an answer on at least 4
       curves; it does not prescribe which answer.

G9c is deliberately not a bet on the outcome. Both answers are results, and
pre-registering one of them would be pre-registering a conclusion.

Run:  python scripts/gate_n7.py
Exit code is 0 only if all three pass.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dst_fit import assess, fit_model, thickening_branch  # noqa: E402
from dst_fit.fitting import synthetic_curve  # noqa: E402
from dst_fit.io_literature import (  # noqa: E402
    GUY_FIG1B_PHI,
    GUY_PMMA_3770,
    RICHARDS_CORNSTARCH,
    read_guy_flow_curves,
    read_richards_series,
)
from dst_fit.models import WYART_CATES  # noqa: E402
from dst_fit.profiling import (  # noqa: E402
    block_bootstrap,
    interval_intersection,
    profile_likelihood,
)

ROOT = Path(__file__).resolve().parents[1]
FIG1B = (ROOT / "Towards a unified description of the rheology of hard-particle suspensions"
         / "Fig 1" / "Fig 1" / "Fig 1b" / "Fig1b.xlsx")
RICHARDS_DIR = (ROOT / "Competing Timescales Lead to Oscillations in Shear-Thickening Suspensions"
                / "Figure1" / "Plotting")

N_TRIALS = 100
COVERAGE_MIN = 80  # G9a, G9b. Same bar as G1b.
MIN_CURVES = 4  # G9c
N_RESAMPLES = 400
LEVEL = 0.95

# The synthetic truth. phi sits close to phi_m so the curve is near jamming,
# which is the regime where the Wald error was suspected of failing. Stresses
# bracket the transition so sigma* is measured rather than bounded.
TRUTH = (0.560, 0.580, 5.0)
SIGMA_GRID = np.logspace(-1, 3, 50)
NOISE = 0.03


def coverage_run() -> dict:
    """Wald, profile and bootstrap coverage of the truth over N_TRIALS curves."""
    truth_phi_m = TRUTH[1]
    hits = {"wald": 0, "profile": 0, "bootstrap": 0}
    widths = {"wald": [], "profile": [], "bootstrap": []}
    asym = {"profile": [], "bootstrap": []}
    n_ok = 0

    for trial in range(N_TRIALS):
        s, eta = synthetic_curve(WYART_CATES, TRUTH, SIGMA_GRID, noise=NOISE, seed=1000 + trial)
        try:
            fit = fit_model(WYART_CATES, s, eta, n_starts=8, seed=trial)
            pr = profile_likelihood(WYART_CATES, s, eta, "phi_m", fit=fit, level=LEVEL)
            bs = block_bootstrap(WYART_CATES, s, eta, "phi_m", fit=fit, level=LEVEL,
                                 n_resamples=150, seed=trial)
        except (RuntimeError, ValueError):
            continue
        n_ok += 1

        i = WYART_CATES.param_names.index("phi_m")
        half = 1.959963984540054 * fit.stderr[i]
        if abs(fit.params[i] - truth_phi_m) <= half:
            hits["wald"] += 1
        widths["wald"].append(2 * half)

        for name, obj in (("profile", pr), ("bootstrap", bs)):
            if obj.contains(truth_phi_m):
                hits[name] += 1
            widths[name].append(obj.width)
            asym[name].append(obj.asymmetry)

    return {
        "n_ok": n_ok,
        "hits": hits,
        "median_width": {k: float(np.median(v)) if v else np.nan for k, v in widths.items()},
        "median_asym": {k: float(np.median(v)) if v else np.nan for k, v in asym.items()},
    }


def real_series(name, sweeps, constants, labels):
    """Wald / profile / bootstrap for phi_m on every curve of a real series."""
    model = WYART_CATES.with_fixed(**constants.as_fixed())
    i = model.param_names.index("phi_m")
    rows = []
    print(f"\n{name}  ({model.describe_fixed()})")
    for s, label in zip(sweeps, labels):
        br = thickening_branch(s.viscosity)
        n_branch = len(range(*br.indices(s.viscosity.size)))
        if n_branch < 10:
            print(f"  {label:>7s}  skipped, {n_branch} branch points")
            continue
        try:
            fit = fit_model(model, s.stress, s.viscosity, branch=br, n_starts=16)
            pr = profile_likelihood(model, s.stress, s.viscosity, "phi_m",
                                    branch=br, fit=fit, level=LEVEL)
            bs = block_bootstrap(model, s.stress, s.viscosity, "phi_m", branch=br,
                                 fit=fit, level=LEVEL, n_resamples=N_RESAMPLES, seed=7)
            iid = block_bootstrap(model, s.stress, s.viscosity, "phi_m", branch=br,
                                  fit=fit, level=LEVEL, n_resamples=N_RESAMPLES,
                                  block_length=1, seed=7)
        except (RuntimeError, ValueError) as exc:
            print(f"  {label:>7s}  refused: {exc}")
            continue
        a = assess(fit)
        rows.append({"label": label, "sweep": s, "fit": fit, "profile": pr,
                     "boot": bs, "iid": iid})
        half = 1.959963984540054 * fit.stderr[i]
        print(
            f"  {label:>7s}  phi_m {fit.params[i]:.4f}   "
            f"Wald [{fit.params[i] - half:.4f}, {fit.params[i] + half:.4f}]   "
            f"profile [{pr.lower:.4f}, {pr.upper:.4f}]   "
            f"block [{bs.lower:.4f}, {bs.upper:.4f}]   "
            f"iid [{iid.lower:.4f}, {iid.upper:.4f}]"
        )
        print(
            f"           {bs.n_runs:3d} sign runs in {bs.n_points} points -> blocks of "
            f"{bs.block_length}; widening profile {pr.widening:.2f}x, "
            f"block {bs.widening:.1f}x, iid {iid.widening:.1f}x; "
            f"runs p = {a.runs_p:.1e}"
        )
    return rows


def widening_needed(rows, key) -> float:
    """How many times wider every interval would have to be to admit one value.

    Computed rather than estimated, because it is the number that decides
    whether the disagreement is worth explaining. The two intervals that set the
    empty intersection are the binding pair; each has to reach half the gap
    further, so the factor is the larger of the two ratios.
    """
    ivs = [r[key] for r in rows]
    lows = np.array([iv.lower for iv in ivs])
    highs = np.array([iv.upper for iv in ivs])
    gap = float(lows.max() - highs.min())
    if gap <= 0:
        return 1.0
    i_lo, i_hi = int(np.argmax(lows)), int(np.argmin(highs))
    factors = []
    for j in (i_lo, i_hi):
        half = 0.5 * (highs[j] - lows[j])
        factors.append((half + 0.5 * gap) / half if half > 0 else np.inf)
    return float(max(factors))


def verdict(rows, key, label):
    if len(rows) < 2:
        print(f"    {label}: too few curves")
        return None
    inter = interval_intersection([r[key] for r in rows], [r["label"] for r in rows])
    print(f"    {label}: {inter.report()}")
    if not inter.non_empty:
        print(f"    {' ' * len(label)}  every interval would have to be "
              f"{widening_needed(rows, key):.1f}x wider to admit one value")
    return inter


def main() -> int:
    failures: list[str] = []
    t0 = time.time()

    # ----------------------------------------------------------------------
    # G9a, G9b: coverage on data the model actually generated
    # ----------------------------------------------------------------------
    print("G9a/G9b: coverage of the truth over "
          f"{N_TRIALS} synthetic curves at phi={TRUTH[0]}, phi_m={TRUTH[1]}, "
          f"sigma*={TRUTH[2]} Pa, {100 * NOISE:.0f}% noise")
    cov = coverage_run()
    print(f"    {cov['n_ok']} of {N_TRIALS} trials completed")
    for name in ("wald", "profile", "bootstrap"):
        w = cov["median_width"][name]
        extra = ""
        if name in cov["median_asym"]:
            extra = f"   median asymmetry {cov['median_asym'][name]:.2f}"
        print(f"      {name:>9s}: covered {cov['hits'][name]:3d}/{cov['n_ok']}   "
              f"median width {w:.4f}{extra}")
    if cov["hits"]["profile"] < COVERAGE_MIN:
        failures.append(f"G9a: profile covered {cov['hits']['profile']}, "
                        f"needed {COVERAGE_MIN}")
    if cov["hits"]["bootstrap"] < COVERAGE_MIN:
        failures.append(f"G9b: bootstrap covered {cov['hits']['bootstrap']}, "
                        f"needed {COVERAGE_MIN}")
    print("    (a method that covered ~100/100 with a huge width would pass this gate "
          "and be useless; the median widths above are what stops that being invisible)")

    # ----------------------------------------------------------------------
    # G9c: the real series
    # ----------------------------------------------------------------------
    guy = [s for s in read_guy_flow_curves(FIG1B, GUY_FIG1B_PHI) if s.sheet != "19%"]
    guy_rows = real_series("Guy et al. 2015 PMMA, d = 3770 nm", guy, GUY_PMMA_3770,
                           [s.sheet for s in guy])

    print("\n  does one phi_m fit all of them?")
    for key, label in (("profile", "profile intervals   "),
                       ("iid", "bootstrap, blocks of 1"),
                       ("boot", "bootstrap, real blocks")):
        verdict(guy_rows, key, label)
    if len(guy_rows) < MIN_CURVES:
        failures.append(f"G9c: only {len(guy_rows)} Guy curves usable, "
                        f"needed {MIN_CURVES}")

    rich = read_richards_series(RICHARDS_DIR)
    rich_rows = real_series("Richards et al. 2019 corn starch", rich, RICHARDS_CORNSTARCH,
                            [f"{s.phi:.2f}" for s in rich])
    print("\n  does one phi_m fit all of them?")
    for key, label in (("profile", "profile intervals   "),
                       ("iid", "bootstrap, blocks of 1"),
                       ("boot", "bootstrap, real blocks")):
        verdict(rich_rows, key, label)

    # ----------------------------------------------------------------------
    # observations
    # ----------------------------------------------------------------------
    print("\nobservations (not gated):")
    for name, rows in (("Guy PMMA", guy_rows), ("Richards corn starch", rich_rows)):
        if not rows:
            continue
        pw = np.array([r["profile"].widening for r in rows])
        bw = np.array([r["boot"].widening for r in rows])
        iw = np.array([r["iid"].widening for r in rows])
        bl = np.array([r["boot"].block_length for r in rows])
        print(f"    {name}: widening over Wald -- profile {pw.min():.2f} to {pw.max():.2f}, "
              f"iid bootstrap {iw.min():.1f} to {iw.max():.1f}, "
              f"block bootstrap {bw.min():.1f} to {bw.max():.1f} "
              f"(blocks of {bl.min():.0f} to {bl.max():.0f})")

    print(f"\n    elapsed {time.time() - t0:.0f} s")
    print()
    if failures:
        print(f"GATE N7 FAIL ({len(failures)} problem(s)):")
        for f in failures:
            print("  " + f)
        return 1
    print("GATE N7 PASS (G9a, G9b, G9c)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
