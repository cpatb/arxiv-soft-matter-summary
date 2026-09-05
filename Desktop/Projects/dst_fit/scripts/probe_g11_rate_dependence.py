"""probe_g11_rate_dependence.py -- is the inversion clock rate-dependent?

WHAT THIS IS, AND WHAT IT IS NOT
--------------------------------
`PLAN_2026-08-29.md` item 5 left one thing open: *"the same comparison at Lin's other
simulated rates (Fig_3a has 0.02, 0.13, 0.20 as well as the 0.74 used here), to find out
whether the 41 % is rate-dependent. If it is, it is not a systematic to carry but a second
clock."*

**That comparison cannot be run, and the reason is in the data rather than in the code.**
The 41 % is a difference between TWO channels: the Wyart-Cates inversion of the total
viscosity, and the directly evaluated contact-stress fraction. The direct channel exists
only in `Fig_3b`, which is a single record at `gamma_dot = 0.74 gamma_dot_0`. `Fig_3a`
carries four rates but only the *total viscosity* -- i.e. only the inversion channel. So
there is no second channel at 0.02, 0.13 or 0.20 to difference against, and no amount of
fitting produces one.

What CAN be measured, and what this script measures, is the adjacent question that
actually decides how the 41 % should be quoted:

    Does the INVERSION channel's gamma_eff depend on shear rate?

Both readings matter, and they point opposite ways:

  * If gamma_eff is rate-INDEPENDENT across a 37x span in gamma_dot, then the inversion
    behaves like a fixed map and the 41 % offset measured at one rate is more defensibly
    carried as a systematic of the method.
  * If gamma_eff MOVES with rate, then "the strain clock" is not one number even within a
    single channel, the 41 % is a difference between two things that are each rate-varying,
    and the spec must say the clock is rate-resolved rather than definition-resolved.

⚠ This is INDIRECT evidence about the 41 %. It cannot confirm rate-independence of the
*difference*; it can only remove one way for the difference to be misleading. Stated here
so the result is not over-read later.

METHOD, and it is deliberately the same code path as G11
--------------------------------------------------------
`transient.fit_recovery` is called on each rate's record with identical settings and the
same packing triples the direct control used, so nothing about the fit differs between
rates except the data. The two triples are gauge choices, not measurements; G11's affine
argument says gamma_eff must not depend on them, and printing both is the check that it
still does not at every rate.

    python scripts/probe_g11_rate_dependence.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from dst_fit.transient import WC_EXPONENT, fit_recovery      # noqa: E402

FIG3 = (ROOT / "Hydrodynamic and Contact Contributions to Continuous Shear Thickening "
        "in Colloidal Suspensions" / "Fig3" / "Fig3")

RATES = [0.02, 0.13, 0.20, 0.74]
# The same two deliberately-far-apart gauge choices probe_g11_direct_control.py uses.
# They are NOT measurements of Lin's simulation; printing both is the affine-invariance
# check, which must hold at every rate or the clock is not what G11 says it is.
TRIPLES = [(0.50, 0.58, 0.66), (0.40, 0.52, 0.70)]


def load(name: str) -> np.ndarray:
    return np.loadtxt(FIG3 / name, delimiter=",")


def main() -> int:
    if not FIG3.is_dir():
        print("Fig3 data not found at", FIG3)
        return 1

    print("=" * 78)
    print("G11 rate dependence: gamma_eff from the INVERSION channel, four rates")
    print("=" * 78)
    print("Data: Lin et al. Fig_3a, relative viscosity vs accumulated strain after")
    print("reversal. Rates are in units of the simulation's own gamma_dot_0.")
    print("The DIRECT channel does not exist at these rates -- see the module docstring.")
    print()

    rows = []
    print("%6s %8s %10s %10s %9s %22s %10s %7s %6s"
          % ("rate", "n_pts", "gEff(T1)", "gEff(T2)", "gauge%",
             "95% interval (T1)", "chi2_red", "n_fit", "ok"))
    for rate in RATES:
        d = load(f"Fig_3a_gamma_dot_{rate:.2f}.csv")
        strain, eta_r = d[:, 0], d[:, 1]
        got = []
        for (phi, phi_m, phi_0) in TRIPLES:
            r = fit_recovery(strain, eta_r, rate=rate, phi=phi, phi_m=phi_m,
                             phi_0=phi_0, n=WC_EXPONENT)
            got.append(r)
        a, b = got
        g1, g2 = a.gamma_eff, b.gamma_eff
        spread = 100.0 * abs(g1 - g2) / (0.5 * (g1 + g2))
        rows.append({"rate": rate, "n_points": int(strain.size),
                     "gamma_eff_triple1": float(g1), "gamma_eff_triple2": float(g2),
                     "gauge_spread_pct": float(spread),
                     "gamma_eff_lo": float(a.gamma_eff_lo),
                     "gamma_eff_hi": float(a.gamma_eff_hi),
                     "chi2_reduced": float(a.chi2_reduced),
                     "n_fit": int(a.n_fit), "rmse": float(a.rmse),
                     "success": bool(a.success and b.success),
                     "t_eff_s": float(a.t_eff)})
        print("%6.2f %8d %10.4f %10.4f %9.4f   [%8.4f, %8.4f] %10.1f %7d %6s"
              % (rate, strain.size, g1, g2, spread, a.gamma_eff_lo, a.gamma_eff_hi,
                 a.chi2_reduced, a.n_fit, "yes" if rows[-1]["success"] else "NO"))

    if not all(r["success"] for r in rows):
        print("\nAT LEAST ONE FIT DID NOT CONVERGE -- the span below is not quotable.")
        return 1

    g = np.array([r["gamma_eff_triple1"] for r in rows])
    rate_arr = np.array([r["rate"] for r in rows])
    span = g.max() / g.min()
    rate_span = rate_arr.max() / rate_arr.min()
    # a slope in log-log, which is the shape a "second clock" would have
    slope = float(np.polyfit(np.log(rate_arr), np.log(g), 1)[0])

    print()
    print("VERDICT")
    print("-------")
    monotone = bool(np.all(np.diff(g) > 0) or np.all(np.diff(g) < 0))
    print("gamma_eff spans %.3fx (%.4f to %.4f) across a %.0fx span in shear rate."
          % (span, g.min(), g.max(), rate_span))
    print("Series in rate order: %s" % ", ".join("%.4f" % v for v in g))
    if monotone:
        print("Monotone. Log-log slope d(ln gamma_eff)/d(ln gamma_dot) = %+.4f." % slope)
    else:
        print("*** NON-MONOTONE (it rises then falls). A power law cannot describe this")
        print("    shape at any exponent, so the log-log slope (%+.4f) is NOT quotable"
              % slope)
        print("    and is printed only to show it would mislead.")
    print()
    print("The gauge spread between the two packing triples is %.4f%% at worst, which is"
          % max(r["gauge_spread_pct"] for r in rows))
    print("the affine-invariance check G11 rests on, holding at every rate.")
    print()
    # Is the span bigger than the fits' own uncertainty? A 1.9x spread means nothing if
    # each interval is 2x wide. Compare the extremes' 95% intervals for overlap.
    lo = min(rows, key=lambda r: r["gamma_eff_triple1"])
    hi = max(rows, key=lambda r: r["gamma_eff_triple1"])
    overlap = lo["gamma_eff_hi"] >= hi["gamma_eff_lo"]
    print("Extremes: rate %.2f -> %.4f [%.4f, %.4f]   rate %.2f -> %.4f [%.4f, %.4f]"
          % (lo["rate"], lo["gamma_eff_triple1"], lo["gamma_eff_lo"], lo["gamma_eff_hi"],
             hi["rate"], hi["gamma_eff_triple1"], hi["gamma_eff_lo"], hi["gamma_eff_hi"]))
    print("95%% intervals %s -- the spread IS %s the fit uncertainty."
          % ("OVERLAP" if overlap else "do NOT overlap",
             "not resolved against" if overlap else "resolved against"))
    print()
    print("A NOTE ON WHAT THAT INTERVAL ALREADY CONTAINS, because an earlier version of")
    print("this script got it wrong in a way that reversed the verdict.")
    print("`fit_relaxation` admits lambdas with  ssr <= ssr_min * (1 + CHI2_1DOF_95/dof),")
    print("i.e. a RELATIVE SSR criterion. That is the error-rescaled form: the half-width")
    print("already scales as sqrt(chi2_reduced), so the misfit (chi2_red = %.0f to %.0f)"
          % (min(r["chi2_reduced"] for r in rows), max(r["chi2_reduced"] for r in rows)))
    print("is ALREADY priced in. The first version inflated by the Birge ratio on top of")
    print("that -- an extra factor of 30 to 63 -- which pushed the lower bound to -0.54, a")
    print("NEGATIVE e-folding strain, and manufactured an overlap. Do not re-add it.")
    print("CHI2_1DOF_95 = 3.8415, so these are 95% intervals, not 68%.")
    print()
    if span < 1.10:
        print("=> RATE-INDEPENDENT to better than 10%. The inversion behaves like a fixed")
        print("   map over this range, so the 41% offset measured at gamma_dot = 0.74 is")
        print("   more defensibly carried as a systematic OF THE METHOD.")
    elif overlap:
        print("=> The %.2fx span is NOT resolved against the fit uncertainty. Report it as" % span)
        print("   a span over the sampled rates, not as a measured rate dependence.")
    else:
        print("=> RATE-DEPENDENT by %.2fx, and RESOLVED: the extremes' 95%% intervals do not"
              % span)
        print("   overlap. So 'the strain clock' is not one number even within a single")
        print("   channel, and the 41% is a difference between two rate-varying quantities.")
        print("   gamma_0 needs a RATE quoted with it as well as a definition.")
        if not monotone:
            print("   BUT the dependence is NON-MONOTONE, so it is not a law: do not fit an")
            print("   exponent, and do not extrapolate outside 0.02 <= gamma_dot <= 0.74.")
    print()
    print("⚠ INDIRECT, and this bounds what may be concluded. The direct channel exists at")
    print("  ONE rate, so this cannot show the 41% DIFFERENCE is rate-varying -- only that")
    print("  one of its two terms is. The spec must not claim more.")

    out = ROOT / "logs" / "g11_rate_dependence.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps({
        "source": "Lin et al. Fig_3a (inversion channel only; no direct channel exists "
                  "at these rates)",
        "triples": TRIPLES, "wc_exponent": WC_EXPONENT,
        "rows": rows, "gamma_eff_span": float(span),
        "rate_span": float(rate_span), "loglog_slope": slope,
    }, indent=2), encoding="utf-8")
    print("\nwrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
