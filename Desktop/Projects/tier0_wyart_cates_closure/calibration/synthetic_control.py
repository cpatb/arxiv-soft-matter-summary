"""Is the transfer map a property of the data, or of the fitter?

The calibration's central claim is that a Wyart-Cates fit to a real flow curve
reports a volume fraction that is not the one the sample was made at. That claim
is worthless without this control, because a biased optimiser, a bad starting
guess or the branch cut would each produce the same signature.

So: generate flow curves *from the model itself*, on the real stress grids, at
the real volume fractions, with the same fixed constants, add log-normal noise,
and fit them back with the same code. A correct fitter must return slope 1 and
intercept 0. Anything the real data does beyond that is the data talking.

The branch cut is exercised both ways, because it is the most plausible
mechanical explanation for the map: ``thickening_branch`` throws away the
low-stress end, which is the part that pins ``phi`` absolutely.

Run:  python calibration/synthetic_control.py
Exit 0 if the control passes, 1 if the fitter is implicated.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DST_FIT = Path(r"C:\Users\compl\Desktop\Projects\dst_fit")
sys.path.insert(0, str(DST_FIT / "src"))
sys.path.insert(0, str(HERE))

from dst_fit import fit_model, thickening_branch  # noqa: E402
from dst_fit.io_literature import GUY_FIG1B_PHI, GUY_PMMA_3770, read_guy_flow_curves  # noqa: E402
from dst_fit.models import WYART_CATES  # noqa: E402

from calibrate import GUY_DIR, N_STARTS  # noqa: E402
from transfer import transfer_map  # noqa: E402

# Truth for the synthetic material. phi_m is Guy et al.'s own published value and
# sigma_star the median this calibration fits to their d = 3770 nm series, so the
# synthetic curves sit in the same part of parameter space as the real ones.
PHI_M_TRUE = 0.558
SIGMA_STAR_TRUE = 9.91  # Pa

# Noise levels, as standard deviations of log viscosity. 0.03 brackets the
# scatter of the deposited curves about their own fits (rmse_log 0.002 to 0.06).
NOISE_LEVELS = (0.0, 0.01, 0.03, 0.10)

# How far the control is allowed to stray before the fitter is implicated. The
# real data gives slope 1.24 and intercept -0.144; these bars are far inside that
# and are fixed here rather than chosen after looking.
SLOPE_TOL = 0.10
INTERCEPT_TOL = 0.05

SEED = 20260830


def synthetic_series(grids, model, noise, rng, phi_m=PHI_M_TRUE,
                     sigma_star=SIGMA_STAR_TRUE, use_branch=True,
                     noise_per_sheet=None):
    """Fit model-generated curves back and return (true phi, fitted phi).

    ``noise_per_sheet`` overrides ``noise`` curve by curve, so the control can be
    run at each real curve's own measured residual scatter rather than at a
    round number. That matched run is the one that decides the question; the
    fixed-``noise`` runs only map out where the mechanism bites.
    """
    true, fitted = [], []
    for sheet, sigma in grids.items():
        phi = GUY_FIG1B_PHI[sheet]
        if phi >= phi_m:
            continue                      # jammed; no flowing curve to generate
        eta = model(sigma, phi, phi_m, sigma_star)
        ok = np.isfinite(eta) & (eta > 0)
        if ok.sum() < 5:
            continue
        s, e = sigma[ok], eta[ok]
        sigma_log = noise if noise_per_sheet is None else noise_per_sheet.get(sheet, noise)
        if sigma_log > 0:
            e = e * np.exp(rng.normal(0.0, sigma_log, e.size))
        branch = thickening_branch(e) if use_branch else slice(0, e.size)
        try:
            fit = fit_model(model, s, e, branch=branch, n_starts=N_STARTS)
        except (RuntimeError, ValueError):
            continue
        true.append(phi)
        fitted.append(float(fit.params[0]))
    return np.array(true), np.array(fitted)


def main() -> int:
    model = WYART_CATES.with_fixed(**GUY_PMMA_3770.as_fixed())
    real = read_guy_flow_curves(GUY_DIR / "Fig 1b" / "Fig1b.xlsx", GUY_FIG1B_PHI)
    grids = {s.sheet: s.stress for s in real}
    rng = np.random.default_rng(SEED)

    print("Synthetic control: curves generated from the model on the real stress grids,")
    print("fitted back with the same constants. A correct fitter returns (1.0000, 0.0000).")
    print(f"Tolerances fixed before the run: |slope - 1| < {SLOPE_TOL}, "
          f"|intercept| < {INTERCEPT_TOL}\n")

    failures = []
    for use_branch in (True, False):
        tag = "thickening branch cut" if use_branch else "whole curve fitted"
        print(f"  {tag}")
        for noise in NOISE_LEVELS:
            t, f = synthetic_series(grids, model, noise, rng, use_branch=use_branch)
            if t.size < 3:
                print(f"      noise {noise:.2f}: only {t.size} curves, skipped")
                continue
            tm = transfer_map(t, f, n_boot=500, seed=SEED)
            bad = (abs(tm.slope - 1.0) > SLOPE_TOL
                   or abs(tm.intercept) > INTERCEPT_TOL)
            if bad:
                failures.append((tag, noise, tm.slope, tm.intercept))
            print(f"      noise {noise:.2f}: slope {tm.slope:.4f}  "
                  f"intercept {tm.intercept:+.4f}  r2 {tm.r2:.5f}  n={tm.n_curves}"
                  f"{'   <- OUT OF TOLERANCE' if bad else ''}")
        print()

    # the real series, and the matched control that decides the question
    rows, per_sheet = [], {}
    for s in real:
        if s.phi is None:
            continue
        try:
            fit = fit_model(model, s.stress, s.viscosity,
                            branch=thickening_branch(s.viscosity), n_starts=N_STARTS)
        except (RuntimeError, ValueError):
            continue
        if fit.parameters_at_bound():
            continue
        rows.append((s.phi, float(fit.params[0])))
        per_sheet[s.sheet] = float(fit.rmse_log)
    tm = transfer_map([r[0] for r in rows], [r[1] for r in rows], seed=SEED)
    print(f"  real Guy d = 3770 nm: slope {tm.slope:.4f}  intercept {tm.intercept:+.4f}"
          f"  r2 {tm.r2:.5f}  n={tm.n_curves}")
    print("    residual scatter of those fits (rmse of log eta), per curve: "
          + ", ".join(f"{k} {v:.3f}" for k, v in per_sheet.items()))

    print("\n  MATCHED CONTROL: synthetic curves noised at each real curve's own scatter,")
    print("  repeated over independent noise draws. This is the run that decides it.")
    print("  Run at two truths, because a synthetic curve can only be generated at")
    print("  phi < phi_m: at the published phi_m = 0.558 the two highest-phi curves")
    print("  cannot be generated at all, and those are the two noisiest, so the run is")
    print("  repeated at phi_m = 0.60 where all six exist.")
    verdicts = []
    for phi_m in (PHI_M_TRUE, 0.60):
        slopes, intercepts, sizes = [], [], []
        for trial in range(20):
            r = np.random.default_rng(SEED + trial)
            t, f = synthetic_series(grids, model, 0.0, r, phi_m=phi_m,
                                    noise_per_sheet=per_sheet)
            if t.size < 3:
                continue
            m = transfer_map(t, f, n_boot=1, seed=SEED)
            slopes.append(m.slope)
            intercepts.append(m.intercept)
            sizes.append(t.size)
        if not slopes:
            print(f"      phi_m = {phi_m}: no usable draws")
            continue
        slopes, intercepts = np.array(slopes), np.array(intercepts)
        lo_s, hi_s = np.percentile(slopes, [2.5, 97.5])
        lo_i, hi_i = np.percentile(intercepts, [2.5, 97.5])
        outside = tm.slope > hi_s or tm.intercept < lo_i
        verdicts.append(outside)
        print(f"\n      truth phi_m = {phi_m:.3f}, {int(np.median(sizes))} of "
              f"{len(grids)} curves generatable")
        print(f"        slope     median {np.median(slopes):.4f}, 95% of draws "
              f"[{lo_s:.4f}, {hi_s:.4f}]")
        print(f"        intercept median {np.median(intercepts):+.4f}, 95% of draws "
              f"[{lo_i:+.4f}, {hi_i:+.4f}]")
        print(f"        reproduces {(np.median(slopes) - 1.0) / (tm.slope - 1.0):.0%} "
              f"of the real slope excess; real is "
              f"{'OUTSIDE' if outside else 'INSIDE'} this spread")
    real_outside = bool(verdicts) and all(verdicts)

    if failures:
        print("\nFixed-noise runs out of tolerance (where the mechanism bites):")
        for tag, noise, sl, ic in failures:
            print(f"  {tag}, noise {noise:.2f}: slope {sl:.4f}, intercept {ic:+.4f}")
    return 0 if real_outside else 1


if __name__ == "__main__":
    raise SystemExit(main())
