"""What does a calibrated residual bootstrap read against Wald, and why is it not 1?

Gate N7 turns on a comparison between two widening factors -- how much wider the
profile interval is than Wald, and how much wider the bootstrap interval is. The
two were measured against two *different* Wald intervals: the profile against
``2 sqrt(F(1, dof)) se``, the bootstrap against the asymptotic ``2 x 1.95996 se``,
differing by 4.9% at dof 26 and 31% at dof 5. Both now go through
:func:`dst_fit.profiling.wald_width` on the finite-sample ruler.

Unifying the ruler exposes a second question that had been hidden inside it:
**what should a correctly calibrated percentile bootstrap read?** Not 1.00. On
data drawn from the model itself with iid noise -- where the bootstrap's own
assumption is true by construction -- two exact factors separate the two
recipes, and both are arithmetic, not noise:

  z/t          a percentile interval cuts at the sampling distribution's own
               quantiles, which under normality are +/- z sd; the Wald interval
               cuts at +/- t_{dof} se.  z/t = 0.955 at dof 27.
  sqrt(1-k/n)  fitted residuals are shrunk -- E[SSR] = sigma^2 (n-k) -- so a
               resample drawn from them is under-dispersed by sqrt((n-k)/n),
               0.949 at n = 30, k = 3, unless they are rescaled.

Their product, 0.906, is what the bootstrap read before the rescaling was added.
With the rescaling, the prediction is z/t alone.

Two hypotheses are tested here, one of which was in the session plan and is
**wrong**:

  H_cold  the shortfall came from refitting each resample single-start from the
          generic heuristic instead of warm-starting from the full-data optimum.
  H_arith the shortfall is the two arithmetic factors above and nothing else.

Pre-registered, before the run:

  C1  after rescaling, the median widening lands within 0.05 of z/t.
  C2  the *cold*-start and warm-start medians, on the same seeds, differ by less
      than 0.02 -- which would REFUTE H_cold. If they differ by more, H_cold has
      support and C1's interpretation has to be revisited.
  C3  the ratio of the bootstrap distribution's sd to the Wald se lands within
      0.05 of 1.00 after rescaling. This is the factor-free version of C1: it
      removes z/t from the comparison, so it tests the dispersion alone.
  C4  only realisations whose fitted parameter is interior to its bounds count.
      A fit on a bound is not a measurement and its stderr is not a width.

Run:  python scripts/probe_bootstrap_calibration.py
Exit code 0 only if C1, C2 and C3 all hold.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.stats import f as f_dist
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dst_fit.fitting import fit_model, synthetic_curve  # noqa: E402
from dst_fit.models import get_model  # noqa: E402
from dst_fit.profiling import block_bootstrap  # noqa: E402

# Wyart-Cates at a packing well clear of jamming, so phi_m is interior and the
# cost surface is not the pathological one. Truth values are the Richards
# corn-starch scale, so the probe exercises the regime the gates run in.
MODEL = "wyart_cates"
TRUTH = (0.44, 0.52, 3.0)  # phi, phi_m, sigma_star
SIGMA = np.logspace(-0.5, 2.0, 30)
NOISE = 0.03  # sd in log viscosity, ~3%, rheometer reproducibility
N_REAL = 40
PARAM = "phi_m"
N_RESAMPLES = 300
LEVEL = 0.95


def one(seed: int, *, warm: bool, rescale: bool):
    """Return (widening, sd_boot / se_wald) or None if the realisation is unusable."""
    sig, eta = synthetic_curve(MODEL, TRUTH, SIGMA, noise=NOISE, seed=seed)
    if sig.size < 12:
        return None
    try:
        fit = fit_model(MODEL, sig, eta, n_starts=12, seed=0)
    except (RuntimeError, ValueError):
        return None
    mdl = get_model(MODEL)
    i = mdl.param_names.index(PARAM)
    lo, hi = mdl.bounds[0][i], mdl.bounds[1][i]
    if not (lo * 1.001 < fit.params[i] < hi * 0.999):  # C4
        return None
    se = float(fit.stderr[i])
    if not np.isfinite(se) or se <= 0:
        return None
    try:
        bs = block_bootstrap(
            MODEL, sig, eta, PARAM,
            fit=fit, block_length=1, n_resamples=N_RESAMPLES, seed=seed,
            boot_starts=(3 if warm else 1), warm_start=warm,
            rescale_residuals=rescale, level=LEVEL,
        )
    except (RuntimeError, ValueError):
        return None
    return bs.widening, float(np.std(bs.values, ddof=1)) / se, bs.dof


def collect(*, warm: bool, rescale: bool):
    w, r, dofs = [], [], []
    for s in range(N_REAL):
        got = one(s, warm=warm, rescale=rescale)
        if got is not None:
            w.append(got[0])
            r.append(got[1])
            dofs.append(got[2])
    return np.asarray(w), np.asarray(r), dofs


def main() -> int:
    warm_w, warm_r, dofs = collect(warm=True, rescale=True)
    cold_w, _, _ = collect(warm=False, rescale=True)
    raw_w, raw_r, _ = collect(warm=True, rescale=False)

    if warm_w.size < 20:
        print(f"only {warm_w.size} usable realisations; the probe cannot decide")
        return 2

    dof = int(np.median(dofs))
    t = float(np.sqrt(f_dist.ppf(LEVEL, 1, dof)))
    z = float(norm.ppf(0.5 * (1 + LEVEL)))
    n = int(SIGMA.size)
    k = get_model(MODEL).n_params
    zt = z / t
    shrink = float(np.sqrt((n - k) / n))

    print(f"{warm_w.size}/{N_REAL} usable realisations, {PARAM}, iid residual "
          f"bootstrap vs Wald on the shared finite-sample ruler, dof = {dof}\n")
    print(f"  predicted from arithmetic alone")
    print(f"    z/t                        {zt:.4f}   (z = {z:.4f}, t = {t:.4f})")
    print(f"    residual shrinkage         {shrink:.4f}   sqrt((n-k)/n), n = {n}, k = {k}")
    print(f"    product, no rescaling      {zt * shrink:.4f}")
    print(f"    with rescaling             {zt:.4f}\n")
    print(f"  measured, median over realisations")
    print(f"    warm start, rescaled       {np.median(warm_w):.4f}   "
          f"[{np.quantile(warm_w, .25):.3f}, {np.quantile(warm_w, .75):.3f}] IQR")
    print(f"    COLD start, rescaled       {np.median(cold_w):.4f}")
    print(f"    warm start, NOT rescaled   {np.median(raw_w):.4f}\n")
    print(f"    sd(bootstrap) / se(Wald)   {np.median(warm_r):.4f}   "
          f"(not rescaled: {np.median(raw_r):.4f})")

    c1 = abs(float(np.median(warm_w)) - zt) < 0.05
    c2 = abs(float(np.median(warm_w)) - float(np.median(cold_w))) < 0.02
    c3 = abs(float(np.median(warm_r)) - 1.0) < 0.05
    print()
    print(f"  C1  median widening within 0.05 of z/t          {'PASS' if c1 else 'FAIL'}")
    print(f"  C2  cold == warm to 0.02  (H_cold REFUTED)      {'PASS' if c2 else 'FAIL'}")
    print(f"  C3  sd(boot)/se(Wald) within 0.05 of 1          {'PASS' if c3 else 'FAIL'}")
    if c2:
        print("\n  H_cold is refuted: the cold start was not the cause of the shortfall.")
        print("  The session plan asserted it was. It is not; the arithmetic is.")
    return 0 if (c1 and c2 and c3) else 1


if __name__ == "__main__":
    raise SystemExit(main())
