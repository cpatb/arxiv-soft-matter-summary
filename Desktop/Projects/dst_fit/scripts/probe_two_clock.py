"""Is the falling lambda real, or a two-component mixture with rate-independent clocks?

Four models, all fitted to the SAME binned, weighted u(gamma) on the recovery
branch of all four rates at once.  Amplitudes are always eliminated exactly
(VARPRO), so every model is a global optimum on its own grid.

  M1  per-rate single exponential, per-rate lambda      3 x 4 = 12 par
  M2  per-rate single exponential, ONE shared lambda    2 x 4 + 1 =  9 par
  M3  per-rate TWO exponentials, TWO shared lambdas     3 x 4 + 2 = 14 par
  M4  per-rate two exponentials, shared slow lambda,
      per-rate fast lambda                              3 x 4 + 1 + 3 = 16 par

M3 is the reviewer's claim.  Note it has MORE parameters than M1, so "fits
better" is not by itself evidence: the comparison has to be penalised, and it
has to be penalised on the EFFECTIVE sample size, not the bin count.

Provenance: written in the Drop Impact session on 2026-08-28 as a review of gate
G11, and copied here so the model comparison quoted in STRAIN_CLOCK_2026-08-28.md
and NEXT.md is reproducible rather than cited from a message. Only the root path
was changed.

A caveat to carry with the output: M3's two shared clocks come out close together
and with opposite-sign amplitudes, so that fit is partly describing curvature
rather than two cleanly separated physical processes. A second fit on a different
window gave 0.492 and 1.909 for the same pair. Cite the model *comparison* and the
surviving rate trend in M4's fast clock; do not present the two lambdas as measured
hydrodynamic and contact clocks.

Run:  python scripts/probe_two_clock.py     (needs the gitignored Lin deposit)
Prints a table; not a gate.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dst_fit.transient import log_bin, reduced_jamming_distance          # noqa: E402
from dst_fit.transient import d_reduced_jamming_distance                 # noqa: E402
from dst_fit.profiling import block_length_from_runs                     # noqa: E402

D = ROOT / ("Hydrodynamic and Contact Contributions to Continuous Shear "
            "Thickening in Colloidal Suspensions") / "Fig2" / "Fig2"
RATES = [0.5, 1.0, 2.0, 5.0]
NBINS = 300
N = 2.0


def prep(rate):
    a = np.loadtxt(D / f"Fig_2b_gamma_dot_{rate}.csv", delimiter=",", usecols=(0, 1))
    g, e, err = log_bin(a[:, 0], a[:, 1], NBINS)
    i0 = int(np.argmin(e))
    g, e, err = g[i0:], e[i0:], err[i0:]
    u = reduced_jamming_distance(e, n=N)
    su = np.abs(d_reduced_jamming_distance(e, n=N)) * err
    return g, u, su


DATA = [prep(r) for r in RATES]
NTOT = sum(len(g) for g, _, _ in DATA)


def ssr_for(lams_per_rate):
    """lams_per_rate[i] is a LIST of decay strains for record i. Amplitudes exact."""
    total, resid_all = 0.0, []
    for (g, u, su), lams in zip(DATA, lams_per_rate):
        w = 1.0 / su
        cols = [np.ones_like(g)] + [np.exp(-(g - g[0]) / L) for L in lams]
        X = np.column_stack(cols) * w[:, None]
        coef, *_ = np.linalg.lstsq(X, u * w, rcond=None)
        r = u * w - X @ coef
        total += float(r @ r)
        resid_all.append(r)
    return total, resid_all


GRID = np.logspace(np.log10(0.02), np.log10(20.0), 260)

# ---- M1: per-rate single lambda -------------------------------------------
m1_lams, m1_ssr = [], 0.0
for i in range(4):
    g, u, su = DATA[i]
    w = 1.0 / su
    vals = []
    for L in GRID:
        X = np.column_stack([np.ones_like(g), np.exp(-(g - g[0]) / L)]) * w[:, None]
        c, *_ = np.linalg.lstsq(X, u * w, rcond=None)
        r = u * w - X @ c
        vals.append(float(r @ r))
    k = int(np.argmin(vals))
    m1_lams.append(GRID[k])
    m1_ssr += vals[k]
print("M1 per-rate single exponential")
print("   lambda =", np.round(m1_lams, 4), f"  SSR = {m1_ssr:.4f}")
print(f"   spread {max(m1_lams)/min(m1_lams):.2f}x, "
      f"{'DECREASING' if m1_lams[-1] < m1_lams[0] else 'INCREASING'}")

# ---- M2: one shared lambda -------------------------------------------------
v2 = [ssr_for([[L]] * 4)[0] for L in GRID]
k2 = int(np.argmin(v2))
print(f"\nM2 one shared lambda = {GRID[k2]:.4f}   SSR = {v2[k2]:.4f}")

# ---- M3: two shared lambdas ------------------------------------------------
best3 = (np.inf, None)
for ia, La in enumerate(GRID):
    for Lb in GRID[ia + 1:]:
        s, _ = ssr_for([[La, Lb]] * 4)
        if s < best3[0]:
            best3 = (s, (La, Lb))
print(f"M3 two shared lambdas = {best3[1][0]:.4f}, {best3[1][1]:.4f}   SSR = {best3[0]:.4f}")

# ---- M4: shared slow, per-rate fast ---------------------------------------
Lslow = best3[1][1]
m4_fast, m4_ssr = [], 0.0
for i, (g, u, su) in enumerate(DATA):
    w = 1.0 / su
    vals = []
    for L in GRID:
        if L >= Lslow:
            vals.append(np.inf)
            continue
        X = np.column_stack([np.ones_like(g), np.exp(-(g - g[0]) / L),
                             np.exp(-(g - g[0]) / Lslow)]) * w[:, None]
        c, *_ = np.linalg.lstsq(X, u * w, rcond=None)
        r = u * w - X @ c
        vals.append(float(r @ r))
    k = int(np.argmin(vals))
    m4_fast.append(GRID[k])
    m4_ssr += vals[k]
print(f"M4 shared slow {Lslow:.4f}, per-rate fast =", np.round(m4_fast, 4),
      f"  SSR = {m4_ssr:.4f}")
print(f"   fast spread {max(m4_fast)/min(m4_fast):.2f}x, "
      f"{'DECREASING' if m4_fast[-1] < m4_fast[0] else 'INCREASING'}")

# ---- effective sample size, from the sign runs of the best model's residuals
_, res3 = ssr_for([[best3[1][0], best3[1][1]]] * 4)
neff = 0
for r in res3:
    bl, nruns = block_length_from_runs(r)
    neff += nruns
print(f"\neffective sample size from residual sign runs (M3): {neff} "
      f"against {NTOT} bins")


def aic(ssr, k, n):
    return n * np.log(ssr / n) + 2 * k + 2 * k * (k + 1) / max(n - k - 1, 1)


print(f"\n{'model':>4} {'k':>4} {'SSR':>10} {'AICc(n=bins)':>14} {'AICc(n=neff)':>14}")
for nm, ssr, k in [("M1", m1_ssr, 12), ("M2", v2[k2], 9),
                   ("M3", best3[0], 14), ("M4", m4_ssr, 16)]:
    print(f"{nm:>4} {k:>4} {ssr:10.4f} {aic(ssr, k, NTOT):14.2f} "
          f"{aic(ssr*neff/NTOT, k, neff):14.2f}")

print("\nAmplitudes of M3, to test the reviewer's 'fast share tracks f_ss' claim:")
La, Lb = best3[1]
print(f"{'gdot':>6} {'offset':>9} {'amp_fast':>10} {'amp_slow':>10} "
      f"{'fast share':>11} {'f_ss(WC)':>9}")
tot = dict(np.loadtxt(D / "Fig_2a_steady_state.csv", delimiter=",", usecols=(0, 1)))
for r, (g, u, su) in zip(RATES, DATA):
    w = 1.0 / su
    X = np.column_stack([np.ones_like(g), np.exp(-(g - g[0]) / La),
                         np.exp(-(g - g[0]) / Lb)]) * w[:, None]
    c, *_ = np.linalg.lstsq(X, u * w, rcond=None)
    share = abs(c[1]) / (abs(c[1]) + abs(c[2])) if (abs(c[1]) + abs(c[2])) > 0 else np.nan
    fss = np.exp(-8.699 / (tot[r] * r))
    print(f"{r:6.1f} {c[0]:9.4f} {c[1]:10.4f} {c[2]:10.4f} {share:11.3f} {fss:9.3f}")


# ---- which model actually wins, and on which sample size -------------------
# The AICc(n=bins) column above is printed only to show how large the mistake
# would be. Bins are not independent draws -- the residuals run in long same-sign
# blocks -- so an AICc computed on 570 of them badly overstates the evidence for
# extra parameters. The neff column is the one to read.
_rows = [("M1", m1_ssr, 12), ("M2", v2[k2], 9), ("M3", best3[0], 14), ("M4", m4_ssr, 16)]
_scored = sorted((aic(ssr * neff / NTOT, k, neff), nm) for nm, ssr, k in _rows)
print()
print(f"   -> on the effective sample size, {_scored[0][1]} wins, beating "
      f"{_scored[1][1]} by {_scored[1][0] - _scored[0][0]:.1f} AICc.")
print("   -> read the AICc(n=neff) column. AICc(n=bins) is printed only to show the")
print("      size of that error and must not be used to choose a model.")
