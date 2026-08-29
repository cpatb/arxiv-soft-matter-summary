"""Recovery strain after shear reversal, Lin et al. Fig2b (2 um silica).

Fix from the first pass: the GLOBAL argmin lands on a single-point dip at the
very end of the 2.0 and 5.0 records. Restrict the minimum to the interior
(gamma < 2), which is where the reversal dip actually sits (gamma ~ 0.2-0.5).

Question: does the min->steady recovery collapse on accumulated STRAIN
(the gamma_0 term of Richards Eq. 4) or on TIME (the 1/t_r term)?
Collapse metric: spread of the 63% recovery landmark across the four rates,
in strain vs in time. Smaller spread = the clock that governs.
"""
import numpy as np
from pathlib import Path

D = Path(r"C:\Users\compl\Desktop\Projects\dst_fit"
         r"\Hydrodynamic and Contact Contributions to Continuous Shear Thickening"
         r" in Colloidal Suspensions\Fig2\Fig2")
RATES = [0.5, 1.0, 2.0, 5.0]
SMOOTH = 101          # points; records are ~10-25k over 3 decades of strain
INTERIOR_MAX = 2.0    # strain below which the reversal minimum must lie

g63s, t63s = [], []
print(f"{'gdot':>5} {'gam_min':>8} {'eta_min':>8} {'eta_ss':>7} "
      f"{'gam_63':>7} {'t_63 [s]':>9}")
for r in RATES:
    a = np.loadtxt(D / f"Fig_2b_gamma_dot_{r}.csv", delimiter=",")
    g, eta = a[:, 0], a[:, 1]
    # boxcar smooth to defeat single-point noise before landmarking
    k = np.ones(SMOOTH) / SMOOTH
    es = np.convolve(eta, k, mode="same")
    es[:SMOOTH], es[-SMOOTH:] = eta[:SMOOTH], eta[-SMOOTH:]

    interior = np.where(g < INTERIOR_MAX)[0]
    imin = interior[int(np.argmin(es[interior]))]
    ss = np.median(eta[g >= g.max() * 0.9])
    target = es[imin] + (1 - np.exp(-1.0)) * (ss - es[imin])
    after = np.arange(imin, len(es))
    hit = after[es[after] >= target]
    assert hit.size, f"no 63% crossing at gdot={r}"
    g63 = g[hit[0]] - g[imin]
    g63s.append(g63)
    t63s.append(g63 / r)
    print(f"{r:5.1f} {g[imin]:8.4f} {es[imin]:8.3f} {ss:7.3f} "
          f"{g63:7.4f} {g63/r:9.4f}")

g63s, t63s = np.array(g63s), np.array(t63s)
print(f"\ncollapse on STRAIN: gamma_63 = {g63s.mean():.3f} +/- {g63s.std(ddof=1):.3f}"
      f"   spread max/min = {g63s.max()/g63s.min():.2f}x")
print(f"collapse on TIME  : t_63     = {t63s.mean():.3f} +/- {t63s.std(ddof=1):.3f} s"
      f"   spread max/min = {t63s.max()/t63s.min():.2f}x")
print(f"\n(rate spans {max(RATES)/min(RATES):.0f}x, so a pure time clock would give "
      f"a {max(RATES)/min(RATES):.0f}x strain spread and a flat time spread.)")
