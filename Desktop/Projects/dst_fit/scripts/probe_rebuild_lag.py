"""Probe: does the frictional network rebuild as a first-order relaxation from the
reversal minimum (Richards Eq. 4), or after a strain LAG?

Eq. 4, df/dgamma = -(f - f_ss)/gamma_0 at constant rate, makes the recovery steepest
where it starts. Figure 2b of the strain-clock note shows the opposite on Lin's
simulation record: slow, then fast. This probe measures that on all eight experimental
records and on the simulation's DIRECT contact-stress channel.

Two models on the package's own log-binned, inverted f(gamma), both starting at the
binned reversal minimum gamma_m, amplitudes eliminated exactly by weighted linear
least squares (variable projection), nonlinear parameters scanned on grids:

    E  (Eq. 4)   f = a + b exp(-(gamma - gamma_m)/lam)
    D  (lagged)  f = a + b exp(-max(gamma - gamma_m - delta, 0)/lam),  delta >= 0

E is D at delta = 0, so D can only fit better. What is reported is (i) delta itself in
strain units, (ii) delta against the fitted lam, (iii) the SSR ratio E/D. Bins carry
the within-bin scatter as sigma, but consecutive bins of a smoothed record are
correlated, so chi2 and AIC here are DESCRIPTIVE, not tests; the discriminator is the
size of delta against lam and the controls.

Controls, run first and required to pass before anything is printed about the data:
  C1 (negative)  a pure exponential from the minimum, sampled and binned like the real
                 records with multiplicative noise: delta must come back < 0.05 lam.
  C2 (positive)  the same with a true lag delta = 0.4: must come back within 0.05.

Affine invariance: f = A + B u, so delta and lam do not depend on (phi, phi_m, phi_0);
checked here at two triples (C3).

Run: python scripts/probe_rebuild_lag.py   (exit 0 = controls passed)

RESULT 2026-09-27 (logs/probe_rebuild_lag.txt): CONTROLS FAILED, so no delta from this
probe is quotable. C1 returned delta 0.15 for a "pure" exponential -- but the control
was MIS-SPECIFIED: its f is flat before gamma = 0.2, which IS a lag from reversal, and
the binned minimum of a flat noisy trough lands anywhere inside it (0.2 - 0.15 = 0.05).
C2 came back 0.555 = 0.4 + the same 0.155. The deeper defect is the definition: a delay
measured from the minimum of a broad trough is ill-posed. Descriptive only: SSR(E)/SSR(D)
is 1.0-6.3 on the experimental records and 2.9/4.5 on the simulation's two channels, so a
single exponential from the minimum misfits most records. Also note: scalar Wyart-Cates
with Richards Eq. 4 predicts NO reversal transient at constant rate (f_ss depends on
|sigma| only), so any reversal-based clock assumes f is reset at reversal.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import gate_g11 as G  # noqa: E402
import probe_g11_direct_control as P  # noqa: E402
from dst_fit.io_literature import read_lin_transients  # noqa: E402
from dst_fit.transient import (  # noqa: E402
    WC_EXPONENT,
    d_reduced_jamming_distance,
    frictional_fraction_from_viscosity,
    log_bin,
)

LAM_GRID = np.geomspace(0.02, 20.0, 400)
DELTA_GRID = np.concatenate([[0.0], np.geomspace(1e-3, 3.0, 300)])
SEED = 20260927


def binned_f(strain, eta_r, phi, phi_m, phi_0, n=WC_EXPONENT):
    g, e, err = log_bin(strain, eta_r)
    i0 = int(np.argmin(e))
    g, e, err = g[i0:], e[i0:], err[i0:]
    f = frictional_fraction_from_viscosity(e, phi, phi_m, phi_0, n=n)
    fe = phi / (phi_0 - phi_m) * np.abs(d_reduced_jamming_distance(e, n=n)) * err
    ok = np.isfinite(f) & np.isfinite(fe) & (fe > 0)
    return g[ok], f[ok], fe[ok]


def scan(g, y, s):
    """Return (ssr[delta, lam], best delta, best lam). Exact amplitudes per grid point."""
    w = 1.0 / s
    gm = g[0]
    ssr = np.empty((DELTA_GRID.size, LAM_GRID.size))
    for i, d in enumerate(DELTA_GRID):
        x = np.maximum(g - gm - d, 0.0)
        for j, lam in enumerate(LAM_GRID):
            A = np.column_stack([np.ones_like(g), np.exp(-x / lam)]) * w[:, None]
            coef, *_ = np.linalg.lstsq(A, y * w, rcond=None)
            r = A @ coef - y * w
            ssr[i, j] = float(r @ r)
    i, j = np.unravel_index(np.argmin(ssr), ssr.shape)
    return ssr, float(DELTA_GRID[i]), float(LAM_GRID[j])


def measure(g, y, s):
    ssr, d, lam = scan(g, y, s)
    ssr_e = float(ssr[0].min())
    lam_e = float(LAM_GRID[np.argmin(ssr[0])])
    ssr_d = float(ssr.min())
    return dict(delta=d, lam=lam, lam_E=lam_e, ssr_ratio=ssr_e / ssr_d,
                at_edge=d >= DELTA_GRID[-1] or lam in (LAM_GRID[0], LAM_GRID[-1]))


def synthetic(delta_true, rng):
    """Mimic a silica record: uniform strain 0..10, eta from f through WC, noise, bin."""
    phi, phi_m, phi_0 = 0.3913, 0.5083, 0.64
    g = np.arange(1, 20000) * 5.0e-4
    lam, f0, fss = 0.8, 0.0, 0.8
    x = np.maximum(g - 0.2 - delta_true, 0.0)
    f = fss + (f0 - fss) * np.exp(-x / lam)
    f[g < 0.2] = f0 + 0.05 * (0.2 - g[g < 0.2])  # a descending approach to the minimum
    phi_j = f * phi_m + (1 - f) * phi_0
    eta = (1 - phi / phi_j) ** (-WC_EXPONENT)
    eta *= 1 + 0.03 * rng.standard_normal(g.size)
    return g, eta, (phi, phi_m, phi_0)


def main() -> int:
    rng = np.random.default_rng(SEED)
    print(f"seed {SEED}; delta grid {DELTA_GRID[1]:.0e}-{DELTA_GRID[-1]} ({DELTA_GRID.size}), "
          f"lam grid {LAM_GRID[0]}-{LAM_GRID[-1]} ({LAM_GRID.size})")
    ok = True
    g, eta, tri = synthetic(0.0, rng)
    m = measure(*binned_f(g, eta, *tri))
    c1 = m["delta"] < 0.05 * m["lam"]
    print(f"C1 negative control (true delta 0, lam 0.8): delta {m['delta']:.4f} lam {m['lam']:.4f} "
          f"SSR E/D {m['ssr_ratio']:.3f}  {'PASS' if c1 else 'FAIL'}")
    g, eta, tri = synthetic(0.4, rng)
    m = measure(*binned_f(g, eta, *tri))
    c2 = abs(m["delta"] - 0.4) < 0.05
    print(f"C2 positive control (true delta 0.4, lam 0.8): delta {m['delta']:.4f} lam {m['lam']:.4f} "
          f"SSR E/D {m['ssr_ratio']:.3f}  {'PASS' if c2 else 'FAIL'}")
    ok &= c1 and c2

    print("\nexperimental records (Lin 2015 Fig_2b silica, Fig_2d PMMA), phi_0 = 0.64 gauge")
    print(f"  {'panel':7s}{'rate':>6s}{'N bins':>7s}{'g_min':>7s}{'delta':>8s}{'lam':>8s}"
          f"{'delta/lam':>10s}{'lam_E':>8s}{'SSR E/D':>9s}")
    rows = []
    for panel, steady in G.PAIRS:
        _, phi, phi_m, _ = G.steady_fit(steady, G.PHI0_DEFAULT)
        for rate, gg, ee in read_lin_transients(G.LIN_ROOT, panel):
            gb, fb, sb = binned_f(gg, ee, phi, phi_m, G.PHI0_DEFAULT)
            m = measure(gb, fb, sb)
            # C3: gauge check at a triple far off the fitted one
            m2 = measure(*binned_f(gg, ee, 0.30, 0.45, 0.58))
            c3 = abs(m2["delta"] - m["delta"]) < 1e-9 and abs(m2["lam"] - m["lam"]) < 1e-9
            ok &= c3
            rows.append((panel, rate, m))
            print(f"  {panel:7s}{rate:6g}{gb.size:7d}{gb[0]:7.3f}{m['delta']:8.3f}{m['lam']:8.3f}"
                  f"{m['delta'] / m['lam']:10.2f}{m['lam_E']:8.3f}{m['ssr_ratio']:9.2f}"
                  f"{'  EDGE' if m['at_edge'] else ''}{'' if c3 else '  C3 FAIL'}")

    contact = P.load("Fig_3b_contact.csv")
    total = P.load("Fig_3b_total.csv")
    gs = contact[:, 0]
    s = contact[:, 1] / total[:, 1]
    phi, phi_m, phi_0 = P.TRIPLES[0]
    print("\nsimulation record (Fig_3b), two channels on the same run")
    for name, eta in (("direct s", P.wc_forward(s, phi, phi_m, phi_0)), ("inverted", total[:, 1])):
        m = measure(*binned_f(gs, eta, phi, phi_m, phi_0))
        print(f"  {name:9s} delta {m['delta']:.3f} lam {m['lam']:.3f} delta/lam "
              f"{m['delta'] / m['lam']:.2f} lam_E {m['lam_E']:.3f} SSR E/D {m['ssr_ratio']:.2f}")

    incl = [(p, r, m) for p, r, m in rows if not (p == "Fig_2d" and r == 0.05)]
    dl = np.array([m["delta"] / m["lam"] for _, _, m in incl])
    print(f"\n7 in-range records: delta/lam {dl.min():.2f} to {dl.max():.2f}; "
          f"delta > 0.05 lam on {int((dl > 0.05).sum())}/7")
    print("controls " + ("PASSED" if ok else "FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
