"""Probe: how much does the fitted strain clock depend on WHERE the fit starts?

fit_recovery starts at the binned reversal minimum. The trough is broad and noisy
(probe_rebuild_lag.py: a delay measured from it is ill-posed), so the start is uncertain
by ~0.1-0.2 strain. This refits every in-range Lin record with the start moved to each
log bin within +/-WINDOW strain of the minimum and prints lambda against start strain,
flagging grid-edge solutions (lambda at the scan bound = no interior optimum).

Run: python scripts/probe_start_sensitivity.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import gate_g11 as G  # noqa: E402
from dst_fit.io_literature import read_lin_transients  # noqa: E402
from dst_fit.transient import (  # noqa: E402
    LAMBDA_GRID,
    d_reduced_jamming_distance,
    fit_relaxation,
    frictional_fraction_from_viscosity,
    log_bin,
)

WINDOW = 0.2
EXCLUDED = {("Fig_2d", 0.05)}  # plateau in f outside [0, 1]; gate G11 rule


def main() -> None:
    lo_edge, hi_edge = LAMBDA_GRID[0], LAMBDA_GRID[1]
    for panel, steady in G.PAIRS:
        _, phi, phi_m, _ = G.steady_fit(steady, G.PHI0_DEFAULT)
        for rate, gg, ee in read_lin_transients(G.LIN_ROOT, panel):
            if (panel, rate) in EXCLUDED:
                continue
            gb, eb, sb = log_bin(gg, ee)
            i0 = int(np.argmin(eb))
            print(f"{panel} rate {rate:g}: binned min at gamma {gb[i0]:.4f} (bin {i0}), "
                  f"eta_min {eb[i0]:.3f}; {gb.size} bins")
            for j in range(gb.size - 12):
                if abs(gb[j] - gb[i0]) > WINDOW:
                    continue
                f = frictional_fraction_from_viscosity(eb[j:], phi, phi_m, G.PHI0_DEFAULT)
                fe = phi / (G.PHI0_DEFAULT - phi_m) * np.abs(d_reduced_jamming_distance(eb[j:])) * sb[j:]
                s = fit_relaxation(gb[j:], f, fe)
                edge = s.lam <= lo_edge * 1.0001 or s.lam >= hi_edge * 0.9999
                amp = s.y_0 - s.y_inf
                print(f"   start {gb[j]:.4f} ({gb[j] - gb[i0]:+.3f})  lam {s.lam:8.4f} "
                      f"[{s.lo:.4f},{s.hi:.4f}]  amp {amp:+.3f}  chi2r {s.chi2_reduced:7.1f}"
                      f"{'  GRID EDGE' if edge else ''}")


if __name__ == "__main__":
    main()
