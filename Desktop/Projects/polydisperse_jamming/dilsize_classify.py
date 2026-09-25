"""dilsize_classify.py -- class rule for PREREG_DILSIZE_2026-09-25 (Amendment 1), fixed before any cell.

spread_*.dat columns (0-based): 1 t/tau, 4 h/D0, 7 jam_frac, 8 umax (U0 units), 10 zc/D0, 11 Rg_r/D0.
Health (the project's admissibility gate, grade_capvol.py): rc 0, t/tau reaches 5.9, and after T_IMPACT
umax <= 1, zc/D0 <= 0.5, h/D0 <= 1.
Class:
  FLUID      J_max <= J_ON
  TRANSIENT  J_max > J_ON but J_end <= J_ON (jammed, then unjammed: the referee's counterexample)
  BALL       J_end > J_ON, Rg_end >= 0.9 Rg_peak and t_95 <= 0.5 (arrests near its peak and stays)
  HYBRID     J_end > J_ON otherwise (stays jammed yet keeps spreading or retracts)
J_end = mean jam_frac over t/tau >= TEND - 0.5; Rg_end likewise.
"""
from __future__ import annotations

import sys

import numpy as np

J_ON = 0.01
T_IMPACT = 0.5          # grade_capvol.py:17


def classify(path: str) -> dict:
    a = np.loadtxt(path, comments="#", ndmin=2)
    t, h, J, um, zc, rg = a[:, 1], a[:, 4], a[:, 7], a[:, 8], a[:, 10], a[:, 11]
    late = t >= T_IMPACT
    end = t >= t[-1] - 0.5
    adm = bool(t[-1] >= 5.9 and um[late].max() <= 1.0 and zc[late].max() <= 0.5 and h[late].max() <= 1.0)
    jmax, jend = float(J.max()), float(J[end].mean())
    rpk, rend = float(rg.max()), float(rg[end].mean())
    t95 = float(t[np.argmax(rg >= 0.95 * rpk)])
    if jmax <= J_ON:
        c = "FLUID"
    elif jend <= J_ON:
        c = "TRANSIENT"
    elif rend >= 0.9 * rpk and t95 <= 0.5:
        c = "BALL"
    else:
        c = "HYBRID"
    return dict(cls=c, adm=adm, J_max=jmax, J_end=jend, Rg_peak=rpk, Rg_end=rend, t95=t95, t_last=float(t[-1]))


if __name__ == "__main__":
    for p in sys.argv[1:]:
        r = classify(p)
        print(f"{r['cls']:9s} adm {int(r['adm'])}  J_max {r['J_max']:.3f} J_end {r['J_end']:.3f}  Rg_peak {r['Rg_peak']:.3f} "
              f"Rg_end {r['Rg_end']:.3f} t95 {r['t95']:.2f}  <- {p}")
