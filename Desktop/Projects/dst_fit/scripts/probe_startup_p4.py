"""P4 of PREREG_startup_clock_2026-09-27.md (Addendum 1): the jam strain after restart
against the recoil during the preceding rest. Fit-free.

gamma_j = strain at the steepest fall of the RAW rate (max -d(rate)/d(gamma), gamma < 0.3);
R = |net strain of the preceding Relax file|. Registered bar: |gamma_j - R| < 0.03 on
>= 80 % of above-onset steps AND Spearman rho(R, gamma_j) > 0.5.

Result 2026-09-27 ~23:42 (logs/probe_startup_p4.txt): FAIL (57 %, rho = +0.19). The recoil
spans only 1.4x (0.049-0.070), so the design cannot test tracking. Post hoc and labelled
so: gamma_j > R on every step, and gamma_j - R falls with stress.
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import probe_startup_clock as P  # noqa: E402

JAM_WINDOW = 0.3
P4_DIFF = 0.03
P4_FRAC = 0.8
P4_RHO = 0.5


def jam_table():
    files = {int(re.match(r"(\d+)_", p.name).group(1)): p for p in P.D.glob("*_*.csv")}
    sig_ss, _, etar_ss, _ = P.steady()
    plateau = float(np.median(etar_ss[sig_ss <= 1.6]))
    onset = float(sig_ss[np.argmax(etar_ss > 2 * plateau)])
    rows = []
    for k, sig, _ in P.step_files():
        if k == 0:
            continue  # no preceding rest
        rel = files[k - 1]
        assert "Relax" in rel.name, rel
        gr = np.array([float(x["Shear Strain"]) for x in csv.DictReader(open(rel, encoding="utf-8"))])
        g, gd, _, _ = P.load(files[k])
        m = g < JAM_WINDOW
        j = int(np.argmin(np.gradient(gd[m], g[m])))
        rows.append(dict(sigma=sig, above=sig >= onset, R=float(abs(gr[-1] - gr[0])),
                         gamma_j=float(g[m][j])))
    return rows, onset


def jam_times():
    """(sigma, gamma_j, t_j) for above-onset steps, t from the 'Interval Time' column.

    Post-hoc confound check (not registered): if the jam were set by rotor spin-up it would
    collapse in TIME rather than strain. Measured 2026-09-27: spreads 1.90x (strain) vs
    2.07x (time) over a 2.75x stress range -- indistinguishable, so inertia is not excluded.
    """
    files = {int(re.match(r"(\d+)_", p.name).group(1)): p for p in P.D.glob("*_*.csv")}
    rows, onset = jam_table()
    out = []
    for k, sig, _ in P.step_files():
        if k == 0 or sig < onset:
            continue
        r = list(csv.DictReader(open(files[k], encoding="utf-8")))
        g = np.array([float(x["Shear Strain"]) for x in r])
        t = np.array([float(x["Interval Time"]) for x in r])
        gd = np.array([float(x["Shear Rate"]) for x in r])
        m = g < JAM_WINDOW
        j = int(np.argmin(np.gradient(gd[m], g[m])))
        out.append((sig, float(g[m][j]), float(t[m][j])))
    return np.array(out)


def main() -> int:
    rows, onset = jam_table()
    for r in rows:
        print(f"{r['sigma']:5.1f} Pa {'ABOVE' if r['above'] else 'below'}  R {r['R']:.4f}  "
              f"gamma_j {r['gamma_j']:.4f}  gamma_j-R {r['gamma_j'] - r['R']:+.4f}")
    ab = [r for r in rows if r["above"]]
    d = np.array([r["gamma_j"] - r["R"] for r in ab])
    frac = float(np.mean(np.abs(d) < P4_DIFF))
    rho, pv = spearmanr([r["R"] for r in ab], [r["gamma_j"] for r in ab])
    print(f"P4 above onset (N={len(ab)}): |gamma_j-R|<{P4_DIFF} on {100 * frac:.0f}% (bar "
          f"{100 * P4_FRAC:.0f}%), rho {rho:+.3f} (p {pv:.3f}, bar >{P4_RHO}) -> "
          f"{'PASS' if frac >= P4_FRAC and rho > P4_RHO else 'FAIL'}")
    print(f"post hoc: gamma_j > R on {sum(r['gamma_j'] > r['R'] for r in rows)}/{len(rows)} steps")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
