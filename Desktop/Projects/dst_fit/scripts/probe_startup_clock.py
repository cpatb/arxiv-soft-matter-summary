"""Probe: the strain clock from constant-stress START-UP, graded against
PREREG_startup_clock_2026-09-27.md.

Data: Richards et al. 2019 deposit, 4 um silica in 87 wt% glycerol-water, phi = 0.574,
Figure3/Data/silica87glyc/<k>_Stress_<sigma>.csv (TA export: Shear Rate, Shear Stress,
Viscosity [mPa s], Shear Strain), steady curve gw_phi_0p574_flow.txt (CR line endings;
columns include 'avg visc' [mPa s] and 'rel. visc').

Method as registered: start at the rate maximum within gamma < 0.3 (end of rotor
spin-up), log-bin, u = 1/(1 - eta_r^(-1/n)), n = 2, fit_relaxation (variable projection;
affine-invariant so no packing triple is needed). eta_s = eta/eta_r from the steady file.

Run: python scripts/probe_startup_clock.py   -> logs/probe_startup_clock.txt
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dst_fit.transient import d_reduced_jamming_distance, fit_relaxation, log_bin  # noqa: E402

D = ROOT / ("Competing Timescales Lead to Oscillations in Shear-Thickening Suspensions"
            "/Figure3/Data/silica87glyc")
N_EXP = 2.0
START_WINDOW = 0.3
P1_BAR = 0.30
P2_FACTOR = 2.0
P3_TOL = 0.30
REVERSAL_MIN_SILICA = 0.54  # smallest reversal clock on Lin's silica (numbers.tex allLoSi)


def steady():
    txt = (D / "gw_phi_0p574_flow.txt").read_bytes().decode("latin-1")
    rows = [r.split("\t") for r in txt.replace("\r\n", "\r").split("\r") if r.strip()]
    head = rows[0]
    i_s = head.index("Shear Stress (Pa)")
    i_v = head.index("avg visc")
    i_r = head.index("rel. visc")
    data = np.array([[float(r[i_s]), float(r[i_v]), float(r[i_r])] for r in rows[1:]])
    data = data[np.argsort(data[:, 0])]
    sigma, eta_mpas, eta_r = data.T
    eta_s = eta_mpas / eta_r  # mPa s, one value per row -- must be constant
    return sigma, eta_mpas, eta_r, eta_s


def step_files():
    out = []
    for p in D.glob("*_Stress_*.csv"):
        m = re.match(r"(\d+)_Stress_([0-9.]+)\.csv", p.name)
        out.append((int(m.group(1)), float(m.group(2)), p))
    return sorted(out)


def load(p):
    r = list(csv.DictReader(open(p, encoding="utf-8")))
    g = np.array([float(x["Shear Strain"]) for x in r])
    gd = np.array([float(x["Shear Rate"]) for x in r])
    eta = np.array([float(x["Viscosity"]) for x in r])
    s = np.array([float(x["Shear Stress"]) for x in r])
    return g, gd, eta, s


def main() -> int:
    sigma_ss, eta_ss, etar_ss, eta_s_rows = steady()
    eta_s = float(np.median(eta_s_rows))
    print(f"steady curve: N = {sigma_ss.size}, sigma {sigma_ss.min():g}-{sigma_ss.max():g} Pa; "
          f"eta_s = {eta_s:.2f} mPa s (row spread {eta_s_rows.min():.2f}-{eta_s_rows.max():.2f})")
    plateau = float(np.median(etar_ss[sigma_ss <= 1.6]))
    onset = float(sigma_ss[np.argmax(etar_ss > 2 * plateau)])
    print(f"low-stress plateau eta_r = {plateau:.1f}; onset (eta_r > 2x plateau) at {onset:g} Pa")

    rows = []
    for k, sig, p in step_files():
        g, gd, eta, s = load(p)
        n0 = g.size
        ok = np.isfinite(g) & np.isfinite(eta) & (eta > 0) & (np.diff(np.r_[-1.0, g]) > 0)
        g, gd, eta = g[ok], gd[ok], eta[ok]
        win = g < START_WINDOW
        i0 = int(np.argmax(np.where(win, gd, -np.inf)))
        er = eta / eta_s
        gb, eb, sb = log_bin(g[i0:] - 0.0, er[i0:])
        u = 1.0 / (1.0 - eb ** (-1.0 / N_EXP))
        ue = np.abs(d_reduced_jamming_distance(eb, n=N_EXP)) * sb
        good = np.isfinite(u) & np.isfinite(ue) & (ue > 0)
        fit = fit_relaxation(gb[good], u[good], ue[good])
        ss_eta_r = float(np.interp(sig, sigma_ss, etar_ss))
        end_eta_r = float(np.median(er[-200:]))
        rise = float(np.median(er[-200:]) / er[i0])
        rows.append(dict(k=k, sigma=sig, above=sig >= onset, gamma_start=float(g[i0]),
                         lam=fit.lam, lo=fit.lo, hi=fit.hi, chi2=fit.chi2_reduced,
                         rise=rise, end=end_eta_r, ss=ss_eta_r, n=n0, nfit=int(good.sum()),
                         edge=fit.lam in (fit.grid[0], fit.grid[-1])))
        r = rows[-1]
        print(f"{sig:5.1f} Pa {'ABOVE' if r['above'] else 'below'}  N {n0}->{g.size}  start gamma "
              f"{r['gamma_start']:.4f}  clock {fit.lam:.4f} [{fit.lo:.4f},{fit.hi:.4f}]  chi2r "
              f"{fit.chi2_reduced:8.1f}  eta rise x{rise:5.2f}  end eta_r {end_eta_r:7.1f} "
              f"vs steady {ss_eta_r:7.1f}{'  GRID EDGE' if r['edge'] else ''}")

    above = [r for r in rows if r["above"] and not r["edge"]]
    below = [r for r in rows if not r["above"] and not r["edge"]]
    print()
    if above:
        worst = max(r["lam"] for r in above)
        p1 = worst < P1_BAR
        print(f"P1  above onset ({len(above)} steps): clocks {min(r['lam'] for r in above):.3f}-"
              f"{worst:.3f}; bar < {P1_BAR}  -> {'PASS' if p1 else 'FAIL'}; smallest reversal "
              f"clock / largest start-up clock = {REVERSAL_MIN_SILICA / worst:.2f}x")
    if above and below:
        ma = float(np.median([r["lam"] for r in above]))
        mb = float(np.median([r["lam"] for r in below]))
        ratio = max(ma, mb) / min(ma, mb)
        p2 = ratio <= P2_FACTOR
        print(f"P2  median clock below onset {mb:.3f} ({len(below)} steps) vs above {ma:.3f}: "
              f"ratio {ratio:.2f} -> {'WITHIN 2x: start-up transient is NOT primarily f' if p2 else 'differ by >2x'}")
    dev = [abs(r["end"] - r["ss"]) / r["ss"] for r in rows]
    print(f"P3  end-of-step eta_r vs steady: worst {100 * max(dev):.1f}% "
          f"(bar {100 * P3_TOL:.0f}%) -> {'PASS' if max(dev) < P3_TOL else 'FAIL'}; "
          f"failing: {[r['sigma'] for r, d in zip(rows, dev) if d >= P3_TOL]}")
    edges = [r["sigma"] for r in rows if r["edge"]]
    print(f"grid-edge fits excluded: {edges}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
