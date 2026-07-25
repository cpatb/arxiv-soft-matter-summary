"""Recompute the bidet optimum and every design lever under both closures.
Prints a comparison the manuscript rewrite will use, and writes two figures."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bidet_model import (LAMINAR, optimum, calibrate_turbulent, kappa_sinusoid,
                         V0_from_dP, NU, film_tau_numeric, laminar_tauII_direct)

r_max, tau_c = 0.015, 5.0
NT = 1.05                      # developed-turbulent radial decay exponent (Phares)
COMM = dict(d=1.0e-3, dP=0.07e6, H=60e-3, theta=np.deg2rad(43))
TURB = calibrate_turbulent(NT, COMM, r_max, tau_c)

def Lmin(o): return o["Q"] * 6e4

print("="*72)
print("BIDET OPTIMUM — laminar Watson II vs calibrated turbulent closure")
print(f"nominal: r_max=15mm, tau_c=5Pa, dP=0.07MPa; turbulent nt={NT}, C_P={TURB.C_tau:.4g}")
print("="*72)

# --- headline optimum at nominal, idealized standoff ---
oL = optimum(LAMINAR, dP=0.07e6, r_max=r_max, tau_c=tau_c)
oT = optimum(TURB,    dP=0.07e6, r_max=r_max, tau_c=tau_c)
print(f"{'':22s}{'d* (mm)':>10s}{'Q* (L/min)':>12s}{'Re_d':>10s}")
print(f"{'laminar (r^-5)':22s}{oL['d_star']*1e3:>10.2f}{Lmin(oL):>12.3f}{oL['Re_d']:>10.0f}")
print(f"{'turbulent (r^-1.05)':22s}{oT['d_star']*1e3:>10.2f}{Lmin(oT):>12.3f}{oT['Re_d']:>10.0f}")

# --- LEVER 1: supply pressure ---
print("\n-- Pressure lever  Q*(dP)/Q*(0.07MPa) --")
dPs = np.array([0.05, 0.07, 0.10, 0.15]) * 1e6
qL = np.array([optimum(LAMINAR, dP=p, r_max=r_max, tau_c=tau_c)["Q"] for p in dPs])
qT = np.array([optimum(TURB,    dP=p, r_max=r_max, tau_c=tau_c)["Q"] for p in dPs])
for p, a, b in zip(dPs, qL/qL[1], qT/qT[1]):
    print(f"  dP={p/1e6:.2f} MPa   laminar {a:5.3f}   turbulent {b:5.3f}")
pexpL = np.polyfit(np.log(dPs), np.log(qL), 1)[0]
pexpT = np.polyfit(np.log(dPs), np.log(qT), 1)[0]
print(f"  => Q* ~ dP^{pexpL:+.2f} (laminar)   Q* ~ dP^{pexpT:+.2f} (turbulent)")

# --- LEVER 2: standoff ---
print("\n-- Standoff lever  Q*(H)/Q*(H/d=6) --")
d_ref_T = oT["d_star"]; d_ref_L = oL["d_star"]
Hs = np.array([6, 20, 40, 60, 100])
qLh = np.array([optimum(LAMINAR, dP=0.07e6, r_max=r_max, tau_c=tau_c, H=hd*d_ref_L)["Q"] for hd in Hs])
qTh = np.array([optimum(TURB,    dP=0.07e6, r_max=r_max, tau_c=tau_c, H=hd*d_ref_T)["Q"] for hd in Hs])
for hd, a, b in zip(Hs, qLh/qLh[0], qTh/qTh[0]):
    print(f"  H/d={hd:3d}   laminar {a:6.3f} (sign!)   turbulent {b:6.3f}")

# --- LEVER 3: oblique ---
print("\n-- Oblique lever  Q*(theta)/Q*(0) --")
ths = np.deg2rad([0, 15, 30, 43, 53, 60])
qLt = np.array([optimum(LAMINAR, dP=0.07e6, r_max=r_max, tau_c=tau_c, theta=t)["Q"] for t in ths])
qTt = np.array([optimum(TURB,    dP=0.07e6, r_max=r_max, tau_c=tau_c, theta=t)["Q"] for t in ths])
for t, a, b in zip(ths, qLt/qLt[0], qTt/qTt[0]):
    print(f"  theta={np.rad2deg(t):2.0f}deg  laminar {a:5.3f} (=sec: {1/np.cos(t):.3f})   turbulent {b:5.3f}")

# --- LEVER 4: oscillation (kappa-corrected coverage) ---
kap = kappa_sinusoid()
H_osc, alpha = 80e-3, np.deg2rad(15)
A = H_osc * alpha
# stadium coverage: instantaneous radius r_c s.t. pi r_c^2 + 4 A r_c = pi r_max^2, kappa-derated
r_c = (np.sqrt(4*A**2 + np.pi**2 * r_max**2) - 2*A) / np.pi
# effective required instantaneous radius with dwell non-uniformity: r_c_eff = r_c / sqrt(kappa)
r_c_eff = r_c / np.sqrt(kap)
# Q* ~ r_max^{n/(2+p)} at fixed dP; ratio uses (r_c_eff/r_max)
def cover_ratio(closure):
    e = closure.n / (2 + closure.p)
    return (min(r_c_eff, r_max) / r_max)**e
print(f"\n-- Oscillation (alpha=15deg, H=80mm): kappa=2/pi={kap:.3f}, r_c={r_c*1e3:.1f}mm, r_c_eff={r_c_eff*1e3:.1f}mm --")
print(f"  laminar   Q*_osc/Q* = {cover_ratio(LAMINAR):.3f}  (save {100*(1-cover_ratio(LAMINAR)):.0f}%)")
print(f"  turbulent Q*_osc/Q* = {cover_ratio(TURB):.3f}  (save {100*(1-cover_ratio(TURB)):.0f}%)")

# ---------- Figure 1: tau_w(r), laminar vs turbulent, with r_max marked ----------
r = np.linspace(2e-3, 20e-3, 200)
V0 = V0_from_dP(0.07e6); d = 1.0e-3
fig, ax = plt.subplots(figsize=(7, 4.4))
ax.loglog(r*1e3, LAMINAR.tau_w(r, V0, d), label="laminar Watson II  ($r^{-5}$)")
ax.loglog(r*1e3, TURB.tau_w(r, V0, d),    label=f"turbulent, calibrated  ($r^{{-{NT}}}$)")
ax.axhline(tau_c, ls=":", c="k", label=r"$\tau_c=5$ Pa")
ax.axvline(r_max*1e3, ls="--", c="grey", label=r"$r_{max}=15$ mm")
ax.set_xlabel("radius r (mm)"); ax.set_ylabel(r"wall shear $\tau_w$ (Pa)")
ax.set_title("Laminar Region-II vs turbulent wall shear: the r$^{-5}$ law\ncollapses far too fast at the cleaning radius")
ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.25)
fig.tight_layout(); fig.savefig("fig_tau_closures.png", dpi=140)

# ---------- Figure 2: Q* vs decay exponent n (band the postdoc asked for) ----------
ns = np.linspace(1.05, 5, 40)
from bidet_model import Closure, RHO
qn = []
for nn in ns:
    p = -0.5 if nn < 3 else 1.0   # illustrative: turbulent branch low-n, laminar high-n
    # recalibrate C_tau at each n to the commercial point for a fair absolute comparison
    cl = calibrate_turbulent(nn, COMM, r_max, tau_c)
    qn.append(optimum(cl, dP=0.07e6, r_max=r_max, tau_c=tau_c, H=60e-3)["Q"]*6e4)
fig2, ax2 = plt.subplots(figsize=(7, 4.2))
ax2.plot(ns, qn, "o-", ms=3)
ax2.axvline(1.05, ls=":", c="C1", label="turbulent (Phares)")
ax2.axvline(5.0, ls=":", c="C3", label="laminar Watson II")
ax2.set_xlabel("wall-shear decay exponent  n"); ax2.set_ylabel(r"$Q^*$ at $H/d{=}60$ (L/min)")
ax2.set_title("Optimum water use vs the assumed wall-shear decay exponent\n(prefactor recalibrated to the commercial point at each n)")
ax2.legend(fontsize=8); ax2.grid(alpha=0.25)
fig2.tight_layout(); fig2.savefig("fig_Qstar_vs_n.png", dpi=140)

# ---------- persist key numbers ----------
out = dict(
    laminar=dict(d_mm=oL['d_star']*1e3, Q_Lmin=Lmin(oL), Re_d=oL['Re_d']),
    turbulent=dict(d_mm=oT['d_star']*1e3, Q_Lmin=Lmin(oT), Re_d=oT['Re_d'], C_P=TURB.C_tau),
    pressure_exp=dict(laminar=pexpL, turbulent=pexpT),
    standoff_ratio_6to60=dict(laminar=float(qLh[3]/qLh[0]), turbulent=float(qTh[3]/qTh[0])),
    oblique_43=dict(laminar=float(qLt[3]/qLt[0]), turbulent=float(qTt[3]/qTt[0])),
    kappa=kap, r_c_mm=r_c*1e3,
    osc_save=dict(laminar=1-cover_ratio(LAMINAR), turbulent=1-cover_ratio(TURB)),
)
json.dump(out, open("bidet_results.json", "w"), indent=2)
print("\nwrote fig_tau_closures.png, fig_Qstar_vs_n.png, bidet_results.json")
