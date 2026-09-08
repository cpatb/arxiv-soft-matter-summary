"""Recompute the bidet optimum and every design lever under three closures:
laminar Watson II (r^-5), the CFD-fitted turbulent closure (DEFAULT), and the
old commercial-point calibration (kept for comparison behind USE_CFD).

The CFD closure (cfd_openfoam/) replaces the single calibrated datum with a
computed (C_P, n_t). Key results and every lever are recomputed; where the CFD
overturns a headline (the standoff lever), that is reported explicitly.

Writes bidet_results.json and two figures.
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bidet_model import (LAMINAR, optimum, calibrate_turbulent, cfd_closure,
                         kappa_sinusoid, V0_from_dP, NU)

HERE = os.path.dirname(os.path.abspath(__file__))
CFD = json.load(open(os.path.join(HERE, "..", "cfd_openfoam", "cfd_fit_summary.json")))

r_max, tau_c = 0.015, 5.0
USE_CFD = True                              # CFD closure is the default (set False for old calibration)

# --- CFD-computed nominal closure (Re_d=1e4, H/d=6) ---
CP_CFD = CFD["nominal"]["C_P"]              # 1.35
NT_CFD = CFD["nominal"]["n_t"]              # 1.90
TURB = cfd_closure(CP_CFD, NT_CFD)

# --- old calibrated closure (comparison / backwards-compat) ---
NT_OLD = 1.05
COMM = dict(d=1.0e-3, dP=0.07e6, H=60e-3, theta=np.deg2rad(43))
TURB_CALIB = calibrate_turbulent(NT_OLD, COMM, r_max, tau_c)

def Lmin(o): return o["Q"] * 6e4
def tau_w_at(closure, r, V0, d): return closure.tau_w(r, V0, d)

print("=" * 74)
print("BIDET OPTIMUM -- laminar Watson II vs CFD-fitted turbulent closure")
print("nominal: r_max=15mm, tau_c=5Pa, dP=0.07MPa")
print(f"CFD closure (Re_d=1e4,H/d=6): n_t={NT_CFD:.3f}, C_P={CP_CFD:.3f}, r2={CFD['nominal']['r2']:.3f}")
_GI = CFD["grid_independence"]
print(f"  x2 grid SENSITIVITY (converged solves): n_t {_GI['rel_change']['n_t']*100:+.1f}%, "
      f"tau_w(r_max) {_GI['rel_change']['tau15']*100:+.1f}%, "
      f"C_P {_GI['rel_change']['C_P']*100:+.1f}%")
# NOT an error estimate, and the word "grid-independence" was removed from this line on
# 2026-09-08 for that reason. C_P DIVERGES under refinement (|eps21/eps32| = 4.589); n_t
# and tau_w(15 mm) are non-monotone. GRID_STUDY_2026-09-03.md has the whole account, and
# the summary JSON carries the same warning beside the numbers.
print(f"  x4: C_P {_GI['rel_change_x4']['C_P']*100:+.1f}% -- the steps GROW, so this is a "
      "DIVERGENCE, not a convergence: quote as sensitivity, never as an error bar")
print("=" * 74)

# ---------------------------------------------------------------------------
# KEY RESULT 1: wall shear at the cleaning radius for the COMMERCIAL nozzle
# (d=1mm, seated). CFD gives tau_w(r_max) ~ tau_c INDEPENDENTLY of the
# calibration -> commercial designs are near-optimal, not over-spec.
# ---------------------------------------------------------------------------
d_comm, V0_comm = 1.0e-3, V0_from_dP(0.07e6)
Re_comm = V0_comm * d_comm / NU
tau_lam_rmax = LAMINAR.tau_w(r_max, V0_comm, d_comm)
tau_cfd_rmax = TURB.tau_w(r_max, V0_comm, d_comm)
tau_cfd_direct = CFD["nominal"]["tau_at_15mm_Pa"]     # DIRECT CFD value (~5.9 Pa)
print("\n[Result 1] Wall shear at r_max for the commercial 1mm nozzle (Re_d=%.0f):" % Re_comm)
print(f"  laminar Watson II : {tau_lam_rmax:9.3f} Pa   (nominal prefactor; O(100)-uncertain -> value unreliable)")
print(f"  CFD closure       : {tau_cfd_rmax:9.3f} Pa   (fit)")
print(f"  CFD direct value  : {tau_cfd_direct:9.3f} Pa   (~ tau_c=5 Pa) -> commercial ~ just-cleans")
print("  NB: the robust CFD result is the SHAPE (r^-1.9, not r^-5); the absolute laminar")
print("      value carries the paper's O(100) matching-constant uncertainty and is not comparable.")

# ---------------------------------------------------------------------------
# LEVER 1: supply pressure. At fixed nozzle geometry tau_w ~ V0^{2+p}=V0^1.5,
# so the min-water optimum scales Q* ~ dP^{pexp}. Compute the exponent from the
# closure optimum (report the regime caveat -- see [Result 3]).
# ---------------------------------------------------------------------------
print("\n-- LEVER 1  supply pressure  Q*(dP)/Q*(0.07MPa) --")
dPs = np.array([0.05, 0.07, 0.10, 0.15]) * 1e6
qL = np.array([optimum(LAMINAR, dP=p, r_max=r_max, tau_c=tau_c)["Q"] for p in dPs])
qT = np.array([optimum(TURB,    dP=p, r_max=r_max, tau_c=tau_c)["Q"] for p in dPs])
pexpL = np.polyfit(np.log(dPs), np.log(qL), 1)[0]
pexpT = np.polyfit(np.log(dPs), np.log(qT), 1)[0]
for p, a, b in zip(dPs, qL / qL[1], qT / qT[1]):
    print(f"  dP={p/1e6:.2f} MPa   laminar {a:5.3f}   turbulent(CFD) {b:5.3f}")
print(f"  => Q* ~ dP^{pexpL:+.2f} (laminar, ~flat)   Q* ~ dP^{pexpT:+.2f} (CFD turbulent)")
print(f"     (analytic: dP-exponent = 0.5[1 - 2(2+p)/(p+n)] = {0.5*(1 - 2*(2-0.5)/(-0.5+NT_CFD)):+.2f}, "
      f"p=-1/2, n={NT_CFD:.2f}; V0 ~ dP^1/2)")

# ---------------------------------------------------------------------------
# LEVER 2: STANDOFF -- the CORRECTED headline. Use DIRECT CFD tau_w(r_max) at
# each H/d (fixed nozzle). tau_w ~ V0^1.5, Q ~ V0 => flow to reach tau_c scales
# as (tau_w_ref/tau_w)^{1/1.5}. CFD shows tau_w(r_max) is nearly standoff-
# independent -> a WEAK lever, overturning the paper's analytic 4.6x claim.
# ---------------------------------------------------------------------------
print("\n-- LEVER 2  standoff (CFD, fixed nozzle, Re_d=1e4)  [CORRECTS the 4.6x claim] --")
so = CFD["standoff_Re10000"]
hd_keys = sorted(so.keys(), key=lambda k: int(k))
tau_ref = so["6"]["tau15"]
print(f"  {'H/d':>4s}{'tau_w(r_max) [Pa]':>18s}{'flow x (rel H/d=6)':>20s}{'fit r2':>8s}")
for k in hd_keys:
    t = so[k]["tau15"]; mult = (tau_ref / t) ** (1 / 1.5)
    flag = "  (spread jet, unreliable fit)" if so[k]["r2"] < 0.9 else ""
    print(f"  {k:>4s}{t:>18.2f}{mult:>20.3f}{so[k]['r2']:>8.3f}{flag}")
# laminar prediction for contrast (its wrong, strong sign):
qLh = np.array([optimum(LAMINAR, dP=0.07e6, r_max=r_max, tau_c=tau_c, H=hd * 0.8e-3)["Q"] for hd in (6, 20, 60)])
print(f"  laminar r^-5 predicts a STRONG (unphysical) standoff sign: "
      f"Q*(H/d=60)/Q*(6) = {qLh[2]/qLh[0]:.3f}  (paper claimed a 4.6x lever)")
print("  => CFD: coherent-jet standoff lever is WEAK (<~10% over H/d 1-20). Far-field")
print("     wall shear is set by conserved jet momentum flux, not centerline decay.")
print("     Caveat: single-phase CFD; jet breakup at large H (not modelled) could add a penalty.")
standoff_flow_span = max((tau_ref / so[k]["tau15"]) ** (1/1.5) for k in ("1", "6", "20")) / \
                     min((tau_ref / so[k]["tau15"]) ** (1/1.5) for k in ("1", "6", "20"))

# ---------------------------------------------------------------------------
# LEVER 3: oblique impingement. Normal velocity component V0 cos(theta); the
# cleaning uses the normal wall shear. Geometric, survives the closure change.
# ---------------------------------------------------------------------------
print("\n-- LEVER 3  oblique impingement  Q*(theta)/Q*(0) --")
ths = np.deg2rad([0, 15, 30, 43, 53, 60])
qLt = np.array([optimum(LAMINAR, dP=0.07e6, r_max=r_max, tau_c=tau_c, theta=t)["Q"] for t in ths])
qTt = np.array([optimum(TURB,    dP=0.07e6, r_max=r_max, tau_c=tau_c, theta=t)["Q"] for t in ths])
for t, a, b in zip(ths, qLt / qLt[0], qTt / qTt[0]):
    print(f"  theta={np.rad2deg(t):2.0f}deg  laminar {a:5.3f} (sec={1/np.cos(t):.3f})   turbulent(CFD) {b:5.3f}")

# ---------------------------------------------------------------------------
# LEVER 4: oscillation. kappa=2/pi dwell non-uniformity; coverage exponent
# n/(2+p) now uses the CFD n_t. Geometric saving, recomputed.
# ---------------------------------------------------------------------------
kap = kappa_sinusoid()
H_osc, alpha = 80e-3, np.deg2rad(15)
A = H_osc * alpha
r_c = (np.sqrt(4 * A**2 + np.pi**2 * r_max**2) - 2 * A) / np.pi
r_c_eff = r_c / np.sqrt(kap)
def cover_ratio(closure):
    e = closure.n / (2 + closure.p)
    return (min(r_c_eff, r_max) / r_max) ** e
print(f"\n-- LEVER 4  oscillation (alpha=15deg,H=80mm): kappa=2/pi={kap:.3f}, "
      f"r_c={r_c*1e3:.1f}mm, r_c_eff={r_c_eff*1e3:.1f}mm --")
print(f"  laminar        Q*_osc/Q* = {cover_ratio(LAMINAR):.3f}  (save {100*(1-cover_ratio(LAMINAR)):.0f}%)")
print(f"  turbulent(CFD) Q*_osc/Q* = {cover_ratio(TURB):.3f}  (save {100*(1-cover_ratio(TURB)):.0f}%)")

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
# Fig 1: tau_w(r) laminar vs CFD-turbulent, with the CFD data points overlaid.
csv = os.path.join(HERE, "..", "cfd_openfoam", "case_Re10000_H6", "tau_w_r.csv")
cfd_data = np.genfromtxt(csv, delimiter=",", names=True)
r = np.linspace(2e-3, 20e-3, 200)
V0, d = V0_from_dP(0.07e6), 1.0e-3
fig, ax = plt.subplots(figsize=(7, 4.4))
ax.loglog(cfd_data["r_mm"], cfd_data["tau_w_Pa"], "o", ms=3, color="C0", alpha=.7,
          label="CFD data (k-$\\omega$ SST)")
ax.loglog(r * 1e3, TURB.tau_w(r, V0, d), "-", color="C1", lw=2,
          label=f"CFD-fitted turbulent ($r^{{-{NT_CFD:.2f}}}$)")
ax.loglog(r * 1e3, LAMINAR.tau_w(r, V0, d), "--", color="C3",
          label="laminar Watson II ($r^{-5}$)")
ax.axhline(tau_c, ls=":", c="k", label=r"$\tau_c=5$ Pa")
ax.axvline(r_max * 1e3, ls="--", c="grey", label=r"$r_{max}=15$ mm")
ax.set_xlabel("radius r (mm)"); ax.set_ylabel(r"wall shear $\tau_w$ (Pa)")
ax.set_title("CFD wall shear follows a shallow turbulent $r^{-1.9}$ law,\n"
             "not the steep laminar Watson $r^{-5}$ (prefactor $O(100)$-uncertain)")
ax.set_ylim(1e-2, 1e3)
ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.25)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_tau_closures.png"), dpi=140)

# Fig 2: Q* vs assumed decay exponent n, with the CFD point marked.
ns = np.linspace(1.0, 5, 40)
qn = []
for nn in ns:
    cl = calibrate_turbulent(nn, COMM, r_max, tau_c)      # recalibrate for a fair absolute compare
    qn.append(optimum(cl, dP=0.07e6, r_max=r_max, tau_c=tau_c, H=60e-3)["Q"] * 6e4)
fig2, ax2 = plt.subplots(figsize=(7, 4.2))
ax2.plot(ns, qn, "o-", ms=3, color="C0")
ax2.axvline(NT_CFD, ls="-", c="C1", lw=2, label=f"CFD-computed $n_t={NT_CFD:.2f}$")
ax2.axvline(1.05, ls=":", c="C2", label="Phares developed WJ (assumed 1.05)")
ax2.axvline(5.0, ls=":", c="C3", label="laminar Watson II")
ax2.set_xlabel("wall-shear decay exponent  $n$"); ax2.set_ylabel(r"$Q^*$ at $H/d{=}60$ (L/min)")
ax2.set_title("Optimum water use vs assumed decay exponent\n(prefactor recalibrated at each $n$; CFD pins $n_t$)")
ax2.legend(fontsize=8); ax2.grid(alpha=0.25)
fig2.tight_layout(); fig2.savefig(os.path.join(HERE, "fig_Qstar_vs_n.png"), dpi=140)

# ---------------------------------------------------------------------------
# Persist
# ---------------------------------------------------------------------------
out = dict(
    closure_source="CFD" if USE_CFD else "calibrated",
    cfd_closure=dict(C_P=CP_CFD, n_t=NT_CFD, r2=CFD["nominal"]["r2"]),
    tau_w_rmax_commercial=dict(laminar_Pa=float(tau_lam_rmax), cfd_fit_Pa=float(tau_cfd_rmax),
                               cfd_direct_Pa=float(tau_cfd_direct), tau_c_Pa=tau_c,
                               note="CFD tau_w(r_max) ~ tau_c => commercial nozzle near-optimal"),
    pressure_exp=dict(laminar=float(pexpL), cfd_turbulent=float(pexpT)),
    standoff_lever=dict(
        finding="WEAK -- CFD overturns the paper's 4.6x claim",
        flow_span_Hd_1_to_20=float(standoff_flow_span),
        tau_w_rmax_vs_Hd_Re1e4={k: so[k]["tau15"] for k in hd_keys},
        laminar_wrong_ratio_60_over_6=float(qLh[2] / qLh[0]),
        caveat="single-phase CFD; jet breakup at large H not modelled"),
    oblique_43=dict(laminar=float(qLt[3] / qLt[0]), cfd_turbulent=float(qTt[3] / qTt[0])),
    oscillation=dict(kappa=kap, r_c_mm=r_c * 1e3,
                     save_laminar=1 - cover_ratio(LAMINAR), save_cfd=1 - cover_ratio(TURB)),
)
json.dump(out, open(os.path.join(HERE, "bidet_results.json"), "w"), indent=2)
print("\nwrote fig_tau_closures.png, fig_Qstar_vs_n.png, bidet_results.json")
print(f"\nSUMMARY: CFD n_t={NT_CFD:.2f} (reversal from laminar 5 holds). "
      f"Standoff lever WEAK (flow span {standoff_flow_span:.2f}x over H/d 1-20), "
      f"correcting the 4.6x claim. tau_w(r_max)~tau_c: commercial near-optimal.")
