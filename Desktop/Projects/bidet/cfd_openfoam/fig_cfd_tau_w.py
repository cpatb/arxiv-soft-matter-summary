#!/usr/bin/env python3
"""
fig_cfd_tau_w.py -- Figure: CFD tau_w(r) vs the fitted turbulent power law vs the
laminar Watson r^-5 law, at the nominal case (Re_d=1e4, H/d=6).

Reads <case>/tau_w_r.csv and the fitted (C_P, n_t) from sweep_results.jsonl (or
recomputes via fit_closure). Saves fig_cfd_tau_w.png in the project.
"""
import sys, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RHO = 1000.0
D = 1.0e-3
NU = 8.9e-7
C2 = 27.0 / 50.0  # Watson Region-II prefactor (laminar), tau = c2 rho V0^2 Re_a (a/r)^5


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    case = sys.argv[1] if len(sys.argv) > 1 else "case_Re10000_H6"
    cdir = os.path.join(here, case)
    data = np.genfromtxt(os.path.join(cdir, "tau_w_r.csv"), delimiter=",", names=True)
    r_mm = data["r_mm"]; r_over_d = data["r_over_d"]; tau = data["tau_w_Pa"]

    # info
    info = {}
    for line in open(os.path.join(cdir, "case_info.txt")):
        kv = line.strip().split("=")
        if len(kv) == 2:
            try: info[kv[0]] = float(kv[1])
            except ValueError: pass
    V0, Re = info["V0"], info["Re_d"]

    # fitted closure from sweep_results (match case+refine==1)
    C_P = n_t = None
    swp = os.path.join(here, "sweep_results.jsonl")
    if os.path.exists(swp):
        for line in open(swp):
            try: j = json.loads(line)
            except json.JSONDecodeError: continue
            if j.get("case") == case:
                C_P, n_t = j["C_P"], j["n_t"]
    if C_P is None:  # fall back to fitting here
        band = (r_over_d >= 2) & (r_over_d <= 20) & (tau > 0)
        LHS = (tau / (RHO * V0**2)) * np.sqrt(Re)
        slope, inter = np.polyfit(np.log(r_over_d[band]), np.log(LHS[band]), 1)
        n_t, C_P = -slope, float(np.exp(inter))

    r = r_mm * 1e-3
    # fitted turbulent law in Pa:  tau = C_P rho V0^2 Re^{-1/2} (r/d)^{-n_t}
    tau_fit = C_P * RHO * V0**2 * Re**(-0.5) * (r / D)**(-n_t)
    # laminar Watson II in Pa: tau = c2 rho V0^2 Re_a (a/r)^5, a=d/2, Re_a=Re/2
    a = D / 2.0; Re_a = Re / 2.0
    tau_lam = C2 * RHO * V0**2 * Re_a * (a / r)**5

    # restrict the fitted turbulent law to its valid fit band (r/d in [2,20])
    band = (r_over_d >= 2) & (r_over_d <= 20)
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.axvspan(2 * D * 1e3, 20 * D * 1e3, color="C1", alpha=0.06,
               label="power-law fit band")
    ax.loglog(r_mm, tau, "o", ms=3.5, color="C0", label="CFD (k-$\\omega$ SST)")
    ax.loglog(r_mm[band], tau_fit[band], "-", color="C1", lw=2.2,
              label=f"fitted turbulent  $r^{{-{n_t:.2f}}}$  ($C_P={C_P:.3g}$)")
    ax.loglog(r_mm, tau_lam, "--", color="C3", lw=1.5, label="laminar Watson II  $r^{-5}$")
    ax.axvline(15.0, ls=":", color="grey", label="$r_{max}=15$ mm")
    ax.set_xlabel("radius $r$ (mm)")
    ax.set_ylabel(r"wall shear stress $\tau_w$ (Pa)")
    ax.set_title(f"Impinging-jet wall shear: CFD vs turbulent fit vs laminar $r^{{-5}}$\n"
                 f"$Re_d={Re:.0f}$, $H/d={info['H_over_d']:.0f}$")
    ax.set_ylim(1e-3, max(tau.max(), tau_lam[np.argmin(np.abs(r_mm-2))]) * 2)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    out = os.path.join(here, "..", "fig_cfd_tau_w.png")
    fig.savefig(out, dpi=140)
    fig.savefig(os.path.join(here, "fig_cfd_tau_w.png"), dpi=140)
    print(f"wrote fig_cfd_tau_w.png  (n_t={n_t:.3f}, C_P={C_P:.4g})")


if __name__ == "__main__":
    main()
