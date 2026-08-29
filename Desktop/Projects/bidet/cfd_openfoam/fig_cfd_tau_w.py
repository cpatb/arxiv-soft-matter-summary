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
               label="fit band")
    ax.loglog(r_mm, tau, "o", ms=3.5, color="C0", label="CFD (k-$\\omega$ SST)")
    ax.loglog(r_mm[band], tau_fit[band], "-", color="C1", lw=2.2,
              label=f"turbulent fit, $r^{{-{n_t:.2f}}}$")
    ax.loglog(r_mm, tau_lam, "--", color="C3", lw=1.5, label="laminar, $r^{-5}$")
    ax.axvline(15.0, ls=":", color="grey", label="$r_{max}$")
    # Mark where the computation actually ends, and where it used to end. Without the
    # first the plot reads as if the curve simply stopped; without the second a reader
    # cannot see that the tau_c = 1 Pa contour used to lie outside the domain.
    Rdom_mm = info["Rdom"] * 1e3
    ax.axvline(Rdom_mm, ls="-", color="k", lw=1.3, alpha=0.7,
               label=f"$R={Rdom_mm:.0f}$ mm")
    if Rdom_mm > 26.0:
        ax.axvline(25.0, ls=":", color="k", lw=1.1, alpha=0.6,
                   label="earlier $R$")
    # tau_c = 1 Pa: the level the narrow domain could not reach.
    ipk = int(np.argmax(tau))
    ro, to = r_mm[ipk:], tau[ipk:]
    if to.min() <= 1.0:
        root = float(np.interp(1.0, to[::-1], ro[::-1]))
        ax.axhline(1.0, ls="-.", color="C2", lw=1.3,
                   label=rf"$\tau_c=1$ Pa, crossed at {root:.1f} mm")
        ax.plot([root], [1.0], "*", ms=11, color="C2", zorder=5)
    else:
        ax.axhline(1.0, ls="-.", color="C2", lw=1.3,
                   label=r"$\tau_c=1$ Pa (not reached in domain)")
    ax.set_xlabel("radius $r$ (mm)")
    ax.set_ylabel(r"wall shear stress $\tau_w$ (Pa)")
    # No title: the manuscript caption carries Re_d, H/d and the fit parameters, and a
    # figure that needs a paragraph printed on it is a caption pasted into the axes.
    # The log ordinate is kept because tau_w spans ~3 decades (0.455 to 389 Pa); the
    # abscissa spans ~3.5. See fig_cfd_tau_w.md.
    ax.set_ylim(0.2, tau.max() * 2.0)
    ax.set_xlim(r_mm.min() * 0.9, info["Rdom"] * 1e3 * 1.15)
    ax.legend(fontsize=7.5, loc="lower left", frameon=False)
    ax.tick_params(direction="in", which="both", top=True, right=True)
    fig.tight_layout()
    out = os.path.join(here, "..", "fig_cfd_tau_w.png")
    fig.savefig(out, dpi=140)
    fig.savefig(os.path.join(here, "fig_cfd_tau_w.png"), dpi=140)
    print(f"wrote fig_cfd_tau_w.png  (n_t={n_t:.3f}, C_P={C_P:.4g})")


if __name__ == "__main__":
    main()
