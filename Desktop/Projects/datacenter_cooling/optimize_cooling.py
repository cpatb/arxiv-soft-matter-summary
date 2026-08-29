"""
optimize_cooling.py — minimize PUE over (supply temperature, airflow provisioning)
subject to the ASHRAE rack-inlet constraint, on a generic hall. Produces the
headline figures and prints the optimum.

Decision variables:  T_supply [C], provisioning p = supplied/demanded airflow.
Objective:           PUE (total facility power / IT power).
Constraints:         hottest rack inlet <= T_allow ;  p >= p_min.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dc_thermal import generic_hall, Hall


def optimize(hall: Hall, T_bounds=(14.0, 30.0), p_bounds=(0.6, 1.4), x0=(18.0, 1.0),
             objective="pue"):
    """Minimize PUE, or total facility power, subject to the ASHRAE inlet limit.

    The two are not the same problem. PUE is total/IT, and the server-fan ramp adds
    to BOTH, so a design that makes the racks hotter can lower the ratio while
    raising the hall's energy bill. `objective="total"` ranks by energy instead.
    """
    if objective not in ("pue", "total"):
        raise ValueError(f"objective must be 'pue' or 'total', got {objective!r}")

    def obj(x):
        return hall.pue(x[0], x[1]) if objective == "pue" else hall.total_kw(x[0], x[1])

    cons = [
        {"type": "ineq", "fun": lambda x: hall.T_allow - hall.rack_inlet_temps(x[0], x[1]).max()},
        {"type": "ineq", "fun": lambda x: x[1] - hall.p_min},
    ]
    res = minimize(obj, np.array(x0, float), method="SLSQP",
                   bounds=[T_bounds, p_bounds], constraints=cons,
                   options={"ftol": 1e-9, "maxiter": 300})
    return res


def containment_sweep(total_kw=200.0, etas=None):
    """Free-airflow sweep: optimize over both (T_s, p) at each containment level."""
    etas = np.linspace(0.0, 1.0, 11) if etas is None else np.asarray(etas)
    pue, Ts, ps = [], [], []
    for eta in etas:
        h = generic_hall(containment=float(eta), total_kw=total_kw, seed=1)
        r = optimize(h)
        pue.append(r.fun); Ts.append(r.x[0]); ps.append(r.x[1])
    return etas, np.array(pue), np.array(Ts), np.array(ps)


def containment_roi(total_kw=200.0, etas=None, provisioning=1.0):
    """Containment ROI at FIXED airflow -- the curve to report.

    The free-airflow sweep above is degenerate: at its optimum the mass-balance
    floor (1-p) exceeds the geometric term phi(eta)*(1+under), so eta drops out
    of B and the sweep is flat. That is a statement about the airflow axis, not
    about containment. Holding p at the design point isolates the containment
    lever, and keeps the model inside the range it is valid over (see DESIGN.md
    Sec. 8: the p > 1 branch carries no thermal benefit, so the airflow axis is
    one-sided and its optimum is not yet trustworthy).
    """
    etas = np.linspace(0.0, 1.0, 11) if etas is None else np.asarray(etas)
    pue, Ts = [], []
    for eta in etas:
        h = generic_hall(containment=float(eta), total_kw=total_kw, seed=1)
        # raise supply until the hottest rack sits exactly on the ASHRAE limit
        T_s = h.T_allow - h.rack_inlet_temps(0.0, provisioning).max()
        Ts.append(T_s); pue.append(h.pue(T_s, provisioning))
    return etas, np.array(pue), np.array(Ts)


def main():
    hall = generic_hall(n_racks=20, total_kw=200.0, containment=0.3, seed=1)
    base = hall.evaluate(18.0, 1.0)          # a naive "cold and over-blown" baseline
    res = optimize(hall)
    opt = hall.evaluate(res.x[0], res.x[1])

    print("=" * 68)
    print("GENERIC HALL — 20 racks, 200 kW IT, containment eta=0.30")
    print("=" * 68)
    print(f"Baseline (T_s=18.0C, p=1.00):  PUE={base['PUE']:.3f}  "
          f"hottest inlet={base['hottest_C']:.1f}C  cooling={base['cooling_kW']:.1f}kW")
    print(f"Optimized:  T_s={opt['T_supply']:.1f}C  p={opt['provisioning']:.2f}  "
          f"PUE={opt['PUE']:.3f}  hottest inlet={opt['hottest_C']:.1f}C "
          f"(limit {hall.T_allow:.0f}C)")
    saved = (base['cooling_kW'] - opt['cooling_kW']) / base['cooling_kW'] * 100
    print(f"Cooling-energy reduction vs baseline: {saved:.1f}%   "
          f"(COP {base['COP']:.2f} -> {opt['COP']:.2f})")
    print(f"Server-fan draw at the optimum: {opt['server_fan_kW']:.2f} kW "
          f"({100 * opt['server_fan_kW'] / opt['IT_nominal_kW']:.2f}% of nominal IT); "
          f"total hall power {opt['total_kW']:.1f} kW")

    # ---- Is PUE a safe objective here? ----
    # Server-fan power sits inside PUE's denominator, so a design that runs the
    # racks hotter can lower the ratio while raising the bill. Whether that
    # actually bites depends on the inlet limit and on the ramp coefficient, so
    # report it rather than assert it.
    print()
    print("Is PUE a safe objective? (PUE optimum vs energy optimum, extra kW burned)")
    print(f"  {'ramp %/K':>9s} " + "".join(f"{f'T_allow {t:.0f}':>13s}" for t in (27, 32, 35, 40)))
    for ramp in (0.0035, 0.004, 0.005, 0.0075, 0.010):
        cells = []
        for T_allow in (27.0, 32.0, 35.0, 40.0):
            hh = generic_hall(n_racks=20, total_kw=200.0, containment=0.3, seed=1,
                              T_allow=T_allow, fan_ramp=ramp)
            p_opt = hh.evaluate(*optimize(hh, objective="pue").x)
            e_opt = hh.evaluate(*optimize(hh, objective="total").x)
            # The energy optimum minimises total power, so this gap cannot truly be
            # negative; a small negative is SLSQP noise and would print as "-0.00".
            gap = max(0.0, p_opt["total_kW"] - e_opt["total_kW"])
            cells.append(f"{gap:10.2f} kW")
        print(f"  {100 * ramp:9.2f} " + "".join(f"{c:>13s}" for c in cells))
    print("  Negligible at the mid-band 0.40 %/K and the 27 C recommended limit;")
    print("  material at the top of the published range on a relaxed envelope.")

    # ---- Figure 1: rack inlet temperature map, baseline vs optimized ----
    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(hall.n)
    ax.bar(x - 0.2, base["rack_inlet_C"], 0.4, label=f"baseline (PUE {base['PUE']:.2f})")
    ax.bar(x + 0.2, opt["rack_inlet_C"], 0.4, label=f"optimized (PUE {opt['PUE']:.2f})")
    ax.axhline(hall.T_allow, ls="--", color="k", lw=1, label=f"ASHRAE limit {hall.T_allow:.0f}C")
    ax.set_xlabel("rack index"); ax.set_ylabel("inlet temperature [C]")
    ax.set_title("Rack-inlet temperatures: the optimizer runs hot racks right to the limit")
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig("fig1_inlet_map.png", dpi=130)

    # ---- Figure 2: PUE vs containment effectiveness (the design lever) ----
    # Reported at fixed airflow. The free-airflow sweep is overlaid because it is
    # flat, and that flatness is the honest caveat rather than something to hide:
    # it says the airflow axis, not containment, dominates the free optimum.
    # PRL style: no titles, short legend, explanation in fig2_containment_sweep.md.
    # PUE is anchored at its physical floor of 1.0 (zero cooling overhead) rather
    # than at the data range -- the containment lever really is ~1% of PUE, and a
    # cut axis would show it as though it filled the panel.
    etas, pue, Ts = containment_roi()
    _, pue_free, _, _ = containment_sweep()
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.1, 2.9))
    axa.plot(etas, pue, "o-", color="C3", ms=4, label="fixed airflow")
    axa.plot(etas, pue_free, "^:", color="grey", lw=1.2, ms=4, label="free airflow")
    axa.set_xlabel(r"containment effectiveness  $\eta$")
    axa.set_ylabel("PUE")
    axa.set_xlim(0, 1.0)
    axa.set_ylim(1.0, max(pue.max(), pue_free.max()) * 1.06)
    axa.legend(frameon=False, fontsize=8, loc="lower left")
    axa.text(0.03, 0.93, "(a)", transform=axa.transAxes, fontsize=9)

    axb.plot(etas, Ts, "s-", color="C0", ms=4)
    axb.set_xlabel(r"containment effectiveness  $\eta$")
    axb.set_ylabel("supply temperature (°C)")
    axb.set_xlim(0, 1.0)
    axb.set_ylim(0, Ts.max() * 1.06)
    axb.text(0.03, 0.93, "(b)", transform=axb.transAxes, fontsize=9)

    for ax in (axa, axb):
        ax.tick_params(direction="in", top=True, right=True, labelsize=8)
        ax.xaxis.label.set_size(9); ax.yaxis.label.set_size(9)
    fig.tight_layout(); fig.savefig("fig2_containment_sweep.png", dpi=200)

    print(f"\nContainment ROI at fixed airflow p=1: PUE {pue[0]:.4f} (eta=0) -> "
          f"{pue[-1]:.4f} (eta=1), dPUE {pue[0]-pue[-1]:.4f}; "
          f"supply T {Ts[0]:.1f} -> {Ts[-1]:.1f} C")
    print(f"Free-airflow sweep: PUE {pue_free[0]:.4f} -> {pue_free[-1]:.4f} "
          f"(dPUE {pue_free[0]-pue_free[-1]:.4f}) -- flat, mass-balance floor binds")
    print("Figures: fig1_inlet_map.png, fig2_containment_sweep.png")


if __name__ == "__main__":
    main()
