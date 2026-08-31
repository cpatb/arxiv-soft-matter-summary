"""Fit phi(eta, p) from the two-axis CFD sweep, and separate advection from conduction.

WHAT THIS REPLACES
    `calibrate.py` fitted phi(eta) = floor + amp*exp(-eta/tau) to a sweep run at ZERO
    forced airflow, and `dc_thermal` then applies it at DESIGN airflow. The result was
    a model in which over-provisioning cannot help: `under = max(0, 1-p)` makes B
    identical at p = 1 and p = 1.4. DESIGN.md sec.8 calls that gap (1).

    Tonight's sweep (`cfd_provisioning.json`, 48 cells, 0 diverged) puts p on the map.

WHAT THE SWEEP SAYS
    Over-provisioning helps, and by a lot: at eta = 0 the recirculated fraction falls
    from 0.241 at p = 0 to 0.044 at p = 1 and 0.023 at p = 2. The current model gives
    exactly zero benefit over that whole range.

AND WHAT IT SAYS ABOUT THE OLD CALIBRATION
    `CFD_RECIRC` puts phi(eta = 1) at 0.0713 and DESIGN.md builds the "containment ROI
    flattens, ~7% floor" conclusion on it. `run_floor_probe.py` shows that floor is
    thermal DIFFUSION, not recirculated air, by sweeping the tracer's diffusivity: at
    eta = 0.95 the tracer signal scales with its own diffusivity all the way to Le = 1,
    where it reproduces the temperature channel to five decimals (0.06945 vs 0.0694),
    while at eta = 0 the same sweep moves it by 0.4% (ratio 1.000 -> 1.004). An advected
    signal is Lewis-independent; a diffused one is proportional to the diffusivity. The
    floor is entirely the second kind.

    This cell runs at Ra = 4e3 to 1e5. A hall is Ra ~ 1e9 to 1e12, where molecular
    conduction across a cold aisle is nothing. So the floor does not transfer, and this
    script reports both channels rather than picking one quietly.

Run:  python calibrate_provisioning.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
OUT = HERE / "cfd_provisioning_fit.json"


def model(X, phi0, k, tau, floor):
    """phi(eta, p) = [phi0 * exp(-eta/tau) - floor] / (1 + p/k) + floor.

    Three pieces, each with a reason rather than a free shape:

    * `exp(-eta/tau)` is the containment falloff `calibrate.py` already fitted and the
      sweep reproduces along its p = 0 edge, so the eta axis is not re-parameterised.
    * `1 / (1 + p/k)` is a dilution law. The inlet stream is a mixture of supply air and
      whatever the buoyant loop carries back; if the loop's strength is set by buoyancy
      and does not care how much supply air arrives, the recirculated FRACTION goes like
      one over the total, which is `1 + p/k` with `k` the loop strength in units of Q_1.
      It is a one-parameter law, not a spline, so it can fail -- and it does; see below.
    * `floor` is additive and p-independent, because the floor probe says it is
      conduction and conduction does not care about the through-flow.
    """
    eta, p = X
    return (phi0 * np.exp(-eta / tau) - floor) / (1.0 + p / k) + floor


def fit(channel: str, rows: list[dict]):
    eta = np.array([r["eta"] for r in rows])
    p = np.array([r["p"] for r in rows])
    y = np.array([r[channel] for r in rows])
    ok = np.isfinite(y)
    (phi0, k, tau, floor), cov = curve_fit(
        model, (eta[ok], p[ok]), y[ok],
        p0=[0.25, 0.3, 0.8, 0.01],
        bounds=([0.01, 0.01, 0.05, 0.0], [1.0, 10.0, 5.0, 0.2]),
        maxfev=40000,
    )
    pred = model((eta[ok], p[ok]), phi0, k, tau, floor)
    resid = y[ok] - pred
    return dict(phi0=float(phi0), k=float(k), tau=float(tau), floor=float(floor),
                rms=float(np.sqrt(np.mean(resid**2))),
                max_abs=float(np.max(np.abs(resid))),
                n=int(ok.sum()),
                stderr=[float(v) for v in np.sqrt(np.diag(cov))])


def main() -> None:
    rows = [r for r in json.load(open(HERE / "cfd_provisioning.json")) if not r["diverged"]]
    grid_eta = sorted({r["eta"] for r in rows})
    grid_p = sorted({r["p"] for r in rows})
    print(f"{len(rows)} converged cells, Ra = {rows[0]['Ra']:.0e}, N = {rows[0]['N']}, "
          f"Q_1 = {rows[0]['Q1']:.4f}")
    print()

    out = {}
    for channel, label in [("phi", "temperature  (advection + conduction)"),
                           ("phi_tracer", "tracer Le=50 (advection only)")]:
        f = fit(channel, rows)
        out[channel] = f
        keys = ("phi0", "k", "tau", "floor")
        pars = {q: f[q] for q in keys}
        print(label)
        print(f"    phi(eta, p) = [{f['phi0']:.4f} exp(-eta/{f['tau']:.4f}) - "
              f"{f['floor']:.4f}] / (1 + p/{f['k']:.4f}) + {f['floor']:.4f}")
        print(f"    rms residual {f['rms']:.4f}, worst {f['max_abs']:.4f}, "
              f"on {f['n']} cells")
        print(f"    phi(0, 0) = {model((0., 0.), **pars):.4f}   "
              f"phi(0, 1) = {model((0., 1.), **pars):.4f}   "
              f"phi(0.95, 1) = {model((0.95, 1.), **pars):.4f}")
        print()

    # --- why the closed form is not shipped -------------------------------------
    f = out["phi"]
    worst_ratio, worst_cell = 0.0, None
    for r in rows:
        pred = model((r["eta"], r["p"]), f["phi0"], f["k"], f["tau"], f["floor"])
        if r["phi"] > 0 and pred > 0:
            ratio = max(pred / r["phi"], r["phi"] / pred)
            if ratio > worst_ratio:
                worst_ratio, worst_cell = ratio, r

    top_eta0 = [r["phi"] for r in rows if r["eta"] == 0.0 and r["p"] == 4.0][0]
    print("why the closed form is not shipped")
    print(f"    `floor` fitted to {f['floor']:.6f}, which is its lower bound. A parameter")
    print("    sitting on a bound is not a measurement, and here it is being asked to")
    print("    absorb an eta-DEPENDENT conduction path: at eta = 0 there is no partition,")
    print("    so the hot wall conducts straight to the inlet face and the large-p limit")
    print(f"    is {top_eta0:.4f}; at eta = 0.95 the partition blocks it and the same limit")
    print("    is 0.0000. One additive constant cannot be both.")
    if worst_cell is not None:
        pred = model((worst_cell["eta"], worst_cell["p"]),
                     f["phi0"], f["k"], f["tau"], f["floor"])
        print(f"    Worst cell: eta = {worst_cell['eta']}, p = {worst_cell['p']} -- "
              f"measured {worst_cell['phi']:.4f}, predicted {pred:.4f}, a factor "
              f"of {worst_ratio:.0f}.")
    print("    So dc_thermal takes the measured 6 x 8 grid and interpolates it; the law")
    print("    above is kept as the record of what was tried, not as the calibration.")
    print()

    # --- the conduction share, read where the tracer is meaningful ---------------
    print("the conduction share, at p = 0 where the tracer is not flushed to zero")
    print()
    print(f"  {'eta':>5}  {'phi_T':>8}  {'phi_tracer':>11}   phi_tracer / phi_T")
    for eta in grid_eta:
        r = next(x for x in rows if x["eta"] == eta and x["p"] == 0.0)
        print(f"  {eta:>5}  {r['phi']:8.4f}  {r['phi_tracer']:11.4f}   "
              f"{r['phi_tracer'] / r['phi']:>10.3f}")
    print()
    print("  The Le = 50 tracer is a LOWER BOUND on the advective share, not a")
    print("  measurement of it: it under-reads any diffusive path by 50x by")
    print("  construction. That is exactly why run_floor_probe.py sweeps the")
    print("  diffusivity instead of trusting this column.")

    old_floor = 0.007 + 0.228 * np.exp(-1 / 0.790)
    p1_top = next(x["phi"] for x in rows if x["eta"] == 0.95 and x["p"] == 1.0)
    print()
    print(f"the old calibration: phi(eta=1) = {old_floor:.4f}, measured at p = 0 only, and")
    print(f"  applied by dc_thermal at DESIGN airflow, where this sweep measures {p1_top:.4f}.")

    out["grid"] = {
        "eta": grid_eta,
        "p": grid_p,
        "phi": [[next(x["phi"] for x in rows if x["eta"] == e and x["p"] == q)
                 for q in grid_p] for e in grid_eta],
        "phi_tracer": [[next(x["phi_tracer"] for x in rows
                             if x["eta"] == e and x["p"] == q)
                        for q in grid_p] for e in grid_eta],
    }
    out["closed_form_rejected"] = {
        "reason": f"floor pinned on its lower bound; worst cell off by {worst_ratio:.0f}x",
        "worst_cell": ({"eta": worst_cell["eta"], "p": worst_cell["p"],
                        "measured": worst_cell["phi"]} if worst_cell else None),
    }
    out["provenance"] = {
        "source": "cfd_provisioning.json",
        "n_cells": len(rows),
        "Ra": rows[0]["Ra"],
        "N": rows[0]["N"],
        "Q1": rows[0]["Q1"],
        "p_anchor": "flux / Q_1, Q_1 = max|psi| of the buoyancy-only solution at eta = 0",
        "caveat": ("p = 1 is anchored on the cell's own buoyant circulation, not on any "
                   "real fan curve. The SHAPE of phi(p) transfers; the location of p = 1 "
                   "does not."),
    }
    json.dump(out, open(OUT, "w"), indent=2)
    print()
    print(f"written {OUT.name}")


if __name__ == "__main__":
    main()
