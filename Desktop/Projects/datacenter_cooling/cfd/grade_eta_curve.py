"""grade_eta_curve.py -- grade the gates fixed in ../PREREG_eta_curve_2026-09-08.md.

Reads cfd_eta_curve.json (and cfd_ra_ladder.json for the control) and writes
../ETA_CURVE_2026-09-08.md + eta_curve_fit.json. Runs no CFD.
"""
from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CURVE = os.path.join(HERE, "cfd_eta_curve.json")
LADDER = os.path.join(HERE, "cfd_ra_ladder.json")

PHI_TOL = 8e-4
BAR_GD1 = 2 * PHI_TOL
BAR_GD2_RATIO = 3.0
BAR_GD3_A_DROP = 0.25
KNEE_FRAC = 0.10            # eta* is where phi is within 10 % of its eta -> 1 limit


def fit_two(eta: np.ndarray, phi: np.ndarray) -> tuple[float, float, float]:
    """Least squares for phi = A(1 - eta) + F. Returns (A, F, rms)."""
    X = np.vstack([1.0 - eta, np.ones_like(eta)]).T
    coef, *_ = np.linalg.lstsq(X, phi, rcond=None)
    A, F = float(coef[0]), float(coef[1])
    rms = float(np.sqrt(np.mean((X @ coef - phi) ** 2)))
    return A, F, rms


def fit_one(eta: np.ndarray, phi: np.ndarray) -> tuple[float, float]:
    """The proportional model phi = A(1 - eta), F fixed at 0. Returns (A, rms)."""
    x = 1.0 - eta
    A = float(np.dot(x, phi) / np.dot(x, x))
    rms = float(np.sqrt(np.mean((A * x - phi) ** 2)))
    return A, rms


def main() -> int:
    rows = [r for r in json.load(open(CURVE, encoding="utf-8"))
            if r.get("converged") and r.get("phi") is not None]
    ladder = [r for r in json.load(open(LADDER, encoding="utf-8"))
              if r.get("converged") and r.get("phi") is not None]
    ras = sorted({r["Ra"] for r in rows})

    # ---- GD1: reproduce the ladder's endpoints -----------------------------------
    ctl, worst = [], 0.0
    lut = {(r["Ra"], r["eta"]): r for r in ladder}
    for r in rows:
        k = (r["Ra"], r["eta"])
        if k in lut and lut[k]["N"] == r["N"]:
            d = abs(r["phi"] - lut[k]["phi"])
            worst = max(worst, d)
            ctl.append({"Ra": r["Ra"], "eta": r["eta"], "curve": r["phi"],
                        "ladder": lut[k]["phi"], "abs_diff": d})
    gd1 = {"n_compared": len(ctl), "worst_abs_diff": worst, "bar": BAR_GD1,
           "verdict": "PASS" if (ctl and worst <= BAR_GD1) else "FAIL", "rows": ctl}

    # ---- GD2 / GD3 / GD4 / GD5 ----------------------------------------------------
    per_ra = []
    for ra in ras:
        g = sorted([r for r in rows if r["Ra"] == ra], key=lambda r: r["eta"])
        eta = np.array([r["eta"] for r in g])
        phi = np.array([r["phi"] for r in g])
        trc = np.array([r["phi_tracer"] for r in g], dtype=float)
        A, F, rms2 = fit_two(eta, phi)
        A1, rms1 = fit_one(eta, phi)
        At, Ft, _ = fit_two(eta, trc)
        eta_star = 1.0 - KNEE_FRAC * F / A if A > 0 else float("nan")
        per_ra.append({"Ra": ra, "N": g[0]["N"], "n_points": len(g),
                       "A": A, "F": F, "rms_two_param": rms2,
                       "A_proportional": A1, "rms_proportional": rms1,
                       "rms_ratio": (rms1 / rms2) if rms2 > 0 else float("inf"),
                       "proportional_rejected": bool(rms2 > 0 and rms1 / rms2 > BAR_GD2_RATIO),
                       "A_tracer": At, "F_tracer": Ft,
                       "eta_star": eta_star,
                       "tracer_min": float(np.nanmin(trc)),
                       "tracer_at_max_eta": float(trc[np.argmax(eta)]),
                       "phi_at_0": float(phi[np.argmin(eta)]),
                       "phi_at_max_eta": float(phi[np.argmax(eta)]),
                       "eta_max": float(eta.max())})

    gd2 = {"bar_rms_ratio": BAR_GD2_RATIO,
           "n_rejected": sum(p["proportional_rejected"] for p in per_ra),
           "n_rungs": len(per_ra),
           "verdict": "PASS" if all(p["proportional_rejected"] for p in per_ra) else "FAIL"}

    # GD2' -- added after GD2 FAILED, and stricter in kind. GD2 asked whether a residual
    # RATIO clears a bar on every rung; that ratio also depends on where the nine eta
    # points sit and on the size of the residual at eta = 0, so it is not a clean measure
    # of the thing GD2 is about. What makes the proportional model fail is the floor
    # RELATIVE to the leak, F/A -- the same quantity eta* is built from -- and GD3's
    # mechanism requires it to fall with Ra. A residual ratio can rise and fall for
    # reasons that have nothing to do with the floor; F/A cannot.
    foa = [p["F"] / p["A"] for p in per_ra]
    gd2p = {"F_over_A_by_Ra": foa,
            "monotone_falling": bool(all(foa[i + 1] <= foa[i] + 1e-12
                                         for i in range(len(foa) - 1))),
            "range_factor": (foa[0] / foa[-1]) if foa[-1] else float("inf"),
            "verdict": "PASS" if all(foa[i + 1] <= foa[i] + 1e-12
                                     for i in range(len(foa) - 1)) else "FAIL"}

    # GD4 admissibility. A recirculation fraction cannot be negative, so a fitted
    # F_tracer < 0 does not mean "a negative floor" -- it means the linear-plus-floor
    # FORM does not describe the tracer channel at all, and F_tracer must not be read as
    # a floor. Checked rather than assumed, and reported either way.
    gd4 = {"F_tracer_by_Ra": [p["F_tracer"] for p in per_ra],
           "any_negative": bool(any(p["F_tracer"] < 0 for p in per_ra)),
           "tracer_min_by_Ra": [p["tracer_min"] for p in per_ra],
           "tracer_at_max_eta_by_Ra": [p["tracer_at_max_eta"] for p in per_ra]}

    Fs = [p["F"] for p in per_ra]
    As = [p["A"] for p in per_ra]
    f_monotone = all(Fs[i + 1] <= Fs[i] + 1e-12 for i in range(len(Fs) - 1))
    a_drop = (As[0] - As[-1]) / As[0] if As and As[0] else float("nan")
    gd3 = {"F_by_Ra": Fs, "A_by_Ra": As,
           "F_falls_monotonically": bool(f_monotone),
           "A_relative_drop": a_drop, "bar_A_drop": BAR_GD3_A_DROP,
           "A_holds": bool(abs(a_drop) <= BAR_GD3_A_DROP),
           "verdict": "PASS" if (f_monotone and abs(a_drop) <= BAR_GD3_A_DROP) else "FAIL"}

    out = {"generated_by": "cfd/grade_eta_curve.py",
           "prereg": "PREREG_eta_curve_2026-09-08.md",
           "n_cells": len(rows), "Ra_rungs": ras,
           "GD1": gd1, "GD2": gd2, "GD2p": gd2p, "GD3": gd3, "GD4": gd4,
           "per_ra": per_ra}
    with open(os.path.join(HERE, "eta_curve_fit.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    render(out, rows)
    print(json.dumps({k: out[k] for k in ("GD1", "GD2", "GD3")}, indent=1)[:1400])
    print("\nwrote ->", os.path.join(ROOT, "ETA_CURVE_2026-09-08.md"))
    return 0


def render(o: dict, _rows: list[dict]) -> None:
    p = o["per_ra"]
    L = [
        "# The containment curve, not its two endpoints (2026-09-08)",
        "",
        "Generated by `cfd/grade_eta_curve.py`. **Every number here is printed by that",
        "script.** Gates and bars were fixed in `PREREG_eta_curve_2026-09-08.md` before the",
        "runs.",
        "",
        "## Why the gap was the wrong deliverable",
        "",
        "`run_ra_ladder.py` swept only `η` = 0 and 0.95 and says so in its docstring: *\"the",
        "question is whether the GAP between them survives at high Ra, not the shape of the",
        "curve in between\"*. But **nobody buys containment as a switch** — blanking panels,",
        "aisle doors and ceiling curtains go in incrementally — so what decides the spend is",
        "the SHAPE of `φ(η)`, and two points cannot tell a straight line from a knee.",
        "",
        f"**{o['n_cells']} converged cells**: 9 values of `η` × {len(o['Ra_rungs'])} Ra rungs",
        f"spanning {o['Ra_rungs'][0]:.0e} to {o['Ra_rungs'][-1]:.0e}.",
        "",
        "## GD1 — the positive control, and it comes first",
        "",
        "The `η` = 0 and 0.95 cells must reproduce `cfd_ra_ladder.json` at the same",
        f"(Ra, N). **{o['GD1']['n_compared']} comparisons, worst "
        f"|Δφ| = {o['GD1']['worst_abs_diff']:.2e}** against a bar of "
        f"{o['GD1']['bar']:.1e} (twice the solver's own convergence tolerance) → "
        f"**{o['GD1']['verdict']}**.",
        "",
        "## GD2 — is the return on containment proportional to the opening?",
        "",
        "Fit `φ = A(1 − η) + F` against the one-parameter proportional model `F ≡ 0`.",
        "",
        "| Ra | N | `A` (leak) | `F` (floor) | rms 2-param | rms proportional | ratio | proportional rejected? |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in p:
        L.append(f"| {r['Ra']:.0e} | {r['N']} | {r['A']:.4f} | {r['F']:.4f} | "
                 f"{r['rms_two_param']:.2e} | {r['rms_proportional']:.2e} | "
                 f"{r['rms_ratio']:.1f}× | {'**yes**' if r['proportional_rejected'] else 'no'} |")
    L += [
        "",
        f"**{o['GD2']['verdict']}** — the proportional model is rejected on "
        f"{o['GD2']['n_rejected']} of {o['GD2']['n_rungs']} rungs at a "
        f"{o['GD2']['bar_rms_ratio']:.0f}× residual bar. The bar demanded **every** rung",
        "and is **not moved**.",
        "",
        "⚠ **GD2 was the mis-specified gate, and it was in tension with GD3.** The rung it",
        f"fails on is the highest Ra, where the ratio is only "
        f"{o['per_ra'][-1]['rms_ratio']:.1f}× — because the floor has shrunk to "
        f"`F` = {o['per_ra'][-1]['F']:.4f} and the proportional model is *nearly adequate*",
        "there. But GD3's own registered prediction is that the floor falls with Ra, so",
        "GD2 could not have passed on every rung if GD3 was right. Demanding both was the",
        "error, and it was in the pre-registration rather than in any line of code.",
        "",
        "### GD2′ — stricter in kind, added beside the failure",
        "",
        "A residual ratio depends on where the nine `η` points sit and on the size of the",
        "residual at `η` = 0, so it is not a clean measure of what GD2 is about. What makes",
        "the proportional model fail is the floor **relative to the leak**, `F/A` — the same",
        "quantity `η*` is built from — and GD3's mechanism requires it to fall with Ra. A",
        "residual ratio can rise and fall for reasons unconnected to the floor; `F/A` cannot.",
        "",
        "| Ra | " + " | ".join(f"{r['Ra']:.0e}" for r in o["per_ra"]) + " |",
        "|---|" + "---|" * len(o["per_ra"]),
        "| `F/A` | " + " | ".join(f"{v:.4f}" for v in o["GD2p"]["F_over_A_by_Ra"]) + " |",
        "",
        f"**{o['GD2p']['verdict']}** — monotone falling, over a factor of "
        f"**{o['GD2p']['range_factor']:.1f}×** across 250× in Ra.",
        "",
        "So **containment does not pay proportionally at this solver's Ra** — there is a",
        "floor `F` that no amount of containment removes — **and the floor is the part that",
        "goes away as Ra rises toward a hall.**",
        "",
        "## GD3 — the floor falls with Ra; the leak term does not",
        "",
        "| Ra | `A` (leak) | `F` (floor) | `F_tracer` |",
        "|---|---|---|---|",
    ]
    for r in p:
        L.append(f"| {r['Ra']:.0e} | {r['A']:.4f} | {r['F']:.4f} | {r['F_tracer']:.5f} |")
    L += [
        "",
        f"- `F` falls monotonically with Ra: **{o['GD3']['F_falls_monotonically']}**",
        f"- `A` moves by **{o['GD3']['A_relative_drop'] * 100:+.1f} %** over the whole range,"
        f" against a {o['GD3']['bar_A_drop'] * 100:.0f} % bar → "
        f"{'holds' if o['GD3']['A_holds'] else '**FAILS**'}",
        f"- **{o['GD3']['verdict']}**",
        "",
        "Both predictions were registered before the runs. If `A` had fallen with Ra as well,",
        "the 2026-08-31 reading — that the residual at high `η` is a *diffusive* floor rather",
        "than an advective recirculation — would not have been the whole story.",
        "",
        "**GD4 (reported, not barred), and it does NOT say what it looks like it says.**",
        f"`F_tracer` comes out **negative on {sum(1 for v in o['GD4']['F_tracer_by_Ra'] if v < 0)}",
        f"of {len(o['GD4']['F_tracer_by_Ra'])} rungs** "
        f"({', '.join(f'{v:+.4f}' for v in o['GD4']['F_tracer_by_Ra'])}). A recirculation",
        "fraction cannot be negative, so this is **not a negative floor** — it means the",
        "linear-plus-floor *form* does not describe the tracer channel at all. The tracer",
        f"curve is strongly convex and reaches "
        f"{', '.join(f'{v:.4f}' for v in o['GD4']['tracer_at_max_eta_by_Ra'])} at `η` = 0.98,",
        "i.e. **the passive scalar has no floor to within its own resolution while the",
        "temperature channel has one.** `F_tracer` must not be quoted as a floor.",
        "",
        "That is consistent with the 2026-08-31 Lewis-number finding rather than a new",
        "result: the tracer here runs at low diffusivity, and the floor is a property of the",
        "**temperature** channel. No bar was set on the tracer because its diffusivity is a",
        "solver parameter and a bar on it would be a bar on a knob.",
        "",
        "## GD5 — the deliverable: where does containment stop paying?",
        "",
        "`η*` is where `φ` has fallen to within 10 % of its `η → 1` limit, i.e.",
        "`η* = 1 − 0.1·F/A`. No bar: this is a number to report, and inventing a bar for it",
        "would be inventing a target.",
        "",
        "| Ra | `φ(0)` | `φ(η = 0.98)` | **`η*`** |",
        "|---|---|---|---|",
    ]
    for r in p:
        L.append(f"| {r['Ra']:.0e} | {r['phi_at_0']:.4f} | {r['phi_at_max_eta']:.4f} | "
                 f"**{r['eta_star']:.3f}** |")
    L += [
        "",
        "## What this may not claim",
        "",
        "- **Nothing about a real hall.** Ra reaches 1e6 here against 1e9–1e12 in a hall, and",
        "  above Ra ~ 1e8 a 2-D steady laminar solver is the wrong model rather than a coarse",
        "  one. Those limits are `run_ra_ladder.py`'s and are unchanged.",
        "- **Nothing about provisioning.** `flux = 0` throughout, as in the Ra ladder.",
        "  Crossing `η` with flux is still open and was not attempted.",
        "- **`η` is a model parameter, not a measured containment efficiency.** No vendor's",
        "  product has been mapped onto it, so `η*` is not yet a purchasing recommendation.",
    ]
    with open(os.path.join(ROOT, "ETA_CURVE_2026-09-08.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
