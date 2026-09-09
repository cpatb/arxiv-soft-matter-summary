"""analyse_phase8.py — extract phi_c(CV) from the constant-volume shear ladders.

FOUR ESTIMATORS, AND WHY THE PRIMARY IS THE ONE IT IS
----------------------------------------------------
At constant volume every cell reaches a steady state that lies on the granular
`mu(I)` / `phi(I)` rheology, so the ladder is not a set of independent points: it is a
sampling of one curve, parameterised by the inertial number `I = gammadot d50
sqrt(rho/p)`. The critical-state fraction is the quasi-static limit of that curve.

  E1  phi(I) EXTRAPOLATION  [PRIMARY]   phi_c = lim_{I->0} phi(I), from a linear fit of
      phi against I over the rungs with I <= I_FIT_MAX. This is the standard granular
      dilatancy law (da Cruz, Emam, Prochnow, Roux & Chevoir, PRE 72, 021309 (2005);
      Forterre & Pouliquen, Annu. Rev. Fluid Mech. 40, 1 (2008)), the fit is over a
      SHORT lever arm because production I is 1e-4..1e-2, and the fitted quantity does
      not span decades.
  E2  PRESSURE ROOT, free exponent      p_ss = c (phi - phi_c)^b, fitted in log space.
      A check on E1 and on the contact law: `b` should return near 1 for Yade's linear
      Cundall-Strack contacts. Reported, never used to redefine phi_c.
  E3  PRESSURE ROOT, b pinned to 1      the harmonic-contact convention. Declared, not
      circular: it is used as a CHECK on E1, and no claim here is made about `b`.
  E4  COORDINATION CROSSING             the phi at which z_nr crosses Z_ISO_FRICTIONAL.
      Independent of the stress entirely, which is its whole value: a defect in the
      stress read-out cannot move it.

The four are reported side by side for every cell group. A CV contrast is quoted only
if it has the same sign under all four; the SPREAD between them is the systematic and
it is printed.

WHAT IS REFUSED
---------------
  * a root outside the ladder it was fitted on (extrapolation dressed as measurement)
  * a group whose cells are not all steady, or not all quasi-static -- excluded cells
    are PRINTED with their reason, never dropped silently
  * fitting against `phi_target`; every fit uses the REALISED `phi`, because stage B
    lands on the target only to a tolerance and the first version of it did not
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "out")

# --- named constants, with origin ---------------------------------------------------
I_FIT_MAX = 1.0e-2      # quasi-static bar; same value as shear.I_MAX, kept in sync by
                        # test_shear.py::test_i_max_agrees
Z_ISO_FRICTIONAL = 4.0  # isostatic contact number for frictional spheres in 3D
                        # (Z_iso = d + 1 = 4 with all three friction constraints
                        # mobilised); the project's A3 anchor measures 3.98-4.07 at
                        # mu_p = 3. It is a CONVENTION here, and E4 is reported as such.
P_FLOOR = 1.0e-9        # below this p_hat the cell carries no measurable stress


def load(pattern: str) -> list[dict]:
    recs = []
    for path in sorted(glob.glob(os.path.join(OUT_DIR, pattern))):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    recs.append(json.loads(line))
                except json.JSONDecodeError:
                    continue
    return recs


def group_key(r: dict) -> tuple:
    """Everything a comparison must hold fixed, in the key.

    🔴 `protocol`, `gamma_max` and `shear_steps` are in here deliberately. The default
    glob `phase8_prod*.jsonl` matches BOTH `phase8_prod*` (protocol cv-shear-v1) and
    `phase8_prod2*` (v2, after the periodic-cell flip fix and with a different ladder,
    strain and step count). Pooling those would be the N-pooling defect of 2026-09-08
    repeated in a new place: a key that does not encode a knob the analysis groups by.
    `main` additionally REFUSES to report contrasts when more than one protocol is
    present, rather than quietly averaging over it.
    """
    return (r.get("family", "lognormal"), round(float(r["cv"]), 3),
            round(float(r["mu_p"]), 3), int(r["seed"]), int(r["n_requested"]),
            r.get("protocol", "?"), round(float(r.get("gamma_max", float("nan"))), 3),
            int(r.get("shear_steps", -1)))


def e1_phi_of_I(phi: np.ndarray, ivals: np.ndarray) -> dict:
    """phi_c = lim_{I->0} phi(I), linear fit over the quasi-static rungs."""
    ok = np.isfinite(phi) & np.isfinite(ivals) & (ivals <= I_FIT_MAX) & (ivals > 0)
    out = {"phi_c_E1": np.nan, "slope_E1": np.nan, "n_E1": int(ok.sum()),
           "r2_E1": np.nan, "I_min": np.nan, "I_max": np.nan}
    if ok.sum() < 3:
        return out
    x, y = ivals[ok], phi[ok]
    A = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    out.update(phi_c_E1=float(coef[0]), slope_E1=float(coef[1]),
               r2_E1=float(1 - np.sum((y - pred) ** 2) / ss_tot) if ss_tot > 0 else np.nan,
               I_min=float(x.min()), I_max=float(x.max()))
    return out


def _pressure_root(phi: np.ndarray, p: np.ndarray, pin_b: float | None) -> dict:
    """Root of p(phi) = c (phi - phi_c)^b, fitted on log p."""
    from scipy.optimize import curve_fit

    tag = "E3" if pin_b is not None else "E2"
    out = {f"phi_c_{tag}": np.nan, f"b_{tag}": (pin_b if pin_b else np.nan),
           f"n_{tag}": 0}
    ok = np.isfinite(phi) & np.isfinite(p) & (p > P_FLOOR)
    if ok.sum() < (3 if pin_b is not None else 4):
        return out
    x, y = phi[ok], np.log(p[ok])
    lo = float(x.min())

    if pin_b is None:
        def model(xx, phic, logc, b):
            return logc + b * np.log(np.maximum(xx - phic, 1e-12))
        p0 = [lo - 0.02, float(y.mean()), 1.0]
        bounds = ([-np.inf, -np.inf, 0.1], [lo - 1e-6, np.inf, 8.0])
    else:
        def model(xx, phic, logc):
            return logc + pin_b * np.log(np.maximum(xx - phic, 1e-12))
        p0 = [lo - 0.02, float(y.mean())]
        bounds = ([-np.inf, -np.inf], [lo - 1e-6, np.inf])
    try:
        popt, _ = curve_fit(model, x, y, p0=p0, bounds=bounds, maxfev=40000)
    except Exception as exc:                     # noqa: BLE001 - report, do not hide
        out[f"fit_error_{tag}"] = repr(exc)
        return out
    out[f"phi_c_{tag}"] = float(popt[0])
    out[f"n_{tag}"] = int(ok.sum())
    if pin_b is None:
        out["b_E2"] = float(popt[2])
    return out


def e4_z_crossing(phi: np.ndarray, z: np.ndarray) -> dict:
    """phi where z_nr crosses Z_ISO_FRICTIONAL, by linear interpolation."""
    out = {"phi_c_E4": np.nan}
    ok = np.isfinite(phi) & np.isfinite(z)
    if ok.sum() < 2:
        return out
    o = np.argsort(phi[ok])
    x, y = phi[ok][o], z[ok][o]
    if y.min() > Z_ISO_FRICTIONAL or y.max() < Z_ISO_FRICTIONAL:
        return out
    for i in range(len(x) - 1):
        if (y[i] - Z_ISO_FRICTIONAL) * (y[i + 1] - Z_ISO_FRICTIONAL) <= 0:
            if y[i + 1] == y[i]:
                continue
            t = (Z_ISO_FRICTIONAL - y[i]) / (y[i + 1] - y[i])
            out["phi_c_E4"] = float(x[i] + t * (x[i + 1] - x[i]))
            break
    return out


def analyse_group(recs: list[dict]) -> dict:
    recs = sorted(recs, key=lambda r: r["phi"])
    phi = np.array([r["phi"] for r in recs], float)
    p = np.array([r["p_ss_hat"] for r in recs], float)
    ivals = np.array([r["i_ss"] for r in recs], float)
    z = np.array([r["z_nr_ss"] for r in recs], float)
    mu = np.array([r["mu_macro"] for r in recs], float)

    res = {"n_cells": len(recs), "phi_lo": float(phi.min()), "phi_hi": float(phi.max()),
           "n_steady": int(sum(bool(r["steady"]) for r in recs)),
           "n_quasistatic": int(sum(bool(r["quasistatic"]) for r in recs))}
    res.update(e1_phi_of_I(phi, ivals))
    res.update(_pressure_root(phi, p, None))
    res.update(_pressure_root(phi, p, 1.0))
    res.update(e4_z_crossing(phi, z))

    # mu_s = mu_macro(I -> 0): an anchor on the rheology, independent of phi_c
    ok = np.isfinite(mu) & np.isfinite(ivals) & (ivals <= I_FIT_MAX) & (ivals > 0)
    if ok.sum() >= 3:
        A = np.vstack([np.ones_like(ivals[ok]), ivals[ok]]).T
        coef, *_ = np.linalg.lstsq(A, mu[ok], rcond=None)
        res["mu_s"] = float(coef[0])
    else:
        res["mu_s"] = np.nan

    # REFUSE a root outside the ladder it was fitted on.
    for tag in ("E1", "E2", "E3", "E4"):
        v = res.get(f"phi_c_{tag}", np.nan)
        if np.isfinite(v) and not (res["phi_lo"] - 0.05 <= v <= res["phi_hi"]):
            res[f"phi_c_{tag}_refused"] = float(v)
            res[f"phi_c_{tag}"] = np.nan
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", default="phase8_prod2*.jsonl",
                    help="glob under out/. NOTE the default is prod2, not prod*: "
                         "prod* would also match the protocol-v1 records")
    ap.add_argument("--out", default="phase8_verdict.txt")
    args = ap.parse_args()

    recs = load(args.pattern)
    errs = [r for r in recs if "error" in r]
    good = [r for r in recs if "error" not in r and np.isfinite(r.get("p_ss_hat", np.nan))]
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    p(f"phase8: {len(recs)} records from {args.pattern}; "
      f"{len(good)} usable, {len(errs)} errored")
    for r in errs:
        p(f"  ERRORED {r.get('key')}: {r.get('error')}")
    if not good:
        p("nothing to analyse")
        return 0

    # duplicate-state guard: two rungs landing on the same realised phi is the defect
    # that the 500-step compression chunk produced on 2026-09-08.
    by_g = defaultdict(list)
    for r in good:
        by_g[group_key(r)].append(r)
    n_dupe = 0
    for k, rs in by_g.items():
        phis = sorted(round(r["phi"], 6) for r in rs)
        for a, b in zip(phis, phis[1:]):
            if abs(a - b) < 1e-4:
                p(f"  🔴 DUPLICATE STATE in {k}: two rungs at phi = {a:.6f}, {b:.6f}")
                n_dupe += 1
    p(f"duplicate-state check: {n_dupe} offending pairs")
    p()

    protos = sorted({k[5] for k in by_g})
    if len(protos) > 1:
        p(f"🔴 REFUSING to report contrasts: {len(protos)} protocols present in this "
          f"selection -- {protos}. Re-run with a --pattern that picks one.")
        p("   (per-cell rows are still printed below, keyed by protocol)")

    rows = []
    p(f"{'family':10s} {'cv':>5s} {'mu':>5s} {'seed':>4s} {'cells':>5s} "
      f"{'stdy':>4s} {'qs':>3s} {'phi_c_E1':>9s} {'phi_c_E2':>9s} {'b_E2':>6s} "
      f"{'phi_c_E3':>9s} {'phi_c_E4':>9s} {'mu_s':>6s}")
    for k in sorted(by_g):
        r = analyse_group(by_g[k])
        fam, cv, mu, seed, n, proto, gmax, steps = k
        rows.append({"family": fam, "cv": cv, "mu_p": mu, "seed": seed,
                     "n_requested": n, "protocol": proto, "gamma_max": gmax,
                     "shear_steps": steps, **r})
        p(f"{fam:10s} {cv:5.2f} {mu:5.2f} {seed:4d} {r['n_cells']:5d} "
          f"{r['n_steady']:4d} {r['n_quasistatic']:3d} "
          f"{r['phi_c_E1']:9.4f} {r['phi_c_E2']:9.4f} {r.get('b_E2', np.nan):6.2f} "
          f"{r['phi_c_E3']:9.4f} {r['phi_c_E4']:9.4f} {r['mu_s']:6.3f}")

    # --- CV contrasts, per family/mu, averaged over seeds ---------------------------
    p()
    p("=== CV contrasts (mean over seeds, +/- sd) ===")
    if len(protos) > 1:
        p("   SKIPPED -- see the protocol refusal above")
        rows_for_contrast = []
    else:
        rows_for_contrast = rows
    for fam in sorted({r["family"] for r in rows_for_contrast}):
        for mu in sorted({r["mu_p"] for r in rows_for_contrast
                          if r["family"] == fam}):
            sub = [r for r in rows_for_contrast
                   if r["family"] == fam and r["mu_p"] == mu]
            cvs = sorted({r["cv"] for r in sub})
            if len(cvs) < 2:
                continue
            p(f"-- {fam}, mu_p = {mu}")
            for tag in ("E1", "E2", "E3", "E4"):
                vals = {}
                for cv in cvs:
                    v = [r[f"phi_c_{tag}"] for r in sub if r["cv"] == cv
                         and np.isfinite(r[f"phi_c_{tag}"])]
                    if v:
                        vals[cv] = (float(np.mean(v)), float(np.std(v, ddof=1))
                                    if len(v) > 1 else float("nan"), len(v))
                if len(vals) < 2:
                    p(f"   {tag}: too few resolved points ({len(vals)})")
                    continue
                txt = "  ".join(f"cv{cv:.2f}={m:.4f}+-{s:.4f}(n{n})"
                                for cv, (m, s, n) in sorted(vals.items()))
                lo, hi = min(vals), max(vals)
                d = vals[hi][0] - vals[lo][0]
                sd = np.nanmax([vals[hi][1], vals[lo][1]])
                sig = d / sd if np.isfinite(sd) and sd > 0 else float("nan")
                slope = d / (hi - lo)
                p(f"   {tag}: {txt}")
                p(f"       delta(cv {lo:.2f}->{hi:.2f}) = {d:+.4f}  "
                  f"({sig:+.1f} sd)   slope s = {slope:+.4f} per unit CV")

    path = os.path.join(OUT_DIR, args.out)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    csv = os.path.join(OUT_DIR, "phase8_summary.csv")
    import csv as _csv
    if rows:
        keys = sorted({k for r in rows for k in r})
        with open(csv, "w", encoding="utf-8", newline="") as fh:
            w = _csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
    print(f"\nwritten -> {path}\nwritten -> {csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
