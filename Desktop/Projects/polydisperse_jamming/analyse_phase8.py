"""analyse_phase8.py — extract phi_c(CV) from the constant-volume shear ladders.

FOUR ESTIMATORS, AND WHY THE PRIMARY IS THE ONE IT IS
----------------------------------------------------
At constant volume every cell reaches a steady state that lies on the granular
`mu(I)` / `phi(I)` rheology, so the ladder is not a set of independent points: it is a
sampling of one curve, parameterised by the inertial number `I = gammadot d50
sqrt(rho/p)`. The critical-state fraction is the quasi-static limit of that curve.

  E2  PRESSURE ROOT, free exponent  [PRIMARY]   p_ss = c (phi - phi_c)^b, fitted in log
      space. The critical state is where the constant-volume grain pressure VANISHES, so
      extrapolating p -> 0 targets it directly. `b` is reported as a diagnostic of the
      contact law, never used to redefine phi_c.
  E3  PRESSURE ROOT, b pinned to 1      the harmonic-contact convention, as a check.
  E4  COORDINATION CROSSING             the phi at which z_nr crosses Z_ISO_FRICTIONAL.
      Independent of the stress entirely, which is its whole value: a defect in the
      stress read-out cannot move it.

  E1  phi(I) EXTRAPOLATION   🔴 **RETRACTED 2026-09-08 AS A phi_c ESTIMATOR**, having
      been the PRIMARY one. It fits phi against the inertial number and extrapolates to
      I = 0, which is the standard construction WHEN the pressure is imposed and the
      shear rate is varied. This protocol does the opposite: `gdot` is IDENTICAL on every
      rung of a ladder (measured: 3.000e-05 across all five at cv = 0, varying by
      1.0000x), so I = gdot d sqrt(rho/p) moves ONLY through p. Therefore

          I -> 0  is the HIGH-PRESSURE, DENSEST end of the ladder,
                  NOT the critical state.

      The extrapolation ran away from phi_c, not toward it. That is why E1 sat 0.074
      ABOVE E2 -- it was reporting the dense end while E2 reported the vanishing-pressure
      end, and the "systematic between estimators" was two different quantities. E1 is
      still computed and printed, relabelled `phi_dense_E1`, because the dense-end value
      is a real number and its CV contrast happens to agree with E2's; it is no longer
      offered as phi_c.

      ⚠ The CV CONTRAST is barely affected (+0.0578 from E1 vs +0.0623 from E2), because
      both ends of the ladder shift together with CV. P1 and P2 therefore survive the
      retraction; what does not survive is any ABSOLUTE phi_c quoted from E1.

The four are reported side by side for every cell group. A CV contrast is quoted only
if it has the same sign under all four; the SPREAD between them is the systematic and
it is printed.

WHAT IS REFUSED
---------------
  * a root outside the ladder it was fitted on (extrapolation dressed as measurement)
  * 🔴 **NOT WHAT THIS LINE USED TO SAY.** It claimed "a group whose cells are not all
    steady, or not all quasi-static -- excluded cells are PRINTED with their reason,
    never dropped silently". That was FALSE for the primary estimator: `e1_phi_of_I`
    does filter on `I <= I_FIT_MAX`, but `_pressure_root` (E2 and E3) and
    `e4_z_crossing` filter on NOTHING except a pressure floor. The flags were COUNTED
    into `n_steady`/`n_quasistatic` and printed, and never used to exclude.
    Found by a numerics reviewer 2026-09-09. MEASURED IMPACT on the headline, prod2:

        as shipped (no gate applied)    dphi_c = +0.0623
        quasi-static cells only         dphi_c = +0.0444
        steady cells only               dphi_c = +0.0485
        both gates applied              dphi_c = +0.0438   (but n = 1 at cv = 0)

    i.e. **the headline falls ~30 %**. The sign and the qualitative conclusion survive;
    the number does not. The exclusions are also CV-ASYMMETRIC (3/3/1 non-quasi-static
    at cv = 0/0.2/0.5), which is the direction that inflates the contrast.
    ⚠ Applying both gates leaves ONE seed at cv = 0, so the fix is NOT simply to switch
    the gates on -- the run needs more cells at the low-CV end. See NEXT_STEPS_2026-09-09.md
    item 1. Until then, quote +0.044 to +0.062 and say why.
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
# 🔴 The isostatic contact number is NOT one constant: it is 2d = 6 for FRICTIONLESS
# spheres and d + 1 = 4 for fully frictional ones. Hardcoding 4 was right for the mu_p =
# 0.3 production arm and WRONG for the mu_p = 0 anchor run, where the pack is isostatic
# at 6 and E4 would have hunted for a crossing at 4 that happens far below jamming --
# a constant that is correct for one arm and silently wrong for another, which is the
# trap this project keeps meeting. Chosen per cell from mu_p.
Z_ISO_FRICTIONLESS = 6.0
Z_ISO_FRICTIONAL = 4.0  # the project's A2/A3 anchors measure 6.37-6.92 at mu_p = 0 and
                        # 3.98-4.07 at mu_p = 3, which is what pins both values.
Z_ISO_MU_BAR = 1.0e-9   # below this mu_p the pack is treated as frictionless
P_FLOOR = 1.0e-9        # below this p_hat the cell carries no measurable stress
B_MAX = 8.0             # upper bound on E2's free exponent. A fit landing ON it is
                        # reported as bound-limited and its phi_c is refused.
PHI_MATCH_TOL = 5.0e-3  # how close two rungs from different CV ladders must be to count
                        # as "matched phi". The ladder spacing is 0.020, so this admits
                        # only genuine near-coincidences and never two different rungs.
PHI_GROUP_TOL = 1.0e-3  # how close two SEEDS' realised phi must be to be one rung;
                        # stage B lands within ~1e-4, so this is loose by 10x


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
    """phi at the DENSE end, lim_{I->0} phi(I), over the quasi-static rungs.

    🔴 NOT phi_c. See the module docstring: on this protocol gdot is fixed, so I moves
    only through p and I -> 0 is the densest rung. Kept and reported because the value
    is real and its CV contrast is informative; renamed so it cannot be read as phi_c.
    """
    ok = np.isfinite(phi) & np.isfinite(ivals) & (ivals <= I_FIT_MAX) & (ivals > 0)
    out = {"phi_dense_E1": np.nan, "slope_E1": np.nan, "n_E1": int(ok.sum()),
           "r2_E1": np.nan, "I_min": np.nan, "I_max": np.nan}
    if ok.sum() < 3:
        return out
    x, y = ivals[ok], phi[ok]
    A = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    out.update(phi_dense_E1=float(coef[0]), slope_E1=float(coef[1]),
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
        bounds = ([-np.inf, -np.inf, 0.1], [lo - 1e-6, np.inf, B_MAX])
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

    # 🔴 A FIT ON A BOUND IS NOT A MEASUREMENT. Both roots are bounded above by
    # `lo - 1e-6` (the root must lie below the lowest rung). On the 2026-09-08 prod2 run
    # E3 returned EXACTLY that bound on 6 of 9 groups -- 0.5745 at cv 0 and 0.5880 at
    # cv 0.2, i.e. the ladder's own bottom rung, with a seed sd of 0.0000 and an
    # apparent significance of 161 sd. That is the bound reporting itself, not phi_c.
    # Flagged and returned as NaN, with the value kept under a `_at_bound` key.
    if abs(popt[0] - (lo - 1e-6)) < 1e-4:
        out[f"phi_c_{tag}_at_bound"] = float(popt[0])
        out[f"phi_c_{tag}"] = np.nan

    # 🔴 AND THE SAME RULE FOR THE EXPONENT. E2's `b` is bounded to [0.1, 8.0]; a fit
    # that lands ON that bound has not found a power law, it has run out of room, and
    # its phi_c is then whatever the bound forced. Measured on the 2026-09-08 P4
    # frictionless anchor: b = 8.00, 7.33, 8.00 on the three seeds, and the phi_c it
    # returned (0.580) missed the anchor by 0.052 while the isostatic crossing hit it to
    # 0.0015. Frictionless jamming is genuinely much sharper than a low-order power law,
    # so this is the fit reporting its own inadequacy -- which is worth keeping, not
    # hiding.
    if pin_b is None and popt[2] >= B_MAX - 1e-6:
        out["b_E2_at_bound"] = float(popt[2])
        out["phi_c_E2_at_bound"] = float(popt[0])
        out["phi_c_E2"] = np.nan
    return out


def e4_z_crossing(phi: np.ndarray, z: np.ndarray, mu_p: float = 1.0) -> dict:
    """phi where z_nr crosses the isostatic count, by linear interpolation.

    The target depends on friction: 6 for frictionless spheres, 4 for frictional ones.
    Passing `mu_p` is mandatory in spirit -- the default is the frictional case only so
    that older callers keep their previous behaviour rather than silently changing it.
    """
    z_iso = (Z_ISO_FRICTIONLESS if mu_p <= Z_ISO_MU_BAR else Z_ISO_FRICTIONAL)
    out = {"phi_c_E4": np.nan, "z_iso_used": z_iso}
    ok = np.isfinite(phi) & np.isfinite(z)
    if ok.sum() < 2:
        return out
    o = np.argsort(phi[ok])
    x, y = phi[ok][o], z[ok][o]
    if y.min() > z_iso or y.max() < z_iso:
        return out
    for i in range(len(x) - 1):
        if (y[i] - z_iso) * (y[i + 1] - z_iso) <= 0:
            if y[i + 1] == y[i]:
                continue
            t = (z_iso - y[i]) / (y[i + 1] - y[i])
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
    res["p_hat_max"] = float(np.nanmax(p))
    res.update(e1_phi_of_I(phi, ivals))
    res.update(_pressure_root(phi, p, None))
    res.update(_pressure_root(phi, p, 1.0))
    res.update(e4_z_crossing(phi, z, mu_p=float(recs[0].get('mu_p', 1.0))))

    # mu_s = mu_macro(I -> 0): an anchor on the rheology, independent of phi_c
    ok = np.isfinite(mu) & np.isfinite(ivals) & (ivals <= I_FIT_MAX) & (ivals > 0)
    if ok.sum() >= 3:
        A = np.vstack([np.ones_like(ivals[ok]), ivals[ok]]).T
        coef, *_ = np.linalg.lstsq(A, mu[ok], rcond=None)
        res["mu_s"] = float(coef[0])
    else:
        res["mu_s"] = np.nan

    # REFUSE a root outside the ladder it was fitted on.
    for tag in ("E2", "E3", "E4"):
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
      f"{'stdy':>4s} {'qs':>3s} {'phi_c_E2*':>9s} {'b_E2':>6s} "
      f"{'phi_c_E3':>9s} {'phi_c_E4':>9s} {'phiDenseE1':>10s} {'mu_s':>6s} "
      f"{'p_hat_max':>10s}")
    p("   (* E2 is PRIMARY; phiDenseE1 is the retracted E1, relabelled -- see docstring)")
    for k in sorted(by_g):
        r = analyse_group(by_g[k])
        fam, cv, mu, seed, n, proto, gmax, steps = k
        rows.append({"family": fam, "cv": cv, "mu_p": mu, "seed": seed,
                     "n_requested": n, "protocol": proto, "gamma_max": gmax,
                     "shear_steps": steps, **r})
        p(f"{fam:10s} {cv:5.2f} {mu:5.2f} {seed:4d} {r['n_cells']:5d} "
          f"{r['n_steady']:4d} {r['n_quasistatic']:3d} "
          f"{r['phi_c_E2']:9.4f} {r.get('b_E2', np.nan):6.2f} "
          f"{r['phi_c_E3']:9.4f} {r['phi_c_E4']:9.4f} "
          f"{r['phi_dense_E1']:10.4f} {r['mu_s']:6.3f} {r['p_hat_max']:10.3e}")

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
            for tag in ("E2", "E3", "E4", "dense_E1"):
                vals = {}
                for cv in cvs:
                    col = ("phi_dense_E1" if tag == "dense_E1"
                           else f"phi_c_{tag}")
                    v = [r[col] for r in sub if r["cv"] == cv
                         and np.isfinite(r[col])]
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

    # --- P3: the CV contrast at MATCHED phi -----------------------------------------
    # The ladders are centred per CV, so they share no exact phi. Where they OVERLAP the
    # nearest rungs are still a direct test of the prediction that raising CV at fixed
    # volume fraction collapses the undrained grain pressure -- and it needs no fit at
    # all, so no estimator systematic enters.
    #
    # 🔴 EVERY overlapping pair is reported, not the closest one. The first version took
    # the single nearest pair; at this ladder spacing TWO pairs tie to within rounding
    # (both 0.0023 apart in phi) and they give 43x and 5.8x, so which one a min() returns
    # is a floating-point accident. Reporting one of them would have been cherry-picking
    # by tie-break.
    p()
    p("=== P3: the CV contrast at MATCHED phi (no fit, no estimator) ===")
    if len(protos) > 1:
        p("   SKIPPED -- more than one protocol present")
    else:
        good2 = [r for r in good if np.isfinite(r.get("phi", np.nan))]
        cvs = sorted({round(float(r["cv"]), 3) for r in good2})
        pairs = []
        for i_, c1 in enumerate(cvs):
            for c2 in cvs[i_ + 1:]:
                A_ = [r for r in good2 if round(float(r["cv"]), 3) == c1]
                B_ = [r for r in good2 if round(float(r["cv"]), 3) == c2]
                seen = set()
                for ra in A_:
                    for rb in B_:
                        d = abs(ra["phi"] - rb["phi"])
                        if d > PHI_MATCH_TOL:
                            continue
                        key = (c1, c2, round(ra["phi"], 4), round(rb["phi"], 4))
                        if key in seen:
                            continue
                        seen.add(key)
                        pairs.append((c1, c2, ra["phi"], rb["phi"], d))
        if not pairs:
            p(f"   no cross-CV phi pair within {PHI_MATCH_TOL}")
        else:
            p(f"   {len(pairs)} overlapping pairs within {PHI_MATCH_TOL} in phi; "
              f"ALL are reported")
            ratios = []
            for c1, c2, p1, p2, d in sorted(pairs, key=lambda t: t[2]):
                A_ = [r for r in good2 if round(float(r["cv"]), 3) == c1
                      and abs(r["phi"] - p1) < PHI_GROUP_TOL]
                B_ = [r for r in good2 if round(float(r["cv"]), 3) == c2
                      and abs(r["phi"] - p2) < PHI_GROUP_TOL]
                pa = np.array([r["p_ss_hat"] for r in A_], float)
                pb = np.array([r["p_ss_hat"] for r in B_], float)
                za = np.array([r["z_nr_ss"] for r in A_], float)
                zb = np.array([r["z_nr_ss"] for r in B_], float)
                ra_ = np.array([r["rattler_frac_ss"] for r in A_], float)
                rb_ = np.array([r["rattler_frac_ss"] for r in B_], float)
                sd = max(pa.std(ddof=1) if len(pa) > 1 else 0.0,
                         pb.std(ddof=1) if len(pb) > 1 else 0.0)
                ratio = pa.mean() / pb.mean()
                ratios.append(ratio)
                nsd = abs(pa.mean() - pb.mean()) / sd if sd > 0 else float("inf")
                p(f"   phi {p1:.5f}(cv{c1:.2f}, n{len(pa)}) vs {p2:.5f}"
                  f"(cv{c2:.2f}, n{len(pb)}), mismatch {d:.5f}")
                p(f"      p_ss {pa.mean():.4e} -> {pb.mean():.4e}  "
                  f"ratio {ratio:6.1f}x  ({nsd:.0f} sd)   "
                  f"Z_nr {za.mean():.2f} -> {zb.mean():.2f}   "
                  f"rattlers {ra_.mean():.3f} -> {rb_.mean():.3f}")
                if p2 > p1:
                    p("      (the higher-CV cell is the DENSER of the two, which works "
                      "AGAINST the effect: conservative)")
                # ROBUSTNESS: analyse_rheology flags the cv = 0.5 cells as not having
                # settled to +/-10 % within gamma_max. Recompute the ratio from a LATE
                # window (the last 15 % of the strain) and print both, rather than
                # asserting that the flag does not matter.
                def _late(rs):
                    v = []
                    for rr in rs:
                        g = np.asarray(rr.get("gamma", []), float)
                        pp = np.asarray(rr.get("p_hat_series", []), float)
                        if g.size and g.size == pp.size:
                            v.append(float(np.mean(pp[g >= 0.85 * g.max()])))
                    return float(np.mean(v)) if v else float("nan")
                la, lb = _late(A_), _late(B_)
                if np.isfinite(la) and np.isfinite(lb) and lb != 0:
                    p(f"      robustness: on a LATE window (last 15 % of strain) the "
                      f"ratio is {la / lb:.1f}x against {ratio:.1f}x -- "
                      f"a {abs(la / lb / ratio - 1) * 100:.0f} % change")
            if len(ratios) > 1:
                p(f"   spread across pairs: {min(ratios):.1f}x to {max(ratios):.1f}x "
                  f"-- a factor {max(ratios) / min(ratios):.1f} between pairs that "
                  f"differ only in which rung they sit on. Quote the RANGE.")

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
