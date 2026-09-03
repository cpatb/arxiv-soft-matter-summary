"""diag_radius_estimator.py — is the confinement's 10x CV work gap the packing, or `rmax`?

THE QUESTION (REPORT §7g)
-------------------------
After the net-force fix, the confining force still does 10x more work on a CV = 0.5
aggregate than on a CV = 0 one (`W_conf/ke0` 2.906 against 0.289). Two readings, and the
sign of the effect does not separate them:

  (P) PACKING.   The CV = 0.5 body launches at Z = 2.33 against 4.17, so it is genuinely
                 looser, deforms more under the same confining pressure, and W = int F.v dt
                 is larger because the grains simply move more. The work gap is then a
                 CONSEQUENCE of the packing result, not a contaminant of it.

  (E) ESTIMATOR. The force is `f_tot = p_cap * 4*pi*rmax^2` with `rmax = max|r_i - r_cm|`,
                 a SINGLE EXTREMAL grain. A loose polydisperse shell with one outlying grain
                 inflates rmax, hence f_tot, hence the work -- and inflates it more at high
                 CV, because a broad size distribution puts more grains far from the centre.

WHY THE GEOMETRIC HALF IS ENOUGH TO TEST (E)
--------------------------------------------
(E) is a claim about the CARVED GEOMETRY ALONE. It needs no impact, no wall and no
dynamics: if `rmax` is inflated relative to a robust radius, that is visible in the static
aggregate the instant it is carved. The force scales as the SQUARE of whatever radius is
used, so the quantity that matters is

    inflation = (rmax / r_robust)^2

and the discriminating number is how that inflation DIFFERS between CV = 0 and CV = 0.5.
This runs in one carve per cell instead of a carve plus a 100-probe impact loop.

PRE-REGISTERED READING, fixed before the numbers were seen
----------------------------------------------------------
Let  G = inflation(CV = 0.5) / inflation(CV = 0),  the factor by which the estimator
favours the polydisperse cell. The observed work gap is 10.06x.

  G >= 3.0     (E) carries most of the gap. The estimator must be replaced before any
               dynamic CV contrast is quoted, and Phase 7 should be re-run.
  G <= 1.5     (E) is dead. The estimator is near-identical at both CVs, so the work gap
               is the packing -- reading (P) -- and it is a consequence of the physics
               rather than an artefact. No code change is warranted.
  otherwise    partial; report the split and run the dynamics arm.

Note the asymmetry that makes this test worth running FIRST: G <= 1.5 settles the question
with no dynamics at all, while G >= 3.0 only tells us to go and measure.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

CVS = (0.0, 0.5)
SEEDS = (1, 2, 3, 4)
MU = 0.3
CARVE_P_HAT = 1e-3
N_PARENT = 800
SHELL_FRACTION = 0.85

# The robust alternatives. p95 is the headline one: high enough to still describe the
# outer surface, low enough that a single stray grain cannot set it.
PERCENTILES = (100.0, 99.0, 95.0, 90.0)


def radius_stats(cv, seed):
    """Carve one aggregate and measure how much its radius estimator depends on outliers."""
    from phase4_impact import carve_aggregate

    t0 = time.time()
    pos, _rad, phi_parent, mom = carve_aggregate(cv, MU, N_PARENT, seed, 1.0, CARVE_P_HAT)
    if pos.shape[0] < 50:
        return {"failed": True, "cv": cv, "seed": seed, "n": int(pos.shape[0])}

    rel = pos - pos.mean(axis=0)
    rr = np.linalg.norm(rel, axis=1)
    out = {"cv": cv, "seed": seed, "n_grains": int(pos.shape[0]),
           "phi_parent": float(phi_parent),
           "cv_realised": float(mom.get("cv", float("nan"))),
           "wall_s": time.time() - t0}
    for q in PERCENTILES:
        out["r_p%g" % q] = float(np.percentile(rr, q))
    rmax = out["r_p100"]
    # the mean radius of the grains the confinement actually acts on
    shell = rr > SHELL_FRACTION * rmax
    out["r_shell_mean"] = float(rr[shell].mean()) if shell.any() else float("nan")
    out["n_shell"] = int(shell.sum())
    # force inflation: the confining force goes as the SQUARE of the radius used
    for q in (99.0, 95.0, 90.0):
        out["inflation_p%g" % q] = float((rmax / out["r_p%g" % q]) ** 2)
    out["inflation_shellmean"] = float((rmax / out["r_shell_mean"]) ** 2)
    print("  cv=%.2f s=%d n=%d rmax=%.3f p95=%.3f infl_p95=%.4f (%.0fs)"
          % (cv, seed, out["n_grains"], rmax, out["r_p95"],
             out["inflation_p95"], out["wall_s"]), flush=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="")
    args = ap.parse_args()
    seeds = tuple(int(x) for x in args.seeds.split(",")) if args.seeds else SEEDS

    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, "diag_radius_estimator.jsonl")
    recs = []
    for cv in CVS:
        for seed in seeds:
            try:
                r = radius_stats(cv, seed)
            except Exception as exc:               # noqa: BLE001
                print("FAILED cv=%s seed=%s: %r" % (cv, seed, exc), flush=True)
                continue
            recs.append(r)
            with open(path, "a", encoding="utf-8") as fh:
                fh.write(json.dumps(r, sort_keys=True) + "\n")

    good = [r for r in recs if not r.get("failed")]
    if not good:
        print("no cells")
        return 1

    def mean_of(g, k):
        return float(np.mean([r[k] for r in g]))

    lines = ["=" * 78,
             "IS THE CONFINEMENT'S CV WORK GAP AN rmax ARTEFACT?",
             "=" * 78, "",
             "Confining force goes as radius^2, so the quantity that matters is the SQUARE",
             "of the ratio between the extremal radius and a robust one.", ""]
    lines.append("  %5s %3s %8s %8s %8s %9s %9s %9s"
                 % ("CV", "n", "rmax", "r_p95", "r_shell",
                    "infl_p99", "infl_p95", "infl_p90"))
    per_cv = {}
    for cv in CVS:
        g = [r for r in good if r["cv"] == cv]
        if not g:
            continue
        per_cv[cv] = {k: mean_of(g, k) for k in
                      ("r_p100", "r_p95", "r_shell_mean",
                       "inflation_p99", "inflation_p95", "inflation_p90")}
        lines.append("  %5.2f %3d %8.3f %8.3f %8.3f %9.4f %9.4f %9.4f"
                     % (cv, len(g), per_cv[cv]["r_p100"], per_cv[cv]["r_p95"],
                        per_cv[cv]["r_shell_mean"], per_cv[cv]["inflation_p99"],
                        per_cv[cv]["inflation_p95"], per_cv[cv]["inflation_p90"]))
    lines.append("")

    tag = "INCONCLUSIVE"
    if 0.0 in per_cv and 0.5 in per_cv:
        lines.append("PRE-REGISTERED DISCRIMINATOR (set before the numbers were seen)")
        lines.append("  G = inflation(CV=0.5) / inflation(CV=0); observed work gap is 10.06x.")
        lines.append("  G >= 3.0 -> estimator carries the gap, replace it and re-run Phase 7")
        lines.append("  G <= 1.5 -> estimator is dead; the gap is the PACKING, no fix needed")
        lines.append("")
        verdicts = []
        for q in ("inflation_p99", "inflation_p95", "inflation_p90"):
            gg = per_cv[0.5][q] / per_cv[0.0][q]
            v = ("ESTIMATOR" if gg >= 3.0 else "PACKING" if gg <= 1.5 else "partial")
            verdicts.append(v)
            lines.append("    %-18s CV0 %7.4f  CV0.5 %7.4f   G = %6.3f  -> %s"
                         % (q, per_cv[0.0][q], per_cv[0.5][q], gg, v))
        lines.append("")
        if all(v == "PACKING" for v in verdicts):
            tag = "PACKING"
            lines.append("  => VERDICT: PACKING. The extremal radius is not meaningfully more")
            lines.append("     inflated at CV = 0.5 than at CV = 0, so it cannot manufacture a")
            lines.append("     10x work gap. Reading (P) stands: the polydisperse aggregate is")
            lines.append("     looser, so the same confining pressure does more work on it.")
            lines.append("     That is a consequence of the packing result, not an artefact of")
            lines.append("     the boundary condition, and no code change is warranted.")
            lines.append("     The dynamic CV contrast stays ORDINAL for the separate reason")
            lines.append("     given in 7g: the boundary is energetically dominant either way.")
        elif any(v == "ESTIMATOR" for v in verdicts):
            tag = "ESTIMATOR"
            lines.append("  => VERDICT: ESTIMATOR. Replace rmax with the p95 radius in")
            lines.append("     phase4/phase7 confine(), and re-run Phase 7 before quoting any")
            lines.append("     dynamic CV contrast.")
        else:
            tag = "PARTIAL"
            lines.append("  => PARTIAL. The geometric screen does not settle it; the dynamics")
            lines.append("     arm is required.")
        lines.append("")
        lines.append("  RADIUS_ESTIMATOR_" + tag)
    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(OUT_DIR, "radius_estimator_verdict.txt"), "w",
              encoding="utf-8") as fh:
        fh.write(txt + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
