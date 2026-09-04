"""analyse_exponent_constrained.py — is the failing exponent anchor a fit degeneracy?

THE STANDING FAILURE
--------------------
REPORT anchor A6: a harmonic contact law predicts `phi = phi_J + c*p^b` with b ~ 1.
Yade's Law2_ScGeom_FrictPhys_CundallStrack IS linear (Fn = kn*delta), so b = 1 is the
right prediction. The free two-parameter fit returns median b = 0.567 over 22 informative
ladders, against a pre-registered bar of [0.6, 1.6]. That has been reported as a FAIL and
the bar has never been moved.

WHY A DEGENERACY IS THE OBVIOUS SUSPECT
---------------------------------------
`phi_J` and `b` are fitted TOGETHER, and they trade off: lowering phi_J raises (phi - phi_J)
at every p, which a smaller b then flattens. The signature of that trade is a correlation
between the fitted b and how far the fitted phi_J has drifted from an INDEPENDENT estimate
of the same quantity. Measured on the existing fits (2026-09-03):

    corr(phi_J_free - phi_J_lin, b_free) = +0.675, p = 0.00056, n = 22
    ladders with b below median:  phi_J_free - phi_J_lin = -0.0365   (20 seed scatters)
    ladders with b above median:  phi_J_free - phi_J_lin = -0.0023   (~1 seed scatter)

So the low-b ladders are exactly the ones whose free phi_J has run far below the linear
extrapolation. That is a degeneracy, not a measurement of the contact law.

WHAT THIS SCRIPT DOES
---------------------
Refits each ladder with `phi_J` FIXED at `phi_J_lin`, the linear-extrapolation estimate,
which is computed WITHOUT the exponent fit and is therefore independent of b. With phi_J
pinned, b comes from a straight line in log-log:

    log(phi - phi_J) = log(c) + b * log(p_hat)

This is a STRICTER test than the one that failed, not a looser one: it removes a free
parameter. A bar is not being moved -- the same [0.6, 1.6] bar is applied to a
better-posed estimator, and both numbers are reported side by side.

PRE-REGISTERED READING (fixed before the constrained fit was run)
-----------------------------------------------------------------
  constrained median b in [0.6, 1.6]  -> the free fit's low b is a DEGENERACY artefact.
                                         NOTE, added after running: this must NOT be read as
                                         confirming b = 1, because phi_J_lin is a linear-fit
                                         intercept and so already assumes b = 1. The correct
                                         conclusion is that the design cannot determine b;
                                         see the conditioning block.
  constrained median b still < 0.6    -> the failure is PHYSICAL, not a fitting artefact.
                                         A6 stays FAILED and the discrepancy with the
                                         harmonic prediction is a real open problem.
  constrained fit unstable / < 3 usable points on most ladders
                                      -> report inconclusive; do not claim either way.

A ladder is admitted on the SAME gate as the free fit (monotone in p, range > 2 seed
scatters). Note that MORE ladders clear that gate than the free fit produced a b for --
the free fit also fails to converge on some of them -- so the headline comparison is
PAIRED, restricted to ladders where BOTH estimators returned a value. The larger
constrained-only set is reported beside it, clearly labelled, because a comparison of two
estimators over two different samples is not a comparison.
"""

from __future__ import annotations

import csv
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

BAR = (0.6, 1.6)
MIN_POINTS = 3
SEED_SCATTER = 0.00179          # median over 106 replicated cells, from analyse_phase2
EXPONENT_MIN_RANGE_IN_SD = 2.0

# Must match analyse_phase2.P_HAT_FIT_MAX. The first version of this file omitted the trim
# and the dedup below, so it refitted a DIFFERENT point set from the one `phi_J_lin` was
# built on -- which is exactly the inconsistency the whole analysis is about. Found by the
# code reviewer 2026-09-03; it moved the conditioning number from 148 to 171-194.
P_HAT_FIT_MAX = 1.1e-3


def family_of(rec):
    k = rec.get("key", "")
    for f in ("uniform", "bidisperse"):
        if f in k:
            return f
    return "lognormal"


def load_cells():
    # DEDUP BY CELL KEY, as analyse_phase2.load_all does. Killed multi-shard launches left
    # the same cell on disk more than once (463 records, 424 unique keys), and np.polyfit
    # would silently double-weight the duplicates.
    by_key = {}
    for fp in glob.glob(os.path.join(OUT, "phase2*.jsonl")):
        with open(fp, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if r.get("failed") or not r.get("converged"):
                    continue
                if r.get("key"):
                    by_key[r["key"]] = r
    return list(by_key.values())


def main() -> int:
    cells = load_cells()
    if not cells:
        print("no phase2 cells")
        return 1

    fits_fp = os.path.join(OUT, "phase2_fits.csv")
    if not os.path.exists(fits_fp):
        print("no phase2_fits.csv -- run analyse_phase2.py first")
        return 1
    fits = list(csv.DictReader(open(fits_fp, encoding="utf-8")))

    def fnum(d, k):
        try:
            return float(d[k])
        except (TypeError, ValueError, KeyError):
            return float("nan")

    # group the raw cells into ladders keyed the same way the fits table is
    ladders: dict[tuple, list] = {}
    for r in cells:
        key = (round(float(r["cv"]), 3), round(float(r["mu_p"]), 3),
               int(r["seed"]), family_of(r), int(r["n_requested"]))
        ladders.setdefault(key, []).append(r)

    rows = []
    for fr in fits:
        # N is not a column in the fits table; it was built from N = 1500 only.
        key = (round(fnum(fr, "cv_nom"), 3), round(fnum(fr, "mu_p"), 3),
               int(float(fr["seed"])), fr["family"], 1500)
        lad = ladders.get(key)
        if not lad:
            continue
        monotone = str(fr.get("phi_monotone", "")).lower() in ("true", "1")
        rng = fnum(fr, "phi_range")
        informative = monotone and rng > EXPONENT_MIN_RANGE_IN_SD * SEED_SCATTER
        phi_J_lin = fnum(fr, "phi_J_lin")
        b_free = fnum(fr, "b_free")
        if not (informative and np.isfinite(phi_J_lin)):
            continue

        lad = [c for c in lad if float(c["p_hat"]) <= P_HAT_FIT_MAX]
        if len(lad) < MIN_POINTS:
            continue
        p = np.array([float(c["p_hat"]) for c in lad])
        phi = np.array([float(c["phi"]) for c in lad])
        o = np.argsort(p)
        p, phi = p[o], phi[o]
        # phi_J is a LOWER bound on phi by construction; keep only points above it
        excess = phi - phi_J_lin
        m = (excess > 0) & (p > 0)
        n_pos = int(m.sum())
        if n_pos < MIN_POINTS:
            rows.append({"cv": key[0], "mu": key[1], "seed": key[2], "family": key[3],
                         "b_free": b_free, "b_con": float("nan"), "n_pos": n_pos})
            continue
        b_con = np.polyfit(np.log(p[m]), np.log(excess[m]), 1)[0]
        # CONDITIONING. The whole question is how hard b is pinned by this ladder, so
        # measure db/dphi_J directly: perturb the assumed phi_J by +/- one seed scatter
        # and refit. If a phi_J uncertainty the size of the seed scatter already moves b
        # by more than the width of the bar, the bar was never testable by this design.
        slopes = []
        for d in (-SEED_SCATTER, SEED_SCATTER):
            e2 = phi - (phi_J_lin + d)
            m2 = (e2 > 0) & (p > 0)
            if int(m2.sum()) >= MIN_POINTS:
                slopes.append(np.polyfit(np.log(p[m2]), np.log(e2[m2]), 1)[0])
        # If the +SEED_SCATTER perturbation drives too many excesses negative, only the
        # low side survives. Dropping those ladders BIASES the median low, because they
        # fail exactly where the excess is smallest relative to noise -- i.e. on the MOST
        # ill-conditioned ladders. Fall back to a one-sided derivative and record which
        # was used, rather than discarding them.
        if len(slopes) == 2:
            db_dphi = (slopes[1] - slopes[0]) / (2 * SEED_SCATTER)
            sided = 2
        elif len(slopes) == 1:
            db_dphi = (b_con - slopes[0]) / SEED_SCATTER
            sided = 1
        else:
            db_dphi, sided = float("nan"), 0
        rows.append({"cv": key[0], "mu": key[1], "seed": key[2], "family": key[3],
                     "b_free": b_free, "b_con": float(b_con), "n_pos": n_pos,
                     "db_dphiJ": float(db_dphi), "sided": sided,
                     "b_con_lo": float(slopes[0]) if slopes else float("nan"),
                     "b_con_hi": (float(slopes[1]) if len(slopes) == 2
                                  else float("nan"))})

    good = [r for r in rows if np.isfinite(r["b_con"])]
    lines = ["=" * 78,
             "IS THE EXPONENT ANCHOR'S FAILURE A FIT DEGENERACY?",
             "=" * 78, "",
             "Same ladders, same [0.6, 1.6] bar. The only change is that phi_J is PINNED to",
             "the independent linear-extrapolation estimate instead of being fitted",
             "alongside b -- one fewer free parameter, so a stricter test, not a looser one.",
             ""]
    lines.append("  ladders admitted: %d; with >= %d points above phi_J_lin: %d"
                 % (len(rows), MIN_POINTS, len(good)))
    if len(good) < 5:
        lines.append("")
        lines.append("  INCONCLUSIVE -- too few ladders survive the positivity requirement.")
        lines.append("  EXPONENT_CONSTRAINED_INCONCLUSIVE")
        txt = "\n".join(lines)
        print(txt)
        with open(os.path.join(OUT, "exponent_constrained_verdict.txt"), "w",
                  encoding="utf-8") as fh:
            fh.write(txt + "\n")
        return 0

    # PAIRED set: both estimators returned a value on the same ladder. Comparing the
    # 22-ladder free sample against the 40-ladder constrained sample would be comparing
    # two different samples, not two estimators.
    paired = [r for r in good if np.isfinite(r["b_free"])]
    bf = np.array([r["b_free"] for r in paired])
    bc_paired = np.array([r["b_con"] for r in paired])
    bc = np.array([r["b_con"] for r in good])
    lines.append("")
    lines.append("  PAIRED -- same ladders, both estimators:")
    lines.append("  %-14s %4s %9s %9s %9s" % ("estimator", "n", "median", "mean", "IQR"))
    lines.append("  %-14s %4d %9.3f %9.3f  [%.3f, %.3f]"
                 % ("free phi_J", len(bf), np.median(bf), bf.mean(),
                    np.percentile(bf, 25), np.percentile(bf, 75)))
    lines.append("  %-14s %4d %9.3f %9.3f  [%.3f, %.3f]"
                 % ("pinned phi_J", len(bc_paired), np.median(bc_paired), bc_paired.mean(),
                    np.percentile(bc_paired, 25), np.percentile(bc_paired, 75)))
    dpair = bc_paired - bf
    lines.append("  per-ladder difference (pinned - free): median %+.3f, %d of %d positive"
                 % (float(np.median(dpair)), int((dpair > 0).sum()), len(dpair)))
    lines.append("")
    lines.append("  CONSTRAINED-ONLY set (ladders the free fit could not fit): n = %d, "
                 "median %.3f" % (len(bc) - len(bc_paired),
                                  float(np.median([r["b_con"] for r in good
                                                   if not np.isfinite(r["b_free"])]))
                                  if len(bc) > len(bc_paired) else float("nan")))
    lines.append("  ALL constrained: n = %d, median %.3f" % (len(bc), float(np.median(bc))))
    lines.append("")
    for fam in sorted({r["family"] for r in good}):
        g = np.array([r["b_con"] for r in good if r["family"] == fam])
        lines.append("    %-11s n=%2d  median b_constrained = %.3f" % (fam, len(g), np.median(g)))
    lines.append("")

    med = float(np.median(bc_paired))   # the paired median is the headline
    passes = BAR[0] <= med <= BAR[1]

    # --- the conditioning, which decides whether ANY of this is a measurement ---
    dd_all = np.array([r.get("db_dphiJ", np.nan) for r in paired], float)
    sided = np.array([r.get("sided", 0) for r in paired], float)
    dd = dd_all[np.isfinite(dd_all)]
    n_one = int(((sided == 1) & np.isfinite(dd_all)).sum())
    n_drop = int((~np.isfinite(dd_all)).sum())
    swing = float(np.median(np.abs(dd))) * SEED_SCATTER if dd.size else float("nan")
    lines.append("")
    lines.append("CONDITIONING -- how hard does this ladder actually pin b?")
    lines.append("  |db/dphi_J| median = %.1f per unit phi  (n = %d; %d one-sided, "
                 "%d unusable)" % (float(np.median(np.abs(dd))), dd.size, n_one, n_drop))
    if n_one:
        two = np.abs(dd_all[(sided == 2) & np.isfinite(dd_all)])
        one = np.abs(dd_all[(sided == 1) & np.isfinite(dd_all)])
        lines.append("    two-sided median %.1f over %d; one-sided median %.1f over %d "
                     "-- the one-sided ladders are the MORE ill-conditioned ones, so "
                     "dropping them would bias this LOW"
                     % (float(np.median(two)) if two.size else float("nan"), two.size,
                        float(np.median(one)) if one.size else float("nan"), one.size))
    lines.append("  => a phi_J error of ONE seed scatter (%.5f) moves b by %.3f"
                 % (SEED_SCATTER, swing))
    lines.append("  => the bar is %.1f wide; one seed scatter of phi_J moves b across "
                 "%.0f%% of it" % (BAR[1] - BAR[0], 100 * swing / (BAR[1] - BAR[0])))
    lines.append("PRE-REGISTERED READING")
    lines.append("  PAIRED constrained median b = %.3f against bar [%.1f, %.1f] -> %s"
                 % (med, BAR[0], BAR[1], "PASS" if passes else "FAIL"))
    lines.append("  (free fit on the SAME ladders: %.3f)" % float(np.median(bf)))
    lines.append("")
    lines.append("")
    lines.append("  !! CIRCULARITY WARNING, and it governs the verdict. phi_J_lin is the")
    lines.append("     INTERCEPT OF A LINEAR FIT, i.e. it already assumes b = 1. Pinning")
    lines.append("     phi_J to it and recovering b ~ 1 is therefore NOT independent")
    lines.append("     evidence for the harmonic law, and must not be reported as a PASS.")
    lines.append("     What the constrained fit does show is the SIZE of the trade-off.")
    if passes:
        lines.append("")
        lines.append("  => VERDICT: THE EXPONENT IS UNDETERMINED BY THIS DESIGN, which is")
        lines.append("     different from either passing or failing. Three facts together:")
        lines.append("       (1) phi_J and b are strongly degenerate, corr = +0.675,")
        lines.append("           p = 0.00056;")
        lines.append("       (2) the phi_J differences driving the b shift are TINY --")
        lines.append("           mostly 0.003-0.015, a few seed scatters -- yet b moves by")
        lines.append("           0.58, so the fit is ill-conditioned, not merely noisy;")
        lines.append("       (3) pinning phi_J to a b=1-assuming estimate returns b ~ 1,")
        lines.append("           which is what circularity predicts and cannot confirm.")
        lines.append("     The honest statement is that this p-ladder cannot determine")
        lines.append("     phi_J and b simultaneously. The original FAIL should be restated")
        lines.append("     as NOT MEASURABLE rather than as evidence against the contact")
        lines.append("     law -- and equally, nothing here vindicates the contact law.")
        lines.append("     Settling it needs a phi_J from an INDEPENDENT observable (e.g.")
        lines.append("     an isostatic contact criterion), not from the same p-ladder.")
    else:
        lines.append("  => The failure is NOT the degeneracy. Removing the free phi_J does")
        lines.append("     not recover b ~ 1, so the discrepancy with the harmonic contact")
        lines.append("     law survives and stays an open problem. A6 REMAINS FAILED.")
    lines.append("")
    # The tag must not say PASS: the constrained estimator is circular, so the finding is
    # that the design cannot measure b, not that b has been measured.
    lines.append("  EXPONENT_CONSTRAINED_"
                 + ("UNDETERMINED" if passes else "FAIL"))

    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(OUT, "exponent_constrained_verdict.txt"), "w",
              encoding="utf-8") as fh:
        fh.write(txt + "\n")
    with open(os.path.join(OUT, "exponent_constrained.csv"), "w", newline="",
              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
