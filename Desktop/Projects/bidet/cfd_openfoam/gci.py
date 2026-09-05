#!/usr/bin/env python3
"""gci.py -- three-grid discretisation-error study for the impinging-jet closure.

WHY THIS EXISTS
---------------
The paper's grid check was two grids and a percentage difference: base vs refine = 2,
giving 0.5 % on `n_t`, 6.5 % on `tau_w(15 mm)` and 12.3 % on `C_P`. A difference between
two grids is not an error estimate -- it has no observed order and no extrapolated value,
so it cannot say whether the finest grid is near the continuum answer or merely near the
next-coarsest one.

That matters for the venue. The ASME *Journal of Fluids Engineering* editorial policy on
numerical accuracy is a hard gate, not a preference: it "will not accept for publication
any paper reporting the numerical solution of a fluids engineering problem that fails to
address the task of systematic truncation error testing and accuracy estimation", and it
names Richardson extrapolation / the Grid Convergence Index as the accepted procedure.
Two grids cannot satisfy it. Three can.

THE PROCEDURE
-------------
Celik, Ghia, Roache, Freitas, Coleman & Raad, "Procedure for estimation and reporting of
uncertainty due to discretization in CFD applications", J. Fluids Eng. 130(7):078001
(2008), which is the policy's own reference. With grids ordered 1 = finest,

    r21 = h2/h1,  r32 = h3/h2
    eps21 = phi2 - phi1,  eps32 = phi3 - phi2,  s = sign(eps32/eps21)
    p  solves  p = |ln|eps32/eps21| + q(p)| / ln(r21),
               q(p) = ln((r21^p - s)/(r32^p - s))            [q = 0 when r21 == r32]
    phi_ext21   = (r21^p * phi1 - phi2) / (r21^p - 1)
    e_a21       = |(phi1 - phi2)/phi1|                        approximate relative error
    e_ext21     = |(phi_ext21 - phi1)/phi_ext21|              extrapolated relative error
    GCI_fine21  = 1.25 * e_a21 / (r21^p - 1)

The refinement is uniform by the same factor in both mesh directions, so the
representative cell size goes as 1/refine exactly. That is ASSERTED against the cell
counts rather than assumed, because a grading change would break it silently.

TWO FAMILIES, and why both
--------------------------
    --family fine     refine 4 / 2 / 1     GCI reported on the 491 520-cell mesh
    --family coarse   refine 2 / 1 / 0.5   GCI reported on the 122 880-cell mesh
                                           (the paper's existing fine grid)

`fine` is the better study and is the one to quote. `coarse` costs minutes rather than
hours and brackets the same fine grid from below, so it is the fallback if the refine = 4
run has not finished -- and running both is itself a check: the observed order should not
depend on which triplet it was measured from, and if it does, the meshes are not in the
asymptotic range and neither GCI means much.

WHAT THE PROCEDURE CANNOT DO, and this is reported rather than hidden
--------------------------------------------------------------------
1. `p` is only meaningful where the three values are MONOTONE. If `s < 0` the solution is
   oscillatory in the grid and Celik et al. exclude the case explicitly: the extrapolation
   is invalid, and the honest output is the oscillation amplitude, not a GCI. This script
   prints that verdict instead of a number when it happens, and flags any `p` outside
   [0.5, 2.5] as outside the asymptotic range.
2. A GCI measures DISCRETISATION error and silently assumes ITERATIVE error is far
   smaller. That assumption is false here unless checked, so the script reads each run's
   final residual per field out of `log.simpleFoam` and prints it beside the GCI, together
   with whether `residualControl` actually fired. A grid that stopped at the iteration cap
   with a residual comparable to the differences being extrapolated does not support a GCI,
   and the reader is told so rather than left to assume.

    python3 gci.py [--family fine|coarse]
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "convergence")

# grid 1 is the FINEST, per Celik et al.'s indexing
FAMILIES = {
    "fine":   [("case_Re10000_H6_r4", 4.0), ("case_Re10000_H6_r2", 2.0),
               ("case_Re10000_H6", 1.0)],
    "coarse": [("case_Re10000_H6_r2", 2.0), ("case_Re10000_H6", 1.0),
               ("case_Re10000_H6_r0", 0.5)],
    # Same triplet, but with the refine = 2 grid taken from the run that was allowed to
    # reach residualControl (extend_r2.sh, converged at iteration 5105) instead of the
    # one that stopped at the 3000-iteration cap. Comparing this family against `coarse`
    # is how much of the "grid effect" was really iteration count.
    "coarse_conv": [("case_Re10000_H6_r2_ext", 2.0), ("case_Re10000_H6", 1.0),
                    ("case_Re10000_H6_r0", 0.5)],
    # The best-resolved family available: refine 4 / 2 / 1, all three meeting
    # residualControl (iterations 11552 / 5105 / 2261).
    #
    # ⚠⚠ ITS p = 2.198 ON C_P IS NOT AN ESTABLISHED OBSERVED ORDER. Three things, all
    # found by a physics-referee pass on 2026-09-03 and all verified:
    #   (1) An earlier comment here explained the coarse family's failure by y+ straddling
    #       1. THAT IS REFUTED. `yplus_profile.py` resolves y+ radially: the 1.52 median is
    #       entirely a near-stagnation effect (r/d 2-5), and at r/d 14-20 -- where
    #       tau_w(15mm) lives and where the non-monotonicity is -- y+ is 0.04 to 0.30 on
    #       EVERY grid including the coarsest. So the difference between `coarse_conv` and
    #       this family is a change of TRIPLET, not the repair of a mechanism.
    #   (2) p = 2.198 exceeds the formal order of the discretisation: fvSchemes has
    #       div(phi,k) and div(phi,omega) as `bounded Gauss upwind`, FIRST order, and nut
    #       (hence the wall shear) depends on both.
    #   (3) p on C_P is 1.703 / 0.315 / 2.198 across `coarse` / `coarse_conv` / this --
    #       a factor of seven from the same handful of grids. A 3-point Richardson fit has
    #       zero residual degrees of freedom and cannot self-diagnose; cross-triplet
    #       comparison is the only check available and it does not support 2.198.
    #
    # ⚠ Separately: watch the refine = 4 run rather than judging it early. Its Uy residual
    # sat flat at ~4.0e-5 from iteration 6400 to 8000 and looked like a stall; it fell
    # through 1e-5 by 11552.
    "fine_conv": [("case_Re10000_H6_r4", 4.0), ("case_Re10000_H6_r2_ext", 2.0),
                  ("case_Re10000_H6", 1.0)],
    # THE INDEPENDENT ORDER CHECK, and it CANNOT CURRENTLY RUN. Ratio sqrt(2) throughout
    # (above Celik et al.'s recommended 1.3), meshes 122880 / 61009 / 30720 -- a triplet
    # disjoint from `fine_conv` except at the ends, which is exactly what is needed to
    # decide whether p = 2.198 means anything.
    #
    # 🔴 case_Re10000_H6_r1_ext ran 20000 iterations and did NOT converge: every field
    # reaches ~1e-9 EXCEPT Uz at 1.8e-2. Uz is the out-of-plane component on a
    # one-cell-thick wedge and should be ~0 (it is 1.6e-6 to 4.9e-6 on the three grids
    # that did converge), so this is a defect of that case, not slowness. Suspect the odd
    # cell counts at refine = 1.414 (nx = 169) against the collapsed wedge axis.
    # Fixing it is the cheapest open question in the study.
    "midfine": [("case_Re10000_H6_r2_ext", 2.0), ("case_Re10000_H6_r1_ext", 1.414),
                ("case_Re10000_H6", 1.0)],
}

# (key in the fit_closure record, label)
QOI = [
    ("n_t", "closure exponent n_t"),
    ("C_P", "closure prefactor C_P"),
    ("tau_at_15mm_Pa", "tau_w at r = 15 mm [Pa]"),
    ("tau_peak_Pa", "peak tau_w [Pa]"),
]

# [UNVERIFIED PROVENANCE, flagged 2026-09-03]. This window is NOT part of Celik et al.'s
# published procedure; it is carried over from this lab's own convergence-study practice on
# the frisbee and drop-impact ladders. It affects only the OK/outside label printed beside a
# GCI, never the GCI value. Resolve it to a source or restate it as a house rule.
#
# ⚠ And a p inside it is not thereby credible: this case's fvSchemes uses `bounded Gauss
# upwind` (first order) for k and omega, on which nut and hence the wall shear depend, so an
# observed order above ~1 is above the formal order of a link in its own causal chain.
ASYMPTOTIC_BAND = (0.5, 2.5)
RESIDUAL_TARGET = 1e-5          # the residualControl set in system/fvSolution


def n_cells(case: str) -> int:
    """Cell count from blockMesh's own log -- not recomputed from gen_case.py."""
    # A restarted case (`_ext`) has no blockMesh log of its own -- it reuses the mesh --
    # so fall back to the sibling it was restarted from. Never guess from gen_case.py's
    # arithmetic: the whole point of reading the log is that the mesh is what blockMesh
    # actually built.
    candidates = [os.path.join(HERE, case, "log.blockMesh")]
    if case.endswith("_ext"):
        candidates.append(os.path.join(HERE, case[:-4], "log.blockMesh"))
    for path in candidates:
        if not os.path.exists(path):
            continue
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                # blockMesh prints "  nCells: 122880". Match the label exactly rather
                # than a bare "cells:", which also hits the polyMesh constructor
                # signature printed above it.
                if line.strip().startswith("nCells:"):
                    return int(line.split(":")[1].strip())
    raise RuntimeError("no nCells line found for %s (tried %s)" % (case, candidates))


def iterative_state(case: str) -> dict:
    """Last initial-residual per field, final iteration, and whether SIMPLE converged."""
    path = os.path.join(HERE, case, "log.simpleFoam")
    last = {}
    n_iter, converged = None, False
    pat = re.compile(r"Solving for (\w+),\s+Initial residual = ([0-9.eE+-]+)")
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = pat.search(line)
            if m:
                last[m.group(1)] = float(m.group(2))
            elif line.startswith("Time = "):
                n_iter = int(float(line.split("=")[1]))
            elif "SIMPLE solution converged" in line:
                converged = True
    worst = max(last.values()) if last else float("nan")
    return {"residuals": last, "iterations": n_iter, "converged": converged,
            "worst_residual": worst}


def fit_record(case: str) -> dict:
    """Re-run fit_closure.py so every record in the study comes from one code path."""
    r = subprocess.run([sys.executable, os.path.join(HERE, "fit_closure.py"), case],
                       capture_output=True, text=True, cwd=HERE)
    if r.returncode != 0:
        raise RuntimeError("fit_closure failed for %s: %s" % (case, r.stderr[-800:]))
    return json.loads(r.stdout.strip().splitlines()[-1])


def observed_order(phi1, phi2, phi3, r21, r32, tol=1e-12, itmax=200):
    """Solve Celik et al. Eq. (2) by fixed-point iteration. Returns (p, s, ok, why)."""
    eps21, eps32 = phi2 - phi1, phi3 - phi2
    if abs(eps21) < tol:
        return float("nan"), 0.0, False, "phi1 == phi2 to machine precision"
    ratio = eps32 / eps21
    s = math.copysign(1.0, ratio)
    if s < 0:
        return float("nan"), s, False, ("non-monotone in the grid (sign of eps32/eps21 is "
                                        "negative); Richardson extrapolation invalid")
    p = 1.0
    for _ in range(itmax):
        q = 0.0 if abs(r21 - r32) < 1e-12 else math.log((r21**p - s) / (r32**p - s))
        p_new = abs(math.log(abs(ratio)) + q) / math.log(r21)
        if abs(p_new - p) < 1e-10:
            p = p_new
            break
        p = p_new
    return p, s, True, ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", choices=sorted(FAMILIES), default="fine")
    args = ap.parse_args()
    cases = FAMILIES[args.family]

    os.makedirs(OUT, exist_ok=True)
    recs, cells, iters = [], [], []
    for case, ref in cases:
        if not os.path.isdir(os.path.join(HERE, case)):
            print("MISSING case dir: %s -- run `bash run_case.sh 10000 6 %g` first"
                  % (case, ref))
            return 1
        rec = fit_record(case)
        if "error" in rec:
            print("case %s has no usable sample: %s" % (case, rec["error"]))
            return 1
        recs.append(rec)
        cells.append(n_cells(case))
        iters.append(iterative_state(case))

    # representative cell size. For this 2-D wedge mesh h ~ sqrt(A/N); the refinement is
    # uniform by the same factor in both directions, so h ~ 1/refine. Check they agree.
    h_from_ref = [1.0 / ref for _, ref in cases]
    h_from_N = [math.sqrt(1.0 / n) for n in cells]
    scale = h_from_ref[0] / h_from_N[0]
    for a, b in zip(h_from_ref, h_from_N):
        assert abs(a - b * scale) / a < 0.02, (
            "cell counts do not scale as refine^2 -- the meshes are not a uniformly "
            "refined family and the GCI procedure does not apply: N = %s" % (cells,))
    h = h_from_ref
    r21, r32 = h[1] / h[0], h[2] / h[1]

    L = []
    L.append("# Grid convergence, impinging jet at Re_d = 1e4, H/d = 6 (`--family %s`)\n"
             % args.family)
    L.append("Three uniformly refined meshes, refinement ratios r21 = %.3f, r32 = %.3f.\n"
             % (r21, r32))
    L.append("| grid | refine | cells | " + " | ".join(q[0] for q in QOI) + " |")
    L.append("|---|---|---|" + "---|" * len(QOI))
    for i, ((case, ref), rec) in enumerate(zip(cases, recs)):
        L.append("| %d %s | %g | %d | " % (i + 1, "(fine)" if i == 0 else "", ref, cells[i])
                 + " | ".join("%.4f" % rec[q[0]] for q in QOI) + " |")
    L.append("")

    # --- iterative convergence, which a GCI silently assumes is negligible -------------
    L.append("## Iterative convergence\n")
    L.append("A GCI measures discretisation error and assumes iterative error is far "
             "smaller. Residual target in `fvSolution` is %g." % RESIDUAL_TARGET)
    L.append("")
    L.append("| grid | iterations | `SIMPLE converged` | worst final residual | fields |")
    L.append("|---|---|---|---|---|")
    for (case, ref), it in zip(cases, iters):
        flds = ", ".join("%s %.2g" % (k, v) for k, v in sorted(it["residuals"].items()))
        L.append("| refine %g | %s | %s | **%.2g**%s | %s |"
                 % (ref, it["iterations"], "yes" if it["converged"] else "**no (cap)**",
                    it["worst_residual"],
                    "" if it["worst_residual"] <= RESIDUAL_TARGET else " ⚠ above target",
                    flds))
    L.append("")

    rows = [("quantity", "phi1_fine", "phi2", "phi3", "p", "phi_extrap",
             "e_a21_%", "e_ext21_%", "GCI_fine_%", "verdict")]
    csv = ["quantity,phi1,phi2,phi3,p,phi_extrap,e_a21,e_ext21,gci_fine,verdict"]
    for key, label in QOI:
        phi1, phi2, phi3 = recs[0][key], recs[1][key], recs[2][key]
        p, _s, ok, why = observed_order(phi1, phi2, phi3, r21, r32)
        if not ok:
            osc = abs(phi2 - phi1) / abs(phi1) * 100.0
            rows.append((label, "%.5g" % phi1, "%.5g" % phi2, "%.5g" % phi3,
                         "--", "--", "%.2f" % osc, "--", "--", "NOT CONVERGENT: " + why))
            csv.append("%s,%.6g,%.6g,%.6g,,,%.4g,,,%s" % (key, phi1, phi2, phi3, osc, why))
            continue
        ext = (r21**p * phi1 - phi2) / (r21**p - 1.0)
        e_a = abs((phi1 - phi2) / phi1)
        e_ext = abs((ext - phi1) / ext) if ext != 0 else float("nan")
        gci = 1.25 * e_a / (r21**p - 1.0)
        lo, hi = ASYMPTOTIC_BAND
        verdict = ("OK" if lo <= p <= hi else
                   "p = %.2f OUTSIDE [%g, %g]: not in the asymptotic range, GCI is "
                   "indicative only" % (p, lo, hi))
        rows.append((label, "%.5g" % phi1, "%.5g" % phi2, "%.5g" % phi3, "%.3f" % p,
                     "%.5g" % ext, "%.2f" % (100 * e_a), "%.2f" % (100 * e_ext),
                     "%.2f" % (100 * gci), verdict))
        csv.append("%s,%.6g,%.6g,%.6g,%.4f,%.6g,%.5g,%.5g,%.5g,%s"
                   % (key, phi1, phi2, phi3, p, ext, e_a, e_ext, gci, verdict))

    L.append("## Discretisation error\n")
    w = [max(len(str(r[i])) for r in rows) for i in range(len(rows[0]))]
    L.append("| " + " | ".join(str(c).ljust(w[i]) for i, c in enumerate(rows[0])) + " |")
    L.append("|" + "|".join("-" * (w[i] + 2) for i in range(len(w))) + "|")
    for r in rows[1:]:
        L.append("| " + " | ".join(str(c).ljust(w[i]) for i, c in enumerate(r)) + " |")
    L.append("")
    L.append("Procedure: Celik et al., J. Fluids Eng. 130(7):078001 (2008); Fs = 1.25 for "
             "three grids.")

    md = "\n".join(L) + "\n"
    stem = os.path.join(OUT, "gci_Re10000_H6_%s" % args.family)
    with open(stem + ".md", "w", encoding="utf-8") as fh:
        fh.write(md)
    with open(stem + ".csv", "w", encoding="utf-8") as fh:
        fh.write("\n".join(csv) + "\n")
    print(md)
    print("wrote %s.md / .csv" % stem)
    return 0


if __name__ == "__main__":
    sys.exit(main())
