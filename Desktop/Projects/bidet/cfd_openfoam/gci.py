#!/usr/bin/env python3
"""gci.py -- three-grid discretisation-error study for the impinging-jet closure.

WHY THIS EXISTS
---------------
The paper's grid check was two grids and a percentage difference: base vs refine = 2,
giving 0.5 % on `n_t`, 6.5 % on `tau_w(15 mm)` and 12.3 % on `C_P`. A difference between
two grids is not an error estimate -- it has no observed order and no extrapolated value,
so it cannot say whether the finest grid is close to the continuum answer or merely close
to the next-coarsest one.

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

Here the refinement is uniform by an integer factor in both mesh directions, so the
representative cell size goes as 1/refine exactly. That is asserted against the cell
counts rather than assumed, because a grading change would break it silently.

WHAT THE PROCEDURE CANNOT DO, and this is reported rather than hidden
--------------------------------------------------------------------
`p` is only meaningful where the three values are MONOTONE. If `s < 0` the solution is
oscillatory in the grid and Celik et al. say so explicitly: the extrapolation is not
valid, and the honest output is the oscillation amplitude, not a GCI. This script prints
that verdict instead of a number when it happens, and it flags any `p` outside [0.5, 2.5]
as outside the asymptotic range -- the same bar the frisbee and drop-impact ladders use.

    python3 gci.py            # reads the three case dirs, writes convergence/
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "convergence")

# grid 1 is the FINEST, per Celik et al.'s indexing
CASES = [
    ("case_Re10000_H6_r4", 4.0),
    ("case_Re10000_H6_r2", 2.0),
    ("case_Re10000_H6", 1.0),
]

# (key in the fit_closure record, label, units)
QOI = [
    ("n_t", "closure exponent n_t", ""),
    ("C_P", "closure prefactor C_P", ""),
    ("tau_at_15mm_Pa", "tau_w at r = 15 mm", "Pa"),
    ("tau_peak_Pa", "peak tau_w", "Pa"),
]

ASYMPTOTIC_BAND = (0.5, 2.5)


def n_cells(case: str) -> int:
    """Cell count from blockMesh's own log -- not recomputed from gen_case.py."""
    path = os.path.join(HERE, case, "log.blockMesh")
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "cells:" in line:
                return int(line.split(":")[1].strip())
    raise RuntimeError("no cell count in " + path)


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
        return float("nan"), s, False, ("non-monotone in the grid (sign of eps32/eps21 "
                                        "is negative); Richardson extrapolation invalid")
    p = 1.0
    for _ in range(itmax):
        q = 0.0 if abs(r21 - r32) < 1e-12 else math.log(
            (r21**p - s) / (r32**p - s))
        p_new = abs(math.log(abs(ratio)) + q) / math.log(r21)
        if abs(p_new - p) < 1e-10:
            p = p_new
            break
        p = p_new
    return p, s, True, ""


def main() -> int:
    os.makedirs(OUT, exist_ok=True)
    recs, cells = [], []
    for case, ref in CASES:
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

    # representative cell size. For this 2-D wedge mesh h ~ sqrt(A/N); the refinement is
    # uniform by an integer factor, so h ~ 1/refine. Check the two agree before using it.
    h_from_ref = [1.0 / ref for _, ref in CASES]
    h_from_N = [math.sqrt(1.0 / n) for n in cells]
    scale = h_from_ref[0] / h_from_N[0]
    for a, b in zip(h_from_ref, h_from_N):
        assert abs(a - b * scale) / a < 0.02, (
            "cell counts do not scale as refine^2 -- the meshes are not a uniformly "
            "refined family and the GCI procedure does not apply: N = %s" % (cells,))
    h = h_from_ref
    r21, r32 = h[1] / h[0], h[2] / h[1]

    lines = []
    lines.append("# Grid convergence, impinging jet at Re_d = 1e4, H/d = 6\n")
    lines.append("Three uniformly refined meshes, ratio r = %.3f / %.3f.\n" % (r21, r32))
    lines.append("| grid | refine | cells | " + " | ".join(q[0] for q in QOI) + " |")
    lines.append("|---|---|---|" + "---|" * len(QOI))
    for i, ((case, ref), rec) in enumerate(zip(CASES, recs)):
        lines.append("| %d %s | %g | %d | " % (i + 1, "(fine)" if i == 0 else "", ref,
                                               cells[i])
                     + " | ".join("%.4f" % rec[q[0]] for q in QOI) + " |")
    lines.append("")

    rows = [("quantity", "phi1_fine", "phi2", "phi3", "p", "phi_extrap",
             "e_a21_%", "e_ext21_%", "GCI_fine_%", "verdict")]
    csv = ["quantity,phi1,phi2,phi3,p,phi_extrap,e_a21,e_ext21,gci_fine,verdict"]
    for key, label, unit in QOI:
        phi1, phi2, phi3 = (recs[0][key], recs[1][key], recs[2][key])
        p, s, ok, why = observed_order(phi1, phi2, phi3, r21, r32)
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
        band = ASYMPTOTIC_BAND
        verdict = ("OK" if band[0] <= p <= band[1] else
                   "p = %.2f OUTSIDE [%g, %g]: not in the asymptotic range, GCI is "
                   "indicative only" % (p, band[0], band[1]))
        rows.append((label, "%.5g" % phi1, "%.5g" % phi2, "%.5g" % phi3, "%.3f" % p,
                     "%.5g" % ext, "%.2f" % (100 * e_a), "%.2f" % (100 * e_ext),
                     "%.2f" % (100 * gci), verdict))
        csv.append("%s,%.6g,%.6g,%.6g,%.4f,%.6g,%.5g,%.5g,%.5g,%s"
                   % (key, phi1, phi2, phi3, p, ext, e_a, e_ext, gci, verdict))

    w = [max(len(str(r[i])) for r in rows) for i in range(len(rows[0]))]
    lines.append("| " + " | ".join(str(c).ljust(w[i]) for i, c in enumerate(rows[0]))
                 + " |")
    lines.append("|" + "|".join("-" * (w[i] + 2) for i in range(len(w))) + "|")
    for r in rows[1:]:
        lines.append("| " + " | ".join(str(c).ljust(w[i]) for i, c in enumerate(r)) + " |")
    lines.append("")
    lines.append("Procedure: Celik et al., J. Fluids Eng. 130(7):078001 (2008), "
                 "Fs = 1.25 for three grids.")

    md = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "gci_Re10000_H6.md"), "w", encoding="utf-8") as fh:
        fh.write(md)
    with open(os.path.join(OUT, "gci_Re10000_H6.csv"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(csv) + "\n")
    print(md)
    print("wrote %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
