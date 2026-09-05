"""closure_sensitivity.py -- what does the CONVERGED closure do to the paper's levers?

WHY THIS EXISTS
---------------
`model/run_bidet.py` reads `C_P` and `n_t` from `cfd_openfoam/cfd_fit_summary.json`,
whose `nominal` block is the **base grid** (refine = 1, 30 720 cells): C_P = 1.3518,
n_t = 1.9023. The three-grid study finished on 2026-09-03 with every grid meeting
`residualControl` and all three at y+ < 1, and it says that number is not converged:

    grid            cells     C_P      n_t
    refine 1        30720     1.3518   1.9023   <- what the paper quotes
    refine 2       122880     1.3928   1.8922
    refine 4       491520     1.5809   1.9338

So the paper's `C_P` is 17 % below the finest grid.

⚠ AND THERE IS NO EXTRAPOLATED VALUE. An earlier version of this file carried a fourth
row, "Richardson extrapolation 1.6333 (p = 2.198, GCI 4.14 %)". That p was the ABSOLUTE
VALUE of -2.198: C_P DIVERGES under refinement -- 1.3518 -> 1.3928 -> 1.5809 with the
steps GROWING, |eps21/eps32| = 4.589 -- so the sequence is running away and there is
nothing to extrapolate to. No quantity in the study has a valid observed order inside the
asymptotic band on a family whose grids all converged. The rows below are three
MEASUREMENTS on three meshes, and the sensitivity they show is a lower bound on how much
the closure moves, not an error bar.

The question this script answers is not "is C_P wrong" -- it is **"does it matter"**.
A closure parameter that moves 21 % is only a problem if a lever the paper quotes moves
with it. The model's optimum has closed form (`bidet_model.optimum`), so the sensitivity
is exact rather than estimated, and it can be read off in one run.

#############################################################################
# READ THIS BEFORE QUOTING ANY Q* FROM THIS SCRIPT (added 2026-09-03 after a
# physics-referee pass, and it is the reason the write-up quotes a DIRECTION
# and not the number).
#
# `optimum()` solves for the nozzle diameter d*, then evaluates the closure at
# the FIXED cleaning radius r_max = 15 mm. Because d* shrinks as C_P grows,
# every row of the table below sits at
#
#       r_max / d*  =  22.95, 23.45, 25.67
#
# against a closure FITTED over r/d in [2, 20] (fit_closure.py's fit_lo/fit_hi)
# and a CFD domain that ends at Rdom = 25 mm, i.e. r/d = 25 at the simulated
# d = 1 mm. So every design point extrapolates past the fit band, and the two
# finest extrapolates past the edge of the computation itself.
#
# The extrapolation is also NOT uniform across the comparison: larger C_P gives
# smaller d* gives worse extrapolation, so the "corrected" rows stand on the
# thinnest ice. The 23.6 % is therefore a sensitivity, not a corrected value.
#
# Note run_bidet.py's own [Result 1] does NOT have this problem -- it holds
# d = 1 mm so r/d = 15, in band. The discipline exists; optimum() escapes it.
# Fixing it properly means extending Rdom to ~50 mm and refitting.
#############################################################################

WHAT IS SWEPT, and what deliberately is not
-------------------------------------------
Only `C_P`, across the three measured values above. `n_t` is held at the base-grid 1.9023 in every
case **on purpose**: it is non-monotone in the grid (1.9023 / 1.8922 / 1.9338, oscillating
within 2.2 %) so there is no converged value to move it to, and varying two parameters
would confound which one carries the change. `n_t`'s own 2.2 % band is reported separately
at the end so the two are not conflated.

    python closure_sensitivity.py
"""

from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "model"))

from bidet_model import RHO, V0_from_dP, cfd_closure, optimum   # noqa: E402

# The three MEASURED values, all from grids that reached residualControl.
#
# ⚠ THE RICHARDSON-EXTRAPOLATION ROW (1.6333) HAS BEEN REMOVED, 2026-09-03.
# It was quoted from an observed order of p = 2.198 which was the ABSOLUTE VALUE of
# -2.198: C_P is DIVERGING under refinement (|eps21/eps32| = 4.589 -- the gap between the
# two finest grids is 4.6x the gap between the two coarsest), so there is nothing for it
# to extrapolate TO. gci.py now refuses that triplet outright. Any "extrapolated C_P" in
# an older artifact is void.
CP_CASES = [
    ("base grid, refine 1 (what the paper quotes)", 1.3518),
    ("refine 2, converged",                          1.3928),
    ("refine 4, converged (finest)",                  1.5809),
]
NT_FIXED = 1.9023
NT_BAND = (1.8922, 1.9338)      # the non-monotone spread, for the note at the end

R_MAX, TAU_C, DP = 0.015, 5.0, 0.07e6

# The closure's own validity window, so the out-of-band check below is not a magic number.
# FIT_LO / FIT_HI are fit_closure.py's defaults (its argv[2], argv[3]); RDOM_MM is
# gen_case.py's RDOM_DEFAULT, i.e. r/d = 25 at the simulated d = 1 mm.
FIT_LO, FIT_HI = 2.0, 20.0
RDOM_MM = 25.0


def main() -> int:
    print("=" * 78)
    print("Closure sensitivity: what the converged C_P does to the paper's levers")
    print("=" * 78)
    print("n_t held at %.4f throughout (it is NON-MONOTONE in the grid, so there is no"
          % NT_FIXED)
    print("converged value to move it to; its own band is reported at the end).")
    print("Nominal: r_max = %.0f mm, tau_c = %.1f Pa, dP = %.2f MPa\n"
          % (1e3 * R_MAX, TAU_C, DP / 1e6))

    base = None
    print("%-46s %8s %9s %9s %9s %11s" % ("C_P from", "C_P", "d* [mm]", "Q* [L/min]",
                                          "vs base", "r_max/d*"))
    rows = []
    for label, cp in CP_CASES:
        o = optimum(cfd_closure(cp, NT_FIXED), DP, R_MAX, TAU_C)
        q_lmin = o["Q_Lmin"]
        if base is None:
            base = q_lmin
        rd = R_MAX / o["d_star"]
        rows.append((label, cp, o["d_star"] * 1e3, q_lmin, q_lmin / base, rd))
        print("%-46s %8.4f %9.4f %9.4f %8.3fx %8.2f %s"
              % (label, cp, o["d_star"] * 1e3, q_lmin, q_lmin / base, rd,
                 "OK" if FIT_LO <= rd <= FIT_HI else "** OUT OF BAND **"))
    if any(not (FIT_LO <= r[5] <= FIT_HI) for r in rows):
        print()
        print("*** EVERY DESIGN POINT IS OUTSIDE THE CLOSURE'S FIT BAND [%.0f, %.0f]."
              % (FIT_LO, FIT_HI))
        print("    The domain also ends at r/d = %.0f (Rdom = %.0f mm at d = 1 mm), so the"
              % (RDOM_MM / 1.0, RDOM_MM))
        print("    finest rows evaluate past the edge of the CFD that produced the fit.")
        print("    Quote the DIRECTION of the change below, not the value.")

    dq = rows[-1][3] / rows[0][3]
    print()
    print("VERDICT")
    print("-------")
    print("Water demand Q* moves %.3fx between the base grid the paper quotes and the" % dq)
    print("FINEST GRID MEASURED -- a %+.1f%% change in the number the whole paper"
          % (100 * (dq - 1)))
    print("optimises. This is a measured spread across meshes, NOT an error bar: C_P is")
    print("diverging under refinement, so the true continuum value is not bracketed by")
    print("these three and could lie further out still.")
    if abs(dq - 1) < 0.05:
        print("That is inside the 5%% band the paper already carries for other sources,")
        print("so the conclusions survive the correction and only the quoted C_P changes.")
    else:
        print("That is NOT a rounding correction. The optimum itself moves and the")
        print("manuscript's quoted Q* has to be recomputed, not just caveated.")

    # why it is so much smaller than the 21% on C_P: the closed form
    print()
    e_d = 1.0 / (NT_FIXED - 0.5)
    e_q = 2.0 / (NT_FIXED - 0.5)
    cp_ratio = CP_CASES[-1][1] / CP_CASES[0][1]
    print("THE SENSITIVITY IS AN EXPONENT, AND IT AMPLIFIES")
    print("------------------------------------------------")
    print("optimum() solves  d* = [tau_c r_max^n nu^p / (C_tau rho V0^(2+p))]^(1/(p+n))")
    print("with p = -1/2 and n = n_t, so")
    print("    d* ~ C_P^(-1/(n_t - 1/2)) = C_P^(-%.4f)      -- damped, exponent < 1" % e_d)
    print("    Q* ~ d*^2 V0 ~ C_P^(-%.4f)                  -- AMPLIFIED, exponent > 1"
          % e_q)
    print("C_P rises %.1f%% between the base grid and the finest measured, and Q* falls"
          % (100 * (cp_ratio - 1)))
    print("%.1f%%: %.4f^(-%.4f) = %.4f. The exponent 2/(n_t - 1/2) = %.3f is the whole"
          % (100 * (1 - cp_ratio ** (-e_q)), cp_ratio, e_q, cp_ratio ** (-e_q), e_q))
    print("story, and because it EXCEEDS 1 the closure error is magnified rather than")
    print("absorbed. It would be worse still at smaller n_t: the exponent diverges as")
    print("n_t -> 1/2, which is the regime the LAMINAR branch (n = 5) is nowhere near but")
    print("a shallower turbulent decay would be. Worth stating in the paper, because it")
    print("says how much closure accuracy the design model actually requires.")

    print()
    print("SEPARATELY, n_t's own band")
    print("--------------------------")
    for nt in (NT_BAND[0], NT_FIXED, NT_BAND[1]):
        o = optimum(cfd_closure(1.5809, nt), DP, R_MAX, TAU_C)
        print("  n_t = %.4f -> Q* = %.4f L/min" % (nt, o["Q_Lmin"]))
    print("  (all at the finest-grid C_P = 1.5809; the n_t spread is non-monotone in the")
    print("   grid so it is an uncertainty band, not a correction.)")
    qs = [optimum(cfd_closure(1.5809, nt), DP, R_MAX, TAU_C)["Q_Lmin"]
          for nt in (NT_BAND[0], NT_BAND[1])]
    print("  -> n_t's 2.2%% grid band alone spans Q* by %.3fx (%.4f to %.4f L/min),"
          % (max(qs) / min(qs), min(qs), max(qs)))
    print("     which is COMPARABLE to the whole C_P correction. So the honest statement")
    print("     is that Q* is uncertain at the ~20-25% level from the closure alone, and")
    print("     the single number the paper currently quotes for it is not supported.")

    out = os.path.join(HERE, "closure_sensitivity.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump({"n_t_fixed": NT_FIXED, "r_max_m": R_MAX, "tau_c_Pa": TAU_C,
                   "dP_Pa": DP, "V0_ms": V0_from_dP(DP), "rho": RHO,
                   "rows": [{"label": r[0], "C_P": r[1], "d_star_mm": r[2],
                             "Q_star_Lmin": r[3], "Q_ratio_vs_base": r[4]}
                            for r in rows],
                   "Q_ratio_base_to_extrapolated": dq}, fh, indent=2)
    print("\nwrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
