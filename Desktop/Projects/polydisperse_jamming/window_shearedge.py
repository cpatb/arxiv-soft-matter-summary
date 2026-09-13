"""window_shearedge.py — the DST window recomputed with a SHEAR-measured lower edge.

THE CLAIM THIS RE-READS
-----------------------
`dst_window.py` measures `W = phi_0 - phi_mu` from two ISOTROPIC ladders and reports

    "IS W CONSTANT IN CV?  dW/dCV = -0.01217 +- 0.00333, t = -3.65 on n = 6 levels.
     RESOLVED as a CHANGE: over CV 0.0 -> 0.5 the window moves by -0.0061 ... -11 %"

and its own docstring already says why that may be the wrong lower edge:

    "The frictionless identification is the safer of the two ...; the frictional one is
     NOT safe -- the project's own shear campaign puts the frictional critical state
     ABOVE the isotropic frictional phi_J."

In Wyart-Cates, `phi_mu` is the jamming point of the SHEARED suspension in the
high-stress (frictional) limit. Until 2026-09-12 this project had no such measurement
that did not read `phi_J` somewhere. It has one now: the rate-matched constant-pressure
ladder returns `phi_c(p*)` at three CV levels from a route with no `phi_J` in it
(`out/constp_ratematch.txt`). This script asks what the window does when that number is
used as the lower edge, and it answers three separate questions that are easy to blur:

  Q1  does the SIGN of dW/dCV survive the change of edge?
  Q2  is either sign RESOLVED against the error the two campaigns actually carry?
  Q3  does the TRAVERSE -- dphi_mu/dCV measured in windows, the project's load-bearing
      number -- survive it?

THE GAUGE PROBLEM, AND HOW IT IS PRICED RATHER THAN ASSERTED
------------------------------------------------------------
`phi_0` is a `p -> 0` intercept (`phi = phi_J + c*p_hat`, `phase2_phiJ.py`). `phi_c(p*)`
is measured AT `p*` = 3e-04. Mixing them puts a `c*p*` bias in the lower edge, and `c`
is NOT CV-independent on the frictional ladders (20.7 -> 12.5 over CV 0 -> 0.5), so the
bias does not cancel out of the contrast. The critical-state compressibility has never
been measured here -- the four-pressure constant-pressure ladder is
`NEXT_STEPS_2026-09-13.md` §2 item 2 -- so the isotropic frictional `c(CV)` is used as a
PROXY and the reading is reported both ways, uncorrected and corrected. Which one is
right is an open measurement, and saying so is the result.

Reads committed artifacts only. Writes out/window_shearedge.txt. No simulation.
"""

from __future__ import annotations

import math
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FITS = os.path.join(OUT, "phase2_fits.csv")

# --- the shear-measured lower edge --------------------------------------------------
# out/constp_ratematch.txt, "blocks": phi_c pooled over both phi starts, 6 cells each,
# rate-matched to gdot = 1.17e-04 (spread 1.001x), N = 1000, p* = 3e-04.
PHI_C_SHEAR = {0.0: 0.591127, 0.2: 0.601989, 0.5: 0.644447}
PHI_C_SHEAR_SD = {0.0: 0.001020, 0.2: 0.001313, 0.5: 0.004336}
PHI_C_SHEAR_N = {0.0: 6, 0.2: 6, 0.5: 6}
P_STAR = 3e-4

N_MATCH = 1500            # the frictionless arm exists at this N only
FAMILY = "lognormal"
CVS = (0.0, 0.2, 0.5)     # the three levels the shear campaign covers


def ladder(mu: float) -> pd.DataFrame:
    """Aggregate the isotropic ladder at this friction, BY CV.

    ⚠ It loads through `dst_window.load_fits` rather than reading the csv directly.
    Reading it directly is what this script did first, and it silently pooled the
    `capfix` re-runs with the originals -- the project's own standing hazard, "the
    HIGHEST cap wins, cell by cell". That gave phi_0(CV 0) = 0.632763 on 6 rows where
    the correct value is 0.632371 on 3, and it moved phi_0(CV 0.5) by 0.0018, which is
    a third of the effect this script is about. Two scripts computing the same quantity
    must be made to agree before either is believed, so this one does not compute it.
    """
    d = pd.read_csv(FITS)
    d = d[(d.n_requested == N_MATCH) & (d.family == FAMILY) & (d.mu_p == mu)]
    d = d[np.isfinite(d.phi_J_lin)]
    # the HIGHEST cap wins, cell by cell -- dst_window.load_fits's rule verbatim, and
    # `check_agrees_with_dst_window` below asserts this reproduces it row for row on the
    # two frictions that function covers. It is copied rather than imported because
    # load_fits hardcodes mu_p in {0, 0.3} and reading D needs mu_p = 1.
    before = len(d)
    d = d.loc[sorted(d.groupby(["family", "cv_nom", "mu_p", "seed"])["max_iter"].idxmax())]
    print(f"   ladder(mu_p = {mu}): {len(d)} rows after the cap rule "
          f"(dropped {before - len(d)}); caps {sorted(d.max_iter.unique())}")
    g = d.groupby("cv_nom").agg(phi_J=("phi_J_lin", "mean"),
                                sd=("phi_J_lin", "std"),
                                c=("slope_lin", "mean"),
                                c_sd=("slope_lin", "std"),
                                n_ladder=("n_ladder", "min"),
                                n=("phi_J_lin", "size"))
    g["se"] = g["sd"] / np.sqrt(g["n"])
    g["c_se"] = g["c_sd"] / np.sqrt(g["n"])
    return g


def check_agrees_with_dst_window() -> None:
    """Reproduce dst_window.load_fits row for row before any new number is believed.

    Reading the csv directly is what this script did first, and it silently pooled the
    `capfix` re-runs with the originals -- the project's own standing hazard, "the
    HIGHEST cap wins, cell by cell". That gave phi_0(CV 0) = 0.632763 on 6 rows where
    the correct value is 0.632371 on 3, and it moved phi_0(CV 0.5) by 0.0018, a third
    of the effect this script is about.
    """
    from dst_window import load_fits
    ref = load_fits()
    for mu in (0.0, 0.3):
        r = ref[(ref.n_requested == N_MATCH) & (ref.family == FAMILY) & (ref.mu_p == mu)]
        mine = ladder(mu)
        g = r.groupby("cv_nom").phi_J_lin.mean()
        for cv in g.index:
            assert abs(g[cv] - mine.phi_J[cv]) < 1e-12, (
                f"mu {mu} cv {cv}: {g[cv]} vs {mine.phi_J[cv]} -- this script's record "
                f"selection has drifted from dst_window's")
    print("   record selection reproduces dst_window.load_fits exactly at mu_p 0.0/0.3")


def wls_slope(x, y, se):
    """Weighted least squares slope and its standard error. n = 3 here, so this has
    ONE residual degree of freedom -- it cannot self-diagnose and the residual is
    printed rather than a goodness-of-fit statistic."""
    x, y, se = np.asarray(x, float), np.asarray(y, float), np.asarray(se, float)
    w = 1.0 / np.maximum(se, 1e-12) ** 2
    xm = np.sum(w * x) / np.sum(w)
    ym = np.sum(w * y) / np.sum(w)
    sxx = np.sum(w * (x - xm) ** 2)
    b = np.sum(w * (x - xm) * (y - ym)) / sxx
    return float(b), float(1.0 / math.sqrt(sxx)), float(ym - b * xm)


def main() -> int:
    lines: list[str] = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("=== window_shearedge.py -- the DST window with a SHEAR-measured lower edge ===")
    say("   Reads out/phase2_fits.csv and the constant-pressure phi_c transcribed in")
    say("   this module's header. No simulation.")
    say()

    check_agrees_with_dst_window()
    fr = ladder(0.0)     # frictionless -> phi_0, already a p -> 0 intercept
    fm = ladder(0.3)     # frictional isotropic -> the OLD lower edge, and c(CV)

    say("=== the three edges at the three CV levels the shear campaign covers ===")
    say("   cv   phi_0 (iso, mu=0)    phi_J (iso, mu=0.3)   phi_c (SHEAR, p*=3e-4)")
    for cv in CVS:
        say(f"  {cv:.1f}   {fr.phi_J[cv]:.6f} +-{fr.se[cv]:.6f}   "
            f"{fm.phi_J[cv]:.6f} +-{fm.se[cv]:.6f}   "
            f"{PHI_C_SHEAR[cv]:.6f} +-{PHI_C_SHEAR_SD[cv] / math.sqrt(PHI_C_SHEAR_N[cv]):.6f}")
    say()
    say("   the shear critical state sits ABOVE the isotropic frictional phi_J by")
    for cv in CVS:
        say(f"      cv {cv}: {PHI_C_SHEAR[cv] - fm.phi_J[cv]:+.4f}")
    say("   which is the reason the shipped lower edge was flagged as NOT safe.")

    # ---------------- the gauge correction, priced ----------------
    say()
    say("=== the pressure gauge, priced ===")
    say("   phi = phi_J + c * p_hat (phase2_phiJ.py), so a fraction measured AT p*")
    say("   carries c*p* above the p -> 0 intercept the frictionless edge is.")
    say("   cv    c(iso, mu=0.3)     c*p*      c(iso, mu=0)")
    for cv in CVS:
        say(f"  {cv:.1f}   {fm.c[cv]:7.3f} +-{fm.c_se[cv]:5.3f}   "
            f"{fm.c[cv] * P_STAR:.5f}   {fr.c[cv]:7.3f} +-{fr.c_se[cv]:5.3f}")
    dc = fm.c[0.0] - fm.c[0.5]
    dc_se = math.hypot(fm.c_se[0.0], fm.c_se[0.5])
    say(f"   c(frictional) FALLS with CV: {fm.c[0.0]:.2f} -> {fm.c[0.5]:.2f}, "
        f"a difference of {dc:.2f} +- {dc_se:.2f} (t = {dc / dc_se:.2f})")
    say("   -> the gauge bias is CV-DEPENDENT and does NOT cancel out of the contrast.")
    say("   ⚠ this is the ISOTROPIC compressibility used as a PROXY for the")
    say("     critical-state one, which has never been measured here.")

    # ---------------- the three readings ----------------
    readings = {}
    say()
    say("=== W = phi_0 - (lower edge), three readings of the same three CV levels ===")

    def report(name, lo, lo_se, note=""):
        W = [fr.phi_J[cv] - lo[cv] for cv in CVS]
        se = [math.hypot(fr.se[cv], lo_se[cv]) for cv in CVS]
        b, b_se, _ = wls_slope(CVS, W, se)
        dW = W[-1] - W[0]
        dW_se = math.hypot(se[-1], se[0])
        readings[name] = (W, se, b, b_se, dW, dW_se)
        say(f"   -- {name} {note}")
        say("      cv      W       +-se")
        for cv, w, s in zip(CVS, W, se):
            say(f"     {cv:.1f}   {w:.4f}   {s:.4f}")
        say(f"      dW/dCV = {b:+.5f} +- {b_se:.5f}  (t = {b / b_se:+.2f}, n = 3, "
            f"1 residual dof)")
        say(f"      W(0.5) - W(0.0) = {dW:+.4f} +- {dW_se:.4f}  "
            f"(t = {dW / dW_se:+.2f}) -> {'RESOLVED' if abs(dW) > 2 * dW_se else 'NOT RESOLVED'}"
            f" at 2 se; that is {100 * dW / W[0]:+.1f} % of W(0)")
        say()

    report("A. isotropic frictional edge (the SHIPPED reading)",
           {cv: fm.phi_J[cv] for cv in CVS}, {cv: fm.se[cv] for cv in CVS},
           "-- dst_window.py's own, restricted to these 3 levels")
    report("B. SHEAR edge, uncorrected (measured at p* = 3e-04)",
           PHI_C_SHEAR,
           {cv: PHI_C_SHEAR_SD[cv] / math.sqrt(PHI_C_SHEAR_N[cv]) for cv in CVS},
           "-- mixes a finite-p edge with a p -> 0 one")
    lo_corr = {cv: PHI_C_SHEAR[cv] - fm.c[cv] * P_STAR for cv in CVS}
    lo_corr_se = {cv: math.hypot(PHI_C_SHEAR_SD[cv] / math.sqrt(PHI_C_SHEAR_N[cv]),
                                 fm.c_se[cv] * P_STAR) for cv in CVS}
    report("C. SHEAR edge, gauge-corrected with the isotropic c(CV) proxy",
           lo_corr, lo_corr_se, "-- both edges now nominally at p -> 0")

    # Reading D: an INDEPENDENT high-friction lower edge that happens to exist at
    # exactly these three CV levels. phi_J saturates in friction (mu_p 1 -> 3 is inside
    # seed scatter, PREREG_DST_WINDOW_2026-09-09 section 1), so the mu_p = 1 isotropic
    # ladder is the closest thing this project has to Wyart-Cates' high-stress phi_mu
    # WITHOUT shearing anything. It was never compared with the shear edge.
    sat = ladder(1.0)
    report("D. isotropic edge at the SATURATED friction mu_p = 1.0",
           {cv: sat.phi_J[cv] for cv in CVS}, {cv: sat.se[cv] for cv in CVS},
           "-- a second high-friction edge, no shear in it")

    # ---------------- Q1-Q3 ----------------
    say("=== Q1  does the SIGN of dW/dCV survive the change of edge? ===")
    signs = {k: ("narrows" if v[2] < 0 else "widens") for k, v in readings.items()}
    for k, s in signs.items():
        say(f"   {k[:1]}: dW/dCV = {readings[k][2]:+.5f} -> the window {s} with CV")
    say(f"   -> the sign {'SURVIVES' if len(set(signs.values())) == 1 else 'REVERSES'} "
        f"between the shipped edge and the shear-measured one.")
    say()
    say("   ⚠ TWO statistics, and they do not agree on the grouping. The WLS slope")
    say("     weights by 1/se^2 and on these three levels that is nearly a fit through")
    say("     CV 0 and 0.2 alone; the endpoint difference uses the axis end to end.")
    say("      reading   dW/dCV (WLS)      dW endpoint      signs agree?")
    for k, v in readings.items():
        say(f"      {k[:1]}        {v[2]:+.5f}         {v[4]:+.4f}          "
            f"{'yes' if v[2] * v[4] > 0 else 'NO'}")
    say("   ⚠ ALL FOUR READINGS SHARE THE UPPER EDGE phi_0, so any grouping among them")
    say("     is a statement about the LOWER edge alone.")
    ends = {k[:1]: v[4] for k, v in readings.items()}
    wls = {k[:1]: v[2] for k, v in readings.items()}
    say(f"     On the ENDPOINT statistic the three high-friction/sheared edges are "
        f"{'ALL POSITIVE' if all(ends[k] > 0 for k in 'BCD') else 'NOT all positive'} "
        f"(B {ends['B']:+.4f}, C {ends['C']:+.4f}, D {ends['D']:+.4f}) against "
        f"A {ends['A']:+.4f}.")
    say(f"     On the WLS slope D is {wls['D']:+.5f}, i.e. the SAME sign as A and "
        f"consistent with zero (t = {wls['D'] / readings['D. isotropic edge at the SATURATED friction mu_p = 1.0'][3]:+.2f}), "
        f"so that grouping does NOT survive the change of statistic.")
    say("     A claim that 'saturated and sheared edges move together' is therefore NOT")
    say("     supported: what is supported is that no edge resolves a width trend, and")
    say("     that the sign of the central value is set by the edge AND the statistic.")

    say()
    say("=== Q2  is any of them RESOLVED against its own error? ===")
    for k, v in readings.items():
        say(f"   {k[:1]}: |dW| = {abs(v[4]):.4f} against 2 se = {2 * v[5]:.4f}  -> "
            f"{'RESOLVED' if abs(v[4]) > 2 * v[5] else 'NOT RESOLVED'}")
    say("   ⚠ reading A's error here is NOT dst_window.py's: that one pairs the two")
    say("     arms per SEED, which cancels the shared realisation variance and is the")
    say("     right error for two ladders sharing seeds. B and C cannot be paired --")
    say("     the two edges come from different campaigns, N and protocol -- so their")
    say("     larger error is a property of the comparison, not a worse analysis.")
    say(f"   ⚠ phi_0(CV 0.5) carries se {fr.se[0.5]:.4f} on a ladder of "
        f"{int(fr.n_ladder[0.5])} rungs (dst_window flags <= 3 as THIN), and it is what "
        f"dominates every error above.")

    say()
    say("=== Q3  the TRAVERSE: dphi_mu/dCV measured in windows ===")
    say("   TWO slopes are printed and they are not the same statistic. The WLS one")
    say("   weights by 1/se^2, and on three levels whose top rung carries 5x the error")
    say("   of the others that is very nearly a two-point fit through CV 0 and 0.2.")
    say("   The ENDPOINT slope (phi_mu(0.5) - phi_mu(0))/0.5 uses the axis under test")
    say("   end to end and is what dst_window.py's 6-level '2.23 windows' compares to.")
    trav = {}
    for name, lo, lo_se in (
            ("A isotropic frictional", {cv: fm.phi_J[cv] for cv in CVS},
             {cv: fm.se[cv] for cv in CVS}),
            ("B shear, uncorrected", PHI_C_SHEAR,
             {cv: PHI_C_SHEAR_SD[cv] / math.sqrt(PHI_C_SHEAR_N[cv]) for cv in CVS}),
            ("C shear, gauge-corrected", lo_corr, lo_corr_se)):
        b, b_se, _ = wls_slope(CVS, [lo[cv] for cv in CVS], [lo_se[cv] for cv in CVS])
        bend = (lo[CVS[-1]] - lo[CVS[0]]) / (CVS[-1] - CVS[0])
        bend_se = math.hypot(lo_se[CVS[-1]], lo_se[CVS[0]]) / (CVS[-1] - CVS[0])
        Wbar = float(np.mean([fr.phi_J[cv] - lo[cv] for cv in CVS]))
        trav[name[0]] = (bend / Wbar, Wbar / bend)
        say(f"   {name}: mean W = {Wbar:.4f}")
        say(f"      WLS      dphi_mu/dCV = {b:+.5f} +- {b_se:.5f} -> "
            f"{b / Wbar:.2f} windows per unit CV (one window at dCV = {Wbar / b:.3f})")
        say(f"      ENDPOINT dphi_mu/dCV = {bend:+.5f} +- {bend_se:.5f} -> "
            f"{bend / Wbar:.2f} windows per unit CV (one window at "
            f"dCV = {Wbar / bend:.3f})")
    tv = [v[0] for v in trav.values()]
    dcv = [v[1] for v in trav.values()]
    say(f"   -> on the endpoint statistic the traverse runs {min(tv):.2f}-{max(tv):.2f} "
        f"windows per unit CV across the three readings (a spread of "
        f"{max(tv) / min(tv):.2f}x), i.e. one full window at "
        f"dCV = {min(dcv):.2f}-{max(dcv):.2f}.")

    say()
    say("=== WHAT THIS DOES AND DOES NOT SETTLE ===")
    say("   SETTLED: the SIGN of the width trend is a property of which lower edge")
    say("     is used. Central values: A narrows, B and C widen. And on the three CV")
    say("     levels the shear campaign covers, NONE of the three readings resolves a")
    say("     change at 2 se -- including A. dst_window.py's 'RESOLVED as a CHANGE'")
    say("     rests on SIX levels and a PER-SEED pairing of two ladders that share")
    say("     their seeds; the shear edge comes from a different campaign, N and")
    say("     protocol, so it can have neither. The honest statement is that the width")
    say("     trend is not established, not that the shipped one is refuted.")
    say(f"   SETTLED: the traverse survives all three readings at "
        f"{min(tv):.2f}-{max(tv):.2f} windows per unit CV (endpoint statistic), so")
    say(f"     'a CV ladder of about {np.mean(dcv):.1f} moves a material across the")
    say("     window' does not depend on the edge definition.")
    say("   NOT SETTLED: the critical-state compressibility. Reading C uses the")
    say("     isotropic proxy. The 4-pressure constant-pressure ladder measures it")
    say("     directly and would also give phi_c(p -> 0) with no proxy at all.")
    say("   NOT SETTLED: phi_0 itself at CV 0.5, whose THIN 3-rung ladder dominates")
    say("     every error bar above. A frictionless constant-pressure ladder would")
    say("     give phi_0 the same way phi_c now comes.")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "window_shearedge.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
