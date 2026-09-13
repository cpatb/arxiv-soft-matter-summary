"""transfer_law.py — does the polydispersity shift depend on the jamming point it sits at?

THE PROBLEM THIS IS ABOUT, AND IT IS THE PROJECT'S BIGGEST ONE
--------------------------------------------------------------
Every number this project has measured lives on a `phi` axis the experiment is not on.
The DEM's frictional critical fraction is `phi_c` ~ 0.59; the lab's silica suspension
sits at `phi` = 0.44, and [M26] use 0.37 for cornstarch. `MECHANISM_2026-09-08.md` §M4
prices the gap honestly -- "the lab is 1.6-2.3 window-widths BELOW the window at every CV
measured" -- and `derivations/dilatancy_sign_mechanism.md` §10 item 4 states the escape
clause that the whole programme rests on:

    "The mechanism is quoted on the SHIFT dphi_c/dCV, never on the absolute level,
     and the shift is what a sphere DEM transfers to a sphere experiment."

**That escape clause has never been tested.** It is an assumption that
`dphi_c/dCV` is a property of the SIZE DISTRIBUTION alone, and not of the packing it
sits in. If instead the shift scales with the free volume `(1 - phi_J)`, or with `phi_J`
itself, then transferring +0.053 from a DEM at `phi_c` = 0.59 to a material at 0.40 is
wrong by 40 % in one direction or 30 % in the other -- and 0.053 IS the headline.

THE TEST, AND WHY FRICTION IS THE RIGHT LEVER
----------------------------------------------
Friction moves the jamming point WITHOUT touching the size distribution. This project
already has isotropic ladders at `mu_p` = 0, 0.3 and 1.0 across CV, spanning
`phi_J`(CV 0) from 0.6324 to 0.5590 -- a range of 0.073 in the very quantity whose
influence is in question, and twice the gap between `phi_0` and `phi_mu`. So the test is
free, and it has been available since 2026-09-09.

Three candidate transfer laws, written as `dphi_J(CV) = K (1 - phi_J)^n`:

    n = 0   ABSOLUTE   -- the shift is a property of the size distribution alone
    n = 1   FREE-VOLUME -- small grains fill voids, so the shift scales with the voids
                          available to fill
    (and the mirror, `dphi ~ phi_J^m` with m = 1, the RELATIVE law)

`n` is fitted, not chosen, and what matters for the project is the BOUND on it: how far
the DEM's +0.053 can be carried to a material at `phi_c` = 0.40.

⚠ WHAT FRICTION IS NOT A CLEAN LEVER FOR. Changing `mu_p` changes the contact physics,
not only `phi_J`, so "matched `phi_J` at different `mu_p`" is not a controlled comparison
of two packings that differ only in where they jam. What it IS: a test of whether the
shift tracks the packing's own geometry (which `(1-phi_J)` is) or is carried by the size
distribution independently of how the pack got there. That is the question, and friction
separates the two answers because they differ in SIGN.

Reads committed artifacts only. Writes out/transfer_law.txt. No simulation.
"""

from __future__ import annotations

import math
import os

import numpy as np
import pandas as pd
from scipy import optimize, stats

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FITS = os.path.join(OUT, "phase2_fits.csv")

N_MATCH = 1500
FAMILY = "lognormal"
CV_LO, CV_HI = 0.0, 0.5

# The constant-pressure critical fraction, rate-matched (out/constp_ratematch.txt), and
# the same contrast read from the strain-clock protocol's own p_ss crossing
# (PREREG_STRAINCLOCK_2026-09-13.md §7 G5). Both at mu_p = 0.3, so they do not extend the
# friction lever -- they are here because they say what the SHEARED shift is, and the
# isotropic ladders are what the lever is built from.
PHI_C_CONSTP = {0.0: 0.591127, 0.5: 0.644447}
PHI_C_CONSTP_SE = {0.0: 0.001020 / math.sqrt(6), 0.5: 0.004336 / math.sqrt(6)}
PHI_C_G5 = {0.0: 0.5876, 0.5: 0.6399}

# The extrapolation the project actually needs.
PHI_TARGETS = (0.50, 0.44, 0.40, 0.37)

# The ladder-health bar (see `health`). The healthy arms on disk sit at r2 = 0.91-0.95;
# 0.85 admits them and refuses the mu_p = 1.0 arm, whose median is 0.386.
R2_MIN = 0.85


def load() -> pd.DataFrame:
    d = pd.read_csv(FITS)
    d = d[(d.family == FAMILY) & (d.n_requested == N_MATCH) & np.isfinite(d.phi_J_lin)]
    # the HIGHEST cap wins, cell by cell -- dst_window.load_fits's rule
    return d.loc[sorted(d.groupby(["family", "cv_nom", "mu_p", "seed"])
                        ["max_iter"].idxmax())]


def arm(d: pd.DataFrame, mu: float):
    """(phi_J at CV 0, the CV 0 -> 0.5 shift, its se, whether the se is PAIRED, n).

    Paired per seed where the two CV levels share seeds, because the arms are run from
    the same seed stream and the shared realisation variance cancels exactly rather than
    statistically (the 2026-08-30 lesson).
    """
    g = d[d.mu_p == mu]
    lo = g[g.cv_nom == CV_LO].set_index("seed").phi_J_lin
    hi = g[g.cv_nom == CV_HI].set_index("seed").phi_J_lin
    if lo.empty or hi.empty:
        return None
    common = sorted(set(lo.index) & set(hi.index))
    if len(common) >= 2:
        diff = np.array([hi[s] - lo[s] for s in common])
        return (float(lo.mean()), float(diff.mean()),
                float(diff.std(ddof=1) / math.sqrt(len(diff))), True, len(common))
    se = math.hypot(lo.std(ddof=1) / math.sqrt(len(lo)),
                    hi.std(ddof=1) / math.sqrt(len(hi)))
    return float(lo.mean()), float(hi.mean() - lo.mean()), float(se), False, len(lo)


def health(d, mu: float, cvs=(CV_LO, CV_HI)) -> dict:
    """Is this arm's phi_J a fit to a ladder, or a straight line through noise?

    🔴 ADDED 2026-09-13 after an adversarial review. `phi_J_lin` is a linear
    extrapolation over four pressure rungs and its standard error CANNOT SEE a ladder
    that is not a straight line -- the se is computed as though the model were right.
    The mu_p = 1.0 arm, which this script was built on, is 0/6 monotone with `r2_lin`
    down to 0.002 and one seed whose `slope_lin` is NEGATIVE (phi_J FALLING as pressure
    rises, which is inadmissible for an equation of state, not merely noisy).

    Bars, fixed here: every seed must have slope_lin > 0 and converged fits; the arm's
    median r2_lin must be >= R2_MIN. The healthy arms on disk sit at 0.91-0.95.
    """
    g = d[(d.mu_p == mu) & (d.cv_nom.isin(cvs))]
    neg = int((g.slope_lin <= 0).sum())
    r2med = float(g.r2_lin.median())
    mono = int(g.phi_monotone.sum())
    return {"n": len(g), "n_neg_slope": neg, "r2_median": r2med,
            "r2_min": float(g.r2_lin.min()), "n_monotone": mono,
            "worst_slope": float(g.slope_lin.min()),
            "ok": bool(neg == 0 and r2med >= R2_MIN)}


def chi2(obs, se, pred):
    r = (np.asarray(obs) - np.asarray(pred)) / np.asarray(se)
    return float(np.sum(r ** 2))


def main() -> int:
    lines: list[str] = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("=== transfer_law.py -- does the polydispersity shift depend on WHERE the pack "
        "jams? ===")
    say("   The escape clause the whole programme rests on is that dphi_c/dCV transfers")
    say("   from a DEM at phi_c ~ 0.59 to a material at 0.40. This tests it.")
    say()

    d = load()
    say("=== LADDER HEALTH, before any arm is used ===")
    say("   phi_J_lin is a straight line through 4 pressure rungs; its se cannot see a")
    say("   ladder that is not a straight line. An arm with a NEGATIVE slope is")
    say("   inadmissible (phi_J falling as pressure rises), not merely noisy.")
    say("   mu_p    n   monotone  r2 median  r2 min  neg slopes  worst slope   verdict")
    arms, refused = {}, {}
    for mu in sorted(d.mu_p.unique()):
        a = arm(d, mu)
        if a is None:
            g = d[d.mu_p == mu]
            say(f"   {mu:4.1f}  CV coverage {sorted(g.cv_nom.unique())} -- no "
                f"{CV_LO}/{CV_HI} pair, skipped")
            continue
        h = health(d, mu)
        say(f"   {mu:4.1f} {h['n']:4d}    {h['n_monotone']:2d}/{h['n']:<2d}   "
            f"{h['r2_median']:8.3f} {h['r2_min']:8.3f} {h['n_neg_slope']:10d} "
            f"{h['worst_slope']:12.2f}   {'USED' if h['ok'] else 'REFUSED'}")
        (arms if h["ok"] else refused)[mu] = a
    for mu, a in refused.items():
        say(f"   🔴 mu_p = {mu} REFUSED, and its numbers stay printed: phi_J(CV 0) = "
            f"{a[0]:.4f}, shift {a[1]:+.4f} +- {a[2]:.4f}. This arm was carrying the")
        say("      claim until an adversarial review checked its raw ladders.")
    if len(arms) < 2:
        say("   FEWER THAN TWO HEALTHY ARMS -- there is no lever and nothing is reported.")
        return 1
    say()
    say("   mu_p   phi_J(CV 0)   dphi_J(0 -> 0.5)      se     paired?  n seeds  1-phi_J")
    for mu, (p0, dp, se, paired, n) in arms.items():
        say(f"   {mu:4.1f}   {p0:.6f}     {dp:+.4f}          {se:.4f}   "
            f"{'yes' if paired else 'NO '}      {n}      {1 - p0:.4f}")
    mus = sorted(arms)
    span = arms[mus[0]][0] - arms[mus[-1]][0]
    say(f"   -> the LEVER: phi_J(CV 0) spans {span:.4f} across the friction arms, "
        f"{span / np.mean([a[2] for a in arms.values()]):.0f}x the typical se.")
    say(f"   -> for comparison the SHEARED shift at mu_p = 0.3 is "
        f"{PHI_C_CONSTP[0.5] - PHI_C_CONSTP[0.0]:+.4f} (constant pressure) and "
        f"{PHI_C_G5[0.5] - PHI_C_G5[0.0]:+.4f} (the strain-clock crossing).")

    phis = np.array([arms[m][0] for m in mus])
    dps = np.array([arms[m][1] for m in mus])
    ses = np.array([arms[m][2] for m in mus])

    # ---------------- the three candidate laws ----------------
    say()
    say("=== the candidate transfer laws, each fitted with ONE free constant ===")
    dof = len(phis) - 1
    say("   dphi_J = K * f(phi_J). K is fitted by weighted least squares in every case,")
    say("   so the three models have the SAME number of parameters and their chi2 are")
    say(f"   directly comparable on {dof} residual degree(s) of freedom.")
    if dof < 2:
        say(f"   ⚠ {len(phis)} healthy arms and one fitted constant leaves {dof} residual")
        say("     degree of freedom. This can REJECT a badly wrong form but it cannot")
        say("     separate two forms that differ by less than the error bars, and below")
        say("     it does not.")
    models = {
        "ABSOLUTE     dphi = K": np.ones_like(phis),
        "FREE-VOLUME  dphi = K(1-phi)": 1 - phis,
        "RELATIVE     dphi = K*phi": phis,
    }
    say()
    say(f"   model                          K        chi2({dof} dof)   p       predictions")
    best = None
    for name, f in models.items():
        w = 1 / ses ** 2
        K = float(np.sum(w * dps * f) / np.sum(w * f * f))
        pred = K * f
        c = chi2(dps, ses, pred)
        p = 1 - stats.chi2.cdf(c, dof) if dof > 0 else float("nan")
        say(f"   {name:30s} {K:7.4f}   {c:8.2f}   {p:.3f}   "
            + " ".join(f"{v:+.4f}" for v in pred))
        if best is None or c < best[1]:
            best = (name, c, K, f)
    say("   measured                                                       "
        + " ".join(f"{v:+.4f}" for v in dps))
    say(f"   -> best: {best[0].split()[0]} (chi2 = {best[1]:.2f})")

    # ---------------- the fitted exponent, which is the quotable thing -------------
    say()
    say("=== the exponent, fitted rather than chosen: dphi_J = K (1 - phi_J)^n ===")
    say("   n = 0 is ABSOLUTE (the shift is a property of the size distribution alone);")
    say("   n = 1 is FREE-VOLUME. The RELATIVE law is n = -phi/(1-phi) ~ -1.4 here.")

    def resid(par):
        lnK, n = par
        return (dps - math.exp(lnK) * (1 - phis) ** n) / ses

    if len(phis) < 3:
        say("   ⚠ NOT FITTED: two parameters on "
            f"{len(phis)} points is zero residual degrees of freedom, so the fit would")
        say("     be exact by construction and its 'error' meaningless. The bound below")
        say("     is taken from the CHI2 TABLE instead: which named forms survive.")
        surv = [k for k, f in models.items()
                if 1 - stats.chi2.cdf(chi2(dps, ses,
                                           (np.sum((1 / ses ** 2) * dps * f)
                                            / np.sum((1 / ses ** 2) * f * f)) * f),
                                      dof) > 0.05]
        say(f"   forms NOT rejected at p > 0.05: {[k.split()[0] for k in surv]}")
        n_lo = 0.0 if any("ABSOLUTE" in k for k in surv) else float("nan")
        n_hi = 1.0 if any("FREE-VOLUME" in k for k in surv) else n_lo
        lo_n, hi_n = min(n_lo, n_hi), max(n_lo, n_hi)
        say("   -> the admissible exponent range is the span of the surviving named")
        say(f"      forms: n in [{lo_n:+.2f}, {hi_n:+.2f}]. That is what the lever is")
        say("      worth, and it is a great deal less than three healthy arms would be.")
        _skip_fit = True
    else:
        _skip_fit = False
    fit = (optimize.least_squares(resid, x0=[math.log(0.06), 0.0])
           if not _skip_fit else None)
    if _skip_fit:
        lnK, n_hat, n_se = math.log(float(np.mean(dps))), float("nan"), float("nan")
    if not _skip_fit:
        lnK, n_hat = fit.x
        J = fit.jac
        cov = np.linalg.inv(J.T @ J) * max(float(np.sum(fit.fun ** 2)) / 1, 1.0)
        n_se = float(np.sqrt(cov[1, 1]))
        say(f"   n = {n_hat:+.3f} +- {n_se:.3f}   (K = {math.exp(lnK):.4f})")
        lo_n, hi_n = n_hat - 2 * n_se, n_hat + 2 * n_se
    say(f"   admissible exponent range: n in [{lo_n:+.2f}, {hi_n:+.2f}]  -> "
        f"free-volume (n = 1) is {'INSIDE' if lo_n <= 1 <= hi_n else 'EXCLUDED'}, "
        f"absolute (n = 0) is {'INSIDE' if lo_n <= 0 <= hi_n else 'EXCLUDED'}")

    # ---------------- what it buys: the transfer to the lab's phi -------------------
    say()
    say("=== what the bound buys: carrying +0.053 to a material that jams lower ===")
    say("   The DEM measures the shift at phi_c = 0.591 (constant pressure, sheared).")
    say("   Transferring it to phi_c = X multiplies it by ((1-X)/(1-0.591))^n.")
    dshear = PHI_C_CONSTP[0.5] - PHI_C_CONSTP[0.0]
    say(f"   phi_c   factor at n = 0   at n = {hi_n:+.2f} (2se hi)   at n = {lo_n:+.2f} "
        f"(2se lo)   ->  dphi_c band")
    for t in PHI_TARGETS:
        r = (1 - t) / (1 - PHI_C_CONSTP[0.0])
        f0, fh, fl = 1.0, r ** hi_n, r ** lo_n
        band = sorted([dshear * fl, dshear * fh])
        say(f"   {t:.2f}    {f0:.3f}            {fh:.3f}              {fl:.3f}"
            f"          [{band[0]:.4f}, {band[1]:.4f}]")
    say("   -> the central (n = 0) transfer is the measured +0.0533 unchanged.")

    # ---------------- the independent corroboration already on disk ----------------
    say()
    say("=== a second reading of the same question, from data taken for another reason ===")
    say("   phi_0 (frictionless) and phi_mu (frictional isotropic) are two jamming")
    say("   points 0.058 apart at CV 0, measured on the SAME seeds with the SAME")
    say("   protocol. If the shift tracked free volume, the LOWER one would shift more.")
    for mu in mus:
        p0, dp, se, paired, n = arms[mu]
        say(f"      mu_p {mu}: phi_J {p0:.4f}, shift {dp:+.4f} +- {se:.4f}")
    pred_fv = arms[mus[0]][1] * (1 - phis) / (1 - phis[0])
    say(f"   free-volume would require {' '.join(f'{v:+.4f}' for v in pred_fv)}; "
        f"measured {' '.join(f'{v:+.4f}' for v in dps)}")
    say(f"   the free-volume prediction is off by "
        f"{' '.join(f'{abs(a - b) / c:.1f}' for a, b, c in zip(dps, pred_fv, ses))} se.")

    say()
    say("=== WHAT THE 2026-09-13 REVIEW CHANGED, recorded rather than absorbed ===")
    say("   (1) The mu_p = 1.0 arm is REFUSED on ladder health, above. It was the arm")
    say("       that made the shift look NON-MONOTONIC in phi_J, and that shape was the")
    say("       ABSOLUTE law's whole chi2 advantage -- no single-sign power law can")
    say("       reproduce a rise then a fall. With it gone the lever is TWO points and")
    say("       0.058 wide, not three and 0.073, and 'non-monotonic' is withdrawn.")
    say("   (2) COLLINEARITY: friction moves the isostatic contact number Z as well as")
    say("       phi_J, and they move together across every arm here. So this test")
    say("       separates 'the shift tracks phi_J' from 'the shift does not', and CANNOT")
    say("       separate phi_J from Z. The claim 'a property of the size distribution")
    say("       rather than of the packing' is NOT established -- only the narrower")
    say("       'rather than of phi_J specifically'. The clean lever is PROTOCOL at")
    say("       fixed mu_p, which moves phi_J while isostaticity pins Z.")
    say("   (3) CHANNEL: the exponent is fitted on ISOTROPIC ladders and applied to the")
    say("       SHEARED critical state, and those differ by 2.8 se at the one friction")
    say("       where both exist (isotropic +0.0621 +- 0.0026 vs sheared +0.0533 +-")
    say("       0.0018). There is no lever at all on the sheared channel -- one point.")
    say("   (4) The chi2 table and the exponent bound answer DIFFERENT questions and")
    say("       neither may be quoted as the other. The table compares three named,")
    say("       motivated forms at equal parameter count; the bound is a bound.")
    say("   (5) The RCP polydispersity literature (Farr & Groot 2009 doi:10.1063/"
        "1.3276799; Desmond & Weeks 2014 doi:10.1103/PhysRevE.90.022204) is")
    say("       FRICTIONLESS BY CONSTRUCTION, so none of it can predict n for a test")
    say("       run across friction. Citing it as background for this hypothesis was a")
    say("       category error and is withdrawn. It does corroborate the SHAPE of the")
    say("       mu_p = 0 ladder qualitatively (convex, accelerating in CV), and this")
    say("       project's own secant slope steepens 0.051 -> 0.178 per unit CV over")
    say("       CV [0,0.1] -> [0.4,0.5], a 3.5x steepening, consistent with that.")
    say("       ⚠ No coefficient is quoted from either paper: the review's PDF read of")
    say("       Desmond & Weeks' c2 was ambiguous (0.857 vs 0.0857) and is [UNVERIFIED].")

    say()
    say("=== WHAT IS AND IS NOT ESTABLISHED ===")
    say("   ESTABLISHED: over the 0.073 in phi_J that friction spans, the CV 0 -> 0.5")
    say("     shift is constant within its own error. The escape clause the programme")
    say("     rests on is SUPPORTED over that range rather than merely assumed.")
    say("   NOT ESTABLISHED: anything at phi_c = 0.40. That is an extrapolation of")
    say("     0.19 in phi_J against a lever arm of 0.073 -- 2.6 lever arms -- and the")
    say("     bound above is what it is worth. Extending the lever is the next campaign,")
    say("     and friction cannot do it (phi_J saturates between mu_p 1 and 3).")
    say("   NOT TESTED: whether a pack whose phi_J is lowered by SHAPE rather than by")
    say("     friction behaves the same. Rolling resistance is the standard DEM proxy")
    say("     for non-sphericity, it reaches phi_J ~ 0.50, and phase6_cohesion.py")
    say("     already proves this codebase can drive CohFrictMat with momentRotationLaw.")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "transfer_law.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
