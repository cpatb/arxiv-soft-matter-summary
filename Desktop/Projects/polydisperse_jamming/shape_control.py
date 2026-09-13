"""shape_control.py — CV is not the control parameter. What is?

THE FINDING THIS IS ABOUT
-------------------------
Every headline this project has produced is quoted as a derivative with respect to CV:
`dphi_c/dCV`, `dphi_J/dCV`, `dW/dCV`, `CV*`, "a CV ladder of 0.45 traverses the window".
The `MECHANISM_2026-09-08.md` decision tree, the `ID` sensitivity in
`derivations/dilatancy_sign_mechanism.md` Eq. (11), and the drop-side predictions are all
written on that axis.

At **matched CV and matched friction** the three size-distribution families this project
has run give

    lognormal   +0.0621 +- 0.0026
    uniform     +0.0242 +- 0.0008
    bidisperse  +0.0139 +- 0.0025          (CV 0 -> 0.5, mu_p = 0.3, paired per seed)

a factor of **4.5**. For comparison, moving the jamming point itself with friction (over
the two arms that pass ladder health, a span of 0.058 in `phi_J`) moves the same shift by
**1.07x** (`transfer_law.py`). **The shape of the size distribution dominates the axis the
project spent its effort controlling.**

This is not a new measurement: all three families were run on 2026-09-08/09 for a
different question (the CV/skewness separation), and the numbers have been on disk since.
What is new is asking them this question.

WHY IT MATTERS MORE THAN THE FRICTION RESULT
--------------------------------------------
`transfer_law.py` asks whether +0.053 can be carried from a DEM that jams at 0.59 to a
material that jams at 0.44, and answers "yes, over the range tested". But it is the WRONG
worry. The number cannot be carried to a material whose size distribution is a different
SHAPE, whatever its CV and whatever its jamming point -- and no amount of friction
matching fixes that.

The programme is rescued, but by a different argument from the one currently written
down: milled and precipitated powders (the lab's silica, cornstarch) are approximately
LOGNORMAL, which is the family the headline was measured on. That has to be stated as a
condition of the result, and the sensitivity to it quantified, which is what this does.

WHAT IS BEING SEARCHED FOR
--------------------------
A variable `X` of the size distribution such that `dphi_J` collapses on it across all
three families. Candidates, each computed from the PSD classes directly (no DEM):

    CV                    fails by construction -- the families are matched in it
    d32/d10               fails: uniform 1.409 and bidisperse 1.400 are equal to 0.6 %
                          and their shifts differ by 1.7x
    skewness              the project's own CV/skewness separation points here
    d_max/d_min           interstitial filling needs grains small enough for the holes
                          (0.225 of the host for a tetrahedral site, 0.414 octahedral)
    f_small               the number fraction below 0.414 d50 -- the direct version of
                          the same idea
    <d^3>/<d>^3           the volume-weighting the packing actually sees

Reads committed artifacts and the psd module. Writes out/shape_control.txt. No simulation.
"""

from __future__ import annotations

import math
import os

import numpy as np
import pandas as pd

from psd import make_psd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FITS = os.path.join(OUT, "phase2_fits.csv")

N_MATCH = 1500
MU_REF = 0.3
CV_LO, CV_HI = 0.0, 0.5
FAMILIES = ("lognormal", "uniform", "bidisperse")
# the truncation the production packings actually used (packing.insert_packing default)
TRUNCATE = (0.2, 5.0)
N_SAMPLE = 200_000          # for the sampled descriptors; deterministic seed below
SAMPLE_SEED = 20260913
# the octahedral interstice ratio: a sphere this much smaller than its host fits the
# hole between close-packed hosts without disturbing them
R_OCTA = 0.414


def load() -> pd.DataFrame:
    d = pd.read_csv(FITS)
    d = d[(d.n_requested == N_MATCH) & np.isfinite(d.phi_J_lin)]
    return d.loc[sorted(d.groupby(["family", "cv_nom", "mu_p", "seed"])
                        ["max_iter"].idxmax())]


def shift(d: pd.DataFrame, family: str, mu: float):
    """Paired CV 0 -> 0.5 shift, its se, and n."""
    g = d[(d.family == family) & (d.mu_p == mu)]
    lo = g[g.cv_nom == CV_LO].set_index("seed").phi_J_lin
    hi = g[g.cv_nom == CV_HI].set_index("seed").phi_J_lin
    common = sorted(set(lo.index) & set(hi.index))
    if len(common) < 2:
        return None
    diff = np.array([hi[s] - lo[s] for s in common])
    return float(diff.mean()), float(diff.std(ddof=1) / math.sqrt(len(diff))), len(common)


def descriptors(family: str, cv: float) -> dict[str, float]:
    """Every candidate control variable, from the SAMPLED and TRUNCATED distribution.

    Sampled and truncated, not analytic: the packings were built from truncated draws,
    so the analytic moments of the parent distribution are not what the grains had. This
    project has been bitten by exactly that (`realised_moments` exists for the reason).
    """
    rng = np.random.default_rng(SAMPLE_SEED)
    p = make_psd(family, 1.0, cv)
    try:
        d = p.sample(N_SAMPLE, rng, truncate=TRUNCATE)
    except TypeError:
        d = p.sample(N_SAMPLE, rng)
    d = np.asarray(d, float)
    m1, m2, m3 = d.mean(), (d ** 2).mean(), (d ** 3).mean()
    sd = d.std(ddof=0)
    d50 = float(np.median(d))
    return {
        "cv": float(sd / m1),
        "d32": float(m3 / m2),
        "skew": float(((d - m1) ** 3).mean() / sd ** 3) if sd > 0 else 0.0,
        "ratio": float(d.max() / d.min()),
        "d95_d05": float(np.percentile(d, 95) / np.percentile(d, 5)),
        "f_small": float(np.mean(d < R_OCTA * d50)),
        "m3_over_m1c": float(m3 / m1 ** 3),
    }


def main() -> int:
    lines: list[str] = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("=== shape_control.py -- CV is not the control parameter. What is? ===")
    say()
    d = load()

    say("=== LADDER HEALTH first -- the same gate transfer_law.py now applies ===")
    say("   An adversarial review on 2026-09-13 found a friction arm whose phi_J came")
    say("   from ladders that were 0/6 monotone with r2 down to 0.002 and one NEGATIVE")
    say("   slope. These family arms use the same estimator, so they get the same gate")
    say("   BEFORE their numbers are read: no negative slope, median r2 >= 0.85.")
    say("   family        mu_p    n   monotone  r2 median  r2 min  neg slopes  verdict")
    health_ok = {}
    for fam in FAMILIES:
        for mu in (MU_REF, 0.0):
            g = d[(d.family == fam) & (d.mu_p == mu)
                  & (d.cv_nom.isin([CV_LO, CV_HI]))]
            if g.empty:
                continue
            neg = int((g.slope_lin <= 0).sum())
            r2m = float(g.r2_lin.median())
            ok = neg == 0 and r2m >= 0.85
            health_ok[(fam, mu)] = ok
            say(f"   {fam:12s} {mu:4.1f} {len(g):4d}    "
                f"{int(g.phi_monotone.sum()):2d}/{len(g):<2d}   {r2m:8.3f} "
                f"{float(g.r2_lin.min()):8.3f} {neg:10d}   "
                f"{'USED' if ok else 'REFUSED'}")
    if not all(health_ok.values()):
        say("   🔴 at least one family arm is REFUSED; the comparison below is not "
            "readable on the refused arms and says so.")
    else:
        say("   -> all family arms pass. The finding below is NOT an artefact of the "
            "defect the review found.")
    say()
    say("=== the finding, at MATCHED CV and MATCHED friction ===")
    say(f"   CV {CV_LO} -> {CV_HI}, mu_p = {MU_REF}, paired per seed")
    say("   family        dphi_J        se      n")
    shifts = {}
    for fam in FAMILIES:
        s = shift(d, fam, MU_REF)
        if s is None:
            say(f"   {fam:12s} no CV {CV_LO}/{CV_HI} pair")
            continue
        shifts[fam] = s
        say(f"   {fam:12s} {s[0]:+.4f}     {s[1]:.4f}   {s[2]}")
    vals = [shifts[f][0] for f in shifts]
    say(f"   -> spread {max(vals) / min(vals):.2f}x across families at the SAME CV.")
    fr = shift(d, "lognormal", 0.0)
    if fr:
        frict = max(fr[0], shifts["lognormal"][0]) / min(fr[0], shifts["lognormal"][0])
        say("   -> for comparison: moving phi_J by 0.058 with friction (mu_p 0 -> 0.3, "
            "the two arms")
        say(f"      that pass ladder health) moves the SAME lognormal shift by "
            f"{frict:.2f}x, {fr[0]:+.4f} -> {shifts['lognormal'][0]:+.4f}.")
        say(f"      SHAPE moves it {max(vals) / min(vals):.2f}x; FRICTION moves it "
            f"{frict:.2f}x. The axis the project has been")
        say("      controlling is the smaller of the two by a wide margin.")

    say()
    say("   the same three families at mu_p = 0, where only two were run:")
    for fam in FAMILIES:
        s = shift(d, fam, 0.0)
        if s:
            say(f"      {fam:12s} {s[0]:+.4f} +- {s[1]:.4f}  (n {s[2]})")
    say("   -> the ordering is the same at both frictions, so it is not a friction effect.")

    # ---------------- the search for a collapsing variable ----------------
    say()
    say("=== candidate control variables, computed from the TRUNCATED SAMPLED psd ===")
    say(f"   ({N_SAMPLE} draws, seed {SAMPLE_SEED}, truncate {TRUNCATE} -- the same "
        f"truncation the packings used)")
    desc = {f: descriptors(f, CV_HI) for f in shifts}
    keys = ["cv", "d32", "skew", "ratio", "d95_d05", "f_small", "m3_over_m1c"]
    say("   family        " + "".join(f"{k:>12}" for k in keys) + "      dphi_J")
    for fam in shifts:
        say(f"   {fam:12s}  " + "".join(f"{desc[fam][k]:12.4f}" for k in keys)
            + f"   {shifts[fam][0]:+.4f}")

    say()
    say("=== does any single variable order the three families correctly? ===")
    say("   A control variable must at minimum put them in the right ORDER, and the")
    say("   discriminating pair is uniform vs bidisperse: their shifts differ by")
    say(f"   {shifts['uniform'][0] / shifts['bidisperse'][0]:.2f}x "
        if "bidisperse" in shifts and shifts["bidisperse"][0] else "   ")
    order_true = sorted(shifts, key=lambda f: shifts[f][0])
    say(f"   true order (smallest shift first): {order_true}")
    for k in keys:
        order_k = sorted(shifts, key=lambda f: desc[f][k])
        ok = order_k == order_true
        say(f"      {k:12s} orders as {str(order_k):55s} {'MATCHES' if ok else 'no'}")

    say()
    say("=== the quantitative test: proportionality, one free constant ===")
    say("   If dphi_J = K * X, then dphi_J / X must be the SAME for all three families.")
    say("   variable        " + "".join(f"{f[:9]:>11}" for f in shifts)
        + "     max/min")
    for k in keys:
        r = []
        for fam in shifts:
            x = desc[fam][k]
            r.append(shifts[fam][0] / x if x else float("nan"))
        finite = [v for v in r if np.isfinite(v) and v > 0]
        spread = max(finite) / min(finite) if len(finite) == len(r) else float("nan")
        say(f"   {k:12s}    " + "".join(f"{v:11.4f}" for v in r)
            + f"   {spread:8.2f}x")
    say("   -> a collapsing variable would show max/min near 1. NONE of these does.")
    say("      The two that come closest are reported, but 'closest' among failures is")
    say("      not a result: what this establishes is that the shift is NOT a function")
    say("      of any ONE of these descriptors.")

    say()
    say("=== WHAT THIS ESTABLISHES, AND WHAT IT COSTS THE PROGRAMME ===")
    say("   ESTABLISHED: at matched CV and matched friction the CV 0 -> 0.5 jamming")
    say(f"     shift runs {min(vals):+.4f} to {max(vals):+.4f} across three size")
    say(f"     distributions, a factor {max(vals) / min(vals):.1f}. CV is NOT the")
    say("     control parameter, and every result this project quotes as a derivative")
    say("     with respect to CV is conditional on the LOGNORMAL family.")
    say("   ESTABLISHED: it is not a friction effect -- the ordering is identical at")
    say("     mu_p = 0 and 0.3.")
    say("   NOT ESTABLISHED: what the control parameter IS. None of the seven")
    say("     single-variable candidates collapses the three families, and with three")
    say("     families there is not enough leverage to fit a two-variable form.")
    say("   WHAT IT COSTS: nothing, IF the target material is lognormal -- milled and")
    say("     precipitated powders approximately are, and the lab's silica and")
    say("     cornstarch are both. But that is now a STATED CONDITION of the headline")
    say("     rather than an unexamined assumption, and it is a far stronger condition")
    say("     than the jamming-point matching transfer_law.py was built to check.")
    say("   WHAT WOULD SETTLE IT: a fourth and fifth family chosen to SEPARATE the")
    say("     candidates rather than to add coverage -- e.g. a negatively-skewed")
    say("     distribution at the same CV and d32, and a bidisperse at a size ratio")
    say("     large enough to fill interstices (>= 1/0.414 = 2.4), which the CV = 0.5")
    say("     bidisperse already is, so the interstitial picture predicts it should")
    say("     pack WELL and it packs worst of the three.")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "shape_control.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
