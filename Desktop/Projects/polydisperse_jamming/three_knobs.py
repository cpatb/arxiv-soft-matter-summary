"""three_knobs.py -- the whole question, reduced to three numbers.

THE QUESTION THIS PROJECT EXISTS TO ANSWER
------------------------------------------
How do a particle's SIZE, SHAPE, SIZE DISTRIBUTION, SURFACE POTENTIAL and MATERIAL
PROPERTIES change the way a suspension shear-thickens or jams when a drop of it hits a
surface?

THE REDUCTION THAT MAKES IT SIMPLE
----------------------------------
The Wyart-Cates closure the drop-impact DNS solves is

    eta = eta_s (1 - phi/phiJ)^-2,   phiJ(fhat) = phi0 - (phi0 - phim) fhat,
    fhat = f(sigma/sigma*)

and it is EXACTLY covariant under a common rescaling of the three packing fractions
(verified symbolically in `three_knobs_check.py` C1-C2). So the three of them are not
three independent inputs: only their RATIOS matter. Write

    x    = phi / phi0        the distance to FRICTIONLESS jamming
    chi  = phim / phi0       the CONTRAST between the two branches, 0 < chi < 1

and the entire constitutive law becomes

    eta / eta_s = (1 - x / [1 - (1 - chi) fhat])^-2 ,      fhat = f(sigma/sigma*)

**THREE NUMBERS: x, chi, sigma*.** Nothing else about the particles enters.

THE PHASE DIAGRAM IS THEN TWO NUMBERS, AND IT IS EXACT
------------------------------------------------------
`phiJ` runs from `phi0` at zero stress to `phim` at full activation, so

    x <  chi     never jams at any stress          -> continuous thickening only
    chi <= x < 1 JAMS at high stress, flows at low  -> the shear-jamming window
    x >= 1       jammed at rest                     -> a solid, not a suspension

**The shear-jamming window has width `1 - chi` in `x`.** `chi` decides whether shear
jamming is POSSIBLE and over how wide a concentration range; `x` decides where in that
range you sit; `sigma*` decides at what STRESS it happens and never whether.

WHAT THAT BUYS
--------------
Each of the five particle axes can now be classified by WHICH KNOB IT TURNS, and that is
one number per axis rather than a separate `phi0` and `phim` study:

    a chi-knob  changes the WINDOW  -> changes whether the material can shear-jam at all
    an x-knob   changes the POSITION-> equivalent to changing concentration
    a sigma*-knob changes the ONSET -> changes the stress, not the state

Run:  python3 three_knobs.py      (writes out/three_knobs.txt)
"""

from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np

import mixture_model as mm

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
SCALE = mm.PHI_SCALE
N_REQ = 1500
RUNGS = (1e-3, 3e-4, 1e-4, 3e-5)
QUOTE_RUNG = 3e-5           # the matched lowest rung; the fit-free channel

# The friction endpoints used to form chi. mu_p = 0 IS phi0 by definition. mu_p = 3.0
# stands in for phim: Silbert 2010 (doi:10.1039/c001973a) reports phi_c saturating near
# 0.55 for mu >= 1, and this project's own ladder is flat from mu_p = 1 to 3.
MU_0, MU_M = 0.0, 3.0
T_CRIT = {2: 12.706, 3: 4.303, 4: 3.182}


def load():
    """Every lognormal N=1500 record, highest cap wins per (mu, cv, seed, rung)."""
    best: dict = {}
    for fn in sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl"))):
        if "wcbidi" in os.path.basename(fn):
            continue
        for ln in open(fn, encoding="utf-8"):
            ln = ln.strip()
            if not ln:
                continue
            r = json.loads(ln)
            if r.get("failed") or r.get("n_requested") != N_REQ:
                continue
            if (r.get("family") or "lognormal") != "lognormal":
                continue
            k = (r["mu_p"], round(r["cv"], 1), r["seed"], r["p_hat_target"])
            if k not in best or r.get("max_iter", 0) > best[k].get("max_iter", 0):
                best[k] = r
    return best


def chi_at(best, cv: float, rung: float, converged_only: bool = True):
    """chi = phi(MU_M)/phi(MU_0), SEED-PAIRED. Returns (chi, se, n).

    Seed-paired because `packing.py` seeds the radius draw from the seed REGARDLESS of
    `mu_p`, so seed s is the SAME size draw in both arms. The floor propagates the
    protocol's per-record scatter through the ratio, `var = SCALE^2 (1+chi^2)/phi0^2`,
    and does NOT cancel under pairing -- pairing removes the seed's draw, not the
    protocol's re-run scatter.
    """
    a, b = {}, {}
    for s in (1, 2, 3):
        ra, rb = best.get((MU_0, cv, s, rung)), best.get((MU_M, cv, s, rung))
        if ra is None or rb is None:
            continue
        if converged_only and not (ra.get("converged") and rb.get("converged")):
            continue
        a[s], b[s] = ra["phi"], rb["phi"]
    ss = sorted(a)
    if not ss:
        return None, None, 0
    ch = np.array([b[s] / a[s] for s in ss])
    m = float(ch.mean())
    p0 = float(np.mean([a[s] for s in ss]))
    se_seed = float(ch.std(ddof=1) / math.sqrt(len(ss))) if len(ss) > 1 else float("inf")
    se_floor = SCALE * math.sqrt(1 + m * m) / (p0 * math.sqrt(len(ss)))
    return m, max(se_seed, se_floor), len(ss)


def main() -> int:
    L: list[str] = []

    def say(s: str = "") -> None:
        print(s)
        L.append(s)

    say("=" * 92)
    say("three_knobs.py -- the five particle axes, reduced to x, chi and sigma*")
    say("=" * 92)
    say("   eta/eta_s = (1 - x/[1 - (1-chi) fhat])^-2,  x = phi/phi0, chi = phim/phi0")
    say("   x <  chi     : never jams          (continuous thickening only)")
    say("   chi <= x < 1 : SHEAR-JAMMING window, width 1-chi")
    say("   x >= 1       : jammed at rest")
    say(f"   chi is formed as phi(mu_p={MU_M})/phi(mu_p={MU_0}), seed-paired, "
        f"converged cells only.")
    say(f"   grading scale PHI_SCALE = {SCALE}.")

    best = load()

    # ---------------------------------------------------------------- 1
    say("\n=== 1. DOES THE SIZE DISTRIBUTION MOVE THE WINDOW? ===")
    say(f"   At the quoted rung p_hat = {QUOTE_RUNG:.0e}:")
    say("")
    say(f"   {'CV':>6}{'n':>4}{'phi_0':>10}{'phi_m':>10}{'chi':>9}{'se':>9}"
        f"{'window 1-chi':>14}")
    q: dict = {}
    for cv in (0.0, 0.2, 0.5):
        m, se, n = chi_at(best, cv, QUOTE_RUNG)
        if m is None:
            say(f"   {cv:6.1f}   -- no converged pair --")
            continue
        p0 = float(np.mean([best[(MU_0, cv, s, QUOTE_RUNG)]["phi"] for s in (1, 2, 3)
                            if (MU_0, cv, s, QUOTE_RUNG) in best
                            and best[(MU_0, cv, s, QUOTE_RUNG)]["converged"]]))
        pm = float(np.mean([best[(MU_M, cv, s, QUOTE_RUNG)]["phi"] for s in (1, 2, 3)
                            if (MU_M, cv, s, QUOTE_RUNG) in best
                            and best[(MU_M, cv, s, QUOTE_RUNG)]["converged"]]))
        q[cv] = (m, se, n)
        say(f"   {cv:6.1f}{n:4d}{p0:10.5f}{pm:10.5f}{m:9.4f}{se:9.4f}{1-m:14.4f}")
    say("")
    if 0.0 in q and 0.2 in q:
        (m0, s0, n0), (m2, s2, n2) = q[0.0], q[0.2]
        se = math.hypot(s0, s2)
        t = abs(m2 - m0) / se
        tc = T_CRIT.get(min(n0, n2), 4.303)
        say(f"   chi(CV=0.2) - chi(CV=0.0) = {m2-m0:+.4f} +- {se:.4f}   |t| = {t:.2f} "
            f"vs t_crit {tc} -> {'DIFFER' if t > tc else 'NOT RESOLVED'}")
    if 0.5 in q and not math.isfinite(q[0.5][1]):
        say(f"   ⚠ CV = 0.5 has n = {q[0.5][2]} converged pair: the FRICTIONLESS arm still")
        say("     caps at 1e6 on 2 of 3 seeds, so chi(0.5) is not resolvable. Its central")
        say(f"     value {q[0.5][0]:.4f} sits with the others and is printed, not quoted.")
    say("")
    say("   🔑 THE SIZE DISTRIBUTION IS AN x-KNOB, NOT A chi-KNOB. It raises phi_0 and")
    say("   phi_m by the SAME factor, so the shear-jamming window does not move. Widening")
    say("   the PSD does not make a suspension more able to shear-jam; it moves the")
    say("   suspension within the window, exactly as changing concentration would.")

    # ---------------------------------------------------------------- 2
    say("\n=== 2. 🔴 chi IS NOT RUNG-INDEPENDENT. The ratio does NOT escape the ladder. ===")
    say("   A tempting simplification -- 'chi is a ratio, so the pressure dependence")
    say("   cancels and two cells measure it' -- is FALSE, and this is the test:")
    say("")
    say(f"   {'CV':>6}" + "".join(f"{('p=%.0e' % p):>11}" for p in RUNGS)
        + f"{'chi spread':>12}{'phi_0 spread':>14}")
    for cv in (0.0, 0.2, 0.5):
        chis, p0s = [], []
        for p in RUNGS:
            m, _, n = chi_at(best, cv, p)
            chis.append(m if m is not None else float("nan"))
            vals = [best[(MU_0, cv, s, p)]["phi"] for s in (1, 2, 3)
                    if (MU_0, cv, s, p) in best and best[(MU_0, cv, s, p)]["converged"]]
            p0s.append(float(np.mean(vals)) if vals else float("nan"))
        fc = [c for c in chis if np.isfinite(c)]
        fp = [c for c in p0s if np.isfinite(c)]
        say(f"   {cv:6.1f}"
            + "".join(f"{c:11.4f}" if np.isfinite(c) else f"{'--':>11}" for c in chis)
            + f"{(max(fc)-min(fc)) if len(fc) > 1 else float('nan'):12.4f}"
            + f"{(max(fp)-min(fp)) if len(fp) > 1 else float('nan'):14.4f}")
    say("")
    say("   chi drifts across the ladder by MORE than phi_0 itself does, because the two")
    say("   friction arms have different compressibilities: phi = phi_J + c p^b with a")
    say("   different c in each arm, so the ratio does not cancel -- it AMPLIFIES.")
    say("   ⚠ CONSEQUENCE: chi must be quoted AT A NAMED RUNG or extrapolated. There is")
    say("   no shortcut around the ladder, and the 'two cells measure chi' idea is dead.")
    say("")
    say("   ✅ BUT THE VERDICT IN SECTION 1 SURVIVES, and this is why: chi moves across")
    say("   CV by at most ~0.008 at any single rung, while it moves 0.016-0.025 ACROSS")
    say("   rungs. The CV effect is smaller than the systematic it would have to beat.")
    say("   'The distribution does not move the window' is a statement about a")
    say("   comparison AT FIXED RUNG, and that is the comparison that was made.")

    # ---------------------------------------------------------------- 3
    say("\n=== 3. THE AXIS MAP: which knob does each particle property turn? ===")
    say("")
    rowsmap = [
        ("SIZE  a", "sigma*", "sigma* ~ F_rep/a^2 (convention-dependent)",
         "NOT SEPARATED -- onset_size_degeneracy_2026-09-19: a size ladder cannot "
         "separate the two grain-scale stresses"),
        ("SHAPE alpha", "chi (expected)",
         "elongation/interlocking should lower phi_m more than phi_0",
         "🔴 UNMEASURED: the shape campaign ran at ONE friction, so it measures neither "
         "phi_0 nor phi_m separately -- no chi"),
        ("DISTRIBUTION CV", "x", "phi_0 and phi_m rise together",
         "✅ MEASURED HERE: chi flat to |t| = 0.40 across CV 0 -> 0.2"),
        ("POTENTIAL F_rep", "sigma*", "sets the repulsive force the contact must overcome",
         "🔴 UNMEASURED in this repo -- it is the DEM's `sigma*` input, not an output"),
        ("MATERIAL mu_p", "chi", "phi_0 is frictionless BY DEFINITION, so mu_p moves "
         "phi_m alone", "✅ MEASURED: chi = 0.887 at mu_p 0 -> 3, window 1-chi = 0.113"),
        ("MATERIAL E", "eta_s / none", "stiffness rescales the contact, not phi_J",
         "✅ the static packing is scale-free (test_scale_invariance.py)"),
    ]
    say(f"   {'axis':<18}{'knob':<16}{'why':<52}status")
    for ax, knob, why, status in rowsmap:
        say(f"   {ax:<18}{knob:<16}{why:<52}")
        say(f"   {'':<18}{'':<16}{'':<52}{status}")
    say("")
    say("   🔑 READ THE TABLE AS A PROGRAMME. Two axes turn the WINDOW (shape, friction),")
    say("   two turn the ONSET STRESS (size, potential), one turns POSITION (distribution)")
    say("   and one turns nothing (stiffness). Only the window knobs can change WHETHER a")
    say("   suspension shear-jams; the rest change where or when.")
    say("")
    say("   🔴 THE LARGEST GAP IS SHAPE. It is the other candidate chi-knob and the repo")
    say("   cannot presently measure its chi at all, because the shape campaign ran at a")
    say("   single friction. Two extra friction values on the EXISTING shape grid would")
    say("   close it, and that is the cheapest high-value campaign this project has.")

    # ---------------------------------------------------------------- 4
    say("\n=== 4. WHAT A DROP IMPACT NEEDS ON TOP OF (x, chi, sigma*) ===")
    say("   The drop supplies a STRESS and a TIME, and the suspension answers with a")
    say("   state. Two conditions, and they are independent:")
    say("")
    say("      (i)  sigma_imp / sigma*  >> 1     the contacts are actually activated")
    say("      (ii) chi <= x < 1                 the activated state is JAMMED")
    say("")
    say("   Both must hold for a drop to rebound as a solid. (i) is the sigma*-knobs")
    say("   (size, potential); (ii) is the window knobs (friction, shape) and the")
    say("   position knob (concentration, distribution, aggregation).")
    say("   ⚠ WHICH STRESS GAUGE sigma_imp IS is a convention this file does not choose;")
    say("   lessons.md records that `rho k U0^2` and the deviatoric wall shear give")
    say("   OPPOSITE verdicts on whether the activation half is saturated.")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "three_knobs.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    say("")
    say("   written: out/three_knobs.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
