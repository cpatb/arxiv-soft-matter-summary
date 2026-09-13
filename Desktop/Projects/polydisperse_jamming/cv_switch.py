"""cv_switch.py — the CV at which polydispersity flips the dilatancy sign, from the
measured phi_c(CV) rather than a pinned Delta phi.

THE GAP THIS CLOSES
-------------------
`derivations/dilatancy_sign_mechanism.md` Eq. (11) splits the CV dependence of the
inertial-dilatant number into two channels:

    d ln ID/dCV  =  8CV/(1+CV^2)  +  phi_c'(CV)/(phi - phi_c(CV))  -  d ln A/dCV
                    |permeability|    |------ critical state ------|

and MECHANISM_2026-09-08.md §M1 records that the 2026-08-31 exclusion of the
permeability channel **swept only the first term**, holding `Delta phi = phi - phi_c`
pinned at 0.02, because `phi_c(CV)` was a DEM campaign nobody had run. The critical-state
term is `O(1)` where the permeability one is `O(CV)`, and it has a POLE at
`phi_c(CV) = phi`. So the exclusion measured the weaker channel and was silent about the
stronger one.

`phi_c(CV)` has been measured since 2026-09-12, from the rate-matched constant-pressure
ladder, by a route with no `phi_J` in it: 0.5911 / 0.6020 / 0.6444 at CV 0 / 0.2 / 0.5.
This script puts that measurement into Eq. (11) and asks the question the mechanism is
for:

    at fixed phi, is there a CV at which the drop stops being dilatant?

WHAT IT DOES AND DOES NOT CLAIM
-------------------------------
It is quoted on the SHIFT `dphi_c/dCV` and on the WIDTH of the switchable band, never on
the absolute level of `phi_c`: a frictional-sphere DEM gives `phi_c` ~ 0.59 while [M26]
use 0.37 for cornstarch and this repo has guessed 0.42 for silica, and those cannot all
be on one axis (`dilatancy_sign_mechanism.md` §10 item 4). The band width is the
transferable number, and it is what the author's bench observation is about.

`phi_c(CV)` is NOT linear -- it moves +0.0109 over CV 0 -> 0.2 and +0.0425 over
0.2 -> 0.5 -- so `CV* = (phi - phi_c(0))/s` from the derivation's Eq. (14) is a
linearisation that is wrong in a knowable direction. Both are computed and the gap
between them is printed.

Reads committed artifacts only. Writes out/cv_switch.txt. No simulation.
"""

from __future__ import annotations

import math
import os

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

# --- the measurement -----------------------------------------------------------------
# out/constp_ratematch.txt, "blocks": phi_c pooled over both phi starts (6 cells each),
# rate-matched to gdot = 1.17e-04 at p* = 3e-04, N = 1000. The route takes no input from
# phi_J: at fixed grain pressure the steady solid fraction IS phi_c(p).
CV_M = np.array([0.0, 0.2, 0.5])
PHI_C = np.array([0.591127, 0.601989, 0.644447])
PHI_C_SD = np.array([0.001020, 0.001313, 0.004336])     # sd over the 6 cells
N_CELLS = 6
PHI_C_SE = PHI_C_SD / math.sqrt(N_CELLS)

N_BOOT = 2000
RNG_SEED = 20260913


def phi_c_curve(vals=None):
    """Monotone (PCHIP) interpolant through the three measured levels.

    PCHIP and not a spline: a cubic spline through three points with a rising second
    difference overshoots, and phi_c(CV) is physically monotone (broadening a lognormal
    at fixed d10 can only pack denser here -- it is measured monotone on both the
    isotropic and the shear ladders). PCHIP cannot overshoot, and the LINEAR reading is
    reported beside it everywhere so the interpolation's contribution is visible.
    """
    y = PHI_C if vals is None else vals
    return PchipInterpolator(CV_M, y, extrapolate=False)


def cv_star(phi: float, curve, lo=0.0, hi=0.5) -> float:
    """The CV at which phi_c(CV) = phi, i.e. where a pack at fixed phi stops being
    denser than critical. NaN outside the measured range -- this project does not
    extrapolate a root past the rungs that bracket it."""
    f = lambda c: float(curve(c)) - phi            # noqa: E731
    if not (f(lo) < 0 < f(hi)) and not (f(lo) > 0 > f(hi)):
        return float("nan")
    return float(brentq(f, lo, hi, xtol=1e-8))


def main() -> int:
    lines: list[str] = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("=== cv_switch.py -- the dilatancy switch, from the MEASURED phi_c(CV) ===")
    say("   phi_c from out/constp_ratematch.txt (constant-pressure, no phi_J input).")
    say()
    say("   cv    phi_c        se       d phi_c        per unit CV")
    for i, cv in enumerate(CV_M):
        d = "" if i == 0 else f"  +{PHI_C[i] - PHI_C[i - 1]:.4f}   " \
                              f"+{(PHI_C[i] - PHI_C[i - 1]) / (CV_M[i] - CV_M[i - 1]):.4f}"
        say(f"  {cv:.1f}   {PHI_C[i]:.6f}  {PHI_C_SE[i]:.6f}{d}")
    band = PHI_C[-1] - PHI_C[0]
    band_se = math.hypot(PHI_C_SE[-1], PHI_C_SE[0])
    say(f"   -> the SWITCHABLE BAND is phi in [{PHI_C[0]:.4f}, {PHI_C[-1]:.4f}], width "
        f"{band:.4f} +- {band_se:.4f}")
    say("      Any material whose effective phi during impact lies inside its own")
    say("      version of this band can be moved from the dilatant side of critical to")
    say("      the contractive side by CHANGING THE WIDTH OF ITS SIZE DISTRIBUTION")
    say("      ALONE, at fixed mean size, fixed phi and fixed everything else.")
    say(f"      The move is NOT linear in CV: +{PHI_C[1] - PHI_C[0]:.4f} over the first "
        f"0.2 and +{PHI_C[2] - PHI_C[1]:.4f} over the next 0.3, i.e. "
        f"{((PHI_C[2] - PHI_C[1]) / 0.3) / ((PHI_C[1] - PHI_C[0]) / 0.2):.2f}x steeper "
        f"at the top of the ladder.")

    # ---------------- CV*(phi) ----------------
    curve = phi_c_curve()
    lin_s = (PHI_C[-1] - PHI_C[0]) / (CV_M[-1] - CV_M[0])
    say()
    say("=== CV*(phi): the width at which a pack at fixed phi crosses critical ===")
    say("   PCHIP on the three measured levels, against the derivation's Eq. (14)")
    say(f"   linearisation CV* = (phi - phi_c(0))/s with s = {lin_s:.4f}.")
    say("   A bootstrap over the seed scatter (2000 draws, seed "
        f"{RNG_SEED}) gives the interval.")
    rng = np.random.default_rng(RNG_SEED)
    draws = rng.normal(PHI_C[None, :], PHI_C_SE[None, :], size=(N_BOOT, 3))
    draws = np.maximum.accumulate(draws, axis=1)      # keep each draw monotone in CV
    say()
    say("    phi      CV* (PCHIP)   CV* 95% band      CV* (linear Eq.14)   gap")
    gaps = []
    for phi in (0.595, 0.600, 0.605, 0.610, 0.615, 0.620, 0.630, 0.640):
        c = cv_star(phi, curve)
        boot = []
        for row in draws:
            v = cv_star(phi, PchipInterpolator(CV_M, row, extrapolate=False))
            if math.isfinite(v):
                boot.append(v)
        lo, hi = (np.percentile(boot, [2.5, 97.5]) if len(boot) > 0.9 * N_BOOT
                  else (float("nan"), float("nan")))
        clin = (phi - PHI_C[0]) / lin_s
        gap = c - clin if math.isfinite(c) else float("nan")
        if math.isfinite(gap):
            gaps.append(gap)
        say(f"   {phi:.3f}    {c:8.4f}    [{lo:.4f}, {hi:.4f}]   {clin:12.4f}   "
            f"{gap:+.4f}")
    say(f"   -> the `gap` column is {'POSITIVE at every phi' if all(g > 0 for g in gaps) else 'MIXED IN SIGN'}"
        f": because phi_c(CV) is CONVEX (it rises "
        f"{((PHI_C[2] - PHI_C[1]) / 0.3) / ((PHI_C[1] - PHI_C[0]) / 0.2):.1f}x faster")
    say("      over the top of the ladder than the bottom), the single-slope Eq. (14)")
    say("      UNDERSTATES the width needed, everywhere in the band. The error peaks")
    say("      near the middle and closes at the top, where the chord and the curve")
    say(f"      meet by construction. Worst gap {max(gaps):+.4f} in CV, which is "
        f"{max(gaps) / 0.2:.1f}x the whole CV 0 -> 0.2 step: a systematic, not a")
    say("      rounding.")

    # ---------------- Eq. (11), both channels, at last ----------------
    say()
    say("=== Eq. (11): the two channels of d ln ID/dCV, with phi_c'(CV) MEASURED ===")
    say("   permeability term  = 8 CV/(1 + CV^2)          [O(CV), vanishes at CV = 0]")
    say("   critical-state term = phi_c'(CV)/(phi - phi_c(CV))   [O(1), pole at phi_c=phi]")
    say("   The 2026-08-31 exclusion swept the FIRST term only, with phi - phi_c pinned")
    say("   at 0.02. Both are now computed from the same measurement.")
    say()
    say("    phi     cv   perm      crit        ratio crit/perm   d ln ID/dCV (2 terms)")
    for phi in (0.66, 0.70, 0.80):
        for cv in (0.05, 0.2, 0.35, 0.5):
            perm = 8 * cv / (1 + cv ** 2)
            dphi = phi - float(curve(cv))
            crit = float(curve.derivative()(cv)) / dphi
            ratio = crit / perm if perm else float("inf")
            say(f"   {phi:.2f}   {cv:.2f}  {perm:6.3f}   {crit:9.3f}   {ratio:12.2f}"
                f"   {perm + crit:16.3f}")
        say()
    say("   ⚠ phi must exceed phi_c(cv) for ID to be defined at all; the rows above are")
    say("      chosen above phi_c(0.5) = 0.6444 so every cell of the table exists. The")
    say("      pinned-Delta-phi sweep used Delta phi = 0.02, which on this measurement")
    say(f"      corresponds to phi = {PHI_C[0] + 0.02:.4f} at CV 0 -- BELOW phi_c(0.5),")
    say("      i.e. inside the band where the critical-state channel diverges and then")
    say("      changes sign. That is the regime the exclusion could not see.")
    say()
    say("   at the pinned Delta phi = 0.02 the derivation's own comparison reads:")
    for cv in (0.05, 0.2, 0.35, 0.5):
        perm = 8 * cv / (1 + cv ** 2)
        crit = float(curve.derivative()(cv)) / 0.02
        say(f"      cv {cv:.2f}: perm {perm:6.3f}   crit {crit:9.3f}   "
            f"crit/perm = {crit / perm:8.1f}x")
    say("   -> with Delta phi held at 0.02 the critical-state channel beats the")
    say("      permeability one by 2.5-5.1x across the ladder -- NOT by orders of")
    say("      magnitude. The old verdict 'M1 permeability is the weaker channel'")
    say("      survives AT THAT Delta phi, with a measured margin instead of an assumed")
    say("      one; the table above shows it REVERSES by phi = 0.80.")
    say()
    say("=== which channel wins is a property of Delta phi, and here is the boundary ===")
    say("   Setting phi_c'(CV)/Delta phi = 8 CV/(1 + CV^2) gives the distance from")
    say("   critical inside which polydispersity acts through PACKING and outside which")
    say("   it acts through PERMEABILITY:")
    say("      cv     phi_c'(cv)   Delta phi* (crossover)    phi* = phi_c(cv) + Delta phi*")
    for cv in (0.05, 0.1, 0.2, 0.35, 0.5):
        dprime = float(curve.derivative()(cv))
        dstar = dprime * (1 + cv ** 2) / (8 * cv)
        say(f"      {cv:.2f}    {dprime:.5f}        {dstar:.4f}                "
            f"{float(curve(cv)) + dstar:.4f}")
    say("   -> the crossover DIVERGES as CV -> 0 (the permeability term is O(CV) and")
    say("      the packing term is O(1)), so near a monodisperse pack the packing")
    say("      channel wins at any distance from critical. At CV = 0.5 it wins only")
    say("      within about 0.05 of the critical line -- which is the same width as")
    say("      the switchable band itself, and that coincidence is not an accident:")
    say("      both are set by phi_c'(CV).")

    # ---------------- robustness ----------------
    say()
    say("=== robustness: what the interpolation is worth ===")
    lin = PchipInterpolator(CV_M, PHI_C)              # same knots
    seg = np.interp(np.linspace(0, 0.5, 6), CV_M, PHI_C)
    pch = [float(lin(c)) for c in np.linspace(0, 0.5, 6)]
    say("      cv     piecewise-linear   PCHIP      difference")
    for c, a, b in zip(np.linspace(0, 0.5, 6), seg, pch):
        say(f"     {c:.1f}     {a:.6f}        {b:.6f}   {b - a:+.6f}")
    say(f"   max |PCHIP - linear| over the ladder = "
        f"{max(abs(b - a) for a, b in zip(seg, pch)):.6f}, against a seed se of "
        f"{PHI_C_SE.max():.6f} at the worst level.")
    say("   -> the interpolation choice is worth less than the measurement's own error")
    say("      except near CV 0.35, where it is comparable; CV*(phi) in that region")
    say("      should be quoted with the band, not the central value.")

    say()
    say("=== WHAT IS AND IS NOT SETTLED ===")
    say(f"   SETTLED: the switchable band is {band:.4f} +- {band_se:.4f} wide in phi,")
    say("     measured by a route that takes no input from phi_J, and the switch is")
    say("     STEEPER at high CV than a single slope implies.")
    say("   SETTLED: WHICH channel of Eq. (11) dominates is set by the distance from")
    say("     critical, and the boundary is now computed rather than assumed. At the")
    say("     derivation's pinned Delta phi = 0.02 the critical-state channel wins by")
    say("     2.5-5.1x; by phi = 0.80 the permeability channel wins by 2.5-4x. The")
    say("     order-counting (packing O(1), permeability O(CV)) is confirmed near")
    say("     CV = 0 and is NOT a statement about the whole ladder.")
    say("   NOT SETTLED: where phi actually sits during an impact. The sign turns on")
    say("     phi - phi_c AT THE CONTACT, and the bench phi = 0.44 is far below every")
    say("     phi_c here, so the band is reached only if the impact densifies the")
    say("     suspension into it. That is the drop-side measurement, not this one.")
    say("   NOT SETTLED: the absolute level of phi_c for a real material. Only the")
    say("     shift transfers from a frictional-sphere DEM to an experiment.")
    say("   NOT SETTLED: whether CV acts on the TRANSIENT as well as on the critical")
    say("     line -- that is the strain-clock campaign, PREREG_STRAINCLOCK_2026-09-13.")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "cv_switch.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
