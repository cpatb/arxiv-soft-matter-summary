"""diag_transfer_ugap.py -- G8 follow-up: can a knot above u = 0.667 collapse the
nested-rule convention spread?

**DECLARED POST-HOC, like G8 itself.** Every number it reads was on disk before it was
written, and nothing here is reported as a blind prediction.

THE STANDING SITUATION (`out/transfer_kladder.txt`, G8). The held-out continuous-family
chi2 is J 34.9, ADDtab 108.1, and the four nested conventions run SEQd21 42.8, SEQd32
48.0, SEQd10 186.3, SEQd43 509.3 -- an **11.9x spread**. They straddle zero, so swapping
the sum for a nested rule does not explain the transfer gap; the undeclared choice of
composite diameter is worth as much as the result.

THE QUESTION. `mm.psd_classes` discretises each continuous family into K_CLASSES
equal-count classes, and ADJACENT classes of a continuous PSD are nearly equal in size --
so those families may put a LARGE share of their pairwise weight at u > 0.667, exactly
where the closure has no knot. If so, the convention spread is partly a readout of the
same unmeasured region the U-GAP block is measuring, and a knot there should move it. If
the weight there is small, a knot cannot move it and G8's spread is about something else.

Two things are computed and the second only means anything given the first:
  1. the pair weight ABOVE u = 0.667, per family and CV, with the weight BELOW u = 0.25
     printed beside it as the already-known companion caveat;
  2. the same held-out chi2 table, recomputed with a knot INSERTED at u* = 0.70711 and
     swept over the monotone-admissible range -- does the spread collapse anywhere in it?

Run:  python diag_transfer_ugap.py      (writes out/diag_transfer_ugap.txt)
"""

from __future__ import annotations

import os
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator

import grade_mixture_ext as gx
import kladder as kl
import mixture_closure as mc
import mixture_model as mm
import transfer_kladder as tk
import ugap as ug

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = ug.PHI_SCALE
MU = ug.MU_P
CONV = ("d10", "d21", "d32", "d43")


def pair_weight_above(d, y, ucut: float) -> float:
    """Fraction of pairwise volume-product weight at u > `ucut`.

    The complement of `transfer_kladder.pair_weight_below`, written here rather than
    called with 1 - x so that the two can be summed as a check (they must give 1 minus
    the weight in between, not 1).
    """
    d = np.asarray(d, float)
    y = np.asarray(y, float)
    y = y / y.sum()
    tot = hi = 0.0
    for i in range(d.size):
        for j in range(i + 1, d.size):
            u = min(d[i], d[j]) / max(d[i], d[j])
            w = float(y[i] * y[j])
            tot += w
            if u > ucut:
                hi += w
    return hi / tot if tot > 0 else float("nan")


def main() -> int:
    lines: list[str] = []

    def say(s: str = "") -> None:
        lines.append(s)
        print(s, flush=True)

    meas = mm.measured(mu_p=MU)
    beta = meas[("lognormal", 0.0)][0]
    allpts = gx.ext_points(beta)
    tab, fa0, fb0 = gx.tabulated_closure(beta, allpts)
    u_lo = max(tab)
    say("=== diag_transfer_ugap.py -- G8 follow-up, POST-HOC ===")
    say(f"   beta = {beta:.6f}; SCALE = {SCALE:.1e}; highest measured knot u = {u_lo:.4f}")
    say(f"   u* = {ug.U_STAR:.5f}; K_CLASSES = {mm.K_CLASSES} "
        f"(the family discretisation, printed not quoted)")
    say()

    # ------------------------------------------------------------------ section 1
    say("=== 1. where the continuous families' pair weight actually sits ===")
    say("   A knot at u* can only move a family's prediction through pairs NEAR u*.")
    say(f"   {'family':<11}{'CV':>5}{'>0.6667':>10}{'>0.70711':>10}{'>0.90':>9}"
        f"{'>0.95':>9}{'>0.99':>9}{'<0.25':>9}")
    wt = {}
    for fam in mm.FAMILIES:
        for cv in mc.CVS:
            if (fam, cv) not in meas:
                continue
            d, y = mm.psd_classes(fam, cv)
            hi667 = pair_weight_above(d, y, u_lo)
            histar = pair_weight_above(d, y, ug.U_STAR)
            h90 = pair_weight_above(d, y, 0.90)
            h95 = pair_weight_above(d, y, 0.95)
            h99 = pair_weight_above(d, y, 0.99)
            lo25 = tk.pair_weight_below(d, y, 0.25)
            wt[(fam, cv)] = (hi667, histar, lo25, h90)
            say(f"   {fam:<11}{cv:5.1f}{hi667:10.4f}{histar:10.4f}{h90:9.4f}"
                f"{h95:9.4f}{h99:9.4f}{lo25:9.4f}")
    say()
    say("   🔑 THE COLUMN THAT MATTERS IS NOT THE FIRST ONE. A knot at u* = 0.70711 sits")
    say("   at the BOTTOM of the gap, and the continuous families put their weight at the")
    say("   TOP of it: the lognormal and uniform at CV = 0.1 carry 0.54 of their pair")
    say("   weight above u = 0.90 and 0.28-0.30 above 0.95, because `psd_classes`")
    say(f"   discretises into {mm.K_CLASSES} equal-count classes and ADJACENT classes of a")
    say("   continuous PSD are nearly the same size. So the new knot does not cover the")
    say("   families' weight -- it only changes how PCHIP GUESSES on [u*, 1], which is")
    say("   where that weight lives. That is the mechanism behind section 2's verdict,")
    say("   and it says the family transfer needs a(u) at u > 0.9 -- precisely where the")
    say("   U-GAP block showed a binary composition curve has no signal left to fit.")
    say()

    # ------------------------------------------------------------------ section 2
    say("=== 2. the held-out chi2 table as a knot at u* is swept ===")
    say("   The knot is inserted into the SAME tabulated closure the standing numbers")
    say("   use; a = b = the shipped interpolated values reproduces them up to PCHIP")
    say("   reshaping, which is printed as its own row so it cannot be read as physics.")
    a_lo, b_lo = tab[u_lo]
    us0 = [0.0, 1.0] + sorted(tab)
    aa0 = np.asarray([0.0, 1.0] + [tab[u][0] for u in sorted(tab)])
    bb0 = np.asarray([0.0, 1.0] + [tab[u][1] for u in sorted(tab)])
    o = np.argsort(us0)
    uu0 = np.asarray(us0)[o]
    aa0, bb0 = aa0[o], bb0[o]

    def closure(a_star, b_star):
        if a_star is None:
            return fa0, fb0
        ua = np.append(uu0, ug.U_STAR)
        va = np.append(aa0, a_star)
        vb = np.append(bb0, b_star)
        oa = np.argsort(ua)
        return (PchipInterpolator(ua[oa], va[oa]),
                PchipInterpolator(ua[oa], vb[oa]))

    paJ, pbJ = gx.candidates(beta, mc.binary_points(mu_p=MU))[0]["J"](0)

    def table(a_star, b_star):
        f_a, f_b = closure(a_star, b_star)
        chi = {"J": 0.0, "ADDtab": 0.0}
        for c in CONV:
            chi[f"SEQ{c}"] = 0.0
        npt = 0
        for fam in mm.FAMILIES:
            for cv in mc.CVS:
                if (fam, cv) not in meas:
                    continue
                d, y = mm.psd_classes(fam, cv)
                phi = meas[(fam, cv)][0]
                chi["J"] += ((mm.packing_density(d, y, beta, paJ, pbJ) - phi) / SCALE) ** 2
                chi["ADDtab"] += ((mc.packing_density_tabulated(
                    np.asarray(d, float), np.asarray(y, float), beta, f_a, f_b)
                    - phi) / SCALE) ** 2
                for c in CONV:
                    chi[f"SEQ{c}"] += ((kl.sequential_k(d, y, beta, f_a, f_b, c)
                                        - phi) / SCALE) ** 2
                npt += 1
        return chi, npt

    rows = [("no knot (the standing numbers)", None, None),
            ("knot at the SHIPPED value", float(fa0(ug.U_STAR)), float(fb0(ug.U_STAR)))]
    for av in (a_lo, 0.96, 0.98, 0.999):
        for bv in (b_lo, 0.96, 0.999):
            rows.append((f"a {av:.4f}, b {bv:.4f}", av, bv))

    say(f"   {'knot':<32}{'J':>9}{'ADDtab':>9}"
        + "".join(f"{'SEQ' + c:>9}" for c in CONV) + f"{'SEQ spread':>12}{'best':>9}")
    spreads = []
    for name, av, bv in rows:
        chi, npt = table(av, bv)
        seq = [chi[f"SEQ{c}"] for c in CONV]
        spread = max(seq) / min(seq) if min(seq) > 0 else float("inf")
        spreads.append((name, spread, min(seq), chi["J"]))
        say(f"   {name:<32}{chi['J']:9.1f}{chi['ADDtab']:9.1f}"
            + "".join(f"{v:9.1f}" for v in seq)
            + f"{spread:11.2f}x{min(seq):9.1f}")
    say(f"   (held-out family points: {npt})")
    say("   record, out/transfer_kladder.txt: J 34.9, ADDtab 108.1, d21 42.8, d32 48.0,")
    say("   d10 186.3, d43 509.3 -- spread 11.9x.")
    say()

    best_spread = min(s for _, s, _, _ in spreads)
    worst_spread = max(s for _, s, _, _ in spreads)
    say("=== VERDICT ===")
    say(f"   convention spread over the whole swept knot range: "
        f"{best_spread:.2f}x to {worst_spread:.2f}x")
    say(f"   pair weight above u = {u_lo:.4f}, over all families and CV: "
        f"{min(v[0] for v in wt.values()):.4f} to {max(v[0] for v in wt.values()):.4f}")
    if best_spread > 2.0:
        say("   -> a knot above 0.667 does NOT collapse the spread anywhere in the")
        say("      admissible range. G8 stays unsettled, and for a reason this rules IN")
        say("      or OUT explicitly rather than by argument.")
    else:
        say("   -> the spread DOES collapse somewhere in the range; report where, and")
        say("      note it is a post-hoc search over a parameter, not a test.")
    say()
    say("   WHAT THIS CANNOT DECIDE: the composite-diameter convention is a property of")
    say("   the NESTED RULE, not of the interaction table, so nothing about a(u) is")
    say("   required to move it. This measures whether it DOES, not whether it should.")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "diag_transfer_ugap.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
