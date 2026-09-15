"""transfer_kladder.py -- does a NON-ADDITIVE combination rule explain the held-out
continuous-family transfer gap?

**DECLARED POST-HOC (PREREG_KLADDER_2026-09-15 G8). NOT A PRE-REGISTERED GATE.** Both
datasets -- the binary calibration set and the three continuous families -- were on disk
before this script was written, so no prediction here was made blind and none is reported
as one.

THE STANDING NUMBER. `grade_mixture_ext.py` transfers the tabulated closure to the
lognormal / uniform / bidisperse families and gets held-out chi2 **108.1** against the
shipped two-exponent J's **34.9**. That is the observation this whole arm exists to
explain: a closure built from binary-MEASURED interactions transfers WORSE than a
smoothed parametric one, and it under-predicts the lognormal monotonically with CV.

THE QUESTION HERE. The ternary block showed the additive combination rule under-fills.
If that is the cause of the transfer gap, replacing the SUM over pairs with the recursive
nested (sequential) rule -- SAME measured interactions, NO new free parameter -- should
close it. If it does not, non-additivity is not the explanation and the gap is something
else.

WHAT IS SWEPT, AND WHY. The sequential rule's only choice is the effective diameter it
gives the coarse composite. The ternary block picked `d32` and `d43` by hand and they
disagreed by rms 3.80 vs 7.19 -- a factor of 1.9 on the headline. Here all four moment
ratios `d10`, `d21`, `d32`, `d43` are evaluated and the SPREAD IS REPORTED BESIDE THE
ANSWER, because a result that moves by more than its own effect when an undeclared
convention is changed is a measurement of the convention.

Two honest limits, printed rather than argued:
  * the continuous families' pair weight extends BELOW u = 0.25, where the tabulated
    a(u), b(u) are interpolated toward the proved u = 0 limit rather than measured. The
    fraction of pair weight below the lowest knot is printed per family.
  * `mm.psd_classes` discretises each family into `mm.K_CLASSES` equal-count classes --
    160, printed by the script rather than quoted here, because the first draft of this
    docstring said 8 and its own output said 160. The sequential rule is order-dependent
    in a way the additive one is not, so the class count is a convention too; it is held
    at the value every other transfer in this project used. Note the BIDISPERSE family
    arrives with K = 2, where the nested rule and the additive one are provably the same
    calculation -- and the table shows all six columns identical on those five rows,
    which is the wiring gate W2 appearing again in the transfer itself.
"""

from __future__ import annotations

import os

import numpy as np

import grade_mixture_ext as gx
import kladder as kl
import mixture_closure as mc
import mixture_model as mm

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = kl.PHI_SCALE
MU = kl.MU_P
CONV = ("d10", "d21", "d32", "d43")


def pair_weight_below(d, y, ucut: float = 0.25) -> float:
    """Fraction of the pairwise volume-product weight sitting at u < `ucut`.

    The linear model's interactions enter weighted by the volume of the classes they
    act between, so y_i * y_j is the natural weight of the pair (i, j). This says how
    much of a family's prediction leans on interpolated rather than measured ground.
    """
    d = np.asarray(d, float)
    y = np.asarray(y, float)
    y = y / y.sum()
    tot = lo = 0.0
    for i in range(d.size):
        for j in range(i + 1, d.size):
            u = min(d[i], d[j]) / max(d[i], d[j])
            w = float(y[i] * y[j])
            tot += w
            if u < ucut:
                lo += w
    return lo / tot if tot > 0 else float("nan")


def main() -> int:
    lines: list[str] = []

    def say(s: str = "") -> None:
        lines.append(s)
        print(s)

    meas = mm.measured(mu_p=MU)
    beta = meas[("lognormal", 0.0)][0]
    pts23 = mc.binary_points(mu_p=MU)
    allpts = gx.ext_points(beta)
    cand, _, _ = gx.candidates(beta, pts23)
    paJ, pbJ = cand["J"](0)
    tab, fa, fb = gx.tabulated_closure(beta, allpts)

    say("=== transfer_kladder.py -- G8, DECLARED POST-HOC, NOT A PRE-REGISTERED GATE ===")
    say(f"   beta = {beta:.6f}; TAB knots at u = "
        f"{', '.join(f'{u:.3f}' for u in sorted(tab))}; J PA {paJ:.4f} PB {pbJ:.4f}")
    say(f"   {mm.K_CLASSES} equal-count classes per family, as every other transfer in "
        "this project uses.")
    say()

    say("=== the held-out family points, ADDITIVE vs the RECURSIVE NESTED rule ===")
    say(f"   {'family':<11}{'CV':>5}{'phi':>9}{'ADDj':>8}{'ADDtab':>8}"
        + "".join(f"{'SEQ' + c:>8}" for c in CONV) + f"{'w(u<.25)':>10}{'K':>4}")
    chi = {"ADDj": 0.0, "ADDtab": 0.0, **{f"SEQ{c}": 0.0 for c in CONV}}
    res = {k: [] for k in chi}
    npt = 0
    for fam in mm.FAMILIES:
        for cv in mc.CVS:
            if (fam, cv) not in meas:
                continue
            d, y = mm.psd_classes(fam, cv)
            phi = meas[(fam, cv)][0]
            v = {"ADDj": float(mm.packing_density(d, y, beta, paJ, pbJ)),
                 "ADDtab": float(mc.packing_density_tabulated(d, y, beta, fa, fb))}
            for c in CONV:
                v[f"SEQ{c}"] = kl.sequential_k(d, y, beta, fa, fb, c)
            for k in chi:
                chi[k] += ((v[k] - phi) / SCALE) ** 2
                res[k].append((v[k] - phi) / SCALE)
            npt += 1
            say(f"   {fam:<11}{cv:5.1f}{phi:9.4f}{v['ADDj']:8.4f}{v['ADDtab']:8.4f}"
                + "".join(f"{v['SEQ' + c]:8.4f}" for c in CONV)
                + f"{pair_weight_below(d, y):10.3f}{d.size:4d}")
    say()
    say(f"=== held-out chi2 over {npt} family points (lower is better) ===")
    say(f"   {'model':>10}{'chi2':>10}{'rms':>8}{'mean signed':>13}{'above':>7}{'below':>7}")
    for k in ("ADDj", "ADDtab", *(f"SEQ{c}" for c in CONV)):
        r = np.array(res[k])
        say(f"   {k:>10}{chi[k]:10.1f}{float(np.sqrt(np.mean(r**2))):8.2f}"
            f"{float(r.mean()):+13.2f}{int((r > 0).sum()):7d}{int((r < 0).sum()):7d}")
    say("   'above' counts points where the MODEL sits above the measurement.")
    say("   record, grade_mixture_ext.txt: J 34.9, tabulated 108.1.")
    say()

    seq = {c: chi[f"SEQ{c}"] for c in CONV}
    best = min(seq, key=lambda c: seq[c])
    say("=== the convention is SWEPT, not chosen -- and this is what decides how much "
        "the answer is worth ===")
    say("   SEQ chi2 across the four composite-diameter conventions: "
        + ", ".join(f"{c} {seq[c]:.1f}" for c in CONV))
    say(f"   best {best} at {seq[best]:.1f}; spread {min(seq.values()):.1f} to "
        f"{max(seq.values()):.1f}, a factor "
        f"{max(seq.values()) / max(min(seq.values()), 1e-9):.2f}")
    say(f"   ADDtab {chi['ADDtab']:.1f} -> best SEQ {seq[best]:.1f}: "
        f"{'closes' if seq[best] < chi['ADDj'] else 'does NOT close'} the gap to J "
        f"({chi['ADDj']:.1f}).")
    gap = chi["ADDtab"] - chi["ADDj"]
    closed = (chi["ADDtab"] - seq[best]) / gap if gap else float("nan")
    say(f"   fraction of the ADDtab - J gap closed by the best convention: {closed:+.2f}")
    say(f"   fraction closed by the WORST convention: "
        f"{(chi['ADDtab'] - max(seq.values())) / gap if gap else float('nan'):+.2f}")
    say("   If those two fractions straddle zero, the convention decides the verdict and")
    say("   the honest statement is that this comparison does not settle the question.")
    say()

    say("=== the lognormal's monotone under-prediction with CV -- the specific symptom ===")
    say(f"   {'CV':>5}{'phi':>9}{'ADDtab res':>12}" + "".join(f"{'SEQ' + c:>9}" for c in CONV))
    for cv in mc.CVS:
        if ("lognormal", cv) not in meas:
            continue
        d, y = mm.psd_classes("lognormal", cv)
        phi = meas[("lognormal", cv)][0]
        rt = (float(mc.packing_density_tabulated(d, y, beta, fa, fb)) - phi) / SCALE
        rs = [(kl.sequential_k(d, y, beta, fa, fb, c) - phi) / SCALE for c in CONV]
        say(f"   {cv:5.1f}{phi:9.4f}{rt:+12.2f}" + "".join(f"{v:+9.2f}" for v in rs))
    say("   (negative = the model sits BELOW the measurement, i.e. under-predicts)")
    say()
    say("=== what this can and cannot decide ===")
    say("   CAN: whether swapping the SUM for a nested rule, with the same measured")
    say("   interactions and no new parameter, moves the held-out chi2 toward J's.")
    say("   CANNOT: which nesting is right -- the convention spread above is the honest")
    say("   width -- and nothing about u < 0.25, where the interactions are interpolated")
    say("   toward a proved limit rather than measured, at the pair weights printed above.")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "transfer_kladder.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
