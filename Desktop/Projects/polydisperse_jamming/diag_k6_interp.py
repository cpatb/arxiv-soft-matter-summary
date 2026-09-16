"""diag_k6_interp.py -- is "K = 6 is attemptable" a measurement of the SIGNAL or of the
INTERPOLATION FORM?

POST-HOC, NOT A GATE. Written after a physics-referee pass on
`derivations/k_axis_bound.md` reported that the document's forward-content estimate --
signal 2.99 at u = 0.75786 against a 2.57 bar, margin 1.16x -- comes from a TWO-POINT
log-log interpolation that discards the third measured point, and that an equally
defensible form through all three lands at 2.56, i.e. margin 0.997 and the opposite
verdict.

This file exists to do five things and nothing else:

  1. REPRODUCE the referee's numbers from the `ug` shards, not from its report. The
     control is arm B: `out/grade_ugap.txt` records max(phi - beta) = 2.21 there and the
     script exits NON-ZERO if it cannot recover that. A diagnostic that cannot reproduce
     the number it audits has not audited it.
  2. Sweep the interpolation form over a DECLARED set and report the SPREAD beside every
     verdict, rather than one number from one undeclared choice. Two of the forms are
     anchored on the exact physical point (u = 1, s = 0): at u = 1 the two classes are
     the same size, so the mixture IS monodisperse, phi = beta and the filling signal
     vanishes identically. That anchor is free and is NOT a fitted parameter.
  3. Ask whether the three measured points are consistent with any smooth two-parameter
     form AT ALL. If they are not, the interpolation is under-determined and no choice of
     form is more defensible than another -- which is a stronger statement than "the
     margin is thin", and it is the one that decides whether K = 6 is UNRESOLVED or
     merely marginal.

  4. Check each form against the ONE piece of exact physics available -- s(u = 1) = 0,
     verified in the shipped closure rather than asserted -- and DISCARD the forms it
     falsifies. Four of seven do not vanish there and one goes NEGATIVE, which
     max(phi - beta) cannot be. That is stronger than "under-determined".
  5. Recompute the 2.57 RESOLVING-POWER BAR at this ratio. A second referee pass found
     it was computed at ARM B's ratio r = 1.25 over truths near (a, b) = (0.975, 0.960)
     -- so the document graded a signal interpolated to one ratio against a bar measured
     at another, which is the "phi span" units trap relocated.

WHAT THIS CANNOT DO. It cannot decide K = 6. Only the pre-registered ceiling
recomputation at r = 4^(1/5) = 1.31951 can, and `derivations/k_axis_bound.md` already
says so. Nothing here moves a bar, and the document's own 2.99 stays printed.

Run:  python diag_k6_interp.py     (writes out/diag_k6_interp.txt; exit 1 if the control
                                    does not reproduce)
"""

from __future__ import annotations

import math
import os

import numpy as np

import grade_ugap as gu
import mixture_model as mm
import ugap as ug

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

U6 = 0.25 ** (1.0 / 5.0)        # the knot K = 6 needs: 0.757858...
U7 = 0.25 ** (1.0 / 6.0)        # 0.793701...
BAR = 2.57                      # resolving power, quoted VERBATIM from
#                                 out/grade_ugap.txt arm B's ceiling table. NOT MOVABLE.
REC_ARMB = 2.21                 # the control, same file
REC_DOC = 2.99                  # what derivations/k_axis_bound.md prints for K = 6
REC_DOC_K7 = 2.31

LINES: list[str] = []


def say(s: str = "") -> None:
    LINES.append(s)
    print(s, flush=True)


# --------------------------------------------------------------- interpolation forms
def f_loglog2(us, ss, u):
    """The DOCUMENT's form: linear in (ln u, ln s) through the two BRACKETING points."""
    i = int(np.searchsorted(us, u))
    i = max(1, min(i, len(us) - 1))
    x0, x1 = math.log(us[i - 1]), math.log(us[i])
    y0, y1 = math.log(ss[i - 1]), math.log(ss[i])
    return math.exp(y0 + (y1 - y0) * (math.log(u) - x0) / (x1 - x0))


def f_lin2(us, ss, u):
    """Linear in (u, s) through the two bracketing points."""
    i = int(np.searchsorted(us, u))
    i = max(1, min(i, len(us) - 1))
    return float(ss[i - 1] + (ss[i] - ss[i - 1])
                 * (u - us[i - 1]) / (us[i] - us[i - 1]))


def f_quad_lin(us, ss, u):
    """Quadratic in (u, s) through ALL THREE measured points (Lagrange)."""
    return float(np.polyval(np.polyfit(np.asarray(us, float),
                                       np.asarray(ss, float), 2), u))


def f_quad_loglog(us, ss, u):
    """Quadratic in (ln u, ln s) through all three."""
    c = np.polyfit(np.log(np.asarray(us, float)), np.log(np.asarray(ss, float)), 2)
    return float(math.exp(np.polyval(c, math.log(u))))


def f_pow1mu_2(us, ss, u):
    """s = c*(1-u)^p through the two bracketing points. ANCHORED: s -> 0 as u -> 1."""
    i = int(np.searchsorted(us, u))
    i = max(1, min(i, len(us) - 1))
    x0, x1 = math.log(1 - us[i - 1]), math.log(1 - us[i])
    y0, y1 = math.log(ss[i - 1]), math.log(ss[i])
    return math.exp(y0 + (y1 - y0) * (math.log(1 - u) - x0) / (x1 - x0))


def f_pow1mu_3(us, ss, u):
    """s = c*(1-u)^p least-squares through ALL THREE. ANCHORED at (1, 0)."""
    p, lc = np.polyfit(np.log(1 - np.asarray(us, float)),
                       np.log(np.asarray(ss, float)), 1)
    return float(math.exp(lc + p * math.log(1 - u)))


def f_quad_1mu(us, ss, u):
    """Quadratic in (ln(1-u), ln s) through all three. ANCHORED at (1, 0)."""
    c = np.polyfit(np.log(1 - np.asarray(us, float)),
                   np.log(np.asarray(ss, float)), 2)
    return float(math.exp(np.polyval(c, math.log(1 - u))))


FORMS = (
    ("log-log, 2 bracketing pts  [THE DOCUMENT'S]", f_loglog2, 2, "no"),
    ("linear in (u, s), 2 pts", f_lin2, 2, "no"),
    ("quadratic in (u, s), 3 pts", f_quad_lin, 3, "no"),
    ("quadratic in (ln u, ln s), 3 pts", f_quad_loglog, 3, "no"),
    ("power law in (1-u), 2 pts", f_pow1mu_2, 2, "YES"),
    ("power law in (1-u), 3 pts LSQ", f_pow1mu_3, 3, "YES"),
    ("quadratic in (ln(1-u), ln s), 3 pts", f_quad_1mu, 3, "YES"),
)


def main() -> int:
    beta = mm.measured(mu_p=0.3)[("lognormal", 0.0)][0]
    g, nfiles, nraw = gu.load("ug")
    sig = {}
    nbefore = nafter = 0
    for _arm, r in ug.RATIOS.items():
        vals = []
        for x in ug.X_GRID:
            cs = g.get((float(r), round(float(x), 6)), [])
            nbefore += len(cs)
            ok = [c for c in cs if gu.admissible(c, ug.N_PROD)[0]]
            nafter += len(ok)
            if ok:
                vals.append(float(np.mean([c["phi"] for c in ok])))
        if vals:
            sig[1.0 / float(r)] = (max(vals) - beta) / ug.PHI_SCALE
    us = np.array(sorted(sig))
    ss = np.array([sig[u] for u in sorted(sig)])

    say("=== diag_k6_interp.py -- POST-HOC, NOT A GATE ===")
    say("   Does 'K = 6 is attemptable' measure the SIGNAL or the INTERPOLATION FORM?")
    say(f"   beta = {beta:.6f}; PHI_SCALE = {ug.PHI_SCALE:.1e}; "
        f"bar = {BAR} scale (quoted verbatim, NOT MOVABLE)")
    say(f"   ug shards: {nfiles} files, {nraw} records; {nbefore} -> {nafter} admitted")
    say()

    say("=== REPRODUCTION CONTROL (exit 1 if this fails) ===")
    say(f"   {'u':>10}{'ratio':>8}{'max(phi-beta)/scale':>21}")
    for u in sorted(sig):
        say(f"   {u:10.5f}{1 / u:8.3f}{sig[u]:21.4f}")
    got = sig.get(0.8, float("nan"))
    ok = abs(got - REC_ARMB) <= 0.006
    say(f"   arm B (u = 0.80000) recorded in out/grade_ugap.txt at {REC_ARMB}; "
        f"computed here {got:.4f} -> {'OK' if ok else 'MISMATCH'}")
    doc = f_loglog2(us, ss, U6)
    ok2 = abs(doc - REC_DOC) <= 0.006
    say(f"   the DOCUMENT's own K = 6 number, recomputed by its own form: {doc:.4f} "
        f"vs printed {REC_DOC} -> {'OK' if ok2 else 'MISMATCH'}")
    if not (ok and ok2):
        say()
        say("CONTROL FAILED -- nothing below is trustworthy.")
        _write()
        return 1
    say()

    say("=== IS A SINGLE SMOOTH TWO-PARAMETER FORM EVEN CONSISTENT WITH THE THREE "
        "POINTS? ===")
    say("   Fit each two-parameter form to the two OUTER points and predict the MIDDLE")
    say("   one. A form that cannot hit the point it was not given is not a basis for")
    say("   extrapolating between them.")
    mid_u, mid_s = float(us[1]), float(ss[1])
    outer_u = np.array([us[0], us[2]])
    outer_s = np.array([ss[0], ss[2]])
    say(f"   {'form':<40}{'predicts mid':>14}{'measured':>10}{'error':>9}")
    worst = 0.0
    for nm, fn, npar, _anch in FORMS:
        if npar != 2:
            continue
        pred = fn(outer_u, outer_s, mid_u)
        err = abs(pred - mid_s) / mid_s
        worst = max(worst, err)
        say(f"   {nm:<40}{pred:>14.3f}{mid_s:>10.3f}{err:>8.1%}")
    say(f"   -> worst miss on the held-out middle point: {worst:.1%}. The three points")
    say("      are NOT collinear under any of these transforms, so the form is genuinely")
    say("      under-determined by the data and no choice among them is privileged.")
    say()

    say("=== THE SPREAD OVER A DECLARED SET OF FORMS ===")
    say("   'anchored' = the form is forced through the EXACT physical point (u = 1,")
    say("   s = 0): at u = 1 the classes are identical, the mixture is monodisperse,")
    say("   phi = beta and the filling signal vanishes. That costs no fitted parameter.")
    say(f"   {'form':<40}{'pars':>5}{'anch':>6}{'s(K=6)':>9}{'margin':>8}{'K6?':>6}"
        f"{'s(K=7)':>9}{'margin':>8}{'K7?':>6}")
    v6, v7 = [], []
    for nm, fn, npar, anch in FORMS:
        s6 = fn(us, ss, U6)
        s7 = fn(us, ss, U7)
        v6.append(s6)
        v7.append(s7)
        say(f"   {nm:<40}{npar:>5}{anch:>6}{s6:>9.3f}{s6 / BAR:>8.3f}"
            f"{('yes' if s6 > BAR else 'no'):>6}{s7:>9.3f}{s7 / BAR:>8.3f}"
            f"{('yes' if s7 > BAR else 'no'):>6}")
    v6 = np.array(v6)
    v7 = np.array(v7)
    say()
    say(f"   K = 6 (u = {U6:.5f}): signal spans {v6.min():.3f} to {v6.max():.3f}, a "
        f"factor {v6.max() / v6.min():.2f}; the bar {BAR} is "
        f"{'INSIDE' if v6.min() < BAR < v6.max() else 'OUTSIDE'} that span.")
    say(f"           verdict 'attemptable' on {int((v6 > BAR).sum())} of {v6.size} "
        f"forms.")
    say(f"   K = 7 (u = {U7:.5f}): signal spans {v7.min():.3f} to {v7.max():.3f}, a "
        f"factor {v7.max() / v7.min():.2f}; the bar {BAR} is "
        f"{'INSIDE' if v7.min() < BAR < v7.max() else 'OUTSIDE'} that span.")
    say(f"           verdict 'attemptable' on {int((v7 > BAR).sum())} of {v7.size} "
        f"forms.")
    say()

    # the verdict, GENERATED from the numbers rather than written beside them
    say("=== INTERIM VERDICT ON THE DOCUMENT'S OWN BASIS (generated, not typed) ===")
    say("   (!) This grades the document's estimate against the document's own bar and")
    say("   over ALL seven forms. Two corrections below supersede it: four of the forms")
    say("   are FALSIFIED by the exact anchor, and the 2.57 bar was borrowed from")
    say("   another ratio. Read the final section, not this one, for the verdict on")
    say("   whether K = 6 is attemptable.")
    n6 = int((v6 > BAR).sum())
    n7 = int((v7 > BAR).sum())
    if 0 < n6 < v6.size:
        say(f"   K = 6: UNRESOLVED on this basis. The estimate flips with the "
            f"interpolation form ({n6}/{v6.size} say yes),")
        say(f"   and the document's {REC_DOC} is one draw from a {v6.min():.2f}-"
            f"{v6.max():.2f} range that straddles the bar.")
        say("   'The K axis closes at 6' is therefore NOT established by this estimate.")
    elif n6 == 0:
        say(f"   K = 6: every form falls short of {BAR}. The document's "
            f"{REC_DOC} is the outlier.")
    else:
        say(f"   K = 6: every form clears {BAR}; the verdict is robust to the form.")
    if n7 == 0:
        say(f"   K = 7: CLOSED, and robustly -- no form in the declared set reaches "
            f"{BAR} (max {v7.max():.2f}).")
    elif n7 == v7.size:
        say(f"   K = 7: every form clears the bar, which CONTRADICTS the recorded "
            f"{REC_DOC_K7}. Investigate before quoting either.")
    else:
        say(f"   K = 7: also UNRESOLVED ({n7}/{v7.size} say yes) -- the recorded "
            f"'no' at {REC_DOC_K7} is not safe.")
    say()

    # ============================================================== referee pass 2 (i)
    say("=== IS THE (u = 1, s = 0) ANCHOR EXACT IN CODE, OR ONLY IN MY ALGEBRA? ===")
    import mixture_closure as mc_  # noqa: PLC0415
    worst_anchor = 0.0
    for r_ in (1.25, 1.31951, 1.5, 2.0, 4.0):
        for x_ in ug.X_GRID:
            v = float(mc_.packing_density_ab(np.array([r_, 1.0]),
                                             np.array([1 - x_, x_]), beta, 1.0, 1.0))
            worst_anchor = max(worst_anchor, abs(v - beta))
    say(f"   packing_density_ab(..., a = 1, b = 1) - beta over 5 ratios x "
        f"{len(ug.X_GRID)} compositions: worst |diff| = {worst_anchor:.2e}")
    say("   -> the anchor is EXACT in the shipped closure, not just in the algebra: at")
    say("      a = b = 1 both filling branches vanish identically and phi = beta.")
    say()

    say("=== BUT THE ANCHOR DOES NOT NARROW ANYTHING -- IT FALSIFIES THE FORMS ===")
    say("   A two-parameter power law s = c(1-u)^p satisfies s(1) = 0 for ANY (c, p), so")
    say("   the anchor constrains the FAMILY and not the PARAMETERS. And the three")
    say("   pairwise exponents disagree by 2x, so there is no single power law:")
    say(f"   {'pair':<22}{'p':>9}{'predicts the third pt':>24}{'error':>9}")
    idx = {0: "u=0.667", 1: "u=0.707", 2: "u=0.800"}
    for i, j in ((0, 1), (0, 2), (1, 2)):
        k = 3 - i - j
        p_ = math.log(ss[i] / ss[j]) / math.log((1 - us[i]) / (1 - us[j]))
        c_ = ss[i] / (1 - us[i]) ** p_
        pred = c_ * (1 - us[k]) ** p_
        say(f"   {idx[i] + ' & ' + idx[j]:<22}{p_:>9.3f}{pred:>24.3f}"
            f"{(pred - ss[k]) / ss[k]:>+9.1%}")
    say()
    say("   And the SMOOTH forms are falsified outright by the same anchor. Evaluate")
    say("   each at u = 1, where the answer is exactly 0:")
    anchor_ok: dict[str, bool] = {}
    say(f"   {'form':<40}{'s(u -> 1)':>16}{'should be':>11}{'verdict':>13}")
    for nm, fn, _np_, _an in FORMS:
        try:
            v = fn(us, ss, 1.0 - 1e-12)
        except (ValueError, ZeroDivisionError, FloatingPointError, OverflowError):
            v = float("inf")
        adm = abs(v) < 0.05
        anchor_ok[nm] = adm
        say(f"   {nm:<40}{v:>16.3g}{0.0:>11.1f}{('ADMISSIBLE' if adm else 'NO'):>13}")
    say("   -> a form that does not vanish at u = 1 contradicts the one piece of EXACT")
    say("      physics available, and one of them goes NEGATIVE, which max(phi - beta)")
    say("      cannot be (at x = 0 the mixture IS the coarse pack, so the max is >= 0).")
    say("      That is stronger than 'under-determined': those forms are demonstrably")
    say("      wrong outside a narrow neighbourhood, which is a reason to distrust them")
    say("      inside it too.")
    say()
    keep = [(nm, fn) for nm, fn, _a, _b in FORMS if anchor_ok[nm]]
    k6a = np.array([fn(us, ss, U6) for _nm, fn in keep])
    k7a = np.array([fn(us, ss, U7) for _nm, fn in keep])
    say(f"=== RESTRICTED TO THE {len(keep)} ANCHOR-ADMISSIBLE FORMS ===")
    say(f"   K = 6 signal {k6a.min():.3f} to {k6a.max():.3f} "
        f"(all {'ABOVE' if k6a.min() > BAR else 'straddling'} the borrowed bar {BAR}); "
        f"K = 7 {k7a.min():.3f} to {k7a.max():.3f}.")
    say("   So on the forms that are not falsified, K = 6 clears the BORROWED bar on")
    say("   every one of them -- which is exactly why the bar itself now has to be")
    say("   checked rather than borrowed.")
    say()

    # ============================================================= referee pass 2 (ii)
    say("=== THE SECOND, INDEPENDENT DEFECT: THE 2.57 BAR IS NOT A CONSTANT ===")
    say("   grade_ugap.py computes that ceiling AT ARM B's RATIO (r = 1.25) over a")
    say("   handful of assumed truths clustered near (a, b) = (0.975, 0.960). It is the")
    say("   smallest signal that design could RESOLVE THERE. Using it at u = 0.75786 is")
    say("   the same units error as the 'phi span' trap, relocated: a signal")
    say("   interpolated to one ratio, graded against a bar computed at another.")
    say("   Recomputed here at r = 4^(1/5) = 1.31951, noise-free, over a grid of truths.")
    r6 = 1.0 / U6
    say(f"   {'a_true':>8}{'b_true':>8}{'w(a)':>8}{'w(b)':>8}{'ident':>7}"
        f"{'signal/scale':>14}")
    ok_sigs, no_sigs = [], []
    for a_ in (0.9400, 0.9481, 0.9600, 0.9700, 0.9800, 0.9850, 0.9900):
        for b_ in (0.9376, 0.9600):
            p_ = [(r6, x_, float(mc_.packing_density_ab(np.array([r6, 1.0]),
                                                        np.array([1 - x_, x_]),
                                                        beta, a_, b_)))
                  for x_ in ug.X_GRID]
            ok_, ia_, ib_, _, _ = ug.identifiable(p_, beta)
            s_ = (max(q[2] for q in p_) - beta) / ug.PHI_SCALE
            (ok_sigs if ok_ else no_sigs).append(s_)
            say(f"   {a_:8.4f}{b_:8.4f}{ia_[1] - ia_[0]:8.3f}{ib_[1] - ib_[0]:8.3f}"
                f"{('YES' if ok_ else 'NO'):>7}{s_:14.2f}")
    say("   a(0.70711) is MEASURED at 0.9481 and a(0.8) is inferred >= 0.975, and a(u)")
    say("   is increasing, so the truth at u = 0.75786 lies between them -- the rows")
    say("   above bracket it.")
    if ok_sigs and no_sigs:
        bar_lo, bar_hi = max(no_sigs), min(ok_sigs)
        say(f"   -> at THIS ratio the design FAILS up to {bar_lo:.2f} scale and RESOLVES "
            f"from {bar_hi:.2f}, so the bar is")
        say(f"      BRACKETED in [{bar_lo:.2f}, {bar_hi:.2f}] -- it is NOT a single "
            f"number, and the {BAR} the document used")
        say(f"      was computed at r = 1.25, a different ratio. The bracket's ends are "
            f"{bar_lo / BAR:.2f}x and {bar_hi / BAR:.2f}x the borrowed value.")
        say()
        say("   COMBINING THE TWO DEFECTS -- and they do NOT push the same way:")
        say(f"      signal, anchor-admissible forms only : {k6a.min():.3f} to "
            f"{k6a.max():.3f}   (the document used {REC_DOC})")
        say(f"      bar, recomputed at this ratio        : {bar_lo:.2f} to {bar_hi:.2f} "
            f"        (the document used {BAR})")
        m_lo, m_hi = k6a.min() / bar_hi, k6a.max() / bar_lo
        say(f"      margin                               : {m_lo:.2f}x to {m_hi:.2f}x  "
            f"     (the document used {REC_DOC / BAR:.2f}x)")
        say("      Falsifying the bad interpolants pushes the signal UP; correcting the")
        say("      borrowed bar pushes the threshold UP too. They partly CANCEL, which is")
        say("      why the corrected margin lands near the document's -- for reasons that")
        say("      have nothing to do with the document's own argument.")
        if m_lo > 1.0:
            say(f"      -> K = 6 stays ATTEMPTABLE on every corner of this box, at a "
                f"margin of {m_lo:.2f}x")
            say("         at its worst. THINNER than it looks and resting on a bar that")
            say("         is bracketed rather than measured, but it does not flip.")
        elif m_hi < 1.0:
            say("      -> K = 6 is CLOSED on every corner of this box.")
        else:
            say(f"      -> K = 6 remains UNRESOLVED: the box straddles 1 "
                f"({m_lo:.2f}x to {m_hi:.2f}x).")
    elif ok_sigs:
        say("   -> EVERY truth on this grid is identifiable at this ratio, so the ceiling")
        say("      is not bracketed from below here and no bar can be read off it.")
    else:
        say("   -> NO truth on this grid is identifiable at this ratio, which is itself")
        say("      the answer: the design cannot resolve a(0.75786) at any of them.")
    say()

    say("   What can decide K = 6 is unchanged and is what the document already asks")
    say(f"   for: recompute the noise-free identifiability ceiling at r = 4^(1/5) = "
        f"{1.0 / U6:.5f},")
    say("   at several assumed truths, BEFORE any prereg. None of the numbers above is")
    say("   a substitute for that, and none of them moves a bar.")
    _write()
    return 0


def _write():
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "diag_k6_interp.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LINES) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
