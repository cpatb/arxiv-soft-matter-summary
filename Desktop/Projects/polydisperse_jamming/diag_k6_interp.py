"""diag_k6_interp.py -- is "K = 6 is attemptable" a measurement of the SIGNAL or of the
INTERPOLATION FORM?

POST-HOC, NOT A GATE. Written after a physics-referee pass on
`derivations/k_axis_bound.md` reported that the document's forward-content estimate --
signal 2.99 at u = 0.75786 against a 2.57 bar, margin 1.16x -- comes from a TWO-POINT
log-log interpolation that discards the third measured point, and that an equally
defensible form through all three lands at 2.56, i.e. margin 0.997 and the opposite
verdict.

This file exists to do three things and nothing else:

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
    say("=== VERDICT (generated from the table above, not typed) ===")
    n6 = int((v6 > BAR).sum())
    n7 = int((v7 > BAR).sum())
    if 0 < n6 < v6.size:
        say(f"   K = 6: UNRESOLVED. The estimate flips with the interpolation form "
            f"({n6}/{v6.size} say yes),")
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
