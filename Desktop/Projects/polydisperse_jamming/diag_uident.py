"""diag_uident.py -- WHY a(u) is unidentifiable at u = 0.8, and what (if anything) fixes it.

`diag_ugap.py` answered the first question: the K = 5 prediction IS sensitive to the
unmeasured region (worst joint excursion 9.59 scale over the monotone-admissible range,
4.7x TAB's own in-sample rms of 2.02, and it carries R(5) - R(4) across zero on 12 of 12
rows). So the region matters and option (b) -- "irrelevant because phi is insensitive" --
is dead.

That leaves (a) measurable with a better-designed block, or (c) genuinely unidentifiable.
This script decides which, BEFORE any cell runs, because the answer is a property of the
MODEL and the DESIGN, not of the data:

  * the profile bar is delta chi2 = 1 PER POINT, so the interval's width is set by the
    CURVATURE of chi2 in the parameter and not by how many points there are. Adding
    points at the same sensitivity adds equally to chi2 and to the bar. That is testable
    by construction rather than by argument (section 3), and if it holds then "more seeds
    and more compositions" cannot shrink anything and the design has to move somewhere
    else;
  * the linear packing model is a MIN over branches, and one of those branches is
    phi = beta identically (section 2). A parameter that the binding branch does not see
    is not measured, however many cells are run -- the same structure that made the
    ternary block's `a(1/4)` degenerate.

Nothing here reads a new number off disk. Every phi in sections 3-5 is SYNTHETIC, drawn
noise-free from the shipped closure, so what is computed is the IDENTIFIABILITY CEILING:
the narrowest interval a perfect experiment at that design could return. A design that
fails noise-free cannot be rescued by seeds.

Run:  python diag_uident.py      (writes out/diag_uident.txt)
"""

from __future__ import annotations

import os
import sys

import numpy as np

import grade_mixture_ext as gx
import mixture_closure as mc
import mixture_model as mm

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = 2.1e-3
MU = 0.3
IDENT_W = 0.25          # the identifiability rule's width bar, verbatim from
                        # grade_mixture_ext.tabulated_closure. NOT MOVABLE.
MIN_PHASE_COUNT = 50
N_PROD = 1500


def branch_tots(x, beta: float, a: float, b: float):
    """The two branch totals of the linear packing model for ONE binary.

    Written out by hand from mixture_closure.packing_density_ab so that it is an ANCHOR
    OUTSIDE the code under test rather than a second call to it:
      coarse-class branch   tot = (1 - a) * y_S          -- sees a only
      fine-class branch     tot = (1 - beta)(1 - b) y_L  -- sees b only
    (the fine branch's 1 - beta + b*beta*(1 - 1/beta) simplifies to (1 - beta)(1 - b)).
    phi = beta / (1 - min(tot)); larger tot means larger phi, so the MIN over branches
    picks the SMALLER tot.
    """
    y_S, y_L = x, 1.0 - x
    return (1.0 - a) * y_S, (1.0 - beta) * (1.0 - b) * y_L


def n_min_frac(r: float, x: float) -> float:
    """Minority number fraction for a binary of size ratio r at fines volume x.

    n_i proportional to y_i / d_i^3 with d_S = 1, d_L = r.
    """
    wS, wL = x, (1.0 - x) / r**3
    return min(wS, wL) / (wS + wL)


def synth(r: float, xs, beta: float, a: float, b: float):
    return [(r, float(x), float(mc.packing_density_ab(np.array([r, 1.0]),
                                                      np.array([1.0 - x, x]),
                                                      beta, a, b)))
            for x in xs]


def phi_fast(a, b, x, beta: float):
    """Vectorised binary closure, from the branch algebra of `branch_tots`.

    This is a SECOND implementation, not a wrapper: `check_fast` asserts it reproduces
    `mixture_closure.packing_density_ab` before any number computed with it is printed.
    It exists only because the shipped profile grid is 201 x 201 re-minimisations per
    parameter and a pure-Python scan of it does not finish.
    """
    tot = np.minimum((1.0 - a) * x, (1.0 - beta) * (1.0 - b) * (1.0 - x))
    return beta / (1.0 - tot)


def check_fast(beta: float, rng) -> float:
    """Worst |phi_fast - packing_density_ab| over random (a, b, x). The anchor."""
    worst = 0.0
    for _ in range(300):
        a, b, x = rng.uniform(1e-4, 1 - 1e-4, 3)
        r = float(rng.uniform(1.05, 6.0))
        ref = mc.packing_density_ab(np.array([r, 1.0]), np.array([1.0 - x, x]),
                                    beta, float(a), float(b))
        worst = max(worst, abs(float(phi_fast(a, b, x, beta)) - ref))
    return worst


def profile_fast(pts, beta: float, which: int, ngrid: int = 201):
    """`mixture_closure.profile_interval`, same grid, same delta chi2 = n bar, vectorised.

    Reproduces the shipped routine's interval; `main` checks that on real data before
    using it, because a faster reimplementation of the thing under test is exactly the
    kind of code that reports a bug as a finding.
    """
    x = np.array([p[1] for p in pts], dtype=float)
    y = np.array([p[2] for p in pts], dtype=float)
    grid = np.linspace(1e-4, 1.0 - 1e-4, ngrid)
    A = grid[:, None, None]
    B = grid[None, :, None]
    chi = (((phi_fast(A, B, x[None, None, :], beta) - y[None, None, :]) / SCALE) ** 2
           ).sum(axis=2)
    prof = chi.min(axis=1) if which == 0 else chi.min(axis=0)
    lo = float(prof.min())
    ok = grid[prof <= lo + len(pts)]
    return (float(ok.min()), float(ok.max())) if ok.size else (0.0, 1.0)


def widths(pts, beta: float):
    la, ha = profile_fast(pts, beta, 0)
    lb, hb = profile_fast(pts, beta, 1)
    return (ha - la), (hb - lb), (la, ha), (lb, hb)


def main() -> int:
    lines: list[str] = []

    def say(s: str = "") -> None:
        lines.append(s)
        print(s, flush=True)

    meas = mm.measured(mu_p=MU)
    beta = meas[("lognormal", 0.0)][0]
    allpts = gx.ext_points(beta)
    _tab, fa0, fb0 = gx.tabulated_closure(beta, allpts)
    ratios = sorted({p[0] for p in allpts})

    say("=== diag_uident.py -- is a(u) MEASURABLE near u = 0.7-0.8, or unidentifiable? ===")
    say(f"   beta = {beta:.6f} (measured, mu_p {MU}); SCALE = {SCALE:.1e}")

    # ---- ANCHORS, before any number below is believed -------------------------------
    w = check_fast(beta, np.random.default_rng(20260916))
    say("   anchor 1 -- the vectorised kernel vs mixture_closure.packing_density_ab over")
    say(f"      300 random (r, a, b, x): worst |diff| = {w:.3e}")
    if w > 1e-12:
        say("      -> ANCHOR FAILED. Nothing below is printed.")
        return 2
    say("   anchor 2 -- the vectorised profile vs mixture_closure.profile_interval on the")
    say("      two ratios that matter, on REAL data:")
    for rr in (1.25, 1.5):
        sub0 = [p for p in allpts if abs(p[0] - rr) < 1e-9 and 0.0 < p[1] < 1.0]
        for wh, nm in ((0, "a"), (1, "b")):
            g_ = profile_fast(sub0, beta, wh)
            s_ = mc.profile_interval(sub0, beta, wh)
            say(f"      r = {rr:4.2f} {nm}: fast [{g_[0]:.4f}, {g_[1]:.4f}]  "
                f"shipped [{s_[0]:.4f}, {s_[1]:.4f}]  "
                f"{'MATCH' if max(abs(g_[0]-s_[0]), abs(g_[1]-s_[1])) < 1e-12 else 'DIFFER'}")
            if max(abs(g_[0] - s_[0]), abs(g_[1] - s_[1])) >= 1e-12:
                say("      -> ANCHOR FAILED. Nothing below is printed.")
                return 2
    say()
    say("   identifiability rule, quoted verbatim from grade_mixture_ext."
        "tabulated_closure:")
    say(f"      both profile intervals narrower than {IDENT_W} AND containing their own "
        f"best fit.")
    say()

    # ------------------------------------------------------------------ section 1
    say("=== 1. the measured state: what each ratio's binary block returns ===")
    say(f"   {'r':>6}{'u':>8}{'n':>4}{'phi span':>10}{'/scale':>8}"
        f"{'a':>8}{'w(a)':>8}{'b':>8}{'w(b)':>8}  ident?")
    state = {}
    for r in ratios:
        sub = [p for p in allpts if abs(p[0] - r) < 1e-9 and 0.0 < p[1] < 1.0]
        a_, b_, _ = mc.fit_ab_direct(sub, beta)
        wa, wb, ia, ib = widths(sub, beta)
        span = max(p[2] for p in sub) - min(p[2] for p in sub)
        ok = (wa < IDENT_W and wb < IDENT_W
              and ia[0] - 1e-9 <= a_ <= ia[1] + 1e-9
              and ib[0] - 1e-9 <= b_ <= ib[1] + 1e-9)
        state[r] = (a_, b_, wa, wb, span, ok)
        say(f"   {r:6.2f}{1.0/r:8.4f}{len(sub):4d}{span:10.4f}{span/SCALE:8.1f}"
            f"{a_:8.4f}{wa:8.3f}{b_:8.4f}{wb:8.3f}  {'YES' if ok else 'NO'}")
    say("   -> u = 0.8 (r = 1.25) is the ONLY ratio that fails, and it fails on BOTH")
    say("      parameters with intervals spanning the entire [0, 1].")
    say()

    # ------------------------------------------------------------------ section 2
    say("=== 2. WHY: the model has an escape branch at b = 1 where phi = beta exactly ===")
    say("   The closure is a MIN over branches. The fine-class branch carries")
    say("   tot = (1 - beta)(1 - b) y_L, so at b = 1 it gives tot = 0 and phi = beta for")
    say("   EVERY composition and EVERY a. If the measured phi never rises far enough")
    say("   above beta, that constant model is inside the bar and a is free.")
    say(f"   {'r':>6}{'u':>8}{'max(phi-beta)':>15}{'/scale':>9}"
        f"{'chi2 at b=1':>13}{'bar = n':>9}{'escape?':>9}")
    for r in ratios:
        sub = [p for p in allpts if abs(p[0] - r) < 1e-9 and 0.0 < p[1] < 1.0]
        # the whole model collapses to the constant beta at b = 1, whatever a is
        chi_esc = min(
            sum(((mc.packing_density_ab(np.array([r, 1.0]), np.array([1.0 - x, x]),
                                        beta, av, 1.0) - phi) / SCALE) ** 2
                for _, x, phi in sub)
            for av in np.linspace(1e-4, 1 - 1e-4, 51))
        dmax = max(p[2] for p in sub) - beta
        say(f"   {r:6.2f}{1.0/r:8.4f}{dmax:15.4f}{dmax/SCALE:9.1f}"
            f"{chi_esc:13.2f}{len(sub):9d}"
            f"{'YES' if chi_esc <= len(sub) else 'no':>9}")
    say("   -> the escape is open at exactly the ratio that fails. It is NOT a noise")
    say("      problem and NOT a sample-size problem: it is the filling SIGNAL, phi -")
    say("      beta, falling to the size of the bar.")
    say()

    # ------------------------------------------------------------------ section 3
    say("=== 3. does adding points or seeds narrow the interval? NO, and by construction ===")
    say("   The bar is delta chi2 = 1 PER POINT (mixture_closure.profile_interval:")
    say("   `prof <= lo + len(pts)`). Duplicating a point set doubles chi2 AND the bar,")
    say("   so the interval cannot move. Demonstrated rather than argued:")
    r_demo = 1.5
    sub = [p for p in allpts if abs(p[0] - r_demo) < 1e-9 and 0.0 < p[1] < 1.0]
    wa1, wb1, _, _ = widths(sub, beta)
    wa2, wb2, _, _ = widths(sub + sub, beta)
    wa4, wb4, _, _ = widths(sub * 4, beta)
    say(f"   r = {r_demo}: n = {len(sub):2d}  w(a) = {wa1:.4f}  w(b) = {wb1:.4f}")
    say(f"   r = {r_demo}: n = {2*len(sub):2d}  w(a) = {wa2:.4f}  w(b) = {wb2:.4f}   "
        f"(same points, twice)")
    say(f"   r = {r_demo}: n = {4*len(sub):2d}  w(a) = {wa4:.4f}  w(b) = {wb4:.4f}   "
        f"(same points, four times)")
    same = abs(wa1 - wa2) < 1e-12 and abs(wa1 - wa4) < 1e-12
    say(f"   -> widths IDENTICAL: {same}. **More compositions at the same sensitivity")
    say("      buy nothing under this rule, and more SEEDS buy less than that -- seeds")
    say("      move the seed-mean phi, not the curvature, and the bar is in SCALE units,")
    say("      not in units of the seed error.** Any design that hopes to shrink these")
    say("      intervals has to move the SENSITIVITY, i.e. WHERE in composition it")
    say("      samples, not HOW MUCH.")
    say()

    # ------------------------------------------------------------------ section 4
    say("=== 4. where the a-signal actually lives: the branch crossover ===")
    say("   The coarse branch (which is the only one that sees a) binds only while")
    say("   (1 - a) y_S < (1 - beta)(1 - b) y_L. Above that composition the fine branch")
    say("   takes over and a stops entering the model at all. Using the SHIPPED")
    say("   interpolated a(u), b(u) -- the blind prediction, not a measurement:")
    say(f"   {'r':>6}{'u':>8}{'a(u)':>8}{'b(u)':>8}{'x*cross':>9}"
        f"{'max a-signal':>13}{'at x':>7}{'n_min@N1500':>12}")
    cand_r = [1.25, 4 / 3, 1.4, 1.4142135623730951, 1.5, 1.6, 2.0]
    xg = np.linspace(0.01, 0.99, 981)
    plan = {}
    for r in cand_r:
        u = 1.0 / r
        av, bv = float(fa0(u)), float(fb0(u))
        tc, tf = branch_tots(xg, beta, av, bv)
        coarse_binds = np.asarray(tc < tf)
        xstar = float(xg[coarse_binds][-1]) if coarse_binds.any() else float("nan")
        # a-signal: how far phi sits above the b=1 constant, where the coarse branch binds
        sig = np.where(coarse_binds, beta / (1 - tc) - beta, 0.0) / SCALE
        i = int(np.argmax(sig))
        plan[r] = (u, av, bv, xstar, float(sig[i]), float(xg[i]))
        say(f"   {r:6.3f}{u:8.4f}{av:8.4f}{bv:8.4f}{xstar:9.3f}"
            f"{sig[i]:13.2f}{xg[i]:7.3f}{n_min_frac(r, xg[i])*N_PROD:12.0f}")
    say("   -> the a-signal PEAKS at the crossover and is zero beyond it. The existing")
    say("      block sampled x in [0.05, 0.60] at every ratio -- at r = 1.25 that is")
    say("      almost entirely BELOW the peak, where the signal is proportional to x and")
    say("      therefore smallest.")
    say()

    # ------------------------------------------------------------------ section 5
    say("=== 5. THE IDENTIFIABILITY CEILING: noise-free synthetic data, per design ===")
    say("   phi is generated FROM the shipped closure with no noise, so chi2_min = 0 and")
    say("   the width returned is the narrowest any experiment at that design could get.")
    say("   A design that fails HERE cannot be rescued by seeds, and one that passes here")
    say("   is an upper bound on what the real block will do, not a promise.")
    designs = {
        "OLD  x in [0.05, 0.60] (the block that ran)":
            (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.60),
        "MID  x in [0.30, 0.75]": (0.30, 0.40, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75),
        "PEAK x straddling the crossover": None,          # filled per ratio
        "WIDE x in [0.05, 0.95]":
            (0.05, 0.15, 0.30, 0.45, 0.60, 0.70, 0.80, 0.90, 0.95),
    }
    say("   n_req / aff = compositions requested / affordable at MIN_PHASE_COUNT; "
        "n_c / n_f = points")
    say("   on the coarse (a-sensitive) and fine (b-sensitive) branch; esc(b=1), esc(a=1) "
        "= chi2 of")
    say("   the two escape models that collapse the closure to the CONSTANT phi = beta, "
        "against a bar of aff.")
    say(f"   {'design':<44}{'r':>6}{'u':>8}{'w(a)':>8}{'w(b)':>8}{'ident?':>8}"
        f"{'n_req':>6}{'aff':>5}{'n_c':>5}{'n_f':>5}{'esc(b=1)':>11}{'esc(a=1)':>11}")
    ceiling: dict[tuple[str, float], tuple[float, float, bool]] = {}
    for name, xs in designs.items():
        for r in cand_r:
            u, av, bv, xstar, _, _ = plan[r]
            if xs is None:
                # eight points straddling the crossover: four below (a-sensitive) and
                # four above (b-sensitive), so BOTH parameters get their own half.
                lo = np.linspace(max(0.05, xstar - 0.30), xstar - 0.02, 4)
                hi = np.linspace(min(0.97, xstar + 0.02), min(0.97, xstar + 0.30), 4)
                use = tuple(np.concatenate([lo, hi]))
            else:
                use = xs
            n_before = len(use)
            use = tuple(x for x in use
                        if n_min_frac(r, x) * N_PROD >= MIN_PHASE_COUNT)
            if len(use) < 4:
                say(f"   {name:<44}{r:6.3f}{u:8.4f}   REFUSED: {len(use)} of {n_before} "
                    f"compositions affordable at N = {N_PROD}")
                continue
            pts = synth(r, use, beta, av, bv)
            wa, wb, ia, ib = widths(pts, beta)
            a_, b_, _ = mc.fit_ab_direct(pts, beta)
            ok = (wa < IDENT_W and wb < IDENT_W
                  and ia[0] - 1e-9 <= a_ <= ia[1] + 1e-9
                  and ib[0] - 1e-9 <= b_ <= ib[1] + 1e-9)
            # how many points each branch actually gets, and whether either escape
            # (a = 1 or b = 1, both of which collapse the closure to the CONSTANT
            # phi = beta) is open. An open escape is what makes a marginal span [0, 1].
            tc, tf = branch_tots(np.asarray(use, float), beta, av, bv)
            n_coarse = int(np.sum(np.asarray(tc < tf)))
            grid_ = np.linspace(1e-4, 1 - 1e-4, 201)
            xv = np.array([q[1] for q in pts])
            yv = np.array([q[2] for q in pts])
            esc_b1 = float(min((((phi_fast(g, 1.0, xv, beta) - yv) / SCALE) ** 2).sum()
                               for g in grid_))
            esc_a1 = float(min((((phi_fast(1.0, g, xv, beta) - yv) / SCALE) ** 2).sum()
                               for g in grid_))
            ceiling[(name, r)] = (wa, wb, ok)
            say(f"   {name:<44}{r:6.3f}{u:8.4f}{wa:8.3f}{wb:8.3f}"
                f"{'YES' if ok else 'NO':>8}{n_before:6d}{len(use):5d}"
                f"{n_coarse:5d}{len(use) - n_coarse:5d}{esc_b1:11.1f}{esc_a1:11.1f}")
        say()

    # a CONTROL the ceiling must pass: the OLD design at r = 2 and 3 is known
    # identifiable on real data, so if it failed here the ceiling machinery is broken.
    say("   CONTROL -- the ceiling must reproduce what the real data already did:")
    for r in (2.0, 1.5):
        got = ceiling.get(("OLD  x in [0.05, 0.60] (the block that ran)", r))
        was = state.get(r)
        if got and was:
            say(f"      r = {r}: ceiling ident {got[2]}, real data ident {was[5]} "
                f"(w(a) ceiling {got[0]:.3f} vs real {was[2]:.3f})")
    say("      and the failing case:")
    got = ceiling.get(("OLD  x in [0.05, 0.60] (the block that ran)", 1.25))
    was = state.get(1.25)
    if got and was:
        say(f"      r = 1.25: ceiling ident {got[2]}, real data ident {was[5]} "
            f"(w(a) ceiling {got[0]:.3f} vs real {was[2]:.3f})")
    say()

    # ------------------------------------------------------------------ section 6
    say("=== 6. the ratio the K-ladder actually needs is NOT 1.25 ===")
    say("   K = 5's four pairs sit at u = 0.70711, i.e. r = 4^(1/4) = 1.41421, NOT at")
    say("   u = 0.8. The ceiling above is identifiable at r = 1.41421 on the OLD design")
    say("   and fails only at r = 1.25, where the filling signal is weakest. So the")
    say("   design question is not 'how do we rescue u = 0.8' but 'does the standard")
    say("   block work at u = 0.7071' -- and the ceiling says yes, AT the shipped a(u).")
    say("   Since the shipped a(u*) is the extrapolation the block exists to replace,")
    say("   the honest next question is how far it can be wrong before the answer flips.")
    say()
    say("   ROBUSTNESS: the ceiling at r = 1.41421 as the TRUTH moves. Rows are the")
    say("   assumed true a(u*), columns the assumed true b(u*). Entry is 'Y' when BOTH")
    say(f"   widths land under {IDENT_W}, else the WIDER of the two.")
    r_star = 4.0**0.25
    a_grid = (0.90, 0.92, 0.9398, 0.96, 0.98, 0.99, 0.995)
    b_grid = (0.86, 0.90, 0.9177, 0.94, 0.96, 0.98, 0.99)
    old_x = (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.60)
    wide_x = (0.05, 0.15, 0.30, 0.45, 0.60, 0.70, 0.80, 0.90, 0.95)
    crit: dict[str, float] = {}
    for nm, xs_ in (("OLD grid", old_x), ("WIDE grid", wide_x)):
        aff = tuple(x for x in xs_
                    if n_min_frac(r_star, x) * N_PROD >= MIN_PHASE_COUNT)
        say(f"   --- {nm}: {len(aff)} of {len(xs_)} compositions affordable at "
            f"N = {N_PROD}: {[round(x, 2) for x in aff]}")
        say("        a\\b " + "".join(f"{bv_:>9.3f}" for bv_ in b_grid))
        worst_ok = 0.0
        for av_ in a_grid:
            row = []
            for bv_ in b_grid:
                pts_ = synth(r_star, aff, beta, av_, bv_)
                wa_, wb_, ia_, ib_ = widths(pts_, beta)
                af_, bf_, _ = mc.fit_ab_direct(pts_, beta)
                ok_ = (wa_ < IDENT_W and wb_ < IDENT_W
                       and ia_[0] - 1e-9 <= af_ <= ia_[1] + 1e-9
                       and ib_[0] - 1e-9 <= bf_ <= ib_[1] + 1e-9)
                if ok_:
                    worst_ok = max(worst_ok, av_)
                row.append("        Y" if ok_ else f"{max(wa_, wb_):9.2f}")
            say(f"      {av_:6.3f}" + "".join(row))
        crit[nm] = worst_ok
        say(f"   -> largest true a(u*) still identifiable anywhere on this grid: "
            f"{worst_ok:.3f}")
        say()

    # ------------------------------------------------------------------ section 7
    say("=== 7. what the r = 1.25 data DOES constrain -- including the COMBINATION ===")
    say("   If a(u) is unidentifiable there, the fallback is to constrain the combination")
    say("   instead. That needs the JOINT region, not the two marginals, because two wide")
    say("   marginals can still enclose a narrow diagonal.")
    sub125 = [p for p in allpts if abs(p[0] - 1.25) < 1e-9 and 0.0 < p[1] < 1.0]
    xv125 = np.array([q[1] for q in sub125])
    yv125 = np.array([q[2] for q in sub125])
    g201 = np.linspace(1e-4, 1 - 1e-4, 201)
    chi2d = (((phi_fast(g201[:, None, None], g201[None, :, None],
                        xv125[None, None, :], beta)
               - yv125[None, None, :]) / SCALE) ** 2).sum(axis=2)
    lo2 = float(chi2d.min())
    inside = chi2d <= lo2 + len(sub125)
    ai, bi = np.where(inside)
    say(f"   n = {len(sub125)} compositions; chi2_min = {lo2:.2f}, "
        f"bar = chi2_min + n = {lo2 + len(sub125):.2f}")
    say(f"   admissible area: {inside.mean():.1%} of the unit square "
        f"({int(inside.sum())} of {inside.size} grid nodes)")
    say(f"   a range inside the joint region: [{g201[ai].min():.4f}, {g201[ai].max():.4f}]")
    say(f"   b range inside the joint region: [{g201[bi].min():.4f}, {g201[bi].max():.4f}]")
    const_ok = bool(inside[-1, :].any() or inside[:, -1].any())
    say(f"   is the CONSTANT model (a = 1 or b = 1, i.e. phi = beta) admissible? "
        f"{'YES' if const_ok else 'no'}")
    say("   -> that is the whole verdict. A pack indistinguishable from monodisperse")
    say("      constrains NO combination of a and b, because every combination that")
    say("      makes the binding branch vanish fits it equally well. u = 0.8 is")
    say("      uninformative about the combination, not only about a alone, so the")
    say("      monotone bracket a(0.667) <= a(0.7071) <= a(0.8) gains nothing at its")
    say("      upper end either.")
    say()

    say("=== what this does NOT decide ===")
    say("   * The ceiling is computed AT the shipped interpolated a(u), b(u) -- the very")
    say("     numbers the block exists to measure. If the truth at u = 0.7 is closer to 1")
    say("     than the PCHIP says, the signal is smaller and the real block does worse")
    say("     than its ceiling; if it is further from 1, better. The ceiling is a")
    say("     DESIGN comparison at a stated working point, not a prediction of the answer.")
    say("   * It says nothing about whether the model's FORM is right at high u.")
    say("   * MIN_PHASE_COUNT is applied per composition at N = 1500; compositions")
    say("     refused above would need N = 3000 and are not silently included.")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "diag_uident.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
