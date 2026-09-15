"""predict_ugap.py -- the BLIND predictions of PREREG_UGAP_2026-09-16.md.

Run and COMMITTED before a single U-GAP cell exists. Everything here is computed from
the shipped closure and the design; nothing reads a `phase12_ug*` file, and grade_ugap.py
refuses to run if one is present when this is invoked with --check.

Every caveat printed here MUST be carried into the claims file. On 2026-09-15
`out/predict_kladder.txt` flagged the u-gap before any cell ran and the claims shipped
without it until a referee asked; this file exists partly so that cannot happen twice,
and grade_ugap.py asserts the CAVEATS block below appears verbatim in the claims.

Run:  python predict_ugap.py      (writes out/predict_ugap.txt)
"""

from __future__ import annotations

import glob
import os
import sys

import numpy as np

import grade_mixture_ext as gx
import kladder as kl
import mixture_closure as mc
import mixture_model as mm
import ugap as ug

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = ug.PHI_SCALE
MU = ug.MU_P

CAVEATS = [
    "C-i   The measured knot enters at ONE friction (mu_p = 0.3), ONE pressure rung",
    "      (p_hat = 1e-3) and ONE N (1500). Nothing here says a(u) is pressure- or",
    "      friction-independent, and the closure has never been tested off that rung.",
    "C-ii  Arm A's knot is measured AT u = 0.70711 exactly, which is K = 5's adjacent",
    "      pair ratio on a span-4 ladder. It does NOT fill the rest of the gap: the",
    "      interval (0.70711, 1) still has no measured point, and PCHIP's shape there",
    "      remains an extrapolation between a measured knot and a proved limit.",
    "C-iii Inserting a knot at u* also moves a(u) and b(u) at NEIGHBOURING u through the",
    "      interpolant's derivatives. out/diag_ugap.txt measured that leakage: the K = 4",
    "      cells, whose pairs all sit at u <= 0.630, move up to 3.49 scale. So any K = 4",
    "      number recomputed after this block is NOT independent of it.",
    "C-iv  The identifiability rule's bar (delta chi2 = 1 per point, width < 0.25) is a",
    "      statement in units of PHI_SCALE, not of the seed error. It therefore does not",
    "      tighten with seeds or with the number of compositions -- demonstrated in",
    "      out/diag_uident.txt section 3 -- and a PASS is not a claim about precision.",
]


def main() -> int:
    lines: list[str] = []

    def say(s: str = "") -> None:
        lines.append(s)
        print(s, flush=True)

    stray = sorted(glob.glob(os.path.join(OUT, "phase12_ug*")))
    if stray:
        print(f"REFUSED: {len(stray)} U-GAP data files already exist; this file is BLIND "
              f"by construction:\n   " + "\n   ".join(stray), file=sys.stderr)
        return 2

    meas = mm.measured(mu_p=MU)
    beta = meas[("lognormal", 0.0)][0]
    allpts = gx.ext_points(beta)
    tab, fa, fb = gx.tabulated_closure(beta, allpts)

    say("=== predict_ugap.py -- BLIND, PREREG_UGAP_2026-09-16.md ===")
    say(f"   {ug.summary()}")
    say(f"   beta = {beta:.6f} (measured, mu_p {MU}); PHI_SCALE = {SCALE:.1e}; "
        f"identifiability width bar {ug.IDENT_WIDTH} (VERBATIM, not movable)")
    say(f"   no U-GAP data file exists: {len(stray)} matched out/phase12_ug*")
    say()

    w = ug.check_kernel(beta, np.random.default_rng(20260916))
    say(f"   ANCHOR -- vectorised kernel vs mixture_closure.packing_density_ab over 300 "
        f"random (r, a, b, x): worst |diff| = {w:.3e}")
    if w > 1e-12:
        say("   -> ANCHOR FAILED; nothing below is printed.")
        return 2
    say()

    say("=== the shipped closure this block is testing ===")
    say("   measured knots at u = " + ", ".join(f"{u:.4f}" for u in sorted(tab))
        + "; proved limits at u = 0 and 1; NOTHING between 0.6667 and 1.")
    say(f"   {'u':>9}{'a':>9}{'b':>9}")
    for u in sorted(tab):
        say(f"   {u:9.4f}{tab[u][0]:9.4f}{tab[u][1]:9.4f}")
    say(f"   interpolated at u* = {ug.U_STAR:.5f}:  a {float(fa(ug.U_STAR)):.4f}   "
        f"b {float(fb(ug.U_STAR)):.4f}   <- the number arm A replaces")
    say()

    # ------------------------------------------------------------------- the ceiling
    say("=== the design's IDENTIFIABILITY CEILING on the registered grid ===")
    say("   Noise-free phi generated from the shipped closure at each arm's own u, then")
    say("   passed through the registered rule. This is the narrowest interval a perfect")
    say("   experiment at this design could return -- an upper bound on the block, not a")
    say("   promise about it.")
    say(f"   {'arm':>4}{'r':>10}{'u':>9}{'n_x':>5}{'x*cross':>9}"
        f"{'n_c':>5}{'n_f':>5}{'w(a)':>8}{'w(b)':>8}  ident?")
    for arm, r in ug.RATIOS.items():
        u = 1.0 / r
        av, bv = float(fa(u)), float(fb(u))
        xs = [c["x"] for c in ug.admitted() if c["arm"] == arm]
        pts = [(r, x, float(mc.packing_density_ab(np.array([r, 1.0]),
                                                  np.array([1 - x, x]), beta, av, bv)))
               for x in xs]
        xc = ug.crossover(beta, av, bv)
        nc = sum(1 for x in xs if x < xc)
        ok, ia, ib, _, _ = ug.identifiable(pts, beta)
        say(f"   {arm:>4}{r:10.5f}{u:9.5f}{len(xs):5d}{xc:9.4f}{nc:5d}{len(xs)-nc:5d}"
            f"{ia[1]-ia[0]:8.3f}{ib[1]-ib[0]:8.3f}  {'YES' if ok else 'NO'}")
    say("   (for comparison, out/diag_uident.txt section 5: the OLD 8-point grid's")
    say("   ceiling is identifiable at r = 1.41421 and 1.5 and NOT at r = 1.25.)")
    say()

    say("=== ROBUSTNESS: how wrong can the shipped a(u*) be before arm A fails? ===")
    say("   Rows: assumed true a(u*). Columns: assumed true b(u*). 'Y' = both widths")
    say(f"   under {ug.IDENT_WIDTH}; otherwise the WIDER of the two is printed.")
    a_grid = (0.90, 0.92, 0.9398, 0.96, 0.97, 0.98, 0.99)
    b_grid = (0.86, 0.90, 0.9177, 0.94, 0.96, 0.98, 0.99)
    xs_a = [c["x"] for c in ug.admitted() if c["arm"] == "A"]
    say("        a\\b " + "".join(f"{v:>9.3f}" for v in b_grid))
    a_limit = 0.0
    for av_ in a_grid:
        row = []
        for bv_ in b_grid:
            pts = [(ug.R_STAR, x,
                    float(mc.packing_density_ab(np.array([ug.R_STAR, 1.0]),
                                                np.array([1 - x, x]), beta, av_, bv_)))
                   for x in xs_a]
            ok, ia, ib, _, _ = ug.identifiable(pts, beta)
            if ok:
                a_limit = max(a_limit, av_)
            row.append("        Y" if ok
                       else f"{max(ia[1]-ia[0], ib[1]-ib[0]):9.2f}")
        say(f"      {av_:6.4f}" + "".join(row))
    say(f"   -> arm A stays identifiable up to a true a(u*) of {a_limit:.4f}.")
    say()

    # ------------------------------------------------ the held-out evidence on direction
    say("=== the HELD-OUT evidence that sets the direction prediction ===")
    say("   r = 1.25 failed identifiability, so `tabulated_closure` never put a knot")
    say("   there: its phi data is genuinely OUT of the closure and is a legal test of")
    say("   the extrapolation. Residuals phi_measured - phi_shipped, in scale units:")
    for r in (1.25, 1.5):
        u = 1.0 / r
        av, bv = float(fa(u)), float(fb(u))
        sub = [p for p in allpts if abs(p[0] - r) < 1e-9 and 0.0 < p[1] < 1.0]
        res = np.array([(phi - float(mc.packing_density_ab(
            np.array([r, 1.0]), np.array([1 - x, x]), beta, av, bv))) / SCALE
            for _, x, phi in sub])
        tag = "HELD OUT" if r == 1.25 else "a FITTED knot (control)"
        say(f"   r = {r:4.2f} (u {u:.4f}, {tag}): n = {res.size}, mean {res.mean():+.2f}, "
            f"rms {np.sqrt((res**2).mean()):.2f}, negative on {int((res < 0).sum())}"
            f"/{res.size}")
        say("      " + " ".join(f"{v:+.2f}" for v in res))
    say("   -> the shipped extrapolation OVER-predicts filling where it is held out, on")
    say("      8 of 8. If that direction persists at u*, the true a(u*) is ABOVE 0.9398,")
    say("      which is the direction that threatens BOTH arm A's identifiability and")
    say("      the K-ladder's saturation claim. That is why it is registered as U2.")
    say()

    # -------------------------------------------------------- the K-ladder consequence
    say("=== the K-ladder consequence the block exists to settle (registered as U3) ===")
    say("   out/diag_ugap.txt swept a(u*), b(u*) over the monotone-admissible range and")
    say("   found R(5) - R(4) can be carried across zero on 12 of 12 rows, so a measured")
    say("   knot decides the saturation verdict rather than confirming it. Blind")
    say("   prediction, evaluated by reinserting a knot at u* and recomputing ADDtab:")
    say(f"   {'N':>6}{'y_S':>7}{'R5-R4 now':>11}{'at a=0.9398':>13}"
        f"{'at a=0.97':>11}{'at a=0.99':>11}")
    import csv
    disk = {}
    with open(os.path.join(OUT, "kladder_cells.csv"), encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            disk[(int(row["N"]), row["label"])] = row
    by_label = {}
    for nn in (kl.N_PROD, kl.N_BIG):
        for c in kl.build_cells(nn)[0]:
            by_label.setdefault(c["label"], c)
    from scipy.interpolate import PchipInterpolator
    us0 = [0.0, 1.0] + sorted(tab)
    aa0 = [0.0, 1.0] + [tab[u][0] for u in sorted(tab)]
    bb0 = [0.0, 1.0] + [tab[u][1] for u in sorted(tab)]
    o = np.argsort(us0)
    uu0 = np.asarray(us0)[o]
    aa0 = np.asarray(aa0)[o]
    bb0 = np.asarray(bb0)[o]
    bshift = float(fb(ug.U_STAR))

    def R_with(n, label, a_star):
        row = disk[(n, label)]
        c = by_label[label]
        bu = float(row["beta_used"])
        ua = np.append(uu0, ug.U_STAR)
        va = np.append(aa0, a_star)
        ub = np.append(uu0, ug.U_STAR)
        vb = np.append(bb0, bshift)
        oa, ob = np.argsort(ua), np.argsort(ub)
        f_a = PchipInterpolator(ua[oa], va[oa])
        f_b = PchipInterpolator(ub[ob], vb[ob])
        add = float(mc.packing_density_tabulated(np.asarray(c["d_rel"], float),
                                                 np.asarray(c["y"], float), bu, f_a, f_b))
        return (float(row["phi"]) - add) / SCALE

    for n in (1500, 3000):
        for ys in kl.YS_BALANCED:
            l4, l5 = f"K4_yS{ys:.2f}", f"K5_yS{ys:.2f}"
            if (n, l4) not in disk or (n, l5) not in disk:
                continue
            now = (float(disk[(n, l5)]["R"]) - float(disk[(n, l4)]["R"]))
            vals = [R_with(n, l5, a) - R_with(n, l4, a) for a in (0.9398, 0.97, 0.99)]
            say(f"   {n:6d}{ys:7.2f}{now:+11.2f}{vals[0]:+13.2f}"
                f"{vals[1]:+11.2f}{vals[2]:+11.2f}")
    say("   (the a = 0.9398 column re-inserts the closure's OWN interpolated value as an")
    say("   explicit knot, so any difference from 'R5-R4 now' is PCHIP reshaping, not")
    say("   physics; it is printed so the two effects cannot be confused.)")
    say()

    say("=== CAVEATS -- these must appear in the claims file ===")
    for c in CAVEATS:
        say("   " + c)
    say()

    say("=== per-cell blind phi predictions, shipped closure ===")
    say(f"   {'arm':>4}{'r':>10}{'x':>7}{'phi_pred':>10}{'branch':>8}{'n_min':>8}")
    for c in ug.admitted():
        r, x = c["ratio"], c["x"]
        u = 1.0 / r
        av, bv = float(fa(u)), float(fb(u))
        p = float(mc.packing_density_ab(np.array([r, 1.0]), np.array([1 - x, x]),
                                        beta, av, bv))
        br = "coarse" if x < ug.crossover(beta, av, bv) else "fine"
        say(f"   {c['arm']:>4}{r:10.5f}{x:7.2f}{p:10.5f}{br:>8}"
            f"{c['n_min_at_judge']:8.0f}")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "predict_ugap.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
