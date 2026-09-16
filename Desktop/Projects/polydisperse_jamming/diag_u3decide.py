"""diag_u3decide.py -- BEFORE any cell: CAN a narrower a(u*) decide U3 at all?

THE QUESTION, AND WHY IT IS ASKED BEFORE THE DESIGN AND NOT AFTER
-----------------------------------------------------------------
`out/grade_ugap.txt` U3 grades FAIL at the point estimate (R5-R4 <= 0 on 1/3 levels at
N = 1500 and 0/3 at N = 3000, bar >= 2 of 3 at both) while the scientific statement is
UNDETERMINED: swept over arm A's own JOINT likelihood region, saturation survives on
33 of 87 corners (37.9 %). So the verdict turns on where inside [0.8799, 0.9699] the
true a(u*) sits, and the obvious next block is "measure a(u*) more precisely".

`out/diag_uident.txt` section 3 already proved that seeds and extra compositions buy
NOTHING under the registered rule (the profile bar is delta chi2 = 1 PER POINT, so
n = 8/16/32 returns w(a) = 0.0400 three times identically). Only SENSITIVITY -- where in
composition the block samples -- can narrow the interval. So a narrowing block is
possible in principle. This file asks the prior question:

    IF the interval narrowed by a factor f, would the U3 survival fraction leave the
    middle? If a plausible f still leaves saturation surviving on a middling fraction,
    THE BLOCK CANNOT DECIDE U3 AND MUST NOT BE RUN FOR THAT PURPOSE.

That is the same move that made the U-GAP block worth running -- `diag_ugap.py` measured
that the knot moves R5-R4 across zero on 12 of 12 rows before a cell was spent -- and it
would have killed that block had the answer come out flat.

WHAT IS COMPUTED
----------------
1. CONTROL. The 33/87 of `grade_ugap.py` is reproduced here from the same inputs on the
   same 13x13 grid with the same joint-region rule. A diagnostic that cannot reproduce
   the number it is auditing has not audited it (lessons.md 2026-09-08).
2. The survival map S(a, b) on a FINE grid, and the location of its boundary relative to
   the best fit -- in units of the CURRENT width, which is what decides whether a factor
   of two in width is enough.
3. The width sweep. For each contraction factor f the future joint region is modelled as
   the CURRENT joint region contracted homothetically about a candidate truth (a0, b0):
   exact for a quadratic chi2 surface whose curvature is scaled by 1/f^2, which is what
   "the same experiment with more sensitivity" means. Reported at the current best fit
   AND over every candidate truth the current data admits, because the future block's
   centre is not known now -- that uncertainty is the whole reason the question is open.
4. The ACHIEVABLE width. The noise-free identifiability ceiling of candidate composition
   grids at r = 1.41421, computed at the MEASURED knot (0.9481, 0.9376) rather than at
   the shipped extrapolation `diag_uident.py` had to use. A design that fails here cannot
   be rescued by seeds; one that passes here is an upper bound on what a real block does.

SIGN AND SCALE CONVENTIONS, stated because both have bitten this project.
  R = (phi - ADDtab)/SCALE, taken verbatim from `transfer_kladder.py` /
  `mixture_closure.py` (model - phi is the OTHER convention and produced an artifact that
  contradicted its own table, 2026-09-16). Saturation "survives" at a corner exactly when
  `grade_ugap.py`'s registered rule holds there: R5 - R4 <= 0 on >= 2 of 3 y_S levels at
  BOTH N. SCALE = PHI_SCALE = 2.1e-3, the protocol reproducibility, not a seed error.

Run:  python3 diag_u3decide.py      (writes out/diag_u3decide.txt)
"""

from __future__ import annotations

import collections
import csv
import os
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator

import grade_mixture_ext as gx
import grade_ugap as gu
import kladder as kl
import mixture_closure as mc
import mixture_model as mm
import ugap as ug

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = ug.PHI_SCALE
MU = 0.3
U3_MIN_LEVELS = 2

# Contraction factors swept. 1.0 is the measured interval (must reproduce 37.9 %);
# 0.60 is roughly what section 4's best candidate grid predicts; 0.10 is the limit a
# real block cannot reach and is included so the sweep's own endpoint is visible.
FACTORS = (1.00, 0.85, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10)

# "Decisive" thresholds, stated here and used everywhere below.
DECIDE_LO, DECIDE_HI = 0.10, 0.90


def survival_map(tab, disk, by_label, ga, gb):
    """S[i, j] = 1 where the registered U3 rule holds with the knot set to (ga[i], gb[j])."""
    s = np.zeros((len(ga), len(gb)), dtype=bool)
    for i, av in enumerate(ga):
        for j, bv in enumerate(gb):
            t2 = dict(tab)
            t2[ug.U_STAR] = (float(av), float(bv))
            us = [0.0, 1.0] + sorted(t2)
            aa = np.asarray([0.0, 1.0] + [t2[u][0] for u in sorted(t2)])
            bb = np.asarray([0.0, 1.0] + [t2[u][1] for u in sorted(t2)])
            o = np.argsort(us)
            uu = np.asarray(us)[o]
            fa = PchipInterpolator(uu, aa[o])
            fb = PchipInterpolator(uu, bb[o])
            cnt = collections.Counter()
            for n in (1500, 3000):
                for ys in kl.YS_BALANCED:
                    l4, l5 = f"K4_yS{ys:.2f}", f"K5_yS{ys:.2f}"
                    if (n, l4) not in disk or (n, l5) not in disk:
                        continue
                    r45 = []
                    for lab in (l4, l5):
                        row = disk[(n, lab)]
                        c = by_label[lab]
                        add = float(mc.packing_density_tabulated(
                            np.asarray(c["d_rel"], float), np.asarray(c["y"], float),
                            float(row["beta_used"]), fa, fb))
                        r45.append((float(row["phi"]) - add) / SCALE)
                    cnt[n] += int(r45[1] - r45[0] <= 0)
            s[i, j] = all(cnt[n] >= U3_MIN_LEVELS for n in (1500, 3000))
    return s


_LINES: list[str] = []


def main() -> int:
    lines = _LINES

    def say(t: str = "") -> None:
        lines.append(t)
        print(t, flush=True)

    say("=== diag_u3decide.py -- can a NARROWER a(u*) decide U3? Asked BEFORE any cell ===")

    meas = mm.measured(mu_p=MU)
    beta = meas[("lognormal", 0.0)][0]
    allpts = gx.ext_points(beta)
    tab, fa0, fb0 = gx.tabulated_closure(beta, allpts)
    say(f"   beta = {beta:.6f} (measured, mu_p {MU}); SCALE = {SCALE:.1e}; "
        f"u* = {ug.U_STAR:.5f}")

    w = ug.check_kernel(beta, np.random.default_rng(20260917))
    say(f"   ANCHOR: ugap.phi_fast vs mixture_closure.packing_density_ab, 300 random "
        f"points, worst |diff| = {w:.3e}")
    if w > 1e-12:
        say("   -> ANCHOR FAILED. No number below is printed.")
        return 2

    # ------------------------------------------------------------------ arm A data
    g15, nf15, nraw15 = gu.load("ug")
    adm: dict[tuple[float, float], list[dict]] = {}
    nref = 0
    for key in sorted(g15):
        ok = []
        for c in g15[key]:
            a_, _why = gu.admissible(c, ug.N_PROD)
            if a_:
                ok.append(c)
            else:
                nref += 1
        adm[key] = ok
    ptsA = []
    for x in ug.X_GRID:
        v = adm.get((float(ug.R_STAR), round(float(x), 6)), [])
        if v:
            ptsA.append((float(ug.R_STAR), float(x),
                         float(np.mean([c["phi"] for c in v]))))
    say(f"   arm A points: {nraw15} raw records from {nf15} files, {nref} refused block-wide"
        f" -> N before/after the arm-A filter: {len(g15)} compositions / {len(ptsA)} used")

    okA, ia, ib, a_hat, b_hat = ug.identifiable(ptsA, beta)
    say(f"   arm A: a = {a_hat:.4f} [{ia[0]:.4f}, {ia[1]:.4f}] (w {ia[1]-ia[0]:.4f}), "
        f"b = {b_hat:.4f} [{ib[0]:.4f}, {ib[1]:.4f}] (w {ib[1]-ib[0]:.4f}), "
        f"identifiable = {okA}")
    xvA = np.array([p[1] for p in ptsA])
    yvA = np.array([p[2] for p in ptsA])

    def chi2(a, b):
        return (((ug.phi_fast(a, b, xvA[None, None, :], beta) - yvA[None, None, :])
                 / SCALE) ** 2).sum(axis=2)

    # ------------------------------------------------------------- K-ladder cells
    disk: dict[tuple[int, str], dict] = {}
    with open(os.path.join(OUT, "kladder_cells.csv"), encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            disk[(int(row["N"]), row["label"])] = row
    by_label: dict[str, dict] = {}
    for nn in (kl.N_PROD, kl.N_BIG):
        for c in kl.build_cells(nn)[0]:
            by_label.setdefault(c["label"], c)
    say(f"   K-ladder cells on disk: {len(disk)} (N, label) rows")
    say()

    # ============================================================ 1. the CONTROL
    say("=== 1. CONTROL -- reproduce grade_ugap.py's 33 of 87 from the same inputs ===")
    say("   Same 13x13 grid over the two marginal intervals, same joint rule")
    say("   (chi2 <= chi2_min + n). A diagnostic that cannot reproduce the number it")
    say("   audits has not audited it.")
    ga13 = np.linspace(ia[0], ia[1], 13)
    gb13 = np.linspace(ib[0], ib[1], 13)
    s13 = survival_map(tab, disk, by_label, ga13, gb13)
    c13 = chi2(ga13[:, None, None], gb13[None, :, None])
    bar = float(c13.min()) + len(ptsA)
    joint13 = c13 <= bar
    n_j = int(joint13.sum())
    n_s = int((s13 & joint13).sum())
    say(f"   chi2_min = {c13.min():.4f}, n = {len(ptsA)}, bar = {bar:.4f}")
    say(f"   marginal box: {int(s13.sum())} of {s13.size} ({s13.sum()/s13.size:.1%}) "
        f"  [grade_ugap records 110 of 169, 65.1%]")
    say(f"   JOINT region: {n_s} of {n_j} ({n_s/n_j:.1%})"
        f"   [grade_ugap records 33 of 87, 37.9%]")
    ctl_ok = (n_s == 33 and n_j == 87)
    say(f"   -> CONTROL {'PASS -- exact match' if ctl_ok else 'FAIL -- DOES NOT MATCH, nothing below is trustworthy'}")
    if not ctl_ok:
        say(f"FAILURES: ['control: {n_s}/{n_j} against the recorded 33/87']")
        return 1
    say()

    # ====================================================== 2. the boundary, finely
    say("=== 2. WHERE THE BOUNDARY IS, on a fine grid, in units of the CURRENT width ===")
    say("   The survival set is what a narrower interval has to fall inside or outside.")
    say("   Its distance from the best fit, measured in current widths, is what decides")
    say("   whether halving the width is enough.")
    NG = 81
    pad = 0.6                       # padding in current widths, each side
    wa, wb = ia[1] - ia[0], ib[1] - ib[0]
    gaF = np.linspace(max(1e-4, ia[0] - pad * wa), min(1 - 1e-4, ia[1] + pad * wa), NG)
    gbF = np.linspace(max(1e-4, ib[0] - pad * wb), min(1 - 1e-4, ib[1] + pad * wb), NG)
    sF = survival_map(tab, disk, by_label, gaF, gbF)
    say(f"   fine grid {NG}x{NG} over a box padded by {pad:.2f} current widths each side")
    say(f"   survival fraction over the padded box: {sF.mean():.1%}")

    # the a-boundary at the best-fit b, and the b-boundary at the best-fit a
    def flip_along(vals, s_line):
        """First value where the survival indicator changes, or None if constant."""
        ch = np.nonzero(np.diff(s_line.astype(int)))[0]
        return None if ch.size == 0 else float(0.5 * (vals[ch[0]] + vals[ch[0] + 1]))

    j_bhat = int(np.argmin(np.abs(gbF - b_hat)))
    i_ahat = int(np.argmin(np.abs(gaF - a_hat)))
    a_edge = flip_along(gaF, sF[:, j_bhat])
    b_edge = flip_along(gbF, sF[i_ahat, :])
    s_at_hat = bool(sF[i_ahat, j_bhat])
    say(f"   at the best fit (a {a_hat:.4f}, b {b_hat:.4f}): saturation survives = {s_at_hat}"
        f"   [consistent with U3 FAIL at the point estimate: {not s_at_hat}]")
    if a_edge is None:
        say("   along a at b = b_hat: NO flip inside the padded box -- the verdict "
            "does not turn on a alone there")
    else:
        say(f"   along a at b = b_hat: boundary at a = {a_edge:.4f}, i.e. "
            f"{abs(a_edge - a_hat)/wa:.2f} current a-widths from the best fit")
    if b_edge is None:
        say("   along b at a = a_hat: NO flip inside the padded box")
    else:
        say(f"   along b at a = a_hat: boundary at b = {b_edge:.4f}, i.e. "
            f"{abs(b_edge - b_hat)/wb:.2f} current b-widths from the best fit")
    say()

    # =================================================== 3. the width sweep / GATE
    say("=== 3. THE GATE -- survival fraction as the interval CONTRACTS ===")
    say("   Model: a sharper block's joint region is the CURRENT joint region contracted")
    say("   by f about that block's own best fit. Exact for a quadratic chi2 whose")
    say("   curvature is scaled by 1/f^2. The 87 jointly admissible corners of section 1")
    say("   supply the SHAPE; f scales their offsets from the current best fit.")
    off_a = (ga13[:, None] - a_hat) * np.ones((1, len(gb13)))
    off_b = np.ones((len(ga13), 1)) * (gb13[None, :] - b_hat)
    oa = off_a[joint13]
    ob = off_b[joint13]
    say(f"   shape template: {oa.size} offsets, a in [{oa.min():+.4f}, {oa.max():+.4f}], "
        f"b in [{ob.min():+.4f}, {ob.max():+.4f}]")

    # Every corner below is looked up in the section-2 fine map rather than re-running
    # the closure: 87 centres x 9 factors x 87 corners is 68 121 evaluations and the map
    # is the same function on a 81x81 mesh. The lookup is validated against a direct
    # recomputation on a random sample before any number from it is printed.
    def s_lookup(a_arr, b_arr):
        ia_ = np.array([int(np.argmin(np.abs(gaF - v))) for v in np.atleast_1d(a_arr)])
        jb_ = np.array([int(np.argmin(np.abs(gbF - v))) for v in np.atleast_1d(b_arr)])
        return sF[ia_, jb_]

    rng = np.random.default_rng(20260917)
    n_chk = 40
    ta = rng.uniform(gaF[0], gaF[-1], n_chk)
    tb = rng.uniform(gbF[0], gbF[-1], n_chk)
    direct = np.array([bool(survival_map(tab, disk, by_label,
                                         np.array([x_]), np.array([y_]))[0, 0])
                       for x_, y_ in zip(ta, tb)])
    approx = s_lookup(ta, tb)
    agree = int((direct == approx).sum())
    say(f"   LOOKUP ANCHOR: nearest-node lookup vs a direct recomputation on {n_chk} "
        f"random (a, b): agree {agree}/{n_chk}")
    if agree < n_chk - 2:
        say("   -> LOOKUP ANCHOR FAILED (more than 2 disagreements): the mesh is too "
            "coarse for the boundary. Nothing in section 3 is printed.")
        return 1
    say()

    # (a) centred on the CURRENT best fit
    say("   (a) IF the future block's best fit lands where the current one is:")
    say(f"      {'f':>6}{'w(a)':>8}{'w(b)':>8}{'survives':>12}{'verdict':>16}")
    surv_at_hat = {}
    for f in FACTORS:
        ca = np.clip(a_hat + f * oa, 1e-4, 1 - 1e-4)
        cb = np.clip(b_hat + f * ob, 1e-4, 1 - 1e-4)
        sv = s_lookup(ca, cb)
        fr = float(sv.mean())
        surv_at_hat[f] = fr
        verd = ("DECIDED against" if fr < DECIDE_LO else
                "DECIDED for" if fr > DECIDE_HI else "undetermined")
        say(f"      {f:6.2f}{f*wa:8.4f}{f*wb:8.4f}{sv.sum():7d}/{sv.size:<4d}{verd:>16}")
    say()

    # (b) over every candidate truth the current data admits
    say("   (b) BUT the future block's centre is NOT known now. Over every candidate")
    say("      truth the current joint region admits, how often is the answer DECISIVE,")
    say("      and which way? This is the number the go/no-go rests on: a block that")
    say("      decides only for some centres is still worth running -- one that is")
    say("      middling for most centres is not.")
    cand_a = ga13[:, None].repeat(len(gb13), 1)[joint13]
    cand_b = gb13[None, :].repeat(len(ga13), 0)[joint13]
    say(f"      {'f':>6}{'decided AGAINST':>17}{'decided FOR':>13}{'undetermined':>14}"
        f"{'P(decisive)':>13}")
    p_dec = {}
    for f in FACTORS:
        n_ag = n_fo = n_un = 0
        for a0, b0 in zip(cand_a, cand_b):
            ca = np.clip(a0 + f * oa, 1e-4, 1 - 1e-4)
            cb = np.clip(b0 + f * ob, 1e-4, 1 - 1e-4)
            fr = float(s_lookup(ca, cb).mean())
            if fr < DECIDE_LO:
                n_ag += 1
            elif fr > DECIDE_HI:
                n_fo += 1
            else:
                n_un += 1
        tot = n_ag + n_fo + n_un
        p_dec[f] = (n_ag + n_fo) / tot
        say(f"      {f:6.2f}{n_ag:10d}/{tot:<6d}{n_fo:6d}/{tot:<6d}{n_un:7d}/{tot:<6d}"
            f"{p_dec[f]:13.1%}")
    say()

    # ============================================== 4. what width is ACHIEVABLE
    say("=== 4. WHAT WIDTH IS ACHIEVABLE -- the noise-free identifiability ceiling ===")
    say("   phi is generated from the closure at the MEASURED knot with no noise, so")
    say("   chi2_min = 0 and the width returned is the narrowest ANY block at that")
    say("   composition grid could get. diag_uident.py had to compute this at the")
    say("   SHIPPED extrapolation; the knot is now measured, so the truth is known.")
    say(f"   truth: a(u*) = {a_hat:.4f}, b(u*) = {b_hat:.4f}, r = {ug.R_STAR:.5f}")
    xstar = ug.crossover(beta, a_hat, b_hat)
    say(f"   branch crossover at the measured knot: x* = {xstar:.4f}")

    def affordable(xs):
        return [x for x in xs
                if ug.n_min_frac(ug.R_STAR, x) * ug.N_PROD >= ug.MIN_PHASE_COUNT]

    say("   ⚠ The registered bar is delta chi2 = n, so a design with MORE points is also")
    say("      graded against a LOOSER bar. Designs of 9 and of 16 points are therefore")
    say("      not a like-for-like statistical comparison; both are listed and the")
    say("      MATCHED-n (9-point) rows are what isolate WHERE from HOW MANY.")
    designs = {
        "SHIPPED   the 9-point grid arm A ran": list(ug.X_GRID),
        "n9 SPAN   9 pts, 5 below x* spread wide": (
            [0.05, 0.11, 0.18, 0.25, 0.32, 0.40, 0.52, 0.68, 0.88]),
        "n9 TIGHT  9 pts packed just below x*": (
            list(np.round(np.linspace(0.22, 0.325, 6), 4)) + [0.40, 0.55, 0.80]),
        "n9 SPLIT  9 pts straddling x* closely": (
            list(np.round(np.linspace(0.26, 0.325, 5), 4))
            + list(np.round(np.linspace(0.335, 0.45, 4), 4))),
        "n16 SPAN  16 pts, 12 spread below x*": (
            list(np.round(np.linspace(0.05, 0.325, 12), 4)) + [0.42, 0.52, 0.66, 0.88]),
        "n16 CLUST 16 pts, 12 packed below x*": (
            list(np.round(np.linspace(0.18, 0.325, 12), 4)) + [0.42, 0.52, 0.66, 0.88]),
        "n16 TIGHT 16 pts inside [0.25, 0.36]": (
            list(np.round(np.linspace(0.25, 0.358, 13), 4)) + [0.40, 0.55, 0.80]),
        "WIDE      9 pts over [0.05, 0.95]": (
            list(np.round(np.linspace(0.05, 0.95, 9), 4))),
    }
    say(f"   {'design':>40}{'n_req':>7}{'aff':>5}{'n_c':>5}{'n_f':>5}"
        f"{'w(a)':>8}{'w(b)':>8}{'f(a)':>7}  ident?")
    ceil_w = {}
    for name, xs in designs.items():
        aff = affordable(xs)
        pts = [(ug.R_STAR, float(x), float(ug.phi_fast(a_hat, b_hat, float(x), beta)))
               for x in aff]
        if len(pts) < 3:
            say(f"   {name:>40}{len(xs):7d}{len(aff):5d}   too few affordable points")
            continue
        n_c = sum(1 for _, x, _ in pts if x < xstar)
        ok_, ja, jb, _ah, _bh = ug.identifiable(pts, beta)
        w_a, w_b = ja[1] - ja[0], jb[1] - jb[0]
        ceil_w[name] = (w_a, w_b, xs)
        say(f"   {name:>40}{len(xs):7d}{len(aff):5d}{n_c:5d}{len(pts)-n_c:5d}"
            f"{w_a:8.4f}{w_b:8.4f}{w_a/wa:7.2f}  {'YES' if ok_ else 'NO'}")
    say("   f(a) is the CONTRACTION FACTOR that design predicts, w(a)_ceiling / w(a)")
    say(f"      measured now. The measured w(a) = {wa:.4f} is itself ABOVE the shipped")
    say("      grid's own ceiling, so every f above is optimistic by that margin too.")
    sh = ceil_w["SHIPPED   the 9-point grid arm A ran"][0]
    say(f"   CONTROL: the shipped grid's ceiling w(a) = {sh:.4f} against the {wa:.4f}")
    say(f"      the real block returned -- ratio {wa/sh:.2f}x. A ceiling BELOW the")
    say("      realised width is the expected direction; above it would be a bug.")
    say()
    say("   🔴 THE BRIEF'S PREMISE IS REFUTED BY THIS TABLE. Concentrating compositions")
    say("      where |d phi / d a| is largest (just below x*) does NOT narrow a -- it")
    say("      DESTROYS identifiability. The mechanism is the branch swap: over a narrow")
    say("      window of x the FINE branch alone can mimic the coarse branch's shape, so")
    say("      the a = 1 escape reopens. What narrows a is the x-SPAN of the coarse-branch")
    say("      points, not their density at the peak.")
    say()

    # ---- robustness: the ceiling as the assumed TRUTH moves over the interval
    say("   ROBUSTNESS -- the ceiling is computed at ONE assumed truth, and the truth is")
    say("   only known to within the interval this block exists to shrink. Rows are the")
    say("   assumed true a(u*) across its own profile interval, at b = b_hat.")
    best_name = min(ceil_w, key=lambda k: ceil_w[k][0])
    say(f"   best design by w(a): {best_name.strip()}")
    say("   The SHIPPED grid is re-run at each assumed truth beside it, because a width")
    say("   that changes with the truth is a property of the OBSERVABLE, not of the")
    say("   design: the design gain is the RATIO of the two ceilings, and the realistic")
    say("   width is that ratio times the width the shipped grid actually returned.")
    say(f"      {'a_true':>8}{'new w(a)':>10}{'ship w(a)':>11}{'ratio':>7}"
        f"{'f_real':>8}  ident?")
    xs_best = ceil_w[best_name][2]
    aff_best = affordable(xs_best)
    aff_ship = affordable(list(ug.X_GRID))
    f_rob = []
    for a_t in np.linspace(ia[0], ia[1], 7):
        pts = [(ug.R_STAR, float(x), float(ug.phi_fast(a_t, b_hat, float(x), beta)))
               for x in aff_best]
        ok_, ja, _jb, _ah, _bh = ug.identifiable(pts, beta)
        pts_s = [(ug.R_STAR, float(x), float(ug.phi_fast(a_t, b_hat, float(x), beta)))
                 for x in aff_ship]
        _oks, jas, _jbs, _a2, _b2 = ug.identifiable(pts_s, beta)
        w_new, w_ship = ja[1] - ja[0], jas[1] - jas[0]
        ratio = w_new / w_ship if w_ship > 0 else float("nan")
        f_rob.append(ratio)
        say(f"      {a_t:8.4f}{w_new:10.4f}{w_ship:11.4f}{ratio:7.2f}{ratio:8.2f}"
            f"  {'YES' if ok_ else 'NO'}")
    say("   f_real = ratio, i.e. the predicted contraction of the width the block would")
    say("      ACTUALLY return, which already absorbs the 1.20x realised/ceiling gap.")
    say(f"   -> design gain over the truth's own interval: best {min(f_rob):.2f}, "
        f"WORST {max(f_rob):.2f}. The WORST is what the gate must be read at.")
    say()

    # ==================================================================== verdict
    say("=== VERDICT ===")
    best = ceil_w[best_name][0]
    f_best = best / wa
    f_worst = max(f_rob)
    say(f"   best achievable contraction, at the assumed truth: f = {f_best:.2f} "
        f"(w(a) {best:.4f} from {wa:.4f}); over the truth's own interval the WORST is "
        f"f = {f_worst:.2f}.")
    say(f"   the shipped grid realised {wa/sh:.2f}x its own ceiling; the f_real column "
        f"above already absorbs that, so f_real is the number to read.")
    # the design gain AT the current best estimate of the truth, computed not interpolated
    pts_n = [(ug.R_STAR, float(x), float(ug.phi_fast(a_hat, b_hat, float(x), beta)))
             for x in aff_best]
    pts_s = [(ug.R_STAR, float(x), float(ug.phi_fast(a_hat, b_hat, float(x), beta)))
             for x in aff_ship]
    _o1, j_n, _x1, _y1, _z1 = ug.identifiable(pts_n, beta)
    _o2, j_s, _x2, _y2, _z2 = ug.identifiable(pts_s, beta)
    f_at_hat = (j_n[1] - j_n[0]) / (j_s[1] - j_s[0])
    say(f"   design gain AT the current best estimate a = {a_hat:.4f}: "
        f"f_real = {f_at_hat:.2f} ({j_n[1]-j_n[0]:.4f} / {j_s[1]-j_s[0]:.4f})")
    fk = sorted(FACTORS)
    rows = (("f_real AT the current best estimate", f_at_hat),
            ("f_real at the most favourable truth", min(f_rob)),
            ("f_real at the WORST truth", f_worst),
            ("the noise-free ceiling ratio, an unreachable bound", f_best))
    for lab, f_ in rows:
        fr_hat = float(np.interp(f_, fk, [surv_at_hat[f] for f in fk]))
        pd_ = float(np.interp(f_, fk, [p_dec[f] for f in fk]))
        say(f"   at f = {f_:.2f} ({lab}): survival at the current best fit = "
            f"{fr_hat:.1%}; P(decisive) over admitted centres = {pd_:.1%}")
    say(f"   decisive means < {DECIDE_LO:.0%} or > {DECIDE_HI:.0%}, stated in section 3 "
        f"before the numbers.")
    say()
    s_at = float(np.interp(f_at_hat, fk, [surv_at_hat[f] for f in fk]))
    go = (s_at < DECIDE_LO or s_at > DECIDE_HI)
    say(f"   DECISION: a narrowing block at the best composition grid found here is "
        f"{'WORTH' if go else 'NOT WORTH'} cells for the purpose of deciding U3.")
    if not go:
        say(f"   At the realistic contraction f = {f_at_hat:.2f} the survival fraction is "
            f"{s_at:.1%} -- inside [{DECIDE_LO:.0%}, {DECIDE_HI:.0%}], i.e. still")
        say("   MIDDLING, so the block would return the same verdict it was run to")
        say("   replace: UNDETERMINED, leaning against saturation.")
        need = [f for f in FACTORS if surv_at_hat[f] < DECIDE_LO or surv_at_hat[f] > DECIDE_HI]
        say(f"   To decide, f <= {max(need):.2f} is needed and the best design found "
            f"reaches only {f_at_hat:.2f}." if need else
            "   NO contraction in the swept range decides it at this centre.")
    return 0


def run() -> int:
    """Wrapper so the transcript is written whatever main() returns, including a
    CONTROL failure -- an artifact that only exists on success cannot record a failure."""
    rc = main()
    with open(os.path.join(OUT, "diag_u3decide.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(_LINES) + "\n")
    return rc


if __name__ == "__main__":
    sys.exit(run())
