"""diag_ugap.py -- BEFORE any new cell: is the K = 5 prediction sensitive to a(u), b(u)
in the unmeasured gap 0.667 < u < 1 at all?

WHY THIS RUNS FIRST
-------------------
`CLAIMS_KLADDER_2026-09-15.md` C1a(ii) carries a caveat the referee extracted: the
tabulated closure has NO measured knot above u = 0.667, and 4 of K = 5's 10 pairs sit at
u = 1/4^(1/4) = 0.70711 -- its most numerous pair type -- so the rung that decides
"and then SATURATES" is evaluated where the closure has no data.

That caveat is only load-bearing IF phi actually moves when a(0.70711) and b(0.70711)
move. The ternary block already produced the mirror case: `a(1/4)` turned out to be
DEGENERATE (phi bit-identical from -3 to +0.6, because the L-S pair term never entered
the binding branch) and the inversion built on it was WITHDRAWN. A degenerate parameter
cannot testify either way -- in EITHER direction. So the honest first move is to ask
what the K = 5 cells' ADDtab actually DOES as those two numbers are swept.

Three outcomes, all informative and all decided by ONE bar (TAB's own in-sample rms,
2.02 scale -- not a new number):

  (A) the total phi excursion over the whole admissible range is < 2.02 scale
      -> the K = 5 prediction is INSENSITIVE to the unmeasured region. C1a(ii)'s caveat
         lifts with no new cells and the saturation claim stands as measured.
  (B) the excursion is >= 2.02 but the SIGN cannot reach the measurement
      -> the caveat is real but bounded; report the bound.
  (C) the excursion is large enough to move the K5-K4 increment across zero
      -> the block is necessary, and PREREG_UGAP is written.

WHAT "ADMISSIBLE RANGE" MEANS, and why two of them are reported
---------------------------------------------------------------
`mixture_model.check_limits` proves a(1) = b(1) = 1 (equal sizes, no separation) and
a(0) = b(0) = 0 (infinite separation, ideal Furnas filling). Both functions are also
required MONOTONE in u by `grade_mixture_ext`'s own printed check. So a value inserted at
u* = 0.70711 is admissible only in [a(0.667), 1] -- the MONOTONE range. That is the
range the physics permits. The UNCONSTRAINED range [0, 1] is swept too and reported
separately, because a bound that holds over [0, 1] needs no monotonicity argument at all
and is the stronger statement when it happens to hold.

WHAT IS HELD FIXED
------------------
Everything except the inserted knot: the same beta, the same measured knots at
u = 0.250/0.333/0.500/0.667, the same proved limits at u = 0 and 1, the same PCHIP
construction, the same cells. Inserting a knot DOES perturb the interpolant on the
neighbouring intervals as well -- that is not a bug, it is exactly what measuring there
would do -- so K = 2, 3, 4 are swept too as a control: their pairs all sit at u <= 0.63,
so a large move there would mean the sweep is leaking rather than measuring.

REPRODUCTION FIRST (lessons.md 2026-09-02, 2026-09-08): every ADDtab recomputed here is
checked against the value on disk in out/kladder_cells.csv before any sweep number is
printed. A diagnostic that cannot reproduce the artifact it audits reports bugs as
findings.

Run:  python diag_ugap.py      (writes out/diag_ugap.txt; exit code is the verdict gate)
"""

from __future__ import annotations

import csv
import os
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator

import grade_mixture_ext as gx
import kladder as kl
import mixture_closure as mc
import mixture_model as mm

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = kl.PHI_SCALE
MU = 0.3

# The bars. BOTH are TAB's own in-sample binary residual statistics, quoted from
# CLAIMS_KLADDER_2026-09-15.md / the ternary block. NOT MOVABLE, NOT NEW.
BAR_RMS = 2.02        # TAB in-sample rms, in scale units
BAR_MAX = 5.74        # TAB in-sample max
SEP_BAR = 1.6         # the project's standing "is this test informative" separation bar

U_GAP = 4.0 ** (-0.25)        # 0.70711 -- K = 5's adjacent-pair ratio on a span-4 ladder
NGRID = 41                    # sweep resolution per parameter


def knots(tab: dict[float, tuple[float, float]]):
    """(u, a, b) arrays of the shipped closure: measured knots + the two proved limits."""
    us = [0.0, 1.0] + sorted(tab)
    aa = [0.0, 1.0] + [tab[u][0] for u in sorted(tab)]
    bb = [0.0, 1.0] + [tab[u][1] for u in sorted(tab)]
    o = np.argsort(us)
    return np.asarray(us)[o], np.asarray(aa)[o], np.asarray(bb)[o]


def interp_with(uu, aa, bb, u_star: float, a_star, b_star):
    """Rebuild (a(u), b(u)) with an extra knot at u_star. None leaves that one alone."""
    ua, va = list(uu), list(aa)
    ub, vb = list(uu), list(bb)
    if a_star is not None:
        ua.append(u_star)
        va.append(float(a_star))
    if b_star is not None:
        ub.append(u_star)
        vb.append(float(b_star))
    oa, ob = np.argsort(ua), np.argsort(ub)
    fa = PchipInterpolator(np.asarray(ua)[oa], np.asarray(va)[oa])
    fb = PchipInterpolator(np.asarray(ub)[ob], np.asarray(vb)[ob])
    return fa, fb


def main() -> int:
    lines: list[str] = []

    def say(s: str = "") -> None:
        lines.append(s)
        print(s, flush=True)

    meas = mm.measured(mu_p=MU)
    beta1500 = meas[("lognormal", 0.0)][0]
    allpts = gx.ext_points(beta1500)
    tab, fa0, fb0 = gx.tabulated_closure(beta1500, allpts)
    uu, aa, bb = knots(tab)

    say("=== diag_ugap.py -- is the K = 5 prediction sensitive to the unmeasured "
        "u > 0.667 region? ===")
    say(f"   beta(N=1500) = {beta1500:.6f} (measured, mu_p {MU}); PHI_SCALE = {SCALE:.1e}")
    say(f"   bars: TAB in-sample rms {BAR_RMS}, max {BAR_MAX} (quoted, NOT new); "
        f"separation bar {SEP_BAR}x")
    say(f"   u* = 4^(-1/4) = {U_GAP:.5f}")
    say()
    say("   the shipped closure's knots (identifiable binary ratios + two PROVED limits):")
    say(f"      {'u':>9}{'a':>10}{'b':>10}   source")
    for u, a, b in zip(uu, aa, bb):
        src = "proved limit" if u in (0.0, 1.0) else f"measured, r = {1.0/u:.2f}"
        say(f"      {u:9.5f}{a:10.4f}{b:10.4f}   {src}")
    say(f"   -> the gap 0.667 < u < 1 contains NO knot. Width in u: "
        f"{1.0 - max(tab):.4f}, the widest interval in the construction.")
    say()

    # ------------------------------------------------------------------ the cells
    cells, _ = kl.build_cells(kl.N_PROD)
    by_label = {c["label"]: c for c in cells}
    cells3k, _ = kl.build_cells(kl.N_BIG)
    for c in cells3k:
        by_label.setdefault(c["label"], c)

    disk: dict[tuple[int, str], dict] = {}
    with open(os.path.join(OUT, "kladder_cells.csv"), encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            disk[(int(row["N"]), row["label"])] = row
    say(f"=== REPRODUCTION: {len(disk)} rows on disk in out/kladder_cells.csv ===")

    worst = 0.0
    nrep = 0
    for (n, label), row in sorted(disk.items()):
        c = by_label.get(label)
        if c is None:
            say(f"   MISSING from build_cells: {label}")
            return 2
        beta = float(row["beta_used"])
        got = float(mc.packing_density_tabulated(
            np.asarray(c["d_rel"], float), np.asarray(c["y"], float), beta, fa0, fb0))
        d = abs(got - float(row["ADDtab"]))
        worst = max(worst, d)
        nrep += 1
    say(f"   recomputed ADDtab on all {nrep} rows; worst |diff| = {worst:.3e} "
        f"({worst/SCALE:.2e} scale)")
    if worst > 5e-5:
        say("   -> REPRODUCTION FAILED. No sweep number is printed; fix the audit first.")
        return 2
    say("   -> reproduces the artifact. The sweep below is of the SHIPPED closure.")
    say()

    # the pair table, so "4 of 10 pairs sit at u*" is printed rather than asserted
    say("=== pair inventory: where each K's pairs actually sit in u ===")
    say(f"   {'K':>3}{'n_pairs':>9}   u values (count)")
    for k in kl.K_LEVELS:
        d = kl.ladder(k)
        us_k = sorted(round(d[i] / d[j], 6)
                      for j in range(len(d)) for i in range(j))
        cnt: dict[float, int] = {}
        for u in us_k:
            cnt[u] = cnt.get(u, 0) + 1
        say(f"   {k:3d}{len(us_k):9d}   "
            + ", ".join(f"{u:.5f} (x{c})" for u, c in sorted(cnt.items())))
    say(f"   -> only K = 5 has any pair above u = {max(tab):.3f}: "
        f"4 of its 10, at u* = {U_GAP:.5f}.")
    say()

    # ------------------------------------------------------- the admissible ranges
    u_lo = max(tab)
    a_lo, b_lo = tab[u_lo]
    say("=== the two ranges swept ===")
    say("   MONOTONE  (what the proved limits + the monotonicity of a, b permit):")
    say(f"      a(u*) in [{a_lo:.4f}, 1.0000]   (a({u_lo:.3f}) = {a_lo:.4f}, a(1) = 1)")
    say(f"      b(u*) in [{b_lo:.4f}, 1.0000]   (b({u_lo:.3f}) = {b_lo:.4f}, b(1) = 1)")
    say("   UNCONSTRAINED  a(u*), b(u*) in [0, 1] -- no monotonicity assumed.")
    say(f"   shipped interpolated values at u*: a {float(fa0(U_GAP)):.4f}  "
        f"b {float(fb0(U_GAP)):.4f}")
    say()

    targets = [(n, f"K{k}_yS{ys:.2f}")
               for n in (1500, 3000) for k in (2, 3, 4, 5) for ys in kl.YS_BALANCED]
    targets = [t for t in targets if t in disk]

    results: dict[str, dict[tuple[int, str], dict]] = {}
    for tag, (alo, ahi, blo, bhi) in (
            ("MONOTONE", (a_lo, 1.0, b_lo, 1.0)),
            ("UNCONSTRAINED", (0.0, 1.0, 0.0, 1.0))):
        ga = np.linspace(alo, ahi, NGRID)
        gb = np.linspace(blo, bhi, NGRID)
        res: dict[tuple[int, str], dict] = {}
        for (n, label) in targets:
            c = by_label[label]
            beta = float(disk[(n, label)]["beta_used"])
            d = np.asarray(c["d_rel"], float)
            y = np.asarray(c["y"], float)
            base = float(disk[(n, label)]["ADDtab"])
            phi_meas = float(disk[(n, label)]["phi"])

            vals_a, vals_b, vals_j = [], [], []
            for va in ga:
                f_a, f_b = interp_with(uu, aa, bb, U_GAP, va, None)
                vals_a.append(mc.packing_density_tabulated(d, y, beta, f_a, f_b))
            for vb in gb:
                f_a, f_b = interp_with(uu, aa, bb, U_GAP, None, vb)
                vals_b.append(mc.packing_density_tabulated(d, y, beta, f_a, f_b))
            for va in ga:
                for vb in gb:
                    f_a, f_b = interp_with(uu, aa, bb, U_GAP, va, vb)
                    vals_j.append(mc.packing_density_tabulated(d, y, beta, f_a, f_b))
            va_ = np.asarray(vals_a, float)
            vb_ = np.asarray(vals_b, float)
            vj_ = np.asarray(vals_j, float)
            res[(n, label)] = {
                "base": base, "phi": phi_meas, "beta": beta,
                "span_a": float(va_.max() - va_.min()) / SCALE,
                "span_b": float(vb_.max() - vb_.min()) / SCALE,
                "span_j": float(vj_.max() - vj_.min()) / SCALE,
                "jmin": float(vj_.min()), "jmax": float(vj_.max()),
                "R_base": (phi_meas - base) / SCALE,
                "R_jmin": (phi_meas - float(vj_.max())) / SCALE,
                "R_jmax": (phi_meas - float(vj_.min())) / SCALE,
            }
        results[tag] = res

        say(f"=== SWEEP -- {tag} range, {NGRID}x{NGRID} grid ===")
        say("   phi excursion in SCALE units as a(u*), b(u*) move over the whole range.")
        say(f"   {'N':>5}{'cell':>13}{'a only':>9}{'b only':>9}{'joint':>9}"
            f"{'R base':>9}{'R range':>18}")
        for (n, label) in targets:
            v = res[(n, label)]
            say(f"   {n:5d}{label:>13}{v['span_a']:9.2f}{v['span_b']:9.2f}"
                f"{v['span_j']:9.2f}{v['R_base']:+9.2f}"
                f"   [{v['R_jmin']:+7.2f}, {v['R_jmax']:+7.2f}]")
        k5 = [v for (n, l), v in res.items() if l.startswith("K5_")]
        k34 = [v for (n, l), v in res.items() if not l.startswith("K5_")]
        say(f"   worst joint excursion, K = 5 cells:      "
            f"{max(v['span_j'] for v in k5):.2f} scale")
        say(f"   worst joint excursion, K = 2/3/4 control: "
            f"{max(v['span_j'] for v in k34):.2f} scale   "
            f"(their pairs all sit at u <= {4.0**(-1/3):.3f})")
        say()

    # ------------------------------------------------------------- the K5-K4 increment
    say("=== the quantity the claim actually rests on: R(5) - R(4) at matched y_S ===")
    say("   C1's saturation statement is that this increment is <= 0. If the whole")
    say("   admissible sweep cannot carry it across zero, the unmeasured region cannot")
    say("   overturn the claim, whatever a(u*) is.")
    verdicts: list[str] = []
    for tag in ("MONOTONE", "UNCONSTRAINED"):
        res = results[tag]
        say(f"   --- {tag} ---")
        say(f"   {'N':>5}{'y_S':>7}{'R5-R4 base':>12}{'min':>9}{'max':>9}"
            f"{'crosses 0?':>12}")
        for n in (1500, 3000):
            for ys in kl.YS_BALANCED:
                k4 = res.get((n, f"K4_yS{ys:.2f}"))
                k5 = res.get((n, f"K5_yS{ys:.2f}"))
                if k4 is None or k5 is None:
                    continue
                base = k5["R_base"] - k4["R_base"]
                # R = (phi - ADD)/scale, so R is DECREASING in ADD; the extreme
                # increments pair the extremes of the two cells' ADD ranges.
                lo = (k5["phi"] - k5["jmax"]) / SCALE - (k4["phi"] - k4["jmin"]) / SCALE
                hi = (k5["phi"] - k5["jmin"]) / SCALE - (k4["phi"] - k4["jmax"]) / SCALE
                cross = "YES" if lo < 0.0 < hi else "no"
                verdicts.append(f"{tag} N{n} yS{ys}: {cross}")
                say(f"   {n:5d}{ys:7.2f}{base:+12.2f}{lo:+9.2f}{hi:+9.2f}{cross:>12}")
        say()

    # --------------------------------------------------------------------- verdict
    say("=== VERDICT ===")
    mono_worst = max(v["span_j"] for (n, l), v in results["MONOTONE"].items()
                     if l.startswith("K5_"))
    unc_worst = max(v["span_j"] for (n, l), v in results["UNCONSTRAINED"].items()
                    if l.startswith("K5_"))
    say(f"   worst K = 5 joint phi excursion: MONOTONE {mono_worst:.2f} scale, "
        f"UNCONSTRAINED {unc_worst:.2f} scale")
    say(f"   bar: TAB's own in-sample rms, {BAR_RMS} scale")

    # candidate separation: how far apart are the two extreme closures' PREDICTIONS,
    # measured against the seed error the cells carry. Below SEP_BAR the sweep itself
    # cannot discriminate and is reported NOT INFORMATIVE rather than as a null.
    ses = []
    for (n, label) in targets:
        if label.startswith("K5_"):
            ses.append(float(disk[(n, label)]["seR"]))
    se_typ = float(np.median(ses))
    sep_mono = mono_worst / se_typ if se_typ else float("inf")
    sep_unc = unc_worst / se_typ if se_typ else float("inf")
    say(f"   candidate separation (excursion / median seR = {se_typ:.2f}): "
        f"MONOTONE {sep_mono:.2f}x, UNCONSTRAINED {sep_unc:.2f}x; bar {SEP_BAR}x")

    rc = 0
    if mono_worst < BAR_RMS:
        say("   -> (A) INSENSITIVE on the MONOTONE range: the whole admissible excursion")
        say("      is below TAB's own in-sample rms. The K = 5 prediction does not")
        say("      depend on the unmeasured region at a level this campaign can see.")
    else:
        say("   -> the MONOTONE excursion EXCEEDS the in-sample rms; the region matters.")
        rc = 1
    if sep_mono < SEP_BAR:
        say(f"   -> and the sweep's own candidate separation is {sep_mono:.2f}x, BELOW the")
        say(f"      {SEP_BAR}x bar: a new binary block at u ~ 0.7-0.8 could not change")
        say("      the K = 5 verdict by more than the seed error already on it. NOT")
        say("      INFORMATIVE to run one for THIS purpose.")
    else:
        say(f"   -> candidate separation {sep_mono:.2f}x >= {SEP_BAR}x: a measurement at")
        say("      u* would be informative.")
    crossings = [v for v in verdicts if v.endswith("YES")]
    say(f"   increments that the sweep can carry across zero: {len(crossings)} of "
        f"{len(verdicts)}")
    for v in crossings:
        say(f"      {v}")
    say()
    say("   WHAT THIS CANNOT DECIDE: it sweeps the closure's PARAMETER, not the")
    say("   closure's FORM. A measurement at u* that fell OUTSIDE [0, 1], or that broke")
    say("   monotonicity, is outside every range swept here -- but both are excluded by")
    say("   mixture_model.check_limits, not by this script. And it says nothing about")
    say("   whether a(u) is measurable there, only whether K = 5 cares.")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "diag_ugap.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return rc


if __name__ == "__main__":
    sys.exit(main())
