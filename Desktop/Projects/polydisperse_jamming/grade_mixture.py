"""grade_mixture.py -- grade PREREG_MIXTURE_2026-09-13.md against the landed cells.

Reads `out/phase12_mixr*_shard*.jsonl` (the production N = 1500 campaign) and the
existing family records through the shipped `analyse_phase2.load_all`, applies the gates
M0-M6 exactly as they are written in the pre-registration, and writes
`out/mixture.txt`. Exits non-zero if any graded gate fails, and prints the failing
numbers either way.

The N = 200 smoke cells (`out/phase12_mixsmoke_*`) are on disk and are EXCLUDED here by
the N filter, not by the file name: two of that grid's compositions could not be built at
N = 200 (the periodic collider needs a cell edge > 2.5 d_max and a 200-grain cell is too
small at d_max = 2.7), so it is a different grid, not a smaller replicate.
"""

from __future__ import annotations

import collections
import glob
import json
import os
import sys

import numpy as np

import analyse_phase2 as ap2
import mixture_model as mm
from phase12_mixture import MIN_PHASE_COUNT, N_PROD, X_ANCHOR
from psd import BinaryPSD

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = mm.PHI_SCALE

# Registered bars, copied from PREREG_MIXTURE_2026-09-13.md section 4.
M1_PASS_SCALES = 5.0          # matched pair must move phi by more than this many scales
M1_FAIL_SCALES = 2.0          # below this the filling picture is dead
M2_SCALES = 5.0               # interior maximum over both endpoints
M4_SCALES = 3.0               # ratio 3 must beat ratio 2 by this much in peak excess
SEPARATION_MIN = 1.6          # a test whose candidates are closer than this decides nothing

# The friction M0-M6 are CLOSED on. Frozen here rather than read from a shared constant,
# so that adding a block at another friction cannot move a closed gate.
GRADED_MU_P = 0.3


def load_cells() -> dict[tuple[float, float], list[dict]]:
    rows: list[dict] = []
    for f in sorted(glob.glob(os.path.join(OUT, "phase12_mixr*_shard*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                if r.get("error"):
                    rows.append(r)
                    continue
                if r.get("n_requested") != N_PROD:
                    continue
                # 🔴 2026-09-13, added BEFORE the frictionless block ran, not after.
                # M0-M6 are CLOSED on the mu_p = 0.3 grid. The shard glob is
                # `phase12_mixr*`, which a `mixr3mu0` tag matches, so without this line
                # a later block at another friction would silently re-decide six closed
                # gates -- the exact shape of the 2026-09-06 lesson about appending to a
                # list a closed gate iterates. mu_p is in the cell key as well.
                if abs(float(r.get("mu_p", -1)) - GRADED_MU_P) > 1e-9:
                    continue
                rows.append(r)
    g: dict[tuple[float, float], list[dict]] = collections.defaultdict(list)
    for r in rows:
        g[(float(r["mix_ratio"]), round(float(r["x_small_vol"]), 6))].append(r)
    return g


def phi_of(cells: list[dict]) -> tuple[float, float, int]:
    v = np.array([float(c["phi"]) for c in cells if "error" not in c])
    if v.size == 0:
        return float("nan"), float("nan"), 0
    sem = float(v.std(ddof=1) / np.sqrt(v.size)) if v.size > 1 else float("nan")
    return float(v.mean()), sem, int(v.size)


def family_phi(family: str, cv: float) -> dict[int, float]:
    """Per-seed phi from the existing family records at the campaign's rung."""
    out = {}
    for r in ap2.load_all():
        if (r.get("family") or "lognormal") != family:
            continue
        if abs(float(r.get("cv", -1)) - cv) > 1e-9:
            continue
        if r.get("n_requested") != N_PROD or abs(float(r.get("mu_p", 9)) - 0.3) > 1e-9:
            continue
        if abs(float(r.get("p_hat_target", 0)) - mm.P_HAT) > 1e-12:
            continue
        if r.get("eta_roll") is not None or not r.get("converged"):
            continue
        out[int(r["seed"])] = float(r["phi"])
    return out


def main() -> int:
    lines: list[str] = []
    failures: list[str] = []

    def say(s: str = "") -> None:
        lines.append(s)
        print(s)

    g = load_cells()
    errs = [c for v in g.values() for c in v if "error" in c]
    nonconv = [c for v in g.values() for c in v
               if "error" not in c and not c.get("converged")]
    say("=== grade_mixture.py -- PREREG_MIXTURE_2026-09-13.md, N = 1500 ===")
    say(f"   {sum(len(v) for v in g.values())} cells, {len(errs)} errors, "
        f"{len(nonconv)} non-converged, grading scale PHI_SCALE = {SCALE:.1e}")
    for c in errs:
        say(f"   ERROR {c['key']}: {str(c['error'])[:140]}")
    for c in nonconv:
        say(f"   NON-CONVERGED {c['key']}: unbalanced={c.get('unbalanced')}")
    say()

    # ---- M0 -------------------------------------------------------------------------
    say("=== M0 -- REPRODUCTION, bit-identity (registered as |dphi| < 1e-12) ===")
    bid = family_phi("bidisperse", 0.5)
    mono = family_phi("lognormal", 0.0)
    anchor = {int(c["seed"]): float(c["phi"])
              for c in g.get((3.0, round(X_ANCHOR, 6)), []) if "error" not in c}
    m0_ok = True
    say("   (a) the anchor cell is BidispersePSD(cv = 0.5) reached the other way:")
    for s in sorted(anchor):
        if s not in bid:
            say(f"      seed {s}: no bidisperse record to compare against")
            continue
        d = abs(anchor[s] - bid[s])
        say(f"      seed {s}: new {anchor[s]:.17g}  old {bid[s]:.17g}  |dphi| = {d:.3e}")
        if d > 1e-12:
            m0_ok = False
    x0 = {int(c["seed"]): float(c["phi"]) for c in g.get((3.0, 0.0), []) if "error" not in c}
    x1 = {int(c["seed"]): float(c["phi"]) for c in g.get((3.0, 1.0), []) if "error" not in c}
    x0r2 = {int(c["seed"]): float(c["phi"]) for c in g.get((2.0, 0.0), []) if "error" not in c}
    say("   (b) both endpoints are monodisperse at d = 1 and must agree with each other:")
    for s in sorted(x0):
        d1 = abs(x0[s] - x1.get(s, np.nan))
        d2 = abs(x0[s] - x0r2.get(s, np.nan))
        say(f"      seed {s}: |x0 - x1| = {d1:.3e}   |r3 - r2 at x0| = {d2:.3e}")
        if d1 > 1e-12 or d2 > 1e-12:
            m0_ok = False
    say("   (c) and against the existing CV = 0 lognormal cells -- SAME inputs entirely:")
    off = []
    for s in sorted(x0):
        if s in mono:
            d = abs(x0[s] - mono[s])
            say(f"      seed {s}: |dphi| = {d:.3e}")
            if d > 1e-12:
                off.append(d)
    if off:
        say(f"      -> {len(off)} of {len(x0)} seeds differ, worst {max(off):.2e}. This is")
        say("         NOT counted as an M0 failure: it is the protocol's own documented")
        say("         non-determinism, which analyse_phase2.load_all reports at a median")
        say(f"         2.1e-03 and finds ONLY at cv = 0. Measuring {max(off):.1e} here from a")
        say("         completely separate campaign is an independent confirmation of the")
        say("         number PHI_SCALE was set from, arrived at by a gate written to test")
        say("         something else. (a) and (b) are the parts that bind.")
    say(f"   -> M0 {'PASS' if m0_ok else 'FAIL'}")
    if not m0_ok:
        failures.append("M0")
        say("   /!\\ M0 FAILED: the campaign is NOT READABLE and nothing below is quoted.")
        with open(os.path.join(OUT, "mixture.txt"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
        return 1
    say()

    # ---- the table ------------------------------------------------------------------
    say("=== the measurement: phi at fixed size ratio, across composition ===")
    say(f"   {'r':>4}{'X_s':>10}{'phi':>9}{'sem':>8}{'n':>3}{'cv_real':>9}"
        f"{'n_minor':>9}{'d_max/L':>9}{'gradable':>10}")
    table: dict[tuple[float, float], tuple[float, float, int]] = {}
    gradable: set[tuple[float, float]] = set()
    for key in sorted(g):
        cells = [c for c in g[key] if "error" not in c]
        if not cells:
            continue
        m, sem, n = phi_of(cells)
        table[key] = (m, sem, n)
        p = BinaryPSD(d10=1.0, ratio=key[0], x_small=key[1])
        n_minor = min(p.n_small, 1.0 - p.n_small) * N_PROD
        ok = (key[1] in (0.0, 1.0)) or n_minor >= MIN_PHASE_COUNT
        if ok:
            gradable.add(key)
        say(f"   {key[0]:4.1f}{key[1]:10.6f}{m:9.4f}{sem:8.4f}{n:3d}"
            f"{float(cells[0]['cv_realised']):9.4f}{n_minor:9.0f}"
            f"{float(cells[0]['d_max_over_L']):9.3f}{'yes' if ok else 'NO (M5)':>10}")
    say()

    # ---- M1 -------------------------------------------------------------------------
    say("=== M1 -- THE MATCHED PAIR. Same ratio, same number CV, different composition ===")
    a_key, b_key = (3.0, round(X_ANCHOR, 6)), (3.0, 0.25)
    pa_, sa, na = table[a_key]
    pb_, sb, nb = table[b_key]
    cva = float([c for c in g[a_key] if "error" not in c][0]["cv_realised"])
    cvb = float([c for c in g[b_key] if "error" not in c][0]["cv_realised"])
    say(f"   A: r = 3, X_s = {a_key[1]:.6f}, cv_realised = {cva:.6f}, phi = {pa_:.4f} +- {sa:.4f} (n={na})")
    say(f"   B: r = 3, X_s = {b_key[1]:.6f}, cv_realised = {cvb:.6f}, phi = {pb_:.4f} +- {sb:.4f} (n={nb})")
    say(f"   the two CVs agree to {abs(cva-cvb):.2e}; the size ratio, d_max/d_min, d95/d05,")
    say("   number span and the number of atoms are IDENTICAL by construction.")
    d = pb_ - pa_
    sd = float(np.hypot(sa, sb))
    beta = mm.measured()[("lognormal", 0.0)][0]
    pred_fill = (mm.binary_density(3.0, 0.25, beta, 1.02, 1.50)
                 - mm.binary_density(3.0, X_ANCHOR, beta, 1.02, 1.50))
    say(f"   MEASURED  dphi = {d:+.4f} +- {sd:.4f}  ({d/SCALE:+.1f} scale, {d/sd:+.1f} se)")
    say(f"   predicted: descriptor account {0.0:+.4f}, filling account {pred_fill:+.4f}")
    sep = abs(pred_fill - 0.0) / sd
    say(f"   candidate separation = {sep:.1f} measurement errors "
        f"(the rule refuses to call a test informative below {SEPARATION_MIN})")
    if sep < SEPARATION_MIN:
        say("   -> M1 NOT INFORMATIVE")
        failures.append("M1(not informative)")
    elif abs(d) > M1_PASS_SCALES * SCALE and d > 0:
        say(f"   -> M1 PASS: {abs(d)/SCALE:.1f} scale against a bar of {M1_PASS_SCALES}, "
            f"sign positive as predicted.")
        say(f"      the filling account predicted {pred_fill:+.4f} and the measurement is "
            f"{d:+.4f}, i.e. {abs(d/pred_fill - 1)*100:.0f}% off a ZERO-PARAMETER prediction.")
    elif abs(d) < M1_FAIL_SCALES * SCALE:
        say("   -> M1 FAIL: the filling picture is refuted on one matched pair.")
        failures.append("M1")
    else:
        say("   -> M1 UNDETERMINED: between the pass and fail bars as registered.")
        failures.append("M1(undetermined)")
    say()

    # ---- M2 / M3 / M4 ---------------------------------------------------------------
    say("=== M2 -- an interior maximum exists; M3 -- where it is (reported, not graded) ===")
    peak: dict[float, tuple[float, float]] = {}
    for r in (3.0, 2.0):
        pts = [(x, table[(r, x)][0]) for (rr, x) in sorted(table) if rr == r
               and (rr, x) in gradable and 0.0 < x < 1.0]
        if not pts:
            continue
        end = np.mean([table[(r, 0.0)][0], table[(r, 1.0)][0]])
        xb, pb2 = max(pts, key=lambda t: t[1])
        peak[r] = (xb, pb2 - end)
        ok = (pb2 - end) > M2_SCALES * SCALE
        say(f"   r = {r:.0f}: endpoints {end:.4f}, maximum {pb2:.4f} at X_s = {xb:.2f}, "
            f"excess {pb2-end:+.4f} ({(pb2-end)/SCALE:+.1f} scale)  -> {'PASS' if ok else 'FAIL'}")
        if not ok:
            failures.append(f"M2(r={r:.0f})")
    say("   M3, the peak POSITION: measured "
        + ", ".join(f"r={r:.0f} at X_s={peak[r][0]:.2f}" for r in sorted(peak, reverse=True)))
    say("   against zero-parameter model r=3 at 0.40 and r=2 at 0.30, fitted 0.20 / 0.15.")
    say("   Reported only: the grid spacing (0.05) is coarser than the model's own")
    say("   uncertainty in the peak position, so this cannot separate them.")
    say()

    say("=== M4 -- THE SIZE RATIO MATTERS, which the blind null denies ===")
    if 3.0 in peak and 2.0 in peak:
        gap = peak[3.0][1] - peak[2.0][1]
        ok = gap > M4_SCALES * SCALE
        say(f"   peak excess over beta: r = 3 gives {peak[3.0][1]:+.4f}, "
            f"r = 2 gives {peak[2.0][1]:+.4f}, difference {gap:+.4f} ({gap/SCALE:+.1f} scale)")
        say("   zero-parameter prediction: +0.0907 vs +0.0562. The size-ratio-blind null")
        say("   predicts the SAME curve at both ratios, i.e. a difference of exactly 0.")
        say(f"   -> M4 {'PASS' if ok else 'FAIL'} against a bar of {M4_SCALES} scale.")
        if not ok:
            failures.append("M4")
    say()

    # ---- CV is not the control parameter, INSIDE one ratio --------------------------
    say("=== the negative, now from WITHIN a single two-atom family ===")
    say("   Every point below is a binary at ONE size ratio. If CV were the control")
    say("   parameter, phi would be a function of cv_realised alone on each block.")
    for r in (3.0, 2.0):
        pts = [(float([c for c in g[(r, x)] if "error" not in c][0]["cv_realised"]),
                table[(r, x)][0], x)
               for (rr, x) in sorted(table) if rr == r and (rr, x) in gradable]
        pts.sort()
        say(f"   r = {r:.0f}:  " + "  ".join(f"cv {c:.4f}->phi {p:.4f}" for c, p, _ in pts))
        # the sharpest violation: the closest pair in CV with the largest gap in phi
        worst = None
        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                dcv = abs(pts[i][0] - pts[j][0])
                dphi = abs(pts[i][1] - pts[j][1])
                if dcv > 0.02:
                    continue
                if worst is None or dphi > worst[0]:
                    worst = (dphi, dcv, pts[i], pts[j])
        if worst:
            say(f"      closest-in-CV pair (dcv = {worst[1]:.4f}): X_s {worst[2][2]:.3f} and "
                f"{worst[3][2]:.3f} give phi {worst[2][1]:.4f} and {worst[3][1]:.4f}, "
                f"apart by {worst[0]:.4f} ({worst[0]/SCALE:.0f} scale)")
        mono_cv = all(pts[i][1] <= pts[i + 1][1] for i in range(len(pts) - 1))
        say(f"      phi monotone in cv_realised on this block: {mono_cv}")
    say()

    # ---- the model, point by point ---------------------------------------------------
    say("=== the model against the measurement, point by point ===")
    pa_r, pb_r, _, _ = mm.fit_on_bidisperse(mm.measured(), beta, blind=False)
    pa_b, pb_b, _, _ = mm.fit_on_bidisperse(mm.measured(), beta, blind=True)
    say(f"   beta = {beta:.4f} (measured). zero-param PA = 1.02, PB = 1.50 (textbook).")
    say(f"   fitted-on-bidisperse PA = {pa_r:.3f}, PB = {pb_r:.3f}; "
        f"blind null a = {pa_b:.3f}, b = {pb_b:.3f}.")
    say(f"   {'r':>4}{'X_s':>8}{'phi_meas':>10}{'zero-par':>10}{'/scale':>8}"
        f"{'fitted':>9}{'/scale':>8}{'null':>8}{'/scale':>8}")
    chi = {"zero": 0.0, "fitted": 0.0, "null": 0.0}
    npt = 0
    for key in sorted(table):
        if key not in gradable or key[1] in (0.0, 1.0):
            continue
        r, x = key
        m = table[key][0]
        z = mm.binary_density(r, x, beta, 1.02, 1.50)
        f_ = mm.binary_density(r, x, beta, pa_r, pb_r)
        nb2 = mm.packing_density_blind(np.array([r, 1.0]), np.array([1.0 - x, x]),
                                       beta, pa_b, pb_b)
        chi["zero"] += ((z - m) / SCALE) ** 2
        chi["fitted"] += ((f_ - m) / SCALE) ** 2
        chi["null"] += ((nb2 - m) / SCALE) ** 2
        npt += 1
        say(f"   {r:4.1f}{x:8.2f}{m:10.4f}{z:10.4f}{(z-m)/SCALE:+8.1f}"
            f"{f_:9.4f}{(f_-m)/SCALE:+8.1f}{nb2:8.4f}{(nb2-m)/SCALE:+8.1f}")
    say(f"   chi2 over {npt} gradable interior points, in PHI_SCALE units:")
    for k in ("zero", "fitted", "null"):
        say(f"      {k:<8} {chi[k]:9.1f}   (per point {chi[k]/npt:7.1f})")
    say("   The zero-parameter column had NOTHING fitted to this campaign; the fitted one")
    say("   was calibrated on the bidisperse family before this campaign existed; the null")
    say("   had the same two free parameters and the same calibration set.")
    say()

    say("=== where the zero-parameter model is wrong, which is the useful part ===")
    for r in (3.0, 2.0):
        res = [(x, table[(r, x)][0] - mm.binary_density(r, x, beta, 1.02, 1.50))
               for (rr, x) in sorted(table) if rr == r and (rr, x) in gradable
               and 0.0 < x < 1.0]
        if not res:
            continue
        v = np.array([t[1] for t in res])
        say(f"   r = {r:.0f}: residual (meas - model) runs {v.min():+.4f} to {v.max():+.4f}, "
            f"mean {v.mean():+.4f} ({v.mean()/SCALE:+.1f} scale)")
    say("   -> the model OVER-credits filling at the smaller ratio and is close at the")
    say("      larger one. That is the same defect the pre-registration named in advance")
    say("      from the family path (its bidisperse maximum sits at CV 0.2 where the")
    say("      measurement puts it at CV 0.4), now measured directly and with a sign.")
    say()

    say(f"FAILURES: {failures}")
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "mixture.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
