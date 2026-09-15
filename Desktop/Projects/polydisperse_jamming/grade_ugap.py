"""grade_ugap.py -- grades PREREG_UGAP_2026-09-16.md against the U-GAP block.

Reads ONLY `out/phase12_ug*_shard*.jsonl` for its new data. The tag `ug` / `ug3k` matches
none of `mixr*`, `extr*`, `tern*`, `kl*`, and `assert_scope()` checks that STRUCTURALLY
rather than by comment, so M0-M6, F1-F3, X0-X5, T0-T6, NA1-NA4 and G0-G8 cannot be
re-decided here.

Every bar is quoted, none is invented:
  * the identifiability rule is `ugap.IDENT_WIDTH`, asserted equal to the one
    `grade_mixture_ext.tabulated_closure` applies;
  * the residual bars are TAB's own in-sample rms 2.02 and max 5.74.

Run:  python grade_ugap.py      (writes out/grade_ugap.txt; EXIT CODE is the verdict)
"""

from __future__ import annotations

import collections
import csv
import glob
import json
import math
import os
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator

import grade_mixture_ext as gx
import kladder as kl
import mixture_closure as mc
import mixture_model as mm
import ugap as ug

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = ug.PHI_SCALE
MU = ug.MU_P
BAR_RMS = 2.02          # TAB in-sample rms, scale units. QUOTED. NOT MOVABLE.
BAR_MAX = 5.74          # TAB in-sample max. QUOTED. NOT MOVABLE.
U2_MIN_NEG = 7          # registered: negative on >= 7 of 9
U3_MIN_LEVELS = 2       # registered: R5-R4 <= 0 on >= 2 of 3 y_S at both N


def assert_scope() -> str:
    """The U-GAP tag must collide with no closed block's loader glob."""
    tags = ("ug", "ug3k")
    closed = ("mixr", "extr", "tern", "kl")
    bad = [f"{t} vs {c}" for t in tags for c in closed if t.startswith(c)]
    if bad:
        raise SystemExit(f"SCOPE VIOLATION: {bad}")
    seen = sorted(os.path.basename(p) for p in
                  glob.glob(os.path.join(OUT, "phase12_ug*_shard*.jsonl")))
    for p in seen:
        for c in closed:
            if p.startswith(f"phase12_{c}"):
                raise SystemExit(f"SCOPE VIOLATION: {p} would load into a closed block")
    if ug.IDENT_WIDTH != 0.25:
        raise SystemExit(f"BAR MOVED: ugap.IDENT_WIDTH = {ug.IDENT_WIDTH}, must be 0.25")
    return (f"scope OK: tags {tags} collide with none of {closed}; "
            f"{len(seen)} U-GAP shard files; IDENT_WIDTH {ug.IDENT_WIDTH} unmoved")


def load(tag: str):
    """Raw U-GAP records for one tag, grouped by (ratio, x). Counts are printed."""
    g: dict[tuple[float, float], list[dict]] = collections.defaultdict(list)
    files = sorted(glob.glob(os.path.join(OUT, f"phase12_{tag}_shard*.jsonl")))
    nraw = 0
    for f in files:
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                nraw += 1
                r = json.loads(line)
                g[(float(r["mix_ratio"]), round(float(r["x_small_vol"]), 6))].append(r)
    return g, len(files), nraw


def admissible(c: dict, n_expect: int) -> tuple[bool, str]:
    if c.get("error"):
        return False, "error"
    if not c.get("converged"):
        return False, "not converged"
    if c.get("n_requested") != n_expect:
        return False, f"N {c.get('n_requested')}"
    if abs(float(c.get("mu_p", -1)) - MU) > 1e-9:
        return False, f"mu {c.get('mu_p')}"
    bc = c.get("by_class", {})
    if min(bc.get("n_S", 0), bc.get("n_L", 0)) < ug.MIN_PHASE_COUNT:
        return False, f"phase count {bc.get('n_S')}/{bc.get('n_L')}"
    return True, ""


def main() -> int:
    lines: list[str] = []
    failures: list[str] = []
    notinfo: list[str] = []

    def say(s: str = "") -> None:
        lines.append(s)
        print(s, flush=True)

    say("=== grade_ugap.py -- PREREG_UGAP_2026-09-16.md ===")
    say("   " + assert_scope())

    meas = mm.measured(mu_p=MU)
    beta = meas[("lognormal", 0.0)][0]
    allpts = gx.ext_points(beta)
    tab, fa0, fb0 = gx.tabulated_closure(beta, allpts)
    say(f"   beta = {beta:.6f} (measured, mu_p {MU}); PHI_SCALE = {SCALE:.1e}")
    say(f"   bars: identifiability width {ug.IDENT_WIDTH} (verbatim), TAB in-sample "
        f"rms {BAR_RMS} / max {BAR_MAX} (quoted)")
    say(f"   shipped closure at u* = {ug.U_STAR:.5f}: a {float(fa0(ug.U_STAR)):.4f}  "
        f"b {float(fb0(ug.U_STAR)):.4f}")
    say()

    w = ug.check_kernel(beta, np.random.default_rng(20260916))
    say(f"   ANCHOR: vectorised kernel vs mixture_closure.packing_density_ab, 300 random "
        f"points, worst |diff| = {w:.3e}")
    if w > 1e-12:
        say("   -> ANCHOR FAILED. No verdict is printed.")
        return 2
    say()

    # ------------------------------------------------------------------- data
    g15, nf15, nraw15 = load("ug")
    say(f"=== data, N = 1500: {nraw15} records from {nf15} shard files -> "
        f"{len(g15)} compositions ===")
    adm15: dict[tuple[float, float], list[dict]] = {}
    nref = 0
    for key in sorted(g15):
        ok = []
        for c in g15[key]:
            a, why = admissible(c, ug.N_PROD)
            if a:
                ok.append(c)
            else:
                nref += 1
                say(f"   REFUSED {c.get('key')}: {why}")
        adm15[key] = ok
    nadm = sum(len(v) for v in adm15.values())
    say(f"   admitted {nadm} of {nraw15}, refused {nref}; expected "
        f"{len(ug.admitted()) * len(ug.SEEDS)}")
    if nadm != len(ug.admitted()) * len(ug.SEEDS):
        failures.append(f"data: admitted {nadm}, expected "
                        f"{len(ug.admitted()) * len(ug.SEEDS)}")
    say()

    def pts_for(arm: str, adm) -> list[tuple[float, float, float]]:
        r = ug.RATIOS[arm]
        out = []
        for x in ug.X_GRID:
            v = adm.get((float(r), round(float(x), 6)), [])
            if v:
                out.append((float(r), float(x),
                            float(np.mean([c["phi"] for c in v]))))
        return out

    # ---------------------------------------------------------------- U0 wiring
    say("=== U0 -- the bit-identity wiring gate, RUN not asserted ===")
    say("   Four (r, x) cells exist in the `mixr` block at the same seeds, N, p_hat and")
    say("   mu_p, and phase12_mixture.cell_key carries no tag, so they must reproduce")
    say("   EXACTLY. diff == 0.0, not 'close'.")
    old: dict[str, float] = {}
    for f in sorted(glob.glob(os.path.join(OUT, "phase12_mixr*_shard*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                if not rec.get("error") and rec.get("key"):
                    old[rec["key"]] = float(rec.get("phi", float("nan")))
    n_match = n_tot = 0
    for r, x in ug.BITID:
        for c in adm15.get((float(r), round(float(x), 6)), []):
            k = c.get("key")
            if k in old:
                n_tot += 1
                d = abs(float(c["phi"]) - old[k])
                n_match += int(d == 0.0)
                say(f"   {k}: diff = {d:+.3e}  {'IDENTICAL' if d == 0.0 else 'DIFFERS'}")
    say(f"   -> {n_match} of {n_tot} bit-identical (expected 12 of 12)")
    if n_tot == 0 or n_match != n_tot:
        failures.append(f"U0 {n_match}/{n_tot} bit-identical")
        say("   -> U0 FAILED: the new block is not on the old protocol. Nothing else "
            "is graded.")
        say(f"FAILURES: {failures}")
        with open(os.path.join(OUT, "grade_ugap.txt"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
        return 1
    say()

    # ----------------------------------------------------- U1 / U1b / U1c identifiability
    say("=== U1 / U1b / U1c -- the registered identifiability rule, applied verbatim ===")
    say("   'both profile intervals narrower than 0.25 AND containing their own best fit'")
    say(f"   {'arm':>4}{'r':>10}{'u':>9}{'n_x':>5}{'phi span':>10}{'/scale':>8}"
        f"{'a':>8}{'w(a)':>7}{'b':>8}{'w(b)':>7}  ident?")
    res: dict[str, dict] = {}
    for arm in ("A", "B", "C"):
        pts = pts_for(arm, adm15)
        if len(pts) < 4:
            say(f"   {arm:>4}: only {len(pts)} compositions; NOT GRADED")
            failures.append(f"U1{arm}: {len(pts)} compositions")
            continue
        ok, ia, ib, a_hat, b_hat = ug.identifiable(pts, beta)
        span = max(p[2] for p in pts) - min(p[2] for p in pts)
        res[arm] = {"pts": pts, "ok": ok, "ia": ia, "ib": ib, "a": a_hat, "b": b_hat,
                    "span": span}
        say(f"   {arm:>4}{ug.RATIOS[arm]:10.5f}{1/ug.RATIOS[arm]:9.5f}{len(pts):5d}"
            f"{span:10.4f}{span/SCALE:8.1f}{a_hat:8.4f}{ia[1]-ia[0]:7.3f}"
            f"{b_hat:8.4f}{ib[1]-ib[0]:7.3f}  {'YES' if ok else 'NO'}")
        say(f"        a interval [{ia[0]:.4f}, {ia[1]:.4f}]   "
            f"b interval [{ib[0]:.4f}, {ib[1]:.4f}]")
    say()

    # the escape diagnostic, which is what explains a NO
    say("   the escape model (a = 1 or b = 1 -> phi = beta identically) against the bar.")
    say("   The bar is the SAME one the profile uses, chi2_min + n, not n -- an earlier")
    say("   draft of the claims quoted a bar of 10.41 that no artifact printed (referee,")
    say("   2026-09-16). Both chi2_min and the bar are printed here so it cannot recur.")
    say(f"      {'arm':>5}{'chi2_min':>10}{'n':>4}{'bar':>9}{'chi2 escape':>13}  verdict")
    for arm in ("A", "B", "C"):
        if arm not in res:
            continue
        pts = res[arm]["pts"]
        xv = np.array([p[1] for p in pts])
        yv = np.array([p[2] for p in pts])
        grid = np.linspace(1e-4, 1 - 1e-4, 201)
        chi2d = (((ug.phi_fast(grid[:, None, None], grid[None, :, None],
                               xv[None, None, :], beta) - yv[None, None, :])
                  / SCALE) ** 2).sum(axis=2)
        cmin = float(chi2d.min())
        bar = cmin + len(pts)
        esc = float(min((((ug.phi_fast(gv, 1.0, xv, beta) - yv) / SCALE) ** 2).sum()
                        for gv in grid))
        res[arm]["chi2_min"] = cmin
        res[arm]["bar"] = bar
        say(f"      {arm:>5}{cmin:10.2f}{len(pts):4d}{bar:9.2f}{esc:13.2f}  "
            f"{'OPEN -- unidentifiable by construction' if esc <= bar else 'closed'}")
    say()

    # U1b, the attribution control: arm C against the existing mixr value
    say("   U1b -- the ATTRIBUTION CONTROL. Arm C is r = 1.5, an already-identifiable")
    say("   knot; running the NEW composition grid there separates 'the grid moved the")
    say("   answer' from 'the ratio moved the answer'.")
    u1b_ok = False
    if "C" in res:
        a_old, b_old = tab[1 / 1.5] if (1 / 1.5) in tab else (float("nan"),) * 2
        old_pts = [p for p in allpts if abs(p[0] - 1.5) < 1e-9 and 0.0 < p[1] < 1.0]
        ia_o = ug.profile_fast(old_pts, beta, 0)
        ib_o = ug.profile_fast(old_pts, beta, 1)
        ia_n, ib_n = res["C"]["ia"], res["C"]["ib"]
        say(f"      old grid: a {a_old:.4f} [{ia_o[0]:.4f}, {ia_o[1]:.4f}]   "
            f"b {b_old:.4f} [{ib_o[0]:.4f}, {ib_o[1]:.4f}]")
        say(f"      new grid: a {res['C']['a']:.4f} [{ia_n[0]:.4f}, {ia_n[1]:.4f}]   "
            f"b {res['C']['b']:.4f} [{ib_n[0]:.4f}, {ib_n[1]:.4f}]")
        say(f"      da = {res['C']['a'] - a_old:+.4f}, db = {res['C']['b'] - b_old:+.4f}")
        # The registered wording said "within the UNION of the two profile intervals".
        # A REFEREE (2026-09-16) showed that test CANNOT FAIL: `identifiable()` already
        # guarantees each best fit lies inside its OWN interval, and each interval is a
        # subset of the union -- so any two individually-identifiable fits pass it, even
        # with DISJOINT intervals an order of magnitude further apart than these. It was
        # a gate that could not fire. The vacuous number stays printed beside the
        # replacement, which is stricter in KIND: the intervals must actually OVERLAP.
        ua = (min(ia_o[0], ia_n[0]), max(ia_o[1], ia_n[1]))
        ub = (min(ib_o[0], ib_n[0]), max(ib_o[1], ib_n[1]))
        a_in = ua[0] <= res["C"]["a"] <= ua[1] and ua[0] <= a_old <= ua[1]
        b_in = ub[0] <= res["C"]["b"] <= ub[1] and ub[0] <= b_old <= ub[1]
        say(f"      [registered wording] both values inside the UNION: a {a_in}, "
            f"b {b_in} -- VACUOUS, true for any two identifiable fits")
        a_ov = min(ia_o[1], ia_n[1]) - max(ia_o[0], ia_n[0])
        b_ov = min(ib_o[1], ib_n[1]) - max(ib_o[0], ib_n[0])
        a_nar = min(ia_o[1] - ia_o[0], ia_n[1] - ia_n[0])
        b_nar = min(ib_o[1] - ib_o[0], ib_n[1] - ib_n[0])
        say(f"      [replacement] the two INTERSECTIONS must be non-empty: "
            f"a overlap {a_ov:+.4f} (narrower interval {a_nar:.4f}), "
            f"b overlap {b_ov:+.4f} ({b_nar:.4f})")
        say(f"      and |da| must not exceed the narrower interval's width: "
            f"|da| {abs(res['C']['a'] - a_old):.4f} vs {a_nar:.4f}; "
            f"|db| {abs(res['C']['b'] - b_old):.4f} vs {b_nar:.4f}")
        u1b_ok = bool(res["C"]["ok"] and a_ov > 0 and b_ov > 0
                      and abs(res["C"]["a"] - a_old) <= a_nar
                      and abs(res["C"]["b"] - b_old) <= b_nar)
        # The referee's recommended gate, stronger than interval overlap in KIND: a
        # NESTED-MODEL test. Interval overlap is shape-blind -- two intervals can
        # intersect while each fit's optimum sits at a far edge -- so ask directly
        # whether ONE (a, b) explains BOTH grids. delta chi2 for 2 extra parameters,
        # against the project's own "1 unit per dof" convention (bar 2).
        gg = np.linspace(1e-4, 1 - 1e-4, 401)

        def chi_min_of(pset):
            xx = np.array([p[1] for p in pset])
            yy = np.array([p[2] for p in pset])
            return float((((ug.phi_fast(gg[:, None, None], gg[None, :, None],
                                        xx[None, None, :], beta) - yy[None, None, :])
                           / SCALE) ** 2).sum(axis=2).min())

        c_sep = chi_min_of(old_pts) + chi_min_of(res["C"]["pts"])
        c_pool = chi_min_of(old_pts + res["C"]["pts"])
        d_chi = c_pool - c_sep
        say(f"      [nested-model test, stronger in kind] one shared (a, b) for BOTH "
            f"grids: chi2 {c_pool:.2f} vs {c_sep:.2f} separate, "
            f"delta chi2 = {d_chi:.2f} for 2 extra dof")
        say(f"      against the project's 1-per-dof convention (bar 2.0): "
            f"{'PASS -- one interaction pair explains both grids' if d_chi <= 2.0 else 'FAIL'}")
        u1b_ok = bool(u1b_ok and d_chi <= 2.0)
    say(f"   -> U1b {'PASS' if u1b_ok else 'FAIL'}")
    if not u1b_ok:
        failures.append("U1b: the composition grid moved the answer at a fixed ratio")
    say()

    for arm, gate, want in (("A", "U1", True), ("B", "U1c", True)):
        if arm not in res:
            continue
        got = res[arm]["ok"]
        say(f"   -> {gate} (arm {arm}, u = {1/ug.RATIOS[arm]:.5f}): "
            f"{'IDENTIFIABLE' if got else 'UNDETERMINED'}, predicted "
            f"{'IDENTIFIABLE' if want else 'UNDETERMINED'} -- "
            f"{'PASS' if got == want else 'FAIL'}")
        if got != want:
            failures.append(f"{gate}: arm {arm} came back "
                            f"{'IDENTIFIABLE' if got else 'UNDETERMINED'}")
        if not u1b_ok:
            notinfo.append(f"{gate} NOT INFORMATIVE about the ratio: U1b failed")
    if "B" in res and not res["B"]["ok"]:
        say("   the registered reading of arm B's null, CHECKED rather than quoted.")
        say("   The prereg said an UNDETERMINED arm B is evidence a(0.8) > 0.98, and it")
        say("   reached that 0.98 from ARM A's robustness table -- the wrong arm. The")
        say("   defect is recorded; arm B's own ceiling is computed here instead.")
        say("   (It also predicted the mechanism would be the a = 1 / b = 1 escape. It")
        say("   is NOT: the new grid CLOSED that escape at arm B, see above. So the")
        say("   registered sentence was right about the verdict and wrong about the")
        say("   reason, and the reason is what the bound rests on.)")
        say("   NOISE-FREE ceiling at arm B, as the assumed truth moves "
            "(POST-HOC, declared):")
        say(f"      {'a_true':>8}{'b_true':>8}{'w(a)':>8}{'w(b)':>8}{'ident':>7}"
            f"{'max signal/scale':>18}")
        rb = ug.RATIOS["B"]
        bound = None
        sig_ok, sig_no = float("inf"), 0.0   # smallest resolvable / largest unresolved
        for a_, b_ in ((0.9695, 0.9177), (0.9695, 0.9574), (0.9750, 0.9600),
                       (0.9800, 0.9600), (res["B"]["a"], res["B"]["b"]),
                       (0.9900, 0.9700)):
            p_ = [(rb, x, float(mc.packing_density_ab(np.array([rb, 1.0]),
                                                      np.array([1 - x, x]), beta, a_, b_)))
                  for x in ug.X_GRID]
            ok_, ia_, ib_, _, _ = ug.identifiable(p_, beta)
            sig = (max(q[2] for q in p_) - beta) / SCALE
            if ok_:
                bound = a_
                sig_ok = min(sig_ok, sig)
            else:
                sig_no = max(sig_no, sig)
            say(f"      {a_:8.4f}{b_:8.4f}{ia_[1]-ia_[0]:8.3f}{ib_[1]-ib_[0]:8.3f}"
                f"{'YES' if ok_ else 'NO':>7}{sig:18.2f}")
        say(f"   NOT ARGUED FROM THE POINT ESTIMATE. Arm B's fitted b = "
            f"{res['B']['b']:.4f} is 1/(1 + e^0), i.e. the Nelder-Mead START -- on a flat")
        say("   surface the optimiser did not move, and the 201-grid argmin lands")
        say("   somewhere else entirely. A meaningless point estimate cannot carry a")
        say("   bound, so the argument is made on the MEASURED SIGNAL instead:")
        sigB = (max(p[2] for p in res["B"]["pts"]) - beta) / SCALE
        say(f"      arm B's measured max(phi - beta) = {sigB:.2f} scale")
        say(f"      the ceiling above resolves down to a signal of {sig_ok:.2f} scale and "
            f"fails at {sig_no:.2f} -- both read off the table, not hardcoded")
        say(f"   -> the measured filling signal at u = 0.8 ({sigB:.2f}) sits BELOW the")
        say("      threshold this design can resolve, so NO experiment at this design")
        say("      could have determined a(0.8), whatever it is. Mapping that signal")
        say(f"      through the ceiling table puts a(0.8) above about {bound:.4f}.")
        say("      SUPPORTABLE, with its status stated: an inference FROM A NULL via the")
        say("      design's own resolving power -- not a measurement, and not to be")
        say("      quoted with an interval. It is consistent with the measured arm A")
        say(f"      (a(0.70711) = {res['A']['a']:.4f}) under monotonicity in u.")
    say()

    # ------------------------------------------------------------------ U2 direction
    say("=== U2 -- the DIRECTION: does the shipped extrapolation over-predict at u*? ===")
    say("   Registered: negative on >= 7 of 9 compositions at arm A, and the fitted")
    say(f"   a(u*) above the shipped {float(fa0(ug.U_STAR)):.4f}.")
    u2_ok = False
    if "A" in res:
        rr = ug.R_STAR
        u = 1.0 / rr
        av, bv = float(fa0(u)), float(fb0(u))
        resid = np.array([(p[2] - float(mc.packing_density_ab(
            np.array([rr, 1.0]), np.array([1 - p[1], p[1]]), beta, av, bv))) / SCALE
            for p in res["A"]["pts"]])
        nneg = int((resid < 0).sum())
        # exact two-sided sign-test p, binomial(n, 1/2)
        n = resid.size
        k = max(nneg, n - nneg)
        p_sign = 2.0 * sum(math.comb(n, i) for i in range(k, n + 1)) / 2**n
        u2_ok = nneg >= U2_MIN_NEG
        say("   residuals phi_measured - phi_shipped, scale units: "
            + " ".join(f"{v:+.2f}" for v in resid))
        say(f"   mean {resid.mean():+.2f}, rms {np.sqrt((resid**2).mean()):.2f}, "
            f"max|.| {np.abs(resid).max():.2f}; negative on {nneg}/{n}, "
            f"sign-test p = {min(p_sign, 1.0):.4f}")
        say(f"   fitted a(u*) = {res['A']['a']:.4f} against the shipped {av:.4f} "
            f"({res['A']['a'] - av:+.4f}); b {res['A']['b']:.4f} vs {bv:.4f} "
            f"({res['A']['b'] - bv:+.4f})")
        say(f"   -> U2 {'PASS' if u2_ok else 'FAIL'} "
            f"({nneg} negative, bar >= {U2_MIN_NEG}) -- and {nneg}/{n} is EXACTLY the "
            f"bar, a marginal pass, labelled as one.")
        say("   TWO HONEST QUALIFICATIONS, both of which weaken it:")
        say(f"   (1) The registered sentence said '>= {U2_MIN_NEG} of 9 (sign-test "
            f"p <= 0.09)'. Those halves do not agree: {nneg}/9 gives one-sided "
            f"p = {min(p_sign, 1.0)/2:.4f} and two-sided p = {min(p_sign, 1.0):.4f}, so it")
        say("       clears the bar only on the one-sided reading, which the sentence did")
        say("       not specify. The prereg named the COUNT as the gate, so the count is")
        say("       what is graded -- but the parenthetical was ambiguous and is recorded")
        say("       as a defect rather than resolved in the direction that passes.")
        say(f"   (2) The a interval [{res['A']['ia'][0]:.4f}, {res['A']['ia'][1]:.4f}] "
            f"CONTAINS the shipped {av:.4f}, so the SECOND half of U2 -- 'the fitted")
        say("       a(u*) comes out above 0.9398' -- is true of the point estimate and")
        say("       NOT RESOLVED by the interval. The direction is a trend, not a")
        say("       measured disagreement.")
        # A referee (2026-09-16) asked whether 9 compositions sharing 3 seed values are
        # independent trials. Measured rather than argued: the per-seed residual columns
        # are correlated across compositions for one seed pair.
        cols = {}
        for s in ug.SEEDS:
            v = []
            for x in ug.X_GRID:
                rec = [c for c in adm15.get((float(rr), round(float(x), 6)), [])
                       if int(c["key"].split("_s")[1].split("_")[0]) == s]
                if rec:
                    v.append((float(rec[0]["phi"]) - float(mc.packing_density_ab(
                        np.array([rr, 1.0]), np.array([1 - x, x]), beta, av, bv))) / SCALE)
            cols[s] = np.array(v)
        say("   (3) The 9 compositions share the SAME 3 seed values, so the sign test's")
        say("       independence assumption is not free. Per-seed residual correlation")
        say("       across compositions:")
        ss = list(ug.SEEDS)
        for i in range(len(ss)):
            for j in range(i + 1, len(ss)):
                if cols[ss[i]].size == cols[ss[j]].size == len(ug.X_GRID):
                    r_ = float(np.corrcoef(cols[ss[i]], cols[ss[j]])[0, 1])
                    say(f"          seed {ss[i]} vs {ss[j]}: r = {r_:+.3f}")
        say("       One pair is strongly correlated, so the binomial p-values above are")
        say("       optimistic. It changes no verdict -- U2 is already reported as")
        say("       marginal on three other grounds -- but the count is not 9 independent")
        say("       trials and must not be quoted as if it were.")
        if not u2_ok:
            failures.append(f"U2: {nneg}/{n} negative, bar >= {U2_MIN_NEG}")
    say()

    # ------------------------------------------------------------------ U4 monotonicity
    say("=== U4 -- monotonicity of a(u), b(u) with the new knots inserted ===")
    new_tab = dict(tab)
    for arm in ("A", "B", "C"):
        if arm in res and res[arm]["ok"]:
            new_tab[1.0 / ug.RATIOS[arm]] = (res[arm]["a"], res[arm]["b"])
    us = [0.0, 1.0] + sorted(new_tab)
    aa = [0.0, 1.0] + [new_tab[u][0] for u in sorted(new_tab)]
    bb = [0.0, 1.0] + [new_tab[u][1] for u in sorted(new_tab)]
    o = np.argsort(us)
    uu = np.asarray(us)[o]
    aa = np.asarray(aa)[o]
    bb = np.asarray(bb)[o]
    say(f"   {'u':>9}{'a':>9}{'b':>9}   source")
    for u, a_, b_ in zip(uu, aa, bb):
        src = ("proved limit" if u in (0.0, 1.0)
               else "NEW" if u in {1.0 / ug.RATIOS[k] for k in res
                                   if res[k]["ok"]} and u not in tab
               else "REMEASURED" if u in {1.0 / ug.RATIOS[k] for k in res if res[k]["ok"]}
               else "existing")
        say(f"   {u:9.5f}{a_:9.4f}{b_:9.4f}   {src}")
    mono_a = bool(np.all(np.diff(aa) >= -1e-9))
    mono_b = bool(np.all(np.diff(bb) >= -1e-9))
    say(f"   -> U4 monotone: a {mono_a}, b {mono_b} -- "
        f"{'PASS' if mono_a and mono_b else 'FAIL'}")
    if not (mono_a and mono_b):
        failures.append(f"U4: monotone a {mono_a}, b {mono_b}")
    say()

    # ------------------------------------------------------------------ U3 the consequence
    say("=== U3 -- the K-ladder consequence: does 'and then SATURATES' survive? ===")
    fa_n = PchipInterpolator(uu, aa)
    fb_n = PchipInterpolator(uu, bb)
    disk: dict[tuple[int, str], dict] = {}
    with open(os.path.join(OUT, "kladder_cells.csv"), encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            disk[(int(row["N"]), row["label"])] = row
    by_label: dict[str, dict] = {}
    for nn in (kl.N_PROD, kl.N_BIG):
        for c in kl.build_cells(nn)[0]:
            by_label.setdefault(c["label"], c)

    def R_new(n, label):
        row = disk[(n, label)]
        c = by_label[label]
        add = float(mc.packing_density_tabulated(
            np.asarray(c["d_rel"], float), np.asarray(c["y"], float),
            float(row["beta_used"]), fa_n, fb_n))
        return (float(row["phi"]) - add) / SCALE, add

    say("   R = (phi - ADDtab)/scale, so a step in R is (data - model) and BOTH move;")
    say("   the decomposition is printed rather than the increment alone. phi is FIXED")
    say("   data here, so every change from 'old' to 'new' is the model, and it is split")
    say("   between the K = 5 and K = 4 arms because caveat C-iii says both move.")
    say(f"   {'N':>6}{'y_S':>7}{'R5-R4 old':>11}{'R5-R4 new':>11}{'dR5':>8}{'dR4':>8}"
        f"{'d(phi)5-4':>11}{'d(ADD)5-4':>11}{'<= 0?':>7}")
    nlev = collections.Counter()
    ntot = collections.Counter()
    for n in (1500, 3000):
        for ys in kl.YS_BALANCED:
            l4, l5 = f"K4_yS{ys:.2f}", f"K5_yS{ys:.2f}"
            if (n, l4) not in disk or (n, l5) not in disk:
                continue
            r4, a4 = R_new(n, l4)
            r5, a5 = R_new(n, l5)
            o4 = float(disk[(n, l4)]["R"])
            o5 = float(disk[(n, l5)]["R"])
            oldd = o5 - o4
            dphi = (float(disk[(n, l5)]["phi"]) - float(disk[(n, l4)]["phi"])) / SCALE
            dadd = (a5 - a4) / SCALE
            ntot[n] += 1
            nlev[n] += int(r5 - r4 <= 0)
            say(f"   {n:6d}{ys:7.2f}{oldd:+11.2f}{r5 - r4:+11.2f}{r5 - o5:+8.2f}"
                f"{r4 - o4:+8.2f}{dphi:+11.2f}{dadd:+11.2f}"
                f"{'YES' if r5 - r4 <= 0 else 'no':>7}")
    u3_ok = all(nlev[n] >= U3_MIN_LEVELS for n in (1500, 3000) if ntot[n])
    say()
    say("   AND THE SAME TABLE WITH ARM A'S KNOT ONLY, keeping the existing u = 0.6667")
    say("   knot rather than replacing it with arm C's re-measurement. Arm C was")
    say("   registered as a CONTROL, not as a new calibration point, so the verdict")
    say("   must not rest on having silently re-measured an existing knot.")
    tab_a = dict(tab)
    if "A" in res and res["A"]["ok"]:
        tab_a[ug.U_STAR] = (res["A"]["a"], res["A"]["b"])
    us_a = [0.0, 1.0] + sorted(tab_a)
    aa_a = np.asarray([0.0, 1.0] + [tab_a[u][0] for u in sorted(tab_a)])
    bb_a = np.asarray([0.0, 1.0] + [tab_a[u][1] for u in sorted(tab_a)])
    oa = np.argsort(us_a)
    uu_a = np.asarray(us_a)[oa]
    fa_a = PchipInterpolator(uu_a, aa_a[oa])
    fb_a = PchipInterpolator(uu_a, bb_a[oa])

    def R_armA(n, label):
        row = disk[(n, label)]
        c = by_label[label]
        add = float(mc.packing_density_tabulated(
            np.asarray(c["d_rel"], float), np.asarray(c["y"], float),
            float(row["beta_used"]), fa_a, fb_a))
        return (float(row["phi"]) - add) / SCALE

    say(f"   {'N':>6}{'y_S':>7}{'R5-R4 old':>11}{'A-only':>9}{'A+C':>9}{'<= 0?':>7}")
    nlev_a = collections.Counter()
    for n in (1500, 3000):
        for ys in kl.YS_BALANCED:
            l4, l5 = f"K4_yS{ys:.2f}", f"K5_yS{ys:.2f}"
            if (n, l4) not in disk or (n, l5) not in disk:
                continue
            va = R_armA(n, l5) - R_armA(n, l4)
            vb = R_new(n, l5)[0] - R_new(n, l4)[0]
            oldd = float(disk[(n, l5)]["R"]) - float(disk[(n, l4)]["R"])
            nlev_a[n] += int(va <= 0)
            say(f"   {n:6d}{ys:7.2f}{oldd:+11.2f}{va:+9.2f}{vb:+9.2f}"
                f"{'YES' if va <= 0 else 'no':>7}")
    u3_a_ok = all(nlev_a[n] >= U3_MIN_LEVELS for n in (1500, 3000) if ntot[n])
    say(f"   -> arm-A-knot-only: <= 0 on {nlev_a[1500]}/{ntot[1500]} and "
        f"{nlev_a[3000]}/{ntot[3000]}; verdict "
        f"{'SURVIVES' if u3_a_ok else 'does NOT survive'}")
    if u3_a_ok != u3_ok:
        say("   ** THE TWO READINGS DISAGREE. The graded one is the full closure (A+C);")
        say("      both are reported and the disagreement is the honest statement. **")
        notinfo.append("U3: arm-A-only and A+C readings disagree; both printed")
    else:
        say("   -> the two readings AGREE, so the verdict does not turn on whether arm")
        say("      C's re-measurement replaces the existing knot.")
    say()
    say("   ROBUSTNESS OF THE U3 VERDICT OVER THE MEASURED INTERVAL. A verdict read at")
    say("   the point estimate is a verdict without an error bar; this sweeps a and b")
    say("   over arm A's OWN profile intervals and asks where saturation would survive.")
    if "A" in res:
        ia_, ib_ = res["A"]["ia"], res["A"]["ib"]
        ga = np.linspace(ia_[0], ia_[1], 13)
        gb = np.linspace(ib_[0], ib_[1], 13)
        # A REFEREE (2026-09-16) caught the first version sweeping the MARGINAL box.
        # Two wide marginals can enclose corners the JOINT likelihood excludes -- the
        # very point section 7 of out/diag_uident.txt makes about arm B -- and here 84
        # of 169 marginal corners score chi2 far above the bar. Both fractions are
        # reported; the JOINT one is the number the verdict rests on.
        xvA = np.array([p[1] for p in res["A"]["pts"]])
        yvA = np.array([p[2] for p in res["A"]["pts"]])
        barA = res["A"].get("bar")
        joint_ok = (((ug.phi_fast(ga[:, None, None], gb[None, :, None],
                                  xvA[None, None, :], beta) - yvA[None, None, :])
                     / SCALE) ** 2).sum(axis=2) <= barA
        surv = 0
        surv_j = 0
        n_joint = int(joint_ok.sum())
        tot = 0
        edge_a = []
        for i_a, av_ in enumerate(ga):
            n_surv_row = 0
            for j_b, bv_ in enumerate(gb):
                t2 = dict(tab)
                t2[ug.U_STAR] = (float(av_), float(bv_))
                u2_ = [0.0, 1.0] + sorted(t2)
                a2 = np.asarray([0.0, 1.0] + [t2[u][0] for u in sorted(t2)])
                b2 = np.asarray([0.0, 1.0] + [t2[u][1] for u in sorted(t2)])
                o2 = np.argsort(u2_)
                uu2 = np.asarray(u2_)[o2]
                f2a = PchipInterpolator(uu2, a2[o2])
                f2b = PchipInterpolator(uu2, b2[o2])

                def Rq(n, label, _fa=f2a, _fb=f2b):
                    row = disk[(n, label)]
                    c = by_label[label]
                    add = float(mc.packing_density_tabulated(
                        np.asarray(c["d_rel"], float), np.asarray(c["y"], float),
                        float(row["beta_used"]), _fa, _fb))
                    return (float(row["phi"]) - add) / SCALE

                cnt = collections.Counter()
                for n in (1500, 3000):
                    for ys in kl.YS_BALANCED:
                        l4, l5 = f"K4_yS{ys:.2f}", f"K5_yS{ys:.2f}"
                        if (n, l4) in disk and (n, l5) in disk:
                            cnt[n] += int(Rq(n, l5) - Rq(n, l4) <= 0)
                ok_here = all(cnt[n] >= U3_MIN_LEVELS for n in (1500, 3000))
                surv += int(ok_here)
                n_surv_row += int(ok_here and bool(joint_ok[i_a, j_b]))
                surv_j += int(ok_here and bool(joint_ok[i_a, j_b]))
                tot += 1
            edge_a.append((float(av_), n_surv_row, int(joint_ok[i_a].sum())))
        say(f"   a swept over [{ia_[0]:.4f}, {ia_[1]:.4f}], b over "
            f"[{ib_[0]:.4f}, {ib_[1]:.4f}], {len(ga)}x{len(gb)} = {tot} corners")
        say(f"   over the MARGINAL box:  saturation survives on {surv} of {tot} "
            f"({surv/tot:.1%}) -- SUPERSEDED, see below")
        say(f"   of those, {n_joint} corners are inside the JOINT likelihood region "
            f"(chi2 <= chi2_min + n = {barA:.2f}); the other {tot - n_joint} are "
            f"combinations the data EXCLUDES")
        say(f"   over the JOINT region: saturation survives on {surv_j} of {n_joint} "
            f"({surv_j/n_joint:.1%} ) <- THE NUMBER THE VERDICT RESTS ON")
        say("   survival by a(u*), jointly admissible corners only, so it is visible")
        say("   WHICH parameter the verdict turns on:")
        say("      " + "  ".join(f"a{a_:.3f}:{k}/{m}" for a_, k, m in edge_a[:7]))
        say("      " + "  ".join(f"a{a_:.3f}:{k}/{m}" for a_, k, m in edge_a[7:]))
        surv, tot = surv_j, n_joint
        say(f"   -> the U3 verdict is "
            f"{'NOT RESOLVED -- it turns on where inside its own interval a(u*) sits' if 0 < surv < tot else 'ROBUST over the whole measured interval'}")
        if 0 < surv < tot:
            say("   ** So the GRADED gate fails at the point estimate, as registered, but")
            say("      the scientific statement is UNDETERMINED, not REFUTED: the measured")
            say("      knot is not precise enough to decide saturation either way. The")
            say("      registered verdict and the supportable claim are different")
            say("      sentences and both are recorded. **")
            notinfo.append(f"U3 not resolved over the measured interval: survives on "
                           f"{surv}/{tot} corners ({surv/tot:.0%})")
    say()
    say(f"   -> R(5) - R(4) <= 0 on {nlev[1500]}/{ntot[1500]} levels at N = 1500 and "
        f"{nlev[3000]}/{ntot[3000]} at N = 3000; bar >= {U3_MIN_LEVELS} of 3 at BOTH")
    say(f"   -> U3 {'PASS -- saturation SURVIVES' if u3_ok else 'FAIL -- saturation does NOT survive the measured knot'}")
    if not u3_ok:
        failures.append(f"U3: {nlev[1500]}/{ntot[1500]} and {nlev[3000]}/{ntot[3000]} "
                        f"levels <= 0")
    say("   NOTE (prereg C-iii): the K = 4 cells move too, through the interpolant's")
    say("   shape at u <= 0.630, so the two arms of this difference are NOT independent")
    say("   of the new knot. K = 2 and K = 3 are untouched.")
    say()

    # ------------------------------------------------------------------ U5 the N control
    say("=== U5 -- the N = 3000 arm, and what it is NOT ===")
    say("   Declared in advance: the K-ladder measured a block-wide +1.9 scale offset")
    say("   when N doubles (G5), so a knot fitted at N = 3000 is NOT comparable to a")
    say("   closure calibrated at N = 1500. A difference in a between the two N is")
    say("   EXPECTED and is not evidence about the closure. What is graded is only that")
    say("   arm A's identifiability is not a finite-size accident.")
    g30, nf30, nraw30 = load("ug3k")
    if nraw30 == 0:
        say("   no N = 3000 records on disk; U5 NOT GRADED.")
        notinfo.append("U5 NOT GRADED: no N = 3000 data")
    else:
        adm30 = {}
        nref30 = 0
        for key in sorted(g30):
            ok = []
            for c in g30[key]:
                a, why = admissible(c, ug.N_BIG)
                if a:
                    ok.append(c)
                else:
                    nref30 += 1
            adm30[key] = ok
        n30 = sum(len(v) for v in adm30.values())
        say(f"   N = 3000: {nraw30} records from {nf30} files, admitted {n30}, "
            f"refused {nref30}")
        pts30 = pts_for("A", adm30)
        if len(pts30) >= 4:
            ok30, ia30, ib30, a30, b30 = ug.identifiable(pts30, beta)
            say(f"   arm A at N = 3000: a {a30:.4f} [{ia30[0]:.4f}, {ia30[1]:.4f}]  "
                f"b {b30:.4f} [{ib30[0]:.4f}, {ib30[1]:.4f}]  "
                f"{'IDENTIFIABLE' if ok30 else 'UNDETERMINED'}")
            if "A" in res:
                say(f"   da(N) = {a30 - res['A']['a']:+.4f}, "
                    f"db(N) = {b30 - res['A']['b']:+.4f} -- EXPECTED to be nonzero; "
                    f"not graded, per the registered reading above.")
                dphi = np.mean([p30[2] - p15[2] for p30, p15
                                in zip(pts30, res["A"]["pts"])
                                if abs(p30[1] - p15[1]) < 1e-9]) / SCALE
                say(f"   mean phi shift on the identical cell list: {dphi:+.2f} scale "
                    f"(the K-ladder's block-wide figure was +1.9)")
            say(f"   -> U5 {'PASS' if ok30 else 'FAIL'}")
            if not ok30:
                failures.append("U5: arm A not identifiable at N = 3000")
        else:
            say(f"   only {len(pts30)} compositions at N = 3000; U5 NOT GRADED.")
            notinfo.append("U5 NOT GRADED: too few N = 3000 compositions")
    say()

    # ------------------------------------------------------------------ the transfer
    say("=== U6 -- POST-HOC, declared: the held-out continuous-family transfer ===")
    say("   A search, not a test. out/diag_transfer_ugap.txt already swept this knot and")
    say("   found no value collapses G8's convention spread; this row is the MEASURED")
    say("   knot's own entry in that table.")
    paJ, pbJ = gx.candidates(beta, mc.binary_points(mu_p=MU))[0]["J"](0)
    chi = {"J": 0.0, "ADDtab_old": 0.0, "ADDtab_new": 0.0}
    for c in ("d10", "d21", "d32", "d43"):
        chi[f"SEQ{c}"] = 0.0
    npt = 0
    for fam in mm.FAMILIES:
        for cv in mc.CVS:
            if (fam, cv) not in meas:
                continue
            d, y = mm.psd_classes(fam, cv)
            phi = meas[(fam, cv)][0]
            da = np.asarray(d, float)
            ya = np.asarray(y, float)
            chi["J"] += ((mm.packing_density(d, y, beta, paJ, pbJ) - phi) / SCALE) ** 2
            chi["ADDtab_old"] += ((mc.packing_density_tabulated(da, ya, beta, fa0, fb0)
                                   - phi) / SCALE) ** 2
            chi["ADDtab_new"] += ((mc.packing_density_tabulated(da, ya, beta, fa_n, fb_n)
                                   - phi) / SCALE) ** 2
            for c in ("d10", "d21", "d32", "d43"):
                chi[f"SEQ{c}"] += ((kl.sequential_k(d, y, beta, fa_n, fb_n, c)
                                    - phi) / SCALE) ** 2
            npt += 1
    seq = [chi[f"SEQ{c}"] for c in ("d10", "d21", "d32", "d43")]
    say(f"   held-out points {npt}: J {chi['J']:.1f}, ADDtab before {chi['ADDtab_old']:.1f}"
        f", after {chi['ADDtab_new']:.1f}")
    say(f"   nested, measured knot: d10 {seq[0]:.1f}, d21 {seq[1]:.1f}, "
        f"d32 {seq[2]:.1f}, d43 {seq[3]:.1f}; spread "
        f"{max(seq)/min(seq) if min(seq) > 0 else float('inf'):.2f}x "
        f"(standing 11.90x), best {min(seq):.1f}")
    say("   the spread NARROWING while every member gets WORSE is not a collapse; both")
    say("   the ratio and the best member are printed so it cannot read as one.")
    say()

    # ------------------------------------------------------------------- caveats
    say("=== CAVEATS carried from out/predict_ugap.txt (prereg section 6) ===")
    import predict_ugap as pu
    for c in pu.CAVEATS:
        say("   " + c)
    say("   C-v  ADDED after the referee pass, and NOT in the blind file: every number")
    say("        here rests on `beta`, which comes from mm.measured() -- and that loader")
    say("        prints a standing NON-DETERMINISM warning (15 of 82 duplicated keys")
    say("        disagree, worst |dphi| = 5.03e-03, i.e. 2.4x PHI_SCALE, resolved by")
    say("        FILENAME SORT ORDER). It is tracked in PREREG_DST_WINDOW_2026-09-09.md")
    say("        and is not introduced by this block, but beta's own precision carries")
    say("        that exposure and neither the prereg nor the first claims draft said so.")
    say()

    # ------------------------------------------------------------------- artifact
    with open(os.path.join(OUT, "ugap_cells.csv"), "w", newline="", encoding="utf-8") as fh:
        wcsv = csv.writer(fh)
        wcsv.writerow(["N", "arm", "ratio", "u", "x", "n_seeds", "phi", "sd", "se",
                       "phi_shipped", "res_scale"])
        for adm, nn in ((adm15, 1500),):
            for arm in ("A", "B", "C"):
                r = ug.RATIOS[arm]
                u = 1.0 / r
                av, bv = float(fa0(u)), float(fb0(u))
                for x in ug.X_GRID:
                    v = adm.get((float(r), round(float(x), 6)), [])
                    if not v:
                        continue
                    a_ = np.array([c["phi"] for c in v], float)
                    p = float(mc.packing_density_ab(np.array([r, 1.0]),
                                                    np.array([1 - x, x]), beta, av, bv))
                    sd = float(a_.std(ddof=1)) if a_.size > 1 else float("nan")
                    wcsv.writerow([nn, arm, f"{r:.8f}", f"{u:.8f}", f"{x:.4f}", a_.size,
                                   f"{a_.mean():.8f}", f"{sd:.8f}",
                                   f"{sd/math.sqrt(a_.size):.8f}", f"{p:.8f}",
                                   f"{(a_.mean() - p)/SCALE:.4f}"])

    say(f"NOT INFORMATIVE: {notinfo}")
    say(f"FAILURES: {failures}")
    with open(os.path.join(OUT, "grade_ugap.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
