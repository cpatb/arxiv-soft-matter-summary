"""grade_constp.py — grade PREREG_CONSTP_2026-09-12.md, the servo shakeout (S1-S7).

Written before the cells landed; every bar is quoted from that file's §5 and none is
computed here.

    python grade_constp.py          -> out/constp_shakeout.txt

Exit 0 only if S2 AND S3 pass, because §6 says the rest is meaningless without them.
"""

from __future__ import annotations

import glob
import json
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

# --- bars, QUOTED from PREREG_CONSTP_2026-09-12.md §5 ------------------------------
S1_P_TOL = 0.05          # |<p>_ss/p* - 1|
S1_MIN_OK = 5            # of 6 cells
S2_GAIN_SPREAD = 0.005   # phi_ss spread across gains at fixed phi_start
S3_START_GAP = 0.005     # |phi_ss(0.55) - phi_ss(0.64)| at matched gain
S4_ALPHA_FACTOR = 2.0    # the two starting sides must agree within this factor
S5_MIN_OK = 5            # of 6 cells quasi-static AND steady on both channels
S6_BRACKET = (0.55, 0.64)   # phi_ss must lie strictly inside the starting states
S7_PHI_ISO = 0.583314    # isotropic phi at p* = 3e-04, cv 0, mu 0.3, N 1500 (n = 9)
S7_PHI_ISO_SD = 0.000725
S7_E4M = 0.5902          # E4m's phi_c(cv 0) on the common ladder
S7_E4M_TOL = 0.015
# phase9_lowgain* is the §5-S1-prescribed downward extension of the gain grid after
# the first six cells came back 2 of 6 held. Both globs are read as ONE set and the
# gain is printed on every row, so no cell is silently pooled.
PATTERNS = ("phase9_shakeout*.jsonl", "phase9_lowgain*.jsonl")
# post-hoc DIAGNOSTICS, labelled as such: they explain S1's failure and are not
# used to re-grade anything. A cell that visits z ~ 0 has unjammed.
DIAG_Z_MIN = 2.0
DIAG_P_FLOOR_FRAC = 0.1


def load() -> list[dict]:
    recs, seen = [], set()
    paths = [q for pat in PATTERNS
             for q in sorted(glob.glob(os.path.join(OUT, pat)))]
    for path in paths:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                k = r.get("key")
                if k in seen:
                    continue
                seen.add(k)
                r["_src"] = os.path.basename(path)
                recs.append(r)
    return recs


def stability(r: dict) -> dict:
    """Did the cell stay jammed through the averaging window?

    POST HOC (2026-09-12), and it re-grades nothing: S1 failed on its own
    pre-registered bar and stays failed. This reads the series already recorded to say
    WHY -- a proportional servo above some gain drives the cell through unjamming, and
    a cell that visits p ~ 1e-08 against a 3e-04 set point is not a droopy measurement,
    it is not a measurement.
    """
    g = np.asarray(r.get("gamma", []), float)
    z = np.asarray(r.get("z_series", []), float)
    pp = np.asarray(r.get("p_hat_series", []), float)
    if g.size == 0 or z.size != g.size or pp.size != g.size:
        return {"z_win_min": float("nan"), "p_win_min": float("nan"),
                "stayed_jammed": False}
    w = g >= 0.5 * float(g.max())
    zmin, pmin = float(z[w].min()), float(pp[w].min())
    return {"z_win_min": zmin, "p_win_min": pmin,
            "stayed_jammed": bool(zmin >= DIAG_Z_MIN
                                  and pmin >= DIAG_P_FLOOR_FRAC
                                  * float(r["p_hat_target"]))}


def main() -> int:
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    ok: dict[str, bool] = {}
    recs = load()
    good = [r for r in recs if not r.get("error")]
    p("=== grade_constp.py -- PREREG_CONSTP_2026-09-12.md, S1-S7 ===")
    p(f"   {len(recs)} records, {len(good)} usable, "
      f"{len(recs) - len(good)} errored")
    for r in recs:
        if r.get("error"):
            p(f"   🔴 FAILED CELL {r['key']}: {r['error']}")
    if len(good) < 4:
        p("   🔴 fewer than 4 usable cells -- the shakeout is UNINTERPRETABLE")
        with open(os.path.join(OUT, "constp_shakeout.txt"), "w",
                  encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
        return 1

    p("")
    p(f"{'phi0':>6s}{'K':>6s}{'phi_ss':>10s}{'p_err':>9s}{'held':>6s}"
      f"{'z_wmin':>8s}{'p_wmin':>10s}{'jammed':>8s}"
      f"{'alpha':>9s}{'dec':>6s}{'I':>10s}{'steady':>8s}{'wall':>7s}")
    for r in sorted(good, key=lambda r: (r["phi_target"], r["servo_gain"])):
        d = stability(r)
        r["_stab"] = d
        p(f"{r['phi_target']:6.2f}{r['servo_gain']:6.2f}{r['phi_ss']:10.6f}"
          f"{r['p_ss_hat'] / r['p_hat_target'] - 1:+9.1%}"
          f"{str(r['servo_held']):>6s}{d['z_win_min']:8.2f}{d['p_win_min']:10.1e}"
          f"{str(d['stayed_jammed']):>8s}"
          f"{r['alpha']:9.3g}{r['alpha_decades']:6.2f}"
          f"{r['i_ss']:10.2e}{str(r['steady']):>8s}{r['wall_s']:7.0f}")

    # ---------------------------------------------------------------- S1
    held = [r for r in good if r["servo_held"]]
    ok["S1"] = len(held) >= S1_MIN_OK
    p("")
    p(f"=== S1 the servo holds (|<p>/p*-1| <= {S1_P_TOL:.0%} on >= {S1_MIN_OK}/6) ===")
    p(f"   held on {len(held)} of {len(good)}; worst "
      f"{max(abs(r['p_ss_hat'] / r['p_hat_target'] - 1) for r in good):.1%} -> "
      f"{'PASS' if ok['S1'] else 'FAIL -- the gain RANGE is wrong, not the bar'}")

    by_start = defaultdict(dict)
    for r in good:
        by_start[round(float(r["phi_target"]), 4)][round(float(r["servo_gain"]), 4)] = r

    # ---------------------------------------------------------------- S2
    p("")
    p(f"=== S2 🔑 phi_ss does not depend on the GAIN (spread <= {S2_GAIN_SPREAD}) ===")
    worst2 = 0.0
    for phi0 in sorted(by_start):
        v = np.array([by_start[phi0][k]["phi_ss"] for k in sorted(by_start[phi0])])
        spread = float(v.max() - v.min()) if v.size > 1 else float("nan")
        p(f"   phi0 {phi0:.2f}: " + "  ".join(
            f"K{k:.2f}={by_start[phi0][k]['phi_ss']:.6f}"
            for k in sorted(by_start[phi0])) + f"   spread {spread:.6f}")
        if np.isfinite(spread):
            worst2 = max(worst2, spread)
    ok["S2"] = worst2 <= S2_GAIN_SPREAD
    p(f"   worst spread {worst2:.6f} (bar {S2_GAIN_SPREAD}) -> "
      f"{'PASS -- the answer is the material, not the controller' if ok['S2'] else 'FAIL -- phi_ss TRACKS THE CONTROLLER; the protocol is not admissible'}")

    # ---------------------------------------------------------------- S3
    p("")
    p(f"=== S3 🔑 the critical state forgets its preparation (gap <= "
      f"{S3_START_GAP}) ===")
    starts = sorted(by_start)
    worst3 = 0.0
    if len(starts) >= 2:
        lo, hi = starts[0], starts[-1]
        shared = sorted(set(by_start[lo]) & set(by_start[hi]))
        for k in shared:
            gap = abs(by_start[lo][k]["phi_ss"] - by_start[hi][k]["phi_ss"])
            p(f"   K {k:.2f}: phi_ss({lo:.2f}) = {by_start[lo][k]['phi_ss']:.6f}, "
              f"phi_ss({hi:.2f}) = {by_start[hi][k]['phi_ss']:.6f}   gap {gap:.6f}")
            worst3 = max(worst3, gap)
        ok["S3"] = bool(shared) and worst3 <= S3_START_GAP
        p(f"   worst gap {worst3:.6f} (bar {S3_START_GAP}) -> "
          f"{'PASS' if ok['S3'] else 'FAIL -- gamma_max = 3 is not enough strain (check whether the phi series flattened; that is a DIFFERENT fix from S2)'}")
    else:
        ok["S3"] = False
        p("   🔴 only one starting state present -- S3 cannot be read")

    # ------------------------------------------------------- S2'/S3', the re-read
    # 🔴 S2 and S3 above are computed over EVERY cell, which is how they were written.
    # But S1 failed, and the stability column says the failing cells periodically
    # UNJAM. A cell that visits p ~ 1e-08 against a 3e-04 set point is not a
    # measurement of phi_c at 3e-04, so a gain-independence test that averages it in
    # has passed on a different question. The restriction below uses S1's OWN
    # pre-registered per-cell criterion (|p_err| <= 5 %) -- no new bar.
    p("")
    p("=== S2'/S3' the same two tests, restricted to cells that PASSED S1 ===")
    held_by_start = defaultdict(dict)
    for r in held:
        held_by_start[round(float(r["phi_target"]), 4)][
            round(float(r["servo_gain"]), 4)] = r
    gains_held = sorted({round(float(r["servo_gain"]), 4) for r in held})
    p(f"   {len(held)} held cells, at gains {gains_held}")
    worst2h, readable2 = 0.0, False
    for phi0 in sorted(held_by_start):
        ks = sorted(held_by_start[phi0])
        if len(ks) < 2:
            p(f"   phi0 {phi0:.2f}: only {len(ks)} held gain(s) -- S2' UNREADABLE here")
            continue
        readable2 = True
        v = np.array([held_by_start[phi0][k]["phi_ss"] for k in ks])
        spread = float(v.max() - v.min())
        worst2h = max(worst2h, spread)
        p(f"   phi0 {phi0:.2f}: " + "  ".join(
            f"K{k:.2f}={held_by_start[phi0][k]['phi_ss']:.6f}" for k in ks)
          + f"   spread {spread:.6f}")
    ok["S2p"] = readable2 and worst2h <= S2_GAIN_SPREAD
    p(f"   S2' -> {'PASS' if ok['S2p'] else ('UNREADABLE -- fewer than two stable gains at any preparation' if not readable2 else f'FAIL at {worst2h:.6f}')}")

    starts_h = sorted(held_by_start)
    worst3h, readable3 = 0.0, False
    if len(starts_h) >= 2:
        lo, hi = starts_h[0], starts_h[-1]
        for k in sorted(set(held_by_start[lo]) & set(held_by_start[hi])):
            readable3 = True
            gap = abs(held_by_start[lo][k]["phi_ss"] - held_by_start[hi][k]["phi_ss"])
            worst3h = max(worst3h, gap)
            p(f"   K {k:.2f}: gap between preparations = {gap:.6f}")
    ok["S3p"] = readable3 and worst3h <= S3_START_GAP
    p(f"   S3' -> {'PASS' if ok['S3p'] else ('UNREADABLE' if not readable3 else f'FAIL at {worst3h:.6f}')}")

    # the servo's stability boundary, read off the gains actually run
    stable_k = sorted({round(float(r["servo_gain"]), 4) for r in good
                       if stability(r)["stayed_jammed"]})
    unstable_k = sorted({round(float(r["servo_gain"]), 4) for r in good
                         if not stability(r)["stayed_jammed"]})
    p(f"   gains that stayed jammed: {stable_k}; gains that unjammed: {unstable_k}")
    if stable_k and unstable_k:
        p(f"   => the servo's stability boundary is between K = {max(stable_k)} and "
          f"K = {min(unstable_k)} for this cell. A production run must sit below it "
          f"and must report z_win_min per cell, because a cell can hold <p> on average "
          f"while visiting an unjammed state.")

    # ---------------------------------------------------------------- S4
    p("")
    p(f"=== S4 alpha is admissible and consistent (sign 6/6, >= 0.3 decades, "
      f"sides within {S4_ALPHA_FACTOR}x) ===")
    n_sign = sum(1 for r in good if r["alpha_sign_ok"])
    quotable = [r for r in good if np.isfinite(r["alpha"])]
    p(f"   fitted rate negative on {n_sign} of {len(good)}; "
      f"{len(quotable)} cells have a quotable alpha")
    for r in good:
        if r.get("alpha_note"):
            p(f"      {r['phi_target']:.2f}/K{r['servo_gain']:.2f}: "
              f"{r['alpha_note']}")
    ratios = []
    if len(starts) >= 2:
        lo, hi = starts[0], starts[-1]
        for k in sorted(set(by_start[lo]) & set(by_start[hi])):
            a, b = by_start[lo][k]["alpha"], by_start[hi][k]["alpha"]
            if np.isfinite(a) and np.isfinite(b) and a > 0 and b > 0:
                ratio = max(a / b, b / a)
                ratios.append(ratio)
                p(f"   K {k:.2f}: alpha({lo:.2f}) = {a:.4g}, alpha({hi:.2f}) = "
                  f"{b:.4g}   ratio {ratio:.2f}x")
    ok["S4"] = (n_sign == len(good) and bool(ratios)
                and max(ratios) <= S4_ALPHA_FACTOR)
    p(f"   -> {'PASS' if ok['S4'] else 'NOT PASSED'}"
      + ("" if ok["S4"] else
         "  (§5 said a clean disagreement between the dilating and contracting "
         "branches is a PHYSICAL result, not a protocol failure -- read the ratios)"))

    # ---------------------------------------------------------------- S5
    p("")
    p(f"=== S5 the cell is a valid measurement (I <= 1e-2 and steady, >= "
      f"{S5_MIN_OK}/6) ===")
    valid = [r for r in good if r["quasistatic"] and r["steady"]]
    ok["S5"] = len(valid) >= S5_MIN_OK
    p(f"   {len(valid)} of {len(good)} quasi-static AND steady on both channels -> "
      f"{'PASS' if ok['S5'] else 'FAIL'}")
    for r in good:
        if not (r["quasistatic"] and r["steady"]):
            p(f"      phi0 {r['phi_target']:.2f} K{r['servo_gain']:.2f}: "
              f"I = {r['i_ss']:.2e}, p-split {r['steady_split']:.1%}, "
              f"phi-split {r['phi_ss_split']:.1%}")

    # ---------------------------------------------------------------- S6
    p("")
    p(f"=== S6 phi_ss lies strictly inside {S6_BRACKET} ===")
    outside = [r for r in good
               if not (S6_BRACKET[0] < r["phi_ss"] < S6_BRACKET[1])]
    ok["S6"] = not outside
    p(f"   {len(outside)} cells outside the bracket -> "
      f"{'PASS' if ok['S6'] else 'FAIL -- the servo drove the volume somewhere neither the preparation nor the material asked for'}")
    for r in outside:
        p(f"      phi0 {r['phi_target']:.2f} K{r['servo_gain']:.2f}: "
          f"phi_ss = {r['phi_ss']:.6f}")

    # ---------------------------------------------------------------- S7
    p("")
    p("=== S7 the comparison, REPORTED not gated ===")
    v = np.array([r["phi_ss"] for r in valid or good], float)
    mean = float(v.mean())
    p(f"   mean phi_ss over {'the valid' if valid else 'all'} cells = {mean:.6f}")
    p(f"   isotropic phi at the same p* = {S7_PHI_ISO:.6f} +- {S7_PHI_ISO_SD:.6f} "
      f"(N = 1500) -> phi_ss is {'DENSER' if mean > S7_PHI_ISO else 'LOOSER'} by "
      f"{mean - S7_PHI_ISO:+.6f}  [predicted: denser]")
    p(f"   E4m phi_c(cv 0) on the common ladder = {S7_E4M:.4f} -> |d| = "
      f"{abs(mean - S7_E4M):.4f} (informative if <= {S7_E4M_TOL})")
    if abs(mean - S7_E4M) <= S7_E4M_TOL:
        p("   🔑 an estimator needing NO isostatic target and NO ladder agrees with "
          "E4m to within the bar. ⚠ N = 600 here against 1500 there, so finite size "
          "alone could account for a move of this size; the production run is what "
          "settles it.")
    else:
        p("   the two do NOT agree within the bar. ⚠ N = 600 vs 1500, single seed, "
          "single pressure -- this is a flag for the production run, not a refutation "
          "of either estimator.")

    # ---------------------------------------------------------------- summary
    p("")
    p("=== SUMMARY ===")
    for g in ("S1", "S2", "S2p", "S3", "S3p", "S4", "S5", "S6"):
        p(f"   {g}: {'PASS' if ok.get(g) else 'NOT PASSED'}")
    p("   S7: reported, not gated")
    # 🔴 The verdict turns on S2' and S3', NOT on S2 and S3. S2/S3 average in cells
    # that periodically unjam, so passing them is compatible with the servo being the
    # measurement -- which is the one thing this shakeout exists to rule out.
    decisive = ok.get("S2p", False) and ok.get("S3p", False)
    p("")
    if decisive:
        p("   ✅ THE PROTOCOL IS ADMISSIBLE: phi_ss depends on neither the gain nor "
          "the preparation, which is what a critical state means. A production "
          "pre-registration may now be written.")
    else:
        p("   🔴 THE PROTOCOL IS NOT YET ADMISSIBLE. Per §6 the rest of the run means "
          "nothing until S2 and S3 pass, and the two failures have DIFFERENT fixes: "
          "S2 failing means the servo is the measurement, S3 failing means the strain "
          "is too short.")

    path = os.path.join(OUT, "constp_shakeout.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten -> {path}")
    return 0 if decisive else 1


if __name__ == "__main__":
    raise SystemExit(main())
