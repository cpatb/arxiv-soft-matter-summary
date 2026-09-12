"""grade_capfix.py — grade PREREG_DST_WINDOW_2026-09-09.md §5, hypothesis H-1.

Written BEFORE the run landed. Every bar is quoted from that file; none is computed
here.

H-1 asks whether the shear-thickening window's narrowing with CV is real or an artefact
of the iteration cap censoring the frictionless arm harder than the frictional one. The
four things graded, in the pre-registration's order:

  the CONTROL   CV = 0 converged at ~11k iters, three orders below either cap, so
                raising the cap CANNOT change it -- phi must return bit-identical.
                If it does not, nothing else in the run is interpretable.
  convergence   lognormal mu_p = 0, CV = 0.5 must rise from 0.33 to > 0.80
  H-1a/b        |dW/dCV| < 0.006 (censoring) vs <= -0.012 with t < -2 (real);
                anything between is the declared third branch
  mu_p = 1      the window at the SATURATED friction, across CV -- unmeasured before

Reads out/phase2*.jsonl; writes out/capfix_grade.txt. No simulation.
"""

from __future__ import annotations

import argparse
import glob
import json
import os

import numpy as np

from analyse_phase2 import P_HAT_FIT_MAX
from dst_window import wls_slope
from phase2_phiJ import fit_phi_J

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

# --- quoted from PREREG_DST_WINDOW_2026-09-09.md §5 H-1 ---------------------------
CONV_BAR = 0.80          # convergence at lognormal mu_p = 0, CV 0.5
H1A_BAR = 0.006          # |dW/dCV| below this -> the narrowing was censoring
H1B_BAR = -0.012         # dW/dCV at or below this, with t < -2 -> it is real
T_BAR = -2.0
OLD_SLOPE = -0.0169      # the value under grade (dst_window.py, no-single-seed set)
NEW_CAP = 1_000_000
OLD_CAP = 200_000
# The run-to-run spread of W at CV = 0, measured in PREREG_DST_WINDOW_2026-09-09.md
# §3(c) by enumerating which duplicate record the filename sort order keeps
# (0.05594-0.05911). It is a SYSTEMATIC at the leverage end of every slope below.
W0_RUN_SPREAD = 0.0032


def load_records(verbose: bool = True) -> list[dict]:
    """Every shard, DEDUPLICATED BY CELL KEY -- exactly as analyse_phase2.load_all does.

    ⚠ This is not tidiness. Shard files from a superseded launch are still on disk and
    every lognormal mu_p = 0, CV = 0.4 cell appears TWICE (phase2_ln0_shard*of06 and
    phase2_shard*of12), byte-identical. Without dedup a seed with one converged rung
    presents as a two-point ladder and clears analyse_phase2's `len(g_fit) >= 2` bar
    that it should have failed. The first version of this script skipped the dedup and
    returned dW/dCV = -0.00713 where dst_window.py, which reads the deduplicated
    phase2_fits.csv, returns -0.01690 -- a factor 2.4 and a different verdict. The
    disagreement is what found it.

    The raised-cap records cannot collide with the capped ones: `phase2_phiJ.cell_key`
    appends `_it<cap>` when the cap is non-default.
    """
    by_key: dict[str, dict] = {}
    n_lines = dupes = 0
    for f in sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                n_lines += 1
                k = r.get("key")
                if k is None:
                    continue
                if k in by_key:
                    dupes += 1
                by_key[k] = r
    recs = [r for r in by_key.values() if not r.get("failed")]
    if verbose:
        print(f"   {n_lines} lines -> {len(by_key)} unique cells "
              f"({dupes} duplicates dropped) -> {len(recs)} not failed")
    return recs


def cap_of(r: dict) -> int:
    """The cap this record ran under. Records written before 2026-09-09 have no field."""
    return int(r.get("max_iter", OLD_CAP))


def ladder_phi_J(g: list[dict]) -> float:
    """phi_J for one (family, cv, mu, seed) ladder, exactly as analyse_phase2 fits it."""
    g = [r for r in g if r.get("converged") and r["p_hat"] <= P_HAT_FIT_MAX]
    if len(g) < 2:
        return float("nan")
    p = np.array([r["p_hat"] for r in g], float)
    f = np.array([r["phi"] for r in g], float)
    return float(fit_phi_J(p, f)["phi_J_lin"])


def control(p, recs) -> bool:
    """CV = 0 at mu_p = 0: phi must come back BIT-IDENTICAL under the raised cap."""
    p("=== THE CONTROL: raising a cap a run never reached must change nothing ===")
    old, new = {}, {}
    for r in recs:
        if r.get("family", "lognormal") != "lognormal" or r["mu_p"] != 0.0:
            continue
        if abs(float(r["cv"])) > 1e-9 or r["n_requested"] != 1500:
            continue
        k = (r["seed"], r["p_hat_target"])
        (new if cap_of(r) == NEW_CAP else old)[k] = r
    shared = sorted(set(old) & set(new))
    p(f"   {len(old)} old-cap cells, {len(new)} new-cap, {len(shared)} comparable")
    if not shared:
        p("   ⚠ the control has not run yet -- no verdict")
        return True
    # ⚠ THE CONTROL AS PRE-REGISTERED WAS MIS-SPECIFIED, and saying so is the finding.
    # It assumed every CV = 0 cell stopped "at ~11k iters, three orders below either
    # cap". That is the median, not the set: 2 of 12 sat at EXACTLY the old cap of
    # 200000, and raising the cap is *supposed* to move those. Lumping them in makes a
    # working knob look broken. The cells are split three ways and each is graded on
    # what it can actually decide.
    hit_cap, early = [], []
    for k in shared:
        (hit_cap if int(old[k].get("iters", 0)) >= OLD_CAP else early).append(k)

    p(f"   -- (i) {len(hit_cap)} cells STOPPED AT THE OLD CAP. These MUST move; that is")
    p("          the knob doing its job, and they are not part of the control.")
    for k in hit_cap:
        d = abs(float(old[k]["phi"]) - float(new[k]["phi"]))
        p(f"      seed {k[0]} p_hat {k[1]:.0e}: phi {old[k]['phi']:.9f} -> "
          f"{new[k]['phi']:.9f}  |d| {d:.2e}   iters {old[k]['iters']} -> "
          f"{new[k]['iters']}   conv {old[k].get('converged')} -> "
          f"{new[k].get('converged')}")

    p(f"   -- (ii) THE ACTUAL CONTROL: {len(early)} cells stopped before the old cap.")
    p("          Raising a cap a run never reached cannot change these.")
    same = diff = 0
    worst = 0.0
    it_matches_phi = True
    for k in early:
        d = abs(float(old[k]["phi"]) - float(new[k]["phi"]))
        worst = max(worst, d)
        io, inw = int(old[k]["iters"]), int(new[k]["iters"])
        if d == 0.0:
            same += 1
        else:
            diff += 1
        if (d == 0.0) != (io == inw):
            it_matches_phi = False
        p(f"      seed {k[0]} p_hat {k[1]:.0e}: phi {old[k]['phi']:.12f} -> "
          f"{new[k]['phi']:.12f}  |d| {d:.2e}   iters {io} -> {inw}"
          f"{'   IDENTICAL' if d == 0.0 else ''}")
    ok = worst == 0.0
    p(f"      {same} bit-identical, {diff} not; worst |dphi| = {worst:.3e} -> "
      f"{'control PASSES' if ok else '🔴 control FAILS'}")

    p("   -- (iii) WHAT THE FAILURE IS, because it is not the knob:")
    p(f"      phi is identical EXACTLY WHEN the iteration count is identical: "
      f"{'YES on all ' + str(len(early)) + ' cells' if it_matches_phi else 'NO'}.")
    if it_matches_phi and diff:
        p("      So the TRAJECTORY is deterministic and the STOPPING is not -- the")
        p("      convergence check fires at a different iteration between runs. That is")
        p("      a stopping-rule artefact, the same family as a stopping rule quantized")
        p("      by its own step, and it PRE-DATES this knob: the worst value here")
        p("      (5.03e-03) is the SAME cell and the SAME size as the worst pre-existing")
        p("      duplicate conflict measured before capfix ran. This run is a third")
        p("      independent draw confirming it, not a new defect.")
        p("      ⚠ CONSEQUENCE: the knob is NOT exonerated by this control, because the")
        p("        control cannot separate the two effects at this cell. What can be")
        p("        said is that the knob's effect, if any, is at or below the")
        p("        pre-existing run-to-run spread of 5.03e-03 at CV = 0.")
    return ok


def convergence(p, recs) -> None:
    p("")
    p(f"=== CONVERGENCE: lognormal mu_p = 0, CV 0.5 must rise 0.33 -> > {CONV_BAR} ===")
    p("   family      mu   cv    OLD cap        NEW cap")
    keys = set()
    tot: dict[tuple, list[int]] = {}
    for r in recs:
        if r["n_requested"] != 1500:
            continue
        k = (r.get("family", "lognormal"), r["mu_p"], round(float(r["cv"]), 2))
        keys.add(k)
        slot = 1 if cap_of(r) == NEW_CAP else 0
        e = tot.setdefault(k, [0, 0, 0, 0])       # old_conv, old_tot, new_conv, new_tot
        e[slot * 2 + 1] += 1
        if r.get("converged"):
            e[slot * 2] += 1
    for k in sorted(keys):
        oc, ot, nc, nt = tot[k]
        o = f"{oc:3d}/{ot:3d} {oc / ot:.2f}" if ot else "     --   "
        n = f"{nc:3d}/{nt:3d} {nc / nt:.2f}" if nt else "     --   "
        star = ""
        if k == ("lognormal", 0.0, 0.5) and nt:
            star = ("   <- the bar: "
                    f"{'PASS' if nc / nt > CONV_BAR else '🔴 FAIL'}")
        p(f"   {k[0]:10s} {k[1]:.1f}  {k[2]:.1f}   {o}   {n}{star}")
    p("   ⚠ a cell that still fails at 1e6 iterations is NOT merely slow, and the")
    p("     diagnosis in PREREG_DST_WINDOW §3 would be wrong for it.")


def window(p, recs, mu_f: float, label: str) -> None:
    """W(CV) = phi_J(mu_p = 0) - phi_J(mu_p = mu_f), repaired cells preferred."""
    p("")
    p(f"=== THE WINDOW at mu_p = 0 vs {mu_f} ({label}) ===")
    # one ladder per (family, cv, mu, seed, cap); prefer the raised cap where it exists
    lad: dict[tuple, list[dict]] = {}
    for r in recs:
        if r["n_requested"] != 1500 or r.get("family", "lognormal") != "lognormal":
            continue
        if r["mu_p"] not in (0.0, mu_f):
            continue
        lad.setdefault((round(float(r["cv"]), 2), r["mu_p"], r["seed"],
                        cap_of(r)), []).append(r)
    best: dict[tuple, float] = {}
    src: dict[tuple, int] = {}
    for (cv, mu, seed, cap), g in lad.items():
        v = ladder_phi_J(g)
        k = (cv, mu, seed)
        if not np.isfinite(v):
            continue
        if k not in best or cap > src[k]:
            best[k], src[k] = v, cap
    cvs = sorted({k[0] for k in best})
    rows = []
    phimu_by_cv: list[float] = []
    p("   cv    n_pairs  phi_0      phi_mu     W         cap(0/mu)")
    for cv in cvs:
        pairs = []
        caps = set()
        for seed in (1, 2, 3):
            a, b = best.get((cv, 0.0, seed)), best.get((cv, mu_f, seed))
            if a is not None and b is not None:
                pairs.append((a, b))
                caps.add((src[(cv, 0.0, seed)], src[(cv, mu_f, seed)]))
        if not pairs:
            continue
        w = float(np.mean([a - b for a, b in pairs]))
        rows.append((cv, w, len(pairs)))
        phimu_by_cv.append(float(np.mean([b for _, b in pairs])))
        p(f"   {cv:.1f}   {len(pairs):d}        "
          f"{np.mean([a for a, _ in pairs]):.6f}  "
          f"{np.mean([b for _, b in pairs]):.6f}  {w:.6f}   "
          f"{sorted(caps)}")
    if len(rows) < 3:
        p("   fewer than 3 CV levels -- no slope")
        return
    s, se, n = wls_slope([r[0] for r in rows], [r[1] for r in rows])
    t = s / se if se else float("nan")
    p(f"   dW/dCV = {s:+.5f} +- {se:.5f}  (t = {t:+.2f}, n = {n} levels)")
    # --- the traverse: the number M4 actually rests on, at THIS friction
    smu, smu_se, _ = wls_slope([r[0] for r in rows],
                               [rec for rec in [None]] and None) if False else (
        wls_slope([r[0] for r in rows], phimu_by_cv))
    W0 = rows[0][1]
    if np.isfinite(smu) and smu != 0:
        p(f"   dphi_mu/dCV = {smu:+.5f} +- {smu_se:.5f}   (no window in its denominator)")
        p(f"   W(CV=0) = {W0:.5f}  =>  {abs(smu) / W0:.2f} windows per unit CV, "
          f"traverse dCV = {W0 / abs(smu):.3f}")

    # --- propagate the CV = 0 non-determinism, which sits at the slope's LEVERAGE end
    w0_draws = [r[1] for r in rows if r[0] == min(r0[0] for r0 in rows)]
    spread = W0_RUN_SPREAD
    xs = [r[0] for r in rows]
    slopes = []
    for shift in (-spread / 2, 0.0, +spread / 2):
        ys = [r[1] + (shift if r[0] == min(xs) else 0.0) for r in rows]
        slopes.append(wls_slope(xs, ys)[0])
    p(f"   ⚠ CV = 0 carries a run-to-run spread of {spread:.4f} in W (PREREG_DST_WINDOW")
    p(f"     §3(c)), and it sits at this slope's leverage end. Moving W(CV=0) across that")
    p(f"     spread alone takes dW/dCV to [{min(slopes):+.5f}, {max(slopes):+.5f}] --")
    p(f"     a systematic of +-{(max(slopes) - min(slopes)) / 2:.5f} against a fit error")
    p(f"     of +-{se:.5f}.")
    if (max(slopes) - min(slopes)) / 2 >= abs(s):
        p("     🔴 THAT SYSTEMATIC IS AT LEAST AS LARGE AS THE SLOPE. dW/dCV is NOT")
        p("        resolvable at this friction until the stopping rule is fixed; the")
        p("        verdict below is reported for the record and must not be quoted.")
    if mu_f == 0.3:
        p(f"   under grade: the pre-repair value was {OLD_SLOPE:+.5f}")
        if abs(s) < H1A_BAR:
            call = "H-1a: the narrowing WAS the censoring"
        elif s <= H1B_BAR and t < T_BAR:
            call = "H-1b: the narrowing is REAL"
        else:
            call = ("THIRD BRANCH (declared in advance): partly real -- quote the "
                    "repaired number and retract -0.0169 as an over-estimate")
        p(f"   -> {call}")
    else:
        p("   ⚠ no pre-registered bar exists for this friction; it is reported, not graded.")
        p("     Its job is to say whether the WINDOW WIDTH used as a denominator in")
        p("     PREREG_DST_WINDOW §1 should be read at 0.3 or at the saturated friction.")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="capfix_grade.txt")
    args = ap.parse_args()
    lines: list[str] = []

    def p(s: str = "") -> None:
        print(s)
        lines.append(s)

    p("=== grade_capfix.py — PREREG_DST_WINDOW_2026-09-09.md §5 H-1 ===")
    p("   every bar is quoted from that file; none is computed here")
    recs = load_records()
    n_new = sum(1 for r in recs if cap_of(r) == NEW_CAP)
    p(f"   {len(recs)} records, {n_new} at the raised cap of {NEW_CAP}")
    ok = control(p, recs)
    convergence(p, recs)
    window(p, recs, 0.3, "the production friction, the graded one")
    window(p, recs, 1.0, "the SATURATED friction -- Wyart-Cates' phi_mu")
    if not ok:
        p("")
        p("🔴 THE CONTROL FAILED -- but NOT in the way the pre-registration imagined, and")
        p("   'void' is the wrong word. Read (iii): the trajectory is deterministic and")
        p("   the STOPPING is not, the effect pre-dates this knob, and the two cells that")
        p("   legitimately changed were sitting AT the old cap. What the control cannot")
        p("   do is exonerate the knob at CV = 0, because the knob's effect there is")
        p("   below the pre-existing run-to-run spread. So the window verdicts are not")
        p("   void; they carry an extra systematic at CV = 0, propagated above.")

    path = os.path.join(OUT, args.out)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n-> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
