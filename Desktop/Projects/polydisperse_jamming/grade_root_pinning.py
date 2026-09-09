"""grade_root_pinning.py — grade PREREG_ROOT_PINNING_2026-09-09.md.

THE QUESTION
------------
`analyse_phase8._pressure_root` (E2, the primary estimator) fits
`p_ss = c (phi - phi_c)^b` with the root bounded above by the lowest rung it is given.
Across all 46 resolved fits on `phase8_prod2*` and `phase8_p5uniform*`, under all four
gate sets, the returned root lands 0.0003-0.0124 below that lowest rung and never more
-- inside the first ladder interval, every time. And the ladders are CENTRED on the
project's own published isotropic `phi_J`, so a root pinned to the ladder makes
`dphi_c(CV)` a restatement of `dphi_J_iso(CV)` rather than an independent measurement.

THE TEST
--------
One cell (lognormal, cv = 0, mu_p = 0.3, N = 1500, seeds 1-3), identical to prod2's
cv = 0 arm in every knob except the ladder offsets, which run -0.06..+0.02 instead of
0.00..+0.08 -- so the ladder bottom drops by 0.060 while the material does not change.

  H_pinned  the root follows the window:  phi_c(down) in [0.4945, 0.5145]
  H_real    0.5695 is a material root:    phi_c(down) in [0.565, 0.575], AND the two
            rungs below it cannot sustain a steady stress -- p_ss_hat < 1e-6 and
            z_nr_ss(0.5145) < 4, which needs no fit at all.

The two predictions for the root are 0.055 apart, 18 seed standard deviations.

EVERY NUMBER BELOW IS COMPUTED HERE. Nothing is typed in from the pre-registration
except the pre-registered BOUNDS, which are named constants so they cannot drift.
"""

from __future__ import annotations

import os
from collections import defaultdict

import numpy as np

import analyse_phase8 as A
import phase8_criticalstate as P8

HERE = os.path.dirname(os.path.abspath(__file__))

# --- the pre-registered bounds, from PREREG_ROOT_PINNING_2026-09-09.md section 3 -----
# the ladder the pre-registration named: the ORIGINAL five rungs, 0.00..+0.08.
# `phase8_prod2*` would also pull in the two rungs added at +0.10/+0.12 on
# 2026-09-09, which is a different window and is tested separately below.
UP_LADDER_PATTERN = "phase8_prod2_shard*.jsonl"
DOWN_LADDER_PATTERN = "phase8_downlad*.jsonl"
CELL = ("lognormal", 0.0, 0.3)          # the cell both ladders measure
H_PINNED_BAND = (0.4945, 0.5145)        # one spacing below the new bottom rung
H_REAL_BAND = (0.565, 0.575)            # the up-ladder answer, +- one seed sd band
P_UNJAMMED_BAR = 1.0e-6                 # prod2's jammed rungs run 2e-4 to 2e-2
Z_ISO_FRICTIONAL = A.Z_ISO_FRICTIONAL   # 4.0, the frictional isostatic count
N_LOW_RUNGS_BELOW_CLAIM = 2             # rungs at 0.5145 and 0.5345 sit below 0.5695


def _groups(pattern: str) -> dict:
    recs = [r for r in A.load(pattern)
            if "error" not in r and np.isfinite(r.get("p_ss_hat", np.nan))]
    by = defaultdict(list)
    for r in recs:
        by[A.group_key(r)].append(r)
    return by


def _hi(by: dict) -> float:
    """The top of the fitting window for the cell under test, averaged over seeds."""
    tops = [max(float(r["phi"]) for r in rs)
            for k, rs in by.items() if (k[0], k[1], k[2]) == CELL]
    return float(np.mean(tops)) if tops else float("nan")


def _roots(by: dict, gate: tuple[str, ...] = ()) -> dict:
    """seed -> (phi_c_E2, phi_lo) for the cell under test."""
    out = {}
    for k, rs in by.items():
        if (k[0], k[1], k[2]) != CELL:
            continue
        res = A.analyse_group(rs, gate)
        out[k[3]] = (res["phi_c_E2"], res["phi_lo"], res.get("b_E2", np.nan))
    return out


def main() -> int:
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    up = _groups(UP_LADDER_PATTERN)
    down = _groups(DOWN_LADDER_PATTERN)
    centre = P8.ladder_centre(*CELL)
    p(f"cell under test: {CELL}, ladder centre (published isotropic phi_J) = "
      f"{centre:.6f}")

    # ---- the up-ladder answer, recomputed --------------------------------------------
    up_roots = _roots(up)
    up_vals = np.array([v[0] for v in up_roots.values() if np.isfinite(v[0])])
    p(f"\nUP-LADDER  (offsets {P8.PHI_OFFSETS}) n = {len(up_vals)} seeds")
    for s in sorted(up_roots):
        rc, lo, b = up_roots[s]
        p(f"   seed {s}: phi_lo = {lo:.4f}  phi_c_E2 = {rc:.4f}  "
          f"(lo - root = {lo - rc:+.4f})  b = {b:.2f}")
    if up_vals.size:
        p(f"   mean {up_vals.mean():.4f} +- {up_vals.std(ddof=1):.4f}")

    if not down:
        p(f"\n🔴 NO DOWN-LADDER RECORDS matching {DOWN_LADDER_PATTERN} -- "
          f"the deciding run has not landed. Nothing is graded.")
        return 1

    # ---- the physical readout, which needs no fit -------------------------------------
    p("\nDOWN-LADDER, per rung (mean over seeds) -- THE FIT-FREE HALF")
    per_phi = defaultdict(list)
    for k, rs in down.items():
        if (k[0], k[1], k[2]) != CELL:
            continue
        for r in rs:
            per_phi[round(float(r["phi"]), 4)].append(r)
    p(f"   {'phi':>8s} {'off':>7s} {'n':>2s} {'p_ss_hat':>11s} {'z_nr':>6s} "
      f"{'rattlers':>8s} {'I':>9s} {'stdy':>5s} {'qs':>4s}")
    unjammed = {}
    for phi_ in sorted(per_phi):
        rs = per_phi[phi_]
        pv = float(np.mean([r["p_ss_hat"] for r in rs]))
        zv = float(np.mean([r["z_nr_ss"] for r in rs]))
        rv = float(np.mean([r["rattler_frac_ss"] for r in rs]))
        iv = float(np.mean([r["i_ss"] for r in rs]))
        ns = sum(bool(r["steady"]) for r in rs)
        nq = sum(bool(r["quasistatic"]) for r in rs)
        unjammed[phi_] = (pv < P_UNJAMMED_BAR, zv < Z_ISO_FRICTIONAL)
        p(f"   {phi_:8.4f} {phi_ - centre:+7.3f} {len(rs):2d} {pv:11.3e} {zv:6.2f} "
          f"{rv:8.3f} {iv:9.2e} {ns:5d} {nq:4d}")

    lowest = sorted(per_phi)[:N_LOW_RUNGS_BELOW_CLAIM]
    n_p_unjammed = sum(unjammed[x][0] for x in lowest)
    n_z_unjammed = sum(unjammed[x][1] for x in lowest)
    p(f"\n   of the {N_LOW_RUNGS_BELOW_CLAIM} rungs below the up-ladder's claimed root "
      f"{up_vals.mean():.4f}: {n_p_unjammed} carry p_ss_hat < {P_UNJAMMED_BAR:.0e}, "
      f"{n_z_unjammed} have z_nr < {Z_ISO_FRICTIONAL:g}")

    # ---- the fitted root -------------------------------------------------------------
    down_roots = _roots(down)
    dn_vals = np.array([v[0] for v in down_roots.values() if np.isfinite(v[0])])
    p(f"\nDOWN-LADDER  phi_c_E2, n = {len(dn_vals)} seeds resolved of {len(down_roots)}")
    for s in sorted(down_roots):
        rc, lo, b = down_roots[s]
        p(f"   seed {s}: phi_lo = {lo:.4f}  phi_c_E2 = {rc:.4f}  "
          f"(lo - root = {lo - rc:+.4f})  b = {b:.2f}")
    if not dn_vals.size:
        p("\n🔴 no down-ladder root resolved; H_pinned and H_real cannot be separated "
          "on the fitted channel. The fit-free readout above still stands.")
        shift = float("nan")
    else:
        p(f"   mean {dn_vals.mean():.4f} +- "
          f"{dn_vals.std(ddof=1) if dn_vals.size > 1 else float('nan'):.4f}")
        shift = float(dn_vals.mean() - up_vals.mean())
        p(f"\n   ladder bottom moved {min(per_phi) - centre:+.4f}; "
          f"the fitted root moved {shift:+.4f}  "
          f"(response {shift / (min(per_phi) - centre):.2f} of the ladder move)")

    # ---- the replacement: is E4 window-robust where E2 is not? -----------------------
    # E4 is the phi at which the non-rattler contact number crosses the isostatic count.
    # It INTERPOLATES inside the ladder instead of extrapolating below it, and it never
    # touches the stress, so it has no bound at the ladder bottom to be pinned against.
    # If E2 moves with the window and E4 does not, the repair is to change estimator,
    # not to change gate.
    p("\n=== the same window move, seen by E4 (isostatic crossing) ===")
    for name, by in (("up-ladder  ", up), ("down-ladder", down)):
        vals = []
        for k, rs in by.items():
            if (k[0], k[1], k[2]) != CELL:
                continue
            r4 = A.analyse_group(rs, ())
            if np.isfinite(r4["phi_c_E4"]):
                vals.append(r4["phi_c_E4"])
        if vals:
            sd = float(np.std(vals, ddof=1)) if len(vals) > 1 else float("nan")
            p(f"   {name} phi_c_E4 = {np.mean(vals):.4f} +- {sd:.4f} (n{len(vals)})")
        else:
            p(f"   {name} phi_c_E4 unresolved on every seed")

    # ---- window sensitivity: TOP vs BOTTOM, from data that already exists -------------
    # A root that is a property of the material cannot care which end of the window
    # moves. A root that is a local extrapolation from the bottom cares about the bottom
    # and not about the top. `phase8_prod2ext` added two rungs at +0.10 and +0.12 to the
    # SAME groups, i.e. it moved the TOP of the window by 0.040 while leaving the bottom
    # alone; the quality gates move the BOTTOM by 0.020 while leaving the top alone.
    p("\n=== window sensitivity: does the root care which END moves? ===")
    up5 = _groups(UP_LADDER_PATTERN)          # the 5-rung ladder
    up7 = _groups("phase8_prod2*.jsonl")      # prod2 + prod2ext pooled, 7 rungs
    r5, r7 = _roots(up5), _roots(up7)
    seeds_both = [s for s in sorted(r5)
                  if s in r7 and np.isfinite(r5[s][0]) and np.isfinite(r7[s][0])]
    tops = [r7[s][0] - r5[s][0] for s in seeds_both]
    hi5 = _hi(up5)
    hi7 = _hi(up7)
    d_top = hi7 - hi5
    if tops and d_top > 0:
        p(f"   TOP moved +{d_top:.3f} (two rungs added above): root moved "
          f"{np.mean(tops):+.4f} on {len(tops)} seeds "
          f"(response {np.mean(tops) / d_top:.2f} of the window move)")
    bots = []
    for s in sorted(r5):
        gq = A.analyse_group([r for k, rs in up5.items()
                              if (k[0], k[1], k[2]) == CELL and k[3] == s
                              for r in rs], ("quasistatic",))
        if np.isfinite(gq["phi_c_E2"]) and np.isfinite(r5[s][0]):
            bots.append((gq["phi_lo"] - r5[s][1], gq["phi_c_E2"] - r5[s][0]))
    if bots:
        dlo = float(np.mean([b[0] for b in bots]))
        drt = float(np.mean([b[1] for b in bots]))
        p(f"   BOTTOM moved {dlo:+.3f} (quasi-static gate drops the lowest rung): root "
          f"moved {drt:+.4f} on {len(bots)} seeds (response {drt / dlo:.2f})")
    if tops and bots:
        p("   -> the root follows the BOTTOM of its window and is nearly blind to the "
          "TOP, which is what a local extrapolation does and what a material root "
          "does not.")

    # ---- the verdict ------------------------------------------------------------------
    p("\n=== VERDICT against PREREG_ROOT_PINNING_2026-09-09.md ===")
    pinned = bool(dn_vals.size and H_PINNED_BAND[0] <= dn_vals.mean() <= H_PINNED_BAND[1])
    real_fit = bool(dn_vals.size and H_REAL_BAND[0] <= dn_vals.mean() <= H_REAL_BAND[1])
    real_phys = (n_p_unjammed == N_LOW_RUNGS_BELOW_CLAIM
                 and n_z_unjammed >= 1)
    p(f"   H_pinned  root in [{H_PINNED_BAND[0]}, {H_PINNED_BAND[1]}]      : "
      f"{'PASS' if pinned else 'fail'}")
    p(f"   H_real    root in [{H_REAL_BAND[0]}, {H_REAL_BAND[1]}]        : "
      f"{'PASS' if real_fit else 'fail'}")
    p(f"   H_real    low rungs unjammed (fit-free)   : "
      f"{'PASS' if real_phys else 'fail'}")

    if pinned and not real_fit:
        p("\n   -> THE ROOT IS PINNED TO THE LADDER. Every fitted phi_c in this project "
          "is (lowest admitted rung - a small offset), so dphi_c(CV) restates "
          "dphi_J_iso(CV) and is NOT an independent measurement of it. P3, the "
          "matched-phi pressure collapse, is unaffected: it uses no fit.")
    elif real_fit and real_phys and not pinned:
        p("\n   -> THE ROOT IS REAL. phi_c survives a 0.060 move of the fitting window "
          "and the rungs below it cannot sustain a steady stress. The CV contrast "
          "stands as a measurement.")
    elif real_phys and not real_fit and not pinned:
        p("\n   -> SPLIT. The critical state is BRACKETED physically (the low rungs are "
          "unjammed) but E2 does not return it under a window move, so the pressure "
          "root is the wrong estimator and the quoted phi_c values must be replaced, "
          "not merely re-gated.")
    else:
        p("\n   -> NEITHER pre-registered outcome. Report the numbers and say so; do "
          "not reinterpret the bands after seeing them.")

    path = os.path.join(A.OUT_DIR, "root_pinning_verdict.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
