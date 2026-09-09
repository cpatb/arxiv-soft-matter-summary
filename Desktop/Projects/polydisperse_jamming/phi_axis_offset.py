"""phi_axis_offset.py — can the lab's silica be on the DILATANT side at all?

THE PROBLEM
-----------
The mechanism in `MECHANISM_2026-09-08.md` turns on the SIGN of `phi - phi_c`. So the
single most important question for applying it to the lab's own experiment is: at the
lab's `phi`, is the material above or below its critical state?

The lab works at `phi = 0.44`. This project has now MEASURED `phi_c` for frictional hard
spheres under constant-volume shear: **0.5695 at mu_p = 0.3, CV = 0**. On its face that
puts the experiment 0.13 BELOW the critical state -- permanently contractive, never
solid-like -- which contradicts the observation the whole project exists to explain.

This module does not resolve that. It BOUNDS it, by asking how far `phi_c` can be pushed
down by the one knob a frictional hard-sphere model has (friction), and then stating what
is left over.

WHY THIS IS NOT JUST THE KNOWN "phi-SCALE OFFSET"
-------------------------------------------------
The Status card already records a phi-axis disagreement (`Drop Impact - Paper 1` §8):
solidification fronts reported at `phi` = 0.42-0.45 here, first solid patches at 0.47 by
Shah et al., a Wyart-Cates `phi_M` of 0.5176 from the repo's own (withdrawn) fit. Those
are three numbers that cannot share an axis. What is new is that a CONTROLLED simulation
now supplies a FLOOR: however much friction you give hard spheres, `phi_J` saturates, so
there is a value below which no frictional hard-sphere pack can have its critical state.
That turns a vague offset into a bound.

    python phi_axis_offset.py
"""

from __future__ import annotations

import csv
import os

import numpy as np

import analyse_phase8 as A

HERE = os.path.dirname(os.path.abspath(__file__))
PHI_LAB = 0.44          # Paper 1 working solid fraction


def friction_ladder() -> list[tuple[float, float, float]]:
    """(mu_p, phi_J, sd) from the COMMITTED isotropic summary, lognormal CV=0, N=1500."""
    path = os.path.join(HERE, "out", "phase2_summary.csv")
    rows = []
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if (r["family"] == "lognormal" and int(r["n_requested"]) == 1500
                    and abs(float(r["cv_nom"])) < 1e-9):
                sd = float(r["phi_J_sd"]) if r["phi_J_sd"] not in ("", "nan") else 0.0
                rows.append((float(r["mu_p"]), float(r["phi_J"]), sd))
    return sorted(rows)


def measured_phi_c_cv0() -> tuple[float, float]:
    """phi_c(mu_p = 0.3, CV = 0) from the constant-volume run, E2 (primary)."""
    recs = [r for r in A.load("phase8_prod2*.jsonl") if "error" not in r]
    vals = []
    for seed in sorted({int(r["seed"]) for r in recs}):
        sel = [r for r in recs
               if int(r["seed"]) == seed and abs(float(r["cv"])) < 1e-9]
        if len(sel) >= 4:
            v = A.analyse_group(sel)["phi_c_E2"]
            if np.isfinite(v):
                vals.append(v)
    return (float(np.mean(vals)),
            float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0)


def main() -> int:
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    phi_c, phi_c_sd = measured_phi_c_cv0()
    lad = friction_ladder()

    p("=== IS THE LAB'S SILICA ABOVE OR BELOW ITS CRITICAL STATE? ===")
    p(f"  lab working point                       phi   = {PHI_LAB:.4f}")
    p(f"  MEASURED critical state (mu_p=0.3, CV=0) phi_c = {phi_c:.4f} +/- {phi_c_sd:.4f}")
    p(f"  phi - phi_c                                    = {PHI_LAB - phi_c:+.4f}")
    p()
    p("  On its face the lab sits 0.13 BELOW the critical state, i.e. permanently")
    p("  CONTRACTIVE -- which would mean the solid-like branch the project exists to")
    p("  explain cannot happen at all. Something in that comparison is wrong, and this")
    p("  module bounds which part.")
    p()
    p("=== HOW FAR CAN FRICTION PUSH phi_c DOWN? (the only knob a hard-sphere model has)")
    p("  Isotropic phi_J against mu_p, from the project's committed A4 ladder:")
    p(f"  {'mu_p':>6s} {'phi_J':>9s} {'sd':>8s} {'drop from mu=0':>15s}")
    base = lad[0][1]
    for mu, pj, sd in lad:
        p(f"  {mu:6.2f} {pj:9.4f} {sd:8.4f} {pj - base:+15.4f}")
    floor = min(pj for _, pj, _ in lad)
    p()
    p(f"  phi_J SATURATES at {floor:.4f}: the step from mu_p = 1 to 3 is "
      f"{lad[-1][1] - lad[-2][1]:+.4f}, i.e. inside the seed scatter. More friction buys")
    p("  essentially nothing. So for FRICTIONAL HARD SPHERES there is a FLOOR.")
    p()
    p(f"  floor - lab phi = {floor - PHI_LAB:+.4f}")
    p("  ⇒ 🔑 NO amount of interparticle friction can bring a hard-sphere critical state")
    p(f"    down to phi = {PHI_LAB}. The gap that friction cannot close is "
      f"{floor - PHI_LAB:.3f} in phi.")
    p()
    p("=== SO ONE OF THESE MUST BE TRUE, and they are testable against each other ===")
    p("  (1) THE LOCAL phi DURING IMPACT IS MUCH HIGHER THAN THE BULK 0.44.")
    need = phi_c / PHI_LAB
    p(f"      Reaching the measured phi_c needs a local densification of "
      f"{need:.3f}x,")
    p(f"      i.e. a volumetric compaction of {100 * (1 - 1 / need):.0f} % at the")
    p("      contact. That is large but it is exactly the Status card's own")
    p("      next-action #2 ('measure the LOCAL phi reached during impact'), and this")
    p("      is the number that measurement has to beat.")
    p()
    p("  (2) THE MATERIAL IS NOT A FRICTIONAL HARD-SPHERE PACK. 800 nm silica at low")
    p("      salt is charge-stabilised: the particles never touch, and a repulsive")
    p("      network is mechanically rigid at far lower solid fractions than a hard")
    p("      sphere pack. This DEM has NO repulsion, so its phi_c is the wrong")
    p("      quantity for that material, and the floor above does not apply to it.")
    p()
    p("  (3) THE EXPERIMENT'S phi IS SYSTEMATICALLY LOW. The Status card already")
    p("      records three irreconcilable numbers on this axis (fronts at 0.42-0.45")
    p("      here, first solid patches at 0.47 in Shah et al., a fitted phi_M of")
    p("      0.5176). A controlled simulation now adds a FLOOR to that list rather")
    p("      than a fourth free number.")
    p()
    p("=== WHAT THIS DOES AND DOES NOT DO TO THE MECHANISM ===")
    p("  DOES NOT TOUCH: the CV CONTRAST. Every claim in MECHANISM_2026-09-08.md is a")
    p("  DIFFERENCE at fixed protocol -- dphi_c/dCV = +0.125, the 43x collapse at")
    p("  matched phi -- and a common offset in the phi axis cancels out of all of them.")
    p()
    p("  DOES TOUCH: any statement about WHICH SIDE the lab's silica is on, including")
    p("  the critical width CV* and the regime map's assignment of the lab to R1. Those")
    p("  need phi_c for THIS material, which this DEM does not measure, and they should")
    p("  be read as conditional on resolution (1), (2) or (3).")
    p()
    p("  🔑 AND IT MAKES ONE CHEAP THING URGENT: option (2) is decidable in this same")
    p("  DEM by adding a repulsive contact law and re-running one ladder. That is a")
    p("  bigger change than a sweep but a much smaller one than a new experiment, and")
    p("  it would say whether a charge-stabilised pack has a critical state anywhere")
    p("  near 0.44.")

    out = os.path.join(HERE, "out", "phi_axis_offset.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
