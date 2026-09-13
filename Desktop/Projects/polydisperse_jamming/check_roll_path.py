"""check_roll_path.py — three gates on the rolling-resistance path, before any campaign.

PRE-REGISTRATION: PREREG_ROLLING_2026-09-13.md (written alongside this).

The arm changes the CONTACT LAW, which is the most invasive change this project has made
to `packing.py`. Three things must hold and each one has failed somewhere in this project
before:

  R-A  THE DEFAULT PATH IS UNTOUCHED. `eta_roll = None` must return bit-identical
       numbers to the pre-patch module on the same cell. A refactor that moves a
       trajectory invalidates every record on disk.

  R-B  THE CONTROL REPRODUCES THE DEFAULT. `CohFrictMat` with cohesion off and
       `etaRoll = 0` (moment always at its zero limit, i.e. NO rolling resistance) is
       the SAME PHYSICS as `FrictMat` through a different code path. It cannot be
       bit-identical -- different law, different integration order -- but it must agree
       on `phi_J` well inside the seed scatter. If it does not, the new path is not a
       controlled extension of the old one and nothing measured on it is comparable.

  R-C  THE KNOB IS NOT INERT. `etaRoll > 0` must MOVE `phi`. This is the gate that
       matters most, because the one way this arm fails silently is
       `always_use_moment_law` defaulting to False (verified by probe_roll_api.py):
       with no bonded contacts the moment law would then never fire, every cell would
       run, and every number would be the cohesionless one. A knob whose non-default
       setting returns the default's answer is inert, not verified.

R-B and R-C are the two halves of the same requirement and must BOTH pass: a broken path
that ignored `etaRoll` entirely would pass R-B and fail R-C; one that got the control
wrong would pass R-C and fail R-B.

Usage (WSL, under scripts/env.sh):
    python3 check_roll_path.py --sidecar packing_pre_roll_src.txt
"""

from __future__ import annotations

import argparse
import importlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TMP_MODULE = "packing_pre_roll"

# One cheap cell, at the pressure the ladders use and a CV the campaign will run.
CELL = dict(cv=0.0, mu_p=0.3, p_hat=3e-4, n=600, seed=1)

# The seed scatter this project measures on phi_J at N = 600-1500 is 0.0007-0.0048
# (out/phase2_summary.csv). R-B's bar is the LOOSE end of that, because a single cell
# carries a single realisation.
R_B_TOL = 0.005
# R-C: etaRoll = 0.3 must move phi by more than R-B's own tolerance, or the two gates
# cannot be told apart.
R_C_MIN = 0.005
ETA_TEST = 0.3

CONTINUOUS = ("phi", "p_hat", "z_all", "z_nr", "rattler_frac", "unbalanced")


def fetch_old(sidecar_name: str) -> str:
    """The pre-patch packing.py, from a sidecar written on the Windows side.

    ~/pj is a file copy, not a checkout, so `git show` is unavailable exactly where Yade
    is. Generate with: git show HEAD:<path>/packing.py > packing_pre_roll_src.txt
    """
    path_in = os.path.join(HERE, sidecar_name)
    if not os.path.exists(path_in):
        raise SystemExit(f"sidecar {path_in} not found -- generate it on the Windows side")
    src = open(path_in, encoding="utf-8").read()
    if "eta_roll" in src:
        raise SystemExit(f"{sidecar_name} already knows about eta_roll -- wrong ref")
    if "def build_and_compress(" not in src:
        raise SystemExit(f"{sidecar_name} does not look like packing.py")
    out = os.path.join(HERE, f"{TMP_MODULE}.py")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(src)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sidecar", default="packing_pre_roll_src.txt")
    args = ap.parse_args()

    path = fetch_old(args.sidecar)
    fails = []
    try:
        old = importlib.import_module(TMP_MODULE)
        import packing as new
        assert not hasattr(old, "interaction_loop"), "the old module must not have it"
        assert hasattr(new, "interaction_loop"), "the new module must have it"
        print(f"cell: {CELL}")
        a = json.loads(old.build_and_compress(**CELL, verbose=True).to_json())
        b = json.loads(new.build_and_compress(**CELL, verbose=True).to_json())
        c = json.loads(new.build_and_compress(**CELL, eta_roll=0.0,
                                              verbose=True).to_json())
        d = json.loads(new.build_and_compress(**CELL, eta_roll=ETA_TEST,
                                              verbose=True).to_json())
    finally:
        if os.path.exists(path):
            os.remove(path)

    # ---------------- R-A ----------------
    print("\n=== R-A  the default path is bit-identical to the pre-patch module ===")
    print(f"{'field':18s}{'pre-patch':>24s}{'post':>24s}  identical")
    bad = []
    for f in CONTINUOUS:
        same = a.get(f) == b.get(f)
        print(f"{f:18s}{a.get(f)!r:>24}{b.get(f)!r:>24}  {'yes' if same else 'NO'}")
        if not same:
            bad.append(f)
    if bad:
        print(f"🔴 R-A FAIL on {bad} -- the patch moved a trajectory on the DEFAULT "
              f"path, so every record on disk is in question")
        fails.append("R-A")
    else:
        print("✅ R-A PASS")

    # ---------------- R-B ----------------
    print("\n=== R-B  CohFrictMat with etaRoll = 0 reproduces FrictMat ===")
    dphi = abs(c["phi"] - b["phi"])
    print(f"   phi: FrictMat {b['phi']:.6f}   CohFrictMat(etaRoll=0) {c['phi']:.6f}   "
          f"|d| = {dphi:.6f}  against {R_B_TOL}")
    print(f"   z_nr: {b['z_nr']:.3f} vs {c['z_nr']:.3f};  "
          f"rattlers {b['rattler_frac']:.3f} vs {c['rattler_frac']:.3f};  "
          f"converged {b['converged']} vs {c['converged']}")
    if dphi > R_B_TOL or not c["converged"]:
        print("🔴 R-B FAIL -- the new path is not a controlled extension of the old one")
        fails.append("R-B")
    else:
        print("✅ R-B PASS")

    # ---------------- R-C ----------------
    print(f"\n=== R-C  etaRoll = {ETA_TEST} is NOT INERT ===")
    move = c["phi"] - d["phi"]
    print(f"   phi: etaRoll=0 {c['phi']:.6f}   etaRoll={ETA_TEST} {d['phi']:.6f}   "
          f"move = {move:+.6f}  (need |move| > {R_C_MIN})")
    print(f"   z_nr {c['z_nr']:.3f} -> {d['z_nr']:.3f};  "
          f"rattlers {c['rattler_frac']:.3f} -> {d['rattler_frac']:.3f};  "
          f"converged {d['converged']}")
    if d["phi"] == c["phi"]:
        print("🔴 R-C FAIL -- IDENTICAL phi: the moment law never fired. Check "
              "always_use_moment_law and momentRotationLaw.")
        fails.append("R-C")
    elif abs(move) < R_C_MIN:
        print(f"🔴 R-C FAIL -- the knob moves phi by only {abs(move):.6f}, inside R-B's "
              f"own tolerance, so the two gates cannot be told apart")
        fails.append("R-C")
    elif move < 0:
        print("⚠ R-C passes on magnitude but the SIGN is wrong: rolling resistance "
              "should make the pack LOOSER (phi lower), not denser. Reported, not "
              "graded -- the sign is a prediction of the pre-registration, not a gate "
              "on the code.")
        print("✅ R-C PASS (magnitude)")
    else:
        print(f"✅ R-C PASS -- rolling resistance lowers phi by {move:.4f}, the "
              f"predicted direction")

    print("\n" + ("FAILURES: " + str(fails) if fails else "ALL THREE GATES PASS"))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
