"""phase10_strainclock.py — the strain clock: frictionless prep, friction on, shear.

PRE-REGISTRATION: PREREG_STRAINCLOCK_2026-09-13.md. Read §1 before changing anything
here; the grid in §2 is fixed and this module is where it is written down.

WHY A NEW DRIVER AND NOT --offsets ON phase8
--------------------------------------------
`phase8_criticalstate.py` places its rungs at UNIFORM offsets from a ladder centre it
looks up in `PHI_J_ISO`, and refuses non-uniform spacing on purpose. This campaign needs
ABSOLUTE fractions, per CV, derived from two DIFFERENT measurements (the constant-pressure
`phi_c(p*)` and the isotropic frictionless `phi_J`), with a per-CV step count so every
cell shears at the same rate. Bending phase8 to that would have meant loosening a guard
that exists for a reason, and it would have mixed `cv-shear-v2` and `cv-shear-v3-fprep`
records under one tag. The two record sets must never pool: the preparation is the
independent variable.

WHAT ONE CELL IS
----------------
  stages A-C  insert + compress to `phi_target` + relax, at mu_p = 0 (FRICTIONLESS)
  ---------   read p_hat_prep, z_nr_prep, rattler_frac_prep; then setContactFriction(0.3)
              with a read-back on every material and every real interaction
  stage D     constant-volume simple shear to gamma_max = 3 at gdot = 1.17e-04,
              sampled 1500 times (Delta gamma = 0.002)

Below `phi_0(CV)` the frictionless pack is unjammed, so `p_hat_prep` ~ 0 and
`z_nr_prep` is NaN with `rattler_frac_prep` = 1.0. That is the OCR -> infinity premise,
NOT a failed cell; G1 in the pre-registration grades it as such.

Usage
-----
    python3 phase10_strainclock.py --smoke                    # one cheap cell
    python3 phase10_strainclock.py --cvs 0.0 --shard 0 --nshards 5 --tag sc
"""

from __future__ import annotations

import argparse
import json
import os
import time

from shear import shear_at_constant_volume

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

# --- the ladder anchors, both MEASURED, both with their source named ----------------
# PHI_MU: constant-pressure phi_c at p* = 3e-04, rate-matched, N = 1000, 6 cells each.
#   out/constp_ratematch.txt, "blocks" section, pooled over both phi starts.
# PHI_0: isotropic frictionless phi_J, lognormal, N = 1500, 3 seeds.
#   out/dst_window.txt, "lognormal: the two isotropic ladders", mu_p = 0.0 rows.
# ⚠ These are in DIFFERENT pressure gauges (finite p* vs p -> 0). See §2 caveat 1.
PHI_MU = {0.0: 0.591127, 0.2: 0.601989, 0.5: 0.644447}
PHI_0 = {0.0: 0.632371, 0.2: 0.645392, 0.5: 0.690401}
PHI_0_SD = {0.0: 0.000722, 0.2: 0.000685, 0.5: 0.004665}

# Rungs in the reduced coordinate x = (phi - phi_mu) / (phi_0 - phi_mu).
# x < 0 is below the frictional critical line (nothing should jam);
# x = 0.85 is close enough to phi_0 that at CV 0.5 (sd 0.0047) the prep MAY be jammed,
# which is why G1 grades only x <= 0.60.
X_RUNGS = (-0.15, 0.10, 0.35, 0.60, 0.85)
X_GRADED_MAX = 0.60

# Per-CV step counts, set from the MEASURED dt so gdot is matched at 1.17e-04 with
# gamma_max = 3. Same numbers as scripts/launch_constp_ratematch.sh, whose realised
# spread was 1.001x. dt varies by seed at CV 0.5, so the realised spread is reported.
SHEAR_STEPS = {0.0: 102_000, 0.2: 200_000, 0.5: 480_000}

GAMMA_MAX = 3.0
N_SAMPLES = 1500          # Delta gamma = 0.002 over gamma_max = 3
MU_SHEAR = 0.3
MU_PREP = 0.0
N_PROD = 1000             # matches the constant-pressure ladder PHI_MU came from
N_SMOKE = 600
SEEDS = (1, 2, 3)
CV_LADDER = (0.0, 0.2, 0.5)
FAMILY = "lognormal"


def rungs_for(cv: float) -> list[float]:
    """The five absolute phi for this CV, rounded to 4 dp (the grid, as registered)."""
    if cv not in PHI_MU:
        raise KeyError(f"no measured phi_mu for CV = {cv}; this grid is pre-registered "
                       f"and a new level needs a measurement, not an interpolation")
    w = PHI_0[cv] - PHI_MU[cv]
    return [round(PHI_MU[cv] + x * w, 4) for x in X_RUNGS]


def cell_key(cv: float, seed: int, phi: float, n: int, steps: int) -> str:
    return (f"sc_{FAMILY}_cv{cv:.3f}_mu{MU_SHEAR:.3f}_prep{MU_PREP:.3f}"
            f"_s{seed}_phi{phi:.4f}_N{n}_g{GAMMA_MAX:.2f}_st{steps}")


def load_done(path: str) -> set[str]:
    done: set[str] = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["key"])
                except Exception:
                    continue
    return done


def build_grid(cvs, seeds, n: int):
    """(cv, seed, phi, x, steps) tuples.

    Ordered CV-OUTERMOST so that round-robin sharding hands each worker one cell from
    each CV: the CV 0.5 cells are 2x the wall time of the CV 0 ones, and a shard that
    drew only CV 0.5 would set the campaign's wall clock on its own.
    """
    grid = []
    for cv in cvs:
        for seed in seeds:
            for x, phi in zip(X_RUNGS, rungs_for(cv)):
                grid.append((cv, seed, phi, x, SHEAR_STEPS[cv], n))
    return grid


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="one cell at N=600 and a short shear, to time the protocol")
    ap.add_argument("--tag", default="")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--n", type=int, default=0)
    ap.add_argument("--cvs", default="")
    ap.add_argument("--seeds", default="")
    ap.add_argument("--n-samples", type=int, default=N_SAMPLES)
    ap.add_argument("--gamma-max", type=float, default=GAMMA_MAX)
    ap.add_argument("--steps-scale", type=float, default=1.0,
                    help="multiply every per-CV step count (smoke/timing only; any "
                         "value but 1.0 changes gdot and is recorded in the key)")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    n = args.n or (N_SMOKE if args.smoke else N_PROD)
    cvs = tuple(float(x) for x in args.cvs.split(",")) if args.cvs else CV_LADDER
    seeds = tuple(int(x) for x in args.seeds.split(",")) if args.seeds else SEEDS
    if args.smoke:
        cvs, seeds = (0.0,), (1,)

    grid = build_grid(cvs, seeds, n)
    if args.smoke:
        grid = [g for g in grid if abs(g[3] - 0.60) < 1e-9]      # one jammed rung
        assert grid, "smoke filter matched no cell -- X_RUNGS moved"

    suffix = ("_smoke" if args.smoke else "") + (f"_{args.tag}" if args.tag else "")
    if args.nshards > 1:
        suffix += f"_shard{args.shard:02d}of{args.nshards:02d}"
    path = os.path.join(OUT_DIR, f"phase10{suffix}.jsonl")
    done = load_done(path)
    mine = sum(1 for i in range(len(grid)) if i % args.nshards == args.shard)
    print(f"phase10 strainclock: {len(grid)} cells, {mine} on shard {args.shard}/"
          f"{args.nshards}, N={n}, gamma_max={args.gamma_max}, "
          f"n_samples={args.n_samples}, {len(done)} already done -> {path}", flush=True)

    t0 = time.time()
    nrun = 0
    for idx, (cv, seed, phi, x, steps, ncell) in enumerate(grid):
        if idx % args.nshards != args.shard:
            continue
        steps = int(round(steps * args.steps_scale))
        key = cell_key(cv, seed, phi, ncell, steps)
        if key in done:
            continue
        try:
            st = shear_at_constant_volume(
                cv=cv, mu_p=MU_SHEAR, phi_target=phi, n=ncell, seed=seed,
                family=FAMILY, gamma_max=args.gamma_max, shear_steps=steps,
                mu_prep=MU_PREP, n_samples=args.n_samples, verbose=True,
            )
            rec = json.loads(st.to_json())
            rec["key"] = key
            rec["x_rung"] = x
            rec["phi_mu_anchor"] = PHI_MU[cv]
            rec["phi_0_anchor"] = PHI_0[cv]
            rec["campaign"] = "strainclock-2026-09-13"
        except Exception as exc:               # noqa: BLE001 - record, never hide
            rec = {"key": key, "cv": cv, "seed": seed, "phi_target": phi,
                   "x_rung": x, "n_requested": ncell, "error": repr(exc),
                   "campaign": "strainclock-2026-09-13"}
            print(f"CELL FAILED {key}: {exc!r}", flush=True)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
        nrun += 1
    print(f"phase10 shard {args.shard}: ran {nrun} cells in {time.time() - t0:.0f} s "
          f"-> {path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
