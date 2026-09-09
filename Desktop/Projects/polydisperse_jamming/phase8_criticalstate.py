"""phase8_criticalstate.py — measure phi_c(CV), the critical-state packing fraction.

WHAT THIS ADDS THAT PHASE 2 CANNOT GIVE
---------------------------------------
`phase2_phiJ.py` measures `phi_J`, the ISOTROPIC jamming point, by extrapolating a
pressure ladder. The mechanism under test (MECHANISM_2026-09-08.md) does not turn on
`phi_J`. It turns on `phi_c`, the CRITICAL-STATE fraction under sustained shear, because
that is the fraction that decides the SIGN of the pore-pressure a sheared saturated
packing generates:

    phi > phi_c   the skeleton must dilate to flow  -> pore fluid into SUCTION
                  -> effective stress rises -> solid-like, jammed, rebounding
    phi < phi_c   the skeleton contracts as it flows -> pore fluid into OVERPRESSURE
                  -> effective stress falls -> liquefied, spreading

`phi_c` appears in the denominator of [M26]'s inertial-dilatant number as `phi - phi_c`,
so a CV dependence of `phi_c` is a CV dependence of ID with a POLE in it, which is a
qualitatively different thing from the smooth `d32^2` factor the permeability channel
supplies. Which of the two dominates is an arithmetic question once `phi_c(CV)` is
measured, and it has never been measured here.

THE MEASUREMENT
---------------
For each (family, CV, mu_p, seed) a LADDER of five volume fractions is sheared at
CONSTANT VOLUME to `gamma_max`, and the steady-state grain pressure `p_ss(phi)` is read.
`phi_c` is the root of `p_ss(phi) = 0`.

The ladder is centred on the ISOTROPIC `phi_J` already published for that cell, with
fixed offsets. That is an efficiency choice, declared before the run: it keeps every
rung inside the near-jamming regime where the harmonic contact law gives `p ~ (phi -
phi_c)^1`, instead of spending half the ladder at pressures where the grains are visibly
soft. It cannot bias the answer toward a particular `phi_c` because the fit is a free
intercept over five rungs and `analyse_phase8.py` REFUSES to report a root that falls
outside the ladder it was fitted on.

Usage
-----
    ~/rl.sh phase8_criticalstate.py --smoke        # one cell, times the protocol
    ~/lp8.sh 12                                    # the production sweep, sharded
"""

from __future__ import annotations

import argparse
import json
import os
import time

from shear import shear_at_constant_volume

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

# --- ladder centres: the project's OWN published isotropic phi_J -------------------
# Source: out/phase2_summary.csv, lognormal family, N = 1500, three seeds, as committed
# on 2026-09-08 after the N-grouping fix. `test_shear.py::test_ladder_centres_match_csv`
# asserts these against that file, so they cannot drift silently.
PHI_J_ISO = {
    ("lognormal", 0.0, 0.0): 0.6331558432715606,
    ("lognormal", 0.0, 0.3): 0.5745332977186427,
    ("lognormal", 0.1, 0.3): 0.5781125394913632,
    ("lognormal", 0.2, 0.3): 0.5880503856512317,
    ("lognormal", 0.3, 0.3): 0.6039853497318549,
    ("lognormal", 0.4, 0.3): 0.6195824881938097,
    ("lognormal", 0.5, 0.3): 0.6368540484855477,
    ("uniform", 0.0, 0.3): 0.5751361139135575,
    ("uniform", 0.5, 0.3): 0.5979189641708998,
}

# Offsets from the ladder centre. Five rungs at 0.020 spacing, 0 to +0.080.
#
# ⚠ MOVED UPWARD 2026-09-08 after the first production run. The ladder was originally
# centred on the isotropic phi_J, on the assumption that the critical state sits near
# it. It does not: at cv = 0 the phi(I -> 0) extrapolation returned phi_c = 0.614-0.616
# against phi_J_iso = 0.5745, i.e. the critical state is ~0.04 ABOVE the isotropic
# jamming point, so five of the seven rungs were unjammed and only two carried signal.
# That is a real result (steady shear densifies a frictional pack past the state
# isotropic compression reaches) and it is also a badly placed sample.
#
# This is a change of SAMPLING, not of hypothesis: the predictions P1-P5 in
# derivations/dilatancy_sign_mechanism.md section 5 are untouched, `analyse_phase8.py`
# still REFUSES a root outside the ladder it was fitted on, and the previous ladder's
# records are kept (they are what located the root). Bracketing a root by moving the
# bracket toward it is a search, not a tuned parameter.
#
# 0.020 spacing is ~4x the largest seed scatter in PHI_J_ISO (0.0048).
PHI_OFFSETS = (0.0, +0.020, +0.040, +0.060, +0.080)

CV_LADDER = (0.0, 0.2, 0.5)
MU_LADDER = (0.3,)
SEEDS = (1, 2, 3)

# N is bounded below by the periodic collider at CV = 0.5 (see packing.insert_packing).
N_PROD = 1500
N_SMOKE = 600

# gamma_max = 1.5 so the averaging window (the last STEADY_FRAC of the strain) is
# gamma in [0.75, 1.5], i.e. most of a unit of strain AFTER the preparation is erased.
# At gamma_max = 1.0 the steadiness gate failed on 4 of 10 locating cells. Reduced from
# 2.0 with the step count cut in proportion once the ladder moved onto the jammed side,
# where the pressure is higher and the inertial number correspondingly smaller.
GAMMA_MAX = 1.5
SHEAR_STEPS = 200_000


def ladder_centre(family: str, cv: float, mu: float) -> float:
    """The published isotropic phi_J for this cell, or the nearest CV in the family.

    Raises rather than guessing: a missing centre means the ladder would be placed on
    an assumption, and this project has been bitten by exactly that before.
    """
    key = (family, round(cv, 3), round(mu, 3))
    if key not in PHI_J_ISO:
        raise KeyError(
            f"no published isotropic phi_J for {key}; add it to PHI_J_ISO from "
            f"out/phase2_summary.csv rather than letting the ladder float"
        )
    return PHI_J_ISO[key]


def cell_key(family: str, cv: float, mu: float, seed: int, phi: float, n: int,
             gamma_max: float, shear_steps: int) -> str:
    return (f"{family}_cv{cv:.3f}_mu{mu:.3f}_s{seed}_phi{phi:.4f}_N{n}"
            f"_g{gamma_max:.2f}_st{shear_steps}")


def load_done(path: str) -> set[str]:
    done = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["key"])
                except Exception:
                    continue
    return done


def build_grid(family: str, cvs, mus, seeds, n: int):
    """(family, cv, mu, seed, phi) tuples, ladder centres resolved."""
    grid = []
    for cv in cvs:
        for mu in mus:
            centre = ladder_centre(family, cv, mu)
            for seed in seeds:
                for off in PHI_OFFSETS:
                    grid.append((family, cv, mu, seed, round(centre + off, 6), n))
    return grid


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="one cheap cell, to time the protocol and check the gates")
    ap.add_argument("--tag", default="")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--family", default="lognormal",
                    choices=["lognormal", "uniform", "bidisperse"])
    ap.add_argument("--n", type=int, default=0)
    ap.add_argument("--cvs", default="")
    ap.add_argument("--mus", default="")
    ap.add_argument("--seeds", default="")
    ap.add_argument("--gamma-max", type=float, default=GAMMA_MAX)
    ap.add_argument("--shear-steps", type=int, default=SHEAR_STEPS)
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    n = args.n or (N_SMOKE if args.smoke else N_PROD)
    cvs = (0.0,) if args.smoke else CV_LADDER
    mus = MU_LADDER
    seeds = (1,) if args.smoke else SEEDS
    if args.cvs:
        cvs = tuple(float(x) for x in args.cvs.split(","))
    if args.mus:
        mus = tuple(float(x) for x in args.mus.split(","))
    if args.seeds:
        seeds = tuple(int(x) for x in args.seeds.split(","))

    grid = build_grid(args.family, cvs, mus, seeds, n)
    if args.smoke:
        # one rung above the centre: a cell that MUST be jammed, so a null result is
        # a protocol failure rather than physics.
        top = round(ladder_centre(args.family, cvs[0], mus[0]) + PHI_OFFSETS[-1], 6)
        grid = [g for g in grid if g[4] == top]
        assert grid, "smoke filter matched no cell -- the ladder rounding moved"

    suffix = ("_smoke" if args.smoke else "") + (f"_{args.tag}" if args.tag else "")
    if args.family != "lognormal":
        suffix += f"_{args.family}"
    if args.nshards > 1:
        suffix += f"_shard{args.shard:02d}of{args.nshards:02d}"
    path = os.path.join(OUT_DIR, f"phase8{suffix}.jsonl")
    done = load_done(path)
    mine = sum(1 for i in range(len(grid)) if i % args.nshards == args.shard)
    print(f"phase8: {len(grid)} cells total, {mine} on shard {args.shard}/"
          f"{args.nshards}, N={n}, gamma_max={args.gamma_max}, "
          f"steps={args.shear_steps}, {len(done)} already done -> {path}", flush=True)

    t0 = time.time()
    nrun = 0
    for idx, (family, cv, mu, seed, phi, ncell) in enumerate(grid):
        if idx % args.nshards != args.shard:
            continue
        key = cell_key(family, cv, mu, seed, phi, ncell, args.gamma_max,
                       args.shear_steps)
        if key in done:
            continue
        try:
            st = shear_at_constant_volume(
                cv=cv, mu_p=mu, phi_target=phi, n=ncell, seed=seed, family=family,
                gamma_max=args.gamma_max, shear_steps=args.shear_steps, verbose=True,
            )
            rec = json.loads(st.to_json())
            rec["key"] = key
            rec["phi_J_iso_centre"] = ladder_centre(family, cv, mu)
        except Exception as exc:                       # noqa: BLE001 - record, do not hide
            rec = {"key": key, "family": family, "cv": cv, "mu_p": mu, "seed": seed,
                   "phi_target": phi, "n_requested": ncell, "error": repr(exc)}
            print(f"CELL FAILED {key}: {exc!r}", flush=True)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
        nrun += 1
    print(f"phase8 shard {args.shard}: ran {nrun} cells in "
          f"{time.time() - t0:.0f} s -> {path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
