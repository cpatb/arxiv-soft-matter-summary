"""phase11_tail.py — is the jamming shift set by the NUMBER width or by the TAIL?

PRE-REGISTRATION: PREREG_TAIL_2026-09-13.md.

THE QUESTION
------------
Every headline this project quotes is a derivative with respect to CV, the NUMBER-weighted
coefficient of variation. At matched CV and matched friction the three size-distribution
families give shifts of +0.0621 / +0.0242 / +0.0139 (lognormal / uniform / bidisperse),
a factor 4.5 (`shape_control.py`). So CV is not the control parameter.

The families differ in many ways at once. This campaign changes ONE thing, inside ONE
family, at ONE friction: the extent of the large-particle TAIL.

  Four truncated lognormals, all with realised NUMBER CV = 0.500, cut at
  d/d10 <= U for U in {2, 3, 5, 8}. The nominal sigma is re-solved per window so the
  realised number CV is matched to 1e-3. Their VOLUME-weighted widths are then
  0.280 / 0.355 / 0.454 / 0.483 -- a factor 1.72 at matched number width.

  NUMBER-WIDTH HYPOTHESIS  phi_J is the same for all four.
  TAIL HYPOTHESIS          phi_J tracks the volume-weighted width, spreading by
                           ~K * 0.203 = 0.016 with K = 0.0787 from the family fit.

The two predictions are 3.2x apart against the registered 0.005 bar, so the test is
informative before it is run.

WHY IT MATTERS: the production truncation (0.2, 5.0) was an IMPLEMENTATION CHOICE, and
the descriptor sweep shows it sits on the rising part of the cv_vol curve (0.280 at U=2,
0.454 at U=5, 0.483 at U=8, saturating only past 8). If the tail hypothesis is right, the
project's headline +0.053 is partly a property of that constant. If the number-width
hypothesis is right, it is not, and CV's failure across families must be something else.

⚠ THE REFEREE'S GATE, pre-registered here because it has already bitten this project.
`phi_J_lin` is a straight line through a 4-rung pressure ladder, and its se cannot see a
ladder that is not a straight line. The mu_p = 1.0 arm on disk is 0/6 monotone with
r2 down to 0.002 and one UNPHYSICAL negative slope, and it was carrying a claim until an
adversarial review found it. Every cell here records `phi_monotone`, `r2_lin` and
`slope_lin`, and the grader REFUSES any arm failing the health bar in section 3 of the
pre-registration.

Usage:
    python3 phase11_tail.py --shard 0 --nshards 12 --tag tail
"""

from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
from scipy.optimize import brentq

from packing import build_and_compress
from psd import make_psd

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

# The pressure ladder phi_J is extrapolated from -- phase2_phiJ.P_HAT_LADDER verbatim,
# so these records sit on the same ladder as every phi_J on disk.
P_HAT_LADDER = (1e-3, 3e-4, 1e-4, 3e-5)

TRUNC_LO = 0.2
UPPERS = (2.0, 3.0, 5.0, 8.0)
TARGET_CV = 0.50
SEEDS = (1, 2, 3)
N_PROD = 1500
MU_P = 0.3
FAMILY = "lognormal"

# The CV = 0 control: truncation cannot move a monodisperse pack, because there is no
# tail to cut. Run at the two extreme windows only.
CONTROL_UPPERS = (2.0, 8.0)

# solver settings for the matched-CV nominal sigma
SOLVE_N = 400_000
SOLVE_SEED = 7
SOLVE_TOL = 1e-4


def realised_cv(nominal_cv: float, upper: float) -> float:
    rng = np.random.default_rng(SOLVE_SEED)
    x = np.asarray(make_psd(FAMILY, 1.0, nominal_cv).sample(
        SOLVE_N, rng, truncate=(TRUNC_LO, upper)), float)
    return float(x.std() / x.mean())


def nominal_for(upper: float, target: float = TARGET_CV) -> float:
    """The nominal CV whose TRUNCATED draw has realised number CV = `target`.

    Truncation removes the large tail, which is where most of the number-CV lives, so the
    nominal must be raised to compensate -- at U = 2 it takes 0.724 to realise 0.500.
    Solved rather than tabulated so the grid is reproducible from the code.
    """
    if target <= 0:
        return 0.0
    return float(brentq(lambda c: realised_cv(c, upper) - target, 0.4, 3.0,
                        xtol=SOLVE_TOL))


def vol_cv(nominal_cv: float, upper: float) -> float:
    rng = np.random.default_rng(SOLVE_SEED)
    x = np.asarray(make_psd(FAMILY, 1.0, nominal_cv).sample(
        SOLVE_N, rng, truncate=(TRUNC_LO, upper)), float)
    w = x ** 3
    w = w / w.sum()
    mv = float((w * x).sum())
    return float(np.sqrt((w * (x - mv) ** 2).sum()) / mv)


def build_grid(n: int):
    """(upper, nominal_cv, target_cv, seed, p_hat) tuples.

    Ordered so round-robin sharding spreads the pressure rungs, which differ in cost
    (the low-pressure rungs take the most iterations).
    """
    grid = []
    for cv_target in (TARGET_CV, 0.0):
        uppers = UPPERS if cv_target > 0 else CONTROL_UPPERS
        for upper in uppers:
            nom = nominal_for(upper, cv_target)
            for seed in SEEDS:
                for p in P_HAT_LADDER:
                    grid.append((upper, nom, cv_target, seed, p, n))
    return grid


def cell_key(upper, nom, cv_t, seed, p, n) -> str:
    return (f"tail_U{upper:.1f}_cvt{cv_t:.2f}_nom{nom:.4f}_s{seed}"
            f"_p{p:.0e}_N{n}_mu{MU_P:.2f}")


def load_done(path: str) -> set[str]:
    done: set[str] = set()
    if os.path.exists(path):
        for line in open(path, encoding="utf-8"):
            try:
                done.add(json.loads(line)["key"])
            except Exception:
                continue
    return done


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="tail")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--n", type=int, default=N_PROD)
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    grid = build_grid(args.n)
    suffix = f"_{args.tag}"
    if args.nshards > 1:
        suffix += f"_shard{args.shard:02d}of{args.nshards:02d}"
    path = os.path.join(OUT_DIR, f"phase11{suffix}.jsonl")
    done = load_done(path)
    mine = sum(1 for i in range(len(grid)) if i % args.nshards == args.shard)
    print(f"phase11 tail: {len(grid)} cells, {mine} on shard {args.shard}/"
          f"{args.nshards}, N={args.n}, {len(done)} done -> {path}", flush=True)

    t0 = time.time()
    nrun = 0
    for idx, (upper, nom, cv_t, seed, p, ncell) in enumerate(grid):
        if idx % args.nshards != args.shard:
            continue
        key = cell_key(upper, nom, cv_t, seed, p, ncell)
        if key in done:
            continue
        try:
            st = build_and_compress(cv=nom, mu_p=MU_P, p_hat=p, n=ncell, seed=seed,
                                    family=FAMILY, truncate=(TRUNC_LO, upper),
                                    verbose=True)
            rec = json.loads(st.to_json())
            rec["key"] = key
            rec["trunc_lo"] = TRUNC_LO
            rec["trunc_hi"] = upper
            rec["cv_nominal"] = nom
            rec["cv_target"] = cv_t
            rec["cv_vol_design"] = vol_cv(nom, upper)
            rec["campaign"] = "tail-2026-09-13"
        except Exception as exc:            # noqa: BLE001 - record, never hide
            rec = {"key": key, "trunc_hi": upper, "cv_nominal": nom,
                   "cv_target": cv_t, "seed": seed, "p_hat_target": p,
                   "error": repr(exc), "campaign": "tail-2026-09-13"}
            print(f"CELL FAILED {key}: {exc!r}", flush=True)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
        nrun += 1
    print(f"phase11 shard {args.shard}: ran {nrun} cells in {time.time() - t0:.0f} s "
          f"-> {path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
