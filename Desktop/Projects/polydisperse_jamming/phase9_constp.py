"""phase9_constp.py — constant-PRESSURE shear runs. Writes out/phase9*.jsonl.

Pre-registration: PREREG_CONSTP_2026-09-12.md. Grader: grade_constp.py.

At constant pressure the steady-state solid fraction IS `phi_c(p)`: no ladder, no root,
no fitting window, no isostatic contact target. See `shear.shear_at_constant_pressure`.

The SHAKEOUT (`--shakeout`) is 6 cells at N = 600 whose whole purpose is to decide
whether the servo is admissible: `phi_ss` must not move with the gain (S2) and must not
remember its preparation (S3). Nothing in it is a measurement of `phi_c`.

    python phase9_constp.py --shakeout --shard 0 --nshards 6
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import time

from shear import SERVO_GAIN, shear_at_constant_pressure

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

# --- the shakeout grid, QUOTED from PREREG_CONSTP_2026-09-12.md §4 -----------------
SHAKE_N = 600
SHAKE_CV = 0.0
SHAKE_MU = 0.3
SHAKE_SEED = 1
SHAKE_P_STAR = 3.0e-4          # a rung of the isotropic ladder, so the answer can be
                               # set beside phi_iso = 0.583314 +- 0.000725 at the same p
SHAKE_PHI_START = (0.55, 0.64)  # one below and one above the expected phi_c
SHAKE_GAINS = (0.25, 0.5, 1.0)
GAMMA_MAX = 3.0                # more strain than the constant-volume protocol needs:
                               # the pack has to travel to phi_c, not just forget its
                               # preparation at fixed volume
SHEAR_STEPS = 200_000          # gdot = 6e-05 at dt = 0.25, so I ~ 3.5e-03 at p*


def cell_key(family, cv, mu, seed, phi0, n, p_star, gain, gamma_max, steps) -> str:
    return (f"{family}_cv{cv:.3f}_mu{mu:.3f}_s{seed}_phi0{phi0:.4f}_N{n}"
            f"_p{p_star:.3e}_K{gain:.2f}_g{gamma_max:.2f}_st{steps}")


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


def shakeout_grid() -> list[tuple]:
    return [("lognormal", SHAKE_CV, SHAKE_MU, SHAKE_SEED, phi0, SHAKE_N, SHAKE_P_STAR,
             k) for phi0, k in itertools.product(SHAKE_PHI_START, SHAKE_GAINS)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shakeout", action="store_true",
                    help="the 6-cell servo shakeout of PREREG_CONSTP §4")
    ap.add_argument("--tag", default="")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--family", default="lognormal",
                    choices=["lognormal", "uniform", "bidisperse"])
    ap.add_argument("--n", type=int, default=0)
    ap.add_argument("--cvs", default="")
    ap.add_argument("--p-stars", default="")
    ap.add_argument("--phi-starts", default="")
    ap.add_argument("--gains", default="")
    ap.add_argument("--seeds", default="")
    ap.add_argument("--mus", default="")
    ap.add_argument("--gamma-max", type=float, default=GAMMA_MAX)
    ap.add_argument("--shear-steps", type=int, default=SHEAR_STEPS)
    ap.add_argument("--dump-dir", default="",
                    help="write per-contact arrays to <dir>/<key>.npz (PREREG_ACTSHAPE)")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    if args.shakeout:
        grid = shakeout_grid()
    else:
        # an explicit grid, so a production run cannot inherit the shakeout's defaults
        need = {"--cvs": args.cvs, "--p-stars": args.p_stars,
                "--phi-starts": args.phi_starts, "--seeds": args.seeds}
        miss = [f for f, v in need.items() if not v]
        if miss:
            raise SystemExit(
                f"a non-shakeout run must state its grid explicitly; missing "
                f"{', '.join(miss)}. The shakeout is N = {SHAKE_N} and is not a "
                f"measurement -- see PREREG_CONSTP_2026-09-12.md §4.")
        cvs = tuple(float(x) for x in args.cvs.split(","))
        mus = tuple(float(x) for x in (args.mus or str(SHAKE_MU)).split(","))
        ps = tuple(float(x) for x in args.p_stars.split(","))
        p0s = tuple(float(x) for x in args.phi_starts.split(","))
        ks = tuple(float(x) for x in (args.gains or str(SERVO_GAIN)).split(","))
        seeds = tuple(int(x) for x in args.seeds.split(","))
        n = args.n or 1500
        grid = [(args.family, cv, mu, s, p0, n, p, k) for cv, mu, s, p0, p, k
                in itertools.product(cvs, mus, seeds, p0s, ps, ks)]

    suffix = ("_shakeout" if args.shakeout else "") + \
        (f"_{args.tag}" if args.tag else "")
    if args.nshards > 1:
        suffix += f"_shard{args.shard:02d}of{args.nshards:02d}"
    path = os.path.join(OUT_DIR, f"phase9{suffix}.jsonl")
    done = load_done(path)
    mine = sum(1 for i in range(len(grid)) if i % args.nshards == args.shard)
    print(f"phase9: {len(grid)} cells total, {mine} on shard {args.shard}/"
          f"{args.nshards}, gamma_max={args.gamma_max}, steps={args.shear_steps}, "
          f"{len(done)} already done -> {path}", flush=True)

    if args.dump_dir:
        os.makedirs(args.dump_dir, exist_ok=True)
    t0 = time.time()
    nrun = 0
    for idx, (family, cv, mu, seed, phi0, ncell, p_star, gain) in enumerate(grid):
        if idx % args.nshards != args.shard:
            continue
        key = cell_key(family, cv, mu, seed, phi0, ncell, p_star, gain,
                       args.gamma_max, args.shear_steps)
        if key in done:
            continue
        try:
            st = shear_at_constant_pressure(
                cv=cv, mu_p=mu, p_hat_target=p_star, n=ncell, seed=seed,
                phi_start_target=phi0, family=family, servo_gain=gain,
                gamma_max=args.gamma_max, shear_steps=args.shear_steps, verbose=True,
                dump_path=(os.path.join(args.dump_dir, key + ".npz")
                           if args.dump_dir else None),
            )
            rec = json.loads(st.to_json())
            rec["key"] = key
            rec["shakeout"] = bool(args.shakeout)
        except Exception as exc:                       # noqa: BLE001 - record, not hide
            rec = {"key": key, "family": family, "cv": cv, "mu_p": mu, "seed": seed,
                   "phi_target": phi0, "n_requested": ncell, "p_hat_target": p_star,
                   "servo_gain": gain, "error": repr(exc)}
            print(f"CELL FAILED {key}: {exc!r}", flush=True)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
        nrun += 1
    print(f"phase9 shard {args.shard}: ran {nrun} cells in "
          f"{time.time() - t0:.0f} s -> {path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
