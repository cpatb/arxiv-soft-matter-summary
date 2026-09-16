"""phase2_shape.py -- the SHAPE campaign runner. PREREG_SHAPE_2026-09-18.md.

A SIBLING of `phase2_phiJ.py`, not an edit to it. `phase2_phiJ.py` is the runner ten
campaigns on disk were written by; adding a second axis to its cell key would change the
key of every future sphere cell and make the resumability of those campaigns depend on a
flag none of them pass. This file imports the fit, the ladders and the resume logic from
it and owns only what is new: the alpha axis, the alpha-aware cell key, and the tag.

WHAT IS NEW ABOUT THE CELL, in one place
-----------------------------------------
`alpha` is passed straight through to `packing.build_and_compress`, which builds a
two-sphere clump of the same VOLUME-EQUIVALENT diameter (so every PSD descriptor is
invariant along the axis), scales the material density so the clump's mass is the union
mass, computes phi from the analytic union volume rather than `utils.porosity()`, and
builds the contact graph over GRAINS rather than member spheres.

alpha = 1.0 builds a SINGLE SPHERE. That is a correction to the prereg forced by
`scripts/probe_clump_s0.sh`, and the reason is in `packing.PackingState.alpha`: a clump
of two COINCIDENT spheres raises four interactions per grain pair, so it is not the
sphere this project has run. The alpha ladder therefore has a SEAM at its first rung and
gate S0 had to be restated; both are registered in the prereg's amendment block.

THE TAG. `shp`, which matches none of `mixr*`, `extr*`, `tern*`, `kl*`, `ug*`, `fab*`,
and the shards land in `out/shape/` rather than the shared `out/phase2*.jsonl` glob that
ten scripts read. `grade_shape.py` asserts both STRUCTURALLY rather than by eye.

Usage (from `scripts/launch_shape.sh`, never by hand)
-----------------------------------------------------
    python3 phase2_shape.py --tag shp --family lognormal --alpha 1.5 \
        --cvs 0.0,0.1,... --seeds 1,2,3 --phats 1e-3,3e-4,1e-4,3e-5 --n 1500 \
        --shard 0 --nshards 20
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import packing
from packing import build_and_compress
from phase2_phiJ import CV_LADDER, OUT_DIR, P_HAT_LADDER, SEEDS, load_done
from shape import ALPHA_LADDER

# The tag every shard of this campaign carries. Asserted structurally in grade_shape.py
# rather than left to the launcher to get right.
SHAPE_TAG = "shp"


def cell_key(cv: float, mu: float, seed: int, p_hat: float, n: int,
             family: str, alpha: float) -> str:
    """Resume key. `alpha` is ALWAYS in it, including at alpha = 1.

    Deliberately NOT backward-compatible with `phase2_phiJ.cell_key`: the alpha = 1
    cells of this campaign are a reproduction of the sphere records, and a key that
    collided with them would let a resume silently serve the sphere record back as
    the shape campaign's own answer -- which is precisely what gate S0a exists to
    test, so it must not be able to pass by aliasing.
    """
    base = f"cv{cv:.3f}_mu{mu:.3f}_s{seed}_p{p_hat:.3e}_N{n}_a{alpha:.3f}"
    key = base if family == "lognormal" else f"{base}_{family}"
    cap = packing.max_iter()
    return key if cap == packing.MAX_ITER else f"{key}_it{cap}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default=SHAPE_TAG)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--family", default="lognormal",
                    choices=["lognormal", "uniform", "bidisperse"])
    ap.add_argument("--n", type=int, default=0, help="override N")
    ap.add_argument("--alphas", default="",
                    help="comma-separated aspect ratios, overriding shape.ALPHA_LADDER")
    ap.add_argument("--mus", default="0.3")
    ap.add_argument("--cvs", default="")
    ap.add_argument("--seeds", default="")
    ap.add_argument("--phats", default="")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    n = args.n or 1500
    alphas = (tuple(float(x) for x in args.alphas.split(','))
              if args.alphas else ALPHA_LADDER)
    cvs = tuple(float(x) for x in args.cvs.split(',')) if args.cvs else CV_LADDER
    mus = tuple(float(x) for x in args.mus.split(','))
    seeds = tuple(int(x) for x in args.seeds.split(',')) if args.seeds else SEEDS
    ladder = (tuple(float(x) for x in args.phats.split(','))
              if args.phats else P_HAT_LADDER)

    suffix = f"_{args.tag}_{args.family}"
    if args.nshards > 1:
        suffix += f"_shard{args.shard:02d}of{args.nshards:02d}"
    path = os.path.join(OUT_DIR, f"phase2{suffix}.jsonl")
    done = load_done(path)
    total = len(alphas) * len(cvs) * len(mus) * len(seeds) * len(ladder)
    mine = (total + args.nshards - 1 - args.shard) // args.nshards
    # N before and after every filter, printed, not inferred.
    print(f"shape: N_cells_total={total} N_on_this_shard={mine} "
          f"N_already_done={len(done)} shard={args.shard}/{args.nshards} "
          f"family={args.family} alphas={alphas} cvs={cvs} seeds={seeds} "
          f"phats={ladder} N_grains={n} -> {path}", flush=True)

    t0 = time.time()
    idx = 0
    n_run = 0
    n_skipped = 0
    n_failed = 0
    for alpha in alphas:
        for cv in cvs:
            for mu in mus:
                for seed in seeds:
                    for p_hat in ladder:
                        idx += 1
                        if (idx - 1) % args.nshards != args.shard:
                            continue
                        key = cell_key(cv, mu, seed, p_hat, n, args.family, alpha)
                        if key in done:
                            n_skipped += 1
                            continue
                        try:
                            st = build_and_compress(cv=cv, mu_p=mu, p_hat=p_hat, n=n,
                                                    seed=seed, verbose=True,
                                                    family=args.family, alpha=alpha)
                            rec = json.loads(st.to_json())
                            n_run += 1
                        except Exception as exc:      # noqa: BLE001
                            print(f"CELL FAILED {key}: {exc!r}", flush=True)
                            rec = {"failed": True, "error": repr(exc), "cv": cv,
                                   "mu_p": mu, "seed": seed, "p_hat_target": p_hat,
                                   "n_requested": n, "family": args.family,
                                   "alpha": alpha}
                            n_failed += 1
                        rec["key"] = key
                        rec["tag"] = args.tag
                        with open(path, "a", encoding="utf-8") as fh:
                            fh.write(json.dumps(rec, sort_keys=True) + "\n")
                        el = time.time() - t0
                        print(f"  [{idx}/{total}] elapsed {el/60:.1f} min", flush=True)
    print(f"shape done in {(time.time()-t0)/60:.1f} min: "
          f"N_run={n_run} N_skipped_already_done={n_skipped} N_failed={n_failed}",
          flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
