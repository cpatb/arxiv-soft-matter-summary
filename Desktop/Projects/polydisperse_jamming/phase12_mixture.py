"""phase12_mixture.py -- the COMPOSITION axis, at fixed size ratio.

PRE-REGISTRATION: PREREG_MIXTURE_2026-09-13.md. Grade with grade_mixture.py.

THE QUESTION
------------
`shape_control.py` established that no single descriptor of the size distribution
collapses the jamming shift, and that the decisive obstruction is the bidisperse
family's NON-MONOTONIC response to CV. `mixture_model.py` proposes the reason: the
bidisperse family is a one-parameter PATH through a two-parameter plane. It is 50:50 by
NUMBER at d10*(1 -/+ CV), so as CV rises its size ratio rises 1.00 -> 3.00 while the
fraction of the SOLID VOLUME carried by the small atom collapses 0.500 -> 0.036. The
family therefore never visits the fines-rich compositions where a pairwise-filling model
puts its maximum, and "ratio 3 packs worst of the three families" is a statement about
composition, not about whether interstitial filling happens.

This campaign puts a grid on the plane the family only ever traced a curve through:

    phi( size ratio r, fines VOLUME fraction X_s ),  r in {2, 3},  X_s in [0, 1]

at the protocol's own p_hat = 1e-3, N = 1500, mu_p = 0.3 rung. No extrapolated phi_J
appears anywhere -- phi is READ, following the COMPANION section of
PREREG_TAIL_2026-09-13.md, which found the fixed-pressure channel is the sound one and
the ladder fit is what was unhealthy.

THE MATCHED PAIR, which is the sharpest single comparison in this project
----------------------------------------------------------------------
At r = 3, X_s = 0.25 the binary has number CV EXACTLY 0.5000 (n_small = 0.9,
d_small = 5/6, d_large = 5/2: mean 1, variance 1/4). The existing `bidisperse` CV = 0.5
cell is r = 3 with number CV exactly 0.5000 as well, at X_s = 0.0357. So the two packs
agree EXACTLY in

    number CV, size ratio, d_max/d_min, d95/d05, number span, and the number of atoms

and differ only in which atom carries the volume. Every descriptor built on those
quantities predicts the same phi. The filling model predicts they differ by ~0.06, which
is ~28x the protocol's own reproducibility. One pair, one measurement, and it cannot be
explained away by a family difference, a friction difference or a truncation choice.

TWO REPRODUCTION GATES, both bit-identity rather than agreement
---------------------------------------------------------------
  * The ANCHOR cell r = 3, X_s = 0.5^3/(0.5^3 + 1.5^3) = 0.035714... is the SAME
    distribution as `BidispersePSD(cv = 0.5)`, reached along the other parameterisation.
    `BinaryPSD.sample` and `BidispersePSD.sample` both build [750 small, 750 large] and
    shuffle with the same generator, so the anchor must reproduce the bidisperse CV = 0.5
    records BIT-IDENTICALLY, not merely within scatter. If it does not, the new family is
    not the same family.
  * Both ENDPOINTS X_s = 0 and X_s = 1 are monodisperse at d = 1 by construction, so they
    must reproduce EACH OTHER and the existing CV = 0 cells bit-identically.
A gate that can only be passed by the code being the same code is worth more than a gate
passed by two numbers being close, and both of these are free.

Usage:
    python3 phase12_mixture.py --shard 0 --nshards 12 --tag mixr3 --ratio 3.0
"""

from __future__ import annotations

import argparse
import json
import os
import time

from packing import build_and_compress
from psd import BinaryPSD

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

P_HAT = 1e-3          # the rung every family record on disk already has
N_PROD = 1500
MU_P = 0.3
SEEDS = (1, 2, 3)

# The composition grid. Chosen to bracket the filling maximum the model puts near
# X_s ~ 0.3 AND to include the matched pair at 0.25, with both monodisperse endpoints.
X_GRID = (0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.60, 1.0)

# The bidisperse family's own composition at CV = 0.5, to the precision Python prints:
# the reproduction anchor. Written as the expression, not a decimal, so it cannot drift.
X_ANCHOR = 0.5 ** 3 / (0.5 ** 3 + 1.5 ** 3)

# Cells with fewer than this many grains of EITHER atom are run and recorded but are
# not gradable: a mixture whose minority phase is a few dozen grains is measuring the
# minority phase's sampling noise. Registered here, applied in grade_mixture.py.
MIN_PHASE_COUNT = 50


def build_grid(ratio: float, n: int, with_anchor: bool):
    xs = list(X_GRID)
    if with_anchor:
        xs.append(X_ANCHOR)
    return [(ratio, x, seed, n) for x in xs for seed in SEEDS]


def cell_key(ratio: float, x: float, seed: int, n: int) -> str:
    return f"mix_r{ratio:.2f}_x{x:.6f}_s{seed}_p{P_HAT:.0e}_N{n}_mu{MU_P:.2f}"


def load_done(path: str) -> set[str]:
    done: set[str] = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["key"])
                except Exception:      # noqa: BLE001 - a truncated last line is not fatal
                    continue
    return done


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="mix")
    ap.add_argument("--ratio", type=float, required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--n", type=int, default=N_PROD)
    ap.add_argument("--no-anchor", action="store_true")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    grid = build_grid(args.ratio, args.n, not args.no_anchor)
    suffix = f"_{args.tag}"
    if args.nshards > 1:
        suffix += f"_shard{args.shard:02d}of{args.nshards:02d}"
    path = os.path.join(OUT_DIR, f"phase12{suffix}.jsonl")
    done = load_done(path)
    mine = sum(1 for i in range(len(grid)) if i % args.nshards == args.shard)
    print(f"phase12 mixture: ratio {args.ratio}, {len(grid)} cells, {mine} on shard "
          f"{args.shard}/{args.nshards}, N={args.n}, {len(done)} done -> {path}",
          flush=True)

    t0 = time.time()
    nrun = 0
    for idx, (ratio, x, seed, ncell) in enumerate(grid):
        if idx % args.nshards != args.shard:
            continue
        key = cell_key(ratio, x, seed, ncell)
        if key in done:
            continue
        # FULL round-trip precision in the family string, not a rounded one. At the
        # anchor composition 0.5^3/(0.5^3+1.5^3) a 6-decimal truncation moves d_small
        # off 0.5 by 1e-6, which is enough to make the packing a different packing and
        # would have destroyed the bit-identity gate the anchor exists to provide.
        fam = f"binary:r{ratio!r}:x{x!r}"
        psd = BinaryPSD(d10=1.0, ratio=ratio, x_small=x)
        try:
            st = build_and_compress(cv=0.0, mu_p=MU_P, p_hat=P_HAT, n=ncell, seed=seed,
                                    family=fam, truncate=None, verbose=True)
            rec = json.loads(st.to_json())
            rec["key"] = key
            rec["mix_ratio"] = ratio
            rec["x_small_vol"] = x
            rec["n_small_design"] = psd.n_small
            rec["d_small_design"] = psd.d_small
            rec["d_large_design"] = psd.d_large
            rec["cv_design"] = psd.cv
            rec["d32_design"] = psd.d32
            rec["is_anchor"] = bool(abs(x - X_ANCHOR) < 1e-12)
            rec["campaign"] = "mixture-2026-09-13"
        except Exception as exc:            # noqa: BLE001 - record, never hide
            rec = {"key": key, "mix_ratio": ratio, "x_small_vol": x, "seed": seed,
                   "p_hat_target": P_HAT, "family": fam, "error": repr(exc),
                   "campaign": "mixture-2026-09-13"}
            print(f"CELL FAILED {key}: {exc!r}", flush=True)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
        nrun += 1
    print(f"phase12 shard {args.shard}: ran {nrun} cells in {time.time() - t0:.0f} s "
          f"-> {path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
