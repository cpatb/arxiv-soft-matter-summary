"""phase12_kladder.py -- the K-LADDER / y_S-WIDENING campaign.

PRE-REGISTRATION: PREREG_KLADDER_2026-09-15.md. Blind predictions:
`predict_kladder.py` -> out/predict_kladder.txt. Grade with `grade_kladder.py`.

This file is `phase12_ternary.py` with the grid lifted out into `kladder.py` and the
record fields renamed. The protocol -- `build_and_compress` at p_hat = 1e-3, mu_p = 0.3,
the same `contact_graph` / `rattler_by_class` structure read off the live scene -- is
byte-for-byte the production one, which is what makes the monodisperse control (G0)
able to land on the SAME beta the binary campaign measured and, at seeds 1-3, on the
BIT-IDENTICAL packing the ternary block's T0 built.

THE GRID IS NOT DEFINED HERE. `kladder.build_cells` is the single source of truth, and
`predict_kladder.py` and `grade_kladder.py` read the same function. lessons.md
2026-09-03: two pieces of code that construct the same identifier by different rules
will eventually disagree.

TAG. `kl` at N = 1500 and `kl3k` at N = 3000, chosen to match NEITHER
`phase12_tern_shard*` NOR `phase12_mixr*` NOR `phase12_extr*`, so T0-T6, NA1-NA4,
M0-M6, F1-F3 and X0-X5 are all CLOSED and this campaign cannot re-decide any of them.

Usage:
    python3 phase12_kladder.py --shard 0 --nshards 20 --tag kl
    python3 phase12_kladder.py --shard 0 --nshards 20 --tag kl3k --n 3000
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time

import kladder as kl
from measure import rattler_by_class
from packing import build_and_compress, contact_graph

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


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
    ap.add_argument("--tag", default="kl")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--n", type=int, default=kl.N_PROD)
    ap.add_argument("--mu", type=float, default=kl.MU_P)
    ap.add_argument("--only-arm", default=None,
                    help="restrict to one arm (T0/K/SPLIT/NA2/KX); default all")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    cells, refused = kl.build_cells(args.n)
    if args.only_arm:
        cells = [c for c in cells if c["arm"] == args.only_arm]
    runs = [(c, s) for c in cells for s in c["seeds"]]
    suffix = f"_{args.tag}"
    if args.nshards > 1:
        suffix += f"_shard{args.shard:02d}of{args.nshards:02d}"
    path = os.path.join(OUT_DIR, f"phase12{suffix}.jsonl")
    done = load_done(path)
    mine = sum(1 for i in range(len(runs)) if i % args.nshards == args.shard)
    print(f"phase12 KLADDER: {len(cells)} cells ({len(refused)} refused), {len(runs)} "
          f"runs, {mine} on shard {args.shard}/{args.nshards}, N={args.n}, "
          f"mu_p={args.mu}, {len(done)} done -> {path}", flush=True)

    t0 = time.time()
    nrun = 0
    for idx, (cell, seed) in enumerate(runs):
        if idx % args.nshards != args.shard:
            continue
        key = kl.cell_key(cell["label"], seed, args.n, args.mu)
        if key in done:
            continue
        fam = cell["family"]
        psd = kl.psd_for(cell["d_rel"], cell["y"])
        try:
            st = build_and_compress(cv=0.0, mu_p=args.mu, p_hat=kl.P_HAT, n=args.n,
                                    seed=seed, family=fam, truncate=None, verbose=True)
            rec = json.loads(st.to_json())
            from yade import O                    # noqa: PLC0415 - WSL only
            counts, neigh, _ = contact_graph(O)
            diam = {b.id: 2.0 * float(b.shape.radius) for b in O.bodies if b is not None}
            dd = psd.diameters
            # Kept for continuity with the binary and ternary runners. For K > 2 the
            # two-way split at the geometric mean of the EXTREMES is ambiguous, which is
            # exactly why `by_design_class` below exists and is the one to read.
            rec["by_class"] = rattler_by_class(
                counts, neigh, diam, d_split=math.sqrt(float(dd[0]) * float(dd[-1])))
            per = {}
            for ci, dc in enumerate(dd):
                ids = [i for i, dv in diam.items() if abs(dv - float(dc)) < 1e-12]
                zs = [counts.get(i, 0) for i in ids]
                per[f"c{ci}"] = {
                    "d": float(dc),
                    "n": len(ids),
                    "z_mean": (float(sum(zs)) / len(zs)) if ids else float("nan"),
                    "rattler_frac": ((sum(1 for z in zs if z < 4) / len(ids))
                                     if ids else float("nan")),
                }
            rec["by_design_class"] = per
            rec["key"] = key
            rec["kl_label"] = cell["label"]
            rec["kl_arm"] = cell["arm"]
            rec["kl_K"] = cell["K"]
            rec["kl_y_S_target"] = cell["y_S_target"]
            rec["kl_q"] = cell["q"]
            rec["d_rel"] = list(cell["d_rel"])
            rec["y_vol_design"] = list(cell["y"])
            rec["n_frac_design"] = [float(v) for v in psd.number_fractions]
            rec["d_design"] = [float(v) for v in psd.diameters]
            rec["cv_design"] = psd.cv
            rec["d32_design"] = psd.d32
            rec["d43_design"] = psd.d43
            rec["size_ratio_design"] = psd.size_ratio
            rec["n_design_counts"] = [int(v) for v in psd.class_counts(args.n)]
            rec["min_phase_count"] = kl.MIN_PHASE_COUNT
            rec["campaign"] = "kladder-2026-09-15"
        except Exception as exc:            # noqa: BLE001 - record, never hide
            rec = {"key": key, "kl_label": cell["label"], "kl_arm": cell["arm"],
                   "kl_K": cell["K"], "kl_y_S_target": cell["y_S_target"],
                   "seed": seed, "p_hat_target": kl.P_HAT, "family": fam,
                   "error": repr(exc), "campaign": "kladder-2026-09-15"}
            print(f"CELL FAILED {key}: {exc!r}", flush=True)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
        nrun += 1
    print(f"phase12 kladder shard {args.shard}: ran {nrun} runs in "
          f"{time.time() - t0:.0f} s -> {path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
