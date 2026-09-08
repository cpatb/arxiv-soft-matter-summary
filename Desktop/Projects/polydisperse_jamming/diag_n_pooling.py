"""diag_n_pooling.py -- is the phase-2 analysis pooling two SYSTEM SIZES into one ladder?

WHAT WAS FOUND. `analyse_phase2.py` de-duplicates on a record key that encodes N, and its
own docstring says so -- but the two keys it then GROUPS by do not contain N:

    ladder fit   groups[(family, cv, mu_p, seed)]                    <- no N
    seed scatter scat[(family, cv, mu_p, p_hat_target)]              <- no N

`out/` holds 36 lognormal cells at N = 3000 beside 427 at N = 1500, so some ladders carry
the same `p_hat` twice at two system sizes, and some "seed scatter" cells are measuring a
system-size difference rather than seed noise. The seed scatter is the DENOMINATOR of the
gate that decides which ladders are informative, and the ladder fit produces `b_free`,
which is the project's pre-registered exponent anchor.

WHAT THIS SCRIPT DOES. It recomputes the two statistics with and without N in the keys and
reports the difference, so the size of the defect is a measurement rather than an argument.

STATUS: `analyse_phase2.py` WAS FIXED on 2026-09-08 (N is now in all three grouping keys,
and anchors A4 and A5 are graded per N). This script is kept as the record of what the
defect was worth, and it still runs: the "POOLED" arm below reconstructs the old behaviour
from the same data. If it ever reports 0 ladders mixing N, the data no longer has two
sizes and the diagnosis is moot rather than passing.

THE BAR IS NOT MINE. `b_free` median in [0.6, 1.6] was pre-registered by this project
before 2026-09-08. Whether separating N moves that verdict is the question, and no bar is
touched either way. The two medians are PRINTED below rather than quoted here, because a
number hardcoded in a docstring is exactly what goes stale when the data grows.

Run:  python diag_n_pooling.py
Out:  out/n_pooling.txt  +  out/n_pooling.json
"""

from __future__ import annotations

import collections
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, HERE)

from phase2_phiJ import fit_phi_J  # noqa: E402

# Copied from analyse_phase2.py so this script states the constants it uses rather than
# importing a module with side effects. If they drift there, they must be changed here.
P_HAT_FIT_MAX = 0.0011
EXPONENT_MIN_RANGE_IN_SD = 2.0
B_BAND = (0.6, 1.6)          # pre-registered by this project BEFORE today. Not touched.


def load() -> list[dict]:
    by_key: dict[str, dict] = {}
    for f in sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                k = r.get("key")
                if k and k not in by_key:
                    by_key[k] = r
    return [r for r in by_key.values()
            if not r.get("failed") and r.get("converged")]


def fam(r: dict) -> str:
    return r.get("family", "lognormal")


def ladders(recs: list[dict], *, split_n: bool) -> dict[tuple, list[dict]]:
    g: dict[tuple, list[dict]] = collections.defaultdict(list)
    for r in recs:
        k = (fam(r), r["cv"], r["mu_p"], r["seed"])
        if split_n:
            k = k + (r["n_requested"],)
        g[k].append(r)
    return g


def seed_sd(recs: list[dict], *, split_n: bool) -> tuple[float, int]:
    s: dict[tuple, list[float]] = collections.defaultdict(list)
    for r in recs:
        k = (fam(r), r["cv"], r["mu_p"], r["p_hat_target"])
        if split_n:
            k = k + (r["n_requested"],)
        s[k].append(r["phi"])
    sds = [float(np.std(v, ddof=1)) for v in s.values() if len(v) > 1]
    return (float(np.median(sds)) if sds else float("nan")), len(sds)


def grade(recs: list[dict], *, split_n: bool) -> dict:
    sd, n_cells = seed_sd(recs, split_n=split_n)
    rows, b_informative, n_no_free, n_nonmono, n_flat = [], [], 0, 0, 0
    for key, g in sorted(ladders(recs, split_n=split_n).items()):
        gf = [r for r in g if r["p_hat"] <= P_HAT_FIT_MAX]
        if len(gf) < 2:
            continue
        p = np.array([r["p_hat"] for r in gf])
        f = np.array([r["phi"] for r in gf])
        fit = fit_phi_J(p, f)
        order = np.argsort(p)
        monotone = bool(np.all(np.diff(f[order]) >= 0) or np.all(np.diff(f[order]) <= 0))
        rng = float(f.max() - f.min())
        rows.append({"key": str(key), "n_points": len(gf), "phi_J_lin": fit["phi_J_lin"],
                     "b_free": fit["b_free"], "monotone": monotone, "phi_range": rng})
        if not np.isfinite(fit["b_free"]):
            n_no_free += 1
            continue
        if not monotone:
            n_nonmono += 1
            continue
        if not (rng > EXPONENT_MIN_RANGE_IN_SD * sd):
            n_flat += 1
            continue
        b_informative.append(fit["b_free"])
    b = np.array(b_informative)
    med = float(np.median(b)) if b.size else float("nan")
    return {"split_n": split_n, "seed_sd": sd, "n_scatter_cells": n_cells,
            "n_ladders": len(rows), "n_informative": int(b.size),
            "n_no_free_fit": n_no_free, "n_nonmonotonic": n_nonmono, "n_flat": n_flat,
            "b_median": med, "b_mean": float(b.mean()) if b.size else float("nan"),
            "b_iqr": [float(np.percentile(b, 25)), float(np.percentile(b, 75))]
                     if b.size else [float("nan")] * 2,
            "verdict": "PASS" if (B_BAND[0] <= med <= B_BAND[1]) else "FAIL",
            "rows": rows}


def main() -> int:
    recs = load()
    lad_pool = ladders(recs, split_n=False)
    mixed = {k: v for k, v in lad_pool.items()
             if len({x["n_requested"] for x in v}) > 1}
    s_pool: dict[tuple, list[dict]] = collections.defaultdict(list)
    for r in recs:
        s_pool[(fam(r), r["cv"], r["mu_p"], r["p_hat_target"])].append(r)
    mixed_scatter = {k: v for k, v in s_pool.items()
                     if len({x["n_requested"] for x in v}) > 1}

    pooled = grade(recs, split_n=False)
    split = grade(recs, split_n=True)

    # per-ladder move, on the ladders that actually mix N
    moves = []
    for k, v in sorted(mixed.items()):
        a = [r for r in v if r["p_hat"] <= P_HAT_FIT_MAX]
        b1500 = [r for r in a if r["n_requested"] == 1500]
        if len(a) < 2 or len(b1500) < 2:
            continue
        fa = fit_phi_J(np.array([r["p_hat"] for r in a]),
                       np.array([r["phi"] for r in a]))
        fb = fit_phi_J(np.array([r["p_hat"] for r in b1500]),
                       np.array([r["phi"] for r in b1500]))
        moves.append({"ladder": str(k), "n_pooled": len(a), "n_1500": len(b1500),
                      "phi_J_pooled": fa["phi_J_lin"], "phi_J_1500": fb["phi_J_lin"],
                      "d_phi_J": fa["phi_J_lin"] - fb["phi_J_lin"],
                      "b_pooled": fa["b_free"], "b_1500": fb["b_free"]})
    d = np.array([abs(m["d_phi_J"]) for m in moves])

    out = {"generated_by": "diag_n_pooling.py",
           "n_converged_cells": len(recs),
           "sizes_present": sorted({r["n_requested"] for r in recs}),
           "n_ladder_groups": len(lad_pool), "n_ladders_mixing_N": len(mixed),
           "n_scatter_cells": len(s_pool), "n_scatter_cells_mixing_N": len(mixed_scatter),
           "pooled": {k: v for k, v in pooled.items() if k != "rows"},
           "split": {k: v for k, v in split.items() if k != "rows"},
           "per_ladder_moves": moves,
           "d_phi_J_median": float(np.median(d)) if d.size else float("nan"),
           "d_phi_J_max": float(np.max(d)) if d.size else float("nan"),
           "d_phi_J_max_in_seed_sd": float(np.max(d) / pooled["seed_sd"]) if d.size else
                                     float("nan"),
           "B_BAND": list(B_BAND)}
    with open(os.path.join(OUT, "n_pooling.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)

    L = ["=" * 78,
         "IS THE PHASE-2 ANALYSIS POOLING TWO SYSTEM SIZES INTO ONE LADDER?",
         "=" * 78, "",
         f"{len(recs)} converged cells, sizes present {out['sizes_present']}.",
         "",
         "analyse_phase2.py de-duplicates on a key that encodes N -- its docstring says",
         "so -- but the two keys it GROUPS by do not:",
         "    ladder fit   groups[(family, cv, mu_p, seed)]",
         "    seed scatter scat[(family, cv, mu_p, p_hat_target)]",
         "",
         f"  ladder groups that mix N : {len(mixed)} of {len(lad_pool)}",
         f"  scatter cells that mix N : {len(mixed_scatter)} of {len(s_pool)}",
         "",
         "-" * 78,
         "1. THE DENOMINATOR: median seed scatter in phi",
         "-" * 78,
         f"  pooled over N : {pooled['seed_sd']:.6f}  ({pooled['n_scatter_cells']} cells)",
         f"  split by N    : {split['seed_sd']:.6f}  ({split['n_scatter_cells']} cells)",
         f"  -> moves by {abs(split['seed_sd'] - pooled['seed_sd']) / pooled['seed_sd'] * 100:.1f} %.",
         "  This is the denominator of the gate that decides which ladders are",
         "  informative, so a small move here is worth reporting even though it is small.",
         "",
         "-" * 78,
         "2. THE LADDERS THAT ACTUALLY MIX N",
         "-" * 78,
         f"  {'ladder':<34}{'n_pool':>7}{'n_1500':>7}{'phiJ_pool':>11}{'phiJ_1500':>11}"
         f"{'delta':>9}{'b_pool':>8}{'b_1500':>8}"]
    for m in moves:
        L.append(f"  {m['ladder']:<34}{m['n_pooled']:>7}{m['n_1500']:>7}"
                 f"{m['phi_J_pooled']:>11.5f}{m['phi_J_1500']:>11.5f}"
                 f"{m['d_phi_J']:>+9.5f}"
                 f"{m['b_pooled']:>8.3f}"
                 f"{(('%8.3f' % m['b_1500']) if np.isfinite(m['b_1500']) else '     n/a')}")
    L += ["",
          f"  |delta phi_J|: median {out['d_phi_J_median']:.5f}, max {out['d_phi_J_max']:.5f}",
          f"  in units of the pooled seed scatter: max {out['d_phi_J_max_in_seed_sd']:.2f}",
          "",
          "  NOTE the direction of the cost: splitting N takes ladders from 8 points to 4",
          "  or from 6 to 3, and fit_phi_J returns b as NaN below 4 points. So separating",
          "  N does not simply change b -- it DISQUALIFIES ladders. That is the honest",
          "  trade and it is why this is a diagnosis, not a patch.",
          "",
          "-" * 78,
          "3. THE PRE-REGISTERED EXPONENT ANCHOR (bar: median b_free in "
          f"[{B_BAND[0]}, {B_BAND[1]}], set by this project BEFORE today)",
          "-" * 78]
    for tag, g in (("POOLED (the behaviour BEFORE 2026-09-08)", pooled), ("SPLIT BY N (what ships now)", split)):
        L.append(f"  {tag}")
        L.append(f"    ladders {g['n_ladders']}, informative {g['n_informative']} "
                 f"(excluded: {g['n_no_free_fit']} no free fit, "
                 f"{g['n_nonmonotonic']} non-monotone, {g['n_flat']} flat)")
        L.append(f"    median b = {g['b_median']:.3f}  mean {g['b_mean']:.3f}  "
                 f"IQR [{g['b_iqr'][0]:.3f}, {g['b_iqr'][1]:.3f}]  -> {g['verdict']}")
    same = pooled["verdict"] == split["verdict"]
    L += ["",
          f"  The verdict is {'UNCHANGED' if same else 'DIFFERENT'} under the split.",
          ""]
    if same:
        L += ["  So the pooling was a real defect in the grouping and it does NOT rescue",
              "  the anchor: the harmonic prediction b = 1 fails either way. That is the",
              "  useful outcome -- a defect that could have explained the failure, ruled",
              "  out by measurement rather than by argument.",
              ""]
    else:
        L += ["  The verdict MOVES, so the shipped number depends on a grouping choice",
              "  that nobody made deliberately. Neither value may be quoted until the",
              "  grouping is fixed and the analysis re-run.",
              ""]
    L += ["-" * 78,
          "WHAT THIS DOES NOT SAY",
          "-" * 78,
          "  Not that N = 3000 is wrong, or that phi_J is wrong. The N ladder itself was",
          "  measured separately (analyse_finite_size.py) and the extrapolated phi_J moved",
          "  0.0033 when N doubled. What is measured HERE is only that two system sizes",
          "  land in one fit and in one scatter estimate without anyone choosing that.",
          ""]
    txt = "\n".join(L) + "\n"
    with open(os.path.join(OUT, "n_pooling.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt)
    print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
