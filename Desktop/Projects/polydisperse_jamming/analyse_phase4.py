"""analyse_phase4.py — does the boundary condition change what an impacting pack does?

Grades the Phase-4 sweep against the scope guard written into SPEC §5 BEFORE the runs:

    At Lambda = 12.5-200 the capillary case is EXPECTED to sit close to the free case.
    If it does, that IS the finding -- capillary confinement is irrelevant during
    impact -- and it must be reported rather than tuned away.

So the test is not "does capillary differ from free" but "WHERE, if anywhere, does the
capillary case depart from free, and is the lab's Lambda range on the free side?"

Every number here inherits the Phase-4 limitation: the model runs at a granular inertial
number I ~ d/D_agg ~ 0.08 against the experiment's 4e-4, so this compares boundary
conditions at matched I and does not predict the experiment's regime boundary.
"""

from __future__ import annotations

import csv
import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

# The lab's own Lambda window, from SPEC §2.1 (We/4 at U0 = 1 to 4 m/s).
LAB_LAMBDA = (12.5, 200.0)

OBS = [("spread", "lateral spread / initial"),
       ("largest_frac_peak", "largest backbone fraction at peak"),
       ("z_peak", "Z at peak contact count"),
       ("rattler_frac_peak", "rattler fraction at peak"),
       ("cor", "coefficient of restitution"),
       ("ke_ratio", "KE retained"),
       ("z_initial", "Z as carved"),
       ("z_launch", "Z at launch")]


def load() -> list[dict]:
    by_key: dict[str, dict] = {}
    n_lines = dropped_failed = dropped_boom = 0
    for f in sorted(glob.glob(os.path.join(OUT, "phase4*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                n_lines += 1
                if r.get("failed"):
                    dropped_failed += 1
                    continue
                if r.get("exploded"):
                    dropped_boom += 1
                    continue
                k = r.get("key")
                if k:
                    by_key[k] = r
    print(f"N: {n_lines} lines -> {len(by_key)} unique cells "
          f"({dropped_failed} failed, {dropped_boom} exploded, dropped)")
    return list(by_key.values())


def lam_plot(r: dict) -> float:
    """Lambda for ordering. free = inf, caged = 0; both are LIMITS, not values."""
    if r["bc"] == "free":
        return math.inf
    if r["bc"] == "caged":
        return 0.0
    return float(r["lam"])


def main() -> int:
    recs = load()
    if not recs:
        print("no phase4 records yet")
        return 1

    groups: dict[tuple, list[dict]] = {}
    for r in recs:
        groups.setdefault((r["cv"], r["bc"], r["lam"]), []).append(r)

    rows = []
    for (cv, bc, lam), g in groups.items():
        row = {"cv": cv, "bc": bc, "lam": lam, "n_seeds": len(g),
               "lam_sort": lam_plot(g[0]),
               "i_inertial": float(np.nanmean([x["i_inertial"] for x in g])),
               "n_grains": float(np.nanmean([x["n_grains"] for x in g]))}
        for key, _ in OBS:
            v = np.array([x.get(key, np.nan) for x in g], float)
            v = v[np.isfinite(v)]
            row[key] = float(v.mean()) if v.size else float("nan")
            row[key + "_sd"] = float(v.std(ddof=1)) if v.size > 1 else float("nan")
        rows.append(row)

    rows.sort(key=lambda r: (r["cv"], r["lam_sort"]))
    fp = os.path.join(OUT, "phase4_summary.csv")
    with open(fp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {fp}")

    lines = ["=" * 78,
             "PHASE 4 -- the drop-like boundary",
             "=" * 78, "",
             "LIMITATION carried by every number below: the model sits at a granular",
             f"inertial number I ~ d/D_agg ~ "
             f"{np.nanmean([r['i_inertial'] for r in rows]):.3f}, against 4e-4 for a",
             "2.5 mm drop of 1 um grains. Boundary conditions are compared at matched I;",
             "no regime boundary here transfers to the experiment.", ""]

    for cv in sorted({r["cv"] for r in rows}):
        g = [r for r in rows if r["cv"] == cv]
        lines.append(f"--- CV = {cv:g}  (N_grains ~ {g[0]['n_grains']:.0f}) ---")
        lines.append(f"  {'bc':<10} {'Lambda':>9} {'n':>3} " +
                     " ".join(f"{lbl[:11]:>12}" for _, lbl in OBS[:5]))
        for r in sorted(g, key=lambda x: x["lam_sort"]):
            lam_s = ("inf" if math.isinf(r["lam_sort"])
                     else f"{r['lam']:.4g}")
            lines.append(f"  {r['bc']:<10} {lam_s:>9} {r['n_seeds']:>3} " +
                         " ".join(f"{r[k]:12.4f}" for k, _ in OBS[:5]))
        lines.append("")

    # --- the pre-registered scope guard --------------------------------------
    lines.append("=" * 78)
    lines.append("PRE-REGISTERED SCOPE GUARD (SPEC §5)")
    lines.append("=" * 78)
    lines.append("'At Lambda = 12.5-200 the capillary case is expected to sit close to")
    lines.append(" the free case. If it does, THAT is the finding.'")
    lines.append("")
    for cv in sorted({r["cv"] for r in rows}):
        g = {(_r["bc"], _r["lam"]): _r for _r in rows if _r["cv"] == cv}
        free = next((v for k, v in g.items() if k[0] == "free"), None)
        caged = next((v for k, v in g.items() if k[0] == "caged"), None)
        if free is None:
            continue
        lines.append(f"CV = {cv:g}:")
        for (bc, lam), r in sorted(g.items(), key=lambda kv: lam_plot(kv[1])):
            if bc != "capillary":
                continue
            in_lab = LAB_LAMBDA[0] <= lam <= LAB_LAMBDA[1]
            parts = []
            for key, _ in ("spread", ""), ("largest_frac_peak", ""):
                fv, cvv = free.get(key, np.nan), r.get(key, np.nan)
                sd = free.get(key + "_sd", np.nan)
                if np.isfinite(fv) and np.isfinite(cvv):
                    d = cvv - fv
                    nsd = abs(d) / sd if (np.isfinite(sd) and sd > 0) else np.inf
                    parts.append(f"{key}: {cvv:.3f} vs free {fv:.3f} "
                                 f"(delta {d:+.3f}, {nsd:.1f} sd)")
            lines.append(f"  Lambda = {lam:7.4g}"
                         f"{'  <-- IN THE LAB WINDOW' if in_lab else ''}")
            for p in parts:
                lines.append(f"      {p}")
        if caged is not None:
            lines.append(f"  caged (Lambda -> 0): spread {caged['spread']:.3f} vs "
                         f"free {free['spread']:.3f}")
        lines.append("")

    txt = "\n".join(lines)
    print(txt)
    fp = os.path.join(OUT, "phase4_verdict.txt")
    with open(fp, "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    print(f"wrote {fp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
