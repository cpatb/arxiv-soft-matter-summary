"""analyse_phase7.py — the impact experiment with the aggregate actually jammed.

Phases 4 and 5 confined the carved aggregate at `p_imp/Lambda`, which was ~100x weaker
than the pressure it was carved at, so it expanded and lost 50-100% of its contacts
before impact (`diag_carve_release.py`). Phase 7 confines at `p_cap = p_carve` and
gates every cell on `z_retention = Z_launch / Z_carved`.

Two things are reported, in this order, because the second is only meaningful if the
first holds:

  1. THE GATE. How many cells were actually jammed when they hit. A cell below
     Z_RETENTION_MIN is reported, not silently dropped, and its numbers are excluded
     from the response statistics.

  2. THE RESPONSE, on the admitted cells only: does a genuinely jammed aggregate
     rebound, and does polydispersity change that?

The comparison against Phase 5 is the point. Same carve pressures, same velocities, same
grains, one design change.
"""

from __future__ import annotations

import csv
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

Z_RETENTION_MIN = 0.80
COR_REBOUND = 0.05        # same bar as Phase 5, set from the C = 0 noise floor
SPREAD_COHERENT = 1.5
COR_ELASTIC_CEILING = 0.8189

# PHYSICAL ADMISSIBILITY. A passive body cannot leave faster than it arrived, so
# COR <= 1 and KE_after <= KE_before. The confining force in Phase 7 is applied every
# timestep by a PyRunner, including during and after the collision, so it can do WORK on
# the aggregate and eject it. The first run of this analysis reported COR = +2.15 with a
# seed sd of 3.05 -- unphysical, and with the two seeds disagreeing by more than the
# mean. Cells violating these bounds are diagnosing that artefact, not a rebound, and
# are reported separately instead of being counted.
COR_MAX_PHYSICAL = 1.0
KE_MAX_PHYSICAL = 1.0


def load(pattern):
    by_key, bad = {}, 0
    for f in sorted(glob.glob(os.path.join(OUT, pattern))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if r.get("failed"):
                    bad += 1
                    continue
                if r.get("key"):
                    by_key[r["key"]] = r
    return list(by_key.values()), bad


def main() -> int:
    recs, bad = load("phase7*.jsonl")
    print(f"N: {len(recs)} cells ({bad} failed)")
    if not recs:
        print("no phase7 records yet")
        return 1

    lines = ["=" * 84,
             "PHASE 7 -- impact on an aggregate that is still jammed at launch",
             "=" * 84, "",
             "Phases 4 and 5 confined at p_imp/Lambda, ~100x below the carve pressure,",
             "so the aggregate expanded and lost 50-100% of its contacts BEFORE impact.",
             "Phase 7 confines at p_cap = p_carve, which makes Lambda = J.", ""]

    # ---- 1. the gate ------------------------------------------------------
    ret = np.array([r["z_retention"] for r in recs], float)
    jam = [r for r in recs if r.get("jammed_at_launch")]
    lines.append("=" * 84)
    lines.append(f"THE GATE: was it jammed at impact?  (Z_launch/Z_carved >= "
                 f"{Z_RETENTION_MIN})")
    lines.append("=" * 84)
    lines.append(f"  retention over {len(recs)} cells: median {np.nanmedian(ret):.3f}, "
                 f"range [{np.nanmin(ret):.3f}, {np.nanmax(ret):.3f}]")
    lines.append(f"  cells passing: {len(jam)} / {len(recs)}")
    if len(jam) < len(recs):
        lines.append("  cells FAILING the gate (reported, and excluded below):")
        for r in sorted((x for x in recs if not x.get("jammed_at_launch")),
                        key=lambda x: x["z_retention"]):
            lines.append(f"      cv={r['cv']:.2f} p_carve={r['carve_p_hat']:.0e} "
                         f"u0={r['u0']:.3f} retention={r['z_retention']:.3f}")
    # the Phase 4/5 comparison
    old, _ = load("phase5*.jsonl")
    if old:
        oret = [r["z_launch"] / r["z_initial"] for r in old
                if r.get("z_initial") and r.get("z_launch") is not None
                and r["z_initial"] > 0]
        if oret:
            lines.append("")
            lines.append(f"  For comparison, Phase 5 (confinement ~100x too weak) had "
                         f"median retention {np.median(oret):.3f} over {len(oret)} cells.")
    lines.append("")

    if not jam:
        lines.append("  No cell passed the gate. Nothing below can be interpreted.")
        txt = "\n".join(lines)
        print(txt)
        with open(os.path.join(OUT, "phase7_verdict.txt"), "w", encoding="utf-8") as fh:
            fh.write(txt + "\n")
        return 1

    # ---- 2. the response, admitted cells only -----------------------------
    groups: dict[tuple, list[dict]] = {}
    for r in jam:
        groups.setdefault((r["cv"], r["carve_p_hat"], r["u0"]), []).append(r)
    rows = []
    for (cv, cp, u0), g in groups.items():
        row = {"cv": cv, "carve_p_hat": cp, "u0": u0, "n": len(g),
               "J": float(np.mean([x["j_impact"] for x in g]))}
        for k in ("cor", "spread", "ke_ratio", "z_carved", "z_launch",
                  "z_retention", "z_peak"):
            v = np.array([x.get(k, np.nan) for x in g], float)
            v = v[np.isfinite(v)]
            row[k] = float(v.mean()) if v.size else float("nan")
            row[k + "_sd"] = float(v.std(ddof=1)) if v.size > 1 else float("nan")
        rows.append(row)
    rows.sort(key=lambda r: (r["cv"], r["J"]))

    fp = os.path.join(OUT, "phase7_summary.csv")
    with open(fp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {fp}")

    lines.append("=" * 84)
    lines.append("THE RESPONSE, on jammed cells only")
    lines.append("=" * 84)
    lines.append(f"  {'CV':>5} {'p_carve':>9} {'u0':>7} {'J':>8} {'n':>3} "
                 f"{'ret':>6} {'COR':>9} {'sd':>8} {'spread':>8} {'Z peak':>7}")
    for r in rows:
        lines.append(f"  {r['cv']:5.2f} {r['carve_p_hat']:9.1e} {r['u0']:7.3f} "
                     f"{r['J']:8.3g} {r['n']:3d} {r['z_retention']:6.2f} "
                     f"{r['cor']:+9.4f} {r['cor_sd']:8.4f} {r['spread']:8.3f} "
                     f"{r['z_peak']:7.2f}")
    lines.append("")

    def admissible(r):
        return (r["cor"] <= COR_MAX_PHYSICAL
                and (not np.isfinite(r.get("ke_ratio", np.nan))
                     or r["ke_ratio"] <= KE_MAX_PHYSICAL))

    inadmissible = [r for r in rows if not admissible(r)]
    reb = [r for r in rows if admissible(r) and r["cor"] > COR_REBOUND
           and r["spread"] < SPREAD_COHERENT]
    lines.append("=" * 84)
    lines.append("DOES A JAMMED AGGREGATE REBOUND?")
    lines.append("=" * 84)
    lines.append(f"  bar: COR > {COR_REBOUND} AND spread < {SPREAD_COHERENT}")
    lines.append(f"  detection ceiling {COR_ELASTIC_CEILING:.4f} (one elastic sphere, "
                 "same damping)")
    if inadmissible:
        lines.append(f"  ⚠ {len(inadmissible)} cell group(s) are PHYSICALLY INADMISSIBLE "
                     f"(COR > {COR_MAX_PHYSICAL} or KE gained) and are excluded:")
        for r in sorted(inadmissible, key=lambda x: -x["cor"]):
            lines.append(f"      CV={r['cv']:.2f} p_carve={r['carve_p_hat']:.0e} "
                         f"J={r['J']:.3g} COR={r['cor']:+.4f} (sd {r['cor_sd']:.4f}) "
                         f"KE={r['ke_ratio']:.3f}")
        lines.append("    The confining PyRunner acts every step, including after impact,")
        lines.append("    so at high p_carve it does work on the aggregate and throws it")
        lines.append("    off the wall. That is a driven body, not a rebounding one.")
        lines.append("")
    lines.append(f"  cells passing: {len(reb)} / {len(rows)} "
                 f"({len(inadmissible)} inadmissible)")
    # Summary statistics must be taken over ADMISSIBLE rows only; the first version
    # reported "COR reaches +2.1515" from a row it had just excluded.
    adm = [r for r in rows if admissible(r)]
    cors = np.array([r["cor"] for r in adm]) if adm else np.array([np.nan])
    # A mean smaller than its own seed scatter is not a measurement.
    reb = [r for r in reb
           if not np.isfinite(r["cor_sd"]) or r["cor"] > 2.0 * r["cor_sd"]]
    unresolved = [r for r in rows if admissible(r) and r["cor"] > COR_REBOUND
                  and r["spread"] < SPREAD_COHERENT and r not in reb]
    if unresolved:
        lines.append(f"  {len(unresolved)} further cell group(s) clear the bar but are "
                     "NOT RESOLVED (mean < 2x seed sd):")
        for r in unresolved:
            lines.append(f"      CV={r['cv']:.2f} J={r['J']:.3g} "
                         f"COR={r['cor']:+.4f} +/- {r['cor_sd']:.4f}")
        lines.append("")
    if reb:
        lines.append(f"  YES. Largest admissible COR is {cors.max():+.4f}.")
        for r in sorted(reb, key=lambda x: -x["cor"]):
            lines.append(f"      CV={r['cv']:.2f} J={r['J']:.3g} COR={r['cor']:+.4f} "
                         f"spread={r['spread']:.3f}")
        lines.append("")
        lines.append("  => Confining the aggregate at the pressure it was carved at is")
        lines.append("     the difference between an impact on a jammed solid and an")
        lines.append("     impact on a loose cloud. Phases 4 and 5 measured the latter.")
    else:
        lines.append(f"  NO. The largest COR on a jammed cell is {cors.max():+.4f}, "
                     f"{COR_ELASTIC_CEILING / max(cors.max(), 1e-9):.0f}x below the")
        lines.append("  elastic ceiling. So the Phase 4/5 conclusion SURVIVES the")
        lines.append("  correction: even a properly confined, still-jammed cohesionless")
        lines.append("  aggregate does not rebound. The confinement error changed the")
        lines.append("  interpretation but not the answer.")

    # does CV matter, among jammed cells?
    lines.append("")
    lines.append("  CV dependence among jammed cells:")
    for cv in sorted({r["cv"] for r in rows}):
        g = [r for r in rows if r["cv"] == cv]
        c = np.array([r["cor"] for r in g])
        s = np.array([r["spread"] for r in g])
        lines.append(f"    CV={cv:.2f}: COR {c.min():+.4f} to {c.max():+.4f}, "
                     f"spread {s.min():.3f} to {s.max():.3f}  (n={len(g)} cells)")

    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(OUT, "phase7_verdict.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
