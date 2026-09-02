"""analyse_phase5.py — the impact phase diagram in (J, CV).

Reads out/phase5*.jsonl and reports, for each boundary condition, how the response
varies with the impact number

    J = p_imp / p_carve = rho u0^2 / p

The question it answers is the one the whole project started from: is there ANY region
of this space where the aggregate behaves as a rebounding solid rather than a spreading
one, and where does the crossover sit?

Two response measures, chosen because they are opposite ends of the same behaviour:
  spread     lateral extent after impact, normalised to before. 1.0 = held together.
  cor        coefficient of restitution of the centre of mass. > 0 = rebounded.
  ke_ratio   kinetic energy retained. The cleanest single number for "elastic or not".
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

# A cell counts as rebounding if the centre of mass leaves the wall with at least this
# fraction of its incoming speed. 0.01 is deliberately generous: the point is to find
# ANY signal, and Phase 4 found |COR| < 0.035 everywhere at a single J.
# A REBOUND requires two things, and the first version only tested one.
#
# Testing COR alone produced 16/34 "rebounds": three caged CV = 0.5 cells at COR
# 0.15-0.17 whose spread was 3.6-4.4 (the aggregate SHATTERED and the centre-of-mass
# recoil was measuring debris leaving), and thirteen more at COR 0.010-0.023, which is
# below the 0.026 noise floor that the C = 0 cells of Phase 6 establish for this
# observable.
#
# So the bar is now STRICTER than the one that produced the false positive, not looser:
# the body must recoil AND still be a body.
COR_REBOUND = 0.05        # was 0.01; the C = 0 cohesion cells give |COR| up to 0.026
SPREAD_COHERENT = 1.5     # spread beyond this is a shattered aggregate, not a rebound

# Was the aggregate still jammed when it hit? Phase 5 confined at p_imp/Lambda, which is
# ~100x weaker than the carve pressure, so some aggregates expanded and lost their
# contacts during the settle (diag_carve_release.py). The retention distribution is
# BIMODAL, not uniformly bad: a spike near 0 and a cluster at 0.85-1.10 (see
# figures/F6_jamming_gate). So the sweep is not wholesale invalid -- it is a mixture, and
# the response statistics must be computed on the jammed subset with the discarded
# fraction reported.
Z_RETENTION_MIN = 0.80

# What COR this pipeline WOULD report for a perfectly elastic body, measured by
# test_cor_anchor.py through the same code path at the production damping of 0.2:
COR_ELASTIC_CEILING = 0.8189


def load():
    by_key, n_lines, bad = {}, 0, 0
    for f in sorted(glob.glob(os.path.join(OUT, "phase5*.jsonl"))):
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
                if r.get("failed") or r.get("exploded"):
                    bad += 1
                    continue
                if r.get("key"):
                    # Legacy records carry cor = -v_z/u0 (a perfect bounce reads -1).
                    # They are identified by the missing convention stamp.
                    if r.get("cor_convention") != "positive_is_rebound":
                        if isinstance(r.get("cor"), (int, float)):
                            r["cor"] = -r["cor"]
                        r["cor_convention"] = "negated_on_read"
                    by_key[r["key"]] = r
    print(f"N: {n_lines} lines -> {len(by_key)} unique cells ({bad} failed/exploded)")
    return list(by_key.values())


def main() -> int:
    recs = load()
    if not recs:
        print("no phase5 records yet")
        return 1

    # --- the jamming gate, applied before anything is interpreted ------------
    def retention(r):
        zc, zl = r.get("z_initial"), r.get("z_launch")
        if zc and zl is not None and zc > 0:
            return zl / zc
        return float("nan")

    for r in recs:
        r["z_retention"] = retention(r)
    ret = np.array([r["z_retention"] for r in recs], float)
    jammed = [r for r in recs if np.isfinite(r["z_retention"])
              and r["z_retention"] >= Z_RETENTION_MIN]
    unjammed = [r for r in recs if r not in jammed]
    print(f"JAMMING GATE: {len(jammed)}/{len(recs)} cells kept "
          f"{Z_RETENTION_MIN:.0%}+ of their as-carved coordination at launch; "
          f"{len(unjammed)} did not and are EXCLUDED from the response statistics "
          f"(median retention overall {np.nanmedian(ret):.3f})")
    recs = jammed
    if not recs:
        print("no cell was jammed at launch; nothing can be interpreted")
        return 1

    groups: dict[tuple, list[dict]] = {}
    for r in recs:
        groups.setdefault((r["bc"], r["carve_p_hat"], r["u0"], r["cv"]), []).append(r)

    rows = []
    for (bc, cp, u0, cv), g in groups.items():
        row = {"bc": bc, "carve_p_hat": cp, "u0": u0, "cv": cv, "n": len(g),
               "J": float(np.mean([x["j_impact"] for x in g]))}
        for k in ("spread", "cor", "ke_ratio", "z_initial", "z_peak",
                  "rattler_frac_peak", "largest_frac_peak"):
            v = np.array([x.get(k, np.nan) for x in g], float)
            v = v[np.isfinite(v)]
            row[k] = float(v.mean()) if v.size else float("nan")
            row[k + "_sd"] = float(v.std(ddof=1)) if v.size > 1 else float("nan")
        rows.append(row)
    rows.sort(key=lambda r: (r["bc"], r["cv"], r["J"]))

    fp = os.path.join(OUT, "phase5_summary.csv")
    with open(fp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {fp}")

    lines = ["=" * 86,
             "PHASE 5 -- the impact phase diagram in (J, CV)",
             "=" * 86, "",
             "J = rho u0^2 / p_carve. J << 1: pack stiff against the impact.",
             "J >> 1: impact overwhelms the confinement.",
             "",
             "Phase 4 sampled a single J = 4 and found no rebound anywhere. This sweep",
             "asks whether that holds across four decades of J.",
             "",
             "GATED: only cells that were still jammed at launch "
             f"(Z_launch/Z_carved >= {Z_RETENTION_MIN}) are counted. The rest expanded",
             "out of contact during the settle because the confinement was ~100x weaker",
             "than the carve pressure; see diag_carve_release.py and F6_jamming_gate.", ""]

    for bc in sorted({r["bc"] for r in rows}):
        lines.append(f"--- boundary: {bc} ---")
        lines.append(f"  {'CV':>5} {'p_carve':>9} {'u0':>7} {'J':>9} {'n':>3} "
                     f"{'spread':>8} {'COR':>9} {'KE kept':>9} {'Z init':>7} {'Z peak':>7}")
        for r in [x for x in rows if x["bc"] == bc]:
            lines.append(f"  {r['cv']:5.2f} {r['carve_p_hat']:9.1e} {r['u0']:7.3f} "
                         f"{r['J']:9.3g} {r['n']:3d} {r['spread']:8.3f} "
                         f"{r['cor']:+9.4f} {r['ke_ratio']:9.4f} "
                         f"{r['z_initial']:7.2f} {r['z_peak']:7.2f}")
        lines.append("")

    # --- does anything rebound, anywhere? ------------------------------------
    lines.append("=" * 86)
    lines.append("DOES ANYTHING REBOUND?")
    lines.append("=" * 86)
    reb = [r for r in rows if r["cor"] > COR_REBOUND
           and r["spread"] < SPREAD_COHERENT]
    recoil_only = [r for r in rows if r["cor"] > COR_REBOUND
                   and r["spread"] >= SPREAD_COHERENT]
    lines.append(f"  bar: COR > {COR_REBOUND} (deliberately generous)")
    lines.append("  detection ceiling: one perfectly elastic sphere measured through")
    lines.append(f"  this same path at the production damping returns "
                 f"COR = {COR_ELASTIC_CEILING:.4f}")
    lines.append("  (test_cor_anchor.py). So a rebound WAS detectable if it occurred.")
    lines.append(f"  and coherence: spread < {SPREAD_COHERENT} "
                 "(otherwise COR measures debris, not a body)")
    lines.append(f"  cells passing BOTH: {len(reb)} / {len(rows)}")
    if recoil_only:
        lines.append(f"  cells with recoil but NO coherence: {len(recoil_only)} "
                     "-- these SHATTERED and are not rebounds:")
        for r in sorted(recoil_only, key=lambda x: -x["cor"]):
            lines.append(f"      {r['bc']:>10} CV={r['cv']:.2f} J={r['J']:.3g} "
                         f"COR={r['cor']:+.4f} spread={r['spread']:.3f}")
    if reb:
        lines.append(f"  {'bc':>10} {'CV':>5} {'J':>9} {'COR':>9} {'spread':>8}")
        for r in sorted(reb, key=lambda x: -x["cor"]):
            lines.append(f"  {r['bc']:>10} {r['cv']:5.2f} {r['J']:9.3g} "
                         f"{r['cor']:+9.4f} {r['spread']:8.3f}")
        js = [r["J"] for r in reb]
        lines.append(f"  -> rebound is found, and it lives at J = "
                     f"{min(js):.3g} to {max(js):.3g}")
    else:
        cors = np.array([r["cor"] for r in rows])
        alljs = np.array([r["J"] for r in rows], float)
        alljs = alljs[np.isfinite(alljs) & (alljs > 0)]
        # State the range ACTUALLY covered, not the range the sweep was designed for.
        dec = (np.log10(alljs.max() / alljs.min()) if alljs.size > 1 else 0.0)
        lines.append(f"  -> NO. COR spans [{cors.min():+.4f}, {cors.max():+.4f}] over "
                     f"{len(rows)} cells covering J = {alljs.min():.3g} to "
                     f"{alljs.max():.3g} ({dec:.1f} decades).")
        lines.append("     A cohesionless frictional aggregate with no interstitial fluid")
        lines.append("     does not store and return impact energy ANYWHERE in this space.")
        lines.append(f"     The largest |COR| seen is {np.abs(cors).max():.4f}, which is "
                     f"{COR_ELASTIC_CEILING / max(np.abs(cors).max(), 1e-9):.0f}x below")
        lines.append("     the elastic ceiling, so this is a measurement, not a null "
                     "observable.")

    # --- where does J actually control the response? --------------------------
    lines.append("")
    lines.append("=" * 86)
    lines.append("DOES J ORGANISE THE RESPONSE AT ALL?")
    lines.append("=" * 86)
    for bc in sorted({r["bc"] for r in rows}):
        for cv in sorted({r["cv"] for r in rows}):
            g = sorted([r for r in rows if r["bc"] == bc and r["cv"] == cv],
                       key=lambda x: x["J"])
            if len(g) < 4:
                continue
            J = np.array([r["J"] for r in g])
            sp = np.array([r["spread"] for r in g])
            ok = np.isfinite(J) & np.isfinite(sp) & (J > 0)
            if ok.sum() < 4:
                continue
            rho = float(np.corrcoef(np.log10(J[ok]), sp[ok])[0, 1])
            lines.append(f"  {bc:>10} CV={cv:.2f}: spread {sp.min():.3f} to "
                         f"{sp.max():.3f} over J {J.min():.2g} to {J.max():.2g}; "
                         f"corr(log J, spread) = {rho:+.3f}")

    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(OUT, "phase5_verdict.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
