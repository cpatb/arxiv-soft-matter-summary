"""analyse_phase6.py — how much tensile strength does a rebound need?

Grades the cohesion ladder against two pre-registered checks:

  DISCRIMINATING NEGATIVE   at C = 0 the cohesive machinery must reproduce the
                            cohesionless result (|COR| <= 0.05). If a bonded-contact
                            model reports a rebound with the bonds switched off, it is
                            broken and nothing else in this file means anything.

  POSITIVE CONTROL          somewhere on the ladder COR must rise clearly above the
                            floor. If it never does, the control has FAILED: either the
                            ladder was too short or the cohesion is not being applied,
                            and "nothing rebounds" cannot then be attributed to the
                            absence of tensile strength.

The number this exists to produce is the threshold: the ratio C = sigma_coh / (rho u0^2)
at which an aggregate of this kind starts to come back.
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

COR_FLOOR = 0.05        # the cohesionless result, |COR| <= 0.035, plus margin
COR_REBOUND = 0.15      # "clearly rebounding", well above the floor and well below
                        # the 0.819 an elastic sphere gives through this damping

# A rebound needs TWO things, and this file originally tested only the first. Phases 5
# and 7 both require recoil AND coherence, because a shattered aggregate whose fragments
# recoil produces a centre-of-mass COR that is not a restitution coefficient. Applying
# the same two-part criterion here is a consistency fix, not a new bar: with six seeds
# the CV = 0.5 / C = 3000 cell reaches COR +0.171, which clears the recoil bar, while its
# spread is 2.649 -- the aggregate is in pieces.
SPREAD_COHERENT = 1.5

# And a mean smaller than twice its own seed scatter is not a measurement. At CV = 0.5
# the COR sd is comparable to the mean at every C above 100.
RESOLVE_FACTOR = 2.0


def load():
    by_key, n_lines, bad = {}, 0, 0
    for f in sorted(glob.glob(os.path.join(OUT, "phase6*.jsonl"))):
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
                    bad += 1
                    continue
                if r.get("key"):
                    by_key[r["key"]] = r
    print(f"N: {n_lines} lines -> {len(by_key)} unique cells ({bad} failed)")
    return list(by_key.values())


def main() -> int:
    recs = load()
    if not recs:
        print("no phase6 records yet")
        return 1

    groups: dict[tuple, list[dict]] = {}
    for r in recs:
        groups.setdefault((r["cv"], r["coh_ratio"]), []).append(r)

    rows = []
    for (cv, C), g in groups.items():
        row = {"cv": cv, "C": C, "n": len(g)}
        for k in ("cor", "spread", "ke_ratio", "z_initial", "z_peak", "n_grains"):
            v = np.array([x.get(k, np.nan) for x in g], float)
            v = v[np.isfinite(v)]
            row[k] = float(v.mean()) if v.size else float("nan")
            row[k + "_sd"] = float(v.std(ddof=1)) if v.size > 1 else float("nan")
        rows.append(row)
    rows.sort(key=lambda r: (r["cv"], r["C"]))

    fp = os.path.join(OUT, "phase6_summary.csv")
    with open(fp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {fp}")

    lines = ["=" * 80,
             "PHASE 6 -- the cohesion ladder",
             "=" * 80, "",
             "C = sigma_coh / (rho u0^2), the bond strength in units of the impact stress.",
             "Reference points: the cohesionless aggregate gives |COR| <= 0.035, and one",
             "perfectly elastic sphere through this same damping gives COR = 0.819.", ""]

    for cv in sorted({r["cv"] for r in rows}):
        lines.append(f"--- CV = {cv:g} ---")
        lines.append(f"  {'C':>7} {'n':>3} {'COR':>9} {'sd':>8} {'spread':>8} "
                     f"{'KE kept':>9} {'Z init':>7} {'Z peak':>7}")
        for r in [x for x in rows if x["cv"] == cv]:
            lines.append(f"  {r['C']:7.4g} {r['n']:3d} {r['cor']:+9.4f} "
                         f"{r['cor_sd']:8.4f} {r['spread']:8.3f} {r['ke_ratio']:9.4f} "
                         f"{r['z_initial']:7.2f} {r['z_peak']:7.2f}")
        lines.append("")

    lines.append("=" * 80)
    lines.append("CHECKS")
    lines.append("=" * 80)
    ok = True

    zero = [r for r in rows if r["C"] == 0.0]
    if zero:
        worst = max(abs(r["cor"]) for r in zero)
        passed = worst <= COR_FLOOR
        ok &= passed
        lines.append(f"  DISCRIMINATING NEGATIVE (C = 0 must not rebound, "
                     f"|COR| <= {COR_FLOOR}):")
        lines.append(f"    worst |COR| at C = 0 is {worst:.4f} -> "
                     f"{'PASS' if passed else 'FAIL -- the bonded model rebounds with '
                        'its bonds off, so nothing here is trustworthy'}")
    else:
        lines.append("  DISCRIMINATING NEGATIVE: no C = 0 cells present")
        ok = False

    def rebounds(r):
        return (r["cor"] >= COR_REBOUND
                and r["spread"] < SPREAD_COHERENT
                and (not np.isfinite(r["cor_sd"])
                     or r["cor"] > RESOLVE_FACTOR * r["cor_sd"]))

    reb = [r for r in rows if rebounds(r)]
    recoil_only = [r for r in rows if r["cor"] >= COR_REBOUND and not rebounds(r)]
    lines.append("")
    lines.append(f"  POSITIVE CONTROL (some C must give COR >= {COR_REBOUND}, "
                 f"spread < {SPREAD_COHERENT}, and mean > {RESOLVE_FACTOR}x its sd):")
    if recoil_only:
        lines.append(f"    {len(recoil_only)} cell group(s) clear the RECOIL bar but fail "
                     "coherence or resolution -- not rebounds:")
        for r in sorted(recoil_only, key=lambda x: -x["cor"]):
            why = []
            if r["spread"] >= SPREAD_COHERENT:
                why.append(f"spread {r['spread']:.2f}")
            if np.isfinite(r["cor_sd"]) and r["cor"] <= RESOLVE_FACTOR * r["cor_sd"]:
                why.append(f"mean {r['cor']:.3f} vs sd {r['cor_sd']:.3f}")
            lines.append(f"      CV={r['cv']:.2f} C={r['C']:g} COR={r['cor']:+.4f} "
                         f"({'; '.join(why)})")
    if reb:
        lo = min(r["C"] for r in reb)
        lines.append(f"    PASS -- rebound appears at C = {lo:g} and above")
        for cv in sorted({r["cv"] for r in rows}):
            g = sorted([r for r in rows if r["cv"] == cv], key=lambda x: x["C"])
            cross = None
            for a, b in zip(g[:-1], g[1:]):
                if not rebounds(a) and rebounds(b):
                    cross = (a["C"], b["C"])
                    break
            lines.append(f"    CV = {cv:g}: threshold bracketed by C = "
                         f"{cross[0]:g} to {cross[1]:g}" if cross else
                         f"    CV = {cv:g}: no crossing inside the ladder")
        lines.append("")
        lines.append("  => The absence of rebound in the cohesionless runs is attributable")
        lines.append("     to the absence of TENSILE STRENGTH, not to the machinery. A")
        lines.append("     jammed drop that bounces needs a bond strength of at least")
        lines.append(f"     ~{lo:g}x the impact stress rho*U0^2, supplied by something this")
        lines.append("     model does not contain: capillary bridges, adhesion, or a")
        lines.append("     stress-activated network with genuine tensile capacity.")
    else:
        best = max(rows, key=lambda r: r["cor"])
        ok = False
        lines.append(f"    FAIL -- the largest COR on the whole ladder is "
                     f"{best['cor']:+.4f} at C = {best['C']:g}.")
        lines.append("    The ladder does not bracket a rebound, so it is EITHER too")
        lines.append("    short OR the cohesion is not reaching the contacts. Until this")
        lines.append("    is resolved, 'nothing rebounds' may NOT be attributed to the")
        lines.append("    absence of cohesion -- it is only established as an observation.")
        lines.append("    Next step: check normalAdhesion on a live interaction, and")
        lines.append("    extend the ladder.")

    lines.append("")
    lines.append("  PHASE6_PASS" if ok else "  PHASE6_FAIL")
    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(OUT, "phase6_verdict.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
