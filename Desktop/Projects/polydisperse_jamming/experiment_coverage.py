"""experiment_coverage.py — can the EXISTING drop-impact data test any of the predictions?

WHY
---
`MECHANISM_2026-09-08.md` §6 proposes five experiments. Before any of them is designed it
is worth asking which axes the lab has ALREADY swept, and whether the ranges already on
disk are wide enough for the predictions to be visible. That is a cheaper question than
any of the experiments and it changes what to build.

This reads the per-drop table **read-only**. Nothing on `E:\\` is written or modified.

⚠ TWO TRAPS, BOTH ON RECORD IN THIS PROJECT, AND BOTH AVOIDED HERE:

  1. The path. `E:\\CPB\\Drop Impact Analyses...` ends in a Unicode private-use character
     (U+F028) that macOS wrote and that cannot be typed. A typed path fails `Test-Path`
     and the parent lists as empty, which reads exactly like deleted data. This module
     GLOBS for the directory and uses the returned name.

  2. The `arrested` column. The project records that this flag "stayed True for 22 drops
     that never jammed after the classification moved". So **this module uses only
     CONTROL variables** -- phi, u0, salt, d0 -- which are set by the experimenter and
     cannot go stale. No derived flag is read, and therefore no conclusion here depends
     on a classification whose definition is uncertain.

    python experiment_coverage.py
"""

from __future__ import annotations

import csv
import glob
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# the shift E-2 must resolve, MEASURED (out/phase8_verdict.txt, estimator E2)
DPHI_C_MEASURED = 0.0623
LADDER_RATIO_MIN = 1.5   # below this a "ladder" is nozzle scatter, not a designed axis


def find_table() -> str | None:
    """Glob for the per-drop table; never type the path (see the docstring)."""
    for pat in (r"E:\CPB\Drop Impact Analyses*\figures_rebuild_*\phase_drops.csv",
                r"E:\CPB\Drop Impact Analyses*\*\phase_drops.csv"):
        hits = sorted(glob.glob(pat))
        if hits:
            return hits[-1]
    return None


def col(rows: list[dict], key: str) -> np.ndarray:
    out = []
    for r in rows:
        try:
            out.append(float(r[key]))
        except (TypeError, ValueError, KeyError):
            out.append(float("nan"))
    return np.array(out, float)


def main() -> int:
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    path = find_table()
    if path is None:
        p("per-drop table not found. ⚠ That is NOT evidence the data is missing: the")
        p("containing directory's name ends in a Unicode private-use character and a")
        p("typed path cannot reach it. Check the drive is mounted, then glob again.")
        return 1

    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    p(f"=== COVERAGE OF THE EXISTING PER-DROP TABLE ({len(rows)} drops) ===")
    p(f"  {os.path.basename(os.path.dirname(path))}/phase_drops.csv, read-only")
    p("  CONTROL variables only -- no derived flag is read, deliberately.")
    p()
    p(f"  {'axis':12s} {'n':>5s} {'min':>9s} {'max':>9s} {'span':>9s} "
      f"{'max/min':>8s} {'levels':>7s}")
    axes = {}
    for name, key in (("phi", "phi"), ("u0 [m/s]", "u0"), ("salt [mM]", "salt"),
                      ("d0 [mm]", "d0_mm")):
        v = col(rows, key)
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue
        ratio = v.max() / v.min() if v.min() > 0 else float("inf")
        axes[name] = v
        p(f"  {name:12s} {v.size:5d} {v.min():9.3f} {v.max():9.3f} "
          f"{v.max() - v.min():9.3f} {ratio:8.2f} "
          f"{len(set(np.round(v, 3))):7d}")

    p()
    p("=== WHAT THIS DATA CAN AND CANNOT TEST ===")

    phi = axes.get("phi")
    if phi is not None:
        span = float(phi.max() - phi.min())
        p("  E-2 (raise phi until the solid branch returns) must resolve a shift of")
        p(f"      {DPHI_C_MEASURED:.4f} in phi. The table's ENTIRE phi span is "
          f"{span:.3f}.")
        if span < DPHI_C_MEASURED:
            p(f"      🔴 {DPHI_C_MEASURED / span:.1f}x TOO SMALL -- the existing ladder "
              f"cannot test E-2")
            p(f"      at all. It needs roughly doubling, to about "
              f"{phi.min():.2f}-{phi.min() + 2.2 * DPHI_C_MEASURED:.2f}.")
        else:
            p("      the span is sufficient.")
        p(f"      (and only {len(set(np.round(phi, 3)))} distinct levels, so the "
          f"threshold in phi")
        p("      could not be LOCATED even if the span were wide enough.)")

    d0 = axes.get("d0 [mm]")
    if d0 is not None:
        ratio = float(d0.max() / d0.min())
        p()
        p("  E-6 (drop-size ladder) needs the three drainage readings to separate,")
        p(f"      which they do over a factor ~2-3 in D0. Here D0 ranges "
          f"{d0.min():.2f}-{d0.max():.2f} mm,")
        p(f"      a ratio of {ratio:.2f}.")
        if ratio < LADDER_RATIO_MIN:
            p("      🔴 That is NOZZLE SCATTER, not a designed axis: the values are")
            p("      near-continuous around a single nominal size. E-6 needs deliberate")
            p("      nozzles, not a re-analysis of what exists.")

    salt = axes.get("salt [mM]")
    if salt is not None:
        lv = sorted(set(np.round(salt, 1)))
        p()
        p(f"  THE SALT AXIS IS A REAL LADDER: {len(lv)} levels, {lv}.")
        p("      🔑 This is the one designed axis the lab already has, and the mechanism")
        p("      makes predictions on it -- the double-layer estimate in")
        p("      phi_axis_offset.py has phi_c_eff moving from 0.23 to 0.56 across")
        p("      0-100 mM. So the CHEAPEST test of the sign mechanism may be a")
        p("      re-analysis of the existing salt series, not a new experiment.")
        p("      ⚠ Testing that needs a regime classification, and the project records")
        p("      that the `arrested` flag went stale once. Establish which criterion is")
        p("      current BEFORE using it -- this module deliberately does not.")

    out = os.path.join(HERE, "out", "experiment_coverage.txt")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
