"""Re-grade the PRE-FIX interval measurement, which the working tree can no longer produce.

Why this script exists
----------------------
`profile_interval_diagnosis.py` reads the committed sidecars. Item 13 then REGENERATED
those sidecars with the corrected profiler, in the same commit that first wrote the
diagnosis artifact -- so from that commit onward the artifact described a state the tree
no longer contained, and re-running the script against today's registry finds **zero**
grid points. That is the right answer for today's registry and it is not the measurement
the diagnosis was making.

The raw ROWS are the measurement and they are preserved verbatim in
`data/profile_interval_prefix_rows.json`, taken from commit `d38bc77`. This script
re-grades those rows and nothing else. It reads no sidecars, refits nothing, and cannot
be affected by any later change to the registry.

What the re-grading changes
---------------------------
The original test asked only whether `m = delta * 5.12e6 / span` is an integer. With
`_REFINEMENTS = 10` the bisection can only return `m` in `[1, 1024]`, so an integer
outside that range is not a grid point -- and in particular `m = 5.12e6` means
`delta == span`, i.e. the endpoint sits AT the box bound. Adding the range test was
forced by NEXT_STEPS item 16, whose G16.1 failed on exactly those rows.

Run:  python scripts/profile_interval_regrade.py
Writes: data/profile_interval_prefix.md and .json
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
ROWS_IN = REPO / "data" / "profile_interval_prefix_rows.json"
OUT_MD = REPO / "data" / "profile_interval_prefix.md"
OUT_JSON = REPO / "data" / "profile_interval_prefix.json"

_STEPS_PER_SPAN = 5000.0
_REFINEMENTS = 10
GRID_DIVISOR = _STEPS_PER_SPAN * (2 ** _REFINEMENTS)   # 5.12e6
GRID_TOL = 1e-6
GRID_M_MIN, GRID_M_MAX = 1, 2 ** _REFINEMENTS          # a bisection cannot exceed this

PROVENANCE_COMMIT = "d38bc77"


def is_grid_point(m: float | None) -> bool:
    if m is None or not np.isfinite(m):
        return False
    return (abs(m - round(m)) < GRID_TOL
            and GRID_M_MIN - GRID_TOL <= m <= GRID_M_MAX + GRID_TOL)


def is_integer_index(m: float | None) -> bool:
    if m is None or not np.isfinite(m):
        return False
    return abs(m - round(m)) < GRID_TOL


def main() -> int:
    raw = json.loads(ROWS_IN.read_text(encoding="utf-8"))
    rows = raw["rows"]
    for r in rows:
        m = r["grid_index"]
        r["on_grid"] = is_grid_point(m)
        r["integer_but_out_of_range"] = is_integer_index(m) and not r["on_grid"]

    closed = [r for r in rows if r["grid_index"] is not None]
    on_grid = [r for r in closed if r["on_grid"]]
    upper = [r for r in closed if r["end"] == "hi"]
    upper_on_grid = [r for r in upper if r["on_grid"]]
    rejected = [r for r in rows if r["integer_but_out_of_range"]]
    worst = max(on_grid, key=lambda r: r["delta_over_value"]) if on_grid else None
    by_slot = Counter((r["model"], r["param"]) for r in on_grid)
    at_bound_like = [r for r in rejected
                     if abs(r["grid_index"] - GRID_DIVISOR) < 1.0]

    summary = {
        "provenance": f"rows as committed at {PROVENANCE_COMMIT}, before the sidecars "
                      "were regenerated with the corrected profiler",
        "n_closed_ends": len(closed),
        "n_on_grid": len(on_grid),
        "n_closed_upper": len(upper),
        "n_upper_on_grid": len(upper_on_grid),
        "n_rejected_out_of_range": len(rejected),
        "n_rejected_at_the_bound": len(at_bound_like),
        "grid_m_range": [GRID_M_MIN, GRID_M_MAX],
        "as_shipped_before_regrade": raw["summary"]["n_on_grid"],
        "as_shipped_upper_before_regrade": raw["summary"]["n_upper_on_grid"],
        "worst": worst,
        "on_grid_by_model_param": {f"{m}.{p}": n for (m, p), n in by_slot.most_common()},
    }
    OUT_JSON.write_text(json.dumps({"summary": summary, "rows": rows}, indent=2),
                        encoding="utf-8")

    lines = [
        "# The pre-fix interval measurement, re-graded",
        "",
        "Written by `scripts/profile_interval_regrade.py` from",
        "`data/profile_interval_prefix_rows.json`, which is the row set committed at",
        f"`{PROVENANCE_COMMIT}`. **The working tree can no longer produce this**: item 13",
        "regenerated every sidecar with the corrected profiler in that same commit, so",
        "`profile_interval_diagnosis.py` run today reads the CORRECTED registry and finds",
        "zero grid points. That is the right answer for today's registry and it is a",
        "different question from the one this file answers.",
        "",
        "## What the re-grade changed",
        "",
        "The original test asked only whether `m = delta * 5.12e6 / span` is an integer.",
        f"With `_REFINEMENTS = {_REFINEMENTS}` a bisection can only return `m` in",
        f"`[{GRID_M_MIN}, {GRID_M_MAX}]`, so an integer outside that range is not a grid",
        "point. The range test was forced by NEXT_STEPS item 16, whose G16.1 failed on",
        "exactly these rows.",
        "",
        f"- rejected as out of range: **{len(rejected)}**",
        "- of those, sitting at `m = 5.12e6`, i.e. `delta == span` and the endpoint AT the",
        f"  box bound rather than at any crossing: **{len(at_bound_like)}**",
        "",
        "## The corrected counts",
        "",
        "| quantity | as first reported | re-graded |",
        "|---|---|---|",
        f"| grid points, all ends | {summary['as_shipped_before_regrade']} | "
        f"**{len(on_grid)}** |",
        f"| closed UPPER ends on the grid | {summary['as_shipped_upper_before_regrade']} of "
        f"{len(upper)} | **{len(upper_on_grid)} of {len(upper)}** |",
        "",
    ]
    if worst is not None:
        lines += [
            f"The worst case is **unchanged**: `{worst['record']}` `{worst['param']}`",
            f"({worst['end']}), whose endpoint sits **{worst['delta_over_value']:.4g}x the",
            f"best-fit value** above it, at a legitimate grid index `m = "
            f"{round(worst['grid_index'])}`.",
            "",
        ]
    lines += [
        "## What this does and does not change",
        "",
        "**It does not change item 13's conclusion or any of its gates.** G13.2 evaluated",
        "the profile statistic AT the shipped endpoints and found a median of 710 times",
        "`f_crit`; that measurement does not go through this detector at all. The",
        "correction moves a COUNT by about one per cent and removes a class of false",
        "positive that would otherwise have been quoted as evidence.",
        "",
        "It does change the number that appears in prose: **201 of 296**, not 203.",
        "",
        "## Which (model, parameter) slots are affected",
        "",
    ]
    for key, n in by_slot.most_common():
        lines.append(f"- `{key}`: {n}")
    lines.append("")
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")

    print(f"re-graded {len(rows)} rows from {PROVENANCE_COMMIT}")
    print(f"  grid points: {summary['as_shipped_before_regrade']} -> {len(on_grid)}")
    print(f"  closed upper on grid: {summary['as_shipped_upper_before_regrade']} -> "
          f"{len(upper_on_grid)} of {len(upper)}")
    print(f"  rejected out of range: {len(rejected)} "
          f"({len(at_bound_like)} of them at delta == span)")
    if worst is not None:
        print(f"  worst delta/value unchanged: {worst['delta_over_value']:.4g} "
              f"(m = {round(worst['grid_index'])})")
    print(f"wrote {OUT_MD.name}, {OUT_JSON.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
