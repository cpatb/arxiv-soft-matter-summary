"""Item 13, step 1: REPRODUCE the diagnosis before changing anything.

``dst_fit.profiling.profile_likelihood`` walks outward from the best fit with

    step = max(wald / 4, span / 5000),   span = |box_bound - best|

and, once the statistic is past ``f_crit``, runs ``_REFINEMENTS = 10`` bisections on
the bracket ``[prev, prev + step]`` and returns the OUTER end. Two consequences:

1. The floor ``span / 5000`` is sized from the parameter's arbitrary declared BOX
   bound, not from anything about the data. For a viscosity bounded at 1e6 that is
   ``step = 200`` on a parameter whose value is 0.0038.
2. Ten bisections of a bracket of width ``step`` leave a resolution of
   ``step / 1024``. When the crossing happens on the FIRST step the returned
   endpoint is therefore ``best + m * span / 5.12e6`` for an integer ``m <= 1024``:
   a point on a bisection GRID, not a likelihood crossing.

This script measures how much of the registry sits on that grid. It reads the
committed sidecars and each model's own declared bounds; it re-fits nothing, so it
cannot be accused of measuring a new fit rather than the shipped one.

**Since item 13 regenerated every sidecar with the corrected profiler, running this
against today's registry is a REGRESSION CHECK and should report zero.** The historical
measurement -- the registry as it was before that regeneration -- is preserved as
``data/profile_interval_prefix_rows.json`` and re-graded by
``scripts/profile_interval_regrade.py``. Anything asking "how bad was it?" must read
that; this file answers "is it still happening?".

The grid signature is exact arithmetic, not a fit: ``(hi - value) * 5.12e6 / span``
is an integer to within ``GRID_TOL`` relative. ``span`` is taken from the model's
own bound for that parameter, so no candidate bound is guessed -- the earlier count
of 117/296 tested only three candidate bounds and was flagged as conservative.

Run:  python scripts/profile_interval_diagnosis.py
Writes: data/profile_interval_diagnosis.md and .json
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
import yaml

from biorheo.constitutive import get_model

REPO = Path(__file__).resolve().parent.parent
RECORDS = REPO / "src" / "biorheo" / "registry" / "records"
OUT_MD = REPO / "data" / "profile_interval_diagnosis.md"
OUT_JSON = REPO / "data" / "profile_interval_diagnosis.json"

# dst_fit.profiling: step floor is span/5000, then _REFINEMENTS = 10 bisections.
_STEPS_PER_SPAN = 5000.0
_REFINEMENTS = 10
GRID_DIVISOR = _STEPS_PER_SPAN * (2 ** _REFINEMENTS)  # 5.12e6

# An integer to one part in 1e-6 of a unit. The quantity is a difference of two
# float64s that were written through YAML at full repr precision, so the residual
# is rounding, not tolerance shopping: donor A's eta_c reproduces m = 1 to 10 digits.
GRID_TOL = 1e-6

# ...AND the integer must be one a bisection can actually produce. With
# `_REFINEMENTS = 10` the returned endpoint is `prev + m * step / 1024`, so `m` lies in
# `[1, 1024]`; `m > 1024` is not expressible by that arithmetic.
#
# Added 2026-09-06 after NEXT_STEPS item 16. Without it the test accepted any `delta/span`
# that happened to be a round fraction, and in particular accepted every end whose `delta`
# EQUALS its `span` -- `m = 5.12e6`, which means the endpoint sits AT the box bound and is
# not a crossing at all. 18 of the 219 previously reported points were of that kind, 5 of
# them at exactly 5120000. The bound is derived from `_REFINEMENTS`, not chosen: it was
# found by a gate (G16.1) failing on those 18 and it would be the right bound regardless.
GRID_M_MIN, GRID_M_MAX = 1, 2 ** _REFINEMENTS


def _sidecars() -> list[Path]:
    return sorted(RECORDS.glob("*.identifiability.yaml"))


def _grid_index(delta: float, span: float) -> float | None:
    """``m`` such that ``delta = m * span / 5.12e6``, or None if not usable."""
    if not np.isfinite(delta) or not np.isfinite(span) or span <= 0 or delta <= 0:
        return None
    return delta * GRID_DIVISOR / span


def _is_grid_point(m: float | None) -> bool:
    """An integer index a 10-fold bisection of one step could actually have produced."""
    if m is None or not np.isfinite(m):
        return False
    return (abs(m - round(m)) < GRID_TOL
            and GRID_M_MIN - GRID_TOL <= m <= GRID_M_MAX + GRID_TOL)


def main() -> int:
    rows = []
    for path in _sidecars():
        side = yaml.safe_load(path.read_text(encoding="utf-8"))
        model_name = side["model"]
        try:
            model = get_model(model_name)
        except Exception as exc:  # a model the registry has but the package renamed
            print(f"SKIP {path.name}: {exc}")
            continue
        lo_b = dict(zip(model.param_names, model.bounds[0]))
        hi_b = dict(zip(model.param_names, model.bounds[1]))
        for pname, p in side["parameters"].items():
            value = float(p["value"])
            for end, at_bound_key, bound in (
                ("hi", "hi_at_bound", hi_b[pname]),
                ("lo", "lo_at_bound", lo_b[pname]),
            ):
                if p.get(at_bound_key):
                    continue  # an open end is not a crossing and has no grid
                endpoint = float(p[end])
                delta = abs(endpoint - value)
                span = abs(float(bound) - value)
                m = _grid_index(delta, span)
                on_grid = _is_grid_point(m)
                # an index that is an integer but OUT OF RANGE: reported, not silently
                # dropped, because it was counted as a grid point until 2026-09-06
                out_of_range = bool(
                    m is not None and np.isfinite(m)
                    and abs(m - round(m)) < GRID_TOL and not on_grid)
                rows.append({
                    "record": path.name.replace(".identifiability.yaml", ""),
                    "model": model_name,
                    "param": pname,
                    "end": end,
                    "value": value,
                    "endpoint": endpoint,
                    "bound": float(bound),
                    "delta": delta,
                    "grid_index": m,
                    "on_grid": bool(on_grid),
                    "integer_but_out_of_range": out_of_range,
                    "delta_over_value": (delta / abs(value)) if value != 0 else float("inf"),
                })

    closed = [r for r in rows if r["grid_index"] is not None]
    on_grid = [r for r in closed if r["on_grid"]]
    upper = [r for r in closed if r["end"] == "hi"]
    upper_on_grid = [r for r in upper if r["on_grid"]]

    # The worst case is the one whose interval is most inflated RELATIVE to the
    # parameter it describes -- a ratio whose denominator is the best fit, which is
    # a declared scale (the shipped value) and not another uncertain quantity.
    worst = max(on_grid, key=lambda r: r["delta_over_value"]) if on_grid else None
    by_model = Counter((r["model"], r["param"]) for r in on_grid)

    n_out_of_range = sum(1 for r in rows if r["integer_but_out_of_range"])
    summary = {
        "n_sidecars": len(_sidecars()),
        "n_integer_but_out_of_range": n_out_of_range,
        "grid_m_range": [GRID_M_MIN, GRID_M_MAX],
        "n_closed_ends": len(closed),
        "n_on_grid": len(on_grid),
        "n_closed_upper": len(upper),
        "n_upper_on_grid": len(upper_on_grid),
        "grid_divisor": GRID_DIVISOR,
        "grid_tol": GRID_TOL,
        "worst": worst,
        "on_grid_by_model_param": {f"{m}.{p}": n for (m, p), n in by_model.most_common()},
    }
    OUT_JSON.write_text(json.dumps({"summary": summary, "rows": rows}, indent=2), encoding="utf-8")

    lines = [
        "# Is the registry still on the bisection grid? (regression check)",
        "",
        "Written by `scripts/profile_interval_diagnosis.py`. Reads the committed sidecars",
        "and each model's own declared bounds; refits nothing.",
        "",
        "**Item 13 regenerated every sidecar with the corrected profiler, so this should",
        "report ZERO.** The historical measurement -- how bad it was before that -- is",
        "`data/profile_interval_prefix.md`, re-graded from the row set committed at",
        "`d38bc77`. Do not quote this file for the size of the defect; it is the check that",
        "the defect is gone.",
        "",
        "## The signature",
        "",
        "`dst_fit.profiling.profile_likelihood` sizes its outward step as",
        "`max(wald/4, span/5000)` with `span = |box_bound - best|`, then bisects the",
        "bracket 10 times and returns the OUTER end. If the statistic crosses on the",
        "first step, the returned endpoint is `best + m * span / 5.12e6` for an integer",
        f"`m <= 1024`. Counted here as `|m - round(m)| < {GRID_TOL:g}` **and**",
        f"`{GRID_M_MIN} <= m <= {GRID_M_MAX}`.",
        "",
        "**The range half of that test was added on 2026-09-06** (NEXT_STEPS item 16).",
        "Without it the test accepted any `delta/span` that happened to be a round",
        "fraction -- including every end whose `delta` EQUALS its `span`, which is",
        f"`m = {GRID_DIVISOR:.3g}` and means the endpoint sits AT the box bound rather than",
        f"at a crossing. **{n_out_of_range}** ends are integer-indexed but out of range and",
        "are now excluded; 5 of them sit at exactly 5120000. The bound is derived from",
        "`_REFINEMENTS = 10`, not chosen to make anything pass.",
        "",
        "## What is on the grid",
        "",
        f"- sidecars read: **{len(_sidecars())}**",
        f"- closed interval ends (an end at a box bound is not a crossing and is excluded): **{len(closed)}**",
        f"- of those, on the grid: **{len(on_grid)}** ({100 * len(on_grid) / max(len(closed), 1):.1f} %)",
        f"- closed UPPER ends: **{len(upper)}**, on the grid: **{len(upper_on_grid)}**",
        "",
    ]
    if worst is not None:
        lines += [
            "## Worst case, by how far the endpoint sits above the value it describes",
            "",
            f"`{worst['record']}` · `{worst['model']}` · `{worst['param']}` ({worst['end']})",
            "",
            f"- best fit: `{worst['value']:.6g}`",
            f"- endpoint: `{worst['endpoint']:.6g}`",
            f"- box bound: `{worst['bound']:.6g}`",
            f"- endpoint - value: `{worst['delta']:.12g}` = "
            f"`{round(worst['grid_index'])} x {worst['bound']:g}/5.12e6` "
            f"(grid index {worst['grid_index']:.10f})",
            f"- that is **{worst['delta_over_value']:.4g}x the best-fit value** above it",
            "",
        ]
    lines += ["## Which (model, parameter) slots are affected", ""]
    for key, n in by_model.most_common():
        lines.append(f"- `{key}`: {n}")
    lines += [
        "",
        "## What this does NOT show",
        "",
        "It does not show the interval is too WIDE by a specific amount, and it is not a",
        "measurement of where the F statistic actually crosses. The grid arithmetic says",
        "only that the returned endpoint is a bisection endpoint. Whether each one is far",
        "from the crossing is the direct check, and that needs the profile re-evaluated at",
        "the returned endpoint -- gate G13.2.",
        "",
    ]
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")

    print(f"closed ends {len(closed)}, on grid {len(on_grid)} "
          f"({100 * len(on_grid) / max(len(closed), 1):.1f} %)")
    print(f"closed upper {len(upper)}, on grid {len(upper_on_grid)}")
    if worst is not None:
        print(f"worst: {worst['record']} {worst['param']} {worst['end']} "
              f"delta/value = {worst['delta_over_value']:.4g}, m = {worst['grid_index']:.10f}")
    print(f"wrote {OUT_MD.name} and {OUT_JSON.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
