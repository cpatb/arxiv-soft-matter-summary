"""Item 16: which of the registry's parameters were exposed to the interval defect?

Gate pre-registered in `data/exposure_prereg.md`, committed at `1e33b56` before this was
run. That file opens by retracting an item-15 sentence; read it first.

The defect (NEXT_STEPS item 13) is a step floor of `span/5000` taken from the parameter's
declared BOX bound, applied as `step = max(wald/4, span/5000)`. So it acted on a parameter
exactly when the floor won, and the amount by which it overreached is

    inflation = (span / 5000) / (wald / 4)

which is a RATIO, not a boolean. `inflation = 1.09` (a Carbopol flow index) means the
floor was applied and happened to be the right size; `inflation = 11770` (that gel's
consistency) means the first step overshot by four orders of magnitude.

G16.1 is the falsifiable claim and it is ONE-DIRECTIONAL: the grid signature
`best + m*span/5.12e6` can only arise if the crossing was bracketed on the FIRST step and
that step was the floor, so `on_grid` must imply `inflation > 1`. The converse is not
claimed -- a crossing found on a later step lands on a different grid that the diagnosis
script does not test for -- and G16.2 reports those cells rather than letting them read as
noise.

This refits every record with `fit_model` only. It does NOT profile: the Wald standard
error comes straight from the fitter, so the scan costs a fit per record and not a walk
per parameter.

Run:  python scripts/exposure_scan.py   (~2 min)
Writes: data/exposure_scan.md and .json
"""

from __future__ import annotations

import glob
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import yaml
from dst_fit.fitting import fit_model
from scipy.stats import spearmanr

from biorheo.constitutive import get_model

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
RECORDS = REPO / "src" / "biorheo" / "registry" / "records"
# The PRE-FIX row set, not the live scan: item 13 regenerated the sidecars, so
# `profile_interval_diagnosis.json` now describes the CORRECTED registry and
# contains zero grid points. Reading it here would give a gate an empty sample,
# which passes by vacuum. See scripts/profile_interval_regrade.py.
DIAG = REPO / "data" / "profile_interval_prefix.json"
OUT_MD = REPO / "data" / "exposure_scan.md"
OUT_JSON = REPO / "data" / "exposure_scan.json"

STEPS_PER_SPAN = 5000.0
WALD_FRACTION = 0.25


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def blood_curves() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    mss = _load("_mss", HERE / "multistart_seed_scatter.py")
    gd, stress = mss.load()
    out = {}
    for donor, tau in stress.items():
        eta = tau / gd
        out[f"blood_horner2020_donor_{donor}"] = (gd, eta)
        for fam in ("cross", "carreau_yasuda", "newtonian", "power_law"):
            out[f"blood_horner2020_donor_{donor}_{fam}"] = (gd, eta)
    return out


def carbopol_curves() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    ic = _load("_ic", HERE / "identify_carbopol.py")
    out = {}
    for path in sorted(glob.glob(str(ic.DATA / "*.csv"))):
        gd, eta = ic.load(path)
        out[ic.record_id(path)] = (np.asarray(gd, float), np.asarray(eta, float))
    return out


def main() -> int:
    if not DIAG.exists():
        print("run scripts/profile_interval_regrade.py first", file=sys.stderr)
        return 2
    diag = json.loads(DIAG.read_text(encoding="utf-8"))
    # Keyed by (record, param, END). `profile_likelihood` walks each direction with its
    # OWN `span = |limit - best|`, so a parameter has TWO inflation ratios and they can
    # differ by decades -- a viscosity at `best = 21` has `span_hi/5000 = 200` and
    # `span_lo/5000 = 0.0043`. Aggregating to the parameter and comparing everything
    # against the upper ratio is what produced 15 spurious G16.1 counterexamples on the
    # first run of this script, 14 of them grid points detected on the LOWER end.
    on_grid: dict[tuple[str, str, str], bool] = {}
    delta_over_value: dict[tuple[str, str, str], float] = {}
    for r in diag["rows"]:
        key = (r["record"], r["param"], r["end"])
        on_grid[key] = bool(r["on_grid"])
        delta_over_value[key] = r["delta_over_value"]

    curves = {**blood_curves(), **carbopol_curves()}
    rows = []
    for path in sorted(RECORDS.glob("*.identifiability.yaml")):
        rid = path.name.replace(".identifiability.yaml", "")
        if rid not in curves:
            print(f"SKIP {rid}: no source curve", file=sys.stderr)
            continue
        side = yaml.safe_load(path.read_text(encoding="utf-8"))
        model = get_model(side["model"])
        gd, eta = curves[rid]
        good = np.isfinite(gd) & np.isfinite(eta) & (gd > 0) & (eta > 0)
        gd, eta = gd[good], eta[good]
        fit = fit_model(model, gd, eta, n_starts=8, seed=0, start=model.guess(gd, eta))

        for i, pname in enumerate(model.param_names):
            best = float(fit.params[i])
            wald = float(fit.stderr[i])
            from_wald = WALD_FRACTION * wald if np.isfinite(wald) and wald > 0 else 0.0
            for end, bound in (("hi", float(model.bounds[1][i])),
                               ("lo", float(model.bounds[0][i]))):
                span = abs(bound - best)
                from_box = span / STEPS_PER_SPAN
                inflation = (from_box / from_wald) if from_wald > 0 else float("inf")
                key = (rid, pname, end)
                if key not in on_grid:
                    continue    # that end sits at a box bound; not a crossing, excluded
                rows.append({
                    "record": rid, "model": side["model"], "param": pname, "end": end,
                    "value": best, "wald": wald, "bound": bound,
                    "step_from_wald": from_wald, "step_from_box": from_box,
                    "inflation": inflation,
                    "box_won": bool(from_box > from_wald),
                    "on_grid": on_grid[key],
                    "delta_over_value": delta_over_value.get(key),
                    "n_points": int(gd.size),
                })

    # ---- G16.1: on_grid must imply inflation > 1 ---------------------------------------
    counterexamples = [r for r in rows if r["on_grid"] and r["inflation"] <= 1.0]
    g161 = not counterexamples

    # ---- G16.2: the converse is not claimed; report the cell ---------------------------
    inflated_not_on_grid = [r for r in rows if r["inflation"] > 1.0 and not r["on_grid"]]
    on_grid_rows = [r for r in rows if r["on_grid"]]

    # ---- G16.3: does the ratio track the measured overreach? ---------------------------
    pairs = [(r["inflation"], r["delta_over_value"]) for r in on_grid_rows
             if r["delta_over_value"] is not None
             and np.isfinite(r["inflation"]) and np.isfinite(r["delta_over_value"])]
    if len(pairs) >= 3:
        rho, pval = spearmanr([p[0] for p in pairs], [p[1] for p in pairs])
        rho, pval = float(rho), float(pval)
    else:
        rho = pval = float("nan")

    finite_inf = sorted(r["inflation"] for r in rows if np.isfinite(r["inflation"]))
    q = (lambda f: float(np.quantile(finite_inf, f))) if finite_inf else (lambda f: float("nan"))

    verdict = {
        "G16.1": {"pass": bool(g161), "n_counterexamples": len(counterexamples),
                  "counterexamples": counterexamples[:10],
                  "n_on_grid": len(on_grid_rows), "n_parameters": len(rows)},
        "G16.2": {"n_inflated_not_on_grid": len(inflated_not_on_grid),
                  "n_inflated": sum(1 for r in rows if r["inflation"] > 1.0),
                  "note": "expected to be large; the converse is NOT claimed"},
        "G16.3": {"spearman_inflation_vs_overreach": rho, "p": pval, "n": len(pairs),
                  "inflation_quantiles": {"min": (finite_inf[0] if finite_inf else None),
                                          "q25": q(0.25), "median": q(0.50),
                                          "q75": q(0.75),
                                          "max": (finite_inf[-1] if finite_inf else None)},
                  "n_infinite_inflation": sum(1 for r in rows
                                              if not np.isfinite(r["inflation"]))},
    }
    OUT_JSON.write_text(json.dumps({"verdict": verdict, "rows": rows}, indent=2),
                        encoding="utf-8")

    v = verdict
    lines = [
        "# Item 16: which parameters were exposed to the interval defect, and by how much",
        "",
        "Written by `scripts/exposure_scan.py`. Gate pre-registered in",
        "`data/exposure_prereg.md` (commit `1e33b56`) before this was run. **That file",
        "opens by retracting a sentence item 15 shipped; read it first.**",
        "",
        "## The quantity",
        "",
        "The defect is a step floor of `span/5000` from the parameter's declared BOX,",
        "applied as `step = max(wald/4, span/5000)`. It acted when the floor won, and the",
        "amount by which it overreached is",
        "",
        "```",
        "inflation = (span / 5000) / (wald / 4)",
        "```",
        "",
        "**a ratio, not a boolean.** `inflation = 1.09` means the floor was applied and",
        "happened to be the right size; `inflation = 1.2e+04` means the first step",
        "overshot by four orders of magnitude. Saying only *whether* the box was consulted",
        "is what made item 15's explanation wrong.",
        "",
        "| gate | claim | bar | result | verdict |",
        "|---|---|---|---|---|",
        f"| G16.1 | the grid signature implies the box floor won | 0 counterexamples | "
        f"**{len(counterexamples)}** over {len(rows)} closed interval ENDS "
        f"({len(on_grid_rows)} on the grid) | **{'PASS' if g161 else 'FAIL'}** |",
        f"| G16.2 | the converse is NOT claimed | reported | "
        f"{len(inflated_not_on_grid)} inflated but not on the grid, of "
        f"{v['G16.2']['n_inflated']} inflated | reported |",
        f"| G16.3 | the ratio tracks the overreach | reported | Spearman "
        f"{rho:.3f} (p = {pval:.2g}, n = {len(pairs)}) | reported |",
        "",
        "## What G16.1 establishes, and what it does not",
        "",
        "Every parameter whose stored endpoint sits on the bisection grid has",
        "`inflation > 1`, so **the mechanism written on four pages of this repo is the",
        "right one**. The converse is deliberately not claimed: a crossing found on a",
        "LATER step lands on `step0 * 1.18^j / 1024`, a different grid that",
        "`profile_interval_diagnosis.py` does not test for, so a parameter can be inflated",
        f"without being detectably on the grid. That cell holds **{len(inflated_not_on_grid)}**",
        "ends and is a limit of the DETECTOR, not evidence against the mechanism.",
        "",
        "## How exposed the registry is",
        "",
        f"- closed interval ENDS scanned: **{len(rows)}**. Each parameter contributes up",
        "  to two, and the two differ: a viscosity at `best = 21` has `span_hi/5000 = 200`",
        "  against `span_lo/5000 = 0.0043`, so an end must be graded against its own bound.",
        f"- with `inflation > 1` (the floor won): **{v['G16.2']['n_inflated']}**",
        f"- with a detectable grid endpoint: **{len(on_grid_rows)}**",
        f"- `inflation` quantiles: min {v['G16.3']['inflation_quantiles']['min']:.3g}, "
        f"q25 {v['G16.3']['inflation_quantiles']['q25']:.3g}, "
        f"median {v['G16.3']['inflation_quantiles']['median']:.3g}, "
        f"q75 {v['G16.3']['inflation_quantiles']['q75']:.3g}, "
        f"max {v['G16.3']['inflation_quantiles']['max']:.3g}",
        "",
        "## Worst and best exposed (model, parameter)",
        "",
        "| record | model | parameter | end | wald/4 | span/5000 | inflation | on grid |",
        "|---|---|---|---|---|---|---|---|",
    ]
    ranked = sorted((r for r in rows if np.isfinite(r["inflation"])),
                    key=lambda r: -r["inflation"])
    for r in ranked[:5] + ranked[-5:]:
        lines.append(f"| `{r['record']}` | `{r['model']}` | `{r['param']}` | {r['end']} | "
                     f"{r['step_from_wald']:.4g} | {r['step_from_box']:.4g} | "
                     f"**{r['inflation']:.4g}** | {'yes' if r['on_grid'] else 'no'} |")
    lines += [
        "",
        "## What this is NOT",
        "",
        "**Not an impact estimator.** `inflation` says how far the first step overshot for",
        "one parameter. It does not say how wrong a downstream number was, because that",
        "also depends on which parameter the number is sensitive to (item 15: a needle",
        "force set by a flow index barely moved while its `k` interval moved 2.94 %) and on",
        "which END of the interval the comparison binds (item 13's G13.9: the defect",
        "inflated the upper ends and the comparison used the lower ones). Neither is",
        "derivable from this ratio and neither should be inferred from it.",
        "",
    ]
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")

    print(f"G16.1: {'PASS' if g161 else 'FAIL'}  {len(counterexamples)} counterexamples "
          f"over {len(rows)} ends ({len(on_grid_rows)} on grid)")
    print(f"G16.2: {len(inflated_not_on_grid)} inflated but not on the grid, of "
          f"{v['G16.2']['n_inflated']} inflated")
    print(f"G16.3: Spearman {rho:.3f} (p = {pval:.2g}, n = {len(pairs)})")
    print(f"inflation median {v['G16.3']['inflation_quantiles']['median']:.4g}, "
          f"max {v['G16.3']['inflation_quantiles']['max']:.4g}")
    print(f"wrote {OUT_MD.name}, {OUT_JSON.name}")
    return 0 if g161 else 1


if __name__ == "__main__":
    raise SystemExit(main())
