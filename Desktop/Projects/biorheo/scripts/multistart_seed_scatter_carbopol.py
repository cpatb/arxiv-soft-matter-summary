"""The seed-scatter question, asked of the 15 Carbopol records (NEXT_STEPS item 10).

`scripts/multistart_seed_scatter.py` measured how far fitted parameters move when only the
optimiser's seed changes, on the 20 BLOOD curves. That left a gap which
`solver.openfoam.check_exportable`'s docstring states explicitly: the 15 Carbopol records
were never in it. They differ from the blood set in every way that could matter -- a
different material, a different instrument, 54 points per curve instead of 16, and, most
importantly, **the shipped parameter values come from a different fitter entirely**
(`carbopol/fit_hb.py`, `curve_fit` on stress residuals), not from `biorheo.identify`.

That last difference changes what is at risk, and the third criterion below is about it.
For a blood record, seed scatter moves the value the record ships. For a Carbopol record
it does not: the value is the csv's. What the seed moves is the SIDECAR -- the interval
and the verdict -- and `write_identifiability` refuses to write a sidecar whose interval
does not contain the shipped csv value. So if that containment holds only at seed 0, the
record's validity is an accident of the seed, and the loader's containment check is
resting on nothing.

PRE-REGISTERED (written before any of this was measured; do not weaken):

  G10.1 Same as G9.1 on independent data. No Carbopol record whose parameters are all
        graded `measurable` may have a seed spread above 1e-3.
  G10.2 Same as G9.3. No Carbopol record that `check_exportable` accepts by default --
        i.e. none without a `degenerate` or `not_identifiable` parameter -- may have a
        seed spread above 1e-3.
  G10.3 THE ONE THAT IS SPECIFIC TO THIS SET. For every seed in 0-4, the profiled
        interval for every parameter must contain the shipped csv value, on all 15
        records. This is what the loader's containment rule assumes and it has only ever
        been checked at seed 0.

Stated in advance so it cannot be rationalised afterwards: the blood set FAILED the
G9.2/G9.3 pair on one Herschel-Bulkley record whose `tau_y` was `bounded_only`. These
records are also Herschel-Bulkley. If G10.1 or G10.2 fails here the same way, that is a
second instance of the same defect and it lands on records that DO ship, which would be a
materially worse result than the blood one -- it must be reported as such and not
softened by pointing at the different fitter.

Run:  python scripts/multistart_seed_scatter_carbopol.py   (a few minutes). Writes
`data/carbopol_multistart_scatter.md` and `.json`; pinned in
`tests/test_multistart_scatter.py`.
"""

from __future__ import annotations

import glob
import io
import json
import os
import warnings
from pathlib import Path

import numpy as np

from biorheo.constitutive import get_model
from biorheo.identify import fit_rate_curve
from biorheo.registry import load_all

HERE = Path(__file__).resolve().parent
DATA = HERE.parents[1] / "carbopol" / "data" / "etd2020_biermann2024"
OUT = HERE.parent / "data"
HB = get_model("herschel_bulkley")

SEEDS = (0, 1, 2, 3, 4)
BAD_VERDICTS = ("degenerate", "not_identifiable")
G101_BAR = 1e-3
G102_BAR = 1e-3


def load(path: str) -> tuple[np.ndarray, np.ndarray]:
    """(recovered rate, mean stress), the same ten lines as `fit_hb.load`."""
    raw = io.open(path, encoding="latin-1").read().splitlines()
    rows = [r.split(";") for r in raw[2:] if r.strip()]

    def num(c: str) -> float:
        return float(c.replace(",", ".")) if c.strip() else np.nan

    a = np.array([[num(c) for c in r] for r in rows])
    tau = a[:, [2, 4, 6]]
    eta = a[:, [3, 5, 7]] * 1e-3
    return np.nanmean(tau / eta, axis=1), np.nanmean(tau, axis=1)


def record_id(path: str) -> str:
    wt, tea = os.path.basename(path).replace("ETD2020-", "").replace(".csv", "").split("-")
    return f"carbopol_etd2020_{int(wt[:2]) / 10:g}wt_{tea[:2]}uL".replace(".", "p")


def main() -> None:
    warnings.filterwarnings("ignore")
    records = load_all(allow_network=False)
    files = sorted(glob.glob(str(DATA / "*.csv")))
    assert len(files) == 15, files

    out: dict[str, dict] = {}
    for f in files:
        g, t = load(f)
        good = np.isfinite(g) & np.isfinite(t) & (g > 0) & (t > 0)
        g, t = g[good], t[good]
        rec = records[record_id(f)]
        shipped = {n: float(v) for n, v in zip(HB.param_names, rec.params_si)}

        fits = [fit_rate_curve(HB, g, t / g, seed=s,
                               temperature_C=rec.envelope.temperature_C,
                               geometry=rec.envelope.geometry) for s in SEEDS]
        vals = np.array([np.asarray(fi.values, dtype=float) for fi in fits])
        cost = np.array([fi.fit.rmse_log for fi in fits], dtype=float)
        mean = np.abs(vals.mean(axis=0))
        rel = np.where(mean > 0, (vals.max(axis=0) - vals.min(axis=0)) /
                       np.where(mean > 0, mean, 1.0), 0.0)

        # G10.3: does every seed's interval still contain the shipped csv value?
        contains = {}
        for si, fi in zip(SEEDS, fits):
            for n in HB.param_names:
                est = fi.params[n]
                lo, hi = float(est.lo), float(est.hi)
                ok = bool(lo <= shipped[n] <= hi)
                contains[f"seed{si}:{n}"] = {"lo": lo, "hi": hi,
                                             "shipped": shipped[n], "contains": ok}
        n_fail = sum(1 for v in contains.values() if not v["contains"])

        verdicts = {n: fits[0].params[n].verdict for n in HB.param_names}
        out[rec.id] = {
            "n_points": int(g.size),
            "param_spread": float(rel.max()),
            "worst_param": HB.param_names[int(np.argmax(rel))],
            "per_param_spread": {n: float(r) for n, r in zip(HB.param_names, rel)},
            "cost_spread": float((cost.max() - cost.min()) / cost.mean()),
            "verdicts": verdicts,
            "verdicts_by_seed": {str(s): {n: fi.params[n].verdict for n in HB.param_names}
                                 for s, fi in zip(SEEDS, fits)},
            "all_measurable": all(v == "measurable" for v in verdicts.values()),
            "has_bad_verdict": any(v in BAD_VERDICTS for v in verdicts.values()),
            "shipped_values": shipped,
            "n_interval_failures": n_fail,
            "containment": contains,
        }
        print(f"{rec.id:34s} N={g.size:3d} spread {rel.max():.3e} ({out[rec.id]['worst_param']}) "
              f"cost {out[rec.id]['cost_spread']:.2e}  "
              f"{''.join(v[0] for v in verdicts.values())}  interval failures {n_fail}")

    rows = list(out.items())
    g101 = [(k, r["param_spread"]) for k, r in rows
            if r["all_measurable"] and r["param_spread"] > G101_BAR]
    g102 = [(k, r["param_spread"]) for k, r in rows
            if not r["has_bad_verdict"] and r["param_spread"] > G102_BAR]
    g103 = [(k, r["n_interval_failures"]) for k, r in rows if r["n_interval_failures"]]
    clean = [r for _, r in rows if r["all_measurable"]]
    accepted = [r for _, r in rows if not r["has_bad_verdict"]]
    # does the VERDICT itself depend on the seed?
    verdict_unstable = [k for k, r in rows
                        if len({tuple(sorted(v.items()))
                                for v in r["verdicts_by_seed"].values()}) > 1]

    L = ["# Do the Carbopol records move when only the fitter's seed changes?", "",
         "Generated by `scripts/multistart_seed_scatter_carbopol.py`; every number "
         "reproduces by running it. NEXT_STEPS item 10, closing the scope gap that "
         "`solver.openfoam.check_exportable`'s docstring names: the blood seed study "
         "never covered these 15 records.", "",
         "These records are independent of the blood set in every way that matters -- "
         "different material, different instrument, 54 points per curve against 16, and "
         "**the shipped values come from a different fitter** (`carbopol/fit_hb.py`, "
         "`curve_fit` on stress residuals). So seed scatter here does not move what the "
         "record ships; it moves the SIDECAR, which is why G10.3 exists.", "",
         "## Pre-registered criteria", "",
         "| criterion | bar | measured | verdict |", "|---|---|---|---|",
         f"| G10.1 all-`measurable` records do not wander | spread <= {G101_BAR:g} | "
         f"{len(clean)}/15 all-measurable, worst "
         f"{max((r['param_spread'] for r in clean), default=float('nan')):.3e} | "
         f"{'PASS' if not g101 else '**FAILED**: ' + str(g101)} |",
         f"| G10.2 nothing the export accepts wanders | spread <= {G102_BAR:g} | "
         f"{len(accepted)}/15 accepted, worst "
         f"{max((r['param_spread'] for r in accepted), default=float('nan')):.3e} | "
         f"{'PASS' if not g102 else '**FAILED**: ' + str(g102)} |",
         f"| G10.3 every seed's interval contains the shipped csv value | 0 failures over "
         f"5 seeds x 3 parameters x 15 records | "
         f"{sum(r['n_interval_failures'] for _, r in rows)} failures | "
         f"{'PASS' if not g103 else '**FAILED**: ' + str(g103)} |", "",
         f"The objective moves by at most "
         f"{max(r['cost_spread'] for _, r in rows):.2e} relative across all 15 records, so "
         "nothing here is a failed search; any parameter movement is flatness.", "",
         f"The VERDICT SET itself is seed-dependent on {len(verdict_unstable)} of 15 "
         "records"
         + (f": {verdict_unstable}. A verdict that changes with the seed is a verdict "
            "about the optimiser, and the sidecar records only seed 0's."
            if verdict_unstable else
            ", so the sidecar's verdicts are a property of the data rather than of the "
            "seed."), "",
         "## Per record", "",
         "| record | N | spread | worst parameter | cost spread | verdicts | interval "
         "failures |", "|---|---|---|---|---|---|---|"]
    for k, r in rows:
        L.append(f"| `{k}` | {r['n_points']} | {r['param_spread']:.3e} | "
                 f"`{r['worst_param']}` | {r['cost_spread']:.2e} | "
                 f"{', '.join(f'{a}={b}' for a, b in r['verdicts'].items())} | "
                 f"{r['n_interval_failures']} |")
    L += ["",
          "**What G10.3 is for.** `write_identifiability` refuses to write a sidecar whose "
          "interval does not contain the shipped csv value, and the loader then refuses a "
          "record whose value falls outside its sidecar interval. Both have only ever been "
          "evaluated at seed 0. If containment held only there, the record's validity "
          "would be an accident of the seed rather than a statement about the fit.", ""]

    (OUT / "carbopol_multistart_scatter.md").write_text("\n".join(L), encoding="utf-8",
                                                        newline="\n")
    (OUT / "carbopol_multistart_scatter.json").write_text(json.dumps(
        {"seeds": list(SEEDS), "G101_bar": G101_BAR, "G102_bar": G102_BAR,
         "G101_violations": g101, "G102_violations": g102, "G103_violations": g103,
         "n_all_measurable": len(clean), "n_export_accepted": len(accepted),
         "worst_spread_all_measurable":
             max((r["param_spread"] for r in clean), default=float("nan")),
         "worst_spread_export_accepted":
             max((r["param_spread"] for r in accepted), default=float("nan")),
         "worst_cost_spread": max(r["cost_spread"] for _, r in rows),
         "total_interval_failures": sum(r["n_interval_failures"] for _, r in rows),
         "verdict_unstable_records": verdict_unstable,
         "records": out}, indent=2), encoding="ascii")
    print("\n".join(L[:22]))


if __name__ == "__main__":
    main()
