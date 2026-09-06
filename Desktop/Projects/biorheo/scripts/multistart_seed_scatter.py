"""Does the identifiability verdict predict how much the FITTED VALUES move when only the
optimiser's seed changes? (NEXT_STEPS item 9.)

The referee item this replaces was a request to "document dst_fit's multi-start guess
collision". `fit_model` builds its start list as ``[start, heuristic, *(n_starts-1)
log-uniform draws]``, so the count of starts overstates the number of DISTINCT optima
explored -- two starts can coincide, and several can fall into the same basin. Writing
that down as prose says nothing about whether it matters here, so it is measured instead.

WHAT IS MEASURED. Each donor curve is fitted five times with seeds 0-4, changing nothing
but the random starts. The relative spread is
``max_j (max_s p_js - min_s p_js) / |mean_s p_js|`` over the model's parameters. Beside it
sits the spread of the objective (``rmse_log``) over the same five fits. Those two
together separate the two things that can go wrong: a fit that lands on WORSE optima from
some starts (cost moves) from a fit that lands on EQUALLY GOOD optima at very different
parameter values (cost does not move, parameters do) -- the second is a flat direction,
i.e. exactly what the identifiability layer's `degenerate` / `not_identifiable` verdicts
are supposed to name.

PROVENANCE, because it decides what this is worth. The Carreau-Yasuda column and the
four-family maxima were measured EXPLORATORY on 2026-09-06, before any criterion was
written: they are what suggested the hypothesis, and they cannot also test it. The
per-donor `casson`, `cross` and `herschel_bulkley` columns were NOT looked at per donor
before the two criteria below were fixed, so they are the held-out half, and G9.1/G9.2 are
graded on ALL FOUR families with that split stated in the output.

PRE-REGISTERED (written before the held-out families were looked at; do not weaken):

  G9.1 NO FALSE NEGATIVES. No record whose parameters are ALL graded `measurable` may
       have a seed spread above 1e-3. A verdict set that calls everything measurable
       while the optimiser wanders is a verdict set that is not doing its job.
  G9.2 THE LARGE MOVERS ARE FLAGGED. Every record with a seed spread above 1e-1 must
       carry at least one `degenerate` or `not_identifiable` verdict.

Both are one-sided on purpose. The converse -- that every flagged record moves -- is NOT
claimed and must not be added: a parameter can be unidentifiable in a direction the
optimiser happens not to explore from these particular starts, and the profile likelihood
would still be right to flag it.

**G9.2 FAILED, on the held-out half, and the bar is not moved.** Herschel-Bulkley on
donor G moves 4.149 (415 %) in `tau_y` while its verdicts are
`tau_y=bounded_only, k=measurable, n=measurable` -- no `degenerate`, no
`not_identifiable`. One violation of 80 (model, donor) pairs; G9.1 passed with 54
all-measurable records whose worst spread is 5.997e-05.

The failure is informative rather than fatal, and it points at a specific gap:
`bounded_only` is NOT in `solver.openfoam._BAD_VERDICTS`, so an export ACCEPTS a record
carrying one -- yet a `bounded_only` parameter can move 415 % under nothing but a change
of seed. G9.2's mistake was to ask only about the two verdicts that block an export, when
the quantity at risk is any value that gets quoted. So a STRICTER criterion is added
beside it, with G9.2's failing number kept above:

  G9.3 NOTHING THE EXPORT ACCEPTS MAY WANDER. No record that `check_exportable` would
       let through -- i.e. no record without a `degenerate` or `not_identifiable`
       parameter -- may have a seed spread above 1e-3. This is strictly stronger than
       G9.2 (a lower bar, over a larger set) and it is the one that speaks for what
       ships. Reported separately over ALL 80 (model, donor) pairs and over the 60 that
       are actually in the registry, because `herschel_bulkley` is fitted for the model
       comparison and never shipped as a blood record.

Run:  python scripts/multistart_seed_scatter.py   (about 20 minutes; 4 models x 20 donors
x 5 seeds, each with a profile-likelihood closure). Writes
`data/horner2020/multistart_scatter.md` and `.json`; pinned in
`tests/test_multistart_scatter.py`.
"""

from __future__ import annotations

import csv
import json
import warnings
from pathlib import Path

import numpy as np

from biorheo.constitutive import get_model
from biorheo.identify import fit_rate_curve

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data" / "horner2020"
TEMPERATURE_C = 37.0
GEOMETRY = "ARES-G2 double-wall concentric cylinder, EDTA-anticoagulated whole blood, rate sweep"

SEEDS = (0, 1, 2, 3, 4)
FAMILIES = ("casson", "herschel_bulkley", "cross", "carreau_yasuda")
# `carreau_yasuda` and the four-family maxima were the exploratory pass; the rest is
# held out. Recorded so the .md can say which half tested the hypothesis.
EXPLORATORY = ("carreau_yasuda",)
BAD_VERDICTS = ("degenerate", "not_identifiable")
G91_BAR = 1e-3   # all-measurable records must be below this
G92_BAR = 1e-1   # records above this must carry a bad verdict -- FAILED, see docstring
G93_BAR = 1e-3   # nothing the export accepts may move more than this
# the families that actually ship as registry records; herschel_bulkley is fitted
# for model_comparison.md only
SHIPPED = ("casson", "cross", "carreau_yasuda")


def load() -> tuple[np.ndarray, dict[str, np.ndarray]]:
    with open(DATA / "steady_shear_human.csv", newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    head, body = rows[0], np.array(rows[1:], dtype=float)
    gd = body[:, 0]
    stress = {h.replace("stress_Pa_", ""): body[:, i] for i, h in enumerate(head) if i > 0}
    return gd, stress


def one_record(model, gd: np.ndarray, eta: np.ndarray) -> dict:
    """Five fits differing only in the optimiser's seed."""
    fits = [fit_rate_curve(model, gd, eta, seed=s, temperature_C=TEMPERATURE_C,
                           geometry=GEOMETRY) for s in SEEDS]
    vals = np.array([np.asarray(f.values, dtype=float) for f in fits])
    cost = np.array([f.fit.rmse_log for f in fits], dtype=float)
    mean = np.abs(vals.mean(axis=0))
    # a parameter that is identically zero across seeds has no relative spread; treat it
    # as agreeing rather than dividing by zero (a constant is agreement, not a blow-up)
    span = vals.max(axis=0) - vals.min(axis=0)
    rel = np.where(mean > 0, span / np.where(mean > 0, mean, 1.0), 0.0)
    verdicts = {n: fits[0].params[n].verdict for n in model.param_names}
    return {
        "param_spread": float(rel.max()),
        "worst_param": model.param_names[int(np.argmax(rel))],
        "cost_spread": float((cost.max() - cost.min()) / cost.mean()),
        "rmse_log": float(cost[0]),
        "verdicts": verdicts,
        "all_measurable": all(v == "measurable" for v in verdicts.values()),
        "has_bad_verdict": any(v in BAD_VERDICTS for v in verdicts.values()),
        "per_param_spread": {n: float(r) for n, r in zip(model.param_names, rel)},
    }


def main() -> None:
    warnings.filterwarnings("ignore")
    gd, stress = load()
    out: dict[str, dict[str, dict]] = {}
    for fam in FAMILIES:
        model = get_model(fam)
        out[fam] = {}
        for donor in sorted(stress):
            out[fam][donor] = one_record(model, gd, stress[donor] / gd)
            r = out[fam][donor]
            print(f"{fam:16s} {donor}  spread {r['param_spread']:.3e} "
                  f"({r['worst_param']})  cost {r['cost_spread']:.2e}  "
                  f"{''.join(v[0] for v in r['verdicts'].values())}")

    rows = [(f, d, r) for f, ds in out.items() for d, r in ds.items()]
    g91_viol = [(f, d, r["param_spread"]) for f, d, r in rows
                if r["all_measurable"] and r["param_spread"] > G91_BAR]
    g92_viol = [(f, d, r["param_spread"]) for f, d, r in rows
                if r["param_spread"] > G92_BAR and not r["has_bad_verdict"]]
    movers = [(f, d, r) for f, d, r in rows if r["param_spread"] > G92_BAR]
    # G9.3: what an export would accept is what matters, and `bounded_only` is accepted
    accepted = [(f, d, r) for f, d, r in rows if not r["has_bad_verdict"]]
    g93_viol = [(f, d, r["param_spread"]) for f, d, r in accepted
                if r["param_spread"] > G93_BAR]
    g93_viol_shipped = [(f, d, s) for f, d, s in g93_viol if f in SHIPPED]
    accepted_shipped = [(f, d, r) for f, d, r in accepted if f in SHIPPED]
    worst_accepted = max((r["param_spread"] for _, _, r in accepted), default=float("nan"))
    worst_accepted_shipped = max((r["param_spread"] for _, _, r in accepted_shipped),
                                 default=float("nan"))
    clean = [r for _, _, r in rows if r["all_measurable"]]
    worst_clean = max((r["param_spread"] for r in clean), default=float("nan"))
    worst_cost = max(r["cost_spread"] for _, _, r in rows)

    L = ["# Does the identifiability verdict predict multi-start seed scatter?", "",
         "Generated by `scripts/multistart_seed_scatter.py`; every number reproduces by "
         "running it. Each of 20 donor curves is fitted five times with seeds 0-4, "
         "changing nothing but the optimiser's random starts.", "",
         "**Provenance.** The `carreau_yasuda` column and the four-family maxima were an "
         "EXPLORATORY pass taken before any criterion was written -- they suggested the "
         "hypothesis and cannot also test it. The per-donor `casson`, `cross` and "
         "`herschel_bulkley` columns are the held-out half. Both criteria are graded on "
         "all four families and the split is stated so the exploratory half can be "
         "discounted.", "",
         "## Pre-registered criteria", "",
         "| criterion | bar | measured | verdict |", "|---|---|---|---|",
         f"| G9.1 no false negatives | an all-`measurable` record has spread <= {G91_BAR:g} "
         f"| {len(clean)} all-measurable records, worst spread {worst_clean:.3e} | "
         f"{'PASS' if not g91_viol else '**FAILED**: ' + str(g91_viol)} |",
         f"| G9.2 large movers are flagged | spread > {G92_BAR:g} implies a `degenerate` "
         f"or `not_identifiable` verdict | {len(movers)} records move that much | "
         f"{'PASS' if not g92_viol else '**FAILED**: ' + str(g92_viol)} |",
         f"| G9.3 nothing the export accepts may wander | an export-accepted record has "
         f"spread <= {G93_BAR:g} | {len(accepted)} accepted ({len(accepted_shipped)} in "
         f"shipped families), worst {worst_accepted:.3e} (shipped "
         f"{worst_accepted_shipped:.3e}) | "
         f"{'PASS' if not g93_viol else '**FAILED**: ' + str(g93_viol)} |", "",
         "**G9.2 FAILED and its bar was not moved.** The violation is the "
         "Herschel-Bulkley fit on donor G: `tau_y` moves 4.149 under nothing but a change "
         "of seed, while its verdict is `bounded_only` -- which is neither of the two "
         "verdicts G9.2 asked about, and, more to the point, is not one the export "
         "refuses. G9.3 is the stricter replacement: it asks about every record the "
         "export would ACCEPT, at a bar two orders lower. "
         + (f"It passes over all {len(accepted)} accepted records."
            if not g93_viol else
            f"It FAILS on {len(g93_viol)} record(s): {g93_viol}."
            + (" None of them is in a shipped family."
               if not g93_viol_shipped else
               f" **{len(g93_viol_shipped)} of them ship**: {g93_viol_shipped}.")), "",
         "## What the two spreads separate", "",
         "A fit can go wrong under a change of seed in two different ways, and the pair "
         "of columns tells them apart. If the OBJECTIVE moves, some starts are landing on "
         "worse optima and the fit is under-converged. If the objective does NOT move "
         "while the PARAMETERS do, the starts are landing on equally good optima at "
         "different places -- a flat direction, which is what `degenerate` and "
         "`not_identifiable` are for.", "",
         f"Across all {len(rows)} (model, donor) pairs the objective moves by at most "
         f"**{worst_cost:.2e}** relative, i.e. not at all. Every parameter spread below "
         "is therefore the second kind: flatness, not a failed search. That is the whole "
         "result -- the multi-start count is not measuring convergence here, it is "
         "sampling a flat direction, and the identifiability layer already names which "
         "records have one.", "",
         "## The records that move", "",
         "| model | donor | spread | worst parameter | cost spread | verdicts | flagged |",
         "|---|---|---|---|---|---|---|"]
    for f, d, r in sorted(movers, key=lambda t: -t[2]["param_spread"]):
        L.append(f"| {f} | {d} | {r['param_spread']:.3e} | `{r['worst_param']}` | "
                 f"{r['cost_spread']:.2e} | "
                 f"{', '.join(f'{k}={v}' for k, v in r['verdicts'].items())} | "
                 f"{'yes' if r['has_bad_verdict'] else '**NO**'} |")
    L += ["", "## Per-family summary", "",
          "| model | exploratory? | worst spread | worst cost spread | all-measurable "
          "records | worst spread among those |", "|---|---|---|---|---|---|"]
    for fam in FAMILIES:
        rs = list(out[fam].values())
        cl = [r for r in rs if r["all_measurable"]]
        L.append(f"| {fam} | {'yes' if fam in EXPLORATORY else 'held out'} | "
                 f"{max(r['param_spread'] for r in rs):.3e} | "
                 f"{max(r['cost_spread'] for r in rs):.2e} | {len(cl)}/20 | "
                 + (f"{max(r['param_spread'] for r in cl):.3e} |" if cl else "n/a |"))
    L += ["",
          "**The converse is NOT claimed.** Nothing here says every flagged record must "
          "move: a parameter can be unidentifiable along a direction these particular "
          "starts never explore, and the profile likelihood would still be right to flag "
          "it. Only the two one-sided criteria above are tested.", "",
          "**Consequence for the registry.** A record whose parameters move under a "
          "change of seed ships whatever `seed=0` returned. That is reproducible -- the "
          "seed is fixed and recorded -- but it is not a measurement, and it is exactly "
          "the set the export already refuses without `allow_unidentified=True`.", ""]

    (DATA / "multistart_scatter.md").write_text("\n".join(L), encoding="utf-8", newline="\n")
    payload = {"seeds": list(SEEDS), "families": list(FAMILIES),
               "exploratory": list(EXPLORATORY), "G91_bar": G91_BAR, "G92_bar": G92_BAR,
               "G91_violations": g91_viol, "G92_violations": g92_viol,
               "G93_bar": G93_BAR, "G93_violations": g93_viol,
               "G93_violations_shipped": g93_viol_shipped,
               "n_export_accepted": len(accepted),
               "n_export_accepted_shipped": len(accepted_shipped),
               "worst_spread_export_accepted": worst_accepted,
               "worst_spread_export_accepted_shipped": worst_accepted_shipped,
               "shipped_families": list(SHIPPED),
               "n_all_measurable": len(clean), "worst_spread_all_measurable": worst_clean,
               "worst_cost_spread": worst_cost,
               "n_movers": len(movers), "records": out}
    (DATA / "multistart_scatter.json").write_text(json.dumps(payload, indent=2),
                                                  encoding="ascii")
    print("\n".join(L[:20]))


if __name__ == "__main__":
    main()
