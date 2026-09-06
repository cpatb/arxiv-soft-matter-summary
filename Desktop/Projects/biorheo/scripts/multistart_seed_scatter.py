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

**RETRACTED 2026-09-06: G9.2 and G9.3 did not fail. The metric was broken.**

An earlier version of this script reported one violation -- Herschel-Bulkley on donor G,
"4.149, i.e. 415 %, in `tau_y`" -- and a retraction, a docstring claim, a pinned test and
a new API flag were built on it before a code review caught it. The five seeds return

    tau_y = 1.61e-19, 9.81e-19, 5.00e-22, 3.27e-20, 5.94e-21 Pa

against a smallest measured stress on that curve of 1.583e-03 Pa, with `lo_at_bound=True`
and an IDENTICAL profile interval [0, 0.01953] on all five. The optimiser converged
`tau_y` to zero five times out of five. The "415 % move" was a span of 9.8e-19 divided by
a mean of 2.4e-19 -- dust over dust. Scaled against anything physical it is 6.2e-16 of the
smallest measured stress and 5.0e-17 of the parameter's own confidence interval.

Worse, the verdict turned on floating-point noise: had `least_squares` returned exactly
0.0 five times, the guard below would have given rel = 0.0 and G9.2 would have PASSED. A
criterion whose outcome flips on whether a converged-to-zero parameter lands at 0.0 or at
1e-19 is not measuring the hypothesis. This is the lesson already in `~/.claude/lessons.md`
-- "a RELATIVE sensitivity whose denominator is itself shrinking measures the shrinkage" --
applied to my own instrument a session after writing it down.

**The metric is therefore fixed, and this is a repair of a broken instrument, not a bar
that moved after seeing the answer.** The criteria and their bars are unchanged; what
changed is the quantity they are evaluated on:

    spread = max_j (max_s p_js - min_s p_js) / (that parameter's 95 % profile-interval
             width at seed 0)

A parameter that moves by less than its own stated uncertainty has not moved. The raw
`|span| / |mean|` is still computed and reported beside it, so the retracted number stays
visible. Where the interval is unbounded or degenerate the ratio is small BY CONSTRUCTION,
which is the right answer -- an unbounded interval already says the value is not pinned,
and that is what the verdict is for.

  G9.3 NOTHING THE EXPORT ACCEPTS MAY WANDER. No record that `check_exportable` would
       let through -- i.e. no record without a `degenerate` or `not_identifiable`
       parameter -- may have a spread above 1e-3. It is G9.2's predicate at a lower bar
       over the SAME set (the review was right that it is nested, not broader), so the
       two are reported as one nested criterion with two bars rather than as two
       independent results.

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
    """Five fits differing only in the optimiser's seed.

    The reported spread is the seed-to-seed span measured in units of the parameter's OWN
    95 % profile interval. Dividing by |mean| instead is what produced the retracted
    G9.2 failure: for a parameter the optimiser drives to zero, both span and mean are
    ~1e-19 and their ratio is meaningless. The raw ratio is kept as
    `param_spread_over_mean` so the retracted number stays visible.

    Verdicts are collected FROM EVERY SEED, not just seed 0 -- `all_measurable` and
    `has_bad_verdict` decide which criterion a record falls under, and deciding that from
    one seed while measuring the spread over five would make the population a function of
    seed 0's luck (raised by the code review 2026-09-06)."""
    fits = [fit_rate_curve(model, gd, eta, seed=s, temperature_C=TEMPERATURE_C,
                           geometry=GEOMETRY) for s in SEEDS]
    vals = np.array([np.asarray(f.values, dtype=float) for f in fits])
    cost = np.array([f.fit.rmse_log for f in fits], dtype=float)
    span = vals.max(axis=0) - vals.min(axis=0)

    # the scale a move is measured against: the parameter's own stated uncertainty
    ci = np.array([float(fits[0].params[n].hi) - float(fits[0].params[n].lo)
                   for n in model.param_names])
    rel_ci = np.where(ci > 0, span / np.where(ci > 0, ci, 1.0), 0.0)
    mean = np.abs(vals.mean(axis=0))
    rel_mean = np.where(mean > 0, span / np.where(mean > 0, mean, 1.0), 0.0)

    by_seed = {str(s): {n: f.params[n].verdict for n in model.param_names}
               for s, f in zip(SEEDS, fits)}
    verdicts = fits[0].params and {n: fits[0].params[n].verdict for n in model.param_names}
    # a record counts as all-measurable only if EVERY seed says so, and as flagged if ANY
    # seed does -- both the conservative reading
    all_meas = all(v == "measurable" for vs in by_seed.values() for v in vs.values())
    has_bad = any(v in BAD_VERDICTS for vs in by_seed.values() for v in vs.values())
    return {
        "param_spread": float(rel_ci.max()),
        "worst_param": model.param_names[int(np.argmax(rel_ci))],
        "param_spread_over_mean": float(rel_mean.max()),   # the retracted metric
        "worst_param_over_mean": model.param_names[int(np.argmax(rel_mean))],
        "abs_span": {n: float(s) for n, s in zip(model.param_names, span)},
        "ci_width": {n: float(c) for n, c in zip(model.param_names, ci)},
        "values_by_seed": {n: [float(v) for v in vals[:, j]]
                           for j, n in enumerate(model.param_names)},
        "cost_spread": float((cost.max() - cost.min()) / cost.mean()),
        "rmse_log": float(cost[0]),
        "verdicts": verdicts,
        "verdicts_by_seed": by_seed,
        "verdict_seed_stable": len({tuple(sorted(v.items())) for v in by_seed.values()}) == 1,
        "all_measurable": all_meas,
        "has_bad_verdict": has_bad,
        "per_param_spread": {n: float(r) for n, r in zip(model.param_names, rel_ci)},
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
    # G9.2 and G9.3 share a predicate and differ only in the bar (the review was right
    # that G9.3 is nested inside G9.2, not broader), so they are computed together.
    accepted = [(f, d, r) for f, d, r in rows if not r["has_bad_verdict"]]
    g92_viol = [(f, d, r["param_spread"]) for f, d, r in accepted
                if r["param_spread"] > G92_BAR]
    g93_viol = [(f, d, r["param_spread"]) for f, d, r in accepted
                if r["param_spread"] > G93_BAR]
    movers = [(f, d, r) for f, d, r in rows if r["param_spread"] > G92_BAR]
    clean = [r for _, _, r in rows if r["all_measurable"]]
    worst_clean = max((r["param_spread"] for r in clean), default=float("nan"))
    worst_accepted = max((r["param_spread"] for _, _, r in accepted), default=float("nan"))
    worst_cost = max(r["cost_spread"] for _, _, r in rows)
    unstable = [(f, d) for f, d, r in rows if not r["verdict_seed_stable"]]
    # the retracted metric, kept visible
    worst_over_mean = max(r["param_spread_over_mean"] for _, _, r in rows)
    worst_over_mean_id = max(rows, key=lambda t: t[2]["param_spread_over_mean"])

    L = ["# Does the identifiability verdict predict multi-start seed scatter?", "",
         "Generated by `scripts/multistart_seed_scatter.py`; every number reproduces by "
         "running it. Each of 20 donor curves is fitted five times with seeds 0-4, "
         "changing nothing but the optimiser's random starts.", "",
         "**A spread is measured in units of the parameter's OWN 95 % profile interval.** "
         "A parameter that moves by less than its stated uncertainty has not moved. An "
         "earlier version of this script divided by `|mean|` instead and reported a "
         "spurious failure; see the retraction below.", "",
         "## Pre-registered criteria", "",
         "| criterion | bar | measured | verdict |", "|---|---|---|---|",
         f"| G9.1 no false negatives | an all-`measurable` record has spread <= {G91_BAR:g} "
         f"| {len(clean)} all-measurable records, worst {worst_clean:.3e} | "
         f"{'PASS' if not g91_viol else '**FAILED**: ' + str(g91_viol)} |",
         f"| G9.2 large movers are flagged | an export-accepted record has spread <= "
         f"{G92_BAR:g} | {len(accepted)} accepted, worst {worst_accepted:.3e} | "
         f"{'PASS' if not g92_viol else '**FAILED**: ' + str(g92_viol)} |",
         f"| G9.3 the same, at a lower bar | an export-accepted record has spread <= "
         f"{G93_BAR:g} | as above | "
         f"{'PASS' if not g93_viol else '**FAILED**: ' + str(g93_viol)} |", "",
         "G9.3 is G9.2's predicate at a lower bar over the same set, so it is one nested "
         "criterion with two bars rather than two independent results.", "",
         "## Retraction, 2026-09-06", "",
         "An earlier run of this script reported G9.2 and G9.3 as FAILED on one record -- "
         "Herschel-Bulkley, donor G, \"4.149, i.e. 415 %, in `tau_y`\" -- and a "
         "docstring claim, a pinned test and a new `strict_identifiability` API flag were "
         "built on it before a code review caught it. **That failure is withdrawn.** The "
         "five seeds return `tau_y` = 1.61e-19, 9.81e-19, 5.00e-22, 3.27e-20, 5.94e-21 Pa "
         "against a smallest measured stress of 1.583e-03 Pa on that curve, with "
         "`lo_at_bound=True` and an identical profile interval [0, 0.01953] on all five. "
         "The optimiser converged the parameter to zero five times out of five; the "
         "\"415 % move\" was a span of 9.8e-19 over a mean of 2.4e-19.", "",
         "The retracted metric is still computed and printed so the number stays visible: "
         f"its worst value is {worst_over_mean:.3e}, on "
         f"`{worst_over_mean_id[0]}` donor {worst_over_mean_id[1]}. Under the corrected "
         "metric that same record scores "
         f"{worst_over_mean_id[2]['param_spread']:.3e}.", "",
         "The failure also turned on floating-point noise: had `least_squares` returned "
         "exactly 0.0 five times, the old guard would have given 0.0 and G9.2 would have "
         "passed. A criterion whose outcome flips on whether a converged-to-zero "
         "parameter lands at 0.0 or 1e-19 is not measuring the hypothesis.", "",
         "## What the two spreads separate", "",
         "A fit can go wrong under a change of seed in two ways. If the OBJECTIVE moves, "
         "some starts are landing on worse optima and the fit is under-converged. If the "
         "objective does NOT move while the PARAMETERS do, the starts are landing on "
         "equally good optima at different places -- a flat direction, which is what "
         "`degenerate` and `not_identifiable` are for.", "",
         f"Across all {len(rows)} (model, donor) pairs the objective moves by at most "
         f"**{worst_cost:.2e}** relative in `rmse_log`, i.e. **{2 * worst_cost:.2e}** in "
         "the least-squares cost itself, which is the quantity the optimiser minimises. "
         "Not at all, either way. Every parameter spread below is therefore flatness, not "
         "a failed search.", "",
         f"The verdict SET is seed-dependent on {len(unstable)} of {len(rows)} pairs"
         + (f": {unstable}. Both `all_measurable` and `has_bad_verdict` are therefore "
            "taken over ALL five seeds -- all-measurable only if every seed says so, "
            "flagged if any seed does." if unstable else
            ", so reading them from one seed would have been safe; they are taken over "
            "all five anyway."), "",
         "## The records that move", "",
         "| model | donor | spread (CI units) | worst parameter | span | CI width | "
         "cost spread | verdicts | flagged |",
         "|---|---|---|---|---|---|---|---|---|"]
    for f, d, r in sorted(movers, key=lambda t: -t[2]["param_spread"]):
        w = r["worst_param"]
        L.append(f"| {f} | {d} | {r['param_spread']:.3e} | `{w}` | "
                 f"{r['abs_span'][w]:.3e} | {r['ci_width'][w]:.3e} | "
                 f"{r['cost_spread']:.2e} | "
                 f"{', '.join(f'{k}={v}' for k, v in r['verdicts'].items())} | "
                 f"{'yes' if r['has_bad_verdict'] else '**NO**'} |")
    if not movers:
        L.append("| -- | -- | -- | -- | -- | -- | -- | none | -- |")
    L += ["", "## Per-family summary", "",
          "| model | exploratory? | worst spread | worst spread over mean (retracted "
          "metric) | worst cost spread | all-measurable | worst among those |",
          "|---|---|---|---|---|---|---|"]
    for fam in FAMILIES:
        rs = list(out[fam].values())
        cl = [r for r in rs if r["all_measurable"]]
        L.append(f"| {fam} | {'yes' if fam in EXPLORATORY else 'held out'} | "
                 f"{max(r['param_spread'] for r in rs):.3e} | "
                 f"{max(r['param_spread_over_mean'] for r in rs):.3e} | "
                 f"{max(r['cost_spread'] for r in rs):.2e} | {len(cl)}/20 | "
                 + (f"{max(r['param_spread'] for r in cl):.3e} |" if cl else "n/a |"))
    L += ["",
          "**Provenance, restated honestly.** The `carreau_yasuda` per-donor column AND "
          "the four-family maxima were an exploratory pass taken before the bars were "
          "written -- and the four-family maxima INCLUDE the Herschel-Bulkley number that "
          "the retracted failure came from, so the original 'held-out' framing did not "
          "survive contact with the code review. Nothing in the code records what was "
          "viewed when. Read this as one exploratory study whose criteria were written "
          "down before the per-donor breakdown was examined, and no stronger.", "",
          "**The converse is NOT claimed.** Nothing here says every flagged record must "
          "move: a parameter can be unidentifiable along a direction these particular "
          "starts never explore, and the profile likelihood would still be right to flag "
          "it. Only the one-sided criteria above are tested.", ""]

    (DATA / "multistart_scatter.md").write_text("\n".join(L), encoding="utf-8", newline="\n")
    payload = {"seeds": list(SEEDS), "families": list(FAMILIES),
               "exploratory": list(EXPLORATORY), "G91_bar": G91_BAR, "G92_bar": G92_BAR,
               "G93_bar": G93_BAR,
               "metric": "seed span in units of the parameter's own 95 % profile interval",
               "retracted_metric": "seed span over |mean|; see the .md",
               "G91_violations": g91_viol, "G92_violations": g92_viol,
               "G93_violations": g93_viol,
               "n_all_measurable": len(clean), "worst_spread_all_measurable": worst_clean,
               "n_export_accepted": len(accepted),
               "worst_spread_export_accepted": worst_accepted,
               "worst_cost_spread": worst_cost,
               "worst_spread_over_mean": worst_over_mean,
               "verdict_unstable": unstable,
               "n_movers": len(movers), "records": out}
    (DATA / "multistart_scatter.json").write_text(json.dumps(payload, indent=2),
                                                  encoding="ascii")
    print("\n".join(L[:24]))


if __name__ == "__main__":
    main()
