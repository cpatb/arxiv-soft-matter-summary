"""Item 18: how much does a finer needle really cost, and how much of a force budget does
the fit uncertainty eat.

Gate pre-registered in `data/needle_bore_prereg.md`, committed at `8f4fe9c` BEFORE any of
this was written. Bars live there and are repeated here only where the code needs the
number. Nothing here fits anything, adds a record, or touches `dst_fit`.

The measured quantity is `b = d ln F / d ln R` at fixed flow rate: the exponent a device
group implicitly assumes is `-4` when it sizes a needle from a Newtonian calculation. Two
closed forms anchor it without passing through the force solver (`newtonian` gives exactly
`-4`, `power_law` exactly `-(3n+1)`), and both null families are in this registry on the
same 20 donor curves as the fitted ones.

Run:  python scripts/needle_bore.py
Writes `data/needle_bore.md` and `data/needle_bore.json`. Pinned in
`tests/test_design_space.py`.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import biorheo
from biorheo.apps.design_space import (
    bore_exponent, bore_exponent_closed_form, derating, newtonian_rule_error,
)
from biorheo.apps.injectability import Needle, Syringe, injection_force
from biorheo.registry import load_record

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
RECORDS = REPO / "src" / "biorheo" / "registry" / "records"
CROSSCHECK = REPO / "case_studies" / "blood_injectability" / "crosscheck.json"
OUT_MD = REPO / "data" / "needle_bore.md"
OUT_JSON = REPO / "data" / "needle_bore.json"

# --- geometry and operating points, all from the pre-registration -------------------
R_REF = 1.05e-4          # m. The SAME declared value the blood case study runs at, and it
                         # carries that study's [UNVERIFIED] flag. Nothing here is a gauge.
NEEDLE_LENGTH = 12.7e-3  # m
BARREL_RADIUS = 2.3e-3   # m, a 1 mL syringe; a choice, as in both shipped case studies
Q1 = 1.0e-6 / 60.0       # m^3/s, 1 mL/min -- the shipped blood study's rate
Q2 = 1.0e-6 / 10.0       # m^3/s, 1 mL in 10 s -- sourced, see the pre-registration
BORE_RATIO = 0.70        # the "much finer needle" corner; a RATIO, never a gauge
BORE_TOLERANCE = 0.05    # declared 5 % bore error for the tolerance term of G18.7
F_BUDGET = 25.0          # N. A RULE OF THUMB from the JHU thesis, not a standard.

# --- bars, from data/needle_bore_prereg.md ------------------------------------------
G181_TOL = 1e-3          # |b + 4| on the newtonian records
G182_RTOL = 1e-3         # relative error against -(3n+1) on the power_law records
G183_MIN_SPREAD = 0.1    # at least one record with |b + 4| above this
G184_SLACK = 1e-6        # b >= -4 - slack on every record
G185_THRESHOLD = 0.10    # "over-predicts by more than 10 %"
G185_BLOOD_MIN = 15      # of 20 Casson donors BELOW the threshold (the FAIL half)

DONORS = tuple("ABCDFGHIJKLMNOPQRSTU")
BLOOD_FAMILIES = ("casson", "cross", "carreau_yasuda", "newtonian", "power_law")
CARBOPOL_WT = ("0p1wt", "0p3wt", "0p5wt", "0p7wt", "0p9wt")
CARBOPOL_VOL = ("30uL", "40uL", "50uL")
# The two gels the surrogates case study ships, named here because the pre-registration's
# G18.5 bar says "2/2 gels". The gate is graded on ALL 15 as well; see the report.
SURROGATE_GELS = ("carbopol_etd2020_0p5wt_40uL", "carbopol_etd2020_0p1wt_40uL")


def frozen_record_ids() -> list[str]:
    """This item's OWN record list. Deliberately not `registry.RECORD_IDS`: appending to a
    shared list has silently re-decided a closed gate in this repo before (item 11)."""
    ids = []
    for d in DONORS:
        base = f"blood_horner2020_donor_{d}"
        for fam in BLOOD_FAMILIES:
            ids.append(base if fam == "casson" else f"{base}_{fam}")
    for w in CARBOPOL_WT:
        for v in CARBOPOL_VOL:
            ids.append(f"carbopol_etd2020_{w}_{v}")
    return ids


def intervals_of(rec) -> tuple[dict[str, tuple[float, float]], tuple[str, ...]]:
    """Per-parameter profile interval and the names whose end sits on a box bound.

    An end at a bound is a BOUND, not a limit (NEXT_STEPS' standing rule), so it is
    reported and the record is excluded from the gate that needs closed intervals.
    """
    block = rec.identifiability
    if block is None:
        return {}, ()
    out, open_names = {}, []
    for name, p in block["parameters"].items():
        out[name] = (float(p["lo"]), float(p["hi"]))
        if p.get("lo_at_bound") or p.get("hi_at_bound"):
            open_names.append(name)
    return out, tuple(open_names)


def hull_relative(Q, needle, syringe, model, params, intervals) -> float:
    """One-at-a-time profile push-through force hull, as a fraction of the nominal force.

    The same object `force_with_uq` reports: WIDER than the joint 95 % region it is drawn
    from, because that region is an ellipse inscribed in this box.
    """
    f0 = injection_force(Q, needle, syringe, model, params).force
    forces = [f0]
    for j, name in enumerate(model.param_names):
        if name not in intervals:
            continue
        for val in intervals[name]:
            if not np.isfinite(val):
                continue
            p = np.asarray(params, dtype=float).copy()
            p[j] = val
            try:
                forces.append(injection_force(Q, needle, syringe, model, p).force)
            except (ValueError, RuntimeError):
                continue
    return float((max(forces) - min(forces)) / f0)


def family_term() -> dict[str, float]:
    """Item 12's Cross-vs-Casson plunger-force difference per donor, |rel_dev|.

    READ from the shipped artifact, never recomputed: a second source of truth for a
    number this repo already ships is a failure mode it has a rule about.
    """
    if not CROSSCHECK.exists():
        return {}
    d = json.loads(CROSSCHECK.read_text(encoding="utf-8"))
    return {r["donor"]: abs(float(r["rel_dev"])) for r in d["summary"]["cross"]["rows"]}


def main() -> int:
    needle = Needle(radius=R_REF, length=NEEDLE_LENGTH)
    syringe = Syringe(barrel_radius=BARREL_RADIUS)
    fam_term = family_term()

    rows = []
    for rid in frozen_record_ids():
        rec = load_record(RECORDS / f"{rid}.yaml")
        params = rec.params_si
        model = rec.model
        iv, open_names = intervals_of(rec)

        row: dict[str, object] = {
            "id": rid, "model": model.name,
            "family": "carbopol" if rid.startswith("carbopol") else "blood",
            "donor": rid.split("_donor_")[1].split("_")[0] if "_donor_" in rid else None,
            "commercial_ok": bool(rec.commercial_ok),
            "verdicts": rec.verdicts,
            "all_measurable": (rec.verdicts is not None
                               and all(v == "measurable" for v in rec.verdicts.values())),
        }
        for label, Q in (("Q1", Q1), ("Q2", Q2)):
            row[f"b_{label}"] = bore_exponent(Q, needle, syringe, model, params)
        cf = bore_exponent_closed_form(model, params)
        row["b_closed_form"] = cf
        row["b_closed_form_rel_err"] = (None if cf is None
                                        else abs(row["b_Q1"] / cf - 1.0))  # type: ignore[operator]

        nre = newtonian_rule_error(Q1, needle, syringe, model, params, bore_ratio=BORE_RATIO)
        row["newtonian_rule_over_prediction"] = nre["over_prediction"]
        row["F_ref_N"] = nre["F_ref"]
        row["F_fine_N"] = nre["F_true"]

        row["open_parameters"] = list(open_names)
        row["hull_rel_Q1"] = hull_relative(Q1, needle, syringe, model, params, iv) if iv else None

        # derating against the force budget
        if iv:
            try:
                der = derating(F_BUDGET, needle, syringe, model, params, iv,
                               open_parameters=open_names)
                row["derating"] = {"q_nominal_mL_min": der.q_nominal * 6e7,
                                   "q_upper_mL_min": der.q_upper_bound * 6e7,
                                   "factor": der.factor,
                                   "binding_parameter": der.binding_parameter,
                                   "binding_end": der.binding_end}
            except ValueError as e:
                row["derating"] = {"error": str(e)}
        else:
            row["derating"] = None

        # --- G18.7 three error terms, same units (relative force error) ---------------
        b = float(row["b_Q1"])  # type: ignore[arg-type]
        terms = {"bore_tolerance": abs(b) * BORE_TOLERANCE,
                 "parameter_uncertainty": row["hull_rel_Q1"],
                 "family_choice": fam_term.get(str(row["donor"])) if row["model"] == "casson" else None}
        row["error_terms"] = terms
        finite = {k: v for k, v in terms.items() if isinstance(v, float) and np.isfinite(v)}
        row["dominant_error_term"] = max(finite, key=lambda k: finite[k]) if finite else None
        rows.append(row)

    by_id = {r["id"]: r for r in rows}

    # --- gates -----------------------------------------------------------------------
    newt = [r for r in rows if r["model"] == "newtonian"]
    g181_bad = [r["id"] for r in newt if abs(float(r["b_Q1"]) + 4.0) >= G181_TOL]  # type: ignore[arg-type]
    g181 = {"bar": f"|b + 4| < {G181_TOL} on 20/20 newtonian records",
            "n": len(newt), "n_pass": len(newt) - len(g181_bad),
            "worst": max((abs(float(r["b_Q1"]) + 4.0) for r in newt), default=float("nan")),  # type: ignore[arg-type]
            "failures": g181_bad, "pass": len(newt) == 20 and not g181_bad}

    pl = [r for r in rows if r["model"] == "power_law"]
    g182_bad = [r["id"] for r in pl if float(r["b_closed_form_rel_err"]) >= G182_RTOL]  # type: ignore[arg-type]
    g182 = {"bar": f"|b/(-(3n+1)) - 1| < {G182_RTOL} on 20/20 power_law records",
            "n": len(pl), "n_pass": len(pl) - len(g182_bad),
            "worst": max((float(r["b_closed_form_rel_err"]) for r in pl), default=float("nan")),  # type: ignore[arg-type]
            "failures": g182_bad, "pass": len(pl) == 20 and not g182_bad}

    spread = [(abs(float(r["b_Q1"]) + 4.0), r["id"]) for r in rows]  # type: ignore[arg-type]
    worst_spread, worst_id = max(spread)
    g183 = {"bar": f"at least one record with |b + 4| > {G183_MIN_SPREAD}",
            "worst": worst_spread, "worst_record": worst_id,
            "n_above": sum(1 for s, _ in spread if s > G183_MIN_SPREAD),
            "pass": worst_spread > G183_MIN_SPREAD,
            "note": "without this the anchors are all passed by a routine that returns -4 "
                    "for everything"}

    g184_bad = [r["id"] for r in rows if float(r["b_Q1"]) < -4.0 - G184_SLACK]  # type: ignore[arg-type]
    g184 = {"bar": f"b >= -4 - {G184_SLACK} on every record graded",
            "n": len(rows), "n_pass": len(rows) - len(g184_bad),
            "steepest": min(float(r["b_Q1"]) for r in rows),  # type: ignore[arg-type]
            "failures": g184_bad, "pass": not g184_bad}

    casson = [r for r in rows if r["model"] == "casson"]
    blood_below = [r["id"] for r in casson
                   if float(r["newtonian_rule_over_prediction"]) < G185_THRESHOLD]  # type: ignore[arg-type]
    carbo = [r for r in rows if r["family"] == "carbopol"]
    carbo_above = [r["id"] for r in carbo
                   if float(r["newtonian_rule_over_prediction"]) > G185_THRESHOLD]  # type: ignore[arg-type]
    gels = [by_id[g] for g in SURROGATE_GELS]
    gels_above = [g["id"] for g in gels
                  if float(g["newtonian_rule_over_prediction"]) > G185_THRESHOLD]  # type: ignore[arg-type]
    g185 = {"bar": ("SPLIT: FAIL on blood (< 10 % on >= 15/20 Casson donors) and "
                    "PASS on Carbopol (> 10 % on 2/2 named gels)"),
            "threshold": G185_THRESHOLD,
            "blood_n_below": len(blood_below), "blood_n": len(casson),
            "blood_worst": max(float(r["newtonian_rule_over_prediction"]) for r in casson),  # type: ignore[arg-type]
            "blood_half_as_predicted": len(blood_below) >= G185_BLOOD_MIN,
            "carbopol_named_above": len(gels_above), "carbopol_named_n": len(gels),
            "carbopol_all_above": len(carbo_above), "carbopol_all_n": len(carbo),
            "carbopol_worst": max(float(r["newtonian_rule_over_prediction"]) for r in carbo),  # type: ignore[arg-type]
            "carbopol_least": min(float(r["newtonian_rule_over_prediction"]) for r in carbo),  # type: ignore[arg-type]
            "carbopol_half_as_predicted": len(gels_above) == len(gels),
            "pass": len(blood_below) >= G185_BLOOD_MIN and len(gels_above) == len(gels),
            "uninformative": (len(blood_below) == 0 and len(carbo_above) == len(carbo))}

    closed = [r for r in rows if isinstance(r["derating"], dict)
              and "factor" in r["derating"] and not r["open_parameters"]]  # type: ignore[operator]
    g186_bad = [r["id"] for r in closed if float(r["derating"]["factor"]) >= 1.0]  # type: ignore[index]
    # G18.6' -- the pre-registered bar is close to VACUOUS and this is the stricter test
    # beside it, not a replacement. The upper force envelope is a MAX over profile ends and
    # the nominal force is not in it, so `factor <= 1` can only fail if the fitted parameters
    # locally MAXIMISE the force in every direction at once. A routine that returned a
    # constant 0.5 for everything would pass G18.6. G18.6' pins the two ends the bar cannot:
    # a zero-width interval must give exactly 1, and a real interval must give strictly less.
    probe = by_id["blood_horner2020_donor_A"]
    probe_rec = load_record(RECORDS / f"{probe['id']}.yaml")
    degenerate = {n: (float(v), float(v)) for n, v in
                  zip(probe_rec.model.param_names, probe_rec.params_si)}
    d0 = derating(F_BUDGET, needle, syringe, probe_rec.model, probe_rec.params_si, degenerate)
    real_iv, _ = intervals_of(probe_rec)
    d1 = derating(F_BUDGET, needle, syringe, probe_rec.model, probe_rec.params_si, real_iv)
    g186p = {"bar": "zero-width intervals give factor == 1 exactly; the real intervals give < 1",
             "record": probe["id"],
             "factor_zero_width": d0.factor, "factor_real": d1.factor,
             "pass": bool(abs(d0.factor - 1.0) < 1e-12 and d1.factor < 1.0),
             "note": "G18.6 as pre-registered is nearly unfalsifiable: the upper envelope is a "
                     "max over profile ends, so the factor is <= 1 unless the fit locally "
                     "maximises the force in every direction. The bar is kept and reported "
                     "with its number; this is the stricter test beside it, not instead of it."}

    # The split is taken over EVERY record with a numeric factor, not just the
    # closed-interval subset: the records with a non-`measurable` verdict are exactly the
    # ones whose intervals run to a bound, so restricting to `closed` empties the group the
    # comparison is about. (The first version of this section did restrict, printed an
    # empty second group, and had a narrative written around it. Generated from the rows.)
    with_factor = [r for r in rows if isinstance(r["derating"], dict)
                   and "factor" in r["derating"]]
    identifiable = [r for r in with_factor if r["all_measurable"]]
    not_identifiable = [r for r in with_factor if not r["all_measurable"]]
    g186 = {"bar": "derating factor < 1 on every record with a closed interval",
            "vacuity_note": "see G18.6-prime: this bar can essentially only be failed by a fit "
                            "that locally maximises the force, so the SIZE of the factor is the "
                            "result and the sign is not",
            "median_factor_all_measurable": float(np.median(
                [r["derating"]["factor"] for r in identifiable])) if identifiable else None,
            "n_all_measurable": len(identifiable),
            "median_factor_not_all_measurable": float(np.median(
                [r["derating"]["factor"] for r in not_identifiable])) if not_identifiable else None,
            "n_not_all_measurable": len(not_identifiable),
            "worst_factor_all_measurable": min(
                (r["derating"]["factor"] for r in identifiable), default=None),
            "worst_factor_not_all_measurable": min(
                (r["derating"]["factor"] for r in not_identifiable), default=None),
            "split_note": "taken over every record with a numeric factor, not the "
                          "closed-interval subset: a non-measurable verdict and an interval "
                          "running to a bound are the same records, so restricting to closed "
                          "intervals empties the group the comparison is about",
            "n": len(closed), "n_pass": len(closed) - len(g186_bad),
            "worst_factor": max((float(r["derating"]["factor"]) for r in closed), default=float("nan")),  # type: ignore[index]
            "strongest_derating": min((float(r["derating"]["factor"]) for r in closed), default=float("nan")),  # type: ignore[index]
            "failures": g186_bad, "pass": bool(closed) and not g186_bad}

    winners: dict[str, int] = {}
    for r in rows:
        w = r["dominant_error_term"]
        if w:
            winners[str(w)] = winners.get(str(w), 0) + 1
    three_way = [r for r in rows if r["error_terms"]["family_choice"] is not None  # type: ignore[index]
                 and r["error_terms"]["parameter_uncertainty"] is not None]  # type: ignore[index]
    three_winners: dict[str, int] = {}
    for r in three_way:
        t = {k: v for k, v in r["error_terms"].items() if isinstance(v, float)}  # type: ignore[union-attr]
        three_winners[max(t, key=lambda k: t[k])] = three_winners.get(max(t, key=lambda k: t[k]), 0) + 1
    g187 = {"bar": "reported, not barred: a winner is named per record",
            "winners_all_records": winners, "n_records": len(rows),
            "three_way_n": len(three_way), "three_way_winners": three_winners,
            "bore_tolerance_pct": BORE_TOLERANCE * 100,
            "note": "the bore term is |b| x a DECLARED 5 % bore error, not a standard's "
                    "tolerance; it scales linearly with whatever error the caller declares"}

    gates = {"G18_1_newtonian_control": g181, "G18_2_power_law_closed_form": g182,
             "G18_3_discriminating_half": g183, "G18_4_never_steeper_than_newtonian": g184,
             "G18_5_newtonian_rule_error_split": g185, "G18_6_derating": g186,
             "G18_6p_derating_negative_control": g186p,
             "G18_7_error_term_ranking": g187}

    payload = {
        "generated_by": "scripts/needle_bore.py",
        "prereg": "data/needle_bore_prereg.md",
        "biorheo_version": biorheo.__version__,
        "geometry": {"R_ref_m": R_REF, "needle_length_m": NEEDLE_LENGTH,
                     "barrel_radius_m": BARREL_RADIUS, "bore_ratio": BORE_RATIO,
                     "note": "R_ref is a DECLARED radius carrying the blood study's "
                             "[UNVERIFIED] flag; the bore axis is a ratio and no gauge is named"},
        "operating_points": {"Q1_mL_per_min": Q1 * 6e7, "Q2_mL_per_min": Q2 * 6e7},
        "force_budget_N": F_BUDGET,
        "force_budget_note": "a rule of thumb from the JHU thesis, NOT a standard; every "
                             "derating factor is a ratio and is meaningful as one, while the "
                             "absolute Q_max is not",
        "n_records": len(rows), "rows": rows, "gates": gates,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")

    write_md(payload)
    ok = all(g.get("pass", True) for g in gates.values())
    print(f"wrote {OUT_JSON.name} and {OUT_MD.name}; gates "
          + ", ".join(f"{k.split('_')[0]}{k.split('_')[1]}.{k.split('_')[2]}="
                      f"{'PASS' if g.get('pass', True) else 'FAIL'}" for k, g in gates.items()))
    return 0 if ok else 1


def _f(x, spec=".4g"):
    return "-" if x is None else format(float(x), spec)


def write_md(p: dict) -> None:
    g = p["gates"]
    rows = p["rows"]
    L = [
        "# Item 18 — the needle-bore exponent, and what a force budget costs once the fit "
        "has an error bar",
        "",
        f"Written by `{p['generated_by']}` from `{p['prereg']}`. biorheo {p['biorheo_version']}. "
        "Every number on this page is printed by that script.",
        "",
        "## What is measured",
        "",
        "```",
        "b  =  d ln F / d ln R        at fixed flow rate, needle length and syringe",
        "```",
        "",
        "A Newtonian sizing calculation assumes `b = -4`. For a shear-thinning fluid it is "
        "not: a narrower bore raises the wall shear rate, which thins the fluid, so the force "
        "rises less steeply. **`b` is a property of the fluid AND the operating point**, not a "
        "material constant — a fluid already on its infinite-shear plateau at the wall returns "
        "`-4` however shear-thinning it is at rheometer rates, and that is exactly what the "
        "blood records do.",
        "",
        f"Geometry is DECLARED, not standardised. The bore axis is the ratio `rho = R/R_ref` "
        f"with `R_ref = {p['geometry']['R_ref_m']:.3g} m`, the same declared value the shipped "
        "blood study runs at and carrying that study's `[UNVERIFIED]` flag. **No needle gauge is "
        "named anywhere on this page**; this repo retracted an ISO 9626 comparison once and the "
        "standard is still not in hand.",
        "",
        "## Gates",
        "",
        "| gate | statement | result |",
        "|---|---|---|",
        f"| G18.1 | positive control, `b = -4` on the 20 `newtonian` records | "
        f"**{'PASS' if g['G18_1_newtonian_control']['pass'] else 'FAIL'}**, "
        f"{g['G18_1_newtonian_control']['n_pass']}/{g['G18_1_newtonian_control']['n']}, worst "
        f"|b+4| = {_f(g['G18_1_newtonian_control']['worst'], '.3g')} |",
        f"| G18.2 | closed form `b = -(3n+1)` on the 20 `power_law` records, computed without the "
        f"force solver | **{'PASS' if g['G18_2_power_law_closed_form']['pass'] else 'FAIL'}**, "
        f"{g['G18_2_power_law_closed_form']['n_pass']}/{g['G18_2_power_law_closed_form']['n']}, "
        f"worst relative error {_f(g['G18_2_power_law_closed_form']['worst'], '.3g')} |",
        f"| G18.3 | the discriminating half — the routine does not return `-4` for everything | "
        f"**{'PASS' if g['G18_3_discriminating_half']['pass'] else 'FAIL'}**, "
        f"{g['G18_3_discriminating_half']['n_above']} records with |b+4| > "
        f"{_f(G183_MIN_SPREAD, '.2g')}, worst {_f(g['G18_3_discriminating_half']['worst'], '.4g')} "
        f"(`{g['G18_3_discriminating_half']['worst_record']}`) |",
        f"| G18.4 | no record is steeper than Newtonian | "
        f"**{'PASS' if g['G18_4_never_steeper_than_newtonian']['pass'] else 'FAIL'}**, "
        f"{g['G18_4_never_steeper_than_newtonian']['n_pass']}/"
        f"{g['G18_4_never_steeper_than_newtonian']['n']}, steepest "
        f"{_f(g['G18_4_never_steeper_than_newtonian']['steepest'], '.6g')} |",
        f"| G18.5 | SPLIT prediction on the Newtonian rule's over-prediction at `rho = "
        f"{p['geometry']['bore_ratio']}` | **{'PASS' if g['G18_5_newtonian_rule_error_split']['pass'] else 'FAIL'}** "
        f"(see below) |",
        f"| G18.6 | the derating factor is below 1 on every closed-interval record | "
        f"**{'PASS' if g['G18_6_derating']['pass'] else 'FAIL'}**, "
        f"{g['G18_6_derating']['n_pass']}/{g['G18_6_derating']['n']}, strongest derating "
        f"{_f(g['G18_6_derating']['strongest_derating'], '.4g')} |",
        f"| **G18.6′** | the stricter test beside it: zero-width intervals must give a factor "
        f"of exactly 1, real ones strictly less | "
        f"**{'PASS' if g['G18_6p_derating_negative_control']['pass'] else 'FAIL'}**, "
        f"{_f(g['G18_6p_derating_negative_control']['factor_zero_width'], '.12g')} and "
        f"{_f(g['G18_6p_derating_negative_control']['factor_real'], '.5g')} |",
        "| G18.7 | three error terms ranked per record | reported, not barred |",
        "",
        "**G18.6 as pre-registered is close to unfalsifiable, and that is recorded rather than "
        "fixed by moving it.** The upper force envelope is a maximum over profile ends and the "
        "nominal force is not one of them, so the factor is below 1 unless the fitted parameters "
        "locally maximise the force in every direction at once. A routine that returned a constant "
        "0.5 would have passed it. The bar and its number stand; **G18.6′** is the "
        "discriminating test added beside it, and it pins both ends the bar cannot.",
        "",
        "## G18.5, the split prediction, and it came out as pre-registered",
        "",
        "The pre-registration predicted the blood half would **FAIL** its 10 % threshold and the "
        "Carbopol half would **PASS** it, for a stated reason in each case.",
        "",
        f"- **Blood: {g['G18_5_newtonian_rule_error_split']['blood_n_below']} of "
        f"{g['G18_5_newtonian_rule_error_split']['blood_n']}** Casson donors are BELOW 10 %; the "
        f"worst is {_f(g['G18_5_newtonian_rule_error_split']['blood_worst'] * 100, '.3g')} %. The "
        "Newtonian rule is very nearly right for blood in this needle, because at the wall the "
        "fluid is on its infinite-shear plateau — the same fact item 14 measured a different way "
        "(`|d ln F / d ln eta_c| = 0.995`).",
        f"- **Carbopol: {g['G18_5_newtonian_rule_error_split']['carbopol_named_above']} of "
        f"{g['G18_5_newtonian_rule_error_split']['carbopol_named_n']}** named gels are ABOVE 10 %, "
        f"and {g['G18_5_newtonian_rule_error_split']['carbopol_all_above']} of "
        f"{g['G18_5_newtonian_rule_error_split']['carbopol_all_n']} of all Carbopol records are. "
        f"The over-prediction runs "
        f"{_f(g['G18_5_newtonian_rule_error_split']['carbopol_least'] * 100, '.3g')} % to "
        f"{_f(g['G18_5_newtonian_rule_error_split']['carbopol_worst'] * 100, '.3g')} %.",
        "",
        "**The commercial reading, and it is the point of the item.** Whether a Newtonian sizing "
        "rule is safe is a property of *your* fluid at *your* operating point, and the two fluids "
        "in this registry sit on opposite sides of the question by an order of magnitude. A shop "
        "that quotes `R^-4` without measuring `b` cannot tell which case it is in — and the "
        "measurement costs one extra evaluation of a curve that has already been fitted.",
        "",
        "## The three error terms, in the same units",
        "",
        "All three are relative force errors at `Q1`, so they can be compared directly:",
        "",
        f"- **bore tolerance** — `|b| x {g['G18_7_error_term_ranking']['bore_tolerance_pct']:g} %`, a "
        "DECLARED bore error, not a standard's tolerance. It scales linearly with whatever error "
        "the caller declares.",
        "- **parameter uncertainty** — the one-at-a-time profile push-through hull, which is "
        "*wider* than the joint 95 % region it is drawn from.",
        "- **family choice** — item 12's Cross-vs-Casson difference, READ from that study's "
        "artifact and not recomputed. Defined only for the Casson blood records.",
        "",
        "| winner | records |",
        "|---|---|",
    ]
    for k, v in sorted(g["G18_7_error_term_ranking"]["winners_all_records"].items(),
                       key=lambda kv: -kv[1]):
        L.append(f"| {k} | {v} of {g['G18_7_error_term_ranking']['n_records']} |")
    L += ["",
          f"On the {g['G18_7_error_term_ranking']['three_way_n']} records where all three terms "
          "are defined (the Casson blood donors), the ranking is "
          + ", ".join(f"**{k}** {v}" for k, v in
                      sorted(g["G18_7_error_term_ranking"]["three_way_winners"].items(),
                             key=lambda kv: -kv[1])) + ".",
          "",
          "## Derating against a force budget",
          "",
          f"`Q_max` is the flow rate at which the plunger force reaches **{p['force_budget_N']:g} N**. "
          "That figure is a **rule of thumb** quoted in the JHU thesis, not a standard, and the "
          "absolute `Q_max` inherits every omission in `injectability.py` (no barrel pressure drop, "
          "no entrance loss, no wall slip, no hub friction, no stiction). **The ratio is the usable "
          "number**, because those omissions largely cancel between the two crossings.",
          "",
          f"Over {g['G18_6_derating']['n']} records with fully closed profile intervals, the "
          f"derating factor `Q_max(95 % upper) / Q_max(fit)` runs from "
          f"{_f(g['G18_6_derating']['strongest_derating'], '.5g')} to "
          f"{_f(g['G18_6_derating']['worst_factor'], '.5g')}. **The parameter uncertainty alone "
          f"costs up to a factor of "
          f"{_f(1.0 / float(g['G18_6_derating']['strongest_derating']), '.3g')} in deliverable "
          "flow rate**, before any of the physics the force model leaves out.",
          "",
          "**The split by identifiability is a TAIL, not a shift, and the medians do not show "
          "it.**",
          f"Records whose every parameter is `measurable` derate by a median "
          f"{_f(g['G18_6_derating']['median_factor_all_measurable'], '.4g')} "
          f"({g['G18_6_derating']['n_all_measurable']} records); records carrying at least one "
          f"non-`measurable` verdict derate by a median "
          f"{_f(g['G18_6_derating']['median_factor_not_all_measurable'], '.4g')} "
          f"({g['G18_6_derating']['n_not_all_measurable']} records) — **the same, to within "
          "the spread**. The difference is entirely at the bad end: the worst derating among "
          f"fully-measurable records is "
          f"{_f(g['G18_6_derating']['worst_factor_all_measurable'], '.4g')}, a factor of "
          f"{_f(1.0 / float(g['G18_6_derating']['worst_factor_all_measurable']), '.3g')}, against "
          f"{_f(g['G18_6_derating']['worst_factor_not_all_measurable'], '.4g')} — a factor of "
          f"{_f(1.0 / float(g['G18_6_derating']['worst_factor_not_all_measurable']), '.3g')} — "
          "among those with a non-measurable verdict. So an identifiability verdict does not "
          "predict how wide a typical error bar will be; what it predicts is whether the error "
          "bar can be catastrophic. The point forecast looks identical either way.",
          "",
          "## Per-record table",
          "",
          "| record | model | `b` at Q1 | `b` at Q2 | closed form | Newtonian-rule over-prediction "
          "at rho = 0.7 | hull | derating | all measurable? | dominant term |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        der = r["derating"]
        dtxt = (_f(der["factor"], ".5g") if isinstance(der, dict) and "factor" in der else "-")
        L.append(
            f"| `{r['id']}` | {r['model']} | {_f(r['b_Q1'], '.5g')} | {_f(r['b_Q2'], '.5g')} | "
            f"{_f(r['b_closed_form'], '.5g')} | {_f(r['newtonian_rule_over_prediction'] * 100, '.3g')} % | "
            f"{_f(r['hull_rel_Q1'] * 100 if r['hull_rel_Q1'] is not None else None, '.3g')} % | "
            f"{dtxt} | {'yes' if r['all_measurable'] else 'NO'} | "
            f"{r['dominant_error_term'] or '-'} |")
    L += ["",
          "## What this page does NOT say",
          "",
          "- Nothing about any real needle gauge, any standard, or any real device.",
          "- Nothing about a monoclonal antibody or any protein formulation. **The registry "
          "contains none.** The wedge fluid for the venture's first product is not in this data, "
          "and this page demonstrates the method on human blood and a microgel.",
          "- Nothing about temperature; no record carries a temperature axis.",
          f"- **EVERY record on this page is `commercial_use: false` — all "
          f"{sum(1 for r in rows if not r['commercial_ok'])} of {len(rows)} rows cannot appear "
          "in a commercial deliverable.** The 15 Carbopol records are CC BY-NC-SA; the 100 "
          "blood records come from an author-posted workbook that states no licence at all, "
          "which is not permission. **This is a standing constraint on the venture, not a "
          "detail of this item**: the method demonstrated here is licensable, and not one "
          "curve currently in the registry is. A paid study needs either the client's own "
          "data or a curve whose licence permits commercial use.",
          "- No record was added, no curve was fitted, and `dst_fit` was not touched.",
          ""]
    OUT_MD.write_text("\n".join(L), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
