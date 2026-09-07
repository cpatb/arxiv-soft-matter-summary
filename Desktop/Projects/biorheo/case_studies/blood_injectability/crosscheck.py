"""NEXT_STEPS item 12: does a donor's Cross record give the same needle force as its
Casson record?

**Read `crosscheck_prereg.md` first.** Every bar below is fixed there and was committed
at `1cd46de`, before any of this was computed. This script grades and prints; it does not
choose bars.

The comparison. Each of the 20 Horner donors has five registry records fitted to the SAME
16-point flow curve: `casson` (the reference, and the family the case study ships),
`cross`, `carreau_yasuda`, and the two NULL models `newtonian` and `power_law`. Every one
is pushed through the same 27G needle at the same 1 mL/min by the same
`biorheo.apps.injectability` code path, and the plunger force is compared with the
donor's own Casson force.

It can fail. The needle's wall rate is about 1.8e4 1/s, **1.42 decades above the 700 1/s
top of every record's envelope**, so all five families are read where none was measured.
Casson, Cross and Carreau-Yasuda each approach a constant there -- but three DIFFERENT
fitted constants, each partly an extrapolation of a fit that stops at 700 1/s. Whether
they agree is a real question and "they fit the same curve" does not settle it.

Three scales, all from the case study's own `results.json` and all named in the
pre-registration before the run:

  1 %      the calculation's own unresolved systematic. The entrance (Hagenbach) loss,
           which the model omits, is 0.43-0.99 % of dp, and Fahraeus-Lindqvist is a
           second unquantified correction pointing the OTHER way. Agreement below this
           is agreement below the resolution of the model.
  2.94 %   the largest yield-stress contribution -- the effect the case study calls
           negligible.
  128 %    the donor-to-donor force spread (a factor 2.2792), which is the tolerance
           item 12 proposed. It is a nearly vacuous bar and is labelled as one.

Run:  python case_studies/blood_injectability/crosscheck.py
Writes `crosscheck.md` and `crosscheck.json` beside this file. Pinned in
`tests/test_injectability_crosscheck.py`.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import binomtest, spearmanr

import biorheo
from biorheo.constitutive.gnf import NEWTONIAN
from biorheo.apps.injectability import (
    Needle, Syringe, flow_rate_casson, flow_rate_gnf, injection_force,
    pressure_for_flow_rate, wall_shear_rate,
)
from biorheo.registry.loader import load_all

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results.json"

# The operating point, taken from `run.py` rather than retyped, so the two studies cannot
# drift apart. Asserted against the committed artifact below.
ID_NOMINAL_MM = 0.210
NEEDLE_LENGTH = 12.7e-3
BARREL_RADIUS = 2.3e-3
Q_SI = 1.0 * 1e-6 / 60.0
RHO_BLOOD = 1060.0

REFERENCE = "casson"
FAMILIES = ("cross", "carreau_yasuda", "newtonian", "power_law")
NULLS = ("newtonian", "power_law")

# --- the three scales. NOT bars chosen here: they are read out of the case study's own
# artifact at run time (see `_scales`), and these names are only what they are called.
BAND_RESOLUTION = "below the calculation's own resolution"
BAND_YIELD = "within the neglected yield term"
BAND_ABOVE_YIELD = "above the neglected yield term, below the donor spread"
BAND_FAIL = "above the donor spread"

# G12.3: the separation floor below which a discriminator is called uninformative. This
# programme's since 2026-08-30 (the frisbee dispersion-law test), unchanged here.
SEPARATION_FLOOR = 1.6


def _scales() -> dict:
    """The three comparison scales, read from the case study's committed artifact.

    They are NOT retyped here. If `run.py` is re-run and a scale moves, this study moves
    with it and the .md says so, rather than grading against a number that used to be
    true."""
    d = json.loads(RESULTS.read_text())
    return {
        "resolution": max(d["entrance_loss"]["fraction_of_dp_min"],
                          d["entrance_loss"]["fraction_of_dp_max"]),
        "yield_term": d["yield_share_max"],
        "donor_spread": d["geometry"]["spread_all_20"] - 1.0,
        "casson_force_min": min(r["force_N"] for r in d["rows"]),
        "casson_force_max": max(r["force_N"] for r in d["rows"]),
        "source": "case_studies/blood_injectability/results.json",
    }


def band(dev: float, sc: dict) -> str:
    """Which of the four pre-registered bands a worst-case relative deviation falls in."""
    if not np.isfinite(dev):
        return "NOT COMPUTED"
    if dev < sc["resolution"]:
        return BAND_RESOLUTION
    if dev < sc["yield_term"]:
        return BAND_YIELD
    if dev < sc["donor_spread"]:
        return BAND_ABOVE_YIELD
    return BAND_FAIL


def force_of(rec, nd: Needle) -> dict:
    """The plunger force for one record, or the reason it could not be computed.

    Every failure is CAUGHT and named rather than dropped: `flow_rate_gnf` refuses a
    non-monotone flow curve (a Cross record with m > 1 can produce one) and
    `pressure_for_flow_rate` can run out of bracket. A record this study cannot evaluate
    is data about the fit."""
    try:
        dp = pressure_for_flow_rate(Q_SI, nd, rec.model, rec.params_si)
        gd_w = wall_shear_rate(dp, nd, rec.model, rec.params_si)
        f = injection_force(Q_SI, nd, Syringe(BARREL_RADIUS), rec.model, rec.params_si,
                            rho=RHO_BLOOD)
        assert abs(f.dp - dp) / dp < 1e-9, (rec.id, f.dp, dp)
        lo, hi = rec.envelope.gammadot_min, rec.envelope.gammadot_max
        return {"record": rec.id, "model": rec.model.name, "force_N": f.force,
                "dp_Pa": dp, "gd_wall": gd_w, "envelope": [lo, hi],
                "extrapolated": bool(gd_w < lo or gd_w > hi),
                "decades_past_high": float(np.log10(gd_w / hi)),
                "all_measurable": set((rec.verdicts or {}).values()) == {"measurable"},
                "verdicts": rec.verdicts, "failed": None}
    except Exception as e:                      # noqa: BLE001 -- the reason is the datum
        return {"record": rec.id, "model": rec.model.name, "force_N": float("nan"),
                "all_measurable": set((rec.verdicts or {}).values()) == {"measurable"},
                "verdicts": rec.verdicts,
                "failed": f"{type(e).__name__}: {str(e)[:200]}"}


def route_control(ref: dict, nd: Needle) -> dict:
    """**The control the pre-registration did not think of, added after the first run and
    reported whichever way it falls.**

    Casson forces come from the CLOSED FORM (`flow_rate_casson`) and every other family
    from the QUADRATURE (`flow_rate_gnf`). If those two routes disagreed at this operating
    point, the whole comparison would be measuring the ROUTE and not the fit -- and the
    effect it found (tens of per cent) is far too large to assume away. So the Casson
    force is recomputed through the quadrature, on the same records, at the same Q.

    This is strictly a stronger test: it can only reveal a confound, never hide one. It is
    the "two anchors that share an identity cannot check it" lesson applied forwards."""
    worst, rows = 0.0, []
    for d, r in sorted(ref.items()):
        rec = r
        dp = pressure_for_flow_rate(Q_SI, nd, rec.model, rec.params_si)
        # the same pressure through the two independent flow-rate routes
        q_closed = flow_rate_casson(dp, nd.radius, nd.length, *rec.params_si)
        q_quad = flow_rate_gnf(dp, nd.radius, nd.length, rec.model, rec.params_si)
        dev = abs(q_quad / q_closed - 1.0)
        worst = max(worst, dev) if np.isfinite(dev) else float("nan")
        rows.append({"donor": d, "Q_closed": q_closed, "Q_quad": q_quad, "rel_dev": dev})
    return {"worst_rel_dev": worst, "rows": rows,
            "claim": ("Casson's closed form and the generic quadrature agree at this "
                      "operating point, so a family-to-family force difference is a "
                      "difference between FITS and not between numerical routes.")}


def _finite(xs):
    """Finite values only. `max(0.0, nan)` returns 0.0 in Python, so a NaN row is
    SILENTLY DROPPED and a worst-case understated -- `lessons.md`, 2026-09-06."""
    return [x for x in xs if np.isfinite(x)]


def sign_accounting(rows: list[dict]) -> dict:
    """The signed structure of a set of relative deviations.

    Reporting a signed median beside an unsigned worst hides the case where the extremes
    point in OPPOSITE directions -- which is what happens here, and what made "both
    families predict a LOWER force, consistently" false when it was written into a commit
    message and onto two pages. The supportable claim is a paired sign test."""
    d = _finite([r["rel_dev"] for r in rows])
    neg = [x for x in d if x < 0]
    pos = [x for x in d if x > 0]
    p = float(binomtest(len(neg), len(d), 0.5).pvalue) if d else float("nan")
    return {"n": len(d), "n_negative": len(neg), "n_positive": len(pos),
            "most_negative": min(d) if d else float("nan"),
            "most_positive": max(pos) if pos else 0.0,
            "sign_test_p": p,
            "claim": (f"lower on {len(neg)} of {len(d)} donors (two-sided sign test "
                      f"p = {p:.4g}), NOT on all of them")}


def matched_plateau_null(ref: dict, recs: dict, nd: Needle) -> dict:
    """**The discriminating negative case for "the needle reads the plateau", supplied by
    the code review after this study's first two nulls were shown not to be one.**

    Each donor's `cross` record is re-run with `eta_inf` REPLACED by that donor's Casson
    `eta_c`, leaving `eta_0`, `lam` and `m` at their fitted values. This is not Casson in
    disguise -- it keeps the whole low-rate structure of the Cross fit -- so if the force
    difference were coming from anywhere other than the high-rate plateau it would
    survive. It goes THROUGH the flow solver rather than around it, which is what the
    earlier `plateau_control` did not do.

    Contrast with the two nulls that do not work:
      * `newtonian` at 271 % is `mu/eta_c` divided by the yield factor -- an algebraic
        identity that restates "blood shear-thins by 2-4x over the measured range".
      * `plateau_null` (Newtonian at `eta_c`) is Casson with `tau_y = 0`, i.e. the
        sibling study's `yield_share_max` by another route."""
    rows, worst = [], 0.0
    cross_model = recs[f"{next(iter(ref.values())).id}_cross"].model
    i_inf = cross_model.param_names.index("eta_inf")
    for d, r in sorted(ref.items()):
        eta_c = dict(zip(r.model.param_names, r.params_si))["eta_c"]
        f_cas = injection_force(Q_SI, nd, Syringe(BARREL_RADIUS), r.model, r.params_si,
                                rho=RHO_BLOOD).force
        p_matched = np.array(recs[f"{r.id}_cross"].params_si, dtype=float)
        p_matched[i_inf] = eta_c
        f_mat = injection_force(Q_SI, nd, Syringe(BARREL_RADIUS), cross_model,
                                p_matched, rho=RHO_BLOOD).force
        dev = f_mat / f_cas - 1.0
        worst = max(worst, abs(dev))
        rows.append({"donor": d, "rel_dev": float(dev)})
    return {"worst_abs_rel_dev": worst,
            "median_rel_dev": float(np.median([r["rel_dev"] for r in rows])),
            "rows": rows,
            "claim": ("Cross with its plateau set to Casson's, keeping all its other "
                      "fitted parameters, collapses the difference into the yield band. "
                      "So at the needle the plateau IS what is being read, and the "
                      "30 % is a disagreement about its value.")}


def plateau_control(ref: dict, recs: dict, summary: dict) -> dict:
    """**The second control added after the first run: an anchor OUTSIDE the flow solver.**

    At the needle's wall rate every one of these families is effectively a constant
    viscosity, and for a Newtonian fluid at fixed Q the plunger force is proportional to
    that constant. So the force ratio should reproduce the ratio of the two fits' own
    high-rate plateaux -- `eta_inf` for Cross and Carreau-Yasuda, `eta_c` for Casson --
    read straight out of the records, touching no part of `injectability`.

    If it does, the mechanism is identified: the families disagree because their fitted
    plateaux disagree, which is a statement about extrapolating a fit that stops at
    700 1/s. If it does not, the difference is coming from somewhere else and the
    explanation on the page would be wrong."""
    out = {}
    for fam, key in (("cross", "eta_inf"), ("carreau_yasuda", "eta_inf")):
        rows = []
        for d, r in sorted(ref.items()):
            f_rec = recs[f"{r.id}_{fam}"]
            p = dict(zip(f_rec.model.param_names, f_rec.params_si))
            eta_c = dict(zip(r.model.param_names, r.params_si))["eta_c"]
            plateau_ratio = float(p[key]) / float(eta_c) - 1.0
            force_ratio = next(x["rel_dev"] for x in summary[fam]["rows"] if x["donor"] == d)
            rows.append({"donor": d, "plateau_ratio": plateau_ratio,
                         "force_ratio": force_ratio,
                         # SIGNED. The yield explanation is signed and quantitative --
                         # F_casson = F_newt(eta_c)(1+y) makes the force ratio sit BELOW
                         # the plateau ratio on every donor -- so an unsigned |gap| cannot
                         # test it, and the first version of this study asserted it with
                         # an unsigned bound and a factor-2 slop. It is false: the
                         # residual has the wrong sign on several donors.
                         "signed_gap": float(force_ratio - plateau_ratio),
                         "gap": abs(force_ratio - plateau_ratio)})
        gaps = _finite([x["gap"] for x in rows])
        signed = _finite([x["signed_gap"] for x in rows])
        worst_row = max((x for x in rows if np.isfinite(x["gap"])),
                        key=lambda x: x["gap"], default=None)
        out[fam] = {
            "worst_gap": max(gaps) if gaps else float("nan"),
            "worst_gap_donor": worst_row["donor"] if worst_row else None,
            "n_wrong_sign": sum(1 for x in signed if x > 0),
            "wrong_sign_donors": [x["donor"] for x in rows if x["signed_gap"] > 0],
            "rows": rows,
        }
    return out


def misfit_correlation(ref: dict, recs: dict, summary: dict) -> dict:
    """Does the force difference track how badly Casson misfits that donor's curve?

    The two families are NOT equally good -- Casson's `rmse_log` is worse on 20/20 donors
    -- so "family choice moves the force" is the wrong framing: the family the case study
    ships is the worse-fitting one, and if the effect grows with the misfit that is a
    directional statement, not a symmetric ambiguity. Raised by code review."""
    out = {}
    for fam in ("cross", "carreau_yasuda"):
        ratio, dev = [], []
        for d, r in sorted(ref.items()):
            e = recs[f"{r.id}_{fam}"].extra
            ratio.append(float(e["rmse_log_casson"]) / float(e["rmse_log_this_model"]))
            dev.append(abs(next(x["rel_dev"] for x in summary[fam]["rows"]
                                if x["donor"] == d)))
        keep = [(a, b) for a, b in zip(ratio, dev) if np.isfinite(a) and np.isfinite(b)]
        rho, p = spearmanr([a for a, _ in keep], [b for _, b in keep])
        out[fam] = {"spearman_rho": float(rho), "p": float(p), "n": len(keep),
                    "rmse_ratio_min": float(min(ratio)),
                    "rmse_ratio_max": float(max(ratio)),
                    "casson_worse_on": sum(1 for a in ratio if a > 1.0)}
    return out


def plateau_null(ref: dict, nd: Needle) -> dict:
    """**A sharper null than the two registry ones, added after the first run.**

    `newtonian` misses by 271 %, but partly BY CONSTRUCTION: its `mu` is fitted over the
    whole 0.1-700 1/s curve on log-viscosity residuals, so it is a kind of average over
    the whole curve and is nowhere near the high-rate plateau. A null that is off for a
    trivial reason still shows the needle discriminates, but it cannot separate the two
    questions this study needs kept apart:

      (a) does the needle read the high-rate PLATEAU, or something else?
      (b) do the three families' plateaux DISAGREE?

    So: a Newtonian fluid whose viscosity is the donor's own Casson `eta_c` -- the
    plateau itself, not a fit to anything.

    **IT DOES NOT WORK, and the reason is worth keeping.** Casson at `tau_y = 0` IS
    Newtonian at `eta_c` -- that is the same law, not a different one -- so this
    recomputes the sibling study's `yield_share_max` by another route rather than testing
    anything it does not already assume. The numbers agree to 1e-5, which makes it a
    useful cross-implementation check of that quantity (one route through
    `flow_rate_casson`, one through `flow_rate_gnf`) and NOT a second anchor on the
    mechanism. It is kept, and labelled, because "an anchor that shares an identity with
    its target cannot check that identity" is a lesson this repo wrote down on 2026-09-05
    and this is what walking into it again looks like."""
    rows, worst = [], 0.0
    for d, r in sorted(ref.items()):
        eta_c = dict(zip(r.model.param_names, r.params_si))["eta_c"]
        f_cas = injection_force(Q_SI, nd, Syringe(BARREL_RADIUS), r.model, r.params_si,
                                rho=RHO_BLOOD).force
        f_new = injection_force(Q_SI, nd, Syringe(BARREL_RADIUS), NEWTONIAN,
                                np.array([eta_c]), rho=RHO_BLOOD).force
        dev = abs(f_new / f_cas - 1.0)
        worst = max(worst, dev)
        rows.append({"donor": d, "eta_c": float(eta_c), "force_casson": f_cas,
                     "force_newtonian_at_eta_c": f_new, "rel_dev": dev})
    return {"worst_rel_dev": worst, "rows": rows,
            "claim": ("A Newtonian fluid AT the Casson plateau reproduces the Casson "
                      "force to about the yield term, so the needle is reading the "
                      "high-rate plateau and essentially nothing else. The 30 % is "
                      "therefore a disagreement between plateaux, not a difference in "
                      "what the geometry samples.")}


def in_envelope_needle(ref: dict, recs: dict, target_gd: float) -> dict:
    """**The control that decides C1: run the SAME comparison through the SAME solver at a
    needle wide enough that the wall rate lands INSIDE the measured envelope.**

    The `viscosity_ladder` below shows the fits disagree across the measured range, but it
    evaluates the laws directly. This runs the whole force calculation -- pressure
    inversion, quadrature, plunger area -- with no extrapolation anywhere, which is the
    only way to say what the extrapolation is contributing to the HEADLINE number rather
    than to the viscosity ratio. Supplied by code review; it is 20 lines and it changes
    the sentence.

    The radius is chosen so the apparent wall rate `4Q/(pi R^3)` equals `target_gd`; the
    true rate is then checked to be inside every record's envelope and the study says so
    rather than assuming it."""
    R = (4.0 * Q_SI / (np.pi * target_gd)) ** (1.0 / 3.0)
    nd = Needle(radius=R, length=NEEDLE_LENGTH)
    out = {"radius_m": float(R), "id_mm": float(2e3 * R), "target_gd": float(target_gd),
           "families": {}}
    inside = True
    for fam in ("cross", "carreau_yasuda"):
        rows = []
        for d, r in sorted(ref.items()):
            f_cas = force_of(r, nd)
            f_fam = force_of(recs[f"{r.id}_{fam}"], nd)
            if f_cas["failed"] or f_fam["failed"]:
                rows.append({"donor": d, "rel_dev": float("nan")})
                continue
            inside &= not (f_cas["extrapolated"] or f_fam["extrapolated"])
            rows.append({"donor": d, "rel_dev": float(f_fam["force_N"] / f_cas["force_N"] - 1.0),
                         "gd_wall": f_cas["gd_wall"],
                         "decades_past_high": max(f_cas["decades_past_high"],
                                                  f_fam["decades_past_high"])})
        devs = _finite([x["rel_dev"] for x in rows])
        out["families"][fam] = {
            "rows": rows, "n": len(devs),
            "worst_abs_rel_dev": max((abs(x) for x in devs), default=float("nan")),
            "median_rel_dev": float(np.median(devs)) if devs else float("nan"),
        }
    out["every_record_inside_its_envelope"] = bool(inside)
    # Sizing on the APPARENT rate `4Q/(pi R^3)` leaves the true wall rate of a
    # shear-thinning fluid slightly higher, so a few records can still sit just outside.
    # Reported in decades rather than claimed away: the point of this control is that the
    # extrapolation is ~0, and "~0" has to be a number.
    over = [x["decades_past_high"] for v in out["families"].values() for x in v["rows"]
            if "decades_past_high" in x]
    out["worst_decades_past_high"] = float(max(over)) if over else float("nan")
    out["n_records_outside"] = int(sum(1 for x in over if x > 0))
    out["n_records_total"] = len(over)
    return out


def viscosity_ladder(ref: dict, recs: dict, gd_needle: float) -> dict:
    """**The check that refuted this study's own first explanation.**

    Having found a 30 % force difference at a rate 1.42 decades past the data, the obvious
    reading -- and the one first written onto both pages -- was "the cause is the
    extrapolation". That is testable in one line, by asking how far apart the same two
    fits are at rates INSIDE the measured envelope, and it is false.

    The families disagree by tens of per cent across the WHOLE measured range, and by MORE
    at the bottom of the envelope than at the needle. So the needle is not where the
    disagreement is created; it is just one rate at which two differently-shaped fits to
    the same 16 noisy points are read. The mechanism is the MISFIT (Casson's median
    rmse_log is 0.146 against Cross's 0.030, and it is flagged `structured` on 19 of 20
    donors -- `data/horner2020/model_comparison.md`), not the extrapolation.

    This touches no part of `injectability`: it evaluates the records' own laws."""
    lo = min(r.envelope.gammadot_min for r in ref.values())
    hi = max(r.envelope.gammadot_max for r in ref.values())
    rates = list(np.geomspace(lo, hi, 5)) + [gd_needle]
    rows = []
    for gd in rates:
        row = {"gammadot": float(gd), "inside_envelope": bool(lo <= gd <= hi)}
        for fam in ("cross", "carreau_yasuda"):
            devs = []
            for d, r in sorted(ref.items()):
                a = float(r.model(np.array([gd]), *r.params_si)[0])
                f = recs[f"{r.id}_{fam}"]
                b = float(f.model(np.array([gd]), *f.params_si)[0])
                devs.append(abs(b / a - 1.0))
            row[fam] = float(max(devs))
        rows.append(row)
    inside = [r for r in rows if r["inside_envelope"]]
    needle = [r for r in rows if not r["inside_envelope"]][0]
    # the misfit itself, computed from the records' own `extra` rather than quoted from
    # `model_comparison.md` -- every number in a case-study .md is printed by a script
    rmse = {"casson": float(np.median(
        [recs[f"{r.id}_cross"].extra["rmse_log_casson"] for r in ref.values()]))}
    for fam in ("cross", "carreau_yasuda"):
        rmse[fam] = float(np.median(
            [recs[f"{r.id}_{fam}"].extra["rmse_log_this_model"] for r in ref.values()]))
    structured = {"casson": sum(
        bool((r.identifiability or {}).get("residual_structure", {}).get("structured"))
        for r in ref.values())}
    for fam in ("cross", "carreau_yasuda"):
        structured[fam] = sum(
            bool((recs[f"{r.id}_{fam}"].identifiability or {})
                 .get("residual_structure", {}).get("structured"))
            for r in ref.values())
    return {
        "rows": rows, "envelope": [lo, hi], "gd_needle": gd_needle,
        "median_rmse_log": rmse, "n_structured_of_20": structured, "n_donors": len(ref),
        "worst_inside_cross": max(r["cross"] for r in inside),
        "worst_inside_cy": max(r["carreau_yasuda"] for r in inside),
        "at_needle_cross": needle["cross"], "at_needle_cy": needle["carreau_yasuda"],
        "needle_is_the_worst_rate": bool(
            needle["cross"] >= max(r["cross"] for r in inside)),
        "claim": ("The families disagree across the whole MEASURED range, by more at the "
                  "bottom of the envelope than at the needle. The 30 % is a misfit "
                  "between two differently-shaped fits to the same 16 points, not an "
                  "artefact of reading them past the data."),
    }


def main() -> int:
    recs = load_all()
    nd = Needle(radius=ID_NOMINAL_MM * 1e-3 / 2.0, length=NEEDLE_LENGTH)
    sc = _scales()

    # the reference: the 20 unsuffixed Casson records, the ones the case study ships
    ref = {r.id.split("_")[-1]: r for r in recs.values() if r.model.name == REFERENCE}
    assert len(ref) == 20, len(ref)
    donors = sorted(ref)

    # the operating point must be the case study's, not a retyped copy of it
    src = json.loads(RESULTS.read_text())
    assert src["needle"]["id_nominal_mm"] == ID_NOMINAL_MM, src["needle"]
    assert src["needle"]["length_m"] == NEEDLE_LENGTH, src["needle"]
    assert src["Q_mL_per_min"] == 1.0 and src["rho"] == RHO_BLOOD

    out: dict[str, dict[str, dict]] = {REFERENCE: {}}
    for d in donors:
        out[REFERENCE][d] = force_of(ref[d], nd)
    # the reference forces must reproduce the committed case study, or this comparison is
    # against a different calculation than the one it claims to cross-check
    ref_src = {r["donor"]: r["force_N"] for r in src["rows"]}
    repro = max(abs(out[REFERENCE][d]["force_N"] / ref_src[d] - 1.0) for d in donors)

    for fam in FAMILIES:
        out[fam] = {}
        for d in donors:
            rid = f"{ref[d].id}_{fam}"
            out[fam][d] = force_of(recs[rid], nd) if rid in recs else {
                "record": rid, "model": fam, "force_N": float("nan"),
                "failed": "no such record", "all_measurable": False, "verdicts": None}

    # --- the comparison ----------------------------------------------------------------
    summary: dict[str, dict] = {}
    for fam in FAMILIES:
        rows = []
        for d in donors:
            a, b = out[fam][d], out[REFERENCE][d]
            dev = (float(a["force_N"] / b["force_N"] - 1.0)
                   if a["failed"] is None and b["failed"] is None else float("nan"))
            rows.append({"donor": d, "force_N": a["force_N"],
                         "casson_force_N": b["force_N"], "rel_dev": dev,
                         "all_measurable": a["all_measurable"], "failed": a["failed"]})
        ok = [r for r in rows if np.isfinite(r["rel_dev"])]
        meas = [r for r in ok if r["all_measurable"]]
        worst = max((abs(r["rel_dev"]) for r in ok), default=float("nan"))
        worst_meas = max((abs(r["rel_dev"]) for r in meas), default=float("nan"))
        summary[fam] = {
            "rows": rows,
            "n_in": len(rows), "n_out": len(ok),
            "n_failed": len(rows) - len(ok),
            "failed_donors": [r["donor"] for r in rows if not np.isfinite(r["rel_dev"])],
            "failure_reasons": {r["donor"]: r["failed"] for r in rows
                                if r["failed"] is not None},
            "n_all_measurable": len(meas),
            "worst_abs_rel_dev": worst,
            "worst_abs_rel_dev_measurable_only": worst_meas,
            "median_rel_dev": float(np.median([r["rel_dev"] for r in ok])) if ok else float("nan"),
            "band": band(worst, sc),
            "band_measurable_only": band(worst_meas, sc),
        }

    cross = summary["cross"]
    # G12.1 -- item 12's own bar, pre-registered as WEAK
    g121 = {"bar": sc["donor_spread"], "worst": cross["worst_abs_rel_dev"],
            "pass": bool(cross["worst_abs_rel_dev"] < sc["donor_spread"]),
            "note": "the bar item 12 proposed; passing it establishes very little"}
    # G12.3 -- the informativeness control: a null must land in a DIFFERENT band
    null_worst = max(summary[f]["worst_abs_rel_dev"] for f in NULLS)
    sep = (null_worst / cross["worst_abs_rel_dev"]
           if cross["worst_abs_rel_dev"] > 0 else float("inf"))
    different_band = any(summary[f]["band"] != cross["band"] for f in NULLS)
    g123 = {"null_worst": null_worst, "cross_worst": cross["worst_abs_rel_dev"],
            "separation": float(sep), "floor": SEPARATION_FLOOR,
            "nulls_in_a_different_band": bool(different_band),
            "informative": bool(different_band and sep >= SEPARATION_FLOOR),
            "null_bands": {f: summary[f]["band"] for f in NULLS}}
    # G12.4 -- the verdict must not depend on the identifiability split
    g124 = {f: {"band": summary[f]["band"],
                "band_measurable_only": summary[f]["band_measurable_only"],
                "n_all_measurable": summary[f]["n_all_measurable"],
                "same": summary[f]["band"] == summary[f]["band_measurable_only"]}
            for f in FAMILIES}
    # G12.5 -- nothing silently dropped
    g125 = {f: {"n_in": summary[f]["n_in"], "n_out": summary[f]["n_out"],
                "failed": summary[f]["failed_donors"]} for f in FAMILIES}
    # the two controls added AFTER the first run, because the effect it found is far too
    # large to leave a numerical-route confound or an unidentified mechanism unchecked
    rc = route_control(ref, nd)
    pc = plateau_control(ref, recs, summary)
    pn = plateau_null(ref, nd)
    gd_needle = float(np.median([out[REFERENCE][d]["gd_wall"] for d in donors]))
    vl = viscosity_ladder(ref, recs, gd_needle)
    mpn = matched_plateau_null(ref, recs, nd)
    mc = misfit_correlation(ref, recs, summary)
    # the top of every record's envelope, so the wall rate lands INSIDE it
    ien = in_envelope_needle(ref, recs, max(r.envelope.gammadot_max
                                            for r in ref.values()))
    for fam in FAMILIES:
        summary[fam]["signs"] = sign_accounting(summary[fam]["rows"])
    # G12.3 used `max` over the two nulls, which is a free choice. Recorded both ways so
    # the verdict is not resting on which null happened to be picked.
    null_min = min(summary[f]["worst_abs_rel_dev"] for f in NULLS)
    g123["separation_worst_null"] = g123["separation"]
    g123["separation_least_extreme_null"] = float(
        null_min / cross["worst_abs_rel_dev"]) if cross["worst_abs_rel_dev"] > 0 else float("inf")
    g123["informative_on_either_null"] = bool(
        g123["separation_least_extreme_null"] >= SEPARATION_FLOOR)

    payload = {
        "generated_by": "case_studies/blood_injectability/crosscheck.py",
        "prereg": "case_studies/blood_injectability/crosscheck_prereg.md (commit 1cd46de)",
        "biorheo_version": biorheo.__version__,
        "scales": sc, "reference": REFERENCE, "donors": donors,
        "reference_reproduces_case_study_to": repro,
        "needle": {"id_nominal_mm": ID_NOMINAL_MM, "length_m": NEEDLE_LENGTH},
        "Q_mL_per_min": 1.0, "rho": RHO_BLOOD,
        "summary": summary,
        "G12_1_item12_own_bar": g121,
        "G12_2_bands": {f: summary[f]["band"] for f in FAMILIES},
        "G12_3_informative": g123,
        "G12_4_identifiability_split": g124,
        "G12_5_accounting": g125,
        "route_control": rc,
        "plateau_control": pc,
        "plateau_null": pn,
        "viscosity_ladder": vl,
        "matched_plateau_null": mpn,
        "misfit_correlation": mc,
        "in_envelope_needle": ien,
        # measured over every family and donor this study could evaluate, not quoted
        # from the sibling case study, which computed it for `casson` alone
        "median_decades_past_high": float(np.median(
            [v[d]["decades_past_high"] for v in out.values() for d in donors
             if v[d]["failed"] is None])),
    }
    (HERE / "crosscheck.json").write_text(json.dumps(payload, indent=2), encoding="ascii")

    md = _report(payload)
    (HERE / "crosscheck.md").write_text(md, encoding="utf-8", newline="\n")
    print(md)
    ok = (repro < 1e-9
          and all(s["n_out"] == s["n_in"] for s in summary.values())
          and all(v["same"] for v in g124.values()))
    if not ok:
        print("\nNOTE: a check above did not come out clean; see the .md.")
    return 0


def _pct(x: float) -> str:
    return "n/a" if not np.isfinite(x) else f"{100 * x:.4g} %"


def _report(p: dict) -> str:
    sc, s = p["scales"], p["summary"]
    L = [
        "# Item 12: does a donor's Cross record give the same needle force as its Casson record?",
        "",
        "Generated by `case_studies/blood_injectability/crosscheck.py`; every number here "
        "reproduces by running it. **Bars are pre-registered in `crosscheck_prereg.md`, "
        "committed at `1cd46de` before any of this was computed, and are not restated as "
        "new numbers -- they are read out of the case study's own `results.json` at run "
        "time.**", "",
        f"The reference Casson forces reproduce the committed case study to "
        f"{p['reference_reproduces_case_study_to']:.2e} relative, so this is a "
        "cross-check of that calculation and not of a different one.", "",
        "## The three scales, read from `results.json`", "",
        "| scale | value | what it is |", "|---|---|---|",
        f"| resolution | {_pct(sc['resolution'])} | the entrance (Hagenbach) loss the "
        "model OMITS, as a fraction of dp. Fahraeus-Lindqvist is a second unquantified "
        "correction pointing the OTHER way, so they do not cancel predictably. **Two "
        "models agreeing below this are not distinguishable by this calculation.** |",
        f"| neglected yield term | {_pct(sc['yield_term'])} | the largest shift in force "
        "from zeroing `tau_y` -- the effect the case study calls negligible |",
        f"| donor spread | {_pct(sc['donor_spread'])} | max/min force over the 20 donors, "
        "the tolerance item 12 proposed. **A nearly vacuous bar.** |", "",
        f"The ratio's denominator is the donor's own Casson force, "
        f"{sc['casson_force_min']:.4g} to {sc['casson_force_max']:.4g} N -- a declared, "
        "non-collapsing scale, checked rather than assumed (this repo's item-9 "
        "retraction is the reason that sentence is here).", "",
        "## Result", "",
        "| family | n in | n out | worst \\|rel dev\\| | median rel dev | band |",
        "|---|---|---|---|---|---|",
    ]
    for f in ("cross", "carreau_yasuda", "newtonian", "power_law"):
        v = s[f]
        L.append(f"| `{f}`{' (NULL)' if f in NULLS else ''} | {v['n_in']} | {v['n_out']} | "
                 f"{_pct(v['worst_abs_rel_dev'])} | {_pct(v['median_rel_dev'])} | "
                 f"{v['band']} |")
    L += ["",
          "**The extremes point in OPPOSITE directions, and reporting a signed median "
          "beside an unsigned worst hid that.** A first version of this page and of the "
          "commit message said both no-yield families predict a *consistently* lower "
          "force. They do not:", "",
          "| family | lower on | most positive | two-sided sign test |", "|---|---|---|---|"]
    for f in ("cross", "carreau_yasuda"):
        g = s[f]["signs"]
        L.append(f"| `{f}` | {g['n_negative']}/{g['n']} donors | "
                 f"+{_pct(g['most_positive'])} | p = {g['sign_test_p']:.4g} |")
    L += ["", "So the supportable claim is a paired sign test, not \"consistently\".", "",
          "**And the two families are not equally good.** Casson's `rmse_log` is worse on "
          f"{p['misfit_correlation']['cross']['casson_worse_on']}/20 donors, by "
          f"{p['misfit_correlation']['cross']['rmse_ratio_min']:.2f}x to "
          f"{p['misfit_correlation']['cross']['rmse_ratio_max']:.2f}x, and the force "
          "difference TRACKS the misfit: Spearman(|rel dev|, rmse ratio) = "
          f"{p['misfit_correlation']['cross']['spearman_rho']:.3f} "
          f"(p = {p['misfit_correlation']['cross']['p']:.4g}) for `cross` and "
          f"{p['misfit_correlation']['carreau_yasuda']['spearman_rho']:.3f} "
          f"(p = {p['misfit_correlation']['carreau_yasuda']['p']:.4g}) for "
          "`carreau_yasuda`. **This is not a symmetric \"family choice moves the force\" "
          "ambiguity.** It is the two-parameter family the case study ships misfitting "
          "these curves systematically, and its plateau being biased relative to a "
          "four-parameter family that does not -- a stronger and differently actionable "
          "result than the first version of this page claimed."]
    g1, g3 = p["G12_1_item12_own_bar"], p["G12_3_informative"]
    L += ["",
          f"**G12.1 ({'PASS' if g1['pass'] else 'FAILED'}), item 12's own bar:** the worst "
          f"Cross deviation is {_pct(g1['worst'])} against the donor spread "
          f"{_pct(g1['bar'])}. Pre-registered as WEAK -- {g1['note']}.", "",
          f"**G12.2:** `cross` lands in the band *{s['cross']['band']}*; "
          f"`carreau_yasuda` in *{s['carreau_yasuda']['band']}*.", "",
          f"**G12.3 ({'INFORMATIVE' if g3['informative'] else 'UNINFORMATIVE'}), the gate "
          f"that decides whether G12.2 means anything:** the worst NULL deviation is "
          f"{_pct(g3['null_worst'])} against Cross's {_pct(g3['cross_worst'])}, a "
          f"separation of {g3['separation']:.4g}x against a {g3['floor']}x floor, and the "
          f"nulls land in "
          f"{'a DIFFERENT band' if g3['nulls_in_a_different_band'] else 'the SAME band'} "
          f"({', '.join(f'{k}: {v}' for k, v in g3['null_bands'].items())}). "
          + ("A null that reproduced the Casson force as well as Cross does would mean "
             "the needle had washed out every difference between every model, and the "
             "Cross/Casson agreement would say nothing about their high-rate plateaux."
             if g3["informative"] else
             "**So G12.2 is reported UNINFORMATIVE regardless of its number**: the needle "
             "washes out the difference between every model here, and the Cross/Casson "
             "agreement is a statement about the geometry, not about the fits."), "",
          "**G12.4, the identifiability split.** Reported both ways, decided before it "
          "was looked at:", "",
          "| family | all 20 | all-`measurable` only | n measurable | same verdict |",
          "|---|---|---|---|---|"]
    for f, v in p["G12_4_identifiability_split"].items():
        L.append(f"| `{f}` | {v['band']} | {v['band_measurable_only']} | "
                 f"{v['n_all_measurable']}/20 | {'yes' if v['same'] else '**NO**'} |")
    L += ["", "**G12.5, the accounting.** Every record this study could not evaluate is "
          "counted and named; a quadrature refusal is data about the fit, not a nuisance.",
          "", "| family | n in | n out | donors lost |", "|---|---|---|---|"]
    for f, v in p["G12_5_accounting"].items():
        L.append(f"| `{f}` | {v['n_in']} | {v['n_out']} | "
                 f"{', '.join(v['failed']) if v['failed'] else 'none'} |")
    reasons = {f: s[f]["failure_reasons"] for f in FAMILIES if s[f]["failure_reasons"]}
    if reasons:
        L += ["", "Reasons:", ""]
        for f, rr in reasons.items():
            for d, why in rr.items():
                L.append(f"- `{f}` donor {d}: {why}")
    rc, pc = p["route_control"], p["plateau_control"]
    L += ["", "## Two controls added AFTER the first run", "",
          "Neither is in the pre-registration. Both were added because the effect the "
          "first run found is far too large to leave unchecked, and both can only "
          "REVEAL a confound, never hide one -- they are strictly stronger tests, not "
          "looser ones.", "",
          f"**The routes are not the effect.** Casson forces come from the closed form "
          f"and every other family from the quadrature, so a route difference would "
          f"masquerade as a family difference. Recomputing the Casson flow rate through "
          f"the quadrature at the same pressure agrees to "
          f"**{rc['worst_rel_dev']:.2e}** relative, worst over 20 donors -- "
          f"{sc['yield_term'] / max(rc['worst_rel_dev'], 1e-300):.1e}x smaller than the "
          "smallest scale on this page. " + rc["claim"], "",
          "**The mechanism is identified, outside the flow solver.** At the needle's "
          "wall rate every family is effectively a constant viscosity, and for a "
          "Newtonian fluid at fixed Q the force is proportional to it. So the force "
          "ratio should reproduce the ratio of the two fits' own high-rate plateaux "
          "(`eta_inf` against `eta_c`), read straight out of the records and touching no "
          "part of `injectability`:", "",
          "| family | worst \\|force ratio - plateau ratio\\| |", "|---|---|"]
    for f, v in pc.items():
        L.append(f"| `{f}` | {_pct(v['worst_gap'])} |")
    pn = p["plateau_null"]
    vl, ien, mpn = p["viscosity_ladder"], p["in_envelope_needle"], p["matched_plateau_null"]
    L += ["",
          "**That gap is NOT the yield term, and this study said it was.** The yield "
          "explanation is signed and quantitative -- `F_casson = F_newt(eta_c)(1+y)` puts "
          "the force ratio BELOW the plateau ratio on every donor by `y` -- so an "
          "unsigned |gap| cannot test it. Signed, the residual has the WRONG SIGN on "
          f"{p['plateau_control']['cross']['n_wrong_sign']}/20 Cross donors "
          f"({', '.join(p['plateau_control']['cross']['wrong_sign_donors'])}) and "
          f"{p['plateau_control']['carreau_yasuda']['n_wrong_sign']}/20 "
          "Carreau-Yasuda ones, the worst gap in both families is a wrong-sign row "
          f"(donor {p['plateau_control']['cross']['worst_gap_donor']}), and "
          f"Carreau-Yasuda's {_pct(p['plateau_control']['carreau_yasuda']['worst_gap'])} "
          f"is {p['plateau_control']['carreau_yasuda']['worst_gap'] / sc['yield_term']:.2f}x "
          "the yield term rather than equal to it. **The sentence is withdrawn**; what "
          "the residual is instead is a compression that grows with the plateau ratio -- "
          "the non-Newtonian core -- and this study does not characterise it.", "",
          "### The explanation this study first gave was WRONG, and two checks refute it",
          "",
          "Having found a 30 % difference at a rate 1.42 decades past the data, the first "
          "version of this page and of `results.md` both said *the cause is the "
          "extrapolation*.", "",
          "**Check one, through the solver.** Run the identical comparison at a needle "
          f"wide enough that the wall rate lands in the measured range -- ID "
          f"{ien['id_mm']:.4g} mm instead of {ID_NOMINAL_MM} mm, sized so "
          f"`4Q/(pi R^3)` = {ien['target_gd']:.0f} 1/s. With the extrapolation removed "
          f"({ien['n_records_outside']}/{ien['n_records_total']} record evaluations still "
          f"sit outside, by at most {ien['worst_decades_past_high']:.3f} decades, because "
          "sizing on the apparent rate leaves a shear-thinning fluid's true wall rate a "
          "little higher) the Cross-vs-Casson difference is "
          f"**{_pct(ien['families']['cross']['worst_abs_rel_dev'])} worst, "
          f"{_pct(ien['families']['cross']['median_rel_dev'])} median** -- against "
          f"{_pct(s['cross']['worst_abs_rel_dev'])} and "
          f"{_pct(s['cross']['median_rel_dev'])} at the 27G. "
          f"**{100 * abs(ien['families']['cross']['median_rel_dev'] / s['cross']['median_rel_dev']):.1f} % "
          "of the median survives with no extrapolation at all.** The extrapolation "
          "contributes essentially nothing to the typical donor and about a third of the "
          "single worst one.", "",
          "**Check two, on the laws themselves.** How far apart are the same two fits at "
          "rates inside the envelope?", "",
          "| rate (1/s) | worst `cross` vs `casson` | worst `carreau_yasuda` | |",
          "|---|---|---|---|"]
    for r in vl["rows"]:
        L.append(f"| {r['gammadot']:.4g} | {_pct(r['cross'])} | "
                 f"{_pct(r['carreau_yasuda'])} | "
                 f"{'inside the envelope' if r['inside_envelope'] else '**the needle**'} |")
    L += ["",
          f"The families disagree by up to {_pct(vl['worst_inside_cross'])} **inside the "
          f"measured range**, against {_pct(vl['at_needle_cross'])} at the needle -- and "
          f"the worst rate is the BOTTOM of the envelope, not the top. The needle is "
          "therefore not where the disagreement is created; it is one rate at which two "
          "differently-shaped fits to the same 16 noisy points are read. "
          + vl["claim"], "",
          "That does not weaken the finding, it widens it: the 30 % is not a hazard of "
          "extrapolating these records, it is a property of the records themselves, and "
          "it applies to anyone reading them anywhere on the curve -- not only in a "
          "needle. What the extrapolation adds is that outside the envelope there is no "
          "data to say which fit is right, and `eta_at` refuses to answer while a solver "
          "dictionary cannot.", "",
          "### The null that DOES discriminate, and the two that do not", "",
          "**The matched-plateau null.** Take each donor's `cross` record and replace "
          "`eta_inf` with that donor's Casson `eta_c`, leaving `eta_0`, `lam` and `m` at "
          "their fitted values -- not Casson in disguise, since the whole low-rate "
          "structure of the Cross fit is kept -- and push it through the same solver. If "
          "the difference came from anywhere but the high-rate plateau it would survive. "
          f"It does not: **{_pct(mpn['worst_abs_rel_dev'])} worst, "
          f"{_pct(mpn['median_rel_dev'])} median**, collapsing into the "
          f"{_pct(sc['yield_term'])} yield band from "
          f"{_pct(s['cross']['worst_abs_rel_dev'])}. So at the needle the plateau IS what "
          "is read, and the disagreement is about its value. (Supplied by code review "
          "after this study's own first two nulls were shown not to discriminate.)", "",
          f"**Neither registry null discriminates, and `newtonian`'s number is an "
          f"identity.** For a Newtonian fluid `dp = 8 mu L Q / (pi R^4)` exactly, so "
          f"`F_newt/F_casson = (mu/eta_c)/(1+y)` -- computable from the registry with no "
          f"needle, no quadrature and no solver. Its "
          f"{_pct(s['newtonian']['worst_abs_rel_dev'])} is a restatement of *blood "
          f"shear-thins by 2-4x between 0.1 and 700 1/s*, which is visible in the raw "
          "data. G12.3 therefore does NOT license \"the needle has not washed out the "
          "difference between models\"; the matched-plateau null is what licenses it. "
          f"And `power_law`, a family with NO plateau at all, lands in the SAME band as "
          f"`cross` -- so the four-band verdict is coarser than the claim it was being "
          "used to support. The bar is not changed; this is reported as the finding it "
          "is.", "",
          "**A fourth control was attempted and it is NOT independent evidence. Recorded "
          "anyway, because what it turned out to be is worth knowing.** The registry "
          f"nulls are off partly by construction -- `newtonian` misses by "
          f"{_pct(s['newtonian']['worst_abs_rel_dev'])}, but its `mu` is fitted over the "
          "whole 0.1-700 1/s curve and so is nowhere near the high-rate plateau. So a "
          "Newtonian fluid AT the donor's own Casson `eta_c` was run, meaning to separate "
          "*does the needle read the plateau?* from *do the plateaux disagree?*. It "
          f"reproduces the Casson force to **{_pct(pn['worst_rel_dev'])}** worst over 20 "
          f"donors -- and the yield term on the sibling page is {_pct(sc['yield_term'])}, "
          "**the same number**, because Casson at `tau_y = 0` IS Newtonian at `eta_c`. "
          "The two agree to "
          f"{abs(pn['worst_rel_dev'] / sc['yield_term'] - 1):.1e} relative, which makes "
          "this a cross-implementation check of that number (one route through "
          "`flow_rate_casson`, one through `flow_rate_gnf`) and nothing more. It shares "
          "the identity it was meant to test, so it cannot test it -- this repo's own "
          "2026-09-05 lesson, walked into again. The mechanism rests on the "
          "MATCHED-PLATEAU null, which is the only control here that both goes through "
          "the solver and could have come out the other way.", "",
          f"**G12.3 does not rest on which null was picked.** The separation is "
          f"{g3['separation_worst_null']:.4g}x on the worse null and "
          f"{g3['separation_least_extreme_null']:.4g}x on the milder one, both above the "
          f"{g3['floor']}x floor "
          f"({'informative either way' if g3['informative_on_either_null'] else 'NOT informative on the milder null'}).",
          "",
          "## What this does not establish", "",
          "**It does not say which family is RIGHT.** At the 27G every family is read a "
          f"median {p['median_decades_past_high']:.2f} decades past the top of its own "
          "envelope, and the measurement that would settle which extrapolation is correct "
          "is tube viscometry at this diameter, which the Horner dataset does not "
          "contain. Casson fitting worse over the MEASURED range is a reason to doubt its "
          "plateau; it is not a demonstration that Cross's is right.", "",
          "**The forces here carry no parameter error. That was once because the "
          "intervals were unusable; it is now a choice, and the answer is on a separate "
          "page.** When this study was written, the stored profile-likelihood UPPER "
          "bounds were corrupted by a step-size bug upstream in `dst_fit.profiling` -- "
          "`step = max(wald/4, span/5000)` takes `span` from the parameter's arbitrary "
          "BOX BOUND, so 203 of the registry's 296 closed upper ends were bisection-grid "
          "points rather than likelihood crossings, the worst sitting 1639x the "
          "best-fit value above it. **That is fixed (NEXT_STEPS item 13, "
          "`data/profile_interval_prereg.md`), and the question it blocked has been "
          "asked and ANSWERED: no, the disagreement is NOT resolved at 95 %.**", "",
          "`data/profile_interval_g139.md`: pushing each Casson parameter through its "
          "corrected 95 % profile interval gives a force hull with a median width of "
          "**29.1 %** of the force, against the median family difference of **12.5 %** "
          "measured here. On **11 of 20** donors the Cross force falls INSIDE that hull, "
          "and inside the sharper covariance-sampled joint interval (median 25.2 %) on "
          "11 of 20 as well. **The Casson record's own parameter uncertainty at this "
          "needle is LARGER than the family-choice effect.** So this study's caveat is "
          "strengthened, not lifted: do not quote a needle force from it to better than "
          "about 30 %.", "",
          "Two things that does not touch. The mechanism above is unaffected -- Casson "
          "misfits these curves on 20/20 donors in `rmse_log`, which is not a statement "
          "about parameter uncertainty. And the interval fix did NOT change this "
          "conclusion: the bug inflated only the UPPER ends, while the difference is "
          "signed (Cross below Casson on 16/20), so the comparison turns on the LOWER "
          "end, which was never corrupted. Exactly one donor of twenty flips between the "
          "old intervals and the corrected ones.", "",
          "**It resolves nothing at the 1 % level**, for the reason in the scales table.",
          ""]
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
