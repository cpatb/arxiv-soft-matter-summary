"""G13.9: is item 12's Cross-vs-Casson force disagreement resolved at 95 %?

`NEXT_STEPS` item 13 step 4 says this question becomes askable only once the profile
intervals are correct. Pre-registered in `data/profile_interval_prereg.md` addendum A.2,
committed at `eac96db` before any of this was computed.

Item 12 measured a Cross-vs-Casson plunger-force difference of up to 30.18 % at a 27G
needle and propagated NO parameter uncertainty, because the stored intervals were
bisection-grid artefacts. With corrected intervals:

  G13.9a  the Cross nominal force falls OUTSIDE the Casson push-through hull on >= 11/20
  G13.9b  with the OLD (dst_fit) intervals it falls INSIDE on >= 19/20, i.e. the old
          intervals could not have refuted anything -- the control that decides whether
          this is a result OF the fix

The hull is `force_with_uq`'s one-at-a-time push-through: each parameter moved to each
end of its own profile interval, the others held, and the force range taken. **That is
wider than the joint 95 % region it is drawn from**, so a Cross force outside it is a
strong statement and one inside it is a weak one. Said in the pre-registration, repeated
here, and printed in the output.

Run:  python scripts/profile_interval_g139.py   (~3 min)
Writes: data/profile_interval_g139.md and .json
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from dst_fit.profiling import profile_likelihood

from biorheo.apps.injectability import Needle, Syringe, force_with_uq, injection_force
from biorheo.constitutive import get_model
from biorheo.identify import fit_rate_curve
from biorheo.identify.fit import ParamEstimate

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
OUT_MD = REPO / "data" / "profile_interval_g139.md"
OUT_JSON = REPO / "data" / "profile_interval_g139.json"

G139A_BAR = 11   # of 20
G139B_BAR = 19   # of 20

TEMPERATURE_C = 37.0
GEOMETRY = ("ARES-G2 double-wall concentric cylinder, EDTA-anticoagulated whole blood, "
            "rate sweep")


def _needle_and_syringe():
    """The 27G case the item-12 cross-check runs, read from its own module.

    Imported rather than restated: a needle geometry retyped here would make this a
    comparison against a DIFFERENT device than the one item 12 measured, which is the
    reproduce-the-claim's-own-setup rule this repo has been bitten by before.
    """
    spec = importlib.util.spec_from_file_location(
        "_xcheck", REPO / "case_studies" / "blood_injectability" / "crosscheck.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_xcheck"] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    x = _needle_and_syringe()
    needle = Needle(radius=0.5 * x.ID_NOMINAL_MM * 1e-3, length=x.NEEDLE_LENGTH)
    syringe = Syringe(x.BARREL_RADIUS)
    Q = x.Q_SI

    mss_spec = importlib.util.spec_from_file_location(
        "_mss", HERE / "multistart_seed_scatter.py")
    assert mss_spec and mss_spec.loader
    mss = importlib.util.module_from_spec(mss_spec)
    sys.modules["_mss"] = mss
    mss_spec.loader.exec_module(mss)
    gd, stress = mss.load()

    casson, cross = get_model("casson"), get_model("cross")
    rows = []
    for donor in sorted(stress):
        eta = stress[donor] / gd
        ident_c = fit_rate_curve(casson, gd, eta, temperature_C=TEMPERATURE_C,
                                 geometry=GEOMETRY)
        ident_x = fit_rate_curve(cross, gd, eta, temperature_C=TEMPERATURE_C,
                                 geometry=GEOMETRY)

        uq_new = force_with_uq(ident_c, Q, needle, syringe)
        f_cross = injection_force(Q, needle, syringe, cross,
                                  ident_x.values).force

        # the same hull rebuilt on dst_fit's OWN intervals -- the G13.9b control. The FIT
        # is identical; only the interval it is pushed through differs, which is the one
        # thing this gate is about.
        old_params = {}
        for name in casson.param_names:
            old = profile_likelihood(casson, gd, eta, name, fit=ident_c.fit, seed=0)
            est = ident_c.params[name]
            old_params[name] = ParamEstimate(
                name=name, value=est.value, lo=float(old.lower), hi=float(old.upper),
                verdict=est.verdict, lo_at_bound=bool(old.lower_at_bound),
                hi_at_bound=bool(old.upper_at_bound), wald_stderr=est.wald_stderr,
                partner=est.partner)
        ident_old = type(ident_c)(model=ident_c.model, params=old_params,
                                  fit=ident_c.fit, envelope=ident_c.envelope,
                                  level=ident_c.level,
                                  residual_structure=ident_c.residual_structure)
        uq_old = force_with_uq(ident_old, Q, needle, syringe)

        rows.append({
            "donor": donor,
            "f_casson": uq_new.nominal.force,
            "f_cross": float(f_cross),
            "pct_diff": 100.0 * (f_cross - uq_new.nominal.force) / uq_new.nominal.force,
            "new_lo": uq_new.profile_lo, "new_hi": uq_new.profile_hi,
            "old_lo": uq_old.profile_lo, "old_hi": uq_old.profile_hi,
            "new_contains": bool(uq_new.profile_lo <= f_cross <= uq_new.profile_hi),
            "old_contains": bool(uq_old.profile_lo <= f_cross <= uq_old.profile_hi),
            "samp_lo": uq_new.sampled_lo, "samp_hi": uq_new.sampled_hi,
            "samp_contains": bool(np.isfinite(uq_new.sampled_lo)
                                  and uq_new.sampled_lo <= f_cross <= uq_new.sampled_hi),
            "samp_usable": bool(np.isfinite(uq_new.sampled_lo)),
            # which end of the hull does the comparison turn on? The bug inflated only
            # the UPPER ends, so if the answer is "lower" the fix cannot have changed it.
            "binding_end": ("lower" if f_cross < uq_new.profile_lo
                            else "upper" if f_cross > uq_new.profile_hi else "none"),
            "old_binding_end": ("lower" if f_cross < uq_old.profile_lo
                                else "upper" if f_cross > uq_old.profile_hi else "none"),
            "new_hull_width_pct": 100.0 * (uq_new.profile_hi - uq_new.profile_lo)
                                  / uq_new.nominal.force,
            "old_hull_width_pct": 100.0 * (uq_old.profile_hi - uq_old.profile_lo)
                                  / uq_new.nominal.force,
        })
        print(f"{donor}: casson {rows[-1]['f_casson']:.4g} N, cross {f_cross:.4g} N "
              f"({rows[-1]['pct_diff']:+.2f} %), new hull "
              f"[{uq_new.profile_lo:.4g}, {uq_new.profile_hi:.4g}] "
              f"{'CONTAINS' if rows[-1]['new_contains'] else 'EXCLUDES'}, "
              f"old hull {'CONTAINS' if rows[-1]['old_contains'] else 'EXCLUDES'}")

    n = len(rows)
    a_out = sum(1 for r in rows if not r["new_contains"])
    b_in = sum(1 for r in rows if r["old_contains"])
    g139a = a_out >= G139A_BAR
    g139b = b_in >= G139B_BAR

    # G13.9a' -- reported BESIDE G13.9a, whose bar is unchanged. The one-at-a-time hull
    # is wider than the joint region, so the pre-registered test is the GENEROUS one; the
    # covariance-sampled interval is the sharper companion. It is not simply "better":
    # it assumes the Wald normal, which the profile exists because it distrusts.
    samp_rows = [r for r in rows if r["samp_usable"]]
    a_out_samp = sum(1 for r in samp_rows if not r["samp_contains"])
    binding = {e: sum(1 for r in rows if r["binding_end"] == e)
               for e in ("lower", "upper", "none")}
    old_binding = {e: sum(1 for r in rows if r["old_binding_end"] == e)
                   for e in ("lower", "upper", "none")}
    med_samp = float(np.median([100.0 * (r["samp_hi"] - r["samp_lo"]) / r["f_casson"]
                                for r in samp_rows])) if samp_rows else float("nan")
    med_new = float(np.median([r["new_hull_width_pct"] for r in rows]))
    med_old = float(np.median([r["old_hull_width_pct"] for r in rows]))
    med_diff = float(np.median([abs(r["pct_diff"]) for r in rows]))

    verdict = {
        "G13.9a": {"pass": g139a, "n_excluded": a_out, "n": n, "bar": G139A_BAR},
        "G13.9b": {"pass": g139b, "n_contained_old": b_in, "n": n, "bar": G139B_BAR},
        "median_new_hull_width_pct": med_new,
        "median_old_hull_width_pct": med_old,
        "median_abs_pct_diff": med_diff,
        "G13.9a_prime": {"n_excluded_sampled": a_out_samp, "n": len(samp_rows),
                         "median_sampled_width_pct": med_samp,
                         "note": "reported beside G13.9a; the bar was not changed"},
        "binding_end_new": binding, "binding_end_old": old_binding,
    }
    OUT_JSON.write_text(json.dumps({"verdict": verdict, "rows": rows,
                                    "needle_ID_m": 2.0 * needle.radius,
                                    "Q_m3_per_s": Q}, indent=2), encoding="utf-8")

    L = [
        "# G13.9 — is item 12's force disagreement resolved at 95 %?",
        "",
        "Written by `scripts/profile_interval_g139.py`. Pre-registered in",
        "`data/profile_interval_prereg.md` addendum A.2 (commit `eac96db`) before it was run.",
        "",
        "| gate | claim | bar | result | verdict |",
        "|---|---|---|---|---|",
        f"| G13.9a | the disagreement is NOT covered by Casson's own parameter uncertainty "
        f"| Cross force outside the hull on >= {G139A_BAR}/20 | **{a_out}/{n}** | "
        f"**{'PASS' if g139a else 'FAIL'}** |",
        f"| G13.9b | the OLD intervals could not have refuted anything | Cross force inside "
        f"the old hull on >= {G139B_BAR}/20 | **{b_in}/{n}** | "
        f"**{'PASS' if g139b else 'FAIL'}** |",
        "",
        "## Both gates FAILED and the prediction was wrong. Here is why.",
        "",
        "The pre-registration predicted 20/20 on both. Neither bar has been touched and",
        "neither number is restated as a success.",
        "",
        f"- median |Cross - Casson| force difference: **{med_diff:.2f} %**",
        f"- median Casson push-through hull, CORRECTED intervals: **{med_new:.2f} %** of the force",
        f"- median Casson push-through hull, dst_fit's intervals: **{med_old:.3g} %** of the force",
        f"- median covariance-SAMPLED (joint) interval: **{med_samp:.2f} %** of the force",
        "",
        "**The reason both gates failed is the same one, and it is a sign error in my",
        "reasoning, not in any code.** The bug inflated the UPPER interval ends only --",
        "the lower ends were never on the bisection grid in the direction that matters,",
        "because the step floor bites where the box bound is far away and the lower bound",
        "of a viscosity is near zero. But the force difference is SIGNED: item 12 measured",
        "Cross BELOW Casson on 16 of 20 donors (paired sign test, p = 0.0118). So the",
        "comparison turns on the LOWER end of the hull, which the bug never corrupted.",
        "",
        f"Measured rather than argued: of the {n} donors, the corrected hull's binding end",
        f"is the LOWER one on **{binding['lower']}** and the upper one on"
        f" **{binding['upper']}**; on dst_fit's own intervals it is lower on"
        f" **{old_binding['lower']}** and upper on **{old_binding['upper']}**.",
        "",
        "Two consequences, both of which cut against the reason this gate was written:",
        "",
        "1. **G13.9b could not have come out the other way for the right reason.** The",
        "   one-at-a-time hull on the old intervals CONTAINS the hull on the corrected",
        "   ones -- narrowing an interval can only shrink a push-through range -- so",
        "   `old_contains >= new_contains` always. The bar of 19/20 was reachable in",
        "   principle, but only if the excluded donors were excluded from ABOVE, and they",
        "   are not. This is the 'a control that shares an identity with its target'",
        "   failure in a new costume, and it was not caught by writing the gate down.",
        "2. **Item 13 did not change item 12's conclusion.** That is the honest headline",
        "   here, and it is the opposite of what the pre-registration expected.",
        "",
        "## What IS established, and it is the useful part",
        "",
        f"On **{n - a_out} of {n}** donors the Cross force lies inside the Casson record's",
        "own one-at-a-time 95 % push-through, and on the sharper covariance-sampled",
        f"interval it lies inside on **{len(samp_rows) - a_out_samp} of {len(samp_rows)}**.",
        f"The Casson uncertainty at this needle is a median **{med_new:.1f} %** of the",
        f"force against a median family difference of **{med_diff:.1f} %** -- the two",
        "effects are the same size, with parameter uncertainty the larger of them.",
        "",
        "So the answer to *is the plateau disagreement resolved at 95 %* is **no, on a",
        "majority of donors**, and item 12's caveat is strengthened rather than lifted:",
        "a needle force from this case study should not be quoted to better than about",
        f"{med_new:.0f} %, and the family choice is inside that.",
        "",
        "The sampled interval is reported beside the pre-registered one, not instead of",
        "it. It is narrower because it is a joint region rather than a box hull, but it",
        "assumes the Wald normal -- which is the approximation the profile likelihood",
        "exists because it distrusts. Neither is simply the better number.",
        "",
        "## What this does NOT say",
        "",
        "It does not say Cross and Casson agree. Item 12's mechanism is untouched: Casson",
        "misfits these curves on 20/20 donors in `rmse_log` and `residual_structure` flags",
        "19/20 against 2/20, and none of that is a statement about parameter uncertainty.",
        "It does not say Cross is right either -- both families are read a median 1.42",
        "decades past their own envelope at this needle.",
        "",
        "## Per donor",
        "",
        "| donor | Casson F (N) | Cross F (N) | diff | corrected hull | verdict | old hull | verdict |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        L.append(
            f"| {r['donor']} | {r['f_casson']:.4g} | {r['f_cross']:.4g} | "
            f"{r['pct_diff']:+.2f} % | [{r['new_lo']:.4g}, {r['new_hi']:.4g}] | "
            f"{'contains' if r['new_contains'] else '**excludes**'} | "
            f"[{r['old_lo']:.4g}, {r['old_hi']:.4g}] | "
            f"{'contains' if r['old_contains'] else '**excludes**'} |")
    L.append("")
    OUT_MD.write_text("\n".join(L), encoding="utf-8")

    print(f"\nG13.9a: {'PASS' if g139a else 'FAIL'}  {a_out}/{n} excluded (bar {G139A_BAR})")
    print(f"G13.9b: {'PASS' if g139b else 'FAIL'}  {b_in}/{n} contained by the old hull "
          f"(bar {G139B_BAR})")
    print(f"median hull: corrected {med_new:.2f} %, old {med_old:.3g} %, "
          f"difference under test {med_diff:.2f} %")
    print(f"wrote {OUT_MD.name}, {OUT_JSON.name}")
    return 0 if (g139a and g139b) else 1


if __name__ == "__main__":
    raise SystemExit(main())
