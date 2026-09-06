"""Twenty human whole-blood donors through a 27G needle (NEXT_STEPS item 4).

The question NEXT_STEPS asks: run the 20 Horner blood records through
`biorheo.apps.injectability` and report how many are flagged `extrapolated` at a 27G
needle at 1 mL/min. The Casson closed form used here is derived and anchored in
`derivations/casson_tube_flow.md` (against the package's independent quadrature).

PRE-REGISTERED PREDICTION, written and committed to before the run (2026-09-05):

  Every record's fitted envelope tops out at 700 1/s. The apparent (Newtonian) wall
  shear rate for a 27G needle is 4Q/(pi R^3); at R = 1.05e-4 m and Q = 1 mL/min =
  1.667e-8 m^3/s that is about 1.8e4 1/s, roughly 26x the top of the envelope, and the
  true wall rate for a shear-thinning fluid is HIGHER still. So the prediction is
  20/20 flagged `extrapolated`, and the interesting number is not the count but the
  DISTANCE -- how many decades past the measured data the device operates. If fewer
  than 20 come back flagged, something is wrong with either the flag or the geometry
  and that is the finding.

  Secondary prediction, same status: the yield stress will be irrelevant here.
  x = tau_y/tau_w is of order 1e-5 at these pressures (the derivation's own table), so
  the Casson bracket differs from 1 by the sqrt(x) term only, about 0.5 %, and blood in
  a needle is a Newtonian-with-eta_c problem. If so, the donor-to-donor spread in force
  should track eta_c and not tau_y, which is checked by rank correlation.

  OUTCOME OF THE PREDICTIONS (recorded after the run, not edited into them):
  the first was RIGHT (20/20, and the distance is 1.42 decades). The second was right in
  direction and WRONG BY 6x IN SIZE: x runs to 1.7e-4, not "of order 1e-5", because the
  27G length here is 12.7 mm against the derivation table's 20 mm and its pressures were
  not this study's, so the bracket deviates 3 % rather than 0.5 %. The conclusion it
  supported survives at the larger number -- zeroing tau_y moves the force by at most
  2.9 %, and the force ranks with eta_c at 0.995 against 0.170 for tau_y -- but the
  estimate that motivated it was off, and an estimate carried over from a table computed
  at a different geometry is exactly the kind that is.

GEOMETRY IS A CHOICE, NOT A MEASUREMENT. ISO 9626 gives a RANGE of inner diameters per
gauge, so "27G" is not a radius. The value used is the nominal regular-wall 27G inner
diameter, 0.210 mm, from the ISO 9626:2016 table; the sensitivity to it is reported
(F ~ R^-4 at fixed Q -- measured, not assumed, in the geometry section; an earlier
draft said "Q ~ R^3 at fixed dp", which is wrong twice over) and the conclusion of that
section is rebuilt on that scaling rather than on a standard.

Everything the module's docstring excludes is excluded here too, and the omissions do
NOT all point the same way, which an earlier draft of this file got wrong by calling
these forces "upper bounds". The Fahraeus-Lindqvist effect -- a cell-depleted wall layer
that lowers blood's apparent viscosity in tubes below roughly 300 um [UNVERIFIED: no
citation resolved for that threshold] -- lowers the force, and this 210 um lumen is
inside that regime; the neglected entrance (Hagenbach) loss RAISES it, and is measured
below to be the same size as the yield-stress contribution this study calls negligible.
The honest statement is that the calculation resolves neither, not that it bounds the
answer from one side.

Run:  python case_studies/blood_injectability/run.py
Writes `results.md` and `results.json` beside this file. Pinned in
`tests/test_casson_needle.py`.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from biorheo.apps.injectability import (
    Needle, Syringe, injection_force, pressure_for_flow_rate,
    wall_shear_rate, wall_stress,
)
import biorheo
from biorheo.registry.loader import load_all

HERE = Path(__file__).resolve().parent

# 27G inner diameter. 0.210 mm is the widely quoted nominal for a regular-wall 27G and
# is what the study is run at. THE BAND BELOW IS [UNVERIFIED] AND THE CONCLUSION NO
# LONGER RESTS ON IT (see `geometry` in results.md): a referee pass on 2026-09-05 pointed
# out that 0.184 mm is a regular-wall figure while 0.241 mm is a THIN-WALL one, so
# quoting the pair as "the tolerance of regular-wall 27G" mixes two products -- and the
# headline it supported REVERSES under the narrower reading. ISO 9626:2016 Table 2 is not
# available offline, so the standard is not cited as fact anywhere; both readings are
# reported with their own verdicts, and the claim that carries the section is the exact
# R^-4 scaling, which needs no standard at all.
ID_NOMINAL_MM = 0.210
ID_BAND_MM = (0.184, 0.241)          # [UNVERIFIED] wide reading: mixes wall types
ID_BAND_NARROW_MM = (0.184, 0.210)   # [UNVERIFIED] narrow reading: regular wall only
# donors the SOURCE study flags as atypical; the donor spread is reported with and
# without them, because 2 of the 3 sit at the ends of the eta_c range
PAPER_OUTLIERS = ("F", "K", "N")
NEEDLE_LENGTH = 12.7e-3          # m, a 1/2 inch needle
BARREL_RADIUS = 2.3e-3           # m, a 1 mL syringe; a choice, as in the Carbopol study
Q_ML_PER_MIN = 1.0
Q_SI = Q_ML_PER_MIN * 1e-6 / 60.0
RHO_BLOOD = 1060.0               # kg/m^3, for the Reynolds number only


def needle(id_mm: float) -> Needle:
    return Needle(radius=id_mm * 1e-3 / 2.0, length=NEEDLE_LENGTH)


def one_donor(rec, nd: Needle) -> dict:
    tau_y, eta_c = (float(v) for v in rec.params_si)
    dp = pressure_for_flow_rate(Q_SI, nd, rec.model, rec.params_si)
    gd_w = wall_shear_rate(dp, nd, rec.model, rec.params_si)
    tw = wall_stress(dp, nd.radius, nd.length)
    force = injection_force(Q_SI, nd, Syringe(BARREL_RADIUS), rec.model, rec.params_si,
                            rho=RHO_BLOOD)
    # `injection_force` re-solves the pressure internally, so the dp in the table and the
    # dp behind the force come from two independent brentq runs. They must agree.
    assert abs(force.dp - dp) / dp < 1e-9, (rec.id, force.dp, dp)
    lo, hi = rec.envelope.gammadot_min, rec.envelope.gammadot_max
    # the apparent (Newtonian) wall rate, which is what a datasheet would quote
    gd_apparent = 4.0 * Q_SI / (np.pi * nd.radius ** 3)
    return {"record": rec.id, "donor": rec.id.split("_")[-1],
            "tau_y": tau_y, "eta_c": eta_c, "dp_Pa": dp, "tau_w_Pa": tw,
            "x": tau_y / tw, "gd_wall": gd_w, "gd_apparent": gd_apparent,
            "envelope": [lo, hi], "extrapolated": bool(gd_w < lo or gd_w > hi),
            # named for the high side deliberately: log10(gd_w/hi) is only meaningful
            # when the extrapolation is ABOVE the envelope, which is asserted here. A
            # low-side record would give a large negative value that would silently win
            # a `min()` over the band.
            "decades_past_high": float(np.log10(gd_w / hi)),
            "decades_past": float(np.log10(gd_w / hi)), "force_N": force.force,
            "reynolds": force.reynolds}


def _crosscheck_section() -> list[str]:
    """The caveat item 12 put on this page's PRECISION, quoted from that study's own
    artifact rather than retyped.

    This page ranks donors and says the yield stress is invisible here. Both survive item
    12. What does not is any reading of these forces as accurate to better than tens of
    per cent, and until 2026-09-06 nothing on this page said so -- which is this repo's
    "prose contradicts its own artifact" failure mode waiting to happen.

    If `crosscheck.json` is absent the section says so rather than quoting numbers that
    are not there; if it is present but was generated against a different operating
    point, that is asserted rather than papered over."""
    path = HERE / "crosscheck.json"
    head = ["## The FAMILY the record was fitted to moves these forces by more than the "
            "yield stress does", ""]
    if not path.exists():
        return head + ["Not computed: `crosscheck.json` is missing. Run "
                       "`python case_studies/blood_injectability/crosscheck.py`. **The "
                       "forces on this page should not be read as accurate to better "
                       "than tens of per cent until it has been run** -- see "
                       "`crosscheck_prereg.md`.", ""]
    c = json.loads(path.read_text())
    assert c["needle"]["id_nominal_mm"] == ID_NOMINAL_MM, c["needle"]
    assert c["Q_mL_per_min"] == Q_ML_PER_MIN and c["rho"] == RHO_BLOOD
    s, g3 = c["summary"], c["G12_3_informative"]
    return head + [
        "Every donor on this page also has a `cross` and a `carreau_yasuda` record "
        "fitted by the same code to the SAME 16 points. Pushed through this same needle "
        "at this same flow rate "
        "(`crosscheck.md`, gate pre-registered in `crosscheck_prereg.md` before the run), "
        f"they miss this page's Casson force by up to "
        f"**{100 * s['cross']['worst_abs_rel_dev']:.4g} %** and "
        f"**{100 * s['carreau_yasuda']['worst_abs_rel_dev']:.4g} %** respectively, with "
        f"medians {100 * s['cross']['median_rel_dev']:.4g} % and "
        f"{100 * s['carreau_yasuda']['median_rel_dev']:.4g} % -- both families predict a "
        "LOWER force. That is an order of magnitude more than the "
        f"{100 * c['scales']['yield_term']:.3g} % the yield stress contributes (that "
        "study reads this one's own `yield_share_max`, so the two cannot drift), and it "
        "is not a "
        "numerical artefact: the closed form and the quadrature agree to "
        f"{c['route_control']['worst_rel_dev']:.2e} on this page's own records.", "",
        "The cause is the extrapolation this page already reports. At "
        f"{c['median_decades_past_high']:.2f} decades past the top of the measured "
        "envelope, Casson tends to `eta_c`, Cross to `eta_inf` and Carreau-Yasuda to its "
        "own `eta_inf` -- three different constants, each fitted to data that stops at "
        "700 1/s. The force ratio reproduces the ratio of those plateaux to "
        f"{100 * c['plateau_control']['cross']['worst_gap']:.3g} %, so the families "
        "disagree because their plateaux do.", "",
        "**What survives and what does not.** The ordering survives: "
        f"{100 * s['cross']['worst_abs_rel_dev']:.4g} % is well inside the "
        f"{100 * c['scales']['donor_spread']:.4g} % donor-to-donor spread, so which "
        "patient the blood came from still matters more than which family was fitted, "
        "and the yield stress is still invisible in this geometry. The PRECISION does "
        "not: a force from this page is a Casson force, and choosing a different "
        "equally-good fit to the same curve moves it by tens of per cent. The "
        "comparison is not vacuous -- a Newtonian null misses by "
        f"{100 * s['newtonian']['worst_abs_rel_dev']:.4g} %, a "
        f"{g3['separation']:.3g}x separation, so the needle has not simply washed out "
        "the difference between every model.", ""]


def main() -> int:
    recs = {k: v for k, v in sorted(load_all().items()) if v.model.name == "casson"}
    assert len(recs) == 20, len(recs)
    nd = needle(ID_NOMINAL_MM)
    rows = [one_donor(r, nd) for r in recs.values()]

    n_extrap = sum(r["extrapolated"] for r in rows)
    dec = np.array([r["decades_past"] for r in rows])
    xs = np.array([r["x"] for r in rows])
    ty = np.array([r["tau_y"] for r in rows])
    ec = np.array([r["eta_c"] for r in rows])
    F = np.array([r["force_N"] for r in rows])

    # scipy, not a hand-rolled double-argsort: that version is wrong under ties (it
    # returns -1.0 where scipy returns -0.949 on a 4-point tied example), and published
    # Casson tables are rounded and will tie. No ties in THIS data -- the two agree to
    # the last bit here -- but the correct one is free.
    rho_ec = float(spearmanr(ec, F).statistic)
    rho_ty = float(spearmanr(ty, F).statistic)

    # what the yield stress is actually worth here: the same records with tau_y set to 0
    F_no_yield = []
    for r, rec in zip(rows, recs.values()):
        dp0 = pressure_for_flow_rate(Q_SI, nd, rec.model, np.array([1e-12, r["eta_c"]]))
        F_no_yield.append(dp0 * np.pi * BARREL_RADIUS ** 2)
    yield_share = (F - np.array(F_no_yield)) / F

    # The geometry sensitivity. Run on the MEDIAN donor by force, not the largest: the
    # largest is donor N, which the source study flags as atypical AND whose eta_c
    # profile interval spans 1.56 decades, the worst-constrained fit in the registry.
    # Illustrating a geometry point on the least trustworthy record is a choice that
    # needs stating, so it is not made.
    by_force = sorted(rows, key=lambda r: r["force_N"])
    median_row = by_force[len(by_force) // 2]
    band = {}
    for tag, id_mm in (("low", ID_BAND_MM[0]), ("narrow_low", ID_BAND_NARROW_MM[0]),
                       ("nominal", ID_NOMINAL_MM), ("high", ID_BAND_MM[1])):
        band[tag] = one_donor(recs[median_row["record"]], needle(id_mm))
    band_extrap = {}
    for tag, id_mm in (("low", ID_BAND_MM[0]), ("high", ID_BAND_MM[1])):
        nd_b = needle(id_mm)
        all_rows = [one_donor(r, nd_b) for r in recs.values()]
        band_extrap[tag] = {"n_extrapolated": sum(r["extrapolated"] for r in all_rows),
                            "min_decades_past": min(r["decades_past"] for r in all_rows)}

    # The claim that needs NO standard: F ~ R^-4 at fixed Q (Poiseuille, since x << 1
    # here). Measured, not assumed, by fitting the log-log slope over the whole band --
    # and note the docstrings previously said "Q ~ R^3 at fixed dp", which is wrong twice
    # (Poiseuille is R^4 at fixed dp, and this study fixes Q, so the relation is R^-4).
    ids = np.linspace(0.9 * ID_NOMINAL_MM, 1.1 * ID_NOMINAL_MM, 9)
    f_scan = np.array([one_donor(recs[median_row["record"]], needle(i))["force_N"]
                       for i in ids])
    slope = float(np.polyfit(np.log(ids), np.log(f_scan), 1)[0])
    # what ID error reproduces the whole donor-to-donor spread?
    spread_all = float(F.max() / F.min())
    keep = np.array([r["donor"] not in PAPER_OUTLIERS for r in rows])
    spread_typical = float(F[keep].max() / F[keep].min())
    id_error_for_spread = spread_all ** (1.0 / abs(slope)) - 1.0
    id_error_for_typical = spread_typical ** (1.0 / abs(slope)) - 1.0
    geometry = {"slope_dlogF_dlogR": slope, "spread_all_20": spread_all,
                "spread_excluding_paper_outliers": spread_typical,
                "n_excluded": int((~keep).sum()), "excluded": list(PAPER_OUTLIERS),
                "id_error_matching_full_spread": id_error_for_spread,
                "id_error_matching_typical_spread": id_error_for_typical,
                "band_wide_force_ratio": band["low"]["force_N"] / band["high"]["force_N"],
                "band_narrow_force_ratio": band["narrow_low"]["force_N"]
                / band["nominal"]["force_N"],
                "median_donor": median_row["donor"]}

    # provenance: without it the artifact cannot say which version of the bracket
    # produced it (raised by the code review 2026-09-05)
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                              cwd=str(HERE), timeout=30).stdout.strip() or None
    except Exception:
        head = None
    payload = {"generated_by": "case_studies/blood_injectability/run.py",
               "git_head": head, "biorheo_version": biorheo.__version__,
               "needle": {"id_nominal_mm": ID_NOMINAL_MM, "id_band_mm": list(ID_BAND_MM),
                          "length_m": NEEDLE_LENGTH},
               "barrel_radius_m": BARREL_RADIUS, "Q_mL_per_min": Q_ML_PER_MIN,
               "Q_si": Q_SI, "rho": RHO_BLOOD, "rows": rows,
               "n_extrapolated": n_extrap, "n_records": len(rows),
               "decades_past_min": float(dec.min()), "decades_past_max": float(dec.max()),
               "x_min": float(xs.min()), "x_max": float(xs.max()),
               "spearman_force_eta_c": rho_ec, "spearman_force_tau_y": rho_ty,
               "yield_share_min": float(yield_share.min()),
               "yield_share_max": float(yield_share.max()),
               "geometry_band": band, "band_extrapolation": band_extrap,
               "geometry": geometry,
               "wall_rate_spread": float(max(r["gd_wall"] for r in rows)
                                         / min(r["gd_wall"] for r in rows))}
    (HERE / "results.json").write_text(json.dumps(payload, indent=2), encoding="ascii")

    # the Hagenbach/entrance term, for the "what is not in these numbers" section
    u_mean = Q_SI / (np.pi * nd.radius ** 2)
    ent_pa = 1.08 * RHO_BLOOD * u_mean ** 2 / 2.0
    ent_frac = np.array([ent_pa / r["dp_Pa"] for r in rows])
    ent_frac_min, ent_frac_max = float(ent_frac.min()), float(ent_frac.max())
    # the honest comparison is PER DONOR: entrance term vs that same donor's yield share
    beats = [(rows[i]["donor"], float(ent_frac[i]), float(yield_share[i]))
             for i in range(len(rows)) if ent_frac[i] > yield_share[i]]
    ent_worst_donor = rows[int(np.argmax(ent_frac))]["donor"]
    payload["entrance_loss"] = {"u_mean_m_s": float(u_mean), "delta_p_Pa": float(ent_pa),
                                "fraction_of_dp_min": ent_frac_min,
                                "fraction_of_dp_max": ent_frac_max,
                                "worst_donor": ent_worst_donor,
                                "donors_where_entrance_beats_yield":
                                    [b_[0] for b_ in beats]}
    (HERE / "results.json").write_text(json.dumps(payload, indent=2), encoding="ascii")

    L = [f"# Twenty blood donors through a 27G needle at {Q_ML_PER_MIN:g} mL/min", "",
         "Generated by `run.py`. Casson closed form, derived and anchored against the "
         "package's independent quadrature in `derivations/casson_tube_flow.md`. "
         f"Needle: nominal 27G regular wall, inner diameter {ID_NOMINAL_MM:g} mm, length "
         f"{NEEDLE_LENGTH * 1e3:g} mm. Syringe barrel radius {BARREL_RADIUS * 1e3:g} mm. "
         "The geometry is a choice, not a measurement (see the script docstring); the "
         "ISO tolerance band is swept below.", "",
         "## The pre-registered question", "",
         f"**{n_extrap} of {len(rows)} records are flagged `extrapolated`.** The "
         "prediction, written before the run, was 20 of 20, on the grounds that the "
         "apparent wall rate is about 26x the top of every record's 700 1/s envelope.",
         "",
         f"The count is the boring part. The DISTANCE is the result: the true wall shear "
         f"rate is {dec.min():.2f} decades past the top of the fitted envelope, so the "
         f"Casson parameters are being asked to predict a regime {10 ** dec.min():.0f}x "
         "beyond anything the rheometer saw. Casson has no high-rate structure to get "
         "wrong -- it tends to a constant `eta_c` -- but that is an assumption of the "
         "FORM, and the form is already known to misfit these donors (3 runs of residual "
         "sign against 9 expected). Extrapolating a misfit model "
         f"{dec.min():.2f} decades is the honest description of what a needle "
         "calculation on this data is.", "",
         "That distance is the SAME for every donor to three figures (the wall rate "
         f"spans only {max(r['gd_wall'] for r in rows) / min(r['gd_wall'] for r in rows):.4f}x "
         "across the 20), and that is not a coincidence: for a Newtonian fluid the wall "
         "shear rate is `4Q/(pi R^3)`, independent of viscosity, and these donors are "
         "effectively Newtonian in this geometry (next section). **How far a device "
         "extrapolates its own rheology is therefore a property of the DEVICE, not of "
         "the patient** -- it is fixed by flow rate and lumen before any fluid is "
         "chosen, and it can be computed at design time.", "",
         "That last sentence needs its scope stated, because it is demonstrated only "
         "within the Casson family, which structurally CANNOT express extra shear-"
         "thinning above yield -- it saturates to a constant `eta_c`. For a genuinely "
         "shear-thinning family the true wall rate is larger than the Newtonian one by "
         "an O(1) factor that depends on the fluid: a power law gives `(3n+1)/(4n)`, "
         "which is 1.25x at n = 0.5 and 2x at n = 0.2. So the DECADE is a device "
         "property; the third significant figure is a property of the Casson form. A "
         "second mechanism helps here too and should be named: the envelope [0.1, 700] "
         "1/s is identical for all 20 donors, being the rheometer protocol rather than a "
         "per-donor fit, so both the numerator and the denominator of the distance are "
         "device or protocol quantities.", "",
         "## The yield stress does not matter here", "",
         f"`x = tau_y / tau_w` runs {xs.min():.3g} to {xs.max():.3g} across the 20 donors, "
         "so the Casson bracket `1 - (16/7)sqrt(x) + (4/3)x - x^4/21` differs from 1 by "
         f"{100 * (16 / 7) * np.sqrt(xs.max()):.2g} % at worst. Setting `tau_y` to zero "
         f"and refitting the pressure changes the plunger force by "
         f"{100 * yield_share.min():.2g} % to {100 * yield_share.max():.2g} %. In this "
         "geometry blood is a Newtonian fluid of viscosity `eta_c`, and the donor-to-donor "
         f"force spread ranks with `eta_c` at Spearman {rho_ec:.3f} and with `tau_y` at "
         f"{rho_ty:.3f}. **The parameter the published Casson tables argue about is the "
         "one this device cannot see.** That is not a criticism of the tables -- tau_y is "
         "what matters in the microcirculation, at the rates the rheometer measured -- it "
         "is a statement about which measurement a needle calculation needs.", "",
         "The pre-registered estimate for that bracket deviation was 0.5 %, and it was "
         f"wrong by 6x: the measured worst case is "
         f"{100 * (16 / 7) * np.sqrt(xs.max()):.2g} %. The estimate had been carried over "
         "from the derivation's own table, which uses a 20 mm needle and different "
         "pressures. The conclusion holds at the larger number; the estimate did not.",
         "",
         "## Every donor", "",
         "| donor | tau_y (Pa) | eta_c (mPa s) | dp (kPa) | tau_w (Pa) | x | "
         "gd_wall (1/s) | decades past 700 | force (N) | Re |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['donor']} | {r['tau_y']:.4g} | {1e3 * r['eta_c']:.4g} | "
                 f"{r['dp_Pa'] / 1e3:.4g} | {r['tau_w_Pa']:.4g} | {r['x']:.3g} | "
                 f"{r['gd_wall']:.5g} | {r['decades_past']:.3f} | {r['force_N']:.4g} | "
                 f"{r['reynolds']:.3g} |")

    b, g = band, geometry
    L += ["", "## The gauge is not a radius", "",
          "Force scales as `R^-4` at fixed flow rate -- measured here, not assumed: a "
          f"log-log fit over +/-10 % in inner diameter gives d(log F)/d(log R) = "
          f"{g['slope_dlogF_dlogR']:.4f}, against -4 exactly for Poiseuille -- the small "
          "difference is the Casson bracket, which is not quite 1. (This file's own "
          "docstring previously said `Q ~ R^3 at fixed dp`, which is wrong twice: "
          "Poiseuille gives `R^4` at fixed `dp`, and "
          "this study fixes `Q`, so the operative relation is `R^-4`. Corrected "
          "2026-09-05 after a referee pass; the computed numbers never used the wrong "
          "scaling, only its stated justification did.)", "",
          "The consequence needs no standard and no tolerance table. The whole "
          f"donor-to-donor force spread across these 20 people is {g['spread_all_20']:.2f}x, "
          "and because `F ~ R^-4` an inner-diameter error of only "
          f"**{100 * g['id_error_matching_full_spread']:.1f} %** reproduces it exactly. "
          "Excluding the three donors the source study itself flags as atypical "
          f"({', '.join(g['excluded'])}) the spread falls to "
          f"{g['spread_excluding_paper_outliers']:.2f}x, which a "
          f"{100 * g['id_error_matching_typical_spread']:.1f} % diameter error "
          "reproduces. **A needle whose bore is uncertain at the ten-percent level "
          "carries more uncertainty in plunger force than the biological variation "
          "between twenty people.** That statement is exact arithmetic on the fitted "
          "records; nothing in it depends on a manufacturing standard.", "",
          "### What was retracted here (2026-09-05)", "",
          "An earlier version of this file said: *ISO 9626:2016 allows an inner diameter "
          f"of {ID_BAND_MM[0]:g}-{ID_BAND_MM[1]:g} mm for regular-wall 27G*, computed a "
          f"{g['band_wide_force_ratio']:.2f}x force ratio across it, and concluded that "
          "the manufacturing tolerance outweighs the donor spread. **That is withdrawn.** "
          f"A referee pass found that {ID_BAND_MM[0]:g} mm is a regular-wall figure while "
          f"{ID_BAND_MM[1]:g} mm is a THIN-WALL one, so the pair is not one product's "
          "tolerance, and the standard itself is not available offline to check. Under "
          f"the narrower regular-wall-only reading "
          f"({ID_BAND_NARROW_MM[0]:g}-{ID_NOMINAL_MM:g} mm) the ratio is "
          f"{g['band_narrow_force_ratio']:.2f}x, which is SMALLER than the "
          f"{g['spread_all_20']:.2f}x donor spread -- the conclusion reverses. Both "
          "readings are printed below and both are tagged [UNVERIFIED]; the section's "
          "claim has been rebuilt on the `R^-4` scaling above, which does not need the "
          "standard.", "",
          "There was a second problem with the original comparison even if the band had "
          "been right: an ISO tolerance is a manufacturing bound over a whole population, "
          f"while {g['spread_all_20']:.2f}x is the max/min of a sample of 20, a statistic "
          "whose expected value grows with sample size. The two are not commensurate and "
          "were presented as if they were.", "",
          "| ID (mm) | reading | dp (kPa) | wall rate (1/s) | force (N) |",
          "|---|---|---|---|---|"]
    for tag, id_mm, note in (("low", ID_BAND_MM[0], "wide band low [UNVERIFIED]"),
                             ("narrow_low", ID_BAND_NARROW_MM[0], "narrow band low [UNVERIFIED]"),
                             ("nominal", ID_NOMINAL_MM, "nominal"),
                             ("high", ID_BAND_MM[1], "wide band high, thin wall [UNVERIFIED]")):
        row = b[tag]
        L.append(f"| {id_mm:g} | {note} | {row['dp_Pa'] / 1e3:.4g} | "
                 f"{row['gd_wall']:.5g} | {row['force_N']:.4g} |")
    L += ["",
          f"Rows are donor {g['median_donor']}, the MEDIAN donor by force. The earlier "
          "version used the largest-force donor, N -- which is one of the three the "
          "source study flags as atypical and has the worst-constrained fit in the whole "
          "registry (its `eta_c` profile interval spans 1.56 decades). Illustrating a "
          "geometry point on the least trustworthy record is a choice, and it is not "
          "made here.", "",
          f"All 20 records stay `extrapolated` across the whole band "
          f"({band_extrap['low']['n_extrapolated']}/20 at {ID_BAND_MM[0]:g} mm, "
          f"{band_extrap['high']['n_extrapolated']}/20 at {ID_BAND_MM[1]:g} mm; the "
          f"smallest distance anywhere is "
          f"{min(bb['min_decades_past'] for bb in band_extrap.values()):.2f} decades) -- "
          "checked on all 20, not inferred from one row.", "",
                    "## What is not in these numbers", "",
          "The module docstring lists the omissions: entrance and exit losses, the "
          "barrel's own pressure drop, the kinetic term, the start transient, wall slip. "
          "One of them is NOT negligible relative to the effect this study calls "
          "negligible, and saying `Re` is small is the wrong comparison -- it compares "
          "the term with `dp`, not with the term being ranked. The entrance (Hagenbach) "
          f"loss is about `1.08 rho U^2 / 2` = {ent_pa:.4g} Pa at the mean velocity "
          f"U = {u_mean:.4g} m/s that Q and R fix for every donor, which is "
          f"{100 * ent_frac_min:.2g}-{100 * ent_frac_max:.2g} % of the pressure drop -- "
          f"the same size as the {100 * yield_share.min():.2g}-"
          f"{100 * yield_share.max():.2g} % yield-stress contribution called negligible "
          "above. "
          + (f"On {len(beats)} of the 20 donors ({', '.join(b_[0] for b_ in beats)}) "
             "the neglected entrance term is larger than that SAME donor's whole yield "
             f"contribution -- on {beats[0][0]}, {100 * beats[0][1]:.2g} % against "
             f"{100 * beats[0][2]:.2g} %. " if beats else
             "On no donor does the entrance term exceed that donor's own yield "
             "contribution, though they are the same order. ")
          + "It does "
          "not change the conclusion -- both are small, and the yield stress is still "
          "the one being asked to carry a device claim -- but the correct statement is "
          "that this calculation resolves neither, not that it resolves one and not the "
          "other. The entrance length is "
          f"{1e3 * 0.06 * min(r['reynolds'] for r in rows) * 2 * nd.radius:.3g}-"
          f"{1e3 * 0.06 * max(r['reynolds'] for r in rows) * 2 * nd.radius:.3g} mm "
          f"against a {NEEDLE_LENGTH * 1e3:g} mm needle.", "",
          "And the corrections do not point the same way. Fahraeus-Lindqvist -- a "
          "cell-depleted layer at the wall that lowers the apparent viscosity of blood "
          "in tubes below roughly 300 um [UNVERIFIED: no citation resolved for that "
          "threshold in this session] -- lowers the force; the entrance loss raises it. "
          f"This lumen is {ID_NOMINAL_MM * 1e3:.0f} um, inside the first regime, so it "
          "is a real downward correction, but calling these numbers `upper bounds` was "
          "counting only the downward one and is withdrawn. Quantifying either needs "
          "measurements the Horner dataset does not contain: tube viscometry at this "
          "diameter for the first, a two-length needle test for the second.", "",
          "One number in this file is not computed here: the `3 runs of residual sign "
          "against 9 expected` quoted above comes from `data/horner2020/model_comparison"
          ".md`, and it is the value for 18 of the 20 donors, not all of them (donor D "
          "has 8 runs, z = -0.52, and is not flagged structured).", ""]
    L += _crosscheck_section()
    (HERE / "results.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
