"""Item 25: 58 real high-concentration monoclonal antibody formulations through a needle.

**The first time anything in this project has been run on the fluid the venture is about.**
`OFFER.md` §6 names it as gap 3: *"No protein or high-concentration formulation curve, so
nothing here has been demonstrated on the wedge's own fluid."*

Source: Zenodo 14409284, CC BY 4.0, supporting Table S1 of Prasnikar et al., Int. J. Pharm.
674 (2025) 125451. 58 formulations of one model mAb at 159.3-212.2 mg/mL, each with a
measured viscosity and its SD. Provenance and the source's own words:
`data/mab_prasnikar2025/SOURCE.md`. Gates pre-registered in
`data/mab_injectability_prereg.md` BEFORE this file existed.

The fluid is Newtonian by the SOURCE's own verdict over a stated rate range (960-3329 1/s),
not by assumption -- which is also why these cannot become registry records: a one-parameter
model fitted to a single viscosity has zero residual degrees of freedom and no
identifiability sidecar can be written for it. G25.7 demonstrates that rather than asserting
it. The study therefore reads the deposited table directly, which is the arriving-data path
a paid study would actually take.

    python case_studies/mab_injectability/run.py

Writes `results.json` and `results.md` beside this file. Exit non-zero if a gate whose
verdict is not a recorded failure comes out other than predicted.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import openpyxl

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SRC = ROOT / "data" / "mab_prasnikar2025" / "raw" / "Table_S1.xlsx"

from biorheo.apps.injectability import (Needle, Syringe,      # noqa: E402
                                        injection_force, wall_shear_rate)
from biorheo.constitutive import get_model                     # noqa: E402
from biorheo.identify import fit_rate_curve                    # noqa: E402

NEWTONIAN = get_model("newtonian")


def _load_blood_module():
    """The blood study's `run.py`, loaded BY PATH under a distinct name.

    A bare `import run` would import THIS file -- both case studies are called `run.py`, and
    the first one on `sys.path` wins. It imported itself silently until a type checker said
    so, which is the kind of defect that produces a confident wrong answer rather than an
    error: every device constant would have resolved to this module's own, unset names.
    """
    import importlib.util
    path = ROOT / "case_studies" / "blood_injectability" / "run.py"
    spec = importlib.util.spec_from_file_location("blood_injectability_run", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["blood_injectability_run"] = mod
    spec.loader.exec_module(mod)
    return mod


# --- the device: taken VERBATIM from the blood case study ------------------------------
# Reusing it introduces no new unverified constant and makes the two studies directly
# comparable, which is the point -- this repo has never put its wedge fluid and its
# demonstration fluid on the same axis. Read from that module so they cannot drift.
blood = _load_blood_module()

NEEDLE = Needle(radius=blood.ID_NOMINAL_MM * 1e-3 / 2.0, length=blood.NEEDLE_LENGTH)
SYRINGE = Syringe(blood.BARREL_RADIUS)
Q_SI = blood.Q_SI
Q_ML_PER_MIN = blood.Q_ML_PER_MIN

# --- what the SOURCE states, quoted in SOURCE.md and not computed here -----------------
# Methods 2.2.4: "we measured the viscosity of highly concentrated protein solutions at
# shear rates of 960-3329 s^-1 ... categorized them as Newtonian. Consequently, samples
# were measured at 25 C and shear rates within the described range".
NEWTONIAN_RANGE = (960.0, 3329.0)     # 1/s
TEMPERATURE_C = 25.0

G252_TOL = 0.20        # the two series' Mooney slopes must agree to 20 %
G251_R2_BAR = 0.70     # below this, the reduction % column does not reproduce
G254_DECADE_BAR = 1.5  # predicted extrapolation, decades


def read_table() -> list[dict]:
    """The 58 formulations, paired ROW BY ROW from the `Measurements` sheet."""
    wb = openpyxl.load_workbook(SRC, data_only=True)
    rows = list(wb["Measurements"].iter_rows(values_only=True))
    hdr = [("" if c is None else str(c)).strip() for c in rows[1]]
    need = ["Label", "Test compound added", "c (VRA) [mM]", "c (mAb) [mg/ml]",
            "Average viscosity [mPas]", "Viscosity_SD [mPas]", "Viscosity reduction [%]"]
    idx = {}
    for name in need:
        if name not in hdr:
            raise ValueError(f"column {name!r} is not in the sheet; header is {hdr[:10]}")
        idx[name] = hdr.index(name)
    out, skipped = [], 0
    for r in rows[2:]:
        if r[idx["Label"]] is None:
            skipped += 1
            continue
        rec = {k: r[i] for k, i in idx.items()}
        # a row contributes only if BOTH quantities the force needs are present
        if not isinstance(rec["c (mAb) [mg/ml]"], (int, float)) or \
           not isinstance(rec["Average viscosity [mPas]"], (int, float)):
            raise ValueError(f"row {rec['Label']!r} is missing a concentration or a "
                             "viscosity; it cannot be silently dropped")
        out.append({
            "label": str(rec["Label"]).strip(),
            "compound": str(rec["Test compound added"]).strip(),
            "c_vra_mM": rec["c (VRA) [mM]"],
            "c_mab_mg_per_mL": float(rec["c (mAb) [mg/ml]"]),
            "eta_mPa_s": float(rec["Average viscosity [mPas]"]),
            "eta_sd_mPa_s": (float(rec["Viscosity_SD [mPas]"])
                             if isinstance(rec["Viscosity_SD [mPas]"], (int, float))
                             else float("nan")),
            "reduction_pct": (float(rec["Viscosity reduction [%]"])
                              if isinstance(rec["Viscosity reduction [%]"], (int, float))
                              else float("nan")),
            # the label suffix "200 mAb" marks the second, higher-concentration series
            "series": "B_200mAb" if "200 mAb" in str(rec["Label"]) else "A",
        })
    print(f"  read {len(out)} formulations, skipped {skipped} blank rows")
    return out


def mooney_slope(rows) -> tuple[float, float, int]:
    """`b` in `eta_base = a exp(b c)`, back-solved from the source's own reduction column."""
    c = np.array([r["c_mab_mg_per_mL"] for r in rows])
    eta = np.array([r["eta_mPa_s"] for r in rows])
    red = np.array([r["reduction_pct"] for r in rows])
    good = np.isfinite(red) & (red < 100.0)
    base = eta[good] / (1.0 - red[good] / 100.0)
    lb, cc = np.log(base), c[good]
    A = np.vstack([np.ones(cc.size), cc]).T
    coef, *_ = np.linalg.lstsq(A, lb, rcond=None)
    pred = A @ coef
    r2 = float(1.0 - np.sum((lb - pred) ** 2) / np.sum((lb - lb.mean()) ** 2))
    return float(coef[1]), r2, int(good.sum())


def force_of(eta_Pa_s: float) -> float:
    return injection_force(Q_SI, NEEDLE, SYRINGE, NEWTONIAN,
                           np.array([eta_Pa_s])).force


def main() -> int:  # noqa: C901
    print("=== reading the deposited table ===")
    rows = read_table()
    for r in rows:
        r["eta_Pa_s"] = r["eta_mPa_s"] * 1e-3

    # ---- G25.1: does the source's own reduction column reproduce? ----------------------
    print("\n=== G25.1: can the source's `reduction %` column be reproduced? ===")
    _b_all, r2_all, n_all = mooney_slope(rows)
    series = {s: [r for r in rows if r["series"] == s] for s in ("A", "B_200mAb")}
    per = {s: mooney_slope(v) for s, v in series.items()}
    for s, (b, r2, n) in per.items():
        print(f"  series {s}: n = {n}, b = {b:.6g} mL/mg, R2 = {r2:.4f}")
    g251 = {
        "bar": f"a single Mooney exponential must reproduce the column with R2 >= "
               f"{G251_R2_BAR}",
        "r2_all_rows": r2_all, "n_all_rows": n_all,
        "r2_main_series": per["A"][1], "n_main_series": per["A"][2],
        "pass": bool(per["A"][1] >= G251_R2_BAR),
        "predicted_failure": True,
        "consequence": "this study uses the MEASURED viscosities and never the reduction "
                       "column. A number that cannot be reproduced cannot be extended.",
    }

    # ---- G25.2: the slope, which is what G25.5 needs ------------------------------------
    b_a, b_b = per["A"][0], per["B_200mAb"][0]
    b = 0.5 * (b_a + b_b)
    g252 = {
        "bar": f"the two series' Mooney slopes agree to within {G252_TOL:.0%}",
        "b_series_A": b_a, "b_series_B": b_b, "b_mean": b,
        "relative_difference": abs(b_a - b_b) / b,
        "doubling_concentration_mg_per_mL": math.log(2.0) / b,
        "pass": bool(abs(b_a - b_b) / b < G252_TOL),
        "note": "the scatter that sinks G25.1 is in the INTERCEPT, not the slope.",
    }
    print(f"  slopes {b_a:.6g} and {b_b:.6g} mL/mg, mean {b:.6g}; viscosity doubles every "
          f"{g252['doubling_concentration_mg_per_mL']:.1f} mg/mL")

    # ---- G25.3: the force, on the blood study's own device ------------------------------
    print(f"\n=== G25.3: force through the blood study's needle at {Q_ML_PER_MIN:g} mL/min "
          f"===")
    for r in rows:
        r["force_N"] = force_of(r["eta_Pa_s"])
        rel = (r["eta_sd_mPa_s"] / r["eta_mPa_s"]) if np.isfinite(r["eta_sd_mPa_s"]) else 0.0
        # a Newtonian force is exactly proportional to viscosity, so the relative error
        # passes through unchanged; no propagation approximation is involved.
        r["force_sd_N"] = r["force_N"] * rel
        r["wall_rate_s"] = wall_shear_rate(
            r["force_N"] / SYRINGE.area, NEEDLE, NEWTONIAN, np.array([r["eta_Pa_s"]]))
    f = np.array([r["force_N"] for r in rows])
    blood_forces = blood_force_reference()
    g253 = {
        "bar": "report the mAb forces against the blood donors' through the SAME device; "
               "if any mAb formulation needs less force than the median blood donor the "
               "prediction is refuted",
        "device": {"needle_id_mm": blood.ID_NOMINAL_MM,
                   "needle_length_m": blood.NEEDLE_LENGTH,
                   "barrel_radius_m": blood.BARREL_RADIUS,
                   "Q_mL_per_min": Q_ML_PER_MIN},
        "force_N_min": float(f.min()), "force_N_max": float(f.max()),
        "force_N_median": float(np.median(f)),
        "blood_force_N_median": blood_forces["median"],
        "blood_force_N_min": blood_forces["min"], "blood_force_N_max": blood_forces["max"],
        "ratio_median_mab_to_median_blood": float(np.median(f) / blood_forces["median"]),
        "n_mab_below_median_blood": int((f < blood_forces["median"]).sum()),
        "pass": bool((f < blood_forces["median"]).sum() == 0),
        "note": "no external injectable-force limit is used. This project retracted an ISO "
                "9626 comparison once for quoting a standard from memory; a human-factors "
                "force limit would be the same mistake.",
    }
    print(f"  mAb force {f.min():.3g}-{f.max():.3g} N, median {np.median(f):.3g} N")
    print(f"  blood median {blood_forces['median']:.3g} N -> ratio "
          f"{g253['ratio_median_mab_to_median_blood']:.1f}x")

    # ---- G25.4: how far the device extrapolates ------------------------------------------
    wall = np.array([r["wall_rate_s"] for r in rows])
    dec = np.log10(wall / NEWTONIAN_RANGE[1])
    g254 = {
        "bar": f"report the extrapolation in decades past the source's own measured "
               f"{NEWTONIAN_RANGE[0]:g}-{NEWTONIAN_RANGE[1]:g} 1/s, and compare with the "
               "other fluids in this project",
        "measured_range_s": list(NEWTONIAN_RANGE),
        "wall_rate_min_s": float(wall.min()), "wall_rate_max_s": float(wall.max()),
        "decades_past_measured_min": float(dec.min()),
        "decades_past_measured_max": float(dec.max()),
        "comparison_decades": {"blood_horner2020 (item 5)": 1.42,
                               "carbopol_pg_vinutha2024 (item 24, to 1e3 1/s)": None,
                               "dmem_sints2024 (item 22, to 1 1/s)": None},
        "pass": bool(dec.max() < G254_DECADE_BAR),
        "predicted": "under 1.5 decades and the smallest of any fluid here",
    }
    print(f"\n=== G25.4: wall rate {wall.min():.4g}-{wall.max():.4g} 1/s, i.e. "
          f"{dec.min():.3f}-{dec.max():.3f} decades past the measured {NEWTONIAN_RANGE[1]:g} "
          f"1/s ===")

    # ---- G25.5: the headline -- which lever, and the confound controlled -----------------
    print("\n=== G25.5: the excipient lever against the concentration lever ===")
    eta = np.array([r["eta_mPa_s"] for r in rows])
    c = np.array([r["c_mab_mg_per_mL"] for r in rows])
    raw_ratio = float(eta.max() / eta.min())
    # normalise every formulation to a reference concentration through the Mooney slope, or
    # the spread mixes the two levers: viscosity is exponential in c and these formulations
    # span 159-212 mg/mL.
    c_ref = float(np.median(c))
    eta_norm = eta * np.exp(-b * (c - c_ref))
    norm_ratio = float(eta_norm.max() / eta_norm.min())
    for r, e in zip(rows, eta_norm):
        r["eta_normalised_mPa_s"] = float(e)
    g255 = {
        "bar": "no bar; this is the measurement the study exists to produce",
        "reference_concentration_mg_per_mL": c_ref,
        "raw_viscosity_ratio": raw_ratio,
        "concentration_normalised_viscosity_ratio": norm_ratio,
        "mooney_slope_mL_per_mg": b,
        "excipient_lever_in_mg_per_mL_raw": math.log(raw_ratio) / b,
        "excipient_lever_in_mg_per_mL_normalised": math.log(norm_ratio) / b,
        "concentration_span_mg_per_mL": float(c.max() - c.min()),
        "doubling_concentration_mg_per_mL": g252["doubling_concentration_mg_per_mL"],
        "best": min(rows, key=lambda r: r["eta_normalised_mPa_s"])["label"],
        "worst": max(rows, key=lambda r: r["eta_normalised_mPa_s"])["label"],
        "pass": True,
    }
    print(f"  raw viscosity spread {raw_ratio:.2f}x -> "
          f"{g255['excipient_lever_in_mg_per_mL_raw']:.1f} mg/mL of concentration")
    print(f"  at MATCHED concentration ({c_ref:.1f} mg/mL) it is {norm_ratio:.2f}x -> "
          f"{g255['excipient_lever_in_mg_per_mL_normalised']:.1f} mg/mL")
    print(f"  the formulations themselves span {g255['concentration_span_mg_per_mL']:.1f} "
          f"mg/mL of concentration")

    # ---- G25.6: measurement against formulation ------------------------------------------
    sd = np.array([r["force_sd_N"] for r in rows])
    have = np.isfinite(sd) & (sd > 0)
    pairs = list(combinations(range(len(rows)), 2))
    disjoint = 0
    for i, j in pairs:
        if not (have[i] and have[j]):
            continue
        lo_i, hi_i = f[i] - 2 * sd[i], f[i] + 2 * sd[i]
        lo_j, hi_j = f[j] - 2 * sd[j], f[j] + 2 * sd[j]
        if max(lo_i, lo_j) > min(hi_i, hi_j):
            disjoint += 1
    frac = disjoint / len(pairs)
    g256 = {
        "bar": "report the fraction of formulation pairs whose force intervals are disjoint "
               "at 2 SD; item 14's blood answer was 36 of 190 (18.9 %)",
        "n_pairs": len(pairs), "n_disjoint": disjoint, "fraction_disjoint": frac,
        "blood_comparison": {"n_pairs": 190, "n_disjoint": 36, "fraction": 36 / 190},
        "median_relative_force_sd": float(np.median(sd[have] / f[have])),
        "pass": bool(frac > 36 / 190),
        "predicted": "a much higher fraction than blood's 18.9 %",
    }
    print(f"\n=== G25.6: {disjoint}/{len(pairs)} pairs separable at 2 SD ({100*frac:.1f} %), "
          f"against blood's 18.9 % ===")

    # ---- G25.7: demonstrate the schema boundary -------------------------------------------
    print("\n=== G25.7: can one of these become a registry record? ===")
    r0 = rows[0]
    try:
        fit_rate_curve(NEWTONIAN, np.array([2000.0]), np.array([r0["eta_Pa_s"]]),
                       temperature_C=TEMPERATURE_C, geometry="RheoSense m-VROC")
        outcome, detail = "SUCCEEDED", ""
    except Exception as exc:                                    # noqa: BLE001
        outcome, detail = "FAILED", f"{type(exc).__name__}: {exc}"
    print(f"  fitting a one-parameter model to one point: {outcome}")
    if detail:
        print(f"    {detail[:160]}")
    g257 = {
        "bar": "the failure must be SHOWN with its exception, not argued",
        "attempted_on": r0["label"], "n_points": 1,
        "outcome": outcome, "exception": detail,
        "pass": bool(outcome == "FAILED"),
        "note": "a Newtonian model has one parameter and each formulation has one "
                "viscosity: zero residual degrees of freedom, so there is nothing for "
                "write_identifiability to write and the loader would reject the record. "
                "This is also what a CLIENT's data usually looks like.",
    }

    gates = {"G25_1_the_reduction_column_does_not_reproduce": g251,
             "G25_2_the_mooney_slope_is_robust": g252,
             "G25_3_force_against_the_blood_study_same_device": g253,
             "G25_4_how_far_the_device_extrapolates": g254,
             "G25_5_which_lever_excipient_or_concentration": g255,
             "G25_6_measurement_against_formulation": g256,
             "G25_7_these_cannot_be_registry_records": g257}

    head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    payload = {"git_head": head, "n_formulations": len(rows), "rows": rows, "gates": gates,
               "newtonian_range_s": list(NEWTONIAN_RANGE), "temperature_C": TEMPERATURE_C,
               "prereg": "data/mab_injectability_prereg.md"}
    (HERE / "results.json").write_text(json.dumps(payload, indent=1), encoding="utf-8")
    write_md(payload)

    print("\n=== gates ===")
    unexpected = []
    for k, g in gates.items():
        ok = g["pass"]
        pf = g.get("predicted_failure", False)
        print(f"  {k}: {'PASS' if ok else 'FAIL'}"
              + ("  (predicted failure -- this is the result)" if pf else ""))
        if (ok and pf) or (not ok and not pf):
            unexpected.append(k)
    print(f"\nwrote {HERE / 'results.json'} and results.md")
    if unexpected:
        print(f"UNEXPECTED: {unexpected}")
        return 1
    return 0


def blood_force_reference() -> dict:
    """The blood donors' forces through the same device, from that study's own artifact."""
    art = json.loads((ROOT / "case_studies" / "blood_injectability" / "results.json")
                     .read_text(encoding="utf-8"))
    forces = [r["force_N"] for r in art["rows"] if isinstance(r.get("force_N"), (int, float))]
    if len(forces) < 10:
        raise ValueError(f"expected the 20 blood donors' forces, found {len(forces)}")
    a = np.array(forces)
    return {"n": len(forces), "median": float(np.median(a)),
            "min": float(a.min()), "max": float(a.max())}


def write_md(p) -> None:
    g, rows = p["gates"], p["rows"]
    L = []
    A = L.append
    A(f"# {p['n_formulations']} high-concentration monoclonal antibody formulations "
      "through a needle")
    A("")
    A("Written by `case_studies/mab_injectability/run.py`. Pre-registration: "
      "`data/mab_injectability_prereg.md`. Source and its own words: "
      "`data/mab_prasnikar2025/SOURCE.md`. Every number here is printed by that script.")
    A("")
    A("**This is the first result in this project computed on the fluid the venture is "
      "actually about.** Everything else here is human blood, a Carbopol microgel, or cell "
      "culture medium.")
    A("")
    g3, g4, g5 = (g["G25_3_force_against_the_blood_study_same_device"],
                  g["G25_4_how_far_the_device_extrapolates"],
                  g["G25_5_which_lever_excipient_or_concentration"])
    A("## The headline")
    A("")
    A(f"At matched mAb concentration the best and worst formulations differ in viscosity by "
      f"**{g5['concentration_normalised_viscosity_ratio']:.2f}x**. Through the Mooney slope "
      f"that is worth **{g5['excipient_lever_in_mg_per_mL_normalised']:.0f} mg/mL of mAb "
      f"concentration** — i.e. the right excipient buys you that much extra concentration at "
      f"the same viscosity, against formulations that themselves span "
      f"{g5['concentration_span_mg_per_mL']:.0f} mg/mL. Viscosity doubles every "
      f"**{g5['doubling_concentration_mg_per_mL']:.1f} mg/mL**.")
    A("")
    A("| lever | size | in mg/mL of mAb |")
    A("|---|---|---|")
    A(f"| excipient, at matched concentration | "
      f"{g5['concentration_normalised_viscosity_ratio']:.2f}x in viscosity | "
      f"{g5['excipient_lever_in_mg_per_mL_normalised']:.0f} |")
    A(f"| excipient, raw (CONFOUNDED — mixes both levers) | "
      f"{g5['raw_viscosity_ratio']:.2f}x | "
      f"{g5['excipient_lever_in_mg_per_mL_raw']:.0f} |")
    A(f"| the formulations' own concentration span | — | "
      f"{g5['concentration_span_mg_per_mL']:.0f} |")
    A("")
    A("The raw row is reported because it is what a spreadsheet gives you and it is the "
      "wrong number: these formulations are not at one concentration, viscosity is "
      "exponential in concentration, so the raw spread mixes the excipient with the "
      "concentration. The normalised row is the answer.")
    A("")
    A("## G25.3 — the force, on the blood study's own device")
    A("")
    d = g3["device"]
    A(f"Needle inner diameter {d['needle_id_mm']:g} mm, length {1e3*d['needle_length_m']:g} "
      f"mm, barrel radius {1e3*d['barrel_radius_m']:g} mm, {d['Q_mL_per_min']:g} mL/min — "
      "reused verbatim from `case_studies/blood_injectability`, so no new unverified "
      "constant enters and the two fluids sit on one axis for the first time.")
    A("")
    A("| | min | median | max |")
    A("|---|---|---|---|")
    A(f"| mAb formulations (N = {p['n_formulations']}) | {g3['force_N_min']:.3g} N | "
      f"{g3['force_N_median']:.3g} N | {g3['force_N_max']:.3g} N |")
    A(f"| blood donors, same device | {g3['blood_force_N_min']:.3g} N | "
      f"{g3['blood_force_N_median']:.3g} N | {g3['blood_force_N_max']:.3g} N |")
    A("")
    A(f"The median mAb formulation needs **{g3['ratio_median_mab_to_median_blood']:.1f}x** "
      f"the force of the median blood donor; "
      f"{g3['n_mab_below_median_blood']} of {p['n_formulations']} fall below it.")
    A("")
    A("**No injectable-force limit is quoted, from a standard or anywhere else.** This "
      "project retracted an ISO 9626 comparison once for quoting a standard from memory, and "
      "a human-factors force limit would be the same mistake. The scale here is this repo's "
      "own blood study.")
    A("")
    A("## G25.4 — how far the device extrapolates past what was measured")
    A("")
    A(f"The source measured, and verified Newtonian behaviour over, "
      f"**{g4['measured_range_s'][0]:g}–{g4['measured_range_s'][1]:g} 1/s**. At this "
      f"delivery rate the needle wall sees "
      f"**{g4['wall_rate_min_s']:.4g}–{g4['wall_rate_max_s']:.4g} 1/s**, i.e. "
      f"**{g4['decades_past_measured_max']:.2f} decades** past the top of the measured "
      f"range. Verdict: **{'PASS' if g4['pass'] else 'FAIL'}** against the pre-registered "
      f"1.5-decade bar.")
    A("")
    A("| fluid | decades of extrapolation |")
    A("|---|---|")
    A(f"| **this mAb** | **{g4['decades_past_measured_max']:.2f}** |")
    A("| blood (item 5) | 1.42 |")
    A("| Carbopol in propylene glycol (item 24) | measured to 1e3 1/s |")
    A("| DMEM + FBS (item 22) | measured to 1 1/s, four to five decades |")
    A("")
    A("*How far a device extrapolates its own rheology is a property of the DEVICE* — item "
      "5's finding, now applied to the fluid the venture is for.")
    A("")
    g6 = g["G25_6_measurement_against_formulation"]
    A("## G25.6 — the measurement is not the limiting term here, unlike blood")
    A("")
    A(f"Median relative force uncertainty from the measured viscosity SD: "
      f"**{100*g6['median_relative_force_sd']:.2f} %**. "
      f"**{g6['n_disjoint']} of {g6['n_pairs']}** formulation pairs "
      f"({100*g6['fraction_disjoint']:.1f} %) have disjoint force intervals at 2 SD.")
    A("")
    A(f"Item 14 asked the same question of the blood donors and got "
      f"{g6['blood_comparison']['n_disjoint']} of {g6['blood_comparison']['n_pairs']} "
      f"({100*g6['blood_comparison']['fraction']:.1f} %) — there the uncertainty swamped the "
      "differences and two donors had to be 17 places apart to be told apart. **Here the "
      "formulations are separable**, so a formulation decision made on these numbers is "
      "supported by them in a way a patient-ranking decision on the blood data was not.")
    A("")
    g1 = g["G25_1_the_reduction_column_does_not_reproduce"]
    A("## G25.1 — a column of the source that does NOT reproduce, and is therefore not used")
    A("")
    A("The paper defines its `viscosity reduction [%]` as the relative decrease against the "
      "starting formulation *at the same mAb concentration*, computed from an exponential "
      "(Mooney) calibration curve. Back-solving that baseline and fitting one exponential "
      f"gives **R² = {g1['r2_main_series']:.4f}** on the {g1['n_main_series']} rows of the "
      f"main series ({g1['r2_all_rows']:.4f} over all {g1['n_all_rows']}). A single "
      "calibration curve would collapse this to numerical precision.")
    A("")
    A("Either more than one curve was used, or the column carries processing the Methods do "
      "not describe. **The consequence is the point: this study uses the measured "
      "viscosities and never the reduction column.** A number that cannot be reproduced "
      "cannot be extended.")
    A("")
    g2 = g["G25_2_the_mooney_slope_is_robust"]
    A(f"What does survive is the **slope**: {g2['b_series_A']:.5g} mL/mg on the main series "
      f"and {g2['b_series_B']:.5g} on the higher-concentration one, agreeing to "
      f"{100*g2['relative_difference']:.1f} %. The scatter that sinks G25.1 is in the "
      "intercept. That slope is what makes the headline's mg/mL conversion possible.")
    A("")
    g7 = g["G25_7_these_cannot_be_registry_records"]
    A("## G25.7 — why these are NOT registry records, demonstrated rather than asserted")
    A("")
    A(f"Fitting a one-parameter Newtonian model to the single measured viscosity of "
      f"`{g7['attempted_on']}` **{g7['outcome']}**:")
    A("")
    A("```")
    A(g7["exception"] or "(no exception raised)")
    A("```")
    A("")
    A("Zero residual degrees of freedom, so there is no profile interval, nothing for "
      "`write_identifiability` to write, and the loader would reject the record. This is a "
      "real boundary of the schema — and it is worth stating plainly because **it is also "
      "what a client's data usually looks like**: a table of formulation viscosities, not a "
      "flow curve. This study therefore reads the deposited table directly, which is the "
      "arriving-data path a paid study would actually take.")
    A("")
    A("## Every formulation")
    A("")
    A("| formulation | excipient | c(mAb) mg/mL | η mPa·s | ± SD | force N | η at matched c |")
    A("|---|---|---|---|---|---|---|")
    for r in sorted(rows, key=lambda r: r["eta_normalised_mPa_s"]):
        sd = ("%.3f" % r["eta_sd_mPa_s"]) if np.isfinite(r["eta_sd_mPa_s"]) else "—"
        A(f"| {r['label']} | {r['compound']} | {r['c_mab_mg_per_mL']:.1f} | "
          f"{r['eta_mPa_s']:.2f} | {sd} | {r['force_N']:.3g} | "
          f"{r['eta_normalised_mPa_s']:.2f} |")
    A("")
    A("## What this page does NOT say")
    A("")
    A("- **Nothing about whether any of these is injectable by a human.** No force limit is "
      "quoted. The scale is this repo's own blood study on the same device.")
    A("- **Nothing about stability.** The source carries SEC monomer and aggregate columns; "
      "stability is not rheology and this package has nothing to say about it. An excipient "
      "that fixes viscosity and ruins the protein is not a win, and that trade-off is in the "
      "source and out of scope here.")
    A("- **Nothing that uses the source's `viscosity reduction` column.** G25.1.")
    A("- **Nothing about shear-thinning at the needle wall.** The source verified Newtonian "
      f"over {g4['measured_range_s'][0]:g}–{g4['measured_range_s'][1]:g} 1/s and that is the "
      "model used; whether it still holds "
      f"{g4['decades_past_measured_max']:.2f} decades higher is what that number is for, and "
      "it is not assumed either way.")
    (HERE / "results.md").write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
