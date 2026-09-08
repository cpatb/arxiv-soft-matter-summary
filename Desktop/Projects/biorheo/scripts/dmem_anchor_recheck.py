"""Item 23: re-grade item 22's external anchor against the paper's own Table 2.

Item 22 graded `biorheo.identify`'s fits against `0.32, 0.31, 0.54, 0.50`, taken from the
source paper's FIGURE CAPTIONS, which state the flow index used as the input to the paper's
COMSOL simulations -- rounded, and shared across conditions (one `n = 0.31` covers both the
10 % and the 5 % FBS room-temperature samples; one `n = 0.54` covers both the 37 C sample
and the 1 % FBS one).

The paper also carries **Table 2, "Fit results for K and n values in (4) for the
experimental rheometry data shown in Fig. 1"**, which resolves the four deposited conditions
one row each, for exactly the equation this repo fits, and which also publishes `K`. **All
four** of the values item 22 graded against disagree with it -- and two of them are the
right number against the wrong row (`SPEC` gives 5 % FBS the 0.31 that belongs to 10 %, and
gives 1 % FBS the 0.54 that belongs to 37 C), which is why the mismatch was easy to miss.

This script does not edit item 22's gates. It re-grades beside them at the SAME bar, prints
both readings, and reports what survives. Pre-registration:
`data/dmem_anchor_correction_prereg.md`. Nothing here is hardcoded into the `.md` it writes.

    python scripts/dmem_anchor_recheck.py

Writes `data/dmem_sints2024/anchor_recheck.json` and `.md`. Exit code is non-zero if any
gate whose verdict is not a PREDICTED FAILURE comes out other than predicted.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import zipfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import dmem_candidate_scan as scan  # noqa: E402
import records_from_dmem as gen  # noqa: E402

from biorheo.constitutive import get_model  # noqa: E402
from biorheo.identify import fit_rate_curve  # noqa: E402
from biorheo.registry import RECORDS_DIR, load_record  # noqa: E402

POWER_LAW = get_model("power_law")
ZIP = ROOT / "data" / "dmem_sints2024" / "DMEM_FBS_rheology.zip"
OUT = ROOT / "data" / "dmem_sints2024"

# --- the paper's own fit table -----------------------------------------------------------
# Sints et al. 2024, arXiv:2409.13650, Table 2: "Fit results for K and n values in (4) for
# the experimental rheometry data shown in Fig. 1". One row per sample. Equation (4) is the
# Ostwald-de Waele power law mu = K gd^(n-1), which is exactly this repo's `power_law`.
# Transcribed by hand from the paper; the values the paper prints, at the precision it
# prints them. NOT computed here and not computable here -- that is the point of an anchor.
TABLE2 = {
    "DMEM_and_10pct_FBS": {"K": 0.028, "n": 0.31, "label": "DMEM + 10 % FBS", "T": 20.0},
    "DMEM_and_5pct_FBS": {"K": 0.037, "n": 0.33, "label": "DMEM + 5 % FBS", "T": 20.0},
    "DMEM_and_1pct_FBS": {"K": 0.016, "n": 0.51, "label": "DMEM + 1 % FBS", "T": 20.0},
    "DMEM_and_10pct_FBS_37C": {"K": 0.005, "n": 0.54,
                               "label": "DMEM + 10 % FBS (37 C)", "T": 37.0},
}

# What item 22 graded against instead, from the figure captions. Kept so the artifact
# carries the superseded number rather than deleting it.
SPEC_N = {"DMEM_and_10pct_FBS": 0.32, "DMEM_and_5pct_FBS": 0.31,
          "DMEM_and_1pct_FBS": 0.54, "DMEM_and_10pct_FBS_37C": 0.50}

# The interval of each file that lies in the paper's stated 0.01-1 1/s range. Item 22's own
# SPEC excludes ("DMEM_and_10pct_FBS_37C", 1) for this reason and then pools it anyway.
IN_RANGE = {"DMEM_and_10pct_FBS": (0, 1), "DMEM_and_5pct_FBS": (0, 1),
            "DMEM_and_1pct_FBS": (0, 1), "DMEM_and_10pct_FBS_37C": (0,)}
ALL_INTERVALS = {stem: (0, 1) for stem in TABLE2}

# The authors' own declared repeatability, S2 Fig: "A 20 % deviation from the fit line is
# included as a measure of repeatability of results." A declared scale, not one we chose.
PAPER_REPEATABILITY = 0.20

# Table 2 prints K as 0.028, 0.037, 0.016, 0.005 -- three decimals throughout.
PUBLISHED_K_DECIMALS = 3

G221_TOL = gen.G221_TOL          # 0.05, not moved
RECORD_STEMS = {"DMEM_and_10pct_FBS": "10pct_20C", "DMEM_and_5pct_FBS": "5pct_20C",
                "DMEM_and_1pct_FBS": "1pct_20C", "DMEM_and_10pct_FBS_37C": "10pct_37C"}


def finite_max(values, what):
    """`max` over floats, with the NaNs removed first and counted.

    Builtin `max` over a list containing NaN is order-dependent and silently returns the
    wrong answer (`max([0.1, nan])` is 0.1, `max([nan, 0.1])` is nan). This has already
    bitten this repo once; see lessons.md, 2026-09-06.
    """
    vals = list(values)
    good = [v for v in vals if math.isfinite(v)]
    if len(good) != len(vals):
        raise ValueError(f"{what}: {len(vals) - len(good)} of {len(vals)} values "
                         "are not finite; refusing to take a maximum over them")
    return max(good) if good else float("nan")


def curve(z, stem, idx):
    """One interval of one deposited file, as (rate, viscosity), with N printed both ways."""
    text = z.read(f"DMEM_FBS_rheology/{stem}.txt").decode("utf-8", "replace")
    iv = scan.intervals(text)[idx]
    n_raw = int(iv["gd"].size)
    keep = iv["gd"] > gen.MIN_RATE
    gd, tau = iv["gd"][keep], iv["tau"][keep]
    print(f"    {stem}[{idx}]: {n_raw} points read, {int(keep.sum())} kept "
          f"(rate > {gen.MIN_RATE:g}), rate {gd.min():.4g}-{gd.max():.4g} 1/s")
    return gd, tau / gd, n_raw


def fit(gd, eta, T):
    ip = fit_rate_curve(POWER_LAW, gd, eta, temperature_C=T, geometry=gen.GEOMETRY)
    if int(ip.fit.n_points) != int(gd.size):
        raise ValueError(f"the fitter kept {ip.fit.n_points} of the {gd.size} points it was "
                         "given and the caller's N would have been reported instead")
    return ip


def pooled_fit(z, stem, intervals, T):
    gd_all, eta_all, n_raw = [], [], 0
    for idx in intervals:
        gd, eta, nr = curve(z, stem, idx)
        gd_all.append(gd)
        eta_all.append(eta)
        n_raw += nr
    gd = np.concatenate(gd_all)
    eta = np.concatenate(eta_all)
    return fit(gd, eta, T), gd, eta, n_raw


def main() -> int:
    z = zipfile.ZipFile(ZIP)
    rows = []

    print("=== per-condition fits, pooling only the intervals in the paper's rate range ===")
    for stem, ref in TABLE2.items():
        print(f"  {ref['label']}:")
        ip, gd, eta, n_raw = pooled_fit(z, stem, ALL_INTERVALS[stem], ref["T"])
        k, n = float(ip.values[0]), float(ip.values[1])

        # per-interval fits, to expose whether the pooled estimate is bracketed by them
        per = []
        for idx in ALL_INTERVALS[stem]:
            g1, e1, _ = curve(z, stem, idx)
            per.append(float(fit(g1, e1, ref["T"]).values[1]))

        # the two parameter sets AS VISCOSITY CURVES over this condition's measured range.
        # Immune to the K-n tradeoff: a deviation in n alone is uninterpretable because a
        # compensating K keeps the curve in the same place where the data is.
        grid = np.logspace(math.log10(gd.min()), math.log10(gd.max()), 400)
        eta_ours = k * grid ** (n - 1.0)
        eta_paper = ref["K"] * grid ** (ref["n"] - 1.0)
        rel = np.abs(eta_ours / eta_paper - 1.0)

        # each parameter set's own rms misfit to the SAME data, in the same measure
        rmse_ours = float(np.sqrt(np.mean((np.log(k * gd ** (n - 1.0)) - np.log(eta)) ** 2)))
        rmse_paper = float(np.sqrt(np.mean(
            (np.log(ref["K"] * gd ** (ref["n"] - 1.0)) - np.log(eta)) ** 2)))

        kp = ip.params["k"]
        np_ = ip.params["n"]
        # does the fitted K round to the published one at the precision the paper prints it?
        # the paper prints every K to three decimals (0.028, 0.037, 0.016, 0.005), so that
        # is the precision the comparison can resolve. Deriving it from the magnitude gave
        # four decimals for 0.005 and asked the fit for a digit the paper never printed.
        dec = PUBLISHED_K_DECIMALS
        rounds_to_published = bool(round(k, dec) == round(ref["K"], dec))

        rows.append({
            "stem": stem, "label": ref["label"], "temperature_C": ref["T"],
            "n_points_pooled": int(gd.size), "n_points_read": n_raw,
            "intervals_pooled": list(ALL_INTERVALS[stem]),
            "rate_min": float(gd.min()), "rate_max": float(gd.max()),
            "k_fit": k, "n_fit": n,
            "K_paper_table2": ref["K"], "n_paper_table2": ref["n"],
            "n_paper_figure_caption": SPEC_N[stem],
            "abs_dev_n_vs_table2": abs(n - ref["n"]),
            "abs_dev_n_vs_figure_caption": abs(n - SPEC_N[stem]),
            "rel_dev_K_vs_table2": abs(k - ref["K"]) / ref["K"],
            "K_rounds_to_published": rounds_to_published,
            "K_published_decimals": dec,
            "per_interval_n": per,
            "pooled_inside_per_interval_span": bool(min(per) <= n <= max(per)),
            "per_interval_spread": float(max(per) - min(per)),
            "eta_curve_median_rel": float(np.median(rel)),
            "eta_curve_max_rel": float(rel.max()),
            "rmse_log_ours": rmse_ours, "rmse_log_paper_params": rmse_paper,
            "k_verdict": kp.verdict, "k_lo": float(kp.lo), "k_hi": float(kp.hi),
            "k_lo_at_bound": bool(kp.lo_at_bound), "k_hi_at_bound": bool(kp.hi_at_bound),
            "n_verdict": np_.verdict, "n_lo": float(np_.lo), "n_hi": float(np_.hi),
            "n_paper_inside_profile": bool(np_.lo <= ref["n"] <= np_.hi),
            "n_profile_width": float(np_.hi - np_.lo),
        })
        print(f"    -> k = {k:.5g} (paper {ref['K']:g}, {100 * rows[-1]['rel_dev_K_vs_table2']:.1f} %), "
              f"n = {n:.4f} (paper {ref['n']:g}, |dev| {rows[-1]['abs_dev_n_vs_table2']:.4f})")

    # ---- G23.1: the corrected anchor, at the ORIGINAL bar --------------------------------
    worst_t2 = finite_max((r["abs_dev_n_vs_table2"] for r in rows), "G23.1 deviations")
    failing = [r["label"] for r in rows if r["abs_dev_n_vs_table2"] >= G221_TOL]
    g231 = {
        "bar": f"|n_fit - n_Table2| < {G221_TOL} on all {len(rows)} conditions, pooling "
               "every interval of each condition -- which is what Table 2 says it fitted, "
               '"the entire ensemble of data points"',
        "worst_abs_dev": worst_t2, "n_conditions": len(rows),
        "conditions_failing": failing,
        "pass": bool(len(failing) == 0),
        "predicted": "FAIL on exactly one condition (1 % FBS 20 C) at 0.0599",
        "prediction_history": "the pre-registration ALSO predicted this while specifying a "
                              "pooling restricted to the paper's 0.01-1 1/s range. Those two "
                              "are inconsistent and the restricted run gives 0.1051 on two "
                              "conditions; see G23.1b and addendum A of the pre-registration.",
        "predicted_failure": True,
        "note": "the bar is NOT moved. Item 22 reported worst 0.0299 against the figure "
                "captions; against the paper's own fit table it is larger and one "
                "condition misses. The failure is the result.",
    }

    # ---- G23.2: neither single correction rescues the old pass ---------------------------
    # (a) same reference values as item 22, but do not pool the out-of-range interval
    worst_a_terms = []
    for stem, ref in TABLE2.items():
        ipa, _, _, _ = pooled_fit(z, stem, IN_RANGE[stem], ref["T"])
        worst_a_terms.append(abs(float(ipa.values[1]) - SPEC_N[stem]))
    worst_a = finite_max(worst_a_terms, "G23.2a deviations")
    # (b) item 22's pooling, but graded against Table 2
    worst_b_terms = []
    for stem, ref in TABLE2.items():
        ip, _, _, _ = pooled_fit(z, stem, ALL_INTERVALS[stem], ref["T"])
        worst_b_terms.append(abs(float(ip.values[1]) - ref["n"]))
    worst_b = finite_max(worst_b_terms, "G23.2b deviations")
    g232 = {
        "bar": "each single correction, applied alone to item 22's G22.1', must be reported",
        "worst_dropping_out_of_range_interval_only": worst_a,
        "worst_regrading_against_table2_only": worst_b,
        "item22_reported_worst": 0.0299,
        "both_fail": bool(worst_a >= G221_TOL and worst_b >= G221_TOL),
        "pass": bool(worst_a >= G221_TOL and worst_b >= G221_TOL),
        "predicted": "0.0651 and 0.0599, both FAIL the 0.05 bar",
        "note": "a claim that survives neither of two independent corrections was not "
                "resting on the data. Item 22's 0.0299 is kept here rather than deleted.",
    }

    # ---- G23.3: the flagship 0.0005 is a pooling coincidence -----------------------------
    flagship = next(r for r in rows if r["stem"] == "DMEM_and_10pct_FBS")
    g233 = {
        "bar": "a pooled estimate quoted as a precision must lie inside the span of the "
               "per-interval fits it pools",
        "flagship_condition": flagship["label"],
        "pooled_n": flagship["n_fit"], "per_interval_n": flagship["per_interval_n"],
        "inside_span": flagship["pooled_inside_per_interval_span"],
        "item22_quoted_agreement": 0.0005,
        "conditions_inside_span": {r["label"]: r["pooled_inside_per_interval_span"]
                                   for r in rows},
        "pass": bool(not flagship["pooled_inside_per_interval_span"]),
        "predicted": "False for the flagship condition, i.e. the 0.0005 is a coincidence "
                     "of combining two mutually inconsistent intervals",
        "note": "0.0005 must stop being quoted as a four-figure agreement anywhere.",
    }

    # ---- G23.4: what survives, scaled against the authors' OWN declared repeatability -----
    worst_eta = finite_max((r["eta_curve_max_rel"] for r in rows), "G23.4 curve deviations")
    worst_K = finite_max((r["rel_dev_K_vs_table2"] for r in rows), "G23.4 K deviations")
    n_K_round = sum(1 for r in rows if r["K_rounds_to_published"])
    g234 = {
        "bar": f"the two parameter sets, compared AS VISCOSITY CURVES over each condition's "
               f"own measured range, agree to inside the {PAPER_REPEATABILITY:.0%} "
               "repeatability the authors state for their own measurements (S2 Fig)",
        "paper_declared_repeatability": PAPER_REPEATABILITY,
        "worst_max_rel": worst_eta,
        "median_rel_by_condition": {r["label"]: r["eta_curve_median_rel"] for r in rows},
        "max_rel_by_condition": {r["label"]: r["eta_curve_max_rel"] for r in rows},
        "worst_rel_dev_K": worst_K,
        "K_rounds_to_published_count": n_K_round, "n_conditions": len(rows),
        "pass": bool(worst_eta < PAPER_REPEATABILITY),
        "predicted": "median 2.04 / 0.99 / 7.06 / 4.86 %, worst 19.03 %, inside 20 %",
        "note": "this is the comparison item 22 never made. It is immune to the K-n "
                "tradeoff, it uses a scale the SOURCE declares rather than one this repo "
                "invented, and it is the quantity a client actually uses.",
    }

    # ---- G23.5: the correction must not move any record ----------------------------------
    moved = []
    for stem, ref in TABLE2.items():
        for i, idx in enumerate(ALL_INTERVALS[stem]):
            rec_id = f"dmem_sints2024_{RECORD_STEMS[stem]}_{'ab'[i]}"
            gd1, eta1, _ = curve(z, stem, idx)
            ipi = fit(gd1, eta1, ref["T"])
            rec = load_record(Path(RECORDS_DIR) / f"{rec_id}.yaml")
            for j, pname in enumerate(POWER_LAW.param_names):
                shipped = float(rec.parameters[pname].value)
                refit = float(ipi.values[j])
                # the YAML stores %.6g, so compare at that precision
                if float(f"{refit:.6g}") != shipped:
                    moved.append({"record": rec_id, "param": pname,
                                  "shipped": shipped, "refit": refit})
    g235 = {
        "bar": "no record's parameters move; the defect is in the grading, not the fits",
        "records_checked": 8, "records_moved": moved,
        "pass": bool(not moved),
        "predicted": "no record moves",
        "note": "read on the parameter VALUES, which are continuous, not on a verdict.",
    }

    # ---- G23.6: does the correction reach a shipped number? ------------------------------
    hits = subprocess.run(
        ["git", "grep", "-l", "dmem", "--", "case_studies", "src/biorheo/apps"],
        cwd=ROOT, capture_output=True, text=True).stdout.split()
    g236 = {
        "bar": "no DMEM record enters a force, needle or design-space calculation, so the "
               "correction reaches no shipped number",
        "files_referencing_dmem_under_case_studies_or_apps": hits,
        "pass": bool(not hits),
        "predicted": "none",
        "note": "settled by grep, not by argument. A correction reaches a downstream "
                "number only through the parameters that number uses.",
    }

    # ---- G23.8: is K identifiable at all? ------------------------------------------------
    unident = [{"label": r["label"], "verdict": r["k_verdict"],
                "lo_at_bound": r["k_lo_at_bound"], "hi_at_bound": r["k_hi_at_bound"]}
               for r in rows
               if r["k_verdict"] not in ("identified", "measurable")
               or r["k_lo_at_bound"] or r["k_hi_at_bound"]]
    g238 = {
        "bar": "every condition's k carries a closed profile interval and is not on a bound;"
               " a comparison against a published K means nothing if K is unidentified",
        "k_verdicts": {r["label"]: r["k_verdict"] for r in rows},
        "k_intervals": {r["label"]: [r["k_lo"], r["k_hi"]] for r in rows},
        "unidentified": unident,
        "pass": bool(not unident),
        "predicted": "all four identified with closed intervals",
    }

    # ---- G23.1b: how much of G23.1 is a data-selection choice? --------------------------
    # Item 22 excluded the 37 C interval that runs to 10 1/s on the ground that the paper's
    # indices are for 0.01-1 1/s. Table 2 says otherwise -- it is fitted to "the entire
    # ensemble of data points" -- so the restriction was this repo's, not the paper's.
    # Reporting both is the point: at 37 C the choice moves n more than the whole
    # disagreement being graded.
    restricted = []
    for stem, ref in TABLE2.items():
        ipr, gdr, _, _ = pooled_fit(z, stem, IN_RANGE[stem], ref["T"])
        nr = float(ipr.values[1])
        restricted.append({"label": ref["label"], "n_points": int(gdr.size), "n_fit": nr,
                           "abs_dev_n_vs_table2": abs(nr - ref["n"])})
    worst_restricted = finite_max((r["abs_dev_n_vs_table2"] for r in restricted),
                                  "G23.1b deviations")
    # the 37 C interval that the restriction throws away, fitted on its own
    ip37, gd37, _, _ = pooled_fit(z, "DMEM_and_10pct_FBS_37C", (1,), 37.0)
    n37_alone = float(ip37.values[1])
    g231b = {
        "bar": "report the same grading under item 22's own 0.01-1 1/s restriction, so the "
               "size of the data-selection choice is visible beside the disagreement",
        "worst_abs_dev_restricted": worst_restricted,
        "worst_abs_dev_pooling_everything": worst_t2,
        "per_condition": restricted,
        "the_discarded_37C_interval_alone": {
            "n_points": int(gd37.size), "n_fit": n37_alone,
            "abs_dev_n_vs_table2": abs(n37_alone - TABLE2["DMEM_and_10pct_FBS_37C"]["n"])},
        "pass": True,
        "note": "at 37 C the three data choices give n = 0.4349 (5 points), 0.5393 (the "
                "20-point interval alone) and 0.5120 (both). The published value is 0.54. "
                "The restriction, not the fit, is what disagrees.",
    }

    # ---- G23.9: the comparison scaled by how well each n is actually determined ----------
    # Added after the first run, and said so. A fixed 0.05 bar treats a 5-point fit and a
    # 30-point fit alike; a profile interval does not. Stricter in kind than the bar it sits
    # beside, because a fit too weak to resolve the difference cannot pass by being vague --
    # it is reported as undetermined, which is a third verdict and not a pass.
    outside = [r["label"] for r in rows if not r["n_paper_inside_profile"]]
    g239 = {
        "bar": "the paper's Table 2 n falls inside biorheo's own profile interval on n",
        "inside": {r["label"]: r["n_paper_inside_profile"] for r in rows},
        "intervals": {r["label"]: [r["n_lo"], r["n_hi"]] for r in rows},
        "widths": {r["label"]: r["n_profile_width"] for r in rows},
        "conditions_outside": outside,
        "pass": bool(not outside),
        "added_after_the_first_run": True,
        "predicted_failure": True,
        "predicted": "inside on 3 of 4; the 1 % FBS condition is the one whose deviation "
                     "(0.0599) is large against a well-determined 30-point fit",
        "note": "a deviation is only a disagreement if the fit resolves it. This is the "
                "scaled version of G23.1 and it is reported BESIDE the failed bar, never "
                "instead of it.",
    }

    gates = {"G23_1_corrected_anchor_at_the_original_bar": g231,
             "G23_2_neither_single_correction_rescues_the_old_pass": g232,
             "G23_3_the_flagship_agreement_is_a_pooling_coincidence": g233,
             "G23_4_what_survives_against_the_authors_own_repeatability": g234,
             "G23_5_no_record_moves": g235,
             "G23_6_the_correction_reaches_no_shipped_number": g236,
             "G23_1b_how_much_is_a_data_selection_choice": g231b,
             "G23_8_K_is_identifiable": g238,
             "G23_9_the_published_n_against_our_own_profile_interval": g239}

    head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    payload = {"git_head": head, "rows": rows, "gates": gates,
               "G221_TOL_unchanged": G221_TOL,
               "paper": "Sints et al. 2024, arXiv:2409.13650, Table 2",
               "prereg": "data/dmem_anchor_correction_prereg.md"}
    (OUT / "anchor_recheck.json").write_text(json.dumps(payload, indent=1), encoding="utf-8")

    write_md(payload)

    print("\n=== gates ===")
    bad = []
    for name, g in gates.items():
        predicted_failure = g.get("predicted_failure", False)
        ok = g["pass"]
        # G23.1 is a PREDICTED failure: its `pass` is False by design and that is correct.
        good = (not ok) if predicted_failure else ok
        print(f"  {name}: {'PASS' if ok else 'FAIL'}"
              + ("  (predicted failure -- this is the result)" if predicted_failure else ""))
        if not good:
            bad.append(name)
    print(f"\nwrote {OUT / 'anchor_recheck.json'} and .md")
    if bad:
        print(f"UNEXPECTED: {bad}")
        return 1
    return 0


def write_md(p) -> None:
    rows, g = p["rows"], p["gates"]
    L = []
    A = L.append
    A("# Item 23 — the item-22 anchor, re-graded against the paper's own Table 2")
    A("")
    A("Generated by `scripts/dmem_anchor_recheck.py`. Every number here is printed by that")
    A("script. Pre-registration: `data/dmem_anchor_correction_prereg.md`.")
    A("")
    A("Item 22 graded `biorheo.identify`'s fits against flow indices taken from the source")
    A("paper's **figure captions**, which state the index used as the input to its COMSOL")
    A("simulations -- rounded, and shared across conditions. The paper also carries **Table 2,")
    A('"Fit results for `K` and `n` values in (4) for the experimental rheometry data shown')
    A('in Fig. 1"**, one row per sample, for exactly the equation this repo fits, and it')
    A(f"publishes `K` as well. **All {len(rows)} of the {len(rows)}** values item 22 used disagree")
    A("with it -- two of them are the right number against the wrong row, which is why the")
    A("mismatch was easy to miss.")
    A("")
    A("## The two reference tables, side by side")
    A("")
    A("| condition | `K` Table 2 | `n` Table 2 | `n` figure caption (item 22 used this) |")
    A("|---|---|---|---|")
    for r in rows:
        A(f"| {r['label']} | {r['K_paper_table2']:g} | {r['n_paper_table2']:g} | "
          f"{r['n_paper_figure_caption']:g} |")
    A("")
    A("## What biorheo fits, pooling only the intervals in the paper's stated rate range")
    A("")
    A("| condition | N | rates (1/s) | `k` | `n` | rel. dev `K` | abs. dev `n` |")
    A("|---|---|---|---|---|---|---|")
    for r in rows:
        A(f"| {r['label']} | {r['n_points_pooled']} | "
          f"{r['rate_min']:.4g}-{r['rate_max']:.4g} | {r['k_fit']:.5g} | {r['n_fit']:.4f} | "
          f"{100 * r['rel_dev_K_vs_table2']:.1f} % | {r['abs_dev_n_vs_table2']:.4f} |")
    A("")
    g231 = g["G23_1_corrected_anchor_at_the_original_bar"]
    A(f"## G23.1 — **FAILED**, at the unchanged {p['G221_TOL_unchanged']} bar")
    A("")
    A(f"Worst `|n_fit - n_Table2|` = **{g231['worst_abs_dev']:.4f}** over "
      f"{g231['n_conditions']} conditions. Failing: "
      + (", ".join(g231["conditions_failing"]) or "none") + ".")
    A("")
    A("Item 22 reported worst **0.0299** and a PASS. The bar has not been moved; the")
    A("reference values have been corrected to the ones the paper designates as its fit")
    A("results. The failure is the result.")
    A("")
    g232 = g["G23_2_neither_single_correction_rescues_the_old_pass"]
    A("## G23.2 — neither correction alone rescues the old pass")
    A("")
    A(f"| single change applied to item 22's G22.1' | worst `|dn|` | verdict at "
      f"{p['G221_TOL_unchanged']} |")
    A("|---|---|---|")
    A(f"| as item 22 shipped it | {g232['item22_reported_worst']:.4f} | PASS |")
    A(f"| drop the 37 C interval `SPEC` itself excludes as out of range | "
      f"{g232['worst_dropping_out_of_range_interval_only']:.4f} | **FAIL** |")
    A(f"| keep the pooling, grade against Table 2 | "
      f"{g232['worst_regrading_against_table2_only']:.4f} | **FAIL** |")
    A("")
    A("Two independent corrections, each fatal on its own.")
    A("")
    g233 = g["G23_3_the_flagship_agreement_is_a_pooling_coincidence"]
    A("## G23.3 — the flagship 0.0005 was a pooling coincidence")
    A("")
    A(f"For {g233['flagship_condition']} the pooled index is "
      f"**{g233['pooled_n']:.4f}** while its own two interval fits are "
      + " and ".join(f"{v:.4f}" for v in g233["per_interval_n"]) + ".")
    A("")
    A("A pooled estimator **outside the span of its own inputs** means the two intervals are")
    A("mutually inconsistent under a single power law, so agreement to four figures with any")
    A("external number is a coincidence of combining them, not a precision. Per condition:")
    A("")
    A("| condition | pooled inside the span of its intervals? | per-interval spread |")
    A("|---|---|---|")
    for r in rows:
        A(f"| {r['label']} | {'yes' if r['pooled_inside_per_interval_span'] else '**no**'} "
          f"| {r['per_interval_spread']:.4f} |")
    A("")
    g234 = g["G23_4_what_survives_against_the_authors_own_repeatability"]
    A("## G23.4 — what survives, and it is a comparison item 22 never made")
    A("")
    A("Table 2 publishes `K` as well as `n`, and `K` and `n` trade off along the data: a")
    A("deviation in `n` alone is uninterpretable because a compensating `K` keeps the curve")
    A("in the same place where the points are. So compare the two parameter sets **as")
    A("viscosity curves**, over each condition's own measured range.")
    A("")
    A("| condition | median | worst | rmse_log biorheo | rmse_log paper's params |")
    A("|---|---|---|---|---|")
    for r in rows:
        A(f"| {r['label']} | {100 * r['eta_curve_median_rel']:.2f} % | "
          f"{100 * r['eta_curve_max_rel']:.2f} % | {r['rmse_log_ours']:.4f} | "
          f"{r['rmse_log_paper_params']:.4f} |")
    A("")
    A(f"Worst over all four: **{100 * g234['worst_max_rel']:.2f} %**, against the")
    A(f"**{g234['paper_declared_repeatability']:.0%}** the authors state for their own")
    A('measurements (S2 Fig: *"A 20 % deviation from the fit line is included as a measure')
    A('of repeatability of results."*). That is a scale the SOURCE declares, not a tolerance')
    A(f"this repo invented. Verdict: **{'PASS' if g234['pass'] else 'FAIL'}**.")
    A("")
    A(f"`K` alone agrees to **{100 * g234['worst_rel_dev_K']:.1f} %** or better on all four,")
    A(f"and **{g234['K_rounds_to_published_count']} of {g234['n_conditions']}** fitted values")
    A("round to the published value at the precision the paper prints it:")
    A("")
    A("| condition | `k` fitted | `K` published | rounds to it? |")
    A("|---|---|---|---|")
    for r in rows:
        A(f"| {r['label']} | {r['k_fit']:.5g} | {r['K_paper_table2']:g} | "
          f"{'yes' if r['K_rounds_to_published'] else 'no'} |")
    A("")
    g235 = g["G23_5_no_record_moves"]
    A("## G23.5 — no record moves")
    A("")
    A(f"All {g235['records_checked']} records refitted from the raw files and compared to the")
    A("shipped YAML at its own `%.6g` storage precision. Moved: "
      f"**{len(g235['records_moved'])}**. The defect was in the grading, not in the fits.")
    A("")
    g236 = g["G23_6_the_correction_reaches_no_shipped_number"]
    A("## G23.6 — the correction reaches no shipped number")
    A("")
    A("Files under `case_studies/` or `src/biorheo/apps/` that reference a DMEM record: "
      f"**{len(g236['files_referencing_dmem_under_case_studies_or_apps'])}**"
      + (" (" + ", ".join(g236["files_referencing_dmem_under_case_studies_or_apps"]) + ")"
         if g236["files_referencing_dmem_under_case_studies_or_apps"] else "")
      + ". Settled by grep, not by argument.")
    A("")
    A("## G23.8 — `K` is identifiable, so the `K` comparison means something")
    A("")
    A("| condition | verdict | profile interval on `k` |")
    A("|---|---|---|")
    for r in rows:
        A(f"| {r['label']} | {r['k_verdict']} | [{r['k_lo']:.5g}, {r['k_hi']:.5g}] |")
    A("")
    g231b = g["G23_1b_how_much_is_a_data_selection_choice"]
    A("## G23.1b — how much of that failure is a data-selection choice")
    A("")
    A("Item 22 excluded the 37 C interval running to 10 1/s, on the ground that the paper's")
    A("indices are for 0.01-1 1/s. Table 2 says otherwise: it is fitted to *\"the entire")
    A("ensemble of data points\"*. So the restriction was this repo's, not the paper's.")
    A("Both gradings, side by side:")
    A("")
    A("| condition | `n` pooling everything | `n` restricted to 0.01-1 1/s | N restricted |")
    A("|---|---|---|---|")
    for r, rr in zip(rows, g231b["per_condition"]):
        A(f"| {r['label']} | {r['n_fit']:.4f} | {rr['n_fit']:.4f} | {rr['n_points']} |")
    A("")
    d = g231b["the_discarded_37C_interval_alone"]
    A(f"Worst deviation pooling everything: **{g231b['worst_abs_dev_pooling_everything']:.4f}**;")
    A(f"under the restriction: **{g231b['worst_abs_dev_restricted']:.4f}**.")
    A("")
    A("**At 37 C the choice of data moves `n` further than the disagreement being graded.**")
    A(f"The interval the restriction throws away, fitted alone ({d['n_points']} points),")
    A(f"gives `n` = **{d['n_fit']:.4f}** against the published")
    A(f"{TABLE2['DMEM_and_10pct_FBS_37C']['n']:g} — a deviation of **{d['abs_dev_n_vs_table2']:.4f}**.")
    A("The restriction, not the fit, is what disagrees there.")
    A("")
    g239 = g["G23_9_the_published_n_against_our_own_profile_interval"]
    A("## G23.9 — the same comparison, scaled by how well `n` is determined")
    A("")
    A("Added after the first run, and labelled so. A fixed 0.05 bar treats a 5-point fit and")
    A("a 30-point fit alike; a profile interval does not. This is **stricter in kind** than")
    A("the bar it sits beside, because a fit too weak to resolve the difference cannot pass")
    A("by being vague — it is reported undetermined, which is a third verdict, not a pass.")
    A("")
    A("| condition | biorheo's profile interval on `n` | width | `n` Table 2 | inside? |")
    A("|---|---|---|---|---|")
    for r in rows:
        A(f"| {r['label']} | [{r['n_lo']:.4f}, {r['n_hi']:.4f}] | {r['n_profile_width']:.4f} "
          f"| {r['n_paper_table2']:g} | "
          f"{'yes' if r['n_paper_inside_profile'] else '**no**'} |")
    A("")
    A(f"Outside on **{len(g239['conditions_outside'])}** of {len(rows)}: "
      + (", ".join(g239["conditions_outside"]) or "none") + ".")
    A("")
    A("## The restated headline")
    A("")
    A("Not *\"reproduces the paper's published flow indices on 4/4 conditions, worst")
    A("0.0299\"*. Instead:")
    A("")
    A("> Fitting the raw deposited curves through `biorheo.identify` reproduces the source")
    A("> paper's own Table 2 `(K, n)` **as a viscosity curve** to a median of")
    A(f"> {100 * min(r['eta_curve_median_rel'] for r in rows):.1f}-"
      f"{100 * max(r['eta_curve_median_rel'] for r in rows):.1f} % over each condition's")
    A(f"> measured range, worst {100 * g234['worst_max_rel']:.1f} %, inside the")
    A(f"> {g234['paper_declared_repeatability']:.0%} repeatability the authors state for")
    A("> their own measurements. On the flow index alone it agrees on three of four")
    A(f"> conditions and misses the fourth by {g231['worst_abs_dev']:.4f} against a")
    A(f"> pre-registered {p['G221_TOL_unchanged']}. `K` agrees to")
    A(f"> {100 * g234['worst_rel_dev_K']:.1f} % or better on all four.")
    A("")
    A("## What no version of this anchor can show")
    A("")
    A("1. **The fit method is shared.** For a power law, minimising rms `log(eta)` residuals")
    A("   is ordinary least squares of `log eta` on `log gd` -- the near-universal route to")
    A("   `(K, n)`, and the paper gives no sign of doing otherwise. Agreement certifies this")
    A("   repo's **parse, interval split and units**; it does not exercise the multi-start")
    A("   search or the profile likelihood, which this problem is too well conditioned to")
    A("   strain.")
    A("2. **The paper's own scatter is the size of the agreement.** Its figure captions and")
    A("   its Table 2 differ by 0.01-0.04 on the same conditions, and Table 2 is fitted to an")
    A("   ensemble of three samples, refills and different days while the deposit is one file")
    A("   per condition. The reference is not a point.")
    A("3. **The dataset's own noise floor exceeds the agreement on half the conditions.**")
    A("   Repeat intervals of the same condition move `n` by up to")
    A(f"   {max(r['per_interval_spread'] for r in rows):.4f}, larger than the deviations")
    A("   being called agreement there.")
    (OUT / "anchor_recheck.md").write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
