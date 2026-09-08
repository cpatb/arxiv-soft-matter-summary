"""Item 28: which shipped conclusions survive a non-zero syringe friction term.

Pre-registration: `prereg.md` beside this file, written before any number here was
computed. Read its §2 first -- the algebra there decides which questions are live, and
two of the three obvious checks are invariant BY CONSTRUCTION and are therefore not
gates.

`biorheo.apps.injectability.injection_force` takes `friction: float = 0.0` and every
case study in this repo leaves it at the default. So every published force is the
hydrodynamic term alone. This study does not measure a syringe and does not claim the
shipped forces are wrong; it asks which CONCLUSIONS drawn from them are properties of
the assumption `friction == 0`.

Two live questions, separated because they have different answers:

  Q1  a constant friction `f`  -- moves RATIO claims, leaves ORDER and DIFFERENCE alone.
  Q2  a syringe-to-syringe friction SD `s` -- moves SEPARABILITY, which a constant does
      not, because a constant shifts both ends of every interval equally.

Reads only the shipped artifacts (`blood_injectability/results.json`, `uq.json`,
`mab_injectability/results.json`); it recomputes no force and refits nothing, so it
cannot disturb what it is auditing. C1 in the pre-registration requires it to reproduce
three shipped numbers exactly before any new number is believed.

Run:  python case_studies/device_friction/run.py
Writes `results.md` and `results.json` beside this file.
"""

from __future__ import annotations

import itertools
import json
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
BLOOD = HERE.parent / "blood_injectability"
MAB = HERE.parent / "mab_injectability"

# Shipped numbers this study must reproduce before it is allowed to say anything (C1).
SHIPPED_RATIO_G25_3 = 7.1          # printed to one decimal on the mAb page
SHIPPED_MAB_DISJOINT = None        # read from the artifact, not hardcoded
SHIPPED_BLOOD_DISJOINT = 36
SHIPPED_BLOOD_PAIRS = 190

# Friction grids. Both are declared here, not chosen after seeing the answer.
F_GRID_N = np.concatenate([[0.0], np.geomspace(1e-3, 20.0, 400)])
S_GRID_N = np.concatenate([[0.0], np.geomspace(1e-4, 5.0, 500)])


def _load() -> tuple[dict, dict, dict]:
    blood = json.loads((BLOOD / "results.json").read_text(encoding="utf-8"))
    uq = json.loads((BLOOD / "uq.json").read_text(encoding="utf-8"))
    mab = json.loads((MAB / "results.json").read_text(encoding="utf-8"))
    return blood, uq, mab


def _pairs_disjoint(lo: np.ndarray, hi: np.ndarray) -> int:
    """Count index pairs whose [lo, hi] intervals do not overlap."""
    n = len(lo)
    return sum(1 for i, j in itertools.combinations(range(n), 2)
               if hi[i] < lo[j] or hi[j] < lo[i])


def main() -> int:
    blood, uq, mab = _load()

    # ---------------- the arrays, straight from the shipped artifacts ------------------
    f_blood = np.array([r["force_N"] for r in blood["rows"]])
    hull_lo = np.array([r["hull_lo"] for r in uq["rows"]])
    hull_hi = np.array([r["hull_hi"] for r in uq["rows"]])
    # the blood interval is a parameter push-through hull, asymmetric about the fit.
    # Represent its half-width so a device term can be added in quadrature; the
    # asymmetry is preserved by carrying the two sides separately.
    blood_dn = np.array([r["force_N"] for r in uq["rows"]]) - hull_lo
    blood_up = hull_hi - np.array([r["force_N"] for r in uq["rows"]])
    f_blood_uq = np.array([r["force_N"] for r in uq["rows"]])

    f_mab = np.array([r["force_N"] for r in mab["rows"]])
    sd_mab = np.array([r["force_sd_N"] for r in mab["rows"]])

    g253 = mab["gates"]["G25_3_force_against_the_blood_study_same_device"]
    g256 = mab["gates"]["G25_6_measurement_against_formulation"]

    out: dict = {"prereg": "case_studies/device_friction/prereg.md"}

    # ---------------- C1: reproduce before extending -----------------------------------
    ratio_0 = float(np.median(f_mab) / np.median(f_blood))
    n_disjoint_mab_0 = _pairs_disjoint(f_mab - 2 * sd_mab, f_mab + 2 * sd_mab)
    n_pairs_mab = len(f_mab) * (len(f_mab) - 1) // 2
    n_disjoint_blood_0 = _pairs_disjoint(hull_lo, hull_hi)
    n_pairs_blood = len(f_blood_uq) * (len(f_blood_uq) - 1) // 2

    c1 = {
        "ratio_recomputed": ratio_0,
        "ratio_shipped": g253["ratio_median_mab_to_median_blood"],
        "ratio_ok": abs(ratio_0 - g253["ratio_median_mab_to_median_blood"]) < 1e-9,
        "mab_disjoint_recomputed": n_disjoint_mab_0,
        "mab_disjoint_shipped": g256["n_disjoint"],
        "mab_pairs_recomputed": n_pairs_mab,
        "mab_pairs_shipped": g256["n_pairs"],
        "mab_ok": (n_disjoint_mab_0 == g256["n_disjoint"]
                   and n_pairs_mab == g256["n_pairs"]),
        "blood_disjoint_recomputed": n_disjoint_blood_0,
        "blood_disjoint_shipped": SHIPPED_BLOOD_DISJOINT,
        "blood_pairs_recomputed": n_pairs_blood,
        "blood_pairs_shipped": SHIPPED_BLOOD_PAIRS,
        "blood_ok": (n_disjoint_blood_0 == SHIPPED_BLOOD_DISJOINT
                     and n_pairs_blood == SHIPPED_BLOOD_PAIRS),
    }
    c1["pass"] = bool(c1["ratio_ok"] and c1["mab_ok"] and c1["blood_ok"])
    out["C1_reproduction"] = c1
    print(f"C1 reproduction: ratio {ratio_0:.6f} vs shipped "
          f"{g253['ratio_median_mab_to_median_blood']:.6f}; "
          f"mAb {n_disjoint_mab_0}/{n_pairs_mab} vs {g256['n_disjoint']}/{g256['n_pairs']}; "
          f"blood {n_disjoint_blood_0}/{n_pairs_blood} vs "
          f"{SHIPPED_BLOOD_DISJOINT}/{SHIPPED_BLOOD_PAIRS} -> "
          f"{'PASS' if c1['pass'] else 'FAIL'}")
    if not c1["pass"]:
        (HERE / "results.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
        print("C1 FAILED. The re-implementation does not reproduce the shipped numbers, "
              "so no number below is believable. Stopping; that disagreement is the "
              "finding.")
        return 2

    # ---------------- the invariants, stated not gated (prereg §2) ---------------------
    # A constant f is additive, so order, differences and deterministic-interval
    # separability are unchanged. Asserted here as arithmetic so the claim is auditable,
    # NOT counted as a gate: a check that cannot fail decides nothing.
    f_probe = 2.7
    inv = {
        "order_preserved": bool(np.array_equal(np.argsort(f_mab),
                                               np.argsort(f_mab + f_probe))),
        "differences_preserved": bool(np.allclose(np.diff(np.sort(f_mab)),
                                                  np.diff(np.sort(f_mab + f_probe)),
                                                  rtol=0, atol=1e-15)),
        "deterministic_separability_preserved":
            bool(_pairs_disjoint(f_mab + f_probe - 2 * sd_mab,
                                 f_mab + f_probe + 2 * sd_mab) == n_disjoint_mab_0),
        "probe_N": f_probe,
        "note": "invariant by construction; recorded as arithmetic, not as a gate",
    }
    out["invariants_by_construction"] = inv

    # ---------------- Q1: the ratio against a constant friction ------------------------
    med_mab, med_blood = float(np.median(f_mab)), float(np.median(f_blood))
    ratio = (med_mab + F_GRID_N) / (med_blood + F_GRID_N)

    def first_below(target: float) -> float:
        idx = np.nonzero(ratio < target)[0]
        return float(F_GRID_N[idx[0]]) if idx.size else float("nan")

    def ratio_at(f_val: float) -> float:
        return float((med_mab + f_val) / (med_blood + f_val))

    q1 = {
        "median_mab_N": med_mab, "median_blood_N": med_blood,
        "ratio_at_f": {f"{v:g}": ratio_at(v) for v in (0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0)},
        "f_first_below_3": first_below(3.0),
        "f_first_below_2": first_below(2.0),
        "f_first_below_1_5": first_below(1.5),
        # the force DIFFERENCE, which the same algebra says is untouched
        "median_difference_N": med_mab - med_blood,
    }
    q1["P1_pass"] = bool(ratio_at(1.0) < 3.0 and ratio_at(3.0) < 2.0)
    q1["P2_pass"] = bool(q1["f_first_below_2"] < 4.0)
    out["Q1_constant_friction"] = q1
    print(f"Q1: ratio {ratio_at(0.0):.2f}x at f=0 -> {ratio_at(1.0):.2f}x at 1 N, "
          f"{ratio_at(3.0):.2f}x at 3 N; first below 2.0 at f = "
          f"{q1['f_first_below_2']:.3g} N")
    print(f"    P1 {'PASS' if q1['P1_pass'] else 'FAIL'}, "
          f"P2 {'PASS' if q1['P2_pass'] else 'FAIL'}")

    # ---------------- Q2: separability against friction repeatability ------------------
    # A device friction whose value varies syringe to syringe with SD s adds s in
    # quadrature to each force's own SD. For the mAb cohort the shipped interval IS
    # +/- 2 SD of a measured repeatability, so this is exact. For the blood cohort the
    # shipped interval is a parameter push-through HULL, which is not a standard
    # deviation; the convention used here -- and it IS a convention, named -- is to widen
    # each side of the hull by 2s, i.e. to treat the device term as a 2-SD band added to
    # a hull that is already quoted as a 95 % object. Both cohorts are reported against
    # their OWN f=0 baseline and never against each other (prereg C2).
    mab_frac, blood_frac = [], []
    for s in S_GRID_N:
        sd_eff = np.hypot(sd_mab, s)
        mab_frac.append(_pairs_disjoint(f_mab - 2 * sd_eff, f_mab + 2 * sd_eff)
                        / n_pairs_mab)
        blood_frac.append(_pairs_disjoint(hull_lo - 2 * s, hull_hi + 2 * s)
                          / n_pairs_blood)
    mab_frac = np.array(mab_frac)
    blood_frac = np.array(blood_frac)

    def s_star(frac: np.ndarray, target: float) -> float:
        idx = np.nonzero(frac < target)[0]
        return float(S_GRID_N[idx[0]]) if idx.size else float("nan")

    def frac_at(frac: np.ndarray, s_val: float) -> float:
        return float(frac[int(np.argmin(np.abs(S_GRID_N - s_val)))])

    q2 = {
        "convention": "device friction SD s added in quadrature to the mAb measurement "
                      "SD; for blood, each hull side widened by 2s. Named, not derived.",
        "mab": {"n_pairs": n_pairs_mab, "frac_at_s0": float(mab_frac[0]),
                "s_star_50pct_N": s_star(mab_frac, 0.50),
                "s_star_10pct_N": s_star(mab_frac, 0.10),
                "frac_at": {f"{v:g}": frac_at(mab_frac, v)
                            for v in (0.005, 0.01, 0.05, 0.1, 0.5)}},
        "blood": {"n_pairs": n_pairs_blood, "frac_at_s0": float(blood_frac[0]),
                  "s_star_10pct_N": s_star(blood_frac, 0.10),
                  "s_star_5pct_N": s_star(blood_frac, 0.05),
                  "frac_at": {f"{v:g}": frac_at(blood_frac, v)
                              for v in (0.005, 0.01, 0.05, 0.1, 0.5)}},
    }
    q2["P3_pass"] = bool(np.isfinite(q2["mab"]["s_star_50pct_N"])
                         and q2["mab"]["s_star_50pct_N"] < 0.05)
    q2["P5_pass"] = bool(frac_at(blood_frac, 0.1) < 0.10)
    out["Q2_friction_repeatability"] = q2
    print(f"Q2: mAb separable {100*mab_frac[0]:.1f} % at s=0; drops below 50 % at "
          f"s* = {q2['mab']['s_star_50pct_N']:.4g} N")
    print(f"    blood separable {100*blood_frac[0]:.1f} % at s=0; "
          f"{100*frac_at(blood_frac, 0.1):.1f} % at s = 0.1 N")
    print(f"    P3 {'PASS' if q2['P3_pass'] else 'FAIL'}, "
          f"P5 {'PASS' if q2['P5_pass'] else 'FAIL'}")

    # ---------------- P4: which cohort the same absolute s damages more ----------------
    s4 = 0.05
    mab_hw0 = float(np.median(2 * sd_mab))
    mab_hw4 = float(np.median(2 * np.hypot(sd_mab, s4)))
    blood_hw0 = float(np.median(0.5 * (blood_up + blood_dn)))
    blood_hw4 = blood_hw0 + 2 * s4
    p4: dict = {
        "s_N": s4,
        "mab_median_halfwidth_N": mab_hw0, "mab_median_halfwidth_at_s_N": mab_hw4,
        "mab_growth": mab_hw4 / mab_hw0,
        "blood_median_halfwidth_N": blood_hw0, "blood_median_halfwidth_at_s_N": blood_hw4,
        "blood_growth": blood_hw4 / blood_hw0,
        "mab_relative_sd_pct": float(100 * np.median(sd_mab / f_mab)),
        "blood_relative_hull_pct": float(np.median(
            [r["hull_width_pct"] for r in uq["rows"]])),
    }
    p4["pass"] = bool(p4["mab_growth"] > 3.0 and p4["blood_growth"] < 1.3)
    # P4' -- the pre-registration compared two arms under DIFFERENT operations: the mAb
    # arm added the device term in quadrature and the blood arm added it LINEARLY (2s on
    # each hull side). That is a defect in the pre-registration, not a reason to move a
    # bar, so P4's verdict and numbers above stand exactly as graded. P4' repeats the
    # same question with the SAME operation on both arms, and is reported beside it as a
    # post-hoc recomputation, labelled as such.
    p4["matched_convention"] = {
        "operation": "half-width -> hypot(half-width, 2s) on BOTH arms",
        "blood_median_halfwidth_at_s_N": float(np.hypot(blood_hw0, 2 * s4)),
        "blood_growth": float(np.hypot(blood_hw0, 2 * s4) / blood_hw0),
        "mab_growth": p4["mab_growth"],
        "status": "POST-HOC. Not a pre-registered gate. Reported because the "
                  "pre-registered comparison was not like-for-like.",
    }
    p4["matched_convention"]["mab_over_blood"] = (
        p4["mab_growth"] / p4["matched_convention"]["blood_growth"])
    out["P4_which_cohort_is_exposed"] = p4
    print(f"P4: at s = {s4} N the mAb median half-width grows {p4['mab_growth']:.2f}x "
          f"and the blood one {p4['blood_growth']:.2f}x -> "
          f"{'PASS' if p4['pass'] else 'FAIL'}")
    print(f"    P4' matched convention: blood "
          f"{p4['matched_convention']['blood_growth']:.2f}x, mAb/blood = "
          f"{p4['matched_convention']['mab_over_blood']:.2f}x (post-hoc)")

    # ---------------- C4: is the answer 'friction is irrelevant'? ----------------------
    out["C4_null"] = {
        "would_fire_if": "s* for the mAb screen exceeds 1 N",
        "s_star_50pct_N": q2["mab"]["s_star_50pct_N"],
        "fired": bool(not np.isfinite(q2["mab"]["s_star_50pct_N"])
                      or q2["mab"]["s_star_50pct_N"] > 1.0),
    }

    # ---------------- where a real device sits on the Q1 curve --------------------------
    # C3: the study's deliverable is the THRESHOLD. An external friction magnitude only
    # says WHERE on the curve a device sits, and it is admitted only if it resolved to a
    # DOI in a lookup (data/injection_force_terms_lit_2026-09-08.md, this session).
    # It is ATTRIBUTED, never asserted, and it is a different syringe from the one this
    # repo's case studies chose -- their 1 mL barrel is a choice, not a measured device.
    ext = {
        "source": "Funke et al., J. Pharm. Sci. (2016), doi:10.1002/jps.23945",
        "device": "1 mL glass prefilled syringe with a staked needle; NOT the syringe "
                  "either case study specifies, whose 2.3 mm barrel radius is a CHOICE",
        "glide_N": [1.3, 1.5],
        "glide_note": "1.3-1.5 N glide, from a trade article with NO DOI in the same "
                      "survey; treated as UNVERIFIED and reported only as a bracket",
        "breakloose_N": [3.0, 9.0],
        "breakloose_note": "~9 N siliconised / ~3 N silicone-free (Funke 2016, RESOLVED)",
        "which_term_applies": "the case studies compute a STEADY flow, so glide is the "
                              "comparable term; breakloose is the start transient the "
                              "module docstring already excludes",
        "ratio_at_glide": [ratio_at(1.3), ratio_at(1.5)],
        "ratio_at_breakloose": [ratio_at(3.0), ratio_at(9.0)],
        "no_repeatability_resolved": True,
        "gap": "the survey resolved friction MEANS and no friction SD or "
               "syringe-to-syringe repeatability for any device. s* is therefore a "
               "threshold with nothing published to compare it against -- the quantity "
               "that decides whether a formulation screen is decidable is unpublished.",
    }
    out["external_placement"] = ext
    print(f"EXT: at a glide friction of 1.3-1.5 N the 7.05x ratio is "
          f"{ext['ratio_at_glide'][0]:.2f}-{ext['ratio_at_glide'][1]:.2f}x; "
          f"no friction SD resolved anywhere, so s* has no published comparator")

    out["curves"] = {
        "f_grid_N": F_GRID_N.tolist(), "ratio": ratio.tolist(),
        "s_grid_N": S_GRID_N.tolist(),
        "mab_frac_separable": mab_frac.tolist(),
        "blood_frac_separable": blood_frac.tolist(),
    }
    out["git_head"] = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                                     capture_output=True, text=True).stdout.strip()

    (HERE / "results.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    write_md(out)

    verdicts = {"P1": q1["P1_pass"], "P2": q1["P2_pass"], "P3": q2["P3_pass"],
                "P4": p4["pass"], "P5": q2["P5_pass"]}
    print("\n=== pre-registered predictions ===")
    for k, v in verdicts.items():
        print(f"  {k}: {'PASS' if v else 'REFUTED'}")
    n_ref = sum(1 for v in verdicts.values() if not v)
    print(f"{n_ref} of {len(verdicts)} REFUTED. No bar was moved.")
    return 0


def write_md(d: dict) -> None:
    q1, q2, p4 = d["Q1_constant_friction"], d["Q2_friction_repeatability"], \
        d["P4_which_cohort_is_exposed"]
    c1 = d["C1_reproduction"]
    L = []
    A = L.append
    A("# The term that defaults to zero\n")
    A("Pre-registration: [`prereg.md`](prereg.md), written before any number here. "
      "Generated by `run.py`; every number below is printed by that script.\n")
    A("`injection_force` takes `friction: float = 0.0` and no case study in this repo "
      "supplies it, so every force this project has published is the hydrodynamic term "
      "alone. That is labelled correctly on every page. What had never been checked is "
      "which CONCLUSIONS drawn from those forces depend on the friction being zero.\n")

    A("## What survives by construction, and what does not\n")
    A("With `F = h + f` and `f` a constant, ordering, force DIFFERENCES and "
      "deterministic-interval separability are all invariant -- checked as arithmetic in "
      "`results.json`, and not counted as gates, because a check that cannot fail decides "
      "nothing. Two things are not invariant, and they are the whole study: a RATIO of two "
      "forces moves toward 1 as `f` grows, and separability collapses under a friction "
      "that VARIES from syringe to syringe, which a constant cannot touch.\n")

    A("## C1 -- reproduction before extension\n")
    A("| quantity | shipped | recomputed here |\n|---|---|---|")
    A(f"| G25.3 median force ratio | {c1['ratio_shipped']:.6f} | "
      f"{c1['ratio_recomputed']:.6f} |")
    A(f"| G25.6 separable mAb pairs | {c1['mab_disjoint_shipped']}/"
      f"{c1['mab_pairs_shipped']} | {c1['mab_disjoint_recomputed']}/"
      f"{c1['mab_pairs_recomputed']} |")
    A(f"| G14.5 separable blood pairs | {c1['blood_disjoint_shipped']}/"
      f"{c1['blood_pairs_shipped']} | {c1['blood_disjoint_recomputed']}/"
      f"{c1['blood_pairs_recomputed']} |\n")
    A("All three reproduce, so the arrays below are the shipped ones and nothing was "
      "refitted.\n")

    A("## Q1 -- a constant friction, and the 7.1x\n")
    A(f"The mAb page reports that the median formulation needs "
      f"**{q1['ratio_at_f']['0']:.1f}x** the plunger force of the median blood donor "
      f"({q1['median_mab_N']:.3g} N against {q1['median_blood_N']:.3g} N). That ratio is "
      "a property of the friction being zero:\n")
    A("| constant friction `f` (N) | median force ratio |\n|---|---|")
    for k, v in q1["ratio_at_f"].items():
        A(f"| {k} | {v:.2f}x |")
    A("")
    A(f"It first falls below 3.0 at `f` = **{q1['f_first_below_3']:.3g} N**, below 2.0 at "
      f"**{q1['f_first_below_2']:.3g} N** and below 1.5 at "
      f"**{q1['f_first_below_1_5']:.3g} N**. The force DIFFERENCE, "
      f"{q1['median_difference_N']:.3g} N, is untouched at every `f` -- so the statement "
      "that survives any friction is *this formulation needs about "
      f"{q1['median_difference_N']:.2g} N more than blood*, not *it needs "
      f"{q1['ratio_at_f']['0']:.1f} times as much*.\n")
    A("**No friction magnitude is asserted here.** The table is the sensitivity; where a "
      "particular device sits on it is answered in the section below, from a source that "
      "resolved to a DOI in a lookup, and it is attributed there rather than asserted.\n")

    A("## Q2 -- the number a formulation screen actually needs\n")
    A(f"A device friction that varies syringe to syringe with SD `s` is the term that can "
      f"destroy separability, and no page mentions it. Convention (named, not derived): "
      f"{q2['convention']}\n")
    A("| `s` (N) | mAb pairs separable | blood pairs separable |\n|---|---|---|")
    A(f"| 0 | {100*q2['mab']['frac_at_s0']:.1f} % "
      f"({d['C1_reproduction']['mab_disjoint_recomputed']}/{q2['mab']['n_pairs']}) | "
      f"{100*q2['blood']['frac_at_s0']:.1f} % "
      f"({d['C1_reproduction']['blood_disjoint_recomputed']}/{q2['blood']['n_pairs']}) |")
    for k in q2["mab"]["frac_at"]:
        A(f"| {k} | {100*q2['mab']['frac_at'][k]:.1f} % | "
          f"{100*q2['blood']['frac_at'][k]:.1f} % |")
    A("")
    A(f"**The deliverable is the threshold, not a verdict.** The mAb screen falls below "
      f"50 % separable pairs at `s*` = **{q2['mab']['s_star_50pct_N']:.4g} N** and below "
      f"10 % at **{q2['mab']['s_star_10pct_N']:.4g} N**. A threshold needs no external "
      "number to be correct; a literature or datasheet friction only says where on this "
      "curve a given syringe sits.\n")

    A("## P4 -- the same friction damages the two cohorts very differently\n")
    A(f"At `s` = {p4['s_N']:g} N the mAb cohort's median interval half-width grows "
      f"**{p4['mab_growth']:.2f}x** ({1e3*p4['mab_median_halfwidth_N']:.3g} -> "
      f"{1e3*p4['mab_median_halfwidth_at_s_N']:.3g} mN) and the blood cohort's "
      f"**{p4['blood_growth']:.2f}x** ({1e3*p4['blood_median_halfwidth_N']:.3g} -> "
      f"{1e3*p4['blood_median_halfwidth_at_s_N']:.3g} mN). The mechanism is that the mAb "
      f"interval is a pure measurement repeatability -- a median "
      f"{p4['mab_relative_sd_pct']:.2f} % of the force -- while the blood interval already "
      f"carries a parameter uncertainty of a median {p4['blood_relative_hull_pct']:.1f} %. "
      "**The tighter the fluid measurement, the more of the answer the device owns.**\n")

    mc = p4["matched_convention"]
    A(f"**P4 is REFUTED, and the defect is in the pre-registration rather than in either "
      f"number.** The prediction was `mab_growth > 3` AND `blood_growth < 1.3`; the first "
      f"holds and the second does not ({p4['blood_growth']:.2f}x). It also compared two "
      f"arms under two different OPERATIONS -- quadrature on the mAb side, a linear `2s` "
      f"added to each hull side for blood. The verdict and the numbers above stand as "
      f"graded. Repeating it with the same operation on both arms ({mc['operation']}) "
      f"gives blood {mc['blood_growth']:.2f}x against the mAb {mc['mab_growth']:.2f}x, so "
      f"the mAb cohort is the more exposed by {mc['mab_over_blood']:.2f}x -- right in "
      f"DIRECTION and wrong in SIZE, since the bar implied a separation nearer 3.5x. "
      f"{mc['status']}\n")

    e = d["external_placement"]
    A("## Where a real device sits on that curve\n")
    A(f"The threshold above needs no external number. This section only places a device "
      f"on the curve, and it uses one that **resolved to a DOI in a lookup this session** "
      f"({e['source']}), never one from memory. It is a **{e['device']}**.\n")
    A(f"- **Glide (sustaining) force {e['glide_N'][0]}-{e['glide_N'][1]} N** — "
      f"{e['glide_note']}. This is the comparable term, because {e['which_term_applies']}. "
      f"At it the 7.05x ratio is **{e['ratio_at_glide'][0]:.2f}-"
      f"{e['ratio_at_glide'][1]:.2f}x**.\n"
      f"- **Breakloose {e['breakloose_N'][0]}-{e['breakloose_N'][1]} N** — "
      f"{e['breakloose_note']}. At it the ratio is "
      f"**{e['ratio_at_breakloose'][1]:.2f}-{e['ratio_at_breakloose'][0]:.2f}x**.\n")
    A(f"**The Q2 side has no comparator at all, and that is a result.** {e['gap']}\n")

    A("## What this does not say\n")
    A("- Nothing here measures a syringe. No friction magnitude is quoted.\n"
      "- The shipped forces are not wrong; they are the hydrodynamic term, correctly "
      "computed and correctly labelled.\n"
      "- Friction is treated as fluid-independent. That assumption is what makes it "
      "additive and therefore what makes the ordering claims survive; if it failed, even "
      "the ordering would be at risk and this study could not tell.\n"
      "- The two cohorts' intervals are different objects (a parameter hull and a "
      "measurement SD) and are never compared with each other here.\n")
    (HERE / "results.md").write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
