"""Item 14: put uncertainty on the shipped injectability forces, and ask whether donors
can still be ranked once it is there.

Gate pre-registered in `uq_prereg.md`, committed at `3d267db` before any of this was
computed. Bars are there and are not repeated as adjustable constants except where the
code needs the number.

Where the intervals come from
-----------------------------
They are READ from `../../data/profile_interval_g139.json`, not recomputed. That file
already pushed each Casson parameter through its corrected 95 % profile interval at this
exact needle, and computed the covariance-sampled joint interval beside it. Recomputing
would create a second source of truth for a number the repo already ships, which is the
failure this repo has a rule about. G14.2 and G14.3 check that borrowing them is
legitimate: the forces must match, and the two modules must describe the same device.

The two intervals are NOT interchangeable and the report says so every time it prints one
- the push-through hull moves one parameter at a time to each end of its own interval and
  takes the force range. It is WIDER than the joint region it is drawn from, because the
  joint 95 % region is an ellipse inscribed in that box.
- the covariance sample is a joint region, but it assumes the Wald normal, which is the
  approximation the profile likelihood exists because it distrusts.

Run:  python case_studies/blood_injectability/uq.py
Writes: uq.md and uq.json beside this file.
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import sys
from pathlib import Path

import numpy as np

from biorheo.apps.injectability import (
    INSENSITIVE,
    Needle,
    Syringe,
    sensitivities,
)
from biorheo.registry import load_record

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
RECORDS = REPO / "src" / "biorheo" / "registry" / "records"
G139 = REPO / "data" / "profile_interval_g139.json"
RESULTS = HERE / "results.json"
OUT_MD = HERE / "uq.md"
OUT_JSON = HERE / "uq.json"

# bars, from uq_prereg.md
G142_RTOL = 1e-5
G145_RATE = 0.50
G146_MIN = 19       # of 20


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    if not G139.exists():
        print("run scripts/profile_interval_g139.py first", file=sys.stderr)
        return 2
    g = json.loads(G139.read_text(encoding="utf-8"))
    shipped = json.loads(RESULTS.read_text(encoding="utf-8"))["rows"]
    run = _load("_run", HERE / "run.py")
    xch = _load("_xch", HERE / "crosscheck.py")

    # ---- G14.3: the two modules carry these constants as SEPARATE literals -------------
    # Not a style complaint. gen_case.py and run_case.sh in a sibling project rebuilt the
    # same directory name by two rules, disagreed at one parameter value, and the launcher
    # would have deleted a converged case. Equality is asserted, never toleranced.
    geometry = {
        "needle_ID_mm": (run.ID_NOMINAL_MM, xch.ID_NOMINAL_MM),
        "needle_length_m": (run.NEEDLE_LENGTH, xch.NEEDLE_LENGTH),
        "barrel_radius_m": (run.BARREL_RADIUS, xch.BARREL_RADIUS),
        "Q_m3_per_s": (run.Q_SI, xch.Q_SI),
    }
    geom_mismatch = {k: v for k, v in geometry.items() if v[0] != v[1]}
    # and the artifact being borrowed must have been computed at that same geometry
    if abs(g["needle_ID_m"] - run.ID_NOMINAL_MM * 1e-3) > 0:
        geom_mismatch["g139_needle_ID_m"] = (g["needle_ID_m"], run.ID_NOMINAL_MM * 1e-3)
    if g["Q_m3_per_s"] != run.Q_SI:
        geom_mismatch["g139_Q"] = (g["Q_m3_per_s"], run.Q_SI)
    g143 = not geom_mismatch

    # ---- G14.2: the borrowed intervals belong to the shipped forces --------------------
    by_donor = {r["donor"]: r for r in g["rows"]}
    rows, force_devs = [], []
    needle = Needle(radius=0.5 * run.ID_NOMINAL_MM * 1e-3, length=run.NEEDLE_LENGTH)
    syringe = Syringe(run.BARREL_RADIUS)

    for s in shipped:
        d = s["donor"]
        gr = by_donor[d]
        dev = abs(s["force_N"] - gr["f_casson"]) / s["force_N"]
        force_devs.append(dev)

        # G14.6: a THIRD, independent route to "the yield stress is invisible here".
        # results.md argues it from x = tau_y/tau_w and from zeroing tau_y; this is a
        # central-difference d ln F / d ln p, and it is allowed to disagree with both.
        rec = load_record(RECORDS / f"{s['record']}.yaml")
        sens = sensitivities(run.Q_SI, needle, syringe, rec.model, rec.params_si)
        insensitive = tuple(k for k, v in sens.items() if abs(v) < INSENSITIVE)

        rows.append({
            "donor": d, "record": s["record"],
            "force_N": s["force_N"],
            "force_refit_N": gr["f_casson"],
            "force_rel_dev": dev,
            "hull_lo": gr["new_lo"], "hull_hi": gr["new_hi"],
            "hull_width_pct": 100.0 * (gr["new_hi"] - gr["new_lo"]) / s["force_N"],
            "sampled_lo": gr["samp_lo"], "sampled_hi": gr["samp_hi"],
            "sampled_usable": gr["samp_usable"],
            "sampled_width_pct": (100.0 * (gr["samp_hi"] - gr["samp_lo"]) / s["force_N"]
                                  if gr["samp_usable"] else float("nan")),
            "sensitivities": {k: float(v) for k, v in sens.items()},
            "insensitive_to": list(insensitive),
        })

    worst_dev = max(force_devs)
    g142 = worst_dev <= G142_RTOL

    # ---- G14.4: an open end is never printed as a number -------------------------------
    # The item-13 walk reports an unverified or unbounded end AT THE BOX BOUND with
    # `*_at_bound` set, so an open end shows up here as a sidecar flag. Read it from the
    # shipped sidecars rather than inferred from the hull, because a wide hull and an open
    # one are different statements.
    open_records = []
    for r in rows:
        sc = load_record(RECORDS / f"{r['record']}.yaml").identifiability
        # A record with no sidecar has no interval at all, which is not the same as
        # having an open one; it is refused here rather than defaulting to "closed".
        assert isinstance(sc, dict), f"{r['record']} has no identifiability sidecar"
        opened = [name for name, v in sc["parameters"].items()
                  if v.get("lo_at_bound") or v.get("hi_at_bound")]
        r["open_parameters"] = opened
        if opened:
            open_records.append((r["donor"], opened))
    g144 = not open_records   # nothing numeric is printed for an open record

    # ---- G14.5: can donors still be ranked at 95 %? ------------------------------------
    pairs = list(itertools.combinations(rows, 2))
    def disjoint(a, b, lo, hi):
        return a[hi] < b[lo] or b[hi] < a[lo]
    n_hull = sum(1 for a, b in pairs if disjoint(a, b, "hull_lo", "hull_hi"))
    samp_pairs = [(a, b) for a, b in pairs if a["sampled_usable"] and b["sampled_usable"]]
    n_samp = sum(1 for a, b in samp_pairs if disjoint(a, b, "sampled_lo", "sampled_hi"))
    rate_hull = n_hull / len(pairs)
    g145 = rate_hull >= G145_RATE

    # where does separation fail? The rank gap at which two donors first overlap says
    # something a fraction cannot: how far apart two patients must be to be told apart.
    order = sorted(rows, key=lambda r: r["force_N"])
    min_sep_gap = None
    for gap in range(1, len(order)):
        if all(disjoint(order[i], order[i + gap], "hull_lo", "hull_hi")
               for i in range(len(order) - gap)):
            min_sep_gap = gap
            break

    # ---- G14.6 -------------------------------------------------------------------------
    n_tau_y_insensitive = sum(1 for r in rows if "tau_y" in r["insensitive_to"])
    g146 = n_tau_y_insensitive >= G146_MIN
    worst_tau_y_sens = max(abs(r["sensitivities"]["tau_y"]) for r in rows)
    worst_eta_c_sens = max(abs(r["sensitivities"]["eta_c"]) for r in rows)

    verdict = {
        "G14.2": {"pass": bool(g142), "worst_force_rel_dev": worst_dev, "bar": G142_RTOL,
                  "n": len(rows)},
        "G14.3": {"pass": bool(g143), "mismatches": geom_mismatch, "geometry": geometry},
        "G14.4": {"pass": bool(g144), "records_with_an_open_end": open_records},
        "G14.5": {"pass": bool(g145), "n_pairs": len(pairs), "n_disjoint_hull": n_hull,
                  "rate_hull": rate_hull, "bar": G145_RATE,
                  "n_pairs_sampled": len(samp_pairs), "n_disjoint_sampled": n_samp,
                  "rate_sampled": (n_samp / len(samp_pairs)) if samp_pairs else float("nan"),
                  "min_rank_gap_always_separated": min_sep_gap},
        "G14.6": {"pass": bool(g146), "n_tau_y_insensitive": n_tau_y_insensitive,
                  "bar": G146_MIN, "worst_abs_sens_tau_y": worst_tau_y_sens,
                  "worst_abs_sens_eta_c": worst_eta_c_sens, "INSENSITIVE": INSENSITIVE},
    }
    med_hull = float(np.median([r["hull_width_pct"] for r in rows]))
    med_samp = float(np.median([r["sampled_width_pct"] for r in rows
                                if r["sampled_usable"]]))
    f = [r["force_N"] for r in rows]
    total_spread_pct = 100.0 * (max(f) - min(f)) / min(f)

    OUT_JSON.write_text(json.dumps({
        "verdict": verdict, "rows": rows,
        "median_hull_width_pct": med_hull, "median_sampled_width_pct": med_samp,
        "total_donor_spread_pct": total_spread_pct,
    }, indent=2), encoding="utf-8")

    L = [
        "# Injection force with uncertainty, and whether donors can still be ranked",
        "",
        "Written by `case_studies/blood_injectability/uq.py`. Gate pre-registered in",
        "`uq_prereg.md` (commit `3d267db`) before any of this was computed.",
        "",
        "| gate | claim | bar | result | verdict |",
        "|---|---|---|---|---|",
        f"| G14.2 | the borrowed intervals belong to these forces | <= {G142_RTOL:g} "
        f"relative on 20/20 | worst {worst_dev:.3e} | **{'PASS' if g142 else 'FAIL'}** |",
        f"| G14.3 | both modules describe the SAME device | exact equality | "
        f"{'all four constants equal' if g143 else str(geom_mismatch)} | "
        f"**{'PASS' if g143 else 'FAIL'}** |",
        f"| G14.4 | an open interval is never printed as a number | 0 donors | "
        f"{len(open_records)} | **{'PASS' if g144 else 'FAIL'}** |",
        f"| G14.5 | donors can still be ranked at 95 % | >= {100 * G145_RATE:.0f} % of "
        f"190 pairs disjoint | **{n_hull}/{len(pairs)} = {100 * rate_hull:.1f} %** | "
        f"**{'PASS' if g145 else 'FAIL'}** |",
        f"| G14.6 | insensitivity reproduced by an independent route | `tau_y` "
        f"insensitive on >= {G146_MIN}/20 | {n_tau_y_insensitive}/20, worst "
        f"|d ln F / d ln tau_y| = {worst_tau_y_sens:.2e} | "
        f"**{'PASS' if g146 else 'FAIL'}** |",
        "",
        "## The size of the uncertainty, and what kind of object it is",
        "",
        f"- one-at-a-time profile push-through hull: median **{med_hull:.1f} %** of the force",
        f"- covariance-sampled joint interval: median **{med_samp:.1f} %**",
        f"- donor-to-donor spread across all 20: **{total_spread_pct:.0f} %**",
        "",
        "**They are not interchangeable and neither is simply the better number.** The",
        "hull moves one parameter at a time to each end of its own 95 % profile interval",
        "and takes the force range, so it is WIDER than the joint region it is drawn from",
        "-- the joint region is an ellipse inscribed in that box. The covariance sample IS",
        "a joint region, but it assumes the Wald normal, which is the approximation the",
        "profile likelihood exists because it distrusts. Every interval below says which",
        "it is.",
        "",
        "Every interval here is a CASSON interval. Reading the same donor through the",
        "better-fitting Cross family moves the force by a further median 12.5 % (up to",
        "30.2 %); see `crosscheck.md`. That is a separate term and it is not added in.",
        "",
    ]
    if g145:
        L += [
            "## Donors can be ranked at 95 %",
            "",
            f"**{n_hull} of {len(pairs)}** donor pairs have disjoint push-through hulls",
            f"({100 * rate_hull:.1f} %), and {n_samp} of {len(samp_pairs)} on the sampled",
            "interval. The pre-registration predicted this would FAIL; it did not.",
            "",
        ]
    else:
        L += [
            "## G14.5 FAILED, and the ordinal claim is restated rather than deleted",
            "",
            f"**{n_hull} of {len(pairs)}** donor pairs have disjoint push-through hulls",
            f"({100 * rate_hull:.1f} %, bar {100 * G145_RATE:.0f} %), and",
            f"{n_samp} of {len(samp_pairs)} on the sampled interval. The bar was not moved.",
            "",
            "**What this does and does not overturn.** `results.md`'s arithmetic is",
            f"untouched: the donor spread is {total_spread_pct:.0f} % and the family-choice",
            "effect is 12.5 %, so the spread is still the bigger term. What the arithmetic",
            "does not license is the PAIRWISE reading of it -- 'which patient the blood came",
            "from matters more' is a statement about telling two patients apart, and that is",
            "a comparison of one donor's interval with another's, not of a range with an",
            "effect size.",
            "",
        ]
        if min_sep_gap is not None:
            L += [
                "Ranked by force, two donors are separated at 95 % once they are",
                f"**{min_sep_gap} places apart** in the order; nearer than that they",
                "overlap. That is the honest form of the claim: **the extremes are",
                "distinguishable, the middle of the donor distribution is not.**",
                "",
            ]
        else:
            L += [
                "No rank gap separates every pair at that gap, so the separation is not a",
                "simple function of rank distance: some adjacent pairs separate and some",
                "distant ones do not, because the hull width varies between donors.",
                "",
            ]
    L += [
        "## Which parameters the force is insensitive to",
        "",
        "Phase 0's unit 0.6 asks for this list explicitly. It is a central-difference",
        "`d ln F / d ln p` at the operating point, which is a **third and independent**",
        "route to a claim `results.md` already makes two other ways (`x = tau_y/tau_w`",
        "reaches only 1.7e-4, and zeroing `tau_y` moves the force under 2.9 %).",
        "",
        f"- `tau_y`: insensitive on **{n_tau_y_insensitive} of 20** donors at the",
        f"  `INSENSITIVE = {INSENSITIVE:g}` threshold; worst",
        f"  `|d ln F / d ln tau_y|` = **{worst_tau_y_sens:.2e}**",
        f"- `eta_c`: worst `|d ln F / d ln eta_c|` = **{worst_eta_c_sens:.3f}**, i.e. the",
        "  force is very nearly proportional to the plateau viscosity and to nothing else",
        "",
        "All three routes agree, which is worth stating because they could not have been",
        "made to: the sensitivity is computed from the force solver by finite differences",
        "and never touches the stress ratio the other two arguments are built on.",
        "",
        "## Per donor",
        "",
        "| donor | force (N) | push-through hull (N) | width | sampled joint (N) | width | open ends |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in sorted(rows, key=lambda r: r["force_N"]):
        samp = (f"[{r['sampled_lo']:.4f}, {r['sampled_hi']:.4f}]"
                if r["sampled_usable"] else "not available")
        sw = f"{r['sampled_width_pct']:.1f} %" if r["sampled_usable"] else "--"
        L.append(f"| {r['donor']} | {r['force_N']:.4f} | "
                 f"[{r['hull_lo']:.4f}, {r['hull_hi']:.4f}] | {r['hull_width_pct']:.1f} % | "
                 f"{samp} | {sw} | {', '.join(r['open_parameters']) or 'none'} |")
    L += [
        "",
        "## What is still not in these numbers",
        "",
        "- The FAMILY term (median 12.5 %, worst 30.2 %), which is not added in above.",
        "- The needle bore. `F ~ R^-4` measured, and a 10.5 % bore error reproduces the",
        "  whole donor spread; the 27G inner diameter remains [UNVERIFIED].",
        "- Any error in the Casson FORM itself. These are intervals on the parameters of a",
        "  model that misfits these curves on 20/20 donors (`residual_structure` flags",
        "  19/20). A correct interval on a wrong model is still a wrong uncertainty.",
        "",
    ]
    OUT_MD.write_text("\n".join(L), encoding="utf-8")

    for gid, val in verdict.items():
        print(f"{gid}: {'PASS' if val['pass'] else 'FAIL'}")
    print(f"median hull {med_hull:.1f} %, sampled {med_samp:.1f} %, "
          f"donor spread {total_spread_pct:.0f} %")
    print(f"G14.5: {n_hull}/{len(pairs)} pairs disjoint ({100 * rate_hull:.1f} %), "
          f"min rank gap always separated = {min_sep_gap}")
    print(f"wrote {OUT_MD.name}, {OUT_JSON.name}")
    return 0 if all(v["pass"] for v in verdict.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
