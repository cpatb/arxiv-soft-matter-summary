"""analyse_skewness.py — is the effect POLYDISPERSITY, or is it SKEWNESS?

Desmond & Weeks, Phys. Rev. E 90, 022204 (2014), arXiv:1303.4627, show that the SKEWNESS
of a particle size distribution can shift the random-close-packing density MORE than the
polydispersity does, while kurtosis barely matters.

That is a confound in this project's own design, because a lognormal's skewness is fixed
by its CV:

    S_lognormal = (3 + CV^2) * CV        so CV = 0.5 gives S = 1.625

A lognormal ladder therefore moves CV and skewness TOGETHER, and a "CV effect" measured on
it alone is really a joint (CV, skewness) effect along the lognormal locus.

`UniformPSD` breaks the degeneracy: skewness is identically ZERO at every CV. So a
matched-CV comparison of the two families separates the axes:

    same CV, different skewness -> if the two families agree, the effect is CV
                                   if they disagree, skewness is doing the work

This script does that comparison for phi_J and for every structural observable, and
reports the difference in units of the seed scatter.
"""

from __future__ import annotations

import csv
import os
import sys

import numpy as np

from analyse_phase2 import load_all
from analyse_phase3 import OBSERVABLES, flatten
from psd import lognormal_skewness

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

AGREE_IN_SCATTERS = 2.0     # families "agree" if they differ by less than this many sd


def main() -> int:
    recs = [flatten(r) for r in load_all()]
    recs = [r for r in recs if not r.get("failed") and r.get("converged")]
    fams = {r.get("family", "lognormal") for r in recs}
    if "uniform" not in fams:
        print("no uniform-family cells yet -- nothing to compare")
        return 1

    lines = ["=" * 78,
             "SKEWNESS vs POLYDISPERSITY",
             "=" * 78, "",
             "A lognormal's skewness is set by its CV: S = (3 + CV^2)*CV.",
             "A uniform PSD has S = 0 at every CV.",
             "So at MATCHED CV the two families differ only in skewness (and higher",
             "moments), and Desmond & Weeks (PRE 90, 022204) warn that skewness can",
             "outweigh polydispersity for phi_RCP.",
             "",
             f"{'CV':>5} {'S(lognormal)':>13} {'S(uniform)':>11}",
             ]
    for cv in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5):
        lines.append(f"{cv:5.2f} {lognormal_skewness(cv):13.3f} {0.0:11.3f}")
    lines.append("")

    # --- phi_J comparison, read from the phase-2 summary ----------------------
    summ_fp = os.path.join(OUT, "phase2_summary.csv")
    rows = []
    if os.path.exists(summ_fp):
        with open(summ_fp, encoding="utf-8") as fh:
            rows = [r for r in csv.DictReader(fh)]
        for r in rows:
            for k in ("cv_nom", "mu_p", "phi_J", "phi_J_sd", "cv_real"):
                try:
                    r[k] = float(r[k])
                except (TypeError, ValueError):
                    r[k] = float("nan")

    lines.append("phi_J at matched CV and mu_p = 0.30")
    lines.append(f"  {'CV':>5} {'lognormal':>11} {'uniform':>11} {'diff':>9} "
                 f"{'sd':>8} {'n_sd':>6}  verdict")
    out_rows = []
    phiJ_agree: list[bool] = []
    for cv in sorted({r["cv_nom"] for r in rows if r["mu_p"] == 0.3}):
        ln = next((r for r in rows if r["mu_p"] == 0.3 and r["cv_nom"] == cv
                   and r["family"] == "lognormal"), None)
        un = next((r for r in rows if r["mu_p"] == 0.3 and r["cv_nom"] == cv
                   and r["family"] == "uniform"), None)
        if ln is None or un is None:
            continue
        d = un["phi_J"] - ln["phi_J"]
        sd = float(np.nanmean([ln["phi_J_sd"], un["phi_J_sd"]]))
        n_sd = abs(d) / sd if sd > 0 else float("inf")
        agree = n_sd < AGREE_IN_SCATTERS
        phiJ_agree.append(agree)
        lines.append(f"  {cv:5.2f} {ln['phi_J']:11.4f} {un['phi_J']:11.4f} "
                     f"{d:+9.4f} {sd:8.4f} {n_sd:6.2f}  "
                     f"{'agree' if agree else 'DIFFER'}")
        out_rows.append({"cv": cv, "phi_J_lognormal": ln["phi_J"],
                         "phi_J_uniform": un["phi_J"], "diff": d,
                         "sd": sd, "n_sd": n_sd, "agree": agree})
    lines.append("")

    # --- structural observables at matched CV, mu, p_hat ----------------------
    lines.append("structural observables at matched CV, mu_p and p_hat")
    lines.append(f"  {'observable':<34} {'CV':>5} {'p_hat':>8} {'lognorm':>9} "
                 f"{'uniform':>9} {'n_sd':>6}  verdict")
    cells: dict[tuple, list] = {}
    for r in recs:
        k = (r.get("family", "lognormal"), r["cv"], r["mu_p"], r["p_hat_target"])
        cells.setdefault(k, []).append(r)

    struct_disagree = 0
    struct_total = 0
    per_obs: dict[str, list[bool]] = {}
    for (fam, cv, mu, ph), g in sorted(cells.items()):
        if fam != "lognormal" or mu != 0.3:
            continue
        ug = cells.get(("uniform", cv, mu, ph))
        if not ug or cv == 0.0:
            continue
        for key, label in OBSERVABLES:
            a = np.array([x.get(key, np.nan) for x in g], float)
            b = np.array([x.get(key, np.nan) for x in ug], float)
            a, b = a[np.isfinite(a)], b[np.isfinite(b)]
            if a.size < 2 or b.size < 2:
                continue
            sd = float(np.mean([a.std(ddof=1), b.std(ddof=1)]))
            d = float(b.mean() - a.mean())
            # An identical value with zero scatter is agreement, not disagreement.
            # `largest backbone component` is 1.0000 in both families at every CV,
            # which gave 0/0 -> inf -> DIFFER in the first version.
            if d == 0.0:
                n_sd, agree = 0.0, True
            elif sd > 0:
                n_sd = abs(d) / sd
                agree = n_sd < AGREE_IN_SCATTERS
            else:
                n_sd, agree = float("inf"), False
            struct_total += 1
            if not agree:
                struct_disagree += 1
            per_obs.setdefault(key, []).append(agree)
            lines.append(f"  {label:<34} {cv:5.2f} {ph:8.0e} {a.mean():9.4f} "
                         f"{b.mean():9.4f} {n_sd:6.2f}  "
                         f"{'agree' if agree else 'DIFFER'}")
            out_rows.append({"cv": cv, "p_hat": ph, "observable": key,
                             "lognormal": a.mean(), "uniform": b.mean(),
                             "n_sd": n_sd, "agree": agree})

    lines.append("")
    lines.append("=" * 78)
    lines.append("VERDICT")
    lines.append("=" * 78)
    if struct_total:
        lines.append(f"  structural comparisons: {struct_total - struct_disagree}"
                     f"/{struct_total} agree within {AGREE_IN_SCATTERS} sd")
    lines.append("")
    if phiJ_agree:
        lines.append(f"  phi_J comparisons: {sum(phiJ_agree)}/{len(phiJ_agree)} "
                     f"agree within {AGREE_IN_SCATTERS} sd")
    lines.append("")
    lines.append("  PER OBSERVABLE -- this is where the answer actually lives:")
    labels = dict(OBSERVABLES)
    cv_driven, skew_driven = [], []
    for key, agrees in per_obs.items():
        n_ok = sum(agrees)
        frac = n_ok / len(agrees)
        tag = "CV-driven" if frac >= 0.75 else "SKEWNESS-sensitive"
        (cv_driven if frac >= 0.75 else skew_driven).append(labels.get(key, key))
        lines.append(f"    {labels.get(key, key):<38} {n_ok}/{len(agrees)} agree "
                     f"-> {tag}")
    lines.append("")
    lines.append("  => The split is the result, and it is not uniform:")
    if cv_driven:
        lines.append("     CV-DRIVEN (families agree, so skewness is not the cause):")
        for x in cv_driven:
            lines.append(f"       - {x}")
    if skew_driven:
        lines.append("     SKEWNESS-SENSITIVE (families differ at matched CV):")
        for x in skew_driven:
            lines.append(f"       - {x}")
        lines.append("     Desmond & Weeks warned exactly this for packing density,")
        lines.append("     and it is reproduced here for FRICTIONAL packings. Any")
        lines.append("     claim about these observables must name the distribution")
        lines.append("     family, not just the CV.")

    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(OUT, "skewness_verdict.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    if out_rows:
        fp = os.path.join(OUT, "skewness_comparison.csv")
        keys = sorted({k for r in out_rows for k in r})
        with open(fp, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(out_rows)
        print(f"wrote {fp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
