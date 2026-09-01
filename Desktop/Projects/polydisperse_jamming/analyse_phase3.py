"""analyse_phase3.py — isolate mechanism M3: structure at matched mechanical state.

WHY THIS NEEDS NO NEW SIMULATIONS
---------------------------------
SPEC §5 Phase 3 asks for structure at matched `phi/phi_J` AND matched `d32`. The second
condition turns out to be free, and the reason is worth stating because it is also the
cleanest statement of what separates M1 from M3:

    A DEM packing of frictional hard spheres is SCALE-INVARIANT. Multiplying every
    diameter by a constant is a change of units and changes no packing observable --
    not phi_J, not Z, not the rattler fraction, not the fabric, not the force
    distribution. So at fixed distribution SHAPE (fixed CV), matching on d10 and
    matching on d32 are the SAME experiment.

M1 and M3 are therefore separated by construction rather than by matching: M1 lives in
the interstitial fluid, where an absolute length enters through the permeability
`kappa ~ d32^2`, and M3 lives in the packing, which cannot see an absolute length at
all. Everything the DEM reports is, by definition, free of M1.

That leaves ONE matching condition, the mechanical one. Cells at equal reduced pressure
`p_hat = p/E` are at equal distance from jamming in the only sense the contact law
knows, so the Phase-2 sweep already contains the matched-condition comparison: fix
`p_hat` and `mu_p`, vary CV, read the structure.

THE PRE-REGISTERED TEST (SPEC §5, unchanged)
--------------------------------------------
    A residual CV effect after matching is M3. If the residual is smaller than the
    seed-to-seed scatter, M3 is REFUTED and the honest answer is that the observation
    is not structural.

Graded here per observable, with the effect expressed in units of seed scatter. An
observable is called MOVED only if |effect| > 2 x pooled seed sd across the CV range.
"""

from __future__ import annotations

import csv
import os
import sys

import numpy as np

from analyse_phase2 import load_all

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

# Observables tested for a residual CV dependence. Each is (key, human label).
OBSERVABLES = [
    ("z_all", "Z over all grains"),
    ("rattler_frac_k2", "rattler fraction (k=2, true rattlers)"),
    ("rattler_frac_k3", "rattler fraction (k=3)"),
    ("z_nr_k2", "Z over backbone (k=2)"),
    ("z_nr_k3", "Z over backbone (k=3)"),
    ("largest_frac_k2", "largest backbone component (k=2)"),
    ("phi", "packing fraction at fixed p_hat"),
    ("force_cv", "contact-force CV"),
    ("force_tail_slope", "contact-force tail slope"),
]

MOVED_IN_SCATTERS = 2.0     # |effect| must exceed this many pooled seed sd


def flatten(r: dict) -> dict:
    """Lift the nested rattler/force dicts into flat keys."""
    out = dict(r)
    for k, v in (r.get("rattler") or {}).items():
        out[k] = v
    f = r.get("force") or {}
    out["force_cv"] = f.get("cv", float("nan"))
    out["force_tail_slope"] = f.get("tail_slope", float("nan"))
    out["force_n"] = f.get("n", 0)
    return out


def main() -> int:
    recs = [flatten(r) for r in load_all()]
    n0 = len(recs)
    recs = [r for r in recs if not r.get("failed") and r.get("converged")]
    print(f"N: {n0} -> {len(recs)} converged, non-failed")
    if not recs:
        return 1

    # group by the matched condition: (mu_p, p_hat_target)
    cells: dict[tuple, list[dict]] = {}
    for r in recs:
        cells.setdefault((r.get("family", "lognormal"), r["mu_p"],
                          r["p_hat_target"]), []).append(r)

    lines = ["=" * 78,
             "PHASE 3 -- structure at matched mechanical state (M3)",
             "=" * 78, "",
             "Matching: equal p_hat = p/E and equal mu_p. Matching on d32 is automatic,",
             "because a DEM packing is scale-invariant -- see the module docstring.",
             f"An observable is MOVED only if |effect| > {MOVED_IN_SCATTERS} x pooled "
             "seed sd.", ""]

    rows = []
    verdicts: dict[str, list[str]] = {}

    for (fam, mu, p_hat), g in sorted(cells.items()):
        by_cv: dict[float, list[dict]] = {}
        for r in g:
            by_cv.setdefault(r["cv"], []).append(r)
        cvs = sorted(by_cv)
        if len(cvs) < 3:
            continue
        n_seeds = {c: len(by_cv[c]) for c in cvs}
        if max(n_seeds.values()) < 2:
            continue

        lines.append(f"--- {fam}, mu_p = {mu:.2f}, p_hat = {p_hat:.1e} "
                     f"(CV values {cvs}, seeds {[n_seeds[c] for c in cvs]}) ---")
        lines.append(f"  {'observable':<38} {'CV=lo':>9} {'CV=hi':>9} "
                     f"{'effect':>9} {'seed sd':>9} {'n_sd':>6}  verdict")

        for key, label in OBSERVABLES:
            means, sds = [], []
            for c in cvs:
                v = np.array([x.get(key, np.nan) for x in by_cv[c]], float)
                v = v[np.isfinite(v)]
                if v.size == 0:
                    means.append(np.nan)
                    sds.append(np.nan)
                    continue
                means.append(float(v.mean()))
                sds.append(float(v.std(ddof=1)) if v.size > 1 else np.nan)
            means = np.array(means)
            sds = np.array(sds)
            if not np.isfinite(means).all():
                continue
            pooled = float(np.nanmean(sds)) if np.isfinite(sds).any() else np.nan
            effect = float(means[-1] - means[0])
            # An observable that does not change is NOT 'moved', however small its
            # scatter. `largest backbone component` is identically 1.0000 at every CV
            # with sd 0, which gave effect/sd = 0/0 -> inf -> 'MOVED' in the first
            # version. A constant is a finding, but it is the opposite finding.
            if effect == 0.0:
                n_sd = 0.0
                moved = False
            elif pooled and np.isfinite(pooled) and pooled > 0:
                n_sd = abs(effect) / pooled
                moved = n_sd > MOVED_IN_SCATTERS
            else:
                # nonzero effect, zero scatter: real but the sd is unestimable
                n_sd = float('inf')
                moved = True
            verdict = "MOVED" if moved else (
                "CONSTANT" if effect == 0.0 else "not resolved")
            # monotonicity is a stronger claim than the endpoints
            mono = (np.all(np.diff(means) >= -1e-12)
                    or np.all(np.diff(means) <= 1e-12))
            if moved and mono:
                verdict += ", monotone"
            lines.append(f"  {label:<38} {means[0]:9.4f} {means[-1]:9.4f} "
                         f"{effect:+9.4f} {pooled:9.4f} {n_sd:6.2f}  {verdict}")
            verdicts.setdefault(key, []).append(verdict)
            rows.append({"family": fam, "mu_p": mu, "p_hat": p_hat,
                         "observable": key,
                         "cv_lo": cvs[0], "cv_hi": cvs[-1],
                         "mean_lo": means[0], "mean_hi": means[-1],
                         "effect": effect, "seed_sd": pooled, "n_sd": n_sd,
                         "monotone": bool(mono), "verdict": verdict,
                         "n_seeds_min": min(n_seeds.values())})
        lines.append("")

    if not rows:
        lines.append("No matched cell had >= 3 CV values with >= 2 seeds yet.")
    else:
        fp = os.path.join(OUT, "phase3_structure.csv")
        with open(fp, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=sorted(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {fp}")

        lines.append("=" * 78)
        lines.append("M3 VERDICT, pooled over matched cells")
        lines.append("=" * 78)
        for key, label in OBSERVABLES:
            vs = verdicts.get(key, [])
            if not vs:
                continue
            n_moved = sum(1 for v in vs if v.startswith("MOVED"))
            lines.append(f"  {label:<38} MOVED in {n_moved}/{len(vs)} matched cells")
        any_moved = any(v.startswith("MOVED") for vs in verdicts.values() for v in vs)
        lines.append("")
        if any_moved:
            lines.append("=> A residual CV dependence SURVIVES matching on the")
            lines.append("   mechanical state. M3 is not refuted; see which observables")
            lines.append("   moved and whether they moved monotonically.")
        else:
            lines.append("=> NO observable moved by more than "
                         f"{MOVED_IN_SCATTERS}x seed scatter.")
            lines.append("   M3 is REFUTED on this evidence: at matched p_hat and mu_p,")
            lines.append("   polydispersity does not restructure the contact network.")
            lines.append("   That is a publishable negative and must be reported as one.")

    txt = "\n".join(lines)
    print(txt)
    fp = os.path.join(OUT, "phase3_verdict.txt")
    with open(fp, "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    print(f"wrote {fp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
