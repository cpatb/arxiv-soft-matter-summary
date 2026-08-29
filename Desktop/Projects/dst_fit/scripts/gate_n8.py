"""Gate N8: does the phi_m drift survive restricting to the authors' own steady curves?

The strategic question for this package, and the one that decides what the
product claims. The headline is that a single flow curve does not determine
`phi_m`, and part of the evidence is that `phi_m` drifts monotonically across
three published phi-series. If that drift lives in curves the source labs
themselves say are not steady flow curves, the evidence is weaker than written.

The admission rule is theirs, not ours
--------------------------------------
This is the whole point. An exclusion chosen after seeing which curves are
inconvenient would be choosing the answer. Both rules below come from the
sources' own artifacts, and are quoted here so they can be checked:

**Richards et al. 2019.** Their deposited `Figure1/Plotting/Figure1.py` carries

    colours = ['tab:blue','tab:blue','tab:blue','tab:blue','tab:red','tab:red','tab:red']
    labelList = ['0.38','0.40','0.42','0.44','0.46','0.48','0.50']

and shades the background red over the region above the `phi_M` flow curve, i.e.
the shear-jammed regime. So **0.46, 0.48 and 0.50 are marked by the authors as
sitting in the jammed regime** and 0.38 to 0.44 as flowing. The paper's subject is
that this material *oscillates* there, so those three are time averages of an
oscillating rate rather than steady flow curves.

**Guy, Hermes & Poon 2015.** The paper states the flow "was unsteady for
phi >= 0.56 in both (a) and (c), and points represent temporal averages", with
massive wall slip above `phi_m`. So `phi >= 0.56` is excluded, on both series.

What was fixed before this ran
------------------------------
Being explicit, because this gate is a confirmation and not a discovery. The
decision rule below was committed in `PLAN_2026-08-29.md` before the gate existed:

> **Kill criterion:** if `phi_m` consistency largely recovers on the steady
> subsets, the "not identifiable" claim is *weaker than currently written* and
> the spec must say so.
> **What survives regardless:** the degeneracy mechanism itself. `phi` and
> `phi_m` enter only through `phi/phi_J`; `r(phi, phi_m)` runs to +0.99 and
> `phi_m - phi` collapses to 0.0013.

The reviewing session had also already reported that the restriction takes
Richards from reduced chi2 221 to 11.0 and 0/7 to 3/4. So the numbers below are
being *checked*, not found.

  N8a  The admission rule reproduces the sources. Richards' excluded set must be
       exactly the red-coloured curves parsed out of their own plotting script,
       not a list typed here. If that file stops saying what it says, this fails
       rather than silently using a stale rule.
  N8b  The consistency test is re-run on the steady subset with the *same*
       statistic and the same 2-sigma rule the failing gates use, and the outcome
       is reported against the kill criterion above. This condition passes if the
       restriction is *decisive either way* -- if the reduced chi-squared moves by
       more than 2x -- because a rule that changes nothing is not worth carrying.
  N8c  The degeneracy mechanism survives the restriction. On the steady subset
       the worst `r(phi, phi_m)` must still exceed 0.95 and the smallest
       `phi_m - phi` must still be below 0.02. This is the part of the claim that
       is algebra plus a measurement, and it must not depend on which curves are
       admitted. If N8c fails, the headline is in real trouble.

Run:  python scripts/gate_n8.py
Exit 0 only if all three pass.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dst_fit import (  # noqa: E402
    assess,
    consistency_across_curves,
    fit_model,
    thickening_branch,
)
from dst_fit.io_literature import (  # noqa: E402
    GUY_FIG1A_PHI,
    GUY_FIG1B_PHI,
    GUY_PMMA_404,
    GUY_PMMA_3770,
    RICHARDS_CORNSTARCH,
    read_guy_flow_curves,
    read_richards_series,
)
from dst_fit.models import WYART_CATES  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RICHARDS = ROOT / "Competing Timescales Lead to Oscillations in Shear-Thickening Suspensions"
RICHARDS_DATA = RICHARDS / "Figure1" / "Plotting"
RICHARDS_SCRIPT = RICHARDS / "Figure1" / "Plotting" / "Figure1.py"
GUY = ROOT / "Towards a unified description of the rheology of hard-particle suspensions" / "Fig 1" / "Fig 1"

# Guy et al.'s own stated threshold. Quoted in the docstring above.
GUY_UNSTEADY_AT = 0.56

N_SIGMA = 2.0
CHI2_MOVE_MIN = 2.0  # N8b
CORR_MIN = 0.95  # N8c
GAP_MAX = 0.02  # N8c


def richards_regimes() -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Parse the authors' own colour coding. Returns (flowing, jammed) labels."""
    text = RICHARDS_SCRIPT.read_text(encoding="utf-8", errors="replace")
    colours = re.search(r"colours\s*=\s*\[([^\]]+)\]", text)
    labels = re.search(r"labelList\s*=\s*\[([^\]]+)\]", text)
    if not colours or not labels:
        raise ValueError(f"{RICHARDS_SCRIPT.name}: cannot find colours/labelList")
    cols = [c.strip().strip("'\"") for c in colours.group(1).split(",")]
    labs = [x.strip().strip("'\"") for x in labels.group(1).split(",")]
    if len(cols) != len(labs):
        raise ValueError(
            f"{RICHARDS_SCRIPT.name}: {len(cols)} colours for {len(labs)} labels"
        )
    flowing = tuple(lab for lab, c in zip(labs, cols) if c == "tab:blue")
    jammed = tuple(lab for lab, c in zip(labs, cols) if c == "tab:red")
    return flowing, jammed


def fit_all(sweeps, model, labeller):
    """Fit each sweep; return rows of (label, phi, phi_m, err, corr, gap)."""
    rows = []
    for s in sweeps:
        br = thickening_branch(s.viscosity)
        try:
            fit = fit_model(model, s.stress, s.viscosity, branch=br, n_starts=16)
        except Exception:
            continue
        phi, phi_m = float(fit.params[0]), float(fit.params[1])
        corr = float(np.abs(np.asarray(fit.corr)[0, 1]))
        rows.append(
            {
                "label": labeller(s),
                "phi": phi,
                "phi_m": phi_m,
                "err": float(fit.stderr[1]),
                "corr": corr,
                "gap": phi_m - phi,
                "well_specified": assess(fit).well_specified,
                "at_bound": fit.parameters_at_bound(),
            }
        )
    return rows


def consistency(rows):
    vals = np.array([r["phi_m"] for r in rows])
    errs = np.array([r["err"] for r in rows])
    return consistency_across_curves(vals, errs, n_sigma=N_SIGMA)


def report(name, rows, keep):
    """Print the full and restricted consistency, return (full, restricted)."""
    # Two admission rules, and both are properties of the source or of the fit
    # rather than of the answer. The second was added after the first run exposed
    # Guy's 19% curve fitting to phi_m = 0.3000 +/- 2.3768, exactly the lower
    # bound: at phi ~ 0.10 that suspension barely thickens, so the curve carries
    # no information about phi_m and the fitter returned the bound. That is not
    # unsteadiness, so the authors' rule does not catch it; it is a constrained
    # optimum, and counting one as a measurement is what N8 exists to stop.
    sub = [r for r in rows if r["label"] in keep and not r["at_bound"]]
    print(f"\n=== {name}")
    print(f"    {'curve':>10} {'phi':>8} {'phi_m':>8} {'+/-':>8} "
          f"{'r':>7} {'phi_m-phi':>10}  admitted")
    for r in rows:
        if r["at_bound"]:
            mark = f"NO ({', '.join(r['at_bound'])} at a bound)"
        elif r["label"] not in keep:
            mark = "NO (authors: not steady)"
        else:
            mark = "yes"
        print(f"    {r['label']:>10} {r['phi']:8.4f} {r['phi_m']:8.4f} {r['err']:8.4f} "
              f"{r['corr']:7.3f} {r['gap']:10.4f}  {mark}")
    if len(rows) < 2 or len(sub) < 2:
        print("    too few curves to test consistency")
        return None, None, []
    full, restricted = consistency(rows), consistency(sub)
    print(f"    all {full.n_total} curves : {full.n_consistent}/{full.n_total} within "
          f"{N_SIGMA:.0f} sigma, reduced chi2 {full.reduced_chi2:.1f}, "
          f"phi_m range {full.values.min():.4f} to {full.values.max():.4f}")
    print(f"    admitted {restricted.n_total} : {restricted.n_consistent}/"
          f"{restricted.n_total} within {N_SIGMA:.0f} sigma, reduced chi2 "
          f"{restricted.reduced_chi2:.1f}, phi_m range "
          f"{restricted.values.min():.4f} to {restricted.values.max():.4f}")
    return full, restricted, sub


def main() -> int:
    failures: list[str] = []
    results = []

    # ---------------- N8a ----------------
    print("N8a  the admission rule, read from the sources rather than typed here")
    if not RICHARDS_SCRIPT.is_file():
        print(f"     {RICHARDS_SCRIPT} missing; cannot verify the rule")
        return 1
    flowing, jammed = richards_regimes()
    print(f"     Richards Figure1.py: flowing (blue) {flowing}, jammed (red) {jammed}")
    if jammed != ("0.46", "0.48", "0.50"):
        failures.append(
            f"N8a: Richards' red-coded set is {jammed}, not the ('0.46','0.48','0.50') "
            f"this gate was written against. The rule has moved; re-read the source."
        )
    print(f"     Guy et al., stated in the paper: unsteady at phi >= {GUY_UNSTEADY_AT}")
    print(f"     {'PASS' if not failures else 'FAIL'}")

    # ---------------- Richards ----------------
    if RICHARDS_DATA.is_dir():
        model = WYART_CATES.with_fixed(**RICHARDS_CORNSTARCH.as_fixed())
        sweeps = read_richards_series(RICHARDS_DATA)
        rows = fit_all(sweeps, model, lambda s: s.sheet.replace("wt%", ""))
        # Richards' labels are '38'...'50'; the script's are '0.38'...'0.50'.
        keep = {lab.lstrip("0.").rstrip() for lab in flowing}
        keep = {f"{int(round(float(lab) * 100))}" for lab in flowing}
        full, restricted, sub = report("Richards et al. 2019 corn starch", rows, keep)
        if full and restricted:
            results.append(("Richards corn starch", full, restricted, sub))

    # ---------------- Guy ----------------
    for name, path, phi_map, constants in (
        ("Guy PMMA d = 404 nm", GUY / "Fig 1a" / "Fig1a.xls", GUY_FIG1A_PHI, GUY_PMMA_404),
        ("Guy PMMA d = 3770 nm", GUY / "Fig 1b" / "Fig1b.xlsx", GUY_FIG1B_PHI, GUY_PMMA_3770),
    ):
        if not path.is_file():
            continue
        model = WYART_CATES.with_fixed(**constants.as_fixed())
        sweeps = read_guy_flow_curves(path, phi_map)
        rows = fit_all(sweeps, model, lambda s: s.sheet)
        keep = {
            s.sheet for s in sweeps
            if s.phi is not None and s.phi < GUY_UNSTEADY_AT
        }
        full, restricted, sub = report(name, rows, keep)
        if full and restricted:
            results.append((name, full, restricted, sub))

    if not results:
        print("\nno deposits present; nothing to test")
        return 1

    # ---------------- N8b ----------------
    print("\nN8b  is the restriction decisive either way?")
    for name, full, restricted, _sub in results:
        move = full.reduced_chi2 / restricted.reduced_chi2 if restricted.reduced_chi2 else np.inf
        verdict = "PASS" if move >= CHI2_MOVE_MIN or move <= 1 / CHI2_MOVE_MIN else "FAIL"
        print(f"     {name:<24} reduced chi2 {full.reduced_chi2:8.1f} -> "
              f"{restricted.reduced_chi2:7.1f}  ({move:.1f}x)   {verdict}")
        if verdict == "FAIL":
            failures.append(
                f"N8b[{name}]: the restriction moves reduced chi2 by only {move:.2f}x"
            )

    # ---------------- N8c ----------------
    print("\nN8c  does the degeneracy mechanism survive on the steady subset?")
    print("     (this is the half of the claim that must not depend on admission)")
    for name, _full, _restricted, sub in results:
        worst_r = max(r["corr"] for r in sub)
        min_gap = min(r["gap"] for r in sub)
        ok = worst_r > CORR_MIN and min_gap < GAP_MAX
        print(f"     {name:<24} worst r(phi, phi_m) {worst_r:6.3f} (needs > {CORR_MIN}), "
              f"smallest phi_m - phi {min_gap:7.4f} (needs < {GAP_MAX})   "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            failures.append(
                f"N8c[{name}]: on the steady subset the degeneracy is weaker than "
                f"claimed -- worst r = {worst_r:.3f}, smallest gap = {min_gap:.4f}. The "
                f"mechanism behind the headline does not survive the authors' own "
                f"admission rule."
            )

    print("")
    print("--- what this does to the claim ---")
    for name, full, restricted, _sub in results:
        recovered = restricted.n_consistent == restricted.n_total
        print(f"  {name}: on the steady curves {restricted.n_consistent} of "
              f"{restricted.n_total} agree on a single phi_m "
              f"({'CONSISTENT' if recovered else 'still inconsistent'}), reduced chi2 "
              f"{full.reduced_chi2:.1f} -> {restricted.reduced_chi2:.1f}")

    if failures:
        print("\nFAILURES")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("\nall N8 conditions pass")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
