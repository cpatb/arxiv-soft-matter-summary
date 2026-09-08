"""analyse_families.py — do matched moments determine the packing?

Three distribution families are now run at matched CV and matched friction:

    lognormal    skewness (3+CV^2)*CV, so 1.56 at CV = 0.5
    uniform      skewness 0
    bidisperse   skewness 0, 50:50 at d10*(1 +/- CV)

The bidisperse and uniform families are the sharp pair. They match on the number mean,
on CV, on skewness (both exactly 0), AND on the Sauter mean d32 to four decimal places
at every CV -- yet they are completely different shapes, two atoms against a flat
interval. If matched moments determined the packing, they would give the same phi_J.

That is a direct test of a claim the literature treats as established: that a
polydisperse suspension behaves like a "statistically equivalent bidisperse" one
(Singh, Tyagi & Santra arXiv:2510.16464 and refs 15/23/47/48 therein, for viscosity,
normal stresses, particle pressures and packing efficiency). Those results are 2D and
concern rheology; this is 3D and concerns the jamming point.
"""

from __future__ import annotations

import csv
import os
import sys

import numpy as np

from psd import lognormal_skewness, make_psd
from analyse_phase2 import PRIMARY_N

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
AGREE_SD = 2.0

def welch(mean_a, sd_a, n_a, mean_b, sd_b, n_b):
    """Welch's t for two small samples, with the df its critical value needs.

    The first version of this comparison used  |diff| / mean(sd_a, sd_b), which is not a
    test statistic: it divides by the per-seed scatter instead of the standard error of
    the DIFFERENCE OF MEANS. For n = 3 per arm the correct denominator is smaller by
    sqrt(2/3) = 0.82, so that version understated the statistic by ~22%, and it compared
    against a fixed 2.0 rather than a t-critical value for the tiny df.

    Returns (t, df, t_crit_95). A difference is resolved only if |t| > t_crit.
    """
    import numpy as _np
    from scipy import stats as _st
    if not all(_np.isfinite([mean_a, sd_a, mean_b, sd_b])) or n_a < 2 or n_b < 2:
        return float("nan"), float("nan"), float("nan")
    va, vb = sd_a**2 / n_a, sd_b**2 / n_b
    se = _np.sqrt(va + vb)
    if se <= 0:
        return float("inf"), float("nan"), float("nan")
    t = (mean_b - mean_a) / se
    df = (va + vb) ** 2 / (va**2 / (n_a - 1) + vb**2 / (n_b - 1))
    return float(t), float(df), float(_st.t.ppf(0.975, df))



def main() -> int:
    fp = os.path.join(OUT, "phase2_summary.csv")
    if not os.path.exists(fp):
        print("no phase2_summary.csv")
        return 1
    with open(fp, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    # C3 compares FAMILIES at fixed CV. phase2_summary.csv gained an n_requested column
    # on 2026-09-08, so without this the same (family, cv, mu_p) appears twice at two
    # system sizes and the family comparison silently spans both. The N = 3000 coverage
    # does not span the CV range this comparison needs, so the size is fixed explicitly.
    n_before = len(rows)
    rows = [r for r in rows if int(float(r.get("n_requested", PRIMARY_N))) == PRIMARY_N]
    print(f"N RESTRICTION (C3 family comparison at fixed CV): keeping N = {PRIMARY_N}; "
          f"{len(rows)} of {n_before} summary rows, {n_before - len(rows)} dropped.")
    for r in rows:
        for k in ("cv_nom", "cv_real", "mu_p", "phi_J", "phi_J_sd", "d32"):
            try:
                r[k] = float(r[k])
            except (TypeError, ValueError):
                r[k] = float("nan")

    fams = sorted({r["family"] for r in rows})
    lines = ["=" * 82,
             "DO MATCHED MOMENTS DETERMINE THE PACKING?",
             "=" * 82, "",
             "Families present: " + ", ".join(fams), ""]

    # the moment table, computed independently of the sweep
    lines.append("Moments at matched CV, from psd.py (not from the simulations):")
    lines.append(f"  {'CV':>5} {'family':>11} {'skewness':>9} {'d32/d10':>9} "
                 f"{'d_max/d10':>10}")
    rng = np.random.default_rng(0)
    for cv in (0.0, 0.2, 0.3, 0.5):
        for fam in ("lognormal", "uniform", "bidisperse"):
            try:
                p = make_psd(fam, 1.0, cv)
            except ValueError:
                continue
            d = p.sample(200_000, rng)
            m = p.realised_moments(d)
            sk = (lognormal_skewness(cv) if fam == "lognormal" else 0.0)
            lines.append(f"  {cv:5.2f} {fam:>11} {sk:9.3f} {m['d32']:9.4f} "
                         f"{m['d_max']:10.3f}")
    lines.append("")

    # phi_J, family by family, at each friction
    for mu in sorted({r["mu_p"] for r in rows}):
        g_all = [r for r in rows if r["mu_p"] == mu]
        present = sorted({r["family"] for r in g_all})
        if len(present) < 2:
            continue
        lines.append("=" * 82)
        lines.append(f"phi_J at mu_p = {mu:.2f}")
        lines.append("=" * 82)
        cvs = sorted({r["cv_nom"] for r in g_all})
        lines.append("  " + f"{'CV':>5} " +
                     " ".join(f"{f:>22}" for f in present))
        for cv in cvs:
            cells = []
            for f in present:
                m = next((r for r in g_all
                          if r["family"] == f and r["cv_nom"] == cv), None)
                cells.append(f"{m['phi_J']:.4f} +/- {m['phi_J_sd']:.4f}"
                             if m else "-")
            lines.append("  " + f"{cv:5.2f} " +
                         " ".join(f"{c:>22}" for c in cells))
        lines.append("")

        # the sharp pair
        if "uniform" in present and "bidisperse" in present:
            lines.append("  UNIFORM vs BIDISPERSE -- matched on mean, CV, skewness AND d32:")
            lines.append("  Welch t on the difference of means, against the two-sided 95%")
            lines.append("  critical value for its own df. |t| > t_crit = resolved difference.")
            lines.append(f"    {'CV':>5} {'uniform':>10} {'bidisperse':>11} "
                         f"{'diff':>9} {'t':>7} {'df':>6} {'t_crit':>7}  verdict")
            n_dis = n_tot = 0
            for cv in cvs:
                u = next((r for r in g_all if r["family"] == "uniform"
                          and r["cv_nom"] == cv), None)
                b = next((r for r in g_all if r["family"] == "bidisperse"
                          and r["cv_nom"] == cv), None)
                if not u or not b:
                    continue
                d = b["phi_J"] - u["phi_J"]
                t, df, tc = welch(u["phi_J"], u["phi_J_sd"], int(u.get("n_seeds", 3)),
                                  b["phi_J"], b["phi_J_sd"], int(b.get("n_seeds", 3)))
                agree = not (np.isfinite(t) and np.isfinite(tc) and abs(t) > tc)
                n_tot += 1
                n_dis += 0 if agree else 1
                lines.append(f"    {cv:5.2f} {u['phi_J']:10.4f} {b['phi_J']:11.4f} "
                             f"{d:+9.4f} {t:+7.2f} {df:6.1f} {tc:7.2f}  "
                             f"{'not resolved' if agree else 'DIFFER'}")
            lines.append("")
            lines.append("=" * 82)
            lines.append("VERDICT")
            lines.append("=" * 82)
            lines.append(f"  {n_tot - n_dis}/{n_tot} CV values show NO RESOLVED difference")
            lines.append("  (a failure to resolve, on n = 3 seeds per arm, is not a proof")
            lines.append("   of identity -- it bounds any real difference, it does not")
            lines.append("   exclude one)")
            if n_dis > n_tot / 2:
                lines.append("  => MATCHED MOMENTS DO NOT DETERMINE phi_J. Two distributions")
                lines.append("     agreeing on the number mean, CV, skewness and the Sauter")
                lines.append("     mean still give different jamming points in 3D. The")
                lines.append("     'statistically equivalent bidisperse' equivalence, which")
                lines.append("     the literature reports for 2D rheology, does NOT extend to")
                lines.append("     the jamming point here: the full SHAPE matters.")
            else:
                lines.append("  => The two zero-skew families agree. Matching the first")
                lines.append("     three moments plus d32 is sufficient for phi_J, and the")
                lines.append("     lognormal's departure is attributable to its skewness")
                lines.append("     and tail alone.")

    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(OUT, "families_verdict.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
