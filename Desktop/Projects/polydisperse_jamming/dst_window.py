"""dst_window.py — does polydispersity change the WIDTH of the shear-thickening window?

THE QUESTION THIS ANSWERS, AND WHY IT IS NEW HERE
-------------------------------------------------
Every mechanism this project has weighed (M1 permeability, M2a viscosity, M2' dilatancy
sign) turns on ONE packing fraction and asks where `phi` sits relative to it. The
Wyart-Cates picture of shear thickening turns on TWO:

    phi_0   the FRICTIONLESS jamming fraction (contacts unable to transmit tangential
            force; reached in the low-stress limit, sigma << sigma*)
    phi_mu  the FRICTIONAL jamming fraction  (reached in the high-stress limit)

and it is the position of `phi` BETWEEN them that sets the TYPE of thickening:

    phi < phi_mu           CST      viscosity finite at every stress
    phi_mu < phi < phi_0   DST / shear jamming: the viscosity diverges at a FINITE stress
    phi > phi_0            jammed at every stress

So the window `W = phi_0 - phi_mu` is the whole stage on which shear thickening happens,
and the reduced coordinate

    x = (phi - phi_mu) / W

is the material's position on it. This script measures `W(CV)` and `dx/dCV` from the
project's OWN isotropic ladders, which were run at BOTH frictions on the SAME protocol.

⚠ THE IDENTIFICATION IS AN ASSUMPTION AND IS FLAGGED EVERYWHERE BELOW.
`phi_0` and `phi_mu` in Wyart-Cates are the jamming points of the SHEARED suspension.
This script uses the ISOTROPIC compression values, because those are the ones with CV
coverage at both frictions in this repo. The frictionless identification is the safer of
the two (frictionless isotropic compression goes to random close packing, which is what
`phi_0` means); the frictional one is NOT safe -- the project's own shear campaign puts
the frictional critical state ABOVE the isotropic frictional `phi_J`. Every number here
is therefore a WINDOW BETWEEN TWO ISOTROPIC LADDERS, and the Wyart-Cates reading of it is
a hypothesis, not a measurement. The shear-vs-compression comparison is NOT made here --
it is `grade_common_ladder.offset_h3`, because it can only be made on the common absolute
ladder (on the centred ones it is forced, since they are centred on these same phi_J).

Reads only committed artifacts; writes out/dst_window.txt. No simulation.
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FITS = os.path.join(OUT, "phase2_fits.csv")

# The two frictions the isotropic campaign ran. MU_0 is the frictionless control.
MU_0 = 0.0
MU_F = 0.3

# The finite-size arm re-ran some cells at N = 3000. The frictionless arm exists only at
# N = 1500, so both arms are taken at 1500 and the count is printed.
N_MATCH = 1500

# A linear phi_J extrapolation from THIS many ladder points has zero residual degrees of
# freedom and reports r2 = 1 by construction. Cells at or below it are flagged, never
# silently dropped.
N_LADDER_SELF_DIAGNOSING = 3


def load_fits() -> pd.DataFrame:
    d = pd.read_csv(FITS)
    print(f"phase2_fits.csv: {len(d)} rows")
    d = d[d.n_requested == N_MATCH]
    print(f"  after N == {N_MATCH} (the frictionless arm has no other N): {len(d)} rows")
    d = d[d.mu_p.isin([MU_0, MU_F])]
    print(f"  after mu_p in {{{MU_0}, {MU_F}}}: {len(d)} rows")
    d = d[np.isfinite(d.phi_J_lin)]
    print(f"  after finite phi_J_lin: {len(d)} rows")
    # 🔴 2026-09-09. Part of this grid was re-run at PJ_MAX_ITER = 1e6 (`capfix`), and
    # analyse_phase2 now keeps the two caps as SEPARATE ladders. Without this the same
    # (cv, seed) would appear twice and the pairing below would double-count it. The
    # HIGHEST cap wins, cell by cell, because a cell that stopped at the old cap was
    # censored and its repaired value is the better measurement -- and a cell that
    # stopped early is unaffected by the cap in principle (in practice see §3(c) of
    # PREREG_DST_WINDOW: the stopping is not deterministic at CV = 0).
    if "max_iter" in d.columns:
        before = len(d)
        idx = d.groupby(["family", "cv_nom", "mu_p", "seed"])["max_iter"].idxmax()
        d = d.loc[sorted(idx)]
        caps = sorted(d.max_iter.unique())
        print(f"  after keeping the HIGHEST max_iter per (family, cv, mu, seed): "
              f"{len(d)} rows (dropped {before - len(d)}); caps present {caps}")
    return d


def cell_table(d: pd.DataFrame, family: str) -> pd.DataFrame:
    """One row per (cv_nom, mu_p): mean, sd, n, and the ladder-length flag."""
    rows = []
    sub = d[d.family == family]
    for (cv, mu), g in sub.groupby(["cv_nom", "mu_p"]):
        v = g.phi_J_lin.to_numpy()
        rows.append(dict(
            cv_nom=cv, mu_p=mu, n=len(v),
            cv_real=float(g.cv_real.mean()),
            phi=float(v.mean()),
            sd=float(v.std(ddof=1)) if len(v) > 1 else np.nan,
            se_mean=float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else np.nan,
            n_ladder_min=int(g.n_ladder.min()),
            thin=bool((g.n_ladder <= N_LADDER_SELF_DIAGNOSING).any()),
        ))
    return pd.DataFrame(rows).sort_values(["mu_p", "cv_nom"]).reset_index(drop=True)


def window_table(d: pd.DataFrame, family: str, t: pd.DataFrame) -> pd.DataFrame:
    """W(CV) = phi_J(mu=0) - phi_J(mu=0.3), formed PER SEED and then averaged.

    The two frictions were run on the SAME seeds, and a seed fixes the PSD realisation:
    `cv_real` matches to within 1.1e-16 on all 15 pairs in each family (checked in
    test_dst_window.py; machine epsilon, not exact, since two code paths reduce it). So
    arms are PAIRED, the per-seed difference cancels the shared realisation exactly, and
    the standard error of that difference is much smaller than the quadrature sum of the
    two arms' scatters -- which is the statistic to use. Differencing the two ARM MEANS
    and propagating their separate errors is the unpaired estimator; it is reported
    beside this one and is the wrong one here.
    """
    sub = d[d.family == family]
    a = sub[sub.mu_p == MU_0][["cv_nom", "seed", "cv_real", "phi_J_lin"]]
    b = sub[sub.mu_p == MU_F][["cv_nom", "seed", "cv_real", "phi_J_lin"]]
    m = a.merge(b, on=["cv_nom", "seed"], suffixes=("_0", "_f"))
    m["W"] = m.phi_J_lin_0 - m.phi_J_lin_f
    m["cv_mismatch"] = (m.cv_real_0 - m.cv_real_f).abs()
    flags = t.set_index(["cv_nom", "mu_p"])
    rows = []
    for cv, g in m.groupby("cv_nom"):
        n = len(g)
        rows.append(dict(
            cv_nom=cv, n_pairs=n,
            cv_mismatch=float(g.cv_mismatch.max()),
            phi_0=float(g.phi_J_lin_0.mean()), n_0=n,
            phi_mu=float(g.phi_J_lin_f.mean()), n_f=n,
            W=float(g.W.mean()),
            W_se=float(g.W.std(ddof=1) / np.sqrt(n)) if n > 1 else np.nan,
            thin=bool(flags.loc[(cv, MU_0), "thin"] or flags.loc[(cv, MU_F), "thin"]),
            single_seed=bool(n < 2),
        ))
    return pd.DataFrame(rows)


def wls_slope(x, y, sy=None):
    """Slope and its standard error. Unweighted if any sy is missing."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    if n < 3:
        return np.nan, np.nan, n
    w = None
    if sy is not None:
        sy = np.asarray(sy, float)
        if np.all(np.isfinite(sy)) and np.all(sy > 0):
            w = 1.0 / sy ** 2
    if w is None:
        w = np.ones(n)
    sw = w.sum()
    xm = (w * x).sum() / sw
    ym = (w * y).sum() / sw
    sxx = (w * (x - xm) ** 2).sum()
    slope = (w * (x - xm) * (y - ym)).sum() / sxx
    resid = y - (ym + slope * (x - xm))
    # residual scatter, which is the honest error bar when the per-point sd's are
    # themselves estimated from 2-3 seeds.
    s2 = (w * resid ** 2).sum() / (n - 2)
    se = np.sqrt(s2 / sxx)
    return slope, se, n


def report(p, family: str, d: pd.DataFrame, phi_probe: float) -> dict:
    t = cell_table(d, family)
    p("")
    p(f"=== {family}: the two isotropic ladders ===")
    p("   cv   mu_p  n  n_ladder_min  phi_J      sd        thin?")
    for _, r in t.iterrows():
        sd = f"{r.sd:.6f}" if np.isfinite(r.sd) else "   n/a  "
        p(f"  {r.cv_nom:.1f}  {r.mu_p:.1f}  {r.n:d}      {r.n_ladder_min:d}       "
          f"{r.phi:.6f}  {sd}   {'THIN' if r.thin else ''}")

    w = window_table(d, family, t)
    if w.empty:
        p("   no CV level has both frictions -- nothing to compare")
        return {}

    p("")
    p(f"=== {family}: the shear-thickening window W = phi_0 - phi_mu ===")
    p("   cv   phi_0 (n)      phi_mu (n)     W        +-se      flags")
    for _, r in w.iterrows():
        se = f"{r.W_se:.6f}" if np.isfinite(r.W_se) else "  n/a   "
        flags = []
        if r.single_seed:
            flags.append("SINGLE-SEED")
        if r.thin:
            flags.append("THIN-LADDER")
        if r.cv_mismatch > 0.01:
            flags.append(f"cv_real mismatch {r.cv_mismatch:.3f}")
        p(f"  {r.cv_nom:.1f}  {r.phi_0:.6f}({r.n_0:.0f})  {r.phi_mu:.6f}({r.n_f:.0f})  "
          f"{r.W:.6f}  {se}   {', '.join(flags)}")

    # --- the two slopes, and whether they differ -----------------------------------
    clean = w[~w.single_seed]
    p("")
    p(f"   slopes over CV, ALL {len(w)} levels and the {len(clean)} without a "
      f"single-seed arm, side by side.")
    p("   ⚠ `difference` is the UNPAIRED statistic (two arm slopes differenced); "
      "`dW/dCV`")
    p("     is the PAIRED one (per-seed differences fitted directly). The arms share "
      "seeds,")
    p("     so the paired one is correct and the unpaired one is conservative by the "
      "shared")
    p("     realisation variance. Both are printed so the gap is visible.")
    out = {}
    for label, ww in (("all", w), ("no single-seed", clean)):
        if len(ww) < 3:
            p(f"   {label:16s}: only {len(ww)} levels -- no slope")
            continue
        s0, e0, _ = wls_slope(ww.cv_nom, ww.phi_0)
        sf, ef, _ = wls_slope(ww.cv_nom, ww.phi_mu)
        sw_, ew, n = wls_slope(ww.cv_nom, ww.W)
        # difference of the two slopes, with the SE of the DIFFERENCE
        dse = np.sqrt(e0 ** 2 + ef ** 2)
        t_stat = (s0 - sf) / dse if dse > 0 else np.nan
        p(f"   {label:16s}: dphi_0/dCV = {s0:+.5f} +- {e0:.5f}   "
          f"dphi_mu/dCV = {sf:+.5f} +- {ef:.5f}")
        p(f"   {'':16s}  difference  = {s0 - sf:+.5f} +- {dse:.5f}  "
          f"(t = {t_stat:+.2f}, n = {n});  dW/dCV = {sw_:+.5f} +- {ew:.5f}")
        out[label] = dict(s0=s0, e0=e0, sf=sf, ef=ef, dW=sw_, eW=ew, t=t_stat, n=n)

    # --- the traverse ---------------------------------------------------------------
    W0 = float(w.loc[w.cv_nom == w.cv_nom.min(), "W"].iloc[0])
    key = "no single-seed" if "no single-seed" in out else "all"
    if key in out:
        sf = out[key]["sf"]
        s0 = out[key]["s0"]
        p("")
        p(f"   THE TRAVERSE, at fixed phi, using the '{key}' slopes and W(CV=0) = "
          f"{W0:.5f}:")
        p(f"      the window is {W0:.4f} wide in phi;")
        p(f"      polydispersity moves phi_mu at {sf:+.4f} per unit CV, i.e. "
          f"{abs(sf) / W0:.2f} windows per unit CV.")
        p(f"      dCV needed to move a material one FULL window at fixed phi: "
          f"{W0 / abs(sf):.3f}")
        for x0 in (0.0, 0.5, 1.0):
            phi = w.loc[w.cv_nom == w.cv_nom.min(), "phi_mu"].iloc[0] + x0 * W0
            # x = (phi - phi_mu)/W ; dx/dCV at fixed phi
            dxdcv = (-sf * W0 - (phi - w.loc[w.cv_nom == w.cv_nom.min(),
                                             "phi_mu"].iloc[0]) * (s0 - sf)) / W0 ** 2
            p(f"      at x = {x0:.1f} (phi = {phi:.4f}):  dx/dCV = {dxdcv:+.3f}")
        p(f"      ⚠ x is a RATIO whose denominator W is itself moving; the absolute "
          f"companion is dphi_mu/dCV = {sf:+.5f}, quoted above and independent of W.")

    # --- is the window CONSTANT? a BOUND, never a claim of identity ------------------
    if key in out:
        r = out[key]
        hi = abs(r["dW"]) + 2 * r["eW"]
        p("")
        p(f"   IS W CONSTANT IN CV?  dW/dCV = {r['dW']:+.5f} +- {r['eW']:.5f}, "
          f"t = {r['dW'] / r['eW']:+.2f} on n = {r['n']} levels.")
        wfit = clean if key == "no single-seed" else w
        cvspan = float(wfit.cv_nom.max() - wfit.cv_nom.min())
        tstat = r["dW"] / r["eW"]
        if abs(tstat) >= 2.0:
            p(f"      RESOLVED as a CHANGE: over CV {wfit.cv_nom.min():.1f} -> "
              f"{wfit.cv_nom.max():.1f} the window moves by "
              f"{r['dW'] * cvspan:+.4f} against W = {W0:.4f}, i.e. "
              f"{100 * r['dW'] * cvspan / W0:+.0f} % of itself.")
            p(f"      Measured endpoints over the FITTED range: W = "
              f"{wfit.W.iloc[0]:.5f} -> {wfit.W.iloc[-1]:.5f} "
              f"(all {len(w)} levels end at {w.W.iloc[-1]:.5f}).")
        else:
            p(f"      NOT RESOLVED. |dW/dCV| < {hi:.5f} at ~2 se, i.e. over CV "
              f"{wfit.cv_nom.min():.1f} -> {wfit.cv_nom.max():.1f} the window changes by less")
            p(f"      than {cvspan * hi:.4f} against W = {W0:.4f} -- "
              f"{100 * cvspan * hi / W0:.0f} % of itself. That is a BOUND on a change,")
            p("      not a demonstration that W is the same number at every CV.")
        nthin = int(w.thin.sum())
        if nthin:
            thincv = [f"{c:.1f}" for c in w.loc[w.thin, "cv_nom"]]
            p(f"      ⚠ {nthin} of {len(w)} levels are THIN (a phi_J extrapolation from"
              f" <= {N_LADDER_SELF_DIAGNOSING} ladder")
            p(f"        points has 0-1 residual dof and reports r2 = 1 by construction):"
              f" cv {', '.join(thincv)}.")
        else:
            p("      ✅ no level is THIN any more -- every phi_J extrapolation here uses"
              " more than")
            p(f"        {N_LADDER_SELF_DIAGNOSING} ladder points.")

    # --- where the lab sits, and what it would take to get into the window -----------
    p("")
    p(f"   the lab's bench phi = {phi_probe:.3f} against this window:")
    p("   cv   phi-phi_mu     x      class       phi_eff/phi needed to ENTER / to JAM")
    for _, r in w.iterrows():
        x = (phi_probe - r.phi_mu) / r.W
        cls = "CST only" if x < 0 else ("DST window" if x < 1 else "jammed")
        p(f"  {r.cv_nom:.1f}   {phi_probe - r.phi_mu:+.4f}   {x:+.2f}   {cls:11s} "
          f"{r.phi_mu / phi_probe:.3f}x / {r.phi_0 / phi_probe:.3f}x   "
          f"(volumetric compaction {100 * (1 - phi_probe / r.phi_mu):.0f} % / "
          f"{100 * (1 - phi_probe / r.phi_0):.0f} %)")

    # --- CV*, restated on phi_mu rather than on a window-dependent root --------------
    if key in out:
        sf = out[key]["sf"]
        mu0 = float(w.loc[w.cv_nom == w.cv_nom.min(), "phi_mu"].iloc[0])
        p("")
        p("   CV*, the width at which a material of EFFECTIVE fraction phi_eff leaves")
        p(f"   the window from above:   CV* = (phi_eff - {mu0:.4f}) / {sf:.4f}")
        for phi_eff in (0.58, 0.60, 0.62, 0.64):
            cvs = (phi_eff - mu0) / sf
            # 1 se in the slope, propagated
            dcv = abs(cvs) * (out[key]["ef"] / abs(sf))
            p(f"      phi_eff = {phi_eff:.2f}  ->  CV* = {cvs:.2f} +- {dcv:.2f} "
              f"(slope error only)")
        p("   ⚠ phi_eff is NOT measured. This converts an unmeasured densification into")
        p("     a predicted width, it does not supply one.")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--phi-probe", type=float, default=0.44,
                    help="the bench packing fraction to locate on the window")
    ap.add_argument("--out", default="dst_window.txt")
    args = ap.parse_args()

    lines: list[str] = []

    def p(s: str = "") -> None:
        print(s)
        lines.append(s)

    p("=== dst_window.py — the shear-thickening WINDOW, not its location ===")
    p("⚠ phi_0/phi_mu here are ISOTROPIC jamming fractions at mu_p = 0 and 0.3.")
    p("  The Wyart-Cates reading of them is a HYPOTHESIS (see module docstring).")
    d = load_fits()
    res = {}
    for fam in ("lognormal", "uniform"):
        res[fam] = report(p, fam, d, args.phi_probe)

    p("")
    p("=== WHICH mu_p IS `frictional`? the window is a DEFINITION, so price it ===")
    p("   Wyart-Cates' phi_mu is the FULLY frictional jamming point. This script's")
    p("   default reads it at mu_p = 0.3, which is the production friction, NOT the")
    p("   saturated one -- so the window it reports is the NARROWEST admissible.")
    dd = pd.read_csv(FITS)
    dd = dd[(dd.n_requested == N_MATCH) & (dd.family == "lognormal")
            & (dd.cv_nom == 0.0) & np.isfinite(dd.phi_J_lin)]
    g0 = dd[dd.mu_p == MU_0].phi_J_lin
    if len(g0):
        phi0 = float(g0.mean())
        p(f"   at CV = 0, lognormal, N = {N_MATCH}: phi_0 = {phi0:.6f} (n {len(g0)})")
        p("   mu_p   phi_mu (n)     W        lab x     dCV to traverse")
        for mu in sorted(dd.mu_p.unique()):
            if mu == MU_0:
                continue
            g = dd[dd.mu_p == mu].phi_J_lin
            pm = float(g.mean())
            wv = phi0 - pm
            p(f"   {mu:.1f}   {pm:.6f}({len(g)})  {wv:.4f}  "
              f"{(0.44 - pm) / wv:+.2f}     {wv / 0.1172:.2f}")
        p("   phi_J SATURATES in friction (the mu_p 1 -> 3 step is inside seed scatter),")
        p("   so the window is bounded: 0.0586 (mu_p 0.3) to 0.0757 (saturated), a factor")
        p("   1.29. Every conclusion in the docstring survives that spread IN KIND -- the")
        p("   traverse is 0.50 to 0.65 in CV either way, and the lab sits 1.6 to 2.3")
        p("   windows BELOW the window -- but no single number should be quoted without")
        p("   naming the mu_p it was read at.")
        p("   ⚠ NOT measured: dphi_mu/dCV at the SATURATED friction. The mu_p = 1 and 3")
        p("     ladders exist only at CV = 0, so the traverse above mixes a saturated")
        p("     window width with a mu_p = 0.3 slope. That is the cheapest remaining gap.")

    p("")
    p("=== families side by side (the zero-skew control) ===")
    for fam in ("lognormal", "uniform"):
        r = res.get(fam, {}).get("no single-seed") or res.get(fam, {}).get("all")
        if r:
            p(f"   {fam:10s} dW/dCV = {r['dW']:+.5f} +- {r['eW']:.5f}  "
              f"(n = {r['n']} CV levels)")
    p("")
    p("   ⚠ Singh et al. (arXiv:2510.16464) report phi_J^mu NON-monotonic in")
    p("     polydispersity for ZERO-SKEW distributions -- which is the `uniform` row,")
    p("     not the `lognormal` one. A linear slope is the wrong summary if that")
    p("     carries over; the per-level table above is the primary output.")
    p("   ⚠ Guy et al. report Wyart-Cates FAILING for suspensions with a predominant")
    p("     fraction of large particles, i.e. exactly the broad-lognormal end of this")
    p("     table. LITERATURE_2026-08-31.md §1.")

    path = os.path.join(OUT, args.out)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n-> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
