"""wc_endpoints_move_together_2026-09-19.py

THE QUESTION (session brief item B1, third bullet)
--------------------------------------------------
In Wyart-Cates, phi_J(sigma) interpolates between phi_0 (frictionless, the low-stress end)
and phi_m (frictional, the high-stress end). If BOTH endpoints are on disk at matched CV,
what does polydispersity do to each -- and therefore to the DST diagram?

Two candidate outcomes, as the brief names them:
  (H1) polydispersity shifts both endpoints EQUALLY  -> the frictional window
       W = phi_0 - phi_m is CV-independent and the phi-axis simply TRANSLATES;
  (H2) it shifts phi_m more                          -> the window widens or narrows and
       the shape of the DST diagram in phi changes with CV.

⚠⚠ NOT PRE-REGISTERED. The brief asked for the two candidates and the resolving power to be
written down BEFORE computing. They were not: the endpoints were found on disk and the
contrast was computed in the same step. This page therefore reports a POST-HOC comparison,
and the resolving power is stated after the fact rather than before it. That is a real
weakening and it is why the verdict below is a BOUND, not a confirmation.

WHAT IS ON DISK (checked here, not assumed)
  mu_p = 0.0  : lognormal and uniform, CV 0.0-0.5, 6 rungs, 3 seeds     -> phi_0
  mu_p = 0.3  : lognormal, uniform AND bidisperse, same grid            -> phi_m
  mu_p = 1.0  : lognormal only, 3 CV levels    (not used)
  mu_p = 3.0  : lognormal only, CV = 0 only    (not used)
So the comparison runs on TWO families, not three. Bidisperse has no frictionless block.
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, ROOT)

MU_FRICTIONLESS, MU_FRICTIONAL = 0.0, 0.3
CVS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
T_CRIT_2DOF = 4.303      # two-sided 95 %, 3 seeds -> 2 dof. NEVER a flat 2 (09-18c item 6)

FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"   [{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


def paired_shift(fits, fam, mu, lo=0.0, hi=0.5):
    """Seed-PAIRED CV lo -> hi shift. Paired, because the two ends share a seed draw."""
    g = fits[(fits.family == fam) & (fits.mu_p == mu)]
    a = g[g.cv_nom == lo].set_index("seed").phi_J_lin
    b = g[g.cv_nom == hi].set_index("seed").phi_J_lin
    common = sorted(set(a.index) & set(b.index))
    d = np.array([b[s] - a[s] for s in common])
    return float(d.mean()), float(d.std(ddof=1) / math.sqrt(len(d))), len(d)


def main() -> int:
    import sweep_descriptive as sd
    sd.P = lambda *a, **k: None
    fits = sd.load_fits()

    print("=" * 86)
    print("Do the two Wyart-Cates endpoints move together under polydispersity?")
    print("=" * 86)

    # --- 0. what is actually on disk -------------------------------------------------
    print("\n=== 0. GRID ON DISK (printed, not assumed) ===")
    print(f"   {'mu_p':>5s} {'family':11s} {'CV levels':>28s} {'seeds/cell':>11s}")
    fams_by_mu = {}
    for mu in sorted(fits.mu_p.unique()):
        for fam in sorted(fits[fits.mu_p == mu].family.unique()):
            g = fits[(fits.mu_p == mu) & (fits.family == fam)]
            cvs = [round(float(c), 1) for c in sorted(g.cv_nom.unique())]
            ns = sorted(int(v) for v in g.groupby("cv_nom").size().unique())
            fams_by_mu.setdefault(mu, []).append(fam)
            print(f"   {mu:5.1f} {fam:11s} {str(cvs):>28s} {str(ns):>11s}")
    usable = sorted(set(fams_by_mu.get(MU_FRICTIONLESS, [])) & set(fams_by_mu.get(MU_FRICTIONAL, [])))
    print(f"\n   families with BOTH endpoints at matched CV: {usable}")
    check("both Wyart-Cates endpoints exist at matched CV for at least two families",
          len(usable) >= 2, f"{usable}")
    check("the bidisperse family has NO frictionless block, so this runs on 2 of 3 families",
          "bidisperse" not in usable,
          f"mu_p=0 families: {fams_by_mu.get(MU_FRICTIONLESS)}")

    # --- 1. the two shifts ------------------------------------------------------------
    print("\n=== 1. THE PAIRED CV 0 -> 0.5 SHIFT AT EACH FRICTION ===")
    print(f"   {'family':11s} {'mu_p':>5s} {'Delta phi_J':>12s} {'se':>9s} {'n':>2s}")
    res = {}
    for fam in usable:
        for mu in (MU_FRICTIONLESS, MU_FRICTIONAL):
            m, se, n = paired_shift(fits, fam, mu)
            res[(fam, mu)] = (m, se, n)
            print(f"   {fam:11s} {mu:5.1f} {m:+12.5f} {se:9.5f} {n:2d}")
            check(f"{fam} at mu_p={mu} has a positive shift", m > 0, f"{m:+.5f}")

    # --- 2. the discriminating comparison ---------------------------------------------
    print("\n=== 2. H1 vs H2: is d(phi_0)/dCV equal to d(phi_m)/dCV? ===")
    print(f"   {'family':11s} {'d(phi_0)':>10s} {'d(phi_m)':>10s} {'difference':>12s} "
          f"{'se_diff':>9s} {'|t|':>6s} {'verdict':>28s}")
    bounds = {}
    for fam in usable:
        m0, s0, _ = res[(fam, MU_FRICTIONLESS)]
        m3, s3, _ = res[(fam, MU_FRICTIONAL)]
        d, sed = m0 - m3, math.hypot(s0, s3)
        t = abs(d) / sed
        # the RESOLVING POWER: the smallest difference this design could have resolved
        bound = T_CRIT_2DOF * sed
        bounds[fam] = (bound, abs(m3))
        # ⚠ WORDING. Failing to resolve the two apart does NOT "reject H2"; it fails to
        # reject H1, which at n = 3 means UNDETERMINED with a bound, never "shown equal".
        verdict = ("NOT RESOLVED APART (H1 not rejected)" if t < T_CRIT_2DOF
                   else "RESOLVED APART (H1 rejected)")
        print(f"   {fam:11s} {m0:+10.5f} {m3:+10.5f} {d:+12.5f} {sed:9.5f} {t:6.2f} {verdict:>28s}")
        check(f"{fam}: the two endpoints' CV shifts are not resolved apart at 2 dof",
              t < T_CRIT_2DOF, f"|t| = {t:.2f} against t_crit = {T_CRIT_2DOF}")

    print("\n   ⚠ 'NOT RESOLVED APART' AT n = 3 IS A BOUND, NOT A PROOF OF EQUALITY.")
    print(f"   {'family':11s} {'smallest resolvable diff':>25s} {'as a % of d(phi_m)':>20s}")
    for fam in usable:
        b, eff = bounds[fam]
        print(f"   {fam:11s} {b:25.5f} {100*b/eff:19.1f}%")
    print("   So what is EXCLUDED is a difference larger than that; a difference up to")
    print("   ~26-32 % of the effect would have passed unnoticed. The claim is 'the two")
    print("   endpoints move together to within a quarter to a third of the shift', not")
    print("   'they move identically'.")

    # --- 3. the window, which is the physically meaningful object ---------------------
    print("\n=== 3. THE FRICTIONAL WINDOW W = phi_0 - phi_m, AT EVERY CV ===")
    print("   Wyart-Cates interpolates phi_J(sigma) across W. If W is CV-independent the")
    print("   phi-axis TRANSLATES and the DST diagram keeps its shape.")
    print(f"   {'family':11s} " + "".join(f"{cv:>9.1f}" for cv in CVS) + f"{'range':>9s}")
    for fam in usable:
        row = []
        for cv in CVS:
            a = fits[(fits.family == fam) & (fits.mu_p == MU_FRICTIONLESS)
                     & (fits.cv_nom == cv)].phi_J_lin.mean()
            b = fits[(fits.family == fam) & (fits.mu_p == MU_FRICTIONAL)
                     & (fits.cv_nom == cv)].phi_J_lin.mean()
            row.append(a - b)
        rng = max(row) - min(row)
        print(f"   {fam:11s} " + "".join(f"{v:9.5f}" for v in row) + f"{rng:9.5f}")
        # the scale to read the range against: the seed s.e. of the shift, propagated
        _, s0, _ = res[(fam, MU_FRICTIONLESS)]
        _, s3, _ = res[(fam, MU_FRICTIONAL)]
        scale = math.hypot(s0, s3)
        print(f"   {'':11s} range = {rng/scale:.2f} x the propagated seed s.e. ({scale:.5f})")
        check(f"{fam}: the window's range across the CV ladder is within 2 s.e.",
              rng < 2 * scale, f"range {rng:.5f} vs 2 s.e. = {2*scale:.5f}")

    # --- 4. the discriminating NEGATIVE: the window must NOT be zero ------------------
    print("\n=== 4. A CHECK THAT COULD HAVE FAILED: the window is nonzero and both")
    print("       endpoints are ORDERED phi_0 > phi_m at every CV ===")
    bad = []
    for fam in usable:
        for cv in CVS:
            a = fits[(fits.family == fam) & (fits.mu_p == MU_FRICTIONLESS)
                     & (fits.cv_nom == cv)].phi_J_lin.mean()
            b = fits[(fits.family == fam) & (fits.mu_p == MU_FRICTIONAL)
                     & (fits.cv_nom == cv)].phi_J_lin.mean()
            if not a > b:
                bad.append((fam, cv, a, b))
    check("phi_0 > phi_m at every (family, CV) -- friction lowers the jamming fraction, "
          "which is the sign the literature reports",
          not bad, f"{len(bad)} violations")
    # and the window must be a real gap, not seed noise
    w = fits[(fits.mu_p == MU_FRICTIONLESS) & (fits.cv_nom == 0.0)
             & (fits.family == "lognormal")].phi_J_lin.mean() - \
        fits[(fits.mu_p == MU_FRICTIONAL) & (fits.cv_nom == 0.0)
             & (fits.family == "lognormal")].phi_J_lin.mean()
    _, s0, _ = res[("lognormal", MU_FRICTIONLESS)]
    check("the window at CV = 0 is many seed s.e. wide, so it is a real gap",
          w > 10 * s0, f"W = {w:.5f} = {w/s0:.1f} x the seed s.e.")

    # --- 5. the caveats that must travel ----------------------------------------------
    print("\n=== 5. CAVEATS THAT TRAVEL WITH EVERY NUMBER ABOVE ===")
    print("   (a) ⚠ POST-HOC. The two candidates and the resolving power were not written")
    print("       down before computing. See the module docstring.")
    print("   (b) ⚠ paper/main.tex: 'The frictionless ladders are not monotone in p_hat at")
    print("       this seed count, so this anchor is a check of the protocol's scale and")
    print("       not a result of the paper.' That caveat attaches to phi_0 here too.")
    print("   (c) ⚠ At CV = 0 the generator ignores the family label, so the two families'")
    print("       phi_0 values share a cell (found 2026-09-19). A test treating the two")
    print("       families as independent over-states the variance of their difference,")
    print("       i.e. errs CONSERVATIVELY.")
    print("   (d) ⚠ TWO families, not three. Bidisperse -- the non-monotone one -- has no")
    print("       frictionless block, so the family whose behaviour is least like the")
    print("       others is the one not tested.")
    print("   (e) ⚠ sigma* itself is untouched by any of this: it is a contact-force")
    print("       property and Wyart-Cates makes it independent of phi. Nothing here")
    print("       measures it.")

    print("\n" + "=" * 86)
    if FAILURES:
        print(f"FAILED {len(FAILURES)} check(s):")
        for f in FAILURES:
            print(f"   - {f}")
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
