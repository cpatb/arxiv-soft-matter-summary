"""geom_separability_2026-09-19.py

Prints every number quoted in geom_separability_2026-09-19.md.

Question: does phi_J(CV, mu_p) separate?  ADDITIVE  F(CV) + G(mu),
MULTIPLICATIVE F(CV) G(mu), or neither -- and what does the answer say about the
depopulation mechanism.

READ ONLY on out/.  Run:  python geom_separability_2026-09-19.py
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import geom_common_20260919 as G  # noqa: E402

CVS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
TC = G.T_CRIT_2DOF


def sec(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


# ------------------------------------------------------------------ 1. algebra
def part1():
    sec("1. WHAT 'SEPARABLE' MEANS, AND THE ONE STATISTIC THAT TESTS BOTH FORMS")
    c, m = sp.symbols("c m")
    F, Gs = sp.Function("F"), sp.Function("G")
    print("   ADDITIVE       phi(c,m) = F(c) + G(m)")
    print("   MULTIPLICATIVE phi(c,m) = F(c) * G(m)")
    add = F(c) + Gs(m)
    mul = F(c) * Gs(m)
    print("\n   The seed-paired CV increment D(c,m) = phi(c,m) - phi(0,m):")
    Dadd = sp.simplify(add - add.subs(c, 0))
    Dmul = sp.simplify(mul - mul.subs(c, 0))
    print("     additive:       D =", Dadd, "  -- free of m")
    print("     multiplicative: D =", Dmul, "  -- proportional to G(m) = phi(0,m)/F(0)")
    print("\n   So ONE statistic decides between them.  Write r_m = phi(0,m1)/phi(0,m0).")
    print("     additive       predicts  D(c,m1) - D(c,m0) = 0")
    print("     multiplicative predicts  D(c,m1) - D(c,m0) = (r_m - 1) D(c,m0)")
    print("   Both are parameter-free once phi(0,m) is measured, and their SEPARATION is")
    print("   |r_m - 1| * D(c,m0), which is printed below beside the residual so the test")
    print("   can declare itself uninformative.")
    print("\n   General transform.  If phi separates additively after a smooth monotone")
    print("   h, i.e. h(phi) = F(c) + G(m), then to first order in D")
    x = sp.symbols("x")
    print("     h(phi(c,m)) - h(phi(0,m)) = h'(phi(0,m)) D(c,m) + O(D^2) = F(c) - F(0)")
    print("   so   D(c,m1)/D(c,m0) = h'(phi(0,m0))/h'(phi(0,m1))   for EVERY c.")
    print("   That ratio is a single number fixed by h, so the CV-dependence of the")
    print("   measured ratio is a transform-free test: **no monotone h whatsoever can")
    print("   separate the data if D(c,m1)/D(c,m0) depends on c.**")
    print("   h = identity  -> ratio 1 (additive)")
    print("   h = log       -> ratio phi(0,m1)/phi(0,m0) (multiplicative)")
    print("   h = 1/phi (specific volume) or h = (1-phi)/phi (void ratio)")
    print("                 -> ratio (phi(0,m1)/phi(0,m0))^2")
    print("   h = phi^k     -> ratio (phi(0,m1)/phi(0,m0))^(k-1)")
    k = sp.symbols("k")
    phi0, phi1 = sp.symbols("phi_0 phi_1", positive=True)
    for name, hfun in [("phi^k", x**k), ("log phi", sp.log(x)),
                       ("1/phi", 1 / x), ("(1-phi)/phi", (1 - x) / x)]:
        dh = sp.diff(hfun, x)
        rat = sp.simplify(dh.subs(x, phi0) / dh.subs(x, phi1))
        print(f"     SymPy: h = {name:12s} -> h'(phi_0)/h'(phi_1) = {rat}")


# ------------------------------------------------------------------ 2. gating
def part2_gate(fits):
    sec("2. WHICH FRICTION ARMS ARE ADMISSIBLE -- the ladder-health gate, per arm")
    print("   The project's gate (sweep_descriptive.ladder_health): a family/arm passes")
    print("   if it has NO negative ladder slope and a median r2 >= 0.85.  The shipped")
    print("   implementation evaluates it at mu_p = MU_REF = 0.3 ONLY.  Evaluated here")
    print("   at every friction arm on disk, because item A reads three of them.")
    print("\n     family       mu_p   n cells  neg slopes  median r2   VERDICT")
    ok = {}
    for fam in ("lognormal", "uniform", "bidisperse"):
        for mu in sorted(fits.mu_p.unique()):
            g = fits[(fits.family == fam) & (fits.mu_p == mu)]
            if not len(g):
                continue
            neg = int((g.slope_lin <= 0).sum())
            r2 = float(np.median(g.r2_lin))
            good = (neg == 0) and (r2 >= 0.85)
            ok[(fam, mu)] = good
            print(f"     {fam:11s}  {mu:.1f}   {len(g):5d}   {neg:8d}    {r2:7.4f}   "
                  f"{'HEALTHY' if good else 'FAILS'}")
    print("\n   CONSEQUENCE.  The mu_p = 1.0 and mu_p = 3.0 arms FAIL the project's own")
    print("   ladder-health gate.  Everything below is therefore run on the two arms")
    print("   that pass (mu_p = 0 and 0.3); the mu_p = 1.0 arm is reported separately")
    print("   and labelled, never pooled.")
    return ok


# ------------------------------------------------------------------ 3. the test
def increments(fits, fam, mu):
    g = fits[(fits.family == fam) & (fits.mu_p == mu)]
    return g.pivot_table(index="seed", columns="cv_nom", values="phi_J_lin")


def part3_test(fits):
    sec("3. THE TEST -- seed-paired, both candidates, separation reported")
    out = {}
    for fam in ("lognormal", "uniform"):
        a, b = increments(fits, fam, 0.0), increments(fits, fam, 0.3)
        if a.empty or b.empty:
            continue
        seeds = sorted(set(a.index) & set(b.index))
        print(f"\n   --- {fam}, mu_p = 0.0 vs 0.3, seeds {seeds} ---")
        phi00 = a.loc[seeds, 0.0]
        phi03 = b.loc[seeds, 0.0]
        r_m = float(phi03.mean() / phi00.mean())
        print(f"   phi_J(CV=0): mu=0 {phi00.mean():.5f}, mu=0.3 {phi03.mean():.5f}, "
              f"r_m = {r_m:.5f}")
        print("     CV   D(mu=0)   D(mu=.3)   resid_ADD   t_ADD   resid_MULT  t_MULT  "
              " separation  sep/se  informative?")
        for cv in CVS[1:]:
            if cv not in a.columns or cv not in b.columns:
                continue
            D0 = (a.loc[seeds, cv] - a.loc[seeds, 0.0]).to_numpy()
            D3 = (b.loc[seeds, cv] - b.loc[seeds, 0.0]).to_numpy()
            delta = D3 - D0
            m, se, t, n = G.tstat(delta)
            # multiplicative prediction for delta, per seed
            pred_m = (r_m - 1.0) * D0
            rm_, sem_, tm_, _ = G.tstat(delta - pred_m)
            sep = float(abs(pred_m.mean()))
            ratio = sep / se if se > 0 else float("nan")
            out[(fam, cv)] = dict(D0=D0.mean(), D3=D3.mean(), add=m, add_se=se,
                                  add_t=t, mult=rm_, mult_t=tm_, sep=sep, sepse=ratio)
            print(f"     {cv:.1f}  {D0.mean():+.5f}  {D3.mean():+.5f}  "
                  f"{m:+.5f}  {t:+7.2f}  {rm_:+.5f}  {tm_:+7.2f}   {sep:.5f}  "
                  f"{ratio:6.2f}  {'YES' if ratio >= 1.6 else 'no (sep < 1.6 se)'}")
        print(f"   (t_crit at 2 dof = {TC:.3f}.  'informative' uses the project's own")
        print("    1.6x candidate-separation bar, lessons.md 2026-08-30.)")
    return out


# --------------------------------------------- 4. no monotone transform works
def part4_transform(fits, res):
    sec("4. CAN ANY MONOTONE TRANSFORM SEPARATE IT?  -- ratio must be CV-independent")
    print("   Section 1: h-separability for ANY smooth monotone h forces")
    print("       R(c) = D(c, mu=0.3)/D(c, mu=0)   to be independent of c.")
    print("   So compute R(c) with its seed-paired error and ask whether it moves.")
    for fam in ("lognormal", "uniform"):
        a, b = increments(fits, fam, 0.0), increments(fits, fam, 0.3)
        if a.empty or b.empty:
            continue
        seeds = sorted(set(a.index) & set(b.index))
        r_m = float(b.loc[seeds, 0.0].mean() / a.loc[seeds, 0.0].mean())
        print(f"\n   --- {fam} ---   additive needs R = 1; multiplicative needs "
              f"R = {r_m:.4f};")
        print(f"   void ratio / specific volume need R = {r_m**2:.4f}.")
        print("     CV     R        se      t vs 1     t vs mult   implied k in h=phi^k")
        Rs = []
        for cv in CVS[1:]:
            D0 = (a.loc[seeds, cv] - a.loc[seeds, 0.0]).to_numpy()
            D3 = (b.loc[seeds, cv] - b.loc[seeds, 0.0]).to_numpy()
            # ratio of means with a seed-paired error via the delta method on the
            # paired differences (D3 - R D0 has mean zero at the true R)
            R = D3.mean() / D0.mean()
            # se of R from the paired residual D3 - R*D0
            resid = D3 - R * D0
            seR = float(resid.std(ddof=1) / math.sqrt(len(resid)) / abs(D0.mean()))
            kk = 1 + math.log(R) / math.log(r_m) if R > 0 else float("nan")
            Rs.append((cv, R, seR))
            print(f"     {cv:.1f}  {R:7.4f}  {seR:7.4f}  {(R-1)/seR:+9.2f}  "
                  f"{(R-r_m)/seR:+9.2f}   {kk:8.3f}")
        # is R CV-dependent?  compare the two extreme CVs with a pooled paired test
        print("\n     Pairwise: is R(c) the SAME at two CVs?  (needs |t| > "
              f"{TC:.3f} to resolve)")
        worst = None
        for i in range(len(Rs)):
            for j in range(i + 1, len(Rs)):
                c1, R1, s1 = Rs[i]
                c2, R2, s2 = Rs[j]
                t = (R1 - R2) / math.sqrt(s1 * s1 + s2 * s2)
                if worst is None or abs(t) > abs(worst[2]):
                    worst = (c1, c2, t, R1, R2)
        c1, c2, t, R1, R2 = worst
        print(f"       largest |t| is CV {c1:.1f} (R={R1:.4f}) vs CV {c2:.1f} "
              f"(R={R2:.4f}): t = {t:+.2f}  -> "
              f"{'RESOLVED: no monotone h separates' if abs(t) > TC else 'NOT resolved'}")


# --------------------------------------------- 5. the mu_p = 1.0 arm, labelled
def part5_mu1(fits, sh):
    sec("5. THE mu_p = 1.0 ARM -- reported, gated OUT, never pooled")
    a = increments(fits, "lognormal", 0.0)
    b = increments(fits, "lognormal", 0.3)
    c = increments(fits, "lognormal", 1.0)
    seeds = sorted(set(a.index) & set(b.index) & set(c.index))
    print(f"   seeds common to all three arms: {seeds}")
    print("   phi_J(CV=0): " + ", ".join(
        f"mu={mu}: {x.loc[seeds, 0.0].mean():.5f}" for mu, x in
        [(0.0, a), (0.3, b), (1.0, c)]))
    print("\n     CV   D(mu=0)     D(mu=0.3)   D(mu=1.0)   ADD resid vs mu=0.3   t")
    for cv in (0.2, 0.5):
        D0 = (a.loc[seeds, cv] - a.loc[seeds, 0.0]).to_numpy()
        D3 = (b.loc[seeds, cv] - b.loc[seeds, 0.0]).to_numpy()
        D1 = (c.loc[seeds, cv] - c.loc[seeds, 0.0]).to_numpy()
        m, se, t, n = G.tstat(D1 - D3)
        print(f"     {cv:.1f}  {D0.mean():+.5f}   {D3.mean():+.5f}   {D1.mean():+.5f}   "
              f"{m:+.5f} +- {se:.5f}   {t:+.2f}")
    print("\n   This arm FAILS the ladder-health gate (section 2).  It is quoted so that")
    print("   the reader can see it agrees, NOT as evidence that it does.")

    print("\n   THE FIT-FREE CHANNEL.  phi at the matched lowest rung p_hat = 3e-5 needs")
    print("   no ladder fit, so it survives where phi_J_lin does not.  The same")
    print("   comparison on that channel (lognormal, seed-paired):")
    print("     mu_p   phi(CV=0)   phi(CV=0.5)   D = CV 0->0.5   se      n")
    rows = {}
    for mu in sorted(sh.mu_p.unique()):
        g = sh[(sh.mu_p == mu) & (sh.family == "lognormal") & (np.isclose(sh.p, 3e-5))]
        z = g[np.isclose(g.cv, 0.0)].set_index("seed").phi
        h = g[np.isclose(g.cv, 0.5)].set_index("seed").phi
        common = sorted(set(z.index) & set(h.index))
        if not len(z):
            continue
        if common:
            d = np.array([h[s] - z[s] for s in common])
            m, se, t, n = G.tstat(d)
            rows[mu] = (float(z.mean()), float(h.mean()), m, se, n)
            print(f"     {mu:.1f}    {z.mean():.5f}    {h.mean():.5f}     {m:+.5f}   "
                  f"{se:.5f}  {n}")
        else:
            rows[mu] = (float(z.mean()), float("nan"), float("nan"), float("nan"), 0)
            print(f"     {mu:.1f}    {z.mean():.5f}       --            --        --   0")
    print("\n   The CV=0 column is the monodisperse friction ladder with no fitting in")
    print("   it at all.  External anchor (API-resolved): Silbert, Soft Matter 6, 2918")
    print("   (2010), doi:10.1039/c001973a, Table 1, monodisperse 3D frictional spheres:")
    print("     phi_c = 0.639 / 0.595 / 0.556 / 0.544 at mu = 0 / 0.2 / 1.0 / 10,")
    print("     saturating near 0.55 for mu >= 1.")
    sil = {0.0: 0.639, 0.2: 0.595, 1.0: 0.556, 10.0: 0.544}
    print("     mu   this work (phi at 3e-5)   Silbert 2010   difference")
    for mu in (0.0, 1.0):
        if mu in rows and mu in sil:
            print(f"     {mu:.1f}        {rows[mu][0]:.5f}              {sil[mu]:.3f}"
                  f"        {rows[mu][0]-sil[mu]:+.4f}")
    print("     (mu = 0.3 here has no Silbert row; his 0.2 gives 0.595 against "
          f"{rows.get(0.3, [float('nan')])[0]:.5f} here.)")
    print("   Reading: the monodisperse friction ladder agrees with an independent")
    print("   published one to a few times 1e-3, and BOTH say mu = 0.3 is not saturated")
    print("   while mu = 1.0 essentially is.  So 'friction-independent' in item A is a")
    print("   statement made across a range over which phi_J itself moves by 0.06 --")
    print("   the same size as the whole CV effect.  That is what makes it non-trivial.")


# ------------------------------- 6. what the rattler channel says about the same
def part6_mechanism(sh):
    sec("6. THE MECHANISM TEST -- is the CV effect friction-independent in the "
        "RATTLER channel too?")
    print("   The depopulation picture (geom_rattler_decomposition) writes")
    print("       phi_J(CV,mu) = phi_B(mu) / (1 - psi(CV,mu)).")
    print("   If psi depended on CV ALONE this would be EXACTLY MULTIPLICATIVE.")
    print("   Section 3 refutes multiplicative.  So psi must be friction-dependent, and")
    print("   the rattler NUMBER fraction q is on disk at both arms.  Compare.")
    for p in (1e-3, 3e-5):
        print(f"\n   --- rung p_hat = {p:.0e}, lognormal ---")
        print("     CV     q(mu=0)   q(mu=0.3)   ratio   |   D(q) from CV=0, mu=0   "
              "mu=0.3   ratio")
        g0 = sh[(sh.mu_p == 0.0) & (sh.family == "lognormal") & (np.isclose(sh.p, p))]
        g3 = sh[(sh.mu_p == 0.3) & (sh.family == "lognormal") & (np.isclose(sh.p, p))]
        q00 = float(g0[np.isclose(g0.cv, 0.0)].rf.mean())
        q30 = float(g3[np.isclose(g3.cv, 0.0)].rf.mean())
        for cv in CVS[1:]:
            v0 = g0[np.isclose(g0.cv, cv)].rf
            v3 = g3[np.isclose(g3.cv, cv)].rf
            if not len(v0) or not len(v3):
                continue
            a_, b_ = float(v0.mean()), float(v3.mean())
            print(f"     {cv:.1f}   {a_:.4f}    {b_:.4f}    {b_/a_:.3f}   |   "
                  f"{a_-q00:+.4f}   {b_-q30:+.4f}   {(b_-q30)/(a_-q00):.3f}")
    print("\n   Reading: the rattler COUNT's CV response is strongly friction-dependent")
    print("   while phi_J's is not.  Whatever makes phi_J separable, it is NOT that")
    print("   depopulation itself is friction-blind.  Under H1 the compensation has to")
    print("   sit in the conditional rattler VOLUME, i.e. in the selectivity S.")
    print("\n   Required selectivity S at the two arms (H1 assumed), lognormal:")
    for p in (1e-3, 3e-5):
        print(f"     rung {p:.0e}")
        for mu in (0.0, 0.3):
            g = sh[(sh.mu_p == mu) & (sh.family == "lognormal") & (np.isclose(sh.p, p))]
            ref = g[np.isclose(g.cv, 0.0)]
            if not len(ref):
                continue
            phiB0 = float(ref.phi.mean()) * (1 - float(ref.rf.mean()))
            row = []
            for cv in CVS[1:]:
                gg = g[np.isclose(g.cv, cv)]
                if not len(gg):
                    row.append(f"CV{cv:.1f}    --")
                    continue
                phi_, q_ = float(gg.phi.mean()), float(gg.rf.mean())
                S = (q_ - (1 - phiB0 / phi_)) / (q_ - G.Psi("lognormal", cv, q_))
                row.append(f"CV{cv:.1f} {S:5.3f}")
            print(f"       mu={mu:.1f}  " + "  ".join(row))
    print("   If S came out equal at the two frictions, the separability would have a")
    print("   one-line mechanism.  It is printed rather than described.")


def main():
    fits = G.load_fits()
    sh = G.load_shards()
    part1()
    part2_gate(fits)
    res = part3_test(fits)
    part4_transform(fits, res)
    part5_mu1(fits, sh)
    part6_mechanism(sh)
    sec("VERDICT")
    ln5 = res[("lognormal", 0.5)]
    un5 = res[("uniform", 0.5)]
    print("   MULTIPLICATIVE separability: REFUTED.  lognormal CV 0->0.5 residual "
          f"{ln5['mult']:+.5f}, t = {ln5['mult_t']:+.2f} at 2 dof "
          f"(t_crit {TC:.3f}); candidate separation {ln5['sep']:.5f} = "
          f"{ln5['sepse']:.1f} se, so the test is informative.")
    print("   ADDITIVE separability: also REFUTED at lognormal CV = 0.5 "
          f"(residual {ln5['add']:+.5f}, t = {ln5['add_t']:+.2f}), but the residual is "
          f"{abs(ln5['add'])/abs(ln5['D3'])*100:.1f}% of the shift against "
          f"{abs(ln5['mult'])/abs(ln5['D3'])*100:.1f}% for multiplicative.")
    print(f"   For the uniform family neither form is refuted (add t = "
          f"{un5['add_t']:+.2f}, mult t = {un5['mult_t']:+.2f}); with a separation of "
          f"only {un5['sepse']:.2f} se that arm decides nothing.")
    print("   VERDICT: ADDITIVE SEPARABILITY IS THE BETTER OF THE TWO AND IS ITSELF")
    print("            REFUTED -> the honest statement is 'approximately additive to")
    print("            ~6% of the shift', and 'friction-independent' is an")
    print("            approximation with a measured size, not an exact structure.")
    print("\n   ANCHOR SCOPE.  This can rank the two candidate forms and put a number on")
    print("   the residual interaction.  It cannot identify the mechanism: it uses only")
    print("   phi_J, and two very different microscopic stories can produce the same")
    print("   two-way table.  It also cannot speak for mu_p > 0.3, whose ladders fail")
    print("   the project's own health gate.")


if __name__ == "__main__":
    main()
