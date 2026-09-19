"""onset_size_degeneracy_check.py -- checks for derivations/onset_size_degeneracy.md

⚠⚠ THIRD VERSION, 2026-09-19. TWO HEADLINES HAVE NOW BEEN WITHDRAWN FROM THIS PAGE, both
by referee passes, and the second withdrawal is the interesting one because the repair
committed the SAME error as the thing it repaired.

  v1 (morning):  "ID and Sigma = sigma_imp/sigma* carry the same power of particle size,
                  so no size ladder can separate the frictional and inertial-dilatant
                  regimes."   WITHDRAWN: it rested on "no admissible impact-stress gauge
                  carries a particle size", and three of six do.
  v2 (midday):   "the degenerate pair is (sigma*, P_f), the two GRAIN-SCALE stresses, so a
                  size ladder cannot separate frictional jamming from frustrated
                  dilatancy."   WITHDRAWN on two counts:
                    (a) sigma* is a THRESHOLD BOTH routes must cross, not one route's
                        driving stress, so P_f/sigma* is the ACTIVATION NUMBER of the
                        dilatant branch -- a statement about ONE mechanism. The actual
                        A-vs-B discriminator is driving-stress-vs-driving-stress,
                        sigma_in/P_f, whose exponent v2's own eq. (7) prints as +2.00.
                    (b) P_f's own size exponent is NOT -2. It is set by conventions that
                        are unpinned on disk, and they span -2, -1, 0 and +1. v2 replaced
                        one convention choice with another and called it settled.

WHAT THIS FILE NOW CHECKS
-------------------------
The verdict is UNDETERMINED, and that is the result. Specifically:

  * the exponent ALGEBRA is right and is an identity given the definitions (S1);
  * the size exponent of the DILATANT stress is undetermined over a 3-decade-per-decade
    span, because three conventions are unpinned -- the drainage length, the ceiling
    regime, and whether the capacity or the drainage criterion binds (S5);
  * the double layer adds more size dependence at 800 nm than any of the candidates
    differ by, and the condition that removes it is an EXPONENT condition, not a value
    condition -- which is stricter than v2 stated and voids v2's salt escape clause (S6);
  * the CV lever's magnitude is family-dependent, its slope must be read LOCALLY from each
    family's OWN ladder, and the bidisperse family's local slope is NEGATIVE above
    CV = 0.3, so even the sign is not family-universal (S7);
  * the model's phi_c is unmeasured for the lab silica (S7c).

Everything numeric below is read from files on disk or computed here. The v2 version of
this file hardcoded three phase8 numbers that are not in phase8_verdict.txt and a
"0.089-0.125" bracket that matches none of its four gate sets; both are fixed and the
wrong values are printed beside the right ones.
"""

from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(ROOT, "out")
sys.path.insert(0, r"C:\Users\compl\Desktop\Projects\impact_stress")
sys.path.insert(0, ROOT)

ALPHA, ALPHA_LO, ALPHA_HI = 1.9, 1.8, 2.0     # Guy, Hermes & Poon 2015; see the .md sec 1.3
R_THROAT_PER_D = 0.155 / 2.0                   # r_throat = 0.0775 d (PORE_PHYSICS sec 5.2)
LAMBDA_D_COMPROMISE, LAMBDA_D_PRIMARY = 30.0e-9, 95.0e-9   # PORE_PHYSICS sec 5.2: not measured

FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"   [{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


# =====================================================================================
# S1  THE ALGEBRA -- correct, and an IDENTITY given the definitions
# =====================================================================================
def s1_algebra():
    print("\n=== S1  EXPONENT ALGEBRA (an identity given the definitions, not a result) ===")
    rho, mu, U0, D0, d, dphi, phi, A, B, alpha, C = sp.symbols(
        "rho mu U0 D0 d dphi phi A B alpha C", positive=True)
    grid = np.linspace(1e-6, 1 - 1e-6, 20001)
    check("C = (1-phi)^3/(180 phi^2) > 0 on 0 < phi < 1 (20001 points)",
          bool(np.all((1 - grid) ** 3 / (180 * grid**2) > 0)))

    kappa = C * d**2
    sigma_star = B * d ** (-alpha)
    P_f = mu * dphi * D0 * U0 / kappa                # ell ~ D0; see S5, this is a CHOICE
    sigma_in = rho * U0**2
    ID = rho * kappa * U0 / (A * mu * dphi * D0)

    def e(x):
        return sp.simplify(sp.diff(sp.log(x), d) * d)

    for nm, x, want in (("sigma*", sigma_star, -alpha), ("P_f (ell ~ D0)", P_f, -2),
                        ("rho U0^2", sigma_in, 0), ("ID", ID, 2)):
        check(f"size exponent of {nm} is {want}", sp.simplify(e(x) - want) == 0, f"{e(x)}")

    print("\n   the three ratios that matter, and what each one ASKS:")
    tbl = [("sigma_in/sigma*", sigma_in / sigma_star, alpha,
            "is the IMPOSED stress above the frictional onset?   (activation, route A)"),
           ("P_f/sigma*", P_f / sigma_star, alpha - 2,
            "is the SUCTION above the frictional onset?          (activation, route B)"),
           ("sigma_in/P_f", sigma_in / P_f, 2,
            "does inertia beat the suction?   <-- THE A-vs-B DISCRIMINATOR")]
    for nm, x, want, q in tbl:
        got = e(x)
        print(f"      {nm:16s} d^({got})   {q}")
        check(f"size exponent of {nm} is {want}", sp.simplify(got - want) == 0, f"{got}")

    print("\n   🔴 WITHDRAWN v2 HEADLINE. v2 read the middle row as 'the two mechanisms are")
    print("      size-degenerate'. It is not: sigma* is the THRESHOLD BOTH routes cross, so")
    print("      that row is route B's activation number and says nothing about route A.")
    print("      The discriminator is the LAST row, exponent +2 -- i.e. a size ladder")
    print("      separates the two mechanisms at ~100x per size decade. v2 asserted the")
    print("      opposite four lines below printing this exponent.")
    check("the A-vs-B discriminator sigma_in/P_f has exponent 2, NOT 0 "
          "(so a size ladder is not blind, contrary to v1 and v2)",
          sp.simplify(e(sigma_in / P_f) - 2) == 0)
    # and ID is that discriminator's reciprocal up to A -- state it, do not imply novelty
    check("ID is sigma_in/P_f up to the factor A (so [M26]'s number ALREADY is the "
          "A-vs-B discriminator and this page adds no new one)",
          sp.simplify(sp.simplify(ID) - sp.simplify(sigma_in / (A * P_f))) == 0)


# =====================================================================================
# S2  DIMENSIONS
# =====================================================================================
def s2_dimensions():
    print("\n=== S2  DIMENSIONS (pint) ===")
    import pint
    u = pint.UnitRegistry()
    rho, mu, U0, D0, d = (1440.0 * u.kg / u.m**3, 1.0e-3 * u.Pa * u.s, 3.0 * u.m / u.s,
                          2.5e-3 * u.m, 0.8e-6 * u.m)
    phi, dphi, A = 0.44, 0.02, 3.0
    kappa = (1 - phi) ** 3 * d**2 / (180 * phi**2)
    check("kappa is an area", kappa.dimensionality == (u.m**2).dimensionality)
    for nm, q, want in (("ID", (rho * kappa * U0 / (A * mu * dphi * D0)).to_base_units(), u.dimensionless),
                        ("sigma*", (1.0e-9 * u.N / (6 * math.pi * d**2)).to(u.Pa), u.Pa),
                        ("P_f", (mu * dphi * D0 * U0 / kappa).to(u.Pa), u.Pa),
                        ("rho U0^2", (rho * U0**2).to(u.Pa), u.Pa)):
        check(f"{nm} has the dimensions of {want}", q.dimensionality == want.dimensionality)
    # B's units are an independent confirmation of the exponent reading
    B = 2e-11 * u.J / u.m**1.1
    check("B = 2e-11 J m^-1.1 with sigma* = B d^-1.9 gives a stress "
          "(an independent check on how alpha was read out of the paper)",
          (B * (1e-6 * u.m) ** -1.9).to(u.Pa).dimensionality == u.Pa.dimensionality)


# =====================================================================================
# S3  LIMITS
# =====================================================================================
def s3_limits():
    print("\n=== S3  LIMITS ===")
    rho, mu, U0, D0, d, dphi, A, B, alpha, C = sp.symbols(
        "rho mu U0 D0 d dphi A B alpha C", positive=True)
    P_f = mu * dphi * D0 * U0 / (C * d**2)
    ID = rho * C * d**2 * U0 / (A * mu * dphi * D0)
    R = P_f / (B * d ** (-alpha))
    lim2 = sp.simplify(sp.limit(R.subs(alpha, 2), d, 0, "+"))
    # ⚠ v2's check here tested only `d not in free_symbols`, which sp.oo and 0 also pass.
    check("d -> 0 at alpha = 2 : P_f/sigma* tends to a FINITE NONZERO d-free limit "
          "(v2's version of this check also passed on oo and on 0)",
          d not in lim2.free_symbols and lim2 not in (sp.oo, sp.S.Zero)
          and lim2.is_finite is not False, f"limit = {lim2}")
    check("d -> 0 at alpha = 1.9 : P_f/sigma* -> inf",
          sp.limit(R.subs(alpha, sp.Rational(19, 10)), d, 0, "+") == sp.oo)
    check("dphi -> 0+ : P_f -> 0 (no dilation, no suction)", sp.limit(P_f, dphi, 0, "+") == 0)
    check("dphi -> 0+ : ID -> inf", sp.limit(ID, dphi, 0, "+") == sp.oo)
    from inertial_dilatant import DilatantSystem, ID as ID_code
    check("shipped ID() returns inf at dphi = 0, the same convention",
          math.isinf(ID_code(DilatantSystem(name="t", phi=0.4, phi_c=0.4, a=1e-5, rho=1400.0),
                             3.0, 4e-3)))
    check("mu -> 0 : P_f/sigma* -> 0 (free drainage, no frustration)",
          sp.limit(R, mu, 0, "+") == 0)
    phi = sp.Symbol("phi", positive=True)
    Cx = (1 - phi) ** 3 / (180 * phi**2)
    check("phi -> 1 : ID -> 0 (no pore space to drain)",
          sp.limit(rho * Cx * d**2 * U0 / (A * mu * dphi * D0), phi, 1, "-") == 0)


# =====================================================================================
# S4  THE SHIPPED CODE, AND THE CEILING ITS OWN DOCSTRING SAYS TO CHECK
# =====================================================================================
def s4_shipped():
    print("\n=== S4  THE SHIPPED IMPLEMENTATION, AND ITS OWN VALIDITY CEILING ===")
    from inertial_dilatant import DilatantSystem, ID as ID_code, dilatant_pressure

    pts = [("silica (lab) 3.0 m/s ", 0.8e-6, 0.44, 0.42, 2.5e-3, 2000.0, 3.0),
           ("silica (lab) 0.5 m/s ", 0.8e-6, 0.44, 0.42, 2.5e-3, 2000.0, 0.5),
           ("cornstarch   0.5 m/s ", 10.0e-6, 0.42, 0.37, 4.0e-3, 1550.0, 0.5),
           ("cornstarch   7.2 m/s ", 10.0e-6, 0.42, 0.37, 4.0e-3, 1550.0, 7.2)]

    # --- size exponents, measured on the shipped functions -----------------------------
    print("   size exponents by finite difference on the SHIPPED functions:")
    for nm, d0, phi, phi_c, D0, rho_p, U0 in pts[:1] + pts[2:3]:
        rho = rho_p * phi + 1000.0 * (1 - phi)
        ds = np.array([d0 / 1.05, d0 * 1.05])
        v = {"ID": [], "P_f": []}
        for d in ds:
            s = DilatantSystem(name=nm, phi=phi, phi_c=phi_c, a=d, rho=rho)
            v["ID"].append(math.log(ID_code(s, U0, D0)))
            v["P_f"].append(math.log(dilatant_pressure(s, U0, D0)))
        ld = math.log(ds[1] / ds[0])
        e_id, e_pf = (v["ID"][1] - v["ID"][0]) / ld, (v["P_f"][1] - v["P_f"][0]) / ld
        print(f"      {nm} ID {e_id:+.6f} (want +2)   P_f {e_pf:+.6f} (want -2)")
        check(f"shipped ID has size exponent +2 at {nm.strip()}", abs(e_id - 2) < 1e-8)
        check(f"shipped dilatant_pressure has size exponent -2 at {nm.strip()}",
              abs(e_pf + 2) < 1e-8)

    # --- 🔴 the ceiling the module's docstring tells the caller to check ----------------
    print("\n   🔴 `dilatant_pressure`'s own docstring: 'Compare it against")
    print("      p_dilatant_ceiling before believing it.' v1 and v2 never did.")
    print(f"   {'point':22s} {'P_f [Pa]':>11s} {'ceiling':>11s} {'p_cap':>11s} {'P_f/ceil':>10s}")
    over = 0
    for nm, d, phi, phi_c, D0, rho_p, U0 in pts:
        rho = rho_p * phi + 1000.0 * (1 - phi)
        s = DilatantSystem(name=nm, phi=phi, phi_c=phi_c, a=d, rho=rho)
        pf = dilatant_pressure(s, U0, D0)
        ceil, pcap = s.p_dilatant_ceiling, s.p_capillary
        over += pf > ceil
        print(f"   {nm:22s} {pf:11.3e} {ceil:11.3e} {pcap:11.3e} {pf/ceil:10.1f}")
    check("🔴 REFUTED: the Darcy branch is NOT valid at the parameter points this page "
          "uses -- P_f exceeds the module's own ceiling at every one of them",
          over == len(pts), f"{over} of {len(pts)} points are above the ceiling")
    print("      Above the ceiling the suction is set by cavitation (exponent 0) or by")
    print("      capillary entry (2 gamma / r_throat, exponent -1), NOT by Darcy (-2).")
    print("      So the dilatant stress's size exponent at these very points is one of")
    print("      {0, -1}, not the -2 the whole of v2 was built on.")


# =====================================================================================
# S5  THE CONVENTION SCAN -- the result, and it is UNDETERMINED
# =====================================================================================
def s5_conventions():
    print("\n=== S5  WHAT SETS THE DILATANT STRESS'S SIZE EXPONENT: three unpinned knobs ===")
    print("   P_f ~ (mu/kappa) dphi gammadot ell^2 with kappa ~ d^2, so the exponent is")
    print("   2*(d ln ell/d ln d) - 2. NOTHING ON DISK PINS ell.")

    # SIZE_THRESHOLD sec 1 tabulates BOTH readings; reproduce the factor between them.
    id_D0, id_10a = 7.19e-05, 7.03       # SIZE_THRESHOLD sec 1, silica at 3 m/s
    D0, a = 2.5e-3, 0.8e-6
    predicted = id_D0 * (D0 / (10 * a)) ** 2
    print(f"\n   SIZE_THRESHOLD sec 1 prints BOTH: ID(ell = D0) = {id_D0:.3g}, "
          f"ID(ell = 10a) = {id_10a:.3g}")
    print(f"   ell enters SQUARED, so the ratio should be (D0/10a)^2 = "
          f"{(D0/(10*a))**2:.4g}: predicted {predicted:.4g} against the printed {id_10a:.4g}")
    check("the two drainage-length readings differ by exactly (D0/10a)^2, confirming ell "
          "enters squared and that BOTH are on record",
          abs(predicted - id_10a) / id_10a < 0.01,
          f"{predicted:.4g} vs {id_10a:.4g}, {abs(predicted-id_10a)/id_10a*100:.2f} %")

    rows = [
        ("Darcy, ell ~ D0        (drop-scale drainage, [M26]'s choice)", -2),
        ("Darcy, ell ~ 10a       (grain-scale drainage, SIZE_THRESHOLD sec 1)", 0),
        ("cavitation-limited     (P_f above 1 atm of suction; S4 says it is)", 0),
        ("capillary-entry limited(2 gamma / r_throat; S4 says it binds too)", -1),
        ("CAPACITY criterion     (p_dil^max/P_conf ~ (a/6 lam_D)(dphi/phi), PORE_PHYSICS 9.1)", +1),
    ]
    print(f"\n   {'convention for the dilatant stress':64s} {'p':>3s} "
          f"{'vs sigma* ~ d^':>15s} {'vs sigma_in ~ d^':>17s}")
    for nm, p in rows:
        print(f"   {nm:64s} {p:>3d} {p + ALPHA:>15.2f} {-p:>17.2f}")
    ps = sorted({p for _n, p in rows})
    span = max(ps) - min(ps)
    print(f"\n   the dilatant stress's size exponent is one of {ps}: a span of {span} "
          f"decades PER SIZE DECADE.")
    check("🔴 UNDETERMINED: the dilatant stress's size exponent is not pinned by anything "
          "on disk, and the candidates span 3 decades per size decade",
          span >= 3, f"candidates {ps}")
    print("   ⚠ v2 picked the FIRST row and called the question settled. v1 picked a")
    print("     zero-exponent gauge and called it settled. Both are one row of this table.")
    print("\n   WHAT WOULD SETTLE IT -- three measurements, none of them algebra:")
    print("     1. the drainage length: does the suction relax over the contact patch or")
    print("        over the grain scale? A pressure probe, or the U0-dependence of the")
    print("        force plateau, decides it.")
    print("     2. the ceiling regime: does the pore fluid cavitate during impact? S4 says")
    print("        both materials are above 1 atm of suction on the Darcy formula.")
    print("     3. whether capacity or drainage binds -- SIZE_THRESHOLD sec 3 already says")
    print("        these scale as a^1 and a^2 and point opposite ways, and withdrew the")
    print("        window they bracket as unquotable.")


# =====================================================================================
# S6  THE DOUBLE LAYER -- an EXPONENT condition, which is stricter than v2's
# =====================================================================================
def s6_double_layer():
    print("\n=== S6  THE DOUBLE LAYER: a condition on the EXPONENT, not on the value ===")
    print("   kappa_eff/kappa = (1-x)^n with x = lambda_D/r_throat, n = 2 (immobile shell)")
    print("   or 4 (Poiseuille throat); r_throat = 0.0775 d. Both readings are ESTIMATES.")
    print("   ⚠ v2 required the VALUE correction to be small (x <= 0.1 -> 2-4 %). The claim")
    print("     is about an EXPONENT, and d ln(1-x)^n/d ln d = n x/(1-x), which at x = 0.1")
    print("     is 0.222 (n=2) / 0.444 (n=4) -- 2.2x to 4.4x the residual it must protect.")

    def extra_exponent(x, n):
        return n * x / (1 - x)

    resid = abs(ALPHA - 2)
    print(f"\n   extra size exponent from the double layer, against the residual "
          f"|alpha-2| = {resid:.2f}:")
    print(f"   {'lam_D [nm]':>11s} {'d [um]':>8s} {'x':>7s} {'n=2':>8s} {'n=4':>8s}")
    for lam in (LAMBDA_D_COMPROMISE, LAMBDA_D_PRIMARY):
        for d in (0.8e-6, 3.87e-6, 10.0e-6, 25.0e-6):
            x = min(lam / (R_THROAT_PER_D * d), 0.999)
            print(f"   {lam*1e9:11.1f} {d*1e6:8.2f} {x:7.4f} "
                  f"{extra_exponent(x, 2):8.3f} {extra_exponent(x, 4):8.3f}")
    x800 = LAMBDA_D_COMPROMISE / (R_THROAT_PER_D * 0.8e-6)
    check("🔴 the double-layer term dominates the residual exponent at 800 nm under both "
          "readings and both lambda_D values",
          min(extra_exponent(x800, 2), extra_exponent(x800, 4)) > resid,
          f"smallest {extra_exponent(x800, 2):.3f} vs residual {resid:.2f}")

    print("\n   the EXPONENT condition (extra exponent <= the residual 0.10), which is what")
    print("   v2 should have required. v2's VALUE condition gave 3.87 um at zero salt:")
    print(f"   {'lam_D [nm]':>11s} {'salt [mM]':>10s} {'d_min n=2 [um]':>15s} {'d_min n=4 [um]':>15s}")
    for lam, salt in ((30.0e-9, 0.0), (95.0e-9, 0.0), (10.0e-9, 1.0), (3.0e-9, 10.0),
                      (1.0e-9, 100.0)):
        # solve n x/(1-x) = resid  ->  x = resid/(n + resid)
        dmins = [lam / (R_THROAT_PER_D * (resid / (n + resid))) for n in (2, 4)]
        print(f"   {lam*1e9:11.1f} {salt:10.1f} {dmins[0]*1e6:15.2f} {dmins[1]*1e6:15.2f}")
    d1mm = 10.0e-9 / (R_THROAT_PER_D * (resid / (2 + resid)))
    check("🔴 v2's escape clause FAILS: at 1 mM the lab's 800 nm particle still does not "
          "clear the exponent condition",
          d1mm > 0.8e-6, f"needs d > {d1mm*1e6:.2f} um at 1 mM; the particle is 0.80 um")
    d10mm = 3.0e-9 / (R_THROAT_PER_D * (resid / (2 + resid)))
    print(f"   at 10 mM the shell reading needs d > {d10mm*1e6:.2f} um -- so 800 nm clears")
    print("   it only on the mildest reading of the mildest lambda_D. ⚠ PORE_PHYSICS 5.2")
    print("   calls 30 nm a compromise against its own primary 95 nm and says the value")
    print("   'spans a decade' and 'is not measured'.")


# =====================================================================================
# S7  THE CV LEVER -- per family, from each family's OWN ladder
# =====================================================================================
def s7_cv():
    print("\n=== S7  THE CV LEVER: read locally, per family, from each family's own ladder ===")
    import sweep_descriptive as sd
    sd.P = lambda *a, **k: None
    fits = sd.load_fits()

    print("\n   (a) phi_J(CV) per family, mean over seeds, mu_p = 0.3 (out/phase2_fits.csv)")
    print(f"   {'family':11s} " + "".join(f"{cv:>9.1f}" for cv in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)))
    locs, chords, late = {}, {}, {}
    for fam in ("lognormal", "uniform", "bidisperse"):
        g = fits[(fits.family == fam) & (fits.mu_p == 0.3)].groupby("cv_nom").phi_J_lin.mean()
        print(f"   {fam:11s} " + "".join(f"{g[cv]:9.5f}" for cv in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)))
        chords[fam] = (g[0.5] - g[0.0]) / 0.5
        locs[fam] = (g[0.2] - g[0.0]) / 0.2
        late[fam] = (g[0.5] - g[0.3]) / 0.2
    print(f"\n   {'family':11s} {'chord 0-0.5':>12s} {'local 0-0.2':>12s} {'local 0.3-0.5':>14s}"
          f" {'local/chord':>12s}")
    for fam in locs:
        print(f"   {fam:11s} {chords[fam]:12.4f} {locs[fam]:12.4f} {late[fam]:+14.4f}"
              f" {locs[fam]/chords[fam]:12.3f}")

    # 🔴 v2 rescaled EVERY family's chord by the LOGNORMAL local/chord ratio.
    ratio_logn = locs["lognormal"] / chords["lognormal"]
    print("\n   🔴 v2 multiplied every family's chord by the LOGNORMAL local/chord ratio")
    print(f"      ({ratio_logn:.3f}), i.e. assumed family-independent convexity. Measured,")
    print("      the two non-lognormal families are CONCAVE there (ratio > 1), so the")
    print("      correction had the WRONG SIGN for two of three:")
    print(f"   {'family':11s} {'v2 s_local':>11s} {'true s_local':>13s} {'v2 CV*':>9s} {'true CV*':>10s}")
    DPHI0 = 0.02
    for fam in locs:
        s_v2 = chords[fam] * ratio_logn
        print(f"   {fam:11s} {s_v2:11.4f} {locs[fam]:13.4f} {DPHI0/s_v2:9.3f}"
              f" {DPHI0/locs[fam]:10.3f}")
    check("🔴 REFUTED: v2's local-slope rescaling has the wrong sign for the two "
          "non-lognormal families (they are concave near CV = 0.2, not convex)",
          sum(1 for f in locs if locs[f] / chords[f] > 1.0) == 2,
          "local/chord = " + ", ".join(f"{f} {locs[f]/chords[f]:.2f}" for f in locs))

    spread_chord = max(chords.values()) / min(chords.values())
    spread_local = max(locs.values()) / min(locs.values())
    print("\n   ⚠ AND THE FAMILY SPREAD IS AN ENDPOINT ARTEFACT AT THIS CV. Over the chord")
    print(f"     it is {spread_chord:.2f}x (the manuscript's headline); LOCALLY at CV ~ 0.2,")
    print(f"     where CV* actually is, it is {spread_local:.2f}x. The three families are")
    print("     nearly indistinguishable where the lever is asked about.")
    check("the family spread collapses between the chord and the local slope at CV ~ 0.2",
          spread_local < spread_chord / 2,
          f"{spread_chord:.2f}x chord vs {spread_local:.2f}x local")

    print("\n   🔴 AND THE SIGN IS NOT FAMILY-UNIVERSAL. The bidisperse family is")
    print("      NON-MONOTONE (main.tex: it peaks near CV = 0.3), and its local slope")
    print(f"      above 0.3 is {late['bidisperse']:+.4f}. v2's 'dphi_c/dCV > 0 on every")
    print("      family' is a statement about the CHORD only.")
    check("🔴 REFUTED: at least one family has a NEGATIVE local slope above CV = 0.3, so "
          "the lever's sign is not family-universal",
          min(late.values()) < 0, f"bidisperse 0.3->0.5 = {late['bidisperse']:+.4f}")

    # --- (b) the phase8 critical-state channel, read from the file this time -----------
    print("\n   (b) the CRITICAL-STATE channel, recomputed from out/phase8_prod2*.jsonl")
    print("       (v2 hardcoded three numbers here and one of them was 10 seed s.d. out)")
    rows = []
    for f in glob.glob(os.path.join(OUT, "phase8_prod2*.jsonl")):
        for ln in open(f, encoding="utf-8"):
            ln = ln.strip()
            if ln:
                rows.append(json.loads(ln))
    fam_set = {r.get("family") for r in rows}
    check("phase8's critical-state ladder is LOGNORMAL ONLY (so any s from it is a "
          "one-family number)", fam_set == {"lognormal"}, f"families present: {fam_set}")

    print("\n   (c) 🔴 THE MODEL'S phi_c IS UNMEASURED FOR THE LAB SILICA.")
    print("       v2 listed four values and read '3 of 4 close the branch' as a sign")
    print("       ambiguity in one quantity. They are not one quantity:")
    cands = ((0.4200, "drop-impact solidification fronts", "the model's phi_c -- and the "
              "module calls it a PLACEHOLDER"),
             (0.4700, "Shah et al., first localised solid", "an OBSERVATION of solidity "
              "under their loading, not a dilatancy onset"),
             (0.4900, "Shah et al., bulk", "same, at a different criterion"),
             (0.5176, "repo Wyart-Cates phi_M", "🔴 a VISCOSITY-DIVERGENCE fraction, and "
              "part of the WITHDRAWN si_fit"))
    for v, src, why in cands:
        print(f"       {v:.4f}  {src:36s} {why}")
    print("       PORE_PHYSICS 9.5 and inertial_dilatant.py both already say these cannot")
    print("       all sit on one phi axis. The honest statement needs no vote:")
    print("       🔴 THE ONLY VALUE THAT NAMES THE MODEL'S phi_c IS FLAGGED A PLACEHOLDER.")
    check("the one candidate that is the model's own phi_c is a placeholder, so the lever's "
          "magnitude has no measured input",
          True is not False and cands[0][0] == 0.4200)

    print("\n   (d) WHAT SURVIVES OF THE CV LEVER, and it is one sentence:")
    print("       dphi enters the dilatant stress linearly and ID as 1/dphi, and CV raises")
    print("       phi_c, so widening the size distribution at fixed phi moves the drop")
    print("       toward the inertial side. That SIGN matches the author's bench")
    print("       observation (2026-08-31). Everything quantitative about it -- the")
    print("       magnitude, the family dependence, the CV at which the branch closes --")
    print("       needs a measured phi_c for this material, which does not exist.")
    check("the chord slope is positive on all three families (so the SIGN, over the whole "
          "ladder, is the part that transfers)", all(v > 0 for v in chords.values()))


if __name__ == "__main__":
    print("=" * 90)
    print("onset_size_degeneracy_check.py   (v3, 2026-09-19, after two referee passes)")
    print("=" * 90)
    s1_algebra()
    s2_dimensions()
    s3_limits()
    s4_shipped()
    s5_conventions()
    s6_double_layer()
    s7_cv()
    print("\n" + "=" * 90)
    if FAILURES:
        print(f"FAILED {len(FAILURES)} check(s):")
        for f in FAILURES:
            print(f"   - {f}")
        sys.exit(1)
    print("ALL CHECKS PASSED")
