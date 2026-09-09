"""size_threshold.py — is there a critical PARTICLE SIZE for solid-like -> inertial?

THE HYPOTHESIS
--------------
The author's, 2026-09-08: *maybe there is simply a size threshold `a*`, above which an
impacting dense suspension stops behaving solid-like and starts behaving as the
inertial-dilatant (re-entrant liquid) mode.*

A size threshold is what a critical value of a dimensionless group IS, so the hypothesis
is not one hypothesis: it is a family, one member per group whose value depends on `a`.
This module enumerates every such group in the lab's problem, computes the `a*` each one
implies at the lab's own operating point, and — the part that decides anything — reports
**how each `a*` moves with the knobs the bench can turn.** Two thresholds that sit at the
same size are still distinguishable if one moves with salt and the other with velocity.

Nothing here is fitted. Every formula is either an identity, a group already in the
repo's own modules, or a textbook definition named at the point of use, and the
convention is stated wherever one exists (a Péclet number and a Stokes number both have
several).

    python size_threshold.py            # the table
    python size_threshold.py --sweep    # how each a* moves with U0, salt, phi, D0

WHAT THIS CAN AND CANNOT SETTLE
-------------------------------
It CANNOT tell you which threshold is the real one. It CAN tell you (i) where each
candidate puts `a*`, (ii) which candidates the lab's two materials already straddle, and
(iii) which pair of knobs separates any two candidates — which is the experiment.
"""

from __future__ import annotations

import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "impact_stress")))

from inertial_dilatant import ID, DilatantSystem  # noqa: E402

# --- the lab's operating point, each value with its origin ---------------------------
PHI = 0.44             # Paper 1 working solid fraction
PHI_C = 0.42           # PLACEHOLDER (impact_stress.lab_silica_dilatant default)
U0 = 3.0               # m/s, bench tower
D0 = 2.5e-3            # m, drop diameter
MU_W = 1.0e-3          # Pa s, water
RHO_P = 2000.0         # kg/m^3, the silica-in-water density model in impact_stress
GAMMA_LV = 0.072       # N/m, suspension-air surface tension
K_B = 1.380649e-23     # J/K, exact (SI 2019)
T_LAB = 293.15         # K
LAMBDA_D = 30.0e-9     # m, Debye length at NOMINAL zero added salt (PORE_PHYSICS 6.3);
                       # the salt series spans 30 nm (0 mM) to 1 nm (100 mM)
R_THROAT_COEF = 0.155 / 2.0   # r_throat = 0.155 * R_particle = 0.0775 * a
                              # (impact_stress.inertial_dilatant.r_throat, corrected
                              # there after review; using `a` understates it 12.9x)

# The ID value at the regime boundary. ⚠ THE KNOWN LANDMINE: [M26] draw the boundary at
# "ID ~ 1", but evaluating their own Eq. 6 at their own Regime-2/3 transition returns
# 6.6e-3 (Status card section 8). Both are carried and both are printed, because a* goes
# as sqrt(ID_boundary) and the two differ by 12.3x in a*.
ID_BOUNDARY = (1.0, 6.6e-3)

# The lab's two materials, for "which side is each on".
MATERIALS = (("silica (lab)", 0.8e-6), ("cornstarch [M26]", 10.0e-6))


def rho_susp(phi: float = PHI) -> float:
    return RHO_P * phi + 1000.0 * (1 - phi)


# =====================================================================================
# the candidate thresholds
# =====================================================================================
def a_drainage_macroscopic(id_boundary: float, u0: float = U0, d0: float = D0,
                           phi: float = PHI, phi_c: float = PHI_C) -> float:
    """T1  ID(a*) = id_boundary, with the drainage length taken as the DROP.

    ID goes as a^2 (through Carman-Kozeny), so a* = a_ref * sqrt(id_boundary/ID(a_ref)).
    This is [M26]'s own reading and it is the one that produces the famous a^2 scaling.
    """
    a_ref = 1.0e-6
    s = DilatantSystem(name="ref", phi=phi, phi_c=phi_c, a=a_ref, rho=rho_susp(phi))
    return a_ref * math.sqrt(id_boundary / ID(s, u0, d0))


def a_drainage_grainscale() -> float:
    """T2  the same criterion with the drainage length PROPORTIONAL TO the grain.

    Then kappa/ell^2 carries no particle size and ID is size-independent, so there is NO
    threshold at all. Returned as NaN, which is the honest value: this candidate's
    prediction is "no transition on the size axis", and that is a falsifiable statement,
    not a missing number.
    """
    return float("nan")


def a_dilatant_strength(lam_d: float = LAMBDA_D, phi: float = PHI,
                        phi_c: float = PHI_C) -> float:
    """T3  frustrated dilation just able to double the confining stress.

    PORE_PHYSICS.md Eq. (9.1): p_dil_max / P_conf = (a / 6 lambda_D) * (dphi/phi), with
    P_conf cancelling, so the criterion is ANCHOR-FREE -- no sigma*, no K_p, no
    permeability. Setting it to 1,

        a* = 6 lambda_D phi / (phi - phi_c)

    Below a*, frustrated dilation cannot raise the particle pressure above what already
    drives it, so it is not available as a strengthening mechanism at all.
    """
    return 6.0 * lam_d * phi / (phi - phi_c)


def a_brownian(u0: float = U0, d0: float = D0, mu: float = MU_W,
               t: float = T_LAB) -> float:
    """T4  Peclet = 1: the colloidal-to-granular divide.

    CONVENTION, and it matters: Pe = 6 pi mu a^3 gammadot / (k_B T), the rotary/shear
    Peclet number for a sphere written with the particle RADIUS as the length. Other
    conventions differ by O(1)-O(10) and therefore by up to ~2x in a*; this one is
    stated so the number can be recomputed under another.

        a* = ( k_B T / (6 pi mu gammadot) )^(1/3) * 2      [diameter, not radius]
    """
    gdot = u0 / d0
    r_star = (K_B * t / (6.0 * math.pi * mu * gdot)) ** (1.0 / 3.0)
    return 2.0 * r_star


def a_stokes(u0: float = U0, d0: float = D0, mu: float = MU_W,
             rho_p: float = RHO_P) -> float:
    """T6  St = 1: grains stop following the fluid.

    St = rho_p a^2 U0 / (9 mu D0)  (the convention used in impact_stress and in the
    derivation's section-1 table), so a* = sqrt(9 mu D0 / (rho_p U0)).
    """
    return math.sqrt(9.0 * mu * d0 / (rho_p * u0))


def a_capillary(u0: float = U0, gamma: float = GAMMA_LV, phi: float = PHI) -> float:
    """T8  the impact stress able to push air through a pore throat -> granulation.

    2 gamma / r_throat = rho U0^2 with r_throat = 0.0775 a, so
        a* = 2 gamma / (0.0775 rho U0^2).
    Above a*, the interface can no longer hold the pack together at this speed.
    """
    return 2.0 * gamma / (R_THROAT_COEF * rho_susp(phi) * u0**2)


def a_grain_inertia(d0: float = D0, phi: float = PHI,
                    rho_p: float = RHO_P) -> float:
    """T7  inertial number I = 1 at the impact stress scale.

    I = gammadot a sqrt(rho_p / P) with gammadot = U0/D0 and P ~ rho U0^2 gives
    I = (a/D0) sqrt(rho_p/rho), independent of U0, so a* = D0 sqrt(rho/rho_p).
    """
    return d0 * math.sqrt(rho_susp(phi) / rho_p)


# =====================================================================================
# how each threshold MOVES -- this is the part that discriminates
# =====================================================================================
SCALINGS = {
    "T1 drainage (ell = D0)":      dict(U0="U0^-1/2", salt="—", dphi="dphi^+1/2",
                                        D0="D0^+1/2", CV="through d32 only"),
    "T2 drainage (ell ~ a)":       dict(U0="—", salt="—", dphi="—", D0="—",
                                        CV="—  (no threshold at all)"),
    "T3 dilatant strength":        dict(U0="NONE", salt="~ lambda_D", dphi="dphi^-1",
                                        D0="NONE", CV="through dphi(CV)"),
    "T4 Brownian (Pe = 1)":        dict(U0="U0^-1/3", salt="—", dphi="—",
                                        D0="D0^+1/3", CV="—"),
    "T6 Stokes (St = 1)":          dict(U0="U0^-1/2", salt="—", dphi="—",
                                        D0="D0^+1/2", CV="—"),
    "T8 capillary / granulation":  dict(U0="U0^-2", salt="—", dphi="—", D0="NONE",
                                        CV="—"),
    "T7 grain inertia (I = 1)":    dict(U0="NONE", salt="—", dphi="—", D0="~ D0",
                                        CV="—"),
}


# Which side of each threshold is the INERTIAL / liquid-like one. This is not a detail:
# three of the seven candidates have the OPPOSITE polarity to the hypothesis and are
# therefore excluded by direction alone, whatever their numbers say.
#   "big"   = a > a* is inertial   (the hypothesis's direction)
#   "small" = a > a* is MORE solid-like  (the reverse)
POLARITY = {
    "T1 drainage (ell = D0)":     ("big",
                                   "ID ~ a^2 and ID >> 1 is the re-entrant liquid"),
    "T2 drainage (ell ~ a)":      ("none", "ID carries no particle size"),
    "T3 dilatant strength":       ("small",
                                   "above a* frustrated dilation CAN strengthen, so "
                                   "big grains are MORE solid-like -- the reverse"),
    "T4 Brownian (Pe = 1)":       ("big",
                                   "above a* the pack is non-Brownian/granular"),
    "T6 Stokes (St = 1)":         ("big", "above a* grains decouple from the fluid"),
    "T8 capillary / granulation": ("big", "above a* the interface cannot hold the pack"),
    "T7 grain inertia (I = 1)":   ("big", "above a* grain inertia matters"),
}


def table() -> list[dict]:
    rows = [
        {"name": f"T1 drainage (ell = D0), ID_b = {ID_BOUNDARY[0]:g}",
         "a_star": a_drainage_macroscopic(ID_BOUNDARY[0]), "key": "T1 drainage (ell = D0)"},
        {"name": f"T1 drainage (ell = D0), ID_b = {ID_BOUNDARY[1]:g}",
         "a_star": a_drainage_macroscopic(ID_BOUNDARY[1]), "key": "T1 drainage (ell = D0)"},
        {"name": "T2 drainage (ell ~ a)", "a_star": a_drainage_grainscale(),
         "key": "T2 drainage (ell ~ a)"},
        {"name": "T3 dilatant strength, 0 mM (lambda_D 30 nm)",
         "a_star": a_dilatant_strength(30e-9), "key": "T3 dilatant strength"},
        {"name": "T3 dilatant strength, 1 mM (lambda_D 10 nm)",
         "a_star": a_dilatant_strength(10e-9), "key": "T3 dilatant strength"},
        {"name": "T3 dilatant strength, 100 mM (lambda_D 1 nm)",
         "a_star": a_dilatant_strength(1e-9), "key": "T3 dilatant strength"},
        {"name": "T4 Brownian (Pe = 1)", "a_star": a_brownian(),
         "key": "T4 Brownian (Pe = 1)"},
        {"name": "T6 Stokes (St = 1)", "a_star": a_stokes(), "key": "T6 Stokes (St = 1)"},
        {"name": "T8 capillary / granulation", "a_star": a_capillary(),
         "key": "T8 capillary / granulation"},
        {"name": "T7 grain inertia (I = 1)", "a_star": a_grain_inertia(),
         "key": "T7 grain inertia (I = 1)"},
    ]
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", action="store_true",
                    help="also print how each a* moves when a knob is turned")
    args = ap.parse_args()

    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    p("=== CANDIDATE SIZE THRESHOLDS at the lab's operating point ===")
    p(f"  phi = {PHI}, phi_c = {PHI_C} (PLACEHOLDER), U0 = {U0} m/s, "
      f"D0 = {D0 * 1e3:.1f} mm, T = {T_LAB} K")
    p(f"  {'threshold':46s} {'a* [um]':>10s}   " + "  ".join(
        f"{n:>18s}" for n, _ in MATERIALS))
    for r in table():
        a = r["a_star"]
        if not math.isfinite(a):
            side = "  ".join(f"{'no threshold':>18s}" for _ in MATERIALS)
            p(f"  {r['name']:46s} {'--':>10s}   {side}")
            continue
        side = "  ".join(
            f"{('ABOVE' if am > a else 'below') + f' ({am / a:.2f}x)':>18s}"
            for _, am in MATERIALS)
        p(f"  {r['name']:46s} {a * 1e6:10.3f}   {side}")

    p()
    p("=== WHICH CANDIDATES SURVIVE: polarity AND separation ===")
    p("  A candidate must do TWO things: put silica below and cornstarch above (so it")
    p("  separates the observation), AND have the right POLARITY (a > a* must be the")
    p("  INERTIAL side, not the solid-like one).")
    p()
    p(f"  {'threshold':30s} {'a*[um]':>9s} {'separates?':>11s} {'polarity':>9s} "
      f"{'VERDICT':>10s}")
    survivors = []
    a_si, a_cs = MATERIALS[0][1], MATERIALS[1][1]
    for key, (pol, why) in POLARITY.items():
        cands = [r for r in table() if r["key"] == key and math.isfinite(r["a_star"])]
        if not cands:
            p(f"  {key:30s} {'--':>9s} {'no':>11s} {pol:>9s} {'EXCLUDED':>10s}")
            p(f"      ({why})")
            continue
        seps = [r for r in cands if a_si < r["a_star"] < a_cs]
        a_txt = ("/".join(f"{r['a_star'] * 1e6:.2f}" for r in cands))[:9]
        sep = "YES" if seps else "no"
        ok = bool(seps) and pol == "big"
        p(f"  {key:30s} {a_txt:>9s} {sep:>11s} {pol:>9s} "
          f"{('SURVIVES' if ok else 'EXCLUDED'):>10s}")
        if not ok:
            p(f"      ({why})" if pol != "big"
              else "      (does not separate: both materials on the same side)")
        if ok:
            survivors.extend(seps)
    p()
    p(f"  🔑 {len(survivors)} of {len(POLARITY)} candidates survive both filters:")
    for r in survivors:
        p(f"     - {r['name']}  ->  a* = {r['a_star'] * 1e6:.3f} um")

    # --- the synthesis: a WINDOW, not a half-line -----------------------------------
    p()
    p("=== THE HYPOTHESIS SHARPENED: a WINDOW in size, not a threshold ===")
    p("  T1 and T3 both grow with particle size and they are NOT redundant, because")
    p("  they grow at different rates and point opposite ways:")
    p("     T3  strengthening capacity   ~ a^1   below it dilatancy CANNOT strengthen")
    p("     T1  drainage                 ~ a^2   above it the pore fluid keeps up")
    p("  so frustrated dilation can hold a drop solid only in the band between them:")
    p()
    p(f"  {'salt':>8s} {'U0[m/s]':>8s} {'lower edge T3':>14s} {'upper edge T1':>14s} "
      f"{'window':>22s}")
    for lam_nm, salt in ((30.0, "0 mM"), (10.0, "1 mM"), (1.0, "100 mM")):
        for u in (0.5, 3.0):
            lo_a = a_dilatant_strength(lam_nm * 1e-9)
            hi_a = a_drainage_macroscopic(6.6e-3, u0=u)
            inside = [n for n, am in MATERIALS if lo_a < am < hi_a]
            txt = ", ".join(inside) if inside else "(neither material)"
            p(f"  {salt:>8s} {u:8.2f} {lo_a * 1e6:14.3f} {hi_a * 1e6:14.3f} "
              f"{txt:>22s}")
    p()
    p("  🔑 THREE things fall out, and each is checkable on an axis already in use.")
    p("  (a) CORNSTARCH sits INSIDE the window at 0.5 m/s and EXITS THE TOP by 3 m/s.")
    p("      That IS [M26]'s re-entrant transition, recovered as a window edge moving")
    p("      down past a fixed particle size as U0 rises.")
    p("  (b) 800 nm SILICA sits BELOW the window at 0 and 1 mM -- too small for")
    p("      frustrated dilation to raise the particle pressure above what already")
    p("      drives it. So silica's solid-like branch, which is OBSERVED, is NOT")
    p("      frustrated dilation; it has to be frictional/repulsive jamming. ⚠ That")
    p("      is a substantive claim about the lab's own material and it says silica")
    p("      and cornstarch are solid for DIFFERENT REASONS.")
    p("  (c) Adding salt lowers the bottom edge past 800 nm, so the SAME silica should")
    p("      enter the window inside the salt series Paper 1 already runs -- no new")
    p("      particles, no new apparatus.")

    # What the hypothesis, if true, says about the ID boundary calibration.
    p()
    p("=== WHAT THE HYPOTHESIS IMPLIES ABOUT THE ID BOUNDARY (the card's landmine) ===")
    a_at_1 = a_drainage_macroscopic(1.0)
    lo = (a_si / a_at_1) ** 2
    hi = (a_cs / a_at_1) ** 2
    p(f"  a* goes as sqrt(ID_boundary), and a*(ID_b = 1) = {a_at_1 * 1e6:.2f} um.")
    p("  For the threshold to sit BETWEEN 800 nm silica and 10 um cornstarch,")
    p(f"      ID_boundary must lie in [{lo:.3g}, {hi:.3g}]")
    p("  which BRACKETS the 6.6e-3 that [M26]'s own Eq. 6 returns at their own")
    p(f"  Regime-2/3 transition, and EXCLUDES ID_boundary = 1 by "
      f"{1.0 / hi:.0f}x.")
    p("  ⚠ This is an inference FROM the hypothesis, not independent support for it:")
    p("  it says the hypothesis is only self-consistent on the smaller calibration.")

    # What CV would have to be for a pure size threshold to explain the CV observation.
    p()
    p("=== CAN A PURE SIZE THRESHOLD EXPLAIN THE CV OBSERVATION? ===")
    p("  On a pure size threshold, CV acts ONLY by moving the moment the threshold")
    p("  uses. T1 uses Carman-Kozeny, whose length is the SAUTER mean, and for a")
    p("  lognormal at fixed number mean d32/d10 = (1 + CV^2)^2. So crossing a* from a")
    p("  starting d10 needs a specific CV -- which is a testable number, not a story.")
    p(f"  {'d10 [um]':>9s} {'a* [um]':>9s} {'d32/d10 needed':>15s} {'CV needed':>10s}")
    for a_star_um, tag in ((a_drainage_macroscopic(6.6e-3) * 1e6, "T1, ID_b=6.6e-3"),):
        for d10_um in (0.8, 1.5, 3.0):
            need = a_star_um / d10_um
            cv_need = math.sqrt(math.sqrt(need) - 1.0) if need >= 1 else 0.0
            p(f"  {d10_um:9.2f} {a_star_um:9.2f} {need:15.2f} {cv_need:10.2f}"
              f"   ({tag})")
    p()
    p("  🔑 At the lab's 800 nm number-mean a pure size threshold needs CV ~ 1.5 to")
    p("     reach a*, and CV ~ 1.1 even at a 1.5 um mean. Those are far broader than")
    p("     the CV 0.2-0.5 this project has simulated. SO: if the author's two samples")
    p("     have CV below ~1, a pure size threshold acting through d32 CANNOT carry")
    p("     them across T1, and something else must be moving as well. If they are")
    p("     broader than that, it can -- which makes the raw size distributions the")
    p("     measurement that decides it (Q1 on the Status card, still open).")

    p()
    p("=== HOW EACH THRESHOLD MOVES -- this is what discriminates them ===")
    p(f"  {'threshold':30s} {'U0':>10s} {'salt':>14s} {'dphi':>12s} {'D0':>10s}")
    for k, s in SCALINGS.items():
        p(f"  {k:30s} {s['U0']:>10s} {s['salt']:>14s} {s['dphi']:>12s} "
          f"{s['D0']:>10s}")
    p()
    p("  🔑 T3 is the ONLY candidate that moves with SALT and the ONLY one that does")
    p("     NOT move with U0. Every other finite threshold falls with U0 and is")
    p("     salt-blind. So a size ladder crossed with a VELOCITY ladder and a SALT")
    p("     ladder identifies which threshold is operating, with no new apparatus.")

    if args.sweep:
        p()
        p("=== the sweeps, explicitly ===")
        p(f"  {'U0 [m/s]':>10s} {'T1(ID_b=1)':>12s} {'T1(6.6e-3)':>12s} "
          f"{'T3 (0 mM)':>11s} {'T4 Brown':>10s} {'T6 Stokes':>10s} {'T8 cap':>10s}")
        for u in (0.5, 1.0, 2.0, 3.0, 5.0):
            p(f"  {u:10.2f} {a_drainage_macroscopic(1.0, u0=u) * 1e6:12.3f} "
              f"{a_drainage_macroscopic(6.6e-3, u0=u) * 1e6:12.3f} "
              f"{a_dilatant_strength() * 1e6:11.3f} "
              f"{a_brownian(u0=u) * 1e6:10.3f} {a_stokes(u0=u) * 1e6:10.3f} "
              f"{a_capillary(u0=u) * 1e6:10.3f}")
        p()
        p(f"  {'lambda_D[nm]':>12s} {'salt':>10s} {'T3 [um]':>10s}   "
          f"is 800 nm silica above it?")
        for lam_nm, salt in ((30.0, "0 mM"), (10.0, "1 mM"), (3.0, "10 mM"),
                             (1.0, "100 mM")):
            a3 = a_dilatant_strength(lam_nm * 1e-9)
            p(f"  {lam_nm:12.1f} {salt:>10s} {a3 * 1e6:10.3f}   "
              f"{'YES -- dilatant strengthening AVAILABLE' if 0.8e-6 > a3 else 'no'}")
        p()
        p("  🔑 T3 crosses the lab's own particle size INSIDE its own salt series.")
        p("     At nominal zero salt the 800 nm silica is BELOW the threshold and")
        p("     frustrated dilation is not available; by 10 mM it is ABOVE it. So the")
        p("     salt axis Paper 1 already sweeps should switch the solid-like branch")
        p("     on and off AT FIXED PARTICLE SIZE -- which is a size-threshold")
        p("     prediction testable without making a single new particle.")
        p()
        p(f"  {'phi':>6s} {'dphi':>8s} {'T3 [um]':>10s} {'T1(6.6e-3) [um]':>16s}")
        for phi in (0.44, 0.46, 0.48, 0.50):
            p(f"  {phi:6.2f} {phi - PHI_C:8.3f} "
              f"{a_dilatant_strength(phi=phi) * 1e6:10.3f} "
              f"{a_drainage_macroscopic(6.6e-3, phi=phi) * 1e6:16.3f}")

    out = os.path.join(HERE, "out", "size_threshold.txt")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
