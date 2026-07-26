"""
Quantitative validation of the 6-DOF integrator against Hummel (2003).

This is the referee-facing validation for P1. It compares our integrator to two
INDEPENDENTLY PUBLISHED, fully-specified results in the same source that
supplied the aerodynamic coefficients, so agreement tests the integrator and the
frame conventions rather than the coefficient fit.

CASE A -- "common throw", Hummel Sec. 2.4.3 and Fig. 2-9 (pp. 24-26)
    Release conditions are given exactly, both in the text and in the thesis
    Appendix A comment block:
        theta0 = 0.192 rad (11.0 deg)   pitch angle of the disc
        beta0  = 0.105 rad ( 6.0 deg)   flight-path angle above horizontal
        alpha0 = theta0 - beta0 = 5.0 deg
        speed  = 14 m/s,  spin = 50 rad/s,  roll = 0,  release height 1 m
    Published outcomes, read from the Fig. 2-9 caption and surrounding text:
        flight distance       ~25   m
        maximum height         2.55 m
        time to ground         3.5  s
        terminal speed        < 4   m/s ("nearly 1/3 its original speed")
        lateral behaviour      curves RIGHT first, then back LEFT after ~3 s
        alpha                  starts near 6 deg, stays < 9 deg for ~1.25 s
        pitch angle theta      varies by less than 2 deg over the whole flight

CASE B -- steady downward glide, Hummel Sec. 2.4.2 (p. 23-24)
    An analytic trim condition, not a simulation: at the angle of attack where
    the pitching moment vanishes the disc neither precesses nor pitches, and
    lift plus drag balance weight along a straight descending path.
        alpha_trim = 9 deg      (CM = 0)
        glide angle gamma = -19.3 deg
        CL = 0.67, CD ~ 0.22-0.24
        speed = 9.1 m/s
    We check the trim angle and the resulting glide angle against a closed-form
    solution, which validates the coefficient evaluation independently of the
    time integration.

Run:  python validation_hummel.py
"""

from __future__ import annotations

import numpy as np

import disc_6dof as d6
from coefficients import (
    DEG, HUMMEL_LONG, HUMMEL_SHORT, POTTS_CROWTHER, TabulatedCoefficients,
    ULTRASTAR,
)

# --- Published reference values (Hummel 2003). Do not edit without a citation.
REF_RANGE_M = 25.0
REF_APEX_M = 2.55
REF_TIME_S = 3.5
# The thesis gives two slightly inconsistent statements for the terminal speed:
# "less than 4 m/sec" and "nearly 1/3 its original speed" (14/3 = 4.67 m/s).
# We test against the ratio form, which is the quantitative one, with the same
# tolerance as the other metrics, and report the absolute value alongside.
REF_TERMINAL_SPEED_MAX = 4.0
REF_TERMINAL_SPEED_RATIO = 1.0 / 3.0
REF_ALPHA_TRIM_DEG = 9.0
REF_GLIDE_ANGLE_DEG = -19.3
REF_TRIM_CL = 0.67

# Tolerances. These are deliberately loose: the reference numbers were read off
# a printed figure, and Hummel's own coefficient sets disagree with each other
# by more than this. A 15% band on range/apex/time is the standard "does the
# trajectory model reproduce the published flight" criterion.
TOL_FRAC = 0.15


def hummel_common_throw(coeffs, geom=ULTRASTAR):
    """Integrate the Hummel Sec. 2.4.3 'common release conditions' flight."""
    speed = 14.0
    beta_deg = 6.0            # flight-path angle above horizontal
    theta_deg = 11.0          # disc pitch angle
    alpha_deg = theta_deg - beta_deg     # = 5.0 by construction
    spin = 50.0
    launch = np.array([0.0, 0.0, 1.0])

    r0, v0, q0, w0 = d6.build_initial_conditions(
        launch_pos=launch,
        speed=speed,
        heading_deg=0.0,
        elevation_deg=beta_deg,
        release_alpha_deg=alpha_deg,
        roll_deg=0.0,
        spin_rad_s=spin,
        handedness="right_backhand",
    )
    res = d6.simulate(r0, v0, q0, w0, coeffs=coeffs, geom=geom,
                      wind=d6.still_air, t_max=10.0)
    return res


def summarize(res, coeffs):
    """Extract the quantities Hummel reports, plus an alpha-range audit."""
    speeds = np.linalg.norm(res.v, axis=1)
    alpha_deg = np.array([d.get("alpha_deg", np.nan) for d in res.diag])
    lo, hi = coeffs.alpha_valid_deg
    frac_out = float(np.mean((alpha_deg < lo) | (alpha_deg > hi)))

    # Lateral excursion: +y is LEFT of the throw for a heading along +x, so a
    # right-hand backhand that "curves right" first shows y < 0 early on.
    y = res.r[:, 1]
    t = res.t
    early = t < 0.6 * t[-1]
    # Hummel Fig. 2-9a: the disc curves right throughout the first ~3 s, then
    # reverses and comes back left. In our z-up, +y-left frame "right" is
    # y decreasing, so the signature is a MINIMUM of y before touchdown.
    i_extreme = int(np.argmin(y))
    reverses = i_extreme < len(y) - 1
    return dict(
        range_m=res.range_m,
        apex_m=float(res.r[:, 2].max()),
        time_s=res.flight_time,
        terminal_speed=float(speeds[-1]),
        initial_speed=float(speeds[0]),
        y_min=float(y.min()),
        y_max=float(y.max()),
        y_final=float(y[-1]),
        y_early_extreme=float(y[early].min() if early.any() else np.nan),
        turn_reverses=bool(reverses),
        t_reversal=float(t[i_extreme]),
        lateral_excursion_m=float(-y[i_extreme]),
        alpha_min=float(np.nanmin(alpha_deg)),
        alpha_max=float(np.nanmax(alpha_deg)),
        alpha_frac_out_of_range=frac_out,
    )


def trim_glide(coeffs, geom=ULTRASTAR):
    """Closed-form steady glide at the zero-pitching-moment angle of attack.

    Solves CM(alpha, q=0) = 0 for alpha, then the straight-glide force balance
        mg cos(gamma) = 0.5 rho v^2 A CL
        mg sin(gamma) = 0.5 rho v^2 A CD
    which gives gamma = -atan(CD/CL) and v from the resultant.
    """
    # CM is monotone in alpha over the physical range for every set here, so a
    # bisection on [-10 deg, 30 deg] is safe and needs no derivative.
    from scipy.optimize import brentq

    f = lambda a: coeffs.CM(a, 0.0)
    a_lo, a_hi = -10.0 * DEG, 30.0 * DEG
    if f(a_lo) * f(a_hi) > 0:
        return None
    alpha_trim = brentq(f, a_lo, a_hi, xtol=1e-10)

    CL = coeffs.CL(alpha_trim)
    CD = coeffs.CD(alpha_trim)
    gamma = -np.arctan2(CD, CL)
    weight = geom.mass * geom.g
    resultant = np.hypot(CL, CD)
    v = np.sqrt(2.0 * weight / (geom.rho_air * geom.area * resultant))
    return dict(alpha_trim_deg=alpha_trim / DEG, CL=CL, CD=CD,
                gamma_deg=gamma / DEG, v_trim=v)


FIG_PATH = "validation_hummel.png"


def _save_figure(sets, geom):
    """Plot the reconstructed flight against Hummel's published markers."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    for label, cf in sets:
        res = hummel_common_throw(cf, geom)
        alpha = np.array([d.get("alpha_deg", np.nan) for d in res.diag])
        style = dict(lw=2.2, zorder=3) if "tabulated" in label else dict(
            lw=1.2, alpha=0.75)
        axes[0].plot(res.r[:, 0], res.r[:, 2], label=label, **style)
        axes[1].plot(res.r[:, 0], res.r[:, 1], label=label, **style)
        axes[2].plot(res.t, alpha, label=label, **style)

    # Published reference values, Hummel (2003) Sec. 2.4.3 / Fig. 2-9.
    axes[0].plot([REF_RANGE_M], [0.0], "k*", ms=15, zorder=5,
                 label="Hummel 2003 (published)")
    axes[0].axhline(REF_APEX_M, color="k", ls=":", lw=1)
    axes[0].annotate(f"published apex {REF_APEX_M} m", (0.4, REF_APEX_M),
                     fontsize=8, va="bottom")
    axes[0].set(xlabel="downrange x [m]", ylabel="height z [m]",
                title="Hummel Fig. 2-9 flight: profile")
    axes[0].legend(fontsize=7)

    axes[1].axhline(0, color="0.6", lw=0.8)
    axes[1].set(xlabel="downrange x [m]", ylabel="lateral y [m]  (+y = left)",
                title="Turn (right) then fade (left)")

    lo, hi = sets[-1][1].alpha_valid_deg
    axes[2].axhspan(lo, hi, color="tab:green", alpha=0.10)
    axes[2].annotate("tabulated data range", (0.05, hi), fontsize=8,
                     va="bottom")
    axes[2].axhline(REF_ALPHA_TRIM_DEG, color="k", ls=":", lw=1)
    axes[2].set(xlabel="time [s]", ylabel=r"angle of attack $\alpha$ [deg]",
                title=r"$\alpha$ excursion vs. validity range")

    for ax in axes:
        ax.grid(alpha=0.3)
    fig.suptitle("P1 6-DOF validation against Hummel (2003) published flight",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG_PATH, dpi=150)
    plt.close(fig)
    print(f"\n  [figure written: {FIG_PATH}]")


def main():
    geom = ULTRASTAR
    print("=" * 78)
    print("P1 VALIDATION AGAINST HUMMEL (2003) -- published reference flights")
    print("=" * 78)
    print(f"\nDisc geometry (Hummel Table 3-1): m={geom.mass} kg, "
          f"d={geom.diameter} m, A={geom.area} m^2")
    print(f"  Is (spin axis) = {geom.Is:.5f} kg m^2, "
          f"Ip (transverse) = {geom.Ip:.5f} kg m^2")
    thin_Is = 0.5 * geom.mass * geom.radius ** 2
    thin_Ip = 0.25 * geom.mass * geom.radius ** 2
    print(f"  (a thin uniform disc would give Is={thin_Is:.5f}, "
          f"Ip={thin_Ip:.5f} -- low by "
          f"{100*(1-thin_Is/geom.Is):.0f}% and "
          f"{100*(1-thin_Ip/geom.Ip):.0f}%)")

    # ---------------- CASE B first: it is analytic and cheap ---------------
    print("\n" + "-" * 78)
    print("CASE B -- steady downward glide (Hummel Sec. 2.4.2)")
    print("-" * 78)
    print(f"  published: alpha_trim = {REF_ALPHA_TRIM_DEG} deg, "
          f"gamma = {REF_GLIDE_ANGLE_DEG} deg, CL = {REF_TRIM_CL}, "
          f"v = 9.1 m/s")
    print(f"  {'coefficient set':<34s} {'a_trim':>8s} {'CL':>7s} {'CD':>7s} "
          f"{'gamma':>8s} {'v':>7s}")
    sets_b = [("hummel_short", HUMMEL_SHORT), ("hummel_long", HUMMEL_LONG),
              ("potts_crowther", POTTS_CROWTHER),
              ("potts_tabulated", TabulatedCoefficients())]
    for label, cf in sets_b:
        tg = trim_glide(cf, geom)
        if tg is None:
            print(f"  {label:<34s}   (no zero-pitching-moment angle in range)")
            continue
        print(f"  {label:<34s} {tg['alpha_trim_deg']:>7.2f}d "
              f"{tg['CL']:>7.3f} {tg['CD']:>7.3f} "
              f"{tg['gamma_deg']:>7.2f}d {tg['v_trim']:>6.2f}")

    # ---------------- CASE A: the full trajectory --------------------------
    print("\n" + "-" * 78)
    print("CASE A -- 'common throw' 25 m flight (Hummel Sec. 2.4.3, Fig. 2-9)")
    print("-" * 78)
    print("  release: v=14 m/s, spin=50 rad/s, theta=11 deg, beta=6 deg, "
          "alpha=5 deg, z0=1 m")
    print(f"  published: range~{REF_RANGE_M} m, apex={REF_APEX_M} m, "
          f"t={REF_TIME_S} s, terminal speed<{REF_TERMINAL_SPEED_MAX} m/s")
    print()
    hdr = (f"  {'coefficient set':<34s} {'range':>7s} {'apex':>6s} "
           f"{'time':>6s} {'v_end':>6s} {'a_min':>6s} {'a_max':>6s} "
           f"{'%out':>5s}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    results = {}
    for label, cf in sets_b:
        res = hummel_common_throw(cf, geom)
        s = summarize(res, cf)
        results[label] = s
        print(f"  {label:<34s} {s['range_m']:>6.2f}m {s['apex_m']:>5.2f}m "
              f"{s['time_s']:>5.2f}s {s['terminal_speed']:>5.2f} "
              f"{s['alpha_min']:>5.1f} {s['alpha_max']:>5.1f} "
              f"{100*s['alpha_frac_out_of_range']:>4.0f}%")

    # ---------------- Pass/fail against the published numbers --------------
    print("\n" + "-" * 78)
    print("AGREEMENT WITH PUBLISHED VALUES (+/- %d%%)" % (100 * TOL_FRAC))
    print("-" * 78)
    checks = [("range_m", REF_RANGE_M, "m"),
              ("apex_m", REF_APEX_M, "m"),
              ("time_s", REF_TIME_S, "s")]
    any_fail = False
    for label, _ in sets_b:
        s = results[label]
        parts = []
        ok_all = True
        for key, ref, unit in checks:
            err = abs(s[key] - ref) / ref
            ok = err <= TOL_FRAC
            ok_all &= ok
            parts.append(f"{key.split('_')[0]}={s[key]:.2f}{unit} "
                         f"({100*err:+.0f}%){'' if ok else ' FAIL'}")
        v_ref = REF_TERMINAL_SPEED_RATIO * s["initial_speed"]
        v_err = abs(s["terminal_speed"] - v_ref) / v_ref
        vend_ok = v_err <= TOL_FRAC
        ok_all &= vend_ok
        parts.append(f"v_end={s['terminal_speed']:.2f} m/s "
                     f"({100*v_err:+.0f}% vs v0/3){'' if vend_ok else ' FAIL'}")
        print(f"  {label:<20s} {'PASS' if ok_all else 'CHECK'}: "
              + "; ".join(parts))
        any_fail |= not ok_all

    print("\n" + "-" * 78)
    print("LATERAL SIGNATURE (Hummel Fig. 2-9a: right first, back left ~3 s)")
    print("-" * 78)
    for label, _ in sets_b:
        s = results[label]
        verdict = (f"turns right {s['lateral_excursion_m']:.2f} m, "
                   f"reverses at t={s['t_reversal']:.2f} s"
                   if s["turn_reverses"] else
                   "turns right monotonically, NO reversal before landing")
        print(f"  {label:<20s} {verdict}")

    _save_figure(sets_b, geom)

    print("\nNotes:")
    print("  * '%out' is the fraction of integration samples whose angle of")
    print("    attack fell outside the coefficient set's stated validity range")
    print("    -- the linear sets are only defensible inside it.")
    print("  * Lateral turn-then-fade is reported by test_disc_6dof.py, which")
    print("    asserts the sign of the early and late lateral excursion.")
    return 1 if any_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
