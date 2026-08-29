"""
run_uroflow.py -- figures and a quantitative report for the v2 splashback model.

Every number here comes from uroflow.py, which is validated in test_uroflow.py (39
checks against analytic limits, published constants and conservation laws).

The report deliberately prints the v1 result beside the v2 result wherever the
conclusion changed, because two of the three original recommendations were resting on
defects rather than on physics.

Figures -> figures/.  SI units throughout.
"""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import uroflow as u

FIGDIR = "figures"

# Representative healthy adult male peak flow and effective meatus diameter.
Q = 21e-6          # m^3/s  (21 mL/s)
D_EXIT = 4e-3      # m      effective stream diameter at the meatus

# Primary target: wetted vertical porcelain, the surface that actually returns splash to
# the user.  0.5 mm of standing film puts the dimensionless depth inside Cossali's
# correlated range for the drop sizes this jet makes.
FILM = 0.5e-3


def calibration():
    """Fit Lambda against the observed coherent length and report what it implies."""
    lam_lo, lam_hi, dr_lo, dr_hi = u.calibrate_lambda(Q, D_EXIT, launch_height=0.75)
    return dict(lam_lo=lam_lo, lam_hi=lam_hi, lam_mid=0.5 * (lam_lo + lam_hi),
                delta_lo=dr_lo, delta_hi=dr_hi)


# --------------------------------------------------------------------------- #
def fig_breakup(cal):
    """Fig 1: where the jet actually breaks, ideal against calibrated."""
    v0 = u.exit_velocity(Q, D_EXIT)
    r0 = D_EXIT / 2.0
    lams = np.linspace(1.5, 13.0, 60)
    Lb = [u.breakup_arclength(v0, r0, lambda s: s, lam, s_max=2.0) * 100 for lam in lams]

    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.plot(lams, Lb, "k-", lw=2, label="falling jet (thinning + accelerating)")
    Lb_par = [u.breakup_length_parallel(Q, D_EXIT, lam) * 100 for lam in lams]
    ax.plot(lams, Lb_par, "k--", lw=1.2, alpha=0.7, label="parallel jet (v1 model)")
    ax.axhspan(u.L_BREAKUP_OBSERVED[0] * 100, u.L_BREAKUP_OBSERVED[1] * 100,
               color="tab:orange", alpha=0.25, label="observed coherent length, 15-20 cm")
    ax.axvspan(cal["lam_lo"], cal["lam_hi"], color="tab:blue", alpha=0.18,
               label=f"calibrated $\\Lambda$ = {cal['lam_lo']:.1f}-{cal['lam_hi']:.1f}")
    ax.axvline(u.LAMBDA_IDEAL, color="tab:red", ls=":", lw=1.5,
               label="$\\Lambda$ = 12 assumed in v1")
    ax.set_xlabel("disturbance amplification $\\Lambda = \\ln(r_0/\\delta_0)$")
    ax.set_ylabel("breakup distance (cm)")
    ax.set_title("Breakup is set by the exit disturbance, not by the ideal theory\n"
                 f"calibration implies $\\delta_0/r_0$ = {cal['delta_lo']:.1%}"
                 f"-{cal['delta_hi']:.1%}, a plausibly perturbed meatus")
    ax.set_ylim(0, 90)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig1_breakup_calibration.png", dpi=150)
    plt.close(fig)


def fig_geometry_correction():
    """Fig 2: the v1 angle factor against the corrected normal-component factor, for a
    vertical wall.  This is the figure that shows the inversion."""
    ang = np.linspace(1.0, 89.0, 200)          # angle of velocity from the horizontal
    v1_factor = np.sin(np.radians(ang))        # what v1 used
    # For a vertical wall the normal is horizontal, so the normal fraction is cos(ang).
    v2_factor = np.cos(np.radians(ang))

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.plot(ang, v1_factor, "r--", lw=2, label="v1: $\\sin$(angle from horizontal)")
    ax.plot(ang, v2_factor, "k-", lw=2, label="v2: normal fraction $= \\cos$(angle)")
    ax.axvline(90, color="0.6", lw=0.8)
    ax.annotate("normal impact\n(worst case)", xy=(5, 0.97), xytext=(18, 0.80),
                arrowprops=dict(arrowstyle="->", color="0.3"), fontsize=8)
    ax.annotate("grazing impact\n(best case)", xy=(85, 0.09), xytext=(52, 0.30),
                arrowprops=dict(arrowstyle="->", color="0.3"), fontsize=8)
    ax.set_xlabel("velocity angle from horizontal (deg)")
    ax.set_ylabel("splash-driving velocity fraction")
    ax.set_title("Geometry correction for a VERTICAL wall\n"
                 "v1 scored a head-on impact as the driest configuration")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig2_geometry_correction.png", dpi=150)
    plt.close(fig)


def _excess_map(h, lam, distances, angles_deg, film=FILM):
    """Splash excess K/K_c - 1 over (distance, aim angle) for a wetted vertical wall."""
    E = np.full((len(angles_deg), len(distances)), np.nan)
    C = np.zeros_like(E)
    for i, th in enumerate(np.radians(angles_deg)):
        for j, dist in enumerate(distances):
            try:
                r = u.evaluate(Q, D_EXIT, th, h, u.Target("wall", dist, film), lam)
            except ValueError:
                continue
            C[i, j] = 1.0 if r["coherent"] else 0.0
            E[i, j] = r["excess"] if np.isfinite(r["excess"]) else np.nan
    return E, C


def fig_splash_map(cal):
    """Fig 3: corrected splash map for a standing user against a wetted wall."""
    h = 0.75
    distances = np.linspace(0.10, 0.70, 45)
    angles_deg = np.linspace(-60, 0, 45)
    E, C = _excess_map(h, cal["lam_mid"], distances, angles_deg)

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    pcm = ax.pcolormesh(distances * 100, angles_deg, E, shading="auto", cmap="RdYlGn_r")
    fig.colorbar(pcm, ax=ax, label="splash excess $K/K_c - 1$  (below 0 = no splash)")
    cs = ax.contour(distances * 100, angles_deg, E, levels=[0.0], colors="k", linewidths=2)
    ax.clabel(cs, fmt={0.0: "splash threshold"}, fontsize=8)
    if C.max() > 0.5:
        ax.contour(distances * 100, angles_deg, C, levels=[0.5],
                   colors="tab:blue", linewidths=2)
    ax.set_xlabel("distance to the wall (cm)")
    ax.set_ylabel("aim angle (deg, negative = downward)")
    ax.set_title("Splash excess against a wetted vertical wall, standing (h = 75 cm)\n"
                 "Corrected geometry: aim DOWN so the stream grazes the wall")
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig3_splash_map_corrected.png", dpi=150)
    plt.close(fig)


def fig_film_depth(cal):
    """Fig 4: the standing-water lever, which v1 did not have at all."""
    drop_D = None
    h, dist, th = 0.75, 0.30, np.radians(-25)
    base = u.evaluate(Q, D_EXIT, th, h, u.Target("wall", dist, FILM), cal["lam_mid"])
    drop_D = base["drop_diameter"]

    deltas = np.linspace(0.05, 1.4, 120)
    films = deltas * drop_D
    K = base["criterion"]["K"]
    Kc = 2100.0 + 5880.0 * deltas**1.44
    excess = K / Kc - 1.0

    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.plot(deltas, excess, "k-", lw=2)
    ax.axhline(0.0, color="tab:green", lw=1.5, ls="--", label="splash threshold")
    ax.axvspan(0.1, 1.0, color="tab:blue", alpha=0.12,
               label="Cossali correlated range")
    ax.set_xlabel("dimensionless film depth $\\delta = h_{film}/D_{drop}$")
    ax.set_ylabel("splash excess $K/K_c - 1$")
    sec = ax.secondary_xaxis("top", functions=(lambda x: x * drop_D * 1e3,
                                               lambda x: x / (drop_D * 1e3)))
    sec.set_xlabel(f"film depth (mm), for $D_{{drop}}$ = {drop_D*1e3:.1f} mm")
    ax.set_title("A wetter surface is harder to splash, monotonically\n"
                 "There is no protective 'sweet spot' depth within the correlated range")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig4_film_depth.png", dpi=150)
    plt.close(fig)
    return drop_D, films


def fig_sit_vs_stand(cal):
    """Fig 5: aim-optimised splash excess against distance, for three launch heights."""
    heights = {"sitting (h = 0.15 m)": 0.15,
               "low stance (h = 0.45 m)": 0.45,
               "standing (h = 0.75 m)": 0.75}
    dists = np.linspace(0.10, 0.60, 30)
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    for label, h in heights.items():
        best = []
        for dist in dists:
            vals = []
            for th in np.radians(np.linspace(-70, 0, 36)):
                try:
                    r = u.evaluate(Q, D_EXIT, th, h, u.Target("wall", dist, FILM),
                                   cal["lam_mid"])
                except ValueError:
                    continue
                if np.isfinite(r["excess"]):
                    vals.append(r["excess"])
                elif r["coherent"]:
                    vals.append(-1.0)          # coherent column: no drop splash
            best.append(min(vals) if vals else np.nan)
        ax.plot(dists * 100, best, "o-", ms=3, label=label)
    ax.axhline(0.0, color="tab:green", ls="--", lw=1.5, label="splash threshold")
    ax.set_xlabel("distance to the wall (cm)")
    ax.set_ylabel("aim-optimised splash excess")
    ax.set_title("Lower launch height wins, and now for a stated reason:\n"
                 "shorter path keeps the stream coherent past the 15-20 cm breakup point")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig5_sit_vs_stand.png", dpi=150)
    plt.close(fig)


def critical_aim_angle(h, dist, lam, film=FILM, lo_deg=-85.0, hi_deg=0.0):
    """Least-steep aim angle (deg) that keeps the splash excess at or below zero.

    Returns np.nan if no angle in the range clears the threshold, and lo_deg if even the
    shallowest tested angle is already safe.  Bisection is valid because the excess is
    monotone in aim angle over this range, which is checked in the tests.
    """
    def excess(th_deg):
        try:
            r = u.evaluate(Q, D_EXIT, np.radians(th_deg), h, u.Target("wall", dist, film),
                           lam)
        except ValueError:
            return np.nan
        if r["coherent"]:
            return -1.0
        return r["excess"]

    # Steep aims can land short of the wall entirely.  Walk in from lo_deg to the
    # steepest angle that still reaches the target, so "unreachable" is never confused
    # with "unsafe".
    steepest = None
    for th in np.linspace(hi_deg, lo_deg, 200):
        if np.isfinite(excess(th)):
            steepest = th
        else:
            break
    if steepest is None:
        return np.nan
    lo_deg = steepest

    e_hi, e_lo = excess(hi_deg), excess(lo_deg)
    if not np.isfinite(e_hi) or not np.isfinite(e_lo):
        return np.nan
    if e_hi <= 0.0:
        return hi_deg            # even a horizontal aim is below threshold
    if e_lo > 0.0:
        return np.nan            # no reachable aim clears the threshold
    a, b = lo_deg, hi_deg
    for _ in range(40):
        m = 0.5 * (a + b)
        if excess(m) > 0.0:
            b = m
        else:
            a = m
    return 0.5 * (a + b)


def fig_critical_angle(cal):
    """Fig 6: the practical deliverable, and the cleanest result in the model.

    Left panel, vertical wall.  The wall normal is horizontal, so the splash-driving
    velocity is v0 cos(theta), which does NOT depend on the launch height.  The three
    height curves therefore coincide exactly wherever they exist; height enters only
    through reachability (a steep aim from a low start lands short) and through
    coherence.  For a wall, aim is the whole lever.

    Right panel, horizontal pool.  The pool normal is vertical, so the splash-driving
    velocity is the fall speed, which grows as sqrt(2 g h).  Here height is the whole
    lever and aim barely matters.

    That contrast is exactly the confusion v1 fell into: it reasoned about a pool and
    computed a wall.
    """
    dists = np.linspace(0.10, 0.60, 26)
    heights = {"sitting (h = 0.15 m)": 0.15,
               "low stance (h = 0.45 m)": 0.45,
               "standing (h = 0.75 m)": 0.75}

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.6))

    ax = axes[0]
    for (label, h), style in zip(heights.items(), ["o-", "s--", "^:"]):
        crit = np.array([critical_aim_angle(h, dd, cal["lam_mid"]) for dd in dists])
        ax.plot(dists * 100, crit, style, ms=4, lw=1.6, alpha=0.85, label=label)
        gone = np.isnan(crit)
        if gone.any():
            ax.plot(dists[gone] * 100, np.full(gone.sum(), 2.0), "x", ms=4, alpha=0.6)
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.set_xlabel("distance to the wetted wall (cm)")
    ax.set_ylabel("shallowest safe aim (deg below horizontal)")
    ax.set_title("Vertical wall: aim is the lever\n"
                 "curves coincide because $v_{normal} = v_0\\cos\\theta$ is height-free;\n"
                 "crosses at the top mark distances with no safe reachable aim")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)

    ax = axes[1]
    hs = np.linspace(0.05, 1.20, 40)
    for film_mm, style in zip([0.3, 0.5, 1.0], ["o-", "s--", "^:"]):
        exc = []
        for h in hs:
            try:
                r = u.evaluate(Q, D_EXIT, np.radians(-80.0), h,
                               u.Target("pool", 0.0, film_mm * 1e-3), cal["lam_mid"])
                exc.append(-1.0 if r["coherent"] else r["excess"])
            except u.StreamMissesTarget:
                exc.append(np.nan)
        ax.plot(hs * 100, exc, style, ms=4, lw=1.6, alpha=0.85,
                label=f"film = {film_mm:.1f} mm")
    ax.axhline(0.0, color="tab:green", ls="--", lw=1.5, label="splash threshold")
    ax.set_xlabel("fall height to the pool (cm)")
    ax.set_ylabel("splash excess $K/K_c - 1$")
    ax.set_title("Horizontal pool: height is the lever\n"
                 "$v_{normal}$ is the fall speed, so aim hardly enters")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig6_critical_aim_angle.png", dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def report(cal):
    v0 = u.exit_velocity(Q, D_EXIT)
    print("=" * 74)
    print("UROFLOW SPLASHBACK  --  v2 model summary")
    print("=" * 74)
    print(f"Exit: Q = {Q*1e6:.0f} mL/s, d = {D_EXIT*1e3:.0f} mm  ->  v0 = {v0:.2f} m/s")
    print(f"      We = {u.weber(v0, D_EXIT):.0f}, Re = {u.reynolds(v0, D_EXIT):.0f}, "
          f"Oh = {u.ohnesorge(D_EXIT):.4f}")

    print("\n-- Jet thinning (absent from v1) " + "-" * 40)
    for h in (0.15, 0.45, 0.75):
        v = u.jet_speed_after_drop(v0, h)
        d_loc = 2 * u.jet_radius(D_EXIT / 2, v0, v)
        print(f"   after {h*100:4.0f} cm fall: v = {v:.2f} m/s, d = {d_loc*1e3:.2f} mm "
              f"({d_loc/D_EXIT:.0%} of exit)")
    print("   -> v1 used the 4.00 mm exit diameter at impact, over-stating We by "
          f"{D_EXIT/ (2*u.jet_radius(D_EXIT/2, v0, u.jet_speed_after_drop(v0, 0.75))):.2f}x")

    print("\n-- Breakup calibration " + "-" * 50)
    print(f"   v1 ideal (Lambda = 12, parallel jet): "
          f"{u.breakup_length_parallel(Q, D_EXIT, 12.0)*100:.0f} cm")
    print(f"   observed (Splash Lab): {u.L_BREAKUP_OBSERVED[0]*100:.0f}"
          f"-{u.L_BREAKUP_OBSERVED[1]*100:.0f} cm")
    print(f"   calibrated Lambda = {cal['lam_lo']:.2f}-{cal['lam_hi']:.2f}, implying an "
          f"exit disturbance delta0/r0 = {cal['delta_lo']:.1%}-{cal['delta_hi']:.1%}")
    print("   -> a falsifiable prediction: a high-speed clip of the first 2 cm should")
    print("      show a radius ripple of a few per cent, not a smooth column.")

    print("\n-- Post-breakup drop size (v1 used the jet diameter) " + "-" * 20)
    d_break = 2 * u.jet_radius(D_EXIT / 2, v0, u.jet_speed_after_drop(v0, 0.175))
    D_drop = u.rayleigh_drop_diameter(d_break)
    print(f"   jet diameter at breakup: {d_break*1e3:.2f} mm")
    print(f"   Rayleigh drop diameter : {D_drop*1e3:.2f} mm  ({D_drop/d_break:.2f} x d_jet)")

    print("\n-- Geometry correction " + "-" * 50)
    wall = u.Target("wall", plane=0.30, film=FILM)
    for label, (vx, vz) in {"head-on  (horizontal arrival)": (2.0, 0.0),
                            "45 deg   (mixed)": (2.0, -2.0),
                            "grazing  (steep arrival)": (0.3, -3.0)}.items():
        v_n, v_t = u.resolve_impact(vx, vz, wall)
        print(f"   {label}: v_normal = {v_n:.2f} m/s, v_tangential = {v_t:.2f} m/s, "
              f"v1 factor = {np.sin(np.arctan2(abs(vz), abs(vx))):.2f}")
    print("   -> v1's factor is smallest for the head-on case, which is the WORST case")
    print("      for a vertical wall.  v1's 'aim shallow' advice was the wall-geometry")
    print("      inversion of its own horizontal-target reasoning.")

    print("\n-- Optimised configurations (wetted wall, film = "
          f"{FILM*1e3:.1f} mm) " + "-" * 8)
    rows = []
    for h in (0.15, 0.45, 0.75):
        for dist in np.linspace(0.12, 0.55, 12):
            for th in np.radians(np.linspace(-70, 0, 30)):
                try:
                    r = u.evaluate(Q, D_EXIT, th, h, u.Target("wall", dist, FILM),
                                   cal["lam_mid"])
                except ValueError:
                    continue
                score = -1.0 if r["coherent"] else r["excess"]
                if np.isfinite(score):
                    rows.append((score, h, dist, np.degrees(th), r))
    rows.sort(key=lambda t: t[0])
    for tag, row in (("DRIEST", rows[0]), ("WETTEST", rows[-1])):
        score, h, dist, thdeg, r = row
        print(f"   {tag}: h = {h*100:.0f} cm, distance = {dist*100:.0f} cm, "
              f"aim = {thdeg:+.0f} deg")
        print(f"      regime = {r['regime']}, path = {r['path']*100:.0f} cm, "
              f"incidence = {r['incidence_deg']:.0f} deg from normal")
        print(f"      v_normal = {r['v_normal']:.2f} m/s, splash excess = {score:+.2f}")
        if r["regime"] == "fragmented":
            print(f"      drop D = {r['drop_diameter']*1e3:.2f} mm, "
                  f"drag speed deficit = "
                  f"{r['criterion'].get('drag_speed_deficit', 0.0):.1%}")

    print("\n-- The practical answer: how steeply must you aim? " + "-" * 22)
    print("   Shallowest aim (deg below horizontal) that stays under the threshold:")
    print("      distance:      15 cm    30 cm    45 cm")
    for label, h in {"sitting  (0.15 m)": 0.15, "standing (0.75 m)": 0.75}.items():
        vals = [critical_aim_angle(h, dd, cal["lam_mid"]) for dd in (0.15, 0.30, 0.45)]
        cells = "  ".join("  none " if not np.isfinite(v) else f"{v:6.0f} " for v in vals)
        print(f"      {label}: {cells}")
    print("   'none' means no reachable aim clears the threshold at that distance, so")
    print("   the only remaining levers are distance and flow rate.")
    print("   The two heights give the SAME angle wherever both are reachable, and that")
    print("   is not a coincidence: at a vertical wall the splash-driving component is")
    print("   v0*cos(theta), which contains no height.  Height enters only through")
    print("   reachability and coherence.  At a horizontal pool the roles swap entirely,")
    print("   because there the normal component is the fall speed.  Conflating the two")
    print("   geometries is precisely the error v1 made.")

    print("\n-- Recommendations, v2 " + "-" * 50)
    print("   1. Minimise the NORMAL velocity component.  Against a vertical wall that")
    print("      means arriving steeply, not shallowly.  v1 had this backwards, because")
    print("      it reasoned about a horizontal target while integrating to a vertical")
    print("      one.  The correction reverses the sign of the aim-angle advice.")
    print("   2. Height is not the lever; path length and incidence are.  A short path")
    print("      helps only because it can stay inside the 15-20 cm coherent length,")
    print("      and a long path helps only because gravity has steepened the arrival.")
    print("      The worst case in the whole scan is a SHORT path aimed horizontally.")
    print("   3. A wetted surface splashes LESS, monotonically over the correlated")
    print("      range, because Cossali's threshold rises with film depth.  Drying the")
    print("      target is counterproductive.")
    print("   4. Drag on the fragments is small but not zero; the model now computes it")
    print("      rather than asserting it (a fraction of a per cent at these sizes).")

    print("\n-- Honest limits " + "-" * 55)
    print("   * Cossali's correlation is fitted for 0.1 <= delta <= 1.  A deep bowl is")
    print("     delta >> 1 and the model flags, but does not correct, the extrapolation.")
    print("   * The meatus is a slit; the effective circular diameter lumps in stream")
    print("     splitting, which is the largest unmodelled geometric effect.")
    print("   * Lambda is calibrated against one reported observation range, so the")
    print("     breakup point inherits that anchor's uncertainty.")
    print("   * Gravitational stretching stabilises a falling jet and is not modelled,")
    print("     so the ideal breakup length is a lower bound.")


def main():
    cal = calibration()
    fig_breakup(cal)
    fig_geometry_correction()
    fig_splash_map(cal)
    fig_film_depth(cal)
    fig_sit_vs_stand(cal)
    fig_critical_angle(cal)
    report(cal)
    print("\nFigures written to figures/.")


if __name__ == "__main__":
    main()
