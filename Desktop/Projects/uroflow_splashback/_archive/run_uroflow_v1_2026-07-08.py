"""
run_uroflow.py -- Splashback maps and recommendations for a urine stream.

Reproduces (qualitatively) the Splash Lab (Hurd & Truscott, APS-DFD 2013) findings
from first principles and turns them into a quantitative optimisation:

  * coherent (pre-breakup) impact splashes far less than a drop train,
  * grazing impact angle beats normal impact,
  * lower launch height (sitting) shortens the path and slows the impact.

Figures -> figures/.  All units SI.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import uroflow as u

FIGDIR = "figures"

# Representative healthy adult male peak flow and effective meatus diameter.
Q = 21e-6          # m^3/s  (21 mL/s)
D = 4e-3           # m      effective stream diameter


def fig_breakup_vs_flow():
    """Fig 1: breakup length vs flow rate & diameter, with the Splash Lab
    empirical anchor (real streams break sooner than ideal theory)."""
    Qs = np.linspace(8e-6, 28e-6, 60)
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for d in [3e-3, 4e-3, 5e-3]:
        Lb = np.array([u.breakup_length(q, d) for q in Qs]) * 100
        ax.plot(Qs * 1e6, Lb, label=f"d = {d*1e3:.0f} mm")
    ax.axhspan(15, 20, color="orange", alpha=0.2,
               label="Splash Lab observed breakup ~15-20 cm")
    ax.set_xlabel("peak flow rate Q (mL/s)")
    ax.set_ylabel("ideal breakup length $L_b$ (cm)")
    ax.set_title("Plateau-Rayleigh breakup length (Weber theory, ln$(r/\\delta_0)$=12)\n"
                 "Ideal theory over-predicts $L_b$; real perturbed streams break sooner")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig1_breakup_length.png", dpi=150)
    plt.close(fig)


def fig_splash_map():
    """Fig 2: splash index over (standing distance, aim angle) for a standing
    user. Shows the coherent 'sweet spot' close to the target."""
    h = 0.75                                   # stream launch height above target impact zone (m)
    distances = np.linspace(0.10, 0.70, 60)    # horizontal distance to target wall (m)
    angles_deg = np.linspace(-45, 5, 60)       # launch angle (deg); negative = aimed downward
    S = np.zeros((len(angles_deg), len(distances)))
    COH = np.zeros_like(S)
    for i, th in enumerate(np.radians(angles_deg)):
        for j, dist in enumerate(distances):
            r = u.splash_index(Q, D, th, h, dist)
            S[i, j] = r["index"] if not r["hit_floor"] else np.nan
            COH[i, j] = 1.0 if r["coherent"] else 0.0

    fig, ax = plt.subplots(figsize=(8, 5))
    pcm = ax.pcolormesh(distances * 100, angles_deg, S, shading="auto", cmap="RdYlGn_r")
    fig.colorbar(pcm, ax=ax, label="splash index (lower = drier)")
    # Overlay the pre-breakup (coherent) boundary.
    ax.contour(distances * 100, angles_deg, COH, levels=[0.5], colors="blue", linewidths=2)
    ax.text(0.30 * 100, -40, "coherent stream\n(pre-breakup)", color="blue", fontsize=9)
    ax.set_xlabel("standing distance to target (cm)")
    ax.set_ylabel("aim angle (deg, negative = downward)")
    ax.set_title("Splashback map (standing).  Blue curve = breakup boundary.\n"
                 "Drier: stand close (stay coherent) and aim shallow/downward")
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig2_splash_map.png", dpi=150)
    plt.close(fig)


def fig_sit_vs_stand():
    """Fig 3: splash index vs standing distance for several launch heights
    (proxy for sitting vs standing)."""
    heights = {"sitting  (h=0.15 m)": 0.15, "low stance (h=0.45 m)": 0.45,
               "standing (h=0.75 m)": 0.75}
    dists = np.linspace(0.10, 0.60, 60)
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for label, h in heights.items():
        # Aim so the stream lands at the target: choose angle to hit (dist, 0).
        idx = []
        for dist in dists:
            # pick the aim angle minimising splash at this distance (optimise aim)
            best = min(
                (u.splash_index(Q, D, th, h, dist)
                 for th in np.radians(np.linspace(-60, 0, 25))),
                key=lambda rr: rr["index"] if not rr["hit_floor"] else 1e9)
            idx.append(best["index"] if not best["hit_floor"] else np.nan)
        ax.plot(dists * 100, idx, "o-", ms=3, label=label)
    ax.set_xlabel("horizontal distance to target (cm)")
    ax.set_ylabel("min splash index (aim-optimised)")
    ax.set_title("Lower launch height (sitting) reduces splashback at every distance")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig3_sit_vs_stand.png", dpi=150)
    plt.close(fig)


def print_recommendations():
    print("=" * 68)
    print("UROFLOW SPLASHBACK -- model summary")
    print("=" * 68)
    v = u.exit_velocity(Q, D)
    print(f"Exit: Q={Q*1e6:.0f} mL/s, d={D*1e3:.0f} mm -> v={v:.2f} m/s, "
          f"We={u.weber(v, D):.0f}, Re={u.reynolds(v, D):.0f}, Oh={u.ohnesorge(D):.4f}")
    print(f"Ideal breakup length L_b = {u.breakup_length(Q, D)*100:.0f} cm "
          f"(real streams break sooner; Splash Lab ~15-20 cm)")

    # Scan configs and report the driest.
    configs = []
    for h_ in [0.15, 0.45, 0.75]:
        for dist in np.linspace(0.12, 0.55, 12):
            for th in np.radians(np.linspace(-55, 0, 20)):
                r = u.splash_index(Q, D, th, h_, dist)
                if not r["hit_floor"]:
                    configs.append((r["index"], h_, dist, np.degrees(th), r))
    configs.sort(key=lambda c: c[0])
    print("\nDriest configuration found:")
    idx, h_, dist, thdeg, r = configs[0]
    print(f"  height {h_*100:.0f} cm, distance {dist*100:.0f} cm, aim {thdeg:+.0f} deg")
    print(f"  -> path {r['path']*100:.0f} cm ({'coherent' if r['coherent'] else 'BROKEN into drops'}), "
          f"impact angle {r['angle_deg']:.0f} deg from horizontal, splash index {idx:.1f}")
    print("\nWettest (worst) configuration found:")
    idx, h_, dist, thdeg, r = configs[-1]
    print(f"  height {h_*100:.0f} cm, distance {dist*100:.0f} cm, aim {thdeg:+.0f} deg")
    print(f"  -> path {r['path']*100:.0f} cm ({'coherent' if r['coherent'] else 'BROKEN into drops'}), "
          f"impact angle {r['angle_deg']:.0f} deg, splash index {idx:.1f}")
    print("\nRecommendations (physics-backed, matching Splash Lab):")
    print("  1. Stand close  -> impact before Rayleigh breakup (coherent << drops).")
    print("  2. Aim shallow  -> grazing impact, sin(angle) small, less normal momentum.")
    print("  3. Sit / lower  -> shorter path + slower impact + easier shallow aim.")


def main():
    fig_breakup_vs_flow()
    fig_splash_map()
    fig_sit_vs_stand()
    print_recommendations()
    print("\nFigures written to figures/.")


if __name__ == "__main__":
    main()
