"""
run_experiments.py -- Validate the ballistic region-of-origin method and quantify
the height bias of the classical straight-line ("stringing") method.

Produces four figures in figures/ and prints a numeric report.

Headline question: by how much does the standard forensic straight-line method
overestimate the height of origin, and does the gravity+drag back-projection fix it?
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

import bpa

FIGDIR = "figures"
rng_global = np.random.default_rng(42)


def scene_figure():
    """Fig 1: one scene -- floor stains, straight-line back-projections (converge
    high) vs ballistic back-trajectories (converge at truth)."""
    origin = np.array([0.0, 0.0, 1.15])
    stains = bpa.emit_and_land(origin, n_drops=250, speed_range=(3.5, 7.0),
                               d0_range=(2.8e-3, 3.2e-3),
                               cone_half_angle=np.radians(30), seed=7, drag=True)
    x0, y0 = bpa.convergence_2d(stains)
    z_sl, z_sl_sd, _ = bpa.straight_line_height(stains, x0, y0)
    rec, dists, trajs = bpa.ballistic_origin(stains, drag=True)

    fig = plt.figure(figsize=(11, 5))
    # Left: 3D.  Clip everything to a readable box around the true origin so the
    # (physically-harmless) far tails of the back-integrated paths don't blow up
    # the axes.
    zmax = 2.2 * origin[2]
    rmax = 1.0                                       # +/- 1 m horizontal window
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    xs = [s.x for s in stains if abs(s.x) < rmax and abs(s.y) < rmax]
    ys = [s.y for s in stains if abs(s.x) < rmax and abs(s.y) < rmax]
    ax.scatter(xs, ys, [0] * len(xs), s=8, c="crimson", alpha=0.5, label="stains (floor)")
    for s in stains[:60]:
        d = np.hypot(s.x - x0, s.y - y0)
        zt = min(d * np.tan(s.alpha), zmax)          # clip straight-line ray
        f = zt / max(d * np.tan(s.alpha), 1e-9)
        xe = s.x + f * (x0 - s.x); ye = s.y + f * (y0 - s.y)
        ax.plot([s.x, xe], [s.y, ye], [0, zt], color="gray", lw=0.4, alpha=0.4)
    for tr in trajs[:60]:
        m = (tr[:, 2] >= 0) & (tr[:, 2] <= zmax) & \
            (np.abs(tr[:, 0]) < rmax) & (np.abs(tr[:, 1]) < rmax)
        ax.plot(tr[m, 0], tr[m, 1], tr[m, 2], color="steelblue", lw=0.4, alpha=0.6)
    ax.scatter(*origin, c="black", s=80, marker="*", label="true origin")
    ax.scatter(rec[0], rec[1], rec[2], c="blue", s=40, marker="o", label="ballistic est.")
    ax.set_xlim(-rmax, rmax); ax.set_ylim(-rmax, rmax); ax.set_zlim(0, zmax)
    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)"); ax.set_zlabel("z (m)")
    ax.set_title("Back-projection: straight-line (gray) vs ballistic (blue)")
    ax.legend(fontsize=7, loc="upper left")
    ax.view_init(elev=12, azim=-60)

    # Right: height comparison bar
    ax2 = fig.add_subplot(1, 2, 2)
    labels = ["true", "straight-line", "ballistic"]
    vals = [origin[2], z_sl, rec[2]]
    errs = [0, z_sl_sd, float(np.std([t[np.argmin(np.abs(t[:, 2] - rec[2])), 2] for t in trajs]))]
    colors = ["black", "gray", "steelblue"]
    ax2.bar(labels, vals, yerr=errs, color=colors, alpha=0.8, capsize=5)
    ax2.set_ylabel("height of origin z (m)")
    ax2.axhline(origin[2], color="black", ls="--", lw=1)
    ax2.set_title(f"straight-line overestimates by "
                  f"{100*(z_sl-origin[2])/origin[2]:+.0f}%\n"
                  f"ballistic error {100*(rec[2]-origin[2])/origin[2]:+.1f}%")
    ax2.set_ylim(0, max(vals) * 1.25)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig1_scene.png", dpi=150)
    plt.close(fig)
    return origin[2], z_sl, rec[2]


def bias_sweep():
    """Fig 2: straight-line vs ballistic height error across a sweep of true
    heights and mean launch speeds. This is the quantitative headline."""
    heights = [0.6, 0.9, 1.2, 1.5, 1.8]
    speeds = [(2.5, 4.0), (4.0, 6.0), (6.0, 9.0)]
    speed_labels = ["slow 2.5-4", "med 4-6", "fast 6-9"]
    results = np.zeros((len(speeds), len(heights), 2))  # [., ., (straight, ballistic)] rel error %

    for i, sr in enumerate(speeds):
        for j, h in enumerate(heights):
            origin = np.array([0.0, 0.0, h])
            stains = bpa.emit_and_land(origin, n_drops=350, speed_range=sr,
                                       d0_range=(2.8e-3, 3.2e-3),
                                       cone_half_angle=np.radians(32),
                                       seed=100 + 7 * i + j, drag=True)
            if len(stains) < 20:
                results[i, j] = np.nan
                continue
            x0, y0 = bpa.convergence_2d(stains)
            z_sl, _, _ = bpa.straight_line_height(stains, x0, y0)
            rec, _, _ = bpa.ballistic_origin(stains, drag=True)
            results[i, j, 0] = 100 * (z_sl - h) / h
            results[i, j, 1] = 100 * (rec[2] - h) / h

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True)
    for i, lab in enumerate(speed_labels):
        axes[0].plot(heights, results[i, :, 0], "o-", label=lab)
        axes[1].plot(heights, results[i, :, 1], "o-", label=lab)
    axes[0].set_title("Straight-line (stringing) height error")
    axes[1].set_title("Ballistic back-projection height error")
    for ax in axes:
        ax.axhline(0, color="k", lw=0.8, ls="--")
        ax.set_xlabel("true height of origin (m)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, title="launch speed (m/s)")
    axes[0].set_ylabel("height error (%)")
    fig.suptitle("Straight-line stringing systematically overestimates height; "
                 "ballistic method removes the gravity bias", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig2_bias_sweep.png", dpi=150)
    plt.close(fig)
    return heights, speed_labels, results


def spread_figure():
    """Fig 3: the drop-impact bridge -- spread factor and spine count vs Weber."""
    d0 = 3.0e-3
    v_n = np.linspace(0.5, 12, 200)
    we = bpa.RHO_BLOOD * d0 * v_n**2 / bpa.SIGMA_BLOOD
    re = bpa.RHO_BLOOD * d0 * v_n / bpa.MU_BLOOD
    xi = bpa.spread_factor(we, re)
    spines = np.array([bpa.n_spines(w) for w in we])

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(we, xi, color="crimson")
    ax[0].set_xlabel("normal Weber number  We$_n$")
    ax[0].set_ylabel(r"spread factor  $\xi = D_{stain}/D_0$")
    ax[0].set_title("Stain size encodes impact speed\n(invertible: measure $\\xi$ $\\to$ We $\\to$ $v_n$)")
    ax[0].grid(alpha=0.3)
    ax2 = ax[0].twiny()
    ax2.set_xlim(v_n.min(), v_n.max()); ax2.set_xlabel("normal impact speed $v_n$ (m/s)")

    ax[1].plot(we, spines, color="steelblue")
    ax[1].set_xlabel("normal Weber number  We$_n$")
    ax[1].set_ylabel("expected rim spines  $n \\approx 1.14\\,$We$^{1/2}$")
    ax[1].set_title("Spine count: a size-independent\nWeber estimator when $D_0$ unknown")
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig3_spread_bridge.png", dpi=150)
    plt.close(fig)


def noise_robustness():
    """Fig 4: reconstruction error under realistic measurement noise on the
    ellipse (W, L), azimuth, and parent drop size D_0."""
    origin = np.array([0.0, 0.0, 1.2])
    base = bpa.emit_and_land(origin, n_drops=220, speed_range=(4.0, 7.0),
                             d0_range=(2.8e-3, 3.2e-3),
                             cone_half_angle=np.radians(32), seed=11, drag=True)
    rng = np.random.default_rng(2024)
    n_trials = 25
    sl_err, ba_err = [], []
    for _ in range(n_trials):
        noisy = []
        for s in base:
            W = s.W * (1 + 0.05 * rng.standard_normal())      # 5% size noise
            L = s.L * (1 + 0.05 * rng.standard_normal())
            L = max(L, W * 1.001)                             # keep W<=L
            az = s.azimuth + np.radians(3) * rng.standard_normal()  # 3 deg azimuth noise
            d0 = s.d0 * (1 + 0.10 * rng.standard_normal())    # 10% drop-size uncertainty
            noisy.append(bpa.Stain(s.x, s.y, s.z, L, W, az, max(d0, 1e-3)))
        x0, y0 = bpa.convergence_2d(noisy)
        z_sl, _, _ = bpa.straight_line_height(noisy, x0, y0)
        rec, _, _ = bpa.ballistic_origin(noisy, drag=True)
        sl_err.append(z_sl - origin[2])
        ba_err.append(rec[2] - origin[2])
    sl_err = np.array(sl_err); ba_err = np.array(ba_err)

    fig, ax = plt.subplots(figsize=(7, 4.3))
    ax.hist(sl_err * 100, bins=15, alpha=0.6, color="gray", label="straight-line")
    ax.hist(ba_err * 100, bins=15, alpha=0.6, color="steelblue", label="ballistic")
    ax.axvline(0, color="k", ls="--", lw=1)
    ax.set_xlabel("height-of-origin error (cm)")
    ax.set_ylabel("count over %d noisy trials" % n_trials)
    ax.set_title("Under 5%% size / 3deg azimuth / 10%% drop-size noise:\n"
                 "straight-line bias %+.1f+-%.1f cm  vs  ballistic %+.1f+-%.1f cm"
                 % (100*sl_err.mean(), 100*sl_err.std(),
                    100*ba_err.mean(), 100*ba_err.std()))
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/fig4_noise.png", dpi=150)
    plt.close(fig)
    return sl_err, ba_err


def main():
    print("=" * 70)
    print("BLOODSTAIN PATTERN ANALYSIS -- ballistic region-of-origin")
    print("=" * 70)
    zt, zsl, zba = scene_figure()
    print(f"\n[Fig 1] single scene, true height {zt:.3f} m:")
    print(f"    straight-line estimate : {zsl:.3f} m  ({100*(zsl-zt)/zt:+.1f}%)")
    print(f"    ballistic estimate     : {zba:.3f} m  ({100*(zba-zt)/zt:+.1f}%)")

    heights, labels, res = bias_sweep()
    print("\n[Fig 2] height-error sweep (mean over speed regimes):")
    sl_mean = np.nanmean(res[:, :, 0])
    ba_mean = np.nanmean(res[:, :, 1])
    print(f"    straight-line mean error : {sl_mean:+.1f}%")
    print(f"    ballistic mean error     : {ba_mean:+.1f}%")
    print(f"    -> straight-line bias reduced ~{abs(sl_mean/ba_mean):.0f}x")

    spread_figure()
    print("\n[Fig 3] spread-factor bridge written.")

    sl_err, ba_err = noise_robustness()
    print("\n[Fig 4] noise robustness:")
    print(f"    straight-line : {100*sl_err.mean():+.1f} +- {100*sl_err.std():.1f} cm")
    print(f"    ballistic     : {100*ba_err.mean():+.1f} +- {100*ba_err.std():.1f} cm")
    print("\nAll figures in figures/. Done.")


if __name__ == "__main__":
    main()
