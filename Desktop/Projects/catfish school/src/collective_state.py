"""Where these schools sit in the standard collective-state diagram.

Tunstrom, Katz, Ioannou, Huepe, Lutz and Couzin (PLoS Comput Biol 9:e1002915,
2013) showed that schooling fish occupy at least three dynamically stable
states, separated by two order parameters:

* **Polarisation** Phi = |<v_hat>|, the alignment of swimming directions.
* **Rotation** M = |<r_hat x v_hat>|, the coherence of circulation about the
  group's centre, where r is measured from the centroid.

    swarm      low Phi, low M    -- cohesive but internally disordered
    milling    low Phi, high M   -- a rotating torus; individuals aligned with
                                   their neighbours but not with the group
    polarised  high Phi, low M   -- everyone travelling the same way

Polarisation alone cannot tell a swarm from a mill: both have low Phi and they
are entirely different behaviours. Measuring only Phi would have described a
rotating fry ball as "disordered", which is wrong.

Also computed here: whether the rind of the school behaves differently from its
core, which is the border-versus-interior question that recurs in the flocking
literature.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import cv2
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from catfish import flow, frames, plotstyle, segment  # noqa: E402

RESULTS = frames.ROOT / "results"
FIGS = frames.ROOT / "figures"

N_FRAMES = 70
SEED = 20260808

#: Depth defining the "core", in body lengths, measured inward from the boundary.
CORE_DEPTH_BL = 1.5

#: Speeds below this fraction of the frame's median are treated as directionless
#: and excluded from the angular order parameters.
SPEED_FLOOR_FRAC = 0.1


def rotational_order(y, x, u, v, cy, cz) -> float:
    """|<r_hat x v_hat>| about the centroid. 1 = perfectly circulating.

    Image row index increases downward. That flips the sign of the cross
    product relative to a right-handed frame, which is irrelevant here because
    the magnitude is taken, but it would matter if a direction of rotation were
    ever reported.
    """
    ry, rx = y - cy, x - cz
    rn = np.hypot(ry, rx)
    sp = np.hypot(u, v)
    ok = (rn > 1e-9) & (sp > 1e-9)
    if ok.sum() < 5:
        return float("nan")
    ry, rx, uu, vv, rn, sp = ry[ok], rx[ok], u[ok], v[ok], rn[ok], sp[ok]
    # z-component of r_hat x v_hat
    cross = (rx / rn) * (vv / sp) - (ry / rn) * (uu / sp)
    return float(abs(np.mean(cross)))


def main() -> None:
    plotstyle.apply()
    FIGS.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    rows = []
    for name in frames.CLIPS:
        mp, fp = RESULTS / f"{name}_meta.json", RESULTS / f"{name}_frames.csv"
        if not (mp.exists() and fp.exists()):
            print(f"skip {name}: pipeline output missing")
            continue
        with open(mp) as fh:
            work = json.load(fh)["working_body_px"]
        clip = frames.load_clip(name)
        fr = pd.read_csv(fp)
        present = fr.loc[fr["present"], "frame"].to_numpy()
        present = present[present < len(clip) - 1]
        if len(present) < 5:
            continue
        picks = rng.choice(present, size=min(N_FRAMES, len(present)), replace=False)

        for i in sorted(picks):
            i = int(i)
            f1 = segment.compute_fields(clip.rgb(i), work)
            g1 = np.clip(f1.L * 255 / 100, 0, 255).astype(np.uint8)
            f2 = segment.compute_fields(clip.rgb(i + 1), work)
            g2 = np.clip(f2.L * 255 / 100, 0, 255).astype(np.uint8)
            env = segment.school_envelope(fields=f1)
            if not env.present:
                continue
            cam = flow.estimate_camera_motion(g1, g2, exclude=env.mask)
            if not cam.ok:
                continue
            fl = flow.dense_flow(g1, g2)
            s = flow.sample_school_velocity(fl, env.mask, work, camera=cam)
            if len(s["u"]) < 60:
                continue

            y, x, u, v = s["y"], s["x"], s["u"], s["v"]
            sp = np.hypot(u, v)
            floor = max(0.05, SPEED_FLOOR_FRAC * float(np.median(sp)))
            keep = sp > floor
            if keep.sum() < 40:
                continue

            cy, cx = float(y[keep].mean()), float(x[keep].mean())
            phi = flow.order_parameters(u[keep], v[keep])["polarisation"]
            M = rotational_order(y[keep], x[keep], u[keep], v[keep], cy, cx)

            # rind vs core
            r_core = max(1, int(round(CORE_DEPTH_BL * work)))
            se = cv2.getStructuringElement(
                cv2.MORPH_ELLIPSE, (2 * r_core + 1, 2 * r_core + 1)
            )
            core_mask = cv2.erode(env.mask.astype(np.uint8), se).astype(bool)
            yi = np.clip(y.astype(int), 0, core_mask.shape[0] - 1)
            xi = np.clip(x.astype(int), 0, core_mask.shape[1] - 1)
            in_core = core_mask[yi, xi] & keep
            in_rind = (~core_mask[yi, xi]) & keep

            def sub(mask):
                if mask.sum() < 25:
                    return np.nan, np.nan
                o = flow.order_parameters(u[mask], v[mask])
                return o["polarisation"], float(np.median(np.hypot(u[mask], v[mask])))

            pol_core, spd_core = sub(in_core)
            pol_rind, spd_rind = sub(in_rind)

            rows.append(
                {
                    "clip": name, "frame": i,
                    "polarisation": phi, "rotation": M,
                    "pol_core": pol_core, "pol_rind": pol_rind,
                    "speed_core_px": spd_core, "speed_rind_px": spd_rind,
                    "n_core": int(in_core.sum()), "n_rind": int(in_rind.sum()),
                }
            )
        print(f"{name}: {sum(r['clip']==name for r in rows)} frames", flush=True)

    df = pd.DataFrame(rows)
    if df.empty:
        print("no data")
        return
    df.to_csv(RESULTS / "collective_state.csv", index=False)

    # --- figure ----------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.0))
    for name, g in df.groupby("clip"):
        c = plotstyle.CLIP_COLOUR[name]
        axes[0].scatter(g["polarisation"], g["rotation"], s=11, alpha=0.6,
                        color=c, linewidths=0, label=name.replace("IMG_", "clip "))
        axes[1].scatter(g["pol_rind"], g["pol_core"], s=11, alpha=0.6,
                        color=c, linewidths=0)
        axes[2].scatter(g["speed_rind_px"], g["speed_core_px"], s=11, alpha=0.6,
                        color=c, linewidths=0)
    axes[0].set_xlabel(r"polarisation  $\Phi$")
    axes[0].set_ylabel(r"rotation  $M$")
    axes[0].set_xlim(0, 1)
    axes[0].set_ylim(0, 1)
    axes[0].legend(loc="upper right")

    # Polarisation is bounded on [0, 1]; speed is not, and a couple of frames
    # caught mid-pan reach 40 px/frame, which would squash every real point into
    # the corner. Clip the speed axes to a high percentile of the data instead of
    # letting two outliers set the scale.
    axes[1].set_xlim(0, 1)
    axes[1].set_ylim(0, 1)
    axes[1].plot([0, 1], [0, 1], color="0.5", ls="--", lw=1.0)

    spd = pd.concat([df["speed_rind_px"], df["speed_core_px"]]).dropna()
    hi = float(np.nanpercentile(spd, 98)) if len(spd) else 1.0
    axes[2].set_xlim(0, hi)
    axes[2].set_ylim(0, hi)
    axes[2].plot([0, hi], [0, hi], color="0.5", ls="--", lw=1.0)
    axes[1].set_xlabel(r"$\Phi$  rind")
    axes[1].set_ylabel(r"$\Phi$  core")
    axes[2].set_xlabel("speed rind")
    axes[2].set_ylabel("speed core")

    for ax, L in zip(axes, "abc"):
        plotstyle.panel_label(ax, L)
    fig.tight_layout()
    fig.savefig(FIGS / "fig9_collective_state.png")
    plt.close(fig)

    summary = {}
    for name, g in df.groupby("clip"):
        phi, M = g["polarisation"].median(), g["rotation"].median()
        # Classification thresholds follow the qualitative scheme in the
        # literature; they are conventions, not measurements, and the raw Phi
        # and M values are reported alongside so a reader can reclassify.
        if phi > 0.65:
            state = "polarised"
        elif M > 0.45:
            state = "milling"
        else:
            state = "swarm"
        summary[name] = {
            "n_frames": int(len(g)),
            "polarisation_median": float(phi),
            "polarisation_iqr": [float(g["polarisation"].quantile(.25)),
                                 float(g["polarisation"].quantile(.75))],
            "rotation_median": float(M),
            "rotation_iqr": [float(g["rotation"].quantile(.25)),
                             float(g["rotation"].quantile(.75))],
            "state": state,
            "pol_core_median": float(g["pol_core"].median()),
            "pol_rind_median": float(g["pol_rind"].median()),
            "speed_core_over_rind": float(
                g["speed_core_px"].median() / g["speed_rind_px"].median()
            ) if g["speed_rind_px"].median() > 0 else np.nan,
        }
    with open(RESULTS / "collective_state.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
