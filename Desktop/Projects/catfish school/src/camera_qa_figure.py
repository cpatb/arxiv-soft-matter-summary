"""Why camera-motion removal is not optional.

The video is handheld. Apparent motion of the school is the sum of the fish
swimming and the phone moving, and the phone's contribution is coherent across
the whole field of view -- which is exactly what an order parameter is built to
detect. Left in, a pan masquerades as perfect collective alignment.

This figure compares polarisation and speed measured with and without the
correction, across many frames, so the size of the artefact is visible rather
than asserted.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from catfish import flow, frames, plotstyle, segment  # noqa: E402

RESULTS = frames.ROOT / "results"
FIGS = frames.ROOT / "figures"

N_FRAMES = 60
SEED = 20260808


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
            raw = flow.sample_school_velocity(fl, env.mask, work, camera=None)
            cor = flow.sample_school_velocity(fl, env.mask, work, camera=cam)
            o_raw = flow.order_parameters(raw["u"], raw["v"])
            o_cor = flow.order_parameters(cor["u"], cor["v"])
            c = cam.velocity_at(np.array([g1.shape[0] / 2]),
                                np.array([g1.shape[1] / 2]))[0]
            rows.append(
                {
                    "clip": name,
                    "frame": i,
                    "cam_shift_px": float(np.hypot(c[0], c[1])),
                    "pol_raw": o_raw["polarisation"],
                    "pol_cor": o_cor["polarisation"],
                    "speed_raw": o_raw["median_speed_px"],
                    "speed_cor": o_cor["median_speed_px"],
                }
            )
        print(f"{name}: {sum(r['clip']==name for r in rows)} frames", flush=True)

    df = pd.DataFrame(rows)
    if df.empty:
        print("no data")
        return
    df.to_csv(RESULTS / "camera_qa.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.0))

    for name, g in df.groupby("clip"):
        c = plotstyle.CLIP_COLOUR[name]
        axes[0].scatter(g["pol_raw"], g["pol_cor"], s=12, alpha=0.7, color=c,
                        linewidths=0, label=name.replace("IMG_", "clip "))
        axes[1].scatter(g["cam_shift_px"], g["pol_raw"] - g["pol_cor"], s=12,
                        alpha=0.7, color=c, linewidths=0)
    axes[0].plot([0, 1], [0, 1], color="0.5", ls="--", lw=1.0)
    axes[0].set_xlabel(r"$\Phi$  uncorrected")
    axes[0].set_ylabel(r"$\Phi$  camera removed")
    axes[0].set_xlim(0, 1.02)
    axes[0].set_ylim(0, 1.02)
    axes[0].legend(loc="lower right")

    axes[1].axhline(0, color="0.5", lw=0.8, ls="--")
    axes[1].set_xlabel(r"camera shift  (px frame$^{-1}$)")
    axes[1].set_ylabel(r"$\Delta\Phi$  (artefact)")
    axes[1].set_xscale("symlog", linthresh=1.0)

    axes[2].hist(df["pol_raw"], bins=25, range=(0, 1), histtype="step",
                 color="0.5", density=True, label="uncorrected")
    axes[2].hist(df["pol_cor"], bins=25, range=(0, 1), histtype="step",
                 color="#0072B2", density=True, label="camera removed")
    axes[2].set_xlabel(r"polarisation  $\Phi$")
    axes[2].set_ylabel("density")
    axes[2].legend(loc="upper left")

    for ax, L in zip(axes, "abc"):
        plotstyle.panel_label(ax, L)
    fig.tight_layout()
    fig.savefig(FIGS / "figS2_camera_qa.png")
    plt.close(fig)

    worst = df.loc[df["pol_raw"].sub(df["pol_cor"]).idxmax()]
    summary = {
        "n_frames": int(len(df)),
        "pol_raw_median": float(df["pol_raw"].median()),
        "pol_corrected_median": float(df["pol_cor"].median()),
        "median_inflation": float((df["pol_raw"] - df["pol_cor"]).median()),
        "max_inflation": float((df["pol_raw"] - df["pol_cor"]).max()),
        "worst_case": {
            "clip": worst["clip"], "frame": int(worst["frame"]),
            "cam_shift_px": float(worst["cam_shift_px"]),
            "pol_raw": float(worst["pol_raw"]), "pol_cor": float(worst["pol_cor"]),
        },
        "frac_frames_inflated_by_over_0.2": float(
            ((df["pol_raw"] - df["pol_cor"]) > 0.2).mean()
        ),
        "camera_shift_px_median": float(df["cam_shift_px"].median()),
        "camera_shift_px_p95": float(df["cam_shift_px"].quantile(0.95)),
    }
    with open(RESULTS / "camera_qa.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
