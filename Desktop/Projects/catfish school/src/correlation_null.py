"""Does the measured correlation length reflect real spatial structure?

Subtracting the school's mean velocity constrains the fluctuations to sum to
zero over the sampled points. That constraint alone forces C(r) negative
somewhere, so *any* field -- including pure noise -- yields a finite zero
crossing. A correlation length is therefore not self-evidently meaningful, and
quoting one without a null is how a measurement of nothing gets published as a
measurement of something.

The null here is a permutation test. Velocity vectors are shuffled among the
sampled positions, which destroys all spatial arrangement while preserving
exactly the set of velocities, the number of points, the school's outline, and
the mean-subtraction constraint. Everything that could produce a spurious
correlation length survives the shuffle; only genuine spatial structure does
not.

If the observed xi sits inside the shuffled distribution, there is no evidence
for spatial correlation beyond the trivial constraint. If it sits far outside,
the correlation is real.
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

#: Frames sampled per clip.
N_FRAMES = 22

#: Shuffles per frame.
N_SHUFFLE = 40

#: Common separation grid in body lengths for averaging C(r).
GRID = np.linspace(0.0, 25.0, 100)

SEED = 20260808


def measure_frame(clip: frames.Clip, i: int, work_px: float):
    """Return (r, C) for the real field and a list of (r, C) for shuffles."""
    f1 = segment.compute_fields(clip.rgb(i), body_px=work_px)
    env = segment.school_envelope(fields=f1)
    if not env.present:
        return None
    g1 = np.clip(f1.L * 255 / 100, 0, 255).astype(np.uint8)
    f2 = segment.compute_fields(clip.rgb(i + 1), body_px=work_px)
    g2 = np.clip(f2.L * 255 / 100, 0, 255).astype(np.uint8)

    cam = flow.estimate_camera_motion(g1, g2, exclude=env.mask)
    fl = flow.dense_flow(g1, g2)
    s = flow.sample_school_velocity(fl, env.mask, work_px, camera=cam)
    if len(s["u"]) < 80:
        return None

    real = flow.correlation_function(s["y"], s["x"], s["u"], s["v"], work_px, n_bins=40)
    if not len(real["r_bl"]):
        return None

    rng = np.random.default_rng(SEED + i)
    shuffles = []
    n = len(s["u"])
    for _ in range(N_SHUFFLE):
        p = rng.permutation(n)
        c = flow.correlation_function(
            s["y"], s["x"], s["u"][p], s["v"][p], work_px, n_bins=40
        )
        if len(c["r_bl"]):
            shuffles.append(c)
    return real, shuffles


def interp(c) -> np.ndarray:
    return np.interp(GRID, c["r_bl"], c["C"], left=np.nan, right=np.nan)


def main() -> None:
    plotstyle.apply()
    FIGS.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    bl_path = RESULTS / "bodylength_summary.json"
    bl_summary = json.load(open(bl_path)) if bl_path.exists() else {}

    out: dict = {}
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.1))

    real_stack_all, null_stack_all = [], []

    for name in frames.CLIPS:
        mp = RESULTS / f"{name}_meta.json"
        fp = RESULTS / f"{name}_frames.csv"
        if not (mp.exists() and fp.exists()):
            print(f"skip {name}: pipeline output missing")
            continue
        with open(mp) as fh:
            work = json.load(fh)["working_body_px"]
        # correlation_function is called with body_px = work, so every xi it
        # returns is in *working* body lengths. Convert to the measured BL that
        # every other result is expressed in; leaving the two in different units
        # would put a wrong number in the report while every ratio still looked
        # right, since ratios are unaffected by the conversion.
        bl_meas = bl_summary.get(name, {}).get("body_px", work)
        to_measured_bl = work / bl_meas
        clip = frames.load_clip(name)
        fr = pd.read_csv(fp)
        present = fr.loc[fr["present"], "frame"].to_numpy()
        present = present[present < len(clip) - 1]
        if len(present) < 5:
            continue
        picks = rng.choice(present, size=min(N_FRAMES, len(present)), replace=False)

        real_curves, null_curves = [], []
        real_xi, null_xi = [], []
        for i in sorted(picks):
            got = measure_frame(clip, int(i), work)
            if got is None:
                continue
            real, shuf = got
            real_curves.append(interp(real))
            real_xi.append(flow.correlation_decay_length(real["r_bl"], real["C"]))
            for c in shuf:
                null_curves.append(interp(c))
                null_xi.append(flow.correlation_decay_length(c["r_bl"], c["C"]))

        if not real_curves:
            continue

        real_mean = np.nanmean(np.array(real_curves), axis=0)
        null_mean = np.nanmean(np.array(null_curves), axis=0)
        xi_real = flow.correlation_decay_length(GRID, real_mean) * to_measured_bl
        xi_null = flow.correlation_decay_length(GRID, null_mean) * to_measured_bl

        real_stack_all.extend(real_curves)
        null_stack_all.extend(null_curves)

        rx = np.array(real_xi, float)
        nx = np.array(null_xi, float)
        rx, nx = rx[np.isfinite(rx)] * to_measured_bl, nx[np.isfinite(nx)] * to_measured_bl
        # One-sided permutation p-value on the per-frame xi distributions.
        p = float((nx >= np.median(rx)).mean()) if len(nx) and len(rx) else np.nan

        out[name] = {
            "n_frames": len(real_curves),
            "n_shuffles": len(null_curves),
            "xi_bl_observed": float(xi_real),
            "xi_bl_shuffled": float(xi_null),
            "ratio_observed_over_shuffled": float(xi_real / xi_null)
            if np.isfinite(xi_null) and xi_null > 0 else np.nan,
            "per_frame_xi_median": float(np.median(rx)) if len(rx) else np.nan,
            "shuffled_xi_median": float(np.median(nx)) if len(nx) else np.nan,
            "shuffled_xi_p95": float(np.percentile(nx, 95)) if len(nx) else np.nan,
            "permutation_p": p,
            "flow_floor_bl": 25.0 / bl_meas,
            "body_px_measured": bl_meas,
        }
        print(f"{name}: {out[name]}", flush=True)

        c = plotstyle.CLIP_COLOUR[name]
        axes[0].plot(GRID, real_mean, color=c, label=name.replace("IMG_", "clip "))
        axes[0].plot(GRID, null_mean, color=c, ls=":", lw=1.0)

    if real_stack_all:
        axes[0].axhline(0, color="0.4", lw=0.7, ls=":")
        axes[0].set_xlabel(r"separation  $r$  (BL)")
        axes[0].set_ylabel(r"$C(r)$")
        axes[0].legend(loc="upper right")
        plotstyle.panel_label(axes[0], "a")

        obs = [v["xi_bl_observed"] for v in out.values()]
        nul = [v["xi_bl_shuffled"] for v in out.values()]
        xpos = np.arange(len(out))
        axes[1].bar(xpos - 0.2, obs, width=0.38, color="#0072B2")
        axes[1].bar(xpos + 0.2, nul, width=0.38, color="0.7")
        axes[1].set_xticks(xpos)
        axes[1].set_xticklabels([k.replace("IMG_0", "") for k in out])
        axes[1].set_ylabel(r"$\xi$  (BL)")
        plotstyle.panel_label(axes[1], "b")

    fig.tight_layout()
    fig.savefig(FIGS / "fig8_correlation_null.png")
    plt.close(fig)

    with open(RESULTS / "correlation_null.json", "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    print("\nwrote", RESULTS / "correlation_null.json")


if __name__ == "__main__":
    main()
