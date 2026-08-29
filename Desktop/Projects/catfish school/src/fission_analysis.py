"""Splitting and rejoining: how a fry ball comes apart and puts itself back together.

Connected-component counting says *whether* the school is in one piece. It says
nothing about the approach to a split, which is the interesting part. A school
does not separate instantaneously: it elongates, develops a waist, and the waist
thins until it breaks. Two envelope descriptors track that without needing the
split to complete:

* **solidity**, the occupied area divided by its convex hull. A waisted outline
  leaves a large crescent of hull unoccupied, so solidity falls as the neck
  forms and recovers when it heals.
* **largest_fraction**, the share of total school area in the biggest piece. It
  is 1 while intact and falls towards 0.5 as a school divides in two.

The script builds event-triggered averages of both around every detected
fission, measures how long separations last, and measures how evenly the school
divides.
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

from catfish import frames, plotstyle  # noqa: E402

RESULTS = frames.ROOT / "results"
FIGS = frames.ROOT / "figures"

#: Window either side of an event for the triggered average, in frames.
WINDOW = 45

#: Separations shorter than this are treated as boundary flicker rather than
#: real fission -- a one-frame break in a thin neck is a segmentation event, not
#: a behavioural one.
MIN_SPLIT_FRAMES = 5


def separation_runs(d: pd.DataFrame) -> list[dict]:
    """Contiguous runs of frames with more than one component, with direction.

    A multi-component episode is not necessarily a split. In IMG_0546 the
    longest one begins at the very first frame with two schools already apart
    and ends when they merge -- a fusion. Direction is read from the ends of the
    run rather than assumed, because calling every separation a "split" would
    describe a merger backwards.
    """
    fr = d["frame"].to_numpy()
    sep = (d["n_components"] > 1).to_numpy()
    runs, start = [], None
    for k, s in enumerate(sep):
        if s and start is None:
            start = k
        elif not s and start is not None:
            runs.append((start, k - 1))
            start = None
    if start is not None:
        runs.append((start, len(sep) - 1))

    out = []
    for i0, i1 in runs:
        before_one = i0 > 0 and not sep[i0 - 1]
        after_one = i1 < len(sep) - 1 and not sep[i1 + 1]
        if before_one and after_one:
            kind = "fission_then_fusion"
        elif before_one:
            kind = "fission"
        elif after_one:
            kind = "fusion"
        else:
            kind = "already_separated"
        out.append(
            {
                "i0": i0, "i1": i1,
                "f0": int(fr[i0]), "f1": int(fr[i1]),
                "n_frames": i1 - i0 + 1,
                "kind": kind,
                # Only a run that begins with a single school has a genuine
                # onset to trigger the event-average on.
                "has_onset": before_one,
            }
        )
    return out


def main() -> None:
    plotstyle.apply()
    FIGS.mkdir(parents=True, exist_ok=True)

    out: dict = {}
    trig_sol, trig_frac = [], []

    for name in frames.CLIPS:
        p = RESULTS / f"{name}_frames.csv"
        if not p.exists():
            continue
        df = pd.read_csv(p)
        d = df[df["present"]].reset_index(drop=True)
        if d.empty or "n_components" not in d:
            continue
        with open(RESULTS / f"{name}_meta.json") as fh:
            fps = json.load(fh)["fps"]

        runs = separation_runs(d)
        real = [r for r in runs if r["n_frames"] >= MIN_SPLIT_FRAMES]
        durations = [r["n_frames"] / fps for r in real]

        # Event-triggered average, only for episodes that begin from a single
        # school: a run already separated at the first frame has no onset, and
        # including it would average over an arbitrary alignment point.
        by_frame = d.set_index("frame")
        for r in real:
            if not r["has_onset"]:
                continue
            idx = np.arange(r["f0"] - WINDOW, r["f0"] + WINDOW + 1)
            sol = np.array([by_frame["solidity"].get(i, np.nan) for i in idx])
            frac = np.array([by_frame["largest_fraction"].get(i, np.nan) for i in idx])
            if np.isfinite(sol).sum() > WINDOW:
                trig_sol.append(sol)
                trig_frac.append(frac)

        # How evenly does it divide? largest_fraction at the deepest point of
        # each separation: 0.5 = an even split, ~0.95 = a few strays shed.
        evenness = []
        for r in real:
            seg = d[(d["frame"] >= r["f0"]) & (d["frame"] <= r["f1"])]
            if len(seg):
                evenness.append(float(seg["largest_fraction"].min()))

        kinds: dict[str, int] = {}
        for r in real:
            kinds[r["kind"]] = kinds.get(r["kind"], 0) + 1

        out[name] = {
            "n_runs_raw": len(runs),
            "n_episodes_real": len(real),
            "episode_kinds": kinds,
            "frac_frames_separated": float((d["n_components"] > 1).mean()),
            "episode_duration_s_median": float(np.median(durations)) if durations else np.nan,
            "episode_duration_s_range": [float(np.min(durations)), float(np.max(durations))]
            if durations else None,
            "fraction_ending_in_reunion": float(
                sum(1 for r in real if r["kind"] in ("fusion", "fission_then_fusion"))
                / len(real)
            ) if real else np.nan,
            "largest_fraction_min_median": float(np.median(evenness)) if evenness else np.nan,
            "max_components": int(d["n_components"].max()),
            "solidity_median": float(d["solidity"].median()),
        }
        print(f"{name}: {out[name]}", flush=True)

    # --- figure ----------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.0))

    if trig_sol:
        t = np.arange(-WINDOW, WINDOW + 1) / 29.97
        S = np.array(trig_sol)
        F = np.array(trig_frac)
        for ax, A, lab in [(axes[0], S, "solidity"),
                           (axes[1], F, "largest fraction")]:
            med = np.nanmedian(A, axis=0)
            q1 = np.nanpercentile(A, 25, axis=0)
            q3 = np.nanpercentile(A, 75, axis=0)
            ax.plot(t, med, color="#0072B2")
            ax.fill_between(t, q1, q3, color="#0072B2", alpha=0.18, linewidth=0)
            ax.axvline(0, color="0.45", lw=0.9, ls="--")
            ax.set_xlabel("time from split (s)")
            ax.set_ylabel(lab)
        axes[0].text(0.02, 0.04, f"n = {len(S)}", transform=axes[0].transAxes,
                     fontsize=7, color="0.35")

    for name in out:
        p = RESULTS / f"{name}_frames.csv"
        d = pd.read_csv(p)
        d = d[d["present"]]
        axes[2].plot(d["t_s"], d["solidity"], lw=0.8, alpha=0.85,
                     color=plotstyle.CLIP_COLOUR[name],
                     label=name.replace("IMG_", "clip "))
    axes[2].set_xlabel("time (s)")
    axes[2].set_ylabel("solidity")
    axes[2].legend(loc="lower left", ncol=2)

    for ax, L in zip(axes, "abc"):
        plotstyle.panel_label(ax, L)
    fig.tight_layout()
    fig.savefig(FIGS / "fig10_fission.png")
    plt.close(fig)

    # Does solidity actually fall before a split? Compare the pre-event baseline
    # against the value at the moment of separation.
    if trig_sol:
        S = np.array(trig_sol)
        base = np.nanmedian(S[:, : WINDOW // 2], axis=1)   # earliest quarter
        atsplit = np.nanmedian(S[:, WINDOW - 3: WINDOW + 1], axis=1)
        ok = np.isfinite(base) & np.isfinite(atsplit)
        if ok.sum() >= 5:
            drop = base[ok] - atsplit[ok]
            from scipy import stats

            w = stats.wilcoxon(base[ok], atsplit[ok]) if ok.sum() >= 6 else None
            out["_necking"] = {
                "n_events": int(ok.sum()),
                "solidity_baseline_median": float(np.median(base[ok])),
                "solidity_at_split_median": float(np.median(atsplit[ok])),
                "median_drop": float(np.median(drop)),
                "frac_events_with_drop": float((drop > 0).mean()),
                "wilcoxon_p": float(w.pvalue) if w is not None else None,
            }
            print("\nnecking test:", json.dumps(out["_necking"], indent=2))

    with open(RESULTS / "fission_analysis.json", "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    print("\nwrote", RESULTS / "fission_analysis.json")


if __name__ == "__main__":
    main()
