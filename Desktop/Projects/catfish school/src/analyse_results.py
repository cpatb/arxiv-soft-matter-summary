"""Normalise the pixel-space pipeline output to body lengths and make figures.

The pipeline stores geometry in pixels. This script builds the body-length
series BL(t) for each clip from the isolated-silhouette samples, converts every
length to body lengths, and produces the figures and the summary table.

Why normalisation happens here and not in the pipeline: BL is itself a
measurement drawn from the same frames, so folding it in during the expensive
pass would bake one particular estimate into the results permanently. Kept
separate, the normalisation can be revised -- or its sensitivity tested --
without recomputing anything.
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

#: Half-width of the rolling window used to smooth the body-length series, in
#: seconds. Body length changes only through camera drift, which is slow; a wide
#: window buys precision without hiding anything real.
BL_WINDOW_S = 4.0

#: Minimum silhouettes inside the window for a BL value to be trusted.
BL_MIN_SAMPLES = 25


def body_length_series(clip_name: str, t_grid: np.ndarray) -> tuple[np.ndarray, dict]:
    """Rolling-median body length in px, evaluated on t_grid.

    Falls back to the clip-wide median where a window is too sparse, so the
    series never contains gaps that would silently drop frames downstream.
    """
    # Prefer the dedicated re-measurement (scale-free criteria). The batch pass
    # collected silhouettes with an acceptance window scaled to the assumed body
    # length, which biased the estimate downward; measure_bodylength.py redoes it.
    v2 = RESULTS / f"{clip_name}_bodylength_v2.csv"
    path = v2 if v2.exists() else RESULTS / f"{clip_name}_bodylength.csv"
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError(f"{path}: no body-length samples")

    t = df["t_s"].to_numpy()
    mj = df["major_px"].to_numpy()
    overall = float(np.median(mj))

    out = np.full(len(t_grid), np.nan)
    n_in = np.zeros(len(t_grid), dtype=int)
    for i, tc in enumerate(t_grid):
        sel = np.abs(t - tc) <= BL_WINDOW_S
        n_in[i] = int(sel.sum())
        if n_in[i] >= BL_MIN_SAMPLES:
            out[i] = float(np.median(mj[sel]))
    n_filled = int(np.isnan(out).sum())
    out = np.where(np.isnan(out), overall, out)

    info = {
        "clip": clip_name,
        "n_samples": int(len(mj)),
        "median_px": overall,
        "iqr_px": [float(np.percentile(mj, 25)), float(np.percentile(mj, 75))],
        "drift_pct": float(100 * (np.nanmax(out) - np.nanmin(out)) / overall),
        "n_frames_using_fallback": n_filled,
    }
    return out, info


def load_clip_results(clip_name: str) -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(RESULTS / f"{clip_name}_frames.csv")
    with open(RESULTS / f"{clip_name}_meta.json") as fh:
        meta = json.load(fh)

    bl, bl_info = body_length_series(clip_name, df["t_s"].to_numpy())
    df["body_px"] = bl

    # --- convert pixels to body lengths -----------------------------------
    df["area_bl2"] = df["area_px2"] / df["body_px"] ** 2
    df["total_area_bl2"] = df["total_area_px2"] / df["body_px"] ** 2
    df["perimeter_bl"] = df["perimeter_px"] / df["body_px"]
    df["gyration_bl"] = df["gyration_px"] / df["body_px"]

    # Equivalent diameter of a disc with the school's area: the single number
    # that best stands for "how big is the school".
    df["size_bl"] = 2.0 * np.sqrt(df["area_bl2"] / np.pi)

    # Speeds: px/frame -> body lengths per second.
    fps = meta["fps"]
    for col in ("mean_speed_px", "median_speed_px"):
        if col in df:
            df[col.replace("_px", "_bl_s")] = df[col] / df["body_px"] * fps

    # xi was computed in units of the working scale, not the measured BL.
    if "xi_bl_work" in df:
        df["xi_bl"] = df["xi_bl_work"] * meta["working_body_px"] / df["body_px"]

    meta["body_length"] = bl_info

    # Persist it beside the working scale. The measured body length is the unit
    # every downstream result is expressed in, and computing it here and keeping
    # it only in memory is what left fringe_analysis.py with no measured length
    # to find -- it silently fell back to the working scale, ~1.4-1.6x larger.
    with open(RESULTS / f"{clip_name}_meta.json", "w") as fh:
        json.dump(meta, fh, indent=2, default=float)

    return df, meta


def summarise(df: pd.DataFrame, meta: dict) -> dict:
    """Clip-level summary over frames where a school is present."""
    d = df[df["present"]].copy() if "present" in df else df.copy()
    if d.empty:
        return {"clip": meta["clip"], "n_present": 0}

    def q(col, f=np.nanmedian):
        return float(f(d[col])) if col in d else float("nan")

    # Camera-motion QA: how much apparent motion was removed, relative to the
    # fish motion that remains. If this is large the kinematics are suspect.
    cam = d["camera_shift_px"] if "camera_shift_px" in d else pd.Series(dtype=float)
    spd = d["median_speed_px"] if "median_speed_px" in d else pd.Series(dtype=float)

    return {
        "clip": meta["clip"],
        "duration_s": round(float(d["t_s"].max() - d["t_s"].min()), 1),
        "n_frames_present": int(len(d)),
        "frac_present": round(float(df["present"].mean()), 3),
        "body_px_median": round(meta["body_length"]["median_px"], 1),
        "body_px_drift_pct": round(meta["body_length"]["drift_pct"], 1),
        "n_bl_samples": meta["body_length"]["n_samples"],
        "size_bl_median": round(q("size_bl"), 1),
        "size_bl_range": [round(float(np.nanpercentile(d["size_bl"], 5)), 1),
                          round(float(np.nanpercentile(d["size_bl"], 95)), 1)],
        "area_bl2_median": round(q("area_bl2"), 1),
        "aspect_ratio_median": round(q("aspect_ratio"), 2),
        "circularity_median": round(q("circularity"), 3),
        "solidity_median": round(q("solidity"), 3),
        "boundary_roughness_median": round(q("boundary_roughness"), 3),
        "packing_fraction_median": round(q("packing_fraction"), 3),
        "n_components_max": int(np.nanmax(d["n_components"])) if "n_components" in d else 0,
        "frac_frames_split": round(float((d["n_components"] > 1).mean()), 3)
        if "n_components" in d else float("nan"),
        "polarisation_median": round(q("polarisation"), 3),
        "speed_bl_s_median": round(q("median_speed_bl_s"), 3),
        "xi_bl_median": round(q("xi_bl"), 2),
        "camera_shift_px_median": round(float(np.nanmedian(cam)), 2) if len(cam) else np.nan,
        "camera_ok_frac": round(float(d["camera_ok"].mean()), 3) if "camera_ok" in d else np.nan,
        "camera_to_fish_motion_ratio": round(
            float(np.nanmedian(cam) / np.nanmedian(spd)), 2
        ) if len(cam) and len(spd) and np.nanmedian(spd) > 0 else np.nan,
    }


# ---------------------------------------------------------------- figures ---


def fig_shape_timeseries(data: dict[str, pd.DataFrame]) -> None:
    fig, axes = plt.subplots(4, 1, figsize=(6.6, 7.4), sharex=True)
    panels = [
        ("size_bl", r"school size  $D$  (BL)"),
        ("aspect_ratio", r"aspect ratio"),
        ("circularity", r"circularity  $4\pi A/P^2$"),
        ("packing_fraction", r"packing fraction  $\phi$"),
    ]
    for ax, (col, lab) in zip(axes, panels):
        for name, df in data.items():
            d = df[df["present"]]
            if col in d:
                ax.plot(d["t_s"], d[col], color=plotstyle.CLIP_COLOUR[name],
                        label=name.replace("IMG_", "clip "), alpha=0.85)
        ax.set_ylabel(lab)
    axes[0].legend(ncol=4, loc="upper right")
    axes[-1].set_xlabel("time (s)")
    for ax, L in zip(axes, "abcd"):
        plotstyle.panel_label(ax, L)
    fig.savefig(FIGS / "fig1_shape_timeseries.png")
    plt.close(fig)


#: Common separation grid, in body lengths, that individual C(r) curves are
#: interpolated onto before averaging.
CORR_GRID_BL = np.linspace(0.0, 25.0, 100)

#: Number of school-size bins used for the scale-free test.
N_SIZE_BINS = 5


def _load_curves(name: str, df: pd.DataFrame, meta_work_px: float) -> list[dict]:
    """Read saved C(r) curves and re-express r in measured body lengths.

    The pipeline computed r in units of the per-clip *working* scale. Converting
    to the measured, time-varying body length has to happen here, where BL(t) is
    known.
    """
    p = RESULTS / f"{name}_correlation.json"
    if not p.exists():
        return []
    with open(p) as fh:
        curves = json.load(fh)

    by_frame = df.set_index("frame")
    out = []
    for c in curves:
        if not c["r_bl"]:
            continue
        f = c["frame"]
        if f not in by_frame.index:
            continue
        row = by_frame.loc[f]
        bl = float(row["body_px"])
        if not np.isfinite(bl) or bl <= 0:
            continue
        scale = meta_work_px / bl
        out.append(
            {
                "frame": f,
                "r_bl": np.asarray(c["r_bl"], float) * scale,
                "C": np.asarray(c["C"], float),
                "size_bl": float(row["size_bl"]),
                "clip": name,
            }
        )
    return out


def _average_C(curves: list[dict]) -> tuple[np.ndarray, np.ndarray, float]:
    """Mean C(r) on the common grid, plus its zero crossing.

    Averaging the curves and *then* taking the crossing is deliberate. The
    single-frame estimator has a spread comparable to its own value (verified in
    test_numerics.py), so averaging per-frame xi would be dominated by noise in
    the tail of that distribution.
    """
    stack = []
    for c in curves:
        if len(c["r_bl"]) < 4:
            continue
        # Only interpolate within each curve's measured range; beyond it, NaN.
        vals = np.interp(CORR_GRID_BL, c["r_bl"], c["C"],
                         left=np.nan, right=np.nan)
        stack.append(vals)
    if not stack:
        return CORR_GRID_BL, np.full_like(CORR_GRID_BL, np.nan), float("nan")
    arr = np.array(stack)
    # Curves reach different maximum separations, so the far tail of the average
    # can rest on one or two frames. Require a quorum before trusting a point.
    n_contrib = np.isfinite(arr).sum(axis=0)
    with np.errstate(invalid="ignore"):
        mean_C = np.nanmean(arr, axis=0)
    mean_C = np.where(n_contrib >= max(5, 0.25 * len(arr)), mean_C, np.nan)
    return CORR_GRID_BL, mean_C, decay_length(CORR_GRID_BL, mean_C)


def decay_length(r: np.ndarray, C: np.ndarray) -> float:
    """Correlation length as the 1/e decay of C(r) -- see catfish.flow."""
    from catfish.flow import correlation_decay_length

    return correlation_decay_length(r, C)


def zero_crossing(r: np.ndarray, C: np.ndarray) -> float:
    from catfish.flow import _first_zero_crossing

    return _first_zero_crossing(r, C)


def fig_correlation(data: dict[str, pd.DataFrame], metas: list[dict]) -> dict:
    """Averaged C(r) per clip, and the correlation-length-versus-size test."""
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1))

    work = {m["clip"]: m["working_body_px"] for m in metas}
    all_curves: list[dict] = []
    per_clip: dict = {}

    for name, df in data.items():
        curves = _load_curves(name, df, work[name])
        if not curves:
            continue
        all_curves.extend(curves)
        r, C, xi = _average_C(curves)
        col = plotstyle.CLIP_COLOUR[name]
        axes[0].plot(r, C, color=col, label=name.replace("IMG_", "clip "))
        # Largest separation actually sampled: xi close to this is a lower
        # bound, not a measurement.
        r_max = float(np.nanmax([c["r_bl"].max() for c in curves]))
        r_reliable = float(np.nanmax(r[np.isfinite(C)])) if np.isfinite(C).any() else np.nan
        zc = zero_crossing(r, C)
        per_clip[name] = {
            "n_curves": len(curves),
            "xi_decay_bl": float(xi),          # 1/e decay length -- the reported xi
            "zero_crossing_bl": float(zc),     # NaN means C stays positive throughout
            "C_stays_positive": bool(not np.isfinite(zc)),
            "r_max_bl": r_max,
            "r_reliable_bl": r_reliable,
            "xi_over_rmax": float(xi / r_reliable)
            if np.isfinite(xi) and r_reliable > 0 else np.nan,
            "saturated": bool(np.isfinite(xi) and np.isfinite(r_reliable)
                              and xi > 0.8 * r_reliable),
            "median_size_bl": float(np.nanmedian([c["size_bl"] for c in curves])),
        }

    axes[0].axhline(0, color="0.4", lw=0.7, ls=":")
    axes[0].set_xlabel(r"separation  $r$  (BL)")
    axes[0].set_ylabel(r"$C(r)$")
    axes[0].legend(loc="upper right")
    plotstyle.panel_label(axes[0], "a")

    # --- scale-free test: bin by school size, average within bin ----------
    stats: dict = {"per_clip": per_clip}
    sizes = np.array([c["size_bl"] for c in all_curves])
    ok = np.isfinite(sizes)
    if ok.sum() > 40:
        edges = np.nanpercentile(sizes[ok], np.linspace(0, 100, N_SIZE_BINS + 1))
        edges = np.unique(edges)
        bin_L, bin_xi, bin_n = [], [], []
        for i in range(len(edges) - 1):
            sel = [c for c in all_curves
                   if edges[i] <= c["size_bl"] <= edges[i + 1]]
            if len(sel) < 10:
                continue
            _, Cb, xib = _average_C(sel)
            if not np.isfinite(xib):
                continue
            bin_L.append(float(np.median([c["size_bl"] for c in sel])))
            bin_xi.append(float(xib))
            bin_n.append(len(sel))

        if len(bin_L) >= 3:
            bin_L = np.array(bin_L)
            bin_xi = np.array(bin_xi)
            axes[1].scatter(bin_L, bin_xi, s=34, color="k", zorder=3)
            slope, intercept = np.polyfit(np.log10(bin_L), np.log10(bin_xi), 1)
            r_p = float(np.corrcoef(np.log10(bin_L), np.log10(bin_xi))[0, 1])
            xs = np.linspace(bin_L.min() * 0.9, bin_L.max() * 1.1, 40)
            axes[1].plot(xs, 10**intercept * xs**slope, color="k", lw=1.1)
            # Reference line of slope 1 (xi proportional to L), anchored at the
            # median, which is what scale-free correlation would look like.
            anchor = np.median(bin_xi / bin_L)
            axes[1].plot(xs, anchor * xs, color="0.55", ls="--", lw=1.0)
            stats.update(
                {
                    "n_size_bins": len(bin_L),
                    "n_curves_per_bin": bin_n,
                    "bin_size_bl": bin_L.tolist(),
                    "bin_xi_bl": bin_xi.tolist(),
                    "loglog_slope": float(slope),
                    "loglog_r": r_p,
                    "xi_over_L_median": float(np.median(bin_xi / bin_L)),
                    "size_range_bl": [float(bin_L.min()), float(bin_L.max())],
                    "size_dynamic_range": float(bin_L.max() / bin_L.min()),
                }
            )
        axes[1].set_xscale("log")
        axes[1].set_yscale("log")
    axes[1].set_xlabel(r"school size  $D$  (BL)")
    axes[1].set_ylabel(r"correlation length  $\xi$  (BL)")
    plotstyle.panel_label(axes[1], "b")

    fig.tight_layout()
    fig.savefig(FIGS / "fig3_correlation.png")
    plt.close(fig)
    return stats


def fig_packing_profile(data: dict[str, pd.DataFrame]) -> dict:
    """Packing fraction from the school's rind inward to its core."""
    fig, ax = plt.subplots(figsize=(3.6, 3.0))
    out = {}
    for name in data:
        p = RESULTS / f"{name}_profiles.json"
        if not p.exists():
            continue
        with open(p) as fh:
            profs = json.load(fh)
        if not profs:
            continue
        # Profiles have different lengths; interpolate each onto a common grid
        # in *relative* depth so schools of different size can be averaged.
        grid = np.linspace(0.0, 1.0, 20)
        stack = []
        for pr in profs:
            d = np.asarray(pr["depth_bl"], dtype=float)
            v = np.asarray(pr["packing"], dtype=float)
            if len(d) < 4 or d.max() <= 0:
                continue
            stack.append(np.interp(grid, d / d.max(), v))
        if not stack:
            continue
        stack = np.array(stack)
        med = np.nanmedian(stack, axis=0)
        lo = np.nanpercentile(stack, 25, axis=0)
        hi = np.nanpercentile(stack, 75, axis=0)
        c = plotstyle.CLIP_COLOUR[name]
        ax.plot(grid, med, color=c, label=name.replace("IMG_", "clip "))
        ax.fill_between(grid, lo, hi, color=c, alpha=0.15, linewidth=0)
        out[name] = {
            "edge_packing": float(med[1]),
            "core_packing": float(med[-1]),
            "core_minus_edge": float(med[-1] - med[1]),
            "n_profiles": int(len(stack)),
        }
    ax.set_xlabel("relative depth from boundary")
    ax.set_ylabel(r"packing fraction  $\phi$")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(FIGS / "fig4_packing_profile.png")
    plt.close(fig)
    return out


def fig_kinematics(data: dict[str, pd.DataFrame]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(8.4, 2.9))
    for name, df in data.items():
        d = df[df["present"]]
        c = plotstyle.CLIP_COLOUR[name]
        if "polarisation" in d:
            v = d["polarisation"].dropna()
            if len(v) > 10:
                axes[0].hist(v, bins=30, range=(0, 1), histtype="step",
                             density=True, color=c, label=name.replace("IMG_", "clip "))
        if "median_speed_bl_s" in d:
            v = d["median_speed_bl_s"].dropna()
            v = v[np.isfinite(v)]
            if len(v) > 10:
                axes[1].hist(v, bins=40, histtype="step", density=True, color=c)
        if "n_components" in d:
            axes[2].plot(d["t_s"], d["n_components"], color=c, lw=0.9, alpha=0.85)

    axes[0].set_xlabel(r"polarisation  $\Phi$")
    axes[0].set_ylabel("density")
    axes[0].legend(loc="upper left")
    axes[1].set_xlabel(r"speed  (BL s$^{-1}$)")
    axes[1].set_ylabel("density")
    axes[2].set_xlabel("time (s)")
    axes[2].set_ylabel("connected components")
    for ax, L in zip(axes, "abc"):
        plotstyle.panel_label(ax, L)
    fig.tight_layout()
    fig.savefig(FIGS / "fig2_kinematics.png")
    plt.close(fig)


def fig_bodylength_qa(metas: list[dict]) -> None:
    """Quality control on the one measurement everything else is divided by."""
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9))
    for meta in metas:
        name = meta["clip"]
        df = pd.read_csv(RESULTS / f"{name}_bodylength.csv")
        if df.empty:
            continue
        c = plotstyle.CLIP_COLOUR[name]
        axes[0].hist(df["major_px"], bins=40, histtype="step", density=True,
                     color=c, label=name.replace("IMG_", "clip "))
        t = np.linspace(df["t_s"].min(), df["t_s"].max(), 120)
        bl, _ = body_length_series(name, t)
        axes[1].plot(t, bl, color=c)
    axes[0].set_xlabel("silhouette major axis (px)")
    axes[0].set_ylabel("density")
    axes[0].legend(loc="upper right")
    axes[1].set_xlabel("time (s)")
    axes[1].set_ylabel("body length (px)")
    for ax, L in zip(axes, "ab"):
        plotstyle.panel_label(ax, L)
    fig.tight_layout()
    fig.savefig(FIGS / "figS1_bodylength_qa.png")
    plt.close(fig)


def main() -> None:
    plotstyle.apply()
    FIGS.mkdir(parents=True, exist_ok=True)

    data, metas = {}, []
    for name in frames.CLIPS:
        if not (RESULTS / f"{name}_frames.csv").exists():
            print(f"skip {name}: no results")
            continue
        df, meta = load_clip_results(name)
        data[name] = df
        metas.append(meta)

    if not data:
        raise SystemExit("no results to analyse")

    summaries = [summarise(data[m["clip"]], m) for m in metas]
    sdf = pd.DataFrame(summaries)
    sdf.to_csv(RESULTS / "summary_by_clip.csv", index=False)

    fig_shape_timeseries(data)
    fig_kinematics(data)
    corr_stats = fig_correlation(data, metas)
    pack_stats = fig_packing_profile(data)
    fig_bodylength_qa(metas)

    events = {}
    for name in data:
        p = RESULTS / f"{name}_events.json"
        if p.exists():
            with open(p) as fh:
                ev = json.load(fh)
            events[name] = {
                "n_fission_frames": sum(1 for e in ev if e["type"] == "fission"),
                "n_fusion_frames": sum(1 for e in ev if e["type"] == "fusion"),
            }

    out = {
        "per_clip": summaries,
        "correlation_scaling": corr_stats,
        "packing_profile": pack_stats,
        "events": events,
    }
    with open(RESULTS / "analysis_summary.json", "w") as fh:
        json.dump(out, fh, indent=2, default=float)

    pd.set_option("display.width", 200, "display.max_columns", 50)
    print(sdf.to_string(index=False))
    print("\ncorrelation scaling:", json.dumps(corr_stats, indent=2, default=float))
    print("packing profile:", json.dumps(pack_stats, indent=2, default=float))
    print("events:", json.dumps(events, indent=2))


if __name__ == "__main__":
    main()
