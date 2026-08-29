"""
Wrinkle visualization: time series, distributions, orientations, superimposed ROIs.

Usage:
    python visualize.py <root_folder>

Reads results.csv produced by analyze.py.
Writes figures to <root_folder>/figures/.
For aligned superimposition, run align.py on a face folder first to produce offsets.json.
"""

import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import tifffile
from pathlib import Path


# ── data helpers ──────────────────────────────────────────────────────────────

def load_results(root: Path) -> pd.DataFrame:
    df = pd.read_csv(root / "results.csv")
    missing = [c for c in ("image_h", "image_w", "n_forks", "roi_file") if c not in df.columns]
    if missing:
        sys.exit(f"results.csv is missing columns {missing} — re-run analyze.py first.")
    return df


def per_image(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse per-segment rows to one row per image."""
    key = ["condition", "sample", "face", "timepoint_min"]

    def _agg(g):
        area = int(g["image_h"].iloc[0]) * int(g["image_w"].iloc[0])
        return pd.Series({
            "n_segments":      len(g),
            "total_length_px": g["length_px"].sum(),
            "mean_length_px":  g["length_px"].mean(),
            "std_length_px":   g["length_px"].std(),
            "areal_density":   g["length_px"].sum() / area,
            "n_forks":         int(g["n_forks"].iloc[0]),
            "roi_file":        g["roi_file"].iloc[0],
        })

    return df.groupby(key).apply(_agg).reset_index()


# ── 1. density + fork count over time ─────────────────────────────────────────

def plot_timeseries(summary: pd.DataFrame, out: Path) -> None:
    tps = sorted(summary["timepoint_min"].unique())
    groups = list(summary.groupby(["condition", "face"]))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("Wrinkle evolution over time")
    colors = cm.tab10(np.linspace(0, 1, max(len(groups), 1)))

    for (cond, face), grp in groups:
        agg = (
            grp.groupby("timepoint_min")
            .agg(
                density_mean=("areal_density", "mean"),
                density_std=("areal_density", "std"),
                forks_mean=("n_forks", "mean"),
                forks_std=("n_forks", "std"),
            )
            .reindex(tps)
        )
        lbl = f"{cond} / {face}"
        axes[0].errorbar(agg.index, agg["density_mean"], yerr=agg["density_std"],
                         marker="o", capsize=3, label=lbl)
        axes[1].errorbar(agg.index, agg["forks_mean"], yerr=agg["forks_std"],
                         marker="o", capsize=3, label=lbl)

    for ax, title, ylabel in zip(
        axes,
        ["Areal density over time", "Fork count over time"],
        ["Total skeleton px / image px", "Branch points per image"],
    ):
        ax.set_xlabel("Time (min)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(out / "timeseries.png", dpi=150)
    plt.close(fig)
    print("Saved timeseries.png")


# ── 2. length histograms per timepoint ────────────────────────────────────────

def plot_length_histograms(df: pd.DataFrame, out: Path) -> None:
    tps = sorted(df["timepoint_min"].unique())
    n = len(tps)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4), sharey=True)
    if n == 1:
        axes = [axes]
    fig.suptitle("Wrinkle length distribution per timepoint")

    for ax, tp in zip(axes, tps):
        lengths = df[df["timepoint_min"] == tp]["length_px"]
        ax.hist(lengths, bins=30, edgecolor="white", linewidth=0.3, color="steelblue")
        ax.axvline(lengths.median(), color="tomato", linestyle="--", linewidth=1.2,
                   label=f"median={lengths.median():.0f}")
        ax.set_title(f"{tp} min  (n={len(lengths)})")
        ax.set_xlabel("Length (px)")
        ax.legend(fontsize=7)
        if ax is axes[0]:
            ax.set_ylabel("Count")

    fig.tight_layout()
    fig.savefig(out / "length_histograms.png", dpi=150)
    plt.close(fig)
    print("Saved length_histograms.png")


# ── 3. length violin over time ────────────────────────────────────────────────

def plot_length_violin(df: pd.DataFrame, out: Path) -> None:
    tps = sorted(df["timepoint_min"].unique())
    data = [df[df["timepoint_min"] == tp]["length_px"].values for tp in tps]

    fig, ax = plt.subplots(figsize=(max(6, len(tps) * 1.5), 5))
    parts = ax.violinplot(data, positions=range(len(tps)), showmedians=True, showextrema=True)
    for pc in parts["bodies"]:
        pc.set_facecolor("steelblue")
        pc.set_alpha(0.6)
    ax.set_xticks(range(len(tps)))
    ax.set_xticklabels([f"{t} min" for t in tps])
    ax.set_xlabel("Timepoint")
    ax.set_ylabel("Segment length (px)")
    ax.set_title("Wrinkle length distribution over time")
    fig.tight_layout()
    fig.savefig(out / "length_violin.png", dpi=150)
    plt.close(fig)
    print("Saved length_violin.png")


# ── 4. orientation rose plots ─────────────────────────────────────────────────

def plot_orientation_roses(df: pd.DataFrame, out: Path) -> None:
    tps = sorted(df["timepoint_min"].unique())
    n = len(tps)
    fig, axes = plt.subplots(1, n, figsize=(3 * n, 3),
                             subplot_kw={"projection": "polar"})
    if n == 1:
        axes = [axes]
    fig.suptitle("Wrinkle orientation distribution", y=1.02)

    bins = np.linspace(0, np.pi, 19)          # 18 × 10° bins
    centers = (bins[:-1] + bins[1:]) / 2

    for ax, tp in zip(axes, tps):
        angles = np.deg2rad(df[df["timepoint_min"] == tp]["orientation_deg"])
        counts, _ = np.histogram(angles % np.pi, bins=bins)
        ax.bar(centers,       counts, width=np.pi / 18, alpha=0.75, color="steelblue")
        ax.bar(centers + np.pi, counts, width=np.pi / 18, alpha=0.75, color="steelblue")
        ax.set_title(f"{tp} min", pad=12)
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_yticklabels([])

    fig.tight_layout()
    fig.savefig(out / "orientation_roses.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved orientation_roses.png")


# ── 5. density vs fork count scatter ─────────────────────────────────────────

def plot_density_vs_forks(summary: pd.DataFrame, out: Path) -> None:
    tps = sorted(summary["timepoint_min"].unique())
    norm = plt.Normalize(min(tps), max(tps))

    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter(
        summary["areal_density"], summary["n_forks"],
        c=summary["timepoint_min"], cmap="viridis", norm=norm,
        s=50, alpha=0.75, edgecolors="none",
    )
    plt.colorbar(sc, ax=ax, label="Time (min)")
    ax.set_xlabel("Areal density (skeleton px / image px)")
    ax.set_ylabel("Fork count")
    ax.set_title("Wrinkle complexity: density vs branching")
    fig.tight_layout()
    fig.savefig(out / "density_vs_forks.png", dpi=150)
    plt.close(fig)
    print("Saved density_vs_forks.png")


# ── 6. top vs bottom wrinkle counts per sample ────────────────────────────────

def plot_topbottom_counts(summary: pd.DataFrame, out: Path) -> None:
    """One panel per sample: integer wrinkle count vs time for top and bottom."""
    samples = list(summary.groupby(["condition", "sample"]))
    n = len(samples)
    ncols = min(2, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4 * nrows),
                             squeeze=False)

    for ax, ((cond, samp), grp) in zip(axes.flat, samples):
        for face, sub in grp.groupby("face"):
            sub = sub.sort_values("timepoint_min")
            ax.plot(sub["timepoint_min"], sub["n_segments"],
                    marker="o", label=str(face))
        ax.set_title(f"{cond} / sample {samp}")
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Wrinkle count")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    for ax in axes.flat[n:]:
        ax.axis("off")

    fig.suptitle("Wrinkle count over time: top vs bottom")
    fig.tight_layout()
    fig.savefig(out / "topbottom_counts.png", dpi=150)
    plt.close(fig)
    print("Saved topbottom_counts.png")


# ── 7. superimposed ROIs ──────────────────────────────────────────────────────

def _load_offsets(face_dir: Path) -> dict:
    p = face_dir / "offsets.json"
    if p.exists():
        with open(p) as f:
            return json.load(f)
    return {}


def plot_superimposed_rois(df: pd.DataFrame, root: Path, out: Path) -> None:
    summary = per_image(df)

    for (cond, samp, face), grp in summary.groupby(["condition", "sample", "face"]):
        grp = grp.sort_values("timepoint_min")
        tps   = grp["timepoint_min"].tolist()
        files = [root / f for f in grp["roi_file"].tolist()]

        first = tifffile.imread(files[0])
        h, w  = first.shape[-2], first.shape[-1]

        face_dir = files[0].parent
        offsets  = _load_offsets(face_dir)

        cmap_colors = cm.plasma(np.linspace(0.1, 0.9, len(tps)))
        composite   = np.zeros((h, w, 3), dtype=np.float32)

        for roi_path, color in zip(files, cmap_colors):
            raw  = tifffile.imread(roi_path)
            mask = (raw > 0).astype(np.float32)
            off  = offsets.get(roi_path.name, {"dx": 0, "dy": 0})
            dy, dx = int(round(off["dy"])), int(round(off["dx"]))

            if dy or dx:
                shifted = np.zeros_like(mask)
                src_y = slice(max(0, -dy), min(h, h - dy))
                dst_y = slice(max(0,  dy), min(h, h + dy))
                src_x = slice(max(0, -dx), min(w, w - dx))
                dst_x = slice(max(0,  dx), min(w, w + dx))
                shifted[dst_y, dst_x] = mask[src_y, src_x]
                mask = shifted

            for c in range(3):
                composite[:, :, c] += mask * color[c]

        composite = np.clip(composite, 0, 1)

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.imshow(composite, origin="upper")
        ax.set_title(f"{cond} / sample {samp} / {face}")
        ax.axis("off")
        sm = plt.cm.ScalarMappable(cmap="plasma",
                                   norm=plt.Normalize(min(tps), max(tps)))
        plt.colorbar(sm, ax=ax, label="Time (min)", shrink=0.5)
        fig.tight_layout()

        fname = f"roi_superimposed_{cond}_s{samp}_{face}.png"
        fig.savefig(out / fname, dpi=150)
        plt.close(fig)
        print(f"Saved {fname}")


# ── main ──────────────────────────────────────────────────────────────────────

def main(root_folder: str) -> None:
    root = Path(root_folder)
    out  = root / "figures"
    out.mkdir(exist_ok=True)

    df      = load_results(root)
    summary = per_image(df)

    plot_timeseries(summary, out)
    plot_length_histograms(df, out)
    plot_length_violin(df, out)
    plot_orientation_roses(df, out)
    plot_density_vs_forks(summary, out)
    plot_topbottom_counts(summary, out)
    plot_superimposed_rois(df, root, out)

    print(f"\nAll figures saved to {out}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: python visualize.py <root_folder>")
    main(sys.argv[1])
