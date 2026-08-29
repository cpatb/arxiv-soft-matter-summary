"""Minimal plotting style: let the data carry the figure.

Axis labels and a scale, nothing else. No titles, no annotation boxes, no
explanatory text printed onto the axes -- the explanation belongs in the
accompanying markdown, where it can be read, revised and disagreed with.
"""

from __future__ import annotations

import matplotlib as mpl

#: Colour-blind-safe qualitative sequence (Okabe-Ito).
PALETTE = [
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#009E73",  # bluish green
    "#CC79A7",  # reddish purple
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#000000",  # black
]

CLIP_COLOUR = {
    "IMG_0544": PALETTE[0],
    "IMG_0545": PALETTE[1],
    "IMG_0546": PALETTE[2],
    "IMG_0547": PALETTE[3],
}


def apply() -> None:
    mpl.rcParams.update(
        {
            "figure.dpi": 130,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "legend.frameon": False,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": False,
            "ytick.right": False,
            "lines.linewidth": 1.3,
            "axes.prop_cycle": mpl.cycler(color=PALETTE),
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def panel_label(ax, letter: str) -> None:
    """Single bold letter, top-left, outside the data area."""
    ax.text(
        -0.16, 1.04, letter, transform=ax.transAxes,
        fontsize=10, fontweight="bold", va="bottom", ha="left",
    )
