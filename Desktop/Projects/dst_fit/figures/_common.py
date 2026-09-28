"""Shared house style and provenance for the strain-clock paper figures.

Spare, PRL-style: axis labels, panel letters, a legend where needed, nothing else.
What a figure means goes in its sibling .md, never on the axes.
"""

from __future__ import annotations

import subprocess
from datetime import date
from pathlib import Path

import matplotlib as mpl

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
DATA_DIR = FIG_DIR / "data"
LIN_ROOT = ROOT / (
    "Hydrodynamic and Contact Contributions to Continuous Shear Thickening "
    "in Colloidal Suspensions"
)

# Okabe-Ito, colourblind-safe
CB = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]

mpl.rcParams.update({
    "font.size": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6.5,
    "legend.frameon": False,
    "axes.linewidth": 0.7,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "lines.linewidth": 1.0,
    "savefig.dpi": 300,
    "figure.dpi": 130,
    "pdf.fonttype": 42,
})


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:  # noqa: BLE001 -- provenance only, never fatal
        return "unknown"


def provenance(script: Path, inputs: str, params: str) -> str:
    # No date/time: the CSV and PDF metadata must regenerate byte-identically.
    return (f"script={script.relative_to(ROOT).as_posix()} commit={git_commit()} "
            f"inputs={inputs} params={params}")


def save(fig, slug: str, prov: str) -> tuple[Path, Path]:
    FIG_DIR.mkdir(exist_ok=True)
    pdf = FIG_DIR / f"fig_{slug}.pdf"
    png = FIG_DIR / f"fig_{slug}.png"
    fig.savefig(pdf, metadata={"Subject": prov, "CreationDate": None, "ModDate": None})
    fig.savefig(png, metadata={"Software": None})
    return pdf, png


def today() -> str:
    return date.today().isoformat()
