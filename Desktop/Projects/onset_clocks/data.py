"""The 75 impacts, loaded from the drop-impact repo's own per-drop table.

Source
------
``Github/Drop Impact/analysis_2026-07-25/c3_master_curve_per_drop.csv``, the
table behind the C3 master-curve analysis. One row per drop, with the arrested
footprint ``beta_p`` (side-view baseline width divided by the drop diameter),
the relaxed footprint ``beta_eq``, the measured impact speed, the measured drop
diameter, and the two knobs the campaign varied: particle volume fraction and
LiCl concentration.

``beta_p`` is a *contact* width, not a silhouette diameter: it starts near zero
rather than at one, and an arrested drop holds a footprint smaller than its own
diameter. Values here run 0.31 to 1.96.

Arrest is recomputed, not read from a flag
------------------------------------------
The repo's own criterion is ``beta_p < 0.9 beta_eq``
(``analysis_2026-07-25/README.md``). It is applied here rather than trusting any
stored boolean, because a stored ``arrested`` column in this programme has
already once outlived the classification that set it -- it stayed True for 22
drops after the rule moved. A criterion in code cannot go stale silently.

What is deliberately NOT used
-----------------------------
``t_relax_s`` and everything derived from it. The un-jamming clock is a
different question with its own extensive analysis in that repo, and mixing it
in here would make the response variable ambiguous. The response is the
arrested footprint and nothing else.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from clocks import Drop

TABLE = Path(r"C:\Users\compl\Github\Drop Impact\analysis_2026-07-25"
             r"\c3_master_curve_per_drop.csv")

# The repo's own arrest criterion, from analysis_2026-07-25/README.md.
ARREST_RATIO = 0.9

# Rows whose beta_p is a bound rather than a point are not measurements of the
# footprint and are excluded, with the count reported.
POINT_KINDS = ("point",)


def load_drops(path: Path = TABLE, points_only: bool = True) -> tuple[list[Drop], dict]:
    """Return the drops and a provenance dict recording every exclusion.

    Nothing is dropped silently: the returned dict carries N before and after
    each filter, which is the rule this lab applies to every join and mask.
    """
    df = pd.read_csv(path)
    prov: dict = {"source": str(path), "n_rows": len(df)}

    needed = ["phi", "salt_mM", "u0", "d0_mm", "beta_p", "beta_eq_used",
              "beta_p_kind", "label"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name}: missing column(s) {missing}")

    finite = df[needed[:-2]].apply(pd.to_numeric, errors="coerce").notna().all(axis=1)
    prov["n_non_finite_dropped"] = int((~finite).sum())
    df = df[finite]

    if points_only:
        keep = df["beta_p_kind"].isin(POINT_KINDS)
        prov["n_bounded_beta_p_dropped"] = int((~keep).sum())
        prov["bounded_kinds"] = sorted(df.loc[~keep, "beta_p_kind"].unique().tolist())
        df = df[keep]
    prov["n_used"] = len(df)

    drops = [
        Drop(label=f"{r.label} {r.salt_mM:g}mM u{r.u0:.2f}",
             phi=float(r.phi), salt_mM=float(r.salt_mM), u0=float(r.u0),
             d0=float(r.d0_mm) * 1e-3, beta_p=float(r.beta_p),
             beta_eq=float(r.beta_eq_used))
        for r in df.itertuples()
    ]
    prov["phi_levels"] = sorted({d.phi for d in drops})
    prov["salt_levels"] = sorted({d.salt_mM for d in drops})
    prov["u0_range"] = (min(d.u0 for d in drops), max(d.u0 for d in drops))
    prov["d0_range_mm"] = (min(d.d0 for d in drops) * 1e3,
                           max(d.d0 for d in drops) * 1e3)
    prov["beta_p_range"] = (min(d.beta_p for d in drops), max(d.beta_p for d in drops))
    prov["beta_p_log_range"] = float(np.log(prov["beta_p_range"][1]
                                            / prov["beta_p_range"][0]))
    prov["n_arrested"] = sum(1 for d in drops if is_arrested(d))
    return drops, prov


def is_arrested(d: Drop) -> bool:
    """The repo's own criterion, recomputed rather than read from a flag."""
    return bool(d.beta_p < ARREST_RATIO * d.beta_eq)
