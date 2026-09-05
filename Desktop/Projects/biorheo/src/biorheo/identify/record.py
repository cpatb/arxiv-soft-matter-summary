"""Attach an identifiability block to a registry record, as a sidecar file.

The record YAML may not contain ``identifiability`` (the loader rejects it), because a
hand-typed verdict is a claim without a fit behind it. This module writes
``<record>.identifiability.yaml`` next to the record from an :class:`IdentifiedModel`,
and the loader picks the sidecar up when it carries the ``written_by`` marker.

Before writing, the record's own stored parameters must fall inside the identify
layer's profile intervals; otherwise the block would describe some other fit than the
numbers the record quotes, and the function refuses.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import yaml

from ..registry.schema import FluidRecord
from .fit import IdentifiedModel

__all__ = ["WRITTEN_BY", "sidecar_path", "write_identifiability"]

WRITTEN_BY = "biorheo.identify"


def sidecar_path(record_path: str | Path) -> Path:
    p = Path(record_path)
    return p.with_name(p.stem + ".identifiability.yaml")


def write_identifiability(record: FluidRecord, ident: IdentifiedModel,
                          data_note: str = "") -> Path:
    """Write the sidecar for ``record`` from ``ident``; returns the path written."""
    if record.path is None:
        raise ValueError(f"{record.id}: record has no path to write a sidecar beside")
    if ident.model.name != record.model.name:
        raise ValueError(f"{record.id}: record model {record.model.name!r} but the fit is "
                         f"{ident.model.name!r}")
    mismatched = []
    for name in record.model.param_names:
        est, stored = ident.params[name], record.parameters[name].si
        if not (est.lo <= stored <= est.hi):
            mismatched.append(f"{name}: record {stored:g} outside fit interval "
                              f"[{est.lo:g}, {est.hi:g}]")
    if mismatched:
        raise ValueError(f"{record.id}: the record's parameters are not the fit's: "
                         + "; ".join(mismatched))

    block = {
        "written_by": WRITTEN_BY,
        "written_on": date.today().isoformat(),
        "record_id": record.id,
        "model": ident.model.name,
        "level": ident.level,
        "fit": {"n_points": int(ident.fit.n_points),
                "rmse_log_eta": float(ident.fit.rmse_log),
                "gammadot_min": ident.envelope.gammadot_min,
                "gammadot_max": ident.envelope.gammadot_max},
        "parameters": {
            name: {"value": est.value, "lo": est.lo, "hi": est.hi, "verdict": est.verdict,
                   "decades": (None if est.decades == float("inf") else est.decades),
                   "partner": est.partner, "lo_at_bound": est.lo_at_bound,
                   "hi_at_bound": est.hi_at_bound}
            for name, est in ident.params.items()},
    }
    if data_note:
        block["data_note"] = data_note
    out = sidecar_path(record.path)
    out.write_text("# written by biorheo.identify; do not edit by hand\n"
                   + yaml.safe_dump(block, sort_keys=False), encoding="utf-8", newline="\n")
    return out
