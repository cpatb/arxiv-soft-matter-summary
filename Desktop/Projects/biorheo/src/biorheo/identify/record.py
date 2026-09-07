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

import math
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



_EXPOSURE_NOTE = (
    "exposure_hi / exposure_lo = (span/5000) / (wald_stderr/4) with span = |bound - value| "
    "for that end's box bound. It is how far dst_fit.profiling's step floor overshot the "
    "Wald-scaled step on the walk that produced this interval: > 1 means the box floor set "
    "the first step, and the magnitude is by how much. Written from 2026-09-06 so a reader "
    "can ask whether an interval came from the defective walk without refitting the curve "
    "(NEXT_STEPS items 13 and 16). Note the two ends differ, often by decades."
)

_STEPS_PER_SPAN = 5000.0
_WALD_FRACTION = 0.25


def _exposure(model, name: str, est) -> dict[str, float | None]:
    """The per-END exposure ratio for one parameter; ``None`` where it is undefined.

    Undefined means the Wald standard error is zero or non-finite -- a flat or degenerate
    direction -- in which case the floor won by construction and no finite ratio describes
    it. ``None`` rather than ``inf`` because YAML round-trips ``.inf`` inconsistently and a
    reader should see an absence, not a number.
    """
    i = model.param_names.index(name)
    wald = float(est.wald_stderr)
    if not math.isfinite(wald) or wald <= 0.0:
        return {"exposure_hi": None, "exposure_lo": None}
    step_wald = _WALD_FRACTION * wald
    out: dict[str, float | None] = {}
    for key, bound in (("exposure_hi", float(model.bounds[1][i])),
                       ("exposure_lo", float(model.bounds[0][i]))):
        out[key] = float(abs(bound - est.value) / _STEPS_PER_SPAN / step_wald)
    return out


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
        "residual_structure": (None if ident.residual_structure is None
                               else ident.residual_structure.as_dict()),
        "parameters": {
            name: {"value": est.value, "lo": est.lo, "hi": est.hi, "verdict": est.verdict,
                   "decades": (None if est.decades == float("inf") else est.decades),
                   "partner": est.partner, "lo_at_bound": est.lo_at_bound,
                   "hi_at_bound": est.hi_at_bound,
                   "wald_stderr": est.wald_stderr,
                   **_exposure(ident.model, name, est)}
            for name, est in ident.params.items()},
        "exposure_note": _EXPOSURE_NOTE,
    }
    if data_note:
        block["data_note"] = data_note
    out = sidecar_path(record.path)
    out.write_text("# written by biorheo.identify; do not edit by hand\n"
                   + yaml.safe_dump(block, sort_keys=False), encoding="utf-8", newline="\n")
    return out
