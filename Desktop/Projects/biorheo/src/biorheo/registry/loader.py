"""Load and validate a fluid record from YAML."""

from __future__ import annotations

import re
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import yaml

from ..constitutive import MODELS, get_model
from ..units import dimensionality, same_dimension, to_si
from .provenance import is_doi_syntax, resolve_doi
from .schema import (
    Envelope,
    EnvelopeError,
    FluidRecord,
    Licence,
    ProvenanceError,
    Source,
    make_parameter,
)

__all__ = ["load_record", "REQUIRED_KEYS", "RECORDS_DIR", "load_all"]

REQUIRED_KEYS = ("id", "fluid", "model", "parameters", "envelope", "sources", "licence")
_OPTIONAL_KEYS = ("notes", "extra")
_FORBIDDEN_KEYS = ("identifiability",)  # only biorheo.identify may write this
RECORDS_DIR = Path(__file__).with_name("records")


def _quantity(name: str, spec: object, expect_dim: str) -> float:
    """SI value of an envelope quantity ``{value, unit}`` with a required dimension."""
    if not isinstance(spec, Mapping) or {"value", "unit"} - set(spec):
        raise EnvelopeError(f"envelope.{name} must be a mapping {{value, unit}}, got {spec!r}")
    unit = str(spec["unit"])
    if dimensionality(unit) != expect_dim:
        raise EnvelopeError(
            f"envelope.{name} has unit {unit!r} of dimension {dimensionality(unit)}, "
            f"expected {expect_dim}")
    return to_si(float(spec["value"]), unit)


def _expected_unit(model_unit: str, si_params: Mapping[str, float]) -> str:
    """Substitute symbolic exponents like ``Pa*s**n`` with the record's own value of n."""
    def sub(m: re.Match) -> str:
        return f"**{si_params[m.group(1)]!r}"
    return re.sub(r"\*\*([A-Za-z_]\w*)", sub, model_unit)


def _check_sources(rec_id: str, sources: object, allow_network: bool) -> tuple[Source, ...]:
    if not isinstance(sources, list) or not sources:
        raise ProvenanceError(f"{rec_id}: 'sources' must be a non-empty list; a record with "
                              "no source is a number from memory")
    out = []
    for i, s in enumerate(sources):
        if not isinstance(s, Mapping) or "doi" not in s:
            raise ProvenanceError(f"{rec_id}: sources[{i}] must be a mapping with a 'doi'")
        doi = str(s["doi"])
        if not is_doi_syntax(doi):
            raise ProvenanceError(f"{rec_id}: sources[{i}].doi {doi!r} is not DOI syntax "
                                  "(bare '10.xxxx/suffix', no prefix)")
        r = resolve_doi(doi, allow_network=allow_network)
        if r is None:
            raise ProvenanceError(
                f"{rec_id}: could not resolve DOI {doi!r} "
                f"({'not in the cache and the network is disallowed' if not allow_network else 'the handle system has no such DOI'})")
        out.append(Source(doi=doi, role=str(s.get("role", "unspecified")),
                          resolved_on=r.resolved_on, label=r.label))
    return tuple(out)


def load_record(path: str | Path, allow_network: bool = False) -> FluidRecord:
    """Parse, validate and return one record; raises on every rule violation.

    ``allow_network=False`` (the default) resolves DOIs from the shipped cache only.
    """
    path = Path(path)
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw, Mapping):
        raise ValueError(f"{path}: top level must be a mapping")

    for k in _FORBIDDEN_KEYS:
        if k in raw:
            raise ProvenanceError(
                f"{path}: '{k}' is present in the file; only biorheo.identify may write it")
    missing = [k for k in REQUIRED_KEYS if k not in raw]
    if missing:
        # name the first missing key so the message says what it is about
        raise ValueError(f"{path}: missing required key(s) {missing}: '{missing[0]}' is required")
    unknown = set(raw) - set(REQUIRED_KEYS) - set(_OPTIONAL_KEYS)
    if unknown:
        raise ValueError(f"{path}: unknown key(s) {sorted(unknown)}")

    rec_id = str(raw["id"])

    # model and parameter set
    model_name = str(raw["model"])
    if model_name not in MODELS:
        raise ValueError(f"{rec_id}: unknown model {model_name!r}; known: {sorted(MODELS)}")
    model = get_model(model_name)
    if not isinstance(raw["parameters"], Mapping):
        raise ValueError(f"{rec_id}: 'parameters' must be a mapping")
    got, want = set(raw["parameters"]), set(model.param_names)
    if got != want:
        raise ValueError(
            f"{rec_id}: parameters {sorted(got)} do not match model {model_name!r}'s "
            f"parameters {sorted(want)}")
    params = {name: make_parameter(name, spec) for name, spec in raw["parameters"].items()}
    si = {name: p.si for name, p in params.items()}
    for name, p in params.items():
        expected = _expected_unit(model.units[name], si)
        if not same_dimension(p.unit, expected):
            raise ValueError(
                f"{rec_id}: parameter {name!r} has unit {p.unit!r} "
                f"({dimensionality(p.unit)}), but {model_name} needs {expected!r} "
                f"({dimensionality(expected)})")
    model.check_params([si[n] for n in model.param_names])

    # envelope
    env = raw["envelope"]
    if not isinstance(env, Mapping):
        raise EnvelopeError(f"{rec_id}: 'envelope' must be a mapping")
    for k in ("gammadot_min", "gammadot_max", "temperature_C", "geometry"):
        if k not in env:
            raise EnvelopeError(f"{rec_id}: envelope is missing {k!r}")
    gmin = _quantity("gammadot_min", env["gammadot_min"], "1 / [time]")
    gmax = _quantity("gammadot_max", env["gammadot_max"], "1 / [time]")
    if not (0 < gmin < gmax):
        raise EnvelopeError(f"{rec_id}: envelope needs 0 < gammadot_min < gammadot_max, "
                            f"got {gmin} and {gmax}")
    envelope = Envelope(gammadot_min=gmin, gammadot_max=gmax,
                        temperature_C=float(env["temperature_C"]), geometry=str(env["geometry"]))

    sources = _check_sources(rec_id, raw["sources"], allow_network)

    lic = raw["licence"]
    if not isinstance(lic, Mapping) or "name" not in lic or "commercial_use" not in lic:
        raise ValueError(f"{rec_id}: licence needs 'name' and an explicit 'commercial_use'")
    if not isinstance(lic["commercial_use"], bool):
        raise ValueError(f"{rec_id}: licence.commercial_use must be true or false, "
                         f"got {lic['commercial_use']!r}")
    licence = Licence(name=str(lic["name"]), commercial_use=bool(lic["commercial_use"]))

    sidecar = _read_sidecar(path, rec_id, model)
    if sidecar is not None:
        _check_sidecar_values(sidecar, np.array([params[n].si for n in model.param_names]),
                              model, path)
    return FluidRecord(
        id=rec_id, fluid=str(raw["fluid"]), model=model, parameters=params,
        envelope=envelope, sources=sources, licence=licence,
        notes=str(raw.get("notes", "")), identifiability=sidecar,
        path=str(path), extra=dict(raw.get("extra") or {}),
    )


_SIDECAR_WRITER = "biorheo.identify"


def _read_sidecar(path: Path, rec_id: str, model) -> dict | None:
    """The identifiability block from ``<record>.identifiability.yaml``, if present.

    Only ``biorheo.identify`` writes that file; a sidecar without its marker, for
    another record id, or for another model is rejected rather than ignored.
    """
    sc = path.with_name(path.stem + ".identifiability.yaml")
    if not sc.exists():
        return None
    block = yaml.safe_load(sc.read_text(encoding="utf-8"))
    if not isinstance(block, Mapping) or block.get("written_by") != _SIDECAR_WRITER:
        raise ProvenanceError(f"{sc}: identifiability sidecar is not marked as written by "
                              f"{_SIDECAR_WRITER}; only that layer may write verdicts")
    if block.get("record_id") != rec_id:
        raise ProvenanceError(f"{sc}: sidecar is for record {block.get('record_id')!r}, "
                              f"not {rec_id!r}")
    if block.get("model") != model.name:
        raise ProvenanceError(f"{sc}: sidecar fitted {block.get('model')!r}, record uses "
                              f"{model.name!r}")
    if set(block.get("parameters", {})) != set(model.param_names):
        raise ProvenanceError(f"{sc}: sidecar parameters do not match the model's")
    return dict(block)


def _check_sidecar_values(sc_block: dict, params_si: np.ndarray, model, path: Path) -> None:
    """Apply at load time the rule ``write_identifiability`` applies at write time: every
    record value must lie inside the sidecar's profile interval. The record may quote
    the source's own fit (the Carbopol records carry ``carbopol/fit_hb.py``'s stress-space
    values, 0.1-4 % from the identify layer's log-viscosity fit) but not a number the
    fit excludes; an edited record does not keep its verdicts."""
    for name, stored in zip(model.param_names, params_si):
        est = sc_block["parameters"][name]
        lo, hi = est.get("lo"), est.get("hi")
        if lo is None or hi is None or not (float(lo) <= float(stored) <= float(hi)):
            raise ProvenanceError(f"{path}: record {name} = {stored:g} lies outside the sidecar's "
                                  f"interval [{lo}, {hi}]; the verdicts are not this record's")


def load_all(directory: str | Path = RECORDS_DIR, allow_network: bool = False) -> dict[str, FluidRecord]:
    """Every ``*.yaml`` record in ``directory``, keyed by id; duplicate ids raise."""
    out: dict[str, FluidRecord] = {}
    for p in sorted(Path(directory).glob("*.yaml")):
        if p.name.endswith(".identifiability.yaml"):
            continue  # sidecars are read by load_record beside their record
        r = load_record(p, allow_network=allow_network)
        if r.id in out:
            raise ValueError(f"duplicate record id {r.id!r}: {out[r.id].path} and {p}")
        out[r.id] = r
    return out
