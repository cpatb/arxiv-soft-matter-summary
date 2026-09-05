"""Unit 0.2 gate: the loader refuses what it should and loads what it should.

Every rule here is a rule the design documents state (03_ARCHITECTURE.md §4 data
contracts): no number without provenance, envelope required, identifiability written by
the identify layer only, out-of-envelope evaluation raises.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pytest
import yaml

from biorheo.registry import (
    EnvelopeError,
    FluidRecord,
    ProvenanceError,
    load_all,
    load_record,
)
from biorheo.registry.provenance import DOI_RE, is_doi_syntax, resolve_doi

# A DOI that is in the shipped cache (Biermann & Mesecke 2024 Carbopol ETD2020 dataset,
# resolved on 2026-09-05). The cache is what makes these tests run offline.
CACHED_DOI = "10.25835/hjn2qcvy"
FAKE_DOI = "10.0000/this-doi-does-not-exist-biorheo-test"


def _write(tmp_path: Path, text: str) -> Path:
    p = tmp_path / "rec.yaml"
    p.write_text(textwrap.dedent(text), encoding="utf-8")
    return p


GOOD = f"""
    id: test_hb_surrogate
    fluid: Carbopol ETD2020 0.5 wt%, 40 uL TEA per 100 g (test copy)
    model: herschel_bulkley
    parameters:
      tau_y: {{value: 23.756, unit: Pa}}
      k: {{value: 8.903, unit: "Pa*s**0.4526"}}
      n: {{value: 0.4526, unit: dimensionless}}
    envelope:
      gammadot_min: {{value: 1.0e-3, unit: 1/s}}
      gammadot_max: {{value: 200.0, unit: 1/s}}
      temperature_C: 20.0
      geometry: smooth 50 mm parallel plates, 0.5 mm gap
    sources:
      - doi: {CACHED_DOI}
        role: dataset
    licence:
      name: CC BY-NC-SA 4.0
      commercial_use: false
    notes: smooth plates, so tau_y is a lower bound (slip uncontrolled)
"""


# --- provenance -------------------------------------------------------------------

def test_doi_syntax():
    assert is_doi_syntax("10.1122/1.5016034")
    assert is_doi_syntax(CACHED_DOI)
    assert not is_doi_syntax("doi:10.1122/1.5016034")  # prefix is not part of a DOI
    assert not is_doi_syntax("https://doi.org/10.1122/1.5016034")
    assert not is_doi_syntax("11.1122/x")
    assert DOI_RE.pattern.startswith("^10\\.")


def test_cached_doi_resolves_offline():
    r = resolve_doi(CACHED_DOI, allow_network=False)
    assert r is not None
    assert r.doi == CACHED_DOI
    assert r.resolved_on  # a date string


def test_fake_doi_does_not_resolve_offline():
    assert resolve_doi(FAKE_DOI, allow_network=False) is None


# --- loader: rejections -------------------------------------------------------------

def test_good_record_loads(tmp_path):
    rec = load_record(_write(tmp_path, GOOD), allow_network=False)
    assert isinstance(rec, FluidRecord)
    assert rec.parameters["tau_y"].si == pytest.approx(23.756)
    assert rec.parameters["n"].si == pytest.approx(0.4526)
    assert rec.commercial_ok is False
    assert rec.identifiability is None


def test_fake_doi_rejected(tmp_path):
    text = GOOD.replace(CACHED_DOI, FAKE_DOI)
    with pytest.raises(ProvenanceError, match="resolve"):
        load_record(_write(tmp_path, text), allow_network=False)


def test_no_sources_rejected(tmp_path):
    d = yaml.safe_load(textwrap.dedent(GOOD))
    d["sources"] = []
    p = tmp_path / "rec.yaml"
    p.write_text(yaml.safe_dump(d), encoding="utf-8")
    with pytest.raises(ProvenanceError, match="source"):
        load_record(p, allow_network=False)


def test_parameter_without_unit_rejected(tmp_path):
    d = yaml.safe_load(textwrap.dedent(GOOD))
    d["parameters"]["tau_y"] = {"value": 23.756}
    p = tmp_path / "rec.yaml"
    p.write_text(yaml.safe_dump(d), encoding="utf-8")
    with pytest.raises(ValueError, match="unit"):
        load_record(p, allow_network=False)


def test_si_field_in_file_rejected(tmp_path):
    # si is computed at load, never stored: a stored one can drift from value+unit.
    d = yaml.safe_load(textwrap.dedent(GOOD))
    d["parameters"]["tau_y"]["si"] = 23.756
    p = tmp_path / "rec.yaml"
    p.write_text(yaml.safe_dump(d), encoding="utf-8")
    with pytest.raises(ValueError, match="si"):
        load_record(p, allow_network=False)


def test_inverted_envelope_rejected(tmp_path):
    d = yaml.safe_load(textwrap.dedent(GOOD))
    d["envelope"]["gammadot_min"], d["envelope"]["gammadot_max"] = (
        d["envelope"]["gammadot_max"], d["envelope"]["gammadot_min"])
    p = tmp_path / "rec.yaml"
    p.write_text(yaml.safe_dump(d), encoding="utf-8")
    with pytest.raises(EnvelopeError):
        load_record(p, allow_network=False)


def test_hand_written_identifiability_rejected(tmp_path):
    d = yaml.safe_load(textwrap.dedent(GOOD))
    d["identifiability"] = {"tau_y": "measurable"}
    p = tmp_path / "rec.yaml"
    p.write_text(yaml.safe_dump(d), encoding="utf-8")
    with pytest.raises(ProvenanceError, match="identify"):
        load_record(p, allow_network=False)


def test_missing_licence_rejected(tmp_path):
    d = yaml.safe_load(textwrap.dedent(GOOD))
    del d["licence"]
    p = tmp_path / "rec.yaml"
    p.write_text(yaml.safe_dump(d), encoding="utf-8")
    with pytest.raises(ValueError, match="licence"):
        load_record(p, allow_network=False)


def test_unknown_model_name_rejected(tmp_path):
    text = GOOD.replace("model: herschel_bulkley", "model: not_a_model")
    with pytest.raises(ValueError, match="model"):
        load_record(_write(tmp_path, text), allow_network=False)


def test_wrong_parameter_set_for_model_rejected(tmp_path):
    d = yaml.safe_load(textwrap.dedent(GOOD))
    del d["parameters"]["n"]
    p = tmp_path / "rec.yaml"
    p.write_text(yaml.safe_dump(d), encoding="utf-8")
    with pytest.raises(ValueError, match="parameters"):
        load_record(p, allow_network=False)


# --- envelope enforcement at evaluation -------------------------------------------

def test_eta_inside_envelope_matches_hand_value(tmp_path):
    rec = load_record(_write(tmp_path, GOOD), allow_network=False)
    # HB in rate space: eta = tau_y/gd + k gd^(n-1). At gd = 10 1/s:
    # 23.756/10 + 8.903 * 10**(0.4526-1) = 2.3756 + 8.903 * 10**(-0.5474)
    # 10**(-0.5474) = 0.28356..., so 8.903*0.28356 = 2.5246; total 4.9002 Pa s.
    assert rec.eta_at(10.0) == pytest.approx(4.9002, rel=2e-4)


def test_eta_outside_envelope_raises(tmp_path):
    rec = load_record(_write(tmp_path, GOOD), allow_network=False)
    with pytest.raises(EnvelopeError):
        rec.eta_at(1000.0)
    with pytest.raises(EnvelopeError):
        rec.eta_at(np.array([1.0, 1e-5]))


# --- the shipped records ------------------------------------------------------------

def test_wrong_unit_dimension_rejected(tmp_path):
    d = yaml.safe_load(textwrap.dedent(GOOD))
    d["parameters"]["tau_y"] = {"value": 23.756, "unit": "Pa*s"}  # a viscosity, not a stress
    p = tmp_path / "rec.yaml"
    p.write_text(yaml.safe_dump(d), encoding="utf-8")
    with pytest.raises(ValueError, match="tau_y"):
        load_record(p, allow_network=False)


def test_shipped_records_load_offline():
    """Every record in the package loads with the network disallowed: all DOIs cached,
    all rules satisfied. The Carbopol set is 15 conditions (5 wt% x 3 TEA doses)."""
    recs = load_all(allow_network=False)
    carbopol = [r for r in recs.values() if r.id.startswith("carbopol_etd2020_")]
    assert len(carbopol) == 15
    for r in carbopol:
        assert r.model.name == "herschel_bulkley"
        assert r.commercial_ok is False  # CC BY-NC-SA
        assert {s.role for s in r.sources} == {"dataset", "paper"}
        # each curve's own recovered rate range: set points 1e-3 and 200 1/s
        assert r.envelope.gammadot_min == pytest.approx(1e-3, rel=0.02)
        assert r.envelope.gammadot_max == pytest.approx(200.0, rel=0.02)
        # the identify sidecar is attached and its interval holds the stored value
        assert r.verdicts == {"tau_y": "measurable", "k": "measurable", "n": "measurable"}
        for name, blk in r.identifiability["parameters"].items():
            assert blk["lo"] <= r.parameters[name].si <= blk["hi"]
        # the sidecar's fitted range is the record's envelope (record stores 6 sig figs)
        assert r.identifiability["fit"]["gammadot_min"] == pytest.approx(r.envelope.gammadot_min, rel=1e-5)
        assert r.identifiability["fit"]["gammadot_max"] == pytest.approx(r.envelope.gammadot_max, rel=1e-5)


def test_shipped_blood_records_carry_their_convention():
    """The 20 Horner 2020 donors (A-U without E) ship as Casson records whose values are
    biorheo's own log-eta fit, with the authors' sqrt-space values beside them in extra.
    Two facts pinned from scripts/records_from_horner.py's 2026-09-05 run: the authors'
    table reproduces to 1e-6 or better, and the Casson form is flagged structured on
    19 of 20 donors (3 runs of residual sign on n = 16 for 18 of them; donor D is the exception)."""
    recs = load_all(allow_network=False)
    blood = [r for r in recs.values() if r.id.startswith("blood_horner2020_donor_")]
    assert sorted(r.extra["donor"] for r in blood) == [c for c in "ABCDFGHIJKLMNOPQRSTU"]
    structured = 0
    for r in blood:
        assert r.model.name == "casson" and r.commercial_ok is False
        assert r.envelope.temperature_C == 37.0
        assert r.envelope.gammadot_min == pytest.approx(0.1) and r.envelope.gammadot_max == pytest.approx(700.0)
        assert float(r.extra["authors_casson_reproduced_to"]) < 1e-6
        # the record is the log-eta fit, never the authors' number: tau_y is lower on every donor
        assert r.parameters["tau_y"].si < r.extra["authors_casson_tau_y_Pa"]
        assert r.identifiability["fit"]["n_points"] == 16
        structured += bool(r.identifiability["residual_structure"]["structured"])
    assert structured == 19
    assert recs["blood_horner2020_donor_D"].identifiability["residual_structure"]["structured"] is False
    assert {r.extra["paper_outlier_flag"] for r in blood if r.extra["donor"] in "FNK"} == {True}


def test_shipped_carbopol_record_matches_fit_table():
    """carbopol/hb_fits_etd2020.csv row '0.5 wt%, 40 uL TEA': tau_y 23.756, K 8.903,
    n 0.4526; hand value at 10 1/s is 4.9002 Pa s (see test above)."""
    r = load_all(allow_network=False)["carbopol_etd2020_0p5wt_40uL"]
    np.testing.assert_allclose(r.params_si, [23.756, 8.903, 0.4526], rtol=0)
    assert r.eta_at(10.0) == pytest.approx(4.9002, rel=2e-4)
    assert r.extra["fit_R2"] == pytest.approx(0.99999)


# --- identifiability sidecars: only biorheo.identify may write verdicts ----------------

def _sidecar(tmp_path: Path, **override) -> Path:
    block = {
        "written_by": "biorheo.identify", "record_id": "test_hb_surrogate",
        "model": "herschel_bulkley", "level": 0.95,
        "fit": {"n_points": 54, "rmse_log_eta": 0.01, "gammadot_min": 1e-3, "gammadot_max": 200.0},
        "parameters": {
            "tau_y": {"value": 23.756, "lo": 23.0, "hi": 24.5, "verdict": "measurable"},
            "k": {"value": 8.903, "lo": 8.5, "hi": 9.3, "verdict": "measurable"},
            "n": {"value": 0.4526, "lo": 0.44, "hi": 0.46, "verdict": "measurable"},
        },
    }
    block.update(override)
    sc = tmp_path / "rec.identifiability.yaml"
    sc.write_text(yaml.safe_dump(block), encoding="utf-8")
    return sc


def test_sidecar_is_attached_when_marked(tmp_path):
    p = _write(tmp_path, GOOD)
    _sidecar(tmp_path)
    r = load_record(p)
    assert r.verdicts == {"tau_y": "measurable", "k": "measurable", "n": "measurable"}
    assert r.identifiability["fit"]["n_points"] == 54


def test_no_sidecar_means_no_verdicts(tmp_path):
    r = load_record(_write(tmp_path, GOOD))
    assert r.identifiability is None and r.verdicts is None


@pytest.mark.parametrize("override, match", [
    ({"written_by": "hand"}, "not marked as written by"),
    ({"record_id": "some_other_record"}, "is for record"),
    ({"model": "casson"}, "sidecar fitted"),
    ({"parameters": {"tau_y": {"verdict": "measurable"}}}, "do not match"),
    # reviewer 2026-09-05: an edited record must not keep its verdicts. k = 8.903 in the
    # record against a sidecar interval that excludes it
    ({"parameters": {
        "tau_y": {"value": 23.756, "lo": 23.0, "hi": 24.5, "verdict": "measurable"},
        "k": {"value": 4.918, "lo": 4.73, "hi": 5.11, "verdict": "measurable"},
        "n": {"value": 0.4526, "lo": 0.44, "hi": 0.46, "verdict": "measurable"}}},
     "outside the sidecar's interval"),
])
def test_sidecar_rejected_not_ignored(tmp_path, override, match):
    p = _write(tmp_path, GOOD)
    _sidecar(tmp_path, **override)
    with pytest.raises(ProvenanceError, match=match):
        load_record(p)


def test_load_all_skips_sidecars_and_reads_them_beside_records(tmp_path):
    _write(tmp_path, GOOD)
    _sidecar(tmp_path)
    recs = load_all(tmp_path)
    assert list(recs) == ["test_hb_surrogate"]
    assert recs["test_hb_surrogate"].verdicts["tau_y"] == "measurable"


def test_write_identifiability_refuses_a_foreign_fit(tmp_path):
    """A sidecar whose interval does not hold the record's own value describes some other
    fit than the numbers the record quotes; the writer refuses rather than attaching it."""
    from biorheo.constitutive import get_model
    from biorheo.identify import fit_rate_curve, sidecar_path, synthetic_rate_curve, write_identifiability

    rec = load_record(_write(tmp_path, GOOD))
    hb = get_model("herschel_bulkley")
    gd = np.logspace(-2, 2, 25)
    foreign = fit_rate_curve(hb, gd, synthetic_rate_curve(hb, [2.0, 1.0, 0.6], gd, 0.03, 0))
    with pytest.raises(ValueError, match="not the fit's"):
        write_identifiability(rec, foreign)
    assert not sidecar_path(rec.path).exists()

    # seed 1: seed 0's 95 % interval for n is one of the expected 1-in-20 misses of the
    # truth (test_identify documents it), and the writer refusing it is correct behaviour
    own = fit_rate_curve(hb, gd, synthetic_rate_curve(hb, rec.params_si, gd, 0.03, 1))
    out = write_identifiability(rec, own, data_note="synthetic self-fit")
    assert out == sidecar_path(rec.path)
    reloaded = load_record(rec.path)
    assert reloaded.verdicts == {"tau_y": "measurable", "k": "measurable", "n": "measurable"}
    assert reloaded.identifiability["data_note"] == "synthetic self-fit"
