"""The L3 bridge for the two families that have NO yield stress: ``cross`` ->
``CrossPowerLaw`` and ``carreau_yasuda`` -> ``BirdCarreau`` (NEXT_STEPS item 6).

Every bar here is pre-registered in ``data/openfoam/nonyield_bridge_prereg.md``, which was
written and committed BEFORE this file existed. Nothing in it may be loosened.

  G1  round trip to 1e-10, the bar the two yield-stress families already meet;
  G2  the export refuses a `degenerate` / `not_identifiable` parameter, on exactly the
      records whose sidecar says so and no others;
  G3  the map a NAME-MATCHER would write for Cross -- biorheo's exponent `m` into
      OpenFOAM's `m`, which is a time -- must FAIL the round trip by more than 1e-3, or
      G1 passing at 1e-10 has discriminated nothing;
  G5  a `BirdCarreau` dictionary with `a` deleted parses to a = 2 (OpenFOAM reads it with
      `getOrDefault`, so a missing entry MEANS 2 and reporting it missing would disagree
      with the solver about what the file says) -- and the export always writes it;
  G8  `cross` and `carreau_yasuda` REFUSE an `eta_max` rather than ignoring it: their
      zero-shear plateau is a fitted parameter, not a regularisation.

G4 (does OpenFOAM honour `a` at all?) and G7 (do these families converge where the
yield-stress ones do not?) need a running solver and live in
``scripts/openfoam_couette_anchor.py``, pinned by ``tests/test_openfoam_couette.py``.

A premise of item 6 was WRONG and is recorded rather than worked around: it said v2606's
``BirdCarreau`` has no ``a`` exponent and the export must refuse ``a != 2``. It has one
(``a_`` is ``getOrDefault("a", dimensionedScalar("a", dimless, 2))``, BirdCarreau.C:81-88,
used in ``calcNu`` line 58), so no such refusal is shipped. G5 tests the hazard that is
actually there, which points the other way.
"""

from __future__ import annotations

import re

import numpy as np
import pytest

from biorheo.registry.loader import load_all
from biorheo.solver.openfoam import (
    ExportRefused, coefficients, export_report, export_transport_properties,
    nu_openfoam, parse_transport_properties,
)

N_RATES = 20
RTOL_GATE = 1e-10        # G1, the existing bar
WRONG_MAP_BAR = 1e-3     # G3, the separation below which G1 would prove nothing
NO_YIELD = ("cross", "carreau_yasuda")
RHO = 1060.0  # kg/m^3; a round test value, not a material property (records carry none)


@pytest.fixture(scope="module")
def records():
    out = {k: r for k, r in load_all().items() if r.model.name in NO_YIELD}
    # 20 blood donors x 2 families, written by scripts/records_from_horner.py
    assert len(out) == 40, len(out)
    return out


def _rates(rec):
    return np.geomspace(rec.envelope.gammadot_min, rec.envelope.gammadot_max, N_RATES)


# --------------------------------------------------------------------------------- G1


def test_g1_round_trip_reproduces_eta_at_on_every_record(records):
    """G1: record -> text -> parsed coefficients -> OpenFOAM's calcNu form reproduces
    rho * eta_at to 1e-10. No cap exists for these families, so nothing may bind."""
    worst = {}
    for rid, rec in records.items():
        rep = export_report(rec, RHO, None, _rates(rec), allow_unidentified=True)
        assert rep["n_bound"] == 0, f"{rid}: a cap bound, but these families have none"
        worst[rid] = rep["max_rel_err_unbound"]
        assert rep["max_rel_err_unbound"] <= RTOL_GATE, (rid, rep["max_rel_err_unbound"])
    # the number, so a regression shows as a change and not just a pass
    assert max(worst.values()) < 1e-14, max(worst.items(), key=lambda kv: kv[1])


def test_g1_transport_model_names_and_coefficient_keys(records):
    want = {"cross": ("CrossPowerLaw", {"nu0", "nuInf", "m", "n"}),
            "carreau_yasuda": ("BirdCarreau", {"nu0", "nuInf", "k", "n", "a"})}
    for rid, rec in records.items():
        text = export_transport_properties(rec, RHO, None, allow_unidentified=True)
        of_name, c = parse_transport_properties(text, expect=rec.model.name)
        assert (of_name, set(c)) == want[rec.model.name], rid


def test_g1_the_cross_map_is_the_swapped_one_not_the_name_matched_one(records):
    """The map itself, read off one record: OpenFOAM's `m` must equal biorheo's `lam`
    (a time, order 10-100 s here) and OpenFOAM's `n` must equal biorheo's `m` (an
    exponent, order 0.5). Checking the round trip alone would not say which is which."""
    rec = records["blood_horner2020_donor_N_cross"]
    p = dict(zip(rec.model.param_names, (float(v) for v in rec.params_si)))
    c = coefficients(rec, RHO, None)
    assert c["m"] == p["lam"] and c["n"] == p["m"]
    assert c["nu0"] == p["eta_0"] / RHO and c["nuInf"] == p["eta_inf"] / RHO


# --------------------------------------------------------------------------------- G2


def test_g2_refuses_exactly_the_records_the_sidecar_calls_unidentified(records):
    """G2: no false positives and no false negatives against the sidecar's own verdicts.
    The counts are the result -- one record could not have shown a refusal that depends
    on identifiability."""
    bad = {"degenerate", "not_identifiable"}
    refused = {"cross": 0, "carreau_yasuda": 0}
    for rid, rec in records.items():
        should = any(v in bad for v in rec.verdicts.values())
        try:
            export_transport_properties(rec, RHO, None)
        except ExportRefused as e:
            assert should, f"{rid}: refused but no parameter is degenerate/unidentifiable"
            # the message must name the offending parameter, not just say no
            assert any(n in str(e) for n in rec.model.param_names if rec.verdicts[n] in bad)
            refused[rec.model.name] += 1
        else:
            assert not should, f"{rid}: exported despite {rec.verdicts}"
        # the override is the only way through
        export_transport_properties(rec, RHO, None, allow_unidentified=True)
    # measured 2026-09-06 over the 20 Horner donors; pinned so a refit that changes the
    # registry's identifiability has to be noticed rather than absorbed
    assert refused == {"cross": 9, "carreau_yasuda": 8}, refused


# --------------------------------------------------------------------------------- G3


def test_g3_the_name_matched_cross_map_fails_the_round_trip(records):
    """G3, the discriminator. Writing biorheo's `m` into OpenFOAM's `m` slot and `lam`
    into `n` is what a map built by matching names would do; it raises nothing anywhere,
    because both slots hold a plain scalar. It must miss by more than 1e-3 on at least
    one rate of every Cross record, or G1's 1e-10 is not evidence about the map."""
    seps = {}
    for rid, rec in records.items():
        if rec.model.name != "cross":
            continue
        gd = _rates(rec)
        c = coefficients(rec, RHO, None)
        with np.errstate(over="ignore"):   # the wrong map exponentiates a time constant
            bad = nu_openfoam("CrossPowerLaw", {**c, "m": c["n"], "n": c["m"]}, gd)
        nu_ref = rec.eta_at(gd) / RHO      # the record, NOT nu_openfoam
        seps[rid] = float(np.max(np.abs(bad / nu_ref - 1.0)))
        assert seps[rid] > WRONG_MAP_BAR, (rid, seps[rid])
    assert len(seps) == 20
    # measured 2026-09-06: the weakest separation is 6.3e-2, 63x the bar
    assert min(seps.values()) > 6e-2, min(seps.items(), key=lambda kv: kv[1])


# --------------------------------------------------------------------------------- G5


def _a_entry(text: str) -> re.Match | None:
    return re.search(r"^\s+a\s+([-+0-9.eE]+);\s*$", text, re.M)


def test_g5_export_always_writes_a_even_when_it_equals_two(records):
    """`a` is optional in OpenFOAM, so omitting it is legal and silently means 2. An
    export that dropped it would ship a different fluid with no error anywhere."""
    n = 0
    for rid, rec in records.items():
        if rec.model.name != "carreau_yasuda":
            continue
        text = export_transport_properties(rec, RHO, None, allow_unidentified=True)
        m = _a_entry(text)
        assert m is not None, f"{rid}: no `a` entry in the exported BirdCarreauCoeffs"
        assert float(m.group(1)) == pytest.approx(
            float(dict(zip(rec.model.param_names, rec.params_si))["a"]), rel=1e-15)
        n += 1
    assert n == 20


def test_g5_a_missing_a_parses_to_two_and_changes_the_law(records):
    """Mirroring `getOrDefault("a", 2)`: a dictionary without `a` must PARSE, not raise,
    and the law it then describes must differ from the record by a margin the anchor's
    A2 bar could see -- otherwise the omission would be undetectable either way."""
    rid = "blood_horner2020_donor_S_carreau_yasuda"
    rec = records[rid]
    text = export_transport_properties(rec, RHO, None)
    stripped = re.sub(r"^\s+a\s+[-+0-9.eE]+;\n", "", text, flags=re.M)
    assert _a_entry(stripped) is None
    of_name, c = parse_transport_properties(stripped, expect="carreau_yasuda")
    assert of_name == "BirdCarreau" and c["a"] == 2.0

    gd = _rates(rec)
    nu_ref = rec.eta_at(gd) / RHO
    sep = float(np.max(np.abs(nu_openfoam("BirdCarreau", c, gd) / nu_ref - 1.0)))
    # measured 2026-09-06: 1.128e-1, i.e. 1e5 times the anchor's BAR_A2 = 1e-6. This is
    # also the separation that lets the anchor decide whether OpenFOAM honours `a` (G4).
    # The pre-registration wrote the ratio the other way up, |nu(a)/nu(a=2) - 1|, which
    # is 1.014e-1 on the same numbers; the bar is 1e-2 and both directions clear it.
    assert sep > 1e-2, sep
    assert sep == pytest.approx(0.11279, rel=5e-3)


# --------------------------------------------------------------------------------- G8


@pytest.mark.parametrize("rid", ["blood_horner2020_donor_N_cross",
                                 "blood_horner2020_donor_S_carreau_yasuda"])
def test_g8_refuses_an_eta_max_rather_than_ignoring_it(records, rid):
    """Accepting a clamp these families cannot use and dropping it is the same
    silent-drop failure the bridge exists to prevent."""
    rec = records[rid]
    with pytest.raises(ValueError, match="eta_max"):
        coefficients(rec, RHO, 1.0)
    with pytest.raises(ValueError, match="eta_max"):
        export_transport_properties(rec, RHO, 1.0, allow_unidentified=True)
    coefficients(rec, RHO, None)  # the only accepted call


def test_g8_the_yield_families_still_require_one():
    """The refusal is per family, not a blanket relaxation: dropping the requirement for
    Herschel-Bulkley or Casson would let an unregularised yield-stress law be exported."""
    recs = load_all()
    for rid in ("carbopol_etd2020_0p1wt_30uL", "blood_horner2020_donor_A"):
        with pytest.raises(ValueError, match="requires eta_max"):
            coefficients(recs[rid], RHO, None)


def test_g8_the_exported_header_says_nu0_is_the_record_not_a_clamp(records):
    text = export_transport_properties(records["blood_horner2020_donor_N_cross"],
                                       RHO, None, allow_unidentified=True)
    assert "no clamp" in text and "fitted parameter" in text
    assert "clamp eta_max" not in text
