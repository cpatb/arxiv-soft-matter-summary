"""OpenFOAM v2606 export round trip (NEXT_STEPS item 3, gate pre-registered there):
record -> transportProperties text -> parsed coefficients -> OpenFOAM's own calcNu form
reproduces ``eta_at / rho`` on 20 rates inside the envelope to a relative 1e-10, and the
export refuses a degenerate / not_identifiable parameter without the override flag.

Anchors outside the module under test: (1) OpenFOAM's ``foamDictionary`` parses the
written file and returns the same numbers (needs WSL + /usr/lib/openfoam/openfoam2606,
skipped otherwise); (2) two hand-evaluated points of the C++ forms, including one where
the Casson clamp binds; (3) ``biorheo.constitutive.regularise(..., "clamp")``, a separate
code path, equals the exported law everywhere including where the cap binds."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from biorheo.constitutive.regularise import regularise
from biorheo.registry.loader import load_all
from biorheo.solver.openfoam import (
    ExportRefused, VSMALL, coefficients, export_report, export_transport_properties,
    nu_openfoam, parse_transport_properties, structured_flag, write_transport_properties,
)

N_RATES = 20
RTOL_GATE = 1e-10
RHO = 1060.0  # kg/m^3; a round test value, not a material property (records carry none)


@pytest.fixture(scope="module")
def records():
    recs = load_all()
    out = {k: r for k, r in recs.items() if r.model.name in ("herschel_bulkley", "casson")}
    assert len(out) == 35, len(out)  # 15 Carbopol HB + 20 blood Casson
    return out


@pytest.fixture(scope="module")
def all_records():
    """Every exportable record, all four families. The `records` fixture above is
    deliberately the two yield-stress families only, because most of this file is about
    the clamp, which the other two do not have."""
    return load_all()


def _rates(rec):
    return np.geomspace(rec.envelope.gammadot_min, rec.envelope.gammadot_max, N_RATES)


def test_round_trip_reproduces_eta_at_on_every_record(records):
    worst = 0.0
    for rec in records.values():
        gd = _rates(rec)
        eta_max = 2.0 * float(rec.eta_at(gd[0]))  # cap above the envelope: never binds
        rep = export_report(rec, RHO, eta_max, gd)
        assert rep["n_bound"] == 0 and rep["rate_where_cap_binds"] is None
        assert rep["max_rel_err_unbound"] <= RTOL_GATE, (rec.id, rep["max_rel_err_unbound"])
        worst = max(worst, rep["max_rel_err_unbound"])
        # the text carries the density and the clamp, so a reader can undo the export
        assert f"rho = {RHO:.10g}" in rep["text"] and f"eta_max = {eta_max:.10g}" in rep["text"]
    assert worst > 0  # the comparison is not trivially exact: %.17g and a division happen


def test_cap_binding_is_reported_and_equals_biorheo_clamp(records):
    for rid in ("carbopol_etd2020_0p1wt_30uL", "blood_horner2020_donor_A"):
        rec = records[rid]
        gd = _rates(rec)
        eta_max = float(rec.eta_at(gd[N_RATES // 2]))  # cap at the middle rate
        rep = export_report(rec, RHO, eta_max, gd)
        assert 0 < rep["n_bound"] < N_RATES
        assert rep["rate_where_cap_binds"] is not None
        assert rep["rate_where_cap_binds"] <= gd[N_RATES // 2]
        assert np.all(rep["nu_openfoam"][rep["cap_binds"]] == eta_max / RHO)
        assert rep["max_rel_err_unbound"] <= RTOL_GATE
        # independent path: biorheo's clamp scheme is the same function
        clamped = regularise(rec.model, "clamp", eta_max=eta_max)
        eta_clamp = clamped(gd, *rec.params_si)
        np.testing.assert_allclose(rep["nu_openfoam"] * RHO, eta_clamp, rtol=RTOL_GATE)


def test_hand_evaluated_points_of_the_cpp_forms():
    # Bingham (n = 1): nu = min(nu0, tau0/sr + k); at sr = 2, tau0 = 1, k = 3, nu0 = 10 -> 3.5
    assert nu_openfoam("HerschelBulkley", {"k": 3.0, "n": 1.0, "tau0": 1.0, "nu0": 10.0}, 2.0) == 3.5
    # same coefficients, sr = 0.1: 1/0.1 + 3 = 13 > nu0 -> 10
    assert nu_openfoam("HerschelBulkley", {"k": 3.0, "n": 1.0, "tau0": 1.0, "nu0": 10.0}, 0.1) == 10.0
    # OpenFOAM's own Casson.H blood example: m 3.934986e-6, tau0 2.9032e-6, nuMax 13.3333e-6,
    # nuMin 3.9047e-6. At sr = 1: (sqrt(2.9032e-6) + sqrt(3.934986e-6))^2
    # = (1.703878e-3 + 1.983680e-3)^2 = 1.359810e-5 > nuMax -> 1.33333e-5 (clamp binds)
    c = {"m": 3.934986e-6, "tau0": 2.9032e-6, "nuMin": 3.9047e-6, "nuMax": 13.3333e-6}
    assert nu_openfoam("Casson", c, 1.0) == 13.3333e-6
    # sr = 100: (sqrt(2.9032e-8) + 1.983680e-3)^2 = (1.703878e-4 + 1.983680e-3)^2 = 4.64000e-6
    # rtol, not a 1e-9 absolute bar: that was a 0.02 % relative tolerance on a value
    # quoted to 6 figures, ~200x looser than the anchor it claims to be
    assert np.isclose(nu_openfoam("Casson", c, 100.0), 4.64000e-6, rtol=2e-6, atol=0)
    # sr = 0 is floored at VSMALL, not a division by zero: HB gives nu0, Casson gives nuMax
    assert nu_openfoam("HerschelBulkley", {"k": 3.0, "n": 0.5, "tau0": 1.0, "nu0": 10.0}, 0.0) == 10.0
    assert nu_openfoam("Casson", c, 0.0) == c["nuMax"]
    assert VSMALL == 1.0e-300


def test_refuses_unidentified_parameters(records):
    rec = records["blood_horner2020_donor_A"]
    ident = dict(rec.identifiability)
    params = {k: dict(v) for k, v in ident["parameters"].items()}
    params["tau_y"]["verdict"] = "degenerate"
    bad = replace(rec, identifiability={**ident, "parameters": params})
    with pytest.raises(ExportRefused, match="degenerate"):
        export_transport_properties(bad, RHO, 1.0)
    text = export_transport_properties(bad, RHO, 1.0, allow_unidentified=True)
    assert "degenerate" in text  # the override is written into the file, not hidden
    nosidecar = replace(rec, identifiability=None)
    with pytest.raises(ExportRefused, match="no identifiability sidecar"):
        export_transport_properties(nosidecar, RHO, 1.0)
    export_transport_properties(nosidecar, RHO, 1.0, allow_unidentified=True)


def test_rejects_bad_density_or_cap(records):
    rec = records["blood_horner2020_donor_A"]
    with pytest.raises(ValueError, match="rho"):
        coefficients(rec, 0.0, 1.0)
    # `cross` used to be the unimplemented example here; it is implemented as of
    # 2026-09-05 (NEXT_STEPS item 6), so the check moved to a family that still is not.
    with pytest.raises(NotImplementedError, match="power_law"):
        from biorheo.constitutive.gnf import POWER_LAW
        coefficients(replace(rec, model=POWER_LAW), RHO, 1.0)


@pytest.mark.parametrize("rid", ["carbopol_etd2020_0p1wt_30uL", "blood_horner2020_donor_A"])
def test_refuses_a_clamp_that_makes_the_law_newtonian(records, rid):
    """A cap at or below the viscosity at the TOP of the envelope binds at every measured
    rate, so the exported file would be a constant-viscosity fluid wearing the record's
    name. Both families are shear-thinning, so eta_at(gammadot_max) is the minimum."""
    rec = records[rid]
    eta_hi = float(rec.eta_at(rec.envelope.gammadot_max))
    with pytest.raises(ValueError, match="Newtonian"):
        coefficients(rec, RHO, eta_hi)
    with pytest.raises(ValueError, match="Newtonian"):
        coefficients(rec, RHO, 0.5 * eta_hi)
    coefficients(rec, RHO, 1.0000001 * eta_hi)             # just above: allowed
    if rec.model.name == "casson":
        # the old Casson-floor branch is now unreachable, and this is why: eta_c is the
        # infinite-rate limit, strictly below the viscosity at any finite rate
        eta_c = float(dict(zip(rec.model.param_names, rec.params_si))["eta_c"])
        assert eta_c < eta_hi


def test_parser_rejects_a_lookalike_block_and_a_partial_one():
    """Two silent-wrong-number bugs found by review on 2026-09-05, pinned.

    (1) The Coeffs block was matched unanchored, so a stale `oldHerschelBulkleyCoeffs`
        kept beside the live one won and its numbers were returned with no error.
    (2) OpenFOAM's own tutorials write the full dimensionedScalar form `k k [dims] val;`.
        That parsed to an EMPTY dict and the function returned success."""
    lookalike = """transportModel HerschelBulkley;
oldHerschelBulkleyCoeffs { k 999; n 9; tau0 9; nu0 9; }
HerschelBulkleyCoeffs    { k 0.0046; n 0.5; tau0 2e-6; nu0 1e-3; }
"""
    of_name, c = parse_transport_properties(lookalike, expect="herschel_bulkley")
    assert of_name == "HerschelBulkley" and c["k"] == 0.0046 and c["nu0"] == 1e-3

    tutorial = """transportModel HerschelBulkley;
HerschelBulkleyCoeffs
{
    k       k [ 0 2 -1 0 0 0 0 ] 0.00016;
    n       n [ 0 0 0 0 0 0 0 ] 0.5;
    tau0    tau0 [ 0 2 -2 0 0 0 0 ] 2e-06;
    nu0     nu0 [ 0 2 -1 0 0 0 0 ] 0.001;
}
"""
    _, c = parse_transport_properties(tutorial, expect="herschel_bulkley")
    assert c == {"k": 0.00016, "n": 0.5, "tau0": 2e-06, "nu0": 0.001}

    # a partial block must raise rather than return a successful, incomplete dict
    with pytest.raises(ValueError, match="missing"):
        parse_transport_properties("transportModel Casson;\nCassonCoeffs { m 1; }",
                                   expect="casson")
    # and a key that merely ends in transportModel must not hijack the model name
    with pytest.raises(ValueError, match="no transportModel"):
        parse_transport_properties("mytransportModel Casson;\nCassonCoeffs { m 1; }")


def test_the_export_carries_the_two_things_a_dictionary_cannot_enforce(records):
    """A solver has no envelope and no goodness-of-fit. Both must at least be written
    into the file, and `strict_form` must be able to refuse a refuted model FORM."""
    rec = records["carbopol_etd2020_0p1wt_30uL"]
    eta_max = 10.0 * float(rec.eta_at(rec.envelope.gammadot_min))
    text = export_transport_properties(rec, RHO, eta_max)
    assert "WARNING: the envelope above is a COMMENT" in text
    rs = structured_flag(rec)
    assert rs is not None and rs["structured"] is True     # this record IS a misfit form
    assert "residual_structure is `structured`" in text and f"{rs['n_runs']} runs" in text
    with pytest.raises(ExportRefused, match="FORM misfits"):
        export_transport_properties(rec, RHO, eta_max, strict_form=True)


def test_parser_accepts_dimensioned_and_commented_entries():
    text = """transportModel HerschelBulkley; // blood
HerschelBulkleyCoeffs
{
    k       [0 2 -1 0 0 0 0] 1e-5;
    n       0.5;   // dimless
    tau0    [0 2 -2 0 0 0 0] 2e-6;
    nu0     1e-3;
}
"""
    of_name, c = parse_transport_properties(text)
    assert of_name == "HerschelBulkley"
    assert c == {"k": 1e-5, "n": 0.5, "tau0": 2e-6, "nu0": 1e-3}


def _wsl_openfoam() -> bool:
    if shutil.which("wsl.exe") is None:
        return False
    r = subprocess.run(["wsl.exe", "-e", "test", "-f", "/usr/lib/openfoam/openfoam2606/etc/bashrc"],
                       capture_output=True)
    return r.returncode == 0


def _foam(cmd: str) -> str:
    r = subprocess.run(["wsl.exe", "-e", "bash", "-c",
                        "source /usr/lib/openfoam/openfoam2606/etc/bashrc >/dev/null 2>&1; " + cmd],
                       capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stderr
    return r.stdout.strip()


@pytest.mark.skipif(not _wsl_openfoam(), reason="needs WSL with OpenFOAM v2606")
@pytest.mark.parametrize("rid", ["carbopol_etd2020_0p1wt_30uL", "blood_horner2020_donor_A",
                                 "blood_horner2020_donor_N_cross",
                                 "blood_horner2020_donor_S_carreau_yasuda"])
def test_openfoam_foamDictionary_reads_the_file_back(all_records, rid, tmp_path):
    """OpenFOAM's own parser must recover every coefficient BIT-EXACTLY, which is
    stricter than the 1e-10 gate on the law. `-precision 17` is required and is not a
    convenience: `foamDictionary` prints through an IOstream whose default write
    precision is 6 significant figures, so the DEFAULT invocation returns
    3.63125e-06 for a stored 3.6312547169811322e-06 -- a printer limit, not a parse
    error, and the two are asserted separately below so the distinction is pinned.
    A user reading an exported coefficient back with a bare `foamDictionary -value`
    therefore sees 6 figures; the solver reads the file, not that output."""
    rec = all_records[rid]
    # the two no-yield families take no clamp; the two yield ones require one
    eta_max = (None if rec.model.name in ("cross", "carreau_yasuda")
               else 2.0 * float(rec.eta_at(rec.envelope.gammadot_min)))
    path = write_transport_properties(rec, RHO, eta_max, tmp_path / "transportProperties")
    of_name, c = parse_transport_properties(path.read_text())
    drive, rest = Path(path).drive[0].lower(), Path(path).as_posix()[2:]
    wsl_path = f"/mnt/{drive}{rest}"
    for key, val in c.items():
        got = _foam(f"foamDictionary -precision 17 -entry {of_name}Coeffs/{key} "
                    f"-value '{wsl_path}'")
        assert float(got) == val, (key, got, repr(val))  # exact double, no tolerance
        # and the default precision is exactly %.6g of the same number: lossy printing,
        # not a different value. If OpenFOAM ever changes its default this fails loudly.
        deflt = _foam(f"foamDictionary -entry {of_name}Coeffs/{key} -value '{wsl_path}'")
        assert float(deflt) == float(f"{val:.6g}"), (key, deflt, f"{val:.6g}")
    assert _foam(f"foamDictionary -entry transportModel -value '{wsl_path}'") == of_name
    # `a` is the one OPTIONAL coefficient (BirdCarreau reads it with getOrDefault), so it
    # is the one most likely to be silently absent. Assert OpenFOAM sees it (item 9 of the
    # code review 2026-09-05: nothing checked that OpenFOAM parses the new families).
    if of_name == "BirdCarreau":
        assert "a" in c and c["a"] != 2.0, c
        got = _foam(f"foamDictionary -precision 17 -entry BirdCarreauCoeffs/a "
                    f"-value '{wsl_path}'")
        assert float(got) == c["a"], (got, c["a"])
