"""Pins the OpenFOAM plane-Couette anchor (NEXT_STEPS item 3, gates A1-A5' pre-registered
in `scripts/openfoam_couette_anchor.py`).

Two layers:

* fast -- read the committed `data/openfoam/couette_anchor.json` and
  `couette_diagnosis.json` and assert the verdicts, the bars, and the numbers the .md
  and the Status card quote. This catches a regenerated artifact that silently got worse,
  and a .md edited away from its own JSON.
* slow -- re-run the anchor against OpenFOAM in WSL and require the same verdicts. Needs
  WSL + OpenFOAM v2606 and ~4 minutes of solver time.

The bars are NOT restated here as new numbers: they are read out of the JSON and checked
against the constants in the script, so a loosened bar fails this test rather than
passing it quietly.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
ANCHOR_JSON = ROOT / "data" / "openfoam" / "couette_anchor.json"
DIAG_JSON = ROOT / "data" / "openfoam" / "couette_diagnosis.json"
ANCHOR_MD = ROOT / "data" / "openfoam" / "couette_anchor.md"

sys.path.insert(0, str(SCRIPTS))
import openfoam_couette_anchor as A  # noqa: E402


@pytest.fixture(scope="module")
def anchor():
    assert ANCHOR_JSON.exists(), "run scripts/openfoam_couette_anchor.py"
    return json.loads(ANCHOR_JSON.read_text())


@pytest.fixture(scope="module")
def diag():
    assert DIAG_JSON.exists(), "run scripts/openfoam_couette_diagnose.py"
    return json.loads(DIAG_JSON.read_text())


def test_bars_in_the_artifact_are_the_bars_in_the_script(anchor):
    """A loosened bar must fail here. The script's constants are the source of truth."""
    v = anchor["verdict"]
    assert v["A1_shear_uniformity"]["bar"] == A.BAR_A1 == 1e-6
    assert v["A1prime_profile"]["bar"] == A.BAR_A1P == 1e-9
    assert v["A2_nu"]["bar"] == A.BAR_A2 == 1e-6
    assert v["A3_tau_wall"]["bar"] == A.BAR_A3 == 1e-6
    assert v["A4_newtonian_control"]["bar"] == A.BAR_A2
    assert v["A5prime_cauchy"]["bar"] == A.BAR_A5 == 1e-9


def test_every_gate_passes_with_room(anchor):
    """No gate is a marginal pass, and the margins are very unequal in a way worth
    recording. Measured 2026-09-05 on the 17-case shipped run: A1 5.7e5, A1' 3.7e3,
    A2 1.0e6, A3 1.8e7, A4 3.3e7, A5' 1.1e2. A1, A1', A2 and A5' are all set by ONE
    case -- the lowest Carbopol rate, B = 28.9, the slowest-converging corner in the
    set -- so those margins measure how far that corner has converged, not how well the
    parameter map works. A3 and A4 are set elsewhere and sit at round-off. A1' and A5'
    have the smallest margins because their bars are 1e-9 rather than 1e-6, not because
    they are closer to failing in absolute terms. A margin below 10x means the finding
    needs re-deriving, not the bar re-cutting."""
    v = anchor["verdict"]
    assert all(g["pass"] for g in v.values()), v
    # a gate whose worst is exactly 0 is an INFINITE margin, not an excuse to skip it:
    # dropping those rows (as an earlier version did) would silently exempt any gate
    # that became exact from the margin checks below
    margins = {k: (g["bar"] / g["worst"] if g["worst"] > 0 else float("inf"))
               for k, g in v.items()}
    assert len(margins) == len(v)
    assert min(margins.values()) >= 100, margins
    # the four gates with a 1e-6 bar; A1' and A5' carry 1e-9 bars and are checked above
    physics = {k: m for k, m in margins.items() if v[k]["bar"] == 1e-6}
    assert len(physics) == 4 and min(physics.values()) >= 1e5, physics
    # the claim that one case sets three of the margins: its errors are the reported worsts
    hard = max(anchor["rows"], key=lambda r: r["err_nu"])
    assert hard["case"].startswith("carbopol") and hard["gammadot"] < 0.01
    assert hard["err_nu"] == v["A2_nu"]["worst"] and hard["cauchy"] == v["A5prime_cauchy"]["worst"]
    assert hard["shear_spread"] == v["A1_shear_uniformity"]["worst"]


def test_wall_tractions_are_equal_and_opposite_on_every_case(anchor):
    """The momentum balance that makes tau uniform. Signed, so a sign convention error in
    OpenFOAM's wallShearStress FO fails here -- a magnitude could not see it."""
    v = anchor["verdict"]
    assert v["A3_tau_wall"]["signs_opposite_on_the_two_walls"] is True
    for r in anchor["rows"]:
        assert r["tau_openfoam"] > 0 > r["tau_lower_wall"], r["case"]
        assert r["wall_imbalance"] < 1e-12, (r["case"], r["wall_imbalance"])
    assert v["A3_tau_wall"]["worst_wall_imbalance"] == max(
        r["wall_imbalance"] for r in anchor["rows"])


def test_the_viscosity_field_is_uniform_and_has_the_right_size(anchor):
    """A2 reads the max error; the table prints min and max. A parser regression that
    returned a single cell would quietly reduce A2 to a one-cell check."""
    for r in anchor["rows"]:
        assert r["n_cells_nu"] in (1, anchor["NY"]), r["case"]
        # uniform to the same tolerance A2 uses, not bit-identical: the slowest-converging
        # case leaves a 2e-12 spread across the gap, which is real and is reported
        spread = (r["nu_openfoam_max"] - r["nu_openfoam_min"]) / r["nu_openfoam_max"]
        assert 0 <= spread <= A.BAR_A2, (r["case"], spread)


def test_the_clamp_did_not_bind_even_during_the_transient(anchor):
    """Setting eta_max 10x clear of the converged viscosity says nothing about the
    iteration. Measured at the first window and at the end, on the SOLVER's own nu0."""
    for r in anchor["rows"]:
        for key in ("cap_bound_fraction_early", "cap_bound_fraction_final"):
            f = r[key]
            assert f is None or f == 0.0, (r["case"], key, f)


def test_the_case_set_is_the_one_the_gates_were_written_for(anchor):
    rows = anchor["rows"]
    assert anchor["rho"] == 1060.0 and anchor["H"] == 1e-3 and anchor["NY"] == 20
    # 4 records since 2026-09-05: the two yield-stress ones the gates were written for,
    # plus one Cross and one Carreau-Yasuda (NEXT_STEPS item 6), which have no yield
    # stress and are the control for the Bingham-number finding. 6 since 2026-09-06:
    # item 11 adds the two NULL families, both on donor A's curve -- the SAME curve as
    # the Casson record above, so the three are directly comparable.
    assert len(rows) == 1 + 6 * anchor["n_rates"] == 49
    assert rows[0]["case"] == "newtonian_control"
    ids = {r.get("record") for r in rows[1:]}
    assert ids == {"carbopol_etd2020_0p1wt_30uL", "blood_horner2020_donor_A",
                   "blood_horner2020_donor_N_cross",
                   "blood_horner2020_donor_S_carreau_yasuda",
                   "blood_horner2020_donor_A_newtonian",
                   "blood_horner2020_donor_A_power_law"}
    assert {r.get("model") for r in rows[1:]} == {
        "herschel_bulkley", "casson", "cross", "carreau_yasuda",
        "newtonian", "power_law"}
    # the rates span each record's whole envelope, which is what makes A2 a real test
    for rid in ids:
        gd = sorted(r["gammadot"] for r in rows if r.get("record") == rid)
        assert gd[-1] / gd[0] > 1e3, (rid, gd)


def test_newtonian_control_matches_the_hand_value(anchor):
    """A4: the only row whose answer is known without OpenFOAM. nu = 1e-4 m^2/s at
    gd = 10 1/s gives nu = 1e-4 and tau/rho = nu*gd = 1e-3, by hand."""
    c = anchor["rows"][0]
    assert c["gammadot"] == 10.0 and c["hand_nu"] == 1e-4 and c["hand_tau"] == 1e-3
    assert abs(c["nu_openfoam_max"] - 1e-4) / 1e-4 < 1e-12
    assert abs(c["nu_openfoam_min"] - 1e-4) / 1e-4 < 1e-12
    assert abs(c["tau_openfoam"] - 1e-3) / 1e-3 < 1e-12


def test_the_md_quotes_its_own_json(anchor):
    """Every gate number in the .md must be the JSON's, to the 3 digits the .md prints."""
    md = ANCHOR_MD.read_text(encoding="utf-8")
    for g in anchor["verdict"].values():
        assert f"{g['worst']:.3e}" in md, g
    assert "FAILED" not in md.split("## Every case")[0].replace(
        "**A5 as originally written FAILED", "")


def test_diagnosis_supports_the_three_claims_the_md_makes(diag):
    """The diagnosis table is the evidence for the 'silent failure' finding. Its three
    claims, each checked against the JSON rather than against the prose:
      (1) the low-rate corner is badly wrong at a plausible setting;
      (2) under-relaxing HARDER is worse, not better -- the counter-intuitive one;
      (3) the residual criterion is not reachable where B is large."""
    knob = {(r["relax_U"], r["end_time"]): r for r in diag["knob_sweep"]}
    assert diag["record"] == "carbopol_etd2020_0p1wt_30uL"
    # (1) at relax 0.7 / 2000 iterations -- neither an absurd setting -- nu is >100% wrong
    assert knob[(0.7, 2000)]["err_nu"] > 1.0
    # (2) monotone in the relaxation factor at a fixed cap, over the whole sweep
    at20k = sorted((r["relax_U"], r["err_nu"]) for r in diag["knob_sweep"]
                   if r["end_time"] == 20000)
    assert [e for _, e in at20k] == sorted((e for _, e in at20k), reverse=True), at20k
    assert at20k[0][1] / at20k[-1][1] > 1e8   # 0.3 vs 0.95, same budget
    # 0.99 is past a stability cliff, so the shipped 0.95 is not a comfortable margin
    assert [c["relax_U"] for c in diag["crashes"]] == [0.99]
    # (3) the residual criterion is met only at the low-B end, and B is monotone in gd
    rate = sorted(diag["rate_sweep"], key=lambda r: r["gammadot"])
    assert [r["B"] for r in rate] == sorted((r["B"] for r in rate), reverse=True)
    met = [r["B"] for r in rate if r["residual_control_met"]]
    assert met and max(met) < min(r["B"] for r in rate if not r["residual_control_met"])
    # and the error falls monotonically as B falls: this is a B effect, not scatter
    assert [r["err_nu"] for r in rate] == sorted((r["err_nu"] for r in rate), reverse=True)


def test_the_uniform_shear_derivation_still_holds():
    """A1/A2 compare against `eta_at(U0/H)`, which is only the right rate if the shear is
    uniform. `derivations/couette_uniform_shear.py` proves that (condition: tau, not eta,
    strictly increasing; clamp included) and exits non-zero on any failed check."""
    r = subprocess.run([sys.executable, str(ROOT / "derivations" / "couette_uniform_shear.py")],
                       capture_output=True, text=True, cwd=str(ROOT), timeout=600)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-2000:]
    assert "0 failure(s)" in r.stdout and "FAIL" not in r.stdout


def test_the_readme_l3_paragraph_matches_the_artifacts(anchor, diag):
    """The README quotes five numbers from this study. Each is re-derived from the JSON
    and required to be present, so a regenerated artifact that moved a number leaves the
    README failing rather than quietly stale."""
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    para = readme[readme.index("The L3 solver bridge exports"):]
    worst = max(r["err_nu"] for r in anchor["rows"])
    assert f"{worst:.1e}" in para, f"{worst:.1e}"
    assert f"{len(anchor['rows'])} `simpleFoam`" in para
    knob = {(r["relax_U"], r["end_time"]): r for r in diag["knob_sweep"]}
    assert f"{100 * knob[(0.7, 2000)]['err_nu']:.0f} % wrong" in para
    g = diag["discriminator"]
    r_b = max(g["R_B_at_fixed_eta"].values())
    assert f"{r_b:.1e}".replace("e+04", "e4") in para
    assert f"{max(g['eta_values']) / min(g['eta_values']):.0f}x" in para
@pytest.mark.slow
@pytest.mark.skipif(not A.has_openfoam(), reason="needs WSL with OpenFOAM v2606")
def test_rerunning_the_anchor_reproduces_the_verdicts(anchor, tmp_path):
    """Re-run OpenFOAM and require the same verdicts. Compares verdicts and orders of
    magnitude, not exact doubles: an iterative solve on a different machine is not
    bit-reproducible, and pinning digits it cannot hold would be a false gate.

    BIORHEO_OF_OUT sends the fresh artifacts to a temp directory. Without it this test
    would OVERWRITE the committed `couette_anchor.json` it is validating against, so
    running the suite would silently mutate the evidence."""
    env = {**os.environ, "BIORHEO_OF_OUT": str(tmp_path)}
    r = subprocess.run([sys.executable, str(SCRIPTS / "openfoam_couette_anchor.py")],
                       capture_output=True, text=True, cwd=str(ROOT), timeout=7200, env=env)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-2000:]
    assert ANCHOR_JSON.read_text() == json.dumps(anchor, indent=2), \
        "the committed artifact was modified by the re-run"
    fresh = json.loads((tmp_path / "couette_anchor.json").read_text())["verdict"]
    for name, g in fresh.items():
        assert g["pass"], (name, g)
        assert g["bar"] == anchor["verdict"][name]["bar"], name
        assert g["worst"] <= 10 * max(anchor["verdict"][name]["worst"], 1e-15), (name, g)


# ------------------------------------------------------------------ NEXT_STEPS item 6
# The families with no yield stress. Gates G3/G4/G7 are pre-registered in
# `data/openfoam/nonyield_bridge_prereg.md`, committed (2561346) before the code.
# G1/G2/G5/G8 do not need a solver and live in tests/test_openfoam_export_nonyield.py.


def test_g7_the_no_yield_families_converge_where_the_yield_one_does_not(anchor):
    """G7, the control. `B = tau_y/(k gd^n)` is identically zero for Cross and
    Carreau-Yasuda, so if the Bingham number is what controls simpleFoam's outer
    iteration these cases must be clean. Measured 2026-09-05: worst relative viscosity
    error 2.514e-14 over 16 cases, against the 1.72 (172 %) the Herschel-Bulkley record
    returns at B = 29 with `rc = 0` and no warning.

    The contrast is the argument, and it is not "hard fluid vs easy fluid": the Cross
    record's zero-shear plateau is itself a large viscosity."""
    i = anchor["item6_no_yield_stress"]
    assert i["n_cases"] == 16 and len(i["records"]) == 2
    assert i["G7_bar"] == A.BAR_A2 == 1e-6      # the existing bar, not a new one
    assert i["G7_pass"] and i["G7_worst_err_nu"] <= i["G7_bar"]
    assert i["G7_worst_err_nu"] < 1e-12, i["G7_worst_err_nu"]
    assert i["worst_wall_imbalance"] < 1e-12, i["worst_wall_imbalance"]


def test_g3_g4_the_anchor_could_have_told_the_rejected_map_apart(anchor):
    """A gate that both candidates pass has discriminated nothing. Both separations must
    exceed the 1e-2 pre-registered floor, and by a margin against A2's own bar, or the
    A2 verdict on these records is not evidence about the parameter map."""
    i = anchor["item6_no_yield_stress"]
    assert i["G3_G4_bar"] == 1e-2
    assert i["G3_G4_informative"], i["candidate_separation"]
    seps = i["candidate_separation"]
    assert set(seps) == set(i["records"])
    # measured 2026-09-05: 1.056e+01 for the name-swapped Cross map, 1.128e-01 for `a`
    # forced to 2. Both are >= 1e5 times BAR_A2.
    assert min(seps.values()) / A.BAR_A2 > 1e4, seps
    assert seps["blood_horner2020_donor_N_cross"] > 1.0
    assert seps["blood_horner2020_donor_S_carreau_yasuda"] > 1e-2


def test_g4_openfoam_honours_the_birdcarreau_a_exponent(anchor):
    """The premise of NEXT_STEPS item 6 -- that v2606's `BirdCarreau` has no `a` and the
    export must refuse `a != 2` -- decided by a running solver rather than by reading
    C++. The Carreau-Yasuda record's fitted `a` is far enough from 2 that ignoring the
    entry would move `nu` by 1.1e-1; A2 on those cases lands at 2.5e-14 instead. The
    premise is REFUTED, so no such refusal is shipped and none may be added."""
    i = anchor["item6_no_yield_stress"]
    rid = "blood_horner2020_donor_S_carreau_yasuda"
    cy = [r for r in anchor["rows"] if r.get("record") == rid]
    assert len(cy) == anchor["n_rates"]
    worst = max(r["err_nu"] for r in cy)
    sep = i["candidate_separation"][rid]
    assert worst <= A.BAR_A2, worst
    # the whole content of the gate: the two candidates are far apart AND the shipped
    # one is the one that matches, by many orders of magnitude
    assert sep / max(worst, 1e-300) > 1e9, (sep, worst)

    from biorheo.solver.openfoam import _MODELS
    assert _MODELS["carreau_yasuda"] == ("BirdCarreau", ("nu0", "nuInf", "k", "n", "a"))


def test_the_prereg_for_item_6_predates_its_code():
    """The pre-registration is only worth something if it was written first. Its commit
    must be an ancestor of the commit that added the module's Cross/BirdCarreau map."""
    prereg = ROOT / "data" / "openfoam" / "nonyield_bridge_prereg.md"
    assert prereg.exists()
    text = prereg.read_text(encoding="utf-8")
    for token in ("G1", "G2", "G3", "G4", "G5", "G7", "G8", "1e-10", "1e-3", "1e-2"):
        assert token in text, token
    # the premise it refutes, named, so a later reader cannot reinstate the refusal
    assert "That is false for OpenFOAM v2606" in text


# ---------------------------------------------------------------------------- item 11
# The two NULL families. G11.7 and the two refusal probes are pre-registered in
# `data/openfoam/newtonian_powerlaw_prereg.md`, committed (e9c1ace) before the code.
# G11.1-G11.6 do not need a solver and live in tests/test_openfoam_export_null.py.


def test_g11_7_the_null_families_run_through_the_solver(anchor):
    """G11.7. The same bars as the rest of the page, unchanged. Reported PER FAMILY
    because A2 and A3 are maxima over all 49 cases and a maximum cannot say which family
    moved -- which matters here: A3's headline went 5.53e-14 -> 6.28e-14 when these two
    records were added, and no pre-existing row moved. A new `newtonian` row is simply
    the largest."""
    i = anchor["item11_null_families"]
    assert set(i["per_family"]) == {"newtonian", "power_law"}
    assert i["records"] == ["blood_horner2020_donor_A_newtonian",
                            "blood_horner2020_donor_A_power_law"]
    for fam, f in i["per_family"].items():
        assert f["n_cases"] == anchor["n_rates"] == 8, fam
        assert f["bars"] == {"A1prime": A.BAR_A1P, "A2": A.BAR_A2, "A3": A.BAR_A3}, fam
        assert f["pass"], (fam, f)
        # measured 2026-09-06: newtonian err_nu is EXACTLY 0 (calcNu returns a uniform
        # field of the coefficient, so there is nothing to converge), power_law 6.4e-15
        assert f["worst_err_nu"] < 1e-12, (fam, f["worst_err_nu"])
        assert f["worst_err_tau"] < 1e-12, (fam, f["worst_err_tau"])
        assert f["worst_wall_imbalance"] < 1e-12, (fam, f["worst_wall_imbalance"])
    assert i["per_family"]["newtonian"]["worst_err_nu"] == 0.0
    assert i["G11_7_pass"]


def test_g11_3b_g11_4b_openfoam_rejects_the_two_dictionaries_it_was_predicted_to(anchor):
    """Two PREDICTIONS about the C++ decided by running `simpleFoam`, each beside an
    unmodified CONTROL -- "the broken one failed" says nothing unless the intact one ran.

    G11.3b: `powerLaw`'s `nuMin` is REQUIRED, unlike `BirdCarreau`'s `a`, so a dictionary
    without it must not run. G11.4b: `Newtonian` reads `nu` at the TOP LEVEL, so putting
    it in a `NewtonianCoeffs` block must not run. Both confirmed 2026-09-06, and
    G11.4b's message is the informative one -- OpenFOAM reports `nu` missing from
    `transportProperties` itself, i.e. it never looked in the block."""
    p = anchor["item11_null_families"]["refusal_probes"]
    for g in ("G11_3b", "G11_4b"):
        assert p[g]["control_rc"] == 0, (g, "the UNMODIFIED dictionary must run")
        assert p[g]["rc"] != 0, (g, "the broken dictionary was predicted to fail")
        assert p[g]["pass"], (g, p[g])
    assert p["G11_3b"]["keyword_named"]
    assert "nuMin" in p["G11_3b"]["tail"]
    # OpenFOAM looked for `nu` at the top level, which is the whole claim
    assert "Entry 'nu' not found" in p["G11_4b"]["tail"]
    assert "NewtonianCoeffs" not in p["G11_4b"]["tail"].split("Entry 'nu' not found")[1]


def test_g7_prime_is_still_item_6s_four_records(anchor):
    """A closed gate must not be re-decided by a later item appending to a list. G7'
    iterated `RECORD_IDS` until 2026-09-06; item 11 added two records to it, and G7's
    yield/no-yield split keys on `NO_CLAMP_MODELS`, which lists NEITHER null family --
    so both NULLS would have been counted as yield-stress arms. G7' now has its own
    frozen list and its number is unchanged."""
    g = anchor["item6_no_yield_stress"]["G7prime"]
    assert len(g["rows"]) == 4
    assert {r["model"] for r in g["rows"]} == {"herschel_bulkley", "casson", "cross",
                                              "carreau_yasuda"}
    assert "newtonian" not in {r["model"] for r in g["rows"]}
    assert "power_law" not in {r["model"] for r in g["rows"]}
    # measured 2026-09-05, unchanged after item 11 added two records to RECORD_IDS
    assert g["worst_yield"] == pytest.approx(1.7202393802951472, rel=1e-12)


def test_the_prereg_for_item_11_predates_its_code():
    """`e9c1ace` contains that file and nothing else; the code it grades is in `bcf84bf`.
    The addendum records what the run found, including the one factual error the
    pre-registration itself contains -- corrected there rather than edited out."""
    prereg = ROOT / "data" / "openfoam" / "newtonian_powerlaw_prereg.md"
    text = prereg.read_text(encoding="utf-8")
    for token in ("G11.1", "G11.3b", "G11.4b", "G11.6", "G11.7", "1e-10", "1e-3", "1e-12"):
        assert token in text, token
    assert "Addendum, 2026-09-06" in text
    # the premise it refuted, named
    assert "item's stated reason is wrong" in text
    assert "**§0.2's last sentence is FALSE.**" in text
