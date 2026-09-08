"""NEXT_STEPS item 11: the L3 OpenFOAM bridge for `newtonian` and `power_law`, and the
explicit refusal for `structural_steady`.

Every gate here was pre-registered in `data/openfoam/newtonian_powerlaw_prereg.md` and
committed at e9c1ace BEFORE any of the code it grades. The bars are G11.1 <= 1e-10,
G11.5 > 1e-3, G11.6 <= 1e-12; they are not to be moved. G11.3b, G11.4b and G11.7 need
OpenFOAM and live in `test_openfoam_couette.py` / the anchor script.
"""
from __future__ import annotations

import re
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

import biorheo.solver.openfoam as openfoam_module
from biorheo.constitutive import get_model
from biorheo.constitutive.structural import STRUCTURAL_STEADY, STRUCTURAL_STEADY_SEPARATE
from biorheo.solver.openfoam import (
    SMALL, ThixotropyNotRepresentable, _envelope_eta_range, coefficients,
    export_report, export_transport_properties, nu_openfoam,
    parse_transport_properties)

RHO = 1060.0
BAR_G111 = 1e-10        # G11.1: round-trip relative error on eta_at
BAR_G115 = 1e-3         # G11.5: a wrong map must miss by MORE than this
BAR_G116 = 1e-12        # G11.6: structural_steady vs Cross(m = 1)
N_RATES = 20


@pytest.fixture(scope="module")
def records(closed_gate_records):
    # CLOSED item-11 gate: graded on the 115-record registry of 2026-09-07, not on
    # whatever the registry holds today. Item 22 added eight `power_law` records and
    # G11.5b's separations are computed over this sample, so an unfrozen fixture
    # would have re-decided a closed gate on a different question. See
    # tests/conftest.py and data/dmem_record_prereg.md §4.
    return closed_gate_records


def _null_records(records, family):
    out = {rid: r for rid, r in records.items() if r.model.name == family}
    assert out, f"no {family} record in the registry; item 11's gate cannot run"
    return out


def _rates(rec):
    gd = np.logspace(np.log10(rec.envelope.gammadot_min),
                     np.log10(rec.envelope.gammadot_max), N_RATES)
    gd[0], gd[-1] = rec.envelope.gammadot_min, rec.envelope.gammadot_max
    return np.clip(gd, rec.envelope.gammadot_min, rec.envelope.gammadot_max)


def _clamps(rec):
    """A cap and a floor that are clear of the record's own envelope by a decade, so
    neither binds and G11.1 grades the LAW rather than the clamp."""
    lo, hi = _envelope_eta_range(rec)
    return 10.0 * hi, 0.1 * lo


# --------------------------------------------------------------------- G11.1 round trip

def test_g11_1_round_trip_newtonian(records):
    worst = 0.0
    for rid, rec in sorted(_null_records(records, "newtonian").items()):
        rep = export_report(rec, RHO, None, _rates(rec))
        assert rep["transportModel"] == "Newtonian"
        assert not rep["cap_binds"].any() and not rep["floor_binds"].any()
        worst = max(worst, rep["max_rel_err_unbound"])
    assert worst <= BAR_G111, f"G11.1 newtonian: {worst:.3e} > {BAR_G111:.0e}"


def test_g11_1_round_trip_power_law(records):
    worst = 0.0
    for rid, rec in sorted(_null_records(records, "power_law").items()):
        eta_max, eta_min = _clamps(rec)
        rep = export_report(rec, RHO, eta_max, _rates(rec), eta_min=eta_min)
        assert rep["transportModel"] == "powerLaw"
        # the whole point of putting the clamps a decade outside: they must not bind, or
        # the round trip would be grading `min`/`max` instead of `k gd^(n-1)`
        assert not rep["cap_binds"].any() and not rep["floor_binds"].any(), rid
        worst = max(worst, rep["max_rel_err_unbound"])
    assert worst <= BAR_G111, f"G11.1 power_law: {worst:.3e} > {BAR_G111:.0e}"


def test_g11_1_the_clamps_do_bind_when_asked_to(records):
    """The negative half of G11.1: a clamp inside the envelope must actually be reported,
    or `cap_binds`/`floor_binds` being all-False above proves nothing."""
    rec = sorted(_null_records(records, "power_law").items())[0][1]
    gd = _rates(rec)
    eta = np.asarray(rec.eta_at(gd), dtype=float)
    lo, hi = _envelope_eta_range(rec)
    # a cap at the median viscosity: binds at the low-rate half, not the high
    rep = export_report(rec, RHO, float(np.median(eta)), gd, eta_min=0.1 * lo)
    assert rep["cap_binds"].any() and not rep["cap_binds"].all()
    assert rep["rate_where_cap_binds"] is not None
    # a floor at the median: the mirror image
    rep = export_report(rec, RHO, 10.0 * hi, gd, eta_min=float(np.median(eta)))
    assert rep["floor_binds"].any() and not rep["floor_binds"].all()
    assert rep["rate_where_floor_binds"] is not None


# ------------------------------------------------------------------- G11.2 the verdicts

def test_g11_2_verdict_refusal_reaches_the_new_families(records):
    """A `degenerate` verdict must stop the export for the new families too.

    Neither shipped family carries one -- all 20 `newtonian` and all 20 `power_law`
    records are all-`measurable`, which is itself worth pinning -- so the gate is graded
    on a record with its sidecar verdict overwritten, and this test says so out loud."""
    for family in ("newtonian", "power_law"):
        recs = _null_records(records, family)
        assert all("degenerate" not in (r.verdicts or {}).values() for r in recs.values())
        assert all(set((r.verdicts or {}).values()) == {"measurable"} for r in recs.values())
        rec = sorted(recs.items())[0][1]
        eta_max, eta_min = (None, None) if family == "newtonian" else _clamps(rec)
        ident = dict(rec.identifiability)
        params = {k: dict(v) for k, v in ident["parameters"].items()}
        params[rec.model.param_names[0]]["verdict"] = "degenerate"
        bad = replace(rec, identifiability={**ident, "parameters": params})
        with pytest.raises(Exception) as e:
            export_transport_properties(bad, RHO, eta_max, eta_min=eta_min)
        assert rec.model.param_names[0] in str(e.value) and "degenerate" in str(e.value)
        # and the override still works
        export_transport_properties(bad, RHO, eta_max, eta_min=eta_min,
                                    allow_unidentified=True)


# ---------------------------------------------------------------- G11.3a the clamp rules

def test_g11_3a_every_clamp_rule_raises(records):
    """The eight refusals of the clamp decision (prereg §1). Accepting and dropping any
    of them silently is the failure G8 exists to prevent."""
    pl = sorted(_null_records(records, "power_law").items())[0][1]
    nw = sorted(_null_records(records, "newtonian").items())[0][1]
    hb = records["carbopol_etd2020_0p1wt_30uL"]
    ca = records["blood_horner2020_donor_A"]
    cr = records["blood_horner2020_donor_N_cross"]
    cy = records["blood_horner2020_donor_S_carreau_yasuda"]
    eta_max, eta_min = _clamps(pl)

    with pytest.raises(ValueError, match="requires eta_min"):
        coefficients(pl, RHO, eta_max)
    with pytest.raises(ValueError, match="requires eta_max"):
        coefficients(pl, RHO, None, eta_min=eta_min)
    for rec, why in ((hb, "no nuMin slot"), (ca, "eta_c"), (cr, "no limiter slots"),
                     (cy, "no limiter slots"), (nw, "no limiter slot")):
        with pytest.raises(ValueError, match="takes no eta_min"):
            coefficients(rec, RHO, 1.0, eta_min=1e-6)
        assert why in str(pytest.raises(
            ValueError, lambda: coefficients(rec, RHO, 1.0, eta_min=1e-6)).value)
    for rec in (nw, cr, cy):
        with pytest.raises(ValueError, match="takes no eta_max"):
            coefficients(rec, RHO, 1.0)


def test_g11_3a_a_floor_at_or_above_the_cap_is_refused(records):
    pl = sorted(_null_records(records, "power_law").items())[0][1]
    eta_max, eta_min = _clamps(pl)
    lo, hi = _envelope_eta_range(pl)
    # a cap INSIDE the envelope, so the "floor at or above the envelope max" check does
    # not fire first and the floor-vs-cap ordering is what is being graded
    mid = float(np.sqrt(lo * hi))
    with pytest.raises(ValueError, match="must be below eta_max"):
        coefficients(pl, RHO, mid, eta_min=mid)
    with pytest.raises(ValueError, match="the floor would bind everywhere"):
        coefficients(pl, RHO, eta_max, eta_min=hi)
    with pytest.raises(ValueError, match="the clamp would bind everywhere"):
        coefficients(pl, RHO, lo, eta_min=eta_min)
    for bad in (0.0, -1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="eta_min must be a positive"):
            coefficients(pl, RHO, eta_max, eta_min=bad)


def test_envelope_eta_range_is_inert_on_the_registry(records):
    """`_envelope_eta_range` replaced a single-point `eta_at(gammadot_max)` on
    2026-09-06. That was correct only for a monotonically shear-thinning law, which
    `power_law` and `herschel_bulkley` do not guarantee (both admit n > 1). The two must
    agree EXACTLY on every record now shipped, or the change moved a shipped refusal
    rather than generalising it."""
    n = 0
    for rid, rec in sorted(records.items()):
        lo, _ = _envelope_eta_range(rec)
        endpoint = float(rec.eta_at(rec.envelope.gammadot_max))
        assert lo == endpoint, (rid, lo, endpoint)
        n += 1
    assert n == 115, n


# ---------------------------------------------------------- G11.4a the Newtonian layout

def test_g11_4a_newtonian_is_written_at_the_top_level(records):
    rec = sorted(_null_records(records, "newtonian").items())[0][1]
    text = export_transport_properties(rec, RHO, None)
    assert "transportModel  Newtonian;" in text
    # the header COMMENT names NewtonianCoeffs deliberately (see the last test in this
    # file), so what must be absent is the BLOCK, not the word
    assert re.search(r"(?<![A-Za-z0-9_])NewtonianCoeffs\s*\{", text) is None
    of_name, c = parse_transport_properties(text, expect="newtonian")
    assert of_name == "Newtonian"
    assert c["nu"] == pytest.approx(float(rec.parameters["mu"].si) / RHO, rel=1e-15)


def test_g11_4a_a_coeffs_block_is_refused_not_ignored(records):
    """OpenFOAM would never READ a NewtonianCoeffs block (Newtonian.C:71 reads the
    top-level `nu`), so a caller who wrote one believes a number reached the solver that
    did not. The parser refuses; ignoring it is the silent-drop class."""
    rec = sorted(_null_records(records, "newtonian").items())[0][1]
    text = export_transport_properties(rec, RHO, None)
    nu = parse_transport_properties(text)[1]["nu"]
    moved = text.replace(f"nu              {nu:.17g};",
                         "NewtonianCoeffs\n{\n" f"    nu      {nu:.17g};\n" "}")
    assert "NewtonianCoeffs" in moved
    with pytest.raises(ValueError, match="NewtonianCoeffs"):
        parse_transport_properties(moved)
    # and a dictionary with neither is an error, not a default: `nu` has no default
    with pytest.raises(ValueError, match="no top-level `nu` entry"):
        parse_transport_properties("transportModel  Newtonian;\n")


NESTED = """transportModel  Newtonian;
nu              1e-06;
somethingElse
{
    inner { a 1; }
    nu   999;
}
"""


def test_the_newtonian_top_level_strip_is_not_fail_open_on_nesting():
    """REGRESSION, found by code review 2026-09-06 and it was a silently WRONG number.

    The first version stripped sub-dictionaries with `re.sub(r"\\w+\\s*\\{[^}]*\\}", ...)`.
    A regex cannot match balanced braces: `[^}]*` stops at the INNER `}`, so everything
    after a NESTED block leaked out as top level. On `NESTED` below it returned
    `nu = 999` -- and returned it with the real entry deleted too, i.e. it succeeded on a
    dictionary OpenFOAM rejects. `_strip_subdicts` counts braces instead."""
    assert parse_transport_properties(NESTED)[1] == {"nu": 1e-06}
    # with the real entry gone the parser must RAISE, not fall back on the buried one
    without = NESTED.replace("nu              1e-06;\n", "")
    with pytest.raises(ValueError, match="no top-level `nu` entry"):
        parse_transport_properties(without)


def test_the_newtonian_branch_reports_which_kind_of_failure(records):
    """One check, three reasons, and the caller needs to tell them apart: absent (the file
    is incomplete), malformed (the file is wrong), and an OpenFOAM construct this parser
    does not implement (the file is fine, this is not OpenFOAM). The Coeffs branch already
    made that distinction; the Newtonian branch did not until review."""
    head = "transportModel  Newtonian;\n"
    with pytest.raises(ValueError, match="no top-level `nu` entry"):
        parse_transport_properties(head)
    with pytest.raises(ValueError, match="does not implement"):
        parse_transport_properties(head + "nu  $nuWater;\n")
    with pytest.raises(ValueError, match="present but unreadable"):
        parse_transport_properties(head + "nu  1,5e-6;\n")
    # a duplicate is not a silent last-one-wins either
    with pytest.raises(ValueError, match="more than one top-level"):
        parse_transport_properties(head + "nu 1e-06;\nnu 5e-06;\n")
    # unbalanced braces are refused rather than guessed
    with pytest.raises(ValueError, match="unbalanced"):
        parse_transport_properties(head + "nu 1e-06;\nblock\n{\n")


def test_the_newtonian_branch_still_accepts_ordinary_dictionaries(records):
    """The mirror of the above, and the reason the fix is a brace counter rather than a
    blanket residue check: a real `transportProperties` carries top-level entries this
    parser has no opinion about, and refusing them would bite exactly the foreign
    dictionaries `expect=None` exists for."""
    rec = sorted(_null_records(records, "newtonian").items())[0][1]
    text = export_transport_properties(rec, RHO, None)
    nu = parse_transport_properties(text)[1]["nu"]
    for extra in ("rho             rho [1 -3 0 0 0 0 0] 1000;\n",
                  "phases          (water air);\n",
                  "sigma           sigma [1 0 -2 0 0 0 0] 0.07;\n",
                  "HerschelBulkleyCoeffs { k 1; n 0.5; tau0 1; nu0 1; }\n"):
        assert parse_transport_properties(text + extra)[1] == {"nu": nu}, extra


def test_the_strip_is_quote_aware(records):
    """A brace inside a double-quoted value is not a brace. `note "a { brace";` is a
    legal free-text FoamFile entry and OpenFOAM's tokenizer is string-aware, so refusing
    it would be a fail-CLOSED check rejecting valid input its own test set never
    contained -- the exact shape of the 2026-09-06 lesson, in a function written that
    day. Found by review before it shipped."""
    head = "transportModel  Newtonian;\nnu              1e-06;\n"
    assert parse_transport_properties(head + 'note "a { brace";\n')[1] == {"nu": 1e-06}
    assert parse_transport_properties(head + 'note "}}}";\n')[1] == {"nu": 1e-06}
    # a quoted brace inside a SUB-DICTIONARY must not leak the block's contents either
    assert parse_transport_properties(
        head + 'blk\n{\n  note "{";\n  nu 999;\n}\n')[1] == {"nu": 1e-06}
    # an ESCAPED quote does not close the string -- OpenFOAM's tokenizer honours `\"`,
    # so refusing this would be the same false-refusal class one layer down
    assert parse_transport_properties(head + 'note "he said \\" hi";\n')[1] == {"nu": 1e-06}
    # a real FoamFile header entry, which is where quoted values actually occur
    assert parse_transport_properties(
        head + 'arch "LSB;label=32;scalar=64";\n')[1] == {"nu": 1e-06}
    # an unterminated string is refused rather than silently swallowing the rest
    with pytest.raises(ValueError, match="unterminated"):
        parse_transport_properties(head + 'note "oops;\n')
    # ... including one left open by a trailing escape
    with pytest.raises(ValueError, match="unterminated"):
        parse_transport_properties(head + 'note "oops \\";\n')
    # and the unbalanced-brace message must not blame the thing that is now handled
    with pytest.raises(ValueError, match="not the cause"):
        parse_transport_properties(head + "blk\n{\n")


def test_the_docstring_does_not_call_foamDictionary_a_layout_check():
    """`foamDictionary -entry nu` reads a top-level key out of a file this module wrote
    with a top-level key, using a GENERIC dictionary reader. It would pass unchanged if
    OpenFOAM's Newtonian actually read `NewtonianCoeffs/nu` -- both candidates pass it,
    so it discriminates nothing about the layout. The docstring claimed it as a second
    independent check until a review said otherwise; this pins the correction."""
    doc = openfoam_module.__doc__ or ""
    assert "only\nDISCRIMINATING check of that layout is gate G11.4b" in doc
    assert "is NOT a second one and must not be cited as such" in doc


def test_g11_4a_expect_still_discriminates_the_model(records):
    rec = sorted(_null_records(records, "newtonian").items())[0][1]
    text = export_transport_properties(rec, RHO, None)
    with pytest.raises(ValueError, match="found transportModel Newtonian"):
        parse_transport_properties(text, expect="power_law")


# ----------------------------------------------------------- G11.5 discriminating maps

def test_g11_5a_newtonian_wrong_map_forgetting_rho(records):
    """`nu = mu` instead of `mu/rho`. It must miss by more than 1e-3, or G11.1 passing at
    1e-10 is evidence about arithmetic and not about the map."""
    seps = []
    for rid, rec in sorted(_null_records(records, "newtonian").items()):
        gd = _rates(rec)
        nu_rec = np.asarray(rec.eta_at(gd), dtype=float) / RHO
        wrong = nu_openfoam("Newtonian", {"nu": float(rec.parameters["mu"].si)}, gd)
        seps.append(float(np.abs(wrong - nu_rec).max() / nu_rec.max()))
    assert len(seps) == 20
    assert min(seps) > BAR_G115, f"G11.5a UNINFORMATIVE: min separation {min(seps):.3e}"
    # it is a pure factor of rho, so it should land at rho - 1 exactly
    assert min(seps) == pytest.approx(RHO - 1.0, rel=1e-12)


def test_g11_5b_power_law_wrong_map_subtracting_one_twice(records):
    """OpenFOAM's powerLaw already computes `pow(sr, n - 1)`, so mapping biorheo's `n` to
    `n - 1` applies the -1 twice. Both slots take a plain dimensionless scalar, so
    nothing raises -- exactly the shape of the CrossPowerLaw name collision."""
    seps = []
    for rid, rec in sorted(_null_records(records, "power_law").items()):
        gd = _rates(rec)
        eta_max, eta_min = _clamps(rec)
        c = dict(coefficients(rec, RHO, eta_max, eta_min=eta_min))
        nu_rec = np.asarray(rec.eta_at(gd), dtype=float) / RHO
        c["n"] = c["n"] - 1.0
        wrong = nu_openfoam("powerLaw", c, gd)
        seps.append(float(np.abs(wrong - nu_rec).max() / nu_rec.max()))
    assert len(seps) == 20
    assert min(seps) > BAR_G115, f"G11.5b UNINFORMATIVE: min separation {min(seps):.3e}"


def test_the_power_law_floor_is_SMALL_not_VSMALL():
    """powerLaw.C:64-68 floors the DIMENSIONLESS product `1[s]*strainRate()` at
    SMALL = 1e-15. Transcribing VSMALL = 1e-300 would change `nu` at zero rate by 285
    decades of argument, which is what a cell of a run started from rest sees.

    A first version of this docstring said "the five other mapped models use VSMALL",
    which is false and a review caught it: only `HerschelBulkley` and `Casson` do (they
    divide by the rate); `CrossPowerLaw`, `BirdCarreau` and `Newtonian` have no floor at
    all, because nothing in them is singular at zero rate."""
    assert SMALL == 1.0e-15
    c = {"k": 1.0, "n": 0.5, "nuMin": 0.0, "nuMax": 1e300}
    at_zero = float(nu_openfoam("powerLaw", c, np.array([0.0]))[0])
    assert at_zero == pytest.approx(SMALL ** -0.5, rel=1e-12)


# ------------------------------------------------- G11.6 the structural_steady refusal

def test_g11_6_structural_steady_equals_cross_with_m_one():
    """The pre-registered numeric claim that decides what the refusal MESSAGE may say.

    Item 11 said `structural_steady` "is thixotropic and has no counterpart there".
    `STRUCTURAL_STEADY.limits` says the opposite. If this passes, the item's stated
    reason is false and the message must give the true one instead."""
    ss, cr = get_model("structural_steady"), get_model("cross")
    gd = np.logspace(-2, 2, N_RATES)
    worst = 0.0
    for eta_inf, beta, alpha_tr in ((3.5e-3, 28.0714286, 0.2), (1.0, 5.0, 10.0),
                                    (1e-4, 1e3, 1e-3)):
        a = ss(gd, eta_inf, beta, alpha_tr)
        b = cr(gd, eta_inf * (1.0 + beta), eta_inf, alpha_tr, 1.0)
        worst = max(worst, float(np.abs(a - b).max() / np.abs(a).max()))
    assert worst <= BAR_G116, f"G11.6: {worst:.3e} > {BAR_G116:.0e}"
    # measured 2026-09-06; recorded so a later change to either model shows up as a move
    # rather than as a still-passing gate
    assert worst < 1e-15, worst


def test_g11_6_the_separate_family_needs_its_own_substitution():
    """`structural_steady_separate` is refused by the same code path, and the refusal
    tells the caller what to write instead -- but it carries `alpha` and `t_r` SEPARATELY
    and multiplies them, so it has no `alpha_tr` to read. A first version of the message
    gave `structural_steady`'s substitution to both families and a review caught it; this
    anchors the one that is actually right, which the message had no test for at all."""
    sss, cr = STRUCTURAL_STEADY_SEPARATE, get_model("cross")
    assert sss.param_names == ("eta_inf", "beta", "alpha", "t_r")
    assert "alpha_tr" not in sss.param_names
    gd = np.logspace(-2, 2, N_RATES)
    worst = 0.0
    for eta_inf, beta, alpha, t_r in ((3.5e-3, 28.0714286, 0.4, 0.5), (1.0, 5.0, 2.0, 5.0),
                                      (1e-4, 1e3, 1e-2, 0.1)):
        a = sss(gd, eta_inf, beta, alpha, t_r)
        b = cr(gd, eta_inf * (1.0 + beta), eta_inf, alpha * t_r, 1.0)
        worst = max(worst, float(np.abs(a - b).max() / np.abs(a).max()))
    assert worst <= BAR_G116, f"G11.6 (separate): {worst:.3e} > {BAR_G116:.0e}"


def test_g11_6_the_refusal_is_explicit_and_says_the_true_reason(records):
    """A bare NotImplementedError was what item 11 asked to replace. The message must
    (a) not claim no counterpart exists, (b) name the structure equation as the reason,
    (c) give the substitution a caller can make deliberately."""
    rec = records["blood_horner2020_donor_A"]
    # the substitution DIFFERS between the two families; asserting one string for both is
    # what certified a wrong instruction before the review
    want = {STRUCTURAL_STEADY: "lam = alpha_tr, m = 1",
            STRUCTURAL_STEADY_SEPARATE: "lam = alpha*t_r"}
    for model, sub in want.items():
        with pytest.raises(ThixotropyNotRepresentable) as e:
            coefficients(replace(rec, model=model), RHO, 1.0)
        msg = str(e.value)
        assert msg.startswith(model.name), msg[:80]
        assert "NOT because no counterpart exists" in msg
        assert "CrossPowerLaw with n = 1" in msg
        assert "STRUCTURE EQUATION" in msg
        assert "silently wrong in every transient" in msg
        assert sub in msg, (model.name, msg)
    # and the two messages must not be interchangeable
    msgs = []
    for model in want:
        with pytest.raises(ThixotropyNotRepresentable) as e:
            coefficients(replace(rec, model=model), RHO, 1.0)
        msgs.append(str(e.value))
    assert msgs[0] != msgs[1]
    # the separate family's message may NAME alpha_tr (it says explicitly that this
    # family has none) but must never INSTRUCT the caller to use one
    assert "lam = alpha_tr," not in msgs[1] and "alpha_tr gd" not in msgs[1]
    assert "has no alpha_tr to read" in msgs[1]
    assert issubclass(ThixotropyNotRepresentable, NotImplementedError)


def test_the_prereg_for_item_11_predates_its_code():
    """The pre-registration is only worth something if it was written first. Its commit
    is `e9c1ace`, which contains that file and NOTHING else; the code it grades landed in
    `bcf84bf`. Every bar this file asserts must appear in it, and so must the two premises
    it went on to refute -- named, so a later reader cannot quietly reinstate either."""
    prereg = (Path(__file__).resolve().parents[1] / "data" / "openfoam"
              / "newtonian_powerlaw_prereg.md").read_text(encoding="utf-8")
    for token in ("G11.1", "G11.2", "G11.3", "G11.4", "G11.5", "G11.6", "G11.7"):
        assert token in prereg, token
    # the bars, in either spelling -- what matters is that the NUMBER this file asserts
    # is the number the pre-registration named, not how it was typed
    for bar in (BAR_G111, BAR_G115, BAR_G116):
        assert any(s in prereg for s in (f"{bar:g}", f"{bar:.0e}".replace("e-0", "e-"),
                                         repr(bar))), bar
    # premise 1, refuted by G11.6: item 11 said structural_steady has no counterpart
    assert "item's stated reason is wrong" in prereg
    # premise 2: the item told us to assume nothing about powerLaw's limiters, and the
    # answer -- that they are REQUIRED, unlike BirdCarreau's `a` -- is recorded there
    assert "REQUIRED, not optional" in prereg
    # and the one factual error the prereg itself contains, which is corrected in the
    # addendum rather than edited out (a pre-registration is not rewritten)
    assert "The other four models this bridge already covers use `VSMALL`" in prereg


def test_the_gate_report_on_disk_matches_what_the_code_computes():
    """`data/openfoam/newtonian_powerlaw_gates.md` is written by
    `scripts/item11_gate_report.py`. This READS it and re-derives the numbers; it does not
    re-run the script, so a failure here cannot overwrite the artifact it is checking --
    three tests in this repo used to do exactly that."""
    md = (Path(__file__).resolve().parents[1] / "data" / "openfoam"
          / "newtonian_powerlaw_gates.md").read_text(encoding="utf-8")
    # the verdicts, and the two constants the report identifies rather than just prints
    assert md.count("| PASS |") == 6, md
    assert "FAILED" not in md and "UNINFORMATIVE |" not in md
    assert f"`rho - 1` = {RHO - 1.0:.6g}" in md
    assert "115 records, of which 20 `newtonian` and 20 `power_law`" in md
    assert f"{BAR_G111:g}" in md and f"{BAR_G116:g}" in md and f"{BAR_G115:g}" in md


def test_the_exported_file_says_the_newtonian_layout_is_deliberate(records):
    """A reader who knows the other five families have a Coeffs block will think a
    missing one is a bug unless the file says otherwise."""
    rec = sorted(_null_records(records, "newtonian").items())[0][1]
    text = export_transport_properties(rec, RHO, None)
    assert "deliberately no NewtonianCoeffs" in text
    assert "Newtonian.C:71" in text
    # and the power-law file must say both limits are required and why
    pl = sorted(_null_records(records, "power_law").items())[0][1]
    eta_max, eta_min = _clamps(pl)
    t = export_transport_properties(pl, RHO, eta_max, eta_min=eta_min)
    assert "BOTH required" in t and "no plateau at either end" in t
    assert "nuMin" in t and "nuMax" in t
