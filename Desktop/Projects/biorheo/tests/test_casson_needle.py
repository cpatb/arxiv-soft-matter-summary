"""Casson needle flow and the blood injectability case study (NEXT_STEPS item 4).

Anchors outside the code under test:
  * `flow_rate_gnf`, the package's numerical quadrature, which inverts the flow curve
    node by node and knows nothing about the closed form;
  * the `tau_y -> 0` Poiseuille limit, computed from `pi R^4 dp / (8 eta L)` directly;
  * hand-evaluated coefficients of the bracket at two values of x;
  * `derivations/casson_tube_flow.py`, which derives the form in SymPy and exits
    non-zero on any failed identity.

Every number the case study's `results.md` states is checked against its own JSON.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from biorheo.apps.injectability import (
    Needle, flow_rate, flow_rate_casson, flow_rate_gnf, wall_stress,
)
from biorheo.constitutive.gnf import CASSON
from biorheo.registry.loader import load_all

ROOT = Path(__file__).resolve().parent.parent
CASE = ROOT / "case_studies" / "blood_injectability"

R, L = 1.0e-4, 2.0e-2
RTOL = 1e-8          # the gate pre-registered in derivations/casson_tube_flow.py


@pytest.fixture(scope="module")
def blood():
    return {k: v for k, v in sorted(load_all().items()) if v.model.name == "casson"}


@pytest.fixture(scope="module")
def case():
    p = CASE / "results.json"
    assert p.exists(), "run case_studies/blood_injectability/run.py"
    return json.loads(p.read_text())


def test_closed_form_matches_the_independent_quadrature(blood):
    """C5: 20 shipped records plus four corners of a (tau_y, eta_c) box."""
    worst = 0.0
    cases = [(float(r.params_si[0]), float(r.params_si[1])) for r in blood.values()]
    cases += [(ty, ec) for ty in (1e-3, 5.0) for ec in (1e-3, 1e-1)]
    for tau_y, eta_c in cases:
        for dp in (1e5, 3e5, 1e6):
            q_c = flow_rate_casson(dp, R, L, tau_y, eta_c)
            q_n = flow_rate_gnf(dp, R, L, CASSON, np.array([tau_y, eta_c]))
            rel = abs(q_c - q_n) / q_n
            worst = max(worst, rel)
    assert worst <= RTOL, worst
    assert worst > 0     # not trivially identical: different algorithms


def test_poiseuille_and_yield_limits():
    eta = 3.5e-3
    dp = 2.0e5
    # tau_y -> 0 is Poiseuille with viscosity eta_c, exactly
    assert np.isclose(flow_rate_casson(dp, R, L, 0.0, eta),
                      np.pi * R ** 4 * dp / (8 * eta * L), rtol=1e-15, atol=0)
    # at yield the flow stops CONTINUOUSLY and CUBICALLY: x = 1 is a TRIPLE root of the
    # bracket (value, slope AND curvature vanish), so Q ~ (1 - x)^3. This test is what
    # corrected the derivation, which first claimed a double root and a quadratic
    # approach: the measured ratio is 8, not 4.
    tw = wall_stress(dp, R, L)
    assert flow_rate_casson(dp, R, L, tw, eta) == 0.0
    assert flow_rate_casson(dp, R, L, 1.5 * tw, eta) == 0.0
    q = np.array([flow_rate_casson(dp, R, L, (1 - e) * tw, eta) for e in (1e-3, 2e-3)])
    # doubling (1 - x) multiplies Q by 8 (cubic), not 4 (quadratic) or 2 (linear)
    assert 7.9 < q[1] / q[0] < 8.1, q[1] / q[0]
    # and the leading coefficient is the one the Taylor expansion gives, -b'''(1)/6 = 1/3
    eps = 1e-4
    lead = flow_rate_casson(dp, R, L, (1 - eps) * tw, eta) / (
        np.pi * R ** 3 * tw / (4 * eta) * eps ** 3)
    assert np.isclose(lead, 1 / 3, rtol=1e-3), lead


def test_hand_evaluated_bracket():
    """The bracket 1 - (16/7)sqrt(x) + (4/3)x - x^4/21 at two hand-worked x."""
    dp, eta_c = 4.0e5, 4.0e-3
    tw = wall_stress(dp, R, L)                 # 4e5 * 1e-4 / 4e-2 = 1000 Pa
    assert tw == 1000.0
    # x = 0.01: sqrt(x) = 0.1 -> 1 - 1.6/7 + 0.013333... - 1e-8/21
    x = 0.01
    hand = 1 - 16 / 7 * 0.1 + 4 / 3 * 0.01 - 0.01 ** 4 / 21
    q = flow_rate_casson(dp, R, L, x * tw, eta_c)
    assert np.isclose(q, np.pi * R ** 3 * tw / (4 * eta_c) * hand, rtol=1e-14)
    # x = 1/4: sqrt(x) = 1/2 -> 1 - 8/7 + 1/3 - 1/(21*256)
    hand4 = 1 - 8 / 7 + 1 / 3 - 1 / (21 * 256)
    q4 = flow_rate_casson(dp, R, L, 0.25 * tw, eta_c)
    assert np.isclose(q4, np.pi * R ** 3 * tw / (4 * eta_c) * hand4, rtol=1e-14)
    # the exact rational: 1 - 8/7 + 1/3 - 1/5376 = 4/21 - 1/5376
    assert np.isclose(hand4, 4 / 21 - 1 / 5376, rtol=0, atol=1e-15)


def test_flow_rate_dispatches_casson_to_the_closed_form(blood):
    """`flow_rate` must route Casson to the closed form, not silently to quadrature."""
    rec = blood["blood_horner2020_donor_A"]
    nd = Needle(radius=R, length=L)
    got = flow_rate(3e5, nd, rec.model, rec.params_si)
    assert got == flow_rate_casson(3e5, R, L, *(float(v) for v in rec.params_si))


def test_the_form_stays_computable_all_the_way_to_yield():
    """The module evaluates the FACTORED bracket. The expanded one is a difference of
    O(1) terms giving an O((1-x)^3) result and is meaningless within about 1e-5 of yield
    (relative error 0.52 there, 9e5 at 1e-7). The factored one has all-positive
    coefficients after the triple root is divided out, so it cannot cancel."""
    dp, eta_c = 4.0e5, 4.0e-3
    tw = wall_stress(dp, R, L)
    pref = np.pi * R ** 3 * tw / (4.0 * eta_c)
    macheps = float(np.finfo(float).eps)
    for eps in (1e-4, 1e-5, 1e-6, 1e-7, 1e-9):
        q = flow_rate_casson(dp, R, L, (1 - eps) * tw, eta_c)
        assert q > 0, eps                                   # never negative, never zero
        # Two error sources, both bounded, neither round-off inside the bracket:
        #   (a) the next term of the series, O(eps) relative;
        #   (b) the CALLER's conditioning -- `(1 - eps) * tw` is not exact, so x carries
        #       ~macheps absolute, 1 - x carries macheps/eps relative, and the bracket
        #       cubes it. No rewriting of the formula can remove (b); it is in the inputs.
        bound = 2 * eps + 4 * macheps / eps
        assert abs(q / (pref * eps ** 3 / 3) - 1) < bound, (eps, q, bound)
        # the expanded form is what this replaced; show it really is unusable there
        if eps <= 1e-5:
            u = np.sqrt(1 - eps)
            expanded = 1 - 16 / 7 * u + 4 / 3 * (1 - eps) - (1 - eps) ** 4 / 21
            assert abs(expanded / (eps ** 3 / 3) - 1) > 0.4, eps
    # monotone in dp all the way down, which pressure_for_flow_rate's brentq relies on
    q = [flow_rate_casson(dp, R, L, (1 - e) * tw, eta_c) for e in (1e-7, 1e-6, 1e-5, 1e-4)]
    assert all(a < b for a, b in zip(q, q[1:]))


def test_rejects_bad_inputs():
    for bad in [(0.0, R, L, 1e-3, 1e-3), (1e5, 0.0, L, 1e-3, 1e-3),
                (1e5, R, 0.0, 1e-3, 1e-3), (1e5, R, L, 1e-3, 0.0),
                (1e5, R, L, -1.0, 1e-3)]:
        with pytest.raises(ValueError):
            flow_rate_casson(*bad)


def test_the_derivation_still_holds():
    r = subprocess.run([sys.executable, str(ROOT / "derivations" / "casson_tube_flow.py")],
                       capture_output=True, text=True, cwd=str(ROOT), timeout=900)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-2000:]
    assert "0 failure(s)" in r.stdout and "FAIL" not in r.stdout


# --------------------------------------------------------- the case study's own claims

def test_hand_computed_anchor_for_the_case_study(case):
    """An anchor for the case study that does not pass through its own code. Donor A is
    effectively Newtonian here (x = 3.06e-5), so Poiseuille gives the pressure to within
    the bracket: dp = 8 eta_c L Q / (pi R^4) / bracket(x). Worked with the shipped
    numbers, this must reproduce results.json's dp for donor A."""
    row = next(r for r in case["rows"] if r["donor"] == "A")
    R_n = case["needle"]["id_nominal_mm"] * 1e-3 / 2.0
    L_n = case["needle"]["length_m"]
    Q = case["Q_si"]
    dp_poiseuille = 8 * row["eta_c"] * L_n * Q / (np.pi * R_n ** 4)
    u_ = np.sqrt(row["x"])
    bracket = ((1 - row["x"]) / (1 + u_)) ** 3 * (
        21 + u_ * (15 + u_ * (10 + u_ * (6 + u_ * (3 + u_))))) / 21
    assert np.isclose(dp_poiseuille / bracket, row["dp_Pa"], rtol=1e-9), (
        dp_poiseuille / bracket, row["dp_Pa"])
    # and the force is that pressure on the plunger area, nothing else
    area = np.pi * case["barrel_radius_m"] ** 2
    assert np.isclose(row["dp_Pa"] * area, row["force_N"], rtol=1e-12)


def test_case_study_answers_the_pre_registered_question(case):
    """20/20 extrapolated at a 27G needle, 1 mL/min. This is the question NEXT_STEPS
    item 4 asked, and the prediction written before the run."""
    assert case["n_records"] == 20
    assert case["n_extrapolated"] == 20
    assert all(r["extrapolated"] for r in case["rows"])
    assert case["decades_past_min"] > 1.4          # 1.42, i.e. 26x past the envelope
    for r in case["rows"]:
        assert r["gd_wall"] > r["envelope"][1], r["donor"]


def test_the_extrapolation_distance_is_a_device_property_not_a_donor_one(case):
    """The wall rate is 4Q/(pi R^3) for a Newtonian fluid, independent of viscosity, and
    these donors are effectively Newtonian here -- so every donor sits the same distance
    past its envelope. That is what makes it computable at design time."""
    assert case["wall_rate_spread"] < 1.01
    gd = np.array([r["gd_wall"] for r in case["rows"]])
    app = np.array([r["gd_apparent"] for r in case["rows"]])
    assert np.all(gd > app)                        # shear-thinning: true rate is higher
    assert np.all((gd - app) / app < 0.05)         # but only just, because x is tiny


def test_the_yield_stress_is_the_parameter_this_device_cannot_see(case):
    assert case["x_max"] < 1e-3
    assert case["yield_share_max"] < 0.05          # zeroing tau_y moves F by < 5 %
    assert case["spearman_force_eta_c"] > 0.95
    assert case["spearman_force_tau_y"] < 0.5
    # the ordering is the claim: eta_c ranks with the force, tau_y does not
    assert case["spearman_force_eta_c"] > 4 * case["spearman_force_tau_y"]


def test_the_geometry_claim_does_not_depend_on_the_unverified_standard(case):
    """RETRACTED 2026-09-05 and replaced. The original claim -- "the ISO 9626 tolerance
    for regular-wall 27G moves the force more than the donor spread" -- rested on a band
    (0.184-0.241 mm) that mixes a regular-wall figure with a THIN-WALL one, and it
    REVERSES under the narrower reading (1.69x against a 2.28x donor spread). ISO
    9626:2016 is not available offline to settle it. The replacement is exact arithmetic
    on the fitted records: F ~ R^-4, so a stated percentage error in bore reproduces a
    stated force spread, with no standard involved. This pins the replacement, and pins
    that the retraction is still on the page."""
    g, rows = case["geometry"], case["rows"]
    # the measured scaling, close to Poiseuille's -4 but not identical (Casson bracket)
    assert -4.0 < g["slope_dlogF_dlogR"] < -3.9

    # the replacement claim, recomputed here from the row forces rather than trusted
    F = np.array([r["force_N"] for r in rows])
    assert np.isclose(g["spread_all_20"], F.max() / F.min(), rtol=1e-12)
    keep = np.array([r["donor"] not in ("F", "K", "N") for r in rows])
    assert np.isclose(g["spread_excluding_paper_outliers"],
                      F[keep].max() / F[keep].min(), rtol=1e-12)
    # excluding the source study's own flagged donors SHRINKS the spread, which makes
    # the geometry point stronger rather than weaker -- that direction is the claim
    assert g["spread_excluding_paper_outliers"] < g["spread_all_20"]
    for key, spread in (("id_error_matching_full_spread", g["spread_all_20"]),
                        ("id_error_matching_typical_spread",
                         g["spread_excluding_paper_outliers"])):
        assert np.isclose((1 + g[key]) ** abs(g["slope_dlogF_dlogR"]), spread, rtol=1e-9)
    # a ten-percent bore error covers the typical donor spread: the headline
    assert g["id_error_matching_typical_spread"] < 0.15

    # the two readings of the band genuinely disagree, which is WHY it was retracted:
    # the wide one beats the donor spread and the narrow one loses to it
    assert g["band_wide_force_ratio"] > g["spread_all_20"] > g["band_narrow_force_ratio"]

    # the illustrative donor is the median by force, not the flagged outlier N
    assert g["median_donor"] not in ("F", "K", "N")
    assert np.isclose([r["force_N"] for r in rows if r["donor"] == g["median_donor"]][0],
                      sorted(F)[len(F) // 2], rtol=1e-12)

    # every record stays extrapolated at BOTH ends -- checked on 20, not inferred from one
    for tag in ("low", "high"):
        assert case["band_extrapolation"][tag]["n_extrapolated"] == 20, tag
        assert case["band_extrapolation"][tag]["min_decades_past"] > 1.0, tag


def test_the_neglected_entrance_loss_is_ranked_against_the_yield_effect(case):
    """"Re is small" compares the entrance term with dp, not with the effect being
    ranked. Measured, the Hagenbach loss is the same ORDER as the yield-stress
    contribution the study calls negligible, so the study resolves neither -- and the
    two corrections have opposite signs, so the "upper bound" framing (which counted
    only the downward Fahraeus-Lindqvist one) is withdrawn."""
    e = case["entrance_loss"]
    assert 0.0 < e["fraction_of_dp_min"] < e["fraction_of_dp_max"] < 0.02
    assert 0.1 < e["fraction_of_dp_max"] / case["yield_share_max"] < 10
    assert 0.1 < e["fraction_of_dp_min"] / case["yield_share_min"] < 10
    md = (CASE / "results.md").read_text(encoding="utf-8")
    assert "do not point the same way" in md and "withdrawn" in md
    assert "[UNVERIFIED" in md          # the Fahraeus-Lindqvist threshold is not cited


def test_the_device_not_patient_claim_carries_its_scope(case):
    """The 3-significant-figure constancy of the extrapolation distance is a property of
    the CASSON FORM, which saturates to a constant eta_c and so cannot shear-thin further
    above yield; a power-law fluid would give an O(1) factor (3n+1)/(4n). The decade-level
    claim survives, the third figure does not, and the report must say so."""
    md = (CASE / "results.md").read_text(encoding="utf-8")
    assert "needs its scope stated" in md and "(3n+1)/(4n)" in md
    # the second mechanism: the envelope is the protocol, identical for all 20 donors
    assert len({tuple(r["envelope"]) for r in case["rows"]}) == 1
    assert "identical for all 20 donors" in md


def test_results_md_quotes_its_own_json(case):
    md = (CASE / "results.md").read_text(encoding="utf-8")
    assert f"**{case['n_extrapolated']} of {case['n_records']} records" in md
    assert f"{case['decades_past_min']:.2f} decades" in md
    assert f"{case['spearman_force_eta_c']:.3f}" in md
    assert f"{case['spearman_force_tau_y']:.3f}" in md
    g = case["geometry"]
    assert f"{g['slope_dlogF_dlogR']:.4f}" in md
    assert f"{100 * g['id_error_matching_full_spread']:.1f} %" in md
    assert f"{g['spread_excluding_paper_outliers']:.2f}x" in md
    assert f"{g['band_narrow_force_ratio']:.2f}x" in md
    # the retraction must stay on the page, not be quietly dropped on a regeneration
    assert "What was retracted here" in md and "That is withdrawn" in md
    # the Fahraeus-Lindqvist caveat is not optional: this lumen is inside that regime
    assert "Fahraeus-Lindqvist" in md
    assert "upper bounds on the real force" not in md   # that framing was withdrawn
