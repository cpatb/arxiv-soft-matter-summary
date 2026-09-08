"""The deliverable page: `biorheo.apps.report`.

This module computes nothing new, so the tests are not about physics. They are about the
three things a page a client pays for has to get right:

1. **The refusals fire, and they change the HEADLINE rather than being appended.** Each
   flag has a test that asserts the wording of the first sentence, because a caveat in a
   footer is a caveat nobody reads.
2. **Nothing is printed that the data does not support** -- an interval that ran into a
   box bound is never worded as a confidence limit, and a record without a sidecar never
   gets an interval at all.
3. **The page cannot drift from the module it quotes.** The omissions list is asserted
   verbatim against `injectability.py`'s own docstring. This repo has shipped a page whose
   prose contradicted its own artifact more than once; the fix each time was a test that
   makes the two impossible to separate.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from biorheo.apps.report import DeviceSpec, OMISSIONS, build, render
from biorheo.registry import load_record

REPO = Path(__file__).resolve().parent.parent
RECORDS = REPO / "src" / "biorheo" / "registry" / "records"
INJECTABILITY = REPO / "src" / "biorheo" / "apps" / "injectability.py"

SPEC = DeviceSpec(needle_radius=1.05e-4, needle_length=12.7e-3, barrel_radius=2.3e-3,
                  flow_rate=1e-6 / 60, density=1059.0, force_budget=25.0,
                  label="test device")

CLOSED = "blood_horner2020_donor_A"              # every interval closed, both measurable
OPEN = "blood_horner2020_donor_A_carreau_yasuda"  # lam, a, n run into a bound


@pytest.fixture(scope="module")
def closed():
    return build(load_record(RECORDS / f"{CLOSED}.yaml"), SPEC)


@pytest.fixture(scope="module")
def open_iv():
    return build(load_record(RECORDS / f"{OPEN}.yaml"), SPEC)


# --- the omissions cannot drift from the module that owns them -----------------------

def test_every_omission_appears_verbatim_in_the_injectability_docstring():
    """If someone edits one module's list the other must fail, not silently disagree."""
    # whitespace-normalised: the docstring line-wraps the kinetic-term phrase across two
    # lines, and a line wrap is not a difference of wording.
    src = " ".join(INJECTABILITY.read_text(encoding="utf-8").split())
    for phrase in OMISSIONS:
        assert phrase in src, phrase


def test_the_rendered_page_lists_every_omission(closed):
    page = render(closed)
    for phrase in OMISSIONS:
        assert phrase in page, phrase


def test_the_page_says_wall_slip_is_not_covered_by_any_flag(closed):
    """Wall slip is a boundary condition, not a constitutive extrapolation, so the
    `extrapolated` flag does not cover it and the page has to say so."""
    page = render(closed)
    assert "wall slip" in page.lower()
    assert "boundary condition" in page


# --- refusals change the headline -----------------------------------------------------

def test_open_interval_is_worded_as_a_bound_not_a_limit(open_iv):
    assert "open_interval" in open_iv.flags
    head = render(open_iv).split("## The answer")[1].split("\n\n")[1]
    assert "BOUND, not a limit" in head
    for name in ("lam", "a", "n"):
        assert name in head, name


def test_a_closed_record_is_worded_as_a_hull_and_names_which_interval(closed):
    assert "open_interval" not in closed.flags
    head = render(closed).split("## The answer")[1].split("\n\n")[1]
    assert "one-at-a-time parameter hull" in head
    assert "BOUND, not a limit" not in head


def test_the_hull_says_it_is_wider_than_the_joint_region(closed):
    """The repo rule: name which interval you are quoting, every time."""
    assert closed.hull["kind"] == "one-at-a-time profile push-through hull"
    assert "inscribed" in closed.hull["wider_than"]


def test_no_sidecar_gives_a_point_estimate_and_says_so(tmp_path):
    """A record with its identifiability block removed must NOT quietly lose its error
    bar: the headline has to say the number is a point estimate."""
    dst = tmp_path / f"{CLOSED}.yaml"
    shutil.copy2(RECORDS / f"{CLOSED}.yaml", dst)      # sidecar deliberately not copied
    d = build(load_record(dst), SPEC)
    assert "no_sidecar" in d.flags
    assert d.hull["relative_width"] == 0.0
    head = render(d).split("## The answer")[1].split("\n\n")[1]
    assert "POINT ESTIMATE" in head and "no uncertainty" in head
    assert "hull" not in head


def test_a_non_commercial_record_is_stamped_before_the_answer(closed):
    """The banner must come BEFORE the headline; a licence problem discovered at the foot
    of a page has already been quoted from."""
    page = render(closed)
    assert "non_commercial" in closed.flags
    assert page.index("NOT LICENSABLE") < page.index("## The answer")


def test_every_registry_record_currently_triggers_the_licence_flag():
    """A standing venture constraint, pinned so it fails the day it stops being true --
    which is the day the venture can ship a paid study off its own registry."""
    for rid in (CLOSED, "carbopol_etd2020_0p5wt_40uL"):
        d = build(load_record(RECORDS / f"{rid}.yaml"), SPEC)
        assert "non_commercial" in d.flags, rid


def test_extrapolation_is_flagged_with_its_distance_in_decades(closed):
    assert "extrapolated" in closed.flags
    assert closed.envelope["decades_past"] == pytest.approx(1.42, abs=0.02)
    page = render(closed)
    assert "decades outside the measured range" in page


def test_a_device_inside_the_envelope_is_not_flagged():
    """The discriminating half: a routine that flagged everything would pass the test
    above and mean nothing."""
    slow = DeviceSpec(needle_radius=1.05e-4, needle_length=12.7e-3, barrel_radius=2.3e-3,
                      flow_rate=5e-11, label="slow")
    d = build(load_record(RECORDS / f"{CLOSED}.yaml"), slow)
    assert "extrapolated" not in d.flags
    assert d.envelope["decades_past"] == 0.0
    assert "operates **inside** the measured range" in render(d)


# --- the content a buyer is paying for ------------------------------------------------

def test_the_page_names_the_parameter_the_force_does_not_depend_on(closed):
    assert closed.sensitivity["insensitive_to"] == ("tau_y",)
    assert closed.sensitivity["dominant"] == "eta_c"
    page = render(closed)
    assert "insensitive to `tau_y`" in page
    assert "will not move this answer" in page


def test_the_insensitivity_claim_has_a_discriminating_half(closed):
    """A sensitivity routine that returned zero for everything would pass the test above.
    The dominant parameter must also be near 1, which such a routine cannot produce."""
    assert abs(closed.sensitivity["d_lnF_d_lnp"]["eta_c"]) > 0.9


def test_the_bore_section_compares_the_two_error_terms(closed):
    page = render(closed)
    assert "d ln F / d ln R" in page
    assert "the parameter uncertainty is the larger term" in page
    assert closed.bore["newtonian_rule_safe_within_10pct"] is True


def test_the_derating_section_quotes_the_ratio_not_the_absolute(closed):
    page = render(closed)
    assert "Derating factor" in page
    assert "**ratio** is the usable number" in page
    assert closed.derating is not None and closed.derating["factor"] < 1.0


def test_no_gauge_is_named_and_the_page_says_why(closed):
    page = render(closed)
    for banned in ("27G", "30G", "27 G", "ISO 9626"):
        assert banned not in page, banned
    assert "A gauge number is a range, not a dimension" in page


def test_provenance_carries_every_doi(closed):
    page = render(closed)
    rec = load_record(RECORDS / f"{CLOSED}.yaml")
    assert len(closed.provenance["sources"]) == len(rec.sources)
    for s in rec.sources:
        assert f"https://doi.org/{s.doi}" in page


UNREACHABLE = "blood_horner2020_donor_B_carreau_yasuda"   # eta_inf at its 1e6 Pa s bound


def test_an_unreachable_interval_end_opens_the_hull_rather_than_being_dropped():
    """A parameter end at which the fluid cannot be delivered at any attainable pressure
    has no force. Silently dropping it would NARROW the interval -- the direction that
    flatters the answer -- so it opens the hull and says so in the HEADLINE.

    Pinned to the one record in the registry that does this rather than guarded by a skip:
    donor B's Carreau-Yasuda fit puts `eta_inf` on its 1e6 Pa s box bound, and at that end
    a 1 mL/min delivery would need more than 1e12 Pa.
    """
    d = build(load_record(RECORDS / f"{UNREACHABLE}.yaml"), SPEC)
    assert "hull_unbounded_above" in d.flags
    assert any("eta_inf" in u for u in d.hull["unreachable"])
    head = render(d).split("## The answer")[1].split("\n\n")[1]
    assert "OPEN" in head and "lower bound on the uncertainty" in head


def test_the_clamp_branch_has_never_fired_on_this_registry(closed, open_iv):
    """**This is a record of untested code, not a passing test of the clamp.**

    `build` clamps an interval end that lies outside the model's own parameter box before
    evaluating a force, because a value one ULP past a bound raises inside the model. No
    record in the registry currently does that, so the clamp's body has never run. It is
    kept because the cost is one `np.clip` and the failure it prevents is an exception in
    a client deliverable -- but it is not verified, and this test says which of those two
    statements is true. If it ever fails, the branch has started firing and needs a real
    test rather than this one.
    """
    for d in (closed, open_iv):
        assert d.hull["clamped"] == (), d.record_id
        assert "interval_end_outside_model_bound" not in d.flags


def test_build_raises_nothing_on_a_record_of_every_family():
    """A record the report cannot fully describe must produce FLAGS, not an exception: a
    partial answer with its limits stated is the deliverable. One record per family here;
    the whole registry is the slow test below."""
    for rid in (CLOSED, OPEN, "blood_horner2020_donor_A_cross",
                "blood_horner2020_donor_A_newtonian", "blood_horner2020_donor_A_power_law",
                "carbopol_etd2020_0p9wt_50uL"):
        d = build(load_record(RECORDS / f"{rid}.yaml"), SPEC)
        assert isinstance(render(d), str), rid


@pytest.mark.slow
def test_build_raises_nothing_on_any_registry_record():
    """The same statement over all 115 records. Slow (~8 min) because every push-through
    end inverts the pressure numerically."""
    for path in sorted(RECORDS.glob("*.yaml")):
        if path.name.endswith(".identifiability.yaml"):
            continue
        d = build(load_record(path), SPEC)
        assert isinstance(render(d), str), path.name
