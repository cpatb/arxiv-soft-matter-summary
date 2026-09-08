"""Item 22: the registry's first commercially licensable records.

Gate pre-registered in `data/dmem_record_prereg.md`. The load-bearing test here is
`test_g22_1p_...`: it is the only place in this repo where a fitted parameter is graded
against a number **published by someone else before this code existed**. Every other
anchor compares biorheo against biorheo, against a closed form, or against a solver.

G22.5 and G22.6 are computed here rather than by the generator, because they are about how
the rest of the package behaves once these records exist, not about writing them.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from biorheo.apps.report import DeviceSpec, build, render
from biorheo.constitutive.gnf import POWER_LAW
from biorheo.constitutive.shift import collapse_shift
from biorheo.registry import load_all, load_record

REPO = Path(__file__).resolve().parent.parent
RECORDS = REPO / "src" / "biorheo" / "registry" / "records"
JSON = REPO / "data" / "dmem_sints2024" / "records.json"
MD = REPO / "data" / "dmem_sints2024" / "records.md"
PREREG = REPO / "data" / "dmem_record_prereg.md"

IDS = (
    "dmem_sints2024_10pct_20C_a", "dmem_sints2024_10pct_20C_b",
    "dmem_sints2024_5pct_20C_a", "dmem_sints2024_5pct_20C_b",
    "dmem_sints2024_1pct_20C_a", "dmem_sints2024_1pct_20C_b",
    "dmem_sints2024_10pct_37C_a", "dmem_sints2024_10pct_37C_b",
)
G221_TOL = 0.05


@pytest.fixture(scope="module")
def art() -> dict:
    assert JSON.exists(), "run scripts/records_from_dmem.py"
    return json.loads(JSON.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def recs() -> dict:
    return {i: load_record(RECORDS / f"{i}.yaml") for i in IDS}


# --- the external anchor ---------------------------------------------------------------

def test_g22_1p_the_pooled_fit_reproduces_the_papers_published_flow_indices(art):
    """**The strongest external anchor in this repo.** The authors of arXiv:2409.13650 fit
    the same Ostwald-de Waele power law and published their flow indices in 2024. This code
    has never seen those numbers except as constants in the generator's SPEC table, and the
    fit is done from the raw deposited files."""
    g = art["gates"]["G22_1p_pooled_fit_reproduces_the_papers_flow_index"]
    assert g["pass"]
    assert g["n_graded"] == 4
    assert g["worst_abs_dev"] < G221_TOL
    # and the discriminating half: a fitter that returned a constant would not land within
    # 0.03 of four DIFFERENT published values spanning 0.31 to 0.54
    ns = [r["n_pooled"] for r in g["rows"]]
    assert max(ns) - min(ns) > 0.2, ns


def test_g22_1_is_recorded_failed_and_its_bar_is_not_moved(art):
    """It compared each INTERVAL against a per-CONDITION published value. The bar and its
    number stay on the page; G22.1' is the stricter test beside it."""
    g = art["gates"]["G22_1_reproduces_the_papers_flow_index"]
    assert g["pass"] is False
    assert g["worst_abs_dev"] > G221_TOL
    # the evidence that it was MIS-specified rather than failed: no fit could have met it
    worst_spread = max(x["spread"] for x in g["per_condition_spread"])
    assert worst_spread > G221_TOL, worst_spread
    assert "FAILED" in MD.read_text(encoding="utf-8")


def test_the_papers_indices_are_carried_on_the_records_themselves(recs):
    """A published value used as a bar has to be visible on the artifact it graded, or a
    later reader cannot check the comparison."""
    for rid, rec in recs.items():
        assert "paper_reported_n" in rec.extra, rid
        if rec.extra["paper_reported_n"] is not None:
            assert 0.3 <= float(rec.extra["paper_reported_n"]) <= 0.55


# --- what these records are ---------------------------------------------------------

def test_there_are_exactly_eight_and_they_are_all_power_law(recs):
    assert len(recs) == 8
    assert all(r.model.name == "power_law" for r in recs.values())


# Item 24 added two more CC BY 4.0 records. This constant exists so that the licensable set
# is asserted EXACTLY rather than as a subset -- a subset assertion would let a future record
# become licensable without anyone noticing, and "which records may be shown in public" is
# the one property of this registry a licence mistake makes expensive.
OTHER_LICENSABLE = ("carbopol_pg_vinutha2024_10C", "carbopol_pg_vinutha2024_50C")


def test_they_are_the_first_licensable_records_and_the_set_is_exact():
    """The venture constraint, restated as it now stands.

    Item 18 found all 115 records non-commercial; the eight DMEM records were the first
    that were not, and item 24 added two more. This test asserted `ok == set(IDS)` -- "the
    first AND ONLY" -- and item 24 broke it loudly, which is the intended behaviour: it is a
    sixth consumer of the registry's contents, and item 24's own pre-registration missed it
    because the grep looked for record COUNTS and for `load_all`, while this one selects on
    `commercial_ok`.
    """
    all_recs = load_all()
    ok = {k for k, v in all_recs.items() if v.commercial_ok}
    assert ok == set(IDS) | set(OTHER_LICENSABLE)
    # 115 frozen + 8 DMEM (item 22) + 2 Carbopol-in-propylene-glycol (item 24)
    assert len(all_recs) == 125
    assert all(load_record(RECORDS / f"{i}.yaml").licence.name == "CC BY 4.0"
               for i in (*IDS, *OTHER_LICENSABLE))
    # and nothing may quietly relicense an existing record
    assert len(all_recs) - len(ok) == 115


def test_the_temperature_is_the_papers_number_not_a_guess(recs):
    """'Room temperature' is 20 C in the source, stated in both Results and Methods. It was
    the hard blocker on building these records at all."""
    temps = sorted({r.envelope.temperature_C for r in recs.values()})
    assert temps == [20.0, 37.0]


def test_every_record_carries_the_pre_shear_protocol(recs):
    """2000 s at 1 1/s before measurement, 160 s per point. These are post-pre-shear
    plateau viscosities and a reader who misses that is reading a different measurement."""
    for rid, rec in recs.items():
        assert "2000 s" in rec.notes and "160 s" in rec.notes, rid
        assert "2000 s pre-shear at 1 1/s" in rec.envelope.geometry, rid


def test_every_record_refuses_to_carry_an_injectability_claim(recs):
    for rid, rec in recs.items():
        assert "CANNOT CARRY AN INJECTABILITY CLAIM" in rec.notes, rid
        assert rec.envelope.gammadot_max <= 10.0, rid


def test_both_source_dois_are_on_every_record(recs):
    for rid, rec in recs.items():
        dois = {s.doi for s in rec.sources}
        assert dois == {"10.5281/zenodo.13828992", "10.48550/arXiv.2409.13650"}, rid


def test_the_single_exclusion_is_recorded_on_the_record_that_has_it(recs):
    """Exactly one point in the whole dataset is dropped, and only the record it came from
    may report a non-zero count."""
    dropped = {rid: rec.extra["n_points_dropped_below_0.001_per_s"]
               for rid, rec in recs.items()}
    assert sum(dropped.values()) == 1, dropped
    assert dropped["dmem_sints2024_10pct_20C_a"] == 1


def test_each_record_can_be_evaluated_across_its_own_envelope(recs):
    for rid, rec in recs.items():
        gd = np.geomspace(rec.envelope.gammadot_min, rec.envelope.gammadot_max, 20)
        eta = rec.eta_at(gd)
        assert np.all(np.isfinite(eta)) and np.all(eta > 0), rid
        assert eta[0] > eta[-1], (rid, "must be shear thinning: n < 1")


# --- G22.5 and G22.6, computed here ------------------------------------------------------

def test_g22_5_a_deliverable_page_is_not_stamped_non_licensable_but_IS_extrapolated(recs):
    """The point of the whole item: a page built on one of these may be shown in public.
    And the second half is what keeps it honest -- the same page must say, above the fold,
    that the device runs four-plus decades outside the data."""
    spec = DeviceSpec(needle_radius=1.05e-4, needle_length=12.7e-3, barrel_radius=2.3e-3,
                      flow_rate=1e-6 / 60, label="G22.5")
    d = build(recs["dmem_sints2024_10pct_20C_b"], spec)
    assert "non_commercial" not in d.flags
    assert "extrapolated" in d.flags
    assert d.envelope["decades_past"] > 4.0, d.envelope
    page = render(d)
    assert "NOT LICENSABLE" not in page
    assert "decades outside the measured range" in page


def test_g22_6_the_two_temperatures_do_NOT_collapse_under_one_shift_factor(recs):
    """**Pre-registered as a predicted FAILURE of the collapse, which is the opposite of
    how a 'does it work' test usually reads.**

    A one-parameter time-temperature shift maps a power law `k -> a^n k` with `n` FIXED
    (item 21's derivation). The paper's own indices are 0.32 at 20 C and 0.5 at 37 C for
    the same 10 % FBS medium, so the two curves differ in a parameter the shift cannot
    touch and the collapse must fail. A SMALL residual here would mean the shift
    implementation or the fitted indices are wrong.
    """
    cold = recs["dmem_sints2024_10pct_20C_b"]
    hot = recs["dmem_sints2024_10pct_37C_a"]
    lo = max(cold.envelope.gammadot_min, hot.envelope.gammadot_min)
    hi = min(cold.envelope.gammadot_max, hot.envelope.gammadot_max)
    gd = np.geomspace(lo, hi, 60)
    out = collapse_shift(gd, cold.eta_at(gd), hot.eta_at(gd))
    assert out["rms_log_residual"] > 0.05, out
    # and the reason, checked rather than asserted: the indices really do differ
    n_cold = float(cold.params_si[POWER_LAW.param_names.index("n")])
    n_hot = float(hot.params_si[POWER_LAW.param_names.index("n")])
    assert abs(n_cold - n_hot) > 0.1, (n_cold, n_hot)


def test_g22_6_has_a_positive_control_on_the_same_machinery(recs):
    """Without this, G22.6 is passed by a `collapse_shift` that returns a large residual
    for everything. A curve against a genuine shift of ITSELF must collapse."""
    from biorheo.constitutive.shift import shift_parameters
    rec = recs["dmem_sints2024_10pct_20C_b"]
    gd = np.geomspace(rec.envelope.gammadot_min, rec.envelope.gammadot_max, 60)
    eta = rec.eta_at(gd)
    shifted = POWER_LAW(gd, *shift_parameters(POWER_LAW, rec.params_si, 2.0))
    out = collapse_shift(gd, eta, shifted)
    assert out["rms_log_residual"] < 1e-6, out
    assert out["a"] == pytest.approx(2.0, rel=1e-3)


# --- the prereg's own discipline ----------------------------------------------------------

def test_the_prereg_enumerated_the_consumers_before_the_change():
    """Item 11's documented failure was appending to a shared list and silently re-deciding
    two closed gates. §4 of this item's pre-registration lists every consumer by grep,
    BEFORE the records were written."""
    t = PREREG.read_text(encoding="utf-8")
    assert "## 4. The consumers that must move, enumerated BEFORE the change" in t
    for f in ("test_offer_claims.py", "test_design_space.py",
              "test_openfoam_export_nonyield.py", "test_openfoam_export_null.py",
              "test_report.py", "test_profiling_wrapper.py"):
        assert f in t, f
    assert "no existing gate may change its VALUE" in t


def test_the_closed_gate_registry_is_frozen_at_115():
    """The mechanism that kept the promise above."""
    from conftest import _closed_gate_records
    assert len(_closed_gate_records()) == 115
