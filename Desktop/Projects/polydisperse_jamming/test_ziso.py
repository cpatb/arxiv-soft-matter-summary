"""test_ziso.py — anchors for the z_iso estimator repaired on 2026-09-12.

The defects this pins are in PREREG_ZISO_2026-09-12.md §2:

  D1  the loader read a SHARED GLOB that two campaigns write into, with no cap filter
      and no dedup by cell key;
  D2  `N_LOW_RUNGS = 5` was documented as rungs and implemented as RECORDS, so with
      three seeds per rung the extrapolation ran on TWO distinct pressures.

Anchors outside the code under test: a synthetic ladder with an intercept that is known
exactly, and a curvature case whose direction is known analytically.
"""

from __future__ import annotations

import numpy as np
import pytest

import measure_z_iso as M

Z_TRUE = 4.3
SLOPE = 2.0
RUNG_TARGETS = (3.0e-5, 1.0e-4, 3.0e-4, 1.0e-3)   # the production ladder
SEEDS_PER_RUNG = 3


def _ladder(targets=RUNG_TARGETS, seeds=SEEDS_PER_RUNG, curvature=0.0, jitter=0.0):
    """z = Z_TRUE + SLOPE*sqrt(p) + curvature*p, exact at every point."""
    rng = np.random.default_rng(20260912)
    recs = []
    for t in targets:
        for _ in range(seeds):
            p = float(t) * (1.0 + jitter * rng.standard_normal())
            recs.append({"p_hat": p, "p_hat_target": float(t),
                         "z_nr": Z_TRUE + SLOPE * np.sqrt(p) + curvature * p})
    return recs


# =====================================================================================
# the anchor: an intercept that is known exactly
# =====================================================================================
@pytest.mark.parametrize("rule", ["r2", "r3", "broken"])
def test_recovers_a_known_intercept(rule):
    """On a ladder that is exactly z_iso + a sqrt(p), every rule must return z_iso."""
    f = M.z_iso_fit(_ladder(), rule=rule)
    assert f["z_iso"] == pytest.approx(Z_TRUE, abs=1e-10)
    assert f["n_pts"] >= 2


def test_a_wrong_intercept_is_rejected():
    """The discriminating negative case: the anchor must not pass for any intercept."""
    f = M.z_iso_fit(_ladder(), rule="r3")
    assert f["z_iso"] != pytest.approx(Z_TRUE + 0.1, abs=1e-3)


# =====================================================================================
# D2 — the rung/record confusion, pinned so it cannot come back silently
# =====================================================================================
def test_shipped_rule_fits_two_distinct_pressures_not_five_rungs():
    """D2. This is the defect, asserted rather than described.

    `N_LOW_RECORDS_BROKEN = 5` with three seeds per rung is 2 rungs, on a ladder that
    offers 4. If this test starts failing because `n_rungs` went up, the shipped rule
    was quietly changed and the historical artifact can no longer be reproduced.
    """
    # a small jitter, because the real records never tie exactly: with exact ties the
    # `x <= sort(x)[4]` bound admits the whole tied rung and keeps 6, which is how the
    # rule behaves on synthetic data and not how it behaved on the ladders.
    f = M.z_iso_fit(_ladder(jitter=0.02), rule="broken")
    assert f["n_pts"] == M.N_LOW_RECORDS_BROKEN
    assert f["n_rungs"] == 2, "the shipped rule is a two-pressure fit"
    tied = M.z_iso_fit(_ladder(), rule="broken")
    assert tied["n_rungs"] == 2 and tied["n_pts"] == 2 * SEEDS_PER_RUNG


def test_rung_correct_rules_use_the_rungs_they_name():
    for rule, want in (("r2", 2), ("r3", 3)):
        f = M.z_iso_fit(_ladder(), rule=rule)
        assert f["n_rungs"] == want
        assert f["n_pts"] == want * SEEDS_PER_RUNG


def test_rungs_are_grouped_by_TARGET_not_by_realised_pressure():
    """The realised p_hat scatters within a rung; sorting on it mixes neighbours.

    With 40 % jitter the realised pressures of adjacent rungs overlap, so a rule that
    grouped on p_hat would not return a clean 2 x 3 selection.
    """
    f = M.z_iso_fit(_ladder(jitter=0.4), rule="r2")
    assert f["n_rungs"] == 2
    assert f["n_pts"] == 2 * SEEDS_PER_RUNG


# =====================================================================================
# the curvature direction, which is what §5's failed directional prediction was about
# =====================================================================================
def test_a_narrower_window_wins_when_the_curvature_is_REAL():
    """Analytic: adding a +curvature*p term biases the intercept upward, more so the
    wider the pressure range. So on data with genuine curvature r2 must beat r3.

    This is why the mu_p = 0 anchor coming out the other way (r3 closer to 6.0) is a
    statement about NOISE and not about curvature: the code responds to curvature in
    the direction predicted, and the data did not.
    """
    recs = _ladder(curvature=500.0)
    e2 = abs(M.z_iso_fit(recs, rule="r2")["z_iso"] - Z_TRUE)
    e3 = abs(M.z_iso_fit(recs, rule="r3")["z_iso"] - Z_TRUE)
    assert e2 < e3
    assert e3 > 1e-6, "with no curvature this test would be vacuous"


# =====================================================================================
# refusals — a fit that cannot be supported must return NaN, not a number
# =====================================================================================
def test_refuses_a_ladder_with_too_few_distinct_rungs():
    two = _ladder(targets=RUNG_TARGETS[:2])
    for rule in ("r2", "r3"):
        assert not np.isfinite(M.z_iso_fit(two, rule=rule)["z_iso"])


def test_rejects_an_unknown_rule():
    with pytest.raises(ValueError):
        M.z_iso_fit(_ladder(), rule="whatever")


def test_rejects_an_unknown_cap():
    with pytest.raises(ValueError):
        M.load_phase2(cap="repared")


# =====================================================================================
# D1 — cap selection and dedup, on the real files
# =====================================================================================
def test_cap_selection_partitions_the_records():
    d, di = M.load_phase2(cap="default")
    r, ri = M.load_phase2(cap="repaired")
    a, ai = M.load_phase2(cap="all")
    assert d and r, "both campaigns must be on disk for this test to mean anything"
    assert all(M._cap_of(x) == M.DEFAULT_MAX_ITER for x in d)
    assert all(M._cap_of(x) != M.DEFAULT_MAX_ITER for x in r)
    # the two caps partition "all": no record is in both, none is lost
    assert len(a) == len(d) + len(r)
    assert ai["n"] == di["n"] + ri["n"]
    # what one pass skipped on the cap is exactly what the other pass read
    assert di["wrong_cap"] == ri["n"] + ri["dupes"]
    assert ri["wrong_cap"] == di["n"] + di["dupes"]
    assert ai["wrong_cap"] == 0


def test_dedup_leaves_no_duplicate_keys_and_reports_the_conflicts():
    recs, info = M.load_phase2(cap="default")
    keys = [r["key"] for r in recs if r.get("key")]
    assert len(keys) == len(set(keys))
    assert info["dupes"] > 0, "the shard duplicates are a fact of this out/ directory"
    # the CV = 0 non-determinism (PREREG_DST_WINDOW_2026-09-09.md §3(c)) must stay
    # DETECTED. A zero here means the detection was removed, not that it was fixed.
    assert len(info["conflicts"]) > 0
    assert all("cv0.000" in k for k, _a, _b in info["conflicts"]), \
        "every known duplicate-state conflict is at CV = 0"


def test_no_dedup_reproduces_the_shipped_over_counting():
    """Kept so the historical artifact can be reproduced ON PURPOSE, not by accident."""
    _, with_ = M.load_phase2(cap="default", dedup=True)
    _, without = M.load_phase2(cap="default", dedup=False)
    assert without["n"] == with_["n"] + with_["dupes"]
    assert without["dupes"] == 0
