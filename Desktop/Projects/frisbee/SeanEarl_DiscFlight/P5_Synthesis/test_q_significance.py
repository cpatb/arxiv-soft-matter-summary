"""
test_q_significance.py -- checks for the Q-ordering significance analysis.

Run from P5_Synthesis:  python -m pytest test_q_significance.py -q
"""
import numpy as np
import pytest

import q_significance as qs

Q_NPZ = "Q.npz"
LANDING_NPZ = "../P4_Coupling/landing_stats.npz"


# --------------------------------------------------------------------------- #
# Holm correction: closed-form checks.
# --------------------------------------------------------------------------- #
def test_holm_matches_hand_computation():
    p = np.array([0.01, 0.04, 0.03])
    # sorted: 0.01 (x3), 0.03 (x2), 0.04 (x1), then enforce monotone non-decreasing
    expected = np.array([0.03, 0.06, 0.06])
    np.testing.assert_allclose(qs.holm(p), expected, rtol=1e-12)


def test_holm_is_monotone_and_bounded():
    rng = np.random.default_rng(0)
    p = rng.random(8)
    adj = qs.holm(p)
    assert np.all(adj <= 1.0)
    assert np.all(adj >= p - 1e-12)
    order = np.argsort(p)
    assert np.all(np.diff(adj[order]) >= -1e-12)


def test_holm_of_a_single_test_is_the_identity():
    np.testing.assert_allclose(qs.holm([0.023]), [0.023], rtol=1e-12)


# --------------------------------------------------------------------------- #
# Bootstrap.
# --------------------------------------------------------------------------- #
def test_bootstrap_ci_brackets_the_point_estimate():
    rng = np.random.default_rng(1)
    a = rng.normal(1.0, 0.1, 40)
    b = rng.normal(0.8, 0.1, 40)
    diff, lo, hi = qs.paired_bootstrap_ci(a, b, n_boot=4000)
    assert lo < diff < hi
    np.testing.assert_allclose(diff, np.mean(a - b), rtol=1e-12)


def test_bootstrap_is_reproducible():
    a, b = np.arange(12.0), np.arange(12.0) * 0.9
    assert qs.paired_bootstrap_ci(a, b, n_boot=2000) == \
        qs.paired_bootstrap_ci(a, b, n_boot=2000)


def test_bootstrap_ci_excludes_zero_for_a_clear_effect():
    rng = np.random.default_rng(2)
    a = rng.normal(1.0, 0.02, 30)
    b = a - 0.5
    _, lo, hi = qs.paired_bootstrap_ci(a, b, n_boot=4000)
    assert lo > 0.0 and hi > 0.0


# --------------------------------------------------------------------------- #
# Throw-subset recovery.  This is the provenance guard.
# --------------------------------------------------------------------------- #
def test_throw_subset_is_recovered_and_excludes_the_hammer():
    """The shipped Q.npz was built with --exclude-throws hammer but does not record it.
    The matcher must recover that rather than silently averaging over all three."""
    sub = qs.resolve_throw_subset(LANDING_NPZ, Q_NPZ)
    assert "hammer" not in sub["kept"]
    assert len(sub["kept"]) == 2
    assert sub["source"] in ("recorded", "inferred by matching")


def test_monte_carlo_error_shrinks_with_more_cells():
    """More contributing (throw, bearing) cells must reduce the error on their mean."""
    two = qs.monte_carlo_error(LANDING_NPZ, [0, 1])["se_success"]
    one = qs.monte_carlo_error(LANDING_NPZ, [0])["se_success"]
    assert np.nanmean(two) < np.nanmean(one)


def test_monte_carlo_error_is_positive_and_finite():
    mc = qs.monte_carlo_error(LANDING_NPZ, [0, 1])
    assert np.all(np.isfinite(mc["se_success"]))
    assert np.all(mc["se_success"] >= 0.0)


# --------------------------------------------------------------------------- #
# The headline conclusions, pinned so a future edit cannot quietly change them.
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def analysis():
    Qd = np.load(Q_NPZ, allow_pickle=True)
    Q, headings = Qd["Q"], Qd["headings"]
    sub = qs.resolve_throw_subset(LANDING_NPZ, Q_NPZ)
    mc = qs.monte_carlo_error(LANDING_NPZ, sub["idx"])
    return dict(Q=Q, headings=headings,
                decomp=qs.variance_decomposition(Q, mc["sr_mean"], mc["se_success"]),
                rows=qs.pairwise_tests(Q, headings))


def test_spatial_spread_dominates_monte_carlo_noise(analysis):
    """If this ever fails, the Q map has become a sampling artefact and the whole
    spatial argument in the paper needs re-examining."""
    for d in analysis["decomp"]:
        assert d["ratio"] > 5.0, d


def test_two_seventy_is_resolved_as_the_worst(analysis):
    rows = {(r["a"], r["b"]): r for r in analysis["rows"]}
    assert rows[(90, 270)]["resolved"]
    assert rows[(210, 270)]["resolved"]
    assert rows[(90, 270)]["mean_diff"] > 0
    assert rows[(210, 270)]["mean_diff"] > 0


def test_ninety_and_two_ten_are_not_separable(analysis):
    """The finding that corrects the synthesis: 'best = 210' is not supported."""
    rows = {(r["a"], r["b"]): r for r in analysis["rows"]}
    r = rows[(90, 210)]
    assert not r["resolved"]
    assert r["ci_lo"] < 0.0 < r["ci_hi"]


def test_verdict_reports_a_worst_but_no_best(analysis):
    v = qs.verdict(analysis["rows"], analysis["headings"])
    assert v["worst_supported"] == 270
    assert v["best_supported"] is None
