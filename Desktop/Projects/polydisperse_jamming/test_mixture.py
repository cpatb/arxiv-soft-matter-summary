"""test_mixture.py -- pin PREREG_MIXTURE_2026-09-13.md against its artifacts.

Fast tests only: everything here reads `out/phase12_mixr*_shard*.jsonl` and the model
modules. No DEM runs.
"""

from __future__ import annotations

import collections
import glob
import json
import os

import numpy as np
import pytest

import mixture_closure as mc
import mixture_model as mm
from phase12_mixture import X_ANCHOR
from psd import BidispersePSD, BinaryPSD, make_psd, parse_binary_family

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
SHARDS = sorted(glob.glob(os.path.join(OUT, "phase12_mixr*_shard*.jsonl")))


def cells(mu_p: float = 0.3):
    """The graded block. Defaults to mu_p = 0.3, which is what M0-M6 are closed on.

    🔴 This filter was MISSING when the frictionless block landed, and two closed-gate
    tests (M1 and M4) immediately went red with frictionless cells pooled into them --
    the same defect the grader had already been guarded against, in the test helper that
    checks the grader. The tests caught it, which is the point, but the lesson is that a
    loader is a loader wherever it lives.
    """
    g = collections.defaultdict(list)
    for f in SHARDS:
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                if r.get("error") or r.get("n_requested") != 1500:
                    continue
                if abs(float(r.get("mu_p", -1)) - mu_p) > 1e-9:
                    continue
                g[(float(r["mix_ratio"]), round(float(r["x_small_vol"]), 6))].append(r)
    return g


# --- the PSD family itself -----------------------------------------------------------

def test_binary_reproduces_bidisperse_bit_identically():
    """The anchor composition IS BidispersePSD(cv=0.5), reached the other way.

    This is what licenses the M0 gate to be stated as bit-identity rather than as
    agreement, so it is pinned here as well as run in the grader.
    """
    b = BinaryPSD(d10=1.0, ratio=3.0, x_small=X_ANCHOR)
    c = BidispersePSD(d10=1.0, cv=0.5)
    r1 = b.radii_for_yade(1500, np.random.default_rng(1))
    r2 = c.radii_for_yade(1500, np.random.default_rng(1))
    assert np.array_equal(r1, r2)
    assert b.d_small == 0.5 and b.d_large == 1.5
    assert b.cv == pytest.approx(0.5, abs=1e-15)


def test_family_string_round_trips_exactly():
    """A 6-decimal family string moves d_small off 0.5 by 1e-6 and would have destroyed
    the bit-identity gate. The runner uses repr; this pins that it must."""
    fam = f"binary:r{3.0!r}:x{X_ANCHOR!r}"
    assert parse_binary_family(fam) == (3.0, X_ANCHOR)
    assert np.array_equal(
        make_psd(fam, 1.0, 0.0).radii_for_yade(1500, np.random.default_rng(1)),
        BidispersePSD(d10=1.0, cv=0.5).radii_for_yade(1500, np.random.default_rng(1)))
    rounded = make_psd("binary:r3.0:x0.035714", 1.0, 0.0)
    assert rounded.d_small != 0.5          # the discriminating negative case


def test_binary_family_rejects_a_nonzero_cv_and_malformed_strings():
    with pytest.raises(ValueError, match="derives its CV"):
        make_psd("binary:r3.0:x0.25", 1.0, 0.5)
    for bad in ("binary:3.0:x0.25", "binary:r3.0", "binary:rX:x0.25", "binary:r3.0:y0.2"):
        with pytest.raises(ValueError):
            make_psd(bad, 1.0, 0.0)
    assert parse_binary_family("lognormal") is None


def test_matched_pair_is_matched_in_every_number_descriptor():
    """M1 rests on the two packs agreeing in CV and ratio. If that ever stops being
    true the gate is measuring something else."""
    a = BinaryPSD(d10=1.0, ratio=3.0, x_small=X_ANCHOR)
    b = BinaryPSD(d10=1.0, ratio=3.0, x_small=0.25)
    assert a.cv == pytest.approx(b.cv, abs=1e-12)
    assert a.cv == pytest.approx(0.5, abs=1e-12)
    assert a.d_large / a.d_small == pytest.approx(b.d_large / b.d_small, rel=1e-12)
    assert a.d32 != pytest.approx(b.d32)   # d32 is NOT matched, and the report says so


# --- the model's limits --------------------------------------------------------------

def test_model_limits_hold():
    assert mm.check_limits() == []


def test_ab_form_reduces_to_the_closed_form():
    assert mc.check_ab_reduces() == []


def test_ratio_one_and_single_class_are_beta():
    for beta in (0.55, 0.5947, 0.64):
        assert mm.packing_density([1.0], [1.0], beta, 1.02, 1.5) == pytest.approx(beta)
        for x in (0.1, 0.5, 0.9):
            assert mm.binary_density(1.0, x, beta, 1.02, 1.5) == pytest.approx(beta)


# --- the landed campaign -------------------------------------------------------------

@pytest.mark.skipif(not SHARDS, reason="campaign not on disk")
def test_the_closed_gates_see_only_the_friction_they_were_closed_on():
    """M0-M6 are closed at mu_p = 0.3 and the shard glob matches any tag, so a block at
    another friction must not reach them. This test is what stops that happening
    silently; it is written to FAIL the day the filter is removed."""
    import grade_mixture as gm
    assert gm.GRADED_MU_P == 0.3
    seen = {float(c["mu_p"]) for v in gm.load_cells().values() for c in v
            if "error" not in c}
    assert seen == {0.3}, seen
    for mu in (0.3, 0.0):
        for r, _x, _phi in mc.binary_points(mu_p=mu):
            assert r > 0          # the call must not raise at either friction
    # and the two loaders must actually disagree when a second friction is on disk
    n03 = len(mc.binary_points(mu_p=0.3))
    n00 = len(mc.binary_points(mu_p=0.0))
    assert n03 > 0
    assert n00 != n03 or n00 == 0


@pytest.mark.skipif(not SHARDS, reason="campaign not on disk")
def test_m0_endpoints_and_anchor_are_bit_identical():
    g = cells()
    for seed in (1, 2, 3):
        x0 = [r for r in g[(3.0, 0.0)] if r["seed"] == seed][0]["phi"]
        x1 = [r for r in g[(3.0, 1.0)] if r["seed"] == seed][0]["phi"]
        r2 = [r for r in g[(2.0, 0.0)] if r["seed"] == seed][0]["phi"]
        assert x0 == x1 == r2


@pytest.mark.skipif(not SHARDS, reason="campaign not on disk")
def test_m1_matched_pair_moves_phi():
    g = cells()
    a = np.mean([r["phi"] for r in g[(3.0, round(X_ANCHOR, 6))]])
    b = np.mean([r["phi"] for r in g[(3.0, 0.25)]])
    assert a == pytest.approx(0.6042, abs=5e-4)
    assert b == pytest.approx(0.6677, abs=5e-4)
    assert b - a == pytest.approx(0.0635, abs=1e-3)
    assert (b - a) / mm.PHI_SCALE > 5.0            # the registered M1 bar
    # and the zero-parameter prediction, which had nothing fitted to this campaign
    beta = 0.5947
    pred = (mm.binary_density(3.0, 0.25, beta, 1.02, 1.50)
            - mm.binary_density(3.0, X_ANCHOR, beta, 1.02, 1.50))
    assert abs(pred / (b - a) - 1.0) < 0.10
    # DISCRIMINATING: a deviation of exactly zero would mean the two are one calculation
    assert abs(pred - (b - a)) > 1e-6


@pytest.mark.skipif(not SHARDS, reason="campaign not on disk")
def test_m4_size_ratio_matters():
    """The blind null predicts the same curve at both ratios, i.e. a difference of 0."""
    g = cells()
    peak = {}
    for r in (2.0, 3.0):
        vals = [(x, float(np.mean([c["phi"] for c in v])))
                for (rr, x), v in g.items() if rr == r and 0.0 < x < 1.0
                and not (rr == 3.0 and x == 0.6)]
        peak[r] = max(v for _, v in vals) - 0.5940
    assert peak[3.0] - peak[2.0] > 3.0 * mm.PHI_SCALE
    assert peak[3.0] == pytest.approx(0.0823, abs=1e-3)
    assert peak[2.0] == pytest.approx(0.0377, abs=1e-3)


@pytest.mark.skipif(not SHARDS, reason="campaign not on disk")
def test_phi_is_not_monotone_in_cv_within_one_ratio():
    """The negative, from inside one two-atom family: CV cannot be the control."""
    g = cells()
    for r in (2.0, 3.0):
        pts = sorted((float(v[0]["cv_realised"]), float(np.mean([c["phi"] for c in v])))
                     for (rr, _x), v in g.items() if rr == r)
        assert not all(pts[i][1] <= pts[i + 1][1] for i in range(len(pts) - 1))


@pytest.mark.skipif(not SHARDS, reason="campaign not on disk")
def test_closure_transfers_better_than_the_textbook_exponents():
    """Calibrated on binaries only, the functional must beat the textbook exponents on
    the 15 held-out continuous-family points -- and beat the CV collapse by a lot."""
    pts = [p for p in mc.binary_points() if p[0] in (2.0, 3.0)]
    meas = mm.measured()
    beta = meas[("lognormal", 0.0)][0]
    pa, pb, _ = mc.fit_on_binaries(pts, beta)
    chi = {"txt": 0.0, "fit": 0.0}
    for fam in mm.FAMILIES:
        for cv in (0.1, 0.2, 0.3, 0.4, 0.5):
            if (fam, cv) not in meas:
                continue
            d, y = mm.psd_classes(fam, cv)
            phi = meas[(fam, cv)][0]
            chi["txt"] += ((mm.packing_density(d, y, beta, 1.02, 1.50) - phi)
                           / mm.PHI_SCALE) ** 2
            chi["fit"] += ((mm.packing_density(d, y, beta, pa, pb) - phi)
                           / mm.PHI_SCALE) ** 2
    assert chi["fit"] < 0.25 * chi["txt"]
    assert chi["fit"] == pytest.approx(63.0, rel=0.15)


@pytest.mark.skipif(not SHARDS, reason="campaign not on disk")
def test_f1_the_composition_effect_is_friction_independent():
    """The matched pair at mu_p = 0 and 0.3, which is the blind half of the campaign.

    beta moves 22 scale between the two blocks and the matched-pair dphi moves 0.4, so
    this asserts the ABSOLUTE agreement, not a ratio. It is written to fail if either
    block is re-run into a different answer.
    """
    import grade_mixture_mu0 as gm0
    g0 = gm0.load(0.0)
    assert g0, "the frictionless block is not on disk"
    a = np.mean([c["phi"] for c in g0[(3.0, round(X_ANCHOR, 6))]])
    b = np.mean([c["phi"] for c in g0[(3.0, 0.25)]])
    d0 = b - a
    assert d0 == pytest.approx(0.0643, abs=1e-3)
    g3 = cells()
    d3 = (np.mean([c["phi"] for c in g3[(3.0, 0.25)]])
          - np.mean([c["phi"] for c in g3[(3.0, round(X_ANCHOR, 6))]]))
    assert abs(d0 - d3) < 1.0 * mm.PHI_SCALE          # friction-independent in ABSOLUTE terms
    # and beta itself moves far more than that, which is what makes the above a claim
    meas = mm.measured(mu_p=0.0), mm.measured(mu_p=0.3)
    dbeta = meas[0][("lognormal", 0.0)][0] - meas[1][("lognormal", 0.0)][0]
    assert abs(dbeta) > 20.0 * mm.PHI_SCALE


@pytest.mark.skipif(not SHARDS, reason="campaign not on disk")
def test_frictionless_endpoints_recover_beta_mu0():
    """The control that the block is at the friction it claims."""
    import grade_mixture_mu0 as gm0
    g0 = gm0.load(0.0)
    beta0 = mm.measured(mu_p=0.0)[("lognormal", 0.0)][0]
    for r in (2.0, 3.0):
        for x in (0.0, 1.0):
            v = float(np.mean([c["phi"] for c in g0[(r, x)]]))
            assert abs(v - beta0) < 0.5 * mm.PHI_SCALE


@pytest.mark.skipif(not SHARDS, reason="campaign not on disk")
def test_ratio_1p25_is_not_identifiable_and_is_refused():
    """The refusal in mixture_closure.py rests on this: at ratio 1.25 the composition
    curve spans under 2 scale, so neither interaction parameter is determined. If that
    ever stops being true the refusal must be revisited rather than left standing."""
    pts = [p for p in mc.binary_points() if p[0] == 1.25]
    assert pts, "the ratio 1.25 block is not on disk"
    span = max(p[2] for p in pts) - min(p[2] for p in pts)
    assert span < 2.0 * mm.PHI_SCALE
    lo, hi = mc.profile_interval(pts, 0.5947, which=1, ngrid=61)
    assert hi - lo > 0.9
