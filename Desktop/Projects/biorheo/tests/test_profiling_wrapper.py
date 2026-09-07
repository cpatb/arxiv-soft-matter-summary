"""Tests for the corrected profile-likelihood interval (NEXT_STEPS item 13).

Three things are pinned here and they are different in kind:

1. **The coupling to `dst_fit`.** The wrapper imports a PRIVATE name so that both
   routines re-minimise the nuisance parameters identically. A private name can move
   without warning, and a silent fallback to a local copy of the objective is how two
   implementations drift apart, so its presence and contract are asserted.
2. **The anchor.** `newtonian` has a closed-form profile interval that passes through
   neither the outward walk nor the root finder. This is the only test here that can
   fail for a reason the module's own self-check would not already catch.
3. **The defect, as a negative control.** `dst_fit`'s endpoint is asserted to sit ON the
   bisection grid and to be far past `f_crit`. If that ever stops being true the wrapper
   has become unnecessary, which is worth failing a test over rather than discovering by
   reading.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from dst_fit.fitting import fit_model
from dst_fit.profiling import profile_likelihood
from scipy.stats import f as f_dist

from biorheo.constitutive import get_model
from biorheo.identify import verdict_for
from biorheo.identify.profiling import (
    ProfileInterval,
    newtonian_closed_form,
    profile_interval,
    profile_statistic,
)

REPO = Path(__file__).resolve().parent.parent
GATES = REPO / "data" / "profile_interval_gates.json"
DIAG = REPO / "data" / "profile_interval_diagnosis.json"     # today: a regression check
PREFIX = REPO / "data" / "profile_interval_prefix.json"       # how bad it was

# dst_fit.profiling: step floor span/5000, then _REFINEMENTS = 10 bisections
GRID_DIVISOR = 5000.0 * 1024


def _curve(model_name: str, params: tuple[float, ...], n: int = 24, noise: float = 0.05,
           seed: int = 12345) -> tuple[np.ndarray, np.ndarray]:
    """A synthetic rate curve with multiplicative log-normal noise.

    Noise is required: `profile_interval` refuses a curve whose rmse is at the level of
    float64 rounding, because a profile against no residual scale is a statement about
    arithmetic precision rather than about data.
    """
    model = get_model(model_name)
    gd = np.logspace(-1, 2.5, n)
    rng = np.random.default_rng(seed)
    eta = np.asarray(model(gd, *params), dtype=float)
    return gd, eta * np.exp(rng.normal(0.0, noise, size=gd.shape))


def _fit(model, gd, eta, seed: int = 0):
    return fit_model(model, gd, eta, n_starts=8, seed=seed, start=model.guess(gd, eta))


# --------------------------------------------------------------------------------------
# 1. the coupling to dst_fit
# --------------------------------------------------------------------------------------

def test_dst_fit_private_nuisance_minimiser_is_present_and_has_its_contract():
    """The wrapper depends on `dst_fit.profiling._nuisance_minimiser(model, x, log_eta, i)`
    returning a callable `(value, warm) -> (ssr, params)`. Asserted rather than assumed:
    if dst_fit renames it the import fails loudly at import time, and if it changes SHAPE
    this test says so."""
    from dst_fit.profiling import _nuisance_minimiser

    model = get_model("casson")
    gd, eta = _curve("casson", (0.01, 0.004))
    run = _nuisance_minimiser(model, gd, np.log(eta), 0)
    ssr, params = run(0.01, np.array([0.01, 0.004]))
    assert isinstance(ssr, float)
    assert np.asarray(params).shape == (model.n_params,)
    assert ssr >= 0.0


# --------------------------------------------------------------------------------------
# 2. the anchor: a closed form that uses none of the module's machinery
# --------------------------------------------------------------------------------------

def test_newtonian_closed_form_matches_the_walk_and_root_finder():
    """For `eta = eta_0` with log-space residuals the profile is exactly
    `SSR(u) - SSR_min = n (u - ybar)^2`, so `F = f_crit` at `u = ybar +/- s sqrt(f_crit/n)`.
    Nothing in that derivation goes through the outward walk or Brent."""
    model = get_model("newtonian")
    gd, eta = _curve("newtonian", (0.004,), n=30, noise=0.2)
    fit = _fit(model, gd, eta)
    got = profile_interval(model, gd, eta, model.param_names[0], fit=fit)
    best, lo, hi = newtonian_closed_form(eta)

    assert got.best == pytest.approx(best, rel=1e-8)
    assert got.lower == pytest.approx(lo, rel=1e-6)
    assert got.upper == pytest.approx(hi, rel=1e-6)
    assert not got.lower_at_bound and not got.upper_at_bound


def test_newtonian_closed_form_is_not_a_restatement_of_the_wald_interval():
    """A discriminating negative case. If the closed form happened to equal the Wald
    interval it would anchor nothing, because the Wald error is what the walk starts
    from. The two must differ: the profile of a log-space location fit is symmetric in
    `log eta` and therefore ASYMMETRIC in `eta`."""
    _, eta = _curve("newtonian", (0.004,), n=30, noise=0.2)
    best, lo, hi = newtonian_closed_form(eta)
    assert (hi - best) / (best - lo) > 1.05


def test_closed_form_refuses_too_few_points():
    with pytest.raises(ValueError, match="at least 2 points"):
        newtonian_closed_form(np.array([1.0]))


# --------------------------------------------------------------------------------------
# 3. the defect, as a negative control on dst_fit
# --------------------------------------------------------------------------------------

def test_dst_fit_endpoint_is_a_bisection_grid_point_and_is_far_past_f_crit():
    """The reason this module exists, asserted rather than described.

    `dst_fit` sizes its step floor as `span/5000` from the parameter's BOX bound and then
    bisects 10 times, so a crossing found on the first step lands on
    `best + m * span / 5.12e6` for an integer `m`. The endpoint is then not a confidence
    limit: `F` there is orders of magnitude past `f_crit`.
    """
    model = get_model("casson")
    gd, eta = _curve("casson", (0.01, 0.004))
    fit = _fit(model, gd, eta)
    old = profile_likelihood(model, gd, eta, "eta_c", fit=fit, seed=0)

    hi_b = float(model.bounds[1][model.param_names.index("eta_c")])
    span = abs(hi_b - old.best)
    m = (old.upper - old.best) * GRID_DIVISOR / span
    assert abs(m - round(m)) < 1e-6, f"grid index {m} is not an integer"

    f_there, f_crit = profile_statistic(model, gd, eta, "eta_c", old.upper, fit=fit)
    assert f_there > 10 * f_crit, f"F = {f_there} against f_crit = {f_crit}"


def test_the_wrapper_endpoint_is_the_crossing_where_dst_fit_s_is_not():
    """Same curve, same fit, same nuisance minimiser: only the walk differs."""
    model = get_model("casson")
    gd, eta = _curve("casson", (0.01, 0.004))
    fit = _fit(model, gd, eta)
    new = profile_interval(model, gd, eta, "eta_c", fit=fit)
    old = profile_likelihood(model, gd, eta, "eta_c", fit=fit, seed=0)

    assert new.crossing_error() < 1e-2
    assert new.upper < old.upper
    f_there, f_crit = profile_statistic(model, gd, eta, "eta_c", new.upper, fit=fit)
    assert f_there == pytest.approx(f_crit, rel=1e-2)


# --------------------------------------------------------------------------------------
# the object's own contract
# --------------------------------------------------------------------------------------

def test_profile_statistic_is_zero_at_the_best_fit():
    model = get_model("casson")
    gd, eta = _curve("casson", (0.01, 0.004))
    fit = _fit(model, gd, eta)
    f_at, _ = profile_statistic(model, gd, eta, "eta_c", float(fit.params[1]), fit=fit)
    assert abs(f_at) < 1e-6


def test_f_crit_is_the_F_distribution_critical_value_and_not_a_chi2_one():
    """A profile interval quoted against the wrong reference distribution is wrong by a
    fixed factor, which no self-check inside the module could see."""
    model = get_model("casson")
    gd, eta = _curve("casson", (0.01, 0.004))
    fit = _fit(model, gd, eta)
    p = profile_interval(model, gd, eta, "eta_c", fit=fit)
    assert p.f_crit == pytest.approx(f_dist.ppf(0.95, 1, gd.size - model.n_params))
    assert p.n_points == gd.size and p.n_params == model.n_params


def test_unknown_parameter_name_raises():
    model = get_model("casson")
    gd, eta = _curve("casson", (0.01, 0.004))
    fit = _fit(model, gd, eta)
    with pytest.raises(ValueError, match="has no parameter"):
        profile_interval(model, gd, eta, "not_a_parameter", fit=fit)
    with pytest.raises(ValueError, match="has no parameter"):
        profile_statistic(model, gd, eta, "not_a_parameter", 1.0, fit=fit)


def test_a_noise_free_curve_is_refused_rather_than_given_an_interval():
    """`rmse` at the level of float64 rounding leaves no residual scale, so an interval
    would describe arithmetic precision. dst_fit refuses this and so does the wrapper."""
    model = get_model("casson")
    gd = np.logspace(-1, 2.5, 24)
    eta = np.asarray(model(gd, 0.01, 0.004), dtype=float)
    fit = _fit(model, gd, eta)
    with pytest.raises(ValueError, match="float64 rounding"):
        profile_interval(model, gd, eta, "eta_c", fit=fit)


def test_verdict_for_accepts_a_ProfileInterval():
    """`verdict_for` is typed for either object; a shim between them would be a place for
    the two to disagree."""
    model = get_model("casson")
    gd, eta = _curve("casson", (0.01, 0.004))
    fit = _fit(model, gd, eta)
    p = profile_interval(model, gd, eta, "eta_c", fit=fit)
    verdict, partner = verdict_for(p, fit.corr[1], 1, model.param_names)
    assert verdict in ("measurable", "bounded_only", "degenerate", "not_identifiable")


def test_an_unverified_end_is_reported_open_and_carries_what_was_rejected():
    """Constructed, not sampled: an end rejected as a discontinuity must be OPEN (so no
    downstream reader quotes it) and must still carry the value and F that were rejected
    (so the rejection is auditable)."""
    p = ProfileInterval(
        param="x", best=1.0, lower=0.5, upper=1e6, lower_at_bound=False,
        upper_at_bound=True, wald_stderr=0.1, f_crit=4.0, n_points=16, n_params=2,
        bounds=(0.0, 1e6), f_at_lower=4.0, f_at_upper=0.2, step_rule="wald", n_evals=10,
        upper_unverified=(123.0, 0.24),
    )
    assert not p.verified
    assert p.upper_at_bound and p.upper == 1e6
    assert p.upper_unverified == (123.0, 0.24)
    # the rejected end is excluded from crossing_error and reported by rejected_error
    assert p.crossing_error() == pytest.approx(0.0)
    assert p.rejected_error() == pytest.approx(abs(0.24 - 4.0) / 4.0)


# --------------------------------------------------------------------------------------
# the committed artifacts
# --------------------------------------------------------------------------------------

def test_gate_artifact_records_all_gates_passing():
    assert GATES.exists(), "run scripts/profile_interval_gates.py"
    g = json.loads(GATES.read_text(encoding="utf-8"))["verdict"]
    for gid in ("G13.2", "G13.3", "G13.4", "G13.5", "G13.8"):
        assert g[gid]["pass"], f"{gid} is recorded as failing: {g[gid]}"


def test_gate_artifact_numbers_are_the_ones_quoted_in_prose():
    """Every number in a `.md` under `data/` is printed by a script; these are the ones
    the Status card and NEXT_STEPS quote, so they are pinned rather than re-derived."""
    g = json.loads(GATES.read_text(encoding="utf-8"))["verdict"]
    assert g["G13.2"]["n"] == 70
    assert g["G13.2"]["rate"] == 1.0
    assert g["G13.2"]["median_ratio"] == pytest.approx(710.559, rel=1e-3)
    assert g["G13.2"]["worst_ratio"] == pytest.approx(8130.12, rel=1e-3)
    assert g["G13.3"]["worst_crossing_error"] < 1e-3
    assert g["G13.4"]["worst_rel"] < 1e-8
    assert g["G13.5"]["max_old_over_new"] == pytest.approx(1227.39, rel=1e-3)


def test_no_profile_beats_the_fit_it_is_centred_on():
    """A free by-product of the corrected walk, and a check on the FITS that nothing else
    in this repo performs. `F < 0` somewhere on a profile means the nuisance
    re-minimisation found a lower SSR than the reported optimum, i.e. the interval is
    drawn around the wrong point. The old walk took one or two steps and could not have
    noticed; this one evaluates the profile densely enough to."""
    g = json.loads(GATES.read_text(encoding="utf-8"))["verdict"]["fit_quality_scan"]
    assert g["n"] == 325
    assert g["n_profiles_beating_their_own_fit"] == 0
    assert g["worst_negative_f"] > -1e-6, g["worst_negative_f"]


def test_the_g139_gates_are_recorded_as_FAILED_and_stay_that_way():
    """Both arms of G13.9 failed against a 20/20 prediction, and the bars were not
    touched. Pinned so a later run cannot quietly turn a recorded failure into a pass by
    moving a bar: if these ever flip, the change must be deliberate and explained."""
    p = REPO / "data" / "profile_interval_g139.json"
    assert p.exists(), "run scripts/profile_interval_g139.py"
    v = json.loads(p.read_text(encoding="utf-8"))["verdict"]
    assert v["G13.9a"]["pass"] is False and v["G13.9a"]["bar"] == 11
    assert v["G13.9b"]["pass"] is False and v["G13.9b"]["bar"] == 19
    # the reason both failed: the comparison binds on the LOWER end, which the bug never
    # touched, so exactly one donor differs between the old and corrected hulls
    assert v["binding_end_new"]["lower"] == 8 and v["binding_end_new"]["upper"] == 1
    assert v["binding_end_old"]["lower"] == 8 and v["binding_end_old"]["upper"] == 0
    # and the salvage: parameter uncertainty is the LARGER of the two effects
    assert v["median_new_hull_width_pct"] > v["median_abs_pct_diff"]


def test_the_prefix_artifact_counts_how_bad_it_was():
    """The HISTORICAL measurement, re-gradable forever from the rows committed at
    `d38bc77`. It cannot be reproduced from the working tree, because item 13 regenerated
    every sidecar in that same commit -- which is why the rows are preserved separately.

    `NEXT_STEPS` first said 117 of 296 and flagged itself conservative; reading each
    model's own declared bound gave 203; requiring the grid index to be one a 10-fold
    bisection can actually produce (item 16) gives 201. All three numbers are pinned in
    the artifact's own prose so the correction cannot quietly revert."""
    assert PREFIX.exists(), "run scripts/profile_interval_regrade.py"
    s = json.loads(PREFIX.read_text(encoding="utf-8"))["summary"]
    assert s["n_closed_upper"] == 296
    assert s["n_upper_on_grid"] == 201
    assert s["as_shipped_upper_before_regrade"] == 203
    assert s["n_rejected_out_of_range"] == 18
    assert s["n_rejected_at_the_bound"] == 5
    assert s["worst"]["delta_over_value"] > 1000
    # a bisection index the arithmetic can produce
    assert 1 <= s["worst"]["grid_index"] <= 1024


def test_the_live_scan_says_the_defect_is_gone():
    """The same scan against TODAY's registry is a regression check and must report zero.
    If this ever fires, some sidecar has been written by the old walk again."""
    assert DIAG.exists(), "run scripts/profile_interval_diagnosis.py"
    s = json.loads(DIAG.read_text(encoding="utf-8"))["summary"]
    assert s["n_closed_ends"] > 500, "the scan found almost nothing; did it read the records?"
    assert s["n_on_grid"] == 0, s["on_grid_by_model_param"]


def test_every_sidecar_carries_the_exposure_ratio():
    """G16.4. `exposure = (span/5000)/(wald/4)` says whether the defective step floor set
    the walk that produced an interval, and by how much it overshot -- so a reader can ask
    that of a shipped number without refitting the curve. Two ends, because the two spans
    differ by decades."""
    import yaml

    records = REPO / "src" / "biorheo" / "registry" / "records"
    sidecars = sorted(records.glob("*.identifiability.yaml"))
    assert len(sidecars) >= 115
    n_defined = 0
    for path in sidecars:
        block = yaml.safe_load(path.read_text(encoding="utf-8"))
        assert "exposure_note" in block, path.name
        for name, est in block["parameters"].items():
            assert "wald_stderr" in est, (path.name, name)
            assert "exposure_hi" in est and "exposure_lo" in est, (path.name, name)
            if est["exposure_hi"] is not None:
                n_defined += 1
                # the two ends are different questions and must not be the same number
                assert est["exposure_lo"] is not None
                assert est["exposure_hi"] != est["exposure_lo"], (path.name, name)
    assert n_defined > 300, n_defined


def test_the_exposure_ratio_reproduces_the_item_15_numbers():
    """A discriminating check: the field must reproduce numbers measured independently by
    `case_studies/injectability_surrogates/uq_check.py`, which computes them from a refit
    and never reads a sidecar. If they agreed only with themselves this would anchor
    nothing."""
    import yaml

    records = REPO / "src" / "biorheo" / "registry" / "records"
    block = yaml.safe_load(
        (records / "carbopol_etd2020_0p5wt_50uL.identifiability.yaml").read_text(
            encoding="utf-8"))
    got = {k: v["exposure_hi"] for k, v in block["parameters"].items()}
    assert got["k"] == pytest.approx(1.177e4, rel=2e-3)
    assert got["tau_y"] == pytest.approx(1755, rel=2e-3)
    # the flow index is the one the defect could barely touch, and it is why that study's
    # forces did not move: its box is O(1), so the floor was the right size anyway
    assert got["n"] == pytest.approx(1.092, rel=2e-3)


def test_every_registry_value_lies_inside_its_corrected_interval():
    """G13.7. Narrowing every interval can only make containment harder, and the Carbopol
    records were the pre-registered risk because their shipped values come from a
    different fitter. `load_all` enforces it; this asserts the population is the whole
    registry and not an empty set."""
    from biorheo.registry import load_all

    recs = load_all()
    assert len(recs) >= 115
