"""Tests for the holes the 2026-08-29 audit found in the coverage.

Every one of these guards a failure mode that has ALREADY happened once in this
project and was caught by a human reading a file listing, not by the suite:

  F5  `landing_stats.py`'s schema guard was never exercised, and no test ever read
      `gust_model` or a `prov_*` key back out of a real artifact. A stale
      `landing_stats.npz` -- eleven hours older than the code that writes it, carrying
      no record of which gust model built it -- reached the manuscript's Q.

  F6  `test_gust_autocorrelation_matches_target_scale` checks the ANALYTIC ensemble
      autocorrelation computed from the deterministic mode weights. It never
      correlates a sampled xi(t), and it runs at a correlation time SHORTER than a
      flight, where the frozen-gust behaviour the temporal model was retired for
      cannot appear at all.

  F7  `validation_hummel.py` prints a 9% trim-speed error against Hummel (2003) and
      nothing asserts on it, so a regression there would be silent.

Plus two tests for the re-graded gate (`regrade_gate_g1.py`).
"""

from __future__ import annotations

import hashlib
import io
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import landing_stats
import synthetic_turbulence as st
from synthetic_turbulence import CorrelatedGusts
from test_coupling import ConstantTurbField, StubField

HERE = Path(__file__).resolve().parent
REAL_ARTIFACT = HERE / "landing_stats.npz"
FIELD = HERE.parent / "P3_CFD" / "flow_fields.npz"


# --------------------------------------------------------------------------------
# F5: the schema guard, and provenance that is actually read back
# --------------------------------------------------------------------------------

def _make_stale_npz(tmp_path):
    """An artifact in the PRE-FIX schema: the keys the July file actually had."""
    p = tmp_path / "stale_landing_stats.npz"
    n = 2
    np.savez_compressed(
        p,
        landing_mean=np.zeros((n, n, 1, 1, 1, 2)),
        landing_cov=np.zeros((n, n, 1, 1, 1, 2, 2)),
        ellipse_axes=np.zeros((n, n, 1, 1, 1, 2)),
        alpha_landing=np.zeros((n, n, 1, 1, 1)),
        success_rate=np.zeros((n, n, 1, 1, 1)),
        blocked_fraction=np.zeros((n, n, 1, 1, 1)),
        launch_x=np.arange(n, dtype=float), launch_y=np.arange(n, dtype=float),
        launch_z=np.ones(n), headings=np.array([90]),
        throws=np.array(["right_backhand"]), bearings=np.array([0.0]),
        n_throws=np.array(1), n_valid=np.zeros((n, n, 1, 1, 1)),
        seed=np.array(1), coeff_set=np.array("potts_tabulated"),
        field_source=np.array("stub"), target_radius_m=np.array(3.0),
        target_range_m=np.array(25.0),
    )
    return p


def _bias_rows(out):
    return [ln for ln in out.splitlines() if ln.strip().startswith("heading ")]


def test_schema_guard_refuses_a_pre_fix_artifact(tmp_path, capsys):
    """The guard must FIRE, and say why, on an artifact without landing_bias.

    Discriminating negative case: this is the artifact shape that silently produced a
    world-frame "mean drift" number with no interpretation.
    """
    p = _make_stale_npz(tmp_path)
    landing_stats.summarize(str(p))
    out = capsys.readouterr().out
    assert "landing_bias not present" in out, \
        f"the schema guard did not fire on a pre-fix artifact; got:\n{out}"
    assert "Re-run monte_carlo.py" in out
    # It must not go on to print bias ROWS. ("downrange" alone is not a discriminator:
    # it appears in the header printed before the guard is reached.)
    assert not _bias_rows(out), \
        "summarize() printed bias rows for an artifact with no landing_bias"


def test_schema_guard_passes_the_real_artifact(capsys):
    """The complement: the guard must NOT fire on the current artifact."""
    if not REAL_ARTIFACT.exists():
        pytest.skip("landing_stats.npz not present")
    landing_stats.summarize(str(REAL_ARTIFACT))
    out = capsys.readouterr().out
    assert "landing_bias not present" not in out
    assert _bias_rows(out), "the real artifact produced no bias rows"


@pytest.mark.skipif(not REAL_ARTIFACT.exists(), reason="landing_stats.npz not present")
def test_artifact_records_which_gust_model_built_it():
    """The question "which gust model is Q built on?" must be answerable from the file.

    This is the exact question that could not be answered on 2026-08-29 and had to be
    settled by comparing file timestamps.
    """
    d = np.load(REAL_ARTIFACT, allow_pickle=False)
    assert "gust_model" in d.files, \
        "the artifact does not record its gust model; a temporal/sweeping switch " \
        "would be undetectable from the file"
    model = str(d["gust_model"])
    assert model in ("temporal", "sweeping"), f"unknown gust model {model!r}"
    assert model == "sweeping", (
        f"landing_stats.npz was built with the {model!r} gust model. The temporal "
        "model overstates turbulence-induced dispersion by ~1.7x (GATE G1); the "
        "manuscript's Q must rest on the sweeping model.")
    assert "gusts_enabled" in d.files and bool(d["gusts_enabled"])


@pytest.mark.skipif(not REAL_ARTIFACT.exists() or not FIELD.exists(),
                    reason="artifacts not present")
def test_artifact_provenance_matches_the_field_it_claims():
    """prov_field_sha256 must be the hash of the field actually on disk.

    Catches a stale artifact carried forward past a field regeneration -- the failure
    the whole prov_ block was added for, but which nothing read back until now.
    """
    d = np.load(REAL_ARTIFACT, allow_pickle=False)
    for key in ("prov_field_sha256", "prov_git_sha", "prov_numpy", "prov_n_modes"):
        assert key in d.files, f"provenance key {key} missing from the artifact"

    h = hashlib.sha256()
    with open(FIELD, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    assert str(d["prov_field_sha256"]) == h.hexdigest(), (
        "landing_stats.npz was computed on a DIFFERENT flow field from the one now "
        f"at {FIELD}. Downstream Q and every manuscript macro are stale.")


@pytest.mark.skipif(not REAL_ARTIFACT.exists(), reason="landing_stats.npz not present")
def test_recorded_gust_dt_is_not_used_by_the_production_path():
    """prov_gust_dt records a constant the production gust class never consumes.

    Not a bug in the results, but a reader auditing reproducibility from the
    provenance block would reasonably conclude this knob affected them. Pinned so that
    if GUST_DT ever does become live, this fails and the note gets corrected.
    """
    src = io.open(st.__file__, encoding="utf-8").read()
    body = src.split("class CorrelatedGusts")[1].split("\nclass ")[0]
    assert "GUST_DT" not in body, (
        "GUST_DT is now referenced inside CorrelatedGusts. It is live; update "
        "monte_carlo._provenance's documentation and delete this test.")


# --------------------------------------------------------------------------------
# F6: the realised gust field, at a correlation scale the site actually has
# --------------------------------------------------------------------------------

class _SiteTurbField(StubField):
    """Uniform k and eps giving a correlation scale the site actually has.

    GATE G1 measured scale_mean = 8.43 (seconds, temporal arm) over the real field,
    spread 4.71-16.25. `ConstantTurbField` gives T_L = k/eps = 1 s, which is SHORTER
    than a flight, so the frozen-gust behaviour the temporal model was retired for
    cannot appear there at all. The distinction only exists at the site's own T_L.
    """

    def __init__(self, r, k_val=3.0, T_L=8.43):
        super().__init__(r)
        self._k = float(k_val)
        self._eps = float(k_val) / float(T_L)

    def k(self, r, t):
        return self._k

    def eps(self, r, t):
        return self._eps


def _variance_partition(model, span, n_real=600, seed=7):
    """Split realised gust variance into a per-flight offset and within-flight wander.

    Returns (total, offset_var, within_var, offset_fraction): `offset_var` is the
    variance across realisations of each realisation's MEAN over the window, and
    `within_var` the mean over realisations of the variance WITHIN the window.
    """
    rng = np.random.default_rng(seed)
    f = _SiteTurbField(np.zeros(3))
    coords = np.linspace(0.0, span, 64)
    rows = []
    for _ in range(n_real):
        g = CorrelatedGusts(f, np.random.default_rng(int(rng.integers(1 << 30))),
                            t_max=8.0, model=model)
        rows.append([g.shape(c)[0] for c in coords])
    a = np.asarray(rows)
    total = float(a.var())
    offset = float(a.mean(axis=1).var())
    within = float(a.var(axis=1).mean())
    return total, offset, within, offset / total


def test_temporal_gusts_are_a_frozen_offset_over_a_flight():
    """The physical content of the WP1 correction, measured on realisations.

    Both models are normalised to unit ENSEMBLE variance and over a flight window both
    realise very nearly that (~0.96), so total variance is NOT what separates them.
    What separates them is how that variance is PARTITIONED.

    At the site's correlation scale (T_L = 8.43 s, from GATE G1) a 3.5 s flight is
    shorter than one correlation time, so the temporal model delivers an almost
    constant offset: ~89% of its variance is a per-flight offset, only ~11% is
    within-flight wander. The sweeping model decorrelates over ~12 m of a ~25 m path
    and splits roughly 55/44. A coherent offset pushes the disc the same way for the
    whole flight; within-flight wander partially cancels. That is why the temporal arm
    disperses landings ~1.7x more, and it is the mechanism behind Phi(N).

    Asserts the ORDERING and bounds the partition -- the reproducible content.
    """
    tot_t, off_t, wit_t, frac_t = _variance_partition("temporal", span=3.5)
    tot_s, off_s, wit_s, frac_s = _variance_partition("sweeping", span=25.0)
    print(f"temporal: total {tot_t:.3f} offset {off_t:.3f} within {wit_t:.3f} "
          f"offset_frac {frac_t:.3f}")
    print(f"sweeping: total {tot_s:.3f} offset {off_s:.3f} within {wit_s:.3f} "
          f"offset_frac {frac_s:.3f}")

    # Both must still realise close to their nominal unit variance over a flight.
    assert 0.85 < tot_t < 1.15, f"temporal total variance {tot_t:.3f}"
    assert 0.85 < tot_s < 1.15, f"sweeping total variance {tot_s:.3f}"

    # The discriminating claim.
    assert frac_t > frac_s, (
        f"the temporal model ({frac_t:.3f}) must be MORE dominated by a frozen "
        f"per-flight offset than the sweeping model ({frac_s:.3f}); if this inverts, "
        "the WP1 correction has lost its mechanism")
    assert frac_t > 0.80, f"temporal offset fraction {frac_t:.3f} (expected ~0.89)"
    assert frac_s < 0.70, f"sweeping offset fraction {frac_s:.3f} (expected ~0.55)"


def test_realised_autocorrelation_matches_the_analytic_one():
    """Sampled xi(t) must actually have the autocorrelation the mode weights claim.

    `autocorrelation()` returns the analytic ensemble ACF from the deterministic
    weights. That is a real check on the Shinozuka-Deodatis discretisation, but it is
    not evidence that a realisation behaves. This closes the loop empirically.
    """
    f = ConstantTurbField(np.zeros(3))
    g0 = CorrelatedGusts(f, np.random.default_rng(11), t_max=8.0, model="temporal")
    scale = g0.scale
    lags = np.linspace(0.0, 2.0 * scale, 9)
    analytic = g0.autocorrelation(lags)

    rng = np.random.default_rng(3)
    span, n_t = 60.0 * scale, 4096
    t = np.linspace(0.0, span, n_t)
    dt = span / (n_t - 1)
    acc = np.zeros_like(lags)
    n_real = 200
    for _ in range(n_real):
        g = CorrelatedGusts(f, np.random.default_rng(int(rng.integers(1 << 30))),
                            t_max=8.0, model="temporal")
        xi = np.array([g.shape(tt)[0] for tt in t])
        xi = xi - xi.mean()
        for i, lag in enumerate(lags):
            shift = int(round(lag / dt))
            if shift == 0:
                acc[i] += float((xi * xi).mean())
            else:
                acc[i] += float((xi[:-shift] * xi[shift:]).mean())
    empirical = acc / n_real
    empirical = empirical / empirical[0]

    np.testing.assert_allclose(
        empirical, analytic / analytic[0], atol=0.12,
        err_msg="the realised autocorrelation departs from the analytic mode-weight "
                "prediction")


# --------------------------------------------------------------------------------
# F7: the trim-speed deviation that was printed but never guarded
# --------------------------------------------------------------------------------

HUMMEL_TRIM_SPEED_MS = 9.1        # Hummel (2003) section 2.4.2, analytic trim glide
TRIM_SPEED_TOLERANCE = 0.12       # see the docstring below


def test_trim_speed_against_hummel_is_bounded():
    """The 6-DOF trim speed sits 9% below Hummel's published 9.1 m/s.

    Every other anchor in the Case B trim glide agrees to better than 0.5% (trim alpha
    9.00 vs 9 deg, glide angle -19.32 vs -19.3 deg, C_L 0.672 vs 0.67), so the
    deviation is specific to the speed and worth bounding rather than ignoring. The
    tolerance is 12%: wide enough to admit the known -9.0%, tight enough that a
    regression of the kind that produced the original five coefficient errors -- each
    of which moved things by tens of percent -- would trip it.

    This is a BOUND on a known, reported deviation, not a claim of agreement. The
    deviation itself is stated in the paper's validation table.
    """
    sys.path.insert(0, str(HERE.parent / "P1_AnalyticalModel"))
    import coefficients as C
    import validation_hummel as VH

    trim = VH.trim_glide(C.get_coefficients("potts_tabulated"))
    assert trim is not None, "trim_glide found no zero-pitching-moment angle"
    v = float(trim["v_trim"])

    # The anchors that DO agree, kept beside the one that does not, so a regression
    # cannot hide by moving them all together.
    assert abs(trim["alpha_trim_deg"] - 9.0) < 0.1, trim
    assert abs(trim["gamma_deg"] - (-19.3)) < 0.1, trim
    assert abs(trim["CL"] - 0.67) < 0.01, trim

    rel = abs(v - HUMMEL_TRIM_SPEED_MS) / HUMMEL_TRIM_SPEED_MS
    print(f"trim speed {v:.3f} m/s vs Hummel {HUMMEL_TRIM_SPEED_MS} "
          f"({100 * (v - HUMMEL_TRIM_SPEED_MS) / HUMMEL_TRIM_SPEED_MS:+.1f}%)")
    assert rel < TRIM_SPEED_TOLERANCE, (
        f"trim speed {v:.2f} m/s is {100 * rel:.1f}% from Hummel's "
        f"{HUMMEL_TRIM_SPEED_MS} m/s, outside the "
        f"{100 * TRIM_SPEED_TOLERANCE:.0f}% bound this deviation was recorded at")


def test_only_the_tabulated_coefficient_set_passes_its_own_band():
    """The paper's numbers rest on ONE of four coefficient sets. Pin that.

    `validation_hummel.py` exits 1 by design, because three of its four sets fail the
    +/-15% band -- that non-zero exit is the script working, not breaking. What must
    hold is that `potts_tabulated`, the set the manuscript uses, is the one that
    passes. If a future change makes another set pass, or this one fail, the choice of
    default has to be revisited rather than silently inherited.
    """
    p1 = HERE.parent / "P1_AnalyticalModel"
    r = subprocess.run([sys.executable, "validation_hummel.py"], cwd=p1,
                       capture_output=True, text=True, timeout=900)
    out = r.stdout
    assert "AGREEMENT WITH PUBLISHED VALUES" in out, \
        "validation_hummel.py no longer scores itself against the published values"

    block = out.split("AGREEMENT WITH PUBLISHED VALUES")[1].split("LATERAL")[0]
    verdicts = {}
    for ln in block.splitlines():
        for name in ("hummel_short", "hummel_long", "potts_crowther",
                     "potts_tabulated"):
            if ln.strip().startswith(name):
                verdicts[name] = "PASS" if "PASS:" in ln else "FAIL"
    assert len(verdicts) == 4, f"expected 4 coefficient sets, parsed {verdicts}"
    assert verdicts["potts_tabulated"] == "PASS", (
        f"the coefficient set the manuscript uses no longer passes: {verdicts}")
    assert [k for k, v in verdicts.items() if v == "PASS"] == ["potts_tabulated"], (
        f"the set of passing coefficient sets changed: {verdicts}. The default is "
        "chosen because it is the only one reproducing the published flights; "
        "re-examine that choice.")
    # Exit code is 1 precisely because the other three fail. Pin the contract.
    assert r.returncode == 1, (
        f"validation_hummel.py returned {r.returncode}; it returns 1 while any "
        "coefficient set fails its band, and three of four do")


# --------------------------------------------------------------------------------
# The re-graded gate
# --------------------------------------------------------------------------------

def test_phi_law_reproduces_the_briefs_own_values():
    """Phi(N) must reproduce the brief's two internally consistent values.

    External anchor: these were computed by the brief's author from the same closed
    form, independently of this implementation. (Its third pair, N = 11.3 with
    sqrt(Phi) = 0.46, is inconsistent with its own formula -- 0.46 corresponds to
    N = 10.9 -- and is recorded as such in regrade_gate_g1._verify_limits.)
    """
    import regrade_gate_g1 as rg
    assert abs(float(np.sqrt(rg.phi(4.0)[0])) - 0.67) < 0.005
    assert abs(float(np.sqrt(rg.phi(5.6)[0])) - 0.60) < 0.005
    # Analytic limits.
    assert abs(float(rg.phi(1e-9)[0]) - 1.0) < 1e-8
    big = 1e4
    assert abs(float(rg.phi(big)[0]) - 8.0 / (3.0 * big)) / (8.0 / (3.0 * big)) < 1e-3
    # Monotone, so the inversion has a unique root.
    grid = np.linspace(1e-3, 50.0, 5000)
    assert np.all(np.diff(rg.phi(grid)) < 0)


def test_regrade_pairs_cells_by_identity_not_by_position():
    """The two arms are stored in worker-completion order and must be realigned.

    Discriminating negative case: shuffling one arm must not change the answer. A
    positional pairing would fail this.
    """
    import regrade_gate_g1 as rg
    src = HERE / "gust_model_comparison.json"
    if not src.exists():
        pytest.skip("gust_model_comparison.json not present")
    d = json.loads(src.read_text())

    a = rg._cells_to_arrays(d["temporal"]["cells"])
    rng = np.random.default_rng(0)
    shuffled = list(d["temporal"]["cells"])
    rng.shuffle(shuffled)
    b = rg._cells_to_arrays(shuffled)

    assert a["key"] == b["key"]
    np.testing.assert_array_equal(a["sum_sq"], b["sum_sq"])

    # And the arms really are stored in different orders in the shipped file, which is
    # what makes the realignment necessary rather than merely tidy.
    raw_t = [rg._cell_key(r) for r in d["temporal"]["cells"]]
    raw_s = [rg._cell_key(r) for r in d["sweeping"]["cells"]]
    assert set(raw_t) == set(raw_s)
    assert raw_t != raw_s, (
        "the arms are now stored in the same order; the realignment in "
        "_cells_to_arrays is still correct but this test no longer discriminates")


def test_realised_N_is_not_any_of_the_briefs_conventions():
    """The reconciliation: N is measured from the run, not chosen from a convention.

    The brief's two contradictory acceptance bands both fixed N by picking an
    integral-scale convention (L = 6.2 / 4.4 / 2.2 m, N = 4.0 / 5.6 / 11.3) and then
    disagreed about which was meant -- a spread the brief itself calls "45% from
    definitional choices alone". The run records what it actually used: a per-cell
    `scale`, and a launch-to-landing chord. Neither matches any convention, so the
    argument was never resolvable by choosing better -- it had to be measured.
    """
    import regrade_gate_g1 as rg
    src = HERE / "gust_model_comparison.json"
    if not (src.exists() and REAL_ARTIFACT.exists()):
        pytest.skip("artifacts not present")
    d = json.loads(src.read_text())
    N, chords, scales = rg.realised_N(d["sweeping"]["cells"], REAL_ARTIFACT)

    assert len(N) == 108
    # The chord must look like a disc throw, not like a target offset.
    assert 20.0 < chords.mean() < 35.0, f"chord mean {chords.mean():.2f} m"
    # The realised correlation length is NOT any of the conventional values.
    for conventional_L in (6.2, 4.4, 2.2):
        assert abs(scales.mean() - conventional_L) > 0.5, (
            f"realised scale {scales.mean():.2f} m now matches the conventional "
            f"{conventional_L} m; the reconciliation's premise has changed")
    assert 2.5 < scales.mean() < 4.5, f"realised scale mean {scales.mean():.2f} m"
    assert 5.0 < N.mean() < 13.0, f"realised N mean {N.mean():.2f}"


def test_phi_predicts_the_isotropic_arm_within_the_briefs_own_tolerance():
    """Phi(N) at the realised N predicts the ISOTROPIC measurement to a few percent.

    Phi is derived for an isotropic exponentially-correlated field, so the isotropic
    arm is the one it is entitled to predict; the anisotropic arm carries a correction
    Phi does not model. Both are checked against the brief's own GATE G2 tolerance
    (line 622, "within ~20%"), and the isotropic arm is additionally required to agree
    far more tightly than that -- which is what makes this a real test rather than a
    band wide enough to admit anything.
    """
    import regrade_gate_g1 as rg
    iso = HERE / "gust_model_comparison_isotropic.json"
    ani = HERE / "gust_model_comparison.json"
    if not (iso.exists() and ani.exists() and REAL_ARTIFACT.exists()):
        pytest.skip("artifacts not present")

    results = {}
    for name, path in (("isotropic", iso), ("anisotropic", ani)):
        d = json.loads(path.read_text())
        R_pred, N, _, _ = rg.pooled_phi_prediction(
            d["temporal"]["cells"], d["sweeping"]["cells"], REAL_ARTIFACT)
        R_meas = d["gate_G1"]["ratio"]
        dev = R_meas / R_pred - 1.0
        results[name] = (R_pred, R_meas, dev)
        print(f"{name:12s} predicted {R_pred:.3f}  measured {R_meas:.3f}  "
              f"deviation {100 * dev:+.1f}%")

    for name, (_, _, dev) in results.items():
        assert abs(dev) <= rg.RECONCILED_TOLERANCE, (
            f"{name}: Phi(N) at the realised N is {100 * dev:+.1f}% from the "
            f"measurement, outside the brief's own +/-20% GATE G2 tolerance")

    # The discriminating part: the isotropic arm is the one Phi should get right.
    assert abs(results["isotropic"][2]) < 0.06, (
        f"Phi(N) no longer predicts the isotropic arm to better than 6% "
        f"({100 * results['isotropic'][2]:+.1f}%); the reconciliation rested on it")
    assert abs(results["isotropic"][2]) < abs(results["anisotropic"][2]), (
        "Phi(N) should fit the isotropic arm BETTER than the anisotropic one, since "
        "it carries no anisotropy correction; that ordering has inverted")
