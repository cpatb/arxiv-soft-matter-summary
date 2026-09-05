"""Pins the result of fitting the structural family to the Horner 2020 step RISE
(scripts/structural_step_horner.py, first run 2026-09-05): the timescale gate G1 passed
on every donor and the shape gate G2 failed, and neither post-hoc variant rescues it.
Counts are pinned exactly so a change in the fitter, the noise estimate or the data
extraction cannot move the verdict silently. The gates themselves are the script's
constants, imported rather than restated."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import structural_step_horner as ssh  # noqa: E402

from biorheo.constitutive.structural import STRUCTURAL_STEADY, lambda_steady  # noqa: E402
from biorheo.identify import fit_rate_curve  # noqa: E402


@pytest.fixture(scope="module")
def per_donor():
    t, step = ssh.read_matrix("step_300_to_1_human.csv", "stress_Pa_")
    rates, steady = ssh.read_matrix("steady_shear_human.csv", "stress_Pa_")
    out = {}
    for d in sorted(step):
        eta_inf, beta, a = fit_rate_curve(STRUCTURAL_STEADY, rates, steady[d] / rates).values
        v = step[d]
        ok = np.isfinite(v)
        vv, tt = v[ok], t[ok]
        imax = int(np.argmax(vv))
        t_rise, s_rise = tt[:imax + 1], vv[:imax + 1]
        tau, lo, hi, verdict, _ = ssh.fit_tau(t_rise, s_rise, eta_inf, beta, a)
        pred = ssh.stress_rise(t_rise, tau, eta_inf, beta, a)
        rmse = float(np.sqrt(np.mean((pred - s_rise) ** 2)))
        tail_t, tail_s = tt[-ssh.NOISE_TAIL:], vv[-ssh.NOISE_TAIL:]
        noise = float(np.std(tail_s - np.polyval(np.polyfit(tail_t, tail_s, 1), tail_t), ddof=2))
        amp = eta_inf * beta * (lambda_steady(1.0, a) - lambda_steady(300.0, a))
        s_inf, s_0 = steady[d][ssh.RATE_INDEX_1], steady[d][14] / ssh.ANCHOR_RATE
        tau_b, rmse_b, signs_b = ssh.fit_tau_anchored(t_rise, s_rise, s_0, s_inf)
        out[d] = dict(tau=tau, verdict=verdict, t95=tau * ssh.LN20, rmse=rmse, noise=noise,
                      ramp=rmse / amp, signs="".join("+" if r > 0 else "-" for r in s_rise - pred),
                      rmse_b=rmse_b, signs_b=signs_b, n_rise=t_rise.size)
    return out


def test_g1_timescale_gate_passes_on_every_donor(per_donor):
    g1 = [ssh.T95_BAND[0] <= r["t95"] <= ssh.T95_BAND[1] and r["verdict"] in ("measurable", "bounded_only")
          for r in per_donor.values()]
    assert sum(g1) == 20 >= ssh.GATE_MIN_DONORS
    assert all(r["verdict"] == "measurable" for r in per_donor.values())
    taus = np.array([r["tau"] for r in per_donor.values()])
    assert 1.1 < taus.min() and taus.max() < 2.3, (taus.min(), taus.max())


def test_g2_shape_gate_fails_and_post_hoc_variants_do_not_rescue_it(per_donor):
    g2 = sum(r["rmse"] <= ssh.SHAPE_FACTOR * r["noise"] for r in per_donor.values())
    assert g2 == 6 < ssh.GATE_MIN_DONORS
    g2_b = sum(r["rmse_b"] <= ssh.SHAPE_FACTOR * r["noise"] for r in per_donor.values())
    assert g2_b == 5 < ssh.GATE_MIN_DONORS
    # the misfit is a few per cent of the rise amplitude, not a gross one
    ramps = np.array([r["ramp"] for r in per_donor.values()])
    assert 0.01 < ramps.min() and ramps.max() < 0.12 and 0.03 < np.median(ramps) < 0.04


def test_residual_sign_pattern_is_systematic(per_donor):
    # steady-anchored fit: data above the model at 1 s, below at 2 and 3 s, on 15/20
    assert sum(r["signs"][0] == "+" and r["signs"][1:3] == "--" for r in per_donor.values()) == 15
    # data-anchored fit: above at 1 s and at the peak, below between, on 13/20
    assert sum(r["signs_b"].startswith("+-") and r["signs_b"].endswith("+") for r in per_donor.values()) == 13
    assert per_donor["D"]["n_rise"] == 16 and all(r["n_rise"] in (4, 5, 6) for d, r in per_donor.items() if d != "D")
