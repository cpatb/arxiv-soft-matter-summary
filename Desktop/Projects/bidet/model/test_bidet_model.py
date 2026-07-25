"""Validation of bidet_model.py — analytic limits, paper cross-check, monotonicity,
and the CFD-fitted turbulent closure (Stage 2)."""
import os, json
import numpy as np
from bidet_model import (LAMINAR, laminar_tauII_direct, optimum, calibrate_turbulent,
                         cfd_closure, film_tau_numeric, kappa_sinusoid, V0_from_dP,
                         flow_rate, CD, NU)

NOM = dict(dP=0.07e6, r_max=0.015, tau_c=5.0)

_CFD_PATH = os.path.join(os.path.dirname(__file__), "..", "cfd_openfoam", "cfd_fit_summary.json")
def _cfd():
    return json.load(open(_CFD_PATH))


def test_common_form_matches_paper_tauII():
    V0, a = 11.48, 0.4e-3
    r = np.array([3e-3, 6e-3, 12e-3])
    direct = laminar_tauII_direct(r, V0, a)
    common = LAMINAR.tau_w(r, V0, 2 * a)   # d = 2a
    assert np.allclose(direct, common, rtol=1e-10), "common-form C_tau must equal c2/64"


def test_reproduces_paper_laminar_optimum():
    o = optimum(LAMINAR, **NOM)
    assert abs(o["d_star"] * 1e3 - 0.80) < 0.02, f"d* should be 0.80 mm, got {o['d_star']*1e3:.3f}"
    Q_noCd = o["Q"] / CD                    # paper drops C_d
    assert abs(Q_noCd * 6e4 - 0.347) < 0.02, f"Q*(no Cd) should be 0.347 L/min, got {Q_noCd*6e4:.3f}"
    assert abs(o["V0"] - 11.48) < 0.05


def test_laminar_pressure_independence():
    q = [optimum(LAMINAR, dP=p, r_max=0.015, tau_c=5.0)["Q"] for p in (0.05e6, 0.10e6, 0.15e6)]
    assert max(q) / min(q) < 1.001, "laminar Q* must be dP-independent (its headline result)"


def test_numeric_film_recovers_r5_slope():
    V0, a = 11.48, 0.4e-3
    r = np.linspace(4e-3, 14e-3, 6)
    tau, r0 = film_tau_numeric(r, V0, a)
    slope = np.polyfit(np.log(r), np.log(tau), 1)[0]
    assert -6.0 < slope < -3.5, f"far-field film shear should steepen toward r^-5, slope={slope:.2f}"


def test_turbulent_calibration_reproduces_commercial():
    comm = dict(d=1.0e-3, dP=0.07e6, H=60e-3, theta=np.deg2rad(43))
    cl = calibrate_turbulent(1.05, comm, r_max=0.015, tau_c=5.0)
    # by construction the commercial nozzle's tau_w at r_max equals tau_c
    V0 = V0_from_dP(comm["dP"]) * np.cos(comm["theta"])
    B = 6.0; Vc = B * V0 * comm["d"] / comm["H"]; d_eff = comm["H"] / B
    Re_d = V0_from_dP(comm["dP"]) * comm["d"] / 1  # placeholder
    tau = cl.C_tau * 1000 * Vc**2 * (V0_from_dP(comm["dP"]) * comm["d"] / 8.9e-7)**cl.p * (d_eff / 0.015)**cl.n
    assert abs(tau - 5.0) / 5.0 < 1e-6, "calibration must make commercial just-clean at tau_c"


def test_standoff_sign_laminar_wrong_turbulent_right():
    # turbulent closure: more standoff -> more water (physical). laminar Region II: wrong sign.
    turb = calibrate_turbulent(1.05, dict(d=1.0e-3, dP=0.07e6, H=60e-3, theta=0.0),
                               r_max=0.015, tau_c=5.0)
    q_close = optimum(turb, dP=0.07e6, r_max=0.015, tau_c=5.0, H=6 * 1.0e-3)["Q"]
    q_far = optimum(turb, dP=0.07e6, r_max=0.015, tau_c=5.0, H=60e-3)["Q"]
    assert q_far > q_close, "turbulent: standoff must COST water (right sign)"
    ql_close = optimum(LAMINAR, dP=0.07e6, r_max=0.015, tau_c=5.0, H=6 * 0.8e-3)["Q"]
    ql_far = optimum(LAMINAR, dP=0.07e6, r_max=0.015, tau_c=5.0, H=60e-3)["Q"]
    assert ql_far < ql_close, "laminar Region II gives the (unphysical) wrong standoff sign"


def test_cfd_closure_form():
    """cfd_closure encodes the Phares form: C_tau=C_P, p=-1/2, n=n_t."""
    cl = cfd_closure(1.35, 1.90)
    assert cl.C_tau == 1.35 and cl.p == -0.5 and cl.n == 1.90


def test_cfd_summary_exists_and_sane():
    j = _cfd()
    n = j["nominal"]
    # computed exponent is far from laminar 5 (reversal holds) but not the assumed 1.05
    assert 1.5 < n["n_t"] < 2.3, f"CFD n_t={n['n_t']} outside expected computed range"
    assert n["r2"] > 0.99, "nominal power-law fit must be good"
    assert n["C_P"] > 0


def test_cfd_pressure_exponent_in_range():
    """Q* ~ dP^pexp with the CFD closure: weak negative, ~ -0.4..-0.7 (CFD-derived)."""
    cl = cfd_closure(_cfd()["nominal"]["C_P"], _cfd()["nominal"]["n_t"])
    dPs = np.array([0.05, 0.07, 0.10, 0.15]) * 1e6
    q = np.array([optimum(cl, dP=p, r_max=0.015, tau_c=5.0)["Q"] for p in dPs])
    pexp = np.polyfit(np.log(dPs), np.log(q), 1)[0]
    assert -0.75 < pexp < -0.40, f"CFD pressure exponent {pexp:.3f} outside CFD-derived range"


def test_cfd_standoff_lever_is_weak():
    """CORRECTED finding: CFD tau_w(r_max) is nearly standoff-independent, so the
    flow to clean r_max varies <~1.3x over H/d 1-20 -- NOT the paper's 4.6x lever."""
    so = _cfd()["standoff_Re10000"]
    tau = [so[k]["tau15"] for k in ("1", "6", "20")]   # reliable-fit standoffs
    mult = [(so["6"]["tau15"] / t) ** (1 / 1.5) for t in tau]
    span = max(mult) / min(mult)
    assert span < 1.3, f"CFD standoff flow span {span:.2f} should be weak (<1.3), got strong"


def test_cfd_tau_rmax_corroborates_tau_c():
    """Independent CFD tau_w(r_max) at the commercial point is within 2x of the
    assumed tau_c=5 Pa -> the commercial-datum calibration is corroborated, not circular."""
    tau_cfd = _cfd()["nominal"]["tau_at_15mm_Pa"]
    assert 0.5 * 5.0 < tau_cfd < 2.0 * 5.0, f"CFD tau_w(r_max)={tau_cfd} not near tau_c=5"


def test_kappa_value():
    assert abs(kappa_sinusoid() - 0.6366) < 1e-3


def test_constraints_satisfied_at_optimum():
    o = optimum(LAMINAR, **NOM)
    assert o["Re_d"] > 3000, "turbulent-impingement constraint C1"
    assert o["We"] > 1, "coherence"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    for fn in fns:
        try:
            fn(); print(f"[PASS] {fn.__name__}"); passed += 1
        except AssertionError as e:
            print(f"[FAIL] {fn.__name__}: {e}")
        except Exception as e:
            print(f"[ERR ] {fn.__name__}: {type(e).__name__}: {e}")
    print(f"\n{passed}/{len(fns)} checks passed")
