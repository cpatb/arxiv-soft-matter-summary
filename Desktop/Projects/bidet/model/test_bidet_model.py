"""Validation of bidet_model.py — analytic limits, paper cross-check, monotonicity."""
import numpy as np
from bidet_model import (LAMINAR, laminar_tauII_direct, optimum, calibrate_turbulent,
                         film_tau_numeric, kappa_sinusoid, V0_from_dP, flow_rate, CD)

NOM = dict(dP=0.07e6, r_max=0.015, tau_c=5.0)


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
