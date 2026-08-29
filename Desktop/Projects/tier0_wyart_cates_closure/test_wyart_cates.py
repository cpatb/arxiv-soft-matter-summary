"""
Invariant tests for the Tier-0 Wyart-Cates closure.

Runs under pytest, or standalone:  python test_wyart_cates.py
Checks the properties the physics must satisfy (guilty until verified).
"""
import numpy as np

import wyart_cates as wc


def test_fss_bounds_and_monotone():
    sig = np.logspace(-2, 4, 200)
    f = wc.frictional_fraction_ss(sig, sigma_star=10.0)
    assert np.all(f >= 0) and np.all(f < 1.0)
    assert np.all(np.diff(f) >= -1e-12)                 # non-decreasing in sigma
    assert wc.frictional_fraction_ss(0.0, 10.0) == 0.0  # safe at sigma=0
    assert wc.frictional_fraction_ss(1e12, 10.0) > 0.99  # -> 1 at high stress


def test_jamming_fraction_endpoints():
    assert np.isclose(wc.jamming_fraction(0.0), wc.PHI0_DEFAULT)
    assert np.isclose(wc.jamming_fraction(1.0), wc.PHI_M_DEFAULT)
    f = np.linspace(0, 1, 50)
    assert np.all(np.diff(wc.jamming_fraction(f)) <= 0)  # decreasing in f


def test_relative_viscosity_diverges_at_jam():
    phi = 0.585  # in the shear-jamming window (phi_m < phi < phi0)
    # frictionless: phi < phi_J -> finite, >= 1
    assert wc.relative_viscosity(phi, f=0.0) >= 1.0
    # fully frictional: phi_J(1) = phi_m = 0.56 < 0.585 -> jammed -> inf
    assert np.isinf(wc.relative_viscosity(phi, f=1.0))
    # monotone increasing with f up to the jam
    fs = np.linspace(0.0, 0.4, 20)
    eta = wc.relative_viscosity(phi, fs)
    assert np.all(np.diff(eta) >= -1e-9)


def test_flow_curve_dst_vs_cst():
    dst = wc.steady_flow_curve(0.585, 10.0)                 # silica-like
    cst = wc.steady_flow_curve(0.520, 8.0, sigma_yield=4.0)  # corn-starch-like
    assert wc.is_discontinuous(dst) is True
    assert wc.is_discontinuous(cst) is False


def test_state_converges_to_fss():
    t = np.linspace(0, 0.5, 500)
    sig = 15.0
    f = wc.relax_state(t, lambda _t: sig, tau=5e-3, sigma_star=10.0)
    np.testing.assert_allclose(f[-1], wc.frictional_fraction_ss(sig, 10.0),
                               rtol=1e-3, atol=1e-3)
    assert np.all(f >= -1e-9) and np.all(f <= 1.0 + 1e-9)


def test_retention_vs_recovery():
    t = np.linspace(0, 60e-3, 1200)
    pulse = 5e-3

    def sig(tt):
        return 1.0e3 if tt < pulse else 0.5

    f_ret = wc.relax_state(t, sig, tau=2e-3, sigma_star=10.0, f_residual=0.6)
    f_rec = wc.relax_state(t, sig, tau=2e-3, sigma_star=10.0, f_residual=0.0)
    idx = np.searchsorted(t, pulse + 10e-3)
    assert f_ret[idx] > 0.5    # stuck pancake retains the network
    assert f_rec[idx] < 0.1    # reversible material recovers
    # irreversibility: retained value must not have decayed to zero by the end
    assert f_ret[-1] > 0.5


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"PASS  {fn.__name__}")
    print(f"\nAll {len(fns)} invariant tests passed.")
