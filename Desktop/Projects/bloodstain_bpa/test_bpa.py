"""Correctness tests for bpa.py -- run with: python -m pytest test_bpa.py -q
(or plain `python test_bpa.py` for a lightweight runner)."""
import numpy as np

import bpa


def test_spread_monotonic_and_invertible():
    d0 = 3e-3
    for v_n in [1.0, 3.0, 6.0, 10.0]:
        we, re = bpa._we_re_from_vn(v_n, d0)
        xi = float(bpa.spread_factor(we, re))
        v_rec = bpa.vn_from_spread(xi, d0)
        assert np.isclose(v_rec, v_n, rtol=1e-3), (v_n, v_rec)


def test_spread_increases_with_speed():
    d0 = 3e-3
    xis = []
    for v_n in np.linspace(1, 12, 8):
        we, re = bpa._we_re_from_vn(v_n, d0)
        xis.append(float(bpa.spread_factor(we, re)))
    assert np.all(np.diff(xis) > 0)


def test_balthazard_angle():
    # W/L = sin(alpha) by construction.
    s = bpa.Stain(x=0, y=0, z=0, L=2e-3, W=1e-3, azimuth=0, d0=3e-3)
    assert np.isclose(np.degrees(s.alpha), 30.0, atol=1e-9)


def test_no_drag_freefall_matches_analytic():
    # Drop released from rest, no drag: z(t) = z0 - 0.5 g t^2 ; check at t=0.3 s.
    t, X, V = bpa.integrate_flight([0, 0, 1.0], [0, 0, 0], 3e-3,
                                   direction=+1, drag=False, t_max=0.3)
    z_analytic = 1.0 - 0.5 * bpa.G * 0.3**2
    assert np.isclose(X[-1, 2], z_analytic, atol=1e-4), (X[-1, 2], z_analytic)


def test_backward_undoes_forward():
    # Forward then backward integration should return to the start (drag on).
    x0 = np.array([0.0, 0.0, 1.2])
    v0 = np.array([4.0, 1.0, -2.0])
    t, X, V = bpa.integrate_flight(x0, v0, 3e-3, direction=+1, drag=True, t_max=0.25)
    xe, ve = X[-1], V[-1]
    tb, Xb, Vb = bpa.integrate_flight(xe, ve, 3e-3, direction=-1, drag=True, t_max=0.25)
    assert np.allclose(Xb[-1], x0, atol=1e-4), (Xb[-1], x0)


def test_roundtrip_recovers_origin():
    origin = np.array([0.0, 0.0, 1.10])
    stains = bpa.emit_and_land(origin, n_drops=400, speed_range=(3.0, 8.0),
                               d0_range=(2.5e-3, 3.5e-3), cone_half_angle=np.radians(35),
                               floor_z=0.0, seed=1, drag=True)
    assert len(stains) > 40
    x0, y0 = bpa.convergence_2d(stains)
    # Horizontal origin recovered well by convergence (gravity barely shifts it).
    assert np.hypot(x0 - origin[0], y0 - origin[1]) < 0.05
    rec, dists, _ = bpa.ballistic_origin(stains, drag=True)
    # Ballistic height within a few cm of truth.
    assert abs(rec[2] - origin[2]) < 0.08, (rec, origin)


if __name__ == "__main__":
    import traceback
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS {fn.__name__}")
            passed += 1
        except Exception:
            print(f"FAIL {fn.__name__}")
            traceback.print_exc()
    print(f"\n{passed}/{len(fns)} passed")
