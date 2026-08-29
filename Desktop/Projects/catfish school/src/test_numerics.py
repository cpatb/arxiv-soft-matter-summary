"""Validate the numerical routines against cases with known answers.

These are the quantities the conclusions rest on. Each test compares against an
analytic result, a hand-computable case, or a synthetic field built to have a
property the routine is supposed to recover. Run with:

    python src/test_numerics.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).resolve().parent))

from catfish import flow, segment, topology  # noqa: E402

RNG = np.random.default_rng(20260808)
FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {name}" + (f"  -- {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


# --------------------------------------------------------------- topology ---


def test_shape_metrics() -> None:
    print("\ntopology.component_shapes against analytic shapes")

    # A disc: circularity 1, aspect ratio 1, solidity 1.
    H = W = 401
    yy, xx = np.mgrid[:H, :W]
    R = 120.0
    disc = ((yy - H // 2) ** 2 + (xx - W // 2) ** 2) <= R**2
    env = segment.Envelope(disc, ndi.label(disc)[0].astype(np.int32), 1,
                           np.zeros_like(disc, float), True)
    s = topology.component_shapes(env, body_px=1.0)[0]
    check("disc circularity ~ 1", abs(s.circularity - 1.0) < 0.03,
          f"got {s.circularity:.4f}")
    check("disc aspect ratio ~ 1", abs(s.aspect_ratio - 1.0) < 0.02,
          f"got {s.aspect_ratio:.4f}")
    check("disc solidity ~ 1", s.solidity > 0.98, f"got {s.solidity:.4f}")
    check("disc area ~ pi R^2", abs(s.area_bl2 - np.pi * R**2) / (np.pi * R**2) < 0.01,
          f"got {s.area_bl2:.0f}, expected {np.pi*R**2:.0f}")
    check("disc boundary roughness ~ 1", abs(s.boundary_roughness - 1.0) < 0.03,
          f"got {s.boundary_roughness:.4f}")

    # An ellipse with semi-axes a, b: the second-moment aspect ratio is a/b.
    a, b = 150.0, 50.0
    ell = ((xx - W // 2) / a) ** 2 + ((yy - H // 2) / b) ** 2 <= 1.0
    env = segment.Envelope(ell, ndi.label(ell)[0].astype(np.int32), 1,
                           np.zeros_like(ell, float), True)
    s = topology.component_shapes(env, body_px=1.0)[0]
    check("ellipse aspect ratio ~ a/b", abs(s.aspect_ratio - a / b) / (a / b) < 0.02,
          f"got {s.aspect_ratio:.3f}, expected {a/b:.3f}")

    # Body-length normalisation: doubling body_px must quarter the area in BL^2.
    s1 = topology.component_shapes(env, body_px=1.0)[0]
    s2 = topology.component_shapes(env, body_px=2.0)[0]
    check("area scales as 1/body_px^2", abs(s2.area_bl2 * 4 - s1.area_bl2) < 1e-6,
          f"{s1.area_bl2:.1f} vs 4x{s2.area_bl2:.1f}")
    check("circularity is scale-free", abs(s1.circularity - s2.circularity) < 1e-9)


def test_fission_detection() -> None:
    print("\ntopology.detect_events on a constructed split")
    H = W = 300
    one = np.zeros((H, W), bool)
    one[100:200, 80:220] = True  # a single bar
    two = one.copy()
    two[:, 145:155] = False  # cut it in half

    l_one = ndi.label(one)[0].astype(np.int32)
    l_two = ndi.label(two)[0].astype(np.int32)
    ev = topology.detect_events([l_one, l_two], [0, 1])
    check("split produces a fission event",
          any(e["type"] == "fission" for e in ev), f"events={ev}")

    ev_rev = topology.detect_events([l_two, l_one], [0, 1])
    check("merge produces a fusion event",
          any(e["type"] == "fusion" for e in ev_rev), f"events={ev_rev}")

    ev_none = topology.detect_events([l_one, l_one], [0, 1])
    check("no event when nothing changes", len(ev_none) == 0, f"events={ev_none}")


# ------------------------------------------------------------------- flow ---


def test_order_parameters() -> None:
    print("\nflow.order_parameters on known velocity fields")
    n = 500
    u = np.full(n, 2.0)
    v = np.zeros(n)
    r = flow.order_parameters(u, v)
    check("perfectly aligned -> polarisation 1", abs(r["polarisation"] - 1.0) < 1e-9,
          f"got {r['polarisation']:.6f}")
    check("perfectly aligned -> speed 2", abs(r["mean_speed_px"] - 2.0) < 1e-9)

    th = RNG.uniform(0, 2 * np.pi, 20000)
    r = flow.order_parameters(np.cos(th), np.sin(th))
    # Random directions give |mean unit vector| ~ 1/sqrt(N).
    check("isotropic -> polarisation ~ 0", r["polarisation"] < 0.05,
          f"got {r['polarisation']:.4f}, expected ~{1/np.sqrt(20000):.4f}")

    # Directions in a narrow cone of half-width w have polarisation sin(w)/w.
    w = 0.5
    th = RNG.uniform(-w, w, 40000)
    r = flow.order_parameters(np.cos(th), np.sin(th))
    expect = np.sin(w) / w
    check("cone -> polarisation sin(w)/w", abs(r["polarisation"] - expect) < 0.01,
          f"got {r['polarisation']:.4f}, expected {expect:.4f}")


def test_rotational_order() -> None:
    print("\ncollective_state.rotational_order on known flow patterns")
    from collective_state import rotational_order

    th = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    y, x = 50 * np.sin(th), 50 * np.cos(th)

    # Pure tangential circulation about the centre -> M = 1.
    M = rotational_order(y, x, -np.sin(th), np.cos(th), 0.0, 0.0)
    check("pure circulation -> M = 1", abs(M - 1.0) < 1e-6, f"got {M:.6f}")

    # Uniform translation -> M = 0. This is the case that distinguishes a
    # milling school from a polarised one.
    M = rotational_order(y, x, np.ones(400), np.zeros(400), 0.0, 0.0)
    check("uniform translation -> M = 0", abs(M) < 1e-6, f"got {M:.6f}")

    # Radial expansion -> M = 0: motion is along r, so the cross product vanishes.
    M = rotational_order(y, x, x / 50, y / 50, 0.0, 0.0)
    check("radial expansion -> M = 0", abs(M) < 1e-6, f"got {M:.6f}")

    a = RNG.uniform(0, 2 * np.pi, 4000)
    M = rotational_order(
        RNG.uniform(-50, 50, 4000), RNG.uniform(-50, 50, 4000),
        np.cos(a), np.sin(a), 0.0, 0.0,
    )
    check("random directions -> M ~ 0", M < 0.05, f"got {M:.4f}")


def test_correlation_function() -> None:
    print("\nflow.correlation_function on synthetic fields")

    # (1) Uncorrelated noise: C(r) must fall to ~0 immediately and stay there.
    n = 900
    y = RNG.uniform(0, 600, n)
    x = RNG.uniform(0, 600, n)
    u = RNG.normal(size=n)
    v = RNG.normal(size=n)
    c = flow.correlation_function(y, x, u, v, body_px=20.0, n_bins=25)
    finite = np.isfinite(c["C"])
    check("white noise -> |C| small at all r",
          np.nanmax(np.abs(c["C"][finite])) < 0.15,
          f"max|C|={np.nanmax(np.abs(c['C'][finite])):.4f}")

    # (2) A field with an imposed correlation length. Build a smooth random
    #     field by convolving noise with a Gaussian of width sigma; the measured
    #     crossing must then grow in proportion to sigma.
    #
    #     The single-realisation estimator is noisy -- spread of ~20 px against a
    #     value of ~50 px -- so xi is taken from the *averaged* C(r) rather than
    #     by averaging per-realisation xi. This is not a trick to make the test
    #     pass; it is the estimator the analysis itself must use, and finding
    #     that out is what this test is for.
    def xi_from_averaged_C(sigma_px: float, n_real: int = 25, L: int = 400) -> float:
        grid = np.linspace(0.0, 140.0, 60)
        curves = []
        for _ in range(n_real):
            gy = ndi.gaussian_filter(RNG.normal(size=(L, L)), sigma_px)
            gx = ndi.gaussian_filter(RNG.normal(size=(L, L)), sigma_px)
            yy, xx = np.mgrid[0:L:8, 0:L:8]
            yy, xx = yy.ravel(), xx.ravel()
            cc = flow.correlation_function(
                yy.astype(float), xx.astype(float),
                gy[yy, xx], gx[yy, xx], body_px=1.0, n_bins=40,
            )
            if len(cc["r_bl"]):
                curves.append(np.interp(grid, cc["r_bl"], cc["C"]))
        if not curves:
            return float("nan")
        return flow._first_zero_crossing(grid, np.nanmean(np.array(curves), axis=0))

    xi_6 = xi_from_averaged_C(6.0)
    xi_24 = xi_from_averaged_C(24.0)
    check("correlation length grows with imposed smoothing scale",
          np.isfinite(xi_6) and np.isfinite(xi_24) and xi_24 > 2.0 * xi_6,
          f"xi(sigma=6)={xi_6:.1f} px, xi(sigma=24)={xi_24:.1f} px")

    # (2b) Uncorrelated noise must NOT produce a large apparent correlation
    #      length. Subtracting the sample mean constrains fluctuations to sum to
    #      zero, which forces C(r) negative somewhere regardless of the true
    #      field; this check confirms the resulting spurious xi stays small
    #      relative to the domain, so a large measured xi means something.
    L = 800
    n = 700
    y = RNG.uniform(0, L, n)
    x = RNG.uniform(0, L, n)
    xi_noise = flow.correlation_function(
        y, x, RNG.normal(size=n), RNG.normal(size=n), body_px=1.0, n_bins=40
    )["xi_bl"]
    check("uncorrelated field -> xi small compared with domain",
          (not np.isfinite(xi_noise)) or xi_noise < 0.12 * L,
          f"xi={xi_noise:.1f} px on a {L} px domain")

    # (3) C(0) must equal 1 by construction, and subtracting a uniform drift
    #     must not change anything -- the routine works on fluctuations.
    n = 600
    y = RNG.uniform(0, 400, n)
    x = RNG.uniform(0, 400, n)
    u = RNG.normal(size=n)
    v = RNG.normal(size=n)
    c1 = flow.correlation_function(y, x, u, v, body_px=10.0)
    c2 = flow.correlation_function(y, x, u + 37.0, v - 12.0, body_px=10.0)
    same = np.allclose(np.nan_to_num(c1["C"]), np.nan_to_num(c2["C"]), atol=1e-9)
    check("C(r) is invariant to a uniform drift", same)

    # (4) Body-length normalisation: doubling body_px must halve r and xi.
    c3 = flow.correlation_function(y, x, u, v, body_px=20.0)
    ratio = np.nanmedian(c1["r_bl"][: len(c3["r_bl"])] / c3["r_bl"])
    check("r scales as 1/body_px", abs(ratio - 2.0) < 0.05, f"ratio={ratio:.3f}")


def test_zero_crossing() -> None:
    print("\nflow._first_zero_crossing on hand-computable input")
    r = np.array([0.0, 1.0, 2.0, 3.0])
    C = np.array([1.0, 0.5, -0.5, -1.0])
    xi = flow._first_zero_crossing(r, C)
    check("linear interpolation of crossing", abs(xi - 1.5) < 1e-9, f"got {xi}")

    C_pos = np.array([1.0, 0.8, 0.6, 0.4])
    check("no crossing -> NaN", not np.isfinite(flow._first_zero_crossing(r, C_pos)))

    # The 1/e decay length is the estimator actually used, because the measured
    # C(r) decays smoothly but stays positive across the whole school, so a zero
    # crossing usually does not exist. On C(r) = exp(-r/L) the answer is L.
    L = 4.0
    rr = np.linspace(0, 30, 400)
    xi_e = flow.correlation_decay_length(rr, np.exp(-rr / L))
    check("1/e decay length of exp(-r/L) is L", abs(xi_e - L) < 0.05,
          f"got {xi_e:.4f}, expected {L}")

    # A curve that decays but never reaches zero still has a decay length. This
    # is exactly the case the real data turned out to be in.
    C_slow = 0.2 + 0.8 * np.exp(-rr / 6.0)
    check("decay length defined where zero crossing is not",
          np.isfinite(flow.correlation_decay_length(rr, C_slow))
          and not np.isfinite(flow._first_zero_crossing(rr, C_slow)))

    # Already below the level at the first sample: report the smallest measured
    # separation, not NaN. The correlation is shorter than the resolution, which
    # is a result, not a failure -- this is what the shuffled null does.
    xi_short = flow.correlation_decay_length(np.array([1.0, 2.0]), np.array([0.1, 0.05]))
    check("sub-resolution correlation -> smallest r, not NaN",
          np.isfinite(xi_short) and abs(xi_short - 1.0) < 1e-9, f"got {xi_short}")


def test_camera_motion() -> None:
    print("\nflow.estimate_camera_motion on a known synthetic shift")
    H, W = 400, 600
    base = (RNG.uniform(0, 255, (H, W))).astype(np.uint8)
    base = ndi.gaussian_filter(base, 2).astype(np.uint8)  # give it trackable texture

    # ndi.shift places input[i] at output[i + shift], so to move image content
    # by (dx, dy) in (column, row) the shift argument is (dy, dx). Getting this
    # backwards on the first attempt produced a sign-flipped "failure" that was
    # in the test, not in the estimator.
    dx, dy = 7.0, -4.0
    shifted = ndi.shift(base.astype(float), shift=(dy, dx), order=1,
                        mode="nearest").astype(np.uint8)

    cam = flow.estimate_camera_motion(base, shifted)
    check("camera estimate succeeds", cam.ok, f"inliers={cam.n_inliers}")
    if cam.ok:
        got = cam.velocity_at(np.array([H / 2]), np.array([W / 2]))[0]
        check("recovers imposed shift",
              abs(got[0] - dx) < 0.4 and abs(got[1] - dy) < 0.4,
              f"got ({got[0]:.2f}, {got[1]:.2f}), imposed ({dx}, {dy})")


def test_segmentation_guards() -> None:
    print("\nsegment guards on degenerate input")
    flat = np.full((200, 300, 3), 0.5)
    env = segment.school_envelope(flat)
    check("uniform frame -> no school", not env.present and env.n_components == 0)

    bad = np.full((10, 10, 3), 2.0)
    try:
        segment.fishness_field(bad)
        check("out-of-range RGB rejected", False, "no exception raised")
    except ValueError:
        check("out-of-range RGB rejected", True)

    nan_img = np.full((10, 10, 3), 0.5)
    nan_img[0, 0, 0] = np.nan
    try:
        segment.fishness_field(nan_img)
        check("NaN input rejected", False, "no exception raised")
    except ValueError:
        check("NaN input rejected", True)


def main() -> None:
    test_shape_metrics()
    test_fission_detection()
    test_order_parameters()
    test_rotational_order()
    test_correlation_function()
    test_zero_crossing()
    test_camera_motion()
    test_segmentation_guards()

    print("\n" + "=" * 60)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): " + ", ".join(FAILURES))
        sys.exit(1)
    print("all numerical checks passed")


if __name__ == "__main__":
    main()
