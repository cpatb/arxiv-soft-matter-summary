"""test_shape.py -- anchors for the clump shape axis (shape.py).

Every anchor here is either a closed-form limit computed by hand, or a NUMERICAL volume
(Monte-Carlo integration of the union of two balls) that does not pass through
`shape.g_of_u` at all. A test that re-derived the union volume from the same formula
would certify the algebra against itself.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from shape import (
    ALPHA_LADDER,
    alpha_to_u,
    clump_centres,
    clump_volume,
    corrugation,
    g_of_u,
    random_axes,
    span_factor,
    sub_sphere_diameter,
)


# ------------------------------------------------------------------ closed-form limits
def test_g_at_u_zero_is_one_sphere():
    """Fully coincident spheres: the union is ONE sphere, so g(0) = pi/6 exactly."""
    assert g_of_u(0.0) == pytest.approx(math.pi / 6.0, rel=0, abs=1e-15)


def test_g_at_u_one_is_two_disjoint_spheres():
    """Tangent spheres: no overlap, so g(1) = 2*(pi/6) = pi/3 exactly."""
    assert g_of_u(1.0) == pytest.approx(math.pi / 3.0, rel=0, abs=1e-15)


def test_g_is_continuous_at_the_tangency_switch():
    assert g_of_u(0.999999) == pytest.approx(g_of_u(1.000001), abs=1e-10)


def test_g_is_monotone_increasing():
    u = np.linspace(0.0, 1.0, 201)
    v = np.array([g_of_u(float(x)) for x in u])
    assert np.all(np.diff(v) > 0), "union volume must grow as the spheres separate"


# ----------------------------------------------- the anchor that is NOT the same algebra
def _mc_union_volume(ds: float, s: float, n: int = 400_000, seed: int = 7) -> float:
    """Monte-Carlo the union of two balls of diameter `ds`, centres `s` apart.

    Independent of shape.g_of_u: it samples points in a bounding box and counts.
    """
    r = ds / 2.0
    rng = np.random.default_rng(seed)
    lx = s / 2.0 + r
    box = np.array([lx, r, r])
    p = rng.uniform(-1.0, 1.0, size=(n, 3)) * box
    c1 = np.array([+s / 2.0, 0.0, 0.0])
    c2 = np.array([-s / 2.0, 0.0, 0.0])
    inside = (np.linalg.norm(p - c1, axis=1) <= r) | (np.linalg.norm(p - c2, axis=1) <= r)
    return float(inside.mean() * np.prod(2 * box))


@pytest.mark.parametrize("alpha", [1.0, 1.25, 1.5, 1.75, 2.0])
def test_construction_preserves_the_volume_equivalent_diameter(alpha):
    """The whole point of the axis: moving alpha must not move the PSD.

    Checked against a Monte-Carlo volume, so a wrong closed form cannot pass.
    """
    d_eq = 1.0
    target = (math.pi / 6.0) * d_eq ** 3
    assert clump_volume(d_eq, alpha) == pytest.approx(target, rel=1e-12)
    ds = sub_sphere_diameter(d_eq, alpha)
    s = alpha_to_u(alpha) * ds
    mc = _mc_union_volume(ds, s)
    assert mc == pytest.approx(target, rel=0.02), (alpha, mc, target)


def test_alpha_one_is_literally_a_sphere_of_the_requested_diameter():
    assert sub_sphere_diameter(2.7, 1.0) == pytest.approx(2.7, rel=1e-14)
    assert corrugation(2.7, 1.0) == pytest.approx(0.0, abs=1e-15)


def test_alpha_two_sub_sphere_is_d_eq_over_cube_root_two():
    assert sub_sphere_diameter(1.0, 2.0) == pytest.approx(2.0 ** (-1.0 / 3.0), rel=1e-14)


# ------------------------------------------------------------------------- the confound
def test_corrugation_is_zero_only_at_alpha_one_and_rises_monotonically():
    """The interlocking half of the knob. A test that only checked the endpoints would
    pass on a routine that returns 0 everywhere except alpha = 2."""
    al = np.linspace(1.0, 2.0, 101)
    eps = np.array([corrugation(1.0, float(a)) for a in al])
    assert eps[0] == 0.0
    assert np.all(np.diff(eps) > 0)
    assert eps[-1] == pytest.approx(2.0 ** (-1.0 / 3.0), rel=1e-12)


def test_corrugation_is_scale_free():
    """It is a FRACTION, so it must not depend on d_eq -- the static DEM is scale-free
    and a shape descriptor that carried a size would break that."""
    for a in ALPHA_LADDER:
        assert corrugation(1.0, a) == pytest.approx(corrugation(37.0, a), rel=1e-12)


# ------------------------------------------------------------------------ geometry bits
def test_clump_centres_are_separated_by_s_and_centred_on_the_origin():
    c = clump_centres(1.0, 1.5, axis=[0.0, 0.0, 3.0])
    assert np.allclose(c.sum(axis=0), 0.0)
    ds = sub_sphere_diameter(1.0, 1.5)
    assert np.linalg.norm(c[0] - c[1]) == pytest.approx(0.5 * ds, rel=1e-12)


def test_clump_centres_reject_a_zero_axis():
    with pytest.raises(ValueError):
        clump_centres(1.0, 1.5, axis=[0.0, 0.0, 0.0])


def test_alpha_outside_one_to_two_is_refused_not_clipped():
    for bad in (0.99, 2.01, 3.0):
        with pytest.raises(ValueError):
            sub_sphere_diameter(1.0, bad)


def test_random_axes_are_uniform_on_the_sphere_not_on_theta_phi():
    """A (theta, phi) sampler clusters at the poles and would put a preferred direction
    into a campaign whose primary observable is contact-orientation anisotropy.

    The discriminating check is the FABRIC of the axes themselves: uniform directions
    give n*aniso^2 ~ 1, a pole-clustered sampler gives a much larger value.
    """
    from measure import fabric_tensor
    v = random_axes(4000, seed=11)
    assert np.allclose(np.linalg.norm(v, axis=1), 1.0)
    t = 4000 * fabric_tensor(v)["aniso"] ** 2
    assert t < 6.0, f"n*aniso^2 = {t:.2f}; these axes are not isotropic"
    # and the test can fail: a theta-phi sampler on the same n scores far higher.
    rng = np.random.default_rng(3)
    th = rng.uniform(0, math.pi, 4000)
    ph = rng.uniform(0, 2 * math.pi, 4000)
    bad = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], axis=1)
    assert 4000 * fabric_tensor(bad)["aniso"] ** 2 > 50.0


def test_random_axes_are_reproducible_from_the_seed():
    assert np.array_equal(random_axes(50, seed=5), random_axes(50, seed=5))
    assert not np.array_equal(random_axes(50, seed=5), random_axes(50, seed=6))


def test_span_factor_says_an_elongated_grain_needs_a_bigger_cell():
    """N = 1500 clears the periodic collider at alpha = 1; the span grows by 1.59x at
    alpha = 2, so the shape campaign cannot simply inherit it."""
    assert span_factor(1.0) == pytest.approx(1.0, rel=1e-14)
    assert span_factor(2.0) == pytest.approx(2.0 ** (2.0 / 3.0), rel=1e-12)
    vals = [span_factor(a) for a in ALPHA_LADDER]
    assert vals == sorted(vals)
