"""shape.py -- aspect ratio as a control variable ORTHOGONAL to CV.

WHY THIS FILE EXISTS, AND THE MISTAKE IT IS BUILT TO AVOID
----------------------------------------------------------
The 2026-09-13 lesson: a one-parameter family is a PATH through a multi-parameter plane,
and that is exactly how this project manufactured its "bidisperse is non-monotone in CV"
negative -- the family was 50:50 BY NUMBER, so its fines carried 0.500 -> 0.036 of the
VOLUME as CV rose and it never visited the composition where the effect peaks. So the
shape axis is written out as a parameterisation, with the statement of what it holds
fixed, BEFORE any cell is run:

    the grain's volume-equivalent diameter  d_eq   -- drawn from the SAME psd.make_psd,
                                                      so CV, d32, skewness and every
                                                      other PSD descriptor are UNCHANGED
                                                      at every aspect ratio
    the aspect ratio                        alpha  -- the new axis
    everything else about the protocol             -- unchanged

That makes (CV, alpha) a GRID and not a path: moving alpha at fixed CV changes the shape
of the grain and nothing about the size distribution, because the construction is
VOLUME-PRESERVING in closed form.

THE CLUMP GEOMETRY, in closed form
-----------------------------------
A grain is a rigid clump of two identical spheres of diameter `ds`, centres a distance
`s` apart. Write u = s/ds, so the aspect ratio is

    alpha = (s + ds)/ds = 1 + u ,        u in [0, 1]

(alpha = 1 is a single sphere, alpha = 2 is a tangent dumbbell). The union volume is

    V = 2*(pi/6)*ds^3 - V_lens ,
    V_lens = (pi/12)*(2*ds + s)*(ds - s)^2     for s <= ds, and 0 for s >= ds

so V = ds^3 * g(u) with

    g(u) = pi/3 - (pi/12)*(2 + u)*(1 - u)^2    for 0 <= u <= 1
    g(u) = pi/3                                for u >= 1

and the volume-preserving sub-sphere diameter for a target d_eq is

    ds = d_eq * ( (pi/6) / g(u) )^(1/3) .

Checks that follow immediately and are asserted in test_shape.py: g(0) = pi/6 exactly
(the two spheres coincide, so ds = d_eq and the grain IS a sphere); g(1) = pi/3 exactly
(tangent, no overlap, ds = d_eq/2^(1/3)); g is continuous and monotone increasing on
[0, 1].

THE CONFOUND, NAMED AND MEASURED RATHER THAN CAVEATED
------------------------------------------------------
A clump of overlapping spheres has a BUMPY surface, so this axis moves elongation AND
surface interlocking together -- a knob that moves two things at once, which is the
failure mode this project has been burned by repeatedly. It is not hidden here; it is
given a number. At the midplane the clump's surface radius is sqrt(rs^2 - (s/2)^2)
against rs at the poles, so the fractional corrugation depth is

    eps(u) = [rs - sqrt(rs^2 - (s/2)^2)] / r_eq
           = (ds/d_eq) * [1 - sqrt(1 - u^2)]            (rs = ds/2, s = u*ds)

which is 0 at u = 0 and rises to (ds/d_eq) = 2^(-1/3) = 0.7937 at u = 1, where the waist
radius is exactly zero -- two tangent spheres have no neck at all. So the confound is not
a detail: at the top of the ladder the grain's surface is as corrugated as a body can be.
`eps` is reported beside `alpha` on every cell, and the LEVELSET control -- a
smooth body at matched alpha and matched volume, for which eps = 0 BY CONSTRUCTION -- is
what turns the caveat into a measurement. Until that control exists, no result from this
axis may be attributed to elongation rather than to interlocking.

No Yade import here: every function is a pure function of numbers, so it runs and is
tested on Windows.
"""

from __future__ import annotations

import math

import numpy as np

# The aspect ladder. alpha = 1 is the sphere the whole project has run so far, and it is
# the FIRST rung deliberately: the shape campaign must reproduce the existing sphere
# results at its own alpha = 1 before anything at alpha > 1 is believed.
ALPHA_LADDER = (1.0, 1.25, 1.5, 1.75, 2.0)


def g_of_u(u: float) -> float:
    """Union volume of the two-sphere clump, in units of the sub-sphere diameter cubed.

    V = ds^3 * g(u), u = s/ds. Closed form; see the module docstring.
    """
    if u < 0:
        raise ValueError(f"u must be >= 0, got {u}")
    if u >= 1.0:
        return math.pi / 3.0
    return math.pi / 3.0 - (math.pi / 12.0) * (2.0 + u) * (1.0 - u) ** 2


def alpha_to_u(alpha: float) -> float:
    """u = alpha - 1, with the admissible range asserted rather than clipped."""
    if not (1.0 <= alpha <= 2.0):
        raise ValueError(
            f"alpha must be in [1, 2] for a two-sphere clump, got {alpha}. "
            "Beyond 2 the spheres separate and the body is no longer connected; a "
            "longer rod needs more sub-spheres, which is a different family."
        )
    return alpha - 1.0


def sub_sphere_diameter(d_eq: float, alpha: float) -> float:
    """The sub-sphere diameter that makes a clump of aspect `alpha` have volume-equivalent
    diameter exactly `d_eq`. This is what keeps the PSD invariant along the shape axis."""
    u = alpha_to_u(alpha)
    return float(d_eq * ((math.pi / 6.0) / g_of_u(u)) ** (1.0 / 3.0))


def clump_volume(d_eq: float, alpha: float) -> float:
    """The union volume of the constructed clump. Must equal (pi/6)*d_eq^3 by design;
    it is computed the long way round so a test can catch a construction that drifts."""
    ds = sub_sphere_diameter(d_eq, alpha)
    return float(ds ** 3 * g_of_u(alpha_to_u(alpha)))


def corrugation(d_eq: float, alpha: float) -> float:
    """Fractional depth of the waist, the INTERLOCKING half of this knob.

    0 at alpha = 1 (a sphere), 2^(-1/3) = 0.7937 at alpha = 2, where the waist radius is
    exactly zero. Reported on every cell so that no result from the clump axis can be
    attributed to elongation without it.
    """
    u = alpha_to_u(alpha)
    ds = sub_sphere_diameter(d_eq, alpha)
    return float((ds / d_eq) * (1.0 - math.sqrt(max(0.0, 1.0 - u * u))))


def clump_centres(d_eq: float, alpha: float, axis) -> np.ndarray:
    """The two sub-sphere centres, relative to the clump centroid, along `axis`.

    `axis` is a 3-vector; it is normalised here. Returned as a (2, 3) array.
    """
    ds = sub_sphere_diameter(d_eq, alpha)
    s = alpha_to_u(alpha) * ds
    a = np.asarray(axis, dtype=float)
    nrm = float(np.linalg.norm(a))
    if nrm == 0.0:
        raise ValueError("axis must be non-zero")
    a = a / nrm
    return np.vstack([+0.5 * s * a, -0.5 * s * a])


def random_axes(n: int, seed: int) -> np.ndarray:
    """`n` orientations drawn UNIFORMLY on the sphere, as an (n, 3) array.

    Uniform on the sphere, not uniform in (theta, phi): the second clusters at the poles
    and would put a preferred direction into a campaign whose primary observable is
    contact-orientation anisotropy. Gaussian-normalise, which is exactly uniform.
    """
    rng = np.random.default_rng(seed)
    v = rng.normal(size=(int(n), 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def span_factor(alpha: float) -> float:
    """Longest extent of the clump divided by d_eq.

    The periodic collider needs cell edge > 2.1 * max body span, and an elongated grain
    is LONGER than the equal-volume sphere, so the N that clears the cell at alpha = 1
    does not clear it at alpha = 2. This is the factor `packing`'s clearance check must
    be multiplied by, and it is why the shape campaign cannot simply inherit N = 1500.
    """
    ds = sub_sphere_diameter(1.0, alpha)
    return float(ds * alpha)
