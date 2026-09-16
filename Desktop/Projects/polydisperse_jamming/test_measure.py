"""test_measure.py — anchors for the structural-observable layer.

Engine-independent: pure functions of a contact graph, so these run on Windows.

Anchors are hand-constructed graphs with a known answer, closed-form values for the
fabric tensor, and an analytically-known distribution for the force tail. No test here
inverts the function it validates.
"""

import numpy as np
import pytest

from measure import (
    backbone,
    fabric_tensor,
    force_moments,
    iterated_rattlers,
    rattler_fraction,
    spanning_cluster,
    z_by_size,
)


def graph_from_edges(edges, isolated=()):
    """Build (counts, neigh) from an undirected edge list."""
    neigh: dict[int, set[int]] = {}
    for a, b in edges:
        neigh.setdefault(a, set()).add(b)
        neigh.setdefault(b, set()).add(a)
    for i in isolated:
        neigh.setdefault(i, set())
    counts = {k: len(v) for k, v in neigh.items()}
    return counts, neigh


# --------------------------------------------------------------------------
# A. Rattlers — hand-built graphs where the ITERATION is what is being tested
# --------------------------------------------------------------------------


def test_isolated_grain_is_a_rattler():
    counts, neigh = graph_from_edges([], isolated=[0])
    assert iterated_rattlers(counts, neigh, min_contacts=4) == {0}
    assert rattler_fraction(counts, neigh, min_contacts=4) == 1.0


def test_fully_connected_six_is_not_rattlers():
    """K6: every grain has 5 contacts >= 4, so nothing is a rattler."""
    edges = [(i, j) for i in range(6) for j in range(i + 1, 6)]
    counts, neigh = graph_from_edges(edges)
    assert iterated_rattlers(counts, neigh, min_contacts=4) == set()
    assert backbone(counts, neigh, min_contacts=4) == set(range(6))


def test_iteration_is_required_not_one_pass():
    """THE test for `iterated_rattlers`, hand-constructed so one pass gives the
    WRONG answer.

    K6 on {0..5}, plus grain 6 attached to four members of the K6, plus grain 7
    attached ONLY to grain 6.

      pass 1: grain 7 has 1 contact  -> rattler. Grain 6 has 5 -> survives.
      after removing 7: grain 6 has 4 -> still survives (4 >= 4).

    Now attach grain 6 to only THREE K6 members instead of four:
      pass 1: grain 7 (1 contact) -> rattler ; grain 6 has 4 -> a ONE-PASS
              algorithm keeps it.
      iterated: removing 7 drops grain 6 to 3 -> grain 6 is a rattler too.

    So the correct answer is {6, 7} and a single pass returns {7}.
    """
    edges = [(i, j) for i in range(6) for j in range(i + 1, 6)]
    edges += [(6, 0), (6, 1), (6, 2), (7, 6)]
    counts, neigh = graph_from_edges(edges)

    # Threshold pinned at 4 explicitly: this test is about the ITERATION, and the
    # module default is now 2 (see measure.RATTLER_THRESHOLDS).
    one_pass = {b for b, c in counts.items() if c < 4}
    assert one_pass == {7}, f"setup wrong: one-pass gave {one_pass}"

    assert iterated_rattlers(counts, neigh, min_contacts=4) == {6, 7}
    assert backbone(counts, neigh, min_contacts=4) == set(range(6))


def test_whole_chain_unravels():
    """A pure chain 0-1-2-...-9: every grain has <= 2 contacts, so ALL are rattlers.
    A discriminating negative -- a chain is connected but carries no load.
    """
    edges = [(i, i + 1) for i in range(9)]
    counts, neigh = graph_from_edges(edges)
    assert iterated_rattlers(counts, neigh, min_contacts=4) == set(range(10))
    assert rattler_fraction(counts, neigh, min_contacts=4) == 1.0
    assert backbone(counts, neigh, min_contacts=4) == set()


def test_min_contacts_threshold_is_honoured():
    """A RING of 10, where every grain has exactly 2 contacts. At min_contacts=2
    nothing is a rattler; at min_contacts=3 the whole ring unravels. Same graph,
    opposite answers -- so the parameter is demonstrably doing something.

    (A chain was tried first and does not discriminate: at min_contacts=2 its free end
    goes, which drops its neighbour to 1, and the cascade removes everything -- the
    same answer as min_contacts=4.)
    """
    ring = [(i, (i + 1) % 10) for i in range(10)]
    counts, neigh = graph_from_edges(ring)
    assert all(c == 2 for c in counts.values())
    assert iterated_rattlers(counts, neigh, min_contacts=2) == set()
    assert iterated_rattlers(counts, neigh, min_contacts=3) == set(range(10))


# --------------------------------------------------------------------------
# B. Connectivity
# --------------------------------------------------------------------------


def test_two_disjoint_blocks_are_two_components():
    e1 = [(i, j) for i in range(6) for j in range(i + 1, 6)]
    e2 = [(i, j) for i in range(10, 16) for j in range(i + 1, 16)]
    counts, neigh = graph_from_edges(e1 + e2)
    bb = backbone(counts, neigh, min_contacts=4)
    sc = spanning_cluster(neigh, bb)
    assert sc["n_components"] == 2
    assert sc["largest_frac"] == pytest.approx(0.5)


def test_single_block_is_one_component():
    edges = [(i, j) for i in range(8) for j in range(i + 1, 8)]
    counts, neigh = graph_from_edges(edges)
    sc = spanning_cluster(neigh, backbone(counts, neigh, min_contacts=4))
    assert sc["n_components"] == 1
    assert sc["largest_frac"] == pytest.approx(1.0)


def test_empty_backbone_gives_zero_not_a_crash():
    edges = [(i, i + 1) for i in range(9)]
    counts, neigh = graph_from_edges(edges)
    sc = spanning_cluster(neigh, backbone(counts, neigh, min_contacts=4))
    assert sc["largest_frac"] == 0.0
    assert sc["n_components"] == 0


# --------------------------------------------------------------------------
# C. Fabric tensor — closed-form values computed by hand
# --------------------------------------------------------------------------


def test_fully_aligned_normals_give_anisotropy_exactly_one():
    """All normals along z:

        F = diag(0, 0, 1)
        dev = F - I/3 = diag(-1/3, -1/3, 2/3)
        aniso = sqrt(1.5 * (1/9 + 1/9 + 4/9)) = sqrt(1.5 * 6/9) = sqrt(1) = 1
    """
    n = np.tile(np.array([0.0, 0.0, 1.0]), (500, 1))
    assert fabric_tensor(n)["aniso"] == pytest.approx(1.0, rel=1e-12)


def test_isotropic_normals_give_near_zero_anisotropy():
    """Uniform points on the sphere -> F = I/3 -> aniso = 0, up to sqrt(N) noise."""
    rng = np.random.default_rng(4)
    n = rng.normal(size=(200_000, 3))
    a = fabric_tensor(n)["aniso"]
    assert a < 0.01, f"isotropic draw gave aniso = {a}"


def test_three_orthogonal_axes_are_isotropic_exactly():
    """Equal numbers along x, y, z give F = I/3 exactly, so aniso = 0 to machine
    precision. A discriminating case: 'aligned' is not the same as 'anisotropic'."""
    n = np.repeat(np.eye(3), 100, axis=0)
    assert fabric_tensor(n)["aniso"] == pytest.approx(0.0, abs=1e-12)


def test_fabric_rejects_degenerate_input():
    assert np.isnan(fabric_tensor(np.zeros((3, 3)))["aniso"])


# --------------------------------------------------------------------------
# D. Force distribution — analytic anchor
# --------------------------------------------------------------------------


def test_exponential_forces_give_tail_slope_minus_one():
    """For f drawn from an exponential with mean <f>, x = f/<f> has pdf e^-x, so
    ln P(x) has slope exactly -1. This is the published DST/granular tail shape
    (anchor A6), and it is checked against a distribution whose answer is known in
    closed form rather than against another run of this code.
    """
    rng = np.random.default_rng(12)
    f = rng.exponential(scale=3.7, size=400_000)
    m = force_moments(f)
    assert m["tail_slope"] == pytest.approx(-1.0, abs=0.06), m
    assert m["cv"] == pytest.approx(1.0, rel=0.02)   # exponential has CV = 1 exactly


def test_narrow_forces_have_a_much_steeper_tail():
    """Discriminating negative: a nearly-monodisperse force set must NOT return the
    exponential slope. If it did, the statistic would be measuring the fitting window
    rather than the distribution."""
    rng = np.random.default_rng(13)
    f = rng.normal(loc=1.0, scale=0.02, size=200_000)
    m = force_moments(f[f > 0])
    assert np.isnan(m["tail_slope"]) or m["tail_slope"] < -5.0, m


def test_too_few_forces_returns_nan_not_a_number():
    m = force_moments(np.array([1.0, 2.0, 3.0]))
    assert np.isnan(m["mean"])
    assert m["n"] == 3


# --------------------------------------------------------------------------
# E. Size-resolved coordination
# --------------------------------------------------------------------------


def test_z_by_size_reports_n_per_bin():
    edges = [(i, j) for i in range(8) for j in range(i + 1, 8)]
    counts, neigh = graph_from_edges(edges)
    diam = {i: 1.0 + 0.1 * i for i in range(8)}
    out = z_by_size(counts, neigh, diam, n_bins=4)
    assert sum(out["n"]) == 8
    assert len(out["d_centre"]) == len(out["z_mean"]) == len(out["n"])


def test_z_by_size_detects_a_real_size_trend():
    """Big grains wired to everything, small grains wired to one big grain only.
    Z must rise with size -- a directional, discriminating assertion.
    """
    big = list(range(4))
    small = list(range(4, 20))
    edges = [(i, j) for i in big for j in big if i < j]
    edges += [(s, big[s % 4]) for s in small]
    counts, neigh = graph_from_edges(edges)
    diam = {**{i: 3.0 for i in big}, **{s: 1.0 for s in small}}
    out = z_by_size(counts, neigh, diam, n_bins=2)
    assert out["z_mean"][-1] > out["z_mean"][0], out


def test_z_by_size_empty_input():
    assert z_by_size({}, {}, {}, n_bins=3)["n"] == []


def test_com_vel_is_mass_weighted_not_number_averaged():
    """A centre-of-mass velocity must be sum(m v)/sum(m), not a mean of velocities.

    `phase4_impact._com_vel` was an unweighted mean until 2026-09-03. At CV = 0 the two
    agree exactly (equal masses), so every monodisperse control passed; at CV = 0.5 the
    masses span ~r^3 and the number-average is dominated by the smallest grains, making the
    error CV-DEPENDENT and putting it directly on the polydispersity contrast.

    Tested here on a stub rather than through Yade, so it runs without the engine: two
    bodies of very different mass moving oppositely. The mass-weighted answer follows the
    heavy one; the number-average returns ~0 and would report a rebound as a null.
    """
    import numpy as np

    class _S:
        def __init__(self, v, m):
            self.vel, self.mass = v, m

    class _B:
        def __init__(self, v, m):
            self.state = _S(v, m)

    class _O:
        def __init__(self, bodies):
            self.bodies = bodies

    # heavy grain moving +1, light grain moving -1; masses 27:1 (a 3x radius ratio)
    sim = _O({0: _B((0.0, 0.0, +1.0), 27.0), 1: _B((0.0, 0.0, -1.0), 1.0)})

    def com_vel(sim_, ids):
        v = np.array([[sim_.bodies[i].state.vel[j] for j in range(3)] for i in ids])
        m = np.array([sim_.bodies[i].state.mass for i in ids], float)
        return (v * m[:, None]).sum(axis=0) / m.sum()

    got = com_vel(sim, [0, 1])[2]
    assert got == pytest.approx((27.0 - 1.0) / 28.0, abs=1e-12)
    # the number-average would have said ~0 -- a rebound reported as a null
    naive = np.array([[sim.bodies[i].state.vel[j] for j in range(3)]
                      for i in [0, 1]]).mean(axis=0)[2]
    assert abs(naive) < 1e-12
    assert abs(got - naive) > 0.9
# --------------------------------------------------------------- backbone-resolved fabric
# Added 2026-09-18 with PREREG_FABRIC_2026-09-18.md. `packing` is engine-independent at
# import time (it imports yade only inside build_and_compress), so these run on Windows.

def test_backbone_mask_requires_BOTH_ends_in_the_backbone():
    """The discriminating case: a contact with one end rattling must be excluded.

    An `or` instead of an `and` passes every other test one would write.
    """
    from packing import backbone_mask
    ids = [[0, 1], [1, 2], [2, 3], [3, 0]]
    m = backbone_mask(ids, backbone={0, 1, 2})
    assert list(m) == [True, True, False, False]


def test_backbone_mask_is_aligned_with_its_normals_row_for_row():
    from packing import backbone_mask
    ids = [[5, 6], [7, 8], [6, 7]]
    m = backbone_mask(ids, backbone={5, 6, 7})
    assert m.shape == (3,)
    assert list(m) == [True, False, True]


def test_backbone_mask_on_an_empty_contact_list():
    from packing import backbone_mask
    m = backbone_mask([], backbone={1, 2})
    assert m.shape == (0,) and m.dtype == bool


def test_backbone_mask_with_an_empty_backbone_selects_nothing():
    from packing import backbone_mask
    m = backbone_mask([[0, 1], [1, 2]], backbone=set())
    assert not m.any()


def test_fabric_null_reproduces_the_closed_form_n_times_aniso_squared_is_one():
    """E[A^2] = 1/n EXACTLY for n independent uniform normals.

    The derivation is in `diag_fabric_null.py`; this is the anchor that does not pass
    through the campaign code, and it is what licenses T = n*A^2 as the observable that
    removes the contact-count confound. PREREG_FABRIC_2026-09-18.md sec. 0(a).
    """
    rng = np.random.default_rng(4242)
    for n in (400, 1600):
        v = rng.normal(size=(600, n, 3))
        v /= np.linalg.norm(v, axis=2, keepdims=True)
        t = np.array([n * fabric_tensor(v[r])["aniso"] ** 2 for r in range(v.shape[0])])
        assert abs(t.mean() - 1.0) < 0.08, (n, t.mean())
        # and the scatter is the chi2_5/5 one, sqrt(2/5) = 0.632 -- a routine returning a
        # constant would pass the mean test and fail this one.
        assert 0.45 < t.std(ddof=1) < 0.85, (n, t.std(ddof=1))


def test_the_isotropic_floor_falls_as_one_over_sqrt_n():
    """Quadrupling n must halve `aniso`. This is the confound the campaign divides out."""
    rng = np.random.default_rng(99)
    def floor(n):
        v = rng.normal(size=(300, n, 3))
        v /= np.linalg.norm(v, axis=2, keepdims=True)
        return float(np.mean([fabric_tensor(v[r])["aniso"] for r in range(v.shape[0])]))
    a1, a4 = floor(500), floor(2000)
    assert 1.8 < a1 / a4 < 2.2, (a1, a4, a1 / a4)
