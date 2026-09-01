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
