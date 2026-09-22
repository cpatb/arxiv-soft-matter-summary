"""test_multiplet.py — anchors for the multiplet (aggregated-suspension) PSD family.

Engine-independent: runs on Windows without Yade.

Same rule as `test_psd.py`: every anchor is something that does NOT pass through the
code under test — a value computed by hand from the definition, an exactly-known
special case, an independent particle-array reduction, or a discriminating negative.

The one place a round trip IS asserted is `test_volume_fraction_round_trip`, and it is
asserted for a stated reason rather than as evidence: `MultipletPSD` hands VOLUME
fractions to `DiscretePSD`, which inverts them back to NUMBER fractions with different
float arithmetic, and the campaign's control variable is the number composition. The
test pins how far that inversion may drift. It is a tolerance, not a result.
"""

import math

import numpy as np
import pytest

from psd import (
    DiscretePSD,
    LogNormalPSD,
    MultipletPSD,
    make_psd,
    parse_multiplet_family,
)

D10 = 0.8e-6


# --------------------------------------------------------------------------
# A. Hand-computed anchors
# --------------------------------------------------------------------------


def test_hand_computed_half_singlets_half_doublets_compact():
    """d_f = 3, n = (1/2, 1/2). Every number below is arithmetic done by hand.

        gamma   = 3/d_f = 1
        d_rel   = (1, 2**(1/3))            = (1, 1.2599210498948732)
        <d>/d1  = (1 + 2**(1/3))/2         = 1.1299605249474366
        <d^2>/d1^2 = (1 + 2**(2/3))/2      = 1.2937005259840998
        CV^2    = <d^2>/<d>^2 - 1          = 0.013228897...
        Lambda  = <k^1>/<k^1>              = 1 exactly
    """
    p = MultipletPSD(d10=D10, d_f=3.0, n_k=(0.5, 0.5))
    c = 2.0 ** (1.0 / 3.0)
    mean = (1.0 + c) / 2.0
    msq = (1.0 + c * c) / 2.0
    cv_hand = math.sqrt(msq / mean**2 - 1.0)

    assert p.gamma == pytest.approx(1.0, abs=0.0)
    assert p.cv == pytest.approx(cv_hand, rel=1e-12)
    # The value spelled out, from 40-digit `decimal` arithmetic on the same expression
    # (not from this module): CV = 0.1150133319511157933436699786706133934539.
    # ⚠ My first transcription of this literal was WRONG (0.11501669) and this assertion
    # is what caught it — which is the reason a spelled-out anchor is kept beside the
    # symbolic one rather than instead of it.
    assert p.cv == pytest.approx(0.11501333195111579, rel=1e-13)
    assert p.lambda_trap == pytest.approx(1.0, abs=1e-15)
    assert p.k_bar == pytest.approx(1.5, rel=1e-15)
    # the monomer diameter inverts d10 = d1 * <k^(gamma/3)>
    assert p.d1 == pytest.approx(D10 / mean, rel=1e-14)
    assert p.size_ratio == pytest.approx(c, rel=1e-12)


def test_lambda_is_exactly_one_at_df3_for_any_composition():
    """A COMPACT multiplet excludes exactly its own solid volume, so it traps no fluid.

    This must hold identically, not approximately, and for every composition — it is
    the statement `k**1 == k`. A family that only satisfied it on average would mean
    the volume bookkeeping was wrong somewhere.
    """
    rng = np.random.default_rng(20260922)
    for _ in range(200):
        n = rng.dirichlet(np.ones(rng.integers(1, 7)))
        p = MultipletPSD(d10=D10, d_f=3.0, n_k=tuple(n))
        assert p.lambda_trap == pytest.approx(1.0, abs=1e-14)


def test_pure_k_mer_has_zero_cv_and_lambda_k_pow_gamma_minus_one():
    """All particles are K-mers: a MONODISPERSE pack of big effective spheres.

    Lambda = K**gamma / K = K**(gamma-1), by hand. This is the discriminating case for
    the two factors being divided in the right order: the inverse would give K**(1-gamma).
    """
    for K, d_f in ((4, 1.8), (7, 2.1), (3, 1.0)):
        gamma = 3.0 / d_f
        p = MultipletPSD(d10=D10, d_f=d_f, n_k=tuple([0.0] * (K - 1) + [1.0]))
        assert p.cv == pytest.approx(0.0, abs=1e-14)
        assert p.k_bar == pytest.approx(float(K), rel=1e-14)
        assert p.lambda_trap == pytest.approx(K ** (gamma - 1.0), rel=1e-12)
        assert p.lambda_trap > 1.0


def test_unaggregated_suspension_is_the_monodisperse_pack():
    """n_k = (1,) is a suspension that has not aggregated: it must be EXACTLY the
    monodisperse pack `LogNormalPSD(cv=0)` produces, at any d_f, because there is
    nothing for d_f to act on."""
    for d_f in (0.9, 1.8, 2.5, 3.0):
        p = MultipletPSD(d10=D10, d_f=d_f, n_k=(1.0,))
        assert p.cv == pytest.approx(0.0, abs=1e-15)
        assert p.lambda_trap == pytest.approx(1.0, abs=1e-15)
        assert p.d1 == pytest.approx(D10, rel=1e-15)
        ref = LogNormalPSD(d10=D10, cv=0.0)
        rng_a, rng_b = np.random.default_rng(3), np.random.default_rng(3)
        np.testing.assert_allclose(p.sample(64, rng_a), ref.sample(64, rng_b), rtol=1e-15)


def test_chain_limit_df1_gives_dk_equals_k_times_d1():
    """d_f = 1 is a straight chain of touching monomers: d_k = k*d1 exactly."""
    p = MultipletPSD(d10=D10, d_f=1.0, n_k=(0.6, 0.3, 0.1))
    assert p.gamma == pytest.approx(3.0, abs=0.0)
    np.testing.assert_allclose(np.asarray(p.d_rel), [1.0, 2.0, 3.0], rtol=1e-14)
    # <k> = 0.6 + 0.6 + 0.3 = 1.5, and d10 = d1 * <k> here
    assert p.k_bar == pytest.approx(1.5, rel=1e-14)
    assert p.d1 == pytest.approx(D10 / 1.5, rel=1e-14)


# --------------------------------------------------------------------------
# B. Independent-route anchors (no moment algebra)
# --------------------------------------------------------------------------


def test_moments_match_an_explicit_particle_array():
    """Rebuild CV, d32/d10, d43/d10 and Lambda ONE PARTICLE AT A TIME.

    This route never forms a moment sum and never touches `k_moment`: it draws an
    aggregation number per particle, converts that particle to a diameter and to a
    solid volume, and reduces with plain means. It is the anchor the class's own
    algebra cannot fake.
    """
    n_k = (0.55, 0.25, 0.13, 0.07)
    for d_f in (1.8, 2.1, 3.0):
        p = MultipletPSD(d10=D10, d_f=d_f, n_k=n_k)
        gamma = 3.0 / d_f
        rng = np.random.default_rng(11)
        kk = rng.choice(np.arange(1.0, 5.0), size=2_000_000, p=np.asarray(n_k))
        d = kk ** (gamma / 3.0)
        v_eff = (np.pi / 6.0) * d**3
        v_sol = kk * (np.pi / 6.0)

        assert p.cv == pytest.approx(d.std(ddof=0) / d.mean(), rel=3e-3)
        assert p.d32 / p.moment(1) == pytest.approx(
            (d**3).mean() / (d**2).mean() / d.mean(), rel=3e-3)
        assert p.d43 / p.moment(1) == pytest.approx(
            (d**4).mean() / (d**3).mean() / d.mean(), rel=3e-3)
        assert p.lambda_trap == pytest.approx(v_eff.sum() / v_sol.sum(), rel=3e-3)


def test_lambda_never_below_one():
    """k**gamma >= k pointwise for gamma >= 1 and k >= 1, so Lambda >= 1 ALWAYS.

    A violation would mean aggregation had been made to REDUCE the excluded volume,
    which no admissible d_f in (0, 3] can do.
    """
    rng = np.random.default_rng(5)
    for _ in range(2000):
        K = int(rng.integers(1, 9))
        n = rng.dirichlet(np.ones(K))
        d_f = float(rng.uniform(0.5, 3.0))
        assert MultipletPSD(d10=D10, d_f=d_f, n_k=tuple(n)).lambda_trap >= 1.0 - 1e-12


def test_lambda_and_cv_both_increase_as_aggregates_open_up():
    """Lowering d_f at fixed composition opens the aggregates: both Lambda and CV rise.

    This is the DISCRIMINATING direction. The competing-lever argument in
    `derivations/multiplet_psd_2026-09-22.md` needs both to move the same way — if CV
    fell as Lambda rose there would be no competition to resolve.
    """
    n_k = (0.60, 0.25, 0.10, 0.05)
    dfs = [3.0, 2.5, 2.1, 1.8, 1.5]
    lam = [MultipletPSD(d10=D10, d_f=f, n_k=n_k).lambda_trap for f in dfs]
    cvs = [MultipletPSD(d10=D10, d_f=f, n_k=n_k).cv for f in dfs]
    assert all(b > a for a, b in zip(lam, lam[1:])), lam
    assert all(b > a for a, b in zip(cvs, cvs[1:])), cvs


# --------------------------------------------------------------------------
# C. It really is the discrete pack it claims to be
# --------------------------------------------------------------------------


def test_volume_fraction_round_trip():
    """TOLERANCE, NOT A RESULT (see the module docstring).

    `MultipletPSD` gives `DiscretePSD` volume fractions; `DiscretePSD` inverts them to
    number fractions with `y * (d_max/d)**3`, which is a different float expression.
    The composition is this family's control variable, so the drift is pinned.
    """
    for d_f in (1.8, 2.4, 3.0):
        n_k = (0.55, 0.25, 0.13, 0.07)
        p = MultipletPSD(d10=D10, d_f=d_f, n_k=n_k)
        np.testing.assert_allclose(p.number_fractions, np.asarray(n_k), rtol=1e-12)


def test_delegates_to_discrete_exactly():
    """A multiplet pack IS the discrete pack built from its own (d_rel, y).

    Same diameters element for element and the same sampled array from the same seed —
    not "the same to within scatter". If this ever fails, the family has grown a second
    sampling path and the anchors of `DiscretePSD` no longer cover it.
    """
    p = MultipletPSD(d10=D10, d_f=1.8, n_k=(0.5, 0.3, 0.2))
    q = DiscretePSD(d10=D10, d_rel=p.d_rel, y=p.y)
    np.testing.assert_array_equal(p.diameters, q.diameters)
    a = p.sample(200, np.random.default_rng(77))
    b = q.sample(200, np.random.default_rng(77))
    np.testing.assert_array_equal(a, b)


def test_class_counts_are_exact_not_binomial():
    n_k = (0.5, 0.3, 0.2)
    p = MultipletPSD(d10=D10, d_f=2.0, n_k=n_k)
    counts = p.class_counts(1000)
    assert counts.sum() == 1000
    assert counts.min() >= 1
    np.testing.assert_allclose(counts / 1000.0, np.asarray(n_k), atol=1e-3)
    # the seed changes POSITIONS only, never the composition
    for seed in (1, 2, 3):
        d = p.sample(1000, np.random.default_rng(seed))
        assert sorted(collections_counter(d).values()) == sorted(counts.tolist())


def collections_counter(a):
    import collections
    return collections.Counter(np.round(a, 18).tolist())


def test_zero_entries_are_dropped_so_a_gap_is_the_pack_it_is():
    """(0.7, 0, 0.3) is singlets and TRIPLETS — a two-class pack with no doublets."""
    p = MultipletPSD(d10=D10, d_f=2.0, n_k=(0.7, 0.0, 0.3))
    assert p.diameters.size == 2
    np.testing.assert_allclose(p.k, [1.0, 3.0])
    assert p.k_bar == pytest.approx(0.7 * 1 + 0.3 * 3, rel=1e-14)


# --------------------------------------------------------------------------
# D. Discriminating negatives
# --------------------------------------------------------------------------


def test_family_string_round_trips_through_make_psd():
    p = make_psd("multiplet:f1.8:n0.70/0.20/0.10", d10=D10, cv=0.0)
    assert isinstance(p, MultipletPSD)
    assert p.d_f == pytest.approx(1.8)
    assert p.n_k == (0.7, 0.2, 0.1)
    assert parse_multiplet_family("multiplet:f1.8:n0.70/0.20/0.10") == (1.8, (0.7, 0.2, 0.1))
    assert parse_multiplet_family("lognormal") is None


def test_malformed_spec_raises_rather_than_defaulting():
    """A silent default would put a pack of the wrong aggregation state into a campaign."""
    for bad in ("multiplet:1.8:n0.7/0.3",        # missing the f
                "multiplet:f1.8:0.7/0.3",        # missing the n
                "multiplet:f1.8",                # too few fields
                "multiplet:fX:n0.7/0.3",         # unparsable d_f
                "multiplet:f1.8:n0.7/oops"):     # unparsable fraction
        with pytest.raises(ValueError):
            parse_multiplet_family(bad)


def test_cv_cannot_be_imposed_from_outside():
    """The whole point of the family is that CV is a READOUT of the aggregation state."""
    with pytest.raises(ValueError, match="derives its CV"):
        make_psd("multiplet:f1.8:n0.7/0.3", d10=D10, cv=0.3)


def test_rejects_inadmissible_fractal_dimension():
    """d_f > 3 would mean a multiplet excluding LESS volume than its own solid."""
    for bad in (3.5, 0.0, -1.0):
        with pytest.raises(ValueError):
            MultipletPSD(d10=D10, d_f=bad, n_k=(0.5, 0.5))
    for bad_n in ((), (0.0, 0.0), (-0.1, 1.0)):
        with pytest.raises(ValueError):
            MultipletPSD(d10=D10, d_f=2.0, n_k=bad_n)


def test_multiplet_cv_is_not_free_to_reach_arbitrary_values():
    """A DISCRIMINATING NEGATIVE, and it is the family's physical content.

    With only singlets and doublets there is a HARD CEILING on the width, because the
    size ratio is fixed at 2**(1/d_f) and only the composition is free. At d_f = 3 that
    ratio is 1.26, and no composition of two classes at ratio 1.26 can reach the
    CV = 0.5 the continuous families are routinely run at. So "raise CV" and "aggregate
    more" are NOT the same experiment, which is the reason this family exists.
    """
    best = max(MultipletPSD(d10=D10, d_f=3.0, n_k=(1 - x, x)).cv
               for x in np.linspace(0.001, 0.999, 999))
    assert best < 0.14, best
    # opening the aggregates up raises the ceiling, but not without limit either
    best_open = max(MultipletPSD(d10=D10, d_f=1.8, n_k=(1 - x, x)).cv
                    for x in np.linspace(0.001, 0.999, 999))
    assert 0.14 < best_open < 0.30, best_open


# --------------------------------------------------------------------------
# E. The family string has to survive becoming a filename
# --------------------------------------------------------------------------


def test_fs_safe_does_not_move_existing_names():
    """DISCRIMINATING: every family ever run through phase2_phiJ must be UNCHANGED.

    If sanitising moved one of these, the next run would write to a new file and
    silently stop resuming the campaign already on disk.
    """
    from phase2_phiJ import fs_safe
    for fam in ("lognormal", "uniform", "bidisperse"):
        assert fs_safe(fam) == fam


def test_fs_safe_removes_every_path_hostile_character():
    from phase2_phiJ import fs_safe
    out = fs_safe("multiplet:f1.8:n0.70/0.20/0.10")
    assert out == "multiplet-f1.8-n0.70_0.20_0.10"
    assert ":" not in out and "/" not in out and "\\" not in out


def test_fs_safe_is_injective_on_the_families_in_play():
    """Two DIFFERENT aggregation states must not collapse onto one filename.

    A collision would pool two campaigns into one shard file, which is the defect
    this repo already records for dedup keys.
    """
    from phase2_phiJ import fs_safe
    fams = [
        "multiplet:f1.8:n0.70/0.20/0.10",
        "multiplet:f1.8:n0.70/0.20/0.11",
        "multiplet:f2.1:n0.70/0.20/0.10",
        "multiplet:f3.0:n0.50/0.50",
        "binary:r3.0:x0.30",
        "lognormal", "uniform", "bidisperse",
    ]
    assert len(set(map(fs_safe, fams))) == len(fams)


def test_the_cell_key_is_NOT_sanitised():
    """The key must keep the raw family string: it is what resume matches on."""
    from phase2_phiJ import cell_key
    fam = "multiplet:f1.8:n0.70/0.20/0.10"
    assert fam in cell_key(0.0, 0.3, 1, 3e-5, 1500, family=fam)
