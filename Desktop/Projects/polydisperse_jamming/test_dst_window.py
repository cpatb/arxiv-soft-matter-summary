"""test_dst_window.py — the statistical choice `dst_window.py` rests on.

The window `W = phi_0 - phi_mu` is a DIFFERENCE between two ladders that were run on the
same seeds. Whether the two arms are PAIRED decides the error bar and therefore the
verdict, so the pairing is asserted here rather than assumed: if a future campaign adds
frictionless cells at seeds the frictional arm does not have, or at a different realised
CV, these tests fail loudly instead of the analysis silently reverting to an unpaired
comparison with a different answer.

`out/phase2_fits.csv` is committed, so no `pytest.skip` guard: if the input is missing
that is a failure, not a reason not to test.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

import dst_window as W

FITS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "out", "phase2_fits.csv")


def test_the_input_is_present():
    assert os.path.exists(FITS), "out/phase2_fits.csv is committed and must be here"


def test_the_two_frictions_share_seeds_AND_realised_CV():
    """The premise of pairing. cv_real matching is the part that could silently break."""
    d = pd.read_csv(FITS)
    d = d[d.n_requested == W.N_MATCH]
    for fam in ("lognormal", "uniform"):
        a = d[(d.family == fam) & (d.mu_p == W.MU_0)][["cv_nom", "seed", "cv_real"]]
        b = d[(d.family == fam) & (d.mu_p == W.MU_F)][["cv_nom", "seed", "cv_real"]]
        m = a.merge(b, on=["cv_nom", "seed"], suffixes=("_0", "_f"))
        assert len(m) > 0, f"{fam}: no paired cells at all"
        worst = float((m.cv_real_0 - m.cv_real_f).abs().max())
        # Machine epsilon, not zero: the same PSD realisation is reduced by two
        # different code paths, so the bar is 1e-12 rather than exact equality. A real
        # unpairing would show up at 1e-3 or above (the CV levels are 0.1 apart).
        assert worst < 1e-12, (
            f"{fam}: the two frictions no longer share PSD realisations "
            f"(worst cv_real mismatch {worst}); the paired difference is no longer valid")
        # ...and assert it is not TRIVIALLY zero because the merge found nothing.
        assert len(m) >= 15, f"{fam}: only {len(m)} paired cells"


def test_pairing_restricts_to_seeds_present_in_BOTH_arms():
    """A seed one arm lacks must not be averaged into the other arm's mean.

    ⚠ REWRITTEN 2026-09-09. The original pinned the discriminating case that existed
    then -- lognormal CV = 0.5 had ONE frictionless seed, so the paired phi_mu had to be
    seed 3's value alone (0.632565) rather than the three-seed mean (0.636854), a 0.0043
    difference. The `capfix` run repaired that arm to three seeds and the case vanished,
    so the test failed for the right reason and is restated as an INVARIANT that holds
    however the campaign grows: `n_pairs` is the size of the seed INTERSECTION, and both
    reported means are taken over exactly those seeds.
    """
    d = W.load_fits()
    t = W.cell_table(d, "lognormal")
    w = W.window_table(d, "lognormal", t)
    sub = d[d.family == "lognormal"]
    for _, r in w.iterrows():
        s0 = set(sub[(sub.mu_p == W.MU_0) & (sub.cv_nom == r.cv_nom)].seed)
        sf = set(sub[(sub.mu_p == W.MU_F) & (sub.cv_nom == r.cv_nom)].seed)
        shared = s0 & sf
        assert r.n_pairs == len(shared), (
            f"cv {r.cv_nom}: n_pairs {r.n_pairs} but the seed intersection is "
            f"{sorted(shared)}")
        for arm, col, seeds in ((W.MU_0, "phi_0", s0), (W.MU_F, "phi_mu", sf)):
            over_shared = sub[(sub.mu_p == arm) & (sub.cv_nom == r.cv_nom)
                              & sub.seed.isin(shared)].phi_J_lin.mean()
            assert np.isclose(getattr(r, col), over_shared, atol=1e-9), (
                f"cv {r.cv_nom} {col}: reported {getattr(r, col):.6f} is not the mean "
                f"over the shared seeds {sorted(shared)} ({over_shared:.6f}); a seed "
                f"outside the intersection leaked in")
        if len(shared) < max(len(s0), len(sf)):
            # the case the original test pinned: assert it is still handled
            wider = sub[(sub.mu_p == (W.MU_0 if len(s0) > len(sf) else W.MU_F))
                        & (sub.cv_nom == r.cv_nom)].phi_J_lin.mean()
            col = "phi_0" if len(s0) > len(sf) else "phi_mu"
            assert not np.isclose(getattr(r, col), wider, atol=1e-9), (
                f"cv {r.cv_nom}: the arm means coincide, so this level cannot "
                f"discriminate paired from unpaired")


def test_the_paired_error_bar_is_SMALLER_than_the_unpaired_one():
    """The discriminating half: a pairing that does not reduce the error is not pairing.

    If these came out equal the arms would be independent and the paired statistic would
    be pointless -- so this asserts the gap exists, not merely that both are finite.
    """
    d = W.load_fits()
    t = W.cell_table(d, "lognormal")
    w = W.window_table(d, "lognormal", t)
    paired = w[np.isfinite(w.W_se)]
    a = t[t.mu_p == W.MU_0].set_index("cv_nom")
    b = t[t.mu_p == W.MU_F].set_index("cv_nom")
    better = 0
    for _, r in paired.iterrows():
        unpaired = np.sqrt(a.loc[r.cv_nom, "se_mean"] ** 2
                           + b.loc[r.cv_nom, "se_mean"] ** 2)
        if np.isfinite(unpaired) and r.W_se < unpaired:
            better += 1
    assert better >= 3, (
        f"the paired se beat the unpaired one on only {better} levels; pairing is "
        f"supposed to cancel the shared PSD realisation")


def test_the_window_is_in_the_range_the_wyart_cates_reading_requires():
    """An external anchor: the literature window for spheres is ~0.58 to ~0.64.

    This is a sanity bar on the IDENTIFICATION, not on the slope. A window of 0.006 or
    0.6 would mean phi_0/phi_mu are not the quantities the docstring claims.
    """
    d = W.load_fits()
    t = W.cell_table(d, "lognormal")
    w = W.window_table(d, "lognormal", t)
    w0 = float(w[w.cv_nom == 0.0].W.iloc[0])
    assert 0.03 < w0 < 0.12, f"W(CV=0) = {w0:.4f} is not a plausible DST window"
    phi0 = float(w[w.cv_nom == 0.0].phi_0.iloc[0])
    assert 0.62 < phi0 < 0.65, (
        f"phi_0 = {phi0:.4f}; frictionless monodisperse RCP is ~0.64 and this is the "
        f"leg the docstring calls the SAFE identification")


def test_thin_ladders_are_flagged_not_dropped():
    """The THIN flag must track `n_ladder_min`, not a belief about which cells are thin.

    ⚠ REWRITTEN 2026-09-09. The original asserted every frictionless cell is THIN, which
    was true until `capfix` repaired CV 0.0 and 0.3 to four-point ladders. Pinning the
    state rather than the RULE is what made it fail; the rule is what the report depends
    on, so that is what is asserted now -- in both directions, so a flag that is always
    True or always False fails.
    """
    d = W.load_fits()
    t = W.cell_table(d, "lognormal")
    for _, r in t.iterrows():
        assert bool(r.thin) == (r.n_ladder_min <= W.N_LADDER_SELF_DIAGNOSING), (
            f"cv {r.cv_nom} mu {r.mu_p}: thin={r.thin} but n_ladder_min="
            f"{r.n_ladder_min} against a bar of {W.N_LADDER_SELF_DIAGNOSING}")
    assert t.thin.any() and not t.thin.all(), (
        "the THIN flag is constant across the table, so it is not discriminating "
        "anything")


def test_a_slope_needs_three_points():
    assert not np.isfinite(W.wls_slope([0.0, 1.0], [0.0, 1.0])[0])
    s, se, n = W.wls_slope([0.0, 1.0, 2.0], [0.0, 1.0, 2.0])
    assert np.isclose(s, 1.0) and n == 3
