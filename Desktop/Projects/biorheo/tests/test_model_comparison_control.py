"""The negative control for the "Cross and Carreau-Yasuda fit these curves 5x better than
Casson" claim (NEXT_STEPS item 5, the one referee item left undone).

The claim is in `data/horner2020/model_comparison.md` and is repeated in the README. The
obvious objection is parameter count: Cross has 4 and Carreau-Yasuda 5, against Casson's
2, so of course they fit 16 points better. The finding-1 referee pass (2026-09-05) ran a
control that answers it and the control was never made a test. It is one now.

THE CONTROL: generate data from a TRUE Casson curve on the real Horner rate grid with the
same multiplicative noise, then fit all three families. If the 5x gap on real blood were
just parameter count, the richer models would beat Casson on true-Casson data too. They
must not: on data a 2-parameter model generates, a 4- or 5-parameter model should fit no
better, and in a log-viscosity least squares on 16 points it should fit WORSE, because the
extra parameters buy nothing and the optimiser has more room to land off the truth.

PRE-REGISTERED GATE (written before the run, 2026-09-05; not loosened afterwards):

  N1  On >= 18 of 20 synthetic true-Casson trials, Casson's rmse_log is <= that of BOTH
      Cross and Carreau-Yasuda. This is the direction that matters: if the richer models
      won here, the real-data gap would be uninterpretable.
  N2  The median ratio rmse(richer) / rmse(Casson) on true-Casson data is >= 1.0, i.e.
      the richer families are NOT better on data whose truth is Casson.
  N3  The same three families on the REAL donor curves give the opposite ordering, with
      Cross and Carreau-Yasuda beating Casson by a factor recomputed here from the
      shipped records rather than read out of the .md.

N1/N2 use synthetic data only and are fast. N3 refits the real records and is marked
slow.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from biorheo.constitutive.gnf import CARREAU_YASUDA, CASSON, CROSS
from biorheo.identify.fit import fit_rate_curve
from biorheo.identify.synthetic import synthetic_rate_curve
from biorheo.registry.loader import load_all

ROOT = Path(__file__).resolve().parent.parent

N_TRIALS = 20
NOISE_CV = 0.01          # the referee's control used 1 %
MIN_CASSON_WINS = 18     # N1
FAMILIES = (CROSS, CARREAU_YASUDA)


@pytest.fixture(scope="module")
def blood():
    return {k: v for k, v in sorted(load_all().items()) if v.model.name == "casson"}


@pytest.fixture(scope="module")
def grid(blood):
    """The real Horner rate grid: 16 rates in [0.1, 700] 1/s. Using the actual grid
    matters -- the gap is a statement about a fit to 16 points at these rates, and a
    denser or wider grid would be a different question."""
    rec = blood["blood_horner2020_donor_A"]
    lo, hi = rec.envelope.gammadot_min, rec.envelope.gammadot_max
    return np.geomspace(lo, hi, 16)


def _rmse_log(model, gd, eta, seed):
    ident = fit_rate_curve(model, gd, eta, seed=seed)
    return float(ident.fit.rmse_log)


def test_richer_models_do_not_beat_casson_on_true_casson_data(blood, grid):
    """N1 and N2. If they did, the 5x gap on real blood would be parameter count."""
    truth = np.array([float(v) for v in blood["blood_horner2020_donor_A"].params_si])
    casson_wins, ratios = 0, {m.name: [] for m in FAMILIES}
    for seed in range(N_TRIALS):
        eta = synthetic_rate_curve(CASSON, truth, grid, NOISE_CV, seed)
        r_casson = _rmse_log(CASSON, grid, eta, seed)
        beaten = False
        for m in FAMILIES:
            r = _rmse_log(m, grid, eta, seed)
            ratios[m.name].append(r / r_casson)
            beaten = beaten or r < r_casson
        casson_wins += not beaten
    assert casson_wins >= MIN_CASSON_WINS, (casson_wins, MIN_CASSON_WINS)
    for name, rs in ratios.items():
        assert float(np.median(rs)) >= 1.0, (name, float(np.median(rs)))


@pytest.mark.slow
def test_the_real_curves_give_the_opposite_ordering(blood):
    """N3: the same three families on real donor curves, recomputed from the shipped
    records rather than read out of `model_comparison.md`. Casson must LOSE here, and by
    a wide margin -- that contrast with the synthetic control is the whole argument."""
    src = ROOT / "data" / "horner2020" / "steady_shear_human.csv"
    if not src.exists():
        pytest.skip("the source CSV is gitignored; regenerate with records_from_horner.py")
    # wide format: one rate column, one stress column per donor
    import csv
    with open(src, newline="") as fh:
        rows = list(csv.DictReader(fh))
    gd_all = np.array([float(r["shear_rate_1_per_s"]) for r in rows])
    wins = 0
    ratios = []
    for rec in list(blood.values())[:6]:                 # six donors is enough to order
        donor = rec.id.split("_")[-1]
        col = f"stress_Pa_{donor}"
        assert col in rows[0], (col, sorted(rows[0])[:3])
        sigma = np.array([float(r[col]) for r in rows])
        gd, eta = gd_all, sigma / gd_all          # the layer fits eta = sigma/gd
        if gd.size < 10:
            continue
        r_casson = _rmse_log(CASSON, gd, eta, 0)
        best = min(_rmse_log(m, gd, eta, 0) for m in FAMILIES)
        ratios.append(r_casson / best)
        wins += best < r_casson
    assert ratios, "no donor curves read"
    assert wins == len(ratios), (wins, len(ratios))
    assert float(np.median(ratios)) > 2.0, float(np.median(ratios))
