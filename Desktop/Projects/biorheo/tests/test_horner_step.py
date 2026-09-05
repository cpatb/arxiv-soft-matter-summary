"""Data-level facts from the Horner 2020 step-down sheet, pinned so that a re-extraction
or a re-reading of the steady sheet cannot silently change them. Numbers from
scripts/step_shear_horner.py, 2026-09-05."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent.parent / "data" / "horner2020"


def _matrix(name: str, prefix: str):
    rows = list(csv.reader(open(DATA / name, newline="", encoding="utf-8")))
    x = np.array([r[0] for r in rows[1:]], dtype=float)
    cols = {h.removeprefix(prefix): np.array([np.nan if r[i] == "" else float(r[i]) for r in rows[1:]])
            for i, h in enumerate(rows[0][1:], start=1)}
    return x, cols


def test_every_donor_overshoots_after_the_step_down():
    t, step = _matrix("step_300_to_1_human.csv", "stress_Pa_")
    assert len(step) == 20 and t[0] == np.float64(1.00037)
    for d, v in step.items():
        vv = v[np.isfinite(v)]
        imax = int(np.argmax(vv))
        assert 0 < imax < vv.size - 1, d
        assert 4.0 <= t[imax] <= 16.5, (d, t[imax])   # samples sit at k + 0.00037 s
        assert vv[-1] / vv[imax] < 0.94, (d, vv[-1] / vv[imax])   # decay of at least 6 % after the peak


def test_steady_one_per_second_point_is_a_sample_of_the_transient():
    t, step = _matrix("step_300_to_1_human.csv", "stress_Pa_")
    rates, steady = _matrix("steady_shear_human.csv", "stress_Pa_")
    assert rates[4] == 1.0
    exact = [d for d in step if np.any(np.isclose(step[d][np.isfinite(step[d])], steady[d][4], rtol=1e-12, atol=0))]
    assert len(exact) == 15 and "C" in exact and "D" not in exact
    near = [d for d in step if d not in exact
            and np.any(np.isclose(step[d][np.isfinite(step[d])], steady[d][4], rtol=1e-4, atol=0))]
    assert sorted(near) == ["A", "B", "F", "J"]
