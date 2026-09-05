"""Regression for ``case_studies/injectability_surrogates/run.py``.

Two anchors that do not pass through the case study's own arithmetic: the share of Q
outside the envelope has a closed form for a Newtonian fluid, and the nominal forces
written to ``results.md`` must be reproduced by ``injection_force`` from a fresh refit.
The full script (with UQ) is not re-run here; its png md5 is recorded in results.md.
"""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import numpy as np
import pytest

from biorheo.apps import injection_force
from biorheo.constitutive import get_model
from biorheo.identify import fit_rate_curve

CASE = Path(__file__).resolve().parents[1] / "case_studies" / "injectability_surrogates"
HB = get_model("herschel_bulkley")


@pytest.fixture(scope="module")
def cs():
    spec = importlib.util.spec_from_file_location("cs_run", CASE / "run.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_share_outside_envelope_newtonian_closed_form(cs):
    # Newtonian (tau_y = 0, n = 1, k = mu): gd(r) = tw r / (R mu), so the share of Q at
    # r > r* is 1 - (r*/R)^4 with r*/R = gd_max mu / tw.  R = 0.1 mm, L = 25 mm:
    # dp = 1e5 Pa -> tw = 1e5 * 1e-4 / 0.05 = 200 Pa; mu = 0.5 -> gd_w = 400 1/s.
    # gd_max = 200 -> r*/R = 0.5 -> share = 1 - 0.0625 = 0.9375
    share = cs.flow_fraction_outside_envelope(1e5, np.array([0.0, 0.5, 1.0]), 200.0)
    assert share == pytest.approx(0.9375, rel=1e-6)
    # gd_max = 100 -> r*/R = 0.25 -> 1 - 0.00390625 = 0.99609375
    assert cs.flow_fraction_outside_envelope(1e5, np.array([0.0, 0.5, 1.0]), 100.0) == pytest.approx(0.99609375, rel=1e-6)
    # wall rate inside the envelope -> nothing is extrapolated
    assert cs.flow_fraction_outside_envelope(1e5, np.array([0.0, 0.5, 1.0]), 400.0) == 0.0


@pytest.mark.skipif(not (CASE / "results.md").exists(), reason="case study not run")
def test_results_md_forces_reproduce(cs):
    text = (CASE / "results.md").read_text(encoding="utf-8")
    sections = re.split(r"^## `(.+?)`$", text, flags=re.M)[1:]
    assert len(sections) == 2 * len(cs.SURROGATES)
    records = cs.load_all(allow_network=False)
    for rid, body in zip(sections[0::2], sections[1::2]):
        rows = re.findall(r"^\| (\d+(?:\.\d+)?) \| [\d.e+-]+ \| ([\d.]+) \|", body, flags=re.M)
        assert [float(t) for t, _ in rows] == list(cs.DOSE_TIMES_S)
        g, t = cs.load(cs.DATA / cs.SURROGATES[rid])
        good = np.isfinite(g) & np.isfinite(t) & (g > 0) & (t > 0)
        rec = records[rid]
        ident = fit_rate_curve(HB, g[good], t[good] / g[good],
                               temperature_C=rec.envelope.temperature_C, geometry=rec.envelope.geometry)
        for t_dose, f_written in rows:
            Q = cs.DOSE_ML * 1e-6 / float(t_dose)
            f = injection_force(Q, cs.NEEDLE, cs.SYRINGE, HB, ident.values, rho=cs.RHO).force
            assert f == pytest.approx(float(f_written), abs=6e-4), (rid, t_dose, f, f_written)
