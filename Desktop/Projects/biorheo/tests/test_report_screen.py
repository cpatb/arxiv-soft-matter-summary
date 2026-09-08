"""`report.build_from_screen`: the deliverable page for SCREEN data (item 25).

One measured viscosity, no flow curve. It exists because item 25 found that the registry
cannot hold such data at all -- a one-parameter model on a single point has zero residual
degrees of freedom -- and refusing to produce a page for it would refuse the shape most
inbound client data actually arrives in.

The tests that matter are the ones checking the page does not quietly look like the fitted
one: the interval must be named a MEASUREMENT interval, and the prose must not say "fitted"
about anything, because nothing was.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from biorheo.apps.report import SCREEN_FLAG, DeviceSpec, build_from_screen, render
from biorheo.registry.schema import Licence, Source

DEVICE = DeviceSpec(
    needle_radius=1.05e-4, needle_length=12.7e-3, barrel_radius=2.3e-3,
    flow_rate=1.0e-6 / 60.0, label="test device",
)
SRC = (Source(doi="10.5281/zenodo.14409284", role="dataset",
              resolved_on="2026-09-08", label="x"),)
CCBY = Licence(name="CC BY 4.0", commercial_use=True)
ETA = 25.343e-3
SD = 0.106448213804752e-3


def make(**kw):
    # annotated so the splat does not collapse every parameter to one union type
    args: dict[str, object] = dict(
        label="s", fluid="f", eta_Pa_s=ETA, spec=DEVICE,
        measured_range=(960.0, 3329.0), temperature_C=25.0,
        geometry="m-VROC", sources=SRC, licence=CCBY)
    args.update(kw)
    return build_from_screen(**args)  # type: ignore[arg-type]


def test_the_flag_replaces_no_sidecar_rather_than_joining_it():
    """`no_sidecar` means somebody forgot to fit; screen data has nothing to fit."""
    d = make(eta_sd_Pa_s=SD)
    assert SCREEN_FLAG in d.flags
    assert "no_sidecar" not in d.flags


def test_the_interval_is_a_measurement_interval_and_says_so():
    d = make(eta_sd_Pa_s=SD)
    assert "MEASUREMENT interval" in d.hull["kind"]
    assert d.hull["measurement_relative_sd"] == pytest.approx(SD / ETA, rel=1e-12)
    # it must NOT claim to be wider than a joint region, the way a parameter hull does
    assert "joint" not in d.hull["wider_than"]
    txt = render(d)
    assert "not** a parameter interval" in txt


def test_the_interval_scales_with_the_sd_multiple_and_is_exact():
    """A Newtonian force is exactly proportional to viscosity, so no propagation
    approximation is involved and the interval must be exactly linear in the SD."""
    d1 = make(eta_sd_Pa_s=SD, sd_multiple=1.0)
    d2 = make(eta_sd_Pa_s=SD, sd_multiple=2.0)
    F = d1.nominal["force_N"]
    assert (d2.hull["hi"] - F) == pytest.approx(2 * (d1.hull["hi"] - F), rel=1e-12)
    assert d1.hull["relative_width"] == pytest.approx(2 * SD / ETA, rel=1e-12)


def test_no_sd_means_no_interval_at_all_and_the_headline_says_so():
    d = make()
    assert "no_replicate_sd" in d.flags
    assert d.hull["lo"] == d.hull["hi"] == d.nominal["force_N"]
    assert d.hull["relative_width"] == 0.0
    assert "no uncertainty of any kind" in render(d)


def test_the_force_is_exactly_proportional_to_the_viscosity():
    """An anchor that needs no reference implementation: Hagen-Poiseuille is linear in eta,
    so doubling the viscosity must double the force to machine precision. A screen path
    that silently changed the physics would fail here."""
    a = make(eta_Pa_s=ETA)
    b = make(eta_Pa_s=2 * ETA)
    assert b.nominal["force_N"] == pytest.approx(2 * a.nominal["force_N"], rel=1e-12)


def test_the_bore_exponent_is_exactly_minus_four():
    """Item 18's closed form for a Newtonian fluid. Not approximately."""
    d = make(eta_sd_Pa_s=SD)
    assert d.bore["exponent"] == pytest.approx(-4.0, abs=1e-9)


def test_the_page_never_calls_anything_fitted():
    """Nothing was fitted. A page that said otherwise would be this repo's recurring
    prose-contradicts-its-own-artifact defect, in the client-facing document."""
    txt = render(make(eta_sd_Pa_s=SD))
    for banned in ("fitted law", "fitted envelope", "identifiability verdict",
                   "profile push-through"):
        assert banned not in txt, banned
    assert "NOT a registry record" in txt
    assert "no flow curve" in txt


def test_the_extrapolation_distance_is_measured_from_the_callers_stated_range():
    d = make(eta_sd_Pa_s=SD, measured_range=(960.0, 3329.0))
    gd_w = d.nominal["wall_shear_rate_1_s"]
    assert d.envelope["decades_past"] == pytest.approx(math.log10(gd_w / 3329.0), rel=1e-12)
    assert "extrapolated" in d.flags
    # a wider stated range must reduce it -- the number is the caller's claim, not ours
    wide = make(eta_sd_Pa_s=SD, measured_range=(960.0, 3.0e4))
    assert wide.envelope["decades_past"] < d.envelope["decades_past"]


def test_it_refuses_input_it_cannot_describe():
    with pytest.raises(ValueError, match="viscosity must be positive"):
        make(eta_Pa_s=0.0)
    with pytest.raises(ValueError, match="SD must be non-negative"):
        make(eta_sd_Pa_s=-1.0)
    with pytest.raises(ValueError, match="increasing positive pair"):
        make(measured_range=(3329.0, 960.0))
    with pytest.raises(ValueError, match="increasing positive pair"):
        make(measured_range=(0.0, 3329.0))


def test_the_licence_comes_from_the_caller_and_is_honoured():
    """Whether a client's own data may appear in a public artifact is the client's answer."""
    nc = build_from_screen(
        label="s", fluid="f", eta_Pa_s=ETA, spec=DEVICE, measured_range=(960.0, 3329.0),
        temperature_C=25.0, geometry="g", sources=SRC,
        licence=Licence(name="proprietary", commercial_use=False))
    assert "non_commercial" in nc.flags
    assert "NOT LICENSABLE" in render(nc)
    assert "NOT LICENSABLE" not in render(make(eta_sd_Pa_s=SD))


def test_the_wall_rate_does_not_depend_on_the_viscosity():
    """4Q/(pi R^3). If a screen page's rate moved with the fluid, the device model is wrong."""
    rates = [make(eta_Pa_s=e).nominal["wall_shear_rate_1_s"]
             for e in (5e-3, 25e-3, 50e-3)]
    assert np.allclose(rates, rates[0], rtol=1e-9)
