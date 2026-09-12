"""test_cap.py — the PJ_MAX_ITER iteration-cap override, exercised at BOTH settings.

The repo's standing rule is that a knob's gate must exercise the knob's NON-default
setting on the channel the claims are made from. The claim this knob supports is that a
re-run at a higher cap can be told apart from the capped record it replaces, so the
tests below check the KEY (which is what analyse_phase2's dedup and phase2's resume both
read), not just that the function returns a number.

No Yade: `packing.max_iter()` and `phase2_phiJ.cell_key()` are both engine-independent.
"""

from __future__ import annotations

import os

import pytest

import packing
import phase2_phiJ


@pytest.fixture(autouse=True)
def _clean_env():
    """PJ_MAX_ITER must not leak between tests, or one test sets another's answer."""
    old = os.environ.pop("PJ_MAX_ITER", None)
    yield
    os.environ.pop("PJ_MAX_ITER", None)
    if old is not None:
        os.environ["PJ_MAX_ITER"] = old


def test_default_cap_is_unchanged():
    assert packing.max_iter() == packing.MAX_ITER == 200_000


def test_override_is_read_at_call_time_not_import_time():
    os.environ["PJ_MAX_ITER"] = "750000"
    assert packing.max_iter() == 750_000
    os.environ["PJ_MAX_ITER"] = "300000"
    assert packing.max_iter() == 300_000, "the cap was cached at import"


def test_empty_string_falls_back_to_the_default():
    # `export PJ_MAX_ITER=` in a launcher must not become int("") -> ValueError.
    os.environ["PJ_MAX_ITER"] = ""
    assert packing.max_iter() == packing.MAX_ITER


def test_default_key_is_byte_identical_to_the_keys_on_disk():
    # This exact string was read out of out/phase2_*.jsonl on 2026-09-09. If the key
    # format drifts, every record ever written stops being resumable and dedup breaks.
    assert (phase2_phiJ.cell_key(0.5, 0.0, 1, 1e-3, 1500)
            == "cv0.500_mu0.000_s1_p1.000e-03_N1500")


def test_a_non_default_cap_makes_a_DIFFERENT_key():
    """The discriminating half: without this the re-run would collide with its target."""
    default = phase2_phiJ.cell_key(0.5, 0.0, 1, 1e-3, 1500)
    os.environ["PJ_MAX_ITER"] = "1000000"
    capped = phase2_phiJ.cell_key(0.5, 0.0, 1, 1e-3, 1500)
    assert capped != default
    assert capped.endswith("_it1000000")
    assert capped.startswith(default)


def test_setting_the_cap_to_its_own_default_leaves_the_key_alone():
    """A launcher that exports the default value must not orphan every existing key."""
    os.environ["PJ_MAX_ITER"] = str(packing.MAX_ITER)
    assert (phase2_phiJ.cell_key(0.5, 0.0, 1, 1e-3, 1500)
            == "cv0.500_mu0.000_s1_p1.000e-03_N1500")


def test_the_family_suffix_still_comes_before_the_cap_suffix():
    os.environ["PJ_MAX_ITER"] = "1000000"
    k = phase2_phiJ.cell_key(0.5, 0.0, 1, 1e-3, 1500, "uniform")
    assert k == "cv0.500_mu0.000_s1_p1.000e-03_N1500_uniform_it1000000"


def test_state_records_the_cap_and_the_tolerance_flag():
    """A cap that is not written into the record is a knob nobody can audit later."""
    st = packing.PackingState(n_requested=1500, cv=0.0, mu_p=0.0,
                              p_hat_target=1e-3, seed=1)
    assert st.max_iter == packing.MAX_ITER
    assert st.tol_met is False
    import json
    rec = json.loads(st.to_json())
    assert "max_iter" in rec and "tol_met" in rec
