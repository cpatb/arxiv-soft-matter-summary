"""Make test_pipeline.py runnable under pytest as well as standalone.

test_pipeline.py was written for its own `main()` runner, which hands every
temp-dir test a freshly made directory as a positional argument named `tmp`.
pytest reads that same parameter as a fixture request, finds no fixture called
`tmp`, and errors out 13 of the tests at setup -- they were never failing, they
were never running. This supplies the missing fixture so `pytest` and
`python test_pipeline.py` agree.
"""
from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture()
def tmp(tmp_path: Path) -> Path:
    """The per-test scratch directory test_pipeline.py's own runner provides."""
    return tmp_path
