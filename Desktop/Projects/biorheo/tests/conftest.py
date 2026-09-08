"""Shared fixtures, and the frozen registry that CLOSED gates must be graded on.

Item 22 added the registry's first commercially licensable records (DMEM + FBS, CC BY 4.0).
Item 11 had already shown what happens when the registry grows underneath a gate that was
decided on an earlier version of it: two closed gates were silently re-decided, and the
suite stayed green throughout, because every verdict was still true **on a different
question**.

`closed_gate_records` is the answer. A gate that was closed against the 115-record registry
of 2026-09-07 keeps being graded on those 115 records, whatever is added later. A gate that
is *about* the registry as it stands today should call `load_all()` directly and say so.

The rule from `data/dmem_record_prereg.md` §4: **no existing gate may change its VALUE**.
Any gate whose sample would grow gets this list instead.
"""

from __future__ import annotations

import pytest

from biorheo.registry import load_all

# The two families in the registry as of 2026-09-07, before item 22. A family added later
# is excluded by default, which is the safe direction for a gate that is already closed.
CLOSED_GATE_FAMILIES = ("blood_horner2020_", "carbopol_etd2020_")
N_CLOSED_GATE_RECORDS = 115


def _closed_gate_records() -> dict:
    recs = {k: v for k, v in load_all().items()
            if k.startswith(CLOSED_GATE_FAMILIES)}
    assert len(recs) == N_CLOSED_GATE_RECORDS, (
        f"the frozen registry for closed gates has {len(recs)} records, not "
        f"{N_CLOSED_GATE_RECORDS}. Something was removed or renamed in a family a closed "
        "gate was decided on; that is a change to the gate, not to the fixture."
    )
    return recs


@pytest.fixture(scope="session")
def closed_gate_records() -> dict:
    """The 115 records of 2026-09-07: 100 Horner blood, 15 Carbopol ETD2020."""
    return _closed_gate_records()
