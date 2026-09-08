"""Item 17's proposal: `data/dstfit_upstream_proposal.md`.

The proposal's recommendation rests on one factual claim — that **no shipped `dst_fit`
script reaches the profile walk except `gate_n7.py`, and only for `phi_m`** — and on one
arithmetic claim, that `phi_m`'s exposure is below 1 at both ends. Both are recomputed here
rather than read off the page, because a decision document that quietly stops being true is
worse than one that was never written.

This file asserts nothing about whether the fix SHOULD happen. That is the founder's
decision (D5) and the proposal says so.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
DSTFIT = REPO.parent / "dst_fit"
DOC = REPO / "data" / "dstfit_upstream_proposal.md"


def _hits(pattern: str, sub: str) -> list[str]:
    rx = re.compile(pattern)
    out = []
    for f in sorted((DSTFIT / sub).rglob("*.py")):
        for k, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            if not line.lstrip().startswith("#") and rx.search(line):
                out.append(f"{f.name}:{k}")
    return out


@pytest.fixture(scope="module")
def doc() -> str:
    assert DOC.exists(), "run scripts/dstfit_upstream_proposal.py"
    return DOC.read_text(encoding="utf-8")


def test_no_shipped_dstfit_script_passes_the_profile_kwarg():
    """The claim the whole recommendation turns on. If this ever fails, a shipped
    `dst_fit` result has started going through the defective walk and item 17 must be
    re-decided, not re-read."""
    assert _hits(r"profile=", "scripts") == []


def test_the_only_shipped_profile_caller_is_gate_n7_on_phi_m():
    hits = _hits(r"profile_likelihood\(", "scripts")
    assert hits and all(h.startswith("gate_n7.py") for h in hits), hits
    src = (DSTFIT / "scripts" / "gate_n7.py").read_text(encoding="utf-8")
    assert re.findall(r'profile_likelihood\([^)]*?"([a-z_]+)"', src, re.S) == ["phi_m", "phi_m"]


def test_phi_m_exposure_is_below_one_at_both_ends():
    """`exposure = (span/5000)/(wald/4)`; the box floor wins only above 1."""
    from dst_fit.models import WYART_CATES
    i = WYART_CATES.param_names.index("phi_m")
    lo = float(WYART_CATES.bounds[0][i])
    hi = float(WYART_CATES.bounds[1][i])
    best, wald = 0.5876, 0.0013           # gate_n7.py's own quoted figure
    for span in (abs(best - lo), abs(hi - best)):
        assert (span / 5000) / (wald / 4) < 1.0, span


def test_the_quoted_phi_m_precision_still_appears_in_gate_n7():
    """The exposure arithmetic uses `phi_m = 0.5876 +/- 0.0013` from `gate_n7.py`'s own
    docstring. If that figure changes, the table in Step 2 is about a different fit."""
    src = (DSTFIT / "scripts" / "gate_n7.py").read_text(encoding="utf-8")
    assert "0.5876" in src and "0.0013" in src


def test_sigma_star_has_the_wide_box_the_proposal_warns_about():
    """The discriminating half: if every box were O(1) the proposal's Step 3 risk would not
    exist, and the recommendation would change."""
    from dst_fit.models import WYART_CATES
    i = WYART_CATES.param_names.index("sigma_star")
    assert float(WYART_CATES.bounds[1][i]) >= 1e4


def test_the_proposal_states_it_changes_nothing(doc):
    assert "This document changes nothing" in doc
    assert "PROPOSAL" in doc
    assert "do not merge" in doc.lower() or "not a patch" in doc.lower()


def test_the_proposal_carries_both_sides_and_a_recommendation(doc):
    assert "## The case FOR fixing it upstream" in doc
    assert "## The case AGAINST, which is not weak" in doc
    assert "## Recommendation" in doc
    assert "If it is rejected, that is a defensible position" in doc


def test_the_proposal_does_not_claim_a_shipped_number_is_wrong(doc):
    """Step 2's finding is that the defect is INERT upstream. A later edit that overstated
    it would make the document argue harder than its own evidence supports."""
    assert "the defect is inert for the only number `dst_fit` currently publishes" in doc
    assert "not because anything is currently wrong" in doc
