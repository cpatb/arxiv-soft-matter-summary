"""The authors' own steadiness rules, read from the sources rather than chosen here.

A calibration that quietly drops the curves it does not like is choosing its
answer. So the admission rule comes from the depositing labs:

- **Richards et al. 2019** colour-code their own deposited ``Figure1.py``:
  ``colours = ['tab:blue', ...x4, 'tab:red', 'tab:red', 'tab:red']`` against
  ``labelList``, with red the shear-jammed regime their paper is *about*
  oscillating in. The red curves are time averages of an oscillating rate. The
  list is parsed out of their script rather than hard-coded, so if the deposit
  ever stops saying it, this raises instead of silently keeping a stale rule.
- **Guy, Hermes & Poon 2015** state it in words: "the flow in both (a) and (c)
  was unsteady for phi >= 0.56, and points represent temporal averages", with
  massive wall slip above ``phi_m``.

This duplicates the rule ``dst_fit``'s gate N8 applies. It is repeated here
rather than imported because N8 is a gate script, not a library, and because a
calibration whose admission rule lives in another project's ``scripts/``
directory is one refactor away from being applied silently and wrongly.
"""
from __future__ import annotations

import re
from pathlib import Path

DST_FIT = Path(r"C:\Users\compl\Desktop\Projects\dst_fit")
RICHARDS_SCRIPT = (DST_FIT / "Competing Timescales Lead to Oscillations in Shear-Thickening Suspensions"
                   / "Figure1" / "Plotting" / "Figure1.py")

# Guy, Hermes & Poon, PRL 115, 088304 (2015), stated in the Figure 1 caption.
GUY_UNSTEADY_AT = 0.56


def richards_jammed_labels() -> tuple[str, ...]:
    """Weight-fraction labels the authors' own script colours as shear-jammed."""
    text = RICHARDS_SCRIPT.read_text(encoding="utf-8", errors="replace")
    colours = re.search(r"colours\s*=\s*\[([^\]]+)\]", text)
    labels = re.search(r"labelList\s*=\s*\[([^\]]+)\]", text)
    if not colours or not labels:
        raise ValueError(f"{RICHARDS_SCRIPT.name}: cannot find colours/labelList")
    cols = [c.strip().strip("'\"") for c in colours.group(1).split(",")]
    labs = [x.strip().strip("'\"") for x in labels.group(1).split(",")]
    if len(cols) != len(labs):
        raise ValueError(
            f"{RICHARDS_SCRIPT.name}: {len(cols)} colours for {len(labs)} labels")
    return tuple(lab for lab, c in zip(labs, cols) if c == "tab:red")


def is_steady(dataset: str, phi: float) -> bool:
    """Would the source lab call a curve at this packing a steady flow curve?"""
    if dataset.startswith("guy"):
        return float(phi) < GUY_UNSTEADY_AT
    if dataset == "richards_cornstarch":
        jammed = tuple(float(x) for x in richards_jammed_labels())
        if not jammed:
            raise ValueError("Richards' script no longer colours any curve red")
        return float(phi) < min(jammed)
    raise KeyError(f"no steadiness rule recorded for {dataset!r}")
