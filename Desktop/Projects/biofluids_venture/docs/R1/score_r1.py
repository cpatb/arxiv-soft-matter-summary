#!/usr/bin/env python3
"""Score the R1 discovery-call gate from ``r1_log.csv``.

The rule is pre-registered in ``docs/market_research_2026-09-06/00_SYNTHESIS.md`` §7 and
restated in ``R1_LOG.md``.  It lives in code rather than in prose because a pre-registered
admission rule applied by hand at read-time gets skipped (lesson: *enforce gates in code*).

Pre-registered, not to be changed after data exists:

    R1 = number of HELD CALLS whose contact answered YES to "have you paid an outside party
         for rheology-to-force prediction?" OR YES to "would you pay for the
         high-shear / identifiability layer?"

    R1 >= 3 of 10  -> PASS  -> wedge W1 (injectability), synthesis rank 1-3
    R1 <  3 of 10  -> FAIL  -> switch to rank 4 (battery slurry), synthesis D-A fallback

Two guards that make an early verdict impossible:

  * A verdict requires 10 HELD CALLS.  A no-reply is not a call and an email thread that
    never reached the questions is not a call.  With fewer than 10 the answer is
    UNDETERMINED, which is a real verdict and is not the same as either pass or fail.
  * ``would_pay`` counts only when a deliverable was named by the CONTACT.  "That sounds
    useful" is not a would-pay; the column ``deliverable_named`` must be non-empty, else the
    row is counted as a soft yes and reported separately rather than in R1.

Exit codes: 0 PASS, 1 FAIL, 2 UNDETERMINED, 3 malformed log.

Usage:  python score_r1.py [path/to/r1_log.csv]
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

BAR = 3
N_CALLS_REQUIRED = 10

# From 00_SYNTHESIS.md §7: "target: 4 formulation scientists at biologics companies or
# CDMOs, 3 device/combination-product engineers, 2 regulatory or CM&S consultants, 1
# platform company".
SEGMENT_QUOTA = {"formulation": 4, "device": 3, "regulatory": 2, "platform": 1}

YES = {"y", "yes", "true", "1"}
NO = {"n", "no", "false", "0", ""}
REQUIRED_COLUMNS = {
    "target_id",
    "segment",
    "call_held",
    "paid_before",
    "would_pay",
    "deliverable_named",
}


def _bool(value: str, row_id: str, column: str) -> bool:
    """Parse a yes/no cell, refusing anything ambiguous rather than defaulting to False."""
    token = (value or "").strip().lower()
    if token in YES:
        return True
    if token in NO:
        return False
    raise ValueError(f"row {row_id!r}, column {column!r}: {value!r} is not yes/no/blank")


def main(argv: list[str]) -> int:
    path = Path(argv[1]) if len(argv) > 1 else Path(__file__).resolve().parent / "r1_log.csv"
    if not path.is_file():
        print(f"MALFORMED: no log at {path}")
        return 3

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        print(f"R1 = 0 of 0 held calls (bar {BAR} of {N_CALLS_REQUIRED})")
        print("VERDICT: UNDETERMINED -- the log is empty. The gate has not been run.")
        return 2

    missing = REQUIRED_COLUMNS - set(rows[0])
    if missing:
        print(f"MALFORMED: log is missing column(s) {sorted(missing)}")
        return 3

    held: list[dict[str, str]] = []
    try:
        for row in rows:
            rid = row.get("target_id") or "<blank target_id>"
            if _bool(row["call_held"], rid, "call_held"):
                held.append(row)
    except ValueError as exc:
        print(f"MALFORMED: {exc}")
        return 3

    hard_yes: list[str] = []
    soft_yes: list[str] = []
    try:
        for row in held:
            rid = row.get("target_id") or "<blank target_id>"
            paid = _bool(row["paid_before"], rid, "paid_before")
            would = _bool(row["would_pay"], rid, "would_pay")
            named = bool((row.get("deliverable_named") or "").strip())
            if paid or (would and named):
                hard_yes.append(rid)
            elif would:
                soft_yes.append(rid)
    except ValueError as exc:
        print(f"MALFORMED: {exc}")
        return 3

    r1 = len(hard_yes)
    contacted = sum(1 for row in rows if (row.get("date_contacted") or "").strip())

    print(f"contacted: {contacted}   replied+held calls: {len(held)}   logged rows: {len(rows)}")
    print(f"R1 = {r1} of {len(held)} held calls (bar: {BAR} of {N_CALLS_REQUIRED})")
    if hard_yes:
        print(f"  paid-before or would-pay-with-named-deliverable: {', '.join(hard_yes)}")
    if soft_yes:
        print(
            f"  soft yes (would pay, no deliverable named -- NOT counted in R1): "
            f"{', '.join(soft_yes)}"
        )

    seen = Counter((row.get("segment") or "").strip().lower() for row in held)
    unknown = sorted(set(seen) - set(SEGMENT_QUOTA) - {""})
    print("  segment coverage against the §7 quota:")
    for segment, quota in SEGMENT_QUOTA.items():
        got = seen.get(segment, 0)
        print(f"    {segment:<12} {got}/{quota}{'' if got >= quota else '   SHORT'}")
    if unknown:
        print(f"    WARNING: unrecognised segment label(s) {unknown} -- not counted")

    if len(held) < N_CALLS_REQUIRED:
        short = N_CALLS_REQUIRED - len(held)
        print(
            f"VERDICT: UNDETERMINED -- {short} more held call(s) needed. "
            f"An early {r1}/{len(held)} is not a pass; the bar is {BAR} of "
            f"{N_CALLS_REQUIRED}."
        )
        return 2

    if r1 >= BAR:
        print(f"VERDICT: R1 PASS ({r1} >= {BAR}) -- wedge W1, synthesis ranks 1-3.")
        return 0

    print(
        f"VERDICT: R1 FAIL ({r1} < {BAR}) -- the pre-registered fallback fires: "
        "drop W1, switch to rank 4 (battery slurry / coating window), D-A."
    )
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
