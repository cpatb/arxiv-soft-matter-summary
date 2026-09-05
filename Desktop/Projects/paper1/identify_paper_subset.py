"""identify_paper_subset.py -- which sessions did the missingness audit actually use?

WHY
---
dec100_coverage.py used to PRINT that restricting to "the six PAPER sessions (the three
salt cells plus 0.42 / 0.44 / 0.45 VF) and to Silica/ cells" reproduces the audit exactly
at n = 160, dec100 = 121, seven whole cells with none. Applying the stated restriction in
code gives 17 whole cells and 11 within-cell omissions, so the printed description is not
the restriction the audit used -- or the tree has changed since.

This script does not guess. It enumerates the sessions actually on disk, then searches
subsets for ones reproducing the audit triple, and REPORTS HOW MANY MATCH. If more than
one subset matches, the subset is not identified by those three numbers and saying so is
the result; a single match is an identification.

Note on the risk: searching subsets to hit a target is overfitting if the target is a
physical result. Here the target is a bookkeeping figure from a prior audit and the
question is purely "which files did it count", so the search is legitimate -- but the
multiplicity is reported precisely because a unique match is what makes it an answer.

    python identify_paper_subset.py
"""

from __future__ import annotations

import csv
import itertools
import sys
from collections import defaultdict
from pathlib import Path

CSV = Path(__file__).resolve().parent / "dec100_coverage.csv"

TARGET_N = 160
TARGET_DEC100 = 121
TARGET_WHOLE = 7

# Only silica cells are in the manuscript; the audit's own description says so and this
# part of the restriction is not in doubt.
def is_silica(cell: str) -> bool:
    return "silica" in cell.lower() and "newtonian_control" not in cell.lower()


def main() -> int:
    if not CSV.is_file():
        print("run dec100_coverage.py first to produce", CSV)
        return 1

    rows = [r for r in csv.DictReader(CSV.open(encoding="utf-8"))
            if is_silica(r["cell"])]
    for r in rows:
        r["has_dec100"] = int(r["has_dec100"])

    bysess = defaultdict(list)
    for r in rows:
        bysess[r["session"]].append(r)

    print("Silica recordings per session (%d sessions, %d recordings total)\n"
          % (len(bysess), len(rows)))
    print("%-46s %6s %8s %7s" % ("session", "n", "dec100", "cells"))
    for s in sorted(bysess):
        rs = bysess[s]
        ncell = len({r["cell"] for r in rs})
        print("%-46s %6d %8d %7d"
              % (s, len(rs), sum(r["has_dec100"] for r in rs), ncell))

    def stats(sessions):
        sub = [r for s in sessions for r in bysess[s]]
        bycell = defaultdict(list)
        for r in sub:
            bycell[(r["session"], r["cell"])].append(r)
        whole = sum(1 for rs in bycell.values()
                    if sum(x["has_dec100"] for x in rs) == 0)
        return len(sub), sum(r["has_dec100"] for r in sub), whole

    names = sorted(bysess)
    matches = []
    # the audit describes SIX sessions; search 4..8 so a miscount in the prose cannot
    # hide a real match
    for k in range(4, 9):
        for combo in itertools.combinations(names, k):
            if stats(combo) == (TARGET_N, TARGET_DEC100, TARGET_WHOLE):
                matches.append(combo)

    print("\n" + "=" * 74)
    print("Subsets reproducing the audit (n=%d, dec100=%d, %d whole cells missing)"
          % (TARGET_N, TARGET_DEC100, TARGET_WHOLE))
    print("=" * 74)
    if not matches:
        print("NONE. No subset of 4-8 sessions reproduces the audit triple, so either")
        print("the tree has changed since the audit or the audit counted something this")
        print("enumeration does not (e.g. it included _bak_reprocess, or counted files")
        print("rather than recordings). The audit's numbers should not be quoted as")
        print("reproducible until that is resolved.")
        return 1

    for combo in matches:
        print("  [%d sessions] %s" % (len(combo), "; ".join(combo)))
    if len(matches) == 1:
        print("\nUNIQUE -- the subset is identified. Put this session list in")
        print("dec100_coverage.py as an explicit tuple, not a prose description.")
        return 0
    print("\n%d subsets match, so these three numbers do NOT identify the subset."
          % len(matches))
    print("Pin it from the Status card's session list instead of from the counts.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
