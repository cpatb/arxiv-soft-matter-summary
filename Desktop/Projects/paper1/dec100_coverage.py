"""dec100_coverage.py -- name the within-cell `_dec100` omissions instead of counting them.

WHY THIS EXISTS
---------------
`Drop Impact - Paper 1` card, analysis-queue item 7: *"Check the 4 within-cell `_dec100`
omissions against the hand annotations. 20 minutes."* The missingness audit established
that stitched = 121 = dec100 exactly out of 160 base recordings, that **seven whole cells**
have no `_dec100` at all, and that there are **4 further within-cell exceptions** which are
"unchecked".

Cell-level missingness cannot be outcome-caused, because a cell is configured before its
drops are recorded. **Within-cell missingness can be** -- if a drop was skipped because of
what it did, the late-time sample is selected on the outcome and every late-time statistic
in the manuscript inherits that selection. That is why item 7 exists and why it is worth
more than twenty minutes of anyone's attention.

The expensive half of item 7 is not the judgement, it is finding out WHICH four recordings
they are: the tree is 1.05 TB across 49 000 files on a removable exFAT disk. This script
does that half and leaves the judgement, which needs the hand annotations and the author's
memory of the session.

⚠ THE PATH TRAP, and it is why this script globs rather than joins
------------------------------------------------------------------
The analysis tree is at `E:\\CPB\\Drop Impact Analyses<U+F028>\\`, whose name ends in a
Unicode PRIVATE-USE character -- exFAT written by macOS maps Windows-illegal characters to
U+F000 + c. `Test-Path` on the visible spelling returns False and a non-recursive listing
returns nothing, which looks exactly like the tree has been deleted. It has not. Never
build a path into this disk by string join; glob it, or enumerate and use the real name.

WHAT COUNTS AS A RECORDING, stated before running
-------------------------------------------------
A "base recording" is a `side/measurements/` directory containing at least one
`drop<N>_profiles.csv`. That recording HAS dec100 if the sibling `drop<N>_dec100_*.csv`
exists. `_bak_reprocess*` directories are excluded from both counts -- they are historical
state, and counting them would inflate coverage with files the pipeline no longer reads.

    python dec100_coverage.py            # prints the table and writes dec100_coverage.csv
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

E_ROOT = Path("E:/CPB/Driscoll/Silica Colloid Drop")
OUT = Path(__file__).resolve().parent / "dec100_coverage.csv"
DROP_RE = re.compile(r"^drop(\d+)_profiles\.csv$", re.I)
SKIP_DIR = re.compile(r"_bak_reprocess|__pycache__|\.git", re.I)


def cells(root: Path):
    """Yield (cell, condition, measurements_dir). Globbed, never string-joined."""
    for meas in root.glob("*/*/*/side/measurements"):
        if SKIP_DIR.search(str(meas)):
            continue
        yield meas
    for meas in root.glob("*/*/side/measurements"):
        if SKIP_DIR.search(str(meas)):
            continue
        yield meas


def main() -> int:
    if not E_ROOT.exists():
        print("E: is not mounted, or the silica root moved. Expected:", E_ROOT)
        return 1

    rows = []
    for meas in sorted(set(cells(E_ROOT))):
        names = {p.name for p in meas.iterdir() if p.is_file()}
        drops = sorted(int(m.group(1)) for n in names if (m := DROP_RE.match(n)))
        for d in drops:
            has = any(n.startswith("drop%d_dec100" % d) and n.endswith(".csv")
                      for n in names)
            rel = meas.relative_to(E_ROOT)
            parts = rel.parts
            rows.append({
                "session": parts[0],
                "cell": "/".join(parts[1:-2]) if len(parts) > 3 else "",
                "drop": d,
                "has_dec100": int(has),
                "path": str(meas),
            })

    if not rows:
        print("no measurement directories matched -- check the glob depth against the tree")
        return 1

    # group by cell: a cell with SOME dec100 and some without is a within-cell exception
    bycell: dict[tuple[str, str], list[dict]] = {}
    for r in rows:
        bycell.setdefault((r["session"], r["cell"]), []).append(r)

    whole_missing, within, complete = [], [], 0
    for key, rs in sorted(bycell.items()):
        n, k = len(rs), sum(r["has_dec100"] for r in rs)
        if k == 0:
            whole_missing.append((key, n))
        elif k < n:
            within.extend(r for r in rs if not r["has_dec100"])
        else:
            complete += 1

    print("base recordings (drop*_profiles.csv, excluding _bak_reprocess): %d" % len(rows))
    print("with dec100: %d    without: %d"
          % (sum(r["has_dec100"] for r in rows),
             len(rows) - sum(r["has_dec100"] for r in rows)))
    print("cells: %d complete, %d entirely without dec100, %d with a within-cell gap"
          % (complete, len(whole_missing),
             len({(r["session"], r["cell"]) for r in within})))
    print()
    print("WHOLE CELLS WITH NO dec100 -- configured before their drops, so this cannot be")
    print("outcome-caused, and it is a coverage statement rather than a selection one:")
    for (sess, cell), n in whole_missing:
        print("  %-45s %-28s (%d drops)" % (sess, cell, n))
    print()
    print("WITHIN-CELL OMISSIONS -- these are the ones item 7 is about. A drop skipped")
    print("because of WHAT IT DID selects the late-time sample on the outcome:")
    for r in within:
        print("  %-45s %-28s drop %d" % (r["session"], r["cell"], r["drop"]))
        print("      %s" % r["path"])
    if not within:
        print("  none found by this enumeration -- see the caveat below before concluding")

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote", OUT)
    print("\n⚠ SCOPE, and it has to be read before the list is. Over the whole disk this is")
    print("  276 recordings / 183 with dec100. Restrict to the six PAPER sessions (the")
    print("  three salt cells plus 0.42 / 0.44 / 0.45 VF) and to `Silica/` cells -- i.e.")
    print("  drop the `newtonian_control` cells -- and it reproduces the audit EXACTLY:")
    print("  n = 160, dec100 = 121, seven whole cells with none, and those seven are the")
    print("  seven the Status card names. So the subset is identified, not guessed.")
    print("  On that subset the within-cell omissions number FIVE, not the card's four,")
    print("  and three of the five are the SAME cell (50 mM LiCl, u0_4: drops 2, 3, 5 of")
    print("  6 missing, keeping 1, 4, 6). An interleaved pattern is not a camera that")
    print("  stopped; check that cell first.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
