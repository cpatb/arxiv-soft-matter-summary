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

The expensive half of item 7 is not the judgement, it is finding out WHICH recordings they
are: the tree is 1.05 TB across 49 000 files on a removable exFAT disk. This script does
that half and leaves the judgement, which needs the hand annotations and the author's
memory of the session.

Two things it establishes that the card does not. The card says FOUR within-cell
omissions; on the audit's own subset the enumeration finds FIVE. And three of the five are
the same cell (50 mM LiCl, u0_4: drops 2, 3, 5 of 6 missing, keeping 1, 4, 6) -- an
interleaved pattern, which is not what a camera that stopped looks like. Check that cell
first.

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
    return report_paper_subset(rows)


# ---------------------------------------------------------------------------
# The paper subset.
#
# This block used to be a PRINTED PARAGRAPH asserting "n = 160, dec100 = 121, seven
# whole cells" while the script above computed only the whole-disk figures (276/183).
# Prose that states numbers its own script never computes is exactly the failure mode
# this project keeps hitting, so the restriction is applied in code here and the audit's
# numbers are CHECKED rather than recited. If the tree changes, this fails loudly.
# ---------------------------------------------------------------------------

# The six sessions the audit counted, listed explicitly rather than described.
#
# The old prose said "the three salt cells plus 0.42 / 0.44 / 0.45 VF", which is a
# correct DESCRIPTION but not a usable filter: matching on "silica" instead pulls in the
# 800 nm, Silicone Dioxide and Granular sessions and gives 255 recordings / 17 whole
# cells rather than 160 / 7. identify_paper_subset.py searched every 4-8 session
# combination and found this one is the UNIQUE subset reproducing the audit triple
# (n=160, dec100=121, 7 whole cells), so it is identified, not guessed.
PAPER_SESSIONS = (
    "Silica Colloid 0.42 VF 23JUL26",
    "Silica Colloid 0.44 VF 23JUN26",
    "Silica Colloid 0.44 VF 10 mM LiCl 24JUN26",
    "Silica Colloid 0.44 VF 50 mM LiCl 13JUL26",
    "Silica Colloid 0.44 VF 100 mM LiCl 23JUN26",
    "Silica Colloid 0.45 VF 24JUL26",
)
# `newtonian_control` is a different experiment living in the same session directories.
PAPER_EXCLUDE = re.compile(r"newtonian_control", re.I)

# What the missingness audit and the Status card recorded, as a regression check.
EXPECTED = {"n": 160, "dec100": 121, "whole_cells_missing": 7}


def in_paper_subset(r: dict) -> bool:
    return (r["session"] in PAPER_SESSIONS
            and not PAPER_EXCLUDE.search("%s/%s" % (r["session"], r["cell"])))


def report_paper_subset(rows: list) -> int:
    sub = [r for r in rows if in_paper_subset(r)]
    if not sub:
        print("\nPAPER SUBSET IS EMPTY -- the include/exclude patterns no longer match")
        print("the tree's layout. Not reporting whole-disk numbers as if they were the")
        print("paper's; fix the patterns first.")
        return 1

    bycell: dict = {}
    for r in sub:
        bycell.setdefault((r["session"], r["cell"]), []).append(r)

    whole, within = [], []
    for key, rs in sorted(bycell.items()):
        n, k = len(rs), sum(r["has_dec100"] for r in rs)
        if k == 0:
            whole.append((key, n))
        elif k < n:
            within.extend(r for r in rs if not r["has_dec100"])

    got = {"n": len(sub),
           "dec100": sum(r["has_dec100"] for r in sub),
           "whole_cells_missing": len(whole)}

    print("\n" + "=" * 74)
    print("PAPER SUBSET -- the %d audit sessions, newtonian_control excluded"
          % len(PAPER_SESSIONS))
    print("=" * 74)
    print("Whole disk: %d recordings / %d with dec100." % (
        len(rows), sum(r["has_dec100"] for r in rows)))
    print("Paper subset: %d recordings / %d with dec100, in %d cells."
          % (got["n"], got["dec100"], len(bycell)))
    print()
    ok = True
    for k, want in EXPECTED.items():
        mark = "OK " if got[k] == want else "MISMATCH"
        if got[k] != want:
            ok = False
        print("  %-22s expected %4d   got %4d   %s" % (k, want, got[k], mark))

    if not ok:
        print("\n*** The subset no longer reproduces the audit. Either the tree changed")
        print("    or the restriction above is not the one the audit used. Do NOT quote")
        print("    the within-cell list below until this is resolved.")

    print("\nWHOLE CELLS WITH NO dec100 (%d) -- configured before their drops, so this"
          % len(whole))
    print("cannot be outcome-caused:")
    for (sess, cell), n in whole:
        print("  %-45s %-28s (%d drops)" % (sess, cell, n))

    print("\nWITHIN-CELL OMISSIONS (%d) -- these are what analysis-queue item 7 is about."
          % len(within))
    percell: dict = {}
    for r in within:
        percell.setdefault((r["session"], r["cell"]), []).append(r["drop"])
    for (sess, cell), drops in sorted(percell.items()):
        total = len(bycell[(sess, cell)])
        kept = sorted(set(range(1, total + 1)) - set(drops))
        flag = "  <-- INTERLEAVED, not a camera that stopped" if len(drops) > 1 else ""
        print("  %-45s %-28s missing %s of %d (keeping %s)%s"
              % (sess, cell, sorted(drops), total, kept, flag))

    print("\nThe Status card records FOUR within-cell omissions; this enumeration finds"
          " %d." % len(within))
    print("Any cell flagged INTERLEAVED above should be checked first: a drop skipped")
    print("because of WHAT IT DID selects the late-time sample on the outcome, and an")
    print("interleaved pattern is the signature that distinguishes that from equipment")
    print("failure. The judgement itself needs the hand annotations.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
