"""Refuse to call the manuscript postable while a placeholder is still in the PIXELS.

This project has been declared "one blocker from arXiv" three times, and each time
the blocker list in the Status card was wrong -- because it was maintained by hand,
against the .tex, while what gets posted is the PDF. The sibling frisbee project lost
a withdrawn claim the same way: `fig4` printed the word "best" into an image, where no
LaTeX macro guard could reach it.

So this checks the COMPILED TEXT, not the source. A placeholder that compiles cleanly
and prints "[INSERT REPOSITORY URL/DOI]" on page 9 is exactly the failure mode a
`latexmk` exit code cannot see.

Run:  python check_ready.py
Exit code 0 only if the PDF carries no placeholder and the build is clean.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pdfplumber

DOC = Path(__file__).resolve().parent / "doc"
PDF = DOC / "chicago_slime_mold.REVIEWED.pdf"
LOG = DOC / "chicago_slime_mold.REVIEWED.log"

# Patterns that must not survive into a posted PDF. Each is paired with what to do
# about it, because a checklist that only says "no" costs a round trip.
PLACEHOLDERS = [
    (r"\[INSERT[^\]]*\]", "fill the bracketed placeholder in the .tex"),
    (r"Prepared for\b", "the author block is still the pre-2026-07-22 placeholder"),
    (r"\bTODO\b|\bFIXME\b|\bXXX\b", "an editing marker reached the PDF"),
    (r"\?\?", r"an unresolved \ref or \cite printed as ??"),
    (r"Lorem ipsum", "filler text"),
]

# Not a placeholder, but not postable either: the affiliation carries no department,
# and this tree names three different ones for the same author. Only the author can
# pick. Detected structurally rather than by pattern.
DEPARTMENTS = ("Department", "Dept.")


def main() -> int:
    if not PDF.exists():
        print(f"FAIL  no PDF at {PDF}; run latexmk -pdf first")
        return 2

    with pdfplumber.open(PDF) as pdf:
        pages = [p.extract_text() or "" for p in pdf.pages]
    text = "\n".join(pages)
    print(f"{PDF.name}: {len(pages)} pages, {len(text)} characters of extracted text")
    if len(text) < 5000:
        print("FAIL  the PDF extracted almost no text; this check would pass vacuously")
        return 2

    bad = 0
    for pat, advice in PLACEHOLDERS:
        for i, page in enumerate(pages, 1):
            for m in re.finditer(pat, page):
                lo = max(0, m.start() - 40)
                print(f"FAIL  p.{i}  {m.group(0)!r}  -- {advice}")
                print(f"        ...{page[lo:m.end() + 40].strip()}...")
                bad += 1

    # The affiliation lives in the first 600 characters of page 1, under the title.
    head = pages[0][:600] if pages else ""
    if not any(d in head for d in DEPARTMENTS):
        print("FAIL  p.1  the affiliation names no department -- and this tree gives three "
              "for the same author (Physics and Astronomy, Mechanical Engineering, "
              "Chemical and Biomolecular Engineering). Only the author can choose.")
        bad += 1

    if LOG.exists():
        log = LOG.read_text(encoding="utf-8", errors="replace")
        for token, label in [("Undefined control sequence", "undefined control sequence"),
                             ("There were undefined references", "undefined reference"),
                             ("Citation", "possible undefined citation")]:
            n = log.count(token)
            if token == "Citation":
                n = log.count("Warning: Citation")
            if n:
                print(f"FAIL  build log: {n} x {label}")
                bad += 1

    if bad == 0:
        print("\nPASS  no placeholder reaches the PDF, and the build log is clean.")
        print("      Remaining blockers are author decisions, not defects:")
        print("        - confirm the author block; it was written by an agent on "
              "2026-07-22, not by the author (see the comment above \author)")
        print("        - the four journal-stage items flagged at .tex lines 27-35")
        return 0
    print(f"\nFAIL  {bad} item(s). The manuscript is not postable.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
