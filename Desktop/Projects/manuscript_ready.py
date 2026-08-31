"""Refuse to call a manuscript postable while a placeholder is still in the PIXELS.

Three projects in this tree have now been declared "one blocker from posting" while a
literal placeholder sat in the compiled PDF, and each time the blocker list was kept by
hand against the .tex while what gets posted is the PDF:

  * slime mold: `[INSERT REPOSITORY URL/DOI]` on page 8 of the REVIEWED build, since
    July, while the Status card said the author line was the only thing left.
  * frisbee: `fig4_throwability` printed the word "best" INTO THE IMAGE, where no LaTeX
    macro guard could reach it, and a withdrawn claim survived three corrections.
  * the drop-impact PRL: a caption describing a panel (b) that does not exist.

So this reads the compiled text, not the source. A `latexmk` exit code cannot see any of
the above.

It also answers the second question these directories raise, which is *which file*.
`bidet/` alone holds a dozen `.tex` variants -- `_1`, `.EDITED`, `.REVIEWED`,
`.JOURNAL`, `_pre-review` -- and the slime-mold Status card pointed at the wrong one for
six weeks. The report lists every sibling with its date and page count so the live one
has to be chosen deliberately.

Usage:
    python manuscript_ready.py <file.pdf> [<file.pdf> ...]
    python manuscript_ready.py --dir <project-dir>     # every PDF with a sibling .tex

Exit code 0 only if every PDF checked is clean.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import re
import sys
from pathlib import Path

import pdfplumber

# Each pattern is paired with what to do about it: a checklist that only says "no"
# costs a round trip.
PLACEHOLDERS: list[tuple[str, str]] = [
    (r"\[INSERT[^\]]*\]", "fill the bracketed placeholder in the .tex"),
    (r"\[TODO[^\]]*\]|\bTODO\b|\bFIXME\b", "an editing marker reached the PDF"),
    (r"\bXXX+\b", "an editing marker reached the PDF"),
    (r"Prepared for\b", "the author block is still a placeholder"),
    (r"\bTBD\b|\bTBA\b", "an unresolved value"),
    (r"Lorem ipsum", "filler text"),
    (r"\?\?(?!\?)", r"an unresolved \ref or \cite printed as ??"),
    (r"\bCITE\b|\\cite\{", "an unexpanded citation command"),
    (r"MISSINGNUM", "a generated macro fell back to its missing-value sentinel"),
    (r"your name here|Author Name|First Last", "a template author line"),
]

# Words that are fine in a body but are a red flag inside a Data Availability or
# acknowledgement sentence, where a real URL or DOI belongs.
AVAILABILITY = re.compile(
    r"(availab\w+|deposit\w+|repositor\w+)[^.]{0,200}\.", re.IGNORECASE)
LOOKS_LIKE_LINK = re.compile(r"https?://|doi\.org|10\.\d{4,}/|zenodo|github", re.IGNORECASE)


def check(pdf_path: Path) -> int:
    if not pdf_path.exists():
        print(f"FAIL  {pdf_path}: no such file")
        return 1
    with pdfplumber.open(pdf_path) as pdf:
        pages = [p.extract_text() or "" for p in pdf.pages]
    text = "\n".join(pages)
    stamp = _dt.datetime.fromtimestamp(pdf_path.stat().st_mtime).strftime("%Y-%m-%d")
    print(f"\n{pdf_path.name}  --  {len(pages)} pages, {len(text)} chars, built {stamp}")
    if len(text) < 2000:
        print("  FAIL  almost no extractable text; this check would pass vacuously "
              "(scanned or image-only?)")
        return 1

    bad = 0
    for pat, advice in PLACEHOLDERS:
        for i, page in enumerate(pages, 1):
            for m in re.finditer(pat, page):
                snippet = " ".join(page[max(0, m.start() - 50):m.end() + 50].split())
                print(f"  FAIL  p.{i}  {m.group(0)!r}  -- {advice}")
                print(f"          ...{snippet}...")
                bad += 1

    # A data-availability sentence that names no link is not a placeholder, but it is
    # not a statement either.
    for i, page in enumerate(pages, 1):
        for m in AVAILABILITY.finditer(page):
            sent = " ".join(m.group(0).split())
            # Word boundaries matter here: "metadata" and "available compute" between
            # them produced a false positive on the frisbee manuscript's p.9.
            if re.search(r"\b(data|code|scripts?|dataset)\b", sent, re.IGNORECASE):
                if not LOOKS_LIKE_LINK.search(sent):
                    print(f"  WARN  p.{i}  availability sentence names no URL or DOI:")
                    print(f"          {sent[:180]}")
    return 1 if bad else 0


def siblings(pdf_path: Path) -> None:
    """List the other manuscript variants next to this one, with dates."""
    stem = re.split(r"[._]", pdf_path.stem)[0]
    others = sorted(p for p in pdf_path.parent.glob(f"{stem}*.tex"))
    if len(others) <= 1:
        return
    print(f"  {len(others)} .tex variants share the stem {stem!r} -- "
          "pick the live one deliberately:")
    for o in others:
        d = _dt.datetime.fromtimestamp(o.stat().st_mtime).strftime("%Y-%m-%d")
        print(f"      {d}  {o.name}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdfs", nargs="*", type=Path)
    ap.add_argument("--dir", type=Path, help="check every PDF that has a sibling .tex")
    a = ap.parse_args()

    targets: list[Path] = list(a.pdfs)
    if a.dir:
        targets += [p for p in sorted(a.dir.rglob("*.pdf"))
                    if p.with_suffix(".tex").exists()]
    if not targets:
        ap.error("nothing to check")

    rc = 0
    for p in targets:
        r = check(p)
        siblings(p)
        rc |= r
    print("\n" + ("PASS  no placeholder reaches any PDF checked." if rc == 0
                  else "FAIL  at least one PDF is not postable."))
    return rc


if __name__ == "__main__":
    sys.exit(main())
