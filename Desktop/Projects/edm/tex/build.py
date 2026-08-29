#!/usr/bin/env python3
"""
Build the LaTeX edition of the EDM / Chicago house guide.

The markdown files in ../THEORY, ../FL_STUDIO, ../PRACTICE and ../LISTENING are
the single source of truth. This script converts each one to a LaTeX body
fragment with pandoc, then emits two products from those same fragments:

  tex/bodies/<slug>.tex   the converted body, \\section-level, no preamble
  tex/<slug>.tex          a standalone article that \\inputs its body
  tex/edm_guide.tex       the bound volume, book class, one \\chapter per source

Both products read the same fragment, so the prose exists once on disk. Editing
a markdown file and re-running this script updates the standalone chapter and
the bound volume together; there is no second copy to fall out of step.

Requires: pandoc (tested on 3.10), XeLaTeX for the subsequent compile.

Usage:
    python build.py            convert only
    python build.py --pdf      convert, then compile the bound volume
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

TEX = Path(__file__).resolve().parent
ROOT = TEX.parent
BODIES = TEX / "bodies"

AUTHOR = "Compiled for C. Burton"
DATE = "August 2026"
RUNNING_TITLE = "EDM \\& Chicago House"

# (source path relative to ROOT, slug, part).  Order is reading order, which is
# also the order of the bound volume.  A file absent from this list is not built,
# so adding a chapter means adding a line here.
MANIFEST = [
    ("README.md",                        "00_how_to_use",     None),

    ("THEORY/01_what_edm_is.md",         "th01_what_edm_is",  "Theory"),
    ("THEORY/02_chicago_house.md",       "th02_chicago",      "Theory"),
    ("THEORY/03_rhythm.md",              "th03_rhythm",       "Theory"),
    ("THEORY/04_harmony.md",             "th04_harmony",      "Theory"),
    ("THEORY/05_sound_design.md",        "th05_sound_design", "Theory"),
    ("THEORY/06_arrangement.md",         "th06_arrangement",  "Theory"),
    ("THEORY/07_mixing.md",              "th07_mixing",       "Theory"),
    ("THEORY/08_genre_map.md",           "th08_genre_map",    "Theory"),
    ("THEORY/09_glossary.md",            "th09_glossary",     "Theory"),

    ("FL_STUDIO/00_setup.md",            "fl00_setup",        "FL Studio"),
    ("FL_STUDIO/01_interface.md",        "fl01_interface",    "FL Studio"),
    ("FL_STUDIO/02_first_beat.md",       "fl02_first_beat",   "FL Studio"),
    ("FL_STUDIO/03_first_track.md",      "fl03_first_track",  "FL Studio"),
    ("FL_STUDIO/04_shortcuts.md",        "fl04_shortcuts",    "FL Studio"),
    ("FL_STUDIO/05_plugins.md",          "fl05_plugins",      "FL Studio"),
    ("FL_STUDIO/06_troubleshooting.md",  "fl06_trouble",      "FL Studio"),

    ("PRACTICE/8_week_plan.md",          "pr01_eight_weeks",  "Practice"),
    ("PRACTICE/exercises.md",            "pr02_exercises",    "Practice"),

    ("LISTENING/chicago_canon.md",       "li01_canon",        "Listening"),
    ("LISTENING/chicago_scene_today.md", "li02_scene",        "Listening"),
]

# A markdown blockquote whose first bold run matches one of these leads becomes
# the corresponding callout box.  Keeping the trigger in the markdown means the
# plain-text edition shows the same emphasis, just without the frame.
CALLOUTS = [
    (("why this works", "the idea", "the point", "in practice"), "intuitionbox"),
    (("watch out", "the trap", "common mistake", "the catch"),   "pitfallbox"),
    (("in short", "summary", "what you have"),                   "reviewbox"),
    (("listen for", "go and hear it", "hear it"),                "listenbox"),
]

QUOTE_RE = re.compile(r"\\begin\{quote\}(.*?)\\end\{quote\}", re.S)
LEAD_RE = re.compile(r"\s*\\textbf\{([^}]*)\}")

LONGTABLE_RE = re.compile(
    r"(\\begin\{longtable\}\[\]\{@\{\})([lcr]+)(@\{\}\})(.*?)(\\end\{longtable\})", re.S)

# A row of source text wider than this many characters will not fit the text
# block at 11pt, so its table needs wrapping columns rather than natural ones.
NATURAL_WIDTH_LIMIT = 74
MIN_COLUMN_FRACTION = 0.07


def split_title(md: str) -> tuple[str, str]:
    """Peel the leading H1 off a markdown file and return (title, remainder)."""
    lines = md.splitlines()
    for i, line in enumerate(lines):
        if line.startswith("# "):
            return line[2:].strip(), "\n".join(lines[i + 1:]).lstrip("\n")
    raise ValueError("no H1 title found")


def apply_callouts(tex: str) -> tuple[str, int]:
    """Promote tagged blockquotes to mdframed callout boxes."""
    count = 0

    def repl(m: re.Match) -> str:
        nonlocal count
        body = m.group(1)
        lead = LEAD_RE.match(body)
        if not lead:
            return m.group(0)
        key = lead.group(1).strip().rstrip(".").lower()
        for keys, env in CALLOUTS:
            if key in keys:
                count += 1
                return f"\\begin{{{env}}}\n{body.strip()}\n\\end{{{env}}}"
        return m.group(0)

    return QUOTE_RE.sub(repl, tex), count


def prefix_anchors(tex: str, slug: str) -> str:
    """Namespace pandoc's heading anchors by chapter.

    pandoc derives \\hypertarget and \\label names from the heading text, so
    "Sampling" in the sound-design chapter and "Sampling" in the glossary
    collide in the bound volume. Prefixing with the slug makes them unique
    without giving up the ability to \\ref them later.
    """
    tex = re.sub(r"\\hypertarget\{([^}]*)\}",
                 lambda m: f"\\hypertarget{{{slug}-{m.group(1)}}}", tex)
    return re.sub(r"\\label\{([^}]*)\}",
                  lambda m: f"\\label{{{slug}-{m.group(1)}}}", tex)


def _visible_len(cell: str) -> int:
    """Approximate the printed width of a table cell in characters."""
    s = re.sub(r"\\[a-zA-Z]+\s*", "", cell)
    return len(s.replace("{", "").replace("}", "").replace("\\", "").strip())


def fix_longtable_widths(tex: str) -> tuple[str, int]:
    """Give wide tables wrapping columns sized to their content.

    pandoc emits `l` columns for pipe tables, which never wrap. That is fine
    for a table of note frequencies and disastrous for one whose last column is
    a sentence: the row simply runs off the page. Any table whose widest row
    exceeds the text block is converted to \\raggedright p-columns with widths
    proportional to the longest cell in each column.
    """
    count = 0

    def repl(m: re.Match) -> str:
        nonlocal count
        pre, spec, post, inner, end = m.groups()
        ncol = len(spec)
        if ncol < 2:
            return m.group(0)

        rows = []
        for line in inner.splitlines():
            line = line.strip()
            if not line.endswith("\\\\"):
                continue
            cells = re.split(r"(?<!\\)&", line[:-2])
            if len(cells) == ncol:
                rows.append([_visible_len(c) for c in cells])
        if not rows:
            return m.group(0)

        widths = [max(r[i] for r in rows) for i in range(ncol)]
        if sum(widths) + 3 * ncol <= NATURAL_WIDTH_LIMIT:
            return m.group(0)          # fits as-is; natural columns look better

        total = sum(widths)
        frac = [max(w / total, MIN_COLUMN_FRACTION) for w in widths]
        norm = sum(frac)
        frac = [f / norm for f in frac]

        # Built by concatenation rather than an f-string: the target text is
        # brace-dense LaTeX and doubling every brace makes it unreadable.
        newspec = "".join(
            r">{\raggedright\arraybackslash}p{(\linewidth - "
            + str(2 * ncol) + r"\tabcolsep) * \real{" + "%.4f" % f + "}}"
            for f in frac)
        count += 1
        return pre + newspec + post + inner + end

    return LONGTABLE_RE.sub(repl, tex), count


def latex_escape(s: str) -> str:
    """Escape a plain-text title for use inside a LaTeX argument."""
    for a, b in (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"),
                 ("$", r"\$"), ("#", r"\#"), ("_", r"\_"), ("{", r"\{"),
                 ("}", r"\}"), ("~", r"\textasciitilde{}"),
                 ("^", r"\textasciicircum{}")):
        s = s.replace(a, b)
    return s


def convert(src: Path, slug: str) -> tuple[str, str]:
    """Return (title, latex body) for one markdown source."""
    md = src.read_text(encoding="utf-8")
    title, rest = split_title(md)

    proc = subprocess.run(
        ["pandoc", "-f", "gfm", "-t", "latex",
         "--no-highlight", "--wrap=preserve"],
        input=rest, capture_output=True, text=True, encoding="utf-8",
    )
    if proc.returncode != 0:
        raise RuntimeError(f"pandoc failed on {src}:\n{proc.stderr}")

    body, n_call = apply_callouts(proc.stdout)
    body = prefix_anchors(body, slug)
    body, n_tab = fix_longtable_widths(body)

    notes = []
    if n_call:
        notes.append(f"{n_call} callout{'s' if n_call != 1 else ''}")
    if n_tab:
        notes.append(f"{n_tab} table{'s' if n_tab != 1 else ''} rewidened")
    if notes:
        print("    " + ", ".join(notes))
    return title, body


STANDALONE = """\
% ---------------------------------------------------------------------------
%  {title}
%
%  Standalone edition of {relsrc}.
%  Generated by tex/build.py -- edit the markdown, not this file.
%  Engine: XeLaTeX.
% ---------------------------------------------------------------------------
\\documentclass[11pt,letterpaper]{{article}}
\\input{{preamble}}

\\title{{\\vspace{{-2.5em}}\\textbf{{{title}}}}}
\\author{{\\normalsize {author}}}
\\date{{\\normalsize {date}}}

\\begin{{document}}
\\maketitle
\\vspace{{-1.5em}}
\\markboth{{{title}}}{{{title}}}

\\input{{bodies/{slug}}}

\\end{{document}}
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", action="store_true",
                    help="compile the bound volume after converting")
    args = ap.parse_args()

    BODIES.mkdir(exist_ok=True)
    built = []

    for rel, slug, part in MANIFEST:
        src = ROOT / rel
        if not src.exists():
            print(f"  MISSING {rel}", file=sys.stderr)
            return 1
        print(f"  {rel}")
        title, body = convert(src, slug)

        (BODIES / f"{slug}.tex").write_text(body, encoding="utf-8")
        (TEX / f"{slug}.tex").write_text(
            STANDALONE.format(title=latex_escape(title), relsrc=rel,
                              slug=slug, author=AUTHOR, date=DATE),
            encoding="utf-8")
        built.append((slug, title, part))

    # --- the bound volume --------------------------------------------------
    out = [
        "% " + "-" * 73,
        "%  Electronic Dance Music and Chicago House -- bound volume",
        "%",
        "%  Generated by tex/build.py from the markdown sources one directory up.",
        "%  Edit the markdown, not this file. Engine: XeLaTeX.",
        "% " + "-" * 73,
        "\\documentclass[11pt,letterpaper,oneside]{book}",
        "\\input{preamble}",
        "",
        "\\title{\\textbf{Electronic Dance Music and Chicago House}\\\\[0.35em]",
        "  \\large A Standing Start: Theory, FL Studio, and the Local Tradition}",
        f"\\author{{\\normalsize {AUTHOR}}}",
        f"\\date{{\\normalsize {DATE}}}",
        "",
        "\\begin{document}",
        "\\maketitle",
        "\\tableofcontents",
        "",
    ]
    current_part = "\x00"
    for slug, title, part in built:
        if part != current_part and part is not None:
            out.append(f"\\part{{{latex_escape(part)}}}")
            current_part = part
        out.append(f"\\chapter{{{latex_escape(title)}}}")
        out.append(f"\\input{{bodies/{slug}}}")
        out.append("")
    out.append("\\end{document}")
    (TEX / "edm_guide.tex").write_text("\n".join(out) + "\n", encoding="utf-8")

    print(f"\n{len(built)} chapters -> {len(built)} standalone .tex "
          f"+ bodies/ + edm_guide.tex")

    if args.pdf:
        print("\ncompiling edm_guide.tex ...")
        r = subprocess.run(
            ["latexmk", "-xelatex", "-interaction=nonstopmode",
             "-halt-on-error", "edm_guide.tex"],
            cwd=TEX, capture_output=True, text=True, encoding="utf-8",
            errors="replace")
        if r.returncode != 0:
            tail = "\n".join((r.stdout or "").splitlines()[-40:])
            print(tail, file=sys.stderr)
            print("COMPILE FAILED", file=sys.stderr)
            return 1
        print("edm_guide.pdf written")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
