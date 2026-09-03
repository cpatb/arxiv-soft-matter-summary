r"""wp4b_excise.py -- WP4b item 2: move the planning memo out of the paper.

Run once, from the frisbee project root:

    python wp4b_excise.py

WHY A SCRIPT AND NOT A ONE-LINER
--------------------------------
`~/.claude/rules/latex-writing.md`: never edit .tex with sed/perl one-liners, because
backslash escapes get mangled. That is not hypothetical here -- the first attempt at this
edit went through a bash heredoc and Python turned `\S\ref{sec:limits}` into
`\Sef{sec:limits}`, because `\r` inside a non-raw string literal is a carriage return.
Every literal in this file is a RAW string for that reason, and the file is written to
disk rather than piped.

WHAT IT DOES
------------
Removes `\section{Feasibility check}` and `\section{Project plan}` -- the 14-week work
breakdown, the risks-and-mitigations table, all of it -- into
PLANNING_MEMO_EXCISED_2026-09-03.tex, rewrites the "How to read this document" block that
described them, and repairs the four cross-references that pointed into them.

The brief's reason is right: a document that opens with its own Gantt chart is a project
report, not a paper. The earlier defence -- a status box declaring the paper sections
authoritative -- does not survive a referee, because the planning text is still in the PDF
and still says things the executed work contradicts (thirteen CFD directions, different
building heights, a different site description).

It is idempotent: re-running after a successful pass exits without changing anything.
"""

from __future__ import annotations

import io
import sys
from pathlib import Path

PAPER = Path("Project1_Frisbee_LoyolaWind_PAPER.REVIEWED.tex")
OUT = Path("PLANNING_MEMO_EXCISED_2026-09-03.tex")

EXCISED_HEADER = r"""% PLANNING_MEMO_EXCISED_2026-09-03.tex
%
% WP4b item 2. The feasibility memo and the project plan, including the 14-week work
% breakdown and the risks-and-mitigations table, removed from
% Project1_Frisbee_LoyolaWind_PAPER.REVIEWED.tex on 2026-09-03 and kept here rather than
% deleted.
%
% ONE THING WAS LOAD-BEARING AND HAS BEEN CARRIED FORWARD RATHER THAN EXCISED WITH THE
% REST: the GATE G5 correction of 2026-09-02, which replaced the feasibility table's claim
% that coupling a CFD field into a disc trajectory is "not done in published literature".
% That correction now lives in the Introduction, where it belongs, and in
% GATE_G5_PRIOR_ART_2026-09-02.md. If this file is ever re-read, do NOT restore the
% table's original wording -- it is prior art twice over.
%
% Nothing else here is cited by the paper. The cross-references that pointed into it
% (sec:feas, sec:plan, sec:feasibility-verdict) were rewritten as prose in the same commit.

"""

NEW_HOWTO = r"""% WP4b item 2 [2026-09-03]: this block described a three-layer document whose first two
% layers are now in PLANNING_MEMO_EXCISED_2026-09-03.tex. Rewritten to describe what is
% actually here; the "status of this version" box went with the sections it was hedging.
\section*{How to read this document}

This is a draft paper, not a project report. \S\ref{sec:limits} states the delivered scope
explicitly and should be read alongside the results: several things the study set out to do
were not done, and it says which. Every number quoted comes from a converged run or a
passing test; nothing is projected.

"""

OLD_CONTRIB1 = (r"The coupling itself is established for other bodies (see "
                r"\S\ref{sec:feasibility-verdict}); the disc and the distribution are "
                r"what is new.")
NEW_CONTRIB1 = (r"The coupling itself is established for other bodies --- Kakimpa "
                r"\emph{et al.} (2012) for windborne debris, and US 11{,}052{,}285\,B2 "
                r"(2016) for a ball in a stadium --- and neither the disc nor the "
                r"distribution is.")

OLD_CONTRIB3 = (r"  \item A case study identifying the lakefront internal boundary layer "
                r"as the dominant local control on Ultimate-disc flight at the chosen site.")
NEW_CONTRIB3 = r"""% WP4b item 1 [2026-09-03]: RETRACTED. This contribution named the lakefront internal
% boundary layer as the dominant control -- the premise the study FALSIFIED -- and so
% contradicted the abstract and the conclusions in the same PDF. Replaced by what was found.
  \item A case study in which that premise \emph{fails}: the study was designed around the
  lakefront internal boundary layer as the dominant local control, and finds instead that a
  45\,m building nine metres off the sideline sets the flow in the disc band. The negative
  result is the site-specific contribution, and it is a magnitude statement both meshes
  agree on."""

OLD_PLANREF = (r"thirteen-direction, $15^\circ$-increment sweep described in the project plan"
               "\n" r"(\S\ref{sec:plan}) \textbf{was not run}.")
NEW_PLANREF = (r"thirteen-direction, $15^\circ$-increment sweep described in the original"
               "\n" r"planning memo (excised from this document; see"
               "\n" r"\texttt{PLANNING\_MEMO\_EXCISED\_2026-09-03.tex}) \textbf{was not run}.")


def main() -> int:
    s = io.open(PAPER, encoding="utf-8").read()
    if r"\section{Feasibility check}" not in s:
        print("already excised -- nothing to do")
        return 0

    a = s.index(r"\section{Feasibility check}")
    b = s.index(r"\section{Abstract (draft)}")
    io.open(OUT, "w", encoding="utf-8").write(EXCISED_HEADER + s[a:b])
    print("wrote %s (%d chars)" % (OUT, b - a))
    s = s[:a] + s[b:]

    how_a = s.index(r"\section*{How to read this document}")
    how_b = s.index(r"\section{Abstract (draft)}")
    s = s[:how_a] + NEW_HOWTO + s[how_b:]

    for old, new, label in ((OLD_CONTRIB1, NEW_CONTRIB1, "contribution 1 dead ref"),
                            (OLD_CONTRIB3, NEW_CONTRIB3, "contribution 3 retraction"),
                            (OLD_PLANREF, NEW_PLANREF, "limitations sec:plan ref")):
        if old not in s:
            print("FAILED to find: %s" % label)
            return 1
        s = s.replace(old, new, 1)
        print("ok: %s" % label)

    for dangling in ("sec:feas", "sec:plan", "sec:feasibility-verdict"):
        if dangling in s:
            print("STILL DANGLING: %s -- refusing to write" % dangling)
            return 1

    io.open(PAPER, "w", encoding="utf-8").write(s)
    print("paper now %d chars" % len(s))
    return 0


if __name__ == "__main__":
    sys.exit(main())
