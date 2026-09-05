"""wp4b_item4_roughness.py -- withdraw the "sigma_u genuinely rose" half of the
"shelter is not protection" claim, which rests on an unmatched control.

WHY
---
Both surviving statements of the claim assert that the empty-fetch case is a
counterfactual with "the identical inflow profile over identical ground". It is not.
Verified directly from case_info.json (check_control_matching.py, 2026-09-04):

    case             z0_inlet   z0_ground     u_star   k_inlet_ref   matched?
    dir_090_empty *     0.001        0.03     0.2226        0.1651   (control)
    dir_090             0.001        0.03     0.2226        0.1651   YES
    dir_210             0.5          0.03     0.6733        1.5113   NO
    dir_270             0.5          0.03     0.6733        1.5113   NO

The GROUND roughness is matched; the INLET roughness is not. dir_210 and dir_270 are
fed inflow carrying sqrt(1.5113/0.1651) = 3.03x the fluctuation of the control before
the flow reaches a single building. Dividing that out:

    heading   memo sigma_u ratio   approach-flow ratio   buildings' own effect
      90              0.91                1.00                 0.91   (clean)
     210              1.51                3.03                 0.50
     270              1.57                3.03                 0.52

So against MATCHED approach flow the buildings REDUCE sigma_u in all three directions,
and the "raised by half again at 210 and 270" claim is entirely the inlet condition.
The only roughness-clean heading is 90 deg -- which is the one where the paper already
concedes the claim does not hold.

WHAT SURVIVES, AND IT IS THE MORE INTERESTING HALF
--------------------------------------------------
The TI argument is untouched and is roughness-clean, because it lives entirely inside
the 90 deg column: TI rises 2.08x while sigma_u FALLS to 0.91, and sigma-ratio over
U-ratio (0.91/0.44 = 2.07) reproduces the measured TI ratio (0.25/0.12 = 2.08) to 0.5%.
That is an arithmetic identity within a single matched pair of cases and it does not
depend on the control being matched to the other two headings.

This was flagged in CLAUDE_CODE_FRISBEE_HARDENING.md section C2 and, as the physics
referee noted on 2026-09-04, the 2026-09-03 retraction pass absorbed the TI half of that
finding and left the roughness half unaddressed.

No sed on .tex (CLAUDE.md); raw strings only.
"""

import re
import sys
from pathlib import Path

TEX = Path(__file__).resolve().parent / "Project1_Frisbee_LoyolaWind_PAPER.REVIEWED.tex"

OLD_A = r"""\paragraph{Shelter is not protection --- when it raises $\sigma_u$.}
The empty-fetch case of \S\ref{sec:cfdvalid} doubles as a no-buildings
counterfactual: the identical inflow profile over identical ground, without the
campus. Comparing the two at $2$~m over the pitch separates what the buildings
took away from what they added. They removed mean wind, cutting it to
\VertUratioRange{} of the unobstructed value. Measured as
$\sigma_u = \sqrt{2k/3}$, they left the fluctuation at \VertSigRatioRange{} of the
control: \emph{raised} by about half at $210^\circ$ and $270^\circ$
(\VertSigRatioHTwoOneZero, \VertSigRatioHTwoSevenZero), and \emph{lowered} at
$90^\circ$ (\VertSigRatioHNineZero)."""

NEW_A = r"""\paragraph{Turbulence intensity is the wrong diagnostic for a thrown disc.}
The empty-fetch case of \S\ref{sec:cfdvalid} doubles as a no-buildings
counterfactual, and it is a valid one \emph{for the $90^\circ$ heading only}. The
control and the $90^\circ$ production case share an inlet roughness of
$z_0 = 0.001$~m; the $210^\circ$ and $270^\circ$ cases were run at $z_0 = 0.5$~m and
so meet inflow already carrying $\sqrt{1.5113/0.1651} = 3.03$ times the control's
fluctuation before reaching a building. Ground roughness is matched across all four
($z_0 = 0.03$~m); the inlet is not.

% WP4b item 4 AMENDED [2026-09-04]: the "sigma_u genuinely rose at 210/270" half is
% WITHDRAWN. Verified from case_info.json by check_control_matching.py -- see the
% \emph{Withdrawn} note below. CLAUDE_CODE_FRISBEE_HARDENING.md C2 raised this and the
% 2026-09-03 pass absorbed only the TI half of it.
\emph{\textbf{Withdrawn.}} An earlier version of this paragraph reported that the
buildings \emph{raised} $\sigma_u$ by about half at $210^\circ$ and $270^\circ$
(ratios \VertSigRatioHTwoOneZero{} and \VertSigRatioHTwoSevenZero{}) and lowered it
at $90^\circ$ (\VertSigRatioHNineZero). Those two rises are the inlet condition, not
the buildings: dividing out the $3.03$ approach-flow ratio leaves $0.50$ and $0.52$,
i.e.\ against matched inflow the campus \emph{reduces} $\sigma_u$ in all three
directions. The comparison is only roughness-clean at $90^\circ$, where the ratio is
\VertSigRatioHNineZero{}. \textbf{No claim that shelter raises $\sigma_u$ can be made
from these cases}; testing it needs an empty-domain control run at
$z_0^{\text{inlet}} = 0.5$~m, which is one additional case.

What follows does not depend on that, because it lives entirely inside the matched
$90^\circ$ pair."""

OLD_A_TAIL = r"""manufactures a factor of two where the physical disturbance the disc responds to has, in
that direction, slightly \emph{decreased}. A disc has never heard of $U$; the group it
responds to is $\sigma_u/V_d$. \textbf{So the claim holds at $210^\circ$ and $270^\circ$,
where $\sigma_u$ genuinely rose, and does not hold at $90^\circ$ --- and no version of it
should be argued from TI.}"""

NEW_A_TAIL = r"""manufactures a factor of two where the physical disturbance the disc responds to has, in
that direction, slightly \emph{decreased}. A disc has never heard of $U$; the group it
responds to is $\sigma_u/V_d$. \textbf{This is the result that survives: turbulence
intensity cannot tell you whether shelter helped or hurt a thrown object, because a
collapsing mean inflates it either way.} Whether shelter ever \emph{raises} $\sigma_u$ is
left open here --- the two headings that appeared to show it are not roughness-matched to
the control, and the one that is shows $\sigma_u$ falling."""

OLD_B = r"""The counterfactual of \S\ref{sec:vertical} sharpens the point above into
something that should hold beyond this site, provided it is stated with the
condition attached. The buildings cut mean wind over the pitch to
\VertUratioRange{} of its unobstructed value while leaving $\sigma_u$ at
\VertSigRatioRange{} of the control. In two of three directions they raised
$\sigma_u$ outright, and there shelter genuinely converted a bias a thrower can aim
off into a variance they cannot. In the third ($90^\circ$) they lowered it, and the
apparent penalty is an artefact of normalising by the collapsed mean.
\textbf{The transferable claim is therefore conditional: shelter is not protection
\emph{when it raises $\sigma_u$}, and turbulence intensity cannot tell you whether it
did.}"""

NEW_B = r"""The counterfactual of \S\ref{sec:vertical} sharpens the point above, within the
limit set there: only the $90^\circ$ pair is roughness-matched, so only it can be read
as a buildings-versus-no-buildings comparison. In that pair the buildings cut mean wind
over the pitch to \VertUratioHNineZero{} of its unobstructed value while leaving
$\sigma_u$ at \VertSigRatioHNineZero{} --- both down, yet turbulence intensity rises by
a factor of $2.08$, purely because the mean in its denominator collapsed faster than
the fluctuation.

\textbf{The transferable claim is therefore about the diagnostic, not about shelter:
turbulence intensity cannot tell you whether shelter helped or hurt a thrown object, and
the quantity that can is $\sigma_u/V_d$.} Whether shelter is ever \emph{anti}-protective
--- whether a building array can raise $\sigma_u$ outright over the ground it shelters
--- is a well-posed question this study does not answer, because the two cases that
appeared to answer it differ from the control in inlet roughness as well as in
buildings. It is worth one matched case to settle."""


def expand_inputs(src: str, base: Path, _depth: int = 0) -> str:
    if _depth > 5:
        return src

    def sub(m):
        p = (base / m.group(1)).with_suffix(".tex")
        if not p.is_file():
            return ""
        return expand_inputs(p.read_text(encoding="utf-8"), p.parent, _depth + 1)

    return re.sub(r"\\(?:input|include)\{([^}]+)\}", sub, src)


def main() -> int:
    src = TEX.read_text(encoding="utf-8")
    orig = src
    patches = [(OLD_A, NEW_A, "shelter paragraph head"),
               (OLD_A_TAIL, NEW_A_TAIL, "shelter paragraph tail"),
               (OLD_B, NEW_B, "transferable-claim paragraph")]
    for old, new, label in patches:
        if old not in src:
            if "WP4b item 4 AMENDED [2026-09-04]" in src:
                print("already applied: %s" % label)
                continue
            print("ANCHOR NOT FOUND (%s) -- refusing to write" % label)
            return 1
        src = src.replace(old, new, 1)
        print("patched: %s" % label)

    if src == orig:
        print("no change")
        return 0

    full = expand_inputs(src, TEX.parent)
    labels = set(re.findall(r"\\label\{([^}]+)\}", full))
    refs = set(re.findall(r"\\(?:ref|eqref)\{([^}]+)\}", full))
    dangling = sorted(refs - labels)
    if dangling:
        print("DANGLING REFERENCES, refusing to write:", dangling)
        return 1

    # every \Vert* macro used must be defined in the generated file
    used = set(re.findall(r"\\(Vert[A-Za-z]+)", src))
    defined = set(re.findall(r"\\newcommand\{\\([A-Za-z]+)\}", full))
    missing = sorted(m for m in used if m not in defined)
    if missing:
        print("UNDEFINED MACROS, refusing to write:", missing)
        return 1

    TEX.write_text(src, encoding="utf-8")
    print("wrote %s (%d -> %d chars); %d Vert macros used, all defined"
          % (TEX.name, len(orig), len(src), len(used)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
