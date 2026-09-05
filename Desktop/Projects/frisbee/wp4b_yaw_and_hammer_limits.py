"""wp4b_yaw_and_hammer_limits.py -- two corrections the earlier passes missed.

(1) YAW STIFFNESS. The manuscript says the absent rim side-force coefficient
    `C_S(alpha)` "removes yaw stiffness, which is one-sided: with it the disc would
    weathercock and shed sideslip faster". CLAUDE_CODE_FRISBEE_HARDENING.md section (a)
    contradicts this, and it is right on two independent grounds, neither of which needs
    a coefficient or a citation:

      * BY SYMMETRY. Yaw stiffness means the aerodynamic yawing moment depends on the
        sideslip angle beta. For a body of revolution the geometry is identical at every
        azimuth, so no such dependence can exist. An alpha-dependent side force C_S(alpha)
        is not a beta-dependent yawing moment, so adding it does not create yaw stiffness.
        The model does not OMIT yaw stiffness; an axisymmetric disc HAS none.

      * BY GYROSCOPICS, which kills it even if a yaw moment did exist. Weathercocking is
        the behaviour of a body that rotates to align with the relative wind. A disc
        released at 90-110 rad/s with I_s = 2.35e-3 kg m^2 carries angular momentum
        L = I_s*omega ~ 0.21-0.26 kg m^2 s^-1, and its response to an applied moment is
        PRECESSION PERPENDICULAR to that moment, not rotation towards it. A rapidly
        spinning disc cannot weathercock.

    So the "one-sided" half of the argument is withdrawn. This does not weaken the
    paragraph's conclusion -- it simplifies it. BOTH missing terms are unsigned, so
    "uncertain in sign, not bounded below" now follows from the symmetry of the problem
    rather than from a race between a one-sided term and a two-sided one.

(2) THE HAMMER LIMITATION SENTENCE, left stale by the item-6 retraction earlier today.
    The limitations section still says the hammer "spends most of its flight outside the
    angle-of-attack range over which the tabulated coefficients exist". Direct execution
    (see wp4b_item6_retract.py) shows alpha stays in [-2.65, +6.41] deg and NEVER leaves
    the tabulated band. The real limitation is the unvalidated release parameterisation.

No sed on .tex (CLAUDE.md); raw strings only.
"""

import re
import sys
from pathlib import Path

TEX = Path(__file__).resolve().parent / "Project1_Frisbee_LoyolaWind_PAPER.REVIEWED.tex"

OLD_LIM_HEAD = r"""\paragraph{Known limitations of the disc model.} The model is axisymmetric and
carries no rim side-force coefficient $C_S(\alpha)$, hence no yaw stiffness."""

NEW_LIM_HEAD = r"""\paragraph{Known limitations of the disc model.} The model is axisymmetric and
carries no rim side-force coefficient $C_S(\alpha)$. It therefore has no yaw
stiffness --- though, as set out below, neither does a real axisymmetric disc."""

OLD_YAW = r"""$C_S(\alpha)$ removes \emph{yaw stiffness}, which is one-sided: with it the disc would
weathercock and shed sideslip faster. But the model also omits the Robins--Magnus side
force, which goes as $\boldsymbol{\omega} \times \mathbf{v}$ and therefore \emph{reverses
with the spin direction} --- a backhand and a forehand spin opposite ways, so the same
crosswind produces a lateral force of opposite sign for the two throws. Its contribution
is handedness- and throw-type-dependent, so it can add to the modelled drift or subtract
from it."""

NEW_YAW = r"""% WP4b item 5 AMENDED [2026-09-04]: the "C_S is one-sided / the disc would weathercock"
% half is WITHDRAWN. An axisymmetric spinning disc has no yaw stiffness BY SYMMETRY, so
% C_S(alpha) cannot supply any; and at 90-110 rad/s the response to a moment is
% precession, not alignment. Raised in CLAUDE_CODE_FRISBEE_HARDENING.md (a); the
% 2026-09-03 pass took the sign argument from it and left this part standing.
$C_S(\alpha)$ was previously described here as removing \emph{yaw stiffness}, one-sidedly,
on the grounds that with it the disc would weathercock and shed sideslip faster.
\textbf{That is withdrawn.} Yaw stiffness means the yawing moment depends on sideslip
angle; for a body of revolution the geometry is identical at every azimuth, so no such
dependence exists. An $\alpha$-dependent side force is not a $\beta$-dependent yawing
moment, and the model does not \emph{omit} yaw stiffness --- an axisymmetric disc has
none. Independently, weathercocking is the behaviour of a body free to rotate into the
relative wind, and a disc spinning at $90$--$110$~rad\,s$^{-1}$ with
$I_s = 2.35\times10^{-3}$~kg\,m$^2$ carries $L \approx 0.21$--$0.26$~kg\,m$^2$\,s$^{-1}$;
its response to an applied moment is precession \emph{perpendicular} to that moment, not
rotation towards it. The model also omits the Robins--Magnus side
force, which goes as $\boldsymbol{\omega} \times \mathbf{v}$ and therefore \emph{reverses
with the spin direction} --- a backhand and a forehand spin opposite ways, so the same
crosswind produces a lateral force of opposite sign for the two throws. Its contribution
is handedness- and throw-type-dependent, so it can add to the modelled drift or subtract
from it. \textbf{Neither missing term is signed}, which is why the conclusion below
follows from the symmetry of the problem rather than from any coefficient."""

OLD_LIMS = r"""The disc model has neither a rim side-force term
$C_S(\alpha)$ nor a Robins--Magnus side force; the first removes yaw stiffness (one-sided)
but the second reverses with spin direction, so \textbf{predicted crosswind drift is
uncertain in sign, not a lower bound} --- see \S\ref{sec:analytical}. Ground is flat and
no tree canopy is included. The hammer throw is exploratory: it spends most of its flight outside
the angle-of-attack range over which the tabulated coefficients exist."""

NEW_LIMS = r"""The disc model has neither a rim side-force term
$C_S(\alpha)$ nor a Robins--Magnus side force. Neither omission is signed --- an
axisymmetric spinning disc has no yaw stiffness to omit, and the Robins--Magnus term
reverses with spin direction --- so \textbf{predicted crosswind drift is
uncertain in sign, not a lower bound} --- see \S\ref{sec:analytical}. Ground is flat and
no tree canopy is included. The hammer throw is exploratory because its release
parameterisation is unvalidated (speed, and a $+80^\circ$ roll that is a guess), not
because its coefficients are extrapolated: its angle of attack stays within
$[-2.7^\circ, +6.4^\circ]$ and never leaves the tabulated range."""


def expand_inputs(src: str, base: Path, _depth: int = 0) -> str:
    if _depth > 5:
        return src

    def sub(m):
        p = (base / m.group(1)).with_suffix(".tex")
        return (expand_inputs(p.read_text(encoding="utf-8"), p.parent, _depth + 1)
                if p.is_file() else "")

    return re.sub(r"\\(?:input|include)\{([^}]+)\}", sub, src)


def main() -> int:
    src = TEX.read_text(encoding="utf-8")
    orig = src
    for old, new, label in ((OLD_LIM_HEAD, NEW_LIM_HEAD, "limitations head"),
                            (OLD_YAW, NEW_YAW, "yaw-stiffness argument"),
                            (OLD_LIMS, NEW_LIMS, "limitations list")):
        if old not in src:
            if "WP4b item 5 AMENDED [2026-09-04]" in src:
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
    dangling = sorted(set(re.findall(r"\\(?:ref|eqref)\{([^}]+)\}", full))
                      - set(re.findall(r"\\label\{([^}]+)\}", full)))
    if dangling:
        print("DANGLING REFERENCES, refusing to write:", dangling)
        return 1

    # the stale hammer claim must not survive anywhere
    stale = [s for s in ("spends most of its flight outside",
                         "spends most of its flight well")
             if s in src]
    if stale:
        print("STALE HAMMER CLAIM STILL PRESENT, refusing to write:", stale)
        return 1

    TEX.write_text(src, encoding="utf-8")
    print("wrote %s (%d -> %d chars)" % (TEX.name, len(orig), len(src)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
