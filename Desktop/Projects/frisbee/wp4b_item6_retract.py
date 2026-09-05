"""wp4b_item6_retract.py -- retract the hammer "clamp artefact" claim.

WHY
---
WP4b item 6 (written 2026-09-03) asserted that the model's 12.1 m still-air hammer
range is "a clamp artefact": that the Potts & Crowther table is edge-clamped at
+30 deg, that a hammer released at alpha_0 = +40 deg therefore "generates no lift at
all", and that 12.1 m reproduces v^2/g = 12.3 m to 2%, so the disc is an unpowered
point mass.

A physics referee challenged it and direct execution refutes every link (2026-09-04):

  alpha over the whole flight        [-2.65, +6.41] deg   -- never leaves the table
  fraction of samples outside band   0.0000
  CL                                 0.0157 .. 0.5313     -- lift IS produced
  peak lift  qbarA*CL                2.254 N  vs  mg = 1.714 N

and the ballistic comparison was against the wrong formula. v^2/g ignores BOTH the
40 deg launch elevation and the 2.0 m launch height. Integrating the same initial
condition with no aerodynamics at all gives 14.208 m, so the 6-DOF hammer flies 15%
SHORT of ballistic (drag), not equal to it. The "agreement to 2%" was a coincidence
between a model result and a formula that does not describe this launch.

The real mechanism is ROLL, and it is monotone (verified by sweep):

    roll_deg   range_m   alpha range
       0.0      2.063    [-0.60, +72.86]   <- THIS is where the table is exceeded
      20.0      5.542    [-1.14, +37.16]
      40.0      8.692    [-1.67, +29.50]
      60.0     10.835    [-2.20, +18.75]
      80.0     12.067    [-2.65,  +6.41]   <- the configured hammer
      90.0     12.356    [-3.27,  +0.00]

`trajectories.hammer()` sets roll_deg = +80, which lays the disc almost on its side.
The lift vector is then nearly horizontal, so it does little work against gravity --
the disc is not lift-LESS, its lift is mis-DIRECTED. Range falls to near-ballistic
because of release geometry, which is a modelling choice, not a table artefact.

Consequence for the author-decisions section: the PREVIOUS default ("retain but report
as exploratory") was withdrawn on the strength of the artefact claim. Since that claim
is refuted, the withdrawal is itself withdrawn and the default is restored -- but its
stated ground ("spends most of its flight outside the band") is ALSO false as the
hammer is actually parameterised, and is corrected here. Option (b)'s "extend the
table past stall" is likewise not the fix, because the table is never exceeded.

Raw strings throughout; no sed on .tex (CLAUDE.md).
"""
import re
import sys
from pathlib import Path

TEX = Path(__file__).resolve().parent / "Project1_Frisbee_LoyolaWind_PAPER.REVIEWED.tex"

OLD_RESULTS = r"""    % WP4b item 6 [2026-09-03]: and the 12.1 m is itself an artefact, which is worth
    % stating here rather than only in the author-decisions section.
    \textbf{And 12.1 m is itself a clamp artefact}: with the hammer released at
    $\mathcal{N}(11, 1.5)$ m\,s$^{-1}$, $v^2/g = 12.3$ m, so the modelled hammer is an
    unpowered point mass to within 2\% --- the coefficient table is edge-clamped at
    $+30^\circ$ and the hammer is released at $+40^\circ$, so it generates no lift at all.
    The hammer
"""

NEW_RESULTS = r"""    % WP4b item 6 RETRACTED [2026-09-04] after a physics referee challenge and direct
    % execution of trajectories.hammer(): alpha stays in [-2.65, +6.41] deg for the whole
    % flight, 0.0000 of samples leave the tabulated band, and CL runs 0.0157..0.5313.
    % The v^2/g comparison was also against the wrong formula -- it ignores the 40 deg
    % elevation and the 2.0 m launch height. True no-aerodynamics range from the same
    % initial condition is 14.208 m, so the disc flies 15% SHORT of ballistic, not equal.
    The short range is a release-geometry effect rather than a table artefact.
    \texttt{trajectories.hammer()} sets a roll of $+80^\circ$, which lays the disc nearly
    on its side; the lift vector is then close to horizontal and does little work against
    gravity. Lift is produced --- $C_L$ reaches 0.53 and the peak lift force, 2.25 N,
    exceeds the disc's weight of 1.71 N --- but it is mis-directed rather than absent, and
    the angle of attack stays within $[-2.7^\circ, +6.4^\circ]$, comfortably inside the
    tabulated band. Sweeping the release roll from $0^\circ$ to $90^\circ$ moves the range
    monotonically from 2.1 m to 12.4 m; only the $0^\circ$ end drives $\alpha$ past
    $+30^\circ$ and out of the table. The hammer
"""

OLD_DECISION = r"""% WP4b item 6 [2026-09-03]: this default is REJECTED. "Exploratory" implies a weak
% result; what the model has for the hammer is no aerodynamic result at all, and the
% still-air range is the evidence. See the amendment.
  \item \textbf{Hammer throw.} \emph{Previous default: retained but reported as
    \textbf{exploratory}, on the grounds that a hammer spends most of its flight well
    outside the $-10^\circ$ to $+30^\circ$ band over which the tabulated Potts \&
    Crowther coefficients exist, so its trajectories rely on edge-clamped values.}

    \textbf{That default is withdrawn, because ``exploratory'' overstates what is
    there.} The model's still-air hammer range is 12.1 m. The hammer release speed is
    $\mathcal{N}(11, 1.5)$ m\,s$^{-1}$, and $v^2/g = 11^2/9.81 = 12.3$ m. \textbf{The
    modelled hammer is, to within 2\%, an unpowered point mass.} That is what
    edge-clamping the coefficient table at $+30^\circ$ produces for a throw released at
    $\alpha_0 = +40^\circ$: no lift, no pitching moment, no aerodynamics --- so 12.1 m is
    a clamp artefact and not a model result, and it is also why the success rate is
    structurally zero (\S\ref{sec:results}) whatever the wind does. A number that
    reproduces $v^2/g$ is not evidence about a disc.

    \textbf{Two admissible options, and ``report it as exploratory'' is not one of
    them.} \emph{(a)} Drop the hammer from the paper entirely; the backhand and forehand
    results do not depend on it. \emph{(b)} Keep it only after extending the coefficient
    table past stall \emph{and} re-parameterising the release --- a real hammer leaves the
    hand at roughly 18--21 m\,s$^{-1}$, not 11 --- and then re-running. \emph{\% FLAG FOR
    AUTHOR: (a) costs nothing and is the default we recommend; (b) is a genuine piece of
    work and would make the hammer the most interesting throw in the study rather than the
    weakest.}"""

NEW_DECISION = r"""% WP4b item 6 was REJECTED here on 2026-09-03 and that rejection is itself RETRACTED
% on 2026-09-04. The "clamp artefact" evidence did not survive execution; see below.
  \item \textbf{Hammer throw.} \emph{Default: retained but reported as
    \textbf{exploratory}.}

    \textbf{A stronger claim was entered here on 2026-09-03 and is now withdrawn.} It
    held that the model's 12.1 m still-air hammer range is a clamp artefact --- that a
    release at $\alpha_0 = +40^\circ$ against a table edge-clamped at $+30^\circ$ leaves
    the disc with ``no lift at all'', and that 12.1 m reproduces $v^2/g = 12.3$ m to
    2\%, making the hammer an unpowered point mass. Direct execution of
    \texttt{trajectories.hammer()} refutes every link. The angle of attack stays in
    $[-2.65^\circ, +6.41^\circ]$ for the entire flight and \emph{no} sample leaves the
    tabulated band, so no edge-clamping occurs; $C_L$ runs 0.0157 to 0.5313 and the peak
    lift force of 2.25 N exceeds the 1.71 N weight, so lift is produced. The ballistic
    comparison was also made against the wrong formula: $v^2/g$ ignores both the
    $40^\circ$ launch elevation and the 2.0 m launch height. Integrating the same initial
    condition with no aerodynamics gives 14.21 m, so the modelled hammer flies 15\%
    \emph{short} of ballistic --- the ``agreement to 2\%'' was a coincidence between a
    model result and a formula that does not describe this launch.

    \textbf{What is actually true is milder, and it restores the original default.} The
    short range is release geometry: \texttt{hammer()} rolls the disc $+80^\circ$, so the
    lift vector is nearly horizontal and does little work against gravity. Lift is
    mis-directed, not absent. Sweeping roll $0^\circ \to 90^\circ$ moves the range
    monotonically 2.1 m $\to$ 12.4 m. Note this also corrects the previous default's
    \emph{stated ground}: the hammer as parameterised does \emph{not} spend most of its
    flight outside the $-10^\circ$ to $+30^\circ$ band --- it never leaves it. It is
    exploratory because the release parameters are unvalidated (a real hammer leaves the
    hand at roughly 18--21 m\,s$^{-1}$, not 11, and the $+80^\circ$ roll is a guess), not
    because the coefficients are extrapolated. \emph{\% FLAG FOR AUTHOR: extending the
    coefficient table past stall is therefore \textbf{not} the fix --- the table is never
    exceeded. If the hammer is to carry weight, the release parameterisation is what needs
    measuring. Dropping it remains admissible and costs nothing, since the backhand and
    forehand results do not depend on it.}"""


def expand_inputs(src: str, base: Path, _depth: int = 0) -> str:
    """Splice in every \\input{...} so label resolution sees the generated files.

    The paper \\input{}s three files produced by make_results_tex.py and friends; the
    Q-by-heading table and its label live in results_generated.tex, not in the
    manuscript. Resolving refs without these gives false positives.
    """
    if _depth > 5:
        return src

    def sub(m: "re.Match[str]") -> str:
        p = (base / m.group(1)).with_suffix(".tex")
        if not p.is_file():
            print("  note: \\input target not found, ignoring: %s" % p.name)
            return ""
        return expand_inputs(p.read_text(encoding="utf-8"), p.parent, _depth + 1)

    return re.sub(r"\\(?:input|include)\{([^}]+)\}", sub, src)


def labels_of(src: str) -> set:
    return set(re.findall(r"\\label\{([^}]+)\}",
                          expand_inputs(src, TEX.parent)))


def refs_of(src: str) -> set:
    return set(re.findall(r"\\(?:ref|eqref)\{([^}]+)\}",
                          expand_inputs(src, TEX.parent)))


def main() -> int:
    src = TEX.read_text(encoding="utf-8")
    orig = src
    for old, new, label in ((OLD_RESULTS, NEW_RESULTS, "results"),
                            (OLD_DECISION, NEW_DECISION, "author-decisions")):
        if old not in src:
            # idempotence: recognise an already-applied patch by a line unique to it
            if "RETRACTED [2026-09-04]" in src and label == "results":
                print("already applied: %s" % label)
                continue
            if "itself RETRACTED" in src and label == "author-decisions":
                print("already applied: %s" % label)
                continue
            print("ANCHOR NOT FOUND (%s) -- refusing to write" % label)
            return 1
        src = src.replace(old, new, 1)
        print("patched: %s" % label)

    if src == orig:
        print("no change")
        return 0

    # Refuse to write if the edit introduced a dangling \ref.
    #
    # This MUST follow \input. The first version scanned only the manuscript file and
    # reported tab:Qbyheading as dangling, which blocked a correct patch: that label is
    # defined in SeanEarl_DiscFlight/P5_Synthesis/results_generated.tex, one of three
    # generated files the paper \input{}s. A guard that cannot see included files
    # rejects every legitimate edit, so it was the guard that was wrong, not the edit.
    dangling = sorted(refs_of(src) - labels_of(src))
    if dangling:
        print("DANGLING REFERENCES, refusing to write:", dangling)
        return 1

    TEX.write_text(src, encoding="utf-8")
    print("wrote %s (%d -> %d chars); %d labels, %d refs, none dangling"
          % (TEX.name, len(orig), len(src),
             len(labels_of(src)), len(refs_of(src))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
