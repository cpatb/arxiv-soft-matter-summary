"""python -m studio <command> <episode.yaml>

    lint      check the spec: fact IDs exist, numbers are sourced, footage exists, chapters
    script    print a teleprompter (markdown) for recording the narration yourself
    full      render the full 16:9 video + SRT + description
    shorts    render every Short (or --only s01)
    thumb     render the thumbnail
    all       lint, full, shorts, thumb

Colorization is separate: python -m studio.colorize --help
"""

from __future__ import annotations

import argparse
import sys

from . import episode as epm
from . import render as rd
from . import timeline as tl
from . import voice as vo


def cmd_lint(ep) -> int:
    dur = rd.film_duration(ep.color) if ep.color.exists() else None
    problems = epm.lint(ep, dur)
    for path in (ep.bw, ep.color):
        if not path.exists():
            problems.append(f"missing master {path}")
    for p in problems:
        print("WARN", p)
    print(f"{len(ep.cues)} full cues, {len(ep.shorts)} shorts, {len(problems)} warning(s)")
    return 0


def cmd_script(ep) -> int:
    rate = 150 / 60     # words per second, a typical documentary read
    print(f"# Narration script: {ep.title}\n")
    print("Record each line as its own file in `recordings/` named by its ID "
          "(e.g. `recordings/open01.wav`). Any line without a recording uses the scratch "
          "TTS voice. Leave ~0.3 s of room tone at each end; the renderer adds the pauses.\n")
    groups = [("Full video", ep.cues)] + [(f"Short {s.id}: {s.title}", s.cues) for s in ep.shorts]
    for name, cues in groups:
        words = sum(len(c.say.split()) for c in cues)
        print(f"## {name}  (~{words / rate / 60:.1f} min at 150 wpm)\n")
        for c in cues:
            if not c.say:
                continue
            have = "recorded" if vo.find_recording(ep.recordings, c.id) else "TODO"
            src = f"  _[{', '.join(c.facts)}]_" if c.facts else ""
            print(f"**{c.id}** ({have}, ~{len(c.say.split()) / rate:.0f}s){src}\n\n> {c.say}\n")
    return 0


def cmd_full(ep, burn: bool) -> int:
    final, placed, probs = rd.render_full(ep, burn_captions=burn)
    print(f"full: {final}  ({tl.total(placed) / 60:.1f} min, {len(placed)} cues)")
    for p in probs:
        print("WARN chapters:", p)
    return 0


def cmd_shorts(ep, only: str | None) -> int:
    for s in ep.shorts:
        if only and s.id != only:
            continue
        final, length = rd.render_short(ep, s)
        flag = "" if length >= 60 else "  (under 60 s: below TikTok's monetization minimum)"
        print(f"short {s.id}: {final}  ({length:.1f}s){flag}")
    return 0


def cmd_thumb(ep) -> int:
    print("thumb:", rd.thumbnail(ep, ep.output / f"{ep.id}_thumb.jpg"))
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="python -m studio", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("command", choices=["lint", "script", "full", "shorts", "thumb", "all"])
    p.add_argument("episode")
    p.add_argument("--only", help="shorts: render just this short id")
    p.add_argument("--burn-captions", action="store_true",
                   help="full: burn captions in (default: separate SRT only)")
    a = p.parse_args(argv)
    ep = epm.load(a.episode)
    if a.command == "lint":
        return cmd_lint(ep)
    if a.command == "script":
        return cmd_script(ep)
    if a.command == "full":
        return cmd_full(ep, a.burn_captions)
    if a.command == "shorts":
        return cmd_shorts(ep, a.only)
    if a.command == "thumb":
        return cmd_thumb(ep)
    cmd_lint(ep)
    cmd_full(ep, a.burn_captions)
    cmd_shorts(ep, a.only)
    return cmd_thumb(ep)


if __name__ == "__main__":
    sys.exit(main())
